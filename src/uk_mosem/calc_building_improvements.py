"""
Script deriving techno-economic parameters for building efficiency measures


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

import pandas as pd
import numpy as np
import camelot as cl

import utils

import logging

logger = logging.getLogger(__name__)
utils.adjust_logger()

if __name__ == "__main__":

    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_building_improvements",
                                         **df_params.iloc[0].to_dict())
        
    # calculate potential for domestic building efficiency measures
    
    # load stock data, including characteristics of wall, roof, etc.
    stock = pd.read_csv(snakemake.input.path_dw_stock_ede)
    stock["TYPE"] = stock["EFF"].str.split("_",expand=True)[0]
    stock = stock.set_index(["LSOA11CD","PROPERTY_TYPE","TYPE","EFF"])
    
    dfl = list()
    for t in stock.index.get_level_values("TYPE").unique():
        st = stock.xs(t,level=2)
        
        # take into account NAs
        # if non-NA exists for LSOA and property type, scale up non NA
        # if only NA exists calc fraction based on average of LA (happening below)
        
        # scale if NA and non-NAs present
        st = st.unstack(fill_value=0).droplevel(0,axis=1)
        
        sts = st.loc[st.sum(axis=1)>st.loc[:,t+"_NA"]]
        sts = sts.divide(1-(sts.loc[:,t+"_NA"]
                         /sts.sum(axis=1)).fillna(0),axis=0)
        sts.loc[:,t+"_NA"] = 0
        st.loc[st.sum(axis=1)>st.loc[:,t+"_NA"]] = sts
        

       
        if st.loc[:,t+"_NA"].sum()>0:
            logger.warning("One or more property types in one or more LSOAs"
                           f" do not have any data to calculate {t}"
                           " characteristics. Using LAD average for property"
                           " type.")
            
            # calculate mean proportions for LADs
            LAD_mean = st.loc[(st[t+"_NA"]==0)
                              &(st.sum(axis=1)>0)].drop(t+"_NA",axis=1)
            LAD_mean = utils.groupby_LAD(LAD_mean).sum()
            LAD_mean = LAD_mean.divide(LAD_mean.sum(axis=1),axis=0)
            
            LAD_mean = LAD_mean.stack()
            LAD_mean.index = LAD_mean.index.rename("EFF",level=2)
            LAD_mean.name = "NUMBER_OF_PROPERTIES"
            
            st = st.stack()
            st.index = st.index.rename("EFF",level=2)
            st.name = "NUMBER_OF_PROPERTIES"
            # merge LAD means on stock data
            lsoas = utils.get_entity_lookup(["LSOA","LAD"]).drop("LAD23NM",axis=1)
            
            st = st.reset_index().merge(right=lsoas,on="LSOA11CD",how="left")
            #hsyss[t] = hsyss[t].loc[hsyss[t]["HEATING_SYSTEM"]!="NA"]
            st = st.merge(right=LAD_mean,
                          on=["LAD23CD","PROPERTY_TYPE","EFF"],
                          how="left",suffixes = ["","_LAD"])
            st = st.drop("LAD23CD",axis=1)
            st = st.set_index(["LSOA11CD","PROPERTY_TYPE","EFF"])
            st = st.unstack()
            
            # replace LSOA values with LAD averages if only NA values exist
            for c in st.columns.get_level_values(1).unique():
                if c == t+"_NA":
                    continue
                st.loc[st[("NUMBER_OF_PROPERTIES",t+"_NA")]>0,
                         ("NUMBER_OF_PROPERTIES",c)] = (st.loc[st[("NUMBER_OF_PROPERTIES",t+"_NA")]>0,
                                                                 ("NUMBER_OF_PROPERTIES_LAD",c)]
                                                              *
                                                              st.loc[st[("NUMBER_OF_PROPERTIES",t+"_NA")]>0,
                                                                     ("NUMBER_OF_PROPERTIES",t+"_NA")])
            st = st.drop("NUMBER_OF_PROPERTIES_LAD",axis=1)
            st = st.drop(t+"_NA",axis=1,level=1)
            st = st.stack(dropna=False)
            
            if st["NUMBER_OF_PROPERTIES"].isna().any():
                logger.warning("Was not able to use LAD average in all instances.")
                
            dfl.append(st)

    stock = pd.concat(dfl)
    
    stock = stock.reset_index()
    # dropping uninsulated solid floors in flats following assumption 
    # in data on efficiency measures no such insulation is available
    stock = stock.loc[~((stock["PROPERTY_TYPE"].str.startswith("Flats"))&
                        (stock["EFF"]=="floor_solid_uninsulated"))]
    #stock = stock.set_index(["LSOA11CD","PROPERTY_TYPE","EFF"])
    #stock = stock.divide(stock.groupby("EFF").sum())
    
    # rename property characteristics using connected efficiency measures, while
    # adding separate entries for different efficiency measures linked to same
    # property characteristic (assuming same distribution across LSOAs) and
    # dropping irrelevant ones
    #stock = stock.reset_index()
    stock = pd.concat([stock,stock[stock["EFF"]=="wall_solid_uninsulated"].assign(EFF="wall_solid_insulation_ex"),
                      stock[stock["EFF"]=="wall_cavity_uninsulated"].assign(EFF="wall_cavity_insulation_eh"),
                      stock[stock["EFF"]=="roof_pitched_partially_insulated"].assign(EFF="roof_insulation_add_eh"),
                      stock[stock["EFF"]=="wall_cavity_uninsulated"].assign(EFF="wall_cavity_insulation_h"),
                      stock[stock["EFF"]=="roof_pitched_partially_insulated"].assign(EFF="roof_insulation_add_h"),
                      stock[stock["EFF"]=="roof_pitched_uninsulated"].assign(EFF="roof_insulation_eh"),
                      stock[stock["EFF"]=="roof_pitched_uninsulated"].assign(EFF="roof_insulation_h")])                                                                       
   

    stock = stock.replace({"wall_solid_uninsulated":"wall_solid_insulation_in",
                           "wall_cavity_uninsulated":"wall_cavity_insulation_e",
                           "roof_pitched_partially_insulated":"roof_insulation_add_e",
                           "roof_pitched_uninsulated":"roof_insulation_e",
                           "windows_double_glazed_old":np.nan,
                           "windows_double_glazed_new":np.nan,
                           "roof_pitched_fully_insulated":np.nan,
                           "windows_single_glazed":"windows_double_from_G",
                           "windows_secondary_glazing":np.nan,
                           "windows_triple_glazed":np.nan,
                           "wall_cavity_insulated":np.nan,
                           "wall_solid_insulated":np.nan,
                           "floor_solid_uninsulated":"floor_solid_insulation",
                           "floor_suspended_uninsulated":"floor_suspended_insulation",
                           "floor_solid_insulated":np.nan,
                           "floor_suspended_insulated":np.nan,
                           "floor_no_floor":np.nan,
                           "roof_no_roof":np.nan,
                           "roof_flat_uninsulated":np.nan,
                           "roof_flat_insulated":np.nan,
                           }).dropna()
    
    stock = stock.set_index(["LSOA11CD","PROPERTY_TYPE","EFF"])
    
    
    # load national (GB) building characteristics/potential for measures
    # for main, and then other measures, clean data
    meas = cl.read_pdf(snakemake.input.path_dw_building_num,
                        pages="247",process_background=False)[0].df
    meas = meas.replace(r'^\s*$', np.nan, regex=True)
    meas = meas.dropna(axis=1,thresh=6).dropna(axis=0,thresh=2)
    meas.iloc[1,2:] = meas.iloc[0,2:]
    meas = meas.dropna(axis=0,thresh=8)
    for i in range(1,len(meas.columns),2):
        meas.iloc[:,i] = meas.iloc[:,i].str.replace(',', '').astype(int)
        
    meas.iloc[:,1] = meas.iloc[:,1]+meas.iloc[:,3]+meas.iloc[:,5]
    frac_NI = meas.iloc[:,7].sum()/(meas.iloc[:,1].sum()+meas.iloc[:,7].sum())
    meas = meas.iloc[:,0:2]
    meas.columns = ["measures","stock"]
    meas = meas.set_index("measures")["stock"]
    
    # load and clean list of potential for other measures (given as percentage
    # of stock for different size categories)
    omeas = pd.read_excel(snakemake.input.path_building_eff,
                          sheet_name="Energy efficiency",
                          skiprows=51,header=1,usecols=[2,3,4],index_col=[0,1])
    
    # calculating the average across sizes with no weighting, adding fraction
    # of single glazed windows (assumed equal to number of non-double glazed),
    # multiply with number of properties derived from data on main measures
    omeas = omeas.groupby("Measure").mean()
    omeas.index.name="measures"
    omeas.columns = ["percentage"]
    
    # load raw data set for non-double glazed fraction
    raw_data = pd.read_csv(snakemake.input.path_set_bd_eff)
    
    sgfrac = raw_data[["TECHNOLOGY","UNIT",
                         "VALUE"]][raw_data["VARIABLE"]=="NonDoubleGlazingFrac"]
    omeas.loc["windows_double_from_G","percentage"] =  sgfrac["VALUE"].values
    
    omeas["stock"] = omeas["percentage"]*meas.iloc[0:7].sum()
    meas = pd.concat([meas,omeas["stock"]])
    
    # rename measures, drop characteristics not associated with measures
    meas = meas.reset_index()
    meas = meas.replace({"Cavity uninsulated extra-hard-to-treat": "wall_cavity_insulation_eh",
                         "Cavity uninsulated hard-to-treat": "wall_cavity_insulation_h",
                         "Cavity uninsulated easy-to-treat": "wall_cavity_insulation_e",
                         "Cavity insulated": np.nan,
                         "Solid uninsulated (internal)": "wall_solid_insulation_in",
                         "Solid uninsulated (external)": "wall_solid_insulation_ex",
                         "Solid insulated": np.nan,
                         "No roof": np.nan,
                         "Less than 100mm extra-hard-to-treat": "roof_insulation_eh",
                         "Less than 100mm hard-to-treat": "roof_insulation_h",
                         "Less than 100mm easy-to-treat": "roof_insulation_e",
                         "100-199mm extra-hard-to-treat": "roof_insulation_add_eh",
                         "100-199mm hard-to-treat": "roof_insulation_add_h",
                         "100-199mm easy-to-treat": "roof_insulation_add_e",
                         "200mm or more": np.nan,
                         "Solid uninsulated": "floor_solid_insulation",
                         "Suspended uninsulated": "floor_suspended_insulation",
                         "Suspended insulated": np.nan,
                         "None": np.nan,
                         "Draught proofing (draught stripping)": "other_draught_proofing",
                         "HW Tank insulation": "other_hw_tank_insulation",
                         "Insulated doors": "other_door_insulation",
                         "Reduced infiltration (foam, strips, sealent use)":"other_reduced_infiltration",
                         "Solid floor insulation": np.nan,
                         "Suspended timber floor insulation": np.nan
                         }).dropna()
    
    meas = meas.rename(columns={"measures":"EFF"})
    meas = meas.set_index("EFF")
    
    # calculate number of each measures available in each LSOA, including
    # heritage constraints
    
    # loading fraction of properties  in conservation area per LSOA      
    careas = pd.read_csv(snakemake.input.path_conserv_areas,
                              index_col="LSOA11CD")
    # apply data on suitability in conservation areas
    hsf = raw_data.loc[raw_data["VARIABLE"]=="HeritageSuitableFraction",
                       ["TECHNOLOGY","VALUE"]]
    hsf= hsf.assign(TECHNOLOGY=hsf['TECHNOLOGY'].str.split(',')).explode('TECHNOLOGY')
    hsf = hsf.rename(columns={"TECHNOLOGY":"EFF"})
    hsf = hsf.set_index("EFF")

    constr = stock.multiply(1-hsf["VALUE"],level=2,axis=0)
    constr = constr.multiply(careas["cafrac"],level=0,axis=0)
    constr = constr.fillna(0)
    
    stock = stock - constr
    
    # divide by sum of available measures to get fraction of overall available 
    # in each property type in each LSOA
    stock = stock.divide(stock.groupby("EFF").sum())
    
    stock = stock["NUMBER_OF_PROPERTIES"] * meas["stock"]
    
    # load cost data for measures
    cost = pd.read_excel(snakemake.input.path_building_eff,
                      sheet_name="Measures and costs",
                      skiprows=5,nrows=19,header=[0,1])
    cost = cost.iloc[:,3:28].set_index("Unnamed: 3_level_0")
    cost.index = [i[0].strip(" *") for i in cost.index]
    cost = cost.apply(pd.to_numeric,errors="coerce")
     
    costsupp = pd.read_excel(snakemake.input.path_building_eff,
                            sheet_name="Measures and costs",
                            skiprows=33,nrows=6,header=0,usecols="D:G")
    costsupp = costsupp[["Cost Type","Average"]]
    costsupp = costsupp.replace({"Construction Overheads - Flat LR": "Flat LR",
                         "Construction Overheads - Flat HR": "Flat HR",
                         "Construction Overheads - Terrace": "Terraced",
                         "Construction Overheads - Semi Detached": "Semi-detached",
                         "Construction Overheads - Detached": "Detached",
                         "Design & Plannning": "Design & Planning"})
    costsupp = costsupp.set_index("Cost Type")
    costsupp.index.name = "EFF"
    costsupp.loc["Flats"] = costsupp.loc["Flat LR":"Flat HR"].mean()
    costsupp = costsupp.drop(["Flat LR","Flat HR"])

    # calculate cost for appropriate building types
    
    # take the cost for medium size properties
    cost = cost.xs("M",level=1, axis=1)
    
    cost["Flats"] = cost["Flat LR"]
    cost["Terraced"] = cost[["Mid Terrace","End Terrace"]].mean(axis=1)
    cost["Semi-detached"] = cost["Semi Detached"]
    cost = cost[["Flats","Semi-detached","Detached","Terraced"]]
    cost = cost.reset_index()
    
    # rename measures, drop measures not incorporated
    cost["index"] = cost["index"].replace({
            "Hard to treat cavities (Unfilled)": "wall_cavity_insulation_h|"
                                                 "wall_cavity_insulation_eh",
            "Easy to treat cavities (Unfilled)": "wall_cavity_insulation_e",
            "Solid wall insulation (Internal)": "wall_solid_insulation_in",
            "Solid wall insulation (External)": "wall_solid_insulation_ex",
            "Loft insulation (HTT)": "roof_insulation_h|"
                                     "roof_insulation_add_h|"
                                     "roof_insulation_eh|"
                                     "roof_insulation_add_eh",
            "Loft insulation": "roof_insulation_e|"
                               "roof_insulation_add_e",
            "Solid floor insulation": "floor_solid_insulation",
            "Suspended timber floor insulation": "floor_suspended_insulation",
            "Draught proofing (draught stripping)": "other_draught_proofing",
            "HW Tank insulation": "other_hw_tank_insulation",
            "Insulated doors": "other_door_insulation",
            "Reduced infiltration (foam, strips, sealent use)":"other_reduced_infiltration",
            "Solid Wall insulation (Thin)":np.nan,
            "Partially filled cavities":np.nan,
            "Secondary glazing": np.nan,
            "Double glazing": "windows_double_from_G",
            "Double Glazing- Slim profile": np.nan,
            "Triple Glazing": np.nan,
            "Ventilation": np.nan
            }).str.split("|")
    cost = cost.dropna(subset=["index"]).explode("index").set_index("index")
    cost.index.name = "EFF"

    # add supplementary cost  
    cost.loc["wall_solid_insulation_ex",:] = (cost.loc["wall_solid_insulation_ex",:]
                                              + costsupp["Average"]
                                              + costsupp.loc["Design & Planning","Average"])
    cost.loc["wall_solid_insulation_in",:] = (cost.loc["wall_solid_insulation_in",:]
                                              + costsupp.loc["Design & Planning","Average"])    
    
    # rearrange cost and adjust base year
    cost = cost.stack()
    cost.index.names = ["EFF","PROPERTY_TYPE"]
    cost = cost.reorder_levels(["PROPERTY_TYPE","EFF"]).to_frame()
    cost["UNIT"] = "2019£"
    cost.columns = ["VALUE","UNIT"]
    cost = utils.adjust_monetary_values(cost, 2015)
    cost = cost.drop("UNIT",axis=1)
    
    if snakemake.params.dic["scen_hh_disagg"]=="T":
        cost = utils.explode_hh(cost, ["OO","RP","RS"])
        
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        cost = utils.explode_hh(cost,
                                  [t+str(i) for t in ["O","P","S"]
                                   for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)]
                                  )
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )    
    # load energy savings data for measures
    savi = pd.read_excel(snakemake.input.path_building_eff,
                      sheet_name="Energy savings",
                      skiprows=35,nrows=23,header=[0,1])
    savi = savi.iloc[:,3:28].set_index("Unnamed: 3_level_0")
    savi.index = [i[0].strip(" *") for i in savi.index]    
    savi = savi.apply(pd.to_numeric,errors="coerce")    

    # calculate savings for appropriate building types
    
    # take the savings for medium size properties
    savi = savi.xs("M",level=1, axis=1)
    
    savi["Flats"] = savi[["Flat CNV",
                         "Flat LR",
                         "Flat HR"]].mean(axis=1)
    savi["Terraced"] = savi[["Mid Terrace","End Terrace"]].mean(axis=1)
    savi["Semi-detached"] = savi["Semi Detached"]
    
    savi = savi[["Flats","Semi-detached","Detached","Terraced"]]
    savi = savi.reset_index()
    
    # rename measures
    savi["index"] = savi["index"].replace({
            "Hard to treat cavities (Unfilled)": "wall_cavity_insulation_h|"
                                                 "wall_cavity_insulation_eh",
            "Easy to treat cavities (Unfilled)": "wall_cavity_insulation_e",
            "Internal (solid) wall insulation": "wall_solid_insulation_in",
            "External (solid) wall insulation": "wall_solid_insulation_ex",
            "Loft insulation (HTT) \n(see loft potential factor in below section)": "roof_insulation_h|"
                                     "roof_insulation_add_h|"
                                     "roof_insulation_eh|"
                                     "roof_insulation_add_eh",
            "Loft insulation\n(see loft potential factor in below section)": "roof_insulation_e|"
                               "roof_insulation_add_e",
            "Solid floor insulation": "floor_solid_insulation",
            "Suspended timber floor insulation": "floor_suspended_insulation",
            "Draught proofing (draught stripping)": "other_draught_proofing",
            "HW Tank insulation": "other_hw_tank_insulation",
            "Insulated doors": "other_door_insulation",
            "Reduced infiltration (foam, strips, sealent use)":"other_reduced_infiltration",
            "Thin internal (solid) wall insulation":np.nan,
            "Partially filled cavities":np.nan,
            "Secondary glazing (From Band G)": np.nan,
            "Secondary glazing (From Band E)": np.nan,
            "Double glazing (From Band G)": "windows_double_from_G",
            "Double glazing (From Band E)": np.nan,
            "Double Glazing- Slim profile (From Band G)": np.nan,
            "Double Glazing- Slim profile (From Band E)": np.nan,
            "Triple Glazing (From Band G)": np.nan,
            "Triple Glazing (From Band E)": np.nan,
            "Ventilation": np.nan
            }).str.split("|")
    savi = savi.dropna(subset=["index"]).explode("index").set_index("index")
    savi.index.name = "EFF"

    # load loft potential factor
    lpfac = pd.read_excel(snakemake.input.path_building_eff,
                      sheet_name="Energy savings",
                      skiprows=66,nrows=6,header=[0])
    lpfac = lpfac .iloc[:,2:4].set_index("Disaggregated Loft Attribute Type (based on level of existing insulation)")
    lpfac.index.name = "Type"   
    lpfac  = lpfac .apply(pd.to_numeric,errors="coerce")  
    
    # apply loft potential factor on savings
    savi.loc[["roof_insulation_e",
              "roof_insulation_h",
              "roof_insulation_eh"]] = (savi.loc[["roof_insulation_e",
                                                 "roof_insulation_h",
                                                "roof_insulation_eh"]]
                                        * lpfac.loc["Less than 100mm ETT","Loft Potential Factor"])
    savi.loc[["roof_insulation_add_e",
              "roof_insulation_add_h",
              "roof_insulation_add_eh"]] = (savi.loc[["roof_insulation_add_e",
                                                 "roof_insulation_add_h",
                                                "roof_insulation_add_eh"]]
                                        * lpfac.loc["100 -199mm ETT","Loft Potential Factor"])                                             
                                                      
    # load in use factor and closure data
    closure = raw_data[["TECHNOLOGY","UNIT",
                         "VALUE"]][(raw_data["VARIABLE"]=="InUseFactorClosure")]
    
    iuf = pd.read_excel(snakemake.input.path_building_eff,
                      sheet_name="Energy efficiency",
                      skiprows=25,nrows=8,header=[0])
    iuf = iuf.iloc[:,2:4].set_index("Measure")
    iuf.columns = ["IUF"]
    iuf.index.name = "EFF"   
    iuf  = iuf.apply(pd.to_numeric,errors="coerce") 
    
    iuf = iuf.reset_index()
    iuf["EFF"] = iuf["EFF"].replace({
            "Cavity wall": "wall_cavity_insulation_h|"
                            "wall_cavity_insulation_eh|"
                            "wall_cavity_insulation_e",
            "Solid wall": "wall_solid_insulation_in|"
                            "wall_solid_insulation_ex",
            "Roof": "roof_insulation_h|"
                                     "roof_insulation_add_h|"
                                     "roof_insulation_eh|"
                                     "roof_insulation_add_eh|"
                                     "roof_insulation_e|"
                                     "roof_insulation_add_e",
            "Floor": "floor_solid_insulation|"
                    "floor_suspended_insulation",
            "Other": "other_draught_proofing|"
                    "other_hw_tank_insulation|"
                    "other_door_insulation|"
                    "other_reduced_infiltration",
            "Glazing": "windows_double_from_G",
            "Behavioural": np.nan,
            }).str.split("|")  
    iuf = iuf.dropna(subset=["EFF"]).explode("EFF").set_index("EFF")
    iuf["ClosureFactor"] = ((1-iuf["IUF"])+iuf["IUF"]*closure["VALUE"].values)/(1-iuf["IUF"])

    # apply closure factor to savings data
    savi = savi.multiply(iuf["ClosureFactor"], axis=0)
    savi = savi.stack()
    savi.index.names = ["EFF","PROPERTY_TYPE"]
    savi = savi.reorder_levels(["PROPERTY_TYPE","EFF"])

    if snakemake.params.dic["scen_hh_disagg"]=="T":
        savi = utils.explode_hh(savi.to_frame(), ["OO","RP","RS"]).squeeze(axis=1)
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        savi = utils.explode_hh(savi.to_frame(),
                                  [t+str(i) for t in ["O","P","S"]
                                   for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)]
                                  ).squeeze(axis=1)
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )   
    # load and process deployment constraints
    dc = pd.read_excel(snakemake.input.path_building_eff,
                      sheet_name="Scenario assumptions",
                      skiprows=42,nrows=3,header=[0],usecols=[3,4,5,6])    
    dc.columns = ["EFF","2025","2030","2035"]
    dc["EFF"] = dc["EFF"].replace({
            "Energy efficiency – cavity walls": "wall_cavity_insulation_h|"
                            "wall_cavity_insulation_eh|"
                            "wall_cavity_insulation_e",
            "Energy efficiency – solid walls": "wall_solid_insulation_in|"
                            "wall_solid_insulation_ex",
            "Energy efficiency – lofts": "roof_insulation_h|"
                                     "roof_insulation_add_h|"
                                     "roof_insulation_eh|"
                                     "roof_insulation_add_eh|"
                                     "roof_insulation_e|"
                                     "roof_insulation_add_e"
  
            }).str.split("|")  
    dc = dc.dropna(subset=["EFF"]).explode("EFF").set_index("EFF")
    dc = dc.stack()
    
    # apply factor to reduce for fraction of stock in NI
    dc = dc*(1-frac_NI)
    
    m = meas.copy()
    m["group"]=m.index.str.split("_",2).str[0:2].str.join("_")
    m = m.reset_index().merge(m.groupby("group").sum(), on="group", how="left",
                                           suffixes=["_",""])
    m = m.set_index("EFF")["stock"]
    
    dc = dc.to_frame().reset_index().merge(m, on=["EFF"],how="left")
    dc["rate"] = dc[0]/dc["stock"]
    dc = dc.set_index(["EFF","level_1"])["rate"]    
    dc = dc.unstack(1)
    
    # load lifetimes for measures
    lt = pd.read_excel(snakemake.input.path_building_eff,
                      sheet_name="Measures and costs",
                      skiprows=70,nrows=19,header=[0],usecols=[3,4])
    lt.columns =["EFF","lifetime"]
    lt = lt.set_index("EFF")
    lt = lt.rename(index=lambda x: x.strip(" *"))
    lt = lt.apply(pd.to_numeric,errors="coerce")
    
    lt = lt.reset_index()
    
    # rename measures
    lt["EFF"] = lt["EFF"].replace({
            "Hard to treat cavities (Unfilled)": "wall_cavity_insulation_h|"
                                                 "wall_cavity_insulation_eh",
            "Easy to treat cavities (Unfilled)": "wall_cavity_insulation_e",
            "Solid wall insulation (Internal)": "wall_solid_insulation_in",
            "Solid wall insulation (External)": "wall_solid_insulation_ex",
            "Loft insulation (HTT)": "roof_insulation_h|"
                                     "roof_insulation_add_h|"
                                     "roof_insulation_eh|"
                                     "roof_insulation_add_eh",
            "Loft insulation": "roof_insulation_e|"
                               "roof_insulation_add_e",
            "Solid floor insulation": "floor_solid_insulation",
            "Suspended timber floor insulation": "floor_suspended_insulation",
            "Draught proofing (draught stripping)": "other_draught_proofing",
            "HW Tank insulation": "other_hw_tank_insulation",
            "Insulated doors": "other_door_insulation",
            "Reduced infiltration (foam, strips, sealent use)":"other_reduced_infiltration",
            "Solid Wall insulation (Thin)":np.nan,
            "Partially filled cavities":np.nan,
            "Secondary glazing": np.nan,
            "Double glazing": "windows_double_from_G",
            "Double Glazing- Slim profile": np.nan,
            "Triple Glazing": np.nan,
            "Ventilation": np.nan
            }).str.split("|")
    lt = lt.dropna(subset=["EFF"]).explode("EFF").set_index("EFF")
  
    

    # load number of properties for each LSOA and type as well as peak demand
    # in order to calculate demand covered by measures
    dw_stock_ex =  pd.read_csv(snakemake.input.path_dw_stock_ex,
                               index_col=["LSOA11CD","PROPERTY_TYPE"],
                               usecols=["LSOA11CD","PROPERTY_TYPE","2019"])
    
    # aggregate, concat number of properties, cap measures to number
    # of properties and reallocate
    slagg = pd.read_csv(snakemake.input.path_sublocal_agg,
                            index_col=["LSOA11CD"]).squeeze()
    
    dw_stock_ex = utils.groupby_LAD(dw_stock_ex,disagg=slagg).sum()
    stock = utils.groupby_LAD(stock,disagg=slagg).sum()
    stock = stock.reset_index().merge(right=dw_stock_ex["2019"],
                                    on=[slagg.name,"PROPERTY_TYPE"],
                                    how="left")
    stock = stock.set_index([slagg.name,"PROPERTY_TYPE","EFF"])
    stock.columns = ["measure","properties"]
    # distribute measure across others if they exceed the number of buildings
    capped = stock[stock["measure"]>stock["properties"]]
    dist = stock[stock["measure"]<stock["properties"]]

    dist["measure"] = (dist["measure"]
                       +((dist["properties"]-dist["measure"])
                       /(dist["properties"]-dist["measure"]).groupby("EFF").sum()
                       *(capped["measure"]-capped["properties"]).groupby("EFF").sum()).fillna(0))
    stock[stock["measure"]<stock["properties"]] = dist
  
    stock["measure_adj"] = stock[["measure","properties"]].min(axis=1)
    
    stock["savings_frac"] = stock["measure_adj"].multiply(savi,
                                                          axis=0).reorder_levels([slagg.name,"PROPERTY_TYPE","EFF"])/stock["properties"]
    
    peakd = pd.read_csv(snakemake.input.path_peakc_lsoa,
                        index_col=["LSOA11CD","PROPERTY_TYPE","TECHNOLOGY"],
                        usecols=["LSOA11CD","PROPERTY_TYPE",
                                 "TECHNOLOGY","2019"])
    # assume low cost eff pack, this is basically irrelevant to substantive
    # results
    peakd = peakd.xs("BELO",level=2)
    peakd = utils.groupby_LAD(peakd,disagg=slagg).sum()
    
    peakd.columns=["peak_demand"]
    stock = stock.reset_index().merge(right=peakd["peak_demand"],
                                    on=[slagg.name,"PROPERTY_TYPE"],
                                    how="left")
    stock = stock.set_index([slagg.name,"PROPERTY_TYPE","EFF"])
    stock["pot_peak_demand_red"] = stock["peak_demand"]*stock["savings_frac"]
    
    # calculate aggregate cost per LAD, property and measure type
    # unit: million [Baseyear]£/GW (= [Baseyear]£/kW)
    stock["agg_cost"] =  stock["measure_adj"].multiply(cost["VALUE"],
                                               axis=0).reorder_levels([slagg.name,
                                                                       "PROPERTY_TYPE",
                                                                       "EFF"])/1000000    
    
    # merge lifetime
    stock = stock.reset_index().merge(right=lt,
                                    on=["EFF"],
                                    how="left")
    stock = stock.set_index([slagg.name,"PROPERTY_TYPE","EFF"])    
  
    # merge and process deployment constraints
    
    stock = stock.reset_index().merge(right=dc,
                                    on=["EFF"],
                                    how="left")
    stock = stock.set_index([slagg.name,"PROPERTY_TYPE","EFF"]) 
    # apply mean to all measures not covered by data
    dcm = stock[["2025","2030","2035"]].mean()
    stock[["2025","2030","2035"]] = stock[["2025",
                                           "2030",
                                           "2035"]].fillna(dcm)
    stock[["2025","2030","2035"]] = stock[["2025","2030","2035"]].multiply(
                                     stock["pot_peak_demand_red"],axis=0)
    
    # aggregate measures to efficiency improvement packages

    # load and process package formulation
    pf = pd.read_excel(snakemake.input.path_building_eff,
                      sheet_name="Energy efficiency",
                      skiprows=4,nrows=16,header=[0],usecols=[3,6,7,8])    
    pf = pf.set_index("Measure")
    pf.columns = ["BELO","BEME","BEHI"]
    pf.index.name = "EFF"
    pf = pf.replace({"*":"1"})
    pf = pf.apply(pd.to_numeric,errors="coerce")
    pf = pf.reset_index()
    pf = pf.drop_duplicates()
    pf["EFF"] = pf["EFF"].replace({
            "Hard to treat cavities (Unfilled)": "wall_cavity_insulation_h|"
                                                 "wall_cavity_insulation_eh",
            "Easy to treat cavities (Unfilled)": "wall_cavity_insulation_e",
            "Internal (solid) wall insulation": "wall_solid_insulation_in",
            "External (solid) wall insulation": "wall_solid_insulation_ex",
            "Loft insulation (HTT)": "roof_insulation_h|"
                                     "roof_insulation_add_h|"
                                     "roof_insulation_eh|"
                                     "roof_insulation_add_eh",
            "Loft insulation": "roof_insulation_e|"
                               "roof_insulation_add_e",
            "Solid floor insulation": "floor_solid_insulation",
            "Suspended timber floor insulation": "floor_suspended_insulation",
            "Draught proofing (draught stripping)": "other_draught_proofing",
            "HW Tank insulation": "other_hw_tank_insulation",
            "Insulated doors": "other_door_insulation",
            "Thin internal (solid) wall insulation":np.nan,
            "Double glazing (From Band G)": "windows_double_from_G",
            "Double glazing (From Band E)": np.nan
            }).str.split("|")
    pf = pf.dropna(subset=["EFF"]).explode("EFF").set_index("EFF").fillna(0)
    pf["BELO"] = pf["BEHI"]+pf["BEME"]+pf["BELO"]
    pf = pf["BELO"].replace({0:np.nan,
                             1:"BEHI",
                             2:"BEME",
                             3:"BELO"})
    # add windows manually to high package
    pf["windows_double_from_G"] = "BEHI"

    # aggregate lifetimes (using weighted average based on agg_cost), cost, 
    # savings based on package formulation
    
    # remove rows where number of measures is 0
    stock = stock[stock["measure_adj"]!=0]

    wm = lambda x: np.average(x, weights=stock.loc[x.index, "agg_cost"])

    stock = stock.reset_index()
    stock["EFF"] = stock["EFF"].replace(pf.to_dict())
    stock = stock.groupby([slagg.name,"PROPERTY_TYPE","EFF"]).agg(agg_cost=("agg_cost", "sum"),
                                                         savings_frac=("savings_frac", "sum"),  
                                                         pot_peak_demand_red=("pot_peak_demand_red", "sum"),  
                                                         lifetime=("lifetime", wm),
                                                         _2025=("2025", "sum"),
                                                         _2030=("2030", "sum"),
                                                         _2035=("2035", "sum")) 
    stock.columns = stock.columns.str.strip("_")
    
    # calculate cost per capacity reduction                                                                
    # unit: million [Baseyear]£/GW (= [Baseyear]£/kW)                                                                 
    stock["cost"] = stock["agg_cost"]/stock["pot_peak_demand_red"]
    
    
        
    ### calculate non-domestic efficiency measure cost and savings
    
    # load non-domestic stock data
    ndstock = pd.read_csv(snakemake.input.path_nd_stock_ex,
                          #index_col=["LSOA11CD","PROPERTY_TYPE"],
                          usecols=["LSOA11CD","PROPERTY_TYPE","2019"])
    ndstock["PROPERTY_TYPE"]="Non-domestic"
    ndstock = ndstock.groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
    
    # load cost and potential data
    costsav = pd.read_excel(snakemake.input.path_nd_building_eff,
                      sheet_name="Table 4.4",
                      skiprows=6,nrows=19,header=0)
    costsav = costsav.dropna(thresh=4).dropna(axis=1,thresh=4)
    costsav = costsav.rename(columns={"Unnamed: 1": "Sector"})
    costsav = costsav.set_index("Sector")
    
    # select space heat measure categories and sum up
    costsav = costsav.loc[["Building fabric"],:].sum(axis=0)
    # costsav = costsav.loc[["Building fabric",
    #                         "Building instrumentation and control"],:].sum(axis=0)
    
    # adjust base year for cost
    c = costsav.copy().to_frame()
    c["UNIT"] = "2015£"
    c=c.loc["Total capital cost of measure (£ thousands)",:].to_frame().transpose()
    c.columns = ["VALUE","UNIT"]
    c = utils.adjust_monetary_values(c, 2015)
    costsav["Total capital cost of measure (£ thousands)"] = c["VALUE"]

    # derive scale up factor based on floor area to incorporate Scotland
    sfac = 1/(1-ndstock[ndstock.index.get_level_values(0).str.startswith("S")].sum()
            /ndstock.sum())[0]
    
    # get uplift factor to take into account different base year, etc. see 6CB
    ndeffup = 1+raw_data.loc[raw_data["VARIABLE"]=="NDEfficiencyPotentialUplift","VALUE"].values[0]                         
        
    # allocate cost across LSOAs based on floor area, while applying scale up
    # factor to adjust for inclusion of Scotland, and heritage exclusions
    
    constr = ndstock.multiply(1-hsf.loc["non-domestic","VALUE"],axis=0)
    constr = constr.multiply(careas["cafrac"],level=0,axis=0)
    constr = constr.fillna(0)
    
    sfrac = ((utils.groupby_LAD(ndstock,disagg=slagg).sum()
              -utils.groupby_LAD(constr,disagg=slagg).sum())
             /utils.groupby_LAD(ndstock,disagg=slagg).sum()).fillna(0)
    ndstock = (ndstock-constr)
    
    ndstock = ndstock/ndstock.sum()
    ndstock.columns = ["areafrac"]
    
    ndstock["agg_cost"] = (ndstock["areafrac"]
                           *costsav["Total capital cost of measure (£ thousands)"]
                           /1000
                           *sfac*ndeffup)
    
    # load overall energy demand for space heat

    rescap = pd.read_csv(snakemake.input.path_res_caps_ht,
                        index_col=["LSOA11CD","PROPERTY_TYPE","TECHNOLOGY"])
    rescap = rescap[~rescap.index.get_level_values(2).str.startswith(("RAUP",
                                                                      "WDIS"))]
    rescap = rescap.reset_index().merge(right=rescap.groupby(["LSOA11CD",
                                                              "PROPERTY_TYPE"]).sum()["2020"],
                          on=["LSOA11CD", "PROPERTY_TYPE"], how="left",
                          suffixes=["","_total"]).set_index(["LSOA11CD",
                                                             "PROPERTY_TYPE",
                                                             "TECHNOLOGY"])                    
    rescap["2020_frac"] = rescap["2020"]/rescap["2020_total"]
    rescap["2020_frac"] = rescap["2020_frac"].fillna(0)
    
    
    anndem_SH = pd.read_csv(snakemake.input.path_nd_ann_dem_sh,
                        usecols=["LSOA11CD","PROPERTY_TYPE","2020"])
    anndem_SH["PROPERTY_TYPE"]="Non-domestic"
    anndem_SH = anndem_SH.groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
    
    techs = [t for t in rescap.index.get_level_values(2).unique()
             if ("DNDO" in t)]
    
    # load eff data set 
    eff_data = pd.read_csv(snakemake.input.path_set_eff)
    
    
    # calculate final energy demand for space heating in GWh
    anndem_SH_fin = anndem_SH.multiply(
                     pd.concat([(rescap.loc[rescap.index.get_level_values(2)==tech,
                                  "2020_frac"].droplevel(2)
                       /eff_data.loc[(eff_data["TECHNOLOGY"]
                                      ==tech)&(
                                      eff_data["MODE_OF_OPERATION"]
                                                     ==1), "VALUE"].values[0]).fillna(0)
                         for tech in techs]).groupby(["LSOA11CD",
                                                      "PROPERTY_TYPE"]).sum()
                                                      ,axis=0).sum()/3600
                                 

    # allocate savings, merge on peak demand before diversity and calculate 
    # potential reduction
    ndstock = utils.groupby_LAD(ndstock,disagg=slagg).sum()
    ndstock["savings_frac"] = (costsav["Total annual energy savings (GWh)"]
                               *sfac*ndeffup/anndem_SH_fin.values[0])*sfrac
    # if demand zero, this creates inf, replace with 0
    ndstock["savings_frac"] = ndstock["savings_frac"].replace([np.inf, -np.inf], 0)
    ndstock = ndstock.reset_index().merge(right=peakd["peak_demand"],
                                    on=[slagg.name,"PROPERTY_TYPE"],
                                    how="left")
    ndstock["EFF"] = "BEST"
    ndstock = ndstock.set_index([slagg.name,"PROPERTY_TYPE","EFF"])
    
    ndstock["pot_peak_demand_red"] = (ndstock["peak_demand"]
                                      *ndstock["savings_frac"])

    # calculate cost per capacity reduction, base year already 2015                                                                
    # unit: million [Baseyear]£/GW (= [Baseyear]£/kW)                                                                 
    ndstock["cost"] = ndstock["agg_cost"]/ndstock["pot_peak_demand_red"]
    # if demand zero, this creates inf, replace with 0
    ndstock["cost"] = ndstock["cost"].replace([np.inf, -np.inf], 0)
    # avoid nan in cost where there is no efficiency potential
    ndstock["cost"] = ndstock["cost"].fillna(0)
    # load lifetime and add
    life_data = pd.read_csv(snakemake.input.path_set_lt)
    
    ndstock["lifetime"] = life_data.loc[(life_data["VARIABLE"]=="OperationalLife")
                                       &(life_data["TECHNOLOGY"]=="BESTDNDO00"),
                                         "VALUE"].values[0]
    
    # set deployment constraint, assuming average from domestic calculation
    ndstock[["2025","2030","2035"]] = dcm
    ndstock[["2025","2030","2035"]] = ndstock[["2025","2030","2035"]].multiply(
                                         ndstock["pot_peak_demand_red"],axis=0)
    
    dndstock = pd.concat([stock,ndstock]).reset_index()
    
    # rename property types, taking into account household disaggregation
    if snakemake.params.dic["scen_hh_disagg"].startswith("T"):
        dndstock["PROPERTY_TYPE"] = ["DDDE"+x.split("|")[1]
                                     if x.startswith("Detached") else
                                     "DDSD"+x.split("|")[1]
                                     if x.startswith("Semi-detached") else
                                     "DDTE"+x.split("|")[1]
                                     if x.startswith("Terraced") else
                                     "DDFL"+x.split("|")[1]
                                     if x.startswith("Flats") else
                                     "DNDO00"
                                     if x=="Non-domestic" else 0 for x
                                     in dndstock["PROPERTY_TYPE"]]
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        dndstock["PROPERTY_TYPE"] = dndstock["PROPERTY_TYPE"].replace({"Flats":"DDFL00",
                                                                       "Detached":"DDDE00",
                                                                       "Semi-detached":"DDSD00",
                                                                       "Terraced":"DDTE00",
                                                                       "Non-domestic":"DNDO00"})
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )
    dndstock["PROPERTY_TYPE"] = dndstock["EFF"] + dndstock["PROPERTY_TYPE"]
    dndstock = dndstock.set_index([slagg.name,"PROPERTY_TYPE","EFF"])    
    
    # create and rearrange dataframe for cost (domestic and non-domestic)
    cost_cap = dndstock["cost"].reset_index()
    cost_cap["YEAR"] = ":*"
    cost_cap["MODE_OF_OPERATION"] = ""
    cost_cap["UNIT"] = "million 2015£/GW"
    cost_cap["VARIABLE"] = "CapitalCost"
        
    cost_cap = cost_cap.drop("EFF", axis=1)
    
    cost_cap = cost_cap.rename(columns={slagg.name:"REGION",
                                        "PROPERTY_TYPE":"TECHNOLOGY",
                                        "cost":"VALUE"})
    cost_cap = cost_cap.set_index([col for col 
                                   in cost_cap.columns
                                   if col!="VALUE"])    
    # save to file
    # unit: million [Baseyear]£/GW (= [Baseyear]£/kW) 
    cost_cap.to_csv(snakemake.output.path_building_measures_cost)
    
    # create and rearrange dataframe for constraint (domestic and non-domestic)
    cap_constr = dndstock["pot_peak_demand_red"].reset_index()
    cap_constr["YEAR"] = "2021"
    cap_constr["UNIT"] = "GW"

    cap_constr = cap_constr.drop("EFF", axis=1)
    
    cap_constr["VARIABLE"] = "TotalAnnualMaxCapacity"
    cap_constr = cap_constr.rename(columns={slagg.name:"REGION",
                                            "PROPERTY_TYPE":"TECHNOLOGY",
                                            "pot_peak_demand_red":"VALUE"})
    cap_constr = cap_constr.set_index([col for col 
                                   in cap_constr.columns
                                   if col!="VALUE"]) 
    # constraint for historical years (further interpolated when creating input
    # data set)
    init_cap_constr =  cap_constr.copy()
    init_cap_constr = init_cap_constr.rename({"2021":"2020"})
    init_cap_constr["VALUE"] = 0
    cap_constr = pd.concat([init_cap_constr,cap_constr])
    
    # save to file
    # unit: million [Baseyear]£/GW (= [Baseyear]£/kW)
    cap_constr.to_csv(snakemake.output.path_building_measures_con)
    
    # create and rearrange dataframe for lifetime (domestic and non-domestic)
    lifetime = dndstock["lifetime"].reset_index()
    lifetime["UNIT"] = "a"
    lifetime["VARIABLE"] = "OperationalLife"
    # round to int
    lifetime["lifetime"] = lifetime["lifetime"].round(0).astype(int)
    # round to int
 
    lifetime = lifetime.drop("EFF", axis=1)
    lifetime = lifetime.rename(columns={slagg.name:"REGION",
                                        "PROPERTY_TYPE":"TECHNOLOGY",
                                        "lifetime":"VALUE"})
    lifetime = lifetime.set_index([col for col 
                                   in lifetime.columns
                                   if col!="VALUE"])      
    
    # save to file
    # unit: a
    lifetime.to_csv(snakemake.output.path_building_measures_lt)
    
    # create and rearrange dataframe for deployment constraint (domestic and non-domestic)
    depcon = dndstock[["2025","2030","2035"]].stack().reset_index()
    depcon["UNIT"] = "GW/a"
    depcon["VARIABLE"] = "TotalAnnualMaxCapacityInvestment"

   
    depcon = depcon.drop("EFF", axis=1)
    depcon = depcon.rename(columns={slagg.name:"REGION",
                                        "PROPERTY_TYPE":"TECHNOLOGY",
                                        0:"VALUE",
                                        "level_3":"YEAR"})
    depcon = depcon.set_index([col for col 
                                   in depcon.columns
                                   if col!="VALUE"])      
    
    # save to file
    # unit: GW/a
    depcon.to_csv(snakemake.output.path_building_measures_dc)
