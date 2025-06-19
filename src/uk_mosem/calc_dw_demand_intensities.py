"""
Script deriving demand intensities for different domestic property types


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

import pandas as pd
import numpy as np

import utils


if __name__ == "__main__":
    
    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_dw_demand_intensities",
                                         **df_params.iloc[43].to_dict())
    # load list of areas/area lookup (LAD23CD)
    lsoas = utils.get_entity_lookup(["LSOA","LAD","CTRY"])
    
    # calculate intensity based on LA-level NEED data
    # this is currently not used but can in future be implemented as an option
    """
    # load NEED data per LAD (EW: LAD21CD, single value used for S merged to 
    # LAD23CD) data for Isles of Scilly are missing, City of London has 
    # no (values for) detached and semi-detached properties
    need_ew_data = pd.read_excel(snakemake.input.path_need_con_ew,
                                 sheet_name=["LA5","LA6"], header=[5,7]) 
    # header=[4,6] for 2017 data
    
    need_s_data = pd.read_excel(snakemake.input.path_need_con_s,
                     sheet_name=["Table A3","Table A4"],
                     usecols=[5],skiprows=137,nrows=10)
    # dataskiprows=105 for 2017 data
    
    # loop through gas and electricity data for both Scotland and EW
    need = dict()
    for need_s,need_ew,f in zip(need_s_data.values(),
                                need_ew_data.values(),["gas","el"]):
        # process the few Scotland values
        
        # name and reshape data, drop irrelevant
        need_s = need_s.dropna(axis=0, how="all")
        need_s.index = ["All","Houses","Detached","Semi-detached","Terraced",
                       "Flats","Unknown"]
        need_s = need_s.loc[["Detached","Semi-detached","Terraced",
                       "Flats"]].T
        
        # add country code as index, data will be applied through out Scotland
        need_s.index=["S92000003"]
        need_s.index.name="CTRY11CD"
        
        # expand data to cover all LAD codes in Scotland, drop country code
        need_s = need_s.merge(right=lsoas[["LAD23CD","CTRY11CD"]],
                            on="CTRY11CD").drop_duplicates()
        need_s = need_s.set_index(["LAD23CD"])
        need_s = need_s.drop("CTRY11CD", axis=1)
    
    
        # process EW data
        
        # drop empty rows/columns and arrange data
        need_ew.iloc[0,:] = need_ew.iloc[0,:].str.replace("\n","")
        need_ew.iloc[1,:] = need_ew.iloc[1,:].str.replace("\n","")
        need_ew = need_ew.dropna(axis=0, how="all")
        need_ew = need_ew.dropna(axis=1, how="all")
        need_ew = need_ew.drop(columns=["Unnamed: 0_level_0","Unnamed: 1_level_0"],
                               level=0)
        need_ew = need_ew.set_index(("Unnamed: 2_level_0","Unnamed: 2_level_1"),
                                    append=False)
        need_ew.index = need_ew.index.rename("LAD21CD")
        
        need_ew.iloc[0,:] = need_ew.iloc[0,:].fillna("All dwellings")
    
        idx = pd.MultiIndex.from_product([["number",f+"_mean",f+"_median"],
                                         need_ew.iloc[0,0:8]],
                                         names=['var', 'type'])
        need_ew.columns = idx
        need_ew = need_ew.iloc[1:,:]
        
        # convert to numeric and calculate aggregates (flats and semi-det. 
        # - also including bungalows), etc.
        need_ew = need_ew.apply(pd.to_numeric,errors='coerce')      
        need_ew[f+"_mean_p","Flats"] = ((need_ew[f+"_mean","Converted flat"]
                                      *need_ew["number","Converted flat"]
                                      +need_ew[f+"_mean","Purpose built flat"]
                                      *need_ew["number","Purpose built flat"])
                                      /(need_ew["number","Purpose built flat"]
                                        +need_ew["number","Converted flat"]))
        need_ew[f+"_mean_p","Terraced"] = ((need_ew[f+"_mean","Mid terrace"]
                                      *need_ew["number","Mid terrace"]
                                      +need_ew[f+"_mean","End terrace"]
                                      *need_ew["number","End terrace"])
                                      /(need_ew["number","End terrace"]
                                        +need_ew["number","Mid terrace"]))
        need_ew[f+"_mean_p","Semi-detached"] = ((need_ew[f+"_mean","Semi detached"]
                                      *need_ew["number","Semi detached"]
                                      +need_ew[f+"_mean","Bungalow"]
                                      *need_ew["number","Bungalow"])
                                      /(need_ew["number","Semi detached"]
                                        +need_ew["number","Bungalow"]))
        need_ew[f+"_mean_p","Detached"] = need_ew[f+"_mean","Detached"]

        # drop irrelvant data 
        need_ew = need_ew[f+"_mean_p"]
        
        # update to new LADCD version, if necessary
        need_ew = utils.update_LADCD(need_ew,from_CD="LAD21CD", how="mean")
        
        # concat S and EW data
        need[f] = pd.concat([need_ew,need_s],sort=True)
        need[f].index.name = "LAD23CD"
        need[f].columns.name = "PROPERTY_TYPE"
        need[f] = need[f].stack(dropna=False).to_frame()
        need[f].columns = [f+"_intensity"]
        
        # fill in missing values gas intensity values
        # assume Westminster values for City of London
        need[f].loc[(need[f].index.get_level_values(0)== "E09000001") &
                    (need[f][f+"_intensity"].isna()),f+"_intensity"]= pd.concat([need[f].loc["E09000033",f+"_intensity"]],keys=['E09000001'], names=['LAD23CD'])
        if f == "gas":  
        # assume Cornwall values for Isles of Scilly
            need[f].loc[(need[f].index.get_level_values(0)== "E06000053") &
                        (need[f][f+"_intensity"].isna()),f+"_intensity"]= pd.concat([need[f].loc["E06000052",f+"_intensity"]],keys=['E06000053'], names=['LAD23CD'])
            
        # drop LAD codes and arrange final dataframe structure
        need[f] = lsoas[["LSOA11CD","LAD23CD","LAD23NM"]].merge(right=need[f].reset_index(),
                            on="LAD23CD",how="left").drop_duplicates()
        need[f] = need[f].drop(columns=["LAD23CD","LAD23NM"]) 
        
        # disaggregate households if necessary
        if snakemake.params.dic["scen_hh_disagg"]=="-":
            pass
        elif snakemake.params.dic["scen_hh_disagg"]=="T":
            need[f] = utils.explode_hh(need[f],["OO","RP","RS"])
        elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
            need[f] = utils.explode_hh(need[f],
                                      [t+str(i) for t in ["O","P","S"]
                                       for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)]
                                      )    
        else:
            raise NotImplementedError(
                    f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                    " not a valid option of the 'scen_hh_disagg' parameter."
                    " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
                )
        need[f] = need[f].set_index(["LSOA11CD","PROPERTY_TYPE"])
        



        
    
    # scale intensities to match LSOA gas consumption, not used
    # calibration is included for annual demand
    
    # # load dwelling stock data
    # dw_stock = pd.read_csv(PATHP+PATH_DW_STOCK_EX,usecols=("LSOA11CD",
    #                                                         "PROPERTY_TYPE",
    #                                                         "2019"))
    # dw_stock = dw_stock.merge(right=dw_stock.groupby(["LSOA11CD"]).sum(),
    #                           on=["LSOA11CD"], how="left", suffixes=["",
    #                                                                   "_total"]) 
    # # loop through fuels, i.e., gas and electricity
    # for dw_con_path,f in zip([PATH_DW_CON_G,PATH_DW_CON_E],
    #                             ["gas","el"]):
        
    #     # load consumption data
    #     dw_con = pd.read_csv(PATHR+dw_con_path)
    #     dw_con = dw_con.rename(columns={"LSOACode":"LSOA11CD"})
    #     dw_con.loc[:,"METERS":] = dw_con.loc[:,"METERS":].apply(pd.to_numeric,
    #                                     errors='coerce').fillna(0)
    #     # if f == "gas":
    #     #     dw_con["METERS"] = (dw_con["METERS"]
    #     #                             - dw_con["NON_CONSUMING_METERS"]) 
    #     dw_con = dw_con[dw_con["YEAR"]==2019]
        
    
    #     # merge LAD-based intensity data and LSOA consumption data
    #     dw_stock_f = dw_stock.merge(right=need[f],on=["LSOA11CD",
    #                                                   "PROPERTY_TYPE"],
    #                               how="left")
    #     dw_stock_f = dw_stock_f.merge(right=dw_con,on=["LSOA11CD"],how="left")    
        
    #     # fill gas data with national sum (i.e., averages for intensity
    #     # scaling) where no gas data given
    #     dw_stock_f[["METERS",
    #                 "KWH"]] = dw_stock_f[["METERS",
    #                                       "KWH"]].fillna(dw_stock_f[["METERS",
    #                                                         "KWH"]].sum())
    #     # fill stock data for two Scottish LSOAs/data zones with no properties
    #     # to nevertheless derive intensities
    #     dw_stock_f.loc[dw_stock_f["2019_total"]==0,
    #                     "2019"] = dw_stock_f[["PROPERTY_TYPE",
    #                                           "2019"]].groupby("PROPERTY_TYPE").sum()["2019"].to_list()*2
    #     dw_stock_f.loc[dw_stock_f["2019_total"]==0,
    #                     "2019_total"] = dw_stock_f.loc[dw_stock_f["2019_total"]> 0,"2019"].sum()
        
    #     # scale intensity data
    #     dw_stock_f["stock_"+f] = (dw_stock_f["2019"] / dw_stock_f["2019_total"]
    #                               * dw_stock_f["METERS"])
    #     dw_stock_f["sic"] = dw_stock_f["stock_"+f] * dw_stock_f[f+"_intensity"]
        
    #     dw_stock_f = dw_stock_f.merge(right=dw_stock_f.groupby(["LSOA11CD"]).sum()["sic"],
    #                               on=["LSOA11CD"], how="left", suffixes=["",
    #                                                                       "_total"])
    #     dw_stock_f[f+"_intensity"] = (dw_stock_f[f+"_intensity"]) 
    #                                   # * dw_stock_f["KWH"]
    #                                   # / dw_stock_f["sic_total"])
        
    #     # set meter number again to zero where no gas data given and calculate
    #     # fraction of properties not connected to gas grid
    #     dw_stock_f.loc[dw_stock_f["LAName"].isna(),"METERS"] = 0   
    #     dw_stock_f["offgrid"] = 1-(dw_stock_f["METERS"]
    #                                 /dw_stock_f["2019_total"]).clip(0,1)
        
    #     need[f] =  dw_stock_f.set_index(["LSOA11CD",
    #                                       "PROPERTY_TYPE"])[[f+"_intensity",
    #                                                         "offgrid"]]
  
        
        
    # calculate space and hot water demand intensity per dwelling based on gas
    # consumption, fraction of gas used for space/water heating in domestic
    # sector, and efficiency of gas boilers
    # unit: kWh/a/dwelling    
  
    # load raw data set
    dem_raw_data = pd.read_csv(snakemake.input.path_set_dw_int)
    eff_raw_data = pd.read_csv(snakemake.input.path_set_eff)
      
    need["gas"]["SH_intensity"] = (need["gas"]["gas_intensity"]
                           * dem_raw_data["VALUE"][dem_raw_data["VARIABLE"]
                                            =="NGDomFracSpaceHeating"].values
                           * eff_raw_data[["TECHNOLOGY","FUEL_IN","UNIT",
                                                "VALUE"]][(eff_raw_data["VARIABLE"]=="Efficiency")&
                                                          (eff_raw_data["FUEL_IN"]=="NGD000")&
                                                          (eff_raw_data["MODE_OF_OPERATION"]==1)&
                                                          (eff_raw_data["TECHNOLOGY"]=="NGBODDDE00,NGBODDSD00,NGBODDTE00,NGBODDFL00")]["VALUE"].squeeze())
    need["gas"]["HW_intensity"] = (need["gas"]["gas_intensity"]
                           * dem_raw_data["VALUE"][dem_raw_data["VARIABLE"]
                                            =="NGDomFracWaterHeating"].values
                           * eff_raw_data[["TECHNOLOGY","FUEL_IN","UNIT",
                                                "VALUE"]][(eff_raw_data["VARIABLE"]=="Efficiency")&
                                                          (eff_raw_data["FUEL_IN"]=="NGD000")&
                                                          (eff_raw_data["MODE_OF_OPERATION"]==3)&
                                                          (eff_raw_data["TECHNOLOGY"]=="NGBODDDE00,NGBODDSD00,NGBODDTE00,NGBODDFL00")]["VALUE"].squeeze())
 
    # calculate space heat demand intensity for new builts assuming a certain
    # reduction in comparison to existing buildings
    need["gas"]["SH_intensity_newbuilt"] = (need["gas"]["SH_intensity"] 
                                            * dem_raw_data["VALUE"][
                                            dem_raw_data["VARIABLE"]
                                            =="DemandFracNewbuilt"].values)
    

    # load residual capacity fractions
    rescap = pd.read_csv(snakemake.input.path_res_caps_ht_f,
                        index_col=["LSOA11CD","PROPERTY_TYPE","TECHNOLOGY"])
    rescap = rescap.drop(["Non-domestic"],level=1)
    rescap = rescap.drop(["RAUP","WDIS"],level=2)
    
    hsys_dw = pd.read_csv(snakemake.input.path_dw_stock_hde, index_col=["LSOA11CD",
                                                              "PROPERTY_TYPE",
                                                              "HEATING_SYSTEM"])
    hsys_dw = hsys_dw.loc[hsys_dw.index.get_level_values("HEATING_SYSTEM").str.startswith("MH")]
    hsys_dw.columns = ["NUMBER"]
    hsys_dw.index  = hsys_dw.index.set_names("TECHNOLOGY",level=2)
    hsys_dw = hsys_dw.reset_index()
    hsys_dw["TECHNOLOGY"] = hsys_dw["TECHNOLOGY"].str[3:]
    hsys_dw = hsys_dw.set_index(["LSOA11CD","PROPERTY_TYPE","TECHNOLOGY"])
    hsys_dw = hsys_dw.unstack(fill_value=0).unstack(fill_value=0).stack().stack()
    
    # calculate non-heat gas demand
    need["gas"]["NHG_intensity"] = (need["gas"]["gas_intensity"]
                           * (1-dem_raw_data.loc[dem_raw_data["VARIABLE"]
                                            =="NGDomFracSpaceHeating",
                                            "VALUE"].values
                              -dem_raw_data.loc[dem_raw_data["VARIABLE"]
                                            =="NGDomFracWaterHeating",
                                            "VALUE"].values)
                           * rescap.loc[rescap.index.get_level_values(2).str.contains("NGBO"),"fraction"].droplevel(2))
    
    # FIXME: this can be deleted - start of calculations to use ECUK data to
    # derive average NHE
    # hsys_dw = pd.read_csv(snakemake.input.path_dw_stock_hde, index_col=["LSOA11CD",
    #                                                           "PROPERTY_TYPE",
    #                                                           "HEATING_SYSTEM"])
    # hsys_dw = hsys_dw.loc[hsys_dw.index.get_level_values("HEATING_SYSTEM").str.startswith("MH")]

    # hsys_dw = hsys_dw.groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
    
    # need["el"]["NHE_intensity"] = 
    
    # avg_nhe =  ((need["el"]["el_intensity"]
    #            *hsys_dw["NUMBER_OF_PROPERTIES"])
    #            .groupby(["PROPERTY_TYPE"]).sum()
    #            /hsys_dw["NUMBER_OF_PROPERTIES"]
    #            .groupby(["PROPERTY_TYPE"]).sum())
    
    # ECUK = pd.read_excel(snakemake.input.path_ecuk_enduse,
    #                      sheet_name="Table U2",
    #                      usecols=(0,1,3,4,7),skiprows=4,nrows=6)
    # ECUK = ECUK.set_index(["Sector","End use","Year"]).sort_index()
    # ECUK_nhe_frac = (ECUK.loc[("Domestic",
    #                            ["Cooking/catering",
    #                             "Lighting and appliances"],
    #                            2022),"Electricity"].sum()
    #                  /ECUK.loc[("Domestic","Overall total",2022),"Electricity"])
    
    # noncapf = (ECUK.loc["Industry",
    #                         ["High temperature process",
    #                               "Low temperature process",
    #                               "Drying/separation",
    #                               "Other",
    #                               "Motors",
    #                               "Compressed air"
    #                               ],
    #                         slice("2017","2022")
    #                        ].groupby("Year").sum()
    #                /
    #                ECUK.loc[slice("Industry","Services"),
    #                        "Overall total",
    #                         slice("2017","2022"),
    #                        ].groupby("Year").sum()
    #                )
    
    
    # calculate average residual fraction per LAD/country for Scotland
    # for each LSOA (otherwise
    # LAD NEED intensity data will be used in conjunction with LSOA
    # residual fractions - better to use LAD averages for both)
    
    # calculate LAD average used along with NEED data for demand intensities
    rescap_avg = rescap.reset_index()
    rescap_avg = rescap_avg.merge(right=hsys_dw.groupby(["LSOA11CD",
                                                     "PROPERTY_TYPE"]).sum(),
                                 on=["LSOA11CD","PROPERTY_TYPE"],
                                 how="left")
    rescap_avg["fraction"] = (rescap_avg["fraction"]*rescap_avg["NUMBER"]).fillna(0)
    
    # groupby all Scottish values to get overall Scottish average (choosing
    # arbitrary Scottish LSOA to aggregate)
    # rescap_avg = rescap_avg.reset_index()
    rescap_avg.loc[:,"LSOA11CD"] = [l if not l.startswith("S")
                                  else "S01008445"
                                  for l in rescap_avg.loc[:,"LSOA11CD"]]
    rescap_avg = rescap_avg.groupby(["LSOA11CD","PROPERTY_TYPE","TECHNOLOGY"]).sum()
    # groupby LADs
    rescap_avg = utils.groupby_LAD(rescap_avg).sum()
    #calculate weighted average
    rescap_avg["fraction"] = (rescap_avg["fraction"]/rescap_avg["NUMBER"]).fillna(0)
    rescap_avg = rescap_avg.drop("NUMBER",axis=1)
    
    # extend Scottish LAD to all Scottish LADs
    rescap_avg = rescap_avg.rename(index={"S12000036":",".join([lad for lad 
                                                       in lsoas["LAD23CD"].unique()
                                                       if lad.startswith("S")])})
    rescap_avg = rescap_avg.reset_index()
    rescap_avg = rescap_avg.assign(LAD23CD=rescap_avg['LAD23CD'].str.split(',')).explode('LAD23CD')
    rescap_avg = rescap_avg.set_index(["LAD23CD","PROPERTY_TYPE","TECHNOLOGY"])

    rescap_avg = lsoas[["LSOA11CD","LAD23CD"]].merge(right=rescap_avg.reset_index(),
                                                     on=["LAD23CD"],
                                                     how="left")
    rescap_avg = rescap_avg.drop("LAD23CD",axis=1)
    rescap_avg = rescap_avg.set_index(["LSOA11CD","PROPERTY_TYPE","TECHNOLOGY"])
    
    
    # calculate non-heat electricity demand
    need["el"]["NHE_intensity"] = (need["el"]["el_intensity"]
                                    -(rescap_avg.loc[rescap_avg.index.get_level_values(2).str.contains("ASHP"),"fraction"].droplevel(2)
                                     *(need["gas"]["SH_intensity"]
                                       /eff_raw_data.loc[(eff_raw_data["TECHNOLOGY"]
                                                      ==("ASHPDDDE00,ASHPDDSD00,"+
                                                         "ASHPDDTE00,ASHPDDFL00"))&(
                                                      eff_raw_data["YEAR"]
                                                                     =="2020")&
                                                   (eff_raw_data["MODE_OF_OPERATION"]==1), "VALUE"].values
                                          +need["gas"]["HW_intensity"]
                                          /eff_raw_data.loc[(eff_raw_data["TECHNOLOGY"]
                                                         ==("ASHPDDDE00,ASHPDDSD00,"+
                                                            "ASHPDDTE00,ASHPDDFL00"))&(
                                                         eff_raw_data["YEAR"]
                                                                        =="2020")&
                                                      (eff_raw_data["MODE_OF_OPERATION"]==2), "VALUE"].values
                                          )
                                    ).fillna(0)
                                    -(rescap_avg.loc[rescap_avg.index.get_level_values(2).str.contains("GSHP"),"fraction"].droplevel(2)
                                    *(need["gas"]["SH_intensity"]
                                      /eff_raw_data.loc[(eff_raw_data["TECHNOLOGY"]
                                                     ==("GSHPDDDE00,GSHPDDSD00,"+
                                              "GSHPDDTE00,GSHPDDFL00"))&(
                                                     eff_raw_data["YEAR"]
                                                                    =="2020")&
                                                  (eff_raw_data["MODE_OF_OPERATION"]==1), "VALUE"].values
                                          +need["gas"]["HW_intensity"]
                                    /eff_raw_data.loc[(eff_raw_data["TECHNOLOGY"]
                                                   ==("GSHPDDDE00,GSHPDDSD00,"+
                                            "GSHPDDTE00,GSHPDDFL00))&(
                                                   eff_raw_data["YEAR"]
                                                                  =="2020")&
                                                (eff_raw_data["MODE_OF_OPERATION"]==2), "VALUE"].values)).fillna(0)
                                    -(rescap_avg.loc[rescap_avg.index.get_level_values(2).str.contains("ELRE"),"fraction"].droplevel(2)
                                    *(need["gas"]["SH_intensity"]
                                      /eff_raw_data.loc[(eff_raw_data["TECHNOLOGY"]
                                                     ==("ELREDDDE00,ELREDDSD00,"+
                                              "ELREDDTE00,ELREDDFL00"))&
                                                        (eff_raw_data["MODE_OF_OPERATION"]==1),"VALUE"].values
                                          +need["gas"]["HW_intensity"]
                                    /eff_raw_data.loc[(eff_raw_data["TECHNOLOGY"]
                                                   ==("ELREDDDE00,ELREDDSD00,"+
                                            "ELREDDTE00,ELREDDFL00"))&
                                                      (eff_raw_data["MODE_OF_OPERATION"]==2),"VALUE"].values)).fillna(0)
                                    -(rescap_avg.loc[rescap_avg.index.get_level_values(2).str.contains("ELST"),"fraction"].droplevel(2)
                                    *(need["gas"]["SH_intensity"]
                                      /eff_raw_data.loc[(eff_raw_data["TECHNOLOGY"]
                                                     ==("ELSTDDDE00,ELSTDDSD00,"+
                                              "ELSTDDTE00,ELSTDDFL00"))&
                                                        (eff_raw_data["MODE_OF_OPERATION"]==1),"VALUE"].values
                                          +need["gas"]["HW_intensity"]
                                    /eff_raw_data.loc[(eff_raw_data["TECHNOLOGY"]
                                                   ==("ELSTDDDE00,ELSTDDSD00,"+
                                            "ELSTDDTE00,ELSTDDFL00"))&
                                                      (eff_raw_data["MODE_OF_OPERATION"]==2),"VALUE"].values)).fillna(0))
    
    # values falling below the first percentile are adjusted including for a 
    # small number (around 100) of LSOA which are below 0), replaced with the 
    # first percentile
    perc = need["el"]["NHE_intensity"].groupby("PROPERTY_TYPE").quantile(0.01)
    need["el"] = need["el"].reset_index().merge(right=perc, how="left",on="PROPERTY_TYPE",
                                suffixes=["","_percentile"])
    need["el"] = need["el"].set_index(["LSOA11CD","PROPERTY_TYPE"])
    need["el"].loc[:,"NHE_intensity"] = need["el"].loc[:,["NHE_intensity","NHE_intensity_percentile"]].max(axis=1)
    
    # need["el"].loc[need["el"]["NHE_intensity"]<0,"NHE_intensity"] = (
    #     (need["el"].loc[need["el"]["NHE_intensity"]<0,"NHE_intensity"]*0).add(
    #         need["el"].loc[need["el"]["NHE_intensity"]>0,"NHE_intensity"].groupby("PROPERTY_TYPE").mean(),level=1))                                            
    """
    
    # calculate intensity based on record-level NEED data
    # "if" implemented to later change if done as option
    if True:
        
        # load and process NEED consumption record data
        epcc = pd.read_csv(snakemake.input.path_need_conr_ew,
                           usecols=["PROP_TYPE","EPC","REGION",
                                    "MAIN_HEAT_FUEL",
                                    "Gcons2019","Econs2019"])
        # only use gas heated buildings (to avoid issue with identifying NHE)
        epcc = epcc.loc[epcc["MAIN_HEAT_FUEL"]==1]
        # remove entries with no consumption data
        epcc = epcc.dropna()
        epcc["PROP_TYPE"] = epcc["PROP_TYPE"].replace({"Bungalow":"Semi-detached",
                                                       "Semi detached":"Semi-detached",
                                                       "Mid terrace":"Terraced",
                                                       "End terrace":"Terraced",
                                                       "Flat":"Flats"})
        # calculate average across EPC to use for HW and NHG
        epcc_avg = epcc.groupby(["PROP_TYPE"]).mean()
        epcc_avg.index.names = ["PROPERTY_TYPE"]
        
        # calculate mean across respective entries
        epcc = epcc.groupby(["REGION","PROP_TYPE","EPC"]).mean()
        
        # # calculate relative consumption for each EPC rating
        # # with regard to average
        # epcc = epcc/epcc_avg
        
        epcc.index.names = ["RGN21CD","PROPERTY_TYPE","EPC_RATING"]
        # apply North East to Scotland
        epcc = pd.concat([epcc,
                          epcc.loc[("E12000001",
                                    slice(None),
                                    slice(None))].rename(index={"E12000001":"S92000003"})])
        
        # merge on LSOAs
        rgn = utils.get_entity_lookup(["LSOA","RGN"])
        
        epcc = rgn.merge(epcc.reset_index(),how="left",on="RGN21CD")
        epcc = epcc.drop("RGN21CD",axis=1)
        epcc = epcc.set_index(["LSOA11CD","PROPERTY_TYPE","EPC_RATING"])
        epcc = epcc.drop("No EPC",level="EPC_RATING")
        
        # disaggregate households if necessary
        if snakemake.params.dic["scen_hh_disagg"]=="-":
            pass
        elif snakemake.params.dic["scen_hh_disagg"]=="T":
            epcc = utils.explode_hh(epcc,["OO","RP","RS"])
            epcc_avg = utils.explode_hh(epcc_avg,["OO","RP","RS"])
        elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
            epcc = utils.explode_hh(epcc,
                                      [t+str(i) for t in ["O","P","S"]
                                       for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)]
                                      )
            epcc_avg = utils.explode_hh(epcc_avg,
                                      [t+str(i) for t in ["O","P","S"]
                                       for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)]
                                      )  
        else:
            raise NotImplementedError(
                    f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                    " not a valid option of the 'scen_hh_disagg' parameter."
                    " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
                )
        
        # load EPC stock data
        stock = pd.read_csv(snakemake.input.path_dw_stock_epc)
        stockn = stock.groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
        stock = stock.replace({"A":"A/B",
                               "B":"A/B",
                               "F":"F/G",
                               "G":"F/G"})
        stock = stock.loc[stock["EPC_RATING"]!='INVALID!']
        
        # calculate fraction of EPC ratings
        stock = stock.groupby(["LSOA11CD","PROPERTY_TYPE","EPC_RATING"]).sum()
        stock = stock/stock.groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
        stock = stock.fillna(0)
        
        ### calc gas/el intensities
        
        # get consumption of A/B to apply to new builds
        epcc_new = epcc.xs("A/B",level="EPC_RATING")
        
        # calculate average consumption taking into account consumption
        # for the respective rating and number of properties of each
        epcc = epcc.merge(stock,left_index=True,right_index=True,how="left")
        epcc.loc[:,"NUMBER_OF_PROPERTIES"] = epcc.loc[:,"NUMBER_OF_PROPERTIES"].fillna(0)
        
        # for property types not existent in 2019 assume any potential new built
        # (in the later base years) has A/B rating.
        epcc.loc[(epcc.groupby(["LSOA11CD",
                                "PROPERTY_TYPE"])
                  ["NUMBER_OF_PROPERTIES"].transform("sum")==0)&
                 (epcc.index.get_level_values("EPC_RATING")=="A/B"),
                 "NUMBER_OF_PROPERTIES"] = 1
        
        epcc.loc[:,"Gcons2019"] = epcc.loc[:,"Gcons2019"]*epcc.loc[:,"NUMBER_OF_PROPERTIES"]
        epcc.loc[:,"Econs2019"] = epcc.loc[:,"Econs2019"]*epcc.loc[:,"NUMBER_OF_PROPERTIES"]
        
        # calculate average consumption
        epcc = epcc.groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()

        
        # FIXME: potential implement scaling to match LAD NEED data
        # epcc = epcc.merge(stockn,left_index=True,right_index=True,how="left")


        ### calculate useful heat and non-heat intensities
        
        # set national averages for HW and NHG consumption per property type
        # load raw data set
        dem_raw_data = pd.read_csv(snakemake.input.path_set_dw_int)
        eff_raw_data = pd.read_csv(snakemake.input.path_set_eff)
        
        
        epcc_avg["Gcons2019_HW"] = (epcc_avg["Gcons2019"]* dem_raw_data["VALUE"][dem_raw_data["VARIABLE"]
                         =="NGDomFracWaterHeating"].values)
        epcc_avg["Gcons2019_NHG"] = (epcc_avg["Gcons2019"]* (1-dem_raw_data.loc[dem_raw_data["VARIABLE"]
                         =="NGDomFracSpaceHeating",
                         "VALUE"].values
           -dem_raw_data.loc[dem_raw_data["VARIABLE"]
                         =="NGDomFracWaterHeating",
                         "VALUE"].values))
        epcc_avg = epcc_avg.drop(["Gcons2019","Econs2019"],axis=1)
        epcc = epcc.reset_index().merge(epcc_avg,on="PROPERTY_TYPE",how="left")
        epcc = epcc.set_index(["LSOA11CD","PROPERTY_TYPE"])
        
        need = dict()
        need["gas"] = ((epcc.loc[:,"Gcons2019"]-epcc["Gcons2019_HW"]
                        -epcc["Gcons2019_NHG"])
                               * eff_raw_data[["TECHNOLOGY","FUEL_IN","UNIT",
                                                    "VALUE"]][(eff_raw_data["VARIABLE"]=="Efficiency")&
                                                              (eff_raw_data["FUEL_IN"]=="NGD000")&
                                                              (eff_raw_data["MODE_OF_OPERATION"]==1)&
                                                              (eff_raw_data["TECHNOLOGY"]=="NGBODDDE00,NGBODDSD00,NGBODDTE00,NGBODDFL00")]["VALUE"].squeeze())
        
        need["gas"] = need["gas"].to_frame()
        need["gas"].columns = ["SH_intensity"]                                            
        need["gas"]["HW_intensity"] = (epcc.loc[:,"Gcons2019_HW"]
                               * eff_raw_data[["TECHNOLOGY","FUEL_IN","UNIT",
                                                    "VALUE"]][(eff_raw_data["VARIABLE"]=="Efficiency")&
                                                              (eff_raw_data["FUEL_IN"]=="NGD000")&
                                                              (eff_raw_data["MODE_OF_OPERATION"]==3)&
                                                              (eff_raw_data["TECHNOLOGY"]=="NGBODDDE00,NGBODDSD00,NGBODDTE00,NGBODDFL00")]["VALUE"].squeeze())
     
        # calculate space heat demand intensity for new builts
        need["gas"]["SH_intensity_newbuilt"] = ((epcc_new.loc[:,"Gcons2019"]-epcc["Gcons2019_HW"]
                                                -epcc["Gcons2019_NHG"]).dropna()
                                                * eff_raw_data[["TECHNOLOGY","FUEL_IN","UNIT",
                                                                     "VALUE"]][(eff_raw_data["VARIABLE"]=="Efficiency")&
                                                                               (eff_raw_data["FUEL_IN"]=="NGD000")&
                                                                               (eff_raw_data["MODE_OF_OPERATION"]==1)&
                                                                               (eff_raw_data["TECHNOLOGY"]=="NGBODDDE00,NGBODDSD00,NGBODDTE00,NGBODDFL00")]["VALUE"].squeeze())
                         
        
        # load residual capacity fractions
        rescap = pd.read_csv(snakemake.input.path_res_caps_ht_f,
                            index_col=["LSOA11CD","PROPERTY_TYPE","TECHNOLOGY"])
        rescap = rescap.drop(["Non-domestic"],level=1)
        rescap = rescap.drop(["RAUP","WDIS"],level=2)
        
        hsys_dw = pd.read_csv(snakemake.input.path_dw_stock_hde, index_col=["LSOA11CD",
                                                                  "PROPERTY_TYPE",
                                                                  "HEATING_SYSTEM"])
        hsys_dw = hsys_dw.loc[hsys_dw.index.get_level_values("HEATING_SYSTEM").str.startswith("MH")]
        hsys_dw.columns = ["NUMBER"]
        hsys_dw.index  = hsys_dw.index.set_names("TECHNOLOGY",level=2)
        hsys_dw = hsys_dw.reset_index()
        hsys_dw["TECHNOLOGY"] = hsys_dw["TECHNOLOGY"].str[3:]
        hsys_dw = hsys_dw.set_index(["LSOA11CD","PROPERTY_TYPE","TECHNOLOGY"])
        hsys_dw = hsys_dw.unstack(fill_value=0).unstack(fill_value=0).stack().stack()
        
        # calculate non-heat gas intensity
        need["gas"]["NHG_intensity"] = (epcc.loc[:,"Gcons2019_NHG"]
                               * rescap.loc[rescap.index.get_level_values(2).str.contains("NGBO"),"fraction"].droplevel(2))
        

        # calculate non-heat electricity intensity
        need["el"] = epcc.loc[:,"Econs2019"]
                              
        need["el"] = need["el"].to_frame()
        need["el"].columns = ["NHE_intensity"]  
        
        # # calc adjustment factor
        # f = stock.merge(epcc,left_index=True,right_index=True,how="left")
        # f["ADJ"] = f["NUMBER_OF_PROPERTIES"]*f["Gcons2019"]
        # f = f.groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
        
        # # apply to gas intensity
        # need["gas"]["gas_intensity"] = need["gas"]["gas_intensity"].multiply(f["ADJ"],fill_value=1)
        # need["gas"].merge(f["ADJ"],left_index=True,right_index=True,how="left")
        # need["gas"]["gas_intensity"].multiply(f["ADJ"],fill_value=1)

    # create separate dataframe and project into future modelling years
    demand_sh_ex = need["gas"]["SH_intensity"].rename("2015").to_frame()
    for i in range(2016,2061):
        demand_sh_ex[str(i)] = demand_sh_ex["2015"]

    demand_hw = need["gas"]["HW_intensity"].rename("2015").to_frame()
    for i in range(2016,2061):
        demand_hw[str(i)] = demand_hw["2015"]
        
    demand_sh_nb = need["gas"]["SH_intensity_newbuilt"].rename("2015").to_frame()
    for i in range(2016,2061):
        demand_sh_nb[str(i)] = demand_sh_nb["2015"]
    
    # assume linear decline in NHG usage towards 2050
    demand_nhg = need["gas"]["NHG_intensity"].rename("2015").to_frame()
    for i in range(2016,2061):

        demand_nhg[str(i)] = demand_nhg["2015"]

    demand_nhgc = demand_nhg.copy()
    
    # overwrite NHG intensity with trajectory based on run
    if snakemake.params.dic["scen_nhg"] != "-":
        
        sets = pd.read_csv(snakemake.input.path_set_sets)
        years = sets.loc[(sets["SET"]=="YEAR"),
                         "VALUE"].to_frame() 
        dev = years.copy()
        dev = dev.rename(columns={"VALUE":"YEAR"})
        dev["VALUE"] = np.nan
        vys = snakemake.params.dic["scen_nhg"].split("a")
        
        for  vy in vys:
            value, year = vy.split("y")
            dev.loc[dev["YEAR"]==year,"VALUE"] = value
            
        dev["VALUE"] = dev["VALUE"].astype(float)
        
        # interpolate
        dev["VALUE"] = dev["VALUE"].interpolate(limit_direction="both")
        
        for i in dev["YEAR"]:
            demand_nhg.loc[(slice(None),slice(None)),str(i)] = (
                demand_nhg.loc[(slice(None),slice(None)),str(i)]
                * dev.loc[dev["YEAR"]==i,"VALUE"].squeeze())
        
            
    if ((snakemake.params.dic["scen_local_gov"] != "-") and
        (snakemake.params.dic["scen_local_gov"] != "A")):
        et = pd.read_csv(snakemake.input.path_loc_gov_em)
        sets = pd.read_csv(snakemake.input.path_set_sets)
        years = sets.loc[(sets["SET"]=="YEAR"),
                         "VALUE"].to_frame()
        years["VALUE"] = years["VALUE"].astype(int)
        
        et = et[["REGION","YEAR","VALUE"]].set_index(["REGION","YEAR"]).unstack().droplevel(0,axis=1)
        for y in years["VALUE"]:
            if y not in et.columns:
                et.loc[:,y] = np.nan
        et = et.sort_values(by="YEAR",axis=1)
        
        et.loc[:,2015:2022]=1
        et = et.interpolate(axis=1)
        et = lsoas.merge(right=et.reset_index(),left_on="LAD23CD",
                         right_on="REGION",how="left").set_index("LSOA11CD")
        et = et.drop(["LAD23CD","REGION","LAD23NM","CTRY11CD"],axis=1)
        if snakemake.params.dic["scen_nhg"] != "-":
            for y in dev["YEAR"]:
                et.loc[:,int(y)] = et.loc[:,int(y)].fillna(dev.loc[dev["YEAR"]==y,"VALUE"].squeeze())
        else:
            et = et.fillna(1)
        
        et.columns=et.columns.astype(str)

        demand_nhg = demand_nhgc*et
                

    
    demand_nhe = need["el"]["NHE_intensity"].rename("2015").to_frame()
    for i in range(2016,2061):
        demand_nhe[str(i)] = demand_nhe["2015"]
        
    # overwrite NHE intensity with trajectory based on run
    if snakemake.params.dic["scen_nhe"] != "-":
        
        sets = pd.read_csv(snakemake.input.path_set_sets)
        years = sets.loc[(sets["SET"]=="YEAR"),
                         "VALUE"].to_frame() 
        dev = years.copy()
        dev = dev.rename(columns={"VALUE":"YEAR"})
        dev["VALUE"] = np.nan
        vys = snakemake.params.dic["scen_nhe"].split("a")
        
        for  vy in vys:
            value, year = vy.split("y")
            dev.loc[dev["YEAR"]==year,"VALUE"] = value
            
        dev["VALUE"] = dev["VALUE"].astype(float)
        
        # interpolate
        dev["VALUE"] = dev["VALUE"].interpolate(limit_direction="both")
        
        for i in dev["YEAR"]:
            demand_nhe.loc[(slice(None),slice(None)),str(i)] = (
                demand_nhe.loc[(slice(None),slice(None)),str(i)]
                * dev.loc[dev["YEAR"]==i,"VALUE"].squeeze())  
    
    # save to files
    # unit: kWh/a/dwelling
    demand_sh_ex.to_csv(snakemake.output.path_dw_int_sh)
    demand_sh_nb.to_csv(snakemake.output.path_dw_int_sh_nb)
    demand_hw.to_csv(snakemake.output.path_dw_int_hw)
    demand_nhg.to_csv(snakemake.output.path_dw_int_nhg)
    demand_nhe.to_csv(snakemake.output.path_dw_int_nhe)
