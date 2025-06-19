"""
Script deriving constraints for building heating technologies


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

import pandas as pd
import numpy as np

import utils

import logging

logger = logging.getLogger(__name__)
utils.adjust_logger()

if __name__ == "__main__":

    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_heat_techs",
                                         **df_params.iloc[0].to_dict())    

#%%
# biomass boiler constraint

    hsys = pd.read_csv(snakemake.input.path_res_caps_ht_f, index_col=["LSOA11CD",
                                                           "PROPERTY_TYPE",
                                                           "TECHNOLOGY"]) 
    # take only NGBO fraction (i.e., on-grid fraction) and calculate
    # offgrid fraction
    hsys = hsys[hsys.index.get_level_values("TECHNOLOGY")=="NGBO"]
    hsys = 1-hsys
    

    # load peak heat demands/technology capacity
    hpeaks = pd.read_csv(snakemake.input.path_peakc_lsoa,
                            index_col=["LSOA11CD",
                                       "PROPERTY_TYPE",
                                       "TECHNOLOGY"])
    
    #lsoas = utils.get_entity_lookup(["LSOA"])
    hsys = hpeaks.index.to_frame(index=False)[["LSOA11CD","PROPERTY_TYPE"]].drop_duplicates().merge(right=hsys.reset_index(),
                                                  on=["LSOA11CD","PROPERTY_TYPE"], how="outer")
    hsys["fraction"] = hsys["fraction"].fillna(1)
    hsys["TECHNOLOGY"] = hsys["TECHNOLOGY"].fillna("NGBO")
    hsys = hsys.set_index(["LSOA11CD","PROPERTY_TYPE","TECHNOLOGY"])
    
    hpeaks = hpeaks[hpeaks.index.get_level_values("TECHNOLOGY")=="BMBO"]
    
    if snakemake.params.dic["scen_hh_disagg"].startswith("T"):
            hpeaks.index = pd.MultiIndex.from_tuples([ (x[0],x[1],x[2]+"DDDE"+x[1].split("|")[1])
                                                  if x[1].startswith("Detached") else
                                                  (x[0],x[1],x[2]+"DDSD"+x[1].split("|")[1])
                                                  if x[1].startswith("Semi-detached") else
                                                  (x[0],x[1],x[2]+"DDTE"+x[1].split("|")[1])
                                                  if x[1].startswith("Terraced") else
                                                  (x[0],x[1],x[2]+"DDFL"+x[1].split("|")[1])
                                                  if x[1].startswith("Flats") else
                                                  (x[0],x[1],x[2]+"DNDO00")
                                                  if x[1]=="Non-domestic" else 0 for x
                                                  in hpeaks.index])
            hsys.index = pd.MultiIndex.from_tuples([ (x[0],x[1],x[2]+"DDDE"+x[1].split("|")[1])
                                                  if x[1].startswith("Detached") else
                                                  (x[0],x[1],x[2]+"DDSD"+x[1].split("|")[1])
                                                  if x[1].startswith("Semi-detached") else
                                                  (x[0],x[1],x[2]+"DDTE"+x[1].split("|")[1])
                                                  if x[1].startswith("Terraced") else
                                                  (x[0],x[1],x[2]+"DDFL"+x[1].split("|")[1])
                                                  if x[1].startswith("Flats") else
                                                  (x[0],x[1],x[2]+"DNDO00")
                                                  if x[1]=="Non-domestic" else 0 for x
                                                  in hsys.index])
    else:              
        hpeaks.index = pd.MultiIndex.from_tuples([ (x[0],x[1],x[2]+"DDDE00")
                                                         if x[1].startswith("Detached") else
                                                         (x[0],x[1],x[2]+"DDSD00")
                                                         if x[1].startswith("Semi-detached") else
                                                         (x[0],x[1],x[2]+"DDTE00")
                                                         if x[1].startswith("Terraced") else
                                                         (x[0],x[1],x[2]+"DDFL00")
                                                         if x[1].startswith("Flats") else
                                                         (x[0],x[1],x[2]+"DNDO00")
                                                         if x[1].startswith("Non-domestic") else 0 for x
                                                         in hpeaks.index])
        hsys.index = pd.MultiIndex.from_tuples([ (x[0],x[1],x[2]+"DDDE00")
                                                         if x[1].startswith("Detached") else
                                                         (x[0],x[1],x[2]+"DDSD00")
                                                         if x[1].startswith("Semi-detached") else
                                                         (x[0],x[1],x[2]+"DDTE00")
                                                         if x[1].startswith("Terraced") else
                                                         (x[0],x[1],x[2]+"DDFL00")
                                                         if x[1].startswith("Flats") else
                                                         (x[0],x[1],x[2]+"DNDO00")
                                                         if x[1].startswith("Non-domestic") else 0 for x
                                                         in hsys.index])
    hpeaks.index.names = ["LSOA11CD","PROPERTY_TYPE","TECHNOLOGY"]   
    hsys.index.names = ["LSOA11CD","PROPERTY_TYPE","TECHNOLOGY"] 
    
    hsys.index = hsys.index.droplevel("TECHNOLOGY")
    
    # load rural urban classification
    ruc = pd.read_csv(snakemake.input.path_ruc, usecols=[0,1,3])
    ruc = ruc[ruc["nation"]!="N"]
    ruc = ruc.rename(columns={"lsoa":"LSOA11CD"})
    ruc = ruc.drop("nation",axis=1)
    #ruc.loc[:,"ukruc-2"] = ruc.loc[:,"ukruc-2"].astype(bool)
    ruc = ruc.set_index("LSOA11CD")
    # ruc = ruc["ukruc-2"]
    
    # load suitability fractions
    raw_data = pd.read_csv(snakemake.input.path_set_htech)
    
    bsf = raw_data.loc[raw_data["VARIABLE"]=="BioBoilerSuitableFraction",
                       ["TECHNOLOGY","VALUE"]]
    bsf = bsf.set_index("TECHNOLOGY")
    
    # explode technologies if household disaggregated
    if snakemake.params.dic["scen_hh_disagg"]=="T":    
        bsf = utils.explode_hh(bsf,
                                  ["OO","RP","RS"],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"): 
        bsf = utils.explode_hh(bsf,
                                  [t+str(i) for t in ["O","P","S"]
                                   for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )
        
    bm_cap_constr = hpeaks.multiply(pd.concat([ruc["ukruc-2"].to_frame(name=str(i))
                                        for i in range(2015,2061)],
                                       axis=1),axis=0)
    bm_cap_constr = bm_cap_constr.multiply(pd.concat([bsf["VALUE"].to_frame(name=str(i))
                                                      for i in range(2015,2061)],
                                                     axis=1),axis=0)
    bm_cap_constr = bm_cap_constr.multiply(pd.concat([hsys["fraction"].to_frame(name=str(i))
                                        for i in range(2015,2061)],
                                       axis=1),axis=0).dropna()
   
    # FIXME: delete, now done below for all cap constraints
    # # load residual capacity for past years, assume constant into future,
    # # and set constraint to at least residual cap
    # rescap = pd.read_csv(snakemake.input.path_res_caps_ht,
    #                                             index_col=["LSOA11CD","PROPERTY_TYPE",
    #                                           "TECHNOLOGY"])
    
    # rescap = rescap[rescap.index.get_level_values("TECHNOLOGY").str.contains("BMBO")]
    # rescap = rescap.reindex_like(bm_cap_constr).fillna(method="ffill",axis=1)
    # rescap = rescap.fillna(0)
    
    # #bm_cap_constr[bm_cap_constr<rescap] = rescap
    # bm_cap_constr = bm_cap_constr.where(cond=~(bm_cap_constr<rescap),other=rescap)
    
#%%
# heritage constraint
  
    careas = pd.read_csv(snakemake.input.path_conserv_areas,
                              index_col="LSOA11CD")
    
    # load peak heat demands/technology capacity
    hpeaks = pd.read_csv(snakemake.input.path_peakc_lsoa,
                            index_col=["LSOA11CD",
                                       "PROPERTY_TYPE",
                                       "TECHNOLOGY"]) 
    if snakemake.params.dic["scen_hh_disagg"].startswith("T"):
            hpeaks.index = pd.MultiIndex.from_tuples([ (x[0],x[1],x[2]+"DDDE"+x[1].split("|")[1])
                                                  if x[1].startswith("Detached") else
                                                  (x[0],x[1],x[2]+"DDSD"+x[1].split("|")[1])
                                                  if x[1].startswith("Semi-detached") else
                                                  (x[0],x[1],x[2]+"DDTE"+x[1].split("|")[1])
                                                  if x[1].startswith("Terraced") else
                                                  (x[0],x[1],x[2]+"DDFL"+x[1].split("|")[1])
                                                  if x[1].startswith("Flats") else
                                                  (x[0],x[1],x[2]+"DNDO00")
                                                  if x[1]=="Non-domestic" else 0 for x
                                                  in hpeaks.index])
    else:              
        hpeaks.index = pd.MultiIndex.from_tuples([ (x[0],x[1],x[2]+"DDDE00")
                                                         if x[1].startswith("Detached") else
                                                         (x[0],x[1],x[2]+"DDSD00")
                                                         if x[1].startswith("Semi-detached") else
                                                         (x[0],x[1],x[2]+"DDTE00")
                                                         if x[1].startswith("Terraced") else
                                                         (x[0],x[1],x[2]+"DDFL00")
                                                         if x[1].startswith("Flats") else
                                                         (x[0],x[1],x[2]+"DNDO00")
                                                         if x[1].startswith("Non-domestic") else 0 for x
                                                         in hpeaks.index])
    hpeaks.index.names = ["LSOA11CD","PROPERTY_TYPE","TECHNOLOGY"]   

    # load suitability fractions
    raw_data = pd.read_csv(snakemake.input.path_set_htech)
    
    hsf = raw_data.loc[raw_data["VARIABLE"]=="HeritageSuitableFraction",
                       ["TECHNOLOGY","VALUE"]]
    hsf = hsf.assign(TECHNOLOGY=hsf['TECHNOLOGY'].str.split(',')).explode('TECHNOLOGY')
    hsf = hsf.set_index("TECHNOLOGY")
    
    # calculate fraction that is unsuitable
    hsf = 1-hsf
    # explode technologies if household disaggregated
    if snakemake.params.dic["scen_hh_disagg"]=="T":
        hsf = utils.explode_hh(hsf,
                                  ["OO","RP","RS"],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        
        hsf = utils.explode_hh(hsf,
                                  [t+str(i) for t in ["O","P","S"]
                                   for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )
        
    constr = hpeaks.multiply(pd.concat([hsf["VALUE"].to_frame(name=str(i))
                                        for i in range(2015,2061)],
                                       axis=1),axis=0)
    constr = constr.multiply(pd.concat([careas["cafrac"].to_frame(name=str(i))
                                        for i in range(2015,2061)],
                                       axis=1),axis=0)
    
    h_cap_constr = constr.dropna()
    
    # cap_constr = cap_constr.droplevel("PROPERTY_TYPE")
    # cap_constr = utils.groupby_LAD(cap_constr).sum()
    
    # cap_constr.columns.name = "YEAR"
    # cap_constr = cap_constr.stack()
  
    # # rearrange dataframe
    # cap_constr = cap_constr.reset_index()
    # cap_constr["UNIT"] = "GW"
    # cap_constr["VARIABLE"] = "TotalAnnualMaxCapacity"
    # cap_constr = cap_constr.rename(columns={"LAD21CD":"REGION",
    #                                         0:"VALUE"})
    # cap_constr = cap_constr.set_index([col for col 
    #                                    in cap_constr.columns
    #                                    if col!="VALUE"])
    
    # cap_constr.to_csv()
    

   
    
    # buildings_split = np.array_split(buildings,256)

    # buildings_split = [[buildings_split[i],carea,"union"]
    #                 for i in range(len(buildings_split))]

    # pool = mp.Pool(min(len(buildings_split),32),maxtasksperchild=2)
    # results = pool.starmap(gpd.overlay,buildings_split)
    # pool.terminate()

    # lsoa_road_length = pd.concat(results)

 

#%%
# space constraint

    stock = pd.read_csv(snakemake.input.path_dw_stock_con)
    stock["CONSTRAINTS"] = stock["CONSTRAINTS"].fillna(False)
    stock["CONSTRAINTS"] = stock["CONSTRAINTS"].astype(bool)
    stock = stock.set_index(["LSOA11CD","PROPERTY_TYPE","CONSTRAINTS"])
    
    stock["scfrac"] = (stock.xs(True,level=2)["NUMBER_OF_PROPERTIES"]
                       /stock.groupby(["LSOA11CD",
                                       "PROPERTY_TYPE"]).sum()["NUMBER_OF_PROPERTIES"])
    stock = stock.dropna().xs(True,level=2)
    
    
    scsf = raw_data.loc[raw_data["VARIABLE"]=="SpaceConstrainedSuitableFraction",
                       ["TECHNOLOGY","VALUE"]]
    scsf= scsf.assign(TECHNOLOGY=scsf['TECHNOLOGY'].str.split(',')).explode('TECHNOLOGY')
    scsf = scsf.set_index("TECHNOLOGY")
    
    # calc fraction that are unsuitable
    scsf = 1 - scsf
    
    # explode technologies if household disaggregated
    if snakemake.params.dic["scen_hh_disagg"]=="T":
        scsf = utils.explode_hh(scsf,
                                  ["OO","RP","RS"],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        
        scsf = utils.explode_hh(scsf,
                                  [t+str(i) for t in ["O","P","S"]
                                   for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )

    constr = hpeaks.multiply(pd.concat([scsf["VALUE"].to_frame(name=str(i))
                                        for i in range(2015,2061)],
                                       axis=1),axis=0)
    constr = constr.multiply(pd.concat([stock["scfrac"].to_frame(name=str(i))
                                        for i in range(2015,2061)],
                                       axis=1),axis=0)
    
    # applying (1-conservation area) fraction to avoid double counting
    constr = constr.multiply(pd.concat([(1-careas["cafrac"]).to_frame(name=str(i))
                                        for i in range(2015,2061)],
                                       axis=1),axis=0)
    
    s_cap_constr = constr.dropna()
    
    cap_constr = pd.concat([s_cap_constr.reset_index(),
                           h_cap_constr.reset_index()])

    
    # rename techs if household disaggregated
    # if snakemake.params.dic["scen_hh_disagg"]=="T":
    #     cap_constr["TECHNOLOGY"] = (cap_constr["TECHNOLOGY"].str[:-2]
    #                                 +cap_constr["PROPERTY_TYPE"].str.split("|")[1])
        
    cap_constr = cap_constr.groupby(["LSOA11CD",
                                    "PROPERTY_TYPE",
                                    "TECHNOLOGY"]).sum()
    # calculate capacity limit
    cap_constr = (hpeaks - cap_constr).dropna()
    
    cap_constr = pd.concat([cap_constr.reset_index(),
                           bm_cap_constr.reset_index()])
    
    cap_constr = cap_constr.set_index(["LSOA11CD",
                                    "PROPERTY_TYPE",
                                    "TECHNOLOGY"])
    
    # process
    cap_constr = cap_constr.reset_index()
    cap_constr = cap_constr.drop("PROPERTY_TYPE",axis=1)
    cap_constr = cap_constr.set_index(["LSOA11CD","TECHNOLOGY"])
    slagg = pd.read_csv(snakemake.input.path_sublocal_agg,
                            index_col=["LSOA11CD"]).squeeze()
    
    cap_constr = utils.groupby_LAD(cap_constr,disagg=slagg).sum()
    
    # load residual capacity for past years, assume constant into future,
    # and set constraint to at least residual cap
    rescap = pd.read_csv(snakemake.input.path_res_caps_ht)
    rescap = rescap.drop("PROPERTY_TYPE",axis=1)
    rescap = rescap.set_index(["LSOA11CD","TECHNOLOGY"])
    rescap = utils.groupby_LAD(rescap,disagg=slagg).sum()
    rescap = rescap.reindex_like(cap_constr).fillna(method="ffill",axis=1)
    rescap = rescap.fillna(0)
    
    if (cap_constr<rescap).any().any():
        
        logger.info("Some constraints for heating technologies are contradicted"
                    " by residual capacities and were adjusted based on them.")
        logger.info((cap_constr<rescap).any(axis=1)[(cap_constr<rescap).any(axis=1)==True])
        
        
    cap_constr = cap_constr.where(cond=~(cap_constr<rescap),other=rescap)
    
    
    cap_constr.columns.name = "YEAR"
    cap_constr = cap_constr.stack()
    
 
    
    # a=cap_constr.merge(rescap,right_index=True,left_index=True,how="left")
    # print((a["2020_x"]<a["2020_y"]).describe())
    # print(a)
    

    
    # rearrange dataframe
    cap_constr = cap_constr.reset_index()
    cap_constr["UNIT"] = "GW"
    cap_constr["VARIABLE"] = "TotalAnnualMaxCapacity"
    cap_constr = cap_constr.rename(columns={slagg.name:"REGION",
                                            0:"VALUE"})
    cap_constr = cap_constr.set_index([col for col 
                                       in cap_constr.columns
                                       if col!="VALUE"])
    # save to file
    # unit: GW
    cap_constr.to_csv(snakemake.output.path_heat_tech_con)

 



 
#%%
# deployment constraint (if applied)
    
    # add heat technology deployment constraint if triggered
    sets = pd.read_csv(snakemake.input.path_set_sets)
    techl = sets.loc[(sets["CATEGORY"]=="building_heat_supply")&
                     (sets["SET"]=="TECHNOLOGY"),
                     ["VALUE","PROPERTY_TYPE"]]
    techl.columns = ["TECHNOLOGY","PROPERTY_TYPE"]
    years = sets.loc[(sets["SET"]=="YEAR"),
                     "VALUE"].to_frame()
    slagg = pd.read_csv(snakemake.input.path_sublocal_agg,
                            index_col=["LSOA11CD"]).squeeze()
    
    if snakemake.params.dic["scen_hh_disagg"]=="T":
        techl = utils.explode_hh(techl,
                              ["OO","RP","RS"],
                              col="TECHNOLOGY")
        techl.loc[techl["TECHNOLOGY"].str.endswith(("OO","RS","RP")),
                  "PROPERTY_TYPE"] = (techl.loc[techl["TECHNOLOGY"].str.endswith(("OO","RS","RP")),
                                                "PROPERTY_TYPE"]
                                                +"|"+
                                                techl.loc[techl["TECHNOLOGY"].str.endswith(("OO","RS","RP")),
                                                          "TECHNOLOGY"].str[-2:])
                                                
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        hh = [t+str(i) for t in ["O","P","S"]
              for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)]
       
        techl = utils.explode_hh(techl,
                              hh,
                              col="TECHNOLOGY")
        techl.loc[techl["TECHNOLOGY"].str.endswith(tuple(hh)),
                  "PROPERTY_TYPE"] = (techl.loc[techl["TECHNOLOGY"].str.endswith(tuple(hh)),
                                                "PROPERTY_TYPE"]
                                                +"|"+
                                                techl.loc[techl["TECHNOLOGY"].str.endswith(tuple(hh)),
                                                          "TECHNOLOGY"].str[-2:])
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )
        
    # FIXME: properly check for right parameter syntax and raise error otherwise
    tag = pd.DataFrame(columns=["TYPE","REGION","TECHNOLOGY",
                                "MODE_OF_OPERATION",
                                "UDC","YEAR","VALUE"])
    limit = pd.DataFrame(columns=["REGION","UDC","YEAR","VALUE"])
    UDCTag = pd.DataFrame(columns=["REGION","UDC","VALUE"])
        
        
    if snakemake.params.dic["scen_htd_con"]=="-":                             
        pass
    else:
        # read technology capacities (only top 125 rows with average values)
        techcaps = pd.read_csv(snakemake.input.path_dw_peakd, 
                               index_col=["LSOA11CD","PROPERTY_TYPE","TECHNOLOGY"],
                               nrows=125)
        
        # load demand and property numbers
        dw_anndem = (pd.read_csv(snakemake.input.path_dw_ann_dem_sh,
                                 index_col=["LSOA11CD","PROPERTY_TYPE"])
                     + pd.read_csv(snakemake.input.path_dw_ann_dem_hw,
                                   index_col=["LSOA11CD","PROPERTY_TYPE"])
                     + pd.read_csv(snakemake.input.path_dw_ann_dem_sh_nb,
                                   index_col=["LSOA11CD","PROPERTY_TYPE"]))
        nd_anndem = (pd.read_csv(snakemake.input.path_nd_ann_dem_sh,
                      index_col=["LSOA11CD","PROPERTY_TYPE"])
                     + pd.read_csv(snakemake.input.path_nd_ann_dem_hw,
                                   index_col=["LSOA11CD","PROPERTY_TYPE"])
                     + pd.read_csv(snakemake.input.path_nd_ann_dem_sh_nb,
                                   index_col=["LSOA11CD","PROPERTY_TYPE"]))
        nd_anndem = nd_anndem.reset_index()
        nd_anndem.loc[:,"PROPERTY_TYPE"] = "Non-domestic"
        nd_anndem = nd_anndem.groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
                
        dwel = (pd.read_csv(snakemake.input.path_dw_stock_ex,
                              index_col=["LSOA11CD","PROPERTY_TYPE"])
                   +pd.read_csv(snakemake.input.path_dw_stock_nb,
                              index_col=["LSOA11CD","PROPERTY_TYPE"]))
        nd_stock =  (pd.read_csv(snakemake.input.path_nd_stock_ex,
                                   index_col=["LSOA11CD","PROPERTY_TYPE"])    
                        +pd.read_csv(snakemake.input.path_nd_stock_nb,
                                   index_col=["LSOA11CD","PROPERTY_TYPE"]))
        ndf = pd.read_csv(snakemake.input.path_set_nd_prop)
        ndf = ndf.loc[ndf["VARIABLE"]=="NonDomesticAverageFloorarea","VALUE"].squeeze()
        
        dwel = dwel.groupby("PROPERTY_TYPE").sum()
        nd_stock["PROPERTY_TYPE"] = "Non-domestic"
        nd_stock = nd_stock.set_index("PROPERTY_TYPE")
        nd_stock = nd_stock.groupby(["PROPERTY_TYPE"]).sum()
        nd_stock = nd_stock/ndf
        
        stock = pd.concat([dwel,nd_stock])
        stock = stock.mean(axis=1)/stock.mean(axis=1).sum()
        stock.name = "VALUE"
        techf = techl.merge(stock.to_frame(), on="PROPERTY_TYPE", how="left")
        techf = techf.loc[:,["TECHNOLOGY","VALUE"]] 
        
        hsys = pd.read_csv(snakemake.input.path_res_caps_ht_f, index_col=["LSOA11CD",
                                                               "PROPERTY_TYPE",
                                                               "TECHNOLOGY"]) 
        # hsys.xs("NGBO",level="TECHNOLOGY")
        
        cons = snakemake.params.dic["scen_htd_con"].split("|")
        cons = [c.replace("HPDNAH","HPDN-DN-ci-u-:*-60000y2028a190000y2035a190000y2060") for c in cons]
        cons = [c.replace("HPDNA","HPDN-DN-ci-u-:*-20000y2023a20000y2030a60000y2050a60000y2060") for c in cons]
        cons = cons + [c.replace("ciw","ci")  for c in cons if c.split("-")[2]=="ciw"]
        # these are also used below for the hydrogen boiler constraint
        al = ";".join(utils.get_entity_lookup(["LAD"])["LAD23CD"].to_list())
        uk = ",".join(utils.get_entity_lookup(["LAD"])["LAD23CD"].to_list())
        ic = "S12000014;S12000040;S12000036;E06000002;E06000004;E06000003;E08000011;E08000013;E08000014;E08000015;E08000012;E06000013;E06000012;E06000010;E06000011;E06000045;E07000086;E07000087;E07000088;E06000044;W06000009;W06000011;W06000012;W06000014;W06000015;W06000022;W06000019"
        cons = [c.replace("-IC-","-"+ic+"-") for c in cons]
        cons = [c.replace("-AL-","-"+al+"-") for c in cons]
        cons = [c.replace("-UK-","-"+uk+"-") for c in cons]
        cons = [c.split("-")[0:4]+[r]+[c.split("-")[5]] for c in cons for r in c.split("-")[4].split(";")]
        cons = ["-".join(c) for c in cons]
        for ic, c in enumerate(cons):
            logger.info("Working on constraint " + str(ic))
            tech,sec,ty,ul,rg,vys = c.split("-")
            techs = pd.DataFrame([t for t in techl["TECHNOLOGY"]
                                  if tech in t],columns=["TECHNOLOGY"])
            techs["UDC"] = "HTC"+str(ic)
            # add region, explode, and disaggregate to sublocal agg if needed
            techs["LAD23CD"] = rg
            techs = techs.assign(LAD23CD=techs["LAD23CD"].str.split(
                                        ',')).explode('LAD23CD')
            techs = techs.set_index(["LAD23CD","TECHNOLOGY","UDC"])
            if rg != ':*':
                techs["VALUE"] = 0
                techs = utils.groupby_LAD(techs,disagg=slagg).mean()
                techs = techs.drop("VALUE",axis=1)
            techs = techs.reset_index()
            techs = techs.rename(columns={slagg.name:"REGION",
                                          "LAD23CD":"REGION"})
            
            if ty.startswith("ci"):
                techs["TYPE"] = "NewCapacity"
                techs["MODE_OF_OPERATION"] =""
                
                if ty == "ciw":
                    t = techcaps.xs("AVERAGE").mean(axis=1)
                    t = t.reset_index().merge(stock,
                                              on="PROPERTY_TYPE",
                                              how="left")
                    t["TECH"] = t.loc[:,"TECHNOLOGY"].str[0:4]
                    t["VALUE"] = t["VALUE"]*t[0]
                    t = t.loc[t["TECHNOLOGY"].str.len()==10]
                    
                    t.loc[:,"VALUE"] = t.loc[:,["TECH",
                                                  "VALUE"]].groupby("TECH").transform("sum")
                    t = t.loc[:,["TECHNOLOGY","VALUE"]]
                    t.loc[:,"VALUE"] = 10**6/t.loc[:,"VALUE"]
                    techs = (techs.merge(t,
                                         how="left",on="TECHNOLOGY"))
                    
                else:
                    # multiplier as number of installations per GW
                    techs = (techs.merge(10**6/techcaps.xs("AVERAGE")
                                         .mean(axis=1).to_frame(),
                                         how="left",on="TECHNOLOGY"))
                    techs = techs.rename(columns={0:"VALUE"})
            elif ty == "ct":
                techs["TYPE"] = "TotalCapacity"
                techs["MODE_OF_OPERATION"]=""
                techs = (techs.merge(10**6/techcaps.xs("AVERAGE")
                                     .mean(axis=1).to_frame(),
                                     how="left",on="TECHNOLOGY"))
                techs = techs.rename(columns={0:"VALUE"})
            elif ty.startswith("a"):
                techs["TYPE"] = "Activity"
                techs["MODE_OF_OPERATION"] = ",".join([str(c) 
                                                       for c in list(ty[1:])])
                techs["VALUE"] = 1
                
            tag = pd.concat([tag,techs])
            
            # region in UDCTag does not matter, this will be replaced by
            # the multi-scale framework
            if ul == "u":
                UDCTag.loc[len(UDCTag)] = ["UK","HTC"+str(ic),0]
            elif ul == "l":
                UDCTag.loc[len(UDCTag)] = ["UK","HTC"+str(ic),1]
               
            lim = years.copy()
            lim["UDC"] = "HTC"+str(ic)
            lim = lim.rename(columns={"VALUE":"YEAR"})
            lim["VALUE"] = np.nan
            vys = vys.split("a")
            for  vy in vys:
                value, year = vy.split("y")
                lim.loc[lim["YEAR"]==year,"VALUE"]=value
            lim["VALUE"] = lim["VALUE"].astype(float)
            # interpolate
            if ul == "u":
                lim["VALUE"] = lim["VALUE"].interpolate(limit_area="inside")
            elif ul == "l":
                lim["VALUE"] = lim["VALUE"].interpolate(limit_area="inside")
            lim = lim.set_index(["YEAR","UDC"]).dropna()
            
            if ty.startswith("ci"):  
                # techcap = techcaps.loc[slice(None),slice(None),reftech]
                # m = {"D":"Domestic","N":"Non-domestic","DN":"AVERAGE"}
                # techcap = techcap.xs(m[sec],level=1).loc["AVERAGE"]
                # techcap.name="VALUE"
                # techcap = techcap.to_frame()
                # techcap.index.name="YEAR"
                # lim["VALUE"] = lim["VALUE"]*techcap["VALUE"]/10**6
                lim = lim.reset_index()
                limit = pd.concat([limit,lim])
            elif ty.startswith("a"):
                # take correct demand
                if sec == "D":
                    anndem = dw_anndem.copy()
                elif sec == "N":
                    anndem = nd_anndem.copy()
                elif sec == "DN":
                    anndem = pd.concat([dw_anndem,nd_anndem])
                    
                anndem = anndem.mul(hsys.xs("NGBO",level="TECHNOLOGY")["fraction"],axis=0)
                # for i in range(2015,2061):
                #     anndem[str(i)] = anndem[str(i)]*hsys.xs("NGBO",level="TECHNOLOGY")["fraction"]                 
                    
                if rg != ":*":
                    anndem = utils.groupby_LAD(anndem,disagg=slagg).sum()
                    anndem = anndem.groupby(slagg.name).sum()
                    anndem = anndem.loc[anndem.index.get_level_values(
                                        slagg.name).isin(techs["REGION"].unique())]
                
                # calculate total demand, convert to TJ, and multiply
                tdem = anndem.sum() / 1000
                tdem.index.name="YEAR"
                lim["VALUE"] = lim["VALUE"]*tdem
                lim = lim.reset_index()
                limit = pd.concat([limit,lim])
                
        tag["YEAR"] = ":*"
        tag = tag.assign(MODE_OF_OPERATION=tag["MODE_OF_OPERATION"].str.split(
                                    ',')).explode('MODE_OF_OPERATION')

        tag = tag[["TYPE","REGION","TECHNOLOGY","MODE_OF_OPERATION",
                   "UDC","YEAR","VALUE"]]
        tag = tag.dropna()
        limit["REGION"] = "UK"
        limit = limit.reset_index()[["REGION","UDC","YEAR","VALUE"]]
        
    # save to file
    # unit: -
    UDCTag.to_csv(snakemake.output.path_heat_tech_udc_tag,index=False)
        
    # save to file
    # unit: -
    tag.to_csv(snakemake.output.path_heat_tech_udc_tagtech,index=False)
    
    # save to file
    # unit: GW
    limit.to_csv(snakemake.output.path_heat_tech_udc_invcon,index=False)
    
#%%
# constraint for historical years

    # constrain any heat tech addition during historical years where residual
    # capacities are set
    sets = pd.read_csv(snakemake.input.path_set_sets)
    techs = sets.loc[(sets["CATEGORY"]=="building_heat_supply"),
                     "VALUE"].to_frame()
    years = sets.loc[(sets["SET"]=="YEAR"),
                     "VALUE"].to_frame()   
    
    # FIXME: now done below only for base year, but check if necessary
    # | allowing ELRE to balance, could potentially be removed
    # techs = techs[~techs["VALUE"].str.contains("ELRE")]
    
    if snakemake.params.dic["scen_hh_disagg"]=="T":
        techs = utils.explode_hh(techs,
                              ["OO","RP","RS"],
                              col="VALUE")
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        techs = utils.explode_hh(techs,
                                  [t+str(i) for t in ["O","P","S"]
                                   for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                  col="VALUE"
                                  )
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )
        
    
    inv_constr = techs
    inv_constr["YEAR"] = ",".join(years["VALUE"].to_list())
    inv_constr["REGION"] = ":*"
    inv_constr["VARIABLE"] = "TotalAnnualMaxCapacityInvestment"
    inv_constr = inv_constr.rename(columns={"VALUE":"TECHNOLOGY"})
    inv_constr["VALUE"] = -1

    inv_constr = inv_constr.assign(YEAR=inv_constr['YEAR'].str.split(',')).explode('YEAR')
    
    rescap = pd.read_csv(snakemake.input.path_res_caps_ht)
    inv_constr.loc[inv_constr["YEAR"].isin(list(rescap.columns[3:])),
                   "VALUE"] = 0
    # FIXME: see FIXME above
    inv_constr.loc[(inv_constr["YEAR"]=="2015")&
                    inv_constr["TECHNOLOGY"].str.startswith("OIBO"),
                    "VALUE"] = -1
    # inv_constr.loc[(inv_constr["YEAR"]=="2021")&
    #                 inv_constr["TECHNOLOGY"].str.startswith("OIBO"),
    #                 "VALUE"] = -1

    inv_constr = inv_constr.set_index([col for col 
                                       in inv_constr.columns
                                       if col!="VALUE"])
    
    # save to file
    # unit: GW
    inv_constr.to_csv(snakemake.output.path_heat_tech_invcon)
