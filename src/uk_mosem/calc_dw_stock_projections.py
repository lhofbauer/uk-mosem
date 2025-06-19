"""
Script deriving domestic building stock projections


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

from zipfile import ZipFile

import pandas as pd
import numpy as np

import utils

import logging

logger = logging.getLogger(__name__)
utils.adjust_logger()
 

if __name__ == "__main__":
#%%
    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_dw_stock_projections",
                                         **df_params.iloc[43].to_dict())    
  
    # load area lookup (LAD23CD)
    lsoas = utils.get_entity_lookup(["LSOA","LAD"])

    # load postcode lookup
    pcd = utils.get_entity_lookup(levels=["PCD","LSOA"])
    pcds = pcd[["PCD","LSOA11CD"]].set_index("PCD")["LSOA11CD"] 
    # load EPCs and arrange data (incl. aggregation with respect to LSOAs)
    logger.info("Loading EPC data")
    # set columns from DECs/EPCs to be retained
    cols = ["POSTCODE","BUILDING_REFERENCE_NUMBER",
            "INSPECTION_DATE","PROPERTY_TYPE", "BUILT_FORM","FLOOR_AREA",
            "CURRENT_ENERGY_RATING",
            #"OPERATIONAL_RATING_BAND", # additional params residual cap.:
            "HOTWATER_DESCRIPTION",
            "MAINHEAT_DESCRIPTION",
            "MAIN_HEATING_FUEL",
            "SECONDHEAT_DESCRIPTION",
            # additional parameters for efficiency improvements calc:
            "GLAZED_TYPE",
            "FLOOR_DESCRIPTION",
            "WINDOWS_DESCRIPTION",
            "WALLS_DESCRIPTION",
            "ROOF_DESCRIPTION",
            # additional parameters:
            "NUMBER_HABITABLE_ROOMS",
            "TENURE"
            ]
        
    # set columns which will be grouped by, i.e., LSOAs and categories
    cols_g = ["LSOA11CD",
              #"OPERATIONAL_RATING_BAND",
              "PROPERTY_TYPE", "BUILT_FORM",# additional params residual cap.:
              "HOTWATER_DESCRIPTION",
              "MAINHEAT_DESCRIPTION",
              "MAIN_HEATING_FUEL",
              "CURRENT_ENERGY_RATING",
              "SECONDHEAT_DESCRIPTION",# additional parameters for eff. improv. calc:
              "GLAZED_TYPE",
              "FLOOR_DESCRIPTION",
              "WINDOWS_DESCRIPTION",
              "WALLS_DESCRIPTION",
              "ROOF_DESCRIPTION",
              "SPACE_CONSTRAINED",
              "TENURE"]
        
    # go through ZIP files and load data
    stock_c = dict()
    stock_c = list()
    for zp in [snakemake.input.path_dw_stock_epc_ew,
               snakemake.input.path_dw_stock_epc_s]:
        
        zip_file = ZipFile(zp)
        
        # skip first row for the Scottish EPC file
        sr=None
        if zp==snakemake.input.path_dw_stock_epc_s:
            sr=1
        
        for file in zip_file.infolist():
            if file.filename.endswith(("certificates.csv",
                                       "_0721.csv")):
                data = pd.read_csv(zip_file.open(file.filename),skiprows=sr,
                                   encoding="ISO-8859-1")
                # rename columns for them to match across zip files
                data = data.rename(columns={"MAIN_FUEL":"MAIN_HEATING_FUEL",
                                            "TOTAL_FLOOR_AREA":"FLOOR_AREA",
                                            #"CURRENT_ENERGY_RATING":"OPERATIONAL_RATING_BAND",
                                            "WALL_DESCRIPTION":"WALLS_DESCRIPTION"})
                # only use the latest certificate for each building, discard
                # the rest
                data = data[cols].sort_values("INSPECTION_DATE")
                data = data[~data[["BUILDING_REFERENCE_NUMBER"]].duplicated(keep="last")]
                # calculate floor area per habitable room
                data.loc[:,"SPACE_CONSTRAINED"] = (data["FLOOR_AREA"]
                                                   /data["NUMBER_HABITABLE_ROOMS"])<16
                data = data.drop("NUMBER_HABITABLE_ROOMS",axis=1)

                # aggregate across LSOA while adding number of properties
                data.loc[:,"NUMBER_OF_PROPERTIES"] = 1
                data = data.drop(["INSPECTION_DATE",
                                  "BUILDING_REFERENCE_NUMBER"],
                                 axis=1)

                

                data["POSTCODE"] = data["POSTCODE"].map(pcds)

                data = data.rename(columns={"POSTCODE":"LSOA11CD"})

                # groupby, this also drops rows with (erroneous) postcodes that 
                # could not be matched
                data = data.groupby(cols_g).sum()
                
                stock_c.append(data)
                # if zp not in stock_c.keys():
                #     stock_c[zp] = data
                # else:
                #     stock_c[zp] = pd.concat([stock_c[zp],data])

    # concatenate data and merge with LSOAs
    # dw_stock = pd.concat(stock_c.values())
    dw_stock = pd.concat(stock_c)
    dw_stock = dw_stock.reset_index()
    dw_stock = dw_stock.groupby(cols_g).sum()
    dw_stock = lsoas["LSOA11CD"].to_frame().merge(right=dw_stock.reset_index(),
                                                  on="LSOA11CD", how="left")
    
    del stock_c
    logger.info("Processing EPC data")
    # process property type and built form parameter to categories used
    # by model
    dw_stock.loc[:,["PROPERTY_TYPE",
                    "BUILT_FORM"]] = dw_stock[["PROPERTY_TYPE",
                                               "BUILT_FORM"]].replace(
                                            {"Bungalow":"Semi-detached",
                                            "Flat":"Flats",
                                            "Maisonette":"Flats",
                                            "Park home":"House",
                                            "Mid-Terrace":"Terraced",
                                            "End-Terrace":"Terraced",
                                            "Semi-Detached":"Semi-detached",
                                            "Enclosed End-Terrace":"Terraced",
                                            "Enclosed Mid-Terrace":"Terraced",
                                            "NO DATA!":"Terraced"})
                                                        
    dw_stock["PROPERTY_TYPE"] = [y if x=="House" else x
                                 for x, y in zip(dw_stock["PROPERTY_TYPE"],
                                                 dw_stock["BUILT_FORM"])]
    dw_stock = dw_stock.drop(["BUILT_FORM"],axis=1)
    
    # process/clean property characteristics parameter 
    #dw_stock.TENURE.value_counts().to_csv("tenure_counts.csv")
    dw_stock[["TENURE","NUMBER_OF_PROPERTIES"]].groupby("TENURE").sum().to_csv("tenure.csv")
    
    # load lookup tables to process EPC data
    mhlu = pd.read_csv(snakemake.input.path_dw_mh_lu,index_col="MAINHEAT_DESCRIPTION",
                        usecols=[0,1])
    hwlu = pd.read_csv(snakemake.input.path_dw_hw_lu,index_col="HOTWATER_DESCRIPTION",
                        usecols=[0,1])
    shlu = pd.read_csv(snakemake.input.path_dw_sh_lu,index_col="SECONDHEAT_DESCRIPTION",
                        usecols=[0,1])
    
    dw_stock.loc[:,"MAINHEAT_DESCRIPTION"] = dw_stock["MAINHEAT_DESCRIPTION"].map(mhlu["MATCH"].to_dict())
    dw_stock.loc[:,"SECONDHEAT_DESCRIPTION"] = dw_stock["SECONDHEAT_DESCRIPTION"].map(shlu["MATCH"].to_dict())
    ht_lu = {"ASHP":"ASHP",
            "GSHP":"GSHP",
            "WSHP":"GSHP",
            "gas":"NGBO",
            "gCHP":"NGBO",
            "biomass":"BMBO",
            "electric_storage":"ELST",
            "electric":"ELRE",
            "oil":"OIBO",
            "oilCHP":"OIBO",
            "LPG":"OIBO",
            "coal":"BMBO",
            "electric_boiler":"ELRE",
            "communal":"HIUM-C",
            "gas_communal":"HIUM-G",
            "oil_communal":"HIUM-B",
            "LPG_communal":"HIUM-B",
            "biomass_communal":"HIUM-B",
            "WH_communal":"HIUM-G",
            "HP_communal":"HIUM-H",
            "coal_communal":"HIUM-B",
            "waste_communal":"HIUM-B",
            "geothermal_communal":"HIUM-H",
            "bCHP_communal":"HIUM-B",
            "gCHP_communal":"HIUM-C",
            "wCHP_communal":"HIUM-B",
            "CHP_communal":"HIUM-C"}
    
    # extract data on heat distribution systems within buildings
    dw_stock["HEATDIST_DESCRIPTION"] = dw_stock["MAINHEAT_DESCRIPTION"].str.split("-",expand=True)[1]
    dw_stock["MAINHEAT_DESCRIPTION"] = dw_stock["MAINHEAT_DESCRIPTION"].str.split("-",expand=True)[0]
    dw_stock["SECONDHEAT_DESCRIPTION"] = dw_stock["SECONDHEAT_DESCRIPTION"].str.split("-",expand=True)[0]
    

    
    # rename heating techs
    dw_stock["MAINHEAT_DESCRIPTION"] = dw_stock["MAINHEAT_DESCRIPTION"].replace(ht_lu)
    dw_stock["SECONDHEAT_DESCRIPTION"] = dw_stock["SECONDHEAT_DESCRIPTION"].replace(ht_lu)
    dw_stock["HEATDIST_DESCRIPTION"] = dw_stock["HEATDIST_DESCRIPTION"].replace({
                                        "underfloor":"RAUP",
                                        "wet":"WDIS",
                                        "air":"ADIS",
                                        "room":"NONE"})
    
    # remove errors with gas boilers on islands without gas grid
    il = ["S12000027","S12000023","E06000053"]
    illsoas = lsoas[lsoas["LAD23CD"].isin(il)]
    # set to oil boiler instead
    dw_stock.loc[dw_stock["LSOA11CD"].isin(illsoas["LSOA11CD"])&(
                 dw_stock["MAINHEAT_DESCRIPTION"]=="NGBO"),
                 "MAINHEAT_DESCRIPTION"] = "OIBO"
    
    # reassign heating systems assumed not available in flats
    # GSHP now assumed possible
    # dw_stock.loc[(dw_stock["PROPERTY_TYPE"]=="Flats")&(
    #              dw_stock["MAINHEAT_DESCRIPTION"]=="GSHP"),
    #              "MAINHEAT_DESCRIPTION"] = "ASHP"
    
    dw_stock.loc[(dw_stock["PROPERTY_TYPE"]=="Flats")&(
                 dw_stock["MAINHEAT_DESCRIPTION"]=="BMBO"),
                 "MAINHEAT_DESCRIPTION"] = "OIBO"

    # process hotwater data
    dw_stock["HOTWATER_DESCRIPTION"] = dw_stock["HOTWATER_DESCRIPTION"].map(hwlu["MATCH"].to_dict())
    # ignore second hot water system, i.e., solar thermal
    dw_stock["HOTWATER_DESCRIPTION"] = dw_stock["HOTWATER_DESCRIPTION"].str.split("|",expand=True)[0]
    # rename hotwater techs
    dw_stock["HOTWATER_DESCRIPTION"] = dw_stock["HOTWATER_DESCRIPTION"].replace({
                                        "electric_immersion":"ELIM",
                                        "main":"main",
                                        "secondary":"secondary",
                                        "gas":"NGBO",
                                        "communal":"HIUM",
                                        "electric_instant":"ELPU",
                                        "oil":"OIBO",
                                        "solid":"BMBO",
                                        "HP":"ASHP"})
    
    # process data with respect to heat system combination, currently not
    # really necessary because HW, second heat system is disregarded but could 
    # potentially later effect choice of overall SHW system if more SWH with 
    # same SH but different HW components are modelled (not the case currently)
    
    # dw_stock["HOTWATER_DESCRIPTION"] = [mh+"_SHW" if hw=="main" else
    #                                     mh+"_SHW" if ((hw=="secondary") and
    #                                                   (mh == sh)) else
    #                                     sh if hw=="secondary" else
    #                                     hw                                        
    #                                     for mh, sh, hw in zip(dw_stock["MAINHEAT_DESCRIPTION"],
    #                                              dw_stock["SECONDHEAT_DESCRIPTION"],
    #                                              dw_stock["HOTWATER_DESCRIPTION"])]    
    
    # set all heat distributions systems to wet if heating system is modelled
    # as wet and is not set as low temperature wet system
    dw_stock.loc[dw_stock["MAINHEAT_DESCRIPTION"].isin(["NGBO","OIBO","BMBO",
                                                       "ASHP","GSHP","HIUM-G",
                                                       "HIUM-C","HIUM-B",
                                                       "HIUM-H"])&(
                 dw_stock["HEATDIST_DESCRIPTION"] !="RAUP"),"HEATDIST_DESCRIPTION"] = "WDIS"
    # set all heat distributions systems to low temperature wet if a HP is
    # installed
    dw_stock.loc[dw_stock["MAINHEAT_DESCRIPTION"].isin(["ASHP","GSHP"]),
                 "HEATDIST_DESCRIPTION"] = "RAUP"      
    # set all remaining air distributions systems (currently not modelled)
    # to room
    dw_stock.loc[dw_stock["HEATDIST_DESCRIPTION"]=="ADIS",
                 "HEATDIST_DESCRIPTION"] = "NONE"      


    
    
    # process/clean property characteristics parameter    
    

    # load lookup tables to process EPC data
    walu = pd.read_csv(snakemake.input.path_dw_wa_lu,index_col="WALLS_DESCRIPTION",
                        usecols=[0,1])
    rolu = pd.read_csv(snakemake.input.path_dw_ro_lu,index_col="ROOF_DESCRIPTION",
                        usecols=[0,1])
    fllu = pd.read_csv(snakemake.input.path_dw_fl_lu,index_col="FLOOR_DESCRIPTION",
                        usecols=[0,1])
    wilu = pd.read_csv(snakemake.input.path_dw_wi_lu,
                       index_col=["WINDOWS_DESCRIPTION","GLAZED_TYPE"],
                        usecols=[0,1,2])
    wilu.index= wilu.index.to_flat_index()
    
    
    dw_stock.loc[:,"WALLS_DESCRIPTION"] = dw_stock["WALLS_DESCRIPTION"].map(walu["MATCH"].to_dict())
    dw_stock.loc[:,"ROOF_DESCRIPTION"] = dw_stock["ROOF_DESCRIPTION"].map(rolu["MATCH"].to_dict())
    dw_stock.loc[:,"FLOOR_DESCRIPTION"] = dw_stock["FLOOR_DESCRIPTION"].map(fllu["MATCH"].to_dict())
    
    dw_stock["WINDOWS_DESCRIPTION"] = list(zip(dw_stock["WINDOWS_DESCRIPTION"], 
                                                      dw_stock["GLAZED_TYPE"]))
    
    dw_stock.loc[:,"WINDOWS_DESCRIPTION"] = dw_stock["WINDOWS_DESCRIPTION"].map(wilu["MATCH"].to_dict())
    
    
    # Code below is not used anymore but left to potentially be used for 
    # (pre-) processing of future EPC downloads to generate lookups
    
    # the lookup created from this includes an empty/NaN row mapped to a value
    # - should be updated or row deleted manually if used again
    
    # delete all except first entry (before first separator |), process
    # dw_stock["WALLS_DESCRIPTION"] = dw_stock["WALLS_DESCRIPTION"].str.split("|",
    #                                     expand=True)[0].str.strip().fillna("")
    # dw_stock["WALLS_DESCRIPTION"] = ["cavity_insulated" 
    #                                   if (("Cavity" in e
    #                                       or "Waliau ceudod" in e)
    #                                       and ("filled" in e 
    #                                       or "insulated" in e
    #                                       or "with" in e
    #                                       or "lenwi" in e
    #                                       or "wedi?u" in e))
    #                                   else "cavity_uninsulated" 
    #                                   if ("Cavity" in e)
    #                                   else "solid_insulated"  
    #                                   if (("Solid" in e
    #                                       or "Sandstone"  in e
    #                                       or "Granite"  in e
    #                                       or "System"  in e
    #                                       or "system"  in e
    #                                       or "Cob"  in e
    #                                       or "Park home"  in e
    #                                       or "Timber"  in e
    #                                       or "Gwenithfaen" in e
    #                                       or "Briciau solet" in e
    #                                       or "Tywodfaen" in e
    #                                       or "bren") 
    #                                       and ("insulated" in e
    #                                           or "with" in e
    #                                           or "wediâu" in e
    #                                           or "wediGÃÃu" in e
    #                                           or "gydag" in e
    #                                           or "wedi?u" in e))
    #                                   else "solid_uninsulated"
    #                                   if ("Solid" in e
    #                                       or "Sandstone"  in e
    #                                       or "Granite"  in e
    #                                       or "System"  in e
    #                                       or "system"  in e
    #                                       or "Cob"  in e
    #                                       or "Park home"  in e
    #                                       or "Timber"  in e
    #                                       or "Gwenithfaen" in e
    #                                       or "Briciau solet" in e
    #                                       or "Tywodfaen" in e
    #                                       or "bren") 
    #                                   else np.nan
    #                                   for e in dw_stock["WALLS_DESCRIPTION"]]

    # delete all except first entry (before first separator |), process   
    # dw_stock["ROOF_DESCRIPTION"] = dw_stock["ROOF_DESCRIPTION"].str.split("|",
    #                                                 expand=True)[0].str.strip().fillna("")    
    # dw_stock["ROOF_DESCRIPTION"] = ["pitched_fully_insulated" 
    #                                   if (("Pitched" in e
    #                                       or "Roof room" in e
    #                                       or "Thatched" in e
    #                                       or "Ar oleddf" in e
    #                                       or "To gwellt" in e
    #                                       or "Ystafell" in e)
    #                                       and ("200 mm" in e 
    #                                       or "200mm" in e 
    #                                       or "250+ mm" in e
    #                                       or "250mm" in e
    #                                       or "270 mm" in e
    #                                       or "300 mm" in e
    #                                       or "300mm" in e
    #                                       or "300+ mm" in e
    #                                       or "350 mm" in e
    #                                       or "400 mm" in e
    #                                       or "400+ mm" in e
    #                                       or "400+  mm" in e
    #                                       or ("insulated" in e and "rafters" not in e)
    #                                       or "Thatched" in e
    #                                       or "thatched" in e
    #                                       or "To gwellt" in e
    #                                       or ("wediâi" in e and "trawstiau" not in e)
    #                                       or ("wediGÃÃi" in e and "trawstiau" not in e)
    #                                       or "wedi?i" in e))
    #                                   else "pitched_partially_insulated" 
    #                                   if (("Pitched" in e
    #                                       or "Roof room" in e
    #                                       or "Thatched" in e
    #                                       or "Ar oleddf" in e
    #                                       or "To gwellt" in e
    #                                       or "Ystafell" in e)
    #                                       and ("100 mm" in e
    #                                     or "150 mm" in e))
    #                                   else "pitched_uninsulated"
    #                                   if ("Pitched" in e
    #                                       or "Roof room" in e
    #                                       or "Ar oleddf" in e
    #                                       or "Ystafell" in e)
    #                                   else "flat_insulated"  
    #                                   if (("Flat" in e
    #                                       or "Yn wastad"  in e) 
    #                                       and ("insulated" in e
    #                                           or "wediâi" in e
    #                                           or "wediGÃÃu" in e
    #                                           or "wedi?u" in e))
    #                                   else "flat_uninsulated"
    #                                   if ("Flat" in e
    #                                       or "Yn wastad"  in e)
    #                                   else "no_roof"
    #                                   if ("dwelling above" in e
    #                                   or "annedd arall uwchben" in e)
    #                                   else np.nan
    #                                   for e in dw_stock["ROOF_DESCRIPTION"]] 

    # delete all except first entry (before first separator |), process
    # dw_stock["FLOOR_DESCRIPTION"] = dw_stock["FLOOR_DESCRIPTION"].str.split("|",
    #                                                 expand=True)[0].str.strip().fillna("") 
    # dw_stock["FLOOR_DESCRIPTION"] = ["suspended_insulated" 
    #                                   if (("Suspended" in e
    #                                       #or "To unheated space" in e
    #                                       #or "To external air" in e
    #                                       #or "I ofod heb ei wresogi" in e
    #                                       or "Crog" in e)
    #                                       and ("insulated" in e
    #                                       or "wediâu" in e
    #                                       or "wediGÃÃu" in e
    #                                       or "wedi?u" in e))
    #                                   else "suspended_uninsulated" 
    #                                   if ("Suspended" in e)
    #                                   else "solid_insulated"  
    #                                   if (("Solid" in e
    #                                       or "Solet"  in e) 
    #                                       and ("insulated" in e
    #                                           or "wediâu" in e
    #                                           or "wediGÃÃu" in e
    #                                           or "wedi?u" in e
    #                                           or "wedi'i" in e))
    #                                   else "solid_uninsulated"
    #                                   if ("Solid" in e
    #                                       or "Solet"  in e) 
    #                                   else "no_floor"
    #                                   if ("below" in e
    #                                       or "anheddiad arall islaw" in e
    #                                       or "eiddo arall islaw" in e)
    #                                   else np.nan
    #                                   for e in dw_stock["FLOOR_DESCRIPTION"]]    

    # # delete all except first entry (before first separator |), process
    # dw_stock["WINDOWS_DESCRIPTION"] = dw_stock["WINDOWS_DESCRIPTION"].str.replace('\s{2,}',
    #                                                                           ' ').fillna("") 
    # dw_stock["WINDOWS_DESCRIPTION"] = ["single_glazed"  
    #                                     if ("ingle" in e
    #                                         or (("Some" in e or "some" in e) 
    #                                             and ("ultiple" in e
    #                                             or "ouble" in e
    #                                             or "econdary" in e))
    #                                             or "Gwydrau sengl" in e
    #                                             or "Rhai gwydrau dwbl" in e)
    #                                     else "double_glazed_new" 
    #                                     if (("ouble glaz" in e
    #                                       or "ultiple glaz" in e
    #                                       or "complex" in e
    #                                       or "Some triple glazing" in e
    #                                       or "Gwydrau dwbl llawn" in e
    #                                       or "Gwydrau dwbl gan mwyaf" in e
    #                                       or "Gwydrau dwbl rhannol" in e
    #                                       or "Gwydrau lluosog ym mhobman" in e)
    #                                     and "after 2002" in t)
    #                                   else "double_glazed_old" 
    #                                     if ("ouble glaz" in e
    #                                       or "ultiple glaz" in e
    #                                       or "complex" in e
    #                                       or "Some triple glazing" in e
    #                                       or "Gwydrau dwbl llawn" in e
    #                                       or "Gwydrau dwbl gan mwyaf" in e
    #                                       or "Gwydrau dwbl rhannol" in e
    #                                       or "Gwydrau lluosog ym mhobman" in e)
    #                                   else "secondary_glazing"
    #                                   if ("econdary glaz" in e
    #                                       or "Gwydrau eilaidd llawn" in e)
    #                                   else "triple_glazed"
    #                                   if ("riple glaz" in e
    #                                       or "High performance glazing" in e
    #                                       or "Gwydrau triphlyg llawn" in e)
    #                                   else np.nan
    #                                   for e, t in zip(dw_stock["WINDOWS_DESCRIPTION"],
    #                                               dw_stock["GLAZED_TYPE"])] 
       
    # process tenure and add to property type if required
    # dw_stock.loc[ :,
    #              ["TENURE","NUMBER_OF_PROPERTIES"]].groupby("TENURE",
    #                                                     dropna=False).sum().to_csv("tenure_scaled.csv") 
    # pd.Series(dw_stock["TENURE"].unique()).to_csv("tenure_types.csv")
    
    dw_stock.loc[dw_stock["TENURE"].isna(),"TENURE"] = "OO"
    dw_stock.loc[dw_stock["TENURE"].str.contains("private"),"TENURE"] = "RP"
    dw_stock.loc[dw_stock["TENURE"].str.contains("social"),"TENURE"] = "RS"
    dw_stock.loc[(~dw_stock["TENURE"].str.contains("RP"))
                  & (~dw_stock["TENURE"].str.contains("RS")),"TENURE"] = "OO"
    
    
    # dw_stock[["TENURE","NUMBER_OF_PROPERTIES"]].groupby("TENURE").sum().to_csv("tenure_proc_scaled.csv")
    
    
    # rearrange data, without deleting nan rows
    dw_stock = dw_stock.drop(["GLAZED_TYPE","HOTWATER_DESCRIPTION",
                              "SECONDHEAT_DESCRIPTION",
                              "MAIN_HEATING_FUEL"],axis=1)
    dw_stock = dw_stock.groupby(["LSOA11CD","PROPERTY_TYPE",
                                 "CURRENT_ENERGY_RATING",
                                 "HEATDIST_DESCRIPTION",
                                 "MAINHEAT_DESCRIPTION",
                                   "FLOOR_DESCRIPTION",
                                   "WALLS_DESCRIPTION",
                                   "WINDOWS_DESCRIPTION",
                                   "ROOF_DESCRIPTION",
                                   "SPACE_CONSTRAINED",
                                   "TENURE"],dropna=False).sum()
    dw_stock = dw_stock.reset_index()
    
    # dw_stock = dw_stock.iloc[::5, :]
    
    # load dwelling type numbers per LSOA/MSOA from VOA data
    logger.info("Loading VOA data")
    dwel_ew = dict()
    dwel_s = dict()
    dwel = dict()
    
    # load totals (without property type) for Scotland to derive
    # numbers for year 2018-2022
    dwel_s_totals = pd.read_csv(snakemake.input.paths_dw_stock_s_totals,
                            encoding="utf-8", skiprows=7,
                            usecols=[0,14,15,16,17,18,19])
    dwel_s_totals = dwel_s_totals.rename(columns = {
    "http://purl.org/linked-data/sdmx/2009/dimension#refArea":"LSOA11CD"})
    dwel_s_totals["LSOA11CD"] = dwel_s_totals["LSOA11CD"].str.split("/").str[-1]
    # convert to numeric and fill blanks/etc. with 0
    dwel_s_totals.iloc[:,1:] = dwel_s_totals.iloc[:,1:].apply(pd.to_numeric,
                                    errors='coerce').fillna(0).astype(int)
    dwel_s_totals = dwel_s_totals.replace(",","", regex=True)
    
    for i,y in enumerate(range(2015,2015+len(snakemake.input.paths_dw_stock_ew))):
        dwel_ew[i] = pd.read_csv(snakemake.input.paths_dw_stock_ew[i],
                                 encoding="ISO-8859-1")

        if y<=2017:
            dwel_s[i] = pd.read_csv(snakemake.input.paths_dw_stock_s[i],
                                    encoding="utf-8", skiprows=7,
                                    usecols=[0,2,3,4,5,6,7])
            
        # for 2018 onward no Scottish data (with property type) available,
        # scale 2017 data based on 2018-2022 data on totals
        else:
            dwel_s[i] = pd.read_csv(snakemake.input.paths_dw_stock_s[2],
                                    encoding="utf-8", skiprows=7,
                                    usecols=[0,2,3,4,5,6,7])

        # create LSOACD column
        dwel_s[i] = dwel_s[i].rename(columns = {
        "http://purl.org/linked-data/sdmx/2009/dimension#refArea":"LSOA11CD"})
        dwel_s[i]["LSOA11CD"] = dwel_s[i]["LSOA11CD"].str.split("/").str[-1]
    
        # convert to numeric and fill blanks/etc. with 0
        dwel_s[i].iloc[:,1:] = dwel_s[i].iloc[:,1:].apply(pd.to_numeric,
                                        errors='coerce').fillna(0).astype(int)
        dwel_ew[i] = dwel_ew[i].replace(",","", regex=True)
        dwel_ew[i].iloc[:,5:] = dwel_ew[i].iloc[:,5:].apply(pd.to_numeric,
                                     errors='coerce').fillna(0).astype(int)
    
        # pick rows with all council bands
        dwel_ew[i] = dwel_ew[i].rename(columns={"band":"BAND"})
        dwel_ew[i] = dwel_ew[i][dwel_ew[i]["BAND"]=="All"]
    
        # add bungalow number to semi-detached houses (similar characterstics)
        dwel_ew[i]["house_semi_total"] = (dwel_ew[i]["house_semi_total"]
                                              +dwel_ew[i]["bungalow_total"])
        
        # deducting mobile homes (including boats) from total
        dwel_ew[i]["all_properties"] = (dwel_ew[i]["all_properties"]
                                    -dwel_ew[i]["caravan_houseboat_mobilehome"])
        
        # discard irrelevant columns and rename rest
        dwel_ew[i] = dwel_ew[i][["ecode","flat_mais_total",
                                 "house_terraced_total","house_semi_total",
                                 "house_detached_total","all_properties"]]  
    
        dwel_ew[i] = dwel_ew[i].rename(columns = {"ecode":"LSOA11CD",
                                        "flat_mais_total":"Flats",
                                        "house_terraced_total":"Terraced",
                                        "house_semi_total":"Semi-detached",
                                        "house_detached_total":"Detached",
                                        "all_properties":"All"})
        
        dwel_s[i] = dwel_s[i][["LSOA11CD","All","Flats","Terraced",
                            "Semi-Detached","Detached"]]
        dwel_s[i] = dwel_s[i].rename(columns = {"Semi-Detached":"Semi-detached"})
        
        # if 2018 or later, adjust based on development of total number
        # this removes also some non current datazone geographies
        if y>2017:
            dwel_s_totals["VALUE"] = (dwel_s_totals[str(y)]/dwel_s_totals["2017"]).fillna(0)
            mult = dwel_s_totals.set_index("LSOA11CD")["VALUE"]
            dwel_s[i] = dwel_s[i].set_index("LSOA11CD")
            dwel_s[i] = dwel_s[i].multiply(mult,axis=0).dropna()
            dwel_s[i] = dwel_s[i].reset_index()
            
        # concat EW and S data and add second column level
        dwel[i] = pd.concat([dwel_ew[i],dwel_s[i]], sort=True)
        dwel[i] = dwel[i].set_index("LSOA11CD")
        
        # scale archetype numbers based on all properties number (and thus
        # including unknows and others)
        dwel[i]["Sum"] = dwel[i][["Flats","Terraced",
                                  "Semi-detached","Detached"]].sum(axis=1)
        
        # replace 0 values for "All" with sum values to avoid downscaling in
        # case of missing data (alternatively, the following line could be
        # uncommented to only scale up, but not down based on the "All" column)
        dwel[i].loc[dwel[i]["All"]==0,"All"] = dwel[i]["Sum"]
        # dwel[i]["All"] = dwel[i][["All","Sum"]].max(axis=1)
        dwel[i] = dwel[i].multiply(dwel[i]["All"],axis=0).divide(
                                            dwel[i]["Sum"],axis=0).fillna(0)
        
        dwel[i] = dwel[i][["Flats","Terraced",
                            "Semi-detached","Detached"]]
        

    # pull data together in one dataframe
    dwel = pd.concat([df.stack() for df in dwel.values()], axis=1)
    dwel.columns = [str(i) for i in range(2015,2023)]
    
    ews_total = dwel.loc[["K04000001","S92000003"],:].groupby(level=1).sum()
    ews_total.index.name = "PROPERTY_TYPE"

    #dwel = dwel.drop(["All"])
    dwel = lsoas["LSOA11CD"].to_frame().merge(right=dwel.reset_index(1),
                                              on="LSOA11CD", how="left")
    dwel = dwel.rename(columns={"level_1":"PROPERTY_TYPE"})
    dwel = dwel.set_index(["LSOA11CD",
                           "PROPERTY_TYPE"])
    
    # scale up data to match GB-wide totals (for each property type, which has
    # before been scaled to total as for all LSOAs)
    dwel = (dwel
            *(ews_total
            /dwel.groupby("PROPERTY_TYPE").sum()))
    
    # round while keeping overall sum per year the same
    dwel = dwel.apply(utils.round_col, axis=0)
    
    
    # merge VOA with EPC data
    logger.info("Scaling EPC data")
    
    if ((snakemake.params.dic["scen_hh_disagg"]=="T") or
        (snakemake.params.dic["scen_hh_disagg"].startswith("TI"))):
        # if tenure is used, scale based on local VOA property type counts and 
        # national data on tenure across property types data
        
        dw_stock = dwel.stack().to_frame().reset_index().merge(right=dw_stock,
                                                                on=["LSOA11CD",
                                                                    "PROPERTY_TYPE"],
                                                                how="left")

        # logger.info(dwel.stack().sum())
        # logger.info(dw_stock["NUMBER_OF_PROPERTIES"].sum())
        
        dw_stock.loc[dw_stock["NUMBER_OF_PROPERTIES"].isna(),
                     "NUMBER_OF_PROPERTIES"] = dw_stock.loc[dw_stock["NUMBER_OF_PROPERTIES"].isna(),0]/3
        # assume mix if no certificate for LSOA and property type
        dw_stock.loc[dw_stock["TENURE"].isna(),"TENURE"] = "OO,RP,RS"
        dw_stock = dw_stock.assign(TENURE=dw_stock["TENURE"].str.split(',')).explode("TENURE")
        
        # dw_stock = dw_stock.set_index([c for c in dw_stock.columns
        #                                if c != "NUMBER_OF_PROPERTIES"])
        dw_stock = dw_stock.merge(right=dw_stock.groupby(["LSOA11CD",
                                "PROPERTY_TYPE","level_2"]).sum()["NUMBER_OF_PROPERTIES"],
                                  on=["LSOA11CD","PROPERTY_TYPE","level_2"], how="left",
                                  suffixes=("","_LPsum"))
        dw_stock["NUMBER_OF_PROPERTIES"] = (dw_stock["NUMBER_OF_PROPERTIES"]
                                            *dw_stock[0]
                                            /dw_stock["NUMBER_OF_PROPERTIES_LPsum"])
                                            
        # load tenure data      
        ehs = pd.read_excel(snakemake.input.path_dw_ehs,
                            sheet_name="AT1_6",
                            usecols="B:J",
                            skiprows=5,
                            nrows=24,
                            #header=None,
                            index_col=0)
        ehs = ehs.loc[~ehs.index.isna()]
        ehs = ehs.rename(index={"all terraced houses":"Terraced",
                                 "semi-detached house":"Semi-detached",
                                 "detached house": "Detached",
                                 "bungalow":"Semi-detached",
                                 "converted flat":"Flats",
                                 "purpose built flat, low rise":"Flats",
                                 "purpose built flat, high rise":"Flats"})
        ehs = ehs.rename(columns={"owner occupied":"OO",
                                 "private rented":"RP",
                                 "all social sector": "RS"})
        
        ehs = ehs.iloc[13:]
        ehs = ehs.loc[:,["OO","RP","RS"]]
        ehs.index.name="PROPERTY_TYPE"
        ehs.columns.name = "TENURE"
        ehs = ehs.groupby("PROPERTY_TYPE").sum()
        ehs = ehs.stack()
        ehs.name = "VALUE"
        ehs = ehs/ehs.groupby("PROPERTY_TYPE").sum()
        
        ten = ehs.reset_index()
        
        
   
        # scale percentage tenure data based on dwelling numbers
        ten = ten.merge(dwel.stack().reset_index().groupby(["PROPERTY_TYPE","level_2"])[0].sum().reset_index(),
                        on="PROPERTY_TYPE",
                        how="left")
        ten["TenTotal"] = ten["VALUE"] * ten[0]
        ten = ten.loc[:,["PROPERTY_TYPE","TENURE","level_2",
                         "TenTotal"]].set_index(["PROPERTY_TYPE","TENURE","level_2"])
        
        dw_stock = dw_stock.merge(right=ten,
                                  on=["PROPERTY_TYPE","TENURE","level_2"], how="left")
        # logger.info(dwel.stack().sum())
        # logger.info(dw_stock["NUMBER_OF_PROPERTIES"].sum())
                                    
        for i in range(15):
            
            if "NUMBER_OF_PROPERTIES_PTsum" in dw_stock.columns:
                dw_stock = dw_stock.drop("NUMBER_OF_PROPERTIES_PTsum",axis=1)
            # scale to match national tenure-property type numbers    
            dw_stock = dw_stock.merge(right=dw_stock.groupby(["level_2",
                                "PROPERTY_TYPE","TENURE"])["NUMBER_OF_PROPERTIES"].sum(),
                                  on=["level_2","PROPERTY_TYPE","TENURE"], how="left",
                                  suffixes=("","_PTsum"))
            
            dw_stock["NUMBER_OF_PROPERTIES"] = (dw_stock["NUMBER_OF_PROPERTIES"]
                                                *dw_stock["TenTotal"]
                                                /dw_stock["NUMBER_OF_PROPERTIES_PTsum"])
            
            # scale to match actual number of properties per type in each LSOA
            dw_stock = dw_stock.drop("NUMBER_OF_PROPERTIES_LPsum",axis=1)
            dw_stock = dw_stock.merge(right=dw_stock.groupby(["LSOA11CD",
                                    "PROPERTY_TYPE","level_2"]).sum()["NUMBER_OF_PROPERTIES"],
                                      on=["LSOA11CD","PROPERTY_TYPE","level_2"], how="left",
                                      suffixes=("","_LPsum"))
            
            dw_stock["NUMBER_OF_PROPERTIES"] = (dw_stock["NUMBER_OF_PROPERTIES"]
                                                *dw_stock[0]
                                                /dw_stock["NUMBER_OF_PROPERTIES_LPsum"])
            
            # calculate measure for convergence towards national statistics
            # (local alignment is ensured as scaled last)
            ch = (dw_stock.groupby(["PROPERTY_TYPE","TENURE","level_2"])["NUMBER_OF_PROPERTIES"].sum()
                  -ten["TenTotal"]).abs()
            chp = ch/ten["TenTotal"]
            # logger.info(ch.describe())
            # logger.info(dwel.stack().sum())
            # logger.info(dw_stock["NUMBER_OF_PROPERTIES"].sum()) 
            
            if ch.abs().max()<100:
                break
        logger.info("Accuracy of scaling (difference on national level):")
        logger.info(ch.describe())
        logger.info(chp.describe())
        dw_stock["NUMBER_OF_PROPERTIES"]=dw_stock["NUMBER_OF_PROPERTIES"].fillna(0)
        # to check   
        # (dw_stock.groupby(["LSOA11CD",
        #                     "PROPERTY_TYPE",
        #                     "level_2"])["NUMBER_OF_PROPERTIES"].sum()
        #   .xs("2015",level="level_2").fillna(100)
        #   -dwel["2015"]).abs().describe()
        

        utils.groupby_LAD(dw_stock[["LSOA11CD","PROPERTY_TYPE",
                                   "TENURE","level_2",
                                   "NUMBER_OF_PROPERTIES"]].set_index(["LSOA11CD",
                                                                      "PROPERTY_TYPE",
                                                                      "TENURE",
                                                                      "level_2"])).sum().to_csv("ladten.csv")
    
    elif snakemake.params.dic["scen_hh_disagg"]=="-":
        # if tenure is not used, scale based on VOA data on number
        # of property types
        
        dw_stock = dw_stock.merge(right=dw_stock.groupby(["LSOA11CD",
                                "PROPERTY_TYPE"]).sum()["NUMBER_OF_PROPERTIES"],
                                  on=["LSOA11CD","PROPERTY_TYPE"], how="left")
        dw_stock = dwel.stack().to_frame().reset_index().merge(right=dw_stock,
                                                                on=["LSOA11CD",
                                                                    "PROPERTY_TYPE"],
                                                                how="left")

    
        del pcd
        # del s
        
        # scale number of properties and floor space based on VOA data
        l = list(range(0, len(dw_stock), int(len(dw_stock)/1000)))+[len(dw_stock)]
        li = list()         
        for i in range(0,len(l)-1):
            li.extend([(z/y*x) if (z!=0 and y!=0 and pd.notna(y))
                       #else 0 if x==0
                       else z 
                       for x, y, z
                       in zip(dw_stock["NUMBER_OF_PROPERTIES_x"][l[i]:l[i+1]],
                              dw_stock["NUMBER_OF_PROPERTIES_y"][l[i]:l[i+1]],
                              dw_stock[0][l[i]:l[i+1]])])
        dw_stock["NUMBER_OF_PROPERTIES"] = li
        
        li = list()         
        for i in range(0,len(l)-1):
            li.extend([(z/y*x) if (z!=0 and y!=0 and pd.notna(y))
                       else 0 if x==0
                       else pd.NA 
                       for x, y, z
                       in zip(dw_stock["FLOOR_AREA"][l[i]:l[i+1]],
                              dw_stock["NUMBER_OF_PROPERTIES_y"][l[i]:l[i+1]],
                              dw_stock[0][l[i]:l[i+1]])] )
        dw_stock["FLOOR_AREA"] = li
        del li 
        
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
                )
                
    logger.info("Rearranging and saving detailed base year data")
    
    # rearrange data, drop unnecessary columns
    dw_stock = dw_stock.rename(columns={"level_2":"YEAR"})
    dw_stock = dw_stock[["LSOA11CD","PROPERTY_TYPE",
                         "YEAR","NUMBER_OF_PROPERTIES",
                         "CURRENT_ENERGY_RATING",
                         "HEATDIST_DESCRIPTION",
                         "MAINHEAT_DESCRIPTION",
                         "FLOOR_DESCRIPTION",
                         "WALLS_DESCRIPTION",
                         "WINDOWS_DESCRIPTION",
                         "ROOF_DESCRIPTION",
                         "SPACE_CONSTRAINED",
                         "TENURE"]]
    dw_stock.loc[:,["CURRENT_ENERGY_RATING",
                         "HEATDIST_DESCRIPTION",
                         "MAINHEAT_DESCRIPTION",
                         "FLOOR_DESCRIPTION",
                         "WALLS_DESCRIPTION",
                         "WINDOWS_DESCRIPTION",
                         "ROOF_DESCRIPTION",
                         "SPACE_CONSTRAINED"]] = dw_stock.loc[:,["CURRENT_ENERGY_RATING",
                         "HEATDIST_DESCRIPTION",
                         "MAINHEAT_DESCRIPTION",
                         "FLOOR_DESCRIPTION",
                         "WALLS_DESCRIPTION",
                         "WINDOWS_DESCRIPTION",
                         "ROOF_DESCRIPTION",
                         "SPACE_CONSTRAINED"]].astype(str)
    
  
    # cast as categories to reduce memory usage
    cat = ["YEAR",
           "CURRENT_ENERGY_RATING",
           "HEATDIST_DESCRIPTION",
           "MAINHEAT_DESCRIPTION",
           "FLOOR_DESCRIPTION",
           "WALLS_DESCRIPTION",
           "WINDOWS_DESCRIPTION",
           "ROOF_DESCRIPTION",
           "SPACE_CONSTRAINED",
           "TENURE"]
    for c in cat:
        dw_stock[c] = dw_stock[c].astype("category")
        
    dw_stock["NUMBER_OF_PROPERTIES"] = pd.to_numeric(dw_stock["NUMBER_OF_PROPERTIES"],
                                                     downcast="float")    

    # if required add tenure and income to property type
    if snakemake.params.dic["scen_hh_disagg"]=="T":
        dw_stock["PROPERTY_TYPE"] = (dw_stock["PROPERTY_TYPE"] +"|"
                                      +dw_stock["TENURE"])
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        
        # logger.info(dw_stock["NUMBER_OF_PROPERTIES"].sum())
        # logger.info(len(dw_stock["NUMBER_OF_PROPERTIES"]))

        # split out Wales and Scotland which are not disaggregated
        dw_stock_WS = dw_stock.loc[~dw_stock["LSOA11CD"].str.startswith("E")]
        dw_stock = dw_stock.loc[dw_stock["LSOA11CD"].str.startswith("E")]
                                  
        # get fraction of properties in income band that is in specific
        # tenure for each LSOA
        ehs = pd.read_excel(snakemake.input.path_ec_ghhi_ten,
                            sheet_name="AT1.3 ",
                            usecols="L:S",
                            skiprows=45,
                            nrows=5,
                            header=None,
                            index_col=0)
        ehs.columns = ["_","_","OO","RP","_","_","RS"]
        ehs = ehs.loc[:,["OO","RP","RS"]]
        ehs.index.name = "INCOME"
        ehs.columns.name = "TENURE"
        
        tcount = dw_stock.loc[:,
                              ["YEAR","TENURE",
                                "NUMBER_OF_PROPERTIES"]].groupby(["YEAR","TENURE"]).sum()
        
        ehs = ehs/100*tcount["NUMBER_OF_PROPERTIES"]

        # groupby income groups used in the model
        # FIXME: make flexible to adjust to chosen income bands
        ehs.index = ["0","2","3","3","3"]
        ehs.index.name = "INCOME"
        ehs = ehs.groupby("INCOME").sum()
        ehs.loc["1",:] = ehs.loc["0",:]/2
        ehs.loc["0",:] = ehs.loc["0",:]/2
        ehs = ehs.sort_index()
        ehs = ehs.stack([0,1])
        ehs.name = "ITtotal"
    
        #ehs = ehs.divide(ehs.groupby("YEAR",axis=1).sum(),axis=0)

        # get number of properties in each income band in each LSOA                                              
        ghhi = pd.read_csv(snakemake.input.path_ec_ghhi_dist,
                                     index_col="LSOA11CD")
        ghhi.columns = [str(i) for i in range(len(ghhi.columns))]
        #ghhi = ghhi.loc[ghhi.index.get_level_values("LSOA11CD").str.startswith("E")]
        ghhi.columns.name = "INCOME"
        tp = dw_stock.loc[:,
                              ["LSOA11CD","YEAR",
                               "NUMBER_OF_PROPERTIES"]].groupby(["LSOA11CD",
                                                                 "YEAR"]).sum()
        tp = tp["NUMBER_OF_PROPERTIES"].multiply(ghhi.stack(),axis=0)
        
        tp.name = "LItotal"

        
        # set entry totals
        # dw_stock = dw_stock.merge(right=dw_stock.groupby([c for c in dw_stock.columns
        #                                                   if ~c.startswith("NUMBER_OF")&
        #                                                   ~c.endswith("total")]).sum()["NUMBER_OF_PROPERTIES"],
        #                           on=[c for c in dw_stock.columns
        #                               if ~c.startswith("NUMBER_OF")&
        #                               ~c.endswith("total")],
        #                           how="left",
        #                           suffixes=("","_LPTXtotal"))
        dw_stock["LPTXtotal"] = dw_stock["NUMBER_OF_PROPERTIES"]
        
        # disaggregate dataframe based on income and set initial income profile
        dw_stock_WS.loc[:,"INCOME"] = ghhi.columns[-1]
        
        dw_stock = dw_stock.loc[dw_stock.index.repeat(len(ghhi.columns))]
        dw_stock.loc[:,"INCOME"] = [r for r in ghhi.columns] * int((len(dw_stock)/len(ghhi.columns)))

        # dw_stock = dw_stock.assign(INCOME=dw_stock["INCOME"].str.split(',')).explode("INCOME")
        
        # set number of properties for all tenures/property types across incomes
        # based on the LSOA income data, as initial values to be scaled below
        ghhi = ghhi.stack()
        ghhi.name = "GHHI"

        dw_stock = dw_stock.merge(ghhi.reset_index(),on=["LSOA11CD","INCOME"],
                                  how="left")

        dw_stock.loc[:,"NUMBER_OF_PROPERTIES"] = (dw_stock.loc[:,"NUMBER_OF_PROPERTIES"]
                                            *dw_stock.loc[:,"GHHI"])
        dw_stock = dw_stock.drop("GHHI",axis=1)

        # merge totals for scaling
        dw_stock = dw_stock.merge(right=ehs,
                                  on=["TENURE","INCOME","YEAR"], how="left")
        dw_stock = dw_stock.merge(right=tp,
                                  on=["LSOA11CD","INCOME","YEAR"], how="left")
        
      
        # logger.info(dw_stock["NUMBER_OF_PROPERTIES"].sum())
        # logger.info(len(dw_stock["NUMBER_OF_PROPERTIES"]))

        dw_stock = dw_stock.merge(right=dw_stock.groupby(["LSOA11CD",
                                "INCOME","YEAR"]).sum()["NUMBER_OF_PROPERTIES"],
                                  on=["LSOA11CD","INCOME","YEAR"], how="left",
                                  suffixes=("","_LIsum"))
        
#%%     # scale to match national income by tenure data, as well as LSOA
        # income split (while keeping LSOA property type and tenure
        # number constant)                                       
        for i in range(7):
            

            # if "NUMBER_OF_PROPERTIES_LPTXsum" in dw_stock.columns:
            #     dw_stock = dw_stock.drop("NUMBER_OF_PROPERTIES_LPTXsum",axis=1)
                
            # dw_stock = dw_stock.merge(right=dw_stock.groupby([c for c in dw_stock.columns
            #                                                   if  (not c.startswith("NUMBER_OF")) &
            #                                                    (not c.endswith("total"))
            #                                                    & (c!="INCOME")],
            #                                                  observed=True,
            #                                                  dropna=False)["NUMBER_OF_PROPERTIES"].sum(),
            #                           on=[c for c in dw_stock.columns
            #                             if  (not c.startswith("NUMBER_OF")) &
            #                              (not c.endswith("total"))
            #                              & (c!="INCOME")],
            #                           how="left",
            #                           suffixes=("","_LPTXsum"))
            
            # dw_stock["NUMBER_OF_PROPERTIES"] = (dw_stock["NUMBER_OF_PROPERTIES"]
            #                                     *dw_stock["LPTXtotal"]
            #                                     /dw_stock["NUMBER_OF_PROPERTIES_LPTXsum"]).fillna(0)
            
            if "NUMBER_OF_PROPERTIES_ITsum" in dw_stock.columns:
                dw_stock = dw_stock.drop("NUMBER_OF_PROPERTIES_ITsum",axis=1)
   
            dw_stock = dw_stock.merge(right=dw_stock.groupby(["TENURE",
                                                              "INCOME",
                                                              "YEAR"])["NUMBER_OF_PROPERTIES"].sum(),
                                      on=["TENURE","INCOME","YEAR"], how="left",
                                      suffixes=("","_ITsum"))

            dw_stock["NUMBER_OF_PROPERTIES"] = (dw_stock["NUMBER_OF_PROPERTIES"]
                                                *dw_stock["ITtotal"]
                                                /dw_stock["NUMBER_OF_PROPERTIES_ITsum"]
                                                *dw_stock["LItotal"]
                                                /dw_stock["NUMBER_OF_PROPERTIES_LIsum"]
                                                ).fillna(0)

            if "NUMBER_OF_PROPERTIES_LIsum" in dw_stock.columns:
                dw_stock = dw_stock.drop("NUMBER_OF_PROPERTIES_LIsum",axis=1)

            dw_stock = dw_stock.merge(right=dw_stock.groupby(["LSOA11CD",
                                    "INCOME","YEAR"]).sum()["NUMBER_OF_PROPERTIES"],
                                      on=["LSOA11CD","INCOME","YEAR"], how="left",
                                      suffixes=("","_LIsum"))

            dw_stock["NUMBER_OF_PROPERTIES"] = (dw_stock["NUMBER_OF_PROPERTIES"]
                                                *dw_stock["LItotal"]
                                                /dw_stock["NUMBER_OF_PROPERTIES_LIsum"]).fillna(0)
            
            
            if "NUMBER_OF_PROPERTIES_LPTXsum" in dw_stock.columns:
                dw_stock = dw_stock.drop("NUMBER_OF_PROPERTIES_LPTXsum",axis=1)
   
            dw_stock = dw_stock.merge(right=dw_stock.groupby([c for c in dw_stock.columns
                                                              if  (not c.startswith("NUMBER_OF")) &
                                                               (not c.endswith("total"))
                                                               & (c!="INCOME")],
                                                             observed=True,
                                                             dropna=False)["NUMBER_OF_PROPERTIES"].sum(),
                                      on=[c for c in dw_stock.columns
                                        if  (not c.startswith("NUMBER_OF")) &
                                         (not c.endswith("total"))
                                         & (c!="INCOME")],
                                      how="left",
                                      suffixes=("","_LPTXsum"))

            dw_stock["NUMBER_OF_PROPERTIES"] = (dw_stock["NUMBER_OF_PROPERTIES"]
                                                *dw_stock["LPTXtotal"]
                                                /dw_stock["NUMBER_OF_PROPERTIES_LPTXsum"]).fillna(0)
            
            
            # to check   
            # (dw_stock.groupby(["LSOA11CD",
            #                     "PROPERTY_TYPE",
            #                     "YEAR"])["NUMBER_OF_PROPERTIES"].sum()
            #   .xs("2015",level="YEAR").fillna(0)
            #   -dwel.loc[dwel.index.get_level_values("LSOA11CD").str.startswith("E"),"2015"]).abs().describe()
            # calculate measure for convergence towards national
            # (tenure & income) and local (income) statistics
            chtp = (dw_stock.groupby(["LSOA11CD","YEAR","INCOME"])["NUMBER_OF_PROPERTIES"].sum()
                  -tp).abs()
            chtpp = chtp/tp
            logger.info(chtp.describe())
            logger.info(chtpp.describe())
            chehs = (dw_stock.groupby(["INCOME","YEAR","TENURE"])["NUMBER_OF_PROPERTIES"].sum()
                  -ehs).abs()
            chehsp = chehs/ehs
            logger.info(chehs.describe())
            logger.info(chehsp.describe())

            # logger.info(dwel.stack().sum())
            # logger.info(dw_stock["NUMBER_OF_PROPERTIES"].sum()) 
            
            if (chtp.max()<10) & (chehs.max()<50000):
                break
        
        dw_stock = dw_stock.drop(["LPTXtotal","LItotal", "ITtotal",
                                  "NUMBER_OF_PROPERTIES_ITsum",
                                  "NUMBER_OF_PROPERTIES_LIsum",
                                  "NUMBER_OF_PROPERTIES_LPTXsum"],
                                 axis=1)

        dw_stock = pd.concat([dw_stock,
                              dw_stock_WS])
        dw_stock = dw_stock.reset_index(drop=True)

        dw_stock["PROPERTY_TYPE"] = (dw_stock["PROPERTY_TYPE"] +"|"
                                      +dw_stock["TENURE"].str[1]
                                      +dw_stock["INCOME"])

        dw_stock = dw_stock.drop("INCOME",axis=1)
        
        
        
    for c in cat:
        dw_stock[c] = dw_stock[c].replace("nan",np.nan)
   
        
    # get data for adjusting demand based on EPC rating
    det = dw_stock.drop(["MAINHEAT_DESCRIPTION",
                         "FLOOR_DESCRIPTION",
                        "WALLS_DESCRIPTION",
                        "WINDOWS_DESCRIPTION",
                        "ROOF_DESCRIPTION",
                        "HEATDIST_DESCRIPTION",
                        "SPACE_CONSTRAINED",
                        "TENURE"],axis=1)

    det = det[det["YEAR"]=="2019"]

    det = pd.concat([det.groupby(["LSOA11CD",
            "PROPERTY_TYPE"]+[d],dropna=False,observed=True).sum() 
                     for d in ["CURRENT_ENERGY_RATING"]])
    det.index = det.index.set_names("EPC_RATING",level=2)
    
    det.to_csv(snakemake.output.path_dw_stock_epc)
    
    # get data for DH residual caps
    det = dw_stock.drop(["CURRENT_ENERGY_RATING",
                         "FLOOR_DESCRIPTION",
                        "WALLS_DESCRIPTION",
                        "WINDOWS_DESCRIPTION",
                        "ROOF_DESCRIPTION",
                        "HEATDIST_DESCRIPTION",
                        "SPACE_CONSTRAINED",
                        "TENURE"],axis=1)
    det["MAINHEAT_DESCRIPTION"] = det["MAINHEAT_DESCRIPTION"].cat.add_categories(["MH_NA"])
    det = det[det["MAINHEAT_DESCRIPTION"].fillna("MH_NA").str.startswith("HIUM")]
    det = det[det["YEAR"]=="2022"]
    det = pd.concat([det.groupby(["LSOA11CD",
            "PROPERTY_TYPE"]+[d],dropna=False,observed=True).sum() 
                     for d in ["MAINHEAT_DESCRIPTION"]])
    det.index = det.index.set_names("DH_TYPE",level=2)
    
    #det = det.unstack(fill_value=0).stack()
    # save to file
    det.to_csv(snakemake.output.path_dw_stock_dh)     
    
    dw_stock["MAINHEAT_DESCRIPTION"] = dw_stock["MAINHEAT_DESCRIPTION"].astype("object")
    dw_stock.loc[dw_stock["MAINHEAT_DESCRIPTION"].fillna("MH_NA").str.startswith("HIUM"),
                 "MAINHEAT_DESCRIPTION"] = dw_stock.loc[dw_stock["MAINHEAT_DESCRIPTION"].fillna("MH_NA").str.startswith("HIUM"),
                              "MAINHEAT_DESCRIPTION"].str.split("-",expand=True)[0]
    dw_stock["MAINHEAT_DESCRIPTION"] = dw_stock["MAINHEAT_DESCRIPTION"].astype("category")
    
    # process dataframe for constraints, choosing year 
    # as most recent base year
    
    # rearrange data in a way that lists characteristics separately (while
    # losing information/detail)
    det = dw_stock.drop(["CURRENT_ENERGY_RATING",
                         "FLOOR_DESCRIPTION",
                        "WALLS_DESCRIPTION",
                        "WINDOWS_DESCRIPTION",
                        "ROOF_DESCRIPTION",
                        "MAINHEAT_DESCRIPTION",
                        "HEATDIST_DESCRIPTION",
                        "TENURE"],axis=1)
    det = det[det["YEAR"]=="2022"]
    det = pd.concat([det.groupby(["LSOA11CD",
            "PROPERTY_TYPE"]+[d],dropna=False,observed=True).sum() 
                     for d in ["SPACE_CONSTRAINED"]])
    det.index = det.index.set_names("CONSTRAINTS",level=2)
    
    #det = det.unstack(fill_value=0).stack()
    # save to file
    det.to_csv(snakemake.output.path_dw_stock_con) 

    # process dataframe for socioecomonic parameters, choosing year 
    # as most recent base year

    # rearrange data in a way that lists characteristics separately (while
    # losing information)
    det = dw_stock.drop(["CURRENT_ENERGY_RATING",
                         "FLOOR_DESCRIPTION",
                        "WALLS_DESCRIPTION",
                        "WINDOWS_DESCRIPTION",
                        "ROOF_DESCRIPTION",
                        "MAINHEAT_DESCRIPTION",
                        "HEATDIST_DESCRIPTION",
                        "SPACE_CONSTRAINED"],axis=1)
    det = det[det["YEAR"]=="2022"]
    det = pd.concat([det.groupby(["LSOA11CD",
            "PROPERTY_TYPE"]+[d],dropna=False,observed=True).sum()
                     for d in ["TENURE"]])
    det.index = det.index.set_names("TENURE",level=2)

    # save to file
    det.to_csv(snakemake.output.path_dw_stock_sec)         
    
    # process dataframe for residual heating system capacities, choosing year 
    # as most recent base year

    dw_stock["MAINHEAT_DESCRIPTION"] = dw_stock["MAINHEAT_DESCRIPTION"].cat.add_categories(["MH_"+c for c in dw_stock["MAINHEAT_DESCRIPTION"].cat.categories]+["MH_NA"])
    dw_stock["HEATDIST_DESCRIPTION"] = dw_stock["HEATDIST_DESCRIPTION"].cat.add_categories(["HD_"+c for c in dw_stock["HEATDIST_DESCRIPTION"].cat.categories]+["HD_NA"])

    dw_stock.loc[dw_stock["MAINHEAT_DESCRIPTION"].notna(),
                 "MAINHEAT_DESCRIPTION"]= "MH_" + dw_stock.loc[dw_stock["MAINHEAT_DESCRIPTION"].notna(),
                                                               "MAINHEAT_DESCRIPTION"].astype(str)

    dw_stock.loc[dw_stock["HEATDIST_DESCRIPTION"].notna(),
                 "HEATDIST_DESCRIPTION"]= "HD_" + dw_stock.loc[dw_stock["MAINHEAT_DESCRIPTION"].notna(),
                                                               "HEATDIST_DESCRIPTION"].astype(str)

    dw_stock.loc[:,"MAINHEAT_DESCRIPTION"] = dw_stock["MAINHEAT_DESCRIPTION"].fillna("MH_NA")

    dw_stock.loc[:,"HEATDIST_DESCRIPTION"] = dw_stock["HEATDIST_DESCRIPTION"].fillna("HD_NA")

    # rearrange data in a way that lists characteristics separately (while
    # losing information), groupby removes any potential NA values in group keys
    det = dw_stock.drop(["CURRENT_ENERGY_RATING",
                         "FLOOR_DESCRIPTION",
                        "WALLS_DESCRIPTION",
                        "WINDOWS_DESCRIPTION",
                        "ROOF_DESCRIPTION",
                        "SPACE_CONSTRAINED",
                        "TENURE"],axis=1)

    det = det[det["YEAR"]=="2022"]

    det = pd.concat([det.groupby(["LSOA11CD",
            "PROPERTY_TYPE"]+[d],dropna=False,observed=True).sum()
                     for d in ["MAINHEAT_DESCRIPTION",
                               "HEATDIST_DESCRIPTION"]])

    det.index = det.index.set_names("HEATING_SYSTEM",level=2)


    # save to file
    det.to_csv(snakemake.output.path_dw_stock_hde)  

    # process dataframe for building efficiency calculations, choose year based
    # on year of other building efficiency data
    
    # code here can be exchanged with alternative code below, if saving
    # of large data files is fine, and detailed data are necessary.
    
    # load stock data, including characteristics of wall, roof, etc.
    dw_stock["FLOOR_DESCRIPTION"] = dw_stock["FLOOR_DESCRIPTION"].cat.add_categories(["floor_"+c for c in dw_stock["FLOOR_DESCRIPTION"].cat.categories]+["floor_NA"])
    dw_stock["ROOF_DESCRIPTION"] = dw_stock["ROOF_DESCRIPTION"].cat.add_categories(["roof_"+c for c in dw_stock["ROOF_DESCRIPTION"].cat.categories]+["roof_NA"])
    dw_stock["WALLS_DESCRIPTION"] = dw_stock["WALLS_DESCRIPTION"].cat.add_categories(["wall_"+c for c in dw_stock["WALLS_DESCRIPTION"].cat.categories]+["wall_NA"])
    dw_stock["WINDOWS_DESCRIPTION"] = dw_stock["WINDOWS_DESCRIPTION"].cat.add_categories(["windows_"+c for c in dw_stock["WINDOWS_DESCRIPTION"].cat.categories]+["windows_NA"])

    dw_stock.loc[dw_stock["FLOOR_DESCRIPTION"].notna(),
                 "FLOOR_DESCRIPTION"]= "floor_" + dw_stock["FLOOR_DESCRIPTION"].astype(str)
    dw_stock.loc[dw_stock["WINDOWS_DESCRIPTION"].notna(),
                 "WINDOWS_DESCRIPTION"]= "windows_" + dw_stock["WINDOWS_DESCRIPTION"].astype(str)
    dw_stock.loc[dw_stock["ROOF_DESCRIPTION"].notna(),
                 "ROOF_DESCRIPTION"]= "roof_" + dw_stock["ROOF_DESCRIPTION"].astype(str)
    dw_stock.loc[dw_stock["WALLS_DESCRIPTION"].notna(),
                 "WALLS_DESCRIPTION"]= "wall_" + dw_stock["WALLS_DESCRIPTION"].astype(str)

    dw_stock["WALLS_DESCRIPTION"] = dw_stock["WALLS_DESCRIPTION"].fillna("wall_NA")
    dw_stock["ROOF_DESCRIPTION"] = dw_stock["ROOF_DESCRIPTION"].fillna("roof_NA")
    dw_stock["WINDOWS_DESCRIPTION"] = dw_stock["WINDOWS_DESCRIPTION"].fillna("windows_NA")
    dw_stock["FLOOR_DESCRIPTION"] = dw_stock["FLOOR_DESCRIPTION"].fillna("floor_NA")

    # rearrange data in way that lists characteristics separately (while
    # losing information)
    det = dw_stock.drop(["CURRENT_ENERGY_RATING",
                         "HEATDIST_DESCRIPTION",
                         "MAINHEAT_DESCRIPTION",
                         "SPACE_CONSTRAINED",
                         "TENURE"],axis=1)
    det = det[det["YEAR"]=="2019"]
    det = pd.concat([det.groupby(["LSOA11CD",
            "PROPERTY_TYPE"]+[d],dropna=False,observed=True).sum()
                     for d in ["FLOOR_DESCRIPTION",
                               "WINDOWS_DESCRIPTION",
                               "WALLS_DESCRIPTION",
                               "ROOF_DESCRIPTION"]])
    det.index = det.index.set_names("EFF",level=2)

    # alternative for detailed version:
    # det = dw_stock.drop("MAIN_HEATING_FUEL",axis=1)
    # det = det[det["YEAR"]=="2019"]
    # det = det.groupby(["LSOA11CD","PROPERTY_TYPE",
    #                                "FLOOR_DESCRIPTION",
    #                                "WALLS_DESCRIPTION",
    #                                "WINDOWS_DESCRIPTION",
    #                                "ROOF_DESCRIPTION"],dropna=False).sum()
    
    # save to file
    det.to_csv(snakemake.output.path_dw_stock_ede)
    
    # remove columns not necessary for remaining calculations
    logger.info("Calculating projections")
    
    dw_stock = dw_stock[["LSOA11CD","PROPERTY_TYPE",
                         "YEAR","NUMBER_OF_PROPERTIES"]]

    dw_stock.loc[:,["LSOA11CD","PROPERTY_TYPE",
                   "YEAR"]] =  dw_stock.loc[:,["LSOA11CD","PROPERTY_TYPE",
                                              "YEAR"]].astype("object")
                                              
    dw_stock = dw_stock.groupby(["LSOA11CD","PROPERTY_TYPE",
                                 "YEAR"]).sum()["NUMBER_OF_PROPERTIES"].unstack(["YEAR"])
    

    # load household projections (percentage change)
    hh_proj_pch_ls = pd.read_csv(snakemake.input.path_hh_proj,
                                 index_col="LSOA11CD")
    
    
    # load demolition data (E:LAD17CD, WS: LAD23CD - no recent changes)
    dwel_d_s = pd.read_excel(snakemake.input.path_dw_stock_d_s,
                              sheet_name="tsDemolitionsFinYr",skiprows=3,
                              nrows=33)
    dwel_d_e = pd.read_excel(snakemake.input.path_dw_stock_d_e,
                             sheet_name=None,skiprows=2,usecols=[3,15,22,23,32]) 
    dwel_d_w = pd.read_csv(snakemake.input.path_dw_stock_d_w,nrows=23) 
    
    for y,d in dwel_d_e.items():
        dwel_d_e[y].iloc[0,0]="LAD17CD"
        dwel_d_e[y] = dwel_d_e[y].rename(columns=dwel_d_e[y].iloc[0]).iloc[1:][
                                                    ["LAD17CD","Demolitions"]]
        dwel_d_e[y] = dwel_d_e[y].rename(columns={"Demolitions":y[0:4]})
        dwel_d_e[y] = dwel_d_e[y].dropna(axis=0,how="any")
        dwel_d_e[y] = dwel_d_e[y].set_index("LAD17CD")

    # fix change of code from first year to the rest for Stevenage,
    # Northumberland, Gateshead, and East Hertfordshire
    dwel_d_e["2012-13"] = dwel_d_e["2012-13"].rename(index={"E08000020":"E08000037",
                                                            "E06000048":"E06000057",
                                                            "E07000101":"E07000243",
                                                            "E07000097":"E07000242"})
    # update LAD codes
    for y,d in dwel_d_e.items(): 
        
        if y=="2019-20":
             dwel_d_e[y] = utils.update_LADCD(dwel_d_e[y],from_CD="LAD19CD")
        elif y=="2020-21":
             dwel_d_e[y] = utils.update_LADCD(dwel_d_e[y],from_CD="LAD20CD")
        else:
            dwel_d_e[y] = utils.update_LADCD(dwel_d_e[y],from_CD="LAD17CD")
    
    dwel_d_e = dwel_d_e["2012-13"].join(list(dwel_d_e.values())[1:])
    

    # process Scotland, Wales data, concat with Scotland, and further process
    dwel_d_s = dwel_d_s.set_index("Local Authority")
    dwel_d_s.index.name = "LAD23NM"
    dwel_d_s.columns = dwel_d_s.columns.str[0:4]
    
    
    dwel_d_w = dwel_d_w.set_index("Unnamed: 0")
    dwel_d_w.index.name = "LAD23NM"
    dwel_d_w.columns = dwel_d_w.columns.str[0:4]

    dwel_d_sw = pd.concat([dwel_d_s,dwel_d_w], sort=True)
    dwel_d_sw.index =  dwel_d_sw.index.str.replace("&","and")
    dwel_d_sw = dwel_d_sw.rename(index={"Edinburgh, City of": "City of Edinburgh",
                                        "Orkney":"Orkney Islands",
                                        "Shetland":"Shetland Islands",
                                        "Scottish Borders, The":"Scottish Borders"})
    
    # set LAD codes instead of names as index
    dwel_d_sw.index = dwel_d_sw.index.str.strip().map(lsoas[["LAD23NM",
                "LAD23CD"]].set_index("LAD23NM").drop_duplicates().iloc[:,0])
    dwel_d_sw.index.name = "LAD23CD"
    
    # concat all dataframes
    dwel_d_esw = pd.concat([dwel_d_e,dwel_d_sw], sort=True)
    
    # calculate average demolition rate per LAD/LSOA
    dwel_d = lsoas[["LAD23NM","LAD23CD"]].drop_duplicates()
    dwel_d = dwel_d.merge(right=dwel_d_esw, on="LAD23CD", how="left")
    dwel_d = dwel_d.set_index("LAD23CD")
    dwel_d["avg_rate"] = dwel_d.loc[:,"2012":"2018"].mean(axis=1)

    # calculate demolition rate per LSOA by assuming equal distribution across
    # LSOA and property types

    dw_stock_r = lsoas[["LSOA11CD","LAD23CD"]].merge(right=dw_stock.reset_index(),on="LSOA11CD",how="left")
    dwel_d_rate = dw_stock_r.loc[:,["LAD23CD","LSOA11CD",
                                    "PROPERTY_TYPE","2022"]].merge(
        right=dwel_d["avg_rate"], on="LAD23CD",
        how="left")

    dwel_d_rate = dwel_d_rate.merge(
        right=dwel_d_rate[["LAD23CD","2022"]].groupby(["LAD23CD"]).sum(),
        on="LAD23CD", how="left")
    
    dwel_d_rate["avg_rate"] = (dwel_d_rate["avg_rate"] * dwel_d_rate["2022_x"] /
                               dwel_d_rate["2022_y"])
    
    dwel_d_rate = dwel_d_rate.set_index(["LSOA11CD",
                                         "PROPERTY_TYPE"])["avg_rate"]
    

    # calculate dwelling stock for future years, separately for existing and
    # new

    

    #dw_stock = dw_stock.reset_index().set_index(["LSOA11CD","PROPERTY_TYPE","YEAR"]).unstack(["YEAR"])

    ndw_stock = dw_stock.copy(deep=True)
    ndw_stock.loc[:,"2015":"2022"] = 0
    

    for y in range(2023,2061):
        dw_stock[str(y)] = (dw_stock[str(y-1)] +
                            pd.concat([-dwel_d_rate, ((dw_stock[str(y-1)]
                                                       +ndw_stock[str(y-1)])
                                                      *(hh_proj_pch_ls[str(y)]))],
                                      axis=1).min(axis=1))
                         
        ndw_stock[str(y)] = ((dw_stock[str(y-1)]+ndw_stock[str(y-1)])
                         *(1+hh_proj_pch_ls[str(y)])-dw_stock[str(y)])

    
    # save results
    dw_stock.to_csv(snakemake.output.path_dw_stock_ex)
    ndw_stock.to_csv(snakemake.output.path_dw_stock_nb)
    
    # round projections
    dw_stock = dw_stock.apply(utils.round_col, axis=0)
    ndw_stock = ndw_stock.apply(utils.round_col, axis=0)

    
    # save results
    dw_stock.to_csv(snakemake.output.path_dw_stock_ex)
    ndw_stock.to_csv(snakemake.output.path_dw_stock_nb)
    
    logger.info("Finished processing and saved projections")

