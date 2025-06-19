"""
Script deriving residual capacity fractions of heating technologies


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

import pandas as pd

import utils

import logging

logger = logging.getLogger(__name__)
utils.adjust_logger()

if __name__ == "__main__":
    
    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_residual_frac_heat_tech",
                                         **df_params.iloc[43].to_dict())
    
    # load heating system data per LSOA from EPC data
    hsys_dw = pd.read_csv(snakemake.input.path_dw_stock_hde, index_col=["LSOA11CD",
                                                              "PROPERTY_TYPE",
                                                              "HEATING_SYSTEM"])
    hsys_nd = pd.read_csv(snakemake.input.path_nd_stock_hde, index_col=["LSOA11CD",
                                                              "PROPERTY_TYPE",
                                                              "YEAR",
                                                              "HEATING_SYSTEM"])
    hsys_nd = hsys_nd.xs(2022,level=2)
    hsys_dw.columns = ["NUMBER_AREA"]
    hsys_nd.columns = ["NUMBER_AREA"] 
    
    hsys_nd = hsys_nd.reset_index()
    hsys_nd["PROPERTY_TYPE"] = "Non-domestic"
    hsys_nd = hsys_nd.groupby(["LSOA11CD",
                               "PROPERTY_TYPE",
                               "HEATING_SYSTEM"]).sum()
    hsys = pd.concat([hsys_dw.reset_index(),hsys_nd.reset_index()])
    hsys["HEATING_SYSTEM"] = hsys["HEATING_SYSTEM"].fillna("NA")
    
    hsys = hsys.groupby(["LSOA11CD",
                         "PROPERTY_TYPE",
                         "HEATING_SYSTEM"],dropna=False).sum()
    
    # process heating and heat distribution separately
    hsyss = dict()
    for t in ["MH","HD"]:
        hsyss[t] = hsys.loc[hsys.index.get_level_values("HEATING_SYSTEM").str.startswith(t)]
        
        # take into account NAs
        # if non-NA exists for LSOA and property type, ignore NA and calc fraction
        # if only NA exists calc fraction based on average of ... (happening below)
        hsyss[t] = hsyss[t].unstack(fill_value=0)
        hsyss[t].loc[hsyss[t].sum(axis=1)>hsyss[t].loc[:,("NUMBER_AREA",t+"_NA")],
                 ("NUMBER_AREA",t+"_NA")] = 0  
        hsyss[t] = hsyss[t].stack()
        
        # calculate fractions
        hsyss[t] = hsyss[t].reset_index().merge(right=hsyss[t].groupby(["LSOA11CD",
                                "PROPERTY_TYPE"]).sum()["NUMBER_AREA"],
                                  on=["LSOA11CD","PROPERTY_TYPE"], how="left",
                                  suffixes=["","_total"])
        hsyss[t]["NUMBER_AREA"] = (hsyss[t]["NUMBER_AREA"]/
                               hsyss[t]["NUMBER_AREA_total"]).fillna(0)
        
        hsyss[t] = hsyss[t].set_index(["LSOA11CD","PROPERTY_TYPE",
                        "HEATING_SYSTEM"])["NUMBER_AREA"]
       
        if hsyss[t].xs(t+"_NA",level=2).sum()>0:
            logger.warning("One or more property types in one or more LSOAs"
                           " do not have any data to calculate residual heating"
                           " system fraction. Using LAD average for property type.")
            
            logger.warning(hsyss[t].xs(t+"_NA",level=2)[hsyss[t].xs(t+"_NA",level=2)>0])
            # calculate mean for LADs
            LAD_mean = hsyss[t].unstack(fill_value=0)
            LAD_mean = LAD_mean.loc[(LAD_mean[t+"_NA"]==0)
                                    &(LAD_mean.sum(axis=1)>0)].drop(t+"_NA",axis=1)
            LAD_mean = utils.groupby_LAD(LAD_mean).mean()
            LAD_mean = LAD_mean.stack()
            LAD_mean.index = LAD_mean.index.rename("HEATING_SYSTEM",level=2)
            LAD_mean.name = "NUMBER_AREA"
            
            # merge LAD means on hsys data
            lsoas = utils.get_entity_lookup(["LSOA","LAD"]).drop("LAD23NM",axis=1)
            hsyss[t] = hsyss[t].reset_index().merge(right=lsoas,on="LSOA11CD",how="left")
            #hsyss[t] = hsyss[t].loc[hsyss[t]["HEATING_SYSTEM"]!="NA"]
            hsyss[t] = hsyss[t].merge(right=LAD_mean,
                              on=["LAD23CD","PROPERTY_TYPE","HEATING_SYSTEM"],
                              how="left",suffixes = ["","_LAD"])
            hsyss[t] = hsyss[t].drop("LAD23CD",axis=1)
            hsyss[t] = hsyss[t].set_index(["LSOA11CD","PROPERTY_TYPE","HEATING_SYSTEM"])
            hsyss[t] = hsyss[t].unstack()
            
            # replace LSOA values with LAD averages if only NA values exist
            for c in hsyss[t].columns.get_level_values(1).unique():
                if c == t+"_NA":
                    continue
                hsyss[t].loc[hsyss[t][("NUMBER_AREA",t+"_NA")]>0,
                         ("NUMBER_AREA",c)] = hsyss[t].loc[hsyss[t][("NUMBER_AREA",t+"_NA")]>0,
                                                                 ("NUMBER_AREA_LAD",c)]
            hsyss[t] = hsyss[t].drop("NUMBER_AREA_LAD",axis=1)
            hsyss[t] = hsyss[t].drop(t+"_NA",axis=1,level=1)
            hsyss[t] = hsyss[t].stack(dropna=False)
            
            if hsyss[t]["NUMBER_AREA"].isna().any():
                logger.warning("Was not able to use LAD average in all instances.")
    

    hsys = pd.concat([df for df in hsyss.values()])
    
    if not isinstance(hsys,pd.Series):
        hsys = hsys["NUMBER_AREA"]
    
    # remove no heating system columns and sum up low and normal temperature
    # wet systems
    hsys = hsys.loc[~(hsys.index.get_level_values(2)=="HD_NONE"),:]
    hsys.loc[slice(None),
                slice(None),
                "HD_WDIS"]=(hsys.loc[slice(None),
                                  slice(None),
                                  "HD_WDIS"].add(
                        +hsys.loc[slice(None),
                                  slice(None),
                                  "HD_RAUP"],fill_value=0))
                           
    # rearrange data
    hsys = hsys.reset_index()
    hsys["HEATING_SYSTEM"] = hsys["HEATING_SYSTEM"].str[3:]
    hsys = hsys.rename(columns={"HEATING_SYSTEM":"TECHNOLOGY"})
    hsys = hsys.set_index(["LSOA11CD","PROPERTY_TYPE",
                    "TECHNOLOGY"])
    hsys.columns=["fraction"]
    # save to file
    # unit: -
    hsys.to_csv(snakemake.output.path_res_caps_ht_f)
  

