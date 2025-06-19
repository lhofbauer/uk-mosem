"""
Script deriving residual capacities of heating technologies


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

import pandas as pd

import utils

if __name__ == "__main__":
    
    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_residual_capacities_heat_tech",
                                         **df_params.iloc[0].to_dict())

    # load peak heat demands per LSOA and property type
    hpeaks = pd.read_csv(snakemake.input.path_peakc_lsoa,
                            index_col=["LSOA11CD",
                                       "PROPERTY_TYPE",
                                       "TECHNOLOGY"])   
    
    # aggregate to LADs
    #hpeaks = utils.groupby_LAD(hpeaks).sum()
    
    # load heating system capacity fractions 
    hsys = pd.read_csv(snakemake.input.path_res_caps_ht_f, index_col=["LSOA11CD",
                                                           "PROPERTY_TYPE",
                                                           "TECHNOLOGY"])
    

    # set a linear trajectory for base years
    # for i in range(2016,2023):
    #     hpeaks.loc[:,str(i)] = (hpeaks.loc[:,str(i-1)]
    #                               +(hpeaks.loc[:,"2023"]-hpeaks.loc[:,"2015"])
    #                               /(2023-2015)
    #                               )
    # set 2022 as 2023 – this is to make sure the residual capacity in 2023 is
    # reasonable and not too low (as 2022 has low demand = low capacity)
    # hpeaks.loc[:,"2022"] = hpeaks.loc[:,"2023"]
    # adjust 2022 capacity (see reasoning above/in the documentation)
    hpeaks.loc[:,"2022"] = (hpeaks.loc[:,"2023"]
                              -(hpeaks.loc[:,"2023"]-hpeaks.loc[:,"2015"])
                              /(2023-2015)
                              )

       
    residual_caps = hpeaks.loc[:,"2015":"2022"].merge(right=hsys,
                                                      left_index=True,
                                                      right_index=True,
                                                      how="right")
    
                                                         

    residual_caps = residual_caps.loc[:,"2015":"2022"].multiply(residual_caps["fraction"],axis="index")
    
    
    # translate to actual technology names by flattening index and then
    # renaming
    
    # nondom = ["Community, arts & leisure", "Hospitality", "Industrial",
    #           "Offices", "Others", "Retail", "Education", "Emergency services",
    #           "Health", "Storage"]
    if snakemake.params.dic["scen_hh_disagg"].startswith("T"):
            residual_caps.index = pd.MultiIndex.from_tuples([ (x[0],x[1],x[2]+"DDDE"+x[1].split("|")[1])
                                                  if x[1].startswith("Detached") else
                                                  (x[0],x[1],x[2]+"DDSD"+x[1].split("|")[1])
                                                  if x[1].startswith("Semi-detached") else
                                                  (x[0],x[1],x[2]+"DDTE"+x[1].split("|")[1])
                                                  if x[1].startswith("Terraced") else
                                                  (x[0],x[1],x[2]+"DDFL"+x[1].split("|")[1])
                                                  if x[1].startswith("Flats") else
                                                  (x[0],x[1],x[2]+"DNDO00")
                                                  if x[1]=="Non-domestic" else 0 for x
                                                  in residual_caps.index])
    else:              
        residual_caps.index = pd.MultiIndex.from_tuples([ (x[0],x[1],x[2]+"DDDE00")
                                                         if x[1]=="Detached" else
                                                         (x[0],x[1],x[2]+"DDSD00")
                                                         if x[1]=="Semi-detached" else
                                                         (x[0],x[1],x[2]+"DDTE00")
                                                         if x[1]=="Terraced" else
                                                         (x[0],x[1],x[2]+"DDFL00")
                                                         if x[1]=="Flats" else
                                                         (x[0],x[1],x[2]+"DNDO00")
                                                         if x[1]=="Non-domestic" else 0 for x
                                                         in residual_caps.index])
    residual_caps.index.names = ["LSOA11CD","PROPERTY_TYPE","TECHNOLOGY"]
    
    # save to file
    # unit: GW
    residual_caps.to_csv(snakemake.output.path_res_caps_ht)

