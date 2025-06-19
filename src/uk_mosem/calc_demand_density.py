"""
Script deriving linear heat demand densities for each local area


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

import pandas as pd

import utils

if __name__ == "__main__":

    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_demand_density",
                                         **df_params.iloc[0].to_dict())
        
    # load annual heat demand per LSOA 
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
    anndem = pd.concat([dw_anndem,nd_anndem])
      
    # load road lengths
    rlengths = pd.read_csv(snakemake.input.path_road_lengths,
                            index_col=["LSOA11CD"])

    # calculate linear heat demand density
    anndem = anndem.groupby(["LSOA11CD"]).sum()
    #anndem["mean"] = anndem.mean(axis=1)                           
    
    demdens = anndem["2020"].divide(rlengths['length'],axis=0)
    demdens.name = 'DDensity'
    
    # save to file
    # unit: GJ/m
    demdens.to_csv(snakemake.output.path_dem_density)


 

