"""
Script deriving a classification for sublocal areas (used for aggregation)


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

import pandas as pd

import utils

if __name__ == "__main__":

    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_sublocal_agg",
                                         **df_params.iloc[0].to_dict())    

    # load DH classification scheme
    raw_data = pd.read_csv(snakemake.input.path_set_dh_cl)
    DHclasses = list(raw_data["VALUE"][raw_data["VARIABLE"]
                                       =="SublocalDHClasses"])
    # add class for everything below
    DHclasses.append(-10)
    DHclasses.sort()
    # load demand density
    demdens = pd.read_csv(snakemake.input.path_dem_density,
                            index_col=["LSOA11CD"])

    # derive sublocal aggregation    
    # FIXME: Use this function to implement a more complex aggegration which 
    # also takes gas grid, etc. into account, if required
    def aggregation(x):
        x = "DH_"+str(DHclasses.index(max([l for l in DHclasses
                                           if x["DDensity"]>l])))
        return x
    
    slagg = demdens.apply(aggregation,axis=1)    
    slagg.name = "LSOA11CD_AGG"

    # save to file
    # unit: -
    slagg.to_csv(snakemake.output.path_sublocal_agg)
    

