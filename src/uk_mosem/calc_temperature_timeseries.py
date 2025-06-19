"""
Script deriving climate/temperature timeseries from MetOffice data


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

import datetime

import pandas as pd
import numpy as np
import rasterio as rio
import rasterio.mask as riom
import geopandas as gpd

import utils

import logging

logger = logging.getLogger(__name__)
utils.adjust_logger()

if __name__ == "__main__":
    
    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_temperature_timeseries",
                                         **df_params.iloc[0].to_dict())
        
    files = sorted(snakemake.input.path_climate_dir)
    files = [file for file in files if file.endswith(".nc")]
    
    tas = rio.open(files[0])
    
    lads = gpd.read_file(snakemake.input.path_lad_bounds) 
    
    # calulate mask/raster-LAD overlaps for all LADs
    def get_mask(polygon, data = tas):        
        mask,affine,none = riom.raster_geometry_mask(dataset=tas,
                                                         shapes=[polygon,],
                                                         crop=False,
                                                         all_touched=True,
                                                         invert=True)
        return mask

    logger.info("Creating masks for grid cells.")
    lads["mask"] = lads.geometry.apply(get_mask)
    logger.info("Created masks for grid cells.")
    # create dataframe with respective time index, i.e., hourly with 30 days
    # per month as for the underlying climate data
    ts = pd.date_range(datetime.datetime(1999, 1, 1, 0),
                            periods=30*24, freq='H')
    ts = [pd.date_range(datetime.datetime(1999, m, 1, 0),
                                  periods=30*24, freq='H') 
                    if m != 2 else
                    pd.date_range(datetime.datetime(1999, m, 1, 0),
                                  periods=28*24, freq='H')
                    for m in range(1,13) ]
    
    # read in temperature data, apply mask, i.e., only include grid squares
    # that are relevant to LAD (binary approach - no fractions used) and 
    # create timeseries in geopandas DataFrame
    logger.info("Applying masks for grid cells.")
    for n, f in enumerate(files):
        logger.info(f"Working on temperature data file {n}.")
        tas = rio.open(f)
        for nn,ind in enumerate(tas.indexes):
            if n==1 and nn >= 672:
                break
            t = tas.read(ind)
            lads[ts[n][nn]] = [np.mean(t[mask]) for mask in lads["mask"]]

    logger.info("Applied masks for grid cells.")
    # create new DataFrame with timeseries
    temp = pd.DataFrame(lads.loc[:,datetime.datetime(1999,1,1,0):])
    temp.index = lads["LAD23CD"]
    

    # load list of areas/area lookup (LAD23CD)
    lads = utils.get_entity_lookup(levels=["LAD"])
    
    
    temp = lads.merge(right=temp,on="LAD23CD",how="left")
    temp = temp.set_index("LAD23CD")
    
    # save to file
    logger.info("Saving file.")
    temp.to_csv(snakemake.output.path_temp)
    logger.info("Saved file.")
