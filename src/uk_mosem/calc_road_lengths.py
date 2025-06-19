"""
Script calculating road lengths per LSOA


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

import logging
import multiprocessing as mp

import pandas as pd
import numpy as np
import geopandas as gpd

import utils


logger = logging.getLogger(__name__)
utils.adjust_logger()

if __name__ == "__main__":

    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_road_lengths",
                                         **df_params.iloc[0].to_dict())
    
    # load list of areas/area lookup (LAD21CD)
    lsoas = utils.get_entity_lookup(["LSOA","LAD"])

    # load LSOA boundaries and road shapefiles, both 
    # OSGB36 / British National Grid (EPSG:27700) (meter as unit)
    roads = gpd.read_file(snakemake.input.path_sh_roads)
    lsoas_bounds = gpd.read_file(snakemake.input.path_sh_lsoas)
   
    
    # def overlay(lines, bounds=lsoas_bounds,how="union"):
    #     return gpd.overlay(lines,bounds,how=how)

    # calculate intersection using #threads cores
    roads_split = np.array_split(roads,256)

    roads_split = [[roads_split[i],lsoas_bounds,"union"]
                    for i in range(len(roads_split))]

    pool = mp.Pool(min(len(roads_split),snakemake.threads),maxtasksperchild=2)
    results = pool.starmap(gpd.overlay,roads_split)
    pool.terminate()

    lsoa_road_length = pd.concat(results)
    

    # process road parts that were not allocated    
    lsoa_road_length = lsoa_road_length.reset_index(drop=True).explode().reset_index(drop=True)
    last_coord_x = lsoa_road_length["geometry"].apply(lambda g: g.coords[-1][0])
    last_coord_y = lsoa_road_length["geometry"].apply(lambda g: g.coords[-1][1])
    lsoa_road_length["endpoint"] = gpd.points_from_xy(last_coord_x,last_coord_y)
    first_coord_x = lsoa_road_length["geometry"].apply(lambda g: g.coords[0][0])
    first_coord_y = lsoa_road_length["geometry"].apply(lambda g: g.coords[0][1])
    lsoa_road_length["startpoint"] = gpd.points_from_xy(first_coord_x,first_coord_y)
    
    # iteratively allocate road parts if start-/endpoint is also part of an
    # already allocated road part
    for i in range(10):
        if len(lsoa_road_length[lsoa_road_length["geo_code"].isna()])==0:
            break

        lrl=lsoa_road_length[lsoa_road_length["geo_code"].isna()].reset_index().merge(right=lsoa_road_length.loc[~lsoa_road_length["geo_code"].isna(),["startpoint","geo_code"]],
                                                                    left_on="endpoint",right_on="startpoint",how="left",suffixes=("","_s"))
        lrl = lrl.drop_duplicates("index")

        lrl=lrl.merge(right=lsoa_road_length.loc[~lsoa_road_length["geo_code"].isna(),["endpoint","geo_code"]],
                  left_on="startpoint",right_on="endpoint",how="left",suffixes=("","_e"))
        lrl = lrl.drop_duplicates("index")

        lrl.loc[~lrl["geo_code_s"].isna(), "geo_code"] = lrl["geo_code_s"]

        lrl.loc[~lrl["geo_code_e"].isna(), "geo_code"] = lrl["geo_code_e"]
        lsoa_road_length[lsoa_road_length["geo_code"].isna()] = lrl.set_index("index")[lsoa_road_length.columns]

    else:
        logger.warning(str(len(lsoa_road_length[lsoa_road_length["geo_code"].isna()]))
              +" road parts could not be allocated.")
    
    
    # merge with LSOA DataFrame
    
    lsoa_road_length['length'] = lsoa_road_length['geometry'].length
    lsoa_road_length = lsoa_road_length[["geo_code",
                                  "length"]].groupby(["geo_code"]).sum().reset_index()
    lsoa_road_length = lsoa_road_length.rename(columns = {"geo_code":"LSOA11CD"})
                          
    lsoa_road_length = lsoas.merge(right=lsoa_road_length,on="LSOA11CD",
                                    how="left")
    # save to file
    # unit: m
    lsoa_road_length.to_csv(snakemake.output.path_road_lengths)
