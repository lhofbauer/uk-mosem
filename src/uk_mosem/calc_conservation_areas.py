"""
Script loading and processing conservation areas


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

import pandas as pd
import geopandas as gpd

import utils



if __name__ == "__main__":
    
    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_conservation_areas",
                                         **df_params.iloc[0].to_dict())

    # load LSOAs
    lsoas = utils.get_entity_lookup(["LSOA"])
    
    # load conservation areas
    carea_E = gpd.read_file(snakemake.input.path_conserv_areas_e)
    carea_W = gpd.read_file(snakemake.input.path_conserv_areas_w)
    carea_S = gpd.read_file(snakemake.input.path_conserv_areas_s)

    # load listed buildings areas (not currently used)
    # listedB_E = gpd.read_file(PATHR+PATH_LISTED_BUILD_E)
    # listedB_W = gpd.read_file(PATHR+PATH_LISTED_BUILD_W)
    # listedB_S = gpd.read_file(PATHR+PATH_LISTED_BUILD_S)
    
    # Comment here and the two lines below are to check fraction of listed
    # buildings in conservation areas
    # listedB = pd.concat([listedB_E,listedB_W,listedB_S])
    # listedB["lb"] = True

    buildings = gpd.read_file(snakemake.input.path_build_map)
    buildings["number"] = 1
    lsoas_bounds = gpd.read_file(snakemake.input.path_sh_lsoas)

    
    # remove 'missing district data'
    carea_E = carea_E[carea_E.UID!=1]
    carea = pd.concat([carea_E,carea_W,carea_S])
    carea["CA"] = True
    
    # calculating fraction of listed buildings that is in conservation areas
    # lp = gpd.sjoin(carea,listedB,how = 'left',
    #                            op ='intersects')
    # lbinca = len(lp[lp["lb"]==True])/len(listedB)

    # process conservation areas
    carea.geometry = carea.buffer(0)
    carea_lsoa = gpd.overlay(carea,lsoas_bounds,how="union", keep_geom_type=True)

    # drop geometries that lie outside LSOAS (at least partly due to inaccurate
    # boundaries for conservation areas reaching, e.g., into sea)
    carea_lsoa = carea_lsoa[~carea_lsoa["geo_code"].isna()]
    carea_lsoa = carea_lsoa[["geometry","CA","geo_code"]]
    carea_lsoa["ID"] = range(0,len(carea_lsoa))

    # match building with LSOAs
    carea_lsoa_bui = gpd.sjoin(carea_lsoa,buildings,how = 'left',
                               op ='intersects')
    
    # calculate fraction of building in careas in each LSOA
    carea_lsoa_bui = carea_lsoa_bui[["ID","number"]].groupby("ID").sum()
    carea_lsoa_bui = carea_lsoa[["ID",
                                 "CA",
                                 "geo_code"]].merge(right=carea_lsoa_bui.reset_index(),
                                                    on="ID",
                                                    how="left")

    carea_lsoa_bui = carea_lsoa_bui.drop("ID",axis=1)                                                
    carea_lsoa_bui = carea_lsoa_bui.groupby(["CA",
                                             "geo_code"]).sum().reset_index()
    carea_lsoa_bui = carea_lsoa_bui.merge(right=carea_lsoa_bui.groupby("geo_code").sum()["number"],
                                          on="geo_code",
                                          how="left",
                                          suffixes=("","_total"))
    carea_lsoa_bui = carea_lsoa_bui.loc[carea_lsoa_bui["CA"]==True]

    carea_lsoa_bui["cafrac"] = (carea_lsoa_bui["number"]
                                /carea_lsoa_bui["number_total"])
    carea_lsoa_bui["cafrac"] = carea_lsoa_bui["cafrac"].fillna(0)
    
    
    carea_lsoa_bui = carea_lsoa_bui.rename(columns={"geo_code":"LSOA11CD"})
    carea_lsoas = lsoas.merge(right=carea_lsoa_bui[["LSOA11CD","cafrac"]],
                              on="LSOA11CD", how="left").fillna(0)    
    carea_lsoas = carea_lsoas.set_index("LSOA11CD")


    # save to file
    # unit: -
    carea_lsoas.to_csv(snakemake.output.path_conserv_areas)
  
