"""
Script deriving local renewable potentials and related cost for biomass


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

import pandas as pd
import numpy as np

import geopandas as gpd

import utils

if __name__ == "__main__":

    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_potentials",
                                         **df_params.iloc[0].to_dict())    
 
    # load LADs and land cover shape files
    lads = gpd.read_file(snakemake.input.path_lad_bounds) 
    lc = gpd.read_file(snakemake.input.path_land_cover)
    
    # intersect each LAD with relevant land cover areas
    data = list()
    for _, lad in lads.iterrows():
        for __, lca in lc.iterrows():
            if ((lca['geometry'] is not None)
                and (211 <= int(lca['CODE_18']) <= 313)
                and lca['geometry'].intersects(lad['geometry'])):
                data.append({'LAD23CD':lad['LAD23CD'], 'lccode': lca['CODE_18'],
                             'area':lad['geometry'].intersection(lca['geometry']).area})

    ladslc = pd.DataFrame(data,columns=['LAD23CD', 'lccode',
                                            'area'])
    
    ladslc["lccode"] = ["forest" if 244 <= int(i) <= 313 else
                    "agri" if 211 <= int(i) <= 243 else
                    None for i in ladslc["lccode"]]
    ladslc = ladslc.dropna()
    ladslc = ladslc.groupby(["LAD23CD","lccode"]).sum()

    # load biomass potential
    
    bmp = pd.read_excel(snakemake.input.path_bm_potential,
                          sheet_name="ENER - NUTS2 BioCom E",
                          header=0,index_col=(0,1,2,3,4,5,7))
    
    bmc = pd.read_excel(snakemake.input.path_bm_potential,
                          sheet_name="COST - NUTS2 BioCom",
                          header=0,index_col=(0,1,2,3,4,5,7))
    bmc = bmc.rename_axis(index={"Energy commoditty ":"E-Comm",
                             "Bio Commodity":"B-Com",
                             "NUTS2":"NUST2"})
    bmc = bmc.droplevel("Units")
    
    bmp = bmp.merge(right=bmc,on=["Year","Scenario","NUTS0","NUST2",
                                  "E-Comm","B-Com"],how="left")
    bmp = bmp.drop(columns=[c for c in bmp.columns if "Unnamed" in c])
    
    # considering only forest residues/products and agricultural residues
    # as biomass to be used in the model
    # (other not-energy-crop categories would be:  MINBIOGAS1)                                                                
    bmp = bmp.rename(index={"MINBIOFRSR1":"forest",
                            "MINBIOWOO":"forest",
                            "MINBIOWOOa":"forest",
                            "MINBIOWOOW1":"forest",
                            "MINBIOWOOW1a":"forest",
                            "MINBIOFRSR1a":"agri",
                            "MINBIOAGRW1":"agri"})
    
    bmp = bmp.loc[(slice(None), "ENS_Med", "UK", slice(None),
                   ("forest","agri"), slice(None), slice(None)),
                  :].reset_index(level=["Scenario","NUTS0"], drop=True) 
                                                                    
    # fill agricultural waste with no cost given with average of other waste
    bmp.loc[bmp.index.get_level_values("E-Comm")
            =="agri",
            "NUTS2 Bio Commodity Cost "] = bmp.loc[bmp.index.get_level_values("E-Comm")
                                =="agri","NUTS2 Bio Commodity Cost "].fillna(bmp.loc[bmp.index.get_level_values("E-Comm")
                                                    =="agri","NUTS2 Bio Commodity Cost "].mean())
    # fill one remaining wood residue missing cost value with average across LADs                           
    bmp = bmp.fillna(bmp.xs("Fuelwood Res",level="B-Com").mean())

    
    bmp.columns = ["potential","cost"]
    
    wm = lambda x: np.average(x, weights=bmp.loc[x.index, "potential"]) if bmp.loc[x.index, "potential"].sum()!=0 else 0

    bmp = bmp.groupby([i for i in bmp.index.names
                        if "B-Com" not in i]).agg(potential=("potential", "sum"),  
                                                  cost=("cost", wm))     

    
    bmp = bmp.reset_index()
    bmp["NUST2"] = bmp["NUST2"].str.replace("UK","TL")
    bmp = bmp.rename(columns={"NUST2":"ITL221CD"})
    
    la = utils.get_entity_lookup(["LAD","ITL2","ITL3"])
    
    # rearrange lookup to fit old version of code used in ENSPRESO data (see
    # data file 'NUTS2 conversion' sheet for details)
    la["ITL221CD"] = ["TLM2" if ("M7" in c2) and ("M24" not in c3) else
                      "TLM3" if ("M8" in c2) or ("M9" in c2) or ("M24" in c3) else
                      c2 for c2,c3 in zip(la["ITL221CD"],la["ITL321CD"])]
    la = la.drop("ITL321CD",axis=1)
    la = la[~la.duplicated(keep="first")]
    
    # merge potential and cost on LADs and process,i.e., set values where no
    # potential
    bmp = la.merge(right=bmp, on="ITL221CD", how="left")
    bmp.loc[bmp["Year"].isna(),"potential"] = 0
    bmp.loc[bmp["Year"].isna(),"cost"] = 0
    bmp.loc[bmp["Year"].isna(),"E-Comm"] = "forest"
    bmp.loc[bmp["Year"].isna(),"Year"] = 2010
    
    bmp = bmp.rename(columns={"E-Comm":"lccode"})
    
    bmp = bmp.merge(right=ladslc,how="left",on=["LAD23CD","lccode"])
    
    # where no relevant land cover is present, i.e., value is nan, set to 0
    bmp["area"] = bmp["area"].fillna(0)
    
    # merge with sum of area in each ITL
    bmp = bmp.merge(right=bmp[["ITL221CD",
                               "lccode",
                               "Year",
                               "area"]].groupby(["ITL221CD",
                                                 "lccode",
                                                 "Year"]).sum(),
                    how="left",on=["ITL221CD","lccode","Year"])
                                                 
    # calculate potential in each LAD based on fraction of relevant area in LAD
    # with respect to the ITL                                             
    bmp["lcfrac"] = bmp["area_x"]/bmp["area_y"]
    bmp["lcfrac"] = bmp["lcfrac"].fillna(0)
    bmp = bmp.drop(["area_x","area_y"],axis=1)
    bmp["potential"] = bmp["potential"]*bmp["lcfrac"]
    
    # groupby currently not really necessary, but if different types of biomass
    # are used in the calculations (then this would need to be adapted to pick/
    # groupby the right types for different potentials)
    bmp = bmp[["LAD23CD","Year","potential",
               "cost"]].groupby(["LAD23CD",
                                 "Year"]).agg(potential=("potential", "sum"),  
                                              cost=("cost", wm))     

    
    # convert from PJ to TJ
    bmp["potential"] = bmp["potential"] * 1000
    # convert from €/GJ to €/MJ  (= Mio. €/TJ)
    bmp["cost"] = bmp["cost"] / 1000
    
    # calc costs of pellet production (currently all biomass assumed to be used
    # as pellets); added below
    
    raw_data = pd.read_csv(snakemake.input.path_set_supply)
    
    pcost = raw_data.loc[raw_data["VARIABLE"]=="PelletProductionCost",["UNIT","VALUE"]] 
    pcost = utils.adjust_monetary_values(pcost, 2015)["VALUE"]
    
    # arrange dataframe
    bmpc = bmp["cost"].reset_index()
    bmpc["UNIT"] = "2010€/MJ"
    bmpc["VARIABLE"] = "VariableCost"
    bmpc["TECHNOLOGY"] = "BMEXSEXT00"
    bmpc["MODE_OF_OPERATION"] = "1,2"
    
    bmpc = bmpc.assign(MODE_OF_OPERATION=bmpc["MODE_OF_OPERATION"]
                       .str.split(',')).explode("MODE_OF_OPERATION")
    bmpc["MODE_OF_OPERATION"] = bmpc["MODE_OF_OPERATION"].astype(int)
    
    bmpc = bmpc.rename(columns={"LAD23CD":"REGION",
                              "cost":"VALUE",
                              "Year":"YEAR"})
    
    # add pellet production cost
    bmpc = utils.adjust_monetary_values(bmpc, 2015)
    bmpc["VALUE"] = bmpc["VALUE"] + pcost.squeeze()
    
    
    bmpc = bmpc.set_index([col for col in bmpc.columns if col!="VALUE"])    
    
    # save to file
    # unit: 2015£/TJ
    bmpc.to_csv(snakemake.output.path_lsupply_cost)    
    
    # arrange dataframe
    bmp = bmp["potential"].reset_index()
    bmp["UNIT"] = "TJ"
    bmp["VARIABLE"] = "TotalTechnologyAnnualActivityUpperLimit"
    bmp["TECHNOLOGY"] = "BMEXSEXT00"
    bmp = bmp.rename(columns={"LAD23CD":"REGION",
                              "potential":"VALUE",
                              "Year":"YEAR"})
    bmp = bmp.set_index([col for col in bmp.columns if col!="VALUE"])    
    
    # save to file
    # unit: TJ
    bmp.to_csv(snakemake.output.path_potential_con)
