"""
Script deriving household projections as basis for building stock projections


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

import pandas as pd

import utils


if __name__ == "__main__":

    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_hh_projections",
                                         **df_params.iloc[0].to_dict())    
        
    # load area lookup (LAD23CD)
    lsoas = utils.get_entity_lookup(["LSOA","LAD"])


    # load household projections (E:LAD17CD, S&W: LAD23CD - no recent changes)
    hh_proj_s = pd.read_csv(snakemake.input.path_dw_hh_p_s,skiprows=3,nrows=34,
                            thousands=",",encoding = "ISO-8859-1")
    hh_proj_e = pd.read_excel(snakemake.input.path_dw_hh_p_e,
                             sheet_name="406",skiprows=4) 
    hh_proj_w = pd.read_csv(snakemake.input.path_dw_hh_p_w,skiprows=7,nrows=23) 
    
    # rearrange data
    hh_proj_w = hh_proj_w.drop(columns=["Unnamed: 0"],index=[0])
    hh_proj_w = hh_proj_w.rename(columns={"Unnamed: 1":"LAD23NM"})
    hh_proj_w = hh_proj_w.set_index("LAD23NM")
    hh_proj_w.columns = hh_proj_w.columns.str.strip()
    
    hh_proj_s = hh_proj_s.drop(index=[0,1])
    hh_proj_s = hh_proj_s.rename(columns={"Unnamed: 0":"LAD23NM"})  
    hh_proj_s = hh_proj_s.set_index("LAD23NM")
    hh_proj_s = hh_proj_s.loc[:, :'2043']
    
    hh_proj_e = hh_proj_e.drop(index=[0],columns=["Area name"])
    hh_proj_e = hh_proj_e.rename(columns={"Area code":"LAD17CD"}) 
    hh_proj_e = hh_proj_e.set_index("LAD17CD")


    hh_proj_e = utils.update_LADCD(hh_proj_e,from_CD="LAD17CD",how="sum")
    

    hh_proj_sw = pd.concat([hh_proj_s,hh_proj_w], sort=True)

    # set LAD codes instead of names as index
    hh_proj_sw.index = hh_proj_sw.index.str.strip().map(lsoas[["LAD23NM",
                "LAD23CD"]].set_index("LAD23NM").drop_duplicates().iloc[:,0])
    hh_proj_sw.index.name = "LAD23CD"
    hh_proj_sw.columns = hh_proj_sw.columns.astype(int)
    
    # concat all dataframes
    hh_proj_esw = pd.concat([hh_proj_e,hh_proj_sw], sort=True)
    
    # rearrange and delete aggregate county, etc. data
    hh_proj = lsoas["LAD23CD"].drop_duplicates()
    hh_proj = hh_proj.to_frame().merge(right=hh_proj_esw, on="LAD23CD",
                                       how="left")
    hh_proj = hh_proj.set_index("LAD23CD")
    
    # save household numbers (only for result plotting etc., not used in
    # in further model input calculations)
    hh_proj.to_csv(snakemake.output.path_hh_num)
    
    # calculate household percentage change per LAD/LSOA
    hh_proj_pch = hh_proj.pct_change(axis="columns",periods=1)
    
    hh_proj_pch_ls = lsoas.merge(right=hh_proj_pch,
                                 on="LAD23CD",how="left").set_index("LSOA11CD")
    
    # continue with development until last model year
    hh_proj_pch_ls = hh_proj_pch_ls.assign(**{str(i):hh_proj_pch_ls[list(
                                            range(2039,2044))].mean(axis=1)
                                              for i in range(2044,2061)})
    hh_proj_pch_ls.columns = hh_proj_pch_ls.columns.astype(str)
    hh_proj_pch_ls = hh_proj_pch_ls.drop(["LAD23CD","LAD23NM"],axis=1)
    
    # save household percentage change from year to year for each LSOA
    # unit: -
    hh_proj_pch_ls.to_csv(snakemake.output.path_hh_proj)
    