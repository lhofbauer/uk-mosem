"""
Script deriving potentials for DH heat sources


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

import pandas as pd
import numpy as np
import geopandas as gpd
import rasterio as rio
import utils
from shapely.geometry import LineString

import xarray as xr



if __name__ == "__main__":

    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_heat_sources",
                                         **df_params.iloc[0].to_dict())

    
#%%
# industrial excess heat potential

    # load area lookup (LAD23CD)
    lsoas = utils.get_entity_lookup(["LSOA","LAD"])
    
    # load excess heat data and arrange (in GWh/a)
    indeh = pd.read_csv(snakemake.input.path_ind_eh, usecols=["Country","geom","CityCode",
                                                    "Excess_Heat_Total"],
                        sep=";")
    indeh = indeh[indeh["Country"]=="United Kingdom"].drop("Country",axis=1)
    indeh = indeh[~indeh["geom"].isna()]
    indeh["geometry"] = indeh["geom"].str.split(";",expand=True)[1]
    indeh['geometry'] = gpd.GeoSeries.from_wkt(indeh['geometry'])
    indeh = gpd.GeoDataFrame(indeh, geometry='geometry')
    indeh = indeh.set_crs(epsg=4326)
    indeh = indeh.to_crs(epsg=27700)

    # load LSOAs boundaries
    lsoas_bounds = gpd.read_file(snakemake.input.path_sh_lsoas)
    
    # allocate heat sources to LSOAs
    lsoas_ieh = gpd.sjoin(lsoas_bounds,indeh, how = 'left', op ='intersects')
    
    lsoas_ieh_df = lsoas_ieh[["geo_code",
                              "Excess_Heat_Total"]].groupby("geo_code").sum()
    lsoas_ieh_df.index.name = "LSOA11CD"
    
    # merge to LADs
    ieh = lsoas.merge(right=lsoas_ieh_df.reset_index(),
                      on="LSOA11CD",how="left")
    ieh = ieh[["LAD23CD","Excess_Heat_Total"]].groupby("LAD23CD").sum()
    
    # convert to TJ
    ieh = ieh*3600/1000
    
    # rearrange dataframe
    ieh = ieh.reset_index()
    ieh["TECHNOLOGY"] = "IEHESDIS00"
    ieh["YEAR"] = ":*"
    ieh["VARIABLE"] = "TotalTechnologyAnnualActivityUpperLimit"

    ieh = ieh.rename(columns={"LAD23CD":"REGION",
                              "Excess_Heat_Total":"VALUE"})
    ieh = ieh.reindex(sorted(ieh.columns), axis=1)
    ieh = ieh.set_index([col for col 
                             in ieh.columns
                             if col!="VALUE"])    

#%%
# waste water heat potential

    # load area lookup (LAD23CD)
    lsoas = utils.get_entity_lookup(["LSOA","LAD"])
    
    # load raw data set and parameters
    raw_data = pd.read_csv(snakemake.input.path_set_hsources)
    eff_data = pd.read_csv(snakemake.input.path_set_eff)
    
    
    wwtdiff = raw_data.loc[raw_data["VARIABLE"]=="TemperatureDiffWasteWater","VALUE"].squeeze() 
    wwpp = raw_data.loc[raw_data["VARIABLE"]=="WasteWaterperPerson","VALUE"].squeeze() 
    swhp_eff = eff_data.loc[(eff_data["TECHNOLOGY"]=="SWHPSDIS00"), "VALUE"].values[0]
    
    # load household population projections (E:LAD17CD,
    # WS: LAD23CD - no recent changes)
    hh_proj_s = pd.read_csv(snakemake.input.path_dw_hhp_p_s,skiprows=2,nrows=1122,
                            thousands=",",encoding="ISO-8859-1")
    hh_proj_e = pd.read_excel(snakemake.input.path_dw_hh_p_e,
                             sheet_name="426",skiprows=4) 
    hh_proj_w = pd.read_csv(snakemake.input.path_dw_hhp_p_w,skiprows=7,nrows=23) 
    
    
    hh_proj_w = hh_proj_w.drop(columns=["Unnamed: 0"],index=[0])
    hh_proj_w = hh_proj_w.rename(columns={"Unnamed: 1":"LAD23NM"})
    hh_proj_w = hh_proj_w.set_index("LAD23NM")
    hh_proj_w.columns = hh_proj_w.columns.str.strip()
    
    hh_proj_s = hh_proj_s.drop(["Gender","Age Group"],axis=1)
    hh_proj_s = hh_proj_s.rename(columns={"Local Authority":"LAD23NM"})
    hh_proj_s = hh_proj_s.groupby("LAD23NM").sum()
    hh_proj_s = hh_proj_s.drop("Scotland")
    
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
    
    
    # extend to all model years
    hh_proj = lsoas["LAD23CD"].drop_duplicates()
    hh_proj = hh_proj.to_frame().merge(right=hh_proj_esw, on="LAD23CD",
                                       how="left")
    hh_proj = hh_proj.set_index("LAD23CD")
    # assume similar development to end of modelling period
    # FIXME: get end/start of modelling period from sets in raw data sheet, 
    # not set in code
    hh_proj = hh_proj.assign(**{str(i):(hh_proj[2043]-hh_proj[2038])/5*(i-2043)+hh_proj[2043]
                                for i in range(2044,2061)})
    # fill more recent values for S and W, largely irrelevant given this is 
    # about past years (with no DH HPs installed)
    hh_proj = hh_proj.fillna(method="backfill",axis=1)
    
    # drop years before modeling period
    hh_proj = hh_proj.drop(list(range(2001,2015)),axis=1)
    
    # calculate potential in TJ (multiply with density of water, number of days
    # per year, heat capacity of water, multiply to scale up and include also
    # the heat from the electricity use itself)
    wwpot = (hh_proj * wwpp * 0.997 * 365 * wwtdiff * 4184 / 10**12
             /(1-1/swhp_eff))
    
    # rearrange dataframe
    wwpot.columns.name="YEAR"
    wwpot = wwpot.stack()
    wwpot.name = "VALUE"
    wwpot = wwpot.reset_index()
    wwpot["TECHNOLOGY"] = raw_data.loc[raw_data["VARIABLE"]=="WasteWaterperPerson","TECHNOLOGY"].squeeze()
    wwpot["VARIABLE"] = "TotalTechnologyAnnualActivityUpperLimit"
    wwpot = wwpot.rename(columns={"LAD23CD":"REGION"})
    wwpot = wwpot.reindex(sorted(wwpot.columns), axis=1)
    wwpot = wwpot.set_index([col for col 
                             in wwpot.columns
                             if col!="VALUE"]) 

#%%
# save all

    # save to file
    # unit: TJ/a
    pd.concat([ieh,
               wwpot]).to_csv(snakemake.output.path_heat_pots)
        
 