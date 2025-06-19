"""
Script deriving timeseries from reanalysis ERA5 data 


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

import pandas as pd
import geopandas as gpd

import atlite

import utils

import logging

logger = logging.getLogger(__name__)
utils.adjust_logger()


if __name__ == "__main__":

    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_reanalysis_timeseries",
                                         **df_params.iloc[0].to_dict())
        
    # load reanalysis data
    cutout = atlite.Cutout(snakemake.input.path_era5)


    gb_shp = gpd.read_file(snakemake.input.path_ctry_bounds)
    gb_shp = gb_shp.dissolve()
    gb_shp = gb_shp.to_crs("epsg:4326")
    
    # load shapefiles
    offw_shp_EWNI = gpd.read_file(snakemake.input.path_offw_ewni_bounds)
    offw_shp_S = gpd.read_file(snakemake.input.path_offw_s_bounds)
    offw_shp = pd.concat([offw_shp_EWNI,offw_shp_S])
    offw_shp = offw_shp.dissolve()
    
    lads = gpd.read_file(snakemake.input.path_lad_bounds)
    lads = lads.to_crs("epsg:4326")
    
    # calculate wind (off- and onshore), and PV capacity factor timeseries
    # FIXME: currently this is derived as one UK timeseries, needs to be
    # changed to LAD if power sector is disaggregated (i.e., this should
    # be based on an option)
    wind_gen,wind_cap = cutout.wind(turbine="Vestas_V90_3MW",
                                   shapes=gb_shp,
                                   index=pd.Index(["UK"],name="GEOGRAPHY"),
                                   #index=pd.Index(gb_shp["CTRY22CD"],name="GEOGRAPHY"),
                                   return_capacity=True)
    wind_gen = wind_gen.to_pandas()
    wind_cap = wind_cap.to_pandas()
    wind_on_cf = wind_gen/wind_cap
    
    wind_gen,wind_cap = cutout.wind(turbine="Enercon_E126_7500kW",
                                   shapes=offw_shp,
                                   index=pd.Index(["UK"],name="GEOGRAPHY"),
                                   return_capacity=True)
    wind_gen = wind_gen.to_pandas()
    wind_cap = wind_cap.to_pandas()
    wind_off_cf = wind_gen/wind_cap

    pv_gen,pv_cap = cutout.pv(panel="CSi",
                              shapes=gb_shp,
                              orientation="latitude_optimal",
                              index=pd.Index(["UK"],name="GEOGRAPHY"),
                              return_capacity=True)
    pv_gen = pv_gen.to_pandas()
    pv_cap = pv_cap.to_pandas()
    pv_cf = pv_gen/pv_cap
    
    # these temperature data are currently not used
    # MetOffice data are likely more accurate
    temp = cutout.temperature(shapes=lads,
                              index=pd.Index(lads["LAD23CD"],
                                             name="GEOGRAPHY"),
                              per_unit=True).to_pandas()

    
    ts = pd.concat([wind_on_cf,
                    wind_off_cf,
                    pv_cf,
                    pv_cf,
                    temp],
                   keys=["OWPPSNAT00","FWPPSNAT00","OSPPSNAT00","RSPPSNAT00",
                         "TEMPERATURE"],
                   names=["VARIABLE"],
                   axis=1)
    
    # save to file
    # unit: -, except °C for temperatures
    ts.to_csv(snakemake.output.path_ra_timeseries)
    


                                                    