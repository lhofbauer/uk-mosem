"""
Script deriving hourly demand timeseries


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

import datetime

import pandas as pd
import demandlib.bdew as dl

import utils


if __name__ == "__main__":
    
    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_demand_timeseries",
                                         **df_params.iloc[0].to_dict())
        
    # load temperature timeseries
    temp = pd.read_csv(snakemake.input.path_temp)
    temp = temp.set_index("LAD23CD")
    temp = temp.loc[:,str(datetime.datetime(1999, 1, 1, 0)):]
    
    
    # create dataframe for timeseries  
    ts = pd.date_range(datetime.datetime(1999, 1, 1, 0),
                            periods=30*24, freq='H')
    ts = ts.append([pd.date_range(datetime.datetime(1999, m, 1, 0),
                                  periods=30*24, freq='H') 
                    if m != 2 else
                    pd.date_range(datetime.datetime(1999, m, 1, 0),
                                  periods=28*24, freq='H')
                    for m in range(2,13) ])

    
    # create heat timeseries for each LAD using demandlib    
    
    demand_h = pd.DataFrame(index=ts)    

    for i in range(len(temp)):
        
        hb = dl.HeatBuilding(demand_h.index, temperature=temp.iloc[i],
                                      shlp_type="MFH",building_class=11,
                                      wind_class=1,
                                      name="MFH",
                                      ww_incl=False)
        demand_h[(temp.index[i],"Flats","SH")] = hb.get_normalized_bdew_profile()
        # scaling below necessary to derive HW demand by substracting (below),
        # see demand lib code for equations
        [a, b, c, d] = hb.get_sigmoid_parameters()
        f = hb.get_weekday_parameters()
        h = (a / (1 + (b / (hb.df['temperature_geo'] - 40)) ** c) + d)
        kw = 1.0 / (sum(h * f) / 24)
        demand_h[(temp.index[i],"Flats","SH")]= demand_h[(temp.index[i],"Flats","SH")]/kw    

        hb = dl.HeatBuilding(demand_h.index, temperature=temp.iloc[i],
                                      shlp_type="EFH",building_class=11,
                                      wind_class=1,
                                      name="EFH",
                                      ww_incl=False)
        demand_h[(temp.index[i],"Detached","SH")] = hb.get_normalized_bdew_profile()
        [a, b, c, d] = hb.get_sigmoid_parameters()
        f = hb.get_weekday_parameters()
        h = (a / (1 + (b / (hb.df['temperature_geo'] - 40)) ** c) + d)
        kw = 1.0 / (sum(h * f) / 24)
        demand_h[(temp.index[i],"Detached","SH")]= demand_h[(temp.index[i],"Detached","SH")]/kw          

        demand_h[(temp.index[i],"Semi-detached","SH")] = demand_h[(temp.index[i],"Detached","SH")]
        demand_h[(temp.index[i],"Terraced","SH")] = demand_h[(temp.index[i],"Detached","SH")]

        hb = dl.HeatBuilding(demand_h.index, temperature=temp.iloc[i],
                                      shlp_type="GHD",building_class=0,
                                      wind_class=1,
                                      name="GHD",
                                      ww_incl=False)
        demand_h[(temp.index[i],"Non-domestic","SH")] = hb.get_normalized_bdew_profile()
        [a, b, c, d] = hb.get_sigmoid_parameters()
        f = hb.get_weekday_parameters()
        h = (a / (1 + (b / (hb.df['temperature_geo'] - 40)) ** c) + d)
        kw = 1.0 / (sum(h * f) / 24)
        demand_h[(temp.index[i],"Non-domestic","SH")]= demand_h[(temp.index[i],"Non-domestic","SH")]/kw   
      
        hb = dl.HeatBuilding(demand_h.index, temperature=temp.iloc[i],
                                      shlp_type="MFH",building_class=11,
                                      wind_class=1,
                                      name="MFH",
                                      ww_incl=True)
        demand_h[(temp.index[i],"Flats","HW")] = hb.get_normalized_bdew_profile()
        [a, b, c, d] = hb.get_sigmoid_parameters()
        f = hb.get_weekday_parameters()
        h = (a / (1 + (b / (hb.df['temperature_geo'] - 40)) ** c) + d)
        kw = 1.0 / (sum(h * f) / 24)
        demand_h[(temp.index[i],"Flats","HW")]= (demand_h[(temp.index[i],"Flats","HW")]/kw
                                                 -demand_h[(temp.index[i],"Flats","SH")])
        
        
        hb = dl.HeatBuilding(demand_h.index, temperature=temp.iloc[i],
                                      shlp_type="EFH",building_class=11,
                                      wind_class=1,
                                      name="EFH",
                                      ww_incl=True)
        demand_h[(temp.index[i],"Detached","HW")] = hb.get_normalized_bdew_profile()
        [a, b, c, d] = hb.get_sigmoid_parameters()
        f = hb.get_weekday_parameters()
        h = (a / (1 + (b / (hb.df['temperature_geo'] - 40)) ** c) + d)
        kw = 1.0 / (sum(h * f) / 24)
        demand_h[(temp.index[i],"Detached","HW")]= (demand_h[(temp.index[i],"Detached","HW")]/kw
                                                    -demand_h[(temp.index[i],"Detached","SH")])

        demand_h[(temp.index[i],"Semi-detached","HW")] = demand_h[(temp.index[i],"Detached","HW")]
        demand_h[(temp.index[i],"Terraced","HW")] = demand_h[(temp.index[i],"Detached","HW")]

        hb = dl.HeatBuilding(demand_h.index, temperature=temp.iloc[i],
                                      shlp_type="GHD",building_class=0,
                                      wind_class=1,
                                      name="GHD",
                                      ww_incl=True)
        demand_h[(temp.index[i],"Non-domestic","HW")] = hb.get_normalized_bdew_profile()
        [a, b, c, d] = hb.get_sigmoid_parameters()
        f = hb.get_weekday_parameters()
        h = (a / (1 + (b / (hb.df['temperature_geo'] - 40)) ** c) + d)
        kw = 1.0 / (sum(h * f) / 24)
        demand_h[(temp.index[i],"Non-domestic","HW")]= (demand_h[(temp.index[i],"Non-domestic","HW")]/kw
                                                        -demand_h[(temp.index[i],"Non-domestic","SH")])
        
         
     
    demand_h = demand_h.T
    demand_h.index = pd.MultiIndex.from_tuples(demand_h.index).rename(["LAD23CD",
                                                                   "PROPERTY_TYPE",
                                                                   "DEMAND_TYPE"])
    # normalize
    demand_h = demand_h.divide(demand_h.sum(axis=1),axis=0) *1000
    
    # create NHE timeseries for each LAD using demandlib  
    
    demand_nhe = pd.DataFrame(index=ts)    
    dom_nhe = dl.ElecSlp(year=1999).get_profile({'h0': 1000},
                                                dyn_function_h0=True).resample("H").mean()
    ndom_nhe = dl.ElecSlp(year=1999).get_profile({'g0': 1000},
                                                 dyn_function_h0=True).resample("H").mean()
    for i in range(len(temp)):
        demand_nhe[(temp.index[i],"Flats")] = dom_nhe
        demand_nhe[(temp.index[i],"Detached")] = dom_nhe
        demand_nhe[(temp.index[i],"Semi-detached")] = dom_nhe
        demand_nhe[(temp.index[i],"Terraced")] = dom_nhe
        demand_nhe[(temp.index[i],"Non-domestic")] = ndom_nhe
        
     
    demand_nhe = demand_nhe.T
    demand_nhe.index = pd.MultiIndex.from_tuples(demand_nhe.index).rename(["LAD23CD",
                                                                   "PROPERTY_TYPE"])   
    
    # create NHG timeseries for each LAD assuming constant demand
    
    demand_nhg = pd.DataFrame(index=ts)    

    for i in range(len(temp)):
        demand_nhg[(temp.index[i],"Flats")] = 1/len(ts)*1000
        demand_nhg[(temp.index[i],"Detached")] = 1/len(ts)*1000
        demand_nhg[(temp.index[i],"Semi-detached")] = 1/len(ts)*1000
        demand_nhg[(temp.index[i],"Terraced")] = 1/len(ts)*1000
        demand_nhg[(temp.index[i],"Non-domestic")] = 1/len(ts)*1000
              
        
     
    demand_nhg = demand_nhg.T
    demand_nhg.index = pd.MultiIndex.from_tuples(demand_nhg.index).rename(["LAD23CD",
                                                                   "PROPERTY_TYPE"])       
    
    # save timeseries
    demand_nhg.to_csv(snakemake.output.path_dem_tseries_nhg)
    demand_nhe.to_csv(snakemake.output.path_dem_tseries_nhe)
    demand_h.xs("SH",level=2).to_csv(snakemake.output.path_dem_tseries_sh)
    demand_h.xs("HW",level=2).to_csv(snakemake.output.path_dem_tseries_hw)
    demand_h.groupby(["LAD23CD",
                      "PROPERTY_TYPE"]).sum().to_csv(snakemake.output.path_dem_tseries_shw)
 
