"""
Script deriving aggregated timeseries, i.e., time slice values and lengths


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

import pandas as pd

import tsam.timeseriesaggregation as tsam

import utils


if __name__ == "__main__":
    
    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_aggregated_timeseries",
                                         **df_params.iloc[4].to_dict())
  
        
    # load demand
    demand_sh = pd.read_csv(snakemake.input.path_dem_tseries_sh, index_col=["LAD23CD",
                                                            "PROPERTY_TYPE"])
    demand_hw = pd.read_csv(snakemake.input.path_dem_tseries_hw, index_col=["LAD23CD",
                                                            "PROPERTY_TYPE"])
    demand_shw = pd.read_csv(snakemake.input.path_dem_tseries_shw, index_col=["LAD23CD",
                                                            "PROPERTY_TYPE"])
    demand_nhe = pd.read_csv(snakemake.input.path_dem_tseries_nhe, index_col=["LAD23CD",
                                                            "PROPERTY_TYPE"])
    demand_nhg = pd.read_csv(snakemake.input.path_dem_tseries_nhg, index_col=["LAD23CD",
                                                            "PROPERTY_TYPE"])
    

    
    demand = pd.concat([demand_sh,demand_hw,demand_shw,demand_nhe,demand_nhg],
                       keys=["SH","HW","SHW","NHE","NHG"])
    ts = demand.T
    
    # add one column to represent average/national profile that is used to
    # pick additional peak period/day
    # this could be implemented as weighted average to capture the national
    # peak better, but might then capture local peaks less well
    ts[("_SHW","_UK","_All")] = ts.xs("SHW",
                                      level=0,
                                      axis=1).mean(axis=1)
    
    # add supply timeseries if required and set relevant parameters
    if snakemake.params.dic["scen_supply_imp"]=="agg":
        
        # set weights so aggregation is mainly based on SHW profiles
        weights = {c:(1000/ts.columns.get_level_values(0).value_counts()["SHW"]
                      if (c[0]=="SHW")
                      else 0.01) for c in ts.columns}
      
    elif snakemake.params.dic["scen_supply_imp"].startswith("det"):

        # load reanalysis timeseries
        ra = pd.read_csv(snakemake.input.path_ra_timeseries,
                        index_col=[0],header=[0,1])
        ra = pd.concat([ra.loc[:,["OWPPSNAT00","FWPPSNAT00",
                                  "OSPPSNAT00","RSPPSNAT00"]]],
                       keys=["Supply"],
                       names=[""],
                       axis=1)

        ts = pd.concat([ts,ra],axis=1)
        # The temperature data and thus demand timeseries are based on 30 day 
        # months (except February) while the supply timeseries can have months
        # with 31 days. Thus some nan for 31st's are dropped here.
        ts = ts.dropna()
        # set weights so aggregation is mainly based on SHW profiles
        weights = {c:(20000/ts.columns.get_level_values(0).value_counts()["SHW"]
                      if (c[0]=="SHW") else
                      10/ts.columns.get_level_values(1).value_counts()["OWPPSNAT00"]
                      if (c[1]=="OWPPSNAT00") else
                      10/ts.columns.get_level_values(1).value_counts()["FWPPSNAT00"]
                      if (c[1]=="FWPPSNAT00") else
                      10/ts.columns.get_level_values(1).value_counts()["OSPPSNAT00"]
                      if (c[1]=="OSPPSNAT00") else
                      10/ts.columns.get_level_values(1).value_counts()["RSPPSNAT00"]
                      if (c[1]=="RSPPSNAT00")
                      else 0.01) for c in ts.columns}
                         
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_supply_imp']}' is currently"
                " not a valid option of the 'scen_supply_imp' parameter."
                " Valid options are 'agg' or 'det'."
            )      
        
 
    
    # set relevant parameters for the time aggregation
    pp,seg = snakemake.params.dic["scen_time_agg"].split("|")
    
    if pp[-1]=="p":
        hpp = 24
        po = [0]*58 + [1]*90 + [1]*90 + [1]*90 + [0]*30
        ntp = 2
    elif pp[-1]=="d":
        hpp = 24
        po = None
    elif pp[-1]=="h":
        hpp = 1
        po = None
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_time_agg']}' is currently"
                " not a valid option of the 'scen_time_agg' parameter."
                " Valid options are ending on 'd' or 'h' before the '|'."
            )
    if pp[:-1].isdigit() and pp[-1]!="p":
        ntp = int(pp[:-1])
    elif  pp[-1]=="p":
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_time_agg']}' is currently"
                " not a valid option of the 'scen_time_agg' parameter."
                " The first part of valid options have to start with an integer"
                " before ending with one letter before the '|'."
                )
    if seg.isdigit():
        ns = int(seg)
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_time_agg']}' is currently"
                " not a valid option of the 'scen_time_agg' parameter."
                " The last part of valid options has to be an integer."
                )


    # aggregate timeseries based on pre-defined day order (no clustering
    # because chronological timeslices are necessary if storage is to be used)
    # while adding the peak day
    aggregation = tsam.TimeSeriesAggregation(ts, noTypicalPeriods = ntp,
                                             hoursPerPeriod = hpp,
                                             clusterMethod = 'k_means',
                                             predefClusterOrder=po,
                                             segmentation=True,
                                             weightDict=weights,
                                             noSegments=ns,
                                             extremePeriodMethod="new_cluster_center",
                                             addPeakMax=[("_SHW","_UK","_All")])
    
    typPeriods = aggregation.createTypicalPeriods()
    
    typPeriods = typPeriods.drop(("_SHW","_UK","_All"),axis=1)
    
    # split supply and demand columns
    typPeriods_dem = typPeriods.loc[:,["SH","HW","SHW","NHE","NHG"]]

    # calculate normalized demand per segment (taking into account
    # the number of times each period is present per year as well as the
    # length of each segment)
    typPeriodsnum = aggregation.clusterPeriodNoOccur
    
    se = pd.Series(typPeriodsnum)
    se = se.rename_axis("period")
    
    typPeriods_dem_scaled = (typPeriods_dem.T
                             * typPeriods_dem.index.get_level_values("Segment Duration")
                             * se.repeat(ns).values)
    typPeriods_dem_scaled = typPeriods_dem_scaled.T
    
    typPeriods_dem_norm = typPeriods_dem_scaled / typPeriods_dem_scaled.sum(axis=0)
    typPeriods_dem_norm = typPeriods_dem_norm.reset_index(drop=True).T
    
    # calculate the length of each 'kind of' segment per year
    typPeriodslength = (typPeriods_dem.index.get_level_values("Segment Duration")
                         * se.repeat(ns))
    typPeriodslength_norm = (typPeriodslength.reset_index(drop=True)
                             / typPeriodslength.sum())
    typPeriodslength_norm.name="VALUE"
    
    # rearrange supply capacity factor dataframe, create empty dataframe if
    # not used
    if "Supply" in typPeriods.columns:
        typPeriods_sup =  typPeriods.xs("Supply",level=0,axis=1)
        typPeriods_sup = typPeriods_sup.reset_index(drop=True).T
        typPeriods_sup.columns = ["TS"+ str(i) 
                                  for i in typPeriods_sup.columns.to_list()]
        typPeriods_sup.index.names = ["TECHNOLOGY","REGION"]
        
    else:
        typPeriods_sup = pd.DataFrame()
    
    typPeriods_dem_norm.columns = ["TS"+ str(i) 
                                for i in typPeriods_dem_norm.columns.to_list()]

    typPeriodslength_norm.index = ["TS"+ str(i)
                                   for i 
                                   in typPeriodslength_norm.index.to_list()]
    

    typPeriodslength_norm.index.name="TIMESLICE"
    

    # explode property types if household disaggregated
    if snakemake.params.dic["scen_hh_disagg"]=="T":
        
        typPeriods_dem_norm = utils.explode_hh(typPeriods_dem_norm,
                                           ["OO","RP","RS"])
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        typPeriods_dem_norm = utils.explode_hh(typPeriods_dem_norm,
                                  [t+str(i) for t in ["O","P","S"]
                                   for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)]
                                  )
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )
    # split typperiods dataframe based on demand type to save in separate files
    
    typPeriods_norm_SH = typPeriods_dem_norm.loc[["SH"],:].droplevel(0)
    typPeriods_norm_HW = typPeriods_dem_norm.loc[["HW"],:].droplevel(0)
    typPeriods_norm_SHW = typPeriods_dem_norm.loc[["SHW"],:].droplevel(0)
    typPeriods_norm_NHE = typPeriods_dem_norm.loc[["NHE"],:].droplevel(0)
    typPeriods_norm_NHG = typPeriods_dem_norm.loc[["NHG"],:].droplevel(0)
    
    # FIXME: also create data on structure of timeslices that feeds into
    # conversionld, etc. parameters

    # save to files
    typPeriods_norm_SH.to_csv(snakemake.output.path_tperiods_sh)
    typPeriods_norm_HW.to_csv(snakemake.output.path_tperiods_hw)
    typPeriods_norm_SHW.to_csv(snakemake.output.path_tperiods_shw)
    typPeriods_norm_NHE.to_csv(snakemake.output.path_tperiods_nhe)
    typPeriods_norm_NHG.to_csv(snakemake.output.path_tperiods_nhg)
    
    typPeriods_sup.to_csv(snakemake.output.path_tperiods_supply)
    
    typPeriodslength_norm.to_csv(snakemake.output.path_tperiods_len)

    
    

