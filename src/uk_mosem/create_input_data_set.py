"""
Script processing raw and processed data to create the OSeMOSYS input data set


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

import os
import shutil
import time
import logging

import numpy as np
import pandas as pd
import frictionless as fl

import utils

DECIMALS=7

logger = logging.getLogger(__name__)
utils.adjust_logger()

if __name__ == "__main__":
    
    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("create_input_data_set",
                                         **df_params.iloc[45].to_dict())    

    # aggregating script creating all param and set files
    
    # FIXME: in several instances, consider if techs are only defined at the
    # right scale if it is necessary or an issue to be addressed elsewhere
#%% 
# calculate global levers

    # set multi-year aggregation
    if snakemake.params.dic["scen_year_agg"]=="5y":
        years_map = pd.Series([2015]*6+[2021]*2+[2023]*2+
                              [e for e in list(range(2025,2056,5))
                               for i in range(5)]+
                              [2060]*1,
                              index=range(2015,2061))
        years_map.to_csv(snakemake.output.path_years_map, index=True)
      
    elif snakemake.params.dic["scen_year_agg"]=="-":

        years_map = None                             

    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_year_agg']}' is currently"
                " not a valid option of the 'scen_year_agg' parameter."
                " Valid options are '5y' or '-'."
            )      
        

#%% 
# create directory

    if os.path.exists(snakemake.params.fdir
                      +snakemake.params.dic["name"]):
        shutil.rmtree(snakemake.params.fdir
                      +snakemake.params.dic["name"])
    time.sleep(2)    
    os.makedirs(snakemake.params.fdir
                +snakemake.params.dic["name"]+"/data/")

#%%
# operational life
      
    # simply get operational life times from raw input values and save to csv
    # unit: a
    
    # load raw data set
    raw_data = pd.read_csv(snakemake.input.path_set_lt)
    
    # if list of technologies given for rows, create separate rows
    raw_data = raw_data.assign(TECHNOLOGY=raw_data['TECHNOLOGY'].str.split(',')).explode('TECHNOLOGY')

    
    # load building measure lifetime and concatenate
    bm_lt = pd.read_csv(snakemake.input.path_building_measures_lt)
    supply_lt = pd.read_csv(snakemake.input.path_supply_lt)
    oplife = pd.concat([raw_data,supply_lt,bm_lt])     
    
    oplife = oplife[["REGION","TECHNOLOGY",
                         "VALUE"]][oplife["VARIABLE"]=="OperationalLife"]

  

    oplife["REGION"][oplife["REGION"].isnull()] = ":*"
    

    # explode technologies if household disaggregated
    if snakemake.params.dic["scen_hh_disagg"]=="T":
        
        oplife = utils.explode_hh(oplife,
                                  ["OO","RP","RS"],
                                  col="TECHNOLOGY")
        
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        oplife = utils.explode_hh(oplife,
                                  [t+str(i) for t in ["O","P","S"]
                                   for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )        
    
    oplife["VALUE"] = oplife["VALUE"].round(DECIMALS)
    oplife.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"OperationalLife.csv", index=False)


#%%
# capital cost
    
    # simply get capital costs from raw input values and save to csv
    # unit: million [Baseyear]£/GW (= [Baseyear]£/kW)
    
    # load raw data set
    raw_data = pd.read_csv(snakemake.input.path_set_costs)
    sets = pd.read_csv(snakemake.input.path_set_sets)
    years = sets["VALUE"][sets["SET"]=="YEAR"].astype(int)
    
    # load local and building measures cost data and concatenate
    local_cost = pd.read_csv(snakemake.input.path_local_cost)
    bm_cost = pd.read_csv(snakemake.input.path_building_measures_cost)
    supply_cost = pd.read_csv(snakemake.input.path_supply_cost)
    cost_data = pd.concat([raw_data,local_cost,supply_cost,bm_cost])
    
    
    capcost = cost_data[["REGION","TECHNOLOGY","YEAR","UNIT",
                         "VALUE"]][cost_data["VARIABLE"]=="CapitalCost"]
    capcost.loc[:,"REGION"] = capcost["REGION"].fillna(value=":*")
    
    # adjust cost values for single base year using the producer price index
    capcost = utils.adjust_monetary_values(capcost,2015,"UNIT")
    capcost = capcost.drop("UNIT",axis=1)

    
    # if no year given, assume data is the same for all years, interpolate, if 
    # necessary, remaining values
    capcost.loc[:,"YEAR"]= capcost["YEAR"].fillna(value=":*")
    capcost = utils.interpolate_timeseries(capcost,model_periods=years,
                                           force_exp=True,
                                           dropna=True)
    
    
    # if list of technologies given for rows, create separate rows
    capcost = capcost.assign(TECHNOLOGY=capcost['TECHNOLOGY'].str.split(',')).explode('TECHNOLOGY')
    
    
    # explode technologies if household disaggregated
    if snakemake.params.dic["scen_hh_disagg"]=="T":
        
        capcost = utils.explode_hh(capcost,
                                  ["OO","RP","RS"],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        capcost = utils.explode_hh(capcost,
                                  [t+str(i) for t in ["O","P","S"]
                                   for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )  
        
        
    # adjust cost based on scenario cost lever, if triggered
    capcost = capcost.set_index(["REGION","TECHNOLOGY","YEAR"])

    if snakemake.params.dic["scen_tech_cost"]=="-":
        pass
    
    else:
        techl = sets.loc[(sets["SET"]=="TECHNOLOGY"),
                         "VALUE"]
        techl.name = "TECHNOLOGY"
        techl = techl.to_frame()
        years = sets.loc[(sets["SET"]=="YEAR"),
                         "VALUE"].to_frame()
        slagg = pd.read_csv(snakemake.input.path_sublocal_agg,
                                index_col=["LSOA11CD"]).squeeze()
    
        if snakemake.params.dic["scen_hh_disagg"]=="T":
            techl = utils.explode_hh(techl,
                                  ["OO","RP","RS"],
                                  col="TECHNOLOGY")
        elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
            techl = utils.explode_hh(techl,
                                      [t+str(i) for t in ["O","P","S"]
                                       for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                      col="TECHNOLOGY")    
        elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
            pass
        else:
            raise NotImplementedError(
                    f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                    " not a valid option of the 'scen_hh_disagg' parameter."
                    " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
                )
        
        adj = snakemake.params.dic["scen_tech_cost"]
       
        adj = adj.split("|")
        
        uk = ";".join(utils.get_entity_lookup(["LAD"])["LAD23CD"].to_list())
        ic = "S12000014;S12000040;S12000036;E06000002;E06000004;E06000003;E08000011;E08000013;E08000014;E08000015;E08000012;E06000013;E06000012;E06000010;E06000011;E06000045;E07000086;E07000087;E07000088;E06000044;W06000009;W06000011;W06000012;W06000014;W06000015;W06000022;W06000019"
        adj = [c.replace("-IC-","-"+ic+"-") for c in adj]
        adj = [c.replace("-UK-","-"+uk+"-") for c in adj]

        for ic, a in enumerate(adj):

            tech,ct,m,rg,vys = a.split("-")
            techs = pd.DataFrame([t for t in techl["TECHNOLOGY"]
                                  if tech in t],columns=["TECHNOLOGY"])
            if ct != "C":
                continue
            # FIXME: currently always applies to all regions
            # add region, explode, and disaggregate to sublocal agg if needed
            techs["LAD23CD"] = rg
            techs = techs.assign(LAD23CD=techs["LAD23CD"].str.split(
                                        ';')).explode('LAD23CD')
            techs = techs.set_index(["LAD23CD","TECHNOLOGY"])
            # if rg != ':*':
            #     techs["VALUE"] = 0
            #     techs = utils.groupby_LAD(techs,disagg=slagg).mean()
            #     techs = techs.drop("VALUE",axis=1)
            techs = techs.reset_index()
            # techs = techs.rename(columns={slagg.name:"REGION",
            #                               "LAD23CD":"REGION"})
            
            lim = years.copy()
            lim = lim.rename(columns={"VALUE":"YEAR"})
            lim["VALUE"] = np.nan
            vys = vys.split("a")
            for  vy in vys:
                value, year = vy.split("y")
                lim.loc[lim["YEAR"]==year,"VALUE"]=value
            lim["VALUE"] = lim["VALUE"].astype(float)
            # interpolate
            lim["VALUE"] = lim["VALUE"].interpolate(limit_area="inside")
            lim["YEAR"] = lim["YEAR"].astype("int")
            lim = lim.set_index(["YEAR"]).dropna()
            
            lim.loc[:,"TECHNOLOGY"] = ",".join(techs["TECHNOLOGY"].to_list())
            lim = lim.assign(TECHNOLOGY=lim["TECHNOLOGY"].str.split(
                                        ',')).explode("TECHNOLOGY")
            lim = lim.set_index("TECHNOLOGY",append=True)
            lim = lim.reorder_levels(["TECHNOLOGY","YEAR"])
            if m=="M":
                capcost = capcost.multiply(lim,
                                           fill_value=1.0,
                                           axis=0)
            elif m=="A":
                capcost = capcost["VALUE"].add(lim,
                                           fill_value=0.0,
                                           axis=0)
    
    capcost = capcost.reset_index()
    capcost = capcost.loc[:,["REGION","TECHNOLOGY","YEAR","VALUE"]]
      
     
    # FIXME: add calculations for storage techs
    capcoststo = pd.DataFrame( columns=["REGION","STORAGE","YEAR","VALUE"])
    
    capcost = utils.aggregate_years(capcost, years_map, agg_method="mean")
    capcost["VALUE"] = capcost["VALUE"].round(DECIMALS)
    capcost.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"CapitalCost.csv", index=False)
    capcoststo = utils.aggregate_years(capcoststo, years_map, agg_method="mean")
    capcoststo["VALUE"] = capcoststo["VALUE"].round(DECIMALS)
    capcoststo.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"CapitalCostStorage.csv", index=False)

#%%
# fix cost
        
    # simply get fixed costs from raw input values and save to csv
    # unit: million [Baseyear]£/GW/a  
    
    # load raw data set
    raw_data = pd.read_csv(snakemake.input.path_set_costs)
    sets = pd.read_csv(snakemake.input.path_set_sets)
    years = sets["VALUE"][sets["SET"]=="YEAR"].astype(int)
    
    # load local processed cost data and concatenate
    local_cost = pd.read_csv(snakemake.input.path_local_cost)
    supply_cost = pd.read_csv(snakemake.input.path_supply_cost)
    cost_data = pd.concat([raw_data,supply_cost,local_cost])
    
    
    fixcost = cost_data[["REGION","TECHNOLOGY","YEAR","UNIT",
                         "VALUE"]][cost_data["VARIABLE"]=="FixedCost"]
    fixcost["REGION"][fixcost["REGION"].isnull()] = ":*"
    
    # adjust cost values for single base year using the producer price index
    fixcost = utils.adjust_monetary_values(fixcost,2015,"UNIT")
    fixcost = fixcost.drop("UNIT",axis=1)

    # if no year given, assume data is the same for all years, interpolate, if 
    # necessary, remaining values
    fixcost["YEAR"]= fixcost["YEAR"].fillna(value=":*")
    fixcost = utils.interpolate_timeseries(fixcost,model_periods=years,
                                           dropna=True)
    
    # if list of technologies given for rows, create separate rows
    fixcost = fixcost.assign(TECHNOLOGY=fixcost['TECHNOLOGY'].str.split(',')).explode('TECHNOLOGY')
    
    # explode technologies if household disaggregated
    if snakemake.params.dic["scen_hh_disagg"]=="T":
        
        fixcost = utils.explode_hh(fixcost,
                                  ["OO","RP","RS"],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        fixcost = utils.explode_hh(fixcost,
                                  [t+str(i) for t in ["O","P","S"]
                                   for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )       
    fixcost = utils.aggregate_years(fixcost, years_map, agg_method="sum")
    fixcost["VALUE"] = fixcost["VALUE"].round(DECIMALS)
    fixcost.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"FixedCost.csv", index=False)
    
#%%
# variable cost
          
    # simply get variable costs from raw input values and save to csv
    # unit: million [Baseyear]£/TJ/a

    # load raw data set
    raw_data = pd.read_csv(snakemake.input.path_set_costs)
    sets = pd.read_csv(snakemake.input.path_set_sets)
    years = sets["VALUE"][sets["SET"]=="YEAR"].astype(int)
    
    # load local processed cost data and concatenate
    local_cost = pd.read_csv(snakemake.input.path_local_cost)
    lsupply_cost = pd.read_csv(snakemake.input.path_lsupply_cost)
    supply_cost = pd.read_csv(snakemake.input.path_supply_cost)
    cost_data = pd.concat([raw_data,lsupply_cost,supply_cost,local_cost])  
    

    varcost = cost_data[["REGION","TECHNOLOGY","YEAR","MODE_OF_OPERATION",
                         "UNIT","VALUE"]][cost_data["VARIABLE"]=="VariableCost"]
    varcost["REGION"][varcost["REGION"].isnull()] = ":*"
    

    # adjust cost values for single base year and in GBP using the PPI and exchange rate
    varcost = utils.adjust_monetary_values(varcost,2015,"UNIT")
    varcost = varcost.drop("UNIT",axis=1)    

    # if no year given, assume data is the same for all years, interpolate, if 
    # necessary, remaining values
    varcost["YEAR"]= varcost["YEAR"].fillna(value=":*")
    varcost = utils.interpolate_timeseries(varcost,model_periods=years,
                                           dropna=True)

    
    # if list of technologies given for rows, create separate rows
    varcost = varcost.assign(TECHNOLOGY=varcost['TECHNOLOGY'].str.split(',')).explode('TECHNOLOGY')
    
    # if list of modes given for rows, create separate rows
    # varcost = varcost.assign(TECHNOLOGY=varcost["MODE_OF_OPERATION"].str.split(',')).explode("MODE_OF_OPERATION")
    
    # explode technologies if household disaggregated
    if snakemake.params.dic["scen_hh_disagg"]=="T":
        
        varcost = utils.explode_hh(varcost,
                                  ["OO","RP","RS"],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        varcost = utils.explode_hh(varcost,
                                  [t+str(i) for t in ["O","P","S"]
                                   for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )  

    varcost = varcost.reindex(["REGION","TECHNOLOGY",
                                 "MODE_OF_OPERATION","YEAR","VALUE"],axis=1)
    varcost = utils.aggregate_years(varcost, years_map, agg_method="mean")
    varcost["VALUE"] = varcost["VALUE"].round(DECIMALS)
    varcost.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"VariableCost.csv", index=False)
    
#%%
# input and output ratios

    # get efficiency and derive in- and output activity ratios and save to csv
    # unit: -
    
    # FIXME: make this more streamlined and clear, while being flexible
    
    # load raw data set
    raw_data = pd.read_csv(snakemake.input.path_set_eff)
    raw_data["YEAR"] = raw_data["YEAR"].astype(object)
    
    sets = pd.read_csv(snakemake.input.path_set_sets)
    years = sets["VALUE"][sets["SET"]=="YEAR"].astype(int)
    
    # load local processed efficiency data and concatenate
    local_eff = pd.read_csv(snakemake.input.path_local_eff)
    supply_eff = pd.read_csv(snakemake.input.path_supply_eff)
    eff_data = pd.concat([raw_data,supply_eff,local_eff])  
    
    inratio = eff_data[["REGION","TECHNOLOGY","FUEL_IN","MODE_OF_OPERATION",
                          "YEAR","VALUE"]][(eff_data["VARIABLE"]=="Efficiency")
                                           & (eff_data["FUEL_IN"].notnull())] 
                                           
    # set input activity ration based on efficiency, except efficiency is 0 
    # (in which case the ratio is 0)
    inratio["VALUE"][inratio["VALUE"]!=0] = 1/(inratio["VALUE"])
                                           
    # allow also for input ratio to be given directly                                       
    indratio = eff_data[["REGION","TECHNOLOGY","FUEL_IN","MODE_OF_OPERATION",
                              "YEAR","VALUE"]][(eff_data["VARIABLE"]=="InputActivityRatio")
                                               & (eff_data["FUEL_IN"].notnull())]
    inratio =  pd.concat([indratio,inratio])
    
    # process
    inratio = inratio.rename(columns={"FUEL_IN":"FUEL"})
    
    inratio.loc[inratio["REGION"].isnull(),"REGION"] = ":*"
    
    # FIXME: introduce a methodology for changing efficiency
    inratio.loc[inratio["YEAR"].isnull(),"YEAR"] = ":*"
    
    # if list of technologies given for rows, create separate rows
    inratiot = inratio.loc[~inratio["FUEL"].str.contains(",")].set_index(["REGION","FUEL","MODE_OF_OPERATION","YEAR",
                                  "VALUE"])
    inratiot = inratiot.assign(TECHNOLOGY=inratiot['TECHNOLOGY'].str.split(','))
    inratiot = inratiot.apply(pd.Series.explode).reset_index()  
    
    # if list of technologies & fuels given for rows, create separate rows
    # (must be the same number of techs and fuels in each row)
    inratiotf = inratio.loc[inratio["FUEL"].str.contains(",")].set_index(["REGION","MODE_OF_OPERATION","YEAR",
                                 "VALUE"]) 
    inratiotf = inratiotf.assign(TECHNOLOGY=inratiotf['TECHNOLOGY'].str.split(','),
                               FUEL=inratiotf['FUEL'].str.split(','))
    inratiotf = inratiotf.apply(pd.Series.explode).reset_index() 
    
    inratio =  pd.concat([inratiot,inratiotf])
    # sort columns
    inratio = inratio.reindex(["REGION","TECHNOLOGY","FUEL",
                                 "MODE_OF_OPERATION","YEAR","VALUE"],axis=1)
    
    # drop duplicates, these should only come from technologies that have two
    # or more fuels as output (given the input fuel column is dropped),
    # which is accounted for below for the output ratio. the input ratio is
    # only defined by the first efficiency given, the second will result in a
    # second, adjusted, i.e., in general not =1, output ratio, i.e., capacity
    # cost need to be defined with respect to the output from the first
    # efficiency given. For technologies with multiple input fuels, respective
    # 'efficiencies' need to be defined as reverse of the intended input ratio
    
    inratio = inratio[~inratio[["REGION","TECHNOLOGY","FUEL",
                                 "MODE_OF_OPERATION",
                                 "YEAR"]].duplicated(keep="first")]

    # remove trailing spaces
    inratio = inratio.apply(lambda x: x.astype(str).str.strip() if x.dtypes=="object" else x)
    
    # interpolate if necessary
    inratio = utils.interpolate_timeseries(inratio,model_periods=years,dropna=True)
    
    # explode technologies if household disaggregated
    if snakemake.params.dic["scen_hh_disagg"]=="T":
        
        inratio = utils.explode_hh(inratio,
                                  ["OO","RP","RS"],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        inratio = utils.explode_hh(inratio,
                                  [t+str(i) for t in ["O","P","S"]
                                   for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )  

    inratio = utils.aggregate_years(inratio, years_map, agg_method="mean")
    inratio["VALUE"] = inratio["VALUE"].round(DECIMALS)
    inratio.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"InputActivityRatio.csv", index=False)
    
    
    outratio = eff_data[["REGION","TECHNOLOGY","FUEL_OUT","MODE_OF_OPERATION",
                           "YEAR","VALUE"]][(eff_data["VARIABLE"]=="Efficiency")
                                            & (eff_data["FUEL_OUT"].notnull())]
                                            

    
    outratio = outratio.rename(columns={"FUEL_OUT":"FUEL"})
    
    outratio["REGION"][outratio["REGION"].isnull()] = ":*"
    
    # FIXME: introduce some methodology for changing efficiency
    outratio["YEAR"][outratio["YEAR"].isnull()] = ":*"
    
    # if list of technologies & fuels given for rows, create separate rows
    # (must be the same number of techs and fuels in each row)
    outratio = outratio.set_index(["REGION","MODE_OF_OPERATION","YEAR",
                                 "VALUE"],append=False) 
    outratio = outratio.assign(TECHNOLOGY=outratio['TECHNOLOGY'].str.split(','),
                               FUEL=outratio['FUEL'].str.split(','))
    outratio = outratio.apply(pd.Series.explode).reset_index()
    
    
    
    # remove duplicated rows (different input fuels but same output)
    
    outratio = outratio[~outratio[["REGION","TECHNOLOGY","FUEL",
                                 "MODE_OF_OPERATION",
                                 "YEAR"]].duplicated(keep="first")]
    
    
    # calculate output ratio if more than one output fuel
    outratio = outratio.set_index(["REGION","TECHNOLOGY",
                                   "MODE_OF_OPERATION",
                                   "YEAR"],append=False)
    # get only the rows where technologies have two or more output fuel,
    # as well as only the first rows of those
    outratio_dup = outratio[outratio.index.duplicated(keep=False)]
    outratio_dup_f = outratio_dup[~outratio_dup.index.duplicated(keep='first')]
    
    # calculate output ratios by dividing all values of one technology by the
    # one given first
    outratio_dup["VALUE"] = outratio_dup["VALUE"]/outratio_dup_f["VALUE"]
    
    
    
    
    # set output activity ration to 1 for remaining non-duplicated indices,
    # except efficiency is 0 (in which case the ratio is 0)
    outratio = outratio[~outratio.index.duplicated(keep=False)]
    outratio["VALUE"][outratio["VALUE"]!=0]=1
    
    
    # allow also for output ratio to be given directly                                       
    outdratio = eff_data[["REGION","TECHNOLOGY","FUEL_OUT","MODE_OF_OPERATION",
                              "YEAR","VALUE"]][(eff_data["VARIABLE"]=="OutputActivityRatio")
                                               & (eff_data["FUEL_OUT"].notnull())]
    
    outdratio = outdratio.rename(columns={"FUEL_OUT":"FUEL"})
    
    outdratio["REGION"][outdratio["REGION"].isnull()] = ":*"
    
    # FIXME: introduce some methodology for changing efficiency
    outdratio["YEAR"][outdratio["YEAR"].isnull()] = ":*"
    
    # if list of technologies & fuels given for rows, create separate rows
    # (must be the same number of techs and fuels in each row)
    outdratio = outdratio.set_index(["REGION","MODE_OF_OPERATION","YEAR",
                                 "VALUE"],append=False) 
    outdratio = outdratio.assign(TECHNOLOGY=outdratio['TECHNOLOGY'].str.split(','),
                               FUEL=outdratio['FUEL'].str.split(','))
    outdratio = outdratio.apply(pd.Series.explode).reset_index()   
    outdratio = outdratio.set_index(["REGION","TECHNOLOGY",
                                   "MODE_OF_OPERATION",
                                   "YEAR"],append=False)
    # concat output ratios
    outratio = pd.concat([outratio_dup,outdratio,outratio]).reset_index()
    
    # sort columns
    outratio = outratio.reindex(["REGION","TECHNOLOGY","FUEL",
                                 "MODE_OF_OPERATION","YEAR","VALUE"],axis=1)
    
    # delete accidental trailing spaces
    outratio = outratio.apply(lambda x: x.astype(str).str.strip() if x.dtypes=="object" else x)
    
    # interpolate if necessary
    outratio = utils.interpolate_timeseries(outratio,model_periods=years)  
    outratio = outratio.dropna()
    
    # explode technologies if household disaggregated
    if snakemake.params.dic["scen_hh_disagg"]=="T":
        
        outratio = utils.explode_hh(outratio,
                                  ["OO","RP","RS"],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        outratio = utils.explode_hh(outratio,
                                  [t+str(i) for t in ["O","P","S"]
                                   for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            ) 
    outratio = utils.aggregate_years(outratio, years_map, agg_method="mean")
    outratio["VALUE"] = outratio["VALUE"].round(DECIMALS)
    outratio.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"OutputActivityRatio.csv", index=False)   

#%%
# availability factor


    # simply get availability factor, if given, from raw input values
    # and save to csv
    # unit: -
    
    # load raw data set
    raw_data = pd.read_csv(snakemake.input.path_set_act_con)
    supply_af = pd.read_csv(snakemake.input.path_supply_af)
    avail_data = pd.concat([raw_data,supply_af]) 
    
    
    availfac = avail_data[["REGION","TECHNOLOGY","YEAR",
                         "VALUE"]][avail_data["VARIABLE"]=="AvailibilityFactor"]
    availfac["REGION"][availfac["REGION"].isnull()] = ":*"
    

    
    # FIXME: introduce some methodology for changing availability factor
    availfac["YEAR"][availfac["YEAR"].isnull()] = ":*"
    
    # if list of technologies given for rows, create separate rows
    availfac = availfac.assign(TECHNOLOGY=availfac['TECHNOLOGY'].str.split(',')).explode('TECHNOLOGY')
    
    # explode technologies if household disaggregated
    if snakemake.params.dic["scen_hh_disagg"]=="T":
        
        availfac = utils.explode_hh(availfac,
                                  ["OO","RP","RS"],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        availfac = utils.explode_hh(availfac,
                                  [t+str(i) for t in ["O","P","S"]
                                   for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )  
    availfac = utils.aggregate_years(availfac, years_map, agg_method="mean") 
    availfac["VALUE"] = availfac["VALUE"].round(DECIMALS)
    availfac.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"AvailabilityFactor.csv", index=False)

#%%
# capacity factor

    # process data to derive capacity factors and save to csv
    # unit: -
    
    # load sets
    sets = pd.read_csv(snakemake.input.path_set_sets)
    # explode technology set if household disaggregated
    if snakemake.params.dic["scen_hh_disagg"]=="T":
        
        sets = utils.explode_hh(sets,
                                  ["OO","RP","RS"],
                                  col="VALUE")
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        sets = utils.explode_hh(sets,
                                  [t+str(i) for t in ["O","P","S"]
                                   for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                  col="VALUE")
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )      
    # load annual heat demands per LSOA (currently adding up space and water),
    # demand fraction and length of typical periods, as well as peak
    # heat demand per LA before diversity (i.e., not the peak demand in the
    # model but calculated from the number of properties, temperature, etc.)
    
    # load sublocal aggregation
    slagg = pd.read_csv(snakemake.input.path_sublocal_agg,
                            index_col=["LSOA11CD"]).squeeze()
    
    dw_anndem_SH = (pd.read_csv(snakemake.input.path_dw_ann_dem_sh,
                             index_col=["LSOA11CD","PROPERTY_TYPE"])
                 + pd.read_csv(snakemake.input.path_dw_ann_dem_sh_nb,
                               index_col=["LSOA11CD","PROPERTY_TYPE"]))
    nd_anndem_SH = (pd.read_csv(snakemake.input.path_nd_ann_dem_sh,
                  index_col=["LSOA11CD","PROPERTY_TYPE"])
                 + pd.read_csv(snakemake.input.path_nd_ann_dem_sh_nb,
                               index_col=["LSOA11CD","PROPERTY_TYPE"])).reset_index()
    
    dw_anndem_HW = pd.read_csv(snakemake.input.path_dw_ann_dem_hw,
                               index_col=["LSOA11CD","PROPERTY_TYPE"])
    nd_anndem_HW =  pd.read_csv(snakemake.input.path_nd_ann_dem_hw)
    
    nd_anndem_SH["PROPERTY_TYPE"]="Non-domestic"
    nd_anndem_HW["PROPERTY_TYPE"]="Non-domestic"
    
    nd_anndem_SH = nd_anndem_SH.groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
    nd_anndem_HW = nd_anndem_HW.groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
    anndem_SH = pd.concat([dw_anndem_SH,nd_anndem_SH])
    anndem_HW = pd.concat([dw_anndem_HW,nd_anndem_HW])
    
    typperiods_HW = pd.read_csv(snakemake.input.path_tperiods_hw,
                            index_col=["LAD23CD","PROPERTY_TYPE"])
    typperiods_HW.columns.name = "TS"
    typperiods_HW = utils.groupby_LAD(typperiods_HW,disagg=slagg).mean().squeeze()
    typperiods_SH = pd.read_csv(snakemake.input.path_tperiods_sh,
                            index_col=["LAD23CD","PROPERTY_TYPE"])
    typperiods_SH.columns.name = "TS"
    typperiods_SH = utils.groupby_LAD(typperiods_SH,disagg=slagg).mean().squeeze()
    typperiods_length = pd.read_csv(snakemake.input.path_tperiods_len,
                            index_col=0)
    typperiods_length.index.name = "TS"
    
    peak_cap = pd.read_csv(snakemake.input.path_peakc_lsoa,
                            index_col=["LSOA11CD",
                                       "PROPERTY_TYPE",
                                       "TECHNOLOGY"])
    peak_cap = utils.groupby_LAD(peak_cap,disagg=slagg).sum()
    
    for i in range(2016,2023):
        peak_cap.loc[:,str(i)] = (peak_cap.loc[:,str(i-1)]
                                  +(peak_cap.loc[:,"2023"]-peak_cap.loc[:,"2015"])
                                  /(2023-2015)
                                  )
        
    if snakemake.params.dic["scen_hh_disagg"].startswith("T"):
        peak_cap.index = pd.MultiIndex.from_tuples([ (x[0],x[1],x[2]+"DDDE"+x[1].split("|")[1])
                                                  if x[1].startswith("Detached") else
                                                  (x[0],x[1],x[2]+"DDSD"+x[1].split("|")[1])
                                                  if x[1].startswith("Semi-detached") else
                                                  (x[0],x[1],x[2]+"DDTE"+x[1].split("|")[1])
                                                  if x[1].startswith("Terraced") else
                                                  (x[0],x[1],x[2]+"DDFL"+x[1].split("|")[1])
                                                  if x[1].startswith("Flats") else
                                                  (x[0],x[1],x[2]+"DNDO00")
                                                  if x[1]=="Non-domestic" else 0 for x
                                                  in peak_cap.index])
    else:              
        peak_cap.index = pd.MultiIndex.from_tuples([ (x[0],x[1],x[2]+"DDDE00")
                                                         if x[1]=="Detached" else
                                                         (x[0],x[1],x[2]+"DDSD00")
                                                         if x[1]=="Semi-detached" else
                                                         (x[0],x[1],x[2]+"DDTE00")
                                                         if x[1]=="Terraced" else
                                                         (x[0],x[1],x[2]+"DDFL00")
                                                         if x[1]=="Flats" else
                                                         (x[0],x[1],x[2]+"DNDO00")
                                                         if x[1]=="Non-domestic" else 0 for x
                                                         in peak_cap.index])

    peak_cap.index.names = [slagg.name,"PROPERTY_TYPE","TECHNOLOGY"]
    peak_cap_SH = peak_cap.loc[peak_cap.index.get_level_values("TECHNOLOGY").isin(
                    sets.loc[(sets["CATEGORY"]=="building_heat_distribution")|
                                 (sets["CATEGORY"]=="building_efficiency"),"VALUE"])]
    peak_cap_SHW = peak_cap.loc[peak_cap.index.get_level_values("TECHNOLOGY").isin(
                    sets.loc[(sets["CATEGORY"]=="building_heat_supply"),"VALUE"])]
    peak_cap_SH = peak_cap_SH.stack().to_frame()
    peak_cap_SH.index = peak_cap_SH.index.set_names("YEAR",level=3)
    peak_cap_SH.columns = ["VALUE"]
    peak_cap_SH = peak_cap_SH.dropna()
    peak_cap_SHW = peak_cap_SHW.stack().to_frame()
    peak_cap_SHW.index = peak_cap_SHW.index.set_names("YEAR",level=3)
    peak_cap_SHW.columns = ["VALUE"] 
    peak_cap_SHW = peak_cap_SHW.dropna()
    # calculate heat demand (power) for each timeslice and LAD within model 
    # (i.e., after diversity) (unit: GW)
    
    anndem_SH = utils.groupby_LAD(anndem_SH,disagg=slagg).sum()
    anndem_SH = anndem_SH.stack().to_frame()
    anndem_SH.index.names =[slagg.name,"PROPERTY_TYPE","YEAR"]
    anndem_HW = utils.groupby_LAD(anndem_HW,disagg=slagg).sum()
    anndem_HW = anndem_HW.stack().to_frame()
    anndem_HW.index.names =[slagg.name,"PROPERTY_TYPE","YEAR"]
    demand_capacity_SH = typperiods_SH / typperiods_length.T.values
    demand_capacity_SH = demand_capacity_SH.stack().to_frame()
    demand_capacity_SH.index.names = [slagg.name,"PROPERTY_TYPE","TIMESLICE"]
    demand_capacity_SH = (demand_capacity_SH.multiply(anndem_SH,axis=0)
                       /(8760*60*60))
    demand_capacity_HW = typperiods_HW / typperiods_length.T.values
    demand_capacity_HW = demand_capacity_HW.stack().to_frame()
    demand_capacity_HW.index.names = [slagg.name,"PROPERTY_TYPE","TIMESLICE"]
    demand_capacity_HW = (demand_capacity_HW.multiply(anndem_HW,axis=0)
                       /(8760*60*60))
    demand_capacity_SH.columns = ["VALUE"]
    demand_capacity_HW.columns = ["VALUE"]
    demand_capacity_SHW = demand_capacity_HW + demand_capacity_SH
    demand_capacity_SH = demand_capacity_SH.dropna()
    demand_capacity_SHW = demand_capacity_SHW.dropna()
    
    
    # calculate capacity factors per property type/technology, LAD, year by
    # dividing model peak demand by before diversity heat demand (this ensures
    # 1) not just the most efficient technologies run in non-peak hours, and 
    # 2) the actual before diversity peak capacity of heat techs needs to be
    # built)
    # capfac_SH = peak_cap_SH.stack().reset_index().merge(right=demand_capacity_SH,
    #                         on=[slagg.name,"PROPERTY_TYPE"],
    #                         )
    # capfac_SHW = peak_cap_SHW.stack().reset_index().merge(right=demand_capacity_SHW,
    #                         on=[slagg.name,"PROPERTY_TYPE"],
    #                         )
    
    capfac_SH = demand_capacity_SH/peak_cap_SH
    capfac_SHW = demand_capacity_SHW/peak_cap_SHW
    capfac = pd.concat([capfac_SH,capfac_SHW])
    capfac = capfac.dropna()
    capfac = capfac.droplevel("PROPERTY_TYPE")
    # remove inf where no demand in some years (and division by 0 occurs)
    capfac["VALUE"] = capfac["VALUE"].clip(upper=1)
    capfac = capfac.reset_index()   
    
    capfac = capfac.rename(columns={slagg.name:"REGION","level_3":"YEAR",
                                    "level_4":"TIMESLICE",
                                    0:"VALUE"})

    # sort columns
    capfac = capfac.reindex(["REGION","TECHNOLOGY","TIMESLICE","YEAR","VALUE"],
                            axis=1)
    
    # load local processed capacity factor data and concatenate
    local_capfac = pd.read_csv(snakemake.input.path_local_capf)
    supply_capfac = pd.read_csv(snakemake.input.path_supply_capfac)
    
    add_capfac = pd.concat([local_capfac,supply_capfac])
    add_capfac = add_capfac[["REGION","TECHNOLOGY","TIMESLICE","YEAR",
                     "VALUE"
                     ]][add_capfac["VARIABLE"]=="CapacityFactor"]
    
    capfac = pd.concat([capfac,add_capfac])  
    

    capfac = utils.aggregate_years(capfac, years_map, agg_method="mean")
    capfac["VALUE"] = np.true_divide(np.ceil(capfac["VALUE"] * 10**DECIMALS), 10**DECIMALS)
    
    capfac.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"CapacityFactor.csv", index=False)   

#%%
# capacity to activity ratio

    # simply assuming 31536 for all techs and save to csv
    # unit: TJ/GW
    
    # FIXME: implement this in a cleaner way, the below creates
    # data for techs that do not exist (e.g., national techs for local areas)
    # captoact = pd.DataFrame([[":*2",":*",31536]], columns=["REGION","TECHNOLOGY",
    #                                                    "VALUE"])
    
    # load sets
    sets = pd.read_csv(snakemake.input.path_set_sets)
    tech = sets["VALUE"][sets["SET"]=="TECHNOLOGY"]
    captoact = pd.DataFrame(data={"REGION":":*","TECHNOLOGY":tech,
                                  "VALUE":31536})

    # explode technologies if household disaggregated
    if snakemake.params.dic["scen_hh_disagg"]=="T":
        
        captoact = utils.explode_hh(captoact,
                                  ["OO","RP","RS"],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        captoact = utils.explode_hh(captoact,
                                  [t+str(i) for t in ["O","P","S"]
                                   for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )  
            
    #captoact = pd.DataFrame(columns=["REGION","TECHNOLOGY","VALUE"])
    captoact["VALUE"] = captoact["VALUE"].round(DECIMALS)
    captoact.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"CapacityToActivityUnit.csv", index=False) 
    
    
#%%
# capacity of one tech unit

    # create empty file 
    # unit: GW
    
    capofoneunit = pd.DataFrame( columns=["REGION","TECHNOLOGY","YEAR","VALUE"])
    capofoneunit["VALUE"] = capofoneunit["VALUE"].round(DECIMALS)
    capofoneunit.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"CapacityOfOneTechnologyUnit.csv", index=False) 
    

#%%
# residual capacities

    # get residual caps from processed input data and calculate for future
    # years based on lifetimes
    # unit: GW
    
    # load raw data set
    raw_data = pd.read_csv(snakemake.input.path_set_lt)
    
    # if list of technologies given for rows, create separate rows
    raw_data = raw_data.assign(TECHNOLOGY=raw_data['TECHNOLOGY'].str.split(',')).explode('TECHNOLOGY')

    
    # load building measure lifetime and concatenate
    supply_lt = pd.read_csv(snakemake.input.path_supply_lt)
    oplife = pd.concat([raw_data,supply_lt])     
    
    oplife = oplife[["REGION","TECHNOLOGY",
                         "VALUE"]][oplife["VARIABLE"]=="OperationalLife"]

  

    oplife["REGION"][oplife["REGION"].isnull()] = ":*"
    

    # explode technologies if household disaggregated
    if snakemake.params.dic["scen_hh_disagg"]=="T":
        
        oplife = utils.explode_hh(oplife,
                                  ["OO","RP","RS"],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        oplife = utils.explode_hh(oplife,
                                  [t+str(i) for t in ["O","P","S"]
                                   for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )

    
    
    oplife = oplife.set_index(["TECHNOLOGY"])
        
    # load model years
    
    sets = pd.read_csv(snakemake.input.path_set_sets)
    years = sets["VALUE"][sets["SET"]=="YEAR"].astype(int)
    
    # load sublocal aggregation
    slagg = pd.read_csv(snakemake.input.path_sublocal_agg,
                            index_col=["LSOA11CD"]).squeeze()
    
    
    
    # load residual capacities for first years
    heat_techs = utils.groupby_LAD(pd.read_csv(snakemake.input.path_res_caps_ht,
                                               index_col=["LSOA11CD","PROPERTY_TYPE",
                                             "TECHNOLOGY"]).droplevel(1),
                                                 disagg=slagg).sum()
    heat_techs.index = heat_techs.index.set_names("REGION",level=0)

    dist_techs = pd.read_csv(snakemake.input.path_local_res_caps,
                             index_col=[slagg.name,"TECHNOLOGY"])
    dist_techs.index = dist_techs.index.set_names("REGION",level=0)
    dhgen_techs = pd.read_csv(snakemake.input.path_dhgen_res_caps,
                             index_col=["LAD23CD","TECHNOLOGY"])
    dhgen_techs.index = dhgen_techs.index.set_names("REGION",level=0)
    
    
    supply_techs = pd.read_csv(snakemake.input.path_supply_rescap,
                               index_col=["REGION","TECHNOLOGY"])

    residual_caps = pd.concat([heat_techs,dist_techs,dhgen_techs,supply_techs])
    residual_caps = residual_caps.dropna(axis=0,how="all")
    residual_caps[["LI","R"]] = residual_caps.apply(lambda x: pd.Series([x.last_valid_index(),
                                                                 x[x.last_valid_index()]],
                                                                 index=["LI","R"]),
                                        axis=1)

    
    # merge lifetime on capacities and calculate future years based on
    # lifetime assuming linear decrease
    residual_caps = residual_caps.merge(right=oplife,
                                        left_index=True,
                                        right_index=True)
    years = [str(y) for y in years]
    residual_caps["R"] = residual_caps["R"]/residual_caps["VALUE"]
    residual_caps[[y for y in years if not y in residual_caps.columns]]=pd.NA
    residual_caps=residual_caps.apply(lambda x: pd.Series([x[i] 
                                               if not pd.isnull(x[i])
                                               else (x[x["LI"]]
                                                     -x["R"]
                                                     *(int(i)-int(x["LI"]))
                                                     )
                                               for i in years],
                                            index=years) ,axis=1)
    
    # set negative values to zero
    residual_caps[residual_caps < 0] = 0
    
    # rearrange dataframe
    residual_caps = residual_caps.stack()
    residual_caps = residual_caps.reset_index()
    residual_caps = residual_caps.rename(columns={"level_2":"YEAR",
                                                  0:"VALUE"})
    
    residual_caps = utils.aggregate_years(residual_caps,
                                          years_map, agg_method="mean")
    residual_caps["VALUE"] = np.true_divide(np.ceil(residual_caps["VALUE"]
                                                    * 10**DECIMALS),
                                            10**DECIMALS)
    
    # FIXME: to be seen if this will be needed or implemented differently
    # round up with one decimal less, to avoid issue with new capacity
    # constraint in historical years
    # residual_caps["VALUE"] = np.true_divide(np.ceil(residual_caps["VALUE"]
    #                                                 * 10**(DECIMALS-2)),
    #                                                 10**(DECIMALS-2))
    
    residual_caps.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"ResidualCapacity.csv", index=False)    
    
#%%
# capacity constraints

    # simply get capacity constraints, if given, from raw input values
    # and save to csv's
    # unit: GW
    

    # load raw data set
    raw_data = pd.read_csv(snakemake.input.path_set_cap_con)
    sets = pd.read_csv(snakemake.input.path_set_sets)
    years = sets["VALUE"][sets["SET"]=="YEAR"].astype(int)
    
    # process max capacity constraints
    # load local and building measure capacity constraints and concatenate
    local_con = pd.read_csv(snakemake.input.path_local_con)
    bm_con = pd.read_csv(snakemake.input.path_building_measures_con)
    ht_con =  pd.read_csv(snakemake.input.path_heat_tech_con)
    supply_con = pd.read_csv(snakemake.input.path_supply_capcon)
    con_data = pd.concat([raw_data,local_con,bm_con,ht_con,supply_con])  

    capcon = con_data[["REGION","TECHNOLOGY",
                       "YEAR","VALUE"
                       ]][con_data["VARIABLE"]=="TotalAnnualMaxCapacity"]
    capcon.loc[capcon["REGION"].isnull(),"REGION"] = ":*"
    

    capcon = utils.interpolate_timeseries(capcon,model_periods=years)
    # fill nan's created through interpolation (due to no measures for a few
    # categories in E09000001)
    capcon.loc[capcon["VALUE"].isna(),"VALUE"] = 0
    
    
    # explode technologies if household disaggregated
    if snakemake.params.dic["scen_hh_disagg"]=="T":
        
        capcon = utils.explode_hh(capcon,
                                  ["OO","RP","RS"],
                                  col="TECHNOLOGY")
        sets = utils.explode_hh(sets,
                                  ["OO","RP","RS"],
                                  col="VALUE")
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        capcon = utils.explode_hh(capcon,
                                  [t+str(i) for t in ["O","P","S"]
                                   for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                  col="TECHNOLOGY")
        sets = utils.explode_hh(sets,
                                  [t+str(i) for t in ["O","P","S"]
                                   for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                  col="VALUE")
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )  
    capcon = utils.aggregate_years(capcon, years_map, agg_method="mean")
    capcon.loc[:,"VALUE"] = np.true_divide(np.ceil(capcon["VALUE"] * 10**DECIMALS), 10**DECIMALS)
    
    capcon.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"TotalAnnualMaxCapacity.csv", index=False) 

    
    # process max capacity investment constraints
    # load building measure capacity constraints and concatenate
    bm_icon = pd.read_csv(snakemake.input.path_building_measures_dc)
    ht_icon = pd.read_csv(snakemake.input.path_heat_tech_invcon)
    
    lg_icon = pd.read_csv(snakemake.input.path_loc_gov_invcon)
    lg_icon = lg_icon.set_index([col for col 
                                 in lg_icon.columns
                                 if col!="VALUE"])
    slagg = pd.read_csv(snakemake.input.path_sublocal_agg,
                            index_col=["LSOA11CD"]).squeeze()
    if not lg_icon.empty:
        lg_icon = utils.groupby_LAD(lg_icon,disagg=slagg).mean()
    lg_icon.index = lg_icon.index.set_names("REGION",level=0)
    lg_icon = lg_icon.reset_index()
    
    icon_data = pd.concat([raw_data,ht_icon,bm_icon,lg_icon])  
    
    capicon = icon_data[["REGION","TECHNOLOGY",
                       "YEAR","VALUE"
                       ]][icon_data["VARIABLE"]=="TotalAnnualMaxCapacityInvestment"]
    capicon.loc[capicon["REGION"].isnull(),"REGION"] = ":*"
    
    capicon.loc[:,"VALUE"] = capicon["VALUE"].astype(float)


    capicon = utils.interpolate_timeseries(capicon,model_periods=years)
    capicon.loc[capicon["VALUE"].isna(),"VALUE"] = -1
    
    # FIXME: this should be added above, but currently should not be interpolated
    supicon = supply_con[["REGION","TECHNOLOGY",
                       "YEAR","VALUE"
                       ]][supply_con["VARIABLE"]=="TotalAnnualMaxCapacityInvestment"]
    supicon["REGION"][supicon["REGION"].isnull()] = ":*"
    
    capicon = pd.concat([capicon,supicon])
    
    # add technology investment ban constraints if triggered
    # FIXME: properly check for right parameter syntax and raise error otherwise
    if snakemake.params.dic["scen_tech_bans"]=="-":                             
        pass
    else:
        
        bans = snakemake.params.dic["scen_tech_bans"]
        bans = bans.replace("jt-01",("HPDDFLP-:*-2023|HPDDTEP-:*-2023|"
                                     "HPDDDEP-:*-2023|HPDDSDP-:*-2023|"
                                     "BELODDFLP-:*-2023|BELODDTEP-:*-2023|"
                                     "BELODDDEP-:*-2023|BELODDSDP-:*-2023|"
                                     "BEMEDDFLP-:*-2023|BEMEDDTEP-:*-2023|"
                                     "BEMEDDDEP-:*-2023|BEMEDDSDP-:*-2023|"
                                     "BEHIDDFLP-:*-2023|BEHIDDTEP-:*-2023|"
                                     "BEHIDDDEP-:*-2023|BEHIDDSDP-:*-2023|"
                                     "HPDDFLS-:*-2023|HPDDTES-:*-2023|"
                                     "HPDDDES-:*-2023|HPDDSDS-:*-2023|"
                                     "BELODDFLS-:*-2023|BELODDTES-:*-2023|"
                                     "BELODDDES-:*-2023|BELODDSDS-:*-2023|"
                                     "BEMEDDFLS-:*-2023|BEMEDDTES-:*-2023|"
                                     "BEMEDDDES-:*-2023|BEMEDDSDS-:*-2023|"
                                     "BEHIDDFLS-:*-2023|BEHIDDTES-:*-2023|"
                                     "BEHIDDDES-:*-2023|BEHIDDSDS-:*-2023|"
                                     "HPDDFLO0-:*-2023|HPDDTEO0-:*-2023|"
                                     "HPDDDEO0-:*-2023|HPDDSDO0-:*-2023|"
                                     "BELODDFLO0-:*-2023|BELODDTEO0-:*-2023|"
                                     "BELODDDEO0-:*-2023|BELODDSDO0-:*-2023|"
                                     "BEMEDDFLO0-:*-2023|BEMEDDTEO0-:*-2023|"
                                     "BEMEDDDEO0-:*-2023|BEMEDDSDO0-:*-2023|"
                                     "BEHIDDFLO0-:*-2023|BEHIDDTEO0-:*-2023|"
                                     "BEHIDDDEO0-:*-2023|BEHIDDSDO0-:*-2023|"
                                     "HPDDFLO1-:*-2023|HPDDTEO1-:*-2023|"
                                     "HPDDDEO1-:*-2023|HPDDSDO1-:*-2023|"
                                     "BELODDFLO1-:*-2023|BELODDTEO1-:*-2023|"
                                     "BELODDDEO1-:*-2023|BELODDSDO1-:*-2023|"
                                     "BEMEDDFLO1-:*-2023|BEMEDDTEO1-:*-2023|"
                                     "BEMEDDDEO1-:*-2023|BEMEDDSDO1-:*-2023|"
                                     "BEHIDDFLO1-:*-2023|BEHIDDTEO1-:*-2023|"
                                     "BEHIDDDEO1-:*-2023|BEHIDDSDO1-:*-2023"))
        bans = bans.replace("jt-03",("HPDDFLP-:*-2023|HPDDTEP-:*-2023|"
                                     "HPDDDEP-:*-2023|HPDDSDP-:*-2023|"
                                     "BELODDFLP-:*-2023|BELODDTEP-:*-2023|"
                                     "BELODDDEP-:*-2023|BELODDSDP-:*-2023|"
                                     "BEMEDDFLP-:*-2023|BEMEDDTEP-:*-2023|"
                                     "BEMEDDDEP-:*-2023|BEMEDDSDP-:*-2023|"
                                     "BEHIDDFLP-:*-2023|BEHIDDTEP-:*-2023|"
                                     "BEHIDDDEP-:*-2023|BEHIDDSDP-:*-2023"))#"|HPDDTEO0-:*-2023|HPDDDEO0-:*-2023|HPDDSDO0-:*-2023|BELODDFLO0-:*-2023|BELODDTEO0-:*-2023|BELODDDEO0-:*-2023|BELODDSDO0-:*-2023")
        bans = bans.split("|")

        for b in bans:
            logger.info("Working on ban: "+b)
            tech,reg,year = b.split("-")
            ic = "S12000014;S12000040;S12000036;E06000002;E06000004;E06000003;E08000011;E08000013;E08000014;E08000015;E08000012;E06000013;E06000012;E06000010;E06000011;E06000045;E07000086;E07000087;E07000088;E06000044;W06000009;W06000011;W06000012;W06000014;W06000015;W06000022;W06000019"
            nic = ",".join([l for l in utils.get_entity_lookup(["LAD"])["LAD23CD"].to_list()
                            if l not in ic.split(";")])
            if reg == "NIC":
                reg = nic
            
            techs = [t for t in sets.loc[sets["SET"]=="TECHNOLOGY","VALUE"]
                     if tech in t]

            byears = [y for y in sets.loc[sets["SET"]=="YEAR","VALUE"].astype(int)
                     if y >= int(year)]
 
            df = pd.DataFrame([[reg,t,byears,0] for t in techs],
                              columns=["REGION",
                                       "TECHNOLOGY",
                                       "YEAR",
                                       "VALUE"]).explode("YEAR")
            df["YEAR"] = df["YEAR"].astype(int)
            df = df.assign(REGION=df['REGION'].str.split(',')).explode('REGION')
            
            # set region of to be overwritten entries to :* so these will be
            # replaced when ban is set for :*
            # FIXME: this would still create issues if other constraints are
            # set for :* and the ban for a specific region (it would not
            # overwrite) - currently not the case
            capicon.loc[pd.Index(capicon[["TECHNOLOGY","YEAR"]]).isin(
                                pd.Index(df.loc[df["REGION"]==":*",
                                                ["TECHNOLOGY","YEAR"]])),
                        "REGION"] = ":*"
            capicon = pd.concat([capicon,df])
            # remove duplicates, i.e., overwrite previous constraint
            capicon = capicon.drop_duplicates(["REGION","TECHNOLOGY","YEAR"],
                                              keep="last")
    logger.info("Processed all bans.")
    # print(capicon[capicon["TECHNOLOGY"].str.startswith("NGBO")])  
    # print(capicon[capicon["TECHNOLOGY"].str.startswith("ASHP")&(
    #     capicon["YEAR"]==2035)])
    # explode technologies if household disaggregated
    if snakemake.params.dic["scen_hh_disagg"]=="T":
        
        capicon = utils.explode_hh(capicon,
                                  ["OO","RP","RS"],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        capicon = utils.explode_hh(capicon,
                                  [t+str(i) for t in ["O","P","S"]
                                   for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )  
    capicon = utils.aggregate_years(capicon, years_map, agg_method="sum",
                                    max_constraint_default=-1)
    capicon["VALUE"] = capicon["VALUE"].round(DECIMALS)
    capicon.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"TotalAnnualMaxCapacityInvestment.csv", index=False) 
    
   # process capacity fraction constraint (for electricity sector)

   # load fraction and tag data
    capfrac_data = pd.read_csv(snakemake.input.path_supply_capfrac)
    
 
    capfrac = capfrac_data[["REGION","TECHNOLOGY",
                            "YEAR","VALUE"
                            ]][capfrac_data["VARIABLE"]=="CapacityFraction"]
    capfrac["YEAR"] = capfrac["YEAR"].astype(int)

    techtag = capfrac_data[["REGION","TECHNOLOGY", "VALUE"
                            ]][capfrac_data["VARIABLE"]=="CapacityFractionTagTechnology"]
 
    capfrac = utils.interpolate_timeseries(capfrac,model_periods=years,
                                           direction="forward",dropna=True)
   
    
    capfrac= utils.aggregate_years(capfrac, years_map, agg_method="mean")
    capfrac["VALUE"] = capfrac["VALUE"].round(DECIMALS)
    
    capfrac.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"CapacityFraction.csv", index=False) 
    
    techtag.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"CapacityFractionTagTechnology.csv", index=False) 
    
    
    ### minimum capacity constraint
    sets = pd.read_csv(snakemake.input.path_set_sets)
    if snakemake.params.dic["scen_mcap_con"]=="-":                             
        lims = pd.DataFrame( columns=["REGION","TECHNOLOGY","YEAR","VALUE"])
    else:
        lims = pd.DataFrame( columns=["REGION","TECHNOLOGY","YEAR","VALUE"])
        cons = snakemake.params.dic["scen_mcap_con"]
        cons = cons.replace("jt-02",("BEMEDDFLS-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEMEDDTES-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEMEDDDES-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEMEDDSDS-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEHIDDFLS-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEHIDDTES-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEHIDDDES-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEHIDDSDS-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEMEDDFLP0-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEMEDDTEP0-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEMEDDDEP0-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEMEDDSDP0-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEHIDDFLP0-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEHIDDTEP0-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEHIDDDEP0-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEHIDDSDP0-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEMEDDFLP1-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEMEDDTEP1-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEMEDDDEP1-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEMEDDSDP1-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEHIDDFLP1-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEHIDDTEP1-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEHIDDDEP1-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEHIDDSDP1-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEMEDDFLO0-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEMEDDTEO0-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEMEDDDEO0-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEMEDDSDO0-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEHIDDFLO0-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEHIDDTEO0-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEHIDDDEO0-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEHIDDSDO0-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEMEDDFLO1-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEMEDDTEO1-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEMEDDDEO1-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEMEDDSDO1-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEHIDDFLO1-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEHIDDTEO1-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEHIDDDEO1-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060|"
                                     "BEHIDDSDO1-t-:*-0y2025a0.2y2030a0.75y2050a0.75y2060"))
        cons = cons.replace("jt-028","BEMEDDFLS-r-:*-0y2025a0.9y2030a0.8y2060|BEMEDDTES-r-:*-0y2025a0.8y2030a0.8y2060|BEMEDDDES-r-:*-0y2025a0.8y2030a0.8y2060|BEMEDDSDS-r-:*-0y2025a0.8y2030a0.8y2060|BEHIDDFLS-r-:*-0y2025a0.8y2030a0.8y2060|BEHIDDTES-r-:*-0y2025a0.8y2030a0.8y2060|BEHIDDDES-r-:*-0y2025a0.8y2030a0.8y2060|BEHIDDSDS-r-:*-0y2025a0.8y2030a0.8y2060|BEMEDDFLP0-r-:*-0y2025a0.8y2030a0.8y2060|BEMEDDTEP0-r-:*-0y2025a0.8y2030a0.8y2060|BEMEDDDEP0-r-:*-0y2025a0.8y2030a0.8y2060|BEMEDDSDP0-r-:*-0y2025a0.8y2030a0.8y2060|BEHIDDFLP0-r-:*-0y2025a0.8y2030a0.8y2060|BEHIDDTEP0-r-:*-0y2025a0.8y2030a0.8y2060|BEHIDDDEP0-r-:*-0y2025a0.8y2030a0.8y2060|BEHIDDSDP0-r-:*-0y2025a0.8y2030a0.8y2060")
        cons = cons.replace("jt-026","BEMEDDFLS-r-:*-0y2025a0.6y2030a0.6y2060|BEMEDDTES-r-:*-0y2025a0.6y2030a0.6y2060|BEMEDDDES-r-:*-0y2025a0.6y2030a0.6y2060|BEMEDDSDS-r-:*-0y2025a0.6y2030a0.6y2060|BEHIDDFLS-r-:*-0y2025a0.6y2030a0.6y2060|BEHIDDTES-r-:*-0y2025a0.6y2030a0.6y2060|BEHIDDDES-r-:*-0y2025a0.6y2030a0.6y2060|BEHIDDSDS-r-:*-0y2025a0.6y2030a0.6y2060|BEMEDDFLP0-r-:*-0y2025a0.6y2030a0.6y2060|BEMEDDTEP0-r-:*-0y2025a0.6y2030a0.6y2060|BEMEDDDEP0-r-:*-0y2025a0.6y2030a0.6y2060|BEMEDDSDP0-r-:*-0y2025a0.6y2030a0.6y2060|BEHIDDFLP0-r-:*-0y2025a0.6y2030a0.6y2060|BEHIDDTEP0-r-:*-0y2025a0.6y2030a0.6y2060|BEHIDDDEP0-r-:*-0y2025a0.6y2030a0.6y2060|BEHIDDSDP0-r-:*-0y2025a0.6y2030a0.6y2060")
        cons = cons.split("|")
        
        techl = sets.loc[(sets["SET"]=="TECHNOLOGY"),
                         "VALUE"]
        techl.name = "TECHNOLOGY"
        years = sets["VALUE"][sets["SET"]=="YEAR"].astype(str)
      
        
        cons = [c.split("-")[0:2]+[r]+[c.split("-")[3]] for c in cons for r in c.split("-")[2].split(";")]
        cons = ["-".join(c) for c in cons]
        
        if snakemake.params.dic["scen_hh_disagg"]=="T":
            techl = utils.explode_hh(techl.to_frame(),
                                  ["OO","RP","RS"],
                                  col="TECHNOLOGY")
        elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
            techl = utils.explode_hh(techl.to_frame(),
                                      [t+str(i) for t in ["O","P","S"]
                                       for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                      col="TECHNOLOGY")    
        elif snakemake.params.dic["scen_hh_disagg"]=="-":
            techl = techl.to_frame()                            
        else:
            raise NotImplementedError(
                    f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                    " not a valid option of the 'scen_hh_disagg' parameter."
                    " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
                )

        
        for ic, c in enumerate(cons):

            tech,ty,rg,vys = c.split("-")
            
            techs = pd.DataFrame([t for t in techl["TECHNOLOGY"]
                                  if tech in t],columns=["TECHNOLOGY"])
            logger.info(techs)
            # # FIXME: currently always applies to all regions
            # # add region, explode, and disaggregate to sublocal agg if needed
            # techs["LAD23CD"] = rg
            # techs = techs.assign(LAD23CD=techs["LAD23CD"].str.split(
            #                             ';')).explode('LAD23CD')
            # techs = techs.set_index(["LAD23CD","TECHNOLOGY"])
            # # if rg != ':*':
            # #     techs["VALUE"] = 0
            # #     techs = utils.groupby_LAD(techs,disagg=slagg).mean()
            # #     techs = techs.drop("VALUE",axis=1)
            # techs = techs.reset_index()
            # # techs = techs.rename(columns={slagg.name:"REGION",
            # #                               "LAD23CD":"REGION"})
            
            lim = years.copy().to_frame()
            lim = lim.rename(columns={"VALUE":"YEAR"})
            lim["VALUE"] = np.nan
            vys = vys.split("a")
            for vy in vys:
                value, year = vy.split("y")
                lim.loc[lim["YEAR"]==year,"VALUE"]=value
            lim["VALUE"] = lim["VALUE"].astype(float)
            # interpolate
            lim["VALUE"] = lim["VALUE"].interpolate(limit_area="inside")
            lim["YEAR"] = lim["YEAR"].astype("int")
            lim = lim.set_index(["YEAR"]).dropna()
            
            lim.loc[:,"TECHNOLOGY"] = ",".join(techs["TECHNOLOGY"].to_list())
            lim = lim.assign(TECHNOLOGY=lim["TECHNOLOGY"].str.split(
                                        ',')).explode("TECHNOLOGY")
            lim.loc[:,"REGION"]= rg
            lim = lim.assign(REGION=lim["REGION"].str.split(
                                        ';')).explode('REGION')

            lim = lim.set_index(["REGION","TECHNOLOGY"],append=True)
            lim = lim.reorder_levels(["REGION","TECHNOLOGY","YEAR"])
            lim = lim.reset_index()
            lim = utils.aggregate_years(lim, years_map, agg_method="mean")
            lim["YEAR"] = lim["YEAR"].astype("str")
            lim = lim.set_index(["REGION","TECHNOLOGY","YEAR"])
            
            if ty=="r":

                # FIXME: this would fail if capicon and capcon are one set 
                # for :* and the other for a number of specific areas
                capicon = capicon.loc[capicon["VALUE"]>=0,:]
                capicon.loc[:,"YEAR"] = capicon.loc[:,"YEAR"].astype(str)
                capicon.loc[:,"YEAR"] = capicon.loc[:,"YEAR"].replace(":*",",".join(list(years_map.astype(str).unique())))
                capicon = capicon.assign(YEAR=capicon["YEAR"].str.split(
                                            ',')).explode("YEAR")
                
                capicon.loc[:,"VALUE"] = capicon.groupby(["REGION","TECHNOLOGY"])["VALUE"].cumsum()
                logger.info(capicon)
                logger.info(capcon)
                capcon = capcon.reset_index()
                capcon.loc[:,"YEAR"] = capcon.loc[:,"YEAR"].astype(str)
                capcon = capcon.loc[capcon["VALUE"]>=0,:]
                capcon.loc[:,"YEAR"] = capcon.loc[:,"YEAR"].replace(":*",",".join(list(years_map.astype(str).unique())))
                capcon = capcon.assign(YEAR=capcon["YEAR"].str.split(
                                            ',')).explode("YEAR")
                
                capcon = capcon.merge(capicon,on=["REGION","TECHNOLOGY","YEAR"],
                                      how="outer",suffixes=("_t","_i"))
                capcon.loc[:,"VALUE"] = 0
                logger.info(capcon)
                capcon.loc[capcon["VALUE_t"].isna(),
                           "VALUE"] = capcon.loc[capcon["VALUE_t"].isna(),
                                                 "VALUE_i"]
                capcon.loc[capcon["VALUE_i"].isna(),
                           "VALUE"] = capcon.loc[capcon["VALUE_i"].isna(),
                                                 "VALUE_t"]
                capcon.loc[~capcon["VALUE_i"].isna()&
                           ~capcon["VALUE_t"].isna(),
                           "VALUE"] = capcon.loc[~capcon["VALUE_i"].isna()&
                                      ~capcon["VALUE_t"].isna(),
                                      ["VALUE_t","VALUE_i"]].min(axis=1)
                capcon = capcon.loc[:,["REGION","TECHNOLOGY","YEAR","VALUE"]]
                capcon = capcon.set_index(["REGION","TECHNOLOGY","YEAR"])
                logger.info(capcon)
                logger.info(lim)
                if rg == ":*":
                    lim = lim.droplevel("REGION")
                lim = lim.mul(capcon.clip(0)).dropna()
                logger.info(lim)
            elif ty=="t":

                #capcon = capcon.reset_index()
                capcon.loc[:,"YEAR"] = capcon.loc[:,"YEAR"].astype(str)
                capcon = capcon.loc[capcon["VALUE"]>=0,:]
                capcon.loc[:,"YEAR"] = capcon.loc[:,"YEAR"].replace(":*",",".join(list(years_map.astype(str).unique())))
                capcon = capcon.assign(YEAR=capcon["YEAR"].str.split(
                                            ',')).explode("YEAR")
                
                capcon = capcon.set_index(["REGION","TECHNOLOGY","YEAR"])
                logger.info(capcon)
                logger.info(lim)
                if rg == ":*":
                    lim = lim.droplevel("REGION")
                lim = lim.mul(capcon.clip(0)).dropna()
                capcon = capcon.reset_index()
                logger.info(lim)
            elif ty=="a":
                pass
            
            lim = lim.reset_index()
            lim = lim.loc[:,["REGION","TECHNOLOGY","YEAR","VALUE"]]    
            lims = pd.concat([lims,lim])
        # lims["VALUE"] = lims["VALUE"].round(DECIMALS)    
        lims.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                    +"/data/"+"TotalAnnualMinCapacity.csv", index=False)
    
         
    # FIXME: Add processing for other capacity constraints if required

    
#%% 
# storage parameters

    # load storage tag data and save to file
    storage_data = pd.read_csv(snakemake.input.path_supply_storage)
  

    stotag = storage_data[["REGION","FUEL","VALUE"
                          ]][storage_data["VARIABLE"]=="StorageTagFuel"]


    stotag.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                +"/data/"+"StorageTagFuel.csv", index=False) 
 


#%%
# activity constraints

    # simply get activity constraints, if given, from raw input values
    # and save to csv
    # unit: GJ
    

    
    # load raw data set
    raw_data = pd.read_csv(snakemake.input.path_set_act_con)
    sets = pd.read_csv(snakemake.input.path_set_sets)
    years = sets["VALUE"][sets["SET"]=="YEAR"].astype(int)

    # load other activity constraints and concatenate
    pot_con = pd.read_csv(snakemake.input.path_potential_con)
    pot_heat = pd.read_csv(snakemake.input.path_heat_pots)
    sup_con = pd.read_csv(snakemake.input.path_supply_actcon)
    con_data = pd.concat([raw_data,pot_con,pot_heat,sup_con])  
    
    actucon = con_data[["REGION","TECHNOLOGY",
                       "YEAR","VALUE",
                       ]][con_data["VARIABLE"]=="TotalTechnologyAnnualActivityUpperLimit"]
    actucon["REGION"][actucon["REGION"].isnull()] = ":*"
    actucon = utils.interpolate_timeseries(actucon,model_periods=years,
                                           dropna=True)
    actlcon = con_data[["REGION","TECHNOLOGY",
                       "YEAR","VALUE"
                       ]][con_data["VARIABLE"]=="TotalTechnologyAnnualActivityLowerLimit"]
    actlcon["REGION"][actlcon["REGION"].isnull()] = ":*"
    actlcon = utils.interpolate_timeseries(actlcon,model_periods=years,
                                           dropna=True)
    
    # explode technologies if household disaggregated
    if snakemake.params.dic["scen_hh_disagg"]=="T":
        
        actucon = utils.explode_hh(actucon,
                                  ["OO","RP","RS"],
                                  col="TECHNOLOGY")
        actlcon = utils.explode_hh(actlcon,
                                  ["OO","RP","RS"],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        actlcon = utils.explode_hh(actlcon,
                                  [t+str(i) for t in ["O","P","S"]
                                   for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                  col="TECHNOLOGY")
        actucon = utils.explode_hh(actucon,
                                  [t+str(i) for t in ["O","P","S"]
                                   for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                  col="TECHNOLOGY")        
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )  
    actucon = utils.aggregate_years(actucon, years_map, agg_method="sum")
    actucon["VALUE"] = actucon["VALUE"].round(DECIMALS)
    actucon.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"TotalTechnologyAnnualActivityUpperLimit.csv", index=False) 
    
    actlcon = utils.aggregate_years(actlcon, years_map, agg_method="sum")
    actlcon["VALUE"] = actlcon["VALUE"].round(DECIMALS)
    actlcon.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"TotalTechnologyAnnualActivityLowerLimit.csv", index=False)

#%%
# UDC constraints
    # simply get UDC constraints, if given, from raw input values
    # and save to csv
    # unit: GW/GJ  
    
    # load raw data set
    raw_data = pd.read_csv(snakemake.input.path_set_act_con)
    sets = pd.read_csv(snakemake.input.path_set_sets)
    years = sets["VALUE"][sets["SET"]=="YEAR"].astype(int)
    
    # load UDC heat tech constraints
    tag = pd.read_csv(snakemake.input.path_heat_tech_udc_tagtech)
    UDCtag = pd.read_csv(snakemake.input.path_heat_tech_udc_tag)
    cons = pd.read_csv(snakemake.input.path_heat_tech_udc_invcon)
    
    # any necessary interpolation should be/is done beforehand,
    # given this requires some knowledge about the type of constraint
    
    # explode technologies if household disaggregated
    if snakemake.params.dic["scen_hh_disagg"]=="T":
        
        tag = utils.explode_hh(tag,
                               ["OO","RP","RS"],
                               col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        tag = utils.explode_hh(tag,
                               [t+str(i) for t in ["O","P","S"]
                                for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )  
    cons = utils.aggregate_years(cons, years_map, agg_method="sum",
                                    max_constraint_default=-1)
    tag = utils.aggregate_years(tag, years_map, agg_method="mean",
                                    max_constraint_default=-1)
    cons["VALUE"] = cons["VALUE"].round(DECIMALS)
    
    invtag = tag[["REGION","TECHNOLOGY","UDC",
                       "YEAR","VALUE"
                       ]][tag["TYPE"]=="NewCapacity"]
    tcaptag = tag[["REGION","TECHNOLOGY","UDC",
                       "YEAR","VALUE"
                       ]][tag["TYPE"]=="TotalCapacity"]
    acttag = tag[["REGION","TECHNOLOGY","MODE_OF_OPERATION","UDC",
                       "YEAR","VALUE"
                       ]][tag["TYPE"]=="Activity"]
    
    UDCtag.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"UDCTag.csv", index=False) 
    cons.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"UDCConstant.csv", index=False) 
    invtag.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"UDCMultiplierNewCapacity.csv", index=False)
    tcaptag.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"UDCMultiplierTotalCapacity.csv", index=False) 
    acttag.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"UDCMultiplierActivity.csv", index=False) 
    
    
#%%
# emission activity ratio

    # simply get emission activity ratios, if given, from raw input values
    # and save to csv
    # unit: kt/TJ    

       
    # load raw data set
    raw_data = pd.read_csv(snakemake.input.path_set_em)
    sets = pd.read_csv(snakemake.input.path_set_sets)
    years = sets["VALUE"][sets["SET"]=="YEAR"].astype(int)
    
    supply_emf = pd.read_csv(snakemake.input.path_supply_emf)
    emratio = pd.concat([raw_data,supply_emf])     
    
    emratio = emratio.loc[emratio["VARIABLE"]=="EmissionActivityRatio",["REGION","TECHNOLOGY","EMISSION","MODE_OF_OPERATION",
                          "YEAR",
                         "VALUE"]]
    emratio["REGION"][emratio["REGION"].isnull()] = ":*"
    
    # if no year given, assume data is the same for all years, interpolate, if 
    # necessary, remaining values
    emratio["YEAR"]= emratio["YEAR"].fillna(value=":*")
    emratio = utils.interpolate_timeseries(emratio,model_periods=years)
    
    # if list of technologies given for rows, create separate rows
    emratio = emratio.assign(TECHNOLOGY=emratio['TECHNOLOGY'].str.split(',')).explode('TECHNOLOGY')
    
    # explode technologies if household disaggregated
    if snakemake.params.dic["scen_hh_disagg"]=="T":
        
        emratio = utils.explode_hh(emratio,
                                  ["OO","RP","RS"],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        emratio = utils.explode_hh(emratio,
                                  [t+str(i) for t in ["O","P","S"]
                                   for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                  col="TECHNOLOGY")
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )  
    emratio = utils.aggregate_years(emratio, years_map, agg_method="mean")
    emratio["VALUE"] = emratio["VALUE"].round(DECIMALS)
    emratio.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"EmissionActivityRatio.csv", index=False) 

    
#%%
# emission limits


    # get emission limits directly from raw data set as well as local, if any
    # and save to csv
    # unit: kt/a
    

    # load raw data set
    raw_data = pd.read_csv(snakemake.input.path_set_em)
    sets = pd.read_csv(snakemake.input.path_set_sets)
    years = sets["VALUE"][sets["SET"]=="YEAR"].astype(int)
    
    if snakemake.params.dic["scen_em_limit"]=="-":                             
        raw_data = raw_data.drop("SCENARIO",axis=1).iloc[0:0]
    else:
        raw_data = raw_data.loc[raw_data["SCENARIO"]==snakemake.params.dic["scen_em_limit"]]
        raw_data = raw_data.drop("SCENARIO",axis=1)
        if raw_data.empty:
            raise NotImplementedError(
                    f"'{snakemake.params.dic['scen_em_limit']}' is currently"
                    " not a valid option of the 'scen_em_limit' parameter."
                    " Check the data file for valid options."
                )  
        
    local_data = pd.read_csv(snakemake.input.path_local_gov_em)
  
    
    emdata = pd.concat([raw_data,local_data])
    anemlimit = emdata[["REGION","EMISSION",
                          "YEAR",
                         "VALUE"]][emdata["VARIABLE"]=="AnnualEmissionLimit"]
    anemlimit.loc[anemlimit["REGION"].isnull(),"REGION"] = ":*"
    # local data overwrites other if given for same index
    anemlimit = anemlimit.drop_duplicates(["REGION","EMISSION","YEAR"],
                                          keep="last")
    anemlimit = utils.interpolate_timeseries(anemlimit,model_periods=years,
                                             direction="forward")
    anemlimit = anemlimit.fillna(-1)
    
    # adjust national emission limit if local emission targets set and adjustment triggered
    if ((snakemake.params.dic["scen_local_gov"]=="CPS")
        or (snakemake.params.dic["scen_local_gov"]=="CA")
        or (snakemake.params.dic["scen_local_gov"]=="A")):
        # load annual demand data    
        dw_anndem_SH = (pd.read_csv(snakemake.input.path_dw_ann_dem_sh,
                                 index_col=["LSOA11CD","PROPERTY_TYPE"])
                     + pd.read_csv(snakemake.input.path_dw_ann_dem_sh_nb,
                                   index_col=["LSOA11CD","PROPERTY_TYPE"]))
        nd_anndem_SH = (pd.read_csv(snakemake.input.path_nd_ann_dem_sh,
                      index_col=["LSOA11CD","PROPERTY_TYPE"])
                     + pd.read_csv(snakemake.input.path_nd_ann_dem_sh_nb,
                                   index_col=["LSOA11CD","PROPERTY_TYPE"])).reset_index()
        
        anndem_sh = pd.concat([dw_anndem_SH,nd_anndem_SH]).groupby("LSOA11CD").sum()
        anndem_sh = utils.groupby_LAD(anndem_sh).sum()
        
        et = local_data[["REGION","YEAR","VALUE"]].set_index(["REGION","YEAR"]).unstack().droplevel(0,axis=1)
        et[et==0]=1
        et = et.fillna(0)
        et.index.name="LAD23CD"
        et.columns = et.columns.astype(str)
        
        red = (anndem_sh*et).fillna(0).sum()/anndem_sh.sum()
        red = 1-red
        red.index.name="YEAR"
        red.index=red.index.astype(float)

        al = anemlimit[anemlimit["REGION"]=="UK"].set_index(["REGION","EMISSION","YEAR"])
        al = al["VALUE"].multiply(red,axis=0)
        #al = al.reset_index()
        al = al.to_frame()
        al = al.rename(columns={0:"VALUE"})
        anemlimit = anemlimit.set_index(["REGION","EMISSION","YEAR"])
        anemlimit.loc["UK"]=al
        anemlimit=anemlimit.reset_index()
        
    # FIXME: introduce some methodology for changing emission limit - if
    # necessary
    anemlimit.loc[anemlimit["YEAR"].isnull(),"YEAR"] = ":*"
    
    anemlimit = utils.aggregate_years(anemlimit, years_map, agg_method="sum",
                                      max_constraint_default=-1)
    anemlimit["VALUE"] = anemlimit["VALUE"].round(DECIMALS)
    anemlimit.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"AnnualEmissionLimit.csv", index=False) 
    
    mpemlimit = raw_data[["REGION","EMISSION",
                         "VALUE"]][raw_data["VARIABLE"]=="ModelPeriodEmissionLimit"]
    mpemlimit["REGION"][mpemlimit["REGION"].isnull()] = ":*"
    mpemlimit["VALUE"] = mpemlimit["VALUE"].round(DECIMALS)
    mpemlimit.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"ModelPeriodEmissionLimit.csv", index=False) 
    

#%%
# reserve margin

    # simply get reserve margin and tags, if given, from raw input values
    # and save to csv
    # unit: -
    
    
    # load raw data set
    raw_data = pd.read_csv(snakemake.input.path_set_reservem)    
    
    rmargin = raw_data[["REGION",
                          "YEAR",
                         "VALUE"]][raw_data["VARIABLE"]=="ReserveMargin"]
    rmargin["REGION"][rmargin["REGION"].isnull()] = ":*"
    
    # FIXME: introduce some methodology for changing reserve marg. (necessary!?)
    rmargin["YEAR"][rmargin["YEAR"].isnull()] = ":*"
    
    rmargin = utils.aggregate_years(rmargin, years_map, agg_method="mean")
    rmargin["VALUE"] = rmargin["VALUE"].round(DECIMALS)
    rmargin.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"ReserveMargin.csv", index=False) 
    
    
    rmargintagfuel = raw_data[["REGION","FUEL",
                          "YEAR",
                         "VALUE"]][raw_data["VARIABLE"]=="ReserveMarginTagFuel"]
    rmargintagfuel["REGION"][rmargintagfuel["REGION"].isnull()] = ":*"
    
    # FIXME: introduce some methodology for changing reserve marg. - 
    # necessary
    rmargintagfuel["YEAR"][rmargintagfuel["YEAR"].isnull()] = ":*"
    # FIXME: this aggregation would need to involve some rounding to 1/0 (?)
    rmargintagfuel = utils.aggregate_years(rmargintagfuel,
                                           years_map, agg_method="mean")
    rmargintagfuel.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"ReserveMarginTagFuel.csv", index=False) 
    
 


    rmargintagtech = raw_data[["REGION","TECHNOLOGY",
                          "YEAR",
                         "VALUE"]][raw_data["VARIABLE"]=="ReserveMarginTagTechnology"]
    rmargintagtech["REGION"][rmargintagtech["REGION"].isnull()] = ":*"
    
    # FIXME: introduce some methodology for changing reserve marg.(necessary!?)
    rmargintagtech["YEAR"][rmargintagtech["YEAR"].isnull()] = ":*"
    # FIXME: this aggregation would need to involve some rounding to 1/0 (?)  
    rmargintagtech = utils.aggregate_years(rmargintagtech, years_map, agg_method="mean")
    rmargintagtech["VALUE"] = rmargintagtech["VALUE"].round(DECIMALS)
    rmargintagtech.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"ReserveMarginTagTechnology.csv", index=False) 
    

#%%
# discount rate

    # simply get discount rate, if given, from raw input values
    # and save to csv
    # unit: -    

    # load raw data set
    raw_data = pd.read_csv(snakemake.input.path_set_disc) 
    
    drate = raw_data[["SCENARIO","REGION",
                         "VALUE"]][raw_data["VARIABLE"]=="DiscountRate"]
    drate = drate.loc[drate["SCENARIO"]==snakemake.params.dic["scen_discount_r"]]
    drate = drate.drop("SCENARIO",axis=1)
    drate["REGION"][drate["REGION"].isnull()] = ":*"
    
    drate["VALUE"] = drate["VALUE"].round(DECIMALS)
    drate.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"DiscountRate.csv", index=False)
    

    idrate = raw_data[["SCENARIO","REGION","TECHNOLOGY",
                         "VALUE"]][raw_data["VARIABLE"]=="DiscountRateIdv"]
    
    idrate = idrate.loc[idrate["SCENARIO"]==snakemake.params.dic["scen_discount_r"]]
    idrate = idrate.drop("SCENARIO",axis=1)
    if idrate.empty or drate.empty:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_discount_r']}' is currently"
                " not a valid option of the 'scen_discount_r' parameter."
                " Valid options are 'default' or 'ten_hurdle'."
            )  
        
    idrate["REGION"][idrate["REGION"].isnull()] = ":*"
    
    idrate["VALUE"] = idrate["VALUE"].round(DECIMALS)
    idrate.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"DiscountRateIdv.csv", index=False) 

#%%
# emission penalty

    # simply get emission penalty, if given, from raw input values, this might
    # later be based on other calculations
    # and save to csv
    # unit: million [Baseyear]£/kt/a
    
    # load raw data set
    raw_data = pd.read_csv(snakemake.input.path_set_em)
    
    empenlimit = raw_data[["REGION","EMISSION",
                          "YEAR",
                         "VALUE"]][raw_data["VARIABLE"]=="EmissionsPenalty"]
    empenlimit["REGION"][empenlimit["REGION"].isnull()] = ":*"
    
    # FIXME: introduce some methodology for changing emission penalty (necessary!?)
    empenlimit["YEAR"][empenlimit["YEAR"].isnull()] = ":*"
    empenlimit = utils.aggregate_years(empenlimit, years_map, agg_method="mean")
    empenlimit["VALUE"] = empenlimit["VALUE"].round(DECIMALS)
    empenlimit.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"EmissionsPenalty.csv", index=False) 
         
#%%
# fratoo params

    # simply get fratoo aggregation and disaggregation method for params from  
    # raw input values and save to csv
    # unit: -   

    # load raw data set
    raw_data = pd.read_csv(snakemake.input.path_set_ftpm)
    
    agg = raw_data[["PARAM",
                      "VALUE"]][raw_data["VARIABLE"]=="ft_param_agg"]

    agg.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"ft_param_agg.csv", index=False)
    
    disagg = raw_data[["PARAM",
                      "VALUE"]][raw_data["VARIABLE"]=="ft_param_disagg"]

    disagg.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"ft_param_disagg.csv", index=False) 
     
    
    
#%%
# fratoo entities


    # load area lookup and arrange dataframe for scale and affiliation and
    # save to csv
    # unit: -   
    # FIXME: consider if useful to get the LAD list from somewhere else
    
    # load sublocal aggregation
    slagg = pd.read_csv(snakemake.input.path_sublocal_agg,
                            index_col=["LSOA11CD"]).squeeze()
    
    lookup = utils.get_entity_lookup(levels=["LAD","CTRY"])
    lookup = lookup.set_index(["LAD23CD"])
    lookup = utils.groupby_LAD(lookup,disagg=slagg).max()
    lookup = lookup.reset_index()[["CTRY11CD","LAD23CD",slagg.name]]
    
    regions = list()
    scales = list()
    
    regions.extend(list(lookup[slagg.name]))
    scales.extend([3]*len(lookup[slagg.name]))
    
    regions.extend(list(lookup["LAD23CD"].drop_duplicates()))
    scales.extend([2]*len(lookup["LAD23CD"].drop_duplicates()))

    regions.extend(list(lookup["CTRY11CD"].drop_duplicates()))
    scales.extend([1]*len(lookup["CTRY11CD"].drop_duplicates()))    
    
    regions.extend(["UK"])
    scales.extend([0])

    scale = pd.DataFrame(data={"REGION":regions,"VALUE":scales})
   
    scale.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"ft_scale.csv", index=False)

    regions = list()
    affiliation = list()
    
    regions.extend(list(lookup[slagg.name]))
    affiliation.extend(list(lookup["LAD23CD"])) 
    
    regions.extend(list(lookup[["LAD23CD",
                                "CTRY11CD"]].drop_duplicates()["LAD23CD"]))
    affiliation.extend(list(lookup[["LAD23CD",
                                "CTRY11CD"]].drop_duplicates()["CTRY11CD"]))  
  
    regions.extend(list(lookup["CTRY11CD"].drop_duplicates()))
    affiliation.extend(["UK"]*len(lookup["CTRY11CD"].drop_duplicates()))    
    
    regions.extend(["UK"])
    affiliation.extend([""])

    affil = pd.DataFrame(data={"REGION":regions,"VALUE":affiliation})
   
    affil.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"ft_affiliation.csv", index=False)
    
    


#%%
# demand parameters


    # get demand data, annual and profile, from processed data and save
    # unit: TJ | -
    # FIXME: rearrange so columns are in right order (timeslice before year)
    
    # load sublocal aggregation
    slagg = pd.read_csv(snakemake.input.path_sublocal_agg,
                            index_col=["LSOA11CD"]).squeeze()
    
    # load annual demand data    
    dw_anndem_SH = (pd.read_csv(snakemake.input.path_dw_ann_dem_sh,
                             index_col=["LSOA11CD","PROPERTY_TYPE"])
                 + pd.read_csv(snakemake.input.path_dw_ann_dem_sh_nb,
                               index_col=["LSOA11CD","PROPERTY_TYPE"]))
    nd_anndem_SH = (pd.read_csv(snakemake.input.path_nd_ann_dem_sh,
                  index_col=["LSOA11CD","PROPERTY_TYPE"])
                 + pd.read_csv(snakemake.input.path_nd_ann_dem_sh_nb,
                               index_col=["LSOA11CD","PROPERTY_TYPE"])).reset_index()
    
    dw_anndem_HW = pd.read_csv(snakemake.input.path_dw_ann_dem_hw,
                               index_col=["LSOA11CD","PROPERTY_TYPE"])
    nd_anndem_HW =  pd.read_csv(snakemake.input.path_nd_ann_dem_hw)
    
    nd_anndem_SH["PROPERTY_TYPE"]="Non-domestic"
    nd_anndem_HW["PROPERTY_TYPE"]="Non-domestic"
    
    nd_anndem_SH = nd_anndem_SH.groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
    nd_anndem_HW = nd_anndem_HW.groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
    anndem_sh = pd.concat([dw_anndem_SH,nd_anndem_SH])
    anndem_hw = pd.concat([dw_anndem_HW,nd_anndem_HW])
    
    
    dw_anndem = pd.read_csv(snakemake.input.path_dw_ann_dem_nhg,
                             index_col=["LSOA11CD","PROPERTY_TYPE"])
    nd_anndem = pd.read_csv(snakemake.input.path_nd_ann_dem_nhg)
    nd_anndem["PROPERTY_TYPE"]="Non-domestic"
    nd_anndem = nd_anndem.groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
    anndem_nhg = pd.concat([dw_anndem,nd_anndem])
    
    dw_anndem = pd.read_csv(snakemake.input.path_dw_ann_dem_nhe,
                             index_col=["LSOA11CD","PROPERTY_TYPE"])
    nd_anndem = pd.read_csv(snakemake.input.path_nd_ann_dem_nhe)
    nd_anndem["PROPERTY_TYPE"]="Non-domestic"
    nd_anndem = nd_anndem.groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
    anndem_nhe = pd.concat([dw_anndem,nd_anndem])
    
    # group LSOA data by sublocal aggregation
    anndem_sh = utils.groupby_LAD(anndem_sh,disagg=slagg).sum()
    anndem_hw = utils.groupby_LAD(anndem_hw,disagg=slagg).sum()
    # calculate without disaggregation for weighting profile below
    anndem_nhe_w = utils.groupby_LAD(anndem_nhe,
                                     disagg=slagg).sum()
    anndem_nhg_w = utils.groupby_LAD(anndem_nhg,
                                     disagg=slagg).sum()
    # calculate for demand values
    anndem_nhe = utils.groupby_LAD(anndem_nhe,
                                   disagg=slagg).sum().groupby(slagg.name).sum()
    anndem_nhg = utils.groupby_LAD(anndem_nhg,
                                   disagg=slagg).sum().groupby(slagg.name).sum()

    

    rdict_sh = { pt+hh+ii:"SH"+"DDDE"+hh[1:]+ii if pt.startswith("Detached") else
                "SH"+"DDSD"+hh[1:]+ii if pt.startswith("Semi-detached") else
                "SH"+"DDTE"+hh[1:]+ii if pt.startswith("Terraced") else
                "SH"+"DDFL"+hh[1:]+ii if pt.startswith("Flats") else
                "SH"+"DNDO" if pt.startswith("Non-domestic") else 0
                for pt in ["Detached","Semi-detached",
                           "Terraced","Flats","Non-domestic"]
                for hh in ["","|OO","|RP","|RS","|O","|S","|P"]
                for ii in ([""]+[str(i) for i in range(0,10)])}
    rdict_hw = { pt+hh+ii:"HW"+"DDDE"+hh[1:]+ii if pt.startswith("Detached") else
                "HW"+"DDSD"+hh[1:]+ii if pt.startswith("Semi-detached") else
                "HW"+"DDTE"+hh[1:]+ii if pt.startswith("Terraced") else
                "HW"+"DDFL"+hh[1:]+ii if pt.startswith("Flats") else
                "HW"+"DNDO" if pt.startswith("Non-domestic") else 0
                for pt in ["Detached","Semi-detached",
                           "Terraced","Flats","Non-domestic"]
                for hh in ["","|OO","|RP","|RS","|O","|S","|P"]
                for ii in ([""]+[str(i) for i in range(0,10)])}

    # rename property types to respective heat demand fuels
    anndem_sh = anndem_sh.rename(index=rdict_sh)
    anndem_hw = anndem_hw.rename(index=rdict_hw)
    
    anndem_nhe["FUEL"] = "ELDSER"
    anndem_nhe = anndem_nhe.set_index(["FUEL"], append=True)
    anndem_nhg["FUEL"] = "NGDSER"
    anndem_nhg = anndem_nhg.set_index(["FUEL"], append=True)
    
    # convert from GJ to TJ
    anndem_sh = anndem_sh / 1000
    anndem_hw = anndem_hw / 1000
    anndem_nhe = anndem_nhe / 1000
    anndem_nhg = anndem_nhg / 1000
    
    # rename indices/columns, stack and save to file
    anndem_sh.index.names = ["REGION","FUEL"]
    anndem_sh.columns.name="YEAR"

    anndem_sh = anndem_sh.stack()
    anndem_sh.name="VALUE"
    anndem_sh = anndem_sh.reset_index()

    anndem_hw.index.names = ["REGION","FUEL"]
    anndem_hw.columns.name="YEAR"

    anndem_hw = anndem_hw.stack()
    anndem_hw.name="VALUE"
    anndem_hw = anndem_hw.reset_index()

    anndem_nhe.index.names = ["REGION","FUEL"]
    anndem_nhe.columns.name="YEAR"

    anndem_nhe = anndem_nhe.stack()
    anndem_nhe.name="VALUE"
    anndem_nhe = anndem_nhe.reset_index()
    
    anndem_nhg.index.names = ["REGION","FUEL"]
    anndem_nhg.columns.name="YEAR"

    anndem_nhg = anndem_nhg.stack()
    anndem_nhg.name="VALUE"
    anndem_nhg = anndem_nhg.reset_index()
    
    anndem = pd.concat([anndem_sh,anndem_hw,anndem_nhe,anndem_nhg])
    
    anndem = utils.aggregate_years(anndem, years_map, agg_method="sum")
    anndem["VALUE"] = anndem["VALUE"].round(DECIMALS)
    anndem.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"SpecifiedAnnualDemand.csv", index=False)
    


    # load demand profile data
    demand_pf_sh = pd.read_csv(snakemake.input.path_tperiods_sh,
                              index_col=["LAD23CD","PROPERTY_TYPE"])
    demand_pf_hw = pd.read_csv(snakemake.input.path_tperiods_hw,
                              index_col=["LAD23CD","PROPERTY_TYPE"])
    demand_pf_nhe = pd.read_csv(snakemake.input.path_tperiods_nhe,
                                index_col=["LAD23CD","PROPERTY_TYPE"])
    demand_pf_nhg = pd.read_csv(snakemake.input.path_tperiods_nhg,
                                index_col=["LAD23CD","PROPERTY_TYPE"])    
    # disaggregate based on sublocal aggregation
    demand_pf_sh = utils.groupby_LAD(demand_pf_sh,disagg=slagg).mean().squeeze()
    demand_pf_hw = utils.groupby_LAD(demand_pf_hw,disagg=slagg).mean().squeeze()
    demand_pf_nhe = utils.groupby_LAD(demand_pf_nhe,disagg=slagg).mean().squeeze()
    demand_pf_nhg = utils.groupby_LAD(demand_pf_nhg,disagg=slagg).mean().squeeze()
    
    # rename property types to respective heat demand fuels
    demand_pf_sh = demand_pf_sh.rename(index=rdict_sh)
    demand_pf_hw = demand_pf_hw.rename(index=rdict_hw)  
    
    # FIXME: this weighting would should ideally be not with average over years
    # but calculated and applied for each year separately 
    nhe_w = (anndem_nhe_w.mean(axis=1)
             /anndem_nhe_w.mean(axis=1).groupby(slagg.name).sum())
    nhg_w = (anndem_nhg_w.mean(axis=1)
             /anndem_nhg_w.mean(axis=1).groupby(slagg.name).sum())
    demand_pf_nhe = demand_pf_nhe.multiply(nhe_w,axis=0).groupby(slagg.name).sum()
    demand_pf_nhg = demand_pf_nhg.multiply(nhg_w,axis=0).groupby(slagg.name).mean()
    
    demand_pf_nhe["FUEL"] = "ELDSER"
    demand_pf_nhe = demand_pf_nhe.set_index(["FUEL"], append=True)
    demand_pf_nhg["FUEL"] = "NGDSER"
    demand_pf_nhg = demand_pf_nhg.set_index(["FUEL"], append=True)  
    
    
    # rename indices/columns
    demand_pf_sh.index.names = ["REGION","FUEL"]
    demand_pf_sh.columns.name="TIMESLICE"
    
    demand_pf_hw.index.names = ["REGION","FUEL"]
    demand_pf_hw.columns.name="TIMESLICE"
    
    demand_pf_nhe.index.names = ["REGION","FUEL"]
    demand_pf_nhe.columns.name="TIMESLICE"

    demand_pf_nhg.index.names = ["REGION","FUEL"]
    demand_pf_nhg.columns.name="TIMESLICE"    

    
    # stack, rename and save to file
    demand_pf_sh = demand_pf_sh.stack().to_frame(name="VALUE")
    demand_pf_hw = demand_pf_hw.stack().to_frame(name="VALUE")
    demand_pf_nhe = demand_pf_nhe.stack().to_frame(name="VALUE")
    demand_pf_nhg = demand_pf_nhg.stack().to_frame(name="VALUE")
    
    # add year index level and set assuming no changes over years
    demand_pf_sh = demand_pf_sh.assign(YEAR=":*").set_index("YEAR", append=True) 
    demand_pf_sh = demand_pf_sh.reset_index()
    demand_pf_hw = demand_pf_hw.assign(YEAR=":*").set_index("YEAR", append=True) 
    demand_pf_hw = demand_pf_hw.reset_index()
    demand_pf_nhe = demand_pf_nhe.assign(YEAR=":*").set_index("YEAR", append=True) 
    demand_pf_nhe = demand_pf_nhe.reset_index()
    demand_pf_nhg = demand_pf_nhg.assign(YEAR=":*").set_index("YEAR", append=True) 
    demand_pf_nhg = demand_pf_nhg.reset_index()   
    
    demand_pf = pd.concat([demand_pf_sh,demand_pf_hw,demand_pf_nhe,demand_pf_nhg])
    
    demand_pf = utils.aggregate_years(demand_pf, years_map, agg_method="mean")
    demand_pf["VALUE"] = demand_pf["VALUE"].round(DECIMALS)
    demand_pf.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"SpecifiedDemandProfile.csv", index=False)

#%%
# year split

    # load timeslice length from processed data and save
    # unit: -
    
    # load timeslice length
    ts_length = pd.read_csv(snakemake.input.path_tperiods_len,
                            index_col=["TIMESLICE"])
    
    # add year index level and set assuming no changes over years
    ts_length = ts_length.assign(YEAR=":*").set_index("YEAR", append=True)
    ts_length = ts_length.reset_index()
    
    ts_length = utils.aggregate_years(ts_length, years_map, agg_method="sum")
    ts_length["VALUE"] = ts_length["VALUE"].round(DECIMALS+1)
    # save to file 
    ts_length.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"YearSplit.csv", index=False)  


#%%
# sets

    # get sets that are defined in raw input file and save to files and
    # also create the ones derived through scripts (regions, timeslices)
    # and save to csv's
    # unit: -   

    # load raw data set
    sets = pd.read_csv(snakemake.input.path_set_sets)
    
    # process all sets defined in raw input data file
    setsl = sets["SET"].unique()
    for s in setsl:
        se = sets["VALUE"][sets["SET"]==s]
        #se.name = "VALUE"
        # explode technologies if household disaggregated
        if s=="YEAR":
            if isinstance(years_map,pd.Series):
                se = se.astype(int).map(years_map).drop_duplicates().astype(str)
            
        if s=="TECHNOLOGY":
            if snakemake.params.dic["scen_hh_disagg"]=="T":
                se = utils.explode_hh(se.to_frame(),
                                      ["OO","RP","RS"],
                                      col="VALUE")
            elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
                se = utils.explode_hh(se.to_frame(),
                                          [t+str(i) for t in ["O","P","S"]
                                           for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                          col="VALUE")
            elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
                pass
            else:
                raise NotImplementedError(
                        f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                        " not a valid option of the 'scen_hh_disagg' parameter."
                        " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
                    )
        if s=="FUEL":
            if snakemake.params.dic["scen_hh_disagg"]=="T":
                se = utils.explode_hh(se.to_frame(),
                                      ["OO","RP","RS"],
                                      col="VALUE")
            elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
                se = utils.explode_hh(se.to_frame(),
                                          [t+str(i) for t in ["O","P","S"]
                                           for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                          col="VALUE")
            elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
                pass
            else:
                raise NotImplementedError(
                        f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                        " not a valid option of the 'scen_hh_disagg' parameter."
                        " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
                    )  
        se.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                      +"/data/"+s+".csv", index=False) 
        
    # deduce UDC set from data
    udc = pd.read_csv(snakemake.input.path_heat_tech_udc_tag,
                      usecols=["UDC"],
                      squeeze=True)
    udc.name = "VALUE"
    udc = udc.drop_duplicates()
    udc.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"UDC.csv", index=False)  
    
    # deduce set of timeslices from timeslices length data and save to csv 
    timeslice = pd.read_csv(snakemake.input.path_tperiods_len,
                            usecols=["TIMESLICE"], squeeze=True)
    timeslice.name = "VALUE"
    timeslice.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"TIMESLICE.csv", index=False)  
    
    # use look up to create regions and save to csv
    # FIXME: potentially derive this from elsewhere - regions might not equal
    # lookup
    
    # load sublocal aggregation
    slagg = pd.read_csv(snakemake.input.path_sublocal_agg,
                            index_col=["LSOA11CD"]).squeeze()
    
    lookup = utils.get_entity_lookup(levels=["LAD","CTRY"])
    lookup = lookup.set_index(["LAD23CD"])
    lookup = utils.groupby_LAD(lookup,disagg=slagg).max()
    lookup = lookup.reset_index()[["CTRY11CD","LAD23CD",slagg.name]]

    regions = list()
    regions.extend(list(lookup[slagg.name]))
    regions.extend(list(lookup["LAD23CD"].drop_duplicates()))
    regions.extend(list(lookup["CTRY11CD"].drop_duplicates()))    
    regions.extend(["UK"])

    regions = pd.Series(data=regions, name="VALUE")
   
    regions.to_csv(snakemake.params.fdir+snakemake.params.dic["name"]
                  +"/data/"+"REGION.csv", index=False)    

#%%
# datapackage

    # create frictionless datapackage from model input files by deducing meta
    # data file
    
    # copy data to folder with scripts (fl does not allow to go backwards in
    # directory)
    # FIXME: implement this this in a cleaner way
    tdir = "./"+snakemake.params.dic["name"]+"/data"
    if os.path.exists(tdir):
        shutil.rmtree(tdir)
    shutil.copytree(snakemake.params.fdir+snakemake.params.dic["name"]+"/data",
                    tdir)
    
    # create meta data package file and adjust
    package = fl.describe(tdir+"/*.csv")
    package.title = snakemake.params.dic["name"]
    package.description = ("This is a UK-MOSEM version running "
                           +snakemake.params.dic["name"]+".")
    
    
    # rename resources with right capitalization (based on filenames)
    # and set keys/indices
    # FIXME: do this properly, maybe check how OSeMOSYS toolbox creates its 
    # packages, also involves checking on fratoo and define what is the proper
    # format/requirements for the package
    for r in package.resource_names:
        name = package.get_resource(r).path.split("/")[-1][:-4]
        package.get_resource(r).name = name
        package.get_resource(name).path = "/".join(package.get_resource(name).path.split("/")[1:])
    for r in package.resource_names:    
        ind = list(pd.read_csv(tdir+"/"+r+".csv").columns)
        ind.remove("VALUE")
        package.get_resource(r).schema.primary_key = ind
         
    
    # save file to original input data files' directory
    # FIXME: implement this in a way that works with different PATHF paths and
    # on different operating systems
    package.to_json(snakemake.params.fdir+snakemake.params.dic["name"]
                    +"/datapackage.json")
    
    # delete temporary data
    if os.path.exists(tdir):
        shutil.rmtree(tdir)
#%%
# model equations

    shutil.copy(snakemake.params.rdir+snakemake.params.dic["model_eq"],
                snakemake.params.fdir+snakemake.params.dic["name"]+"/"
                +snakemake.params.dic["model_eq"])


#%%
# built milestone

    # create file signalling end for snakemake model building work flow
    # (given above filepaths do not include wildcards)
    
    pd.DataFrame().to_csv(snakemake.output.path_built_ms)
    
