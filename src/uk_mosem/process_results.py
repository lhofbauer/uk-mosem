"""
Script to process results


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license

"""

import sys
import os
import logging

import pandas as pd


import fratoo as ft
import utils

ft.set_verbosity("INFO")

logger = logging.getLogger(__name__)

def process_results(model=None,params=None,xscale=None,years_map=None,rescapacity=False,
                   discounted=True,
                   technologies=None,demands=None,names=None ):
        
        
        
    logger.info("Processing results")
   
    if params is None:
        return
    
    results = model.results
   
    for i in range(len(results)):

        
        if "extendedcost" in params:
               
            ### calculate basic variables to be used for calculations                                               
            results[i]["UseByTechnology"] = ((results[i]["RateOfActivity"]
                                              *results[i]["InputActivityRatio"]).dropna()
                                             *results[i]["YearSplit"]).groupby(["RUN",
                                                                                "REGION",
                                                                                "TECHNOLOGY",
                                                                                "YEAR",
                                                                                "TIMESLICE",
                                                                                "FUEL"]).sum()                                                                                                                
            results[i]["ProductionByTechnology"] = ((results[i]["RateOfActivity"]
                                                      *results[i]["OutputActivityRatio"]).dropna()
                                                     *results[i]["YearSplit"]).groupby(["RUN",
                                                                                             "REGION",
                                                                                             "TECHNOLOGY",
                                                                                             "YEAR",
                                                                                             "TIMESLICE",
                                                                                             "FUEL"]).sum()            
            results[i]["ProductionByTechnologyByMode"] = ((results[i]["RateOfActivity"]
                                                      *results[i]["OutputActivityRatio"]).dropna()
                                                     *results[i]["YearSplit"]).groupby(["RUN",
                                                                                             "REGION",
                                                                                             "TECHNOLOGY",
                                                                                             "YEAR",
                                                                                             "TIMESLICE",
                                                                                             "MODE_OF_OPERATION",
                                                                                             "FUEL"]).sum()
            
            
            # calculate yearly and model period activity per tech
            yact = ((results[i]["RateOfActivity"]
                    *results[i]["YearSplit"]).groupby(["RUN",
                                                      "REGION",
                                                      "TECHNOLOGY",
                                                      "YEAR"]).sum())
            # remove zeros (that would lead to inf's later)
            yact = yact[yact["VALUE"]>0]
                                                       
            mpact = yact.groupby(["RUN",
                                  "REGION",
                                  "TECHNOLOGY"]).sum()
            # remove zeros (that would lead to inf's later)
            mpact = mpact[mpact["VALUE"]>0]    

                                                    
                                                                                        
                                                                                        
            # set initial cost for energy carrier (ec) outputs                                                                           
            ecoutcost = results[i]["ProductionByTechnologyByMode"].copy()
            ecoutcost = ecoutcost[~(ecoutcost["VALUE"]==0)]
            ecoutcost.loc[:,"VALUE"] = 0
            
            # set output level for which cost are calculated (generally 1,
            # except in case of several output energy carriers in one mode)
            ol = results[i]["ProductionByTechnologyByMode"]
            ol = ol / ol.groupby([c for c in ol.index.names
                                  if c !='FUEL']).transform("sum")
            ol = ol.dropna()
            
            # calculate required activity level
            # *ol could be another vector that splits 'activity', i.e., 
            # also cost, differently across ec outputs of same mode
            # (and tech)
            act = (ol/results[i]["OutputActivityRatio"]*ol).dropna()
            
            # calculate fix cost
            if discounted:
                fixcost = ((results[i]["TotalCapacityAnnual"]
                         *results[i]["FixedCost"]).dropna()
                         /results[i]["DiscountFactorMid"])
            else:
                fixcost = (results[i]["TotalCapacityAnnual"]
                         *results[i]["FixedCost"]).dropna()
                       
            results[i]["CostFixedProcessed"] = fixcost.copy()
            
            # take into account costs that would otherwise be unnaccounted
            # because no activity is present in a certain year (or
            # during the entire operational life). For power technologies
            # distribute cost on other techs, for other techs, save them
            # in separate dataframe for later consideration  
            fixcost = fixcost.reorder_levels([n for n in yact.index.names])
            
            # get unaccounted cost
            unacc = fixcost[~fixcost.index.isin(yact.index)].dropna()
            # get accounted cost
            fixcost = fixcost[fixcost.index.isin(yact.index)]
            
            # get unaccounted power sector costs and groupby year
            unaccpp = unacc[unacc.index.get_level_values("TECHNOLOGY"
                            ).str.contains('|'.join(["PPSNAT","SPSNAT"]))]
            unaccpp= unaccpp.groupby(["RUN","REGION","YEAR"]).sum()
            
            # get accounted power sector cost to calculate distribution
            accpp = fixcost[fixcost.index.get_level_values("TECHNOLOGY"
                            ).str.contains('|'.join(["PPSNAT","SPSNAT"]))]
            # add unaccounted power sector cost to fixcost
            fixcost = fixcost.add(accpp/accpp.groupby(["RUN","REGION",
                                                        "YEAR"]).transform('sum')*unaccpp,
                                  fill_value=0)
            # create dataframe for still unaccounted cost
            fixcost_ua = unacc[~unacc.index.get_level_values("TECHNOLOGY"
                            ).str.contains('|'.join(["PPSNAT","SPSNAT"]))]
            
            # add unaccounted fixcost to fixcost of other years that are
            # accounted for
            # FIXME: improve this so cost are only allocated to relevant years,
            # maybe following a similar approach to capital cost
            fixcost = fixcost.add((fixcost/fixcost.groupby(["RUN",
                                            "REGION",
                                            "TECHNOLOGY"]).sum()
                       * fixcost_ua.groupby(["RUN","REGION","TECHNOLOGY"]).sum()).dropna(),
                                  fill_value=0)
            # remove reallocated cost from unaccounted cost
            fixcost_ua = fixcost_ua.sub(((fixcost/fixcost.groupby(["RUN",
                                            "REGION",
                                            "TECHNOLOGY"]).sum()
                       * fixcost_ua.groupby(["RUN","REGION","TECHNOLOGY"]).sum()).dropna().groupby(["RUN","REGION","TECHNOLOGY"]).sum())
                       *(fixcost_ua/fixcost_ua.groupby(["RUN",
                                                       "REGION",
                                                       "TECHNOLOGY"]).sum()).dropna(),fill_value=0).dropna()
            results[i]["CostFixedProcessed"] = pd.concat([fixcost.copy(),
                                                          fixcost_ua])
            

                
                
            # calculate fixcost per set activity level (relating to
            # one unit of output)
            fixcost = (fixcost / yact).dropna()
            fixcost = (fixcost * act).dropna()
            
            # FIXME: add emission penalty cost                             
            # calculate capital cost
            if discounted:
                capcost = (((results[i]["CapitalCost"]
                           * results[i]["CapitalRecoveryFactor"]
                           * results[i]["PvAnnuity"]
                           * results[i]["NewCapacity"]
                           /results[i]["DiscountFactor"]
                           ).subtract(
                           results[i]["DiscountedSalvageValue"],
                           fill_value=0)))
            else:
                # calculate SalvageCost (without incorporated discounting)
                years = results[i]["YEAR"].copy()
                years.columns=["YEAR"]
                years = years.set_index("YEAR")
                years["VALUE"] = years.index.get_level_values("YEAR")
                ltech = results[i]["OperationalLife"].index.get_level_values("TECHNOLOGY").unique()
                years = pd.concat([years]*len(ltech),names=["TECHNOLOGY"],keys=list(ltech))
                # calculate Pvannuity (without incorporated discounting)
                results[i]["PvAnnuityUndiscounted"] = results[i]["OperationalLife"].copy()
                
                results[i]["SalvageValueUndiscounted"] = (results[i]["CapitalCost"]
                                              * results[i]["CapitalRecoveryFactor"]
                                              * results[i]["PvAnnuityUndiscounted"]
                                              * results[i]["NewCapacity"]
                                              *(1-(years.max().squeeze()-years+1)
                                              /results[i]["OperationalLife"]).clip(lower=0))
                
                    
                capcost = (((results[i]["CapitalCost"]
                           * results[i]["CapitalRecoveryFactor"]
                           * results[i]["PvAnnuityUndiscounted"]
                           * results[i]["NewCapacity"]
                           ).subtract(
                           results[i]["SalvageValueUndiscounted"],
                           fill_value=0)))
                
                # adjust cost of capital and capital cost for gas grid
                if "jt" in params:
                    
                    # calculate weighted average capital cost for GAGR and apply
                    # across all areas, this will then also apply to
                    # residual capacity
                    cc = (results[i]["CapitalCost"].xs("GAGRTDIS00",level="TECHNOLOGY")
                          *results[i]["NewCapacity"].add(results[i]["ResidualCapacity"],fill_value=0).xs("GAGRTDIS00",level="TECHNOLOGY"))
                    cc = (cc.groupby(["RUN","YEAR"]).sum()
                          /results[i]["NewCapacity"].add(results[i]["ResidualCapacity"],fill_value=0).xs("GAGRTDIS00",level="TECHNOLOGY").groupby(["RUN","YEAR"]).sum())
                    cc = cc.fillna(0)
                    cc.loc[:,"TECHNOLOGY"] = "GAGRTDIS00"
                    cc = cc.set_index("TECHNOLOGY",append=True)
                    cc = results[i]["CapitalCost"].reset_index().merge(cc,
                                                         on=["RUN",
                                                             "TECHNOLOGY",
                                                             "YEAR"],
                                                         how="left",
                                                         suffixes=("","_"))
                    logger.info("2")
                    cc = cc.set_index(["RUN","REGION","TECHNOLOGY","YEAR"])
                    cc.loc[~cc["VALUE_"].isna(),"VALUE"] = cc.loc[~cc["VALUE_"].isna(),"VALUE_"]
                    cc = cc.drop("VALUE_",axis=1)
                    results[i]["CapitalCost"] = cc
                    logger.info("Adjusted capital cost")
                    
                    hurdle = pd.read_csv(snakemake.input.path_supply_hurdle,
                                         index_col=["TECHNOLOGY","YEAR"])
                    hurdle = hurdle.groupby(["TECHNOLOGY"]).mean()
                    hurdle.columns = ["VALUE"]
                    h = pd.read_csv(snakemake.input.path_set_supply)
                    h = h.loc[h["VARIABLE"]=="HurdleRate",["TECHNOLOGY",
                                                           "VALUE"]].squeeze()
                    hurdle = pd.concat([hurdle,
                                        h.set_index("TECHNOLOGY")])
                    
                    if "hl" in params:
                        hurdle["VALUE"] = hurdle["VALUE"]*0.8
                    if "hh" in params:
                        hurdle["VALUE"] = hurdle["VALUE"]*1.20
 
                    results[i]["DiscountRateIdv"] = results[i]["DiscountRateIdv"].reset_index().merge(hurdle,on="TECHNOLOGY",how="left",suffixes=("","_"))
                    results[i]["DiscountRateIdv"].loc[~results[i]["DiscountRateIdv"]["VALUE_"].isna(),
                                                      "VALUE"]=results[i]["DiscountRateIdv"].loc[~results[i]["DiscountRateIdv"]["VALUE_"].isna(),
                                                                                        "VALUE_"]
                    results[i]["DiscountRateIdv"] =  results[i]["DiscountRateIdv"].drop("VALUE_",axis=1)                                                                           
                    results[i]["DiscountRateIdv"] = results[i]["DiscountRateIdv"].set_index(["RUN","REGION","TECHNOLOGY"])                                                          
                                                                              
                    results[i]["CapitalRecoveryFactor"] = (1- (1 + results[i]["DiscountRateIdv"])**(-1))/(1-(1+results[i]["DiscountRateIdv"])**(-(results[i]["OperationalLife"])))
                    results[i]["PvAnnuityUndiscounted"] = results[i]["OperationalLife"].copy()
                    
                    results[i]["SalvageValueUndiscounted"] = (results[i]["CapitalCost"]
                                                  * results[i]["CapitalRecoveryFactor"]
                                                  * results[i]["PvAnnuityUndiscounted"]
                                                  * results[i]["NewCapacity"]
                                                  *(1-(years.max().squeeze()-years+1)
                                                  /results[i]["OperationalLife"]).clip(lower=0))    
                    capcost = (((results[i]["CapitalCost"]
                                * results[i]["CapitalRecoveryFactor"]
                                * results[i]["PvAnnuityUndiscounted"]
                                * results[i]["NewCapacity"]
                                ).subtract(
                                results[i]["SalvageValueUndiscounted"],
                                fill_value=0)))
            
                               
            results[i]["CostCapitalProcessed"] = capcost.copy()
            capcost = capcost.reorder_levels(["RUN",
                                                "REGION",
                                                "YEAR",
                                                "TECHNOLOGY"])

            # take into account costs that would otherwise be unnaccounted
            # because no activity is present in a certain year (or
            # during the entire operational life). For power technologies
            # distribute cost on other techs, for other techs, save them
            # in separate dataframe for later consideration

            # get unaccounted cost
            unacc = capcost[~capcost.index.droplevel(2).isin(mpact.index)].dropna()
            # get accounted cost
            capcost = capcost[capcost.index.droplevel(2).isin(mpact.index)]
            
            # get unaccounted power sector costs
            unaccpp = unacc[unacc.index.get_level_values("TECHNOLOGY"
                            ).str.contains('|'.join(["PPSNAT","SPSNAT"]))]
            unaccpp= unaccpp.groupby(["RUN","REGION","YEAR"]).sum()
            
            # get accounted power sector cost to calculate distribution
            accpp = capcost[capcost.index.get_level_values("TECHNOLOGY"
                            ).str.contains('|'.join(["PPSNAT","SPSNAT"]))]
            # add unaccounted power sector cost to capcost
            capcost = capcost.add(accpp/accpp.groupby(["RUN","REGION","YEAR"]).transform('sum')*unaccpp,
                                  fill_value=0)
            # create dataframe for still unaccounted cost
            capcost_ua = unacc[~unacc.index.get_level_values("TECHNOLOGY"
                            ).str.contains('|'.join(["PPSNAT","SPSNAT"]))]
            
            ### distribute capcost across years
            # arrange dataframe to reflect years new built capacity from
            # a certain year is in operation/existent
            capcostdist = results[i]["NewCapacity"].copy()
            capcostdist = capcostdist[capcostdist["VALUE"]>0]
            capcostdist = capcostdist.assign(Y=capcostdist.index.get_level_values("YEAR")).set_index('Y', append=True)
            capcostdist = capcostdist.unstack().droplevel(0,axis=1)
            capcostdist.columns = pd.to_datetime(capcostdist.columns,format='%Y')
            capcostdist = capcostdist.fillna(0)
            capcostdist = capcostdist.apply(lambda x: x.rolling(str(int(results[i]["OperationalLife"].loc[x.name[:-1],"VALUE"])*364)+"D",min_periods=0).sum(),axis=1)
            capcostdist.columns = capcostdist.columns.year
            
            # apply activity to capacity ratio for each year
            results[i]["TotalActivityAnnual"] = ((results[i]["RateOfActivity"]
                                                          *results[i]["YearSplit"]).dropna()
                                                          ).groupby(["RUN","REGION","TECHNOLOGY","YEAR"]).sum()
            # for certain technologies, spread cost equally across years
            # (i.e., set activity as capacity so fraction below will be 1)
            taa = results[i]["TotalActivityAnnual"].copy()
            
            taa.loc[taa.index.get_level_values("TECHNOLOGY").str.startswith(("H2PP",
                                                                              "BSPP",
                                                                              "AESP"))]=results[i]["TotalCapacityAnnual"]
            # capcostdist = capcostdist*(results[i]["TotalActivityAnnual"]
            #              /results[i]["TotalCapacityAnnual"]).dropna().unstack().fillna(0).droplevel(0,axis=1)
            capcostdist = capcostdist*(taa
                          /results[i]["TotalCapacityAnnual"]).unstack().fillna(0).droplevel(0,axis=1)
            # where there is no activity (row filled with 0 above) distributed
            # equally (add 1's)
            re = capcostdist[capcostdist.sum(skipna=False,axis=1)==0].stack()
            re[re.index.get_level_values("YEAR")<=re.index.get_level_values("Y")]=1
            capcostdist[capcostdist.sum(skipna=False,axis=1)==0] = re.unstack()

            capcostdist = capcostdist.divide(capcostdist.sum(axis=1),axis=0)
            capcostdist = capcostdist.dropna(axis = 0, how = 'all')
            capcostdist = capcostdist.reorder_levels(capcost.index.names)
            
            # use capcost distribution to distribute capcost across
            # years
            capcost = capcostdist.multiply(capcost["VALUE"],axis=0)
            capcost = capcost.stack()
            capcost = capcost.groupby(["RUN","REGION","TECHNOLOGY","Y"]).sum()
            capcost.index = capcost.index.rename("YEAR",level=3)
            capcost.name = "VALUE"
            capcost = capcost.to_frame()
            
            # do the same for base capital cost without redistributed capcost
            capcostdistproc = capcostdist
            results[i]["CostCapitalProcessed"] = results[i]["CostCapitalProcessed"].reorder_levels(capcostdist.index.names)
            results[i]["CostCapitalProcessed"] = capcostdist.multiply(results[i]["CostCapitalProcessed"]["VALUE"],axis=0)
            results[i]["CostCapitalProcessed"] = results[i]["CostCapitalProcessed"].stack()
            results[i]["CostCapitalProcessed"] = results[i]["CostCapitalProcessed"].groupby(["RUN","REGION","TECHNOLOGY","Y"]).sum()
            results[i]["CostCapitalProcessed"].index = results[i]["CostCapitalProcessed"].index.rename("YEAR",level=3)
            results[i]["CostCapitalProcessed"].name = "VALUE"
            results[i]["CostCapitalProcessed"] = results[i]["CostCapitalProcessed"].to_frame()
            
            # calculate and add capacity cost for residual capactity
            # if required, this is assumed to be undiscounted in any 
            # case
            # if jt, this takes into account the update CRF and PvAU from above
            if rescapacity:
                rescap = results[i]["ResidualCapacity"].copy()
                rccapcost = (((rescap["VALUE"]*(1/xscale["VALUE"]))
                             .groupby(["RUN",
                                       "REGION",
                                       "TECHNOLOGY"]).sum()
                             /results[i]["OperationalLife"]["VALUE"])
                             *results[i]["CapitalCost"].xs(
                                 min(results[i]["YEAR"]["VALUE"]),
                                 level="YEAR")["VALUE"]
                           * results[i]["CapitalRecoveryFactor"]["VALUE"]
                           * results[i]["PvAnnuityUndiscounted"]["VALUE"]).to_frame()

                
                rccapcost.columns = ["VALUE"]
                
                # calculate distribution
                # TotalActivityAnnual is reindexed and filled with 1 for nan
                # so that techologies with residual capacity that is not
                # used are also included (and separated below)
                rccapcostdist = rescap.copy()
                rccapcostdist = (rccapcostdist*(results[i]["TotalActivityAnnual"].reindex_like(rescap).fillna(results[i]["TotalCapacityAnnual"]*(1/xscale["VALUE"]))
                             /results[i]["TotalCapacityAnnual"]).dropna()).dropna()
                rccapcostdist = rccapcostdist/rccapcostdist.groupby(["RUN","REGION","TECHNOLOGY"]).transform("sum")
                rccapcost = (rccapcost*rccapcostdist).dropna()
                
                # add cost to basic cost with redistribution
                results[i]["CostCapitalProcessed"] = results[i]["CostCapitalProcessed"].add(rccapcost,fill_value=0)
                
                # take into account costs that would otherwise be unnaccounted
                # because no activity is present in a certain year (or
                # during the entire operational life). For power technologies
                # distribute cost on other techs, for other techs, save them
                # in separate dataframe for later consideration

                # get unaccounted cost
                noact = results[i]["TotalActivityAnnual"].reindex_like(rescap).fillna(0).groupby(["RUN","REGION","TECHNOLOGY"]).sum()
                noact = noact[noact["VALUE"]==0]
                rcunacc = rccapcost[rccapcost.index.droplevel(3).isin(noact.index)].dropna()
                # get accounted cost
                rccapcost = rccapcost[~rccapcost.index.droplevel(3).isin(noact.index)]
                # add accounted to capcost to get total based on NewCapacity
                # and ResidualCapacity
                capcost = capcost.add(rccapcost,fill_value=0)
                
                # get unaccounted power sector costs
                rcunaccpp = rcunacc[rcunacc.index.get_level_values("TECHNOLOGY"
                                ).str.contains('|'.join(["PPSNAT","SPSNAT"]))]
                rcunaccpp= rcunaccpp.groupby(["RUN","REGION","YEAR"]).sum()
                
                # get accounted power sector cost to calculate distribution
                accpp = capcost[capcost.index.get_level_values("TECHNOLOGY"
                                ).str.contains('|'.join(["PPSNAT","SPSNAT"]))]
                # add unaccounted power sector cost to capcost
                capcost = capcost.add(accpp/accpp.groupby(["RUN","REGION","YEAR"]).transform('sum')*rcunaccpp,
                                      fill_value=0)
                # create dataframe for still unaccounted cost
                rccapcost_ua = rcunacc[~rcunacc.index.get_level_values("TECHNOLOGY"
                                ).str.contains('|'.join(["PPSNAT","SPSNAT"]))]
                capcost_ua = pd.concat([capcost_ua,rccapcost_ua])
            

                # FIXME: the lines below could be used for a more detailed
                # cost calculation for residual caps
                # rescap = results[i]["ResidualCapacity"].copy()
                # rescap = rescap.unstack().droplevel(0,axis=1)
                # rescap = -rescap.diff(axis=1,periods=-1)
                # rescap = rescap.stack().to_frame()
                # #rescap = rescap.reset_index(level="YEAR")
                # rescap["BY"] = rescap.index.get_level_values("YEAR").values
                # rescap["BY"] = (rescap["BY"]
                #                    - results[i]["OperationalLife"]["VALUE"]
                #                    +1)
            
            # calculate capcost per activity
            capcost = (capcost / yact).dropna()

            # capital cost per output level
            capcost = (capcost * act).dropna() 
            
            # calculate variable cost
            if discounted:
                varcost = (act 
                           *results[i]["VariableCost"]
                           /results[i]["DiscountFactorMid"]).dropna()
                
                results[i]["CostVariableProcessed"] = ((results[i]["RateOfActivity"]
                                                        *results[i]["YearSplit"]
                                                        *results[i]["VariableCost"]
                                                        /results[i]["DiscountFactorMid"]
                                                        ).dropna()
                                                       ).groupby(["RUN",
                                                                  "REGION",
                                                                  "TECHNOLOGY",
                                                                  "YEAR"]).sum()
            else:
                # add other price components otherwise not captured
                if "jt" in params:
                    
                    # wood pellets (add profit based on industry average)

                    raw_data = pd.read_csv(snakemake.input.path_set_supply)
                    
                    pm = raw_data.loc[raw_data["VARIABLE"]=="ProfitMarginEnergyIndustry","VALUE"].squeeze()
                    
                    
                    results[i]["VariableCost"].loc[
                        (results[i]["VariableCost"].index
                        .get_level_values("TECHNOLOGY")=="BMEXSEXT00")] = (
                        results[i]["VariableCost"].loc[
                            (results[i]["VariableCost"].index
                            .get_level_values("TECHNOLOGY")=="BMEXSEXT00")]
                        /(1-pm))
                    results[i]["VariableCost"].loc[
                        (results[i]["VariableCost"].index
                        .get_level_values("TECHNOLOGY")=="BMTRTTRA00")] = (
                        results[i]["VariableCost"].loc[
                            (results[i]["VariableCost"].index
                            .get_level_values("TECHNOLOGY")=="BMTRTTRA00")]
                        /(1-pm))
                    
                    # gas and electricity (add profit and other cost based
                    # on Ofgem price cap methodology)
                    pcap_el = pd.read_excel(snakemake.input.path_ec_pcap,
                                         sheet_name="1b Historical level tables",
                                         skiprows=27,usecols="B:V",
                                         index_col=0,
                                         nrows=13
                                         ).dropna(axis=1,how="all")
                    pcap_ng = pd.read_excel(snakemake.input.path_ec_pcap,
                                         sheet_name="1b Historical level tables",
                                         skiprows=27,usecols="AT:BN",
                                         index_col=0,
                                         nrows=13
                                         ).dropna(axis=1,how="all")
                    pcap_ng.index.names = ["COST"]
                    pcap_ng = pcap_ng.replace({"-":0}).fillna(0)
                    pcap_ng.apply(pd.to_numeric)
                    pcap_el.index.names = ["COST"]
                    pcap_el = pcap_el.replace({"-":0}).fillna(0)
                    pcap_el.apply(pd.to_numeric)

                    # pcap_ng.columns = [c.split(".")[0] for c in pcap_ng.columns]
                    pcap_ng.columns = [c.split(" ")[1] for c in pcap_ng.columns]
                    pcap_el.columns = [c.split(" ")[1] for c in pcap_el.columns]
                    
                    pcap_ng.columns.names = ["YEAR"]
                    pcap_el.columns.names = ["YEAR"]
                    
                    # aggregate to yearly data
                    pcap_ng = pcap_ng.groupby("YEAR",axis=1).mean()
                    pcap_ng.loc[:,"2015"] = pcap_ng.loc[:,"2017"]
                    pcap_ng.loc[:,"2016"] = pcap_ng.loc[:,"2017"]
                    pcap_ng.columns = pcap_ng.columns.astype(int)
                    pcap_el = pcap_el.groupby("YEAR",axis=1).mean()
                    pcap_el.loc[:,"2015"] = pcap_el.loc[:,"2017"]
                    pcap_el.loc[:,"2016"] = pcap_el.loc[:,"2017"]
                    pcap_el.columns = pcap_el.columns.astype(int)
                    
                    con_ng = 12000
                    pcap_ng = pcap_ng.stack()
                    pcap_ng.name="VALUE"
                    pcap_ng = pcap_ng.reset_index()
                    pcap_ng.loc[:,"UNIT"] = "2017£"
                    pcap_ng.loc[pcap_ng["YEAR"]>2016,
                                "UNIT"] = pcap_ng["YEAR"].astype(str)+"£"
                    pcap_ng = utils.adjust_monetary_values(pcap_ng, 2015)
                    
                    
                    inc = ("CM","AA","PC","OC","SMNCC","PAAC","PAP","EBIT","HAP")
                    # FIXME: add climate change levy
                    pcap_ng_b = pcap_ng.loc[pcap_ng["COST"].isin(inc)].copy()
                    # scaling operational cost based on consumption
                    pcap_ng_b.loc[pcap_ng_b["COST"]=="OC",
                                  "VALUE"] = pcap_ng_b.loc[pcap_ng_b["COST"]=="OC",
                                                "VALUE"]*12000/(27778*1000)
                    pcap_ng_b = pcap_ng_b.groupby("YEAR").sum()/con_ng
                                                           
                    pcap_ng = pcap_ng.loc[pcap_ng["COST"].isin(inc)]
                    pcap_ng = pcap_ng.groupby("YEAR").sum()/con_ng
                    
                    con_el = 3100
                    pcap_el = pcap_el.stack()
                    pcap_el.name="VALUE"
                    pcap_el = pcap_el.reset_index()
                    pcap_el.loc[:,"UNIT"] = "2017£"
                    pcap_el.loc[pcap_el["YEAR"]>2016,
                                "UNIT"] = pcap_el["YEAR"].astype(str)+"£"
                    pcap_el = utils.adjust_monetary_values(pcap_el, 2015)
                    
                    # not including "CM" as in power sector costs
                    inc = ("AA","PC","OC","SMNCC","PAAC","PAP","EBIT","HAP")
                    
                    pcap_el_b = pcap_el.loc[pcap_el["COST"].isin(inc)].copy()
                    # scaling operational cost based on consumption
                    pcap_el_b.loc[pcap_el_b["COST"]=="OC",
                                  "VALUE"] = pcap_el_b.loc[pcap_el_b["COST"]=="OC",
                                                "VALUE"]*12000/(20000*1000)
                    pcap_el_b = pcap_el_b.groupby("YEAR").sum()/con_el
                    
                    pcap_el = pcap_el.loc[pcap_el["COST"].isin(inc)]
                    pcap_el = pcap_el.groupby("YEAR").sum()/con_el
                    
                    # recalculate as 10^6 GBP/PJ
                    pcap_ng = pcap_ng/10**6/(3.6/10**6)
                    pcap_el = pcap_el/10**6/(3.6/10**6)
                    pcap_ng_b = pcap_ng_b/10**6/(3.6/10**6)
                    pcap_el_b = pcap_el_b/10**6/(3.6/10**6)
                    
                    # extend 2020 values towards future years
                    
                    for y in range(2020,2061):
                        pcap_ng.loc[y] = pcap_ng.loc[2020]
                        pcap_el.loc[y] = pcap_el.loc[2020]
                        pcap_ng_b.loc[y] = pcap_ng_b.loc[2020]
                        pcap_el_b.loc[y] = pcap_el_b.loc[2020] 
                        
                    if "ll" in params:
                        for y in range(2021,2061):
                            pcap_ng.loc[y] = pcap_ng.loc[2020]-pcap_ng.loc[2020]*0.2*(y-2020)/(2050-2020)
                            pcap_el.loc[y] = pcap_el.loc[2020]-pcap_el.loc[2020]*0.2*(y-2020)/(2050-2020)
                            pcap_ng_b.loc[y] = pcap_ng_b.loc[2020]-pcap_ng_b.loc[2020]*0.2*(y-2020)/(2050-2020)
                            pcap_el_b.loc[y] = pcap_el_b.loc[2020]-pcap_el_b.loc[2020]*0.2*(y-2020)/(2050-2020)
                    if "lh" in params:
                        for y in range(2021,2061):
                             pcap_ng.loc[y] = pcap_ng.loc[2020]+pcap_ng.loc[2020]*0.2*(y-2020)/(2050-2020)
                             pcap_el.loc[y] = pcap_el.loc[2020]+pcap_el.loc[2020]*0.2*(y-2020)/(2050-2020)
                             pcap_ng_b.loc[y] = pcap_ng_b.loc[2020]+pcap_ng_b.loc[2020]*0.2*(y-2020)/(2050-2020)
                             pcap_el_b.loc[y] = pcap_el_b.loc[2020]+pcap_el_b.loc[2020]*0.2*(y-2020)/(2050-2020)
  
                    # aggregate years
                    pcap_ng =utils.aggregate_years(pcap_ng.reset_index(),
                                                   years_map,
                                                   agg_method="mean")
                    pcap_ng = pcap_ng.set_index("YEAR")
                    pcap_el =utils.aggregate_years(pcap_el.reset_index(),
                                                   years_map,
                                                   agg_method="mean")
                    pcap_el = pcap_el.set_index("YEAR")
                    pcap_ng_b =utils.aggregate_years(pcap_ng_b.reset_index(),
                                                   years_map,
                                                   agg_method="mean")
                    pcap_ng_b = pcap_ng_b.set_index("YEAR")
                    pcap_el_b =utils.aggregate_years(pcap_el_b.reset_index(),
                                                   years_map,
                                                   agg_method="mean")
                    pcap_el_b = pcap_el_b.set_index("YEAR")
                    

                    # add to dist network variable cost
                    results[i]["VariableCost"].loc[
                        results[i]["VariableCost"].index
                        .get_level_values("TECHNOLOGY")=="ELGRTDIS00"] = (
                        results[i]["VariableCost"].loc[
                            results[i]["VariableCost"].index
                            .get_level_values("TECHNOLOGY")=="ELGRTDIS00"] + pcap_el)
                    
                    results[i]["VariableCost"].loc[
                        results[i]["VariableCost"].index
                        .get_level_values("TECHNOLOGY")=="GAGRTDIS00"] = (
                        results[i]["VariableCost"].loc[
                            results[i]["VariableCost"].index
                            .get_level_values("TECHNOLOGY")=="GAGRTDIS00"] + pcap_ng)

                    
                    # add non-domestic rates to non-domestic (DH) users
                    results[i]["VariableCost"].loc[
                        results[i]["VariableCost"].index
                        .get_level_values("TECHNOLOGY")=="NGBOSDIS00"] = (
                        results[i]["VariableCost"].loc[
                            results[i]["VariableCost"].index
                            .get_level_values("TECHNOLOGY")=="NGBOSDIS00"].add(
                        results[i]["InputActivityRatio"].loc[
                            results[i]["InputActivityRatio"].index
                            .get_level_values("TECHNOLOGY")=="NGBOSDIS00"].droplevel("FUEL")
                        *pcap_ng_b,fill_value=0))
                                
                    results[i]["VariableCost"].loc[
                        results[i]["VariableCost"].index
                        .get_level_values("TECHNOLOGY")=="NGCHSDIS00"] = (
                        results[i]["VariableCost"].loc[
                            results[i]["VariableCost"].index
                            .get_level_values("TECHNOLOGY")=="NGCHSDIS00"].add(
                        results[i]["InputActivityRatio"].loc[
                            results[i]["InputActivityRatio"].index
                            .get_level_values("TECHNOLOGY")=="NGCHSDIS00"].droplevel("FUEL")
                        *pcap_ng_b,fill_value=0))
                                
                    results[i]["VariableCost"].loc[
                        results[i]["VariableCost"].index
                        .get_level_values("TECHNOLOGY")=="RWHPSDIS00"] = (
                        results[i]["VariableCost"].loc[
                            results[i]["VariableCost"].index
                            .get_level_values("TECHNOLOGY")=="RWHPSDIS00"].add(
                        results[i]["InputActivityRatio"].loc[
                            results[i]["InputActivityRatio"].index
                            .get_level_values("TECHNOLOGY")=="RWHPSDIS00"].droplevel("FUEL")
                        *pcap_el_b,fill_value=0))
                    results[i]["VariableCost"].loc[
                        results[i]["VariableCost"].index
                        .get_level_values("TECHNOLOGY")=="SWHPSDIS00"] = (
                        results[i]["VariableCost"].loc[
                            results[i]["VariableCost"].index
                            .get_level_values("TECHNOLOGY")=="SWHPSDIS00"].add(
                        results[i]["InputActivityRatio"].loc[
                            results[i]["InputActivityRatio"].index
                            .get_level_values("TECHNOLOGY")=="SWHPSDIS00"].droplevel("FUEL")
                        *pcap_el_b,fill_value=0))
                    
                varcost = (act 
                           *results[i]["VariableCost"]).dropna()
                results[i]["CostVariableProcessed"] = ((results[i]["RateOfActivity"]
                                                        *results[i]["YearSplit"]
                                                        *results[i]["VariableCost"]
                                                        ).dropna()
                                                       ).groupby(["RUN",
                                                                  "REGION",
                                                                  "TECHNOLOGY",
                                                                  "YEAR"]).sum()

                    
            results[i]["CostTotalProcessed"] =  (results[i]["CostCapitalProcessed"].add(
                                                results[i]["CostFixedProcessed"], fill_value=0
                                                ).add(
                                                results[i]["CostVariableProcessed"], fill_value=0
                                                ))                                                  
                                                                  
            # rearrange to keep output ec for calculation of upstream cost
            act.index = act.index.set_names("FUEL_OUT",level="FUEL")
            
            # iterate through energy carrier-technology chain until no changes
            # observed (and max 20 iterations)
            for ii in range(0,20):

                # make copy to check if changed in the end of the loop
                cc = ecoutcost.copy() 
                
                # rearrange fuel column, i.e., add fuel region
                ecoutcost = ecoutcost.reset_index()
                ecoutcost.loc[~ecoutcost["FUEL"].str.contains(":"),"FUEL"] = ecoutcost.loc[~ecoutcost["FUEL"].str.contains(":"),"REGION"]+":"+ecoutcost.loc[~ecoutcost["FUEL"].str.contains(":"),"FUEL"]
                ecoutcost = ecoutcost.set_index([c for c in ecoutcost.columns if c !="VALUE"])  


                ### calculate upstream ec cost

                # calculate ec input demand
                ecidem = (act * results[i]["InputActivityRatio"]).dropna()

                # calculate fraction of ec produced by tech & mode
                prodfrac = results[i]["ProductionByTechnologyByMode"].reset_index()
                prodfrac.loc[~prodfrac["FUEL"].str.contains(":"),"FUEL"] = prodfrac.loc[~prodfrac["FUEL"].str.contains(":"),"REGION"]+":"+prodfrac.loc[~prodfrac["FUEL"].str.contains(":"),"FUEL"]
                prodfrac = prodfrac.set_index([c for c in prodfrac.columns if c !="VALUE"])
                #prodfrac = prodfrac.droplevel("REGION")

                # for storage, this is calculated for each timeslice but
                # as fraction to total year, within each timeslice
                prodfracsto = prodfrac / prodfrac.groupby([c for c in prodfrac.index.names
                                                        if ((c !="MODE_OF_OPERATION")
                                                        and (c != "TECHNOLOGY")
                                                        and (c != "REGION")
                                                        and (c != "TIMESLICE"))]).transform('sum')
                prodfrac = prodfrac / prodfrac.groupby([c for c in prodfrac.index.names
                                                        if ((c !="MODE_OF_OPERATION")
                                                        and (c != "TECHNOLOGY")
                                                        and (c != "REGION"))]).transform('sum')
                
                prodfracsto = prodfracsto.dropna()
                prodfrac = prodfrac.dropna()

                
                # calculate ec cost (not tech/mode specific)
                # for storage fuel, average cost across timeslices is 
                # calculated
                fsto = results[i]["StorageTagFuel"].index.get_level_values("FUEL").unique().to_list()
                
                #prodfrac = prodfrac.reorder_levels(order=list(ecoutcost.index.names))
                eccost = ((prodfrac.loc[~prodfrac.index.get_level_values("FUEL").str.endswith(tuple(fsto))]
                           * ecoutcost).dropna()
                          .groupby(["RUN","REGION",
                                    "YEAR","TIMESLICE","FUEL"]).sum())

                eccoststo = ((prodfracsto.loc[prodfracsto.index.get_level_values("FUEL").str.endswith(tuple(fsto))]
                              * ecoutcost).dropna()
                          .groupby(["RUN","REGION",
                                    "YEAR","FUEL"]).sum())

                eccoststo = pd.concat(len(results[i]["TIMESLICE"]["VALUE"].to_list())*[eccoststo],
                                      keys=results[i]["TIMESLICE"]["VALUE"].to_list(),
                                      names=["TIMESLICE"])

                eccoststo = eccoststo.reorder_levels(order=list(eccost.index.names))

                eccost = pd.concat([eccost,eccoststo])

                # rearrange ec demand and cost to accommodate for cross-region
                # input activity ratios
                ecidem = ecidem.reset_index()
                ecidem.loc[~ecidem["FUEL"].str.contains(":"),"FUEL"] = ecidem.loc[~ecidem["FUEL"].str.contains(":"),"REGION"]+":"+ecidem.loc[~ecidem["FUEL"].str.contains(":"),"FUEL"]
                ecidem = ecidem.set_index([c for c in ecidem.columns if c !="VALUE"])
                eccost = eccost.droplevel("REGION")

                # calculate ec input cost
                ecicost= (ecidem * eccost).droplevel("FUEL",axis=0)
                ecicost = ecicost.dropna()
                # summing up if more than one input fuel per mode
                ecicost = ecicost.groupby(ecicost.index.names).sum()
                ecicost.index = ecicost.index.set_names("FUEL",level="FUEL_OUT")

                # calculate overall ec output cost (per tech & mode)
                # ecoutcost = varcost.reset_index().merge(right=ecicost.reset_index(),
                #                                 on=[c for c in ecicost.index.names],
                #                                 how="outer")

                ecoutcost = varcost.add(ecicost, fill_value=0)
                ecoutcost = ecoutcost.add(fixcost, fill_value=0)
                ecoutcost = ecoutcost.add(capcost, fill_value=0)
                #ecoutcost = capcost.add(ecicost, fill_value=0)

                # normalize to output of 1
                ecoutcost = ecoutcost / ol
            
                # break if no changes are observed, i.e., enough steps are done
                # to capture the entire ec tech chain
                if cc.equals(ecoutcost):
                    break
            else:
                logger.warning("Encountered an issue while calculating"
                               " energy carrier cost. Exited after"
                               " calculation 20 steps without finishing.")
                
            
            
            # rearrange ec cost
            eccost = eccost.copy().reset_index()
            eccost["REGION"] = eccost["FUEL"].str.split(":",expand=True)[0]
            eccost["FUEL"] = eccost["FUEL"].str.split(":",expand=True)[1]            
            eccost = eccost.set_index([c for c in eccost.columns if c!="VALUE"])
            results[i]["EnergyCarrierCost"] = eccost
            
            # rearrange ecout cost
            ecocost = ecoutcost.copy().reset_index()
            ecocost["REGION"] = ecocost["FUEL"].str.split(":",expand=True)[0]
            print(ecocost["FUEL"])
            # ecocost["FUEL"] = ecocost["FUEL"].str.split(":",expand=True)[1]            
            ecocost = ecocost.set_index([c for c in ecocost.columns if c!="VALUE"])
            results[i]["EnergyCarrierOutCost"] = ecocost
            
            # calculate cost for meeting demand
            dem = (results[i]["SpecifiedAnnualDemand"]
                   *results[i]["SpecifiedDemandProfile"]).dropna()
            eccost.index = eccost.index.reorder_levels(dem.index.names)
            results[i]["DemandCost"] = (dem*eccost).dropna()


            # calculate cost per demand
            results[i]["CostPerDemand"] = (results[i]["DemandCost"].groupby(["RUN","REGION","FUEL","YEAR"]).sum()
                                           /dem.groupby(["RUN","REGION","FUEL","YEAR"]).sum())


            # save unaccounted costs
            results[i]["UnaccountedCapitalCost"] = capcost_ua
            results[i]["UnaccountedFixedCost"] = fixcost_ua
            
            # FIXME: incorporate more detailed check that calculates 
            # any unaccounted cost that results from fuels being produced but 
            # not used, e.g., due to (potentially wrong) constraints -> simply 
            # take production of each fuel minus use+demand and multiply the 
            # rest times eccoutcost (?)
            
            # check if cost totals match
            if ~((results[i]["DemandCost"].sum()
                +results[i]["UnaccountedCapitalCost"].sum()
                +results[i]["UnaccountedFixedCost"].sum()).round(3)
                == results[i]["CostTotal"].sum().round(3)).squeeze():
                logger.warning("The total cost of the energy carrier cost"
                               " calculation and normal cost calculation do"
                               " not match. This is normal behaviour if the"
                               " cost of residual capacity is included, or"
                               " undiscounted costs are used.")
                logger.warning(str(results[i]["DemandCost"].sum()))
                logger.warning(str(results[i]["UnaccountedCapitalCost"].sum()))
                logger.warning(str(results[i]["UnaccountedFixedCost"].sum()))
                logger.warning(str(results[i]["CostTotal"].sum()))

            # FIXME: create option to choose marginal/shadow price
            # eccost = -results[i]["EnergyBalanceEachTS5"]   
            # calculate energy carrier cost per demand for households
            # (not including cost for building heating technologies itself)
            costpd = dict()
            for t,d,n in zip(technologies,demands,names):
                                                            
                ubt = (results[i]["UseByTechnology"]
                       [results[i]["UseByTechnology"]
                        .index.get_level_values("TECHNOLOGY").isin(t)])                                                                 
                cost= ubt * eccost
                cost = cost[cost["VALUE"]!=0].dropna()
                cost = cost.groupby(["RUN","REGION","YEAR"]).sum()
                dem = (results[i]["SpecifiedAnnualDemand"]
                       [results[i]["SpecifiedAnnualDemand"]
                        .index.get_level_values("FUEL").isin(d)])
                dem = dem.groupby(["RUN","REGION","YEAR"]).sum()
                costpd[n] = cost/dem
            logger.info(str(ii)+"p")
            results[i]["FuelCostPerDemand"] = pd.concat([df for df in costpd.values()],
                                                      keys=[n for n in costpd.keys()],
                                                      names=(["CATEGORY"]
                                                             +[c for c
                                                               in costpd[[n for n in costpd.keys()][0]].index.names]))
            
       
            
        

            results[i].pop("ProductionByTechnologyByMode")
            results[i].pop("ProductionByTechnology")
            results[i].pop("UseByTechnology")
            model.results = results
            #model.save_results("../results/study2/processedtemp/")
        #results[i].pop("CapacityFactor")

    model.results = results
    
    logger.info("Processed results")
    
    return model


if __name__ == "__main__":
    
    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("process_results",
                                          **df_params.iloc[0].to_dict())
        
    # create model and load results
    model = ft.Model()
    model.load_results(snakemake.params.odir
                       +snakemake.params.dic["name"]
                       +"_"
                       +snakemake.params.dic["run_app"]
                       +".zip")
        
    # process results for analysis
    years_map = pd.read_csv(snakemake.input.path_years_map,index_col=0)
    xscale=1/years_map.value_counts().to_frame()
    xscale.index = [y[0] for y in xscale.index]
    xscale.index.names=["YEAR"]
    xscale.columns=["VALUE"]
    
    
    if snakemake.params.dic["run_exp_res"] is not False:
        
        # FIXME: make this flexible based on snakemake parameters
        if snakemake.params.dic["scen_hh_disagg"]=="T":
            
            names=[p+te for p in ["DDDE","DDSD","DDTE","DDFL"]
                   for te in ["OO","RP","RS"]]
            techs = [[t+l  for t in ["ASHP","GSHP","HIUM","NGBO",
                                     "OIBO","ELRE","ELST","BMBO"]]
                     for l in names]
        elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
            names=[p+te for p in ["DDDE","DDSD","DDTE","DDFL"]
                   for te in [t+str(i) for t in ["O","P","S"]
                    for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)]]
            techs = [[t+l  for t in ["ASHP","GSHP","HIUM","NGBO",
                                     "OIBO","ELRE","ELST","BMBO"]]
                     for l in names]
        elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
            names=[p for p in ["DDDE","DDSD","DDTE","DDFL"]]
            techs = [[t+l+"00"  for t in ["ASHP","GSHP","HIUM","NGBO",
                                     "OIBO","ELRE","ELST","BMBO"]]
                     for l in names]
        else:
            raise NotImplementedError(
                    f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                    " not a valid option of the 'scen_hh_disagg' parameter."
                    " Valid options are 'T' or '-'."
                )        
       
        dems = [[d+l  for d in ["SH","HW"]]
                 for l in names]
        
        process_results(params=snakemake.params.dic["run_exp_res"].split("|"),
                        model=model,xscale=xscale,years_map=years_map.squeeze(),rescapacity=True,
                        discounted=False,
                        technologies=techs,demands=dems,names=names)
   
    # save results
    model.save_results(snakemake.params.odir)

    # create milestone file
    pd.DataFrame().to_csv(snakemake.output.path_resprocessed_ms)
