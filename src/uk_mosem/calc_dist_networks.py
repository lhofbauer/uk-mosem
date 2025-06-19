"""
Script deriving techno-economic parameters for distribution networks


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

import pandas as pd
import numpy as np

import utils


if __name__ == "__main__":
    
    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_dist_networks",
                                         **df_params.iloc[0].to_dict())
  

#%%
# cost
    
    # calculate district heat, electricity and gas distribution network
    # capital cost

    # load raw data set and cost parameters
    raw_data = pd.read_csv(snakemake.input.path_set_dist_net)
    eff_data = pd.read_csv(snakemake.input.path_set_eff)
    
    DHncost = raw_data[["TECHNOLOGY","UNIT",
                         "VALUE"]][raw_data["VARIABLE"]=="DHNetworkCost"]    
    DHncost = utils.adjust_monetary_values(DHncost,2015,"UNIT")
    DHncost = DHncost.drop("UNIT",axis=1)
    
    DHccost = raw_data[["TECHNOLOGY","UNIT",
                         "VALUE"]][raw_data["VARIABLE"]=="DHNetworkConnectionCost"]    
    DHccost = utils.adjust_monetary_values(DHccost,2015,"UNIT")
    DHccost = DHccost.drop("UNIT",axis=1)
    DHccost = DHccost.rename(columns={"TECHNOLOGY":"PROPERTY_TYPE"}).set_index("PROPERTY_TYPE")
    
    omfrac = raw_data.loc[raw_data["VARIABLE"]=="O&MCostFraction",
                          "VALUE"].squeeze()
    NGnlen = raw_data[["TECHNOLOGY","UNIT",
                         "VALUE"]][raw_data["VARIABLE"]=="NGResidualNetworkLength"]
    
    # explode technologies if household disaggregated
    if snakemake.params.dic["scen_hh_disagg"]=="T":
        
        DHccost = utils.explode_hh(DHccost,
                                  ["OO","RP","RS"],
                                  col="PROPERTY_TYPE")
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        
        DHccost = utils.explode_hh(DHccost,
                                  [t+str(i) for t in ["O","P","S"]
                                   for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)],
                                  col="PROPERTY_TYPE")
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )     
    DH_HIU_eff = eff_data[["TECHNOLOGY","UNIT",
                         "VALUE"]][(eff_data["VARIABLE"]=="Efficiency") &
                                   (eff_data["TECHNOLOGY"]=="HIUMDNDO00")&
                                   (eff_data["MODE_OF_OPERATION"]==1)]
                                   
    GAncost = raw_data[["TECHNOLOGY","UNIT",
                         "VALUE"]][raw_data["VARIABLE"]=="GANetworkCost"]    
    GAncost = utils.adjust_monetary_values(GAncost,2015,"UNIT")
    GAncost = GAncost.drop("UNIT",axis=1)
    GA_boiler_eff = eff_data[["TECHNOLOGY","FUEL_IN","UNIT",
                         "VALUE"]][(eff_data["VARIABLE"]=="Efficiency")&
                                   (eff_data["FUEL_IN"]=="NGD000")&
                                   (eff_data["MODE_OF_OPERATION"]==2)&
                                   (eff_data["TECHNOLOGY"]=="H2BODNDO00")]
    H2blim = raw_data[["TECHNOLOGY","UNIT",
                       "VALUE"]][raw_data["VARIABLE"]=="H2BlendLimit"]                                
    H2ncost = raw_data[["TECHNOLOGY","UNIT",
                         "VALUE"]][raw_data["VARIABLE"]=="H2UpgradeCost"]    
    H2ncost = utils.adjust_monetary_values(H2ncost,2015,"UNIT")
    H2ncost = H2ncost.drop("UNIT",axis=1)
    H2_boiler_eff = eff_data[["TECHNOLOGY","FUEL_IN","UNIT",
                         "VALUE"]][(eff_data["VARIABLE"]=="Efficiency")&
                                   (eff_data["FUEL_IN"]=="H2D000")&
                                   (eff_data["MODE_OF_OPERATION"]==1)&
                                   (eff_data["TECHNOLOGY"]=="H2BODNDO00")]
                                   
    ELncost = raw_data[["TECHNOLOGY","UNIT",
                         "VALUE"]][raw_data["VARIABLE"]=="ELNetworkCost"]    
    ELncost = utils.adjust_monetary_values(ELncost,2015,"UNIT")
    ELncost = ELncost.drop("UNIT",axis=1)
    ELLVfraction = raw_data[["TECHNOLOGY","UNIT",
                         "VALUE"]][raw_data["VARIABLE"]=="ELNetworkCostLVFraction"]
    
    
    
    # load annual demands per LSOA (adding up space and water),
    # demand fraction and length of typical periods

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
    

    typperiods_sh = pd.read_csv(snakemake.input.path_tperiods_sh,
                               index_col=["LAD23CD","PROPERTY_TYPE"])
    typperiods_hw = pd.read_csv(snakemake.input.path_tperiods_hw,
                               index_col=["LAD23CD","PROPERTY_TYPE"])
    
    typperiods_nhg = pd.read_csv(snakemake.input.path_tperiods_nhg,
                                 index_col=["LAD23CD","PROPERTY_TYPE"])
    typperiods_nhe = pd.read_csv(snakemake.input.path_tperiods_nhe,
                                 index_col=["LAD23CD","PROPERTY_TYPE"])
    
    typperiods_sh.columns.name = "TS"
    typperiods_hw.columns.name = "TS"
    typperiods_nhg.columns.name = "TS"
    typperiods_nhe.columns.name = "TS"
    typperiods_length = pd.read_csv(snakemake.input.path_tperiods_len,
                            index_col=0)
    typperiods_length.index.name = "TS"
    
 
    # load sublocal aggregation
    slagg = pd.read_csv(snakemake.input.path_sublocal_agg,
                            index_col=["LSOA11CD"]).squeeze()
    
    # no aggregation only used to move time periods to LSOA resolution for
    # space heat and hot water
    slagg_na = slagg.copy()
    slagg_na.iloc[:] = slagg_na.index.values


    

    # calculate peak gas and power demand (power) for each LSOA
    # within model (i.e., after diversity) (unit: GW)

    # anndem_sh_agg = utils.groupby_LAD(anndem_sh,disagg=slagg).sum()
    # anndem_hw_agg = utils.groupby_LAD(anndem_hw,disagg=slagg).sum()
    anndem_nhg_agg = utils.groupby_LAD(anndem_nhg,disagg=slagg).sum()
    anndem_nhe_agg = utils.groupby_LAD(anndem_nhe,disagg=slagg).sum()

    tp_sh = typperiods_sh / typperiods_length.T.values
    tp_sh = tp_sh.stack("TS")
    tp_sh = utils.groupby_LAD(tp_sh,disagg=slagg_na).mean().squeeze()
    tp_sh = tp_sh.reset_index()
    tp_sh["LSOA11CD"] = tp_sh[slagg_na.name].str.split("|",expand=True)[1]
    tp_sh = tp_sh.drop(slagg_na.name,axis=1)
    tp_sh = tp_sh.set_index(["LSOA11CD","PROPERTY_TYPE","TS"]).squeeze()

    tp_hw = typperiods_hw / typperiods_length.T.values
    tp_hw = tp_hw.stack("TS")
    tp_hw = utils.groupby_LAD(tp_hw,disagg=slagg_na).mean().squeeze()
    tp_hw = tp_hw.reset_index()
    tp_hw["LSOA11CD"] = tp_hw[slagg_na.name].str.split("|",expand=True)[1]
    tp_hw = tp_hw.drop(slagg_na.name,axis=1)
    tp_hw = tp_hw.set_index(["LSOA11CD","PROPERTY_TYPE","TS"]).squeeze()
    
    tp_nhg = typperiods_nhg / typperiods_length.T.values
    tp_nhg = tp_nhg.stack("TS")
    tp_nhg = utils.groupby_LAD(tp_nhg,disagg=slagg).mean().squeeze()
    tp_nhe = typperiods_nhe / typperiods_length.T.values
    tp_nhe = tp_nhe.stack("TS")
    tp_nhe = utils.groupby_LAD(tp_nhe,disagg=slagg).mean().squeeze()
    
    # calculate HW and SH demand profile for all LSOAs
    anndem_sh_ts = anndem_sh.multiply(pd.concat([tp_sh.to_frame(name=str(i))
                                                    for i in range(2015,2061)],
                                                  axis=1),axis=0)
    anndem_sh_ts = (anndem_sh_ts.groupby(level=["LSOA11CD","PROPERTY_TYPE",
                                                "TS"]).sum()
                          /(8760*60*60)*10**6)
    
    anndem_hw_ts = anndem_hw.multiply(
                                        pd.concat([tp_hw.to_frame(name=str(i))
                                                    for i in range(2015,2061)],
                                                  axis=1),axis=0)
    anndem_hw_ts = (anndem_hw_ts.groupby(level=["LSOA11CD","PROPERTY_TYPE",
                                                "TS"]).sum()
                          /(8760*60*60)*10**6)
    
    
    # calculate NHE and NHG demand profile for aggregate LSOAs
    anndem_nhg_agg_ts = anndem_nhg_agg.multiply(
                                        pd.concat([tp_nhg.to_frame(name=str(i))
                                                   for i in range(2015,2061)],
                                                  axis=1),axis=0)    
    anndem_nhg_agg_ts = (anndem_nhg_agg_ts.groupby(level=[slagg.name,
                                                               "TS"]).sum()
                            /(8760*60*60)*10**6)

    anndem_nhe_agg_ts = anndem_nhe_agg.multiply(
                                        pd.concat([tp_nhe.to_frame(name=str(i))
                                                   for i in range(2015,2061)],
                                                  axis=1),axis=0)
    anndem_nhe_agg_ts = (anndem_nhe_agg_ts.groupby(level=[slagg.name,
                                                               "TS"]).sum()
                            /(8760*60*60)*10**6)
    
    # load necessary efficiency data and other data
    eff_data = pd.read_csv(snakemake.input.path_set_eff)
    raw_data = pd.read_csv(snakemake.input.path_set_dist_net)
                        
                                   
    EL_heating_eff = eff_data[["TECHNOLOGY","UNIT",
                         "VALUE"]][(eff_data["VARIABLE"]=="Efficiency")&
                                   (eff_data["MODE_OF_OPERATION"]==1)&
                                   (eff_data["TECHNOLOGY"]=="ELREDNDO00")]
    DH_HIU_eff = eff_data[["TECHNOLOGY","UNIT",
                         "VALUE"]][(eff_data["VARIABLE"]=="Efficiency") &
                                   (eff_data["TECHNOLOGY"]=="HIUMDNDO00")&
                                   (eff_data["MODE_OF_OPERATION"]==1)]

    sets = pd.read_csv(snakemake.input.path_set_sets)
    years = sets["VALUE"][sets["SET"]=="YEAR"]
    
 
    # get efficiencies for gas boilers
    pt = DHccost.index.get_level_values("PROPERTY_TYPE").dropna().unique()
    NGBO_eff = pd.DataFrame([[0]*2]*len(pt),index=pt,columns=["SH","HW"])
    NGBO_eff.index.name="PROPERTY_TYPE"
    NGBO_eff["SH"] = eff_data.loc[(eff_data["TECHNOLOGY"]
                                         ==("NGBODDDE00,NGBODDSD00,"+
                                            "NGBODDTE00,NGBODDFL00"))&
                                        (eff_data["MODE_OF_OPERATION"]==1),
                                        "VALUE"].squeeze()
    
    NGBO_eff.loc["Non-domestic","SH"] = eff_data.loc[(eff_data["TECHNOLOGY"]
                                              ==("NGBODNDO00"))&
                                             (eff_data["MODE_OF_OPERATION"]==1),
                                             "VALUE"].squeeze()
    NGBO_eff["HW"] = eff_data.loc[(eff_data["TECHNOLOGY"]
                                         ==("NGBODDDE00,NGBODDSD00,"+
                                            "NGBODDTE00,NGBODDFL00"))&
                                        (eff_data["MODE_OF_OPERATION"]==3),
                                        "VALUE"].squeeze()
    
    NGBO_eff.loc["Non-domestic","HW"] = eff_data.loc[(eff_data["TECHNOLOGY"]
                                              ==("NGBODNDO00"))&
                                             (eff_data["MODE_OF_OPERATION"]==3),
                                             "VALUE"].squeeze()
    
    # calculate demand for EL
    anndem_EL_sh_ts = (anndem_sh_ts
                       /1).dropna()
    anndem_EL_hw_ts = (anndem_hw_ts
                       /1).dropna()
    anndem_EL_ts = anndem_EL_hw_ts + anndem_EL_sh_ts
    anndem_EL_ts_agg_pt = utils.groupby_LAD(anndem_EL_ts,disagg=slagg).sum()
    anndem_EL_ts_agg = anndem_EL_ts_agg_pt.groupby([slagg.name,"TS"]).sum()
    
    # calculate EL cap (GW)
    demand_peak_el = (anndem_nhe_agg_ts
                      +anndem_EL_ts_agg).groupby([slagg.name]).max().max(axis=1)/10**6
    demand_peak_el_rel = demand_peak_el/demand_peak_el.sum()
    
    # calculate demand for DH
    anndem_DH_sh_ts = (anndem_sh_ts
                       /DH_HIU_eff["VALUE"].squeeze()).dropna()
    anndem_DH_hw_ts = (anndem_hw_ts
                       /DH_HIU_eff["VALUE"].squeeze()).dropna()
    anndem_DH_ts = anndem_DH_hw_ts + anndem_DH_sh_ts
    anndem_DH_ts_agg_pt = utils.groupby_LAD(anndem_DH_ts,disagg=slagg).sum()
    anndem_DH_ts_agg = anndem_DH_ts_agg_pt.groupby([slagg.name,"TS"]).sum()
    
    # calculate DH cap (GW)
    demand_peak_dh = anndem_DH_ts_agg.groupby([slagg.name]).max().max(axis=1)/10**6
    
    # calculate demand taking into account efficiency of boilers
    anndem_NG_sh_ts = anndem_sh_ts.multiply(1/NGBO_eff["SH"],
                                            level="PROPERTY_TYPE",
                                            axis=0).dropna()
    anndem_NG_hw_ts = anndem_hw_ts.multiply(1/NGBO_eff["HW"],
                                            level="PROPERTY_TYPE",
                                            axis=0).dropna()
    anndem_NG_ts = anndem_NG_hw_ts + anndem_NG_sh_ts
    anndem_NG_ts_agg = utils.groupby_LAD(anndem_NG_ts,disagg=slagg).sum()
    anndem_NG_ts_agg = anndem_NG_ts_agg.groupby([slagg.name,"TS"]).sum()
    
    # calculate NG capacity 
    # convert to GW
    
    demand_peak_ga = (anndem_nhg_agg_ts
                     + anndem_NG_ts_agg).groupby([slagg.name]).max().max(axis=1)/10**6
    
    
    
    # calculate network lengths
    # load road lengths
    road_lengths = pd.read_csv(snakemake.input.path_road_lengths,
                            index_col=["LSOA11CD"],usecols=["LSOA11CD",
                                                            "length"])
    
    # aggregate road lengths per LAD and sublocal class
    rlengths = utils.groupby_LAD(road_lengths,disagg=slagg).sum()
    
    
    # calculate (potential) network length for gas
    
    # load heating tech residual fractions for DH residual caps
    hsys = pd.read_csv(snakemake.input.path_res_caps_ht_f, index_col=["LSOA11CD",
                                                           "PROPERTY_TYPE",
                                                           "TECHNOLOGY"])
    # calculate demand taking into account residual capacity of boilers and
    # efficiency
    anndem_NG_sh_ts = anndem_sh_ts.multiply(hsys.xs("NGBO",
                                            level=2)["fraction"]
                                            /NGBO_eff["SH"]
                                            ,axis=0).dropna()
    anndem_NG_hw_ts = anndem_hw_ts.multiply(hsys.xs("NGBO",
                                            level=2)["fraction"]
                                            /NGBO_eff["HW"]
                                            ,axis=0).dropna()
    anndem_NG_ts = anndem_NG_hw_ts + anndem_NG_sh_ts
    anndem_NG_ts_agg = utils.groupby_LAD(anndem_NG_ts,disagg=slagg).sum()
    anndem_NG_ts_agg = anndem_NG_ts_agg.groupby([slagg.name,"TS"]).sum()
    
    # calculate residual capacity as fraction of peak demand in on-grid areas,
    # convert to GW

    # offgas = utils.groupby_LAD(offgas,disagg=slagg).sum()
    # offgas.columns = pd.MultiIndex.from_tuples(offgas.columns,names=["var","YEAR"])
    # offgas = offgas["noff"]/offgas["ntotal"]
    
    GA_rescap_agg = (anndem_nhg_agg_ts
                     + anndem_NG_ts_agg).groupby([slagg.name]).max().max(axis=1)/10**6
    
    ga_frac = GA_rescap_agg.divide(demand_peak_ga,axis=0)
    
    gnlengths = rlengths * NGnlen["VALUE"].squeeze()/((ga_frac.multiply(rlengths["length"],
                                                                       axis=0)).sum()/1000)

    
    # calculate connection cost
    dw_stock = (pd.read_csv(snakemake.input.path_dw_stock_ex,
                           index_col=["LSOA11CD","PROPERTY_TYPE"])
                +pd.read_csv(snakemake.input.path_dw_stock_nb,
                                       index_col=["LSOA11CD","PROPERTY_TYPE"]))
    nd_stock = (pd.read_csv(snakemake.input.path_nd_stock_ex,
                           index_col=["LSOA11CD","PROPERTY_TYPE"])
                +pd.read_csv(snakemake.input.path_nd_stock_nb,
                                       index_col=["LSOA11CD","PROPERTY_TYPE"])).reset_index()    
    nd_stock["PROPERTY_TYPE"]="Non-domestic"
    nd_stock = nd_stock.groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
    stock = pd.concat([dw_stock,nd_stock])
    
    # assuming average number of properties/floor space over modelling period
    stock = stock.mean(axis=1)
    stock = stock.multiply(DHccost["VALUE"],axis=0)
    stock = stock.groupby("LSOA11CD").sum()
    # aggregate connection/internal pipework cost
    stock = utils.groupby_LAD(stock,disagg=slagg).sum()
    
    # calculate DH cost relative to demand (capacity) for the model
    # unit: 10^6 £/GW (£/kW)
    
    LADDHcost = (rlengths * DHncost["VALUE"].squeeze()).add(stock[0],axis=0)/10**6
    
    # FIXME: Here and below only the non-domestic tech efficiency is taken into
    # account, this could be updated in future
    LADDHcost_cap = LADDHcost.divide(demand_peak_dh,
                                    axis=0)
    LADDHcost_cap = LADDHcost_cap.replace([np.inf, -np.inf],
                                          np.nan).dropna(subset=["length"],
                                                          how="all")
    # calculate gas grid cost relative to demand (capacity) for the model
    # unit: 10^6 £/GW (£/kW)
    LADGAcost = gnlengths * GAncost["VALUE"].squeeze()/10**6
    
    LADGAcost_cap = LADGAcost.divide(demand_peak_ga, axis=0)
    LADGAcost_cap = LADGAcost_cap.replace([np.inf, -np.inf],
                                          np.nan).dropna(subset=["length"],
                                                         how="all")    

    # calculate H2 upgrade cost relative to demand (capacity) for the model,
    # divided through blend limit factor to adjust with respect to residual 
    # capacity that is implemented to allow for blending up to a certain limit 
    # without added cost
    # unit: 10^6 £/GW (£/kW)
    
    H2lim_fac = 1/(GA_boiler_eff["VALUE"].squeeze()
                   /H2_boiler_eff["VALUE"].squeeze()/H2blim["VALUE"].squeeze()
                   -(GA_boiler_eff["VALUE"].squeeze()
                     /H2_boiler_eff["VALUE"].squeeze()-1))
    
    LADH2cost = (gnlengths * H2ncost["VALUE"].squeeze() 
                 / (1-H2lim_fac))/10**6
    
    LADH2cost_cap = LADH2cost.divide(demand_peak_ga,
                                    axis=0)
    LADH2cost_cap = LADH2cost_cap.replace([np.inf, -np.inf],
                                          np.nan).dropna(subset=["length"],
                                                         how="all")
                                                         
                                                         
    # FIXME: this is currently not used and would need some reworking to
    # ensure the adjusted costs are sensible
                                           
    # calculate low voltage locally driven distribution network cost factor
    # used to differentiate cost between, e.g., rural and urban areas, by
    # multiply local grid cost factor with fraction of road length divided by 
    # fraction of peak demand (both with respect to entire modelled region)
    
    # ELLVncost = (ELncost["VALUE"].squeeze()*ELLVfraction["VALUE"].squeeze()
    #              / rlengths.sum().squeeze())
    
    # calculate reinforcement cost per capacity for the model
    # unit: 10^6 £/GW (£/kW)
    # LADELcost_cap = (ELncost["VALUE"].squeeze()
    #                  *(1-ELLVfraction["VALUE"].squeeze())
    #                  +rlengths.divide(demand_peak_el_rel,axis=0)*ELLVncost)
    LADELcost_cap = ELncost["VALUE"].squeeze()*rlengths/rlengths
                        
    
    LADELcost_cap = LADELcost_cap.replace([np.inf, -np.inf],
                                          np.nan).dropna(subset=["length"],
                                                         how="all")                                                             
    
                                                         
                                                         
    # rearrange dataframes
    LADDHcost_cap = LADDHcost_cap.reset_index()
    LADDHcost_cap["YEAR"] = ":*"
    LADDHcost_cap["MODE_OF_OPERATION"] = ""
    LADDHcost_cap["UNIT"] = "2015£/kW"
    LADDHcost_cap["VARIABLE"] = "CapitalCost"
    LADDHcost_cap["TECHNOLOGY"] = DHncost["TECHNOLOGY"].squeeze()
    LADDHcost_cap = LADDHcost_cap.rename(columns={slagg.name:"REGION",
                                                  "length":"VALUE"})
    LADDHcost_cap = LADDHcost_cap.set_index([col for col 
                                             in LADDHcost_cap.columns
                                             if col!="VALUE"])
    # create dataframe for fixed costs
    LADDHcost_fix = LADDHcost_cap.copy().reset_index()
    LADDHcost_fix["VALUE"] = LADDHcost_fix["VALUE"] * omfrac
    LADDHcost_fix["VARIABLE"] = "FixedCost"
    LADDHcost_fix = LADDHcost_fix.set_index([col for col 
                                             in LADDHcost_fix.columns
                                             if col!="VALUE"])
    

    LADELcost_cap = LADELcost_cap.reset_index()
    LADELcost_cap["YEAR"] = ":*"
    LADELcost_cap["MODE_OF_OPERATION"] = ""
    LADELcost_cap["UNIT"] = "2015£/kW"
    LADELcost_cap["VARIABLE"] = "CapitalCost"
    LADELcost_cap["TECHNOLOGY"] = ELncost["TECHNOLOGY"].squeeze()
    LADELcost_cap = LADELcost_cap.rename(columns={slagg.name:"REGION",
                                                  "length":"VALUE"})
    LADELcost_cap = LADELcost_cap.set_index([col for col 
                                             in LADELcost_cap.columns
                                             if col!="VALUE"])
    
    # create dataframe for fixed costs
    LADELcost_fix = LADELcost_cap.copy().reset_index()
    LADELcost_fix["VALUE"] = LADELcost_fix["VALUE"] * omfrac
    LADELcost_fix["VARIABLE"] = "FixedCost"
    LADELcost_fix = LADELcost_fix.set_index([col for col 
                                             in LADELcost_fix.columns
                                             if col!="VALUE"])

    LADGAcost_cap = LADGAcost_cap.reset_index()
    LADGAcost_cap["YEAR"] = ":*"
    LADGAcost_cap["MODE_OF_OPERATION"] = ""
    LADGAcost_cap["UNIT"] = "2015£/kW"
    LADGAcost_cap["VARIABLE"] = "CapitalCost"
    LADGAcost_cap["TECHNOLOGY"] = GAncost["TECHNOLOGY"].squeeze()
    LADGAcost_cap = LADGAcost_cap.rename(columns={slagg.name:"REGION",
                                                  "length":"VALUE"})
    LADGAcost_cap = LADGAcost_cap.set_index([col for col 
                                             in LADGAcost_cap.columns
                                             if col!="VALUE"])
    
    # create dataframe for fixed costs
    LADGAcost_fix = LADGAcost_cap.copy().reset_index()
    LADGAcost_fix["VALUE"] = LADGAcost_fix["VALUE"] * omfrac
    LADGAcost_fix["VARIABLE"] = "FixedCost"
    LADGAcost_fix = LADGAcost_fix.set_index([col for col 
                                             in LADGAcost_fix.columns
                                             if col!="VALUE"])
    
    LADH2cost_cap = LADH2cost_cap.reset_index()
    LADH2cost_cap["YEAR"] = ":*"
    LADH2cost_cap["MODE_OF_OPERATION"] = ""
    LADH2cost_cap["UNIT"] = "2015£/kW"
    LADH2cost_cap["VARIABLE"] = "CapitalCost"
    LADH2cost_cap["TECHNOLOGY"] = H2ncost["TECHNOLOGY"].squeeze()
    LADH2cost_cap = LADH2cost_cap.rename(columns={slagg.name:"REGION",
                                                  "length":"VALUE"})
    LADH2cost_cap = LADH2cost_cap.set_index([col for col 
                                             in LADH2cost_cap.columns
                                             if col!="VALUE"])    
    # save to file
    # unit: 10^6 £/GW (£/kW)
    pd.concat([LADELcost_cap,
               LADDHcost_cap,
               LADGAcost_cap,
               LADH2cost_cap,
               LADELcost_fix,
               LADDHcost_fix,
               LADGAcost_fix]).to_csv(snakemake.output.path_local_cost)
    
        
    
#%%
# efficiencies
    # calculate district heating, gas, and electricity grid efficiency/losses
    # electricity and gas losses currently fixed and not calculated here

    # load loss factor

    raw_data = pd.read_csv(snakemake.input.path_set_dist_net)
    
    DHnlosses = raw_data[["TECHNOLOGY","FUEL_IN","FUEL_OUT","MODE_OF_OPERATION",
                          "VALUE"]][raw_data["VARIABLE"]=="HeatLossFactor"]
    # assume pipelength value is one direction only, multiply for return flow
    avgpl = raw_data["VALUE"][raw_data["VARIABLE"]=="AvgInternalPipeLength"].squeeze()*2

    nd_raw_data = pd.read_csv(snakemake.input.path_set_nd_prop)
    ND_afa = nd_raw_data["VALUE"][nd_raw_data["VARIABLE"]=="NonDomesticAverageFloorarea"].squeeze()
    eff_data = pd.read_csv(snakemake.input.path_set_eff)   
    DH_HIU_eff = eff_data[["TECHNOLOGY","UNIT",
                         "VALUE"]][(eff_data["VARIABLE"]=="Efficiency") &
                                   (eff_data["TECHNOLOGY"]=="HIUMDNDO00")&
                                   (eff_data["MODE_OF_OPERATION"]==1)]
    
    
    # load road lengths
    road_lengths = pd.read_csv(snakemake.input.path_road_lengths,
                            index_col=["LSOA11CD"],usecols=["LSOA11CD",
                                                            "length"])
    
    # load sublocal aggregation
    slagg = pd.read_csv(snakemake.input.path_sublocal_agg,
                            index_col=["LSOA11CD"]).squeeze()
    
    # aggregate road lengths per LAD and sublocal class
    rlengths = utils.groupby_LAD(road_lengths,disagg=slagg).sum()
    
    # calculate internal pipe length
    dw_stock = (pd.read_csv(snakemake.input.path_dw_stock_ex,
                           index_col=["LSOA11CD","PROPERTY_TYPE"])
                +pd.read_csv(snakemake.input.path_dw_stock_nb,
                                       index_col=["LSOA11CD","PROPERTY_TYPE"]))
    nd_stock = (pd.read_csv(snakemake.input.path_nd_stock_ex,
                           index_col=["LSOA11CD","PROPERTY_TYPE"])
                +pd.read_csv(snakemake.input.path_nd_stock_nb,
                                       index_col=["LSOA11CD","PROPERTY_TYPE"])).reset_index()    
    nd_stock["PROPERTY_TYPE"] = "Non-domestic"
    nd_stock = nd_stock.groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
    stock = pd.concat([dw_stock,nd_stock/ND_afa])

    stock = stock.groupby("LSOA11CD").sum()
    
    # aggregate connection/internal pipework length
    stock = utils.groupby_LAD(stock,disagg=slagg).sum()
    stock = stock * avgpl
    
    # calculated lossed for each LAD and class
    # multiplied times 2 to account for return flow
    LADDHlosses = stock.add(rlengths["length"]*2,axis=0) * DHnlosses["VALUE"].squeeze()*8760*60*60/10**9

    # load annual heat demands per LSOA, adding up space and water
    dw_anndem = (pd.read_csv(snakemake.input.path_dw_ann_dem_sh,
                             index_col=["LSOA11CD","PROPERTY_TYPE"])
                 + pd.read_csv(snakemake.input.path_dw_ann_dem_hw,
                               index_col=["LSOA11CD","PROPERTY_TYPE"])
                 + pd.read_csv(snakemake.input.path_dw_ann_dem_sh_nb,
                               index_col=["LSOA11CD","PROPERTY_TYPE"]))
    nd_anndem = (pd.read_csv(snakemake.input.path_nd_ann_dem_sh,
                  index_col=["LSOA11CD","PROPERTY_TYPE"])
                 + pd.read_csv(snakemake.input.path_nd_ann_dem_hw,
                               index_col=["LSOA11CD","PROPERTY_TYPE"])
                 + pd.read_csv(snakemake.input.path_nd_ann_dem_sh_nb,
                               index_col=["LSOA11CD","PROPERTY_TYPE"])).reset_index()
    nd_anndem["PROPERTY_TYPE"]="Non-domestic"
    nd_anndem = nd_anndem.groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
    anndem = pd.concat([dw_anndem,nd_anndem]).divide(DH_HIU_eff["VALUE"].squeeze())


    
    # calculate loss relative to annual demand for each year to come up with
    # efficiency
    # pumping 'loss' is taken into account through electricity IAR, see data set
    anndem_agg = utils.groupby_LAD(anndem,disagg=slagg.squeeze()).sum()
    anndem_agg = anndem_agg.groupby([slagg.name]).sum()
    LADDHlosses_dem =  1-1/anndem_agg.divide(LADDHlosses, axis=0)
    LADDHlosses_dem[LADDHlosses_dem<0] = 0
    

    
    # rearrange dataframe
    LADDHlosses_dem.columns.name="YEAR"
    LADDHlosses_dem = LADDHlosses_dem.stack()
    
    LADDHlosses_dem.name = "VALUE"
    LADDHlosses_dem = LADDHlosses_dem.reset_index()
    LADDHlosses_dem["FUEL_IN"] = DHnlosses["FUEL_IN"].squeeze()
    LADDHlosses_dem["FUEL_OUT"] = DHnlosses["FUEL_OUT"].squeeze()
    LADDHlosses_dem["MODE_OF_OPERATION"] = DHnlosses["MODE_OF_OPERATION"].squeeze()
    LADDHlosses_dem["VARIABLE"] = "Efficiency"
    LADDHlosses_dem["TECHNOLOGY"] = DHnlosses["TECHNOLOGY"].squeeze()
    LADDHlosses_dem = LADDHlosses_dem.rename(columns={slagg.name:"REGION"})
    LADDHlosses_dem = LADDHlosses_dem.set_index([col for col 
                                                 in LADDHlosses_dem.columns
                                                 if col!="VALUE"])
    # save to file
    # unit: -
    LADDHlosses_dem.to_csv(snakemake.output.path_local_eff)

    
#%% constraints

    # calculate district heating and gas network constraints
    # FIXME: this is not necessary anymore given sublocal representation
    # has been added, leaving it here commented out in case of a need to
    # reimplement this
        
    # # load necessary efficiency data
    # eff_data = pd.read_csv(snakemake.input.path_set_eff)
    # DH_HIU_eff = eff_data[["TECHNOLOGY","UNIT",
    #                      "VALUE"]][(eff_data["VARIABLE"]=="Efficiency") &
    #                                (eff_data["TECHNOLOGY"]=="HIUMDNDO00")&
    #                                (eff_data["MODE_OF_OPERATION"]==1)]
    # GA_boiler_eff = eff_data[["TECHNOLOGY","FUEL_IN","UNIT",
    #                      "VALUE"]][(eff_data["VARIABLE"]=="Efficiency")&
    #                                (eff_data["FUEL_IN"]=="NGD000")&
    #                                (eff_data["MODE_OF_OPERATION"]==2)&
    #                                (eff_data["TECHNOLOGY"]=="H2BODNDO00")]
                 
    # # load annual heat and gas demands per LSOA (currently adding up space and 
    # # water), demand fraction and length of typical periods
    
    # dw_anndem = (pd.read_csv(snakemake.input.path_dw_ann_dem_sh,
    #                          index_col=["LSOA11CD","PROPERTY_TYPE"])
    #              + pd.read_csv(snakemake.input.path_dw_ann_dem_hw,
    #                            index_col=["LSOA11CD","PROPERTY_TYPE"])
    #              + pd.read_csv(snakemake.input.path_dw_ann_dem_sh_nb,
    #                            index_col=["LSOA11CD","PROPERTY_TYPE"]))
    # nd_anndem = (pd.read_csv(snakemake.input.path_nd_ann_dem_sh,
    #               index_col=["LSOA11CD","PROPERTY_TYPE"])
    #              + pd.read_csv(snakemake.input.path_nd_ann_dem_hw,
    #                            index_col=["LSOA11CD","PROPERTY_TYPE"])
    #              + pd.read_csv(snakemake.input.path_nd_ann_dem_sh_nb,
    #                            index_col=["LSOA11CD","PROPERTY_TYPE"])).reset_index()
    # nd_anndem["PROPERTY_TYPE"]="Non-domestic"
    # nd_anndem = nd_anndem.groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
    # anndem_h = pd.concat([dw_anndem,nd_anndem])
    
    # dw_anndem = pd.read_csv(snakemake.input.path_dw_ann_dem_nhg,
    #                          index_col=["LSOA11CD","PROPERTY_TYPE"])
    # nd_anndem = pd.read_csv(snakemake.input.path_nd_ann_dem_nhg)
    # nd_anndem["PROPERTY_TYPE"]="Non-domestic"
    # nd_anndem = nd_anndem.groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
    # anndem_nhg = pd.concat([dw_anndem,nd_anndem])
    
    
    # typperiods_shw = pd.read_csv(snakemake.input.path_tperiods_shw,
    #                            index_col=["LAD21CD","PROPERTY_TYPE"])
    # typperiods_nhg = pd.read_csv(snakemake.input.path_tperiods_nhg,
    #                              index_col=["LAD21CD","PROPERTY_TYPE"])
    # typperiods_shw.columns.name = "TS"
    # typperiods_nhg.columns.name = "TS"
    # typperiods_length = pd.read_csv(snakemake.input.path_tperiods_len,
    #                         index_col=0)
    # typperiods_length.index.name = "TS"
    
 
    # # load sublocal aggregation
    # slagg = pd.read_csv(snakemake.input.path_sublocal_agg,
    #                         index_col=["LSOA11CD"])
    
    # # create separate class for non disaggregated, set all SLClass strings
    # # to c_0, i.e., currently no differentiation for electricity and gas
    # slagg_na = slagg.copy()    
    # slagg_na["SLClass"] = "c_0"
    

    
    # # calculate peak heat demand (power) for each LAD within model 
    # # (i.e., after diversity) (unit: GW)
    # anndem_h_agg = utils.groupby_LAD(anndem_h,disagg=slagg.squeeze()).sum()
    # #anndem_h_agg_max = anndem_h_agg.max(axis=1)
    
    # tp_h = typperiods_shw / typperiods_length.T.values
    # tp_h = tp_h.stack("TS")

    # demand_h_agg_ts = anndem_h_agg.reorder_levels([0,2,1]).multiply(
    #                                     pd.concat([tp_h.to_frame(name=str(i))
    #                                                for i in range(2015,2061)],
    #                                               axis=1),axis=0)
    # demand_h_agg_ts = (demand_h_agg_ts.groupby(level=["LAD21CD",
    #                                                   "SLClass",
    #                                                   "TS"]).sum()
    #                       /(8760*60*60))
    
    # demand_peak_h = demand_h_agg_ts.groupby(["LAD21CD",
    #                                          "SLClass"]).max()
    

    # # calculate peak gas demand (power) for each LAD without disagg within
    # # model (i.e., after diversity) (unit: GW)

    # anndem_h_agg_na = utils.groupby_LAD(anndem_h,disagg=slagg_na.squeeze()).sum()
    # anndem_nhg_agg_na = utils.groupby_LAD(anndem_nhg,disagg=slagg_na.squeeze()).sum()
    
    # tp_h = typperiods_shw / typperiods_length.T.values
    # tp_h = tp_h.stack("TS")
    # tp_nhg = typperiods_nhg / typperiods_length.T.values
    # tp_nhg = tp_nhg.stack("TS")   
    
    # anndem_h_agg_na_ts = anndem_h_agg_na.reorder_levels([0,2,1]).multiply(
    #                                     pd.concat([tp_h.to_frame(name=str(i))
    #                                                for i in range(2015,2061)],
    #                                               axis=1),axis=0)
    # anndem_nhg_agg_na_ts = anndem_nhg_agg_na.reorder_levels([0,2,1]).multiply(
    #                                     pd.concat([tp_nhg.to_frame(name=str(i))
    #                                                for i in range(2015,2061)],
    #                                               axis=1),axis=0)    
    # anndem_gas_agg_na_ts = (anndem_nhg_agg_na_ts
    #                      + anndem_h_agg_na_ts.divide(GA_boiler_eff["VALUE"].squeeze()))
    
    # anndem_gas_agg_max_na = anndem_gas_agg_na_ts.groupby(["LAD21CD",
    #                                          "SLClass","TS"]).sum()/(8760*60*60)
    
    # demand_peak_gas_na = anndem_gas_agg_max_na.groupby(["LAD21CD",
    #                                          "SLClass"]).max()
    
    
    # # take model period peak, apply efficiency of heat interface unit
    # # and convert to GW
    # DH_cap_constr = demand_peak_h.divide(DH_HIU_eff["VALUE"].squeeze())
    # DH_cap_constr.columns.name = "YEAR"
    # DH_cap_constr = DH_cap_constr.stack()
    # GA_cap_constr = demand_peak_gas_na
    # GA_cap_constr.columns.name = "YEAR"
    # GA_cap_constr = GA_cap_constr.stack()    
    # # rearrange dataframe
    # # FIXME: solve the issues of having technology names in the code here and
    # # in other places
    # DH_cap_constr = DH_cap_constr.reset_index()
    # DH_cap_constr["UNIT"] = "GW"
    # DH_cap_constr["VARIABLE"] = "TotalAnnualMaxCapacity"
    # DH_cap_constr["SLClass"] = ("DHMTTDIS00"[:-1]
    #                             +DH_cap_constr["SLClass"].str[-1])
    # DH_cap_constr = DH_cap_constr.rename(columns={"LAD21CD":"REGION",
    #                                               "SLClass":"TECHNOLOGY",
    #                                               0:"VALUE"})
    # DH_cap_constr = DH_cap_constr.set_index([col for col 
    #                                          in DH_cap_constr.columns
    #                                          if col!="VALUE"])
    

    
    # # rearrange dataframe
    # GA_cap_constr = GA_cap_constr.reset_index()
    # GA_cap_constr["UNIT"] = "GW"
    # GA_cap_constr["VARIABLE"] = "TotalAnnualMaxCapacity"
    # GA_cap_constr["SLClass"] = ("GAGRTDIS00"[:-1]
    #                             +GA_cap_constr["SLClass"].str[-1])
    # GA_cap_constr = GA_cap_constr.rename(columns={"LAD21CD":"REGION",
    #                                               "SLClass":"TECHNOLOGY",
    #                                               0:"VALUE"})
    # GA_cap_constr = GA_cap_constr.set_index([col for col 
    #                                          in GA_cap_constr.columns
    #                                          if col!="VALUE"])
    
    DH_cap_constr = pd.DataFrame(columns=["REGION","TECHNOLOGY","YEAR",
                                          "VARIBALE","UNIT","VALUE"])
    GA_cap_constr = pd.DataFrame(columns=["REGION","TECHNOLOGY","YEAR",
                                          "VARIBALE","UNIT","VALUE"])
    # save to file
    # unit: GW
    pd.concat([DH_cap_constr,
               GA_cap_constr]).to_csv(snakemake.output.path_local_con)
    
    
    

    
#%% 
# capacity factors
    
    # calculate district heating and gas network capacity factors
    # FIXME: NG and H2 currently not used (and not working),
    # change or delete (see notes)
    
    # calculate capacity factors per technology/sublocal area, LAD, timeslice, 
    # and year by dividing each demand by overall peak demand. This ensures
    # 1) not just the most efficient technologies run in non-peak hours, and 
    # 2) the actual capacity necessary is built (see documentation of DH)

    # load necessary efficiency data
    eff_data = pd.read_csv(snakemake.input.path_set_eff)
    GA_boiler_eff = eff_data[["TECHNOLOGY","FUEL_IN","UNIT",
                         "VALUE"]][(eff_data["VARIABLE"]=="Efficiency")&
                                   (eff_data["FUEL_IN"]=="NGD000")&
                                   (eff_data["MODE_OF_OPERATION"]==2)&
                                   (eff_data["TECHNOLOGY"]=="H2BODNDO00")]
                                   
    # load annual heat and gas demands per LSOA (currently adding up space and 
    # water), demand fraction and length of typical periods
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
    
    
    typperiods_sh = pd.read_csv(snakemake.input.path_tperiods_sh,
                               index_col=["LAD23CD","PROPERTY_TYPE"])
    typperiods_hw = pd.read_csv(snakemake.input.path_tperiods_hw,
                               index_col=["LAD23CD","PROPERTY_TYPE"])
    
    typperiods_nhg = pd.read_csv(snakemake.input.path_tperiods_nhg,
                                 index_col=["LAD23CD","PROPERTY_TYPE"])
    typperiods_sh.columns.name = "TS"
    typperiods_hw.columns.name = "TS"
    typperiods_nhg.columns.name = "TS"
    typperiods_length = pd.read_csv(snakemake.input.path_tperiods_len,
                            index_col=0)
    typperiods_length.index.name = "TS"
    
 
    # load sublocal aggregation
    slagg = pd.read_csv(snakemake.input.path_sublocal_agg,
                            index_col=["LSOA11CD"]).squeeze()

    # calculate peak heat demand (power) for each LAD within model 
    # (i.e., after diversity) (unit: GW)
    
    anndem_sh_agg = utils.groupby_LAD(anndem_sh,disagg=slagg).sum()
    anndem_hw_agg = utils.groupby_LAD(anndem_hw,disagg=slagg).sum()
    anndem_nhg_agg = utils.groupby_LAD(anndem_nhg,disagg=slagg).sum()
    
    tp_sh = typperiods_sh / typperiods_length.T.values
    tp_sh = tp_sh.stack("TS")
    tp_sh = utils.groupby_LAD(tp_sh,disagg=slagg).mean().squeeze()


    tp_hw = typperiods_hw / typperiods_length.T.values
    tp_hw = tp_hw.stack("TS")
    tp_hw = utils.groupby_LAD(tp_hw,disagg=slagg).mean().squeeze()
    
    tp_nhg = typperiods_nhg / typperiods_length.T.values
    tp_nhg = tp_nhg.stack("TS")
    tp_nhg = utils.groupby_LAD(tp_nhg,disagg=slagg).mean().squeeze()
    
    # calculate HW and SH demand profile
    anndem_sh_agg_ts = anndem_sh_agg.multiply(pd.concat([tp_sh.to_frame(name=str(i))
                                                    for i in range(2015,2061)],
                                                  axis=1),axis=0)
    anndem_sh_agg_ts = (anndem_sh_agg_ts.groupby(level=[slagg.name,
                                                "TS"]).sum()
                          /(8760*60*60)*10**6)
    
    anndem_hw_agg_ts = anndem_hw_agg.multiply(
                                        pd.concat([tp_hw.to_frame(name=str(i))
                                                    for i in range(2015,2061)],
                                                  axis=1),axis=0)
    anndem_hw_agg_ts = (anndem_hw_agg_ts.groupby(level=[slagg.name,
                                                "TS"]).sum()
                          /(8760*60*60)*10**6)
    
    anndem_shw_agg_ts = anndem_hw_agg_ts + anndem_sh_agg_ts

    demand_peak_h = anndem_shw_agg_ts.groupby([slagg.name]).max()
    
    # calculate peak gas demand (power) for each sublocal area within
    # model (i.e., after diversity) (unit: GW)

    anndem_nhg_agg_ts = anndem_nhg_agg.multiply(
                                        pd.concat([tp_nhg.to_frame(name=str(i))
                                                   for i in range(2015,2061)],
                                                  axis=1),axis=0) 
    

    anndem_gas_agg_ts = (anndem_nhg_agg_ts
                         + anndem_shw_agg_ts.divide(GA_boiler_eff["VALUE"].squeeze()))
    
    demand_peak_gas = anndem_gas_agg_ts.groupby([slagg.name]).max() 
    
    # calculate capacity factor
    
    anndem_gas_agg_ts.columns.name = "YEAR"
    demand_capacity_gas = anndem_gas_agg_ts.stack()
    # demand_capacity_gas_na = demand_capacity_gas_na.groupby(level=["LAD21CD",
    #                                                                "SLClass",
    #                                                                "YEAR",
    #                                                                "TS"]).sum()       
    DH_capfac =  (anndem_shw_agg_ts)/demand_peak_h
    DH_capfac.columns.name = "YEAR"
    DH_capfac = DH_capfac.stack()
    #DH_capfac = DH_capfac.dropna()

    GA_capfac =  anndem_gas_agg_ts/demand_peak_gas
    GA_capfac.columns.name = "YEAR"
    GA_capfac = GA_capfac.stack()
    GA_capfac = GA_capfac.dropna()
    GA_capfac = DH_capfac.to_frame()[0:0]
    H2_capfac = GA_capfac.copy()

    # rearrange dataframes
    DH_capfac = DH_capfac.reset_index()
    DH_capfac["UNIT"] = "-"
    DH_capfac["VARIABLE"] = "CapacityFactor"
    DH_capfac["TECHNOLOGY"] = "DHMTTDIS00"
    DH_capfac = DH_capfac.rename(columns={slagg.name:"REGION",
                                    "TS":"TIMESLICE",
                                    0:"VALUE"})
    DH_capfac = DH_capfac.set_index([col for col 
                               in DH_capfac.columns
                               if col!="VALUE"])

    GA_capfac = GA_capfac.reset_index()
    GA_capfac["UNIT"] = "-"
    GA_capfac["VARIABLE"] = "CapacityFactor"
    GA_capfac["TECHNOLOGY"] = "GAGRTDIS00"
    GA_capfac = GA_capfac.rename(columns={slagg.name:"REGION",
                                    "TS":"TIMESLICE",
                                    0:"VALUE"})
    GA_capfac = GA_capfac.set_index([col for col 
                               in GA_capfac.columns
                               if col!="VALUE"])
    
    H2_capfac = H2_capfac.reset_index()
    H2_capfac["UNIT"] = "-"
    H2_capfac["VARIABLE"] = "CapacityFactor"
    H2_capfac["TECHNOLOGY"] = "H2UGTDIS00"
    H2_capfac = H2_capfac.rename(columns={slagg.name:"REGION",
                                    "TS":"TIMESLICE",
                                    0:"VALUE"})
    H2_capfac = H2_capfac.set_index([col for col 
                               in H2_capfac.columns
                               if col!="VALUE"])
    
    # save to file
    # unit: -
    pd.concat([DH_capfac,
               GA_capfac,
               H2_capfac]).to_csv(snakemake.output.path_local_capf)
    
#%% 
# residual capacities
    # calculate gas and electricity distribution network residual capacity
    
    # load data on properties off gas grid per LSOA
    # og = pd.read_excel(snakemake.input.path_offgas,sheet_name=["2015","2016","2017",
    #                                                      "2018","2019","2020"],
    #                        header=[0],skiprows=4,usecols=["LSOA code",
    #                                                       "Estimated number \nof properties not \non the gas grid\n[note 3]",
    #                                                       "Number of \ndomestic \nproperties\n[note 1]"],
    #                        index_col="LSOA code")
    

    # offgas = pd.concat(og).unstack(0)
    # offgas = offgas.rename(columns={"Estimated number \nof properties not \non the gas grid\n[note 3]":"noff",
    #                                 "Number of \ndomestic \nproperties\n[note 1]":"ntotal"})
    # offgas.index.name = "LSOA11CD"

    
    # load heating tech residual fractions for DH residual caps
    hsys = pd.read_csv(snakemake.input.path_res_caps_ht_f, index_col=["LSOA11CD",
                                                           "PROPERTY_TYPE",
                                                           "TECHNOLOGY"])
        
    # load necessary efficiency data and other data
    eff_data = pd.read_csv(snakemake.input.path_set_eff)
    raw_data = pd.read_csv(snakemake.input.path_set_dist_net)
    GA_boiler_eff = eff_data[["TECHNOLOGY","FUEL_IN","UNIT",
                         "VALUE"]][(eff_data["VARIABLE"]=="Efficiency")&
                                   (eff_data["FUEL_IN"]=="NGD000")&
                                   (eff_data["MODE_OF_OPERATION"]==2)&
                                   (eff_data["TECHNOLOGY"]=="H2BODNDO00")]                           
    H2blim = raw_data[["TECHNOLOGY","UNIT",
                       "VALUE"]][raw_data["VARIABLE"]=="H2BlendLimit"]
    H2_boiler_eff = eff_data[["TECHNOLOGY","FUEL_IN","UNIT",
                         "VALUE"]][(eff_data["VARIABLE"]=="Efficiency")&
                                   (eff_data["FUEL_IN"]=="H2D000")&
                                   (eff_data["MODE_OF_OPERATION"]==1)&
                                   (eff_data["TECHNOLOGY"]=="H2BODNDO00")]
                                   
    EL_heating_eff = eff_data[["TECHNOLOGY","UNIT",
                         "VALUE"]][(eff_data["VARIABLE"]=="Efficiency")&
                                   (eff_data["MODE_OF_OPERATION"]==1)&
                                   (eff_data["TECHNOLOGY"]=="ELREDNDO00")]
    DH_HIU_eff = eff_data[["TECHNOLOGY","UNIT",
                         "VALUE"]][(eff_data["VARIABLE"]=="Efficiency") &
                                   (eff_data["TECHNOLOGY"]=="HIUMDNDO00")&
                                   (eff_data["MODE_OF_OPERATION"]==1)]
                                   
    life_data = pd.read_csv(snakemake.input.path_set_lt)
    
    GA_rescap_red_start = life_data[["TECHNOLOGY","UNIT",
                         "VALUE"]][(life_data["VARIABLE"]=="ResCapRedStart")&
                                   (life_data["TECHNOLOGY"]=="GAGRTDIS00")]
    EL_rescap_red_start = life_data[["TECHNOLOGY","UNIT",
                         "VALUE"]][(life_data["VARIABLE"]=="ResCapRedStart")&
                                   (life_data["TECHNOLOGY"]=="ELGRTDIS00")]
                                   
    sets = pd.read_csv(snakemake.input.path_set_sets)
    years = sets["VALUE"][sets["SET"]=="YEAR"]
    

                                   
    # load annual heat and gas demands per LSOA, adding up space and 
    # water, demand fraction and length of typical periods
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
    
    # dw_anndem = (pd.read_csv(snakemake.input.path_dw_ann_dem_sh,
    #                          index_col=["LSOA11CD","PROPERTY_TYPE"])
    #              + pd.read_csv(snakemake.input.path_dw_ann_dem_hw,
    #                            index_col=["LSOA11CD","PROPERTY_TYPE"])
    #              + pd.read_csv(snakemake.input.path_dw_ann_dem_sh_nb,
    #                            index_col=["LSOA11CD","PROPERTY_TYPE"]))
    # nd_anndem = (pd.read_csv(snakemake.input.path_nd_ann_dem_sh,
    #               index_col=["LSOA11CD","PROPERTY_TYPE"])
    #              + pd.read_csv(snakemake.input.path_nd_ann_dem_hw,
    #                            index_col=["LSOA11CD","PROPERTY_TYPE"])
    #              + pd.read_csv(snakemake.input.path_nd_ann_dem_sh_nb,
    #                            index_col=["LSOA11CD","PROPERTY_TYPE"])).reset_index()
    # nd_anndem["PROPERTY_TYPE"]="Non-domestic"
    # nd_anndem = nd_anndem.groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
    # anndem_h = pd.concat([dw_anndem,nd_anndem])
    
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
    
    typperiods_hw = pd.read_csv(snakemake.input.path_tperiods_hw,
                               index_col=["LAD23CD","PROPERTY_TYPE"])
    typperiods_sh = pd.read_csv(snakemake.input.path_tperiods_sh,
                               index_col=["LAD23CD","PROPERTY_TYPE"])
    typperiods_nhg = pd.read_csv(snakemake.input.path_tperiods_nhg,
                                 index_col=["LAD23CD","PROPERTY_TYPE"])
    typperiods_nhe = pd.read_csv(snakemake.input.path_tperiods_nhe,
                                 index_col=["LAD23CD","PROPERTY_TYPE"])
    
    typperiods_hw.columns.name = "TS"
    typperiods_sh.columns.name = "TS"
    typperiods_nhg.columns.name = "TS"
    typperiods_nhe.columns.name = "TS"
    
    typperiods_length = pd.read_csv(snakemake.input.path_tperiods_len,
                            index_col=0)
    typperiods_length.index.name = "TS"
    
 
    # load sublocal aggregation
    slagg = pd.read_csv(snakemake.input.path_sublocal_agg,
                            index_col=["LSOA11CD"]).squeeze()
    
    # no aggregation only used to move time periods to LSOA resolution for
    # space heat and hot water
    slagg_na = slagg.copy()
    slagg_na.iloc[:] = slagg_na.index.values
    
    # calculate peak gas and power demand (power) for each LSOA
    # within model (i.e., after diversity) (unit: GW)

    # anndem_sh_agg = utils.groupby_LAD(anndem_sh,disagg=slagg).sum()
    # anndem_hw_agg = utils.groupby_LAD(anndem_hw,disagg=slagg).sum()
    anndem_nhg_agg = utils.groupby_LAD(anndem_nhg,disagg=slagg).sum()
    anndem_nhe_agg = utils.groupby_LAD(anndem_nhe,disagg=slagg).sum()

    tp_sh = typperiods_sh / typperiods_length.T.values
    tp_sh = tp_sh.stack("TS")
    tp_sh = utils.groupby_LAD(tp_sh,disagg=slagg_na).mean().squeeze()
    tp_sh = tp_sh.reset_index()
    tp_sh["LSOA11CD"] = tp_sh[slagg_na.name].str.split("|",expand=True)[1]
    tp_sh = tp_sh.drop(slagg_na.name,axis=1)
    tp_sh = tp_sh.set_index(["LSOA11CD","PROPERTY_TYPE","TS"]).squeeze()

    tp_hw = typperiods_hw / typperiods_length.T.values
    tp_hw = tp_hw.stack("TS")
    tp_hw = utils.groupby_LAD(tp_hw,disagg=slagg_na).mean().squeeze()
    tp_hw = tp_hw.reset_index()
    tp_hw["LSOA11CD"] = tp_hw[slagg_na.name].str.split("|",expand=True)[1]
    tp_hw = tp_hw.drop(slagg_na.name,axis=1)
    tp_hw = tp_hw.set_index(["LSOA11CD","PROPERTY_TYPE","TS"]).squeeze()
    
    tp_nhg = typperiods_nhg / typperiods_length.T.values
    tp_nhg = tp_nhg.stack("TS")
    tp_nhg = utils.groupby_LAD(tp_nhg,disagg=slagg).mean().squeeze()
    tp_nhe = typperiods_nhe / typperiods_length.T.values
    tp_nhe = tp_nhe.stack("TS")
    tp_nhe = utils.groupby_LAD(tp_nhe,disagg=slagg).mean().squeeze()
    
    # calculate HW and SH demand profile for all LSOAs
    anndem_sh_ts = anndem_sh.multiply(pd.concat([tp_sh.to_frame(name=str(i))
                                                    for i in range(2015,2061)],
                                                  axis=1),axis=0)
    anndem_sh_ts = (anndem_sh_ts.groupby(level=["LSOA11CD","PROPERTY_TYPE",
                                                "TS"]).sum()
                          /(8760*60*60)*10**6)
    
    anndem_hw_ts = anndem_hw.multiply(
                                        pd.concat([tp_hw.to_frame(name=str(i))
                                                    for i in range(2015,2061)],
                                                  axis=1),axis=0)
    anndem_hw_ts = (anndem_hw_ts.groupby(level=["LSOA11CD","PROPERTY_TYPE",
                                                "TS"]).sum()
                          /(8760*60*60)*10**6)
    
    
    # calculate NHE and NHG demand profile for aggregate LSOAs
    anndem_nhg_agg_ts = anndem_nhg_agg.multiply(
                                        pd.concat([tp_nhg.to_frame(name=str(i))
                                                   for i in range(2015,2061)],
                                                  axis=1),axis=0)    
    anndem_nhg_agg_ts = (anndem_nhg_agg_ts.groupby(level=[slagg.name,
                                                               "TS"]).sum()
                            /(8760*60*60)*10**6)

    anndem_nhe_agg_ts = anndem_nhe_agg.multiply(
                                        pd.concat([tp_nhe.to_frame(name=str(i))
                                                   for i in range(2015,2061)],
                                                  axis=1),axis=0)
    anndem_nhe_agg_ts = (anndem_nhe_agg_ts.groupby(level=[slagg.name,
                                                               "TS"]).sum()
                            /(8760*60*60)*10**6)
    
    # calculate demand for DH based on residual HIUM fraction
    anndem_DH_sh_ts = (anndem_sh_ts.multiply(hsys.xs("HIUM",
                                            level=2)["fraction"],axis=0)
                       /DH_HIU_eff["VALUE"].squeeze()).dropna()
    anndem_DH_hw_ts = (anndem_hw_ts.multiply(hsys.xs("HIUM",
                                            level=2)["fraction"],axis=0)
                       /DH_HIU_eff["VALUE"].squeeze()).dropna()
    anndem_DH_ts = anndem_DH_hw_ts + anndem_DH_sh_ts
    anndem_DH_ts_agg_pt = utils.groupby_LAD(anndem_DH_ts,disagg=slagg).sum()
    anndem_DH_ts_agg = anndem_DH_ts_agg_pt.groupby([slagg.name,"TS"]).sum()
    
    # calculate DH residual cap (GW)
    DH_rescap_agg = anndem_DH_ts_agg.groupby([slagg.name]).max()/10**6
    
    # calculate residual caps for DH heat generation 
    dhgen = pd.read_csv(snakemake.input.path_dw_stock_dh, index_col=["LSOA11CD",
                                                                     "PROPERTY_TYPE",
                                                                     "DH_TYPE"])
    #dhgen = utils.groupby_LAD(dhgen,disagg=slagg).sum()
    pt_frac = anndem_DH_ts.groupby(["LSOA11CD","PROPERTY_TYPE"]).max()
    # for the fraction, do not use non-domestic data (for capacity it is used)
    pt_frac = pt_frac.drop("Non-domestic",level=1)
    pt_frac = pt_frac/pt_frac.groupby("LSOA11CD").sum()
    dhgen = dhgen/dhgen.groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
    dhgen = dhgen["NUMBER_OF_PROPERTIES"]*pt_frac["2022"]
    dhgen = dhgen.groupby(["LSOA11CD","DH_TYPE"]).sum()
    
    dhgen = dhgen.rename(index={"HIUM-C":"NGCHSDIS00",
                                "HIUM-B":"BMBOSDIS00",
                                "HIUM-G":"NGBOSDIS00",
                                "HIUM-H":"RWHPSDIS00"})
    dhgen.index.names = ["LSOA11CD","TECHNOLOGY"]
    
    #dhgen = utils.groupby_LAD(dhgen).sum()
    
    # DH generation rescap in GW
    anndem_DH_ts_LAD = utils.groupby_LAD(anndem_DH_ts.groupby(["LSOA11CD",
                                                               "TS"]).sum()).sum()
    
    DH_gen_rescap = (anndem_DH_ts.groupby(["LSOA11CD","TS"]).sum()
                         .groupby("LSOA11CD").max()
                         .multiply(dhgen,level="LSOA11CD",axis=0)
                         /10**6)
    DH_gen_rescap_agg = utils.groupby_LAD(DH_gen_rescap).sum()
    

        

    
    # calculate residual cap for electricity distribution grid
    anndem_EL_sh_ts = (anndem_sh_ts.multiply(
                    (hsys.loc[hsys.index.get_level_values(2).str.contains("ASHP")
                                |hsys.index.get_level_values(2).str.contains("AWHP")
                                |hsys.index.get_level_values(2).str.contains("AAHP"),"fraction"].groupby(["LSOA11CD",
                                                                                                          "PROPERTY_TYPE"]).sum()
                                /eff_data.loc[(eff_data["TECHNOLOGY"]
                                              ==("ASHPDDDE00,ASHPDDSD00,"+
                                                  "ASHPDDTE00,ASHPDDFL00"))&(
                                              eff_data["YEAR"]
                                                              =="2020")&
                                            (eff_data["MODE_OF_OPERATION"]==1), "VALUE"].values
                            ).fillna(0).add(
                        (hsys.loc[hsys.index.get_level_values(2).str.contains("GSHP"),"fraction"].droplevel(2)
                              /eff_data.loc[(eff_data["TECHNOLOGY"]
                                              ==("GSHPDDDE00,GSHPDDSD00,"+
                                                "GSHPDDTE00,GSHPDDFL00"))&(
                                              eff_data["YEAR"]
                                                            =="2020")&
                                          (eff_data["MODE_OF_OPERATION"]==1), "VALUE"].values
                                  ).fillna(0),fill_value=0).add(
                          (hsys.loc[hsys.index.get_level_values(2).str.contains("ELRE"),"fraction"].droplevel(2)
                              /eff_data.loc[(eff_data["TECHNOLOGY"]
                                              ==("ELREDDDE00,ELREDDSD00,"+
                                      "ELREDDTE00,ELREDDFL00"))&
                                                (eff_data["MODE_OF_OPERATION"]==1),"VALUE"].values
                                  ).fillna(0),fill_value=0).add(
                          (hsys.loc[hsys.index.get_level_values(2).str.contains("ELST"),"fraction"].droplevel(2)
                              /eff_data.loc[(eff_data["TECHNOLOGY"]
                                              ==("ELSTDDDE00,ELSTDDSD00,"+
                                      "ELSTDDTE00,ELSTDDFL00"))&
                                                (eff_data["MODE_OF_OPERATION"]==1),"VALUE"].values).fillna(0),fill_value=0)
                          ,axis=0)
                        ).fillna(0)
    anndem_EL_hw_ts = (anndem_hw_ts.multiply(
                    (hsys.loc[hsys.index.get_level_values(2).str.contains("ASHP")
                                |hsys.index.get_level_values(2).str.contains("AWHP")
                                |hsys.index.get_level_values(2).str.contains("AAHP"),"fraction"].groupby(["LSOA11CD",
                                                                                                          "PROPERTY_TYPE"]).sum()
                                /eff_data.loc[(eff_data["TECHNOLOGY"]
                                              ==("ASHPDDDE00,ASHPDDSD00,"+
                                                  "ASHPDDTE00,ASHPDDFL00"))&(
                                              eff_data["YEAR"]
                                                              =="2020")&
                                            (eff_data["MODE_OF_OPERATION"]==2), "VALUE"].values
                            ).fillna(0).add(
                        (hsys.loc[hsys.index.get_level_values(2).str.contains("GSHP"),"fraction"].droplevel(2)
                              /eff_data.loc[(eff_data["TECHNOLOGY"]
                                              ==("GSHPDDDE00,GSHPDDSD00,"+
                                                "GSHPDDTE00,GSHPDDFL00"))&(
                                              eff_data["YEAR"]
                                                            =="2020")&
                                          (eff_data["MODE_OF_OPERATION"]==2), "VALUE"].values
                                  ).fillna(0),fill_value=0).add(
                          (hsys.loc[hsys.index.get_level_values(2).str.contains("ELRE"),"fraction"].droplevel(2)
                              /eff_data.loc[(eff_data["TECHNOLOGY"]
                                              ==("ELREDDDE00,ELREDDSD00,"+
                                      "ELREDDTE00,ELREDDFL00"))&
                                                (eff_data["MODE_OF_OPERATION"]==2),"VALUE"].values
                                  ).fillna(0),fill_value=0).add(
                          (hsys.loc[hsys.index.get_level_values(2).str.contains("ELST"),"fraction"].droplevel(2)
                              /eff_data.loc[(eff_data["TECHNOLOGY"]
                                              ==("ELSTDDDE00,ELSTDDSD00,"+
                                      "ELSTDDTE00,ELSTDDFL00"))&
                                                (eff_data["MODE_OF_OPERATION"]==2),"VALUE"].values).fillna(0),fill_value=0)
                          ,axis=0)
                        ).fillna(0)
    anndem_EL_shw_ts =  anndem_EL_sh_ts + anndem_EL_hw_ts     
    anndem_EL_shw_ts = utils.groupby_LAD(anndem_EL_shw_ts,disagg=slagg).sum()
    anndem_EL_shw_ts = anndem_EL_shw_ts.groupby([slagg.name,"TS"]).sum()
  

    # load average headroom
    raw_data = pd.read_csv(snakemake.input.path_set_dist_net)
    EL_spare_cap = raw_data.loc[raw_data["VARIABLE"]=="ELNetworkSpareCapacity",
                                "VALUE"].squeeze()
    
    # adjust spare capacity based on scenario
    if snakemake.params.dic["scen_networks"]=="-":
        pass
    else:
        EL_spare_cap = EL_spare_cap*float(snakemake.params.dic["scen_networks"])
        
    # calculate residual capacity for electricity network (GW), taking into
    # account headroom
    EL_rescap_agg = ((anndem_nhe_agg_ts
                     + anndem_EL_shw_ts).groupby([slagg.name]).max()
                     /(1-EL_spare_cap)
                     /10**6)
    

    # calculate demand for fossil gas based on residual NGBO fraction, for
    # H2 upgrade this is multiplied by blending limit
    
    # get efficiencies for gas boilers
    pt = hsys.index.get_level_values("PROPERTY_TYPE").unique()
    NGBO_eff = pd.DataFrame([[0]*2]*len(pt),index=pt,columns=["SH","HW"])
    NGBO_eff["SH"] = eff_data.loc[(eff_data["TECHNOLOGY"]
                                         ==("NGBODDDE00,NGBODDSD00,"+
                                            "NGBODDTE00,NGBODDFL00"))&
                                        (eff_data["MODE_OF_OPERATION"]==1),
                                        "VALUE"].squeeze()
    
    NGBO_eff.loc["Non-domestic","SH"] = eff_data.loc[(eff_data["TECHNOLOGY"]
                                              ==("NGBODNDO00"))&
                                             (eff_data["MODE_OF_OPERATION"]==1),
                                             "VALUE"].squeeze()
    NGBO_eff["HW"] = eff_data.loc[(eff_data["TECHNOLOGY"]
                                         ==("NGBODDDE00,NGBODDSD00,"+
                                            "NGBODDTE00,NGBODDFL00"))&
                                        (eff_data["MODE_OF_OPERATION"]==3),
                                        "VALUE"].squeeze()
    
    NGBO_eff.loc["Non-domestic","HW"] = eff_data.loc[(eff_data["TECHNOLOGY"]
                                              ==("NGBODNDO00"))&
                                             (eff_data["MODE_OF_OPERATION"]==3),
                                             "VALUE"].squeeze()
    # calculate demand taking into account residual capacity of boilers and
    # efficiency
    anndem_NG_sh_ts = anndem_sh_ts.multiply(hsys.xs("NGBO",
                                            level=2)["fraction"]
                                            /NGBO_eff["SH"]
                                            ,axis=0).dropna()
    anndem_NG_hw_ts = anndem_hw_ts.multiply(hsys.xs("NGBO",
                                            level=2)["fraction"]
                                            /NGBO_eff["HW"]
                                            ,axis=0).dropna()
    anndem_NG_ts = anndem_NG_hw_ts + anndem_NG_sh_ts
    anndem_NG_ts_agg = utils.groupby_LAD(anndem_NG_ts,disagg=slagg).sum()
    anndem_NG_ts_agg = anndem_NG_ts_agg.groupby([slagg.name,"TS"]).sum()
    
    # calculate residual capacity as fraction of peak demand in on-grid areas,
    # convert to GW

    # offgas = utils.groupby_LAD(offgas,disagg=slagg).sum()
    # offgas.columns = pd.MultiIndex.from_tuples(offgas.columns,names=["var","YEAR"])
    # offgas = offgas["noff"]/offgas["ntotal"]
    
    GA_rescap_agg = (anndem_nhg_agg_ts
                     + anndem_NG_ts_agg).groupby([slagg.name]).max()/10**6

    # calculate factor that takes into account input/output ratio adjustment
    # for hydrogen for calculating actual energy-wise blending ratio
    H2lim_fac = 1/(GA_boiler_eff["VALUE"].squeeze()
                   /H2_boiler_eff["VALUE"].squeeze()/H2blim["VALUE"].squeeze()
                   -(GA_boiler_eff["VALUE"].squeeze()
                     /H2_boiler_eff["VALUE"].squeeze()-1))
    H2_rescap_agg = GA_rescap_agg * H2lim_fac    
     
    
    # rearrange dataframes
    py = [str(i) for i in range(2015,2023)]
    DH_rescap_agg = DH_rescap_agg[py].reset_index()
    DH_rescap_agg["TECHNOLOGY"] = "DHMTTDIS00"
     
    DH_rescap_agg = DH_rescap_agg.set_index([slagg.name,"TECHNOLOGY"])
    
    DH_gen_rescap_agg = DH_gen_rescap_agg[py]
     

    ly = py[-1]
    crescapyears = range(int(ly)+1,
                         int(GA_rescap_red_start["VALUE"].squeeze()+1))
    GA_rescap_agg = GA_rescap_agg[py]
    GA_rescap_agg[[str(y) for y in crescapyears]] = 0
    GA_rescap_agg[[str(y) for y in crescapyears]]=(GA_rescap_agg[[ly]
                                                     *len(crescapyears)])
    # limit residual cap to max capacity, necessary for falling
    # demands in certain LAs
    GA_rescap_agg = GA_rescap_agg.clip(upper=demand_peak_gas)
    GA_rescap_agg = GA_rescap_agg.reset_index()
    
    GA_rescap_agg["TECHNOLOGY"] = "GAGRTDIS00"
    GA_rescap_agg = GA_rescap_agg.set_index([slagg.name,"TECHNOLOGY"])


    H2_rescap_agg = H2_rescap_agg.reset_index()
    H2_rescap_agg["TECHNOLOGY"] = "H2UGTDIS00"


    crescapyears = range(int(ly)+1,
                         int(GA_rescap_red_start["VALUE"].squeeze()+1))
    H2_rescap_agg[[str(y) for y in crescapyears]] = 0
    H2_rescap_agg[[str(y) for y in crescapyears]]=(H2_rescap_agg[[ly]
                                                     *len(crescapyears)])  
    H2_rescap_agg = H2_rescap_agg.set_index([slagg.name,"TECHNOLOGY"])

    EL_rescap_agg = EL_rescap_agg[py].reset_index()
    EL_rescap_agg["TECHNOLOGY"] = "ELGRTDIS00"

    crescapyears = range(int(ly)+1,
                         int(EL_rescap_red_start["VALUE"].squeeze()+1))
    EL_rescap_agg[[str(y) for y in crescapyears]] = 0
    EL_rescap_agg[[str(y) for y in crescapyears]]=(EL_rescap_agg[[ly]
                                                     *len(crescapyears)])  
    EL_rescap_agg = EL_rescap_agg.set_index([slagg.name,"TECHNOLOGY"])
   
    # save to file
    # unit: GW
    pd.concat([DH_rescap_agg,
               GA_rescap_agg,
               H2_rescap_agg,
               EL_rescap_agg]).to_csv(snakemake.output.path_local_res_caps)
    
    DH_gen_rescap_agg.to_csv(snakemake.output.path_dhgen_res_caps)
   