"""
Script calculating (before diversity) peak heat demand per LSOA & property type


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

import pandas as pd

import utils


import logging

logger = logging.getLogger(__name__)
utils.adjust_logger()

if __name__ == "__main__":

    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_peak_demand",
                                         **df_params.iloc[0].to_dict())    
    logger.info("Loading load factors and annual demand data.")
    # load load factors for each technology
    loadf = pd.read_excel(snakemake.input.path_building_eff,
                          sheet_name="Efficiencies",
                          skiprows=2,header=1,usecols=[2,3,4,5],index_col=[0],
                          nrows=8)
    loadf.columns=["Balanced","Headwinds","Tailwinds"]
    loadf = loadf.rename(index={"Gas boiler":"NGBO",
                                "Oil boiler":"OIBO",
                                "Bio boiler":"BMBO",
                                "Hydrogen boiler (Hydrogen-only and Hyready boiler)":"H2BO",
                                "Electric resistive heating":"ELRE",
                                "Storage heater":"ELST"})
    loadf.index.name = "TECHNOLOGY"

    # add additional techs
    # this is irrelevant for the model results, only scales the needed capacity
    #loadf.loc["AAHP"] = loadf.loc["ASHP"]
    #loadf.loc["AWHP"] = loadf.loc["ASHP"]
    loadf.loc["HIUM"] = loadf.loc["OIBO"]
    loadf.loc["WDIS"] = loadf.loc["OIBO"]
    loadf.loc["RAUP"] = loadf.loc["OIBO"]
    loadf.loc["BELO"] = loadf.loc["OIBO"]
    loadf.loc["BEME"] = loadf.loc["OIBO"]
    loadf.loc["BEHI"] = loadf.loc["OIBO"]
    loadf.loc["BEST"] = loadf.loc["OIBO"]

    # load demand
    dw_anndem = (pd.read_csv(snakemake.input.path_dw_ann_dem_sh,
                             index_col=["LSOA11CD","PROPERTY_TYPE"])
                 + pd.read_csv(snakemake.input.path_dw_ann_dem_sh_nb,
                               index_col=["LSOA11CD","PROPERTY_TYPE"]))
    nd_anndem = (pd.read_csv(snakemake.input.path_nd_ann_dem_sh,
                  index_col=["LSOA11CD","PROPERTY_TYPE"])
                 + pd.read_csv(snakemake.input.path_nd_ann_dem_sh_nb,
                               index_col=["LSOA11CD","PROPERTY_TYPE"])).reset_index()
    nd_anndem["PROPERTY_TYPE"]="Non-domestic"
    nd_anndem = nd_anndem.groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
    anndem_SH = pd.concat([dw_anndem,nd_anndem])
    
#%%    
    # calculate per dwelling peak demand/unit size, for limited further use,
    # see below
    logger.info("Calculate per dwelling peak demand.")
    nd_raw_data = pd.read_csv(snakemake.input.path_set_nd_prop)
    ND_afa = nd_raw_data["VALUE"][nd_raw_data["VARIABLE"]=="NonDomesticAverageFloorarea"].squeeze()

    # load number of dwellings (domestic) and floor space (non-dom.) per LSOA
    dwel_ex = pd.read_csv(snakemake.input.path_dw_stock_ex,
                          index_col=["LSOA11CD","PROPERTY_TYPE"])
    dwel_nb = pd.read_csv(snakemake.input.path_dw_stock_nb,
                          index_col=["LSOA11CD","PROPERTY_TYPE"])
    nd_stock_ex =  pd.read_csv(snakemake.input.path_nd_stock_ex,
                               index_col=["LSOA11CD","PROPERTY_TYPE"])    
    nd_stock_nb =  pd.read_csv(snakemake.input.path_nd_stock_nb,
                               index_col=["LSOA11CD","PROPERTY_TYPE"])
    nd_stock = nd_stock_ex+nd_stock_nb
    nd_stock = nd_stock.reset_index()
    nd_stock["PROPERTY_TYPE"] = "Non-domestic"
    
    # ND_afa is the average based on the EPC/DEC data
    nd_stock = nd_stock.groupby(["LSOA11CD","PROPERTY_TYPE"]).sum() /ND_afa
    stock = pd.concat([dwel_ex+dwel_nb,nd_stock])
    SH_int = anndem_SH/stock
    SH_int = SH_int.dropna(how="all")
    
    lf = pd.concat([SH_int.index.to_frame(),loadf.T],axis=1)
    lf = lf.set_index(["LSOA11CD","PROPERTY_TYPE"])
    lf.columns.name="TECHNOLOGY"
    lf = lf.stack(dropna=False)
    lf = lf.reset_index().merge(right=loadf,on="TECHNOLOGY").drop(0,axis=1)
    
    caps = SH_int.merge(right=lf,on=["LSOA11CD","PROPERTY_TYPE"])
    
    caps = caps.set_index(["LSOA11CD","PROPERTY_TYPE","TECHNOLOGY"])
    caps.loc[:,"2015":"2060"] = caps.loc[:,"2015":"2060"].divide((caps["Balanced"]
                                                           *8760*3600),axis=0)*10**6
    caps = caps.loc[:,"2015":"2060"]
    #caps = utils.groupby_LAD(caps).mean()
    
    caps_avg = (caps * stock).groupby("TECHNOLOGY").sum()/stock.sum()
    caps_avg["LSOA11CD"] = "AVERAGE"
    caps_avg["PROPERTY_TYPE"] = "AVERAGE"
    caps_avg = caps_avg.reset_index().set_index(["LSOA11CD","PROPERTY_TYPE",
                                                 "TECHNOLOGY"])
    caps_avg_det = (caps * stock).groupby(["PROPERTY_TYPE",
                                           "TECHNOLOGY"]).sum()/stock.groupby("PROPERTY_TYPE").sum()
    caps_avg_det["LSOA11CD"] = "AVERAGE"
    caps_avg_det = caps_avg_det.reset_index().set_index(["LSOA11CD",
                                                         "PROPERTY_TYPE",
                                                         "TECHNOLOGY"])
    
    
    caps_avg_dw = ((caps.drop("Non-domestic",level=1)
                   *stock.drop("Non-domestic",level=1)
                   ).groupby("TECHNOLOGY").sum()
                   /stock.drop("Non-domestic",level=1).sum())
    caps_avg_dw["LSOA11CD"] = "AVERAGE"
    caps_avg_dw["PROPERTY_TYPE"] = "Domestic"
    caps_avg_dw = caps_avg_dw.reset_index().set_index(["LSOA11CD","PROPERTY_TYPE",
                                                 "TECHNOLOGY"])
    dw = ["Detached","Flats","Semi-detached","Terraced"]
    if snakemake.params.dic["scen_hh_disagg"]=="T":
        dw = [p+"|"+t for p in ["Detached","Flats","Semi-detached","Terraced"]
              for t in ["OO","RP","RS"]]
    elif snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        dw = [p+"|"+t for p in ["Detached","Flats","Semi-detached","Terraced"]
              for t in [t+str(i) for t in ["O","P","S"]
              for i in range(len(snakemake.params.dic["scen_hh_disagg"].split("|")[1])+1)]]
    elif snakemake.params.dic["scen_hh_disagg"]=="-":                             
        pass
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_hh_disagg']}' is currently"
                " not a valid option of the 'scen_hh_disagg' parameter."
                " Valid options are 'T', 'TI|[1-9DIGITS]', or '-'."
            )  
    caps.to_csv("caps.csv")
    stock.to_csv("stock.csv")
    SH_int.to_csv("SH_int.csv")
    caps_avg_nd = ((caps.drop(dw,level=1)
                   *stock.drop(dw,level=1)
                   ).groupby("TECHNOLOGY").sum()
                   /stock.drop(dw,level=1).sum())
    caps_avg_nd["LSOA11CD"] = "AVERAGE"
    caps_avg_nd["PROPERTY_TYPE"] = "Non-domestic"
    caps_avg_nd = caps_avg_nd.reset_index().set_index(["LSOA11CD","PROPERTY_TYPE",
                                                 "TECHNOLOGY"])
    
    caps = pd.concat([caps_avg_det,caps])
    
    if snakemake.params.dic["scen_hh_disagg"].startswith("T"):
        caps.index = pd.MultiIndex.from_tuples([ (x[0],x[1],x[2]+"DDDE"+x[1].split("|")[1])
                                                  if x[1].startswith("Detached") else
                                                  (x[0],x[1],x[2]+"DDSD"+x[1].split("|")[1])
                                                  if x[1].startswith("Semi-detached") else
                                                  (x[0],x[1],x[2]+"DDTE"+x[1].split("|")[1])
                                                  if x[1].startswith("Terraced") else
                                                  (x[0],x[1],x[2]+"DDFL"+x[1].split("|")[1])
                                                  if x[1].startswith("Flats") else
                                                  (x[0],x[1],x[2]+"DNDO00")
                                                  if x[1]=="Non-domestic" else 0 for x
                                                  in caps.index])
    else:              
        caps.index = pd.MultiIndex.from_tuples([ (x[0],x[1],x[2]+"DDDE00")
                                                         if x[1]=="Detached" else
                                                         (x[0],x[1],x[2]+"DDSD00")
                                                         if x[1]=="Semi-detached" else
                                                         (x[0],x[1],x[2]+"DDTE00")
                                                         if x[1]=="Terraced" else
                                                         (x[0],x[1],x[2]+"DDFL00")
                                                         if x[1]=="Flats" else
                                                         (x[0],x[1],x[2]+"DNDO00")
                                                         if x[1]=="Non-domestic" else 0 for x
                                                         in caps.index])
        

    caps.index.names = ["LSOA11CD","PROPERTY_TYPE","TECHNOLOGY"]
    
    caps = pd.concat([caps_avg,caps_avg_dw,caps_avg_nd,caps])
    
    
    # save to file (for heat pump constraint - to get from number of HPs to
    # capacity and for reference -  when manually choosing technology 
    # cost brackets)
    # unit: kW
    caps.loc[:,"2015":"2060"].to_csv(snakemake.output.path_dw_peakd)  
#%%    
    # calculate per LSOA (not per dwelling) peak demands
    logger.info("Calculate per LSOA peak demand.")
    
    lf = pd.concat([pd.DataFrame(index=anndem_SH.index).reset_index(),loadf.T],axis=1)
    lf = lf.set_index(["LSOA11CD","PROPERTY_TYPE"])
    lf.columns.name="TECHNOLOGY"
    lf = lf.stack(dropna=False)
    lf = lf.reset_index().merge(right=loadf,on="TECHNOLOGY").drop(0,axis=1)    

    capsa = anndem_SH.merge(right=lf[["LSOA11CD","PROPERTY_TYPE",
                                    "TECHNOLOGY","Balanced"]],
                           on=["LSOA11CD","PROPERTY_TYPE"])
    
    capsa = capsa.set_index(["LSOA11CD","PROPERTY_TYPE","TECHNOLOGY"])
    capsa = capsa.divide(capsa["Balanced"]*8760*3600,axis=0)
    capsa = capsa.drop("Balanced",axis=1)
    
    # save to file
    # unit: GW    
    capsa.to_csv(snakemake.output.path_peakc_lsoa) 
