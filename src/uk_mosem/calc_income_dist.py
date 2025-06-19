"""
Script deriving gross household income distribution for LSOAs


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""
import logging

import pandas as pd
import numpy as np
from scipy.stats import lognorm, norm

import utils

logger = logging.getLogger(__name__)
utils.adjust_logger()

if __name__ == "__main__":

    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_income_dist",
                                         **df_params.iloc[43].to_dict())
        
    # load income data
    dist = pd.read_excel(snakemake.input.path_ec_ghhi_add,
                         sheet_name="Gross Occupied Address LSOA",
                         skiprows=2)
    dist_reg = pd.read_excel(snakemake.input.path_ec_ghhi_add,
                             sheet_name="Gross Occupied Address Region",
                             skiprows=2)
    

    br = dist_reg.loc[dist_reg["Region Code"] == "E92000001"]
    
    ### calculate all percentiles (for results post processing)
    
    pe = br.copy()
    # function to get log normal distributions
    def get_lognorm_from_perc(v_a, perc_a, v_b, perc_b):

        v_a = np.log(v_a)
        v_b = np.log(v_b)
        perc_a_p = norm.ppf(perc_a)
        perc_b_p = norm.ppf(perc_b)
        
        scale = (v_b - v_a) / (perc_b_p - perc_a_p)
        mean = ((v_a * perc_b_p) - (v_b * perc_a_p)) / (perc_b_p - perc_a_p)
        return lognorm(s=scale, scale=np.exp(mean))
    
    pe = pe.drop([c for c in pe.columns if "percentile" not in c],axis=1)
    pe.columns = pe.columns.str[0:2]
    pe = pe.T
    pe = pe.reset_index()
    pe.columns = ["PERCENTILE","VALUE"]
    pe = pe.astype(int)
    
    # estimate distributions
    dists = list()
    for i in range(len(pe.iloc[:-1,:])):
        ln = get_lognorm_from_perc(pe.iloc[i,1],pe.iloc[i,0]/100,
                                    pe.iloc[i+1,1], pe.iloc[i+1,0]/100)
        dists.append(ln)
    # add dist for first and last
    dists.insert(0,dists[0])
    dists.append(dists[-1])
    
    # calculate all percentiles
    pe.loc[9] = [0,0]
    pe = pe.sort_values("PERCENTILE")
    pe = pe.reset_index(drop=True)
    pe.loc[:,"INCOME"] = ""
    for i,p in enumerate(pe["PERCENTILE"]):
        pe.iloc[i,2]= ",".join([str(dists[i].ppf((p+ii)/100)) for ii in range(1,10)])
    pe.loc[:,"INCOME"] = pe.loc[:,"VALUE"].astype(str)+","+pe.loc[:,"INCOME"]
    pe = pe.assign(INCOME=pe["INCOME"].str.split(',')).explode("INCOME")
    pe = pe.reset_index(drop=True)
    pe.loc[100] = dists[9].ppf(0.995)
    pe["PERCENTILE"] = pe.index.values
    pe = pe.drop("VALUE",axis=1)
    
    
    ### calculate local allocation in national brackets
    
    # set national brackets
    if snakemake.params.dic["scen_hh_disagg"].startswith("TI"):
        brackets = snakemake.params.dic["scen_hh_disagg"].split("|")[1]
        br = br.loc[:,[i+'0th percentile (£)' for i in brackets]]
    else: # default option
        br = br.loc[:,['10th percentile (£)',
                       '20th percentile (£)',
                       '40th percentile (£)']]   
    br.columns = [c[0:2]+"p" for c in br.columns]
    br.index = ["GHHI"]
    # br.to_csv(snakemake.output.path_ec_ghhi_dist_br)
    
    # final processing of above dist based on brackets
    pe.loc[:,"BRACKET"] = 0
    brs = [0]+list(br.columns.str[0:2].astype(int))+[100]
    for i in range(0,4):
        pe.loc[(pe["PERCENTILE"]<=brs[i+1])&
               (pe["PERCENTILE"]>brs[i]),"BRACKET"] = i
        
    pe.to_csv(snakemake.output.path_ec_ghhi_dist_perc,index=False)  

    
    # apply national brackets to LSOAs
    dist = dist.set_index("LSOA code")
    dist = dist.drop([c for c in dist.columns if "th perc" not in c],axis=1)
    dist.columns = [c[0:2]+"p" for c in dist.columns]  
    
    # estimate detailed local distributions
    di = pd.DataFrame([],columns=[i for i in range(0,101,1)],dtype=float)
    
    for r in range(len(dist.index)):
        logger.info(str(r))
        pe = dist.iloc[r,:].T.reset_index()
        pe["index"] = pe["index"].str[0:-1].astype(int)
        name = pe.columns[1]
        pe.columns = ["PERCENTILE","VALUE"]
        dists = list()
        for i in range(len(pe.iloc[:-1,:])):
            ln = get_lognorm_from_perc(pe.iloc[i,1],pe.iloc[i,0]/100,
                                        pe.iloc[i+1,1], pe.iloc[i+1,0]/100)
            dists.append(ln)
        # add dist for first and last
        dists.insert(0,dists[0])
        dists.append(dists[-1])
        
        # calculate all percentiles
        pe.loc[9] = [0,0]
        pe = pe.sort_values("PERCENTILE")
        pe = pe.reset_index(drop=True)
        pe.loc[:,"INCOME"] = ""
        for i,p in enumerate(pe["PERCENTILE"]):
            pe.iloc[i,2]= ",".join([str(dists[i].ppf((p+ii)/100)) for ii in range(1,10)])
        pe.loc[:,"INCOME"] = pe.loc[:,"VALUE"].astype(str)+","+pe.loc[:,"INCOME"]
        pe = pe.assign(INCOME=pe["INCOME"].str.split(',')).explode("INCOME")
        pe = pe.reset_index(drop=True)
        pe.loc[100] = dists[9].ppf(0.995)
        pe["PERCENTILE"] = pe.index.values
        pe = pe.drop("VALUE",axis=1)
        pe["INCOME"] = pe["INCOME"].astype(float)
        
        di.loc[name,:] = pe["INCOME"].values
    
    di.columns = [str(c)+"p" for c in di.columns]
    dist = di 
    
    
    dist_base = dist.copy()
    
    # go through each national bracket, calculate LSOA percentiles
    # this uses linear interpolation as unclear how local dist looks like
    for b,bv in zip(br.columns,br.iloc[0]):
        
        dist["perc"] = ((dist_base<=bv).astype(int).sum(axis=1)-1)*0.01
        up = (dist_base>bv).idxmax(axis=1)
        low = ((up.str[0:-1]).astype(int)-1).astype(str)+"p"
        idx, cols = pd.factorize(up)
        dist["up"]=dist_base.reindex(cols, axis=1).to_numpy()[np.arange(len(dist_base)), idx]
        idx, cols = pd.factorize(low)
        dist["low"]=dist_base.reindex(cols, axis=1).to_numpy()[np.arange(len(dist_base)), idx] 
        dist["addperc"] = (bv-dist["low"])/(dist["up"]-dist["low"])*0.01
        dist[b+"n"] = dist["perc"] + dist["addperc"]
    
    dist = dist.loc[:,[c+"n" for c in br.columns]]
    
    dist_base = dist.copy()
    for i,c in enumerate(dist.columns[1:]):
        dist.loc[:,c] = dist_base.loc[:,c]-dist_base.iloc[:,i]
     
    dist.loc[:,"100pn"] = 1-dist.sum(axis=1)
    dist.index.names = ["LSOA11CD"]    
    # save to file
    # unit: -
    dist.to_csv(snakemake.output.path_ec_ghhi_dist)


