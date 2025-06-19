"""
Script deriving local governance related parameters/constraints


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

import pandas as pd
import numpy as np
import utils

if __name__ == "__main__":
    
    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_local_gov",
                                         **df_params.iloc[0].to_dict())
        
    inv_constr = pd.DataFrame(columns=["LAD23CD","TECHNOLOGY", "YEAR",
                                       "VARIABLE",
                                       "VALUE"]).set_index(["LAD23CD",
                                                            "TECHNOLOGY",
                                                            "YEAR",
                                                            "VARIABLE"])
        
    if (snakemake.params.dic["scen_local_gov"].startswith(("C","A")) and
        not snakemake.params.dic["scen_local_gov"].startswith("CP")) :
        
        lads = utils.get_entity_lookup(["LAD"]).drop("LAD23NM",axis=1)
        
        # load LA climate net-zero pledges and analysis of plans
        # 1) only pledges of existing LAs are taken into account, if LAs have
        # been merged previous targets are not counted
        # 2) pledges with 'Unclear' coverage are not taken into account
        plans = pd.read_excel(snakemake.input.path_loc_gov_pa,
                              sheet_name="individual_answers",
                              usecols=("answer_id",
                                        "local-authority-code",
                                        "question_id",
                                        "score"))
        lu = pd.read_csv(snakemake.input.path_loc_gov_lac,
                           usecols=("gss-code",
                                    "local-authority-code"),
                           index_col="local-authority-code")
        plans["LAD23CD"] = plans["local-authority-code"].map(lu["gss-code"])
        
        plans = plans[["LAD23CD","question_id","score"]].set_index(["LAD23CD","question_id"])
        hstr = plans.xs("s2_m&a_q8_sp4",level="question_id")
        hstr = lads.merge(hstr,on="LAD23CD",how="left").set_index("LAD23CD")
        hstr = hstr.astype(bool)
        
        # assume no information as no plan
        hstr = hstr.fillna(0)

        pled = pd.read_csv(snakemake.input.path_loc_gov_pl,usecols=("gss_code",
                                                                    #"council",
                                                                    "scope",
                                                                    "target"))
           
        pled = lads.merge(right=pled,left_on="LAD23CD",right_on="gss_code",
                          how="left").drop("gss_code",axis=1)
        
        pled = pled.merge(hstr,on="LAD23CD",how="left")
        
        
        # get list of LAs with no heat decarbonization plan
        nonam = pled.loc[~pled["score"],"LAD23CD"].unique()
        # nonam = pled.loc[pled["target"].isna()
        #                  |(pled["target"]>=2050),"LAD23CD"].unique()
        
        pled = pled.set_index(["LAD23CD","scope"])
        pled = pled.xs("Whole area",level=1)
        pled = pled.dropna()
        pled["target"] = pled["target"].astype(int) 
        
        # ignore target where no building strategy
        pled = pled.loc[pled["score"],"target"]
        
        # adjust so targets are between 2040 and 2050
        pled = pled.clip(lower=2040)
        pled = pled.astype(str)
        pled = pled.to_frame()
        
        # save to file - for plotting purposes only
        pled.to_csv(snakemake.output.path_local_pledges) 
        
        pled = pled.set_index("target", append=True)
        pled["nz"] = 0
        pled = pled.unstack("target").droplevel(0,axis=1)
        pled[[str(y) for y in range(min(pled.columns.astype(int)),2061) if str(y) not in pled.columns]] = np.nan
        pled = pled.reindex(sorted(pled.columns), axis=1)
        pled = pled.fillna(method="ffill",axis=1)
        pled.columns.name="YEAR"

        # stack, also removing nans
        pled = pled.stack()
        pled.name = "VALUE"
        
        # ban DH in LADs without buildings plan   
        inv_constr = pd.DataFrame(["DHMTTDIS00"],columns=["TECHNOLOGY"])
        inv_constr["YEAR"] = ":*"
        inv_constr["LAD23CD"] = ",".join(list(nonam))
        inv_constr["VARIABLE"] = "TotalAnnualMaxCapacityInvestment"
        inv_constr["VALUE"] = 0
        # disaggregate
        # slagg = pd.read_csv(snakemake.input.path_sublocal_agg,
        #                         index_col=["LSOA11CD"]).squeeze()

        
        inv_constr = inv_constr.assign(LAD23CD=inv_constr['LAD23CD'].str.split(',')).explode("LAD23CD")
        inv_constr = inv_constr[["LAD23CD","TECHNOLOGY","YEAR","VARIABLE","VALUE"]]

        inv_constr = inv_constr.set_index([col for col 
                                           in inv_constr.columns
                                           if col!="VALUE"])
        # inv_constr = utils.groupby_LAD(inv_constr,disagg=slagg).mean()
        # inv_constr.index = inv_constr.index.set_names("REGION",level=0)
        
        if "C" not in snakemake.params.dic["scen_local_gov"]:
            pled = pd.DataFrame(columns=["LAD23CD", "YEAR","VALUE"]).set_index(["LAD23CD",
                                                                        "YEAR"])
            
        if "A" not in snakemake.params.dic["scen_local_gov"]:
            inv_constr = pd.DataFrame(columns=["LAD23CD","TECHNOLOGY", "YEAR",
                                               "VARIABLE",
                                               "VALUE"]).set_index(["LAD23CD",
                                                                    "TECHNOLOGY",
                                                                    "YEAR",
                                                                    "VARIABLE"])


    elif snakemake.params.dic["scen_local_gov"]=="CP" or snakemake.params.dic["scen_local_gov"]=="CPS":
        
        lads = utils.get_entity_lookup(["LAD"]).drop("LAD23NM",axis=1)

        # load LA climate net-zero pledges
        # only pledges of existing LAs are taken into account, if LAs have
        # been merged previous targets are not counted
        pled = pd.read_csv(snakemake.input.path_loc_gov_pl,usecols=("gss_code",
                                                                    #"council",
                                                                    "scope",
                                                                    "target"))
           
        pled = lads.merge(right=pled,left_on="LAD23CD",right_on="gss_code",
                          how="left").drop("gss_code",axis=1)
    
        pled = pled.set_index(["LAD23CD","scope"])
        pled = pled.xs("Whole area",level=1)
        pled = pled.dropna()
        pled = pled.astype(int).astype(str)
        
        # save to file - for plotting purposes only
        pled.to_csv(snakemake.output.path_local_pledges) 
        
        pled = pled.set_index("target",append=True)
        pled["nz"] = 0
        pled = pled.unstack("target").droplevel(0,axis=1)
        pled[[str(y) for y in range(min(pled.columns.astype(int)),2061) if str(y) not in pled.columns]] = np.nan
        pled = pled.reindex(sorted(pled.columns), axis=1)
        pled = pled.fillna(method="ffill",axis=1)
        pled.columns.name="YEAR"

        # stack, also removing nans
        pled = pled.stack()
        pled.name = "VALUE"
    
    elif snakemake.params.dic["scen_local_gov"] == "-":
        
        pled = pd.DataFrame(columns=["LAD23CD", "YEAR","VALUE"]).set_index(["LAD23CD",
                                                                    "YEAR"])
        pled.to_csv(snakemake.output.path_local_pledges) 
        
    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_local_gov']}' is currently"
                " not a valid option of the 'scen_local_gov' parameter."
                " Valid options are 'CP' or '-'."
            )
        
        
    # rearrange dataframe 
    pled = pled.reset_index()
    pled["UNIT"] = "kt/a"
    pled["VARIABLE"] = "AnnualEmissionLimit"
    pled = pled.rename(columns={"LAD23CD": "REGION"})
    pled["EMISSION"] = "CD"

    pled = pled.set_index([col for col 
                                   in pled.columns
                                   if col!="VALUE"])    
    # save emission constraint to file
    # unit: kt/a 
    pled.to_csv(snakemake.output.path_loc_gov_em)  
    
            
    # save investment constraint to file
    # unit: GW
    inv_constr.to_csv(snakemake.output.path_loc_gov_invcon)
