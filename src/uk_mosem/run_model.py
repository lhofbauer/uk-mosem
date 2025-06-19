"""
Script to run model


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license

"""
import argparse
import json

import sys
import os
import pandas as pd


import fratoo as ft

import utils

 
ft.set_verbosity("INFO")

if __name__ == "__main__":
    
    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("run_model",
                                         **df_params.iloc[0].to_dict())

    # parse command line options
    parser = argparse.ArgumentParser(description="Parser for options for script",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-s", "--solver", help="Solver to use")
    parser.add_argument("-m", "--model", help="Model code to use")
    parser.add_argument("-r", "--run", help="Name of run")
    parser.add_argument("-p", "--processes", help="Number of processes")
    args = parser.parse_args()
    config = vars(args)
    
    
    if "snakemake" in globals():
        solver = snakemake.params.dic["run_solver"]
        processes = snakemake.params.dic["run_processes"]
        run = snakemake.params.dic["run_app"]
        model_eq = snakemake.params.dic["model_eq"]
        
    # process options
        

    td = "./temp/"
    eq = snakemake.params.fdir+snakemake.params.dic["name"]+"/"+model_eq
    d = snakemake.params.fdir+snakemake.params.dic["name"]+"/"
    
    
    # load model
    model = ft.Model(data=d,model=eq,process=True,tempdir=td)
    
     
    # define different sets of model runs
    pled = pd.read_csv(snakemake.input.path_agg_local_pledges)
    pled = pled.set_index("LAD23CD")
    ladreg = utils.get_entity_lookup(["LAD","RGN"])
    ladreg = ladreg.set_index("RGN21CD")["LAD23CD"]
    

    if run == "UK|LA|SO":
        reg= [[[r for r in model.get_entities(3)]]]
        names = [snakemake.params.dic["name"]+"_"+run]
        func = None
    elif run == "UK|L2|SO":
        reg= [[[r for r in model.get_entities(2)]]]
        names = [snakemake.params.dic["name"]+"_"+run]
        func = None
    elif run == "UK|L2|MO":
        reg= [[[r] for r in model.get_entities(2)]]
        names = [snakemake.params.dic["name"]+"_"+run]
        func = None
    elif run == "UK|LA|MO":
        en = model.get_entities(2)
        reg=[[model.ms_struct["ft_affiliation"]
         [model.ms_struct["ft_affiliation"]["VALUE"]==lad].index.to_list()
           for lad in en]]
        names = [snakemake.params.dic["name"]+"_"+run]
        func = None
    elif run == "UK|LP|SO":
        srt = list(set([r.split("|")[1] for r in model.get_entities(3)]))
        reg=[[[[r for r in model.get_entities(3) if (t in r) and (pled.loc[r.split("|")[0]].squeeze()==py)]
           for t in srt for py in pled["target"].unique() if
           [r for r in model.get_entities(3) if (t in r) and (pled.loc[r.split("|")[0]].squeeze()==py)]]+
              [[r for r in model.get_entities(2) if (pled.loc[r].squeeze()==py)]
                 for py in pled["target"].unique()]]]
        names = [snakemake.params.dic["name"]+"_"+run]
        func = None
        
    elif run == "UK|UK|SO":
        srt = list(set([r.split("|")[1] for r in model.get_entities(3)]))
        reg=[[
            [[r for r in model.get_entities(3) if t in r]
           for t in srt]+
              [[r for r in model.get_entities(2)]]
             ]]
        names = [snakemake.params.dic["name"]+"_"+run]
        func = None
    elif run == "UK|U2|SO":
        reg=[[["UK"]]]
        names = [snakemake.params.dic["name"]+"_"+run]
        func = None      
        
    elif run == "UK|CY|SO":
        srt = list(set([r.split("|")[1] for r in model.get_entities(3)]))
        reg=[[[[r for r in model.get_entities(3) if (t in r) and r.startswith(cy[0])]
           for t in srt for cy in model.get_entities(1)]+
              [[r for r in model.get_entities(2) if r.startswith(cy[0])]
                 for cy in model.get_entities(1)]]]
        names = [snakemake.params.dic["name"]+"_"+run]
        func = None
    elif run == "UK|C2|SO":
        reg=[[[r for r in model.get_entities(1)]]]
        names = [snakemake.params.dic["name"]+"_"+run]
        func = None
    elif run == "UK|C2|MO":
        reg=[[[r] for r in model.get_entities(1)]]
        names = [snakemake.params.dic["name"]+"_"+run]
        func = None
    elif run == "UK|CY|MO":
        srt = list(set([r.split("|")[1] for r in model.get_entities(3)]))
        reg=[[[[r for r in model.get_entities(3) if (t in r) and r.startswith(cy[0])]
           for t in srt]+
              [[r for r in model.get_entities(2) if r.startswith(cy[0])]]
              for cy in model.get_entities(1)]]
        names = [snakemake.params.dic["name"]+"_"+run]
        func = None
        
    elif run == "UK|RG|SO":
        srt = list(set([r.split("|")[1] for r in model.get_entities(3)]))
        reg=[[[[r for r in model.get_entities(3) if (t in r) and (r.split("|")[0] in list(ladreg.loc[rg]))]
           for t in srt for rg in ladreg.index.unique()]+
              [[r for r in model.get_entities(2) if (r in list(ladreg.loc[rg]))]
                 for rg in ladreg.index.unique()]]]
        with open("regrg.json", "w") as f:
            json.dump(reg, f)
        names = [snakemake.params.dic["name"]+"_"+run]
        func = None
        
    elif run == "UK|R2|SO":
        reg=[[[[r for r in model.get_entities(3) if (r.split("|")[0] in list(ladreg.loc[rg]))]
           for rg in ladreg.index.unique()]+
              [[r for r in model.get_entities(2) if (r in list(ladreg.loc[rg]))]
                 for rg in ladreg.index.unique()]]]
        names = [snakemake.params.dic["name"]+"_"+run]
        with open("regr2.json", "w") as f:
            json.dump(reg, f)
        func = None
    elif run == "UK|R2|MO":
        reg=[[[[r for r in model.get_entities(2) if (r in list(ladreg.loc[rg]))]]
                 for rg in ladreg.index.unique()]]
        reg=[[[[r for r in model.get_entities(3) if (r.split("|")[0] in list(ladreg.loc[rg]))]
           ]+
              [[r for r in model.get_entities(2) if (r in list(ladreg.loc[rg]))]
                 ]
              for rg in ladreg.index.unique()]]
        names = [snakemake.params.dic["name"]+"_"+run]
        with open("regr2m.json", "w") as f:
            json.dump(reg, f)
        func = None
    elif run == "UK|RG|MO":
        srt = list(set([r.split("|")[1] for r in model.get_entities(3)]))
        reg=[[[[r for r in model.get_entities(3) if (t in r) and (r.split("|")[0] in list(ladreg.loc[rg]))]
           for t in srt]+
              [[r for r in model.get_entities(2) if (r in list(ladreg.loc[rg]))]
                 ]
              for rg in ladreg.index.unique()]]
        with open("regrgmo.json", "w") as f:
            json.dump(reg, f)
        names = [snakemake.params.dic["name"]+"_"+run]
        func = None
        
    elif run == "UK|RY|SO":
        ladregd = ladreg.pop("E12000003")
        srt = list(set([r.split("|")[1] for r in model.get_entities(3)]))
        reg=[[[[r for r in model.get_entities(3) if (t in r) and (r.split("|")[0] in list(ladreg.loc[rg]))]
           for t in srt for rg in ladreg.index.unique()]+
              [r for r in model.get_entities(3) if (r.split("|")[0] in list(ladregd))]+
              [[r for r in model.get_entities(2) if (r in list(ladreg.loc[rg]))]
                 for rg in ladreg.index.unique()]]]
        with open("regry.json", "w") as f:
            json.dump(reg, f)
        names = [snakemake.params.dic["name"]+"_"+run]
        func = None
    elif run == "UK|RF|SO":
        ladregd = ladreg.loc[ladreg.isin(["E06000010","E06000011","E06000065",
                                         "E08000019","E08000035"])]
        ladreg = ladreg.loc[~ladreg.isin(["E06000010","E06000011","E06000065",
                                         "E08000019","E08000035"])]
        srt = list(set([r.split("|")[1] for r in model.get_entities(3)]))
        reg=[[[[r for r in model.get_entities(3) if (t in r) and (r.split("|")[0] in list(ladreg.loc[rg]))]
           for t in srt for rg in ladreg.index.unique()]+
              [r for r in model.get_entities(3) if (r.split("|")[0] in list(ladregd))]+
              [[r for r in model.get_entities(2) if (r in list(ladreg.loc[rg]))]
                 for rg in ladreg.index.unique()]]]
        with open("regrf.json", "w") as f:
            json.dump(reg, f)
        names = [snakemake.params.dic["name"]+"_"+run]
        func = None
    elif run.startswith("num"):
        reg= [[[r for r in model.get_entities(2)[0:int(run[3:])]]]]
        names = ["LAs"+run]
        func = None
    elif run.startswith("indnum"):
        reg=[[model.ms_struct["ft_affiliation"]
         [model.ms_struct["ft_affiliation"]["VALUE"]==lad].index.to_list()
           for lad in model.get_entities(2)[0:int(run[6:])]]]
        names = [snakemake.params.dic["name"]]
        func = None

    elif run == "all_ind_agg":
        reg= [[["UK"]],[[r] for r in model.get_entities(2)]]
        names = ["Aggregated","Individual"]
        func = None
    else:
        reg= [[["UK"]]]
        names = ["Aggregated"]
        func = None
    
            
    # run model
    duals = ["AnnualEmissionsLimit"]
    model.perform_runs(names,reg,
                        autoinclude=True, processes=processes,
                        weights="SpecifiedAnnualDemand",
                        join_results=True, overwrite=True,solver=solver,
                        duals=duals, warmstart=False)
    model.expand_results()
    
    # save results
    
    model.save_results(snakemake.params.odir)

    # create milestone file
    if not os.path.exists(os.path.dirname(snakemake.output.path_run_ms)):
        os.makedirs(os.path.dirname(snakemake.output.path_run_ms))

    pd.DataFrame().to_csv(snakemake.output.path_run_ms)
