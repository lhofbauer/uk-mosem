"""
Script deriving annual demands for each local area


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

import pandas as pd
import numpy as np

import utils

if __name__ == "__main__":
    
    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_annual_demand",
                                         **df_params.iloc[0].to_dict())    

    # load number of dwellings (domestic) and floor space (non-dom.) per LSOA
  
    dwel_ex = pd.read_csv(snakemake.input.path_dw_stock_ex,
                          index_col=["LSOA11CD","PROPERTY_TYPE"])
    dwel_nb = pd.read_csv(snakemake.input.path_dw_stock_nb,
                          index_col=["LSOA11CD","PROPERTY_TYPE"])
    nd_stock_ex =  pd.read_csv(snakemake.input.path_nd_stock_ex,
                               index_col=["LSOA11CD","PROPERTY_TYPE"])    
    nd_stock_nb =  pd.read_csv(snakemake.input.path_nd_stock_nb,
                               index_col=["LSOA11CD","PROPERTY_TYPE"])   
                     
    # dwel_existing = dwel_existing.loc[:,"2015":]
    # dwel_new = dwel_new.loc[:,"2015":]
    demands = [ (snakemake.input.path_dw_int_sh,
                 snakemake.output.path_dw_ann_dem_sh,"sh_dw_ex"),
                (snakemake.input.path_dw_int_hw,
                 snakemake.output.path_dw_ann_dem_hw,"hw_dw_ex_nb"),
                (snakemake.input.path_dw_int_sh_nb,
                 snakemake.output.path_dw_ann_dem_sh_nb,"sh_dw_nb"),
                (snakemake.input.path_dw_int_nhe,
                 snakemake.output.path_dw_ann_dem_nhe,"nhe_dw_ex_nb"),
                (snakemake.input.path_dw_int_nhg,
                 snakemake.output.path_dw_ann_dem_nhg,"nhg_dw_ex_nb"),
                (snakemake.input.path_nd_int_sh,
                 snakemake.output.path_nd_ann_dem_sh,"sh_nd_ex"),
                (snakemake.input.path_nd_int_hw,
                 snakemake.output.path_nd_ann_dem_hw,"hw_nd_ex_nb"),
                (snakemake.input.path_nd_int_sh_nb,
                 snakemake.output.path_nd_ann_dem_sh_nb,"sh_nd_nb"),
                (snakemake.input.path_nd_int_nhe,
                 snakemake.output.path_nd_ann_dem_nhe,"nhe_nd_ex_nb"),
                (snakemake.input.path_nd_int_nhg,
                 snakemake.output.path_nd_ann_dem_nhg,"nhg_nd_ex_nb")]

    anndem = dict()
    
    for i,d,b in demands:
        
        # load demand intensity per building type/floor space and LSOA, for
        # domestic or non-domestic properties
        if "dw" in b:
            demint = pd.read_csv(i,
                                 index_col=["LSOA11CD","PROPERTY_TYPE"])
            # multiply number of dwellings by demand intensities and convert to GJ
            if "ex_nb" in b:
                anndem[b] = (dwel_ex+dwel_nb).multiply(demint) * 3.6 / 1000
            elif "ex" in b:
                anndem[b] = dwel_ex.multiply(demint) * 3.6 / 1000
            elif "nb" in b:
                anndem[b] = dwel_nb.multiply(demint) * 3.6 / 1000
        elif "nd" in b:
            demint = pd.read_csv(i,
                                 index_col=["LSOA11CD","PROPERTY_TYPE"])
            # multiply number of dwellings by demand intensities and convert to GJ
            if "ex_nb" in b:
                anndem[b] = (nd_stock_ex+nd_stock_nb).multiply(demint) * 3.6 / 1000
            elif "ex" in b:
                anndem[b] = nd_stock_ex.multiply(demint)* 3.6 / 1000
            elif "nb" in b:
                anndem[b] = nd_stock_nb.multiply(demint)* 3.6 / 1000   
                 
        anndem[b] = anndem[b].fillna(0)
        
    

    if snakemake.params.dic["scen_dem_calib"]=="-":
        pass
    # if required, calibrate annual demands based on electricity and gas
    # consumption statistics (at MSOA level)     
    elif snakemake.params.dic["scen_dem_calib"]=="MSOA":

        
        cdata = dict()
        cdata["gas_dw"] = pd.read_csv(snakemake.input.path_dw_con_g)
        cdata["gas_nd"] = pd.read_csv(snakemake.input.path_nd_con_g)
        cdata["el_dw"] = pd.read_csv(snakemake.input.path_dw_con_e)
        cdata["el_nd"] = pd.read_csv(snakemake.input.path_nd_con_e)
        
        for k in cdata.keys():
            if k == "gas_dw":
                cdata[k] = cdata[k].rename(columns={"LACode":"LAD21CD",
                                      "MSOACode":"MSOA11CD",
                                      "msoa.name":"MSOA11CD",
                                      "year":"YEAR",
                                      "total.kwh":"kWh",
                                      "KWH":"kWh"})
            else:
                cdata[k] = cdata[k].rename(columns={"LACode":"LAD21CD",
                                      "MSOACode":"MSOA11CD",
                                      "msoa.code":"MSOA11CD",
                                      "year":"YEAR",
                                      "total.kwh":"kWh",
                                      "KWH":"kWh"})
                
        # aggregate residential data to MSOA
        for k in cdata.keys():
            if k.endswith("dw"):  
                cdata[k] = cdata[k][["MSOA11CD","YEAR","kWh"]].groupby(["MSOA11CD","YEAR"]).sum()
                cdata[k] = cdata[k].reset_index()
    
        cdata["gas_nd"] = cdata["gas_nd"][["MSOA11CD","YEAR","kWh"]]
        cdata["el_nd"] = cdata["el_nd"][["LAD21CD","MSOA11CD","YEAR","kWh"]]
        
        # reallocate unallocated meters
        for k in ["el_dw","gas_dw","gas_nd"]:
            sf = (cdata[k][cdata[k]["MSOA11CD"].str.contains("Unallocated")][["YEAR","kWh"]].groupby("YEAR").sum()
                  /cdata[k][~cdata[k]["MSOA11CD"].str.contains("Unallocated")][["YEAR","kWh"]].groupby("YEAR").sum())
            cdata[k] = cdata[k].set_index(["MSOA11CD","YEAR"])
            cdata[k] = cdata[k][~cdata[k].index.get_level_values("MSOA11CD").str.contains("Unallocated")].multiply((1+sf),level=1)
            cdata[k] = cdata[k].reset_index()
            
        # allocate half-hourly LAD-level electricity consumption data across MSOAs
        # reallocate overall unallocated values
        sf = (cdata["el_nd"][cdata["el_nd"]["LAD21CD"].str.contains("Unallocated")][["YEAR","kWh"]].groupby("YEAR").sum()
              /cdata["el_nd"][~cdata["el_nd"]["LAD21CD"].str.contains("Unallocated")][["YEAR","kWh"]].groupby("YEAR").sum())
        
        cdata["el_nd"] = cdata["el_nd"].set_index(["LAD21CD","MSOA11CD","YEAR"])
        cdata["el_nd"] = cdata["el_nd"][~cdata["el_nd"].index.get_level_values("LAD21CD").str.contains("Unallocated")].multiply((1+sf),level=2)
        cdata["el_nd"] = cdata["el_nd"].reset_index()
        
        sf = (1-(cdata["el_nd"][cdata["el_nd"]["MSOA11CD"].str.contains("Half-Hourly")].groupby(["LAD21CD","YEAR"]).sum()["kWh"]/
              cdata["el_nd"].groupby(["LAD21CD","YEAR"])["kWh"].sum()).fillna(0))
        sf.name = "sf"
        sf = sf.reset_index()
        
        cdata["el_nd"] = cdata["el_nd"].merge(right=sf, how="left", on=["LAD21CD","YEAR"])
        cdata["el_nd"]["kWh"] = cdata["el_nd"]["kWh"]/cdata["el_nd"]["sf"]
        cdata["el_nd"] = cdata["el_nd"].drop("sf",axis=1)
        
        cdata["el_nd"] = cdata["el_nd"][~cdata["el_nd"]["MSOA11CD"].str.contains("Half-Hourly")]
        
        for k in cdata.keys():
            cdata[k]["YEAR"] = cdata[k]["YEAR"].astype(int).astype(str)
     
        
        # gas = pd.concat([gas_dw,gas_nd]).groupby(["MSOA11CD","YEAR"]).sum()
        # gas = gas.reset_index()
        # el = pd.concat([el_dw,el_nd]).groupby(["MSOA11CD","YEAR"]).sum()
        # el = el.reset_index()
        
        lsoas = utils.get_entity_lookup(["LSOA","MSOA"])
        for k in cdata.keys():
            cdata[k] = lsoas.drop("LSOA11CD",
                                  axis=1).drop_duplicates().merge(right=cdata[k],
                                                                  how="left",
                                                                  on="MSOA11CD").fillna(0) 
        for k in cdata.keys():
            # recalculate as GJ
            cdata[k]["kWh"] = cdata[k]["kWh"]* 3.6/10**3                                                 
            cdata[k] = cdata[k][~cdata[k]["kWh"].isna()]
            cdata[k] = cdata[k][["MSOA11CD","YEAR","kWh"]].set_index(["MSOA11CD","YEAR"])["kWh"].unstack().dropna(axis=1,how="all")
        
        # adjust non-domestic consumption to account for end uses that are
        # not capture in model (industrial processes) but in the consumption
        # statistics
        
        ECUK = pd.read_excel(snakemake.input.path_ecuk_enduse,
                             sheet_name="Table U2",
                             usecols=(0,1,3,4,7),skiprows=4,nrows=400)
        ECUK = ECUK.set_index(["Sector","End use","Year"]).sort_index()

        noncapf = (ECUK.loc["Industry",
                                ["High temperature process",
                                      "Low temperature process",
                                      "Drying/separation",
                                      "Other",
                                      "Motors",
                                      "Compressed air"
                                      ],
                                slice("2017","2022")
                               ].groupby("Year").sum()
                       /
                       ECUK.loc[slice("Industry","Services"),
                               "Overall total",
                                slice("2017","2022"),
                               ].groupby("Year").sum()
                       )
        noncapf.index.name = "YEAR"
        noncapf.index = noncapf.index.astype(str)
        # assume 2017 value for 2015-2016
        noncapf.loc["2015"] = noncapf.loc["2017"]
        noncapf.loc["2016"] = noncapf.loc["2017"]
        noncapf = noncapf.sort_index()
        noncapf.columns = ["gas_nd","el_nd"]
        #MSOA_alloc = nd_stock_ex.xs("Industrial",level=1).loc[:,"2015":"2021"]
        #MSOA_alloc = lsoas.merge(right=MSOA_alloc, how="left", on="LSOA11CD").groupby("MSOA11CD").sum()
        #MSOA_alloc = MSOA_alloc/MSOA_alloc.sum()
        
        #MSOA_alloc = nd_stock_ex.groupby("LSOA11CD").sum().loc[:,"2015":"2021"]*MSOA_alloc
        
        #MSOA_alloc = lsoas.merge(right=MSOA_alloc, how="left", on="LSOA11CD").groupby("MSOA11CD").sum()
        for k in ["gas_nd","el_nd"]:
            MSOA_alloc = cdata[k]
            MSOA_alloc = MSOA_alloc/MSOA_alloc.sum()

            non_alloc_con = MSOA_alloc.multiply((cdata[k].sum()*noncapf[k]).dropna()).dropna(axis=1)
            cdata[k] = cdata[k]-non_alloc_con
    
    
        ad = anndem.copy()
        ad["sh_dw_ex_nb"] = ad["sh_dw_ex"].add(ad["sh_dw_nb"],
                                                        fill_value=0)
        ad["sh_nd_ex_nb"] = ad["sh_nd_ex"].add(ad["sh_nd_nb"],
                                                        fill_value=0)
        del ad["sh_nd_nb"], ad["sh_dw_nb"], ad["sh_nd_ex"], ad["sh_dw_ex"]
        
        for k in ad.keys():
            if "nd" in k:
                ad[k] = ad[k].reset_index()
                ad[k]["PROPERTY_TYPE"] = "Non-domestic"
                ad[k] =  ad[k].groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
       
        
        # anndem_sh = pd.concat([anndem_sh_dw,anndem_sh_nd])
        # anndem_hw = pd.concat([anndem_hw_dw,anndem_hw_nd])
        # anndem_nhe = pd.concat([anndem_nhe_dw,anndem_nhe_nd])
        # anndem_nhg = pd.concat([anndem_nhg_dw,anndem_nhg_nd])
        
        
        # load heating tech residual fractions
        hsys = pd.read_csv(snakemake.input.path_res_caps_ht_f, index_col=["LSOA11CD",
                                                                "PROPERTY_TYPE",
                                                                "TECHNOLOGY"])
        eff_data = pd.read_csv(snakemake.input.path_set_eff)
        
        # calculate gas and electricity consumption for SH and HW
        for k in list(ad.keys()):
            if "sh" in k:
                ad[k+"_el"] = (ad[k].multiply(
                                (hsys.loc[hsys.index.get_level_values(2).str.contains("ASHP")
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
                ad[k+"_ng"] = (ad[k].multiply(
                                (hsys.loc[hsys.index.get_level_values(2).str.contains("NGBO"),
                                          "fraction"].groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
                                            /eff_data.loc[(eff_data["TECHNOLOGY"]
                                                          ==("NGBODDDE00,NGBODDSD00,"+
                                                              "NGBODDTE00,NGBODDFL00"))&
                                                        (eff_data["MODE_OF_OPERATION"]==1), "VALUE"].values
                                        ).fillna(0),axis=0)
                                    ).fillna(0)
                del ad[k]
                
                                                  
            elif "hw" in k:
                ad[k+"_el"] = (ad[k].multiply(
                                (hsys.loc[hsys.index.get_level_values(2).str.contains("ASHP")
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
                                                  
               
                ad[k+"_ng"] = (ad[k].multiply(
                                (hsys.loc[hsys.index.get_level_values(2).str.contains("NGBO"),
                                          "fraction"].groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
                                            /eff_data.loc[(eff_data["TECHNOLOGY"]
                                                          ==("NGBODDDE00,NGBODDSD00,"+
                                                              "NGBODDTE00,NGBODDFL00"))&
                                                        (eff_data["MODE_OF_OPERATION"]==3), "VALUE"].values
                                        ).fillna(0),axis=0)
                                    ).fillna(0)
                del ad[k]
        # aggregate data for SH and HW
        ad["dw_ex_nb_el"] = pd.concat([ad[k] for k in ad.keys()
                                       if ("dw" in k)
                                       and ("el" in k)]).groupby(["LSOA11CD",
                                                                  "PROPERTY_TYPE"]).sum()
        ad["dw_ex_nb_ng"] = pd.concat([ad[k] for k in ad.keys()
                                       if ("dw" in k)
                                       and ("ng" in k)]).groupby(["LSOA11CD",
                                                                  "PROPERTY_TYPE"]).sum()
        ad["nd_ex_nb_el"] = pd.concat([ad[k] for k in ad.keys()
                                       if ("nd" in k)
                                       and ("el" in k)]).groupby(["LSOA11CD",
                                                                  "PROPERTY_TYPE"]).sum()
        ad["nd_ex_nb_ng"] = pd.concat([ad[k] for k in ad.keys()
                                       if ("nd" in k)
                                       and ("ng" in k)]).groupby(["LSOA11CD",
                                                                  "PROPERTY_TYPE"]).sum()
        for k in list(ad.keys()):
            if ("sh" in k) or ("hw" in k):
                del ad[k]
                
        # aggregate to MSOAs    
        for k in ad.keys():
            ad[k] = lsoas.merge(right=ad[k], how="left", on="LSOA11CD").groupby("MSOA11CD").sum()
    
        # calculate scaling factors based on gas/electricity consumption from
        # empirical data and model
    
        years = ["2015","2016","2017","2018","2019","2020","2021"]
        sf = dict()
        sf["dw_gas"] = (cdata["gas_dw"][years].sum(axis=1)
                        /(ad["dw_ex_nb_ng"][years]+ad["nhg_dw_ex_nb"][years]).sum(axis=1)).fillna(1)
        
        
        sf["dw_el"] = ((cdata["el_dw"][years].sum(axis=1)-ad["dw_ex_nb_el"][years].multiply(sf["dw_gas"],axis=0).sum(axis=1))
                        /(ad["nhe_dw_ex_nb"][years]).sum(axis=1)).fillna(1).replace([np.inf, -np.inf],
                                                                       1)
        # FIXME: check if with updated residual caps there are still factors < 0
        # currently these are set to first percentile, i.e., electricity demand
        # in the model (just from heating) will be higher than actual consumption
        # data (while NHE is then set scale to lower bound) - but this only 
        # happens for around 65 (out of 8480 MSOAs) currently
                                                                     
        sf["dw_el"].loc[sf["dw_el"]<
                        sf["dw_el"].quantile(0.01)] = sf["dw_el"].quantile(0.01)  
        #(z - ngfac*a)/b
        # FIXME: this ignores a very few instances where there is actual 
        # consumption but modelled demand is 0 (these are kept at 0)
        sf["nd_gas"] = ((cdata["gas_nd"][years].sum(axis=1)
                        /(ad["nd_ex_nb_ng"][years]
                          +ad["nhg_nd_ex_nb"][years]).sum(axis=1)).fillna(1).replace([np.inf, -np.inf],
                                                                         1))
                                                                                     
        # some factor are very high which can, for example, be due to industrial
        # use of gas not capture in model - cap value based on 99th percentile
        sf["nd_gas"].loc[sf["nd_gas"]>
                        sf["nd_gas"].quantile(0.99)] = sf["nd_gas"].quantile(0.99)  
        
                                                                     
        sf["nd_el"] = (((cdata["el_nd"][years].sum(axis=1)-ad["nd_ex_nb_el"][years].multiply(sf["nd_gas"],axis=0).sum(axis=1))
                        /ad["nhe_nd_ex_nb"][years].sum(axis=1)).fillna(1).replace([np.inf, -np.inf],
                                                                         1))
                                                                                  
        # FIXME: See above - here 211 MSOAs
        sf["nd_el"].loc[sf["nd_el"]>
                        sf["nd_el"].quantile(0.99)] = sf["nd_el"].quantile(0.99)  
        sf["nd_el"].loc[sf["nd_el"]<
                        sf["nd_el"].quantile(0.05)] = sf["nd_el"].quantile(0.05)  
        # comparison for checks                                                                    
        # dw_comp_gas = (cdata["gas_dw"][years]-ad["dw_ex_nb_ng"][years]-ad["nhg_dw_ex_nb"][years]).describe()
        # dw_comp_el = (cdata["el_dw"][years]-ad["dw_ex_nb_el"][years]-ad["nhe_dw_ex_nb"][years]).describe()
        # nd_comp_gas = (cdata["gas_nd"][years]-ad["nd_ex_nb_ng"][years]-ad["nhg_nd_ex_nb"][years]).describe()
        # nd_comp_el = (cdata["el_nd"][years]-ad["nd_ex_nb_el"][years]-ad["nhe_nd_ex_nb"][years]).describe()
        
        for k in sf.keys():
            sf[k].name ="sf"
            sf[k] = sf[k].reset_index()
            sf[k] = lsoas.merge(right=sf[k],on="MSOA11CD",how="right")
            sf[k] = sf[k].set_index("LSOA11CD")
            sf[k] = sf[k].drop("MSOA11CD",axis=1)
            
        for k in anndem.keys():
            if (("sh" in k) or ("hw" in k) or ("nhg" in k)) and ("dw" in k):
                anndem[k] = anndem[k].multiply(sf["dw_gas"]["sf"],
                                               axis=0,
                                               level=0)
            if (("sh" in k) or ("hw" in k) or ("nhg" in k)) and ("nd" in k):
                anndem[k] = anndem[k].multiply(sf["nd_gas"]["sf"],
                                               axis=0,
                                               level=0)
            if ("nhe" in k) and ("dw" in k):
                anndem[k] = anndem[k].multiply(sf["dw_el"]["sf"],
                                               axis=0,
                                               level=0)
            if ("nhe" in k) and ("nd" in k):
                anndem[k] = anndem[k].multiply(sf["nd_el"]["sf"],
                                               axis=0,
                                               level=0)
    
    # if required, calibrate annual demands at GB level based on ECUK data 
    elif snakemake.params.dic["scen_dem_calib"]=="GB" :
        

        raw_data = pd.read_csv(snakemake.input.path_set_nd_prop) 
        NIGDPfrac = raw_data["VALUE"][raw_data["VARIABLE"]=="NIGDPfraction"].values[0]
        NIPOPfrac = raw_data["VALUE"][raw_data["VARIABLE"]=="NIPOPfraction"].values[0]
        
        ECUK = pd.read_excel(snakemake.input.path_ecuk_enduse,
                              sheet_name="Table U2",
                              usecols=(0,1,3,4,7,10),skiprows=4,nrows=300)
        ECUK = ECUK.set_index(["Sector","End use","Year"]).sort_index()
        

        ECUK_data = ECUK.xs(["Domestic","Space heating"],
                            level=[0,1])["Total"].rename("DW_SH").to_frame()
        ECUK_data["DW_HW"] = ECUK.xs(["Domestic","Water heating"],
                                     level=[0,1])["Total"]
        ECUK_data["DW_NHG"] = ECUK.xs(["Domestic","Cooking/catering"],
                                     level=[0,1])["Natural gas"]
        
        ECUK_data["DW_NHE"] = (ECUK.loc[("Domestic",
                                        ("Cooking/catering",
                                         "Lighting and appliances"),
                                        slice(None)),
                                       "Electricity"]
                               .groupby("Year").sum())
        
        ECUK_data["ND_SH"] = (ECUK.loc[(("Services","Industry"),
                                        "Space heating",
                                        slice(None)),
                                       "Total"]
                               .groupby("Year").sum())  
        
        ECUK_data["ND_HW"] = (ECUK.loc[(("Services","Industry"),
                                        "Water heating",
                                        slice(None)),
                                       "Total"]
                               .groupby("Year").sum())  
        # no other NHG end use for Industry assumed to be captured demand
        # (thus this only includes service)
        ECUK_data["ND_NHG"] = (ECUK.loc[(("Services"),
                                        ("Cooking/catering",
                                         "Cooling and ventilation",
                                         "Other"),
                                        slice(None)),
                                       "Natural gas"]
                               .groupby("Year").sum())  
        
        ECUK_data["ND_NHE"] = (pd.concat([
                                        ECUK.loc[(("Services"),
                                        ("Cooking/catering",
                                         "Cooling and ventilation",
                                         "Computing",
                                         "Lighting",
                                         "Other"),
                                        slice(None)),
                                       "Electricity"],
                                        ECUK.loc[(("Industry"),
                                        ("Refrigeration",
                                         "Lighting"),
                                        slice(None)),
                                       "Electricity"]
                                        ])     
                               .groupby("Year").sum()) 
        
        # convert to GJ (from ktoe)
        ECUK_data = ECUK_data * 11.63 * 3600
        
        # apply NI factor to get GB data
        # fraction of NI/UK GDP/population (2015-2021 - 2022 for GDP not 
        # in dataset at the time of implementation)
        ECUK_data.loc[:,[c for c in
                         ECUK_data.columns
                         if "DW" in c]] = ECUK_data.loc[:,[c for c in
                                          ECUK_data.columns
                                          if "DW" in c]] * (1-NIPOPfrac)
                                                           
        ECUK_data.loc[:,[c for c in
                         ECUK_data.columns
                         if "ND" in c]] = ECUK_data.loc[:,[c for c in
                                          ECUK_data.columns
                                          if "ND" in c]] * (1-NIGDPfrac)
        # get relevant demands
        eff_data = pd.read_csv(snakemake.input.path_set_eff)
        hsys = pd.read_csv(snakemake.input.path_res_caps_ht_f, index_col=["LSOA11CD",
                                                               "PROPERTY_TYPE",
                                                               "TECHNOLOGY"])
        # add non-domestic detail, assuming average
        hsys = hsys.reset_index()
        sec = ";".join(['Others', 'Community, arts & leisure', 'Emergency services',
               'Education', 'Health', 'Storage', 'Industrial', 'Hospitality',
               'Offices', 'Retail'])
        hsys = hsys.replace({"Non-domestic":sec})
        hsys = hsys.assign(PROPERTY_TYPE=hsys['PROPERTY_TYPE'].str.split(';')).explode('PROPERTY_TYPE')
        hsys = hsys.set_index(["LSOA11CD","PROPERTY_TYPE","TECHNOLOGY"])
        
        
        ad = pd.DataFrame()
        for d in ECUK_data.columns:
            if "SH" in d:
                end = "_ex"
            else:
                end = "_ex_nb"
                
            astr = (d.split("_")[1]+"_"+d.split("_")[0]+end).lower()
            ys =[str(y) for y in ECUK_data[d].index.to_list()]
            ad[d] = anndem[astr].sum().loc[ys]
                                                  
            # applying efficiency to ECUK data to get from final to useful
            # energy

            if "DW" in d:
                oibo = "OIBODDDE00,OIBODDSD00,OIBODDTE00,OIBODDFL00"
                bmbo = "BMBODDDE00,BMBODDSD00,BMBODDTE00"
                hium = "HIUMDDDE00,HIUMDDSD00,HIUMDDTE00,HIUMDDFL00"
                ngbo = "NGBODDDE00,NGBODDSD00,NGBODDTE00,NGBODDFL00"
                ashp = "ASHPDDDE00,ASHPDDSD00,ASHPDDTE00,ASHPDDFL00"
                gshp = "GSHPDDDE00,GSHPDDSD00,GSHPDDTE00,GSHPDDFL00"
                elre = "ELREDDDE00,ELREDDSD00,ELREDDTE00,ELREDDFL00"
                elst = "ELSTDDDE00,ELSTDDSD00,ELSTDDTE00,ELSTDDFL00"
                
            elif "ND" in d:
                oibo = "OIBODNDO00"
                bmbo = "BMBODNDO00"
                hium = "HIUMDNDO00"
                ngbo = "NGBODNDO00"
                ashp = "ASHPDNDO00"
                gshp = "GSHPDNDO00"
                elre = "ELREDNDO00"
                elst = "ELSTDNDO00"
                
            if "SH" in d:
                mng = 1
                mot = 1
            elif"HW" in d:
                mng = 3
                mot = 2
                
            if "HW" in d or "SH" in d:
            # eff =   eff_data.loc[(eff_data["TECHNOLOGY"]
            #               ==("NGBODDDE00,NGBODDSD00,"+
            #                   "NGBODDTE00,NGBODDFL00"))&
            #             (eff_data["MODE_OF_OPERATION"]==3), "VALUE"].squeeze()
                eff = (anndem[astr].multiply(
                                (hsys.loc[hsys.index.get_level_values(2).str.contains("OIBO"),
                                          "fraction"].groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
                                    *eff_data.loc[(eff_data["TECHNOLOGY"]
                                                  ==oibo)&
                                                (eff_data["MODE_OF_OPERATION"]==mot), "VALUE"].values
                                                                        ).fillna(0).add(
                                hsys.loc[hsys.index.get_level_values(2).str.contains("BMBO"),
                                          "fraction"].groupby(["LSOA11CD","PROPERTY_TYPE"]).sum()
                                    *eff_data.loc[(eff_data["TECHNOLOGY"]
                                                  ==bmbo)&
                                                (eff_data["MODE_OF_OPERATION"]==mot), "VALUE"].values
                                                                        ).fillna(0).add(
                                hsys.loc[hsys.index.get_level_values(2).str.contains("HIUM"),
                                         "fraction"].groupby(["LSOA11CD",
                                                              "PROPERTY_TYPE"]).sum()
                                    *eff_data.loc[(eff_data["TECHNOLOGY"]
                                                  ==hium)&
                                                (eff_data["MODE_OF_OPERATION"]==mot), "VALUE"].values
                                                                        ).fillna(0).add(
                                hsys.loc[hsys.index.get_level_values(2).str.contains("NGBO"),
                                         "fraction"].groupby(["LSOA11CD",
                                                              "PROPERTY_TYPE"]).sum()
                                    *eff_data.loc[(eff_data["TECHNOLOGY"]
                                                  ==ngbo)&
                                                (eff_data["MODE_OF_OPERATION"]==mng), "VALUE"].values
                                                                        ).fillna(0).add(
                                    hsys.loc[hsys.index.get_level_values(2).str.contains("ASHP")
                                            |hsys.index.get_level_values(2).str.contains("AWHP")
                                            |hsys.index.get_level_values(2).str.contains("AAHP"),"fraction"].groupby(["LSOA11CD",
                                                                                                                      "PROPERTY_TYPE"]).sum()
                                            *eff_data.loc[(eff_data["TECHNOLOGY"]
                                                          ==ashp)&(
                                                          eff_data["YEAR"]
                                                                          =="2020")&
                                                        (eff_data["MODE_OF_OPERATION"]==mot), "VALUE"].values
                                        ).fillna(0).add(
                                    (hsys.loc[hsys.index.get_level_values(2).str.contains("GSHP"),"fraction"].droplevel(2)
                                          *eff_data.loc[(eff_data["TECHNOLOGY"]
                                                          ==gshp)&(
                                                          eff_data["YEAR"]
                                                                        =="2020")&
                                                      (eff_data["MODE_OF_OPERATION"]==mot), "VALUE"].values
                                              ).fillna(0),fill_value=0).add(
                                      (hsys.loc[hsys.index.get_level_values(2).str.contains("ELRE"),"fraction"].droplevel(2)
                                          *eff_data.loc[(eff_data["TECHNOLOGY"]
                                                          ==elre)&
                                                            (eff_data["MODE_OF_OPERATION"]==1),"VALUE"].values
                                              ).fillna(0),fill_value=0).add(
                                      (hsys.loc[hsys.index.get_level_values(2).str.contains("ELST"),"fraction"].droplevel(2)
                                          *eff_data.loc[(eff_data["TECHNOLOGY"]
                                                          ==elst)&
                                                            (eff_data["MODE_OF_OPERATION"]==mot),"VALUE"].values).fillna(0),fill_value=0)
                                      ,axis=0)
                                    ).fillna(0).sum()/anndem[astr].sum()
                eff = eff[ys]                           
                eff.index.name = "YEAR"
                eff.index = pd.to_numeric(eff.index,errors='coerce')                               
                                                  
            # elif ("SH" in d or "HW" in d) and "ND" in d:
            #     eff =   eff_data.loc[(eff_data["TECHNOLOGY"]
            #                   =="NGBODNDO00")&
            #                 (eff_data["MODE_OF_OPERATION"]==3), "VALUE"].squeeze()
            else:
                eff = 1
            print(d)
            print(eff)
            ECUK_data[d] = ECUK_data[d]* eff
            
        ad.index = pd.to_numeric(ad.index,errors='coerce')
           

        
        # calculate scaling factor, assume same for 2015 and 2016 as for 2017
        # assume average for future years
        sf = ECUK_data/ad
        sf.loc[2015] = sf.loc[2017]
        sf.loc[2016] = sf.loc[2017]
        for y in [int(y) for y in anndem["sh_dw_ex"].columns
                if int(y) not in sf.index]:
            sf.loc[y] = sf.loc[2017:2022].mean()
        
        sf.index = sf.index.astype(str)
        
        # scale demand
        for k in anndem.keys():
            f = (k.split("_")[1]+"_"+k.split("_")[0]).upper()
            anndem[k] = anndem[k].multiply(sf[f],axis=1)

                          

    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_dem_calib']}' is currently"
                " not a valid option of the 'scen_dem_calib' parameter."
                " Valid options are '-','GB', or 'MSOA'."
            )              
            
    
    
    # save to file
    # unit: GJ/a
    for i,d,b in demands:
        anndem[b].to_csv(d)
    

    
