"""
Script deriving demand intensities for the non-domestic building stock


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

import pandas as pd
import numpy as np

import utils



if __name__ == "__main__":
    
    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_nd_demand_intensities",
                                         **df_params.iloc[4].to_dict())
    
    # load year set
    sets = pd.read_csv(snakemake.input.path_set_sets)
    years = sets.loc[(sets["SET"]=="YEAR"),
                     "VALUE"].to_frame()   
    
    # load list of areas/area lookup
    lsoas = utils.get_entity_lookup(["LSOA"])
    
    # load raw data set
    dem_raw_data = pd.read_csv(snakemake.input.path_set_nd_dem_int)
    efrac = dem_raw_data.loc[dem_raw_data["VARIABLE"]=="NonDomesticElecCookingFraction","VALUE"].values
    eff_data = pd.read_csv(snakemake.input.path_set_eff)
    
    eff = eff_data.loc[(eff_data["TECHNOLOGY"]
                  =="NGBODNDO00")&
                (eff_data["MODE_OF_OPERATION"]==1),"VALUE"].values
    
    # load BEES intensity data
    BEES_int = pd.read_excel(snakemake.input.path_nd_stock_bees_ova,
                             sheet_name="Figure 3.11",skiprows=4,nrows=14)
    BEES_int = BEES_int.dropna()
    
    BEES_int = BEES_int.rename(columns={"Unnamed: 0":"Enduse",
                                        "All sectors":"Others"})
    BEES_int.columns.name="PROPERTY_TYPE"
    BEES_int = BEES_int.set_index("Enduse")
    
    # split catering/cooking in gas and electric
    BEES_int.loc["Catering_EL"] = BEES_int.loc["Catering"]* efrac
    BEES_int.loc["Catering_GA"] = BEES_int.loc["Catering"]* (1-efrac)
    BEES_int = BEES_int.drop("Catering")
    BEES_int = BEES_int.drop(columns="Military")
    BEES_int = BEES_int.rename(index={"Heating":"SH",
                                      "Hot water":"HW",
                                      "Cooling & humidification ":"NHE",
                                      "Cooled storage":"NHE",
                                      "ICT equipment":"NHE",
                                      "Small power":"NHE",
                                      "Fans":"NHE",
                                      "Lighting":"NHE",
                                      "Other":"NHE",
                                      "Catering_EL":"NHE",
                                      "Catering_GA":"NHG"})
    BEES_int = BEES_int.groupby("Enduse").sum()
    
    # get effiencies to get to useful demand for heating
    
    # load non-dom properties for base years
    nd_stock =  pd.read_csv(snakemake.input.path_nd_stock_ex,
                              index_col=["LSOA11CD","PROPERTY_TYPE"]) 
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
    
    oibo = "OIBODNDO00"
    bmbo = "BMBODNDO00"
    hium = "HIUMDNDO00"
    ngbo = "NGBODNDO00"
    ashp = "ASHPDNDO00"
    gshp = "GSHPDNDO00"
    elre = "ELREDNDO00"
    elst = "ELSTDNDO00"
    
    eff = dict()
    for d in ["SH","HW"]:
        if d=="SH":
            mng = 1
            mot = 1
        elif d=="HW":
            mng = 3
            mot = 2

        eff[d] = (nd_stock.multiply(
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
                            ).fillna(0).sum()/nd_stock.sum()
        eff[d] = eff[d].mean()
        
    # convert to useful energy for space and water heating
    BEES_int.loc["HW"] = BEES_int.loc["HW"] * eff["HW"]
    BEES_int.loc["SH"] = BEES_int.loc["SH"] * eff["SH"]
    # derive space heat intensity for new builts as fraction of intensity of
    # existing stock
    BEES_int.loc["SH_NB"] = (BEES_int.loc["SH"]
                                   * dem_raw_data.loc[dem_raw_data["VARIABLE"]=="DemandFracNewbuiltND","VALUE"].values)  
    
    # rearrange and apply across all LSOAs
    nd_int = lsoas.merge(right=BEES_int.T.reset_index(),how="cross")
    nd_int = nd_int.set_index(["LSOA11CD","PROPERTY_TYPE"])
    nd_int.columns.name="Enduse"
    nd_int = nd_int.stack().to_frame()
    nd_int.columns=["2015"]
    
    # extend to all years
    for y in range(2016,2061):
        nd_int[str(y)] = nd_int["2015"]
        
    nd_intc = nd_int.copy()
    
    # overwrite NHE intensity with trajectory based on run
    if snakemake.params.dic["scen_nhe"] != "-":
        dev = years.copy()
        dev = dev.rename(columns={"VALUE":"YEAR"})
        dev["VALUE"] = np.nan
        vys = snakemake.params.dic["scen_nhe"].split("a")
        
        for  vy in vys:
            value, year = vy.split("y")
            dev.loc[dev["YEAR"]==year,"VALUE"] = value
            
        dev["VALUE"] = dev["VALUE"].astype(float)
        
        # interpolate
        dev["VALUE"] = dev["VALUE"].interpolate(limit_direction="both")
        
        for i in dev["YEAR"]:
            nd_int.loc[(slice(None),slice(None),"NHE"),str(i)] = (
                nd_int.loc[(slice(None),slice(None),"NHE"),str(i)]
                * dev.loc[dev["YEAR"]==i,"VALUE"].squeeze())
            
    # overwrite NHG intensity with trajectory based on run
    if snakemake.params.dic["scen_nhg"] != "-":
        dev = years.copy()
        dev = dev.rename(columns={"VALUE":"YEAR"})
        dev["VALUE"] = np.nan
        vys = snakemake.params.dic["scen_nhg"].split("a")
        
        for  vy in vys:
            value, year = vy.split("y")
            dev.loc[dev["YEAR"]==year,"VALUE"] = value
            
        dev["VALUE"] = dev["VALUE"].astype(float)
        
        # interpolate
        dev["VALUE"] = dev["VALUE"].interpolate(limit_direction="both")
        
        for i in dev["YEAR"]:
            nd_int.loc[(slice(None),slice(None),"NHG"),str(i)] = (
                nd_int.loc[(slice(None),slice(None),"NHG"),str(i)]
                * dev.loc[dev["YEAR"]==i,"VALUE"].squeeze())
            
    if ((snakemake.params.dic["scen_local_gov"] != "-") and
        (snakemake.params.dic["scen_local_gov"] != "A")):
        lsoas = utils.get_entity_lookup(["LSOA","LAD"])
        et = pd.read_csv(snakemake.input.path_loc_gov_em)
        sets = pd.read_csv(snakemake.input.path_set_sets)
        years = sets.loc[(sets["SET"]=="YEAR"),
                         "VALUE"].to_frame()
        years["VALUE"] = years["VALUE"].astype(int)
        
        et = et[["REGION","YEAR","VALUE"]].set_index(["REGION","YEAR"]).unstack().droplevel(0,axis=1)
        for y in years["VALUE"]:
            if y not in et.columns:
                et.loc[:,y] = np.nan
        et = et.sort_values(by="YEAR",axis=1)
        
        et.loc[:,2015:2022]=1
        et = et.interpolate(axis=1)
        et = lsoas.merge(right=et.reset_index(),left_on="LAD23CD",
                         right_on="REGION",how="left").set_index("LSOA11CD")
        et = et.drop(["LAD23CD","REGION","LAD23NM"],axis=1)
        if snakemake.params.dic["scen_nhg"] != "-":
            for y in dev["YEAR"]:
                et.loc[:,int(y)] = et.loc[:,int(y)].fillna(dev.loc[dev["YEAR"]==y,"VALUE"].squeeze())
        else:
            et = et.fillna(1)
        
        et.columns=et.columns.astype(str)

        nd_int.loc[(slice(None),slice(None),"NHG"),slice(None)] =(
        nd_intc.loc[(slice(None),slice(None),"NHG"),slice(None)]*et)


    # save to files
    # unit: kWh/a/m^2
    nd_int.xs("SH",level=2).to_csv(snakemake.output.path_nd_int_sh)
    nd_int.xs("SH_NB",level=2).to_csv(snakemake.output.path_nd_int_sh_nb)
    nd_int.xs("HW",level=2).to_csv(snakemake.output.path_nd_int_hw)
    nd_int.xs("NHG",level=2).to_csv(snakemake.output.path_nd_int_nhg)
    nd_int.xs("NHE",level=2).to_csv(snakemake.output.path_nd_int_nhe)

