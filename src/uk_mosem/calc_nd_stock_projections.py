"""
Script deriving non-domestic building stock projections


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

from zipfile import ZipFile

import pandas as pd

import utils


if __name__ == "__main__":
    
    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_nd_stock_projections",
                                         **df_params.iloc[0].to_dict())
        
    # load area lookup (LAD23CD)
    lsoas = utils.get_entity_lookup(["LSOA","LAD"])

    # load postcode  lookup
    pcd = utils.get_entity_lookup(levels=["PCD","LSOA"])
    
    # load DECs, and non-domestic EPCs and arrange data (incl. aggregation
    # with respect to LSOAs)

    # set columns from DECs/EPCs to be retained
    cols = ["POSTCODE","BUILDING_REFERENCE_NUMBER",
            "INSPECTION_DATE",
            #"OPERATIONAL_RATING_BAND",
            #"PROPERTY_TYPE",
            "PROPERTY_TYPE",
            "MAIN_HEATING_FUEL",
            "BUILDING_ENVIRONMENT",
            "FLOOR_AREA"]
        
     # set columns which will be grouped by, i.e., LSOAs and categories
    cols_g = ["LSOA11CD",
              #"OPERATIONAL_RATING_BAND",
              #"PROPERTY_TYPE",
              # "BUILDING_REFERENCE_NUMBER",
              "PROPERTY_TYPE",
              "MAIN_HEATING_FUEL",
              "BUILDING_ENVIRONMENT"]
    
    stock_c=dict()
    for zp in [snakemake.input.path_nd_stock_dec_ew,
               snakemake.input.path_nd_stock_epc_ew,
               snakemake.input.path_nd_stock_epc_s]:
        
        zip_file = ZipFile(zp)
        
        # skip first row for the Scottish EPC file
        sr=None
        if zp==snakemake.input.path_nd_stock_epc_s:
            sr=1
        
        for file in zip_file.infolist():
            if file.filename.endswith(("certificates.csv",
                                       "_0721.csv")):
                data = pd.read_csv(zip_file.open(file.filename),skiprows=sr,
                                   encoding="ISO-8859-1")
                
                # rename columns for them to match across zip files
                data = data.rename(columns={"ASSET_RATING_BAND":"OPERATIONAL_RATING_BAND",
                                            "TOTAL_FLOOR_AREA":"FLOOR_AREA",
                                            "CURRENT_ENERGY_PERFORMANCE_BAND":"OPERATIONAL_RATING_BAND"})
                
                # only use the latest certificate for each building, discard
                # the rest
                data = data[cols].sort_values("INSPECTION_DATE")
                data = data[~data[["BUILDING_REFERENCE_NUMBER"]].duplicated(keep="last")]
                
                # aggregate across LSOA while adding number of properties
                data["NUMBER_OF_PROPERTIES"] = 1
                # data = data.drop(["INSPECTION_DATE","BUILDING_REFERENCE_NUMBER"],
                #                  axis=1)
                s = pcd[["PCD","LSOA11CD"]].set_index('PCD')['LSOA11CD']
                data["POSTCODE"] = data["POSTCODE"].map(s)
                data = data.rename(columns={"POSTCODE":"LSOA11CD"})
                
                # aggregate
                data = data.groupby(cols_g+["BUILDING_REFERENCE_NUMBER",
                                            "INSPECTION_DATE"]).sum()
                
                if zp not in stock_c.keys():
                    stock_c[zp] = data
                else:
                    stock_c[zp] = pd.concat([stock_c[zp],data])
    
    # concatenate data and merge with LSOAs
    nd_stock = pd.concat(stock_c.values())
    nd_stock = nd_stock.reset_index()
    # drop duplicate
    nd_stock = nd_stock.sort_values("INSPECTION_DATE")
    nd_stock = nd_stock[~nd_stock[["BUILDING_REFERENCE_NUMBER"]].duplicated(keep="last")]
    nd_stock = nd_stock.drop(["INSPECTION_DATE","BUILDING_REFERENCE_NUMBER"],
                      axis=1)
    nd_stock = nd_stock.groupby(cols_g).sum()
    
    nd_stock = lsoas.merge(right=nd_stock.reset_index(), on="LSOA11CD", how="left")
    
    ### process property type information
    
    # split floor space and number of properties if more than one value given
    
    # remove ; in the end of strings
    nd_stock.loc[nd_stock["PROPERTY_TYPE"].notna()
                 & nd_stock["PROPERTY_TYPE"].str.endswith(";"),
                 "PROPERTY_TYPE"] = nd_stock.loc[nd_stock["PROPERTY_TYPE"].notna() 
                                                 & nd_stock["PROPERTY_TYPE"].str.endswith(";"),
                                                 "PROPERTY_TYPE"].str[:-1]                                           
    nd_stock = nd_stock.assign(PROPERTY_TYPE=nd_stock["PROPERTY_TYPE"].str.split(";"))
    
    nd_stock["NUMBER_OF_PROPERTIES"] = [n/len(pt) if isinstance(pt,list) else
                                        n for n, pt in zip(nd_stock["NUMBER_OF_PROPERTIES"],
                                               nd_stock["PROPERTY_TYPE"])]
    nd_stock["FLOOR_AREA"] = [f/len(pt) if isinstance(pt,list) else
                              f for f, pt in zip(nd_stock["FLOOR_AREA"],
                                               nd_stock["PROPERTY_TYPE"])]
    nd_stock = nd_stock.explode("PROPERTY_TYPE")
    
    # fill nan with default values (PROPERTY_TYPE, MAIN_HEATING_FUEL is given 
    # whenever FLOOR_AREA is - in current data set-, so no assumption is/needs 
    # to be applied and below is just to fill values where FLOOR_AREA is also nan)
    nd_stock["PROPERTY_TYPE"] = nd_stock["PROPERTY_TYPE"].fillna("Others - Miscellaneous 24hr activities")
    nd_stock["MAIN_HEATING_FUEL"] = nd_stock["MAIN_HEATING_FUEL"].fillna("Other")
    nd_stock["BUILDING_ENVIRONMENT"] = nd_stock["BUILDING_ENVIRONMENT"].fillna("NONE")
    nd_stock["FLOOR_AREA"] = nd_stock["FLOOR_AREA"].fillna(0)
    nd_stock["NUMBER_OF_PROPERTIES"] = nd_stock["NUMBER_OF_PROPERTIES"].fillna(0)
    
    # aggregate property types
    ptlu = pd.read_csv(snakemake.input.path_nd_pt_lu,index_col="PROPERTY_TYPE",
                        usecols=[0,1])
    nd_stock["PROPERTY_TYPE"] = nd_stock["PROPERTY_TYPE"].map(ptlu["MATCH"].to_dict())    
    nd_stock = nd_stock.groupby(cols_g).sum()
    
    # add transport to others, remove domestic and other types not covered
    nd_stock = nd_stock.rename({"Transport":"Others"})
    nd_stock = nd_stock.drop("Domestic",level="PROPERTY_TYPE")
    
    nd_stock = nd_stock.reset_index()
    
    
    
    # process heating system data
    nd_stock = nd_stock.replace({"Natural Gas":"NGBO",
                         "Anthracite":"BMBO",
                         "Smokeless Fuel (inc Coke)":"BMBO",
                         "Grid Displaced Electricity":"ELRE",
                         "Grid Supplied Electricity":"ELRE",
                         "Dual Fuel Appliances (Mineral + Wood)":"BMBO",
                         "Biogas":"BMBO",
                         "LPG":"OIBO",
                         "Oil":"OIBO",
                         "District Heating":"HIUM",
                         "Waste Heat":"HIUM",
                         "Coal":"BMBO",
                         "Other":"ELRE",
                         "Biomass":"BMBO"})

    nd_stock = nd_stock.rename(columns={"MAIN_HEATING_FUEL":"MAINHEAT_DESCRIPTION",
                                        "BUILDING_ENVIRONMENT":"HEATDIST_DESCRIPTION"}) 
    
    nd_stock.loc[(nd_stock["MAINHEAT_DESCRIPTION"]=="ELRE")&(
                 nd_stock["HEATDIST_DESCRIPTION"] =="Air Conditioning"),
                "MAINHEAT_DESCRIPTION"] = "ASHP"
    # set all heat distributions systems to low temperature wet if a HP is
    # installed
    nd_stock.loc[nd_stock["MAINHEAT_DESCRIPTION"]=="ELRE",
                 "HEATDIST_DESCRIPTION"] = "NONE"
    nd_stock.loc[nd_stock["MAINHEAT_DESCRIPTION"]!="ELRE",
                 "HEATDIST_DESCRIPTION"] = "WDIS"   
    nd_stock.loc[nd_stock["MAINHEAT_DESCRIPTION"]=="ASHP",
                 "HEATDIST_DESCRIPTION"] = "RAUP"
    
    nd_stock = nd_stock.groupby(["LSOA11CD",
                                 "PROPERTY_TYPE",
                                 "HEATDIST_DESCRIPTION",
                                 "MAINHEAT_DESCRIPTION"],dropna=False).sum()
    nd_stock = nd_stock.reset_index()
   
    # load VOA data and overwrite EPC where appropriate VOA data exist
    # (Retail, offices, industryŝ)
    
    # load property numbers from VOA data for England & Wales, single Scotland
    # data point from data sheet
    
    zip_file = ZipFile(snakemake.input.path_nd_stock_ndr)
    stock_EW_T = pd.read_csv(zip_file.open("Table FS_OA1.1.csv"))
    stock_EW_R = pd.read_csv(zip_file.open("Table FS_OA2.1.csv"))
    stock_EW_O = pd.read_csv(zip_file.open("Table FS_OA3.1.csv"))
    stock_EW_I = pd.read_csv(zip_file.open("Table FS_OA4.1.csv"))
    stock_EW_R["PROPERTY_TYPE"] = "Retail"
    stock_EW_O["PROPERTY_TYPE"] = "Offices"
    stock_EW_I["PROPERTY_TYPE"] = "Industrial"
    stock_EW_VOA = pd.concat([stock_EW_R,stock_EW_O,stock_EW_I])
    stock_EW_VOA = stock_EW_VOA.rename(columns={"ons_code":"LSOA11CD"})
    stock_EW_VOA = stock_EW_VOA[["LSOA11CD","PROPERTY_TYPE","2014",
                                 "2015","2016","2017","2018","2019","2020",
                                 "2021","2022"]]
    stock_EW_VOA = stock_EW_VOA.set_index(["LSOA11CD","PROPERTY_TYPE"])
    stock_EW_VOA = stock_EW_VOA.apply(pd.to_numeric,errors='coerce').fillna(0).astype(int)
    # convert to m2 from thousands m2
    stock_EW_VOA = stock_EW_VOA * 1000
    
    # merge on EPC stock data and overwrite if appropriate
    # 2014 is chosen here as it aligns with BEES data used below
    nd_stock = nd_stock.merge(right=nd_stock.groupby(["LSOA11CD",
                                                      "PROPERTY_TYPE"]).sum()["FLOOR_AREA"],
                              on=["LSOA11CD","PROPERTY_TYPE"],
                              how="left",suffixes=("","_total"))
    stock_EW_VOA = lsoas.merge(right=stock_EW_VOA.reset_index(),
                                on="LSOA11CD", how="left")
    stock_EW_VOA = stock_EW_VOA.drop(columns=["LAD23CD","LAD23NM"])
    
    # remove Scotland
    stock_EW_VOA = stock_EW_VOA[~stock_EW_VOA["2014"].isna()]
    
    nd_stock = nd_stock.merge(right=stock_EW_VOA[["LSOA11CD",
                                                   "PROPERTY_TYPE",
                                                   "2014"]],
                              on=["LSOA11CD","PROPERTY_TYPE"], how="outer")
    
    # assume gas boiler as default heating
    nd_stock["MAINHEAT_DESCRIPTION"] = nd_stock["MAINHEAT_DESCRIPTION"].fillna("NGBO")
    nd_stock["HEATDIST_DESCRIPTION"] = nd_stock["HEATDIST_DESCRIPTION"].fillna("WDIS")
    # if no data from VOA (v) keep EPC data (f)
    # if no EPC data (f), take VOA data
    # if both data available, scale EPC data (which has heating, etc. detail)
    nd_stock["FLOOR_AREA"] = [f if pd.isna(v) else
                              v if pd.isna(f) else
                              f if ((f<1000)&(v==0)) else
                              f*v/t
                              for f,t,v in zip(nd_stock["FLOOR_AREA"],
                                               nd_stock["FLOOR_AREA_total"],
                                               nd_stock["2014"])]
    nd_stock = nd_stock.drop(columns=["FLOOR_AREA_total",
                                     "2014",
                                     "NUMBER_OF_PROPERTIES"])

    # load overall floor area data for Scotland and England & Wales
    raw_data = pd.read_csv(snakemake.input.path_set_nd_prop)    
    stock_S = raw_data["VALUE"][raw_data["VARIABLE"]=="NonDomesticFloorareaScotland"].values[0]
    stock_S_f = raw_data[["TECHNOLOGY","VALUE"]][raw_data["VARIABLE"]=="NonDomesticFractionScotland"]
    stock_S_f = stock_S_f.rename(columns={"TECHNOLOGY":"PROPERTY_TYPE"})
    stock_S_f = stock_S_f.set_index("PROPERTY_TYPE")
    stock_S_f["VALUE"] = stock_S_f["VALUE"]/stock_S_f["VALUE"].sum()
    
    stock_S = stock_S_f * stock_S
    
    stock_EW = pd.read_excel(snakemake.input.path_nd_stock_bees_app,
                             sheet_name="Table 3.1",
                             usecols=(1,14),skiprows=9,nrows=52)
    stock_EW.columns = ["PROPERTY_TYPE","VALUE"]
    stock_EW["PROPERTY_TYPE"] = stock_EW["PROPERTY_TYPE"].str.strip()
    stock_EW = stock_EW.replace({"Sport":"Community, arts & leisure",
                                 "Leisure":"Community, arts & leisure",
                                 "Community":"Community, arts & leisure",
                                 "Military":"Others",
                                 "Emergency Services":"Emergency services",
                                 "Factory":"Industrial",
                                 "Warehouse":"Storage"})
    stock_EW = stock_EW.dropna()
    stock_EW = stock_EW.groupby("PROPERTY_TYPE").sum()
    
    # scale EPC/DEC data to match the overall national floor area data
    nd_stock_S = nd_stock[nd_stock["LSOA11CD"].str.startswith("S")]
    nd_stock_S = nd_stock_S.merge(right=nd_stock_S.groupby("PROPERTY_TYPE").sum()["FLOOR_AREA"],
                              on="PROPERTY_TYPE", how="left",suffixes=("","_frac"))
    nd_stock_S = nd_stock_S.merge(right=stock_S,on="PROPERTY_TYPE", how="left",suffixes=("","_frac"))
    nd_stock_S["FLOOR_AREA_frac"] = (nd_stock_S["FLOOR_AREA_frac"]
                                     /nd_stock_S["VALUE"])
    nd_stock_S["FLOOR_AREA"] = (nd_stock_S["FLOOR_AREA"]
                                /nd_stock_S["FLOOR_AREA_frac"])
    nd_stock_EW = nd_stock[~nd_stock["LSOA11CD"].str.startswith("S")]
    nd_stock_EW = nd_stock_EW.merge(right=nd_stock_EW.groupby("PROPERTY_TYPE").sum()["FLOOR_AREA"],
                               on="PROPERTY_TYPE", how="left",suffixes=("","_frac"))
    nd_stock_EW = nd_stock_EW.merge(right=stock_EW,on="PROPERTY_TYPE", how="left",suffixes=("","_frac"))
    nd_stock_EW["FLOOR_AREA_frac"] = (nd_stock_EW["FLOOR_AREA_frac"]
                                      /nd_stock_EW["VALUE"])
    nd_stock_EW["FLOOR_AREA"] = (nd_stock_EW["FLOOR_AREA"]
                                 /nd_stock_EW["FLOOR_AREA_frac"])   

    nd_stock = pd.concat([nd_stock_S,nd_stock_EW])
    nd_stock = nd_stock.drop(columns=["FLOOR_AREA_frac",
                                      "VALUE"])

    # extend to all historic years based on average change in EW

    #nd_stock.columns = ["2014"]
    stock_EW_T = stock_EW_T.apply(pd.to_numeric,errors='coerce').fillna(0).astype(int)
    stock_ch =  stock_EW_T[["2014","2015","2016","2017","2018",
                             "2019","2020","2021","2022"]].sum().to_frame()
    stock_ch["EW"] = stock_ch[0]/stock_ch.loc["2014",0]
    stock_ch["S"] = stock_ch[0]/stock_ch.loc["2017",0]
    
    for y in range(2015,2023):
        nd_stock.loc[~nd_stock["LSOA11CD"].str.startswith("S"),str(y)] = nd_stock.loc[~nd_stock["LSOA11CD"].str.startswith("S"),"FLOOR_AREA"] * stock_ch.loc[str(y),"EW"]   
        nd_stock.loc[nd_stock["LSOA11CD"].str.startswith("S"),str(y)] = nd_stock.loc[nd_stock["LSOA11CD"].str.startswith("S"),"FLOOR_AREA"] * stock_ch.loc[str(y),"S"]

    # rearrange data
    nd_stock = nd_stock.drop(columns="FLOOR_AREA")
    nd_stock = nd_stock.set_index(["LSOA11CD",
                                   "PROPERTY_TYPE",
                                   "HEATDIST_DESCRIPTION",
                                   "MAINHEAT_DESCRIPTION"])
    nd_stock.columns.name = "YEAR"
    nd_stock = nd_stock.stack().reset_index()
    nd_stock = nd_stock.rename(columns={0:"FLOOR_AREA"})

    # process dataframe for residual heating system capacities

    nd_stock["MAINHEAT_DESCRIPTION"]= "MH_" + nd_stock["MAINHEAT_DESCRIPTION"].astype(str)
    nd_stock["HEATDIST_DESCRIPTION"]= "HD_" + nd_stock["HEATDIST_DESCRIPTION"].astype(str)

    # rearrange data in way that lists characteristics separately (while
    # losing information/detail on combinations)
    det = pd.concat([nd_stock.groupby(["LSOA11CD",
                                       "PROPERTY_TYPE",
                                       "YEAR"]+[d]).sum() for d in ["MAINHEAT_DESCRIPTION",
                                                  "HEATDIST_DESCRIPTION"]])
    det.index = det.index.set_names("HEATING_SYSTEM",level=3)
    

    # save to file
    det.to_csv(snakemake.output.path_nd_stock_hde)  
    
    nd_stock = nd_stock.set_index(["LSOA11CD", "PROPERTY_TYPE",
                                   "YEAR","MAINHEAT_DESCRIPTION",
                                   "HEATDIST_DESCRIPTION"])["FLOOR_AREA"].unstack("YEAR")    

    nd_stock = nd_stock.groupby(["LSOA11CD", "PROPERTY_TYPE"]).sum()
    

    # load household projections (percentage change)
    hh_proj_pch_ls = pd.read_csv(snakemake.input.path_hh_proj, index_col="LSOA11CD")
    
    # calculate stock for future years, separately for existing and
    # new    

    nnd_stock = nd_stock.copy(deep=True)
    nnd_stock.loc[:,"2015":"2022"] = 0
    
    d_rate = raw_data["VALUE"][raw_data["VARIABLE"]=="NonDomesticDemolitionRate"].values 
    
    
    for y in range(2023,2061):
        
        nd_stock[str(y)] = (nd_stock[str(y-1)]
                            - (nnd_stock[str(y-1)] + nd_stock[str(y-1)])
                            *d_rate)
        nnd_stock[str(y)] = ((nnd_stock[str(y-1)] + nd_stock[str(y-1)])*
                             (1+hh_proj_pch_ls[str(y)])-nd_stock[str(y)])

    ### calibrate based on ECUK data on space heat - unused using factor
    ### taking into account scope of BEES instead (see below), calibration with
    ### ECUK is done for annual demands
    # NIfrac = raw_data["VALUE"][raw_data["VARIABLE"]=="NIGDPfraction"].values[0]
    
    # ECUK = pd.read_excel(snakemake.input.path_ecuk_enduse,
    #                      sheet_name="Table U1",
    #                      usecols=(0,1,3,4,5,6),skiprows=5,nrows=128)
    # ECUK = ECUK.set_index(["End use","Year"]).sort_index()

    # hist_dem_SH = ECUK.loc[("Space heating",
    #                         slice("2015","2022")),
    #                        slice("Industrial","Service")]
    # # sum up and convert to GWh (from ktoe)
    # hist_dem_SH =  hist_dem_SH.sum().sum() * 11.63
    
    # # reduce total by fraction of NI/UK GDP (2015-2021 - 2022 not available at
    # # the time of implementation)
    # hist_dem_SH = hist_dem_SH * (1-NIfrac)
    # demint = pd.read_csv(snakemake.input.path_nd_int_sh,
    #                      index_col=["LSOA11CD","PROPERTY_TYPE"])
    # eff_data = pd.read_csv(snakemake.input.path_set_eff)
    
    # eff = eff_data.loc[(eff_data["TECHNOLOGY"]
    #               =="NGBODNDO00")&
    #             (eff_data["MODE_OF_OPERATION"]==1),"VALUE"].squeeze()
    
    # anndem_SH = nd_stock.multiply(demint) / 10**6 /eff
    # anndem_SH = anndem_SH.sum()["2015":"2022"].sum()
    # f = hist_dem_SH/anndem_SH
    
    # load factor to uplift floor space based on scope of BEES
    f = raw_data["VALUE"][raw_data["VARIABLE"]=="NonDomesticFloorareaUpliftFactor"].squeeze()  
    
    # scale
    nd_stock = f*nd_stock
    nnd_stock = f*nnd_stock
    
    # save results
    nd_stock.to_csv(snakemake.output.path_nd_stock_ex)
    nnd_stock.to_csv(snakemake.output.path_nd_stock_nb)

    
