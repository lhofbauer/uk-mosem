"""
Script deriving techno-economic parameters for supply technologies, including
transmission


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

import pandas as pd
import numpy as np
import geopandas as gpd

import utils


def create_detailed_supply():
    
    # set up fossil fuel supply
    
    cost_data = pd.read_excel(snakemake.input.path_supply_prices,
                              sheet_name="CP1")# nrows=42,header=0,usecols="C:AU"
    cost_data = cost_data.loc[~cost_data.iloc[:,[1,11]].isna().all(axis=1)]
    
    cost = cost_data.iloc[[1,2,7,11],11:]
    cost = cost.set_index("Unnamed: 11")
    cost.columns = cost.iloc[0].astype(int)
    cost = cost.iloc[1:]
    
    cost = cost.rename(index={"Gas price": "NGSUSNAT00",
                              "Oil price": "LFSUSNAT00",
                              "Coal price":"CLSUSNAT00"})
    cost.index.name = "TECHNOLOGY"
    cost.columns.name = "YEAR"
    
    # convert to 2015GBP/MJ (= 10^6 GBP/TJ)
    cost = cost/cost_data.iloc[33,2]*cost_data.iloc[26,10]/3600
    
    
    ### currently oil is adjusted to represent heating oil prices
    
    # load heating oil prices
    fcost_data = pd.read_excel(snakemake.input.path_supply_fprices,
                              sheet_name="4.1.2 (excl VAT)",
                              skiprows=5,
                              usecols=[0,13])

    fcost_data.columns=["YEAR","VALUE"]
    fcost_data = fcost_data.groupby("YEAR").mean()
    fcost_data = fcost_data.drop(2024)
    fcost_data = fcost_data.reset_index()
    data = pd.read_csv(snakemake.input.path_set_supply) 
    d = data.loc[(data["VARIABLE"]=="HeatingOilEnergyDensity"),
                 "VALUE"].squeeze()
    # convert to XGBP/MJ (= 10^6 GBP/TJ)
    fcost_data["VALUE"] = fcost_data["VALUE"]/100/d/1000
    
    # convert to 2015GBP/MJ (= 10^6 GBP/TJ)
    fcost_data["UNIT"] = fcost_data["YEAR"].astype(str)+"£/MJ"
    fcost_data = utils.adjust_monetary_values(fcost_data, 2015,deflator="GDP")
    fcost_data= fcost_data.drop("UNIT",axis=1)
    fcost_data = fcost_data.set_index("YEAR")
    proj = cost.loc["LFSUSNAT00"]/cost.loc["LFSUSNAT00",2023]
    proj.name="VALUE"    
    fcost_data = pd.concat([fcost_data,
                            proj.loc[2024:].to_frame()])
    fcost_data.loc[2024:2050] =  fcost_data.loc[2024:2050] *fcost_data.loc[2023]
    
    # overwrite cost
    cost.loc["LFSUSNAT00",:] = fcost_data.loc[2010:,"VALUE"]    
    
    # create and rearrange dataframe for cost 
    ffvarcost = cost.stack()
    ffvarcost.name = "VALUE"
    ffvarcost = ffvarcost.reset_index()
    ffvarcost["MODE_OF_OPERATION"] = "1"
    ffvarcost["UNIT"] = "million 2015£/TJ"
    ffvarcost["VARIABLE"] = "VariableCost"
    ffvarcost["REGION"] = ":*0"


    # create and save efficiencies dataframe
    ffeff = pd.DataFrame(["NGSUSNAT00",
                          "LFSUSNAT00",
                          "CLSUSNAT00"],columns=["TECHNOLOGY"])

    ffeff["VALUE"]=1
    ffeff["MODE_OF_OPERATION"] = "1"
    ffeff["UNIT"] = "-"
    ffeff["VARIABLE"] = "Efficiency"
    ffeff["REGION"] = ":*0"
    ffeff["YEAR"] = ":*"
    ffeff["FUEL_OUT"] = ["NGSTRA","LFSTRA","CLSTRA"]
    ffeff["FUEL_IN"] = np.nan
  
  
    
    elgparams = pd.read_excel(snakemake.input.path_supply_el,
                              sheet_name="Technical and cost assumptions",
                              header=[0,1],
                              index_col=[0,1])
    elgparams = elgparams.iloc[:,1:]
    elgparams.index.names = ["PARAMETER","SUBCATEGORY"]
    elgparams.columns.names = ["TECHNOLOGY","YEAR"]
   
    elgparams = pd.concat([elgparams,
                           elgparams[["CCGT H Class"]].rename(
                               columns={"CCGT H Class": "H2 CCGT"})], axis=1)
    
    eltechs = {"CCGT H Class":"NGPPSNAT00",
               "H2 CCGT":"H2PPSNAT00",
               "Hydro 516MW":"HYPPSNAT00",
              "Dedicated Biomass":"BMPPSNAT00",
              "Onshore Wind": "OWPPSNAT00",
              "Offshore Wind": "FWPPSNAT00",
              "Large-scale Solar":"OSPPSNAT00",
              "Solar PV 10-50kW":"RSPPSNAT00"
              }
    eltechs_infuel = {"NGPPSNAT00":"NGSTRA",
                      "H2PPSNAT00":"H0SSTO",
                      "HYPPSNAT00":np.nan,
                      "BMPPSNAT00":"BMSTRA",
                      "OWPPSNAT00":np.nan,
                      "FWPPSNAT00":np.nan,
                      "OSPPSNAT00":np.nan,
                      "RSPPSNAT00":np.nan
                      }
    eltechs_outfuel = {"NGPPSNAT00":"E1STRA",
                      "H2PPSNAT00":"E0STRA",
                      "HYPPSNAT00":"E0STRA",
                      "BMPPSNAT00":"E0STRA",
                      "OWPPSNAT00":"E0STRA",
                      "FWPPSNAT00":"E0STRA",
                      "OSPPSNAT00":"E0STRA",
                      "RSPPSNAT00":"E0STRA"
                      }
    
    
    hygparams = pd.read_excel(snakemake.input.path_supply_h2,
                              sheet_name="Technical and cost assumptions",
                              skiprows=1,
                              header=[0,1])
    hygparams.iloc[:,0] = hygparams.iloc[:,0].fillna(method='ffill')
    hygparams = hygparams.set_index([c for c in hygparams.columns[:3]])
    
    hygparams = hygparams.droplevel(1)
    hygparams.index.names = ["PARAMETER","SUBCATEGORY"]
    hygparams.columns.names = ["TECHNOLOGY","YEAR"]
    hygparams = pd.concat([hygparams, hygparams[["Alkaline "]].rename(
                columns={"Alkaline ": "Alkaline_storage"})], axis=1)
    
    hytechs = {"Alkaline ":"AEHPSNAT00",
               "Alkaline_storage":"AESPSNAT00",
              "Steam Methane Reformation with CCUS":"SRHPSNAT00"
              }
    hytechs_infuel = {"AEHPSNAT00":"E0STRA",
                      "AESPSNAT00":"E0STRA",
                      "SRHPSNAT00":"NGSTRA"
                      }
    hytechs_outfuel = {"AEHPSNAT00":"H0STRA",
                       "AESPSNAT00":"H0SSTO",
                      "SRHPSNAT00":"H1STRA"
                      }
    
    # calculate model params

    elgparams = elgparams[(elgparams.index.get_level_values(1)!="Low")&
                        (elgparams.index.get_level_values(1)!="High")]
    elgparams = elgparams.droplevel(1)
     
    elgparams.loc[["Reference plant size"],:] = elgparams.loc[["Reference plant size"],:].replace(0,np.nan)
    
    hygparams = hygparams[(hygparams.index.get_level_values(1)!="Low")&
                        (hygparams.index.get_level_values(1)!="High")]
    hygparams = hygparams.droplevel(1)
    hygparams = hygparams.loc[:,~hygparams.columns.get_level_values(1).to_series().str.contains(".").fillna(False).values]
    hygparams.loc[["Reference plant size"],:] = hygparams.loc[["Reference plant size"],:].replace(0,np.nan)
    

    # costs in GBP per kW or GBP per kWh
    
    # recalculated connection charges (assuming most is connection), cost
    # are small for all but offshore wind anyway, as investment cost in GBP/kW
    data = pd.read_csv(snakemake.input.path_set_supply) 
    r = data.loc[(data["VARIABLE"]=="RoI"),"VALUE"].squeeze()
    ccs = data.loc[(data["VARIABLE"]=="CCS_charge"),"VALUE"].squeeze()
    elgparams.loc[["Connection and use of system charges"],
                  :] = (elgparams.loc[["Connection and use of system charges"],
                                     :].replace(0,np.nan)
                        * [(1-(1/(1+r)**l))/r for l in elgparams.loc["Operating lifetime",:]]
                        /1000)
    elgparams.loc[["Connection and use of system charges"],
                  :] = elgparams.loc[["Connection and use of system charges"],
                                :].fillna(0)
    elgparams.loc["Infrastructure"] = (1000 
                                      * elgparams.loc["Infrastructure"].fillna(0)
                                      /elgparams.loc["Reference plant size",:]
                                      /1000)
    elgparams.loc["CapitalCost"] =  elgparams.loc[["Infrastructure",
                                                "Pre-development cost",
                                                "Construction cost",
                                                "Connection and use of system charges"]].sum()
    elgparams.loc["VariableCost"] =  elgparams.loc["Variable O&M"]/1000/3.6

    elgparams.loc["FixedCost"] =  elgparams.loc[["Fixed O&M",
                                                 "Insurance"]].fillna(0).sum()/1000                                          
    elgparams.loc["Efficiency"] =  elgparams.loc["Average fuel efficiency"].replace("-",1)                                          
    elgparams.loc["OperationalLife"] =  elgparams.loc["Operating lifetime"]
    
    # for flexible technologies this is equivalent to the availability factor,
    # for other techs it incorporates capacity factor, etc. Thus, this is not
    # used for wind and solar technologies that have a separate
    # capacity factor applied but is used to incorporate capacity factor
    # element for hydro
    elgparams.loc["AvailabilityFactor"] = elgparams.loc["Average load factor (net of availablility)"]
    
    hygparams.loc["CapitalCost"] =  hygparams.loc["CAPEX"]
    hygparams.loc["VariableCost"] =  hygparams.loc["Variable OPEX"]/3.6
    hygparams.loc["FixedCost"] =  hygparams.loc["Fixed OPEX"]                                         
    hygparams.loc["Efficiency"] =  1/hygparams.loc["Thermal conversion efficiency\n"
                                                 "(kWh input fuel or heat / kWh H2 HHV)*"].fillna(hygparams.loc["Electrical conversion efficiency \n(kWhe/ kWh H2 HHV)*"])                                         
    hygparams.loc["OperationalLife"] =  hygparams.loc["Operating lifetime"]
    hygparams.loc["AvailabilityFactor"] =  hygparams.loc["Availability/Maximum load factor"]
    
    # add CO2 storage cost
    hygparams.loc["VariableCost",
                  ("Steam Methane Reformation with CCUS",slice(None))] =     hygparams.loc["VariableCost",
                                                                                           ("Steam Methane Reformation with CCUS",slice(None))]+ccs

    # save hurdle rate for use in preprocessing to estimate price
    hurdle = pd.concat([hygparams.loc["Hurdle Rate",hytechs.keys()],
                          elgparams.loc["Hurdle Rate",eltechs.keys()]],axis=0)
    hurdle = hurdle.rename(index={**eltechs,**hytechs})
    
    hurdle.to_csv(snakemake.output.path_supply_hurdle)
    
    # rearrange params
    supparams = pd.concat([hygparams.loc[["CapitalCost","VariableCost",
                                         "FixedCost","Efficiency",
                                         "OperationalLife",
                                         "AvailabilityFactor"],hytechs.keys()],
                          elgparams.loc[["CapitalCost","VariableCost",
                                         "FixedCost","Efficiency",
                                         "OperationalLife",
                                         "AvailabilityFactor"],eltechs.keys()]],axis=1)
    
    # create and rearrange dataframe for capital cost
    capcost = supparams.loc["CapitalCost"]
    capcost = capcost.rename(index={**eltechs,**hytechs})
    capcost.name = "VALUE"
    capcost = capcost.reset_index()
    capcost.loc[capcost["TECHNOLOGY"].str.contains("HP"),
                "UNIT"] = "million 2020£/GW"
    capcost.loc[~capcost["TECHNOLOGY"].str.contains("HP"),
                "UNIT"] = "million 2018£/GW"
    capcost["VARIABLE"] = "CapitalCost"
    capcost["REGION"] = ":*0"

    
    # create and rearrange dataframe for fixcost
    fixcost = supparams.loc["FixedCost"]
    fixcost = fixcost.rename(index={**eltechs,**hytechs})
    fixcost.name = "VALUE"
    fixcost = fixcost.reset_index()
    fixcost.loc[fixcost["TECHNOLOGY"].str.contains("HP"),
                "UNIT"] = "million 2020£/GW"
    fixcost.loc[~fixcost["TECHNOLOGY"].str.contains("HP"),
                "UNIT"] = "million 2018£/GW"
    fixcost["VARIABLE"] = "FixedCost"
    fixcost["REGION"] = ":*0"

    # create and rearrange dataframe for variable cost
    varcost = supparams.loc["VariableCost"]
    varcost = varcost.rename(index={**eltechs,**hytechs})
    varcost.name = "VALUE"
    varcost = varcost.reset_index()
    varcost["MODE_OF_OPERATION"] = "1"
    varcost.loc[varcost["TECHNOLOGY"].str.contains("HP"),
                "UNIT"] = "million 2020£/TJ"
    varcost.loc[~varcost["TECHNOLOGY"].str.contains("HP"),
                "UNIT"] = "million 2018£/TJ"
    varcost["VARIABLE"] = "VariableCost"
    varcost["REGION"] = ":*0"
    
    
    # additional techs (nuclear, coal, etc.) and transmission cost
    sup_data = pd.read_csv(snakemake.input.path_set_supply)  
    
    costadd = sup_data[["TECHNOLOGY","UNIT", "MODE_OF_OPERATION",
                       "VARIABLE","REGION",
                       "YEAR",
                       "VALUE"]][(sup_data["VARIABLE"]=="VariableCost")|
                                 (sup_data["VARIABLE"]=="FixedCost")|
                                 (sup_data["VARIABLE"]=="CapitalCost")] 
    costadd["REGION"] = costadd["REGION"].fillna(":*0")
    
    # cost are saved to file, after adding transmission cost below (see rescap
    # section)

    # create and rearrange dataframe for efficiency
    eff = supparams.loc["Efficiency"]
    eff = eff.rename(index={**eltechs,**hytechs})
    eff.name = "VALUE"
    eff = eff.reset_index()
    eff["MODE_OF_OPERATION"] = "1"
    eff["UNIT"] = "-"
    eff["VARIABLE"] = "Efficiency"
    eff["REGION"] = ":*0"
    eff["FUEL_IN"] = eff["TECHNOLOGY"].replace({**eltechs_infuel,
                                                **hytechs_infuel})
    eff["FUEL_OUT"] = eff["TECHNOLOGY"].replace({**eltechs_outfuel,
                                                **hytechs_outfuel})
    
    # additional techs (nuclear, coal, etc.) and transmission efficiencies
    eff_data = pd.read_csv(snakemake.input.path_set_supply)
    
    effadd = eff_data[["TECHNOLOGY","UNIT", "MODE_OF_OPERATION",
                       "VARIABLE","REGION","FUEL_OUT","FUEL_IN",
                       "YEAR",
                       "VALUE"]][eff_data["VARIABLE"]=="Efficiency"] 
    effadd["REGION"] = effadd["REGION"].fillna(":*0")
    eff = pd.concat([eff,effadd,ffeff])
    
    eff = eff.set_index([col for col 
                         in eff.columns
                         if col!="VALUE"])  
    # save to file
    # unit: -
    eff.to_csv(snakemake.output.path_supply_eff)

    # create and rearrange dataframe for the availability factor
    lifetime = supparams.loc["OperationalLife"].dropna()
    lifetime = lifetime.rename(index={**eltechs,**hytechs})
    # lifetime in data does not change with time (and not implemented in
    # OSeMOSYS), dropping year
    lifetime = lifetime.xs(2025,level=1)
    lifetime.name = "VALUE"
    lifetime = lifetime.reset_index()
    lifetime["UNIT"] = "a"
    lifetime["VARIABLE"] = "OperationalLife"
    lifetime["REGION"] = ":*0"
    
    # additional techs (nuclear) and transmission lifetimes
    sup_data = pd.read_csv(snakemake.input.path_set_supply)
    
    lifeadd = sup_data[["TECHNOLOGY","UNIT",
                       "VARIABLE","REGION",
                       "VALUE"]][(sup_data["VARIABLE"]=="OperationalLife")] 
    lifeadd["REGION"] = lifeadd["REGION"].fillna(":*0")

    lifetime = pd.concat([lifetime,
                          lifeadd])

    lifetime = lifetime.set_index([col for col 
                                   in lifetime.columns
                                   if col!="VALUE"])     
    # save to file
    # unit: a
    lifetime.to_csv(snakemake.output.path_supply_lt)
    
    # create and rearrange dataframe for the availability factor
    availfac = supparams.loc["AvailabilityFactor"].dropna()
    availfac = availfac.rename(index={**eltechs,**hytechs})
    # discard solar and wind factors (added separately below)
    availfac = availfac[~availfac.index.get_level_values(level=0).str.startswith(
                        ("OWPP","FWPP","OSPP","RSPP"))]
    availfac.name = "VALUE"
    availfac = availfac.reset_index()
    availfac["UNIT"] = "-"
    availfac["VARIABLE"] = "AvailabilityFactor"
    availfac["REGION"] = ":*0"
    
    # additional techs (nuclear, etc.) and transmission availability factors
    sup_data = pd.read_csv(snakemake.input.path_set_supply)
    
    addavailfac = sup_data[["TECHNOLOGY","UNIT","YEAR",
                       "VARIABLE","REGION",
                       "VALUE"]][(sup_data["VARIABLE"]=="AvailabilityFactor")] 
    addavailfac["REGION"] = addavailfac["REGION"].fillna(":*0")

    availfac = pd.concat([availfac,
                          addavailfac,
                          ])
    availfac = availfac.set_index([col for col 
                                   in availfac.columns
                                   if col!="VALUE"])     
    # save to file
    # unit: -
    availfac.to_csv(snakemake.output.path_supply_af)

    # create and rearrange dataframe for emission factors
    
    emf_data = pd.read_csv(snakemake.input.path_set_supply)
    
    emfac = emf_data[["TECHNOLOGY","EMISSION","UNIT",
                      "MODE_OF_OPERATION",
                       "VARIABLE","REGION",
                       "YEAR",
                       "VALUE"]][emf_data["VARIABLE"]=="EmissionActivityRatio"] 


    emfac = emfac.set_index([col for col 
                                   in emfac.columns
                                   if col!="VALUE"])    
    # save to file
    # unit: kt/TJ 
    emfac.to_csv(snakemake.output.path_supply_emf)
    
    
    # load capacity factors and rearrange dataframe
    capfac = pd.read_csv(snakemake.input.path_tperiods_supply,
                         index_col=["TECHNOLOGY",
                                    "REGION"])
    capfac.columns.name="TIMESLICE"
    capfac = capfac.stack()
    capfac.name = "VALUE"
    capfac = capfac.reset_index()
    capfac["UNIT"] = "-"
    capfac["YEAR"] = ":*"
    capfac["VARIABLE"] = "CapacityFactor"

    capfac = capfac.set_index([col for col 
                               in capfac.columns
                               if col!="VALUE"])  
    # save to file
    # unit: -
    capfac.to_csv(snakemake.output.path_supply_capfac)
    
      
    
    # create residual capacity
    
    # use dukes power plant list for non-renewable power plants
    res_dukes = pd.read_excel(snakemake.input.path_supply_el_resc,
                              sheet_name="5.11 Full list",
                              skiprows=5)
    

    # load postcode lookup and allocate LAD code where possible
    pcd = utils.get_entity_lookup(levels=["PCD","LAD"]) 
    pcd = pcd[["PCD","LAD23CD"]].set_index("PCD")["LAD23CD"]
    pcd.index = pcd.index.str.replace(" ","")
    res_dukes["Postcode"] = res_dukes["Postcode"].map(pcd)
    res_dukes = res_dukes.rename(columns={"Postcode":"LAD23CD"})


    # allocate remaining to LADs by using coordinates
    res_dukes_xy = res_dukes[res_dukes["LAD23CD"].isna()&
                              res_dukes["X-Coordinate"].notna()]
    

 
    res_dukes_xy["geometry"] = gpd.points_from_xy(res_dukes_xy["X-Coordinate"],
                                                  res_dukes_xy["Y-Coordinate"])
    res_dukes_gdf = gpd.GeoDataFrame(res_dukes_xy,crs="EPSG:27700")
    lads = gpd.read_file(snakemake.input.path_lad_bounds) 
    res_dukes_xy["LAD23CD"] = gpd.sjoin(res_dukes_gdf,
                                        lads,how="left")["LAD23CD_right"].values
    res_dukes_xy = res_dukes_xy[res_dukes_xy["LAD23CD"].notna()]
    res_dukes_xy = res_dukes_xy.drop(columns="geometry")
    
    # insert LAD codes from spatial join using coordinates
    res_dukes = res_dukes.fillna(res_dukes_xy)
    
    # rows with neither coordinates nor postcode are discarded, these are
    # renewable capacities (or in NI - with one minor exception of 20 MW) that  
    # are accounted for through other file
    res_dukes = res_dukes[res_dukes["LAD23CD"].notna()]
    
    # allocate to technologies
    res_dukes["TECHNOLOGY"] = ["OSPPSNAT00" if ty=="PV" else
                                "OWPPSNAT00" if ty=="Onshore" else
                                "FWPPSNAT00" if ty=="Offshore" else
                                "NUPPSNAT00" if te=="Nuclear" else
                                "BMPPSNAT00" if te=="Bioenergy" else
                                "HYPPSNAT00" if te=="Hydro" else
                                "CLPPSNAT00" if (te=="Fossil Fuel")
                                                and (pf=="Coal") else
                                "NGPPSNAT00" if (te=="Fossil Fuel")
                                                and ((pf=="Natural Gas")
                                                     or(pf=="Diesel/Gas Oil")
                                                     or(pf=="Sour Gas"))
                                                else
                                np.nan
                                for te,ty,pf in zip(res_dukes["Technology"],
                                                    res_dukes["Type"],
                                                    res_dukes["Primary Fuel"])]
    # drop renewable capacities (processed above for potential future use)
    res_dukes = res_dukes[res_dukes["TECHNOLOGY"].str.startswith(("NU",
                                                                  "NG",
                                                                  "CL"))]
    # calculate in GW
    res_dukes.loc[:,"VALUE"] = res_dukes["InstalledCapacity (MW)"]/1000
    
    # calculate residual capacities for each year
    res_dukes = res_dukes.merge(right=lifetime,
                                            on = "TECHNOLOGY",
                                            how="left",
                                            suffixes=("","_lt"))
    
    res_dukes.loc[:,"DECOM_YEAR"] = res_dukes["Year Commissioned"]+res_dukes["VALUE_lt"]
    # if calculated decomissioning year is before 2022, assume quarter of
    # operational lifetime left (assuming retrofitting has extended lifetime)
    res_dukes.loc[res_dukes["DECOM_YEAR"]<2022,
                  "DECOM_YEAR"] = (0.25*res_dukes["VALUE_lt"])+2022
    
    for y in range(2015,2061):
        res_dukes[y] = [ v if ((y<dy) and (y>=cy)) else
                        0
                        for v,cy,dy in zip(res_dukes["VALUE"],
                                        res_dukes["Year Commissioned"],
                                        res_dukes["DECOM_YEAR"])]
    
    res_dukes = res_dukes.loc[:,["LAD23CD","TECHNOLOGY"]
                              +list(range(2015,2061))]
    
    res_dukes = res_dukes.groupby(["LAD23CD","TECHNOLOGY"]).sum()
    res_dukes.columns.names = ["YEAR"]
    res_dukes = res_dukes.stack()
    res_dukes.name = "VALUE"
    
   
    
    res_re = dict()
    for y in range(2015,2022):
        res_re[y] = pd.read_excel(snakemake.input.path_supply_el_rescr,
                                  sheet_name="LA - Capacity, "+str(y),
                                  skiprows=5)
        # remove NI and minor fraction of unallocated
        res_re[y] = res_re[y][(res_re[y]["Country"]=="England")|
                              (res_re[y]["Country"]=="Wales")|
                              (res_re[y]["Country"]=="Scotland")]
        res_re[y] = res_re[y][["Local Authority Code [note 1]",
                                "Photovoltaics","Onshore Wind",
                                "Hydro","Offshore Wind",
                                "Plant Biomass"]]
        res_re[y] = res_re[y].rename(
                    columns={"Local Authority Code [note 1]":"LAD21CD",
                                "Photovoltaics":"OSPPSNAT00",
                                "Onshore Wind":"OWPPSNAT00",
                                "Offshore Wind":"FWPPSNAT00",
                                "Hydro":"HYPPSNAT00",
                                "Plant Biomass":"BMPPSNAT00"})
        
        res_re[y] = res_re[y].set_index(["LAD21CD"])
        # convert to GW
        res_re[y] = res_re[y]/1000
        
        lads = utils.get_entity_lookup("LAD")["LAD23CD"]
        if y<2019 and y>2015:
            res_re[y] = utils.update_LADCD(res_re[y],
                                            from_CD="LAD20CD",
                                            how="sum")
        elif y==2015:
            res_re[y] = utils.update_LADCD(res_re[y],
                                            from_CD="LAD17CD",
                                            how="sum")
        else:
            res_re[y] = utils.update_LADCD(res_re[y],
                                            from_CD="LAD21CD",
                                            how="sum")
        res_re[y] = lads.to_frame().merge(right=res_re[y],on="LAD23CD",
                                          how="left")
        res_re[y]["YEAR"] = y
        res_re[y].columns.name = "TECHNOLOGY"
        res_re[y] = res_re[y].set_index(["LAD23CD","YEAR"]).stack()
        
    res_re = pd.concat([res_re[k] for k in res_re.keys()])
    res_re.name = "VALUE"
    
    # concat residual capacities
    
    rescaps = pd.concat([res_dukes.reset_index(),
                          res_re.reset_index()])
    

    rescaps = rescaps.rename(columns= {"LAD23CD":"REGION"})
    
    # currently aggregating to national level (next two lines can be removed for
    # local version)
    rescaps["REGION"] = ":*0"
    rescaps = rescaps.groupby(["REGION","TECHNOLOGY","YEAR"]).sum().reset_index()
    
    
    rescaps = rescaps.set_index(["REGION","TECHNOLOGY","YEAR"]).unstack("YEAR")
    rescaps = rescaps.droplevel(0,axis=1)

    # create and rearrange dataframe for the residual capacity - none
    # for generation but for transmission
    # FIXME: potentially put this into function given it is used here and
    # in other supply calculations
    
    # load residual distribution network capacity
    slagg = pd.read_csv(snakemake.input.path_sublocal_agg,
                            index_col=["LSOA11CD"]).squeeze()
    raw_data = pd.read_csv(snakemake.input.path_set_dist_net)
    EL_spare_cap = raw_data.loc[raw_data["VARIABLE"]=="ELNetworkSpareCapacity",
                                "VALUE"].squeeze()
    dist_rescap = pd.read_csv(snakemake.input.path_local_res_caps,
                              index_col=[slagg.name,"TECHNOLOGY"])
    
    # adjust spare capacity based on scenario
    if snakemake.params.dic["scen_networks"]=="-":
        pass
    else:
        EL_spare_cap = EL_spare_cap*float(snakemake.params.dic["scen_networks"])
    
    # set power transmission residual as current distribution need (not taking
    # into account spare capacity/headroom)
    el_tra_rescap = dist_rescap.xs("ELGRTDIS00",level=1) * (1-EL_spare_cap)
    el_tra_rescap["TECHNOLOGY"] = "ELGRTTRA00"
    
    # set gas transmission residual as current distribution need 
    # applying uplift to avoid any issue with slightly too small grid capacity
    ga_tra_rescap = dist_rescap.xs("GAGRTDIS00",level=1) * 1.05
    
    life_data = pd.read_csv(snakemake.input.path_set_lt)
    GA_rescap_red_start = life_data[["TECHNOLOGY","UNIT",
                         "VALUE"]][(life_data["VARIABLE"]=="ResCapRedStart")&
                                   (life_data["TECHNOLOGY"]=="NGGRTTRA00")]
    GA_rescap_red_end = life_data[["TECHNOLOGY","UNIT",
                         "VALUE"]][(life_data["VARIABLE"]=="ResCapRedEnd")&
                                   (life_data["TECHNOLOGY"]=="NGGRTTRA00")]                            
    py = [str(i) for i in range(2015,2021)]
    ga_tra_rescap = ga_tra_rescap[py]
    ly = py[-1]
    srescapyears = range(int(ly)+1,
                         int(GA_rescap_red_start["VALUE"].squeeze()+1))
    erescapyears = range(int(ly)+1,
                         int(GA_rescap_red_end["VALUE"].squeeze()+1))
    ga_tra_rescap[[str(y) for y in erescapyears]] = pd.NA
    ga_tra_rescap[[str(y) for y in srescapyears]]=(ga_tra_rescap[[ly]
                                                     *len(srescapyears)])
    ga_tra_rescap[str(int(GA_rescap_red_end["VALUE"].squeeze()))] = 0
    ga_tra_rescap = ga_tra_rescap.apply(pd.to_numeric)
    ga_tra_rescap = ga_tra_rescap.interpolate(axis=1)
    ga_tra_rescap["TECHNOLOGY"] = "NGGRTTRA00"
     
    # concat all residual capacities
    tra_rescaps = pd.concat([el_tra_rescap,ga_tra_rescap])
 
    # rearrange and groupby LAD
    tra_rescaps = tra_rescaps.reset_index()
    tra_rescaps["REGION"] = tra_rescaps[slagg.name].str.split("|",expand=True)[0]
    tra_rescaps = tra_rescaps.drop(slagg.name,axis=1)
    tra_rescaps = tra_rescaps.groupby(["REGION","TECHNOLOGY"]).sum()
    tra_rescaps.columns = tra_rescaps.columns.astype(int)

    rescaps = pd.concat([rescaps,tra_rescaps])
 
    # save to file
    # unit: GW
    rescaps.to_csv(snakemake.output.path_supply_rescap)
    
    # calc and add H2 and NG transmission cost
    # this simply calculates cost based on a cost/length value, the current
    # network length and current residual capacity in the model
    tra_data = pd.read_csv(snakemake.input.path_set_supply)
    
    ngcost = tra_data.loc[(tra_data["VARIABLE"]=="NGTraNetworkCost"),"VALUE"].squeeze()
    h2cost = tra_data.loc[(tra_data["VARIABLE"]=="H2TraNetworkCost"),"VALUE"].squeeze()
    nglen = tra_data.loc[(tra_data["VARIABLE"]=="NGTraNetworkLength"),"VALUE"].squeeze()
    omfrac = tra_data.loc[(tra_data["VARIABLE"]=="O&MCostFraction"),"VALUE"].squeeze()
    
    gatracapcost = pd.DataFrame([("NGGRTTRA00",(ngcost/1000
                                             *nglen
                                             /ga_tra_rescap["2015"].sum())),
                              ("H2GRTTRA00",(h2cost/1000
                                             *nglen
                                             /ga_tra_rescap["2015"].sum()))],
                             columns=["TECHNOLOGY","VALUE"])
    gatracapcost["UNIT"] = "million 2016£/GW"
    gatracapcost["VARIABLE"] = "CapitalCost"
    gatracapcost["REGION"] = ":*2"
    gatracapcost["YEAR"] = ":*"
    
    gatrafixcost = gatracapcost.copy()
    gatrafixcost["VARIABLE"] = "FixedCost"
    gatrafixcost["VALUE"] = gatrafixcost["VALUE"]*omfrac
    
    cost = pd.concat([capcost,
                      gatracapcost,
                      fixcost,
                      gatrafixcost,
                      varcost,
                      costadd,
                      ffvarcost])
    cost = cost.set_index([col for col 
                           in cost.columns
                           if col!="VALUE"])  
    # save to file
    # unit (fixed and capital): million [Baseyear]£/GW (= [Baseyear]£/kW)
    # unit (variable cost): million [Baseyear]£/TJ (= [Baseyear]£/TJ) 
    cost.to_csv(snakemake.output.path_supply_cost)

    # create capacity constraints based on potentials
    
    luew = pd.read_csv(snakemake.input.path_ward_lookup_11,usecols=["WD11CD","LAD22CD"])
    luew = luew[~luew["WD11CD"].duplicated(keep="last")]
    luew = luew.set_index("LAD22CD")
    luew["WD11CD"] = luew["WD11CD"] +","
    luew = utils.update_LADCD(luew,from_CD="LAD21CD",how="sum")
    luew["WD11CD"] = luew["WD11CD"].str[:-1]
    luew = luew.assign(WD11CD=luew["WD11CD"].str.split(',')).explode('WD11CD')
    luew = luew.reset_index()
    
    lu = pd.read_csv(snakemake.input.path_ward_lookup_16,index_col="LAD16CD",usecols=["WD16CD","LAD16CD"])
    lu["WD16CD"] = lu["WD16CD"] +","
    lu = utils.update_LADCD(lu,from_CD="LAD17CD",how="sum")
    lu["WD16CD"] = lu["WD16CD"].str[:-1]
    lu = lu.assign(WD16CD=lu["WD16CD"].str.split(',')).explode('WD16CD')
    lu = lu.reset_index()
    #loup.columns = ["WD11CD","WD11NM","LAD22CD","LAD22NM","FID"]
    
    pot = pd.read_csv(snakemake.input.path_supply_pot)
    pot = pot[pot["id"].str.startswith("UK")]
    pot["id"] = pot["id"].str[2:]
    pot = pot.fillna(0)
    
    pot1 = luew.merge(right=pot,left_on="WD11CD",right_on="id",how="right")
    pot2 = lu.merge(right=pot,left_on="WD16CD",right_on="id",how="inner")
    
    pot1 = pot1[["LAD23CD"]+list(pot.columns)].set_index("id")
    pot2 = pot2[["LAD23CD"]+list(pot.columns)].set_index("id")
    
    # manuall allocation for a few wards not captured by the two lookups (no
    # one complete lookup available for the ward code version of the data)
    mall = {"E05008953":"E06000030",
            "E05008959":"E06000030",
            "E05008964":"E06000030",
            "E05008973":"E07000220",
            "E05008978":"E07000220",
            "E05008987":"E07000220",
            "E05009029":"E07000241"}
    
    for k,v in mall.items():
        pot1.loc[k,"LAD23CD"]=v
        
    ladpot = pot1.fillna(pot2).dropna()
    ladpot = ladpot.groupby("LAD23CD").sum()

    
    ladpot = ladpot.rename(columns={"rooftop_pv_mw":"RSPPSNAT00",
                                    "open_field_pv_mw":"OSPPSNAT00",
                                    "onshore_wind_mw":"OWPPSNAT00",
                                    "offshore_wind_mw":"FWPPSNAT00"})
    # convert to GW
    ladpot = ladpot/1000
    
    ladpot.index.name = "REGION"
    
    # currently aggregating to national level (next three lines can be removed 
    # for local version)
    ladpot.index = [":*0"]*len(ladpot)
    ladpot.index.name = "REGION"
    ladpot = ladpot.groupby("REGION").sum()
    
    # add max of residual hydro capacity, i.e., do not allow further expansion
    ladpot["HYPPSNAT00"] = rescaps.xs("HYPPSNAT00",level=1).max(axis=1)
    ladpot.columns.name = "TECHNOLOGY"
    ladpot = ladpot.stack()
    ladpot.name = "VALUE"
    ladpot = ladpot.reset_index()
    
    ladpot["YEAR"] = ":*"
    ladpot["UNIT"] = "GW"
    ladpot["VARIABLE"] = "TotalAnnualMaxCapacity"
    
    # add capacity investment constraint

    sup_data = pd.read_csv(snakemake.input.path_set_supply)
    
    capicon = sup_data[["TECHNOLOGY","UNIT",
                       "VARIABLE","REGION",
                       "YEAR",
                       "VALUE"]][(sup_data["VARIABLE"]=="TotalAnnualMaxCapacityInvestment")] 
    capicon["REGION"] = capicon["REGION"].fillna(":*0")
    
    capcon = pd.concat([ladpot,capicon])
    
    capcon = capcon.set_index([col for col 
                                in capcon.columns
                                if col!="VALUE"])  
    # save to file
    # unit: GW
    capcon.to_csv(snakemake.output.path_supply_capcon)
    
    # add activity constraint - currently empty/not used

    sup_data = pd.read_csv(snakemake.input.path_set_supply)
    
    actcon = sup_data[["TECHNOLOGY","UNIT",
                       "VARIABLE","REGION",
                       "YEAR",
                       "VALUE"]][(sup_data["VARIABLE"]=="TotalTechnologyAnnualActivityLowerLimit")] 
    actcon["REGION"] = actcon["REGION"].fillna(":*0")
    
    
    actcon = actcon.set_index([col for col 
                                in actcon.columns
                                if col!="VALUE"])  
    # save to file
    # unit: GW
    actcon.to_csv(snakemake.output.path_supply_actcon)
    
    
    # set param for simple storage representation
    
    tagsto = pd.DataFrame([[":*0","E0SSTO","StorageTagFuel",1],
                           [":*0","H0SSTO","StorageTagFuel",1],
                           [":*0","H0STRA","StorageTagFuel",1]],
                          columns=["REGION","FUEL","VARIABLE","VALUE"])
    tagsto = tagsto.set_index(["REGION","FUEL","VARIABLE"])
    
    # save to file
    # unit: -
    tagsto.to_csv(snakemake.output.path_supply_storage)
    
    # set capacity constraints based on Future Energy Scenarios (FES) by 
    # National Grid
    
    eltechs = {"ES.E.13":"FWPPSNAT00",
               "ES.E.14":"OWPPSNAT00",
               "ES.E.16":"Solar",
               "ES.E.17":"NGPPSNAT00",
               "ES.E.18":"NGCCS",
               "ES.E.19":"BECCS",
               "ES.E.20":"BMPPSNAT00",
               "ES.E.21":"NUPPSNAT00",
               "ES.E.23":"H2PPSNAT00",
               "ES.E.26":"BSPPSNAT00",
               "FL.13":"AESPSNAT00" 
              }
    
    fes_data = pd.read_excel(snakemake.input.path_supply_prices,
                          sheet_name=list(eltechs.keys()))# nrows=42,header=0,usecols="C:AU"
    
    for k in eltechs.keys():
        if k=="ES.E.16":
            r = 4
            c = 12
        elif (k=="ES.E.17" or k=="ES.E.18" or k=="ES.E.19" or k=="ES.E.20" 
            or k=="ES.E.21" or k=="ES.E.23" or k=="ES.E.26"):
            r = 4
            c = 8
        elif k=="FL.13":
            r = 7
            c = 14
        else:
            r = 5
            c = 13
        
        fes_data[k] = fes_data[k].iloc[r:r+5,c:]
        fes_data[k] = fes_data[k].dropna(how="all",axis=1)
        fes_data[k].columns = (fes_data[k].iloc[0,:].fillna("SCENARIO")
                                                    .replace("\xa0","SCENARIO"))
        fes_data[k] = fes_data[k].iloc[1:,:].dropna(how="all",axis=1)
        fes_data[k] = fes_data[k].set_index("SCENARIO")
        fes_data[k].columns = pd.DatetimeIndex(fes_data[k].columns).year
    
    
    
    # split solar capacity into rooftop and open field based on potential
    frt = (ladpot.loc[ladpot["TECHNOLOGY"]=="RSPPSNAT00","VALUE"].squeeze()
           /ladpot.loc[ladpot["TECHNOLOGY"].str.startswith(("RS","OS")),"VALUE"].sum())
    
    fes_data["ES.E.16.RT"] = fes_data["ES.E.16"] * frt
    fes_data["ES.E.16.OF"] = fes_data["ES.E.16"] * (1-frt)
    del eltechs["ES.E.16"]
    eltechs["ES.E.16.RT"] = "RSPPSNAT00"
    eltechs["ES.E.16.OF"] = "OSPPSNAT00"
    
    # add BECCS to normal biomass
    fes_data["ES.E.20"] = fes_data["ES.E.20"]+fes_data["ES.E.19"]
    del eltechs["ES.E.19"]  
    # add Gas CCS to hydrogen power plant as CCS is not represented in model
    # the chosen scenario (Leading the Way) does not include a lot of Gas CCS
    # in any case
    fes_data["ES.E.23"] = fes_data["ES.E.23"]+fes_data["ES.E.18"]
    del eltechs["ES.E.18"]  
    
    # concat data
    caps = pd.concat([fes_data[k] for k in eltechs.keys()],
                      keys=eltechs.values(),
                      names=["TECHNOLOGY"],
                      axis=0)
    
    # choose scenario
    if snakemake.params.dic["scen_supply_imp"].endswith("ST"):
        caps = caps.xs("System Transformation",level=1)
    else:
        caps = caps.xs("Leading the Way",level=1)
    
    # calculate fractions
    caps = caps/caps.sum()
    
    # rearrange dataframe and save
    caps.columns.name = "YEAR"
    caps = caps.stack()
    caps.name = "VALUE"
    caps = caps.reset_index()
    
    caps["UNIT"] = "-"
    caps["REGION"] = ":*0"
    caps["VARIABLE"] = "CapacityFraction"
    
    # add technology tag param
    techtag = pd.DataFrame(caps["TECHNOLOGY"].unique(),columns=["TECHNOLOGY"])
    techtag["VALUE"] = 1
    techtag["UNIT"] = "-"
    techtag["REGION"] = ":*0"
    techtag["VARIABLE"] = "CapacityFractionTagTechnology"
    
    # do not use natural gas power plants (causes issues given it needs
    # phasing out)
    techtag.loc[techtag["TECHNOLOGY"]=="NGPPSNAT00","VALUE"] = 0
    
    caps = pd.concat([caps,techtag])
    
    
    caps = caps.set_index([col for col 
                                in caps.columns
                                if col!="VALUE"])  
    # save to file
    # unit: -
    caps.to_csv(snakemake.output.path_supply_capfrac)
   

if __name__ == "__main__":

    if "snakemake" not in globals():
        df_params = pd.read_csv("./runs/default_runs.csv")
        snakemake = utils.mock_snakemake("calc_supply_techs",
                                         **df_params.iloc[0].to_dict())
        
    # create aggregate or detailed supply sector
    if snakemake.params.dic["scen_supply_imp"]=="agg":
        
        # currently not implemented
        # create_aggregated_supply()
        pass
    
    elif snakemake.params.dic["scen_supply_imp"].startswith("det"):

        create_detailed_supply()                             

    else:
        raise NotImplementedError(
                f"'{snakemake.params.dic['scen_supply_imp']}' is currently"
                " not a valid option of the 'scen_supply_imp' parameter."
                " Valid options are 'agg' or 'det'."
            )      
        
