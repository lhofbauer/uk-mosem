"""
Module with utility functions


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
includes code with Copyright 2017-2023 The PyPSA-Eur Authors
"""

import sys
import logging

import yaml
import pandas as pd
import numpy as np


# load config, in particular file paths
with open("config.yaml", 'r') as file:
    config = yaml.safe_load(file)


def adjust_logger(level="INFO"):
    root_logger = logging.getLogger()
    
    # remove existing handlers
    if root_logger.hasHandlers():
        for h in root_logger.handlers:
            root_logger.removeHandler(h)
    
    
    console = logging.StreamHandler(stream=sys.stdout)
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
    )
    console.setFormatter(formatter)
    root_logger.addHandler(console)
    root_logger.setLevel(level) 
    

def update_LADCD(df, from_CD="LAD17CD", how="sum"):
    """ Convert data index to LAD23CD
    

    Parameters
    ----------
    df : DataFrame
        Dataframe to be converted. It needs to be indexed over the LAD code.
    from_CD : str
        LAD code version of the input DataFrame. Options are 'LAD17CD',
        'LAD19CD', or 'LAD20CD. The default is 'LAD17CD'.
    how : str, optional
        Function to be applied to aggegrate data. Options are 'sum' or 'mean'.
        The default is 'sum'.

    Returns
    -------
    df : DataFrame
        Dataframe indexed over LAD23CD

    """

    df = df.rename(index={"S12000024":"S12000048",
                          "S12000015":"S12000047",
                          "S12000044":"S12000050",
                          "S12000046":"S12000049",
                          "E06000048":"E06000057",
                          "E08000020":"E08000037"})
    
    df.index.name = "LAD23CD"
    
    if how == "sum" or how == "mean":
        
        if from_CD=="LAD17CD":
            
            # Dorset Council
            df.loc["E06000059"] = getattr(df.loc[["E07000049",
                                                  "E07000050",
                                                  "E07000051",
                                                  "E07000052",
                                                  "E07000053"]],how)(axis=0)
            df = df.drop(["E07000049",
                          "E07000050",
                          "E07000051",
                          "E07000052",
                          "E07000053"])
            # BCP council
            df.loc["E06000058"] = getattr(df.loc[["E06000028",
                                                  "E06000029",
                                                  "E07000048"]],how)(axis=0)
            df = df.drop(["E06000028",
                          "E06000029",
                          "E07000048"])
            
            # West Suffolk council
            df.loc["E07000245"] = getattr(df.loc[["E07000201",
                                                  "E07000204"]],how)(axis=0)
            df = df.drop(["E07000201",
                          "E07000204"])
            
            # East Suffolk council
            df.loc["E07000244"] = getattr(df.loc[["E07000205",
                                                  "E07000206"]],how)(axis=0)
            df = df.drop(["E07000205",
                          "E07000206"])
            
            # Somerset and West Taunton council
            df.loc["E07000246"] = getattr(df.loc[["E07000190",
                                                  "E07000191"]],how)(axis=0)        
            df = df.drop(["E07000190",
                          "E07000191"])
            
        if (from_CD=="LAD17CD") or (from_CD == "LAD19CD"):
            # Buckinghamshire
            df.loc["E06000060"] = getattr(df.loc[["E07000004",
                                                  "E07000005",
                                                  "E07000006",
                                                  "E07000007"]],how)(axis=0)        
            df = df.drop(["E07000004",
                          "E07000005",
                          "E07000006",
                          "E07000007"])
        
        if ((from_CD=="LAD17CD") or (from_CD == "LAD19CD")
            or (from_CD == "LAD20CD")):
            
            # North Northamptshire council
            df.loc["E06000061"] = getattr(df.loc[["E07000150",
                                                  "E07000152",
                                                  "E07000153",
                                                  "E07000156"]],how)(axis=0)
            df = df.drop(["E07000150",
                          "E07000152",
                          "E07000153",
                          "E07000156"])

            # West Northamptshire council
            df.loc["E06000062"] = getattr(df.loc[["E07000151",
                                                  "E07000154",
                                                  "E07000155"]],how)(axis=0)
            df = df.drop(["E07000151",
                          "E07000154",
                          "E07000155"])
            
        if ((from_CD=="LAD17CD") or (from_CD == "LAD19CD")
            or (from_CD == "LAD20CD") or (from_CD == "LAD21CD")):
            
            # Cumberland council
            df.loc["E06000063"] = getattr(df.loc[["E07000026",
                                                  "E07000028",
                                                  "E07000029"]],how)(axis=0)
            df = df.drop(["E07000026",
                          "E07000028",
                          "E07000029"])

            # Westmorland and Furness council
            df.loc["E06000064"] = getattr(df.loc[["E07000027",
                                                  "E07000031",
                                                  "E07000030"]],how)(axis=0)
            df = df.drop(["E07000027",
                          "E07000031",
                          "E07000030"])
            
            
            # North Yorkshire council
            df.loc["E06000065"] = getattr(df.loc[["E07000163",
                                                  "E07000165",
                                                  "E07000164",
                                                  "E07000166",
                                                  "E07000167",
                                                  "E07000168",
                                                  "E07000169"]],how)(axis=0)
            df = df.drop(["E07000163",
                          "E07000165",
                          "E07000164",
                          "E07000166",
                          "E07000167",
                          "E07000168",
                          "E07000169"])
            
            # Somerset council
            df.loc["E06000066"] = getattr(df.loc[["E07000187",
                                                  "E07000188",
                                                  "E07000189",
                                                  "E07000246"]],how)(axis=0)
            df = df.drop(["E07000187",
                          "E07000188",
                          "E07000189",
                          "E07000246"])
            
        
    else:
        raise  ValueError(f"Function '{how}' is not implemented")
        
    return df


def groupby_LAD(df, disagg=None):
    """ Groupby data by local authority district with optional disaggregation
    

    Parameters
    ----------
    df : DataFrame
        Dataframe with data to be grouped by LAD. Needs to be indexed by
        LSOA11CD, MSOA11CD, or LAD23CD.
    disagg : Series
        Series that maps the values of the index of df, e.g., LSOAs, to a
        sublocal class that it will be aggregated to.

    Returns
    -------
    df_gb : DataFrameGroupBy
            Pandas DataFrameGroupBy object to perform aggregating function on

    """
        
    # load list of areas/area lookup (LAD23CD)
    lads = get_entity_lookup()
    
    # groupby and return
    index = df.index.names
    df = df.reset_index()
    df = lads.merge(right=df, on=index[0], how="left")
    
    if disagg is not None:
        df = df.merge(right=disagg, on=disagg.index.name, how="left")
        df.loc[:,disagg.name] = df["LAD23CD"] +"|"+ df[disagg.name]
        df_gb = df.groupby([disagg.name]
                           +[ind for ind in index[1:]])
    else:
        df_gb = df.groupby(["LAD23CD"]
                           +[ind for ind in index[1:]])    
    return df_gb


def get_entity_lookup(levels=["LSOA","MSOA","LAD"]):
    """ Get list/lookup of local areas
    

    Parameters
    ----------
    levels : list of str, optional
        The geographical levels to be included. Options are 'PCD', 'LSOA',
        'MSOA', 'LAD', 'ITL2','ITL3','RGN', and 'CTRY'. The default is 
        ["LSOA","MSOA","LAD"].

    Returns
    -------
    lookup: DataFrame
        Dataframe containing all the areas.

    """
    
    # create dict of columns to be loaded from lookup file along with how they
    # will be renamed
    cols = dict()
    if "PCD" in levels:
        cols["pcds"]="PCD"
    if "LSOA" in levels:
        cols["lsoa11cd"]="LSOA11CD"
    if "MSOA" in levels:
        cols["msoa11cd"]="MSOA11CD"            

    cols["ladcd"]="LAD21CD"
    cols["ladnm"]="LAD21NM"

   
    # load file of area lookup (LAD21CD), drop duplicates, and rename
    lookup = pd.read_csv(config["fp"]["PATHR"]+config["fp"]["PATH_GEO_LOOKUPS"],
                        usecols=cols.keys(),
                        encoding = "ISO-8859-1")
    lookup = lookup.drop_duplicates()
    lookup = lookup.dropna(axis=0, how="all")
    lookup = lookup.rename(columns=cols)
    
    if ("ITL2" in levels) or ("ITL3" in levels):
        lookupITL = pd.read_excel(config["fp"]["PATHR"]+config["fp"]["PATH_GEO_LOOKUPS_ITL"],
                            sheet_name="LAD21_LAU121_ITL21_UK_LU",
                            usecols=["LAD21CD"]+[i+"21CD" for i in levels if
                                                 "ITL" in i])
        lookup = lookup.merge(right=lookupITL,how="left",on="LAD21CD")
        lookup = lookup.drop_duplicates()
    if "RGN" in levels:
        lookupRGN = pd.read_csv(config["fp"]["PATHR"]+config["fp"]["PATH_GEO_LOOKUPS_RGN"],
                            usecols=["LAD21CD","RGN21CD","RGN21NM"])
        lookup = lookup.merge(right=lookupRGN,how="left",on="LAD21CD")
        lookup = lookup.drop_duplicates()
        # add country codes for Wales and Scotland
        lookup.loc[lookup["LAD21CD"].str.startswith("S"),"RGN21CD"] = "S92000003"
        lookup.loc[lookup["LAD21CD"].str.startswith("W"),"RGN21CD"] = "W92000004"
        lookup.loc[lookup["LAD21CD"].str.startswith("S"),"RGN21NM"] = "Scotland"
        lookup.loc[lookup["LAD21CD"].str.startswith("W"),"RGN21NM"] = "Wales"
    
    # update to LAD23CD (from LAD21CD)
    
    lookup = lookup.rename(columns={"LAD21CD":"LAD23CD",
                                    "LAD21NM":"LAD23NM"})  
    dic = {"E07000026":"E06000063",
           "E07000028":"E06000063",
           "E07000029":"E06000063",
           "E07000027":"E06000064",
           "E07000031":"E06000064",
           "E07000030":"E06000064",
           "E07000163":"E06000065",
           "E07000164":"E06000065",
           "E07000165":"E06000065",
           "E07000166":"E06000065",
           "E07000167":"E06000065",
           "E07000168":"E06000065",
           "E07000169":"E06000065",
           "E07000187":"E06000066",
           "E07000188":"E06000066",
           "E07000189":"E06000066",
           "E07000246":"E06000066"}
    lookup = lookup.replace(dic)
    
    lookup.loc[lookup["LAD23CD"]=="E06000063","LAD23NM"]="Cumberland"
    lookup.loc[lookup["LAD23CD"]=="E06000064","LAD23NM"]="Westmorland and Furness"
    lookup.loc[lookup["LAD23CD"]=="E06000065","LAD23NM"]="North Yorkshire"
    lookup.loc[lookup["LAD23CD"]=="E06000066","LAD23NM"]="Somerset"

    lookup = lookup.drop_duplicates()
    
    # drop large consumer and other postcodes (if used), Northern Ireland and
    # pseudo entries
    lookup = lookup.dropna(subset=["LAD23CD"],axis=0)
    lookup = lookup[~lookup["LAD23CD"].str.startswith("N")]
    
    lookup = lookup[lookup["LAD23CD"]!="M99999999"]
    lookup = lookup[lookup["LAD23CD"]!="L99999999"]
    
    # add country codes if necessary
    def CTRY (row):
        if row["LAD23CD"].startswith("E"):
            return "E92000001"
        elif row["LAD23CD"].startswith("S"):
            return "S92000003"
        if row["LAD23CD"].startswith("W"):
            return "W92000004"
        
    if "CTRY" in levels:
       lookup["CTRY11CD"] = lookup.apply(lambda row: CTRY(row), axis=1)
    
    if "LAD" not in levels:
        lookup = lookup.drop(["LAD23CD","LAD23NM"],axis=1).drop_duplicates()
    
    # return dataframe
    return lookup

    
def interpolate_timeseries(df, model_periods, model_periods_column="YEAR",
                           col_to_ipol="VALUE",
                           direction="both",
                           force_exp=False,
                           dropna=False):
    """ Interpolate between model time periods
    

    Parameters
    ----------
    df : DataFrame
        Pandas DataFrame for which values are to be interpolated. Any index is
        disregarded and sequential index is set.
    model_periods : list of int
        List of all model periods.
    model_periods_column : str
        Name of column describing model periods. The default is 'YEAR'.
    direction : str
        In which directions to allow interpolation. Either 'forward',
        'backward', or 'both'. The default is 'both'.
    force_exp : boolean
        If to set all years explicitely if ':*' is given. The default is False.        
    dropna : boolean
        If to drop rows with missing value when re-stacking the index.
        The default is False.

    Returns
    -------
    df_ip : DataFrame
        Dataframe where values have been interpolated.

    """
    
    col = df.columns
    ocol = [col for col in df.columns
            if not (col == model_periods_column or
                    col == col_to_ipol)]

    # only process where not all years are set with abbreviation
    df_tip = df[df[model_periods_column]!=":*"]
    df_tip.loc[:,model_periods_column] = df_tip[model_periods_column].astype(int)
    
    df_ntip = df.loc[df[model_periods_column]==":*",:].copy()
    
    if force_exp:
        df_ntip.loc[:,model_periods_column] = ",".join([str(a) for a in model_periods])
        df_ntip.loc[:,model_periods_column] = df_ntip.loc[:,model_periods_column].str.split(',')
        df_ntip = df_ntip.explode(model_periods_column)
        df_ntip.loc[:,model_periods_column] = df_ntip.loc[:,model_periods_column].astype(int)
        # df_ntip.assign(WD11CD=luew["WD11CD"].str.split(',')).explode('WD11CD')
        
    # return original dataframe if nothing to process
    if df_tip.empty:
        return df
    
    md_years = sorted([y for y in set(model_periods.tolist()
                           +df_tip[model_periods_column].unique().tolist())])

    # create new index with each model year present
    df_ind = df_tip[ocol]
    df_ind = df_ind.set_index(list(df_ind.columns))
    df_ind = df_ind[~df_ind.index.duplicated()]
    df_ind[md_years]=1
    df_ind.columns.name=model_periods_column
    df_ind = df_ind.stack().to_frame()

    # merge in values, and unstack index levels apart from model periods
    # for interpolation
    df_ip = df_ind.merge(right=df_tip,how="left",
                          left_index=True,right_on=ocol+[model_periods_column])[col]
    df_ip = df_ip.set_index(model_periods_column)   
    df_ip = df_ip.set_index(ocol, append=True)
    df_ip = df_ip.unstack(ocol)
    

    # interpolate and rearrange
    df_ip = df_ip.interpolate(method="index",axis="index",
                              limit_direction=direction)
    df_ip = df_ip.stack(ocol,dropna=dropna)
    df_ip = df_ip.reorder_levels([c for c in col
                                  if not c == col_to_ipol]).sort_index()
    df_ip = df_ip.drop([y for y in md_years if (y not in model_periods.tolist())],level=2)
    df_ip = df_ip.reset_index()  
    df = pd.concat([df_ntip,
                    df_ip],axis=0).reset_index(drop=True)

    return df

def aggregate_years(df, years_map, model_periods_column="YEAR",
                    col_to_agg="VALUE",
                    max_constraint_default=None,
                    agg_method="mean"):
    """ Aggregate years
    

    Parameters
    ----------
    df : DataFrame
        Pandas DataFrame for which years are to be aggregated. Any index is
        disregarded and sequential index is set.
    years_map : Series
        Series mapping each year to an (aggregate, multi-year period) year.
    model_periods_column : str
        Name of column describing model periods. The default is 'YEAR'.
    max_constraint_default: float
        Default for a maximum constraint. The default is None.
    agg_method : str
        String naming the aggregation method. The default is 'mean'.
        

    Returns
    -------
    df_agg : DataFrame
        Dataframe where years have been aggregated.

    """
    
    if df.empty or (not isinstance(years_map,pd.Series) and years_map == None):
        return df

    df = df.reset_index(drop=True)
    
    # only process where not all years are set with abbreviation
    if agg_method=="sum" and not df[df[model_periods_column]==":*"].empty:
        df_tagg = df
        df_tagg.loc[:,model_periods_column] = df_tagg[model_periods_column].astype(str)
        df_tagg.loc[:,model_periods_column] = df_tagg[model_periods_column].replace({":*":
                                                                         ",".join(list(years_map.index.astype(str)))})
        df_tagg.loc[:,model_periods_column] = df_tagg[model_periods_column].str.split(",") 
        df_tagg = df_tagg.explode(model_periods_column)                           

    else:   
        df_tagg =  df[df[model_periods_column]!=":*"]

    df_tagg.loc[:,model_periods_column] = df_tagg[model_periods_column].astype(int)
    
    # return original dataframe if nothing to process
    if df_tagg.empty:
        return df
    
    df_tagg.loc[:,model_periods_column] = df_tagg[model_periods_column].map(years_map)
    # if max constraint given, replace values with nan so those are not applied
    # across the aggregated periods (as it results in nan values that are
    # dropped)
    if max_constraint_default is not None:
        td = df_tagg.loc[df_tagg[col_to_agg]==max_constraint_default].drop_duplicates()
        td = td.set_index([c for c in td.columns if c!=col_to_agg])
        # df_tagg.loc[:,col_to_agg] = df_tagg.loc[:,col_to_agg].replace({max_constraint_default:np.nan})
    
    if agg_method=="mean":
        df_agg = df_tagg.groupby([c for c in df_tagg.columns if c!=col_to_agg]).mean()
    if agg_method=="sum":
        df_agg = df_tagg.groupby([c for c in df_tagg.columns if c!=col_to_agg]).sum()    
    
    if max_constraint_default is not None:
        df_agg.loc[df_agg.index.isin(td.index),:] = max_constraint_default
    
    df_agg = df_agg.reset_index()
    
    if agg_method=="sum" and not df[df[model_periods_column]==":*"].empty:
        df = df_agg
    else:
        df = pd.concat([df[df[model_periods_column]==":*"],
                        df_agg],axis=0).reset_index(drop=True)
    
    return   df

def adjust_monetary_values(df, base_year, unit_column="UNIT", deflator="GDP"):
    """ Adjust monetary value to base year £ using a price index, exchange rate
    

    Parameters
    ----------
    df : DataFrame
        Pandas DataFrame for which values are to be processed.
    base_year : int
        Base year to be used.
    unit_column : str
        Name of column describing monetary unit. Values must be in the format
        "[...] YYYY[CurrencySymbol] [...]. Allowed currencies are € and £.
        The default is 'UNIT'.
    deflator : str
        Name of deflator to be used. Allowed are 'GDP' or 'PPI'.
        The default is 'GDP.

    Returns
    -------
    df : DataFrame
        Dataframe where values have been adjusted.

    """
    
    # load and process producer price index
    if deflator=="PPI":
        pi= ("Inputs into Production of Electricity,"
               " Transmission and Distribution Services")    
        ppi = pd.read_csv(config["fp"]["PATHR"]+config["fp"]["PATH_PPI"],
                          index_col=["Title"],
                          usecols=[pi,"Title"])
        ppi = ppi.loc["2009":"2022",:].astype("float")
        ppi.index = ppi.index.astype("int32")
        ppi.index.name='Year'
    elif deflator=="GDP":
        pi= "Gross domestic product at market prices:Implied deflator:SA"    
        ppi = pd.read_csv(config["fp"]["PATHR"]+config["fp"]["PATH_GDP_DEF"],
                          index_col=["Title"],
                          usecols=[pi,"Title"])
        ppi = ppi.loc["2009":"2023",:].astype("float")
        ppi.index = ppi.index.astype("int32")
        ppi.index.name='Year'
    ppi.columns=["EPPI"]
    ppi["multiplier"] = ppi["EPPI"].loc[base_year] / ppi["EPPI"]
    
    # load and process exchange rates
    exr_eur = pd.read_csv(config["fp"]["PATHR"]+config["fp"]["PATH_EX_RATE_EUR"],
                      index_col=["Title"])
    exr_usd = pd.read_csv(config["fp"]["PATHR"]+config["fp"]["PATH_EX_RATE_USD"],
                      index_col=["Title"])
    exr_eur = exr_eur.loc["2009":"2020",:].astype("float")
    exr_usd = exr_usd.loc["2009":"2022",:].astype("float")
    exr_eur.index = exr_eur.index.astype("int32")
    exr_usd.index = exr_usd.index.astype("int32")
    exr_usd.index.name='Year'
    exr_eur.index.name='Year'
    exr_usd.columns=["USD_RATE"]
    exr_eur.columns=["EUR_RATE"]
    exr_usd["USD_RATE"] = 1/exr_usd["USD_RATE"]
    exr_eur["EUR_RATE"] = 1/exr_eur["EUR_RATE"]
    exr = pd.concat([exr_usd,exr_eur],axis=1)
    
    ppiexr = pd.concat([ppi["multiplier"],exr],axis=1)
    
    # merge multiplier on dataframe and process
    df["uyear"] = df[unit_column].str.split(r"£|€|\$").str[0].str[-4:].astype("int32")
    df = df.merge(right=ppiexr,how="left",left_on="uyear",
                   right_index=True)
    # apply exchange rate where necessary
    df.loc[df[unit_column].str.contains("€"),"multiplier"] = (
        df.loc[df[unit_column].str.contains("€"),"multiplier"]
        *df.loc[df[unit_column].str.contains("€"),"EUR_RATE"])
    df.loc[df[unit_column].str.contains("\$"),"multiplier"] = (
        df.loc[df[unit_column].str.contains("\$"),"multiplier"]
        *df.loc[df[unit_column].str.contains("\$"),"USD_RATE"])
    
    df["VALUE"] = df["VALUE"] * df["multiplier"]
    df[unit_column] = (df[unit_column].str.split(r"£|€|\$").str[0].str[:-4]
                       + str(base_year)+"£"
                       + df[unit_column].str.split(r"£|€|\$").str[1])
    df = df.drop(["uyear","multiplier","USD_RATE","EUR_RATE"],axis=1)
    
    return df


def round_col(s):
    """ Round values of a series with (almost) consistent sum
    

    Parameters
    ----------
    s : Series
        Pandas Series with values to be rounded.

    Returns
    -------
    src : Series
        Pandas Series with rounded values.

    """
    sr = s.round()
    nr = round(s.sum()) - sr.sum()
    nrs = np.sign(nr)
    srd = s-sr
    
    # old, much less efficient versions, left for reference only
    # ind = [i for o, (e, i) in enumerate(reversed(sorted((nrs*e,i)
    #                             for i,e in enumerate(srd)))) if o < nrs*nr]
    # ind = [i for (e, i) in sorted(((nrs*e,i)
    #                             for i,e in enumerate(srd)), key=lambda t: t[0], reverse=True)[0:int(nrs*nr)]]   
    ind = (nrs*srd).nlargest(int(nrs*nr)).index
    
    src = [R + nrs*1 if i in ind else R for i,R in enumerate(sr)]

    return src

def get_colour_map(analysis="default",
                   palette="tol-muted",
                   continuous=False):
    """ Get colour mapping for analysis
    

    Parameters
    ----------
    analysis : str, optional
        Name of the analysis. The default is "default".
    palette : str, optional
        Colour palette to be used. The default is "tol-muted".
    continuous : bool, optional
        If a continous colour scale is used. The default is False.

    Returns
    -------
    ca : list, dict
        Colour mapping.

    """
    
    cp = pd.read_csv(config["fp"]["PATHR"]+config["fp"]["PATH_PP_COL_PAL"],
                     index_col=["PALETTE"])
    cp = cp.loc[palette,:]
    cp = cp.sort_values("PC_CODE")
    
    if not continuous:
        ca = pd.read_csv(config["fp"]["PATHR"]+config["fp"]["PATH_PP_COL_ALL"],
                         index_col=["ANALYSIS"]) 
        ca = ca.loc[analysis,:]
        ca = ca.merge(right=cp[["PC_CODE",
                                "COLOUR_CODE"]],
                      how="left",
                      on="PC_CODE")
        ca = ca.drop("PC_CODE",axis=1)
        ca = ca.set_index("ARTEFACT")
        ca = ca["COLOUR_CODE"].to_dict()
    if continuous:
        ca = list()
        ca.append([0,cp.loc[cp["PC_CODE"]==-1,"COLOUR_CODE"].squeeze()])
        colours = cp.loc[cp["PC_CODE"]!=-1,"COLOUR_CODE"].to_list()
        for s, c in zip(np.linspace(0.001,1,len(colours)), colours):
            ca.append([s, c])
            
    return ca


def explode_hh(df,disagg,col="PROPERTY_TYPE"):
    """ Disaggregate households
    
    Parameters
    ----------
    df : DataFrame
        Pandas DataFrame with a 'PROPERTY_TYPE' row or index level.
    disagg : list of str
        List of disaggregation items.
        
    Returns
    -------
    df : DataFrame
        Dataframe with exploded property type column.

    """
    
    if col not in df.columns:
        cols = df.columns
        df = df.reset_index()
        
    if col == "PROPERTY_TYPE":
        pp = ["Detached","Semi-detached","Terraced","Flats"]
        df_ex = df[df[col].str.contains('|'.join(pp))]
        
        df_ex.loc[:,"PROPERTY_TYPE_"]=df_ex["PROPERTY_TYPE"]+"|"+disagg[0]
        for d in disagg[1:]:
            df_ex.loc[:,"PROPERTY_TYPE_"] = df_ex["PROPERTY_TYPE_"]+","+df_ex["PROPERTY_TYPE"]+"|"+d
        df_ex.loc[:,"PROPERTY_TYPE"] = df_ex["PROPERTY_TYPE_"]
        df_ex = df_ex.drop(["PROPERTY_TYPE_"],axis=1)
    
        df_ex = df_ex.assign(PROPERTY_TYPE=df_ex["PROPERTY_TYPE"].str.split(",")).explode('PROPERTY_TYPE')
        df = pd.concat([df[~df[col].str.contains('|'.join(pp))],
                        df_ex],axis=0)
        
    elif col == "TECHNOLOGY":
        pt = ["DDDE00","DDSD00","DDTE00","DDFL00"]
        pf = ["DDDE","DDSD","DDTE","DDFL"]
        df_ex = df[df[col].str.contains('|'.join(pt))]
        
        df_ex.loc[:,col+"_"]=df_ex.loc[:,col].str[:-2]+disagg[0]
        for d in disagg[1:]:
            df_ex.loc[:,col+"_"] = df_ex[col+"_"]+","+df_ex[col].str[:-2]+d
        df_ex.loc[:,col] = df_ex[col+"_"]
        df_ex = df_ex.drop([col+"_"],axis=1)
        df_ex.loc[:,col] = df_ex[col].str.split(",")
        
        if "FUEL" in df.columns:
            bf = df_ex["FUEL"].str.contains('|'.join(pf))
            df_ex.loc[:,"FUEL"+"_"] = df_ex.loc[:,"FUEL"]
            df_ex.loc[bf,"FUEL"+"_"]=df_ex.loc[bf,"FUEL"+"_"]+disagg[0]
            df_ex.loc[~bf,"FUEL"+"_"]=(df_ex.loc[~bf,"FUEL"+"_"]+",").str.repeat(len(disagg)).str[:-1]
            for d in disagg[1:]:
                df_ex.loc[bf,"FUEL"+"_"] = (df_ex.loc[bf,"FUEL"+"_"]
                                            +","+df_ex.loc[bf,"FUEL"]+d)
            df_ex.loc[:,"FUEL"] = df_ex.loc[:,"FUEL"+"_"]
            df_ex = df_ex.drop(["FUEL"+"_"],axis=1)
            df_ex.loc[:,"FUEL"] = df_ex["FUEL"].str.split(",") 
            
        df_ex = df_ex.reset_index(drop=True)   
        df_ex = df_ex.apply(pd.Series.explode)
        #df_ex = df_ex.explode(col)
        
        df = pd.concat([df[~df[col].str.contains('|'.join(pt))],
                        df_ex],axis=0)
    elif col == "VALUE":
        pt = ["DDDE00","DDSD00","DDTE00","DDFL00"]
        pf = ["DDDE","DDSD","DDTE","DDFL"]
        df_ex_t = df[df[col].str.endswith(tuple(pt))]
        df_ex_f = df[df[col].str.endswith(tuple(pf))]
        
        df_ex_t.loc[:,col+"_"]=df_ex_t[col].str[:-2]+disagg[0]
        for d in disagg[1:]:
            df_ex_t.loc[:,col+"_"] = df_ex_t[col+"_"]+","+df_ex_t[col].str[:-2]+d
        df_ex_t.loc[:,col] = df_ex_t[col+"_"]
        df_ex_t = df_ex_t.drop([col+"_"],axis=1)
        df_ex_t.loc[:,col] = df_ex_t[col].str.split(",")
        df_ex_t = df_ex_t.explode(col)
        
        df_ex_f.loc[:,col+"_"]=df_ex_f[col]+disagg[0]
        for d in disagg[1:]:
            df_ex_f.loc[:,col+"_"] = df_ex_f[col+"_"]+","+df_ex_f[col]+d
        df_ex_f.loc[:,col] = df_ex_f[col+"_"]
        df_ex_f = df_ex_f.drop([col+"_"],axis=1)
        df_ex_f.loc[:,col] = df_ex_f[col].str.split(",")        
        df_ex_f = df_ex_f.explode(col)
        
        df = pd.concat([df[~df[col].str.endswith(tuple(pt+pf))],
                        df_ex_t,
                        df_ex_f],axis=0)
        
    if "cols" in locals():
        df = df.set_index([c for c in df.columns if c not in cols])
        
    return df


# The mock_snakemake function is largely adopted from PyPSA-EUR (Copyright
# 2017-2023 The PyPSA-Eur Authors) under an MIT licence.
from pathlib import Path
def mock_snakemake(rulename, **wildcards):
    """
    This function is expected to be executed from the 'scripts'-directory of '
    the snakemake project. It returns a snakemake.script.Snakemake object,
    based on the Snakefile.
    If a rule has wildcards, you have to specify them in **wildcards.
    Parameters
    ----------
    rulename: str
        name of the rule for which the snakemake object should be generated
    **wildcards:
        keyword arguments fixing the wildcards. Only necessary if wildcards are
        needed.
        
    Returns
    ----------
    snakemake : Snakemake
        Mock snakemake
    """
    import os

    import snakemake as sm
    from packaging.version import Version, parse
    from snakemake.script import Snakemake

    script_dir = Path(__file__).parent.resolve()
    assert (
        Path.cwd().resolve() == script_dir
    ), f"mock_snakemake has to be run from the repository scripts directory {script_dir}"
    os.chdir(script_dir)
    for p in sm.SNAKEFILE_CHOICES:
        if os.path.exists(p):
            snakefile = p
            break
    kwargs = dict(rerun_triggers=[]) if parse(sm.__version__) > Version("7.7.0") else {}
    workflow = sm.Workflow(snakefile, overwrite_configfiles=[], **kwargs)
    workflow.include(snakefile)
    workflow.global_resources = {}
    rule = workflow.get_rule(rulename)
    dag = sm.dag.DAG(workflow, rules=[rule])
    wc = dict(wildcards)
    print(wc)
    job = sm.jobs.Job(rule, dag, wc)

    def make_accessable(*ios):
        for io in ios:
            for i in range(len(io)):
                io[i] = os.path.abspath(io[i])

    make_accessable(job.input, job.output, job.log)
    snakemake = Snakemake(
        job.input,
        job.output,
        job.params,
        job.wildcards,
        job.threads,
        job.resources,
        job.log,
        job.dag.workflow.config,
        job.rule.name,
        None,
    )
    # create log and output dir if not existent
    for path in list(snakemake.log) + list(snakemake.output):
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    os.chdir(script_dir)
    return snakemake