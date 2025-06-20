
================
Workflow steps
================

This page describes in some detail the processing steps of the UK-MOSEM workflow, including the underlying assumptions. Yet, the processing setup is complex and not all elements are currently covered in detail here. As the model is further developed, this description could be extended and further detail provided. The workflow is implemented using Snakemake.

It is also important to note that the workflow configuration levers outlined below are not necessarily all robust to different use cases, i.e., consideration and adaptation of the relevant parts of the code can be necessary.

.. ***************
.. Overview
.. ***************

.. add over rulegraph and general explanation + model doc summary from thesis?
.. point out that sources are in meta data files, that source code is available and well documented?
.. note sources can be found with data


************************
Step-by-step description
************************



----------------------------------------------
Household projections (calc_hh_projections.py)
----------------------------------------------

This rule derives household projections as basis for building stock projections.

Main input data are:

* Household projections by ONS, National Records of Scotland, and StatsWales

Configuration options are:

* None

Outputs are:

* Number of households in each LA for each year up to 2060 (only for postprocessing/post-model analysis)
* Percentage change of the number of households in each LSOA for each year of the modelling period


Main calculation steps are:

* Simple calculation to derive percentage change from household projections


Main assumptions are:

* Projection (percentage change) for LAs are applied across LSOAs, i.e., are assumed to be equal across all LSOAs within one LA.
* Percentage change from 2044 (first year not covered by the household projections) onwards is based on the average of the 5 last years of data (2039-2043).


---------------------------------------------------------------
Non-domestic energy intensities (calc_nd_demand_intensities.py)
---------------------------------------------------------------

This rule derives demand intensities for the non-domestic building stock.

Main input data are:

* BEES energy consumption and intensity data for non-domestic buildings


Configuration options are:

* Scenario assumption for the future trajectory for non-heat gas (scen_nhg)  and electricity demand (scen_nhe)  – format: [Fraction of base year value]y[Year]a[Fraction of base year value]y[Year]a [...] - the calculation will interpolate (and back/forward fill) between these values.
* Local emission targets set through (scen_local_gov) are used to adjust non-heat gas intensity in the respective local authorities accordingly, i.e., reaching zero in the target year.


Outputs are:

* Energy intensity for space heat, hot water, non-heat gas, and non-heat electricity for different subsectors for each LSOA and each year


Main calculation steps are:

* The rule largely only loads and aggregates BEES energy intensity data, and applies those across all LSOAs and years (with separate assumptions on the trajectory of the use of non-heat end uses).
* For space and hot water heating, the final energy intensities are re-calculated as useful energy intensities using heating technology efficiencies.


Main assumptions are:

* National BEES data for England and Wales are applied to all local authorities, i.e., intensities are assumed to be the same across GB, including also LAs in Scotland.
* Heating energy intensity for new builds is assumed to decrease by the same percentage from the existing stock than an estimate for the same for domestic buildings based on NEED data for England and Wales (i.e., 30% reduction).
* The calculations assume average intensities across all buildings within each non-domestic sector.
* The calculations assume no changes in intensities, except for non-heat end uses.
* Non-electricity/non-gas fuels are not considered but are only a small fraction of overall final energy consumption for buildings.
* BEES and in turn this rule only considers consumption for building services, not industrial processes.
* BEES intensities for space heating and hot water are adjusted to useful energy using the efficiency of installed heating technologies. This uses currently a simple average across non-domestic sectors and LSOAs, weighted based on floor space.


------------------------------------------------------------
Domestic property projections (calc_dw_stock_projections.py)
------------------------------------------------------------

This rule derives domestic building stock projections.

Main input data are:

* EPC data
* VOA data and Scottish dwelling data
* Data on demolitions
* English Housing Survey data on tenure by property type and income group (used if tenure/income lever activated)
* household projections (output from previous rule)

Configuration options are

* Tenure and income lever (scen_hh_disagg) – if to integrate tenure, or tenure and income in property type (and, hence, it being separately processed throughout the input data processing) – format: "-" if no disaggregation, "T" for disaggregation based on tenure, "TI|[1-9DIGITS]" for disaggregation based on tenure and income, where the [1-9DIGITS] is a number of digits that represent the deciles that income groups will be split into (e.g., 124 means 4 income groups will be represented: -1, 1-2, 2-4, 4-).


Outputs are:

* Total number of different tenure instances and building configurations (heating technologies, building efficiency) for each LSOA and property type for each of the base years
* Projections for domestic properties (existing and new build) for each property type and LSOA for all model years

Main calculations are:

* Use data on dwellings, household projections, and demolitions to derive property projections (if tenure/income are included this is also integrated in the projections, see below)
* Process EPC data, including scaling based on VOA data (if tenure is included, scaling also takes into account national data on tenure by property type, if income and tenure is included it takes into account a local income distribution and national data on income across tenure and property types), to create a data set on socioeconomic and building characteristics for properties.



Main assumptions are:

EPC data related -

* The processing is using the last available certificate per property (even if that is then applied to a previous base year in model).
* Space constrained homes are assumed to be properties with less than 16 m^2 per room (calculated simply from total square meter divided by number of habitable rooms).
* Property types are aggregated to 4 different categories (bungalows to semi-detached based on similar NEED consumption; no data to terraced).
* Heating and heat distribution technologies are processed and aggregated to technologies available in the model (combined communal and district heating, WSHP to GSHP, etc.).
* Other building characteristics are processed and aggregated to match the structure in model.
* When processing tenure data, owner occupied is assumed as default if allocation to private or social rent is not possible.
* Socioeconomic and building characteristics are saved separately for each variable, not as combination (using 2022, except 2019 for building efficency related data, to align with other data).
* Scaling of data based on VOA and EHS data (if tenure/income lever activated) assumes EHS data (fractions) apply to entire GB.
* If tenure/income is integrated, LSOA income data is applied across property types and taking into account national data on income by tenure, through iterative scaling.

VOA data related -

* Property types are aggregated to 4 different categories (bungalows to semi-detached based on similar NEED consumption; unknown distributed across; vans, boats, excluded)
* Scottish data for 2018 to 2022 are scaled 2017 data using totals (not disaggregated by property type) for 2018-2022 due to a lack of available data at the time of implementation.
* The average demolition rate is calculated based on the 2012-2018 period, assumed the same across LSOAs in each LA, and constant over time - if the decrease of property numbers based on the decrease of household number is larger, this number is used as demolition rate.


----------------------------------------------------------------
Non-domestic property projections (calc_nd_stock_projections.py)
----------------------------------------------------------------

This rule derives non-domestic building stock projections.

Main input data are:

* EPC/DEC data
* VOA data and Scottish data on non-domestic floor space
* Non-domestic floor space data for England and Wales from BEES data

Configuration options are

* None


Outputs are:

* Floor space with certain building configurations (heating technologies and building heat distribution) for each LSOA and property type for each of the base years
* Projections for floor space for non-domestic properties (existing and new build) for each property type/non-domestic sector and LSOA for all model years.

Main calculations are:

* EPC data and VOA data (for categories available) are used to derive floor space data for different non-domestic categories and heating characteristics in each LSOA for base years (2014 - England and Wales, 2017 - Scotland).
* These data are scaled to match national totals for each non-domestic category in Scotland, and England & Wales (EW).
* This is then extended to other historical years based on average floor space development in EW, and projected assuming demolition rates and overall change based on household projection.
* An uplift factor is applied to take into account the scope of BEES and the Scottish analysis of non-domestic floor space.


Main assumptions are:

* The processing is using the last available certificate per property/building (even if that is then applied to a previous base year in model) – also as no build date is provided in the certificates.
* The floor space is split equally across different property types if more than one is given for a certificate.
* Property types are aggregated to different categories based on BEES categories.
* Heating and heat distribution technologies are processed and aggregated to technologies available in the model (e.g., combining communal and district heating). Heating technologies for electricity based heating are processed based on simple assumptions (e.g., if air-conditioned, assume heat pump, otherwise electric resistance). If no data on heating technology is given, gas boiler are assumed as the default. The model currently only includes wet distribution systems.
* The uplift factor simply lifts up floor space for all sectors across GB based on the estimated scope (90% of floor area covered) of BEES, which is assumed to be also applicable to the Scottish analysis.


--------------------------------------------------------------------------------------------------
Residual capacity of building heat technologies - Fractions (calc_residual_fractions_heat_tech.py)
--------------------------------------------------------------------------------------------------

This rule derives residual capacity fractions of heating technologies.

Main input data are:

* Data on domestic and non-domestic heating systems from EPC data (output from previous rule)

Configuration options are:

* None

Outputs are:

* Fractions of installed heat technologies, including building heat distribution


Main calculations are:

* This calculates the fraction of properties (domestic)/floorspace (non-domestic) that is heated with a certain technology in the base period.


Main assumptions are:

* This uses data for one base year (2022) as the underlying EPC data are not used in a way that differentiates between base years.
* If no heating technology data are available for a property type in an LSOA (EPC data), the average value of the property type in the LAD is used.

----------------------------------------------------------------------------------------------------
Residual capacity of building heat technologies - Capacities (calc_residual_capacities_heat_tech.py)
----------------------------------------------------------------------------------------------------

This rule derives residual capacities of heating technologies.

Main input data are:

* Fractions of installed heat technologies, including wet heating systems, for the base years
* Peak heat demand/unit sized for each LSOA, property type, year

Configuration options are

* None


Outputs are

* Capacities of installed heat technologies, including wet heating systems, for base years (2015-2022)


Main calculations are

* This calculates capacities by multiplying the fractions with peak demand capacities for each base year (2015-2022)


Main assumptions are

* As before, this assumes constant heating technology residual fractions across all base years as no suitable data set that could capture such changes within base years is available.
* The residual capacity for 2022 is adjusted to ensure that an appropriate level of residual capacity will be projected for 2023 and as the base year values vary due to changing demand (2022 has low demand and would result in less than expected residual capacity in 2023 when demand is again higher).


-----------
Demand intensity for domestic properties (calc_dw_demand_intensities.py)
-----------

This rule derives demand intensities for different domestic property types.

Main input data are:

* NEED data for EW
* Heating technology efficiencies 
* EPC data on energy efficiency bands of properties in each LSOA (output from previous rule)

Configuration options are:

* Scenario assumption for the future trajectory for non-heat gas (scen_nhg)  and electricity demand (scen_nhe)  – format: [Fraction of base year value]y[Year]a[Fraction of base year value]y[Year]a [...] - the calculation will interpolate (and back/forward fill) between these values.
* Local emission targets set through (scen_local_gov) are used to adjust non-heat gas intensity in the respective local authorities accordingly, i.e., reaching zero by the target year.

Outputs are:

* demand intensities for each end use for each property type and LSOA and all model years


Main calculations are:

* NEED record level data are processed to calculate consumption intensities for each property type and energy efficiency band for each region.
* Consumption intensities, together with heating efficiency data and data on energy efficiency bands of properties in each LSOA, are used to calculate demand intensities.


Main assumptions are:

* The calculations assume the intensity is the same across a region for each energy efficiency band, and the same across Scotland as a whole (using value for the North East).
* As above, bungalows are considered as semi-detached properties.
* The calculation uses the model gas boiler efficiency to calculate useful energy demand from gas consumption 
* The calculation assumes non-heat gas consumption only exists in properties that have a gas boiler installed in the base years.
* The calculations assume hot water and non-heat gas demand are independent of efficiency band and are calculated based on averages across efficiency bands.
* For new builds, the demand intensities are assumed equal, except the demand for space heating, which is assumed to be the same as for properties in the A/B efficiency band.


--------------------------------------------------------------------------------
Annual demands for domestic and non-domestic properties (calc_annual_demands.py)
--------------------------------------------------------------------------------

This rule derives annual energy demands for domestic and non-domestic properties for each local area.

Main input data are:

* Intensities (output from previous rule)
* Property projections (output from previous rule)
* Subnational gas and electricity consumption statistics – if MSOA calibration triggered
* Heating technology residual fractions – if MSOA calibration triggered (output from previous rule)
* ECUK data – if GB or MSOA calibration triggered
* Heating technology efficiencies – if GB or MSOA calibration triggered

Configuration options are:

* Parameter defining if and what calibration of demands is performed (scen_dem_calib) – syntax: "-" for none, "GB" for calibration based on ECUK data, or "MSOA" for calibration based on subnational energy consumption statistics


Outputs are:

* Annual demands (SH, HW, NHE, NHG) for each property type, including non-domestic, and LSOA for all years


Main calculations are:

* The annual demands are derived by multiplying intensities with property numbers (domestic) and floor area (non-domestic) projections.
* For GB calibration: A calibration based on ECUK consumption data for different end uses is performed.
* For MSOA calibration: A calibration based on LSOA/MSOA consumption data is performed (experimental). 


Main assumptions are:

* GB calibration: This applies a GDP/population ratio to deduct Northern Ireland from ECUK demands, and uses heating technology stock and efficiencies from the model to convert ECUK data to useful energy for heating demands. It assumes a average scaling factor of past years with data for future years.
* MSOA calibration (experimental): see source code for details.

--------------------------------------------------------------------------------------
Peak demands/unit size for domestic and non-domestic properties (calc_peak_demands.py)
--------------------------------------------------------------------------------------

This rule derives (before diversity) peak heat demand per LSOA, property type, and technology

Main input data are:

* Annual demands (output from previous rule)
* Property projections (output from previous rule)
* Heat technology load factors


Configuration options are:
* None


Outputs are:

* Unit size/peak heat demand (before diversity) for each LSOA and property type (total and per property)


Main calculations are:

* This uses load factors (specifically derived by another analysis for this purpose) to calculate unit size/peak heat demand based on annual space heat demand.


Main assumptions are:

* Load factors for technologies that are not part of the other analysis are set to the one for oil boilers (central value) - this is largely irrelevant and mainly influences, e.g., how capacities are scaled in the results.
* These calculations are for each year, and consider both existing and new builds. Hence, they calculate an average of those.
* An average floor space per non-domestic property is assumed across all LSOAs.
* It is assumed the load factors for sizing also apply to non-domestic properties.



-------------------------------------------------------------------------------
Timeseries for capacity factors and temperature (calc_reanalysis_timeseries.py)
-------------------------------------------------------------------------------

This rule derives timeseries from reanalysis ERA5 data. 

Main input data are:

* ERA5 climate data
* Shape files for UK, LADs, and designated offshore wind areas


Configuration options are:

* None


Outputs are:

* Hourly timeseries for capacity factors (onshore wind, offshore wind, pv) for the UK, and ambient temperature for each LAD for 1999


Main calculations are:

* ERA5 reanalysis climate data are used to derive capacity factors and temperature for the gerographies using the atlite package.


Main assumptions are

* Capacities are built equally across the relevant geopgraphies (PV: land area, wind onshore: land area, wind offshore: designated offshore areas) without considering factors like unsuitable areas (these are taken into account in terms of potentials).
* This only uses data for the year 1999, no multi-year period.

-------------------------------------------------------
Temperature timeseries (calc_temperature_timeseries.py)
-------------------------------------------------------

This rule derives climate/temperature timeseries for each local authority from MetOffice data.

Main input data are:

* MetOffice Local climate projections
* Shape files for LADs

Configuration options are:

* None

Outputs are:

* Hourly timeseries for ambient temperature for each LAD for 1999.


Main calculations are:

* This uses MetOffice data to derive ambient temperature for each LAD, appling a binary mask to gridded data based on LAD shapes.


Main assumptions are:

* This is using MetOffice local projections data as this could allow using projected temperature patterns in future (not currently used).


---------------------------------------------
Demand timeseries (calc_demand_timeseries.py)
---------------------------------------------

This rule derives hourly demand timeseries.

Main input data are:

* BDEW demand profiles (through demandlib, see demandlib docs)
* Temperature timeseries (output from previous rule)


Configuration options are:

* None

Outputs are:

* Hourly timeseries for all demands (SH, HW, NHE, NHG), property types (including non-domestic), and all LADs.


Main calculations are:

* This calculation uses the demandlib library to generate demand profiles taking into account temperature timeseries.


Main assumptions are:

* This applies German standard demand profiles from BDEW, integrated in demandlib (for further demandlib config, see code, and refer to the demandlib documentation).
* The NHE demand profile might include some electric heating, although larger electric consumers, e.g., heat pump are not considered for the standard profile (see BDEW documentation).
* The NHG profile is assumed to be constant across time.
* Profiles are based on 1999 temperature data, and are assumed constant over the modelling horizon.

---------------------------------------------------
Aggregate timeseries (calc_aggregate_timeseries.py)
---------------------------------------------------

This rule derives aggregated timeseries, i.e., time slice values and lengths.

Main input data are:

* Relevant timeseries for demands and capacity factors for power generation technologies (output from previous rule)

Configuration options are:

* Aggregation approach (scen_time_agg) to define how timeseries are to be aggregated – format: NP|S where N is one or more digits of an integer defining the number of typical periods, P is either h (hour), d (day), p(day, predefined order) defining the length of the typical period, and S is one or more digits of an integer defining the number of segments within a typical period.


Outputs are:

* Aggregated timeseries for demands and capacity factors in the form of the fraction of demand in each of the timeslices and the fraction of the year each timeslice represents


Main calculations are:

* Using the tsam library to cluster the timeseries based on a k-means algorithm and following the chosen configuration.


Main assumptions are:

* The clustering assumes a weighting across the different timeseries, where (almost) all weight is giving to the space and hot water demand curve and only neglible weight to other timeseries.
* The peak heat demand is specifically added as an additional cluster center.
* The resulting GB peak demands are slightly lower than the original based on the timeseries. There is a question if to capture UK peak correctly or LA peaks correctly – which can be different due to the structure of the model (e.g., timeslices are the same across LAs, etc.).


-----------------------------------
Road lengths (calc_road_lengths.py)
-----------------------------------

This rule derives road lengths per LSOA.

Main input data are:

* OS openroads - GIS data of the road network in GB from 2020
* Geospatial data of LSOA boundaries

Configuration options are:

* None


Outputs are:

* Total road lengths in each LSOA

Main calculations are:

* This calculates the intersection between the road network layer and LSOA boundaries to calculate the road lengths in each LSOA.


Main assumptions are:

* This includes all road types captured by OS openroads.
* Roads segments that cannot directly be match with an LSOA are iteratively added to LSOAs if their start-/endpoint is part of a segment that is already allocated to an LSOA.


--------------------------------------------
Linear heat density (calc_demand_density.py)
--------------------------------------------

This rule derives linear heat demand densities for each LSOA.

Main input data are:

* Annual demand for space heating and hot water for domestic and non-domestic properties (output from previous rule)
* Total road length for each LSOA (output from previous rule)



Configuration options are:

* None


Outputs are:

* Linear heat density for each LSOA

Main calculations are:

* This calculates the linear heat density by dividing the total heat demand per LSOA by the total road length in each LSOA.


Main assumptions are:

* This assumes the road length as a proxy for the length of a heat network - this is a simplification but a common assumption.


----------------------------------------------
Classification of LSOAs (calc_sublocal_agg.py)
----------------------------------------------

This rule derives a classification for sublocal areas.

Input data are

* Linear heat density of each LSOA (output from previous rule)
* Thresholds of linear heat density for classes



Configuration options are

* None


Outputs are

* Class allocation for each LSOA

Main calculations are

* Each LSOA is allocated to a class of LSOAs based on its heat density


Main assumptions are

* The heat density thresholds are currently defined at the 0.7, 0.8, and 0.9 quantile of the energy demand sorted by heat density, i.e., there are 4 different classes (but this can be updated in the input data set)


------------------------------------------------------------------------
Distribution networks and DH generation residual (calc_dist_networks.py)
------------------------------------------------------------------------

This rule derives techno-economic parameters for distribution networks, as well as residual capacity of DH residual generation capacity (currently linked to network capacity).

Main input data are:

* Techno-economic parameters for distribution networks and other data from the input data set (average floor space per non-domestic building, current gas distribution network length, efficiencies of heat technologies etc.) 
* Road lengths per LSOA (output from previous rule)
* Annual demands for all end-uses (output from previous rule)
* Timeslice demand fractions and lengths (output from previous rule)
* LSOA class allocation (output from previous rule)

Configuration options are:

* None


Outputs are:

* Capital and fixed cost for distribution network technologies for each sublocal area (constant for all years)
* Efficiency of DH networks (i.e., 1 - losses)
* Residual capacities for DH generation

Main calculations for costs are:

For natural gas, H2 retrofit, and district heating network –

* First the total cost for installing the respective network in the sublocal area is calculated by multiply a cost per meter length with the road length in the respective area. The road length for the gas network is scaled based on data for the current network length and residual capacity. For district heating, the building/property connection cost to all buildings is added.
*  To calculate the cost per GW of installed network capacity, the total cost for each network is then divided by the respective peak demand (for DH the SHW peak, for gas and H2 network, the SHW + NHG demand - taking into account the respective building heating technology efficiency for the SHW peak to get to the final energy demand peak).

For electricity –

* The average reinforcement/replacement cost per kW are taken from the literature.

For all –

* Annual fixed cost are calculated as 0.1% of respective capital costs.

Main assumptions are:

* The current implementation aggregates the connection cost, i.e., connection cost that are actually different for different property types (including domestic/non-domestic) are averaged within the model and applied across all property types.
* The cost calculation for DH assumes all properties along the network will be connected and use the network.
* The way the cost are calculated implies several simplifications: 1) the areas are assumed to be homogeneous, single nodes 2) if building retrofit decreases the peak demand, this also reduces the cost for building the networks (this is generally not a large fraction and it might to some extent capture a realistic effect – through using small, cheaper pipes), 3) for DH and gas, the GW capacity of network is not actually the capacity of the grid but a measure for the extent of the network to meet the respective demand.
* For power distribution networks, it is assumed it is already in place and connected to all properties.


Main calculations for efficiency are:

* Total annual losses are calculated by multipling the length of the network (road length + average internal pipe length times number of properties) with the heat loss factor per length of network.
* Relative losses (1-efficiency) are calculated by dividing the total annual losses by the total annual heat demand.

Main assumptions are:

* The network temperature is assumed to be 70°C with a 50°C return temperature.
* This means the relative efficiency/losses of a heat network in the model are independent of the operation of a heat network.

Main calculations for capacity factors are:

* Capacity factors are calculated for the DH network to avoid heat being supplied to other areas in non-peak timeslices. These are calculated as heat demand (power) in each sublocal area, timeslice, and year divided by the peak heat demand in each respective sublocal area and year.
* For electricity this is not necessary as it is assumed all properties are connected. For gas it is not possible to calculate this with the chosen approach as the shape of the demand curve depends on the uptake of gas for heating (given there is a constant NHG demand that is added). Hence this is not implemented but has little impact as to use gas from the network in non-peak timeslice the model would need to build additional capacity that is not useful in the peak timeslice.

Main assumptions are:

* 

Main calculations for residual capacities are:

* Residual capacities for district heating networks is calculated as network capacity that is required to supply heat to the residual capacity of heat interface units.
* Residual capacities for gas networks is calculated as network capacity that is required to supply gas to meet NHG demand and what is required by the residual capacity of gas boilers.
* Residual capacities for the electricity network are calculated as network capacity that is required to supply power to meet NHE demand and what is require by the residual capacity of electricity-based heating. A spare factor is applied to take into account the GB-average relative spare capacity.
* Residual capacity for DH generation is calculated based on EPC data on heat supply to HIUMs and the peak as calculated for the DH network

Main assumptions are:

* Residual capacity are calculate based on above until the year 2022, from when they linearly decrease (DH) or stay constant until the year when they start decreasing (EL/GA).
* A residual capacity for H2 retrofit is added only as a modelling approach to allow for mixing of hydrogen in the existing grid.
* For DH generation residual capacity, if no information given Gas CHP is assumed.

-----------------------------------------------
Conservation areas (calc_conservation_areas.py)
-----------------------------------------------

This rule loads and processes conservation areas.

Main input data are:

* Spatial data on conservation areas in England, Wales, and Scotland
* Spatial data on LSOA borders
* Spatial data on properties in Great Britain


Configuration options are:

* None


Outputs are:

* Fraction of properties in each LSOA that are part of a conservation area


Main calculations are:

* The GIS data on LSOAs and conservation areas are processed to derive the parts of LSOAs protected by conservation areas, if any.
* The number of properties in each of those parts is then divided by the total number of properties in the LSOA to get to the fraction of properties protected by conservation areas.


Main assumptions are:

* The approach does not differentiate between property types (including between domestic and non-domestic) and does only derive a generic fraction for each LSOA (in the current model version this is irrelevant as each LSOA is either completely or not at all covered by a conservation area).
* No change in conversation areas going forward are assumed.



------------------------------------
Energy supply (calc_supply_techs.py)
------------------------------------

The rule derives techno-economic parameters for supply technologies, including transmission.

Main input data are:

* Oil, gas, coal price projections from the Future Energy Scenarios (FES)
* Heating oil price statistics
* Capacities of power generation technologies from FES
* Power and hydrogen generation technology cost data from BEIS and other sources
* Emission factors from BEIS (now DESNZ)
* Residual power sector capacitity from DUKES and BEIS renewable energy per local authority
* Data on renewable potentials
* Capacity factors for renewable technologies (output from previous rule)

Configuration options are

* None ("scen_supply_imp" has been removed but is still present in some parts for potential reimplementation).


Outputs are:

* Fossil fuel price projections
* Characterization of power and hydrogen generation technologies
* Characterization of transmission technologies

Main calculations for import/extraction of fossil fuels are:

* The price for fossil fuels is set based on projections from FES. There is no differentiation between import and extraction of local fossil fuels. Crude oil price is translated to heating oil prices using heating oil price statistics.

Main calculations and assumptions for the power sector:

* The power sector only includes a set of core technologies (nuclear power plants, dam hydro, natural gas-based combined cycle gas turbine power plant, coal power plants, hydrogen-based combined cycle gas turbine power plant, dedicated biomass power plant, utility solar PV, roof-top solar PV, onshore wind, offshore wind, and battery storage).
* The techno-economic data for power generation technologies are mainly based on BEIS (now DESNZ) electricity generation costs assumptions (assuming 'medium' cost values) and where not part of this dataset other sources.
* The characteristics for H2 gas power plants are assumed to be the same as a normal CCGT power plant.
* The entire power sector is modeled at the national level, i.e., even rooftop solar feeds through the transmission grid.
* Natural gas and hydrogen power plants are assumed as 'CCGT H Class', biomass as 'Dedicated Biomass', hydro as 'Hydro 516MW', onshore wind as 'Onshore Wind', offshore wind as 'Offshore Wind', open-field solar as 'Large-Scale Solar', and rooftop solar as 'Solar PV 10-50kW'.
* Emissions factors are taken from DESNZ/BEIS data. These are added to the transmission grid technologies to allow for attribution to local authorities.
* Capacity factor for VRE are based on previous calculations.
* Residual capacities for renewables are taken from BEIS statistics, for other technologies residual capacities are derived from the DUKES list of power plants. This is deriving residual capacities for each LAD, but these are currently aggregated. The data only account for major producers for nuclear and fossil fuel plants, but that only leaves a relatively small amount that is not covered.
* Residual capacity in the dataset is only given for PV, which is all allocated to utility solar. This does not make much of a difference given both are currently treated similarly in the model.
* Variable renewable potentials are applied for solar and wind based on an existing, spatial analysis. No additional hydro potential is assumed.
* A simple storage representation is implemented by assuming production and use of stored fuels (electricity in batteries and hydrogen) only needs to be balanced over an entire year, not in each timeslice.
* The development of new capacities in the power sector is constrained to follow the 'Leading the Way' scenario of the FES. The fraction of new capacities (this only takes into account technologies relevant to the power sector going forward, i.e., zero-emission technologies, while some technologies not represented in the model are added to others for the purpose of calculating the fractions).


Main calculations and assumptions for the transmission technologies:

* Natural gas, Hydrogen, electricity, heating oil, biomass transmission are implemented as technologies that transmit the respective energy carrier from the national to the local level, with no transmission directly between local authorities possible (currently largely irrelevant given the national level supply sector).
* Natural gas and H2 transmission cost are calculated by multiplying a per length cost with the current network length divided by the base year total capacity.
* Residual electricity transmission grid is assumed to be equal to distribution grid capacity without spare capacity. Gas transmission grid is also assumed to be equal to distribution grid capacity in the base year with a small uplift to avoid any issues (capacity here is again less about what energy can be transported as it is not expected any of these grid will need to be extended, but potentially replaced). The decomissioning of the existing grid over time is set based on the literature.


Main calculations and assumptions for the hydrogen sector:

* The data for hydrogen generation technologies are based on BEIS Hydrogen Production Costs assumptions. The model only includes steam methane reforming with CCUS and alkaline electrolysis.
* For technical modelling reasons, two electrolysis technologies exists. One is directly providing hydrogen, the other provides hydrogen to be stored.
* It is assumed there is no residual capacity.

------------------------------------------------------------
Building efficiency measures (calc_building_improvements.py)
------------------------------------------------------------

This rule derives techno-economic parameters for building efficiency measures.

Main input data are:

* Characteristics of the domestic building stock with respect to windows, walls, roof, and floor type for each property type and LSOA (output from previous rule)
* Number of properties in a conservation area in each LSOA (output from previous rule)
* Total number of each of the potential efficiency measures
* Techno-economic parameters of efficiency measures
* Peak space heat demand per property type and LSOA (output from previous rule)

* Total cost and reduction in final demand for space heating for efficiency measures in non-domestic buildings
* Non-domestic property stock for each LSOA (output from previous rule)

Configuration options are:
* None


Outputs are:

* Cost, capacity factors, capacity constraint of three domestic and one non-domestic energy efficiency packages for each LSOA aggregation and property type

Main calculations for domestic properties are:

* This uses property characteristics from the property projections (EPC data) to derive the theoretical potential for building efficiency measures (without considering that measures might not be possible to combine) for each property type and LSOA
* It applies a reduction based on conservation areas reducing the relevant numbers of theoretically possible measures in respective LSOAs.
* The actual possible number of each of the measures is loaded as national total and allocated across properties based on above calculation 
* The cost per measure for each property type and the savings potential for each measure and property type (taking into account loft and in-use factor), the maximal annual installations constraint, and the lifetime for each of the measures are loaded and processed.
* Measures are then aggregated to three different efficiency packages (low, medium, high).


Main assumptions for domestic properties are

* Measures that are not covered by data on maximal installation per year are set to the average of all measures with provided data (i.e., average of percentage of total potential measures that can be installed per year)
* The cost is calculated as cost per measure times the number of measures divided by the product of peak space heat demand and efficiency savings (percent of total demand for a property type in an LSOA).
* Where the number of measures calculated this way exceeds the number of properties, the measures will be reallocated across other LSOAs and property types.

Main calculations for non-domestic properties are

* The calculation is similar as for domestic properties except 1) only one generic efficiency measures exists, 2) no spatial data on distribution of potential exists, so this is disagreggated simply based on floor space (taking into account conservation areas)
* The cost is calculated as cost per measure times the number of measures divided by product of peak space heat demand and efficiency savings (percent of total demand for a property type in a LSOA)


Main assumptions for non-domestic properties are

* Given there is no spatial data this approach simply assumes the potential for efficiency measures/demand reduction per floor space is spread equally across LSOAs (taking into account conservation areas).
* An uplift factor is applied to take into account, among others, new buildings that have been built since the underlying analysis has been undertaken.
* A constraint on the maximal annual installation is set as average rate of domestic measures.
* The total potential for efficiency measures is scaled up based on floor space to also capture Scotland.
* This only takes into account the potential of building fabric measures from the BEES data for now (also as it is difficult to allocate some of the other measures to space heat demand).


---------------------------------------------------------------
Constraints for building heat technologies (calc_heat_techs.py)
---------------------------------------------------------------

This rule derives constraints for building heating technologies.

Main input data are:


* Rural-urban classification for all LSOAs in GB
* Heritage and space-constrained suitability fractions for relevant technologies
* Conservation area data (output from previous rule)
* Space constraint data (output from previous rule)
* Fraction of residual capacities of heat technologies (output from previous rule)
* Heat demand peaks (output from previous rule)

Configuration options are:

* Scenario heat technology deployment constraint ("scen_htd_con") format: T-S-Y-L-R-N, where T is a string that is used to filter technologies the constraint applies to (e.g, "HPD" will constrain all technologies that include "HPD" – all domestic heat pumps), S is the sector (either "D" for domestic, "N" for non-domestic, or "ND" for both – this is only used to pick the annual demand if relevant for the constraint), "Y" is the type of constraint (either "ci" for capacity investment, "ct" for total capacity, or "a" for activity), "L" defines if an upper ("u") or lower ("l") is set, "R" defines the region the constrained is applied to (":*" is used for all regions, several regions can be separate by ";" to aggregate or "," for separate constraints for each region), "N", defines the actual limit in certain years in the format of [Limit]y[Year]a[Limit]y[Year]a[...], e.g., "200000y2023a200000y2050a600000y2060". For capacity constraints the number refers to the number of installations, for activity constraints to the fraction of demand.

Outputs are:

* Total installed capacity constraints for the relevant technologies for all years and LSOAs
* Capacity investment constraint for historic years
* Scenario deployment/activity constraint for specific technologies, years, and regions

Main calculations are:

* The constraint on biomass boilers is calculated by multiplying the peak demand (for biomass boilers) with the urban rural classication (urban: 0, rural: 1), the offgrid (i.e., non-gas boiler) fraction, and the suitability fractions.
* The constraint on heat pumps due to heritage considerations is calculated by multiplying the peak demand with the fraction of properties within the LSOA in conservation areas, and the relevant suitability fractions.
* The constraint on heat pumps due to space constraints is calculated by multiplying the peak demand with the fraction of properties within the LSOA that are space constrained, the relevant suitability fractions, and the fraction of properties in the LSOA not in conservation areas (to avoid double counting).
* The constraint on heat pumps are then added up and subtracted from the total peak demand to get to the constraint.
* A constraint for historical years is calculated that sets the new capacity for all technologies to zero for historic years (except oil boilers that are used to balance any mismatch)
* The flexible constraint is calculated based on user input explained above.

Main assumptions are

* Biomass boilers are assumed to be only suitable in properties in rural areas off the gas grid (i.e., without current gas boiler), with a suitability fraction applied to account for suitability based on space and other requirements.
* A simplified approach is used to avoid double counting heritage and space constraints for heat pumps (see above). This is likely increasing the allowed capacity of HPs properties in conservation areas that would still be suitable to have HPs installed in the model but might actually be space constrained.

-----------------------------------
Biomass supply (calc_potentials.py)
-----------------------------------

This rule derives local renewable potentials and related cost for biomass.

Main input data are:

* LAD boundaries
* Spatial landcover data
* Biomass potential and cost data (NUTS2)
* Pellet production cost

Configuration options are

* None

Outputs are

* Biomass potential and supply cost for each local authority

Main calculations are:

* Biomass potentials are disaggregated to LAD-level by using land cover data.
* The supply cost for each NUTS2 area are similary disaggregated and cost for pellet production cost is added to derive overall pellet supply cost.

Main assumptions are:

* This only considering agricultural residues and forest residues/products, it does not include energy crops.
* Splitting from NUTS2 to LAD relies on forest cover for forest residues and agricultural land for agricultural residues, splitting potential relative to land area for the respective type of land.


-----------------------------------------------------------------
Heat sources for district heating networks (calc_heat_sources.py)
-----------------------------------------------------------------

This rule derives potentials for DH heat sources.

Main input data are:

* Data on excess heat potential (including location)
* LSOA boundaries
* Population projections
* Waste water production per person and temperature difference if used with heat pump
* Efficiency of waster water heat pump

Configuration options are:

* None


Outputs are:

* Potential for excess heat recovery for each LSOA
* Potential for/maximum generation of heat by waste water heat pumps for each LAD

Main calculations are:

* For excess heat, the raw data are loaded and then aggregated to LADs.
* For waste water, the potential is calculated by multiplying the average waste water production per person with the population in each LAD, the temperature difference of the water before/after the heat pump, and a factor to include the electricity that also heats up the water (given the constraint is with respect to the total output).

Main assumptions are:

* Both are applied as totals over the year with no temporal resolution.


------------------------------------
Local governance (calc_local_gov.py)
------------------------------------

This rule derives local governance related parameters/constraints.

Main input data are:

* Local climate pledges
* Local climate plans


Configuration options are:

* Switch to use or not use climate pledges ("scen_local_gov") – format: "-" for no targets, "CA" for targets clipped at 2040 and only in LAs with building-related plan, as well as no DH investments in LAs without building-related plan, "CP" for the targets as in the data, "CPS" as "CP" but additionally adjusts emission constraint to ensure freed up emission budget is not used by others.


Outputs are:

* Emission constraint on local emissions
* Capacity constraint on district heating


Main calculations are:

* This uses the target year from pledges (for whole area, not just council operations) to derive to emission constraint (0 for target year and after).
* For the capacity constraint, information on existing building-related plans is used to constrain district heating investement to 0 if no plan is in place in a local authority.

Main assumptions are:

* This assumes net-zero pledges means no emission from the building sector.


-----------------------------------------------
Model input data set (create_input_data_set.py)
-----------------------------------------------

This rule arranges the OSeMOSYS input data set based on raw and processed data.

Input data are

* Most of the outputs from previous rules


Configuration options are

* Aggregation of years ("scen_year_agg") – format: '5y' for 5 year periods (except initial and final year) or '-' for no aggregation
* Adjustment of capital cost for scenarios ("scen_tech_cost") – format: '-' for no changes, T-Y-M-R-N, where T is a string that is used to filter technologies the adjustment applies to (e.g, "HPD" will affect all technologies that include "HPD" – all domestic heat pumps), "Y" is the type of constraint (currently only 'C' for capital cost is implemented), "M" defines if changes are addition ("A") or multiplication ("M"), "R" defines the region the adjustment is applied to (":*" is used for all regions, several regions can be separate by ";" to aggregate or "," for separate constraints for each region), "N", defines the actual adjustment value in certain years in the format of [Limit]y[Year]a[Limit]y[Year]a[...], e.g., "200000y2023a200000y2050a600000y2060".
* Introducing technology bans ("scen_tech_bans") – format: '-' for none, otherwise T-R-Y|[...], where T is a string to filter for the technologies the ban should be applied to, R are the geographic entities it should be applied to, and Y is the year the ban is implemented.
* Maximum capacity constraint ("scen_mcap_con")

Outputs are:

* Model input data set to run with the multi-scale framework fratoo


Main calculations are

* This script mainly performs minor processing and restructuring of data to save in the format of an OSeMOSYS/fratoo model file.
* The script applies temporal aggregation of years depending on the scenario config.
* The script adjust the base year of monetary values if necessary.
* The script interpolates data across model years if necessary (flat before first and after last value)
* The scipt calculates residual capacity going forward by assuming a equal age distribution of existing capacity and thus a linear decrease of capacity in line with each technologies' lifetime.
* The script calculates capacity factors for heat technologies to ensure accurate capacities for heat technologies need to be built (taking into account sizing factors and peak demand in the model) and proportinally constraining operation in non-peak timeslices to avoid technologies with cheaper running cost supplying beyond the buildings they are installed in.
* The script adjusts the national emission constraint to take into account local targets, if applied, to ensure faster emission reduction in some local authorities does not free up emissions for others to go slower (if triggered).


* Input and output activity ratio are calculated based on the efficiency - where generally the input activity ration is 1/efficiency and the output activity ratio is 1 (see code for details).
* Capacity factor for building heat technologies are calculated based on LSOA peak demand before diversity (required capacity) and after peak diversity (model heat demand).

Main assumptions are

* -


------------------
Run (run_model.py)
------------------

This rule runs the model.

Main input data are:

* Input data set


Configuration options are:

* "run_solver", "model_eq", "run_processes", "run_app" (see code and example run configuration file for details)


Outputs are:

* Run results


Main calculations are:

* This rule simply takes the input data set and performs a run unsing the fratoo framework and based on the provided run parameters for a scenario.

Main assumptions are:

* -


------------------------------------
Process results (process_results.py)
------------------------------------

This processes the results of model runs.

Main input data are:

* Run results


Configuration options are:

* "run_exp_res" (see code and example run configuration file for details)


Outputs are:

* Processed/extended run results


Main calculations are:

* This rule processes the run results. Core element is the calculation of costs of supply of energy carriers, which is calculated by adding up costs across the supply chain through a generic, iterative approach.

Main assumptions are:

* When calculating energy carrier supply costs, investment cost are generall allocated based on the operation of the technology, i.e., costs are allocated to years based on the use of the technology in that year. More details are provided in the code.



