"""
OSeMOSYS Pyomo code, based on the standard OSeMOSYS version (see copyright notice below)


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license
"""

# OSeMOSYS standard version's copyright notice
# ============================================================================
#
#    Copyright [2010-2015] [OSeMOSYS Forum steering committee see: www.osemosys.org]
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
# ============================================================================


from __future__ import division
from pyomo.environ import *
from pyomo.core import *
from pyomo.opt import SolverFactory

model = AbstractModel()


###############
#    Sets     #
###############

model.YEAR = Set()
model.TECHNOLOGY = Set()
model.TIMESLICE = Set()
model.FUEL = Set()
model.EMISSION = Set()
model.MODE_OF_OPERATION = Set()
model.REGION = Set()
model.SEASON = Set()
model.DAYTYPE = Set()
model.DAILYTIMEBRACKET = Set()
model.FLEXIBLEDEMANDTYPE = Set()
model.STORAGE = Set()

model.UDC = Set()

#####################
#    Parameters     #
#####################

########			Global 						#############

model.YearSplit = Param(model.TIMESLICE, model.YEAR)
model.DiscountRate = Param(model.REGION, default=0.05)
model.DiscountRateIdv = Param(model.REGION, model.TECHNOLOGY, default=0.05)
model.DiscountFactor = Param(model.REGION, model.YEAR,mutable=True)
model.DiscountFactorMid = Param(model.REGION, model.YEAR,mutable=True)
model.OperationalLife = Param(model.REGION, model.TECHNOLOGY, default=1)
model.CapitalRecoveryFactor= Param(model.REGION, model.TECHNOLOGY,mutable=True)
model.PvAnnuity= Param(model.REGION, model.TECHNOLOGY,mutable=True)




model.DiscountRateStorage = Param(model.REGION, model.STORAGE, default=0.05)
model.DiscountFactorStorage = Param(model.REGION,model.STORAGE, model.YEAR,mutable=True)

model.DaySplit = Param(model.DAILYTIMEBRACKET, model.YEAR, default=0.00137)
model.Conversionls = Param(model.TIMESLICE, model.SEASON, default=0)
model.Conversionld = Param(model.TIMESLICE, model.DAYTYPE, default=0)
model.Conversionlh = Param(model.TIMESLICE, model.DAILYTIMEBRACKET, default=0)
model.DaysInDayType = Param(model.SEASON, model.DAYTYPE, model.YEAR, default=7)
model.TradeRoute = Param(model.REGION, model.REGION, model.FUEL, model.YEAR, default=0)
model.DepreciationMethod = Param(model.REGION, default=1)


def dfac(model,r,y):
    model.DiscountFactor[r,y] = (1 + model.DiscountRate[r])**(y - min(yy for yy in model.YEAR) + 0.0)
def dfacmid(model,r,y):
    model.DiscountFactorMid[r,y] = (1 + model.DiscountRate[r])**(y - min(yy for yy in model.YEAR) + 0.5)
def crfac(model,r,t):
    model.CapitalRecoveryFactor[r,t] = (1- (1 + model.DiscountRateIdv[r,t])**(-1))/(1-(1+model.DiscountRateIdv[r,t])**(-(model.OperationalLife[r,t])))
    #model.CapitalRecoveryFactor[r,t]=1
def pva(model,r,t):
    model.PvAnnuity[r,t] = (1- (1 + model.DiscountRate[r])**(-(model.OperationalLife[r,t])))*(1+model.DiscountRate[r])/model.DiscountRate[r]
    #model.PvAnnuity[r,t]=1
def dfacs(model,r,s,y):
    model.DiscountFactorStorage[r,s,y] = (1 + model.DiscountRateStorage[r,s])**(y - min(yy for yy in model.YEAR) + 0.0)


model.builtgb0 = BuildAction(model.REGION,model.YEAR,rule=dfac)
model.builtgb1 = BuildAction(model.REGION,model.YEAR,rule=dfacmid)
model.builtgb2 = BuildAction(model.REGION,model.TECHNOLOGY,rule=crfac)
model.builtgb3 = BuildAction(model.REGION,model.TECHNOLOGY,rule=pva)
model.builtgb4 = BuildAction(model.REGION,model.STORAGE,model.YEAR,rule=dfacs)

########			Demands 					#############

model.SpecifiedAnnualDemand = Param(model.REGION, model.FUEL, model.YEAR, default=0)
model.SpecifiedDemandProfile = Param(
    model.REGION, model.FUEL, model.TIMESLICE, model.YEAR, default=0
)
model.AccumulatedAnnualDemand = Param(model.REGION, model.FUEL, model.YEAR, default=0)

#########			Performance					#############

model.CapacityToActivityUnit = Param(model.REGION, model.TECHNOLOGY, default=1)
model.TechWithCapacityNeededToMeetPeakTS = Param(model.REGION,model.TECHNOLOGY,default=0)
model.CapacityFactor = Param(
    model.REGION, model.TECHNOLOGY, model.TIMESLICE, model.YEAR, default=1
)
model.AvailabilityFactor = Param(model.REGION, model.TECHNOLOGY, model.YEAR, default=1)
model.ResidualCapacity = Param(model.REGION, model.TECHNOLOGY, model.YEAR, default=0)
model.InputActivityRatio = Param(
    model.REGION,
    model.TECHNOLOGY,
    model.FUEL,
    model.MODE_OF_OPERATION,
    model.YEAR,
    default=0,
)
model.OutputActivityRatio = Param(
    model.REGION,
    model.TECHNOLOGY,
    model.FUEL,
    model.MODE_OF_OPERATION,
    model.YEAR,
    default=0,
)

#########			Technology Costs			#############

model.CapitalCost = Param(model.REGION, model.TECHNOLOGY, model.YEAR, default=0)
model.VariableCost = Param(
    model.REGION, model.TECHNOLOGY, model.MODE_OF_OPERATION, model.YEAR, default=0.00001
)
model.FixedCost = Param(model.REGION, model.TECHNOLOGY, model.YEAR, default=0)

#########           		Storage                 		#############

model.TechnologyToStorage = Param(
    model.REGION, model.TECHNOLOGY, model.STORAGE, model.MODE_OF_OPERATION, default=0
)
model.TechnologyFromStorage = Param(
    model.REGION, model.TECHNOLOGY, model.STORAGE, model.MODE_OF_OPERATION, default=0
)
model.StorageLevelStart = Param(model.REGION, model.STORAGE, default=999)
model.StorageMaxChargeRate = Param(model.REGION, model.STORAGE, default=99)
model.StorageMaxDischargeRate = Param(model.REGION, model.STORAGE, default=99)
model.MinStorageCharge = Param(model.REGION, model.STORAGE, model.YEAR, default=0)
model.OperationalLifeStorage = Param(model.REGION, model.STORAGE, default=99)
model.CapitalCostStorage = Param(model.REGION, model.STORAGE, model.YEAR, default=0)
model.ResidualStorageCapacity = Param(
    model.REGION, model.STORAGE, model.YEAR, default=0
)

#########			Capacity Constraints		#############

model.CapacityOfOneTechnologyUnit = Param(
    model.REGION, model.TECHNOLOGY, model.YEAR, default=0
)
model.TotalAnnualMaxCapacity = Param(
    model.REGION, model.TECHNOLOGY, model.YEAR, default=-1
)
model.TotalAnnualMinCapacity = Param(
    model.REGION, model.TECHNOLOGY, model.YEAR, default=0
)

#########			Investment Constraints		#############

model.TotalAnnualMaxCapacityInvestment = Param(
    model.REGION, model.TECHNOLOGY, model.YEAR, default=-1
)
model.TotalAnnualMinCapacityInvestment = Param(
    model.REGION, model.TECHNOLOGY, model.YEAR, default=0
)

#########			Activity Constraints		#############

model.TotalTechnologyAnnualActivityUpperLimit = Param(
    model.REGION, model.TECHNOLOGY, model.YEAR, default=-1
)
model.TotalTechnologyAnnualActivityLowerLimit = Param(
    model.REGION, model.TECHNOLOGY, model.YEAR, default=0
)
model.TotalTechnologyModelPeriodActivityUpperLimit = Param(
    model.REGION, model.TECHNOLOGY, default=-1
)
model.TotalTechnologyModelPeriodActivityLowerLimit = Param(
    model.REGION, model.TECHNOLOGY, default=0
)

#########			Reserve Margin				#############

model.ReserveMarginTagTechnology = Param(
    model.REGION, model.TECHNOLOGY, model.YEAR, default=0
)
model.ReserveMarginTagFuel = Param(model.REGION, model.FUEL, model.YEAR, default=0)
model.ReserveMargin = Param(model.REGION, model.YEAR, default=-1)

#########			RE Generation Target		#############

model.RETagTechnology = Param(model.REGION, model.TECHNOLOGY, model.YEAR, default=0)
model.RETagFuel = Param(model.REGION, model.FUEL, model.YEAR, default=0)
model.REMinProductionTarget = Param(model.REGION, model.YEAR, default=0)

#########			Emissions & Penalties		#############

model.EmissionActivityRatio = Param(
    model.REGION,
    model.TECHNOLOGY,
    model.EMISSION,
    model.MODE_OF_OPERATION,
    model.YEAR,
    default=0,
)
model.EmissionsPenalty = Param(model.REGION, model.EMISSION, model.YEAR, default=0)
model.AnnualExogenousEmission = Param(
    model.REGION, model.EMISSION, model.YEAR, default=0
)
model.AnnualEmissionLimit = Param(
    model.REGION, model.EMISSION, model.YEAR, default=-1
)
model.ModelPeriodExogenousEmission = Param(model.REGION, model.EMISSION, default=0)
model.ModelPeriodEmissionLimit = Param(model.REGION, model.EMISSION, default=-1)


# Added parameters

model.UDCMultiplierTotalCapacity = Param(model.REGION, model.TECHNOLOGY, model.UDC, model.YEAR, default=0)
model.UDCMultiplierNewCapacity = Param(model.REGION, model.TECHNOLOGY, model.UDC, model.YEAR, default=0)
model.UDCMultiplierActivity = Param(model.REGION, model.TECHNOLOGY,model.MODE_OF_OPERATION, model.UDC, model.YEAR, default=0)
model.UDCConstant = Param(model.REGION, model.UDC, model.YEAR, default=-1)
model.UDCTag = Param(model.REGION, model.UDC, default=-1)



### Reduced sets
model.MODETECHNOLOGYFUELOUT = Set(within=model.FUEL*model.MODE_OF_OPERATION*model.TECHNOLOGY)
model.MODETECHNOLOGYFUELIN = Set(within=model.FUEL*model.MODE_OF_OPERATION*model.TECHNOLOGY)
model.MODExTECHNOLOGYperStorageto = Set(within=model.STORAGE*model.MODE_OF_OPERATION*model.TECHNOLOGY)
model.MODExTECHNOLOGYperStoragefrom = Set(within=model.STORAGE*model.MODE_OF_OPERATION*model.TECHNOLOGY)
model.TIMESLICEofSEASON = Set(within=model.SEASON*model.TIMESLICE)
model.TIMESLICEofDAYTYPE = Set(within=model.DAYTYPE*model.TIMESLICE)
model.TIMESLICEofDAILYTIMEBRACKET = Set(within=model.DAILYTIMEBRACKET*model.TIMESLICE)
model.TIMESLICEofSDB = Set(within=model.SEASON*model.DAYTYPE
                           *model.DAILYTIMEBRACKET*model.TIMESLICE)
model.MODETECHNOLOGYEMISSION = Set(within=model.EMISSION*model.MODE_OF_OPERATION*model.TECHNOLOGY)


model.MODETECHNOLOGYFUELOUTind = Set(model.FUEL,within=model.MODE_OF_OPERATION*model.TECHNOLOGY)
model.MODETECHNOLOGYFUELINind = Set(model.FUEL,within=model.MODE_OF_OPERATION*model.TECHNOLOGY)
# model.MODExTECHNOLOGYperStorageto = Set(within=model.STORAGE*model.MODE_OF_OPERATION*model.TECHNOLOGY)
# model.MODExTECHNOLOGYperStoragefrom = Set(within=model.STORAGE*model.MODE_OF_OPERATION*model.TECHNOLOGY)
model.MODETECHNOLOGYEMISSIONind = Set(model.EMISSION,within=model.MODE_OF_OPERATION*model.TECHNOLOGY)



model.MODETECHNOLOGY = Set(within=model.TECHNOLOGY*model.MODE_OF_OPERATION)
model.MODETECHNOLOGYind = Set(model.TECHNOLOGY,within=model.MODE_OF_OPERATION)
model.TECHNOLOGYCAPYEARS = Set(model.REGION,model.TECHNOLOGY,
                               model.YEAR,within=model.YEAR)

def mxtfout(model):
    if len(model.MODETECHNOLOGYFUELOUT.ordered_data())>0:
        for f,t,m in model.MODETECHNOLOGYFUELOUT:
                model.MODETECHNOLOGYFUELOUTind[f].add((t,m))
        return
    for f in model.FUEL:
        for m in model.MODE_OF_OPERATION:
            for t in model.TECHNOLOGY:
                for r in model.REGION:
                    for y in model.YEAR:
                        if model.OutputActivityRatio[r, t, f, m, y] != 0:
                            model.MODETECHNOLOGYFUELOUT.add((f,m,t))
                            break
                    else:
                        continue
                    break

          
def mxtfin(model):
    if len(model.MODETECHNOLOGYFUELIN.ordered_data())>0:
        for f,t,m in model.MODETECHNOLOGYFUELIN:
                model.MODETECHNOLOGYFUELINind[f].add((t,m))
        return
    for f in model.FUEL:
        for m in model.MODE_OF_OPERATION:
            for t in model.TECHNOLOGY:
                for r in model.REGION:
                    for y in model.YEAR:
                        if model.InputActivityRatio[r, t, f, m, y] != 0:
                            model.MODETECHNOLOGYFUELIN.add((f,m,t))
                            break
                    else:
                        continue
                    break
    pass
                
def mxtsout(model):
    for s in model.STORAGE:
        for m in model.MODE_OF_OPERATION:
            for t in model.TECHNOLOGY:
                for r in model.REGION:
                    if TechnologyToStorage[r,t,s,m] > 0:
                        model.MODExTECHNOLOGYperStorageto.add((s,m,t))
                        break
                    
def mxtsin(model):
    for s in model.STORAGE:
        for m in model.MODE_OF_OPERATION:
            for t in model.TECHNOLOGY:
                for r in model.REGION:
                    if TechnologyFromStorage[r,t,s,m] > 0:
                        model.MODExTECHNOLOGYperStoragefrom.add((s,m,t))
                        break
def tos(model):
    for ls in model.SEASON:
        for l in model.TIMESLICE:
            if model.Conversionls[l,ls] == 1:
                model.TIMESLICEofSEASON.add((ls,l))
def tod(model):
    for ld in model.DAYTYPE:
        for l in model.TIMESLICE:
            if model.Conversionld[l,ld] == 1:
                model.TIMESLICEofDAYTYPE.add((ld,l))
def todtb(model):
    for lh in model.DAILYTIMEBRACKET:
        for l in model.TIMESLICE:
            if model.Conversionlh[l,lh] == 1:
                model.TIMESLICEofDAILYTIMEBRACKET.add((lh,l))   
             
def mxte(model):
    if len(model.MODETECHNOLOGYEMISSION.ordered_data())>0:
        for e,t,m in model.MODETECHNOLOGYEMISSION:
                model.MODETECHNOLOGYEMISSIONind[e].add((t,m))
        return
    for e in model.EMISSION:
        for m in model.MODE_OF_OPERATION:
            for t in model.TECHNOLOGY:
                for r in model.REGION:
                    for y in model.YEAR:
                        if model.EmissionActivityRatio[r, t, e, m, y] != 0:
                            model.MODETECHNOLOGYEMISSION.add((e,m,t))
                            break
                    else:
                        continue
def tm(model):
    if len(model.MODETECHNOLOGY.ordered_data())>0:
        return
    for t in model.TECHNOLOGY:
        for m in model.MODE_OF_OPERATION:
            
            for f in model.FUEL:
                if (f,m,t) in (model.MODETECHNOLOGYFUELOUT
                              |model.MODETECHNOLOGYFUELIN):
                    model.MODETECHNOLOGY.add((t,m))
                    break
            else:
                continue
            for s in model.STORAGE:
                if (s,m,t) in (model.MODExTECHNOLOGYperSTORAGEto
                              |model.MODExTECHNOLOGYperSTORAGEfrom):
                    model.MODETECHNOLOGY.add((t,m))
                    break
            else:
                continue
            for e in model.EMISSION:
                if (e,m,t) in model.MODETECHNOLOGYEMISSION:
                    model.MODETECHNOLOGY.add((t,m))
                    break
            else:
                continue
 
def tmind(model):
    for t,m in model.MODETECHNOLOGY:
            model.MODETECHNOLOGYind[t].add(m)
         
 
def tcy(model,r,t,y):
    for yy in model.YEAR:
        if (((y-yy)< model.OperationalLife[r,t])and ((y-yy)>=0)):
            model.TECHNOLOGYCAPYEARS[r,t,y].add(yy)



# def tm(model):
#     if sum(len(model.MODETECHNOLOGYsparse[t].ordered_data()) for t in model.TECHNOLOGY)>0:
#         return      
#     for t in model.TECHNOLOGY:
#         for m in model.MODE_OF_OPERATION:
            
#             for f in model.FUEL:
#                 if (f,m,t) in (model.MODETECHNOLOGYFUELOUT
#                               |model.MODETECHNOLOGYFUELIN):
#                     model.MODETECHNOLOGYsparse[t].add((m))
#                     break
#             else:
#                 continue
#             for s in model.STORAGE:
#                 if (s,m,t) in (model.MODExTECHNOLOGYperSTORAGEto
#                               |model.MODExTECHNOLOGYperSTORAGEfrom):
#                     model.MODETECHNOLOGYsparse[t].add((m))
#                     break
#             else:
#                 continue
#             for e in model.EMISSION:
#                 if (e,m,t) in model.MODETECHNOLOGYEMISSION:
#                     model.MODETECHNOLOGYsparse[t].add((m))
#                     break
#             else:
#                 continue
                
model.built1 = BuildAction(rule=mxtfout)
model.built2 = BuildAction(rule=mxtfin)
model.built3 = BuildAction(rule=mxtsout)
model.built4 = BuildAction(rule=mxtsin)
model.built5 = BuildAction(rule=tos)
model.built6 = BuildAction(rule=tod)
model.built7 = BuildAction(rule=todtb)
model.built8 = BuildAction(rule=mxte)
model.built9 = BuildAction(rule=tm)
model.built9a = BuildAction(rule=tmind)
model.built10 = BuildAction(model.REGION,model.TECHNOLOGY, model.YEAR,rule=tcy)

# additions for UK model

model.CapacityFractionTagTechnology = Param(model.REGION,model.TECHNOLOGY,
                                            default=0)
model.CapacityFraction = Param(model.REGION,model.TECHNOLOGY,model.YEAR,
                               default=-1)
# simple storage
model.StorageTagFuel = Param(model.REGION,model.FUEL,
                             default=0)


######################
#   Model Variables  #
######################

########			Demands 					#############

# model.RateOfDemand = Var(
#     model.REGION,
#     model.TIMESLICE,
#     model.FUEL,
#     model.YEAR,
#     domain=NonNegativeReals,
#     initialize=0.0,
# )
# model.Demand = Var(
#     model.REGION,
#     model.TIMESLICE,
#     model.FUEL,
#     model.YEAR,
#     domain=NonNegativeReals,
#     initialize=0.0,
# )

########     		Storage                 		#############

model.NewStorageCapacity = Var(
    model.REGION, model.STORAGE, model.YEAR, domain=NonNegativeReals, initialize=0.0
)

model.SalvageValueStorage = Var(
    model.REGION, model.STORAGE, model.YEAR, domain=NonNegativeReals, initialize=0.0
)

model.StorageLevelYearStart = Var(
    model.REGION, model.STORAGE, model.YEAR, domain=NonNegativeReals, initialize=0.0
)
model.StorageLevelYearFinish = Var(
    model.REGION, model.STORAGE, model.YEAR, domain=NonNegativeReals, initialize=0.0
)
model.StorageLevelSeasonStart = Var(
    model.REGION,
    model.STORAGE,
    model.SEASON,
    model.YEAR,
    domain=NonNegativeReals,
    initialize=0.0,
)
model.StorageLevelDayTypeStart = Var(
    model.REGION,
    model.STORAGE,
    model.SEASON,
    model.DAYTYPE,
    model.YEAR,
    domain=NonNegativeReals,
    initialize=0.0,
)
model.StorageLevelDayTypeFinish = Var(
    model.REGION,
    model.STORAGE,
    model.SEASON,
    model.DAYTYPE,
    model.YEAR,
    domain=NonNegativeReals,
    initialize=0.0,
)

# model.RateOfStorageCharge = Var(
#     model.REGION,
#     model.STORAGE,
#     model.SEASON,
#     model.DAYTYPE,
#     model.DAILYTIMEBRACKET,
#     model.YEAR,
#     initialize=0.0,
# )
# model.RateOfStorageDischarge = Var(
#     model.REGION,
#     model.STORAGE,
#     model.SEASON,
#     model.DAYTYPE,
#     model.DAILYTIMEBRACKET,
#     model.YEAR,
#     initialize=0.0,
# )
# model.NetChargeWithinYear = Var(
#     model.REGION,
#     model.STORAGE,
#     model.SEASON,
#     model.DAYTYPE,
#     model.DAILYTIMEBRACKET,
#     model.YEAR,
#     initialize=0.0,
# )
# model.NetChargeWithinDay = Var(
#     model.REGION,
#     model.STORAGE,
#     model.SEASON,
#     model.DAYTYPE,
#     model.DAILYTIMEBRACKET,
#     model.YEAR,
#     initialize=0.0,
# )

# model.StorageLowerLimit = Var(
#     model.REGION, model.STORAGE, model.YEAR, domain=NonNegativeReals, initialize=0.0
# )
# model.StorageUpperLimit = Var(
#     model.REGION, model.STORAGE, model.YEAR, domain=NonNegativeReals, initialize=0.0
# )
# model.AccumulatedNewStorageCapacity = Var(
#     model.REGION, model.STORAGE, model.YEAR, domain=NonNegativeReals, initialize=0.0
# )

# model.CapitalInvestmentStorage = Var(
#     model.REGION, model.STORAGE, model.YEAR, domain=NonNegativeReals, initialize=0.0
# )
# model.DiscountedCapitalInvestmentStorage = Var(
#     model.REGION, model.STORAGE, model.YEAR, domain=NonNegativeReals, initialize=0.0
# )

# model.DiscountedSalvageValueStorage = Var(
#     model.REGION, model.STORAGE, model.YEAR, domain=NonNegativeReals, initialize=0.0
# )
# model.TotalDiscountedStorageCost = Var(
#     model.REGION, model.STORAGE, model.YEAR, domain=NonNegativeReals, initialize=0.0
# )

#########		    Capacity Variables 			#############
# FIXME: this is actually an domain=NonNegativeInteger, but this causes some solvers not to provide shadow prices
model.NumberOfNewTechnologyUnits = Var(
    model.REGION, model.TECHNOLOGY, model.YEAR, domain=NonNegativeReals, initialize=0
)
model.NewCapacity = Var(
    model.REGION, model.TECHNOLOGY, model.YEAR, domain=NonNegativeReals, initialize=0.0
)
# model.AccumulatedNewCapacity = Var(
#     model.REGION, model.TECHNOLOGY, model.YEAR, domain=NonNegativeReals, initialize=0.0
# )
# model.TotalCapacityAnnual = Var(
#     model.REGION, model.TECHNOLOGY, model.YEAR, domain=NonNegativeReals, initialize=0.0
# )

#########		    Activity Variables 			#############


model.RateOfActivity = Var(
    model.REGION,
    model.TIMESLICE,
    model.MODETECHNOLOGY,
    model.YEAR,
    domain=NonNegativeReals,
    initialize=0.0,
)

model.Trade = Var(
    model.REGION, model.REGION, model.TIMESLICE, model.FUEL, model.YEAR, initialize=0.0
)

model.UseAnnual = Var(
    model.REGION, model.FUEL, model.YEAR, domain=NonNegativeReals, initialize=0.0
)


# model.RateOfTotalActivity = Var(
#     model.REGION,
#     model.TECHNOLOGY,
#     model.TIMESLICE,
#     model.YEAR,
#     domain=NonNegativeReals,
#     initialize=0.0,
# )
# model.TotalTechnologyAnnualActivity = Var(
#     model.REGION, model.TECHNOLOGY, model.YEAR, domain=NonNegativeReals, initialize=0.0
# )
# model.TotalAnnualTechnologyActivityByMode = Var(
#     model.REGION,
#     model.TECHNOLOGY,
#     model.MODE_OF_OPERATION,
#     model.YEAR,
#     domain=NonNegativeReals,
#     initialize=0.0,
# )
# model.RateOfProductionByTechnologyByMode = Var(
#     model.REGION,
#     model.TIMESLICE,
#     model.TECHNOLOGY,
#     model.MODE_OF_OPERATION,
#     model.FUEL,
#     model.YEAR,
#     domain=NonNegativeReals,
#     initialize=0.0,
# )
# model.RateOfProductionByTechnology = Var(
#     model.REGION,
#     model.TIMESLICE,
#     model.TECHNOLOGY,
#     model.FUEL,
#     model.YEAR,
#     domain=NonNegativeReals,
#     initialize=0.0,
# )
# model.ProductionByTechnology = Var(
#     model.REGION,
#     model.TIMESLICE,
#     model.TECHNOLOGY,
#     model.FUEL,
#     model.YEAR,
#     domain=NonNegativeReals,
#     initialize=0.0,
# )
# model.ProductionByTechnologyAnnual = Var(
#     model.REGION,
#     model.TECHNOLOGY,
#     model.FUEL,
#     model.YEAR,
#     domain=NonNegativeReals,
#     initialize=0.0,
# )
# model.RateOfProduction = Var(
#     model.REGION,
#     model.TIMESLICE,
#     model.FUEL,
#     model.YEAR,
#     domain=NonNegativeReals,
#     initialize=0.0,
# )
# model.Production = Var(
#     model.REGION,
#     model.TIMESLICE,
#     model.FUEL,
#     model.YEAR,
#     domain=NonNegativeReals,
#     initialize=0.0,
# )
# model.RateOfUseByTechnologyByMode = Var(
#     model.REGION,
#     model.TIMESLICE,
#     model.TECHNOLOGY,
#     model.MODE_OF_OPERATION,
#     model.FUEL,
#     model.YEAR,
#     domain=NonNegativeReals,
#     initialize=0.0,
# )
# model.RateOfUseByTechnology = Var(
#     model.REGION,
#     model.TIMESLICE,
#     model.TECHNOLOGY,
#     model.FUEL,
#     model.YEAR,
#     domain=NonNegativeReals,
#     initialize=0.0,
# )
# model.UseByTechnologyAnnual = Var(
#     model.REGION,
#     model.TECHNOLOGY,
#     model.FUEL,
#     model.YEAR,
#     domain=NonNegativeReals,
#     initialize=0.0,
# )
# model.RateOfUse = Var(
#     model.REGION,
#     model.TIMESLICE,
#     model.FUEL,
#     model.YEAR,
#     domain=NonNegativeReals,
#     initialize=0.0,
# )
# model.UseByTechnology = Var(
#     model.REGION,
#     model.TIMESLICE,
#     model.TECHNOLOGY,
#     model.FUEL,
#     model.YEAR,
#     domain=NonNegativeReals,
#     initialize=0.0,
# )
# model.Use = Var(
#     model.REGION,
#     model.TIMESLICE,
#     model.FUEL,
#     model.YEAR,
#     domain=NonNegativeReals,
#     initialize=0.0,
# )

# model.TradeAnnual = Var(
#     model.REGION, model.REGION, model.FUEL, model.YEAR, initialize=0.0
# )

# model.ProductionAnnual = Var(
#     model.REGION, model.FUEL, model.YEAR, domain=NonNegativeReals, initialize=0.0
# )


#########		    Costing Variables 			#############
model.VariableOperatingCost = Var(
    model.REGION,
    model.TECHNOLOGY,
    model.TIMESLICE,
    model.YEAR,
    domain=NonNegativeReals,
    initialize=0.0,
)

model.SalvageValue = Var(
    model.REGION, model.TECHNOLOGY, model.YEAR, domain=NonNegativeReals, initialize=0.0
)
model.DiscountedCapitalInvestment = Var(
    model.REGION, model.TECHNOLOGY, model.YEAR, domain=NonNegativeReals, initialize=0.0
)

model.DiscountedSalvageValue = Var(
    model.REGION, model.TECHNOLOGY, model.YEAR, domain=NonNegativeReals, initialize=0.0
)
model.OperatingCost = Var(
    model.REGION, model.TECHNOLOGY, model.YEAR, domain=NonNegativeReals, initialize=0.0
)
# model.CapitalInvestment = Var(
#     model.REGION, model.TECHNOLOGY, model.YEAR, domain=NonNegativeReals, initialize=0.0
# )

# model.DiscountedOperatingCost = Var(
#     model.REGION, model.TECHNOLOGY, model.YEAR, domain=NonNegativeReals, initialize=0.0
# )

# model.AnnualVariableOperatingCost = Var(
#     model.REGION, model.TECHNOLOGY, model.YEAR, domain=NonNegativeReals, initialize=0.0
# )
# model.AnnualFixedOperatingCost = Var(
#     model.REGION, model.TECHNOLOGY, model.YEAR, domain=NonNegativeReals, initialize=0.0
# )


# model.TotalDiscountedCostByTechnology = Var(
#     model.REGION, model.TECHNOLOGY, model.YEAR, domain=NonNegativeReals, initialize=0.0
# )
# model.TotalDiscountedCost = Var(
#     model.REGION, model.YEAR, domain=NonNegativeReals, initialize=0.0
# )

# model.ModelPeriodCostByRegion = Var(
#     model.REGION, domain=NonNegativeReals, initialize=0.0
# )

#########			Reserve Margin				#############

# model.TotalCapacityInReserveMargin = Var(
#     model.REGION, model.YEAR, domain=NonNegativeReals, initialize=0.0
# )
# model.DemandNeedingReserveMargin = Var(
#     model.REGION, model.TIMESLICE, model.YEAR, domain=NonNegativeReals, initialize=0.0
# )

#########			RE Gen Target				#############

# model.TotalREProductionAnnual = Var(model.REGION, model.YEAR, initialize=0.0)
# model.RETotalDemandOfTargetFuelAnnual = Var(model.REGION, model.YEAR, initialize=0.0)

# model.TotalTechnologyModelPeriodActivity = Var(
#     model.REGION, model.TECHNOLOGY, initialize=0.0
# )

#########			Emissions					#############
model.DiscountedTechnologyEmissionsPenalty = Var(
    model.REGION, model.TECHNOLOGY, model.YEAR, domain=NonNegativeReals, initialize=0.0
)

model.ModelPeriodEmissions = Var(
    model.REGION, model.EMISSION, domain=NonNegativeReals, initialize=0.0
)

# model.AnnualTechnologyEmissionByMode = Var(
#     model.REGION,
#     model.TECHNOLOGY,
#     model.EMISSION,
#     model.MODE_OF_OPERATION,
#     model.YEAR,
#     domain=NonNegativeReals,
#     initialize=0.0,
# )
# model.AnnualTechnologyEmission = Var(
#     model.REGION,
#     model.TECHNOLOGY,
#     model.EMISSION,
#     model.YEAR,
#     domain=NonNegativeReals,
#     initialize=0.0,
# )
# model.AnnualTechnologyEmissionPenaltyByEmission = Var(
#     model.REGION,
#     model.TECHNOLOGY,
#     model.EMISSION,
#     model.YEAR,
#     domain=NonNegativeReals,
#     initialize=0.0,
# )
# model.AnnualEmissions = Var(
#     model.REGION, model.EMISSION, model.YEAR, domain=NonNegativeReals, initialize=0.0
# )
# model.AnnualTechnologyEmissionsPenalty = Var(
#     model.REGION, model.TECHNOLOGY, model.YEAR, domain=NonNegativeReals, initialize=0.0
# )


######################
# Objective Function #
######################


# def ObjectiveFunction_rule(model):
#     return (
#         (sum((
#             (sum(model.NewCapacity[r,t,yy] 
#                for yy in model.YEAR if (((y-yy)< model.OperationalLife[r,t])and
#                                         ((y-yy)>=0)))+ model.ResidualCapacity[r,t,y])
#          *model.FixedCost[r,t,y]
#          +
#          sum(model.RateOfActivity[r,l,(t,m),y] * model.YearSplit[l, y]
#              *model.VariableCost[r,t,m,y]
#              for l in model.TIMESLICE for tt,m in model.MODETECHNOLOGY if (tt==t))
#         )/ model.DiscountFactorMid[r,y]
#         +
#         (model.CapitalCost[r,t,y] * model.NewCapacity[r,t,y] 
#         * model.CapitalRecoveryFactor[r,t] * model.PvAnnuity[r,t] 
#         / model.DiscountFactor[r,y] )
#         +
#         model.DiscountedTechnologyEmissionsPenalty[r,t,y] 
#         - model.DiscountedSalvageValue[r,t,y]
#         for t in model.TECHNOLOGY
#         for r in model.REGION
#         for y in model.YEAR)
#         +
#         sum(model.CapitalCostStorage[r,s,y] * model.NewStorageCapacity[r,s,y] 
#             / model.DiscountFactor[r,y]
#             -model.DiscountedSalvageValueStorage[r,s,y]
#             for s in model.STORAGE
#             for r in model.REGION
#             for y in model.YEAR)
#         )
#         )

def ObjectiveFunction_rule(model):
    return (

          sum(model.RateOfActivity[r,l,(t,m),y] * model.YearSplit[l, y]
              *model.VariableCost[r,t,m,y]/ model.DiscountFactorMid[r,y]
            for l in model.TIMESLICE 
            for t,m in model.MODETECHNOLOGY
            for r in model.REGION
            for y in model.YEAR)       
        +
        sum(
            (sum(model.NewCapacity[r,t,yy] 
                for yy in model.TECHNOLOGYCAPYEARS[r,t,y])+ model.ResidualCapacity[r,t,y])
          *model.FixedCost[r,t,y]/ model.DiscountFactorMid[r,y]
          +
        (model.CapitalCost[r,t,y] * model.NewCapacity[r,t,y] 
        * model.CapitalRecoveryFactor[r,t] * model.PvAnnuity[r,t] 
        / model.DiscountFactor[r,y] )
        +
        model.DiscountedTechnologyEmissionsPenalty[r,t,y] 
        - model.DiscountedSalvageValue[r,t,y]
        for t in model.TECHNOLOGY
        for r in model.REGION
        for y in model.YEAR)
        +
        sum(model.CapitalCostStorage[r,s,y] * model.NewStorageCapacity[r,s,y] 
            / model.DiscountFactor[r,y]
            -model.DiscountedSalvageValueStorage[r,s,y]
            for s in model.STORAGE
            for r in model.REGION
            for y in model.YEAR)
        )
        


model.OBJ = Objective(rule=ObjectiveFunction_rule, sense=minimize)


#####################
# Constraints       #
#####################


# def SpecifiedDemand_rule(model, r, f, l, y):
#     return (
#         model.SpecifiedAnnualDemand[r, f, y]
#         * model.SpecifiedDemandProfile[r, f, l, y]
#         / model.YearSplit[l, y]
#         == model.RateOfDemand[r, l, f, y]
#     )


# model.SpecifiedDemand = Constraint(
#     model.REGION, model.FUEL, model.TIMESLICE, model.YEAR, rule=SpecifiedDemand_rule
# )


#########       	Capacity Adequacy A	     	#############


# def TotalNewCapacity_1_rule(model, r, t, y):
#     return model.AccumulatedNewCapacity[r, t, y] == sum(
#         model.NewCapacity[r, t, yy]
#         for yy in model.YEAR
#         if ((y - yy < model.OperationalLife[r, t]) and (y - yy >= 0))
#     )


# model.TotalNewCapacity_1 = Constraint(
#     model.REGION, model.TECHNOLOGY, model.YEAR, rule=TotalNewCapacity_1_rule
# )

# def TotalAnnualCapacity_rule(model, r, t, y):
#     return (
#         model.AccumulatedNewCapacity[r, t, y] + model.ResidualCapacity[r, t, y]
#         == model.TotalCapacityAnnual[r, t, y]
#     )


# model.TotalAnnualCapacity_constraint = Constraint(
#     model.REGION, model.TECHNOLOGY, model.YEAR, rule=TotalAnnualCapacity_rule
# )


# def TotalActivityOfEachTechnology_rule(model, r, t, l, y):
#     return (
#         sum(model.RateOfActivity[r, l, t, m, y] for m in model.MODE_OF_OPERATION)
#         == model.RateOfTotalActivity[r, t, l, y]
#     )


# model.TotalActivityOfEachTechnology = Constraint(
#     model.REGION,
#     model.TECHNOLOGY,
#     model.TIMESLICE,
#     model.YEAR,
#     rule=TotalActivityOfEachTechnology_rule,
# )



def ConstraintCapacity_rule(model, r, l, t, y):
    return (
        sum(model.RateOfActivity[r,l,(t,m),y] 
            for m in model.MODETECHNOLOGYind[t])
        <= ((sum(model.NewCapacity[r,t,yy]
            for yy in model.TECHNOLOGYCAPYEARS[r,t,y])+ model.ResidualCapacity[r,t,y])
            *model.CapacityFactor[r,t,l,y]*model.CapacityToActivityUnit[r,t])    
    )
        


model.ConstraintCapacity = Constraint(
    model.REGION,
    model.TIMESLICE,
    model.TECHNOLOGY,
    model.YEAR,
    rule=ConstraintCapacity_rule,
)

def TotalNewCapacity_2_rule(model,r,t,y):
    if model.CapacityOfOneTechnologyUnit[r,t,y] != 0:
        return model.CapacityOfOneTechnologyUnit[r,t,y]*model.NumberOfNewTechnologyUnits[r,t,y] == model.NewCapacity[r,t,y]
    else:
        return Constraint.Skip
    
model.TotalNewCapacity_2 = Constraint(model.REGION, model.TECHNOLOGY, model.YEAR, rule=TotalNewCapacity_2_rule)



#########       	Capacity Adequacy B		 	#############


def PlannedMaintenance_rule(model, r, t, y):
    return (
        sum(
            model.RateOfActivity[r,l,(t,m),y] * model.YearSplit[l, y]
            for l in model.TIMESLICE for m in model.MODETECHNOLOGYind[t]
        )
        <= (sum(
            (sum(model.NewCapacity[r,t,yy] 
                for yy in model.TECHNOLOGYCAPYEARS[r,t,y])+ model.ResidualCapacity[r,t,y])
            * model.CapacityFactor[r, t, l, y]
            * model.YearSplit[l, y]
            for l in model.TIMESLICE
        )
        * model.AvailabilityFactor[r, t, y]
        * model.CapacityToActivityUnit[r, t]
        )
    )   

model.PlannedMaintenance = Constraint(
    model.REGION, model.TECHNOLOGY, model.YEAR, rule=PlannedMaintenance_rule
)



#########	        Energy Balance A    	 	#############


def EnergyBalanceEachTS4_rule(model, r, rr, l, f, y):
    if model.TradeRoute[r,rr,f,y] !=0:
        return model.Trade[r, rr, l, f, y] + model.Trade[rr, r, l, f, y] == 0
    else:
        return Constraint.Skip

model.EnergyBalanceEachTS4 = Constraint(
    model.REGION,
    model.REGION,
    model.TIMESLICE,
    model.FUEL,
    model.YEAR,
    rule=EnergyBalanceEachTS4_rule,
)

# FIXME: if/else structure below is not existent in standard OSeMOSYS version
def EnergyBalanceEachTS5_rule(model, r, l, f, y):
    if model.StorageTagFuel[r,f]!=1:
        return (sum(model.RateOfActivity[r,l,(t,m),y]
                    *model.OutputActivityRatio[r,t,f,m,y]
                    *model.YearSplit[l,y]
                    for m,t in model.MODETECHNOLOGYFUELOUTind[f])
                    >= (
                    model.SpecifiedAnnualDemand[r,f,y]
                    * model.SpecifiedDemandProfile[r,f,l,y]
                    + sum(model.RateOfActivity[r,l,(t,m),y]
                                *model.InputActivityRatio[r,t,f,m,y]
                                *model.YearSplit[l,y]
                                for m,t in model.MODETECHNOLOGYFUELINind[f])
                    + sum(model.Trade[r, rr, l, f, y] * model.TradeRoute[r, rr, f, y]
                          for rr in model.REGION))
        )
    else:
        return Constraint.Skip


model.EnergyBalanceEachTS5 = Constraint(
    model.REGION,
    model.TIMESLICE,
    model.FUEL,
    model.YEAR,
    rule=EnergyBalanceEachTS5_rule,
)


# def RateOfFuelProduction1_rule(model, r, l, f, t, m, y):
#     if model.OutputActivityRatio[r, t, f, m, y] != 0:
#         return (
#             model.RateOfProductionByTechnologyByMode[r, l, t, m, f, y]
#             == model.RateOfActivity[r, l, t, m, y]
#             * model.OutputActivityRatio[r, t, f, m, y]
#         )
#     else:
#         return model.RateOfProductionByTechnologyByMode[r, l, t, m, f, y] == 0


# model.RateOfFuelProduction1 = Constraint(
#     model.REGION,
#     model.TIMESLICE,
#     model.FUEL,
#     model.TECHNOLOGY,
#     model.MODE_OF_OPERATION,
#     model.YEAR,
#     rule=RateOfFuelProduction1_rule,
# )

# def RateOfFuelProduction1_rule(model,r,l,f,t,m,y):
# return model.RateOfProductionByTechnologyByMode[r,l,t,m,f,y] == model.RateOfActivity[r,l,t,m,y]*model.OutputActivityRatio[r,t,f,m,y]
# model.RateOfFuelProduction1 = Constraint(model.REGION, model.TIMESLICE, model.FUEL, model.TECHNOLOGY, model.MODE_OF_OPERATION, model.YEAR, rule=RateOfFuelProduction1_rule)

# def RateOfFuelProduction2_rule(model,r,l,f,t,y):
# return  model.RateOfProductionByTechnology[r,l,t,f,y] == sum(model.RateOfProductionByTechnologyByMode[r,l,t,m,f,y] for m in model.MODE_OF_OPERATION if model.OutputActivityRatio[r,t,f,m,y] != 0)
# model.RateOfFuelProduction2 = Constraint(model.REGION, model.TIMESLICE, model.FUEL, model.TECHNOLOGY, model.YEAR, rule=RateOfFuelProduction2_rule)


# def RateOfFuelProduction2_rule(model, r, l, f, t, y):
#     return model.RateOfProductionByTechnology[r, l, t, f, y] == sum(
#         model.RateOfProductionByTechnologyByMode[r, l, t, m, f, y]
#         for m in model.MODE_OF_OPERATION
#     )


# model.RateOfFuelProduction2 = Constraint(
#     model.REGION,
#     model.TIMESLICE,
#     model.FUEL,
#     model.TECHNOLOGY,
#     model.YEAR,
#     rule=RateOfFuelProduction2_rule,
# )


# def RateOfFuelProduction3_rule(model, r, l, f, y):
#     return model.RateOfProduction[r, l, f, y] == sum(
#         model.RateOfProductionByTechnology[r, l, t, f, y] for t in model.TECHNOLOGY
#     )


# model.RateOfFuelProduction3 = Constraint(
#     model.REGION,
#     model.TIMESLICE,
#     model.FUEL,
#     model.YEAR,
#     rule=RateOfFuelProduction3_rule,
# )

# def RateOfFuelUse1_rule(model,r,l,f,t,m,y):
# if model.InputActivityRatio[r,t,f,m,y] != 0:
# return model.RateOfActivity[r,l,t,m,y]*model.InputActivityRatio[r,t,f,m,y] == model.RateOfUseByTechnologyByMode[r,l,t,m,f,y]
# else:
# return model.RateOfUseByTechnologyByMode[r,l,t,m,f,y] == 0
# model.RateOfFuelUse1 = Constraint(model.REGION, model.TIMESLICE, model.FUEL, model.TECHNOLOGY, model.MODE_OF_OPERATION, model.YEAR, rule=RateOfFuelUse1_rule)


# def RateOfFuelUse1_rule(model, r, l, f, t, m, y):
#     return (
#         model.RateOfActivity[r, l, t, m, y] * model.InputActivityRatio[r, t, f, m, y]
#         == model.RateOfUseByTechnologyByMode[r, l, t, m, f, y]
#     )


# model.RateOfFuelUse1 = Constraint(
#     model.REGION,
#     model.TIMESLICE,
#     model.FUEL,
#     model.TECHNOLOGY,
#     model.MODE_OF_OPERATION,
#     model.YEAR,
#     rule=RateOfFuelUse1_rule,
# )

# def RateOfFuelUse2_rule(model,r,l,f,t,y):
# return model.RateOfUseByTechnology[r,l,t,f,y] == sum(model.RateOfUseByTechnologyByMode[r,l,t,m,f,y] for m in model.MODE_OF_OPERATION if model.InputActivityRatio[r,t,f,m,y] != 0)
# model.RateOfFuelUse2 = Constraint(model.REGION, model.TIMESLICE, model.FUEL, model.TECHNOLOGY, model.YEAR, rule=RateOfFuelUse2_rule)


# def RateOfFuelUse2_rule(model, r, l, f, t, y):
#     return model.RateOfUseByTechnology[r, l, t, f, y] == sum(
#         model.RateOfUseByTechnologyByMode[r, l, t, m, f, y]
#         for m in model.MODE_OF_OPERATION
#     )


# model.RateOfFuelUse2 = Constraint(
#     model.REGION,
#     model.TIMESLICE,
#     model.FUEL,
#     model.TECHNOLOGY,
#     model.YEAR,
#     rule=RateOfFuelUse2_rule,
# )


# def RateOfFuelUse3_rule(model, r, l, f, y):
#     return (
#         sum(model.RateOfUseByTechnology[r, l, t, f, y] for t in model.TECHNOLOGY)
#         == model.RateOfUse[r, l, f, y]
#     )


# model.RateOfFuelUse3 = Constraint(
#     model.REGION, model.TIMESLICE, model.FUEL, model.YEAR, rule=RateOfFuelUse3_rule
# )


# def EnergyBalanceEachTS1_rule(model, r, l, f, y):
#     return (
#         model.RateOfProduction[r, l, f, y] * model.YearSplit[l, y]
#         == model.Production[r, l, f, y]
#     )


# model.EnergyBalanceEachTS1 = Constraint(
#     model.REGION,
#     model.TIMESLICE,
#     model.FUEL,
#     model.YEAR,
#     rule=EnergyBalanceEachTS1_rule,
# )


# def EnergyBalanceEachTS2_rule(model, r, l, f, y):
#     return model.RateOfUse[r, l, f, y] * model.YearSplit[l, y] == model.Use[r, l, f, y]


# model.EnergyBalanceEachTS2 = Constraint(
#     model.REGION,
#     model.TIMESLICE,
#     model.FUEL,
#     model.YEAR,
#     rule=EnergyBalanceEachTS2_rule,
# )


# def EnergyBalanceEachTS3_rule(model, r, l, f, y):
#     return (
#         model.RateOfDemand[r, l, f, y] * model.YearSplit[l, y]
#         == model.Demand[r, l, f, y]
#     )


# model.EnergyBalanceEachTS3 = Constraint(
#     model.REGION,
#     model.TIMESLICE,
#     model.FUEL,
#     model.YEAR,
#     rule=EnergyBalanceEachTS3_rule,
# )





#########        	Energy Balance B		 	#############

        
# FIXME: 'or ...' in if clause is not existent in standard OSeMOSYS version
def EnergyBalanceEachYear4_rule(model, r, f, y):
    if model.AccumulatedAnnualDemand[r, f, y]!=0 or model.StorageTagFuel[r,f]==1:
        return (sum(model.RateOfActivity[r,l,(t,m),y]
                    *model.OutputActivityRatio[r,t,f,m,y]
                    *model.YearSplit[l,y]
                    for m,t in model.MODETECHNOLOGYFUELOUTind[f]
                    for l in model.TIMESLICE)
            >=  sum(model.RateOfActivity[r,l,(t,m),y]
                        *model.InputActivityRatio[r,t,f,m,y]
                        *model.YearSplit[l,y]
                        for m,t in model.MODETECHNOLOGYFUELINind[f]
                        for l in model.TIMESLICE)
            + sum(
                model.Trade[r, rr,l, f, y] * model.TradeRoute[r, rr, f, y]
                for rr in model.REGION
                for l in model.TIMESLICE
            )
            + model.AccumulatedAnnualDemand[r, f, y]
        )
    else:
        return Constraint.Skip


model.EnergyBalanceEachYear4 = Constraint(
    model.REGION, model.FUEL, model.YEAR, rule=EnergyBalanceEachYear4_rule
)
    
        

# def EnergyBalanceEachYear1_rule(model, r, f, y):
#     return (
#         sum(model.Production[r, l, f, y] for l in model.TIMESLICE)
#         == model.ProductionAnnual[r, f, y]
#     )


# model.EnergyBalanceEachYear1 = Constraint(
#     model.REGION, model.FUEL, model.YEAR, rule=EnergyBalanceEachYear1_rule
# )


# def EnergyBalanceEachYear2_rule(model, r, f, y):
#     return (
#         sum(model.Use[r, l, f, y] for l in model.TIMESLICE) == model.UseAnnual[r, f, y]
#     )


# model.EnergyBalanceEachYear2 = Constraint(
#     model.REGION, model.FUEL, model.YEAR, rule=EnergyBalanceEachYear2_rule
# )


# def EnergyBalanceEachYear3_rule(model, r, rr, f, y):
#     return (
#         sum(model.Trade[r, rr, l, f, y] for l in model.TIMESLICE)
#         == model.TradeAnnual[r, rr, f, y]
#     )


# model.EnergyBalanceEachYear3 = Constraint(
#     model.REGION, model.REGION, model.FUEL, model.YEAR, rule=EnergyBalanceEachYear3_rule
# )




#########        	Accounting Technology Production/Use	#############


# def FuelProductionByTechnology_rule(model, r, l, t, f, y):
#     return (
#         model.RateOfProductionByTechnology[r, l, t, f, y] * model.YearSplit[l, y]
#         == model.ProductionByTechnology[r, l, t, f, y]
#     )


# model.FuelProductionByTechnology = Constraint(
#     model.REGION,
#     model.TIMESLICE,
#     model.TECHNOLOGY,
#     model.FUEL,
#     model.YEAR,
#     rule=FuelProductionByTechnology_rule,
# )


# def FuelUseByTechnology_rule(model, r, l, t, f, y):
#     return (
#         model.RateOfUseByTechnology[r, l, t, f, y] * model.YearSplit[l, y]
#         == model.UseByTechnology[r, l, t, f, y]
#     )


# model.FuelUseByTechnology = Constraint(
#     model.REGION,
#     model.TIMESLICE,
#     model.TECHNOLOGY,
#     model.FUEL,
#     model.YEAR,
#     rule=FuelUseByTechnology_rule,
# )


# def AverageAnnualRateOfActivity_rule(model, r, t, m, y):
#     return (
#         sum(
#             model.RateOfActivity[r, l, t, m, y] * model.YearSplit[l, y]
#             for l in model.TIMESLICE
#         )
#         == model.TotalAnnualTechnologyActivityByMode[r, t, m, y]
#     )


# model.AverageAnnualRateOfActivity = Constraint(
#     model.REGION,
#     model.TECHNOLOGY,
#     model.MODE_OF_OPERATION,
#     model.YEAR,
#     rule=AverageAnnualRateOfActivity_rule,
# )


# def ModelPeriodCostByRegion_rule(model, r):
#     return model.ModelPeriodCostByRegion[r] == sum(
#         model.TotalDiscountedCost[r, y] for y in model.YEAR
#     )


# model.ModelPeriodCostByRegion_constraint = Constraint(
#     model.REGION, rule=ModelPeriodCostByRegion_rule
# )

# def ModelPeriodCost_rule(model):
# return model.ModelPeriodCost == sum(model.ModelPeriodCostByRegion[r] for r in model.REGION)
# model.ModelPeriodCost_Constraint = Constraint(rule=ModelPeriodCost_rule)

#########       	Storage Equations 		     	#############

# FIXME: Implement storage equations and constraints
# s.t. S5_and_S6_StorageLevelYearStart{r in REGION, s in STORAGE, y in YEAR}:
#     if y = min{yy in YEAR} min(yy)
#     then StorageLevelStart[r,s]
#     else StorageLevelYearStart[r,s,y-1] + sum{ls in SEASON, ld in DAYTYPE, lh in DAILYTIMEBRACKET, l in TIMESLICEofSDB[ls,ld,lh]} (sum{(m,t) in MODExTECHNOLOGYperSTORAGEto[s]} (RateOfActivity[r,l,t,m,y-1] * TechnologyToStorage[r,t,s,m]) - (sum{(m,t) in MODExTECHNOLOGYperSTORAGEfrom[s]} RateOfActivity[r,l,t,m,y-1] * TechnologyFromStorage[r,t,s,m])) * YearSplit[l,y-1]
#     = StorageLevelYearStart[r,s,y];
# s.t. S7_and_S8_StorageLevelYearFinish{r in REGION, s in STORAGE, y in YEAR}:
#     if y < max{yy in YEAR} max(yy)
#     then StorageLevelYearStart[r,s,y+1]
#     else StorageLevelYearStart[r,s,y] + sum{ls in SEASON, ld in DAYTYPE, lh in DAILYTIMEBRACKET, l in TIMESLICEofSDB[ls,ld,lh]}  (sum{(m,t) in MODExTECHNOLOGYperSTORAGEto[s]} (RateOfActivity[r,l,t,m,y] * TechnologyToStorage[r,t,s,m]) - (sum{(m,t) in MODExTECHNOLOGYperSTORAGEfrom[s]} RateOfActivity[r,l,t,m,y] * TechnologyFromStorage[r,t,s,m])) * YearSplit[l,y]
#     = StorageLevelYearFinish[r,s,y];
# s.t. S9_and_S10_StorageLevelSeasonStart{r in REGION, s in STORAGE, ls in SEASON, y in YEAR}:
#     if ls = min{lsls in SEASON} min(lsls)
#     then StorageLevelYearStart[r,s,y]
#     else StorageLevelSeasonStart[r,s,ls-1,y] + sum{ld in DAYTYPE, lh in DAILYTIMEBRACKET, l in TIMESLICEofSDB[ls-1,ld,lh]} (sum{(m,t) in MODExTECHNOLOGYperSTORAGEto[s]} (RateOfActivity[r,l,t,m,y] * TechnologyToStorage[r,t,s,m]) - (sum{(m,t) in MODExTECHNOLOGYperSTORAGEfrom[s]} RateOfActivity[r,l,t,m,y] * TechnologyFromStorage[r,t,s,m])) * YearSplit[l,y]
#     = StorageLevelSeasonStart[r,s,ls,y];
# s.t. S11_and_S12_StorageLevelDayTypeStart{r in REGION, s in STORAGE, ls in SEASON, ld in DAYTYPE, y in YEAR}:
#     if ld = min{ldld in DAYTYPE} min(ldld)
#     then StorageLevelSeasonStart[r,s,ls,y]
#     else StorageLevelDayTypeStart[r,s,ls,ld-1,y] + sum{lh in DAILYTIMEBRACKET} (((sum{(m,t) in MODExTECHNOLOGYperSTORAGEto[s], l in TIMESLICEofSDB[ls,ld-1,lh]} RateOfActivity[r,l,t,m,y] * TechnologyToStorage[r,t,s,m]) - (sum{(m,t) in MODExTECHNOLOGYperSTORAGEfrom[s], l in TIMESLICEofSDB[ls,ld-1,lh]} RateOfActivity[r,l,t,m,y] * TechnologyFromStorage[r,t,s,m])) * DaySplit[lh,y]) * DaysInDayType[ls,ld-1,y]
#     = StorageLevelDayTypeStart[r,s,ls,ld,y];
# s.t. S13_and_S14_and_S15_StorageLevelDayTypeFinish{r in REGION, s in STORAGE, ls in SEASON, ld in DAYTYPE, y in YEAR}:
#     if ls = max{lsls in SEASON} max(lsls) && ld = max{ldld in DAYTYPE} max(ldld)
#     then StorageLevelYearFinish[r,s,y]
#     else if ld = max{ldld in DAYTYPE} max(ldld)
#     then StorageLevelSeasonStart[r,s,ls+1,y]
#     else StorageLevelDayTypeFinish[r,s,ls,ld+1,y] - sum{lh in DAILYTIMEBRACKET} (((sum{(m,t) in MODExTECHNOLOGYperSTORAGEto[s], l in TIMESLICEofSDB[ls,ld+1,lh]} RateOfActivity[r,l,t,m,y] * TechnologyToStorage[r,t,s,m]) - (sum{(m,t) in MODExTECHNOLOGYperSTORAGEfrom[s], l in TIMESLICEofSDB[ls,ld+1,lh]} RateOfActivity[r,l,t,m,y] * TechnologyFromStorage[r,t,s,m])) * DaySplit[lh,y]) * DaysInDayType[ls,ld+1,y]
#     = StorageLevelDayTypeFinish[r,s,ls,ld,y];


#########       	Storage Constraints 		     	#############
# s.t. SC1_LowerLimit_BeginningOfDailyTimeBracketOfFirstInstanceOfDayTypeInFirstWeekConstraint{r in REGION, s in STORAGE, ls in SEASON, ld in DAYTYPE, lh in DAILYTIMEBRACKET, y in YEAR}:
#     0 <= (StorageLevelDayTypeStart[r,s,ls,ld,y]+sum{lhlh in DAILYTIMEBRACKET:lh-lhlh>0} (((sum{(m,t) in MODExTECHNOLOGYperSTORAGEto[s], l in TIMESLICEofSDB[ls,ld,lhlh]} RateOfActivity[r,l,t,m,y] * TechnologyToStorage[r,t,s,m]) - (sum{(m,t) in MODExTECHNOLOGYperSTORAGEfrom[s], l in TIMESLICEofSDB[ls,ld,lhlh]} RateOfActivity[r,l,t,m,y] * TechnologyFromStorage[r,t,s,m])) * DaySplit[lhlh,y]))-MinStorageCharge[r,s,y]*(sum{yy in YEAR: y-yy < OperationalLifeStorage[r,s] && y-yy>=0} NewStorageCapacity[r,s,yy]+ResidualStorageCapacity[r,s,y]);
# s.t. SC1_UpperLimit_BeginningOfDailyTimeBracketOfFirstInstanceOfDayTypeInFirstWeekConstraint{r in REGION, s in STORAGE, ls in SEASON, ld in DAYTYPE, lh in DAILYTIMEBRACKET, y in YEAR}:
#     (StorageLevelDayTypeStart[r,s,ls,ld,y]+sum{lhlh in DAILYTIMEBRACKET:lh-lhlh>0} (((sum{(m,t) in MODExTECHNOLOGYperSTORAGEto[s], l in TIMESLICEofSDB[ls,ld,lhlh]} RateOfActivity[r,l,t,m,y] * TechnologyToStorage[r,t,s,m]) - (sum{(m,t) in MODExTECHNOLOGYperSTORAGEfrom[s], l in TIMESLICEofSDB[ls,ld,lhlh]} RateOfActivity[r,l,t,m,y] * TechnologyFromStorage[r,t,s,m])) * DaySplit[lhlh,y]))-(sum{yy in YEAR: y-yy < OperationalLifeStorage[r,s] && y-yy>=0} NewStorageCapacity[r,s,yy]+ResidualStorageCapacity[r,s,y]) <= 0;
# s.t. SC2_LowerLimit_EndOfDailyTimeBracketOfLastInstanceOfDayTypeInFirstWeekConstraint{r in REGION, s in STORAGE, ls in SEASON, ld in DAYTYPE, lh in DAILYTIMEBRACKET, y in YEAR}: 0 <= if ld > min{ldld in DAYTYPE} min(ldld) then (StorageLevelDayTypeStart[r,s,ls,ld,y]-sum{lhlh in DAILYTIMEBRACKET:lh-lhlh<0} (((sum{(m,t) in MODExTECHNOLOGYperSTORAGEto[s], l in TIMESLICEofSDB[ls,ld,lhlh]} RateOfActivity[r,l,t,m,y] * TechnologyToStorage[r,t,s,m]) - (sum{(m,t) in MODExTECHNOLOGYperSTORAGEfrom[s], l in TIMESLICEofSDB[ls,ld,lhlh]} RateOfActivity[r,l,t,m,y] * TechnologyFromStorage[r,t,s,m])) * DaySplit[lhlh,y]))-MinStorageCharge[r,s,y]*(sum{yy in YEAR: y-yy < OperationalLifeStorage[r,s] && y-yy>=0} NewStorageCapacity[r,s,yy]+ResidualStorageCapacity[r,s,y]);
# s.t. SC2_UpperLimit_EndOfDailyTimeBracketOfLastInstanceOfDayTypeInFirstWeekConstraint{r in REGION, s in STORAGE, ls in SEASON, ld in DAYTYPE, lh in DAILYTIMEBRACKET, y in YEAR}: if ld > min{ldld in DAYTYPE} min(ldld) then (StorageLevelDayTypeStart[r,s,ls,ld,y]-sum{lhlh in DAILYTIMEBRACKET:lh-lhlh<0} (((sum{(m,t) in MODExTECHNOLOGYperSTORAGEto[s], l in TIMESLICEofSDB[ls,ld,lhlh]} RateOfActivity[r,l,t,m,y] * TechnologyToStorage[r,t,s,m]) - (sum{(m,t) in MODExTECHNOLOGYperSTORAGEfrom[s], l in TIMESLICEofSDB[ls,ld,lhlh]} RateOfActivity[r,l,t,m,y] * TechnologyFromStorage[r,t,s,m])) * DaySplit[lhlh,y]))-(sum{yy in YEAR: y-yy < OperationalLifeStorage[r,s] && y-yy>=0} NewStorageCapacity[r,s,yy]+ResidualStorageCapacity[r,s,y]) <= 0;
# s.t. SC3_LowerLimit_EndOfDailyTimeBracketOfLastInstanceOfDayTypeInLastWeekConstraint{r in REGION, s in STORAGE, ls in SEASON, ld in DAYTYPE, lh in DAILYTIMEBRACKET, y in YEAR}:  0 <= (StorageLevelDayTypeFinish[r,s,ls,ld,y] - sum{lhlh in DAILYTIMEBRACKET:lh-lhlh<0} (((sum{(m,t) in MODExTECHNOLOGYperSTORAGEto[s], l in TIMESLICEofSDB[ls,ld,lhlh]} RateOfActivity[r,l,t,m,y] * TechnologyToStorage[r,t,s,m]) - (sum{(m,t) in MODExTECHNOLOGYperSTORAGEfrom[s], l in TIMESLICEofSDB[ls,ld,lhlh]} RateOfActivity[r,l,t,m,y] * TechnologyFromStorage[r,t,s,m])) * DaySplit[lhlh,y]))-MinStorageCharge[r,s,y]*(sum{yy in YEAR: y-yy < OperationalLifeStorage[r,s] && y-yy>=0} NewStorageCapacity[r,s,yy]+ResidualStorageCapacity[r,s,y]);
# s.t. SC3_UpperLimit_EndOfDailyTimeBracketOfLastInstanceOfDayTypeInLastWeekConstraint{r in REGION, s in STORAGE, ls in SEASON, ld in DAYTYPE, lh in DAILYTIMEBRACKET, y in YEAR}:  (StorageLevelDayTypeFinish[r,s,ls,ld,y] - sum{lhlh in DAILYTIMEBRACKET:lh-lhlh<0} (((sum{(m,t) in MODExTECHNOLOGYperSTORAGEto[s], l in TIMESLICEofSDB[ls,ld,lhlh]} RateOfActivity[r,l,t,m,y] * TechnologyToStorage[r,t,s,m]) - (sum{(m,t) in MODExTECHNOLOGYperSTORAGEfrom[s], l in TIMESLICEofSDB[ls,ld,lhlh]} RateOfActivity[r,l,t,m,y] * TechnologyFromStorage[r,t,s,m])) * DaySplit[lhlh,y]))-(sum{yy in YEAR: y-yy < OperationalLifeStorage[r,s] && y-yy>=0} NewStorageCapacity[r,s,yy]+ResidualStorageCapacity[r,s,y]) <= 0;
# s.t. SC4_LowerLimit_BeginningOfDailyTimeBracketOfFirstInstanceOfDayTypeInLastWeekConstraint{r in REGION, s in STORAGE, ls in SEASON, ld in DAYTYPE, lh in DAILYTIMEBRACKET, y in YEAR}:         0 <= if ld > min{ldld in DAYTYPE} min(ldld) then (StorageLevelDayTypeFinish[r,s,ls,ld-1,y]+sum{lhlh in DAILYTIMEBRACKET:lh-lhlh>0} (((sum{(m,t) in MODExTECHNOLOGYperSTORAGEto[s], l in TIMESLICEofSDB[ls,ld,lhlh]} RateOfActivity[r,l,t,m,y] * TechnologyToStorage[r,t,s,m]) - (sum{(m,t) in MODExTECHNOLOGYperSTORAGEfrom[s], l in TIMESLICEofSDB[ls,ld,lhlh]} RateOfActivity[r,l,t,m,y] * TechnologyFromStorage[r,t,s,m])) * DaySplit[lhlh,y]))-MinStorageCharge[r,s,y]*(sum{yy in YEAR: y-yy < OperationalLifeStorage[r,s] && y-yy>=0} NewStorageCapacity[r,s,yy]+ResidualStorageCapacity[r,s,y]);
# s.t. SC4_UpperLimit_BeginningOfDailyTimeBracketOfFirstInstanceOfDayTypeInLastWeekConstraint{r in REGION, s in STORAGE, ls in SEASON, ld in DAYTYPE, lh in DAILYTIMEBRACKET, y in YEAR}: if ld > min{ldld in DAYTYPE} min(ldld) then (StorageLevelDayTypeFinish[r,s,ls,ld-1,y]+sum{lhlh in DAILYTIMEBRACKET:lh-lhlh>0} (((sum{(m,t) in MODExTECHNOLOGYperSTORAGEto[s], l in TIMESLICEofSDB[ls,ld,lhlh]} RateOfActivity[r,l,t,m,y] * TechnologyToStorage[r,t,s,m]) - (sum{(m,t) in MODExTECHNOLOGYperSTORAGEfrom[s], l in TIMESLICEofSDB[ls,ld,lhlh]} RateOfActivity[r,l,t,m,y] * TechnologyFromStorage[r,t,s,m])) * DaySplit[lhlh,y]))-(sum{yy in YEAR: y-yy < OperationalLifeStorage[r,s] && y-yy>=0} NewStorageCapacity[r,s,yy]+ResidualStorageCapacity[r,s,y]) <= 0;
# s.t. SC5_MaxChargeConstraint{r in REGION, s in STORAGE, ls in SEASON, ld in DAYTYPE, lh in DAILYTIMEBRACKET, y in YEAR}: sum{(m,t) in MODExTECHNOLOGYperSTORAGEto[s], l in TIMESLICEofSDB[ls,ld,lh]} RateOfActivity[r,l,t,m,y] * TechnologyToStorage[r,t,s,m] <= StorageMaxChargeRate[r,s];
# s.t. SC6_MaxDischargeConstraint{r in REGION, s in STORAGE, ls in SEASON, ld in DAYTYPE, lh in DAILYTIMEBRACKET, y in YEAR}: sum{(m,t) in MODExTECHNOLOGYperSTORAGEfrom[s], l in TIMESLICEofSDB[ls,ld,lh]} RateOfActivity[r,l,t,m,y] * TechnologyFromStorage[r,t,s,m] <= StorageMaxDischargeRate[r,s];
# #

#########       	Capital Costs 		     	#############


# def UndiscountedCapitalInvestment_rule(model, r, t, y):
#     return (
#         model.CapitalCost[r, t, y] * model.NewCapacity[r, t, y]
#         == model.CapitalInvestment[r, t, y]
#     )


# model.UndiscountedCapitalInvestment = Constraint(
#     model.REGION, model.TECHNOLOGY, model.YEAR, rule=UndiscountedCapitalInvestment_rule
# )


# def DiscountedCapitalInvestment_rule(model, r, t, y):
#     return (
#         model.CapitalInvestment[r, t, y]
#         / ((1 + model.DiscountRate[r]) ** (y - min(model.YEAR)))
#         == model.DiscountedCapitalInvestment[r, t, y]
#     )


# model.DiscountedCapitalInvestment_constraint = Constraint(
#     model.REGION, model.TECHNOLOGY, model.YEAR, rule=DiscountedCapitalInvestment_rule
# )


#########        	Operating Costs 		 	#############


# def OperatingCostsVariable_rule(model, r, t, l, y):
#     return (
#         sum(
#             model.TotalAnnualTechnologyActivityByMode[r, t, m, y]
#             * model.VariableCost[r, t, m, y]
#             for m in model.MODE_OF_OPERATION
#         )
#         == model.AnnualVariableOperatingCost[r, t, y]
#     )


# model.OperatingCostsVariable = Constraint(
#     model.REGION,
#     model.TECHNOLOGY,
#     model.TIMESLICE,
#     model.YEAR,
#     rule=OperatingCostsVariable_rule,
# )


# def OperatingCostsFixedAnnual_rule(model, r, t, y):
#     return (
#         model.TotalCapacityAnnual[r, t, y] * model.FixedCost[r, t, y]
#         == model.AnnualFixedOperatingCost[r, t, y]
#     )


# model.OperatingCostsFixedAnnual = Constraint(
#     model.REGION, model.TECHNOLOGY, model.YEAR, rule=OperatingCostsFixedAnnual_rule
# )


# def OperatingCostsTotalAnnual_rule(model, r, t, y):
#     return (
#         model.AnnualFixedOperatingCost[r, t, y]
#         + model.AnnualVariableOperatingCost[r, t, y]
#         == model.OperatingCost[r, t, y]
#     )


# model.OperatingCostsTotalAnnual = Constraint(
#     model.REGION, model.TECHNOLOGY, model.YEAR, rule=OperatingCostsTotalAnnual_rule
# )


# def DiscountedOperatingCostsTotalAnnual_rule(model, r, t, y):
#     return (
#         model.OperatingCost[r, t, y]
#         / ((1 + model.DiscountRate[r]) ** (y - min(model.YEAR) + 0.5))
#         == model.DiscountedOperatingCost[r, t, y]
#     )


# model.DiscountedOperatingCostsTotalAnnual = Constraint(
#     model.REGION,
#     model.TECHNOLOGY,
#     model.YEAR,
#     rule=DiscountedOperatingCostsTotalAnnual_rule,
# )


#########       	Total Discounted Costs	 	#############


# def TotalDiscountedCostByTechnology_rule(model, r, t, y):
#     return (
#         model.DiscountedOperatingCost[r, t, y]
#         + model.DiscountedCapitalInvestment[r, t, y]
#         + model.DiscountedTechnologyEmissionsPenalty[r, t, y]
#         - model.DiscountedSalvageValue[r, t, y]
#         == model.TotalDiscountedCostByTechnology[r, t, y]
#     )


# model.TotalDiscountedCostByTechnology_constraint = Constraint(
#     model.REGION,
#     model.TECHNOLOGY,
#     model.YEAR,
#     rule=TotalDiscountedCostByTechnology_rule,
# )

# def TotalDiscountedCost_rule(model,r,y):
# return sum(model.TotalDiscountedCostByTechnology[r,t,y] for t in model.TECHNOLOGY) + sum(model.TotalDiscountedStorageCost[r,s,y] for s in model.STORAGE) == model.TotalDiscountedCost[r,y]
# model.TotalDiscountedCost_constraint = Constraint(model.REGION, model.YEAR, rule=TotalDiscountedCost_rule)


# def TotalDiscountedCost_rule(model, r, y):
#     return (
#         sum(model.TotalDiscountedCostByTechnology[r, t, y] for t in model.TECHNOLOGY)
#         == model.TotalDiscountedCost[r, y]
#     )


# model.TotalDiscountedCost_constraint = Constraint(
#     model.REGION, model.YEAR, rule=TotalDiscountedCost_rule
# )


#########      		Total Capacity Constraints 	##############


def TotalAnnualMaxCapacityConstraint_rule(model, r, t, y):
    if model.TotalAnnualMaxCapacity[r,t,y]!=-1:
        return ((sum(model.NewCapacity[r,t,yy] 
            for yy in model.TECHNOLOGYCAPYEARS[r,t,y])+ model.ResidualCapacity[r,t,y])
                <= model.TotalAnnualMaxCapacity[r, t, y])
    else:
        return Constraint.Skip

model.TotalAnnualMaxCapacityConstraint = Constraint(
    model.REGION,
    model.TECHNOLOGY,
    model.YEAR,
    rule=TotalAnnualMaxCapacityConstraint_rule,
)


def TotalAnnualMinCapacityConstraint_rule(model, r, t, y):
    if model.TotalAnnualMinCapacity[r,t,y]>0:
        return ((sum(model.NewCapacity[r,t,yy] 
            for yy in model.TECHNOLOGYCAPYEARS[r,t,y])+ model.ResidualCapacity[r,t,y])
                >= model.TotalAnnualMinCapacity[r, t, y])
    else:
        return Constraint.Skip
    

model.TotalAnnualMinCapacityConstraint = Constraint(
    model.REGION,
    model.TECHNOLOGY,
    model.YEAR,
    rule=TotalAnnualMinCapacityConstraint_rule,
)


#########           Salvage Value            	#############

# FIXME: consider removing DepreciationMethod here, as in other version
def SalvageValueAtEndOfPeriod1_rule(model, r, t, y):
    if (
        model.DepreciationMethod[r] == 1
        and ((y + model.OperationalLife[r, t] - 1) > max(model.YEAR))
        and model.DiscountRate[r] > 0
    ):
        return (model.SalvageValue[r, t, y] == model.CapitalCost[r, t, y]
    * model.NewCapacity[r, t, y] * model.CapitalRecoveryFactor[r,t]
        * model.PvAnnuity[r,t]*(
            1
            - (
                ((1 + model.DiscountRate[r]) ** (max(model.YEAR) - y + 1) - 1)
                / ((1 + model.DiscountRate[r]) ** model.OperationalLife[r, t] - 1)
            )
        ))
    elif (
        model.DepreciationMethod[r] == 1
        and ((y + model.OperationalLife[r, t] - 1) > max(model.YEAR))
        and model.DiscountRate[r] == 0
    ) or (
        model.DepreciationMethod[r] == 2
        and (y + model.OperationalLife[r, t] - 1) > (max(model.YEAR))
    ):
        return (model.SalvageValue[r, t, y] == model.CapitalCost[r, t, y]
                * model.NewCapacity[r, t, y] * model.CapitalRecoveryFactor[r,t] 
                * model.PvAnnuity[r,t] * (
            1 - (max(model.YEAR) - y + 1) / model.OperationalLife[r, t]
        ))
    else:
        return model.SalvageValue[r, t, y] == 0
    


model.SalvageValueAtEndOfPeriod1 = Constraint(
    model.REGION, model.TECHNOLOGY, model.YEAR, rule=SalvageValueAtEndOfPeriod1_rule
)


def SalvageValueDiscountedToStartYear_rule(model, r, t, y):
    return model.DiscountedSalvageValue[r, t, y] == model.SalvageValue[r, t, y] / (
        (1 + model.DiscountRate[r]) ** (1 + max(model.YEAR) - min(model.YEAR))
    )


model.SalvageValueDiscountedToStartYear = Constraint(
    model.REGION,
    model.TECHNOLOGY,
    model.YEAR,
    rule=SalvageValueDiscountedToStartYear_rule,
)


#########    		New Capacity Constraints  	##############


def TotalAnnualMaxNewCapacityConstraint_rule(model, r, t, y):
    if model.TotalAnnualMaxCapacityInvestment[r, t, y]!=-1:
        return (model.NewCapacity[r, t, y] 
                <= model.TotalAnnualMaxCapacityInvestment[r, t, y])
    else:
        return Constraint.Skip


model.TotalAnnualMaxNewCapacityConstraint = Constraint(
    model.REGION,
    model.TECHNOLOGY,
    model.YEAR,
    rule=TotalAnnualMaxNewCapacityConstraint_rule,
)


def TotalAnnualMinNewCapacityConstraint_rule(model, r, t, y):
    if model.TotalAnnualMaxCapacityInvestment[r, t, y]>0:
        return (model.NewCapacity[r, t, y] 
                >= model.TotalAnnualMinCapacityInvestment[r, t, y])
    else:
        return Constraint.Skip
        

model.TotalAnnualMinNewCapacityConstraint = Constraint(
    model.REGION,
    model.TECHNOLOGY,
    model.YEAR,
    rule=TotalAnnualMinNewCapacityConstraint_rule,
)



#########   		Annual Activity Constraints	##############


# def TotalAnnualTechnologyActivity_rule(model, r, t, y):
#     return (
#         sum(
#             model.RateOfTotalActivity[r, t, l, y] * model.YearSplit[l, y]
#             for l in model.TIMESLICE
#         )
#         == model.TotalTechnologyAnnualActivity[r, t, y]
#     )


# model.TotalAnnualTechnologyActivity = Constraint(
#     model.REGION, model.TECHNOLOGY, model.YEAR, rule=TotalAnnualTechnologyActivity_rule
# )


def TotalAnnualTechnologyActivityUpperLimit_rule(model, r, t, y):
    if model.TotalTechnologyAnnualActivityUpperLimit[r,t,y]!=-1:
        return (
            sum(
                model.RateOfActivity[r,l,(t,m),y] * model.YearSplit[l, y]
                for l in model.TIMESLICE for m in model.MODETECHNOLOGYind[t]
            )
            <= model.TotalTechnologyAnnualActivityUpperLimit[r, t, y]
        )
    else:
        return Constraint.Skip





model.TotalAnnualTechnologyActivityUpperlimit = Constraint(
    model.REGION,
    model.TECHNOLOGY,
    model.YEAR,
    rule=TotalAnnualTechnologyActivityUpperLimit_rule,
)


def TotalAnnualTechnologyActivityLowerLimit_rule(model, r, t, y):
    if model.TotalTechnologyAnnualActivityLowerLimit[r,t,y]>0:
        return (
            sum(
                model.RateOfActivity[r,l,(t,m),y] * model.YearSplit[l, y]
                for l in model.TIMESLICE for m in model.MODETECHNOLOGYind[t]
            )
            >= model.TotalTechnologyAnnualActivityLowerLimit[r, t, y]
        )
    else:
        return Constraint.Skip
    

model.TotalAnnualTechnologyActivityLowerlimit = Constraint(
    model.REGION,
    model.TECHNOLOGY,
    model.YEAR,
    rule=TotalAnnualTechnologyActivityLowerLimit_rule,
)


#########    		Total Activity Constraints 	##############


# def TotalModelHorizonTechnologyActivity_rule(model, r, t):
#     return (
#         sum(model.TotalTechnologyAnnualActivity[r, t, y] for y in model.YEAR)
#         == model.TotalTechnologyModelPeriodActivity[r, t]
#     )


# model.TotalModelHorizonTechnologyActivity = Constraint(
#     model.REGION, model.TECHNOLOGY, rule=TotalModelHorizonTechnologyActivity_rule
# )


def TotalModelHorizonTechnologyActivityUpperLimit_rule(model, r, t):
    if model.TotalTechnologyModelPeriodActivityUpperLimit[r, t]!=-1:
        return (
            sum(
                model.RateOfActivity[r,l,(t,m),y] * model.YearSplit[l, y]
                for l in model.TIMESLICE
                for m in model.MODETECHNOLOGYind[t]
                for y in model. YEAR
            )
            <= model.TotalTechnologyModelPeriodActivityUpperLimit[r, t]
        )
    else:
        return Constraint.Skip


model.TotalModelHorizonTechnologyActivityUpperLimit = Constraint(
    model.REGION,
    model.TECHNOLOGY,
    rule=TotalModelHorizonTechnologyActivityUpperLimit_rule,
)


def TotalModelHorizonTechnologyActivityLowerLimit_rule(model, r, t):
    if model.TotalTechnologyModelPeriodActivityLowerLimit[r, t]>0:
        return (
            sum(
                model.RateOfActivity[r,l,(t,m),y] * model.YearSplit[l, y]
                for l in model.TIMESLICE
                for m in model.MODETECHNOLOGind[t]
                for y in model. YEAR
            )
            >= model.TotalTechnologyModelPeriodActivityLowerLimit[r, t]
        )
    else:
        return Constraint.Skip


model.TotalModelHorizonTechnologyActivityLowerLimit = Constraint(
    model.REGION,
    model.TECHNOLOGY,
    rule=TotalModelHorizonTechnologyActivityLowerLimit_rule,
)


#########   		RE Production Constraint		##############


# FIXME: added yearsplit on the left of equation, check this is correct
def RE4_EnergyConstraint_rule(model, r, y):
    if model.REMinProductionTarget[r,y]>0:
        return ((model.REMinProductionTarget[r,y]
                *sum(model.RateOfActivity[r,l,(t,m),y]
                    *model.OutputActivityRatio[r,t,f,m,y]
                    *model.YearSplit[l,y]
                    *model.RETagFuel[r,f,y]
                    for m,t in model.MODETECHNOLOGYFUELOUTind[f]
                    for l in model.TIMESLICE
                    for f in model.FUEL))
                    <=
                    (sum(model.RateOfActivity[r,l,(t,m),y]
                        *model.OutputActivityRatio[r,t,f,m,y]
                        *model.YearSplit[l,y]
                        *model.RETagTechnology[r,t,y]
                        for m,t in model.MODETECHNOLOGYFUELOUTind[f]
                        for l in model.TIMESLICE
                        for f in model.FUEL))         
                    )
    else:
        return Constraint.Skip

model.RE4_EnergyConstraint = Constraint(
    model.REGION,
    model.YEAR,
    rule=RE4_EnergyConstraint_rule,
)
#########   		Emissions Accounting		##############


def DiscountedEmissionsPenaltyByTechnology_rule(model, r, t, y):
    return (
        sum(model.RateOfActivity[r,l,(t,m),y]
            *model.EmissionActivityRatio[r,t,e,m,y]
            *model.YearSplit[l,y]
            *model.EmissionsPenalty[r,e,y]
            /model.DiscountFactorMid[r,y]
            for e in model.EMISSION
            for m,tt in model.MODETECHNOLOGYEMISSIONind[e] if (t==tt)
            for l in model.TIMESLICE
            )
        
        == model.DiscountedTechnologyEmissionsPenalty[r, t, y]
    )


model.DiscountedEmissionsPenaltyByTechnology = Constraint(
    model.REGION,
    model.TECHNOLOGY,
    model.YEAR,
    rule=DiscountedEmissionsPenaltyByTechnology_rule,
)


# FIXME: generation of this this and rule below fails under certain
# circumstances when, e.g., emission is not generated by any tech and thus not
# but limit not set to -1 by default and thus this is generated but does not
#  returns not a valid expression
def AnnualEmissionsLimit_rule(model, r, e, y):
    if model.AnnualEmissionLimit[r, e, y] !=-1:
        return (
            (sum(model.RateOfActivity[r,l,(t,m),y]
                *model.EmissionActivityRatio[r,t,e,m,y]
                *model.YearSplit[l,y]
                for m,t in model.MODETECHNOLOGYEMISSIONind[e]
                for l in model.TIMESLICE)
            +model.AnnualExogenousEmission[r,e,y])
            <= model.AnnualEmissionLimit[r, e, y]
        )
    else:
        return Constraint.Skip


model.AnnualEmissionsLimit = Constraint(
    model.REGION, model.EMISSION, model.YEAR, rule=AnnualEmissionsLimit_rule
)

def ModelPeriodEmissionsLimit_rule(model, r, e):
    if model.ModelPeriodEmissionLimit[r, e] !=-1:
        return (
            (sum(model.RateOfActivity[r,l,(t,m),y]
                *model.EmissionActivityRatio[r,t,e,m,y]
                *model.YearSplit[l,y]
                for m,t in model.MODETECHNOLOGYEMISSIONind[e]
                for l in model.TIMESLICE
                for y in model.YEAR)
            + model.ModelPeriodExogenousEmission[r,e])
            <= model.ModelPeriodEmissionLimit[r, e]
        )
    else:
        return Constraint.Skip


model.ModelPeriodEmissionsLimit = Constraint(
    model.REGION, model.EMISSION, rule=ModelPeriodEmissionsLimit_rule
)




# def AnnualEmissionProductionByMode_rule(model, r, t, e, m, y):
#     if model.EmissionActivityRatio[r, t, e, m, y] != 0:
#         return (
#             model.EmissionActivityRatio[r, t, e, m, y]
#             * model.TotalAnnualTechnologyActivityByMode[r, t, m, y]
#             == model.AnnualTechnologyEmissionByMode[r, t, e, m, y]
#         )
#     else:
#         return model.AnnualTechnologyEmissionByMode[r, t, e, m, y] == 0


# model.AnnualEmissionProductionByMode = Constraint(
#     model.REGION,
#     model.TECHNOLOGY,
#     model.EMISSION,
#     model.MODE_OF_OPERATION,
#     model.YEAR,
#     rule=AnnualEmissionProductionByMode_rule,
# )


# def AnnualEmissionProduction_rule(model, r, t, e, y):
#     return (
#         sum(
#             model.AnnualTechnologyEmissionByMode[r, t, e, m, y]
#             for m in model.MODE_OF_OPERATION
#         )
#         == model.AnnualTechnologyEmission[r, t, e, y]
#     )


# model.AnnualEmissionProduction = Constraint(
#     model.REGION,
#     model.TECHNOLOGY,
#     model.EMISSION,
#     model.YEAR,
#     rule=AnnualEmissionProduction_rule,
# )


# def EmissionPenaltyByTechAndEmission_rule(model, r, t, e, y):
#     return (
#         model.AnnualTechnologyEmission[r, t, e, y] * model.EmissionsPenalty[r, e, y]
#         == model.AnnualTechnologyEmissionPenaltyByEmission[r, t, e, y]
#     )


# model.EmissionPenaltyByTechAndEmission = Constraint(
#     model.REGION,
#     model.TECHNOLOGY,
#     model.EMISSION,
#     model.YEAR,
#     rule=EmissionPenaltyByTechAndEmission_rule,
# )


# def EmissionsPenaltyByTechnology_rule(model, r, t, y):
#     return (
#         sum(
#             model.AnnualTechnologyEmissionPenaltyByEmission[r, t, e, y]
#             for e in model.EMISSION
#         )
#         == model.AnnualTechnologyEmissionsPenalty[r, t, y]
#     )


# model.EmissionsPenaltyByTechnology = Constraint(
#     model.REGION, model.TECHNOLOGY, model.YEAR, rule=EmissionsPenaltyByTechnology_rule
# )



# def EmissionsAccounting1_rule(model, r, e, y):
#     return (
#         sum(model.AnnualTechnologyEmission[r, t, e, y] for t in model.TECHNOLOGY)
#         == model.AnnualEmissions[r, e, y]
#     )


# model.EmissionsAccounting1 = Constraint(
#     model.REGION, model.EMISSION, model.YEAR, rule=EmissionsAccounting1_rule
# )


# def EmissionsAccounting2_rule(model, r, e):
#     return (
#         sum(model.AnnualEmissions[r, e, y] for y in model.YEAR)
#         == model.ModelPeriodEmissions[r, e] - model.ModelPeriodExogenousEmission[r, e]
#     )


# model.EmissionsAccounting2 = Constraint(
#     model.REGION, model.EMISSION, rule=EmissionsAccounting2_rule
# )



#########   		Reserve Margin Constraint	############## 


def ReserveMarginConstraint_rule(model, r, l, y):
    if model.ReserveMargin[r,y]!=-1:
        return (
            sum(model.RateOfActivity[r,l,(t,m),y]
                * model.OutputActivityRatio[r,t,f,m,y]
                * model.ReserveMarginTagFuel[r,f,y]
                * model.ReserveMargin[r,y]
                for f in model.FUEL
                for m,t in model.MODETECHNOLOGYFUELOUTind[f]
                )
            <= sum( 
                    (sum(model.NewCapacity[r,t,yy] 
                        for yy in model.TECHNOLOGYCAPYEARS[r,t,y])+ model.ResidualCapacity[r,t,y])
                    * model.ReserveMarginTagTechnology[r,t,y]
                    * model.CapacityToActivityUnit[r,t]
                    for t in model.TECHNOLOGY)
        )
    else:
        return Constraint.Skip


ReserveMarginConstraint = Constraint(
    model.REGION, model.TIMESLICE, model.YEAR, rule=ReserveMarginConstraint_rule
)




# def ReserveMargin_TechnologiesIncluded_rule(model, r, l, y):
#     return (
#         sum(
#             (
#                 model.TotalAnnualCapacity[r, t, y]
#                 * model.ReserveMarginTagTechnology[r, t, y]
#                 * model.CapacityToActivityUnit[r, t]
#             )
#             for t in model.TECHNOLOGY
#         )
#         == model.TotalCapacityInReserveMargin[r, y]
#     )


# ReserveMargin_TechnologiesIncluded = Constraint(
#     model.REGION,
#     model.TIMESLICE,
#     model.YEAR,
#     rule=ReserveMargin_TechnologiesIncluded_rule,
# )


# def ReserveMargin_FuelsIncluded_rule(model, r, l, y):
#     return (
#         sum(
#             (model.RateOfProduction[r, l, f, y] * model.ReserveMarginTagFuel)
#             for f in model.FUEL
#         )
#         == model.DemandNeedingReserveMargin[r, l, y]
#     )


# ReserveMargin_FuelsIncluded = Constraint(
#     model.REGION, model.TIMESLICE, model.YEAR, rule=ReserveMargin_FuelsIncluded_rule
# )

#%%
# added constraints


def UDC1_UserDefinedConstraintInequality_rule(model, r, u, y):
    if (model.UDCTag[r,u]==0) and (model.UDCConstant[r,u,y]!=-1):
        return ((sum(model.UDCMultiplierTotalCapacity[r,t,u,y]
                     *(sum(model.NewCapacity[r,t,yy] 
                           for yy in model.TECHNOLOGYCAPYEARS[r,t,y])+ model.ResidualCapacity[r,t,y])  
                     for t in model.TECHNOLOGY)
                +sum(model.UDCMultiplierNewCapacity[r,t,u,y]
                     *model.NewCapacity[r,t,y] 
                     for t in model.TECHNOLOGY)
                 +sum(model.UDCMultiplierActivity[r,t,m,u,y]
                     *sum(model.RateOfActivity[r,l,(t,m),y] * model.YearSplit[l, y]
                          for l in model.TIMESLICE)
                     for t,m in model.MODETECHNOLOGY))
                <= model.UDCConstant[r,u,y])
    if (model.UDCTag[r,u]==1) and (model.UDCConstant[r,u,y]!=-1):
        return ((sum(model.UDCMultiplierTotalCapacity[r,t,u,y]
                     *(sum(model.NewCapacity[r,t,yy] 
                           for yy in model.TECHNOLOGYCAPYEARS[r,t,y])+ model.ResidualCapacity[r,t,y])  
                     for t in model.TECHNOLOGY)
                +sum(model.UDCMultiplierNewCapacity[r,t,u,y]
                     *model.NewCapacity[r,t,y] 
                     for t in model.TECHNOLOGY)
                 +sum(model.UDCMultiplierActivity[r,t,m,u,y]
                     *sum(model.RateOfActivity[r,l,(t,m),y] * model.YearSplit[l, y]
                          for l in model.TIMESLICE)
                     for t,m in model.MODETECHNOLOGY))
                >= model.UDCConstant[r,u,y])                       
    else:
        return Constraint.Skip

model.UDC1_UserDefinedConstraintInequality = Constraint(
    model.REGION,
    model.UDC,
    model.YEAR,
    rule=UDC1_UserDefinedConstraintInequality_rule,
)


#%%
# additions for UK model

def CapacityFractionConstraint_rule(model, r, t, y):
    if model.CapacityFractionTagTechnology[r,t]!=0:
        return (((sum(model.NewCapacity[r,t,yy] 
            for yy in model.TECHNOLOGYCAPYEARS[r,t,y]))#+ model.ResidualCapacity[rr,t,y])
                    )
                >= (sum((sum(model.NewCapacity[r,tt,yy] 
                    for yy in model.TECHNOLOGYCAPYEARS[r,tt,y]))#+ model.ResidualCapacity[r,tt,y])
                        * model.CapacityFractionTagTechnology[r,tt]
                            for tt in model.TECHNOLOGY)
                    *model.CapacityFraction[r,t,y]
                    *0.95))
    else:
        return Constraint.Skip

model.CapacityFractionConstraint = Constraint(
    model.REGION,
    model.TECHNOLOGY,
    model.YEAR,
    rule=CapacityFractionConstraint_rule,
)



