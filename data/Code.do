*** Restart ***

clear
set pagesize 10000

*** Set Working Directory ***

cd "" /* Update directory */

*** Import Data ***

import excel using "Data.xlsx", firstrow

*** Generate Variables ***

tsset Year, yearly

foreach var in BroadMoneySupply ConsumerPriceIndex EquityPrices AgriculturalOutput IndustrialOutput ServicesOutput Exports Imports GDP GDPUS {

gen `var'Growth = ((`var'/l.`var')-1)*100

}

*** Set Parameters ***

local Impulse Crisis /* Crisis, CrisisUS, ReinhartandRogoff, SchularickandTaylor, Turner, Exogenous, CrisisUnweighted, Crisis85th, Crisis95th, CrisisNoRuns, CrisisScaled */
local Response GDPGrowth /* GDPGrowth, GDPUSGrowth, BroadMoneySupplyGrowth, EquityPricesGrowth, ExportsGrowth, ImportsGrowth, AgriculturalOutputGrowth, IndustrialOutputGrowth, ServicesOutputGrowth */
local Control /* BankRate, GovernmentRevenue, GovernmentSpending, ConsolYield, ConsumerPriceIndexGrowth, EquityPricesGrowth, BankRate GovernmentRevenue GovernmentSpending ConsolYield ConsumerPriceIndexGrowth EquityPricesGrowth */
local Endogenous `Impulse' `Response' `Control' /* `Impulse' `Response' `Control', `Response' `Impulse' `Control' */
local Lags 3 /* 1, 3, 5 */
local Start 1750 /* 1750, 1800, 1826, 1870, 1939 */
local End 1938 /* 1913, 1915, 1938, 2008 */

*** VARs ***

var `Endogenous' if Year>=`Start' & Year<=`End', lags(1/`Lags')

irf create irf, step(5) set(irf, replace)

irf table oirf, impulse(`Impulse') response(`Impulse' `Response') level(90) stderror
