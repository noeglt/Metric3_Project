.libPaths("C:/Users/camil/Documents/ENS/Cours/Econométrie/packages")
setwd("C:/Users/camil/Documents/ENS/Cours/working directory")

library(fredr)       
library(vars)        
library(tseries)     
library(urca)        
library(ggplot2)     
library(dplyr)      
library(tidyr)       
library(lubridate)
library(janitor)

#Data loading and basic cleaning-------------------------------------------------------------------------------

fredr_set_key(Sys.getenv("FRED_API_KEY"))

production = fredr(series_id = "IPG21112S", 
                              observation_start = as.Date("1973-01-01"),
                              observation_end   = as.Date("2026-02-02"))

index = fredr(series_id = "IGREA", 
                     observation_start = as.Date("1973-01-01"),
                     observation_end   = as.Date("2026-02-02"))

price = read.csv("U.S._Crude_Oil_Imported_Acquisition_Cost_by_Refiners.csv") #Not from FRED because we extrapolate for 1973 

#LLM helped for the data cleaning

price = price %>%
  slice(5:n()) %>%
  select(Month = 1, Cost = 2)

price = price %>%
  mutate(
    Month = myd(paste(Month, "01")), 
    Cost = as.numeric(Cost)
  ) %>%
  arrange(Month) %>% 
  rename(date = Month)


df <- left_join(price, index, by = "date")

df <- df %>% rename(index = value)

df <- left_join(df, production, by = "date")

df <- df %>% rename(prod = value)

df <- df %>% select(date, index, prod, Cost)

df <- df %>% mutate(cost = log(Cost)) %>% 
  mutate(growth_prod = (prod - lag(prod)) / lag(prod)
  ) %>% filter(!is.na(growth_prod))

df <- df %>% select(cost, date, index, growth_prod)

#Visualize the series-------------------------------------------------------


#Stationarity tests----------------------------------------------------


# Three complementary tests:
#   PP     (Phillips-Perron): H0 = unit root. Like ADF but uses a
#          nonparametric correction for autocorrelation and
#          heteroskedasticity in the errors rather than adding lags.
#   DF-GLS (Elliott-Rothenberg-Stock): H0 = unit root. Applies a GLS
#          demeaning / detrending step first; typically more powerful than
#          ADF/PP when the root is close to one.
#   KPSS   (Kwiatkowski-Phillips-Schmidt-Shin): H0 = STATIONARY (null
#          flipped!). Reject -> evidence against stationarity.
#
# Rule of thumb: if PP and DF-GLS fail to reject AND KPSS rejects, the evidence
# for a unit root is strong. If they disagree, the picture is mixed.
# Flag it in the diary and sanity-check with differencing.


# 5. Lag order selection -------------------------------------------------------


# 6. Estimate the VAR ----------------------------------------------------------


# 7. Stability & residual diagnostics ------------------------------------------


# 8. Impulse Response Functions ------------------------------------------------

# 9. Forecast Error Variance Decomposition ------------------------------------


# 10. Granger causality (extension) -------------------------------------------

