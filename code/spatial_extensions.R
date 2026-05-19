setwd("C:/Users/taher/OneDrive/Bureau/econometrics 3")
install.packages(c(
  "fredr",
  "vars",
  "tseries",
  "urca",
  "ggplot2",
  "dplyr",
  "tidyr",
  "lubridate",
  "janitor",
  "readr",
  "tidyverse",
  "svars",
  "Cairo"
))

library(fredr)       
library(vars)        
library(tseries)     
library(urca)        
library(ggplot2)     
library(dplyr)      
library(tidyr)       
library(lubridate)
library(janitor)
library(lubridate)
library(readr)
library(tidyverse)
library(svars)
library(Cairo)

#Data loading and basic cleaning-------------------------------------------------------------------------------

fredr_set_key(("ab808af2bcc04e9a32c684991350a0cb"))

kilian = fredr(
  series_id = "IGREA", 
  observation_start = as.Date("1973-01-01"),
  observation_end   = as.Date("2026-02-02"))

df_clean2 <- read_csv("Crude_prod_TS.csv")

df_clean2 <- df_clean2 |>
  select("...5", "...6")

df_clean2 <- df_clean2 %>%
  slice(-(1:4))

df_clean2 <- df_clean2 %>%
  rename(date = "...5", value = "...6") %>%
  mutate(
    date = ym(date),
    value = as.numeric(value)
  ) %>%
  arrange(date)

df_clean2 <- df_clean2 %>%rename(prod = value)


price = read.csv("U.S._Crude_Oil_Imported_Acquisition_Cost_by_Refiners.csv") #Not from FRED because we extrapolate for 1973 

#LLM helped for the data cleaning

price = price |>
  slice(5:n()) |>
  select(Month = 1, Cost = 2)

price = price |>
  mutate(
    Month = myd(paste(Month, "01")), 
    Cost = as.numeric(Cost)
  ) |>
  arrange(Month) |> 
  rename(date = Month)

df_main <- left_join(price, kilian, by = "date") |>
  rename(index = value)  |>
  left_join(df_clean2, by = "date") |>
  select(date, index, prod, Cost) |>
  mutate(
    growth_prod = (log(prod) - log(lag(prod)))*100) 

#Now let us deflate the cost by the CPI index to get the real cost of oil (in which year USD?)

cpi = fredr(
  series_id = "CPIAUCSL",
  observation_start = as.Date("1973-01-01")
)

cpi <- cpi |>
  select(date, cpi = value)

df_main <- df_main |>
  left_join(cpi, by = "date")

df_main <- df_main|>
  mutate( real_price = (Cost / cpi)*100,
          log_real_price = log(real_price))

#FRANCE

indus_prod_fr<- fredr(
  series_id = "FRAPROINDMISMEI", 
  observation_start = as.Date("1973-01-01"))

indus_prod_fr <- indus_prod_fr|> select(date,value) 

df_france <- df_main |> left_join(indus_prod_fr, by = "date")

#CANADA

indus_prod_can <- fredr(
  series_id = "CANPROINDMISMEI",
  observation_start = as.Date("1973-01-01"))

indus_prod_can <- indus_prod_can |> select(date, value)

df_canada <- df_main |> left_join(indus_prod_can, by = "date")

#JAPON 

indus_prod_jpn <- fredr(
  series_id = "JPNPROINDMISMEI",
  observation_start = as.Date("1973-01-01"))

indus_prod_jpn <- indus_prod_jpn |> select(date, value)

df_japan  <- df_main |> left_join(indus_prod_jpn, by = "date")

#NORWAY

indus_prod_nor <- fredr(
  series_id = "NORPROINDMISMEI",
  observation_start = as.Date("1973-01-01"))

indus_prod_nor <- indus_prod_nor |> select(date, value)

df_norway <- df_main |> left_join(indus_prod_nor, by = "date")

build_svar_data <- function(data){
  df <- data %>%
    select(date,log_real_price,growth_prod,index_kilian=index,index_industrial=value) %>%
    filter(date<as.Date("2026-01-01")) %>%
    filter(!is.na(growth_prod)) %>%
    na.omit()
  return(df)}

df_SVAR_fr<- build_svar_data(df_france)
df_SVAR_can <-build_svar_data(df_canada)
df_SVAR_jpn <- build_svar_data(df_japan) 
df_SVAR_nor <- build_svar_data(df_norway)

run_ur_tests<-  function(data){

dfbis <- data %>% select(!date)

ur_tbl <- data.frame(                                  # Empty container
  variable   = colnames(dfbis),
  pp_stat    = NA_real_, pp_pval   = NA_real_,
  dfgls_stat = NA_real_, dfgls_cv5 = NA_real_,
  kpss_stat  = NA_real_, kpss_pval = NA_real_
)
for (j in seq_along(colnames(dfbis))) {               # Loop over 3 variables
  y <- dfbis[, j]                                     # Pick the j-th series
  
  pp  <- pp.test(y, alternative = "stationary")        # Phillips-Perron
  ers <- ur.ers(y, type = "DF-GLS", model = "constant")# DF-GLS (ERS 1996)
  kp  <- kpss.test(y, null = "Level")                  # KPSS (null: stationary)
  
  ur_tbl$pp_stat[j]    <- unname(pp$statistic)         # Store PP stat / p-value
  ur_tbl$pp_pval[j]    <- pp$p.value
  ur_tbl$dfgls_stat[j] <- as.numeric(ers@teststat)     # Uses 5% crit. val.
  ur_tbl$dfgls_cv5[j]  <- ers@cval[, "5pct"]
  ur_tbl$kpss_stat[j]  <- unname(kp$statistic)         # KPSS stat / p-value
  ur_tbl$kpss_pval[j]  <- kp$p.value
}
return(ur_tbl)}

ur_tbl_fr <- run_ur_tests(df_SVAR_fr)
ur_tbl_can <- run_ur_tests(df_SVAR_can)
ur_tbl_jpn <- run_ur_tests(df_SVAR_jpn)
ur_tbl_nor <- run_ur_tests(df_SVAR_nor)

#au cas où 
dfbis_fr <- df_SVAR_fr%>% select(!date)
dfbis_jpn <- df_SVAR_jpn%>% select(!date)
dfbis_nor <- df_SVAR_nor%>% select(!date)
dfbis_can <- df_SVAR_can%>% select(!date)

lag_select <- function(data, lag_max = 24){
  dfbis <- data %>% select(-date)
  lag_sel <- VARselect(dfbis, lag.max = lag_max, type = "const")
  print(lag_sel$selection)
  print(lag_sel$criteria)
  return(lag_sel)}

lag_select(df_SVAR_fr)
lag_select(df_SVAR_jpn)
lag_select(df_SVAR_nor)
lag_select(df_SVAR_can)

var_fit_fr  <- VAR(dfbis_fr,  p = 3, type = "const")
summary(var_fit_fr)

var_fit_jpn <- VAR(dfbis_jpn, p = 3, type = "const")
summary(var_fit_jpn)

var_fit_nor <- VAR(dfbis_nor, p = 3, type = "const")
summary(var_fit_nor)

var_fit_can <- VAR(dfbis_can, p = 3, type = "const")
summary(var_fit_can)

# =========================
# FRANCE
# =========================
roots_mod_fr <- roots(var_fit_fr, modulus = TRUE)
cat("FR max modulus:", round(max(roots_mod_fr), 3), "\n")

sc_test_fr <- serial.test(var_fit_fr, lags.pt = 36, type = "PT.asymptotic")
print(sc_test_fr)

nm_test_fr <- normality.test(var_fit_fr, multivariate.only = TRUE)
print(nm_test_fr$jb.mul)


# =========================
# CANADA
# =========================
roots_mod_can <- roots(var_fit_can, modulus = TRUE)
cat("CAN max modulus:", round(max(roots_mod_can), 3), "\n")

sc_test_can <- serial.test(var_fit_can, lags.pt = 36, type = "PT.asymptotic")
print(sc_test_can)

nm_test_can <- normality.test(var_fit_can, multivariate.only = TRUE)
print(nm_test_can$jb.mul)


# =========================
# JAPON
# =========================
roots_mod_jpn <- roots(var_fit_jpn, modulus = TRUE)
cat("JPN max modulus:", round(max(roots_mod_jpn), 3), "\n")

sc_test_jpn <- serial.test(var_fit_jpn, lags.pt = 36, type = "PT.asymptotic")
print(sc_test_jpn)

nm_test_jpn <- normality.test(var_fit_jpn, multivariate.only = TRUE)
print(nm_test_jpn$jb.mul)


# =========================
# NORVEGE
# =========================
roots_mod_nor <- roots(var_fit_nor, modulus = TRUE)
cat("NOR max modulus:", round(max(roots_mod_nor), 3), "\n")

sc_test_nor <- serial.test(var_fit_nor, lags.pt = 36, type = "PT.asymptotic")
print(sc_test_nor)

nm_test_nor <- normality.test(var_fit_nor, multivariate.only = TRUE)
print(nm_test_nor$jb.mul)

H_fr <- 18
H_can <- 18
H_jpn <- 18
H_nor <- 18

vars <- c("growth_prod","index_kilian","log_real_price","index_industrial")

df_ordered_fr <- dfbis_fr[, c("growth_prod","index_kilian","log_real_price","index_industrial")]
var_ordered_fr <- VAR(df_ordered_fr, p = 3, type = "const")
vars_ordered_fr <- c("growth_prod","index_kilian","log_real_price","index_industrial")

df_ordered_can <- dfbis_can[, c("growth_prod","index_kilian","log_real_price","index_industrial")]
var_ordered_can <- VAR(df_ordered_can, p = 3, type = "const")
vars_ordered_can <- c("growth_prod","index_kilian","log_real_price","index_industrial")

df_ordered_jpn <- dfbis_jpn[, c("growth_prod","index_kilian","log_real_price","index_industrial")]
var_ordered_jpn <- VAR(df_ordered_jpn, p = 3, type = "const")
vars_ordered_jpn <- c("growth_prod","index_kilian","log_real_price","index_industrial")

df_ordered_nor <- dfbis_nor[, c("growth_prod","index_kilian","log_real_price","index_industrial")]
var_ordered_nor <- VAR(df_ordered_nor, p = 3, type = "const")
vars_ordered_nor <- c("growth_prod","index_kilian","log_real_price","index_industrial")

shock_names <- c(
  growth_prod = "Oil supply shock",
  index_kilian = "Aggregate demand shock",
  log_real_price = "Oil-specific demand shock")

ylabs <- c(
  growth_prod = "Oil production",
  index_kilian = "Real activity",
  log_real_price = "Real price of oil",
  index_industrial="industrial production"
)

irf_68_std_fr <- lapply(vars, \(s) irf(var_ordered_fr, impulse = s, response = vars,
                                          n.ahead = H_fr, boot = TRUE, ci = 0.68, runs = 100, cumulative = FALSE))

irf_95_std_fr <- lapply(vars, \(s) irf(var_ordered_fr, impulse = s, response = vars,
                                          n.ahead = H_fr, boot = TRUE, ci = 0.95, runs = 100, cumulative = FALSE))

irf_68_cum_fr <- lapply(vars, \(s) irf(var_ordered_fr, impulse = s, response = vars,
                                          n.ahead = H_fr, boot = TRUE, ci = 0.68, runs = 100, cumulative = TRUE))

irf_95_cum_fr <- lapply(vars, \(s) irf(var_ordered_fr, impulse = s, response = vars,
                                          n.ahead = H_fr, boot = TRUE, ci = 0.95, runs = 100, cumulative = TRUE))

names(irf_68_std_fr) <- vars; names(irf_95_std_fr) <- vars
names(irf_68_cum_fr) <- vars; names(irf_95_cum_fr) <- vars

irf_68_std_can <- lapply(vars, \(s) irf(var_ordered_can, impulse = s, response = vars,
                                            n.ahead = H_can, boot = TRUE, ci = 0.68, runs = 100, cumulative = FALSE))

irf_95_std_can <- lapply(vars, \(s) irf(var_ordered_can, impulse = s, response = vars,
                                            n.ahead = H_can, boot = TRUE, ci = 0.95, runs = 100, cumulative = FALSE))

irf_68_cum_can <- lapply(vars, \(s) irf(var_ordered_can, impulse = s, response = vars,
                                            n.ahead = H_can, boot = TRUE, ci = 0.68, runs = 100, cumulative = TRUE))

irf_95_cum_can <- lapply(vars, \(s) irf(var_ordered_can, impulse = s, response = vars,
                                            n.ahead = H_can, boot = TRUE, ci = 0.95, runs = 100, cumulative = TRUE))

names(irf_68_std_can) <- vars; names(irf_95_std_can) <- vars
names(irf_68_cum_can) <- vars; names(irf_95_cum_can) <- vars

irf_68_std_jpn <- lapply(vars, \(s) irf(var_ordered_jpn, impulse = s, response = vars,
                                            n.ahead = H_jpn, boot = TRUE, ci = 0.68, runs = 100, cumulative = FALSE))

irf_95_std_jpn <- lapply(vars, \(s) irf(var_ordered_jpn, impulse = s, response = vars,
                                            n.ahead = H_jpn, boot = TRUE, ci = 0.95, runs = 100, cumulative = FALSE))

irf_68_cum_jpn <- lapply(vars, \(s) irf(var_ordered_jpn, impulse = s, response = vars,
                                            n.ahead = H_jpn, boot = TRUE, ci = 0.68, runs = 100, cumulative = TRUE))

irf_95_cum_jpn <- lapply(vars, \(s) irf(var_ordered_jpn, impulse = s, response = vars,
                                            n.ahead = H_jpn, boot = TRUE, ci = 0.95, runs = 100, cumulative = TRUE))

names(irf_68_std_jpn) <- vars; names(irf_95_std_jpn) <- vars
names(irf_68_cum_jpn) <- vars; names(irf_95_cum_jpn) <- vars

irf_68_std_nor <- lapply(vars, \(s) irf(var_ordered_nor, impulse = s, response = vars,
                                            n.ahead = H_nor, boot = TRUE, ci = 0.68, runs = 100, cumulative = FALSE))

irf_95_std_nor <- lapply(vars, \(s) irf(var_ordered_nor, impulse = s, response = vars,
                                            n.ahead = H_nor, boot = TRUE, ci = 0.95, runs = 100, cumulative = FALSE))

irf_68_cum_nor <- lapply(vars, \(s) irf(var_ordered_nor, impulse = s, response = vars,
                                            n.ahead = H_nor, boot = TRUE, ci = 0.68, runs = 100, cumulative = TRUE))

irf_95_cum_nor <- lapply(vars, \(s) irf(var_ordered_nor, impulse = s, response = vars,
                                            n.ahead = H_nor, boot = TRUE, ci = 0.95, runs = 100, cumulative = TRUE))

names(irf_68_std_nor) <- vars; names(irf_95_std_nor) <- vars
names(irf_68_cum_nor) <- vars; names(irf_95_cum_nor) <- vars

plot_panel2 <- function(shock,response,country=c("fr","can","jpn","nor")) {
  country <- match.arg(country)
  if(country=="fr"){
    irf_68_std <- irf_68_std_fr
    irf_95_std <- irf_95_std_fr
    irf_68_cum <- irf_68_cum_fr
    irf_95_cum <- irf_95_cum_fr
    H <- H_fr
  }

  if(country=="can"){
    irf_68_std <- irf_68_std_can
    irf_95_std <- irf_95_std_can
    irf_68_cum <- irf_68_cum_can
    irf_95_cum <- irf_95_cum_can
    H <- H_can
  }

  if(country=="jpn"){
    irf_68_std <- irf_68_std_jpn
    irf_95_std <- irf_95_std_jpn
    irf_68_cum <- irf_68_cum_jpn
    irf_95_cum <- irf_95_cum_jpn
    H <- H_jpn
  }

  if(country=="nor"){
    irf_68_std <- irf_68_std_nor
    irf_95_std <- irf_95_std_nor
    irf_68_cum <- irf_68_cum_nor
    irf_95_cum <- irf_95_cum_nor
    H <- H_nor
  }

  if(response=="growth_prod"){
    obj_68 <- irf_68_cum
    obj_95 <- irf_95_cum
  } else {
    obj_68 <- irf_68_std
    obj_95 <- irf_95_std
  }

  y  <- obj_68[[shock]]$irf[[shock]][, response]
  l1 <- obj_68[[shock]]$Lower[[shock]][, response]
  u1 <- obj_68[[shock]]$Upper[[shock]][, response]
  l2 <- obj_95[[shock]]$Lower[[shock]][, response]
  u2 <- obj_95[[shock]]$Upper[[shock]][, response]

  if(shock=="growth_prod"){
    y <- -y
    old_l1 <- l1; l1 <- -u1; u1 <- -old_l1
    old_l2 <- l2; l2 <- -u2; u2 <- -old_l2
  }

  if(response=="log_real_price"){
    y <- 100*y
    l1 <- 100*l1
    u1 <- 100*u1
    l2 <- 100*l2
    u2 <- 100*u2
  }

  x <- 0:H

  plot(x,y,
       type="l",
       lwd=2,
       main=paste0(shock_names[shock]," - ",toupper(country)),
       ylab=unname(ylabs[response]),
       xlab="",
       ylim=range(c(y,l1,u1,l2,u2),na.rm=TRUE))

  abline(h=0,col="gray")

  lines(x,l1,lty=2)
  lines(x,u1,lty=2)

  lines(x,l2,lty=3)
  lines(x,u2,lty=3)
}

par(mfrow = c(4,4), mar = c(3,4,3,1))
for(shock in vars){
  for(response in vars){
    plot_panel2(shock,response,country="fr")
  }
}
par(mfrow = c(1,1))

par(mfrow = c(4,4), mar = c(3,4,3,1))
for(shock in vars){
  for(response in vars){
    plot_panel2(shock,response,country="can")
  }
}
par(mfrow = c(1,1))


par(mfrow = c(4,4), mar = c(3,4,3,1))
for(shock in vars){
  for(response in vars){
    plot_panel2(shock,response,country="jpn")
  }
}
par(mfrow = c(1,1))

par(mfrow = c(4,4), mar = c(3,4,3,1))
for(shock in vars){
  for(response in vars){
    plot_panel2(shock,response,country="nor")
  }
}
par(mfrow = c(1,1))








































