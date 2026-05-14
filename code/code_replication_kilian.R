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


#Data loading and basic cleaning-------------------------------------------------------------------------------

fredr_set_key(Sys.getenv("FRED_API_KEY"))


kilian = fredr(
  series_id = "IGREA", 
  observation_start = as.Date("1973-01-01"),
  observation_end   = as.Date("2026-02-02")
)

df_clean2 <- read_csv("data/Crude_prod_TS.csv")

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


price = read.csv("data/U.S._Crude_Oil_Imported_Acquisition_Cost_by_Refiners.csv") #Not from FRED because we extrapolate for 1973 

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
    growth_prod = (log(prod) - log(lag(prod)))*100,
        indexbis = index / 100) 

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
  mutate( real_price = Cost / cpi,
  log_real_price = log(real_price)) 


df_SVAR <- df_main |>
  select(date, log_real_price, growth_prod, index) %>%  #change 1
  filter(date < as.Date("2008-01-01")) |> 
  filter(!is.na(growth_prod))
#Visualize the series-------------------------------------------------------

ggplot(data = df_SVAR, aes(x = date))+
  geom_line(aes(y = log_real_price)) + 
  geom_line(aes(y = growth_prod))

p <- ggplot(df_SVAR, aes(x = date)) +
  geom_line(aes(y = log_real_price), color = "steelblue") +
  geom_line(aes(y = growth_prod), color = "firebrick")

print(p)

#Stationarity tests----------------------------------------------------

dfbis <- df_SVAR %>% select(!date)

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
print(ur_tbl)

#Due to non-stationarity, we have to difference the series.

df_stat <- df_main |> 
  mutate(
    growth_real_price = (log_real_price - lag(log_real_price))*100
  ) |>
  filter(!is.na(growth_real_price)) |>
  filter(date < as.Date("2008-01-01")) |> 
  select(growth_real_price, date, growth_prod, index)




df_stat_bis <- df_stat %>% select(!date)

ur_tbl_stat <- data.frame(                                  # Empty container
  variable   = colnames(df_stat_bis),
  pp_stat    = NA_real_, pp_pval   = NA_real_,
  dfgls_stat = NA_real_, dfgls_cv5 = NA_real_,
  kpss_stat  = NA_real_, kpss_pval = NA_real_
)
for (j in seq_along(colnames(df_stat_bis))) {               # Loop over 3 variables
  y <- df_stat_bis[, j]                                     # Pick the j-th series
  
  pp  <- pp.test(y, alternative = "stationary")        # Phillips-Perron
  ers <- ur.ers(y, type = "DF-GLS", model = "constant")# DF-GLS (ERS 1996)
  kp  <- kpss.test(y, null = "Level")                  # KPSS (null: stationary)
  
  ur_tbl_stat$pp_stat[j]    <- unname(pp$statistic)         # Store PP stat / p-value
  ur_tbl_stat$pp_pval[j]    <- pp$p.value
  ur_tbl_stat$dfgls_stat[j] <- as.numeric(ers@teststat)     # Uses 5% crit. val.
  ur_tbl_stat$dfgls_cv5[j]  <- ers@cval[, "5pct"]
  ur_tbl_stat$kpss_stat[j]  <- unname(kp$statistic)         # KPSS stat / p-value
  ur_tbl_stat$kpss_pval[j]  <- kp$p.value
}
print(ur_tbl_stat)

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

#We have to start with differencing

# 5. Lag order selection -------------------------------------------------------

lag_sel <- VARselect(df_stat_bis, lag.max = 24, type = "const")
print(lag_sel$selection)
print(lag_sel$criteria)

# 6. Estimate the VAR ----------------------------------------------------------

p_star <- as.integer(lag_sel$selection["HQ(n)"])
cat("VAR lag order selected by HQ: p =", p_star, "\n")
cat("  Other criteria: AIC =", lag_sel$selection["AIC(n)"],
    ", BIC =", lag_sel$selection["SC(n)"],
    ", FPE =", lag_sel$selection["FPE(n)"], "\n")

var_fit <- VAR(df_stat_bis, p = p_star, type = "const")  # Estimate VAR(p) by OLS
summary(var_fit)                                    # Coefficients + res. stats


# 7. Stability & residual diagnostics ------------------------------------------


# (a) Stability: all eigenvalues of the companion matrix should lie inside the unit circle.
roots_mod <- roots(var_fit, modulus = TRUE)
cat("Max modulus of companion roots:", round(max(roots_mod), 3), "\n")

# (b) Serial correlation: Portmanteau / Breusch-Godfrey (H0: no autocorr.)
sc_test <- serial.test(var_fit, lags.pt = 36, type = "PT.asymptotic")
print(sc_test)

# (c) Normality of residuals (Jarque-Bera multivariate, H0: normal errors)
nm_test <- normality.test(var_fit, multivariate.only = TRUE)
print(nm_test$jb.mul)

# 8. Impulse Response Functions ------------------------------------------------


cumulate_irf <- function(irf_obj, impulse_var, vars_to_cumsum) {
  for (v in vars_to_cumsum) {
    irf_obj$irf[[impulse_var]][, v]   <- cumsum(irf_obj$irf[[impulse_var]][, v])
    irf_obj$Lower[[impulse_var]][, v] <- cumsum(irf_obj$Lower[[impulse_var]][, v])
    irf_obj$Upper[[impulse_var]][, v] <- cumsum(irf_obj$Upper[[impulse_var]][, v])
  }
  return(irf_obj)
}

df_ordered = df_stat_bis[, c("growth_prod", "index", "growth_real_price")]
var_ordered <- VAR(df_ordered, p =p_star, type = "const")  # Re-estimate with ordered vars


#IRF positive shock in growth_prod-----------------------------------

irf_supply <- irf(var_ordered, impulse = "growth_prod", response = c("growth_prod", "index", "growth_real_price"),
n.ahead = 18, boot = TRUE, ci = 0.95)
irf_supply <- cumulate_irf(irf_supply, "growth_prod", c("growth_prod", "growth_real_price"))


plot(irf_supply, main = "Impulse Response to a shock in growth_prod", xlab = "Months", ylab = "Response")

# IRF for a negative shock in growth_prod (as in Kilian 2009)---------


# AD shock

var_ordered <- VAR(df_ordered, p = p_star, type = "const")  # Re-estimate with ordered vars

irf_ad <- irf(var_ordered, impulse = "index", response = c("growth_prod", "index", "growth_real_price"),
n.ahead = 18, boot = TRUE, ci = 0.95)
irf_ad <- cumulate_irf(irf_ad, "index", c("growth_prod", "growth_real_price"))

plot(irf_ad, main = "Impulse Response to a shock in index", xlab = "Months", ylab = "Response")


# OIl-specific demand shock

var_ordered <- VAR(df_ordered, p = p_star, type = "const")  # Re-estimate with ordered vars

irf_oild <- irf(var_ordered, impulse = "growth_real_price", response = c("growth_prod", "index", "growth_real_price"),
n.ahead = 18, boot = TRUE, ci = 0.95)
irf_oild <- cumulate_irf(irf_oild, "growth_real_price", c("growth_prod", "growth_real_price"))


plot(irf_oild, main = "Impulse Response to a shock in growth_real_price", xlab = "Months", ylab = "Response")



# 9. Forecast Error Variance Decomposition ------------------------------------

# Reduced-form residuals from the ordered VAR
u_hat <- residuals(var_ordered)
Sigma_u <- crossprod(u_hat) / nrow(u_hat)
A0_inv <- t(chol(Sigma_u))

eps_hat <- t(solve(A0_inv, t(u_hat)))
eps_hat <- as.data.frame(eps_hat)

colnames(eps_hat) <- c(
  "Oil supply shock",
  "Aggregate demand shock",
  "Oil-specific demand shock"
)

# Add dates from df_stat, because df_stat_bis has no date column
eps_hat$date <- tail(df_stat$date, nrow(eps_hat))
eps_hat$date <- as.Date(eps_hat$date)

# Annual averages --------------------------------------------------------

eps_hat_annual <- eps_hat %>%
  mutate(year = lubridate::year(date)) %>%
  group_by(year) %>%
  summarise(
    `Oil supply shock` = mean(`Oil supply shock`, na.rm = TRUE),
    `Aggregate demand shock` = mean(`Aggregate demand shock`, na.rm = TRUE),
    `Oil-specific demand shock` = mean(`Oil-specific demand shock`, na.rm = TRUE),
    .groups = "drop"
  )

# Plot annual structural shocks -----------------------------------------

par(mfrow = c(3, 1), mar = c(3, 4, 3, 1))

plot(eps_hat_annual$year, eps_hat_annual$`Oil supply shock`,
     type = "l", main = "Oil supply shock",
     xlab = "", ylab = "")

abline(h = 0, col = "gray")

plot(eps_hat_annual$year, eps_hat_annual$`Aggregate demand shock`,
     type = "l", main = "Aggregate demand shock",
     xlab = "", ylab = "")

abline(h = 0, col = "gray")

plot(eps_hat_annual$year, eps_hat_annual$`Oil-specific demand shock`,
     type = "l", main = "Oil-specific demand shock",
     xlab = "Year", ylab = "")

abline(h = 0, col = "gray")

par(mfrow = c(1, 1))

# 10. Granger causality (extension) -------------------------------------------



