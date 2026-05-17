setwd("C:/Users/taher/OneDrive/Bureau/econometrics 3")
install.packages(c("fredr","vars","tseries","urca","ggplot2","dplyr","tidyr","lubridate","readr",
                   "tidyverse"))
library(fredr)       
library(vars)        
library(tseries)     
library(urca)        
library(ggplot2)     
library(dplyr)      
library(tidyr)       
library(lubridate)
library(readr)
library(tidyverse)


#Data loading and basic cleaning-------------------------------------------------------------------------------

fredr_set_key("ab808af2bcc04e9a32c684991350a0cb")


kilian = fredr(
  series_id = "IGREA", 
  observation_start = as.Date("1973-01-01"),
  observation_end   = as.Date("2026-02-02")
)

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

ggplot(data = df_main, aes(x = date))+
  geom_line(aes(y = real_price))


df_SVAR <- df_main |>
  select(date, log_real_price, growth_prod, index) %>%  #change 1
  filter(date < as.Date("2026-01-01")) |> 
  filter(!is.na(growth_prod))

df_7385 <- df_SVAR %>%
  filter(date >= as.Date("1973-01-01"),
         date <= as.Date("1985-12-01"))

df_8607 <- df_SVAR %>%
  filter(date >= as.Date("1986-01-01"),
         date <= as.Date("2007-12-01"))

df_0825 <- df_SVAR %>%
  filter(date >= as.Date("2008-01-01")) %>%
  na.omit()

#Visualize the series-------------------------------------------------------

ggplot(data = df_SVAR, aes(x = date))+
  geom_line(aes(y = log_real_price))

ggplot(data = df_7385, aes(x = date))+
  geom_line(aes(y = log_real_price))

ggplot(data = df_8607, aes(x = date))+
  geom_line(aes(y = log_real_price))

ggplot(data = df_0825, aes(x = date))+
  geom_line(aes(y = log_real_price))

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

# =========================
# 1973-1985
# =========================

dfbis_7385 <- df_7385 %>% select(!date)

ur_tbl_7385 <- data.frame(
  variable   = colnames(dfbis_7385),
  pp_stat    = NA_real_, pp_pval   = NA_real_,
  dfgls_stat = NA_real_, dfgls_cv5 = NA_real_,
  kpss_stat  = NA_real_, kpss_pval = NA_real_
)

for (j in seq_along(colnames(dfbis_7385))) {
  
  y <- dfbis_7385[, j]
  
  pp  <- pp.test(y, alternative = "stationary")
  ers <- ur.ers(y, type = "DF-GLS", model = "constant")
  kp  <- kpss.test(y, null = "Level")
  
  ur_tbl_7385$pp_stat[j]    <- unname(pp$statistic)
  ur_tbl_7385$pp_pval[j]    <- pp$p.value
  ur_tbl_7385$dfgls_stat[j] <- as.numeric(ers@teststat)
  ur_tbl_7385$dfgls_cv5[j]  <- ers@cval[, "5pct"]
  ur_tbl_7385$kpss_stat[j]  <- unname(kp$statistic)
  ur_tbl_7385$kpss_pval[j]  <- kp$p.value
}

print(ur_tbl_7385)


# =========================
# 1986-2007
# =========================

dfbis_8607 <- df_8607 %>% select(!date)

ur_tbl_8607 <- data.frame(
  variable   = colnames(dfbis_8607),
  pp_stat    = NA_real_, pp_pval   = NA_real_,
  dfgls_stat = NA_real_, dfgls_cv5 = NA_real_,
  kpss_stat  = NA_real_, kpss_pval = NA_real_
)

for (j in seq_along(colnames(dfbis_8607))) {
  
  y <- dfbis_8607[, j]
  
  pp  <- pp.test(y, alternative = "stationary")
  ers <- ur.ers(y, type = "DF-GLS", model = "constant")
  kp  <- kpss.test(y, null = "Level")
  
  ur_tbl_8607$pp_stat[j]    <- unname(pp$statistic)
  ur_tbl_8607$pp_pval[j]    <- pp$p.value
  ur_tbl_8607$dfgls_stat[j] <- as.numeric(ers@teststat)
  ur_tbl_8607$dfgls_cv5[j]  <- ers@cval[, "5pct"]
  ur_tbl_8607$kpss_stat[j]  <- unname(kp$statistic)
  ur_tbl_8607$kpss_pval[j]  <- kp$p.value
}

print(ur_tbl_8607)

# =========================
# 2008-2025
# =========================

dfbis_0825 <- df_0825 %>% select(!date)

ur_tbl_0825 <- data.frame(
  variable   = colnames(dfbis_0825),
  pp_stat    = NA_real_, pp_pval   = NA_real_,
  dfgls_stat = NA_real_, dfgls_cv5 = NA_real_,
  kpss_stat  = NA_real_, kpss_pval = NA_real_
)

for (j in seq_along(colnames(dfbis_0825))) {
  
y <- dfbis_0825[, j]
  
  pp  <- pp.test(y, alternative = "stationary")
  ers <- ur.ers(y, type = "DF-GLS", model = "constant")
  kp  <- kpss.test(y, null = "Level")
  
  ur_tbl_0825$pp_stat[j]    <- unname(pp$statistic)
  ur_tbl_0825$pp_pval[j]    <- pp$p.value
  ur_tbl_0825$dfgls_stat[j] <- as.numeric(ers@teststat)
  ur_tbl_0825$dfgls_cv5[j]  <- ers@cval[, "5pct"]
  ur_tbl_0825$kpss_stat[j]  <- unname(kp$statistic)
  ur_tbl_0825$kpss_pval[j]  <- kp$p.value
}

print(ur_tbl_0825)


#index et log pas stationnaire, mais growth oui 

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
is.na(dfbis)
lag_sel <- VARselect(dfbis, lag.max = 24, type = "const")
print(lag_sel$selection)
print(lag_sel$criteria)

lag_sel_7385 <- VARselect(dfbis_7385, lag.max = 24, type = "const")
print(lag_sel_7385$selection)
print(lag_sel_7385$criteria)

lag_sel_8607 <- VARselect(dfbis_8607, lag.max = 24, type = "const")
print(lag_sel_8607$selection)
print(lag_sel_8607$criteria)

lag_sel_0825 <- VARselect(dfbis_0825, lag.max = 24, type = "const")
print(lag_sel_0825$selection)
print(lag_sel_0825$criteria)

# 6. Estimate the VAR ----------------------------------------------------------

p_star <- as.integer(lag_sel$selection["HQ(n)"])
cat("VAR lag order selected by HQ: p =", p_star, "\n")
cat("  Other criteria: AIC =", lag_sel$selection["AIC(n)"],
    ", BIC =", lag_sel$selection["SC(n)"],
    ", FPE =", lag_sel$selection["FPE(n)"], "\n")

var_fit <- VAR(dfbis, p = 24, type = "const")  # Estimate VAR(p) by OLS
summary(var_fit)                                    # Coefficients + res. stats

p_star_7385  <- as.integer(lag_sel_7385$selection["HQ(n)"])
cat("VAR lag order selected by HQ: p =", p_star, "\n")
cat("  Other criteria: AIC =", lag_sel_7385$selection["AIC(n)"],
    ", BIC =", lag_sel_7385$selection["SC(n)"],
    ", FPE =", lag_sel_7385$selection["FPE(n)"], "\n")

var_fit_7385 <- VAR(dfbis_7385, p = 24, type = "const")  # Estimate VAR(p) by OLS
summary(var_fit_7385)                                    # Coefficients + res. stats

p_star_8607 <- as.integer(lag_sel_8607$selection["HQ(n)"])
cat("VAR lag order selected by HQ: p =", p_star, "\n")
cat("  Other criteria: AIC =", lag_sel_8607$selection["AIC(n)"],
    ", BIC =", lag_sel_8607$selection["SC(n)"],
    ", FPE =", lag_sel_8607$selection["FPE(n)"], "\n")

var_fit_8607 <- VAR(dfbis_8607, p = 24, type = "const")  # Estimate VAR(p) by OLS
summary(var_fit_8607)                                    # Coefficients + res. stats


p_star_0825 <- as.integer(lag_sel_0825$selection["HQ(n)"])
cat("VAR lag order selected by HQ: p =", p_star, "\n")
cat("  Other criteria: AIC =", lag_sel_0825$selection["AIC(n)"],
    ", BIC =", lag_sel_0825$selection["SC(n)"],
    ", FPE =", lag_sel_0825$selection["FPE(n)"], "\n")

var_fit_0825 <- VAR(dfbis_0825, p = 24, type = "const")  # Estimate VAR(p) by OLS
summary(var_fit_0825)                                    # Coefficients + res. stats


# 7. Stability & residual diagnostics ------------------------------------------


# (a) Stability: all eigenvalues of the companion matrix should lie inside the unit circle.
roots_mod <- roots(var_fit, modulus = TRUE)
cat("Max modulus of companion roots:", round(max(roots_mod), 3), "\n")

roots_mod_7385  <- roots(var_fit_7385, modulus = TRUE)
cat("Max modulus of companion roots:", round(max(roots_mod_7385), 3), "\n")

roots_mod_8607 <- roots(var_fit_8607, modulus = TRUE)
cat("Max modulus of companion roots:", round(max(roots_mod_8607), 3), "\n")

roots_mod_0825 <- roots(var_fit_0825, modulus = TRUE)
cat("Max modulus of companion roots:", round(max(roots_mod_0825), 3), "\n")

# (b) Serial correlation: Portmanteau / Breusch-Godfrey (H0: no autocorr.)
sc_test <- serial.test(var_fit, lags.pt = 36, type = "PT.asymptotic")
print(sc_test)

sc_test_7385   <- serial.test(var_fit_7385, lags.pt = 36, type = "PT.asymptotic")
print(sc_test_7385)

sc_test_8607 <- serial.test(var_fit_8607, lags.pt = 36, type = "PT.asymptotic")
print(sc_test_8607)

sc_test_0825  <- serial.test(var_fit_0825, lags.pt = 36, type = "PT.asymptotic")
print(sc_test_0825)

# (c) Normality of residuals (Jarque-Bera multivariate, H0: normal errors)
nm_test <- normality.test(var_fit, multivariate.only = TRUE)
print(nm_test$jb.mul)

nm_test_7385 <- normality.test(var_fit_7385, multivariate.only = TRUE)
print(nm_test_7385$jb.mul)

nm_test_8607 <- normality.test(var_fit_8607, multivariate.only = TRUE)
print(nm_test_8607$jb.mul)

nm_test_0825 <- normality.test(var_fit_0825, multivariate.only = TRUE)
print(nm_test_0825$jb.mul)

# 8. Impulse Response Functions ------------------------------------------------

df_ordered = dfbis[, c("growth_prod", "index", "log_real_price")]
var_ordered <- VAR(df_ordered, p =24, type = "const")  # Re-estimate with ordered vars

df_ordered_7385 = dfbis_7385[, c("growth_prod", "index", "log_real_price")]
var_ordered_7385 <- VAR(df_ordered_7385, p =24, type = "const")  # Re-estimate with ordered vars

df_ordered_8607 = dfbis_8607[, c("growth_prod", "index", "log_real_price")]
var_ordered_8607 <- VAR(df_ordered_8607, p =24, type = "const")  # Re-estimate with ordered vars

df_ordered_0825  = dfbis_0825[, c("growth_prod", "index", "log_real_price")]
var_ordered_0825 <- VAR(df_ordered_0825, p =24, type = "const")  # Re-estimate with ordered vars

# 3x3 IRFs ---------------------------------------------------

H <- 18
vars <- c("growth_prod", "index", "log_real_price")

shock_names <- c(
  growth_prod = "Oil supply shock",
  index = "Aggregate demand shock",
  log_real_price = "Oil-specific demand shock"
)

ylabs <- c(
  growth_prod = "Oil production",
  index = "Real activity",
  log_real_price = "Real price of oil"
)

# 1973–1985

# Standard IRFs
irf_68_std_7385 <- lapply(vars, \(s) irf(var_ordered_7385, impulse = s, response = vars,
n.ahead = H, boot = TRUE, ci = 0.68, runs = 100, cumulative = FALSE))

irf_95_std_7385 <- lapply(vars, \(s) irf(var_ordered_7385, impulse = s, response = vars,
n.ahead = H, boot = TRUE, ci = 0.95, runs = 100, cumulative = FALSE))

# Cumulative IRFs
irf_68_cum_7385 <- lapply(vars, \(s) irf(var_ordered_7385, impulse = s, response = vars,
n.ahead = H, boot = TRUE, ci = 0.68, runs = 100, cumulative = TRUE))

irf_95_cum_7385 <- lapply(vars, \(s) irf(var_ordered_7385, impulse = s, response = vars,
 n.ahead = H, boot = TRUE, ci = 0.95, runs = 100, cumulative = TRUE))

names(irf_68_std_7385) <- vars; names(irf_95_std_7385) <- vars
names(irf_68_cum_7385) <- vars; names(irf_95_cum_7385) <- vars

# 1986–2007
# Standard IRFs
irf_68_std_8607 <- lapply(vars, \(s) irf(var_ordered_8607, impulse = s, response = vars,
n.ahead = H, boot = TRUE, ci = 0.68, runs = 100, cumulative = FALSE))

irf_95_std_8607 <- lapply(vars, \(s) irf(var_ordered_8607, impulse = s, response = vars,
n.ahead = H, boot = TRUE, ci = 0.95, runs = 100, cumulative = FALSE))

# Cumulative IRFs
irf_68_cum_8607 <- lapply(vars, \(s) irf(var_ordered_8607, impulse = s, response = vars,
n.ahead = H, boot = TRUE, ci = 0.68, runs = 100, cumulative = TRUE))

irf_95_cum_8607 <- lapply(vars, \(s) irf(var_ordered_8607, impulse = s, response = vars,
n.ahead = H, boot = TRUE, ci = 0.95, runs = 100, cumulative = TRUE))

names(irf_68_std_8607) <- vars; names(irf_95_std_8607) <- vars
names(irf_68_cum_8607) <- vars; names(irf_95_cum_8607) <- vars

# 2008–2025

# Standard IRFs
irf_68_std_0825 <- lapply(vars, \(s) irf(var_ordered_0825, impulse = s, response = vars,n.ahead = H, 
boot = TRUE, ci = 0.68, runs = 100, cumulative = FALSE))

irf_95_std_0825 <- lapply(vars, \(s) irf(var_ordered_0825, impulse = s, response = vars,
n.ahead = H, boot = TRUE, ci = 0.95, runs = 100, cumulative = FALSE))

# Cumulative IRFs
irf_68_cum_0825 <- lapply(vars, \(s) irf(var_ordered_0825, impulse = s, response = vars,
n.ahead = H, boot = TRUE, ci = 0.68, runs = 100, cumulative = TRUE))

irf_95_cum_0825 <- lapply(vars, \(s) irf(var_ordered_0825, impulse = s, response = vars,
n.ahead = H, boot = TRUE, ci = 0.95, runs = 100, cumulative = TRUE))

names(irf_68_std_0825) <- vars; names(irf_95_std_0825) <- vars
names(irf_68_cum_0825) <- vars; names(irf_95_cum_0825) <- vars

ylim_fixed <- list(
  growth_prod    = c(-25, 15),
  index          = c(-5,  10),
  log_real_price = c(-7,  12)
)

plot_panel <- function(shock, response,period=c("7385","8607","0825")) {
  if(period =="7385"){
    irf_68<-irf_68_7385
    irf_95<-irf_95_7385}

  if(period =="8607"){
    irf_68<-irf_68_8607
    irf_95<-irf_95_8607}
  
  if(period =="0825"){
    irf_68<-irf_68_0825
    irf_95<-irf_95_0825}

  if (response == "growth_prod") {
    obj_68 <- irf_68_cum
    obj_95 <- irf_95_cum
  } else {
    obj_68 <- irf_68_std
    obj_95 <- irf_95_std
  }
  
  y  <- irf_68[[shock]]$irf[[shock]][, response]
  l1 <- irf_68[[shock]]$Lower[[shock]][, response]
  u1 <- irf_68[[shock]]$Upper[[shock]][, response]
  l2 <- irf_95[[shock]]$Lower[[shock]][, response]
  u2 <- irf_95[[shock]]$Upper[[shock]][, response]

  # Supply shock in Kilian = negative production shock
  if (shock == "growth_prod") {
    y  <- -y
    old_l1 <- l1; l1 <- -u1; u1 <- -old_l1
    old_l2 <- l2; l2 <- -u2; u2 <- -old_l2
  }
  
  if (response == "growth_prod") {
    y  <- cumsum(y)
    l1 <- cumsum(l1)
    u1 <- cumsum(u1)
    l2 <- cumsum(l2)
    u2 <- cumsum(u2)
  }
  

  # Log price response -> percent response
  if (response == "log_real_price") {
    y  <- 100 * y
    l1 <- 100 * l1
    u1 <- 100 * u1
    l2 <- 100 * l2
    u2 <- 100 * u2
  }
  
  x <- 0:H
  
  plot(x, y, type = "l", lwd = 2,
       main = unname(shock_names[shock]),
       ylab = unname(ylabs[response]),
       xlab = "",
       ylim = range(c(y, l1, u1, l2, u2), na.rm = TRUE))
  
  abline(h = 0, col = "gray")
  lines(x, l1, lty = 2)
  lines(x, u1, lty = 2)
  lines(x, l2, lty = 3)
  lines(x, u2, lty = 3)
}
#------------
par(mfrow = c(3, 3), mar = c(3, 4, 3, 1))

for (shock in vars) {
  for (response in vars) {
    plot_panel(shock, response, period = "7385")
  }
}

par(mfrow = c(1, 1))
#--------------------------
par(mfrow = c(3, 3), mar = c(3, 4, 3, 1))

for (shock in vars) {
  for (response in vars) {
    plot_panel(shock, response, period = "8607")
  }
}

par(mfrow = c(1, 1))

#-----------------------

par(mfrow = c(3, 3), mar = c(3, 4, 3, 1))

for (shock in vars) {
  for (response in vars) {
    plot_panel(shock, response, period = "0825")
  }
}

par(mfrow = c(1, 1))

# 9. Forecast Error Variance Decomposition ------------------------------------

# Reduced-form residuals from the ordered VAR
u_hat <- residuals(var_ordered)
Sigma_u <- crossprod(u_hat) / nrow(u_hat)
A0_inv <- t(chol(Sigma_u))

eps_hat <- t(solve(A0_inv, t(u_hat)))
eps_hat <- scale(eps_hat)
eps_hat <- as.data.frame(eps_hat)

colnames(eps_hat) <- c(
  "Oil supply shock",
  "Aggregate demand shock",
  "Oil-specific demand shock"
)

u_hat_7385 <- residuals(var_ordered_7385)
Sigma_u_7385 <- crossprod(u_hat_7385) / nrow(u_hat_7385)
A0_inv_7385 <- t(chol(Sigma_u_7385))

eps_hat_7385 <- t(solve(A0_inv_7385, t(u_hat_7385)))
eps_hat_7385 <- scale(eps_hat_7385)
eps_hat_7385 <- as.data.frame(eps_hat_7385)

colnames(eps_hat_7385) <- c(
  "Oil supply shock",
  "Aggregate demand shock",
  "Oil-specific demand shock"
)

u_hat_8607 <- residuals(var_ordered_8607)
Sigma_u_8607 <- crossprod(u_hat_8607) / nrow(u_hat_8607)
A0_inv_8607 <- t(chol(Sigma_u_8607))

eps_hat_8607 <- t(solve(A0_inv_8607, t(u_hat_8607)))
eps_hat_8607 <- scale(eps_hat_8607)
eps_hat_8607 <- as.data.frame(eps_hat_8607)

colnames(eps_hat_8607) <- c(
  "Oil supply shock",
  "Aggregate demand shock",
  "Oil-specific demand shock"
)

u_hat_0825 <- residuals(var_ordered_0825)
Sigma_u_0825 <- crossprod(u_hat_0825) / nrow(u_hat_0825)
A0_inv_0825 <- t(chol(Sigma_u_0825))

eps_hat_0825 <- t(solve(A0_inv_0825, t(u_hat_0825)))
eps_hat_0825 <- scale(eps_hat_0825)
eps_hat_0825 <- as.data.frame(eps_hat_0825)

colnames(eps_hat_0825) <- c(
  "Oil supply shock",
  "Aggregate demand shock",
  "Oil-specific demand shock"
)


# Add dates from df_stat, because df_stat_bis has no date column
eps_hat$date <- tail(df_stat$date, nrow(eps_hat))
eps_hat$date <- as.Date(eps_hat$date)

eps_hat_7385$date <- tail(df_7385$date, nrow(eps_hat_7385))
eps_hat_7385$date <- as.Date(eps_hat_7385$date)

eps_hat_8607$date <- tail(df_8607$date, nrow(eps_hat_8607))
eps_hat_8607$date <- as.Date(eps_hat_8607$date)

eps_hat_0825$date <- tail(df_0825$date, nrow(eps_hat_0825))
eps_hat_0825$date <- as.Date(eps_hat_0825$date)

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

eps_hat_annual_7385 <- eps_hat_7385 %>%
  mutate(year = lubridate::year(date)) %>%
  group_by(year) %>%
  summarise(
    `Oil supply shock` = mean(`Oil supply shock`, na.rm = TRUE),
    `Aggregate demand shock` = mean(`Aggregate demand shock`, na.rm = TRUE),
    `Oil-specific demand shock` = mean(`Oil-specific demand shock`, na.rm = TRUE),
    .groups = "drop"
  )

eps_hat_annual_8607 <- eps_hat_8607 %>%
  mutate(year = lubridate::year(date)) %>%
  group_by(year) %>%
  summarise(
    `Oil supply shock` = mean(`Oil supply shock`, na.rm = TRUE),
    `Aggregate demand shock` = mean(`Aggregate demand shock`, na.rm = TRUE),
    `Oil-specific demand shock` = mean(`Oil-specific demand shock`, na.rm = TRUE),
    .groups = "drop"
  )

eps_hat_annual_0825 <- eps_hat_0825 %>%
  mutate(year = lubridate::year(date)) %>%
  group_by(year) %>%
  summarise(
    `Oil supply shock` = mean(`Oil supply shock`, na.rm = TRUE),
    `Aggregate demand shock` = mean(`Aggregate demand shock`, na.rm = TRUE),
    `Oil-specific demand shock` = mean(`Oil-specific demand shock`, na.rm = TRUE),
    .groups = "drop"
  )
# Plot annual structural shocks -----------------------------------------

plot_structural_shocks<- function(period = c("7385", "8607", "0825")){

if (period == "7385") {data<-eps_hat_annual_7385 
  x_lim <- c(1975, max(data$year))}
if (period == "8607") {data<-eps_hat_annual_8607 
  x_lim <- c(1986, max(data$year))}
if (period == "0825") {data<-eps_hat_annual_0825 
  x_lim <- c(2008, max(data$year))}

png(filename = paste0("structural_shocks_", period, ".png"),width = 900, height = 700, res = 132)

par(mfrow = c(3, 1), mar = c(2.5, 4, 2.5, 1), oma = c(1, 0, 0, 0))

shocks <- c("Oil supply shock", "Aggregate demand shock", "Oil-specific demand shock")

y_lim <- c(-1, 1)

for (shock in shocks) {
  plot(data$year, data[[shock]],
        type = "n",
        main = paste0(shock, " (", period, ")"),
        xlab = "", ylab = "",
        xlim = x_lim,
        ylim = y_lim,
        xaxs = "i",
        yaxs = "i",
        las  = 1,
        tcl  = -0.3,
        mgp  = c(2, 0.5, 0))
    
    abline(h = 0, col = "gray", lwd = 0.8)
    lines(data$year, data[[shock]], lwd = 1)
  }
  
  par(mfrow = c(1, 1))
  dev.off()
}

plot_structural_shocks("7385")
plot_structural_shocks("8607")
plot_structural_shocks("0825")
