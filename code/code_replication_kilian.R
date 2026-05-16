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


fredr_set_key(Sys.getenv("FRED_API_KEY"))


kilian = fredr(
  series_id = "IGREA", 
  observation_start = as.Date("1973-01-01"),
  observation_end   = as.Date("2026-02-02")
)

plot(kilian$date, kilian$value, type = "l")

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
  growth_prod = (log(prod) - log(lag(prod)))*100*12) 

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
  select(date, log_real_price, growth_prod, index) |>  
  filter(date < as.Date("2008-01-01")) |> 
  filter(!is.na(growth_prod))

df_SVAR <- df_SVAR |> mutate(index = index * (24.08 / sd(index, na.rm = TRUE)))
#Visualize the series-------------------------------------------------------

ggplot(data = df_SVAR, aes(x = date))+
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


# 5. Lag order selection -------------------------------------------------------

lag_sel <- VARselect(dfbis, lag.max = 24, type = "const")
print(lag_sel$selection)
print(lag_sel$criteria)

# 6. Estimate the VAR ----------------------------------------------------------

p_star <- as.integer(lag_sel$selection["HQ(n)"])
cat("VAR lag order selected by HQ: p =", p_star, "\n")
cat("  Other criteria: AIC =", lag_sel$selection["AIC(n)"],
    ", BIC =", lag_sel$selection["SC(n)"],
    ", FPE =", lag_sel$selection["FPE(n)"], "\n")

var_fit <- VAR(dfbis, p = 24, type = "const")  # Estimate VAR(p) by OLS
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

df_ordered = dfbis[, c("growth_prod", "index", "log_real_price")]
var_ordered <- VAR(df_ordered, p =24, type = "const")  # Re-estimate with ordered vars

H <- 18
vars <- c("growth_prod", "index", "log_real_price")
vars_ordered <- c("growth_prod", "index", "log_real_price")

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

# 1. Compute Standard IRFs (for Index and Price)
irf_68_std <- lapply(vars, \(s) irf(var_ordered, impulse = s, response = vars,
                                    n.ahead = H, boot = TRUE, ci = 0.68, runs = 100, cumulative = FALSE))
irf_95_std <- lapply(vars, \(s) irf(var_ordered, impulse = s, response = vars,
                                    n.ahead = H, boot = TRUE, ci = 0.95, runs = 100, cumulative = FALSE))

# 2. Compute Cumulative IRFs (for Production ONLY)
irf_68_cum <- lapply(vars, \(s) irf(var_ordered, impulse = s, response = vars,
                                    n.ahead = H, boot = TRUE, ci = 0.68, runs = 100, cumulative = TRUE))
irf_95_cum <- lapply(vars, \(s) irf(var_ordered, impulse = s, response = vars,
                                    n.ahead = H, boot = TRUE, ci = 0.95, runs = 100, cumulative = TRUE))

names(irf_68_std) <- vars; names(irf_95_std) <- vars
names(irf_68_cum) <- vars; names(irf_95_cum) <- vars

ylim_fixed <- list(
  growth_prod    = c(-25, 15),
  index          = c(-5,  10),
  log_real_price = c(-7,  12)
)

plot_panel <- function(shock, response) {
  
  # Select the correct IRF object based on the response variable
  if (response == "growth_prod") {
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
  
  # Supply shock in Kilian = negative production shock
  if (shock == "growth_prod") {
    y  <- -y
    old_l1 <- l1; l1 <- -u1; u1 <- -old_l1
    old_l2 <- l2; l2 <- -u2; u2 <- -old_l2
  }
  
  # Log price response -> percent response
  if (response == "log_real_price") {
    y  <- 100 * y
    l1 <- 100 * l1
    u1 <- 100 * u1
    l2 <- 100 * l2
    u2 <- 100 * u2
  }
  
  # Notice that the manual cumsum() block has been completely removed.
  
  x <- 0:H
  
  plot(x, y, type = "l", lwd = 2,
       main = unname(shock_names[shock]),
       ylab = unname(ylabs[response]),
       xlab =  if (shock == "log_real_price") "Months" else "",
       xaxs = "i",
       ylim = ylim_fixed[[response]])
  
  abline(h = 0, col = "gray")
  lines(x, l1, lty = 2)
  lines(x, u1, lty = 2)
  lines(x, l2, lty = 3)
  lines(x, u2, lty = 3)
}

png("irf_plots.png", width = 900, height = 700, res = 132) #to make it look as R studio output 
#We may use CairoPNG() instead

par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))

for (shock in vars) {
  for (response in vars) {
    plot_panel(shock, response)
  }
}

par(mfrow = c(1, 1))

dev.off()






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

eps_hat$date <- tail(df_SVAR$date, nrow(eps_hat))
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
