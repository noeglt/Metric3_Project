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
library(httr)
library(jsonlite)


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
  filter(date < as.Date("2008-01-01")) |> 
  filter(!is.na(growth_prod))
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

###



## FIRST ATTEMPT: point estimates and scale are much closer to Killian's results, yet, CIs too wide (specifically for supply shock)
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

# Compute IRFs once per shock
irf_68 <- lapply(vars, \(s) irf(var_ordered, impulse = s, response = vars,
                               n.ahead = H, boot = TRUE, ci = 0.68, runs = 100))
irf_95 <- lapply(vars, \(s) irf(var_ordered, impulse = s, response = vars,
                               n.ahead = H, boot = TRUE, ci = 0.95, runs = 100))

names(irf_68) <- vars
names(irf_95) <- vars

plot_panel <- function(shock, response) {
  
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

par(mfrow = c(3, 3), mar = c(3, 4, 3, 1))

for (shock in vars) {
  for (response in vars) {
    plot_panel(shock, response)
  }
}

par(mfrow = c(1, 1))


## SECOND ATTEMPS: confidence intervals are less wide, yet, scale seems off 

H <- 18
B <- 100
set.seed(123)

vars <- c("growth_prod", "index", "log_real_price")
K <- length(vars)
p <- var_ordered$p

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

transform_irf <- function(z, shock, response) {
  if (shock == "growth_prod") z <- -z
  if (response == "growth_prod") z <- cumsum(z)
  if (response == "log_real_price") z <- 100 * z
  z
}

point_psi <- Psi(var_ordered, nstep = H - 1)
nH <- dim(point_psi)[3]

boot_irfs <- array(
  NA_real_,
  dim = c(K, K, nH, B),
  dimnames = list(
    response = vars,
    shock = vars,
    horizon = 0:(nH - 1),
    draw = NULL
  )
)

Y <- as.matrix(var_ordered$y)
T_full <- nrow(Y)
u_hat <- residuals(var_ordered)
coef_mat <- Bcoef(var_ordered)

for (b in 1:B) {
  
  eta <- sample(c(-1, 1), size = nrow(u_hat), replace = TRUE)
  u_star <- u_hat * eta
  
  Y_star <- Y
  Y_star[1:p, ] <- Y[1:p, ]
  
  for (t in (p + 1):T_full) {
    x_lags <- unlist(lapply(1:p, \(lag) Y_star[t - lag, ]))
    x_t <- c(x_lags, 1)
    Y_star[t, ] <- as.numeric(coef_mat %*% x_t + u_star[t - p, ])
  }
  
  Y_star <- as.data.frame(Y_star)
  colnames(Y_star) <- vars
  
  var_star <- try(VAR(Y_star, p = p, type = "const"), silent = TRUE)
  
  if (!inherits(var_star, "try-error")) {
    boot_irfs[, , , b] <- Psi(var_star, nstep = H - 1)
  }
}

plot_panel <- function(shock, response) {
  
  s <- match(shock, vars)
  r <- match(response, vars)
  
  point <- transform_irf(point_psi[r, s, ], shock, response)
  
  draws <- apply(
    boot_irfs[r, s, , , drop = FALSE],
    4,
    function(z) transform_irf(as.numeric(z), shock, response)
  )
  
  low1 <- apply(draws, 1, quantile, probs = 0.16,  na.rm = TRUE)
  up1  <- apply(draws, 1, quantile, probs = 0.84,  na.rm = TRUE)
  low2 <- apply(draws, 1, quantile, probs = 0.025, na.rm = TRUE)
  up2  <- apply(draws, 1, quantile, probs = 0.975, na.rm = TRUE)
  
  x <- 0:(length(point) - 1)
  
  plot(x, point, type = "l", lwd = 2,
       main = unname(shock_names[shock]),
       ylab = unname(ylabs[response]),
       xlab = "",
       ylim = range(c(point, low1, up1, low2, up2), na.rm = TRUE))
  
  abline(h = 0, col = "gray")
  lines(x, low1, lty = 2)
  lines(x, up1,  lty = 2)
  lines(x, low2, lty = 3)
  lines(x, up2,  lty = 3)
}

par(mfrow = c(3, 3), mar = c(3, 4, 3, 1))

for (shock in vars) {
  for (response in vars) {
    plot_panel(shock, response)
  }
}

par(mfrow = c(1, 1))



### THIRD ATTEMPT with kilianr package => IRFs not much consistent.. 

library(kilianr)

df_ordered <- dfbis[, c("growth_prod", "index", "log_real_price")]
df_ordered <- as.data.frame(df_ordered)

# Keep your dataframe name, but use Kilian-style internal variable names
colnames(df_ordered) <- c("oilsupply", "aggdemand", "rpoil")

H <- 18
nrep <- 100
p <- 24

# Estimate ordered VAR
var_ordered <- vars::VAR(df_ordered, p = p, type = "const")

# Objects needed by kilianr::irfvar()
Ahat <- vars::Bcoef(var_ordered)[, 1:(ncol(df_ordered) * p)]
SIGMAhat <- crossprod(residuals(var_ordered)) / nrow(residuals(var_ordered))
B0inv <- t(chol(SIGMAhat))

# Build ols object in the format expected by kilianr::bootstrap_wild()
sol_update <- kilianr::olsvarc(
  y = as.matrix(df_ordered),
  p = p
)

Ahat <- sol_update$Ahat
B0inv <- t(chol(sol_update$SIGMAhat))

# Point IRFs
irf_update <- kilianr::irfvar(
  Ahat = Ahat,
  B0inv = B0inv,
  p = p,
  h = H,
  var_order = colnames(df_ordered),
  var_cumsum = 1,
  negative_shocks = 1
)
dat_irf <- irf_update$irf_tidy

# 95% CI = +/- 2 standard errors
dat_irf_ci95 <- kilianr::bootstrap_wild(
  olsobj = sol_update,
  irfobj = irf_update,
  nrep = nrep,
  standard_factor = 2.0,
  bootstrap_seed = 676,
  display_progress_bar = TRUE
)

# 68% CI = +/- 1 standard error
dat_irf_ci68 <- kilianr::bootstrap_wild(
  olsobj = sol_update,
  irfobj = irf_update,
  nrep = nrep,
  standard_factor = 1.0,
  bootstrap_seed = 676,
  display_progress_bar = TRUE
)

plot_one_irf <- function(y, title, ylab) {
  
  lo <- paste0(y, "_lo")
  hi <- paste0(y, "_hi")
  
  ylim <- range(
    dat_irf[[y]],
    dat_irf_ci68[[lo]], dat_irf_ci68[[hi]],
    dat_irf_ci95[[lo]], dat_irf_ci95[[hi]],
    na.rm = TRUE
  )
  
  plot(dat_irf$horizon, dat_irf[[y]],
       type = "l", lwd = 2,
       main = title, xlab = "", ylab = ylab,
       ylim = ylim)
  
  abline(h = 0, col = "gray")
  
  lines(dat_irf_ci68$horizon, dat_irf_ci68[[lo]], lty = 2)
  lines(dat_irf_ci68$horizon, dat_irf_ci68[[hi]], lty = 2)
  
  lines(dat_irf_ci95$horizon, dat_irf_ci95[[lo]], lty = 3)
  lines(dat_irf_ci95$horizon, dat_irf_ci95[[hi]], lty = 3)
}

par(mfrow = c(3, 3), mar = c(3, 4, 3, 1))

plot_one_irf("response_shock_oilsupply_oilsupply", "Oil supply shock", "Oil production")
plot_one_irf("response_shock_aggdemand_oilsupply", "Oil supply shock", "Real activity")
plot_one_irf("response_shock_rpoil_oilsupply", "Oil supply shock", "Real price of oil")

plot_one_irf("response_shock_oilsupply_aggdemand", "Aggregate demand shock", "Oil production")
plot_one_irf("response_shock_aggdemand_aggdemand", "Aggregate demand shock", "Real activity")
plot_one_irf("response_shock_rpoil_aggdemand", "Aggregate demand shock", "Real price of oil")

plot_one_irf("response_shock_oilsupply_rpoil", "Oil-specific demand shock", "Oil production")
plot_one_irf("response_shock_aggdemand_rpoil", "Oil-specific demand shock", "Real activity")
plot_one_irf("response_shock_rpoil_rpoil", "Oil-specific demand shock", "Real price of oil")

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

# 10. Granger causality (extension) -------------------------------------------
granger_test_prod <- causality(var_ordered, cause = "growth_prod")
print(granger_test_prod$Granger)

granger_test_index <- causality(var_ordered, cause = "index")
print(granger_test_index$Granger)

granger_test_price <- causality(var_ordered, cause = "growth_real_price")
print(granger_test$Granger)

