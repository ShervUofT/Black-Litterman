rm(list = ls())


# Black-Litterman vs Naive Mean-Variance
# 47 Fama-French Industry Portfolios - Daily
# Uses LN19 formula: mu_E = g V w_M
# Then mu_BL = [ (tau V)^-1 + P' D^-1 P ]^-1 [ (tau V)^-1 mu_E + P' D^-1 Q ]

# Load data

indus <- read.csv("/Users/ibro/Downloads/49_Industry_Portfolios_Daily.csv",
                  skip = 9, check.names = FALSE)
ff    <- read.csv("/Users/ibro/Downloads/F-F_Research_Data_Factors_daily.csv",
                  skip = 4, check.names = FALSE)

# Monthly file that contains Average Firm Size
sizefile <- readLines("/Users/ibro/Downloads/49_Industry_Portfolios.csv")

names(indus)[1] <- "date"
names(ff)[1]    <- "date"

indus <- indus[grepl("^[0-9]{8}$", indus$date), ]
ff    <- ff[grepl("^[0-9]{8}$", ff$date), ]

indus$date <- as.numeric(as.character(indus$date))
ff$date    <- as.numeric(as.character(ff$date))

indus[-1] <- lapply(indus[-1], function(x) as.numeric(as.character(x)))
ff[-1]    <- lapply(ff[-1], function(x) as.numeric(as.character(x)))

# Keep all industries in the 49-industry file
ind_names <- names(indus)[-1]
indus <- indus[, c("date", ind_names)]

# Replace missing codes with NA
for (j in 2:ncol(indus)) {
  indus[[j]][indus[[j]] %in% c(-99.99, -999)] <- NA
}
for (j in 2:ncol(ff)) {
  ff[[j]][ff[[j]] %in% c(-99.99, -999)] <- NA
}

# Keep only factors needed
ff <- ff[, c("date", "Mkt-RF", "RF")]

# -----------------------------
# Read Average Firm Size section
# -----------------------------

size_start <- grep("Average Firm Size", sizefile)

size_raw <- read.csv(
  text = paste(sizefile[(size_start + 1):length(sizefile)], collapse = "\n"),
  check.names = FALSE
)

names(size_raw)[1] <- "date"
size_raw <- size_raw[grepl("^[0-9]{6}$", size_raw$date), ]

size_raw$date <- as.numeric(as.character(size_raw$date))
size_raw[-1]  <- lapply(size_raw[-1], function(x) as.numeric(as.character(x)))

size_raw <- size_raw[, c("date", ind_names)]

for (j in 2:ncol(size_raw)) {
  size_raw[[j]][size_raw[[j]] %in% c(-99.99, -999)] <- NA
}

# Convert monthly size data into normalized benchmark weights
size_mat <- as.matrix(size_raw[, ind_names])

for (i in 1:nrow(size_mat)) {
  good <- !is.na(size_mat[i, ])
  size_mat[i, !good] <- 0
  s <- sum(size_mat[i, ])
  if (s > 0) {
    size_mat[i, ] <- size_mat[i, ] / s
  } else {
    size_mat[i, ] <- rep(1 / ncol(size_mat), ncol(size_mat))
  }
}

size_weights <- data.frame(date = size_raw$date, size_mat)
row.names(size_weights) <- NULL

# -----------------------------
# Merge daily returns and factors
# -----------------------------

allret <- merge(indus, ff, by = "date")

# Start from first date with complete data
ok <- complete.cases(allret)
first_complete <- which(ok)[1]
allret <- allret[first_complete:nrow(allret), ]

# Remove any remaining incomplete rows
allret <- allret[complete.cases(allret), ]
row.names(allret) <- NULL

# Convert percent to decimals
allret[, -1] <- allret[, -1] / 100

# -----------------------------
# Map each daily date to monthly size weights
# -----------------------------

allret$ym <- floor(allret$date / 100)
size_weights$ym <- size_weights$date

wM_daily <- merge(
  allret[, c("date", "ym")],
  size_weights[, c("ym", ind_names)],
  by = "ym",
  all.x = TRUE,
  sort = FALSE
)

# Keep daily dates aligned
wM_daily <- wM_daily[match(allret$date, wM_daily$date), ]
row.names(wM_daily) <- NULL

# If any size row is missing, carry forward last available weights
for (i in 1:nrow(wM_daily)) {
  if (any(is.na(wM_daily[i, ind_names]))) {
    if (i == 1) {
      wM_daily[i, ind_names] <- rep(1 / length(ind_names), length(ind_names))
    } else {
      wM_daily[i, ind_names] <- wM_daily[i - 1, ind_names]
    }
  }
}

# -----------------------------
# Basic setup
# -----------------------------

nasset <- length(ind_names)

nyear <- length(unique(floor(allret$date / 10000)))
dpy   <- nrow(allret) / nyear

ret      <- as.matrix(allret[, ind_names])
rf       <- allret$RF
xrm      <- allret$`Mkt-RF`
indxret  <- ret - rf
dates    <- allret$date

# -----------------------------
# Inputs
# -----------------------------

window      <- 252
view_window <- 63
tau         <- 1 / window
gamma       <- 3
ridge       <- 1e-6

# Relative-view matrix P
P <- matrix(0, nrow = 3, ncol = nasset)
colnames(P) <- ind_names
rownames(P) <- c("Softw-Telcm", "Chips-Util", "Hlth-Steel")

P[1, c("Softw", "Telcm")] <- c(1, -1)
P[2, c("Chips", "Util")]  <- c(1, -1)
P[3, c("Hlth", "Steel")]  <- c(1, -1)

# -----------------------------
# Long-only optimizer:
# maximize w'mu - 0.5*gamma*w'Sigma*w
# s.t. sum(w)=1, w>=0
# -----------------------------
library(quadprog)

mv_weights_longonly <- function(mu, Sigma, gamma = 3, ridge = 1e-6) {
  Sigma <- Sigma + diag(ridge, ncol(Sigma))
  
  Dmat <- gamma * Sigma
  dvec <- mu
  
  Amat <- cbind(rep(1, length(mu)), diag(length(mu)))
  bvec <- c(1, rep(0, length(mu)))
  
  sol <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1)
  as.vector(sol$solution)
}

# -----------------------------
# Storage
# -----------------------------

nroll <- nrow(indxret) - window

wBL   <- matrix(NA, nrow = nroll, ncol = nasset)
wNV   <- matrix(NA, nrow = nroll, ncol = nasset)

muBL_store <- matrix(NA, nrow = nroll, ncol = nasset)
muE_store  <- matrix(NA, nrow = nroll, ncol = nasset)
muN_store  <- matrix(NA, nrow = nroll, ncol = nasset)

retBL_excess <- rep(NA, nroll)
retNV_excess <- rep(NA, nroll)
retBL_raw    <- rep(NA, nroll)
retNV_raw    <- rep(NA, nroll)

backtest_dates <- dates[(window + 1):length(dates)]

# -----------------------------
# Rolling BL and Naive MV
# -----------------------------

for (t in window:(nrow(indxret) - 1)) {
  
  k <- t - window + 1
  
  Rwin  <- indxret[(t - window + 1):t, , drop = FALSE]
  xrm_w <- xrm[(t - window + 1):t]
  
  # Annualized covariance
  V <- cov(Rwin) * dpy
  V <- V + diag(ridge, nasset)
  
  # Naive sample means
  muN <- colMeans(Rwin) * dpy
  
  # Risk aversion
  g <- mean(xrm_w) / var(xrm_w)
  
  # Monthly size-based benchmark weights for day t
  wM_t <- as.numeric(wM_daily[t, ind_names])
  wM_t <- wM_t / sum(wM_t)
  
  # Equilibrium returns from LN19
  muE <- as.vector(g * V %*% wM_t)
  
  # Relative views Q from recent spreads
  recent_mu <- colMeans(Rwin[(window - view_window + 1):window, , drop = FALSE]) * dpy
  
  Q <- matrix(c(
    recent_mu["Softw"] - recent_mu["Telcm"],
    recent_mu["Chips"] - recent_mu["Util"],
    recent_mu["Hlth"]  - recent_mu["Steel"]
  ), ncol = 1)
  
  # Diagonal uncertainty matrix D
  D <- diag(diag(P %*% V %*% t(P)))
  D <- D + diag(ridge, nrow(D))
  
  # BL posterior
  A_BL <- solve(tau * V)
  B_BL <- t(P) %*% solve(D) %*% P
  rhs  <- A_BL %*% muE + t(P) %*% solve(D) %*% Q
  muBL <- as.vector(solve(A_BL + B_BL, rhs))
  
  # Long-only portfolio weights
  w_bl <- mv_weights_longonly(muBL, V, gamma = gamma, ridge = ridge)
  w_nv <- mv_weights_longonly(muN,  V, gamma = gamma, ridge = ridge)
  
  wBL[k, ] <- w_bl
  wNV[k, ] <- w_nv
  
  muBL_store[k, ] <- muBL
  muE_store[k, ]  <- muE
  muN_store[k, ]  <- muN
  
  # Next-day returns
  next_excess <- indxret[t + 1, ]
  next_raw    <- ret[t + 1, ]
  
  retBL_excess[k] <- sum(w_bl * next_excess)
  retNV_excess[k] <- sum(w_nv * next_excess)
  
  retBL_raw[k] <- sum(w_bl * next_raw)
  retNV_raw[k] <- sum(w_nv * next_raw)
}

colnames(wBL) <- ind_names
colnames(wNV) <- ind_names
colnames(muBL_store) <- ind_names
colnames(muE_store)  <- ind_names
colnames(muN_store)  <- ind_names

# -----------------------------
# Summary statistics
# -----------------------------

ann_mean_bl <- mean(retBL_raw, na.rm = TRUE) * dpy
ann_vol_bl  <- sd(retBL_raw, na.rm = TRUE) * sqrt(dpy)
sharpe_bl   <- ann_mean_bl / ann_vol_bl

ann_mean_nv <- mean(retNV_raw, na.rm = TRUE) * dpy
ann_vol_nv  <- sd(retNV_raw, na.rm = TRUE) * sqrt(dpy)
sharpe_nv   <- ann_mean_nv / ann_vol_nv

turnover_bl <- mean(rowSums(abs(wBL[-1, ] - wBL[-nrow(wBL), ])), na.rm = TRUE)
turnover_nv <- mean(rowSums(abs(wNV[-1, ] - wNV[-nrow(wNV), ])), na.rm = TRUE)

conc_bl <- mean(rowSums(wBL^2), na.rm = TRUE)
conc_nv <- mean(rowSums(wNV^2), na.rm = TRUE)

stats <- rbind(
  BL = c(ann_return = ann_mean_bl, ann_vol = ann_vol_bl, sharpe = sharpe_bl,
         turnover = turnover_bl, concentration = conc_bl),
  Naive_MV = c(ann_return = ann_mean_nv, ann_vol = ann_vol_nv, sharpe = sharpe_nv,
               turnover = turnover_nv, concentration = conc_nv)
)

print(round(stats, 4))

# -----------------------------
# CAPM attribution
# -----------------------------

capm_data <- data.frame(
  date = backtest_dates,
  BL_excess = retBL_excess,
  NV_excess = retNV_excess,
  MKT = xrm[(window + 1):length(xrm)]
)

fit_bl <- lm(BL_excess ~ MKT, data = capm_data)
fit_nv <- lm(NV_excess ~ MKT, data = capm_data)

summary(fit_bl)
summary(fit_nv)

# -----------------------------
# Cumulative return plots
# -----------------------------

cum_bl <- cumprod(1 + retBL_raw)
cum_nv <- cumprod(1 + retNV_raw)

plot(backtest_dates, cum_bl,
     type = "l", lwd = 2, col = "red",
     xlab = "Date", ylab = "Cumulative Wealth",
     main = "Cumulative Return: BL vs Naive MV")
lines(backtest_dates, cum_nv, lwd = 2, lty = 2, col = "purple")
legend("topleft",
       legend = c("Black-Litterman", "Naive MV"),
       col = c("red", "purple"),
       lwd = 2, lty = c(1, 2), bty = "n")

# Optional log-scale version
plot(backtest_dates, log(cum_bl),
     type = "l", lwd = 2, col = "red",
     xlab = "Date", ylab = "Log Cumulative Wealth",
     main = "Log Cumulative Return: BL vs Naive MV")
lines(backtest_dates, log(cum_nv), lwd = 2, lty = 2, col = "purple")
legend("topleft",
       legend = c("Black-Litterman", "Naive MV"),
       col = c("red", "purple"),
       lwd = 2, lty = c(1, 2), bty = "n")


# -----------------------------
# Summary of average weights
# -----------------------------

avg_w_bl <- colMeans(wBL, na.rm = TRUE)
avg_w_nv <- colMeans(wNV, na.rm = TRUE)

avg_weights <- data.frame(
  Industry = ind_names,
  BL_AvgWeight = avg_w_bl,
  NaiveMV_AvgWeight = avg_w_nv
)

# Full table
print(avg_weights)

# Top 10 largest average weights in each strategy
top_bl <- avg_weights[order(-avg_weights$BL_AvgWeight), ]
top_nv <- avg_weights[order(-avg_weights$NaiveMV_AvgWeight), ]

cat("\nTop 10 average weights: Black-Litterman\n")
print(head(top_bl, 10))

cat("\nTop 10 average weights: Naive MV\n")
print(head(top_nv, 10))

# Concentration check through max average weight
cat("\nLargest average weight in BL:", round(max(avg_w_bl), 4), "\n")
cat("Largest average weight in Naive MV:", round(max(avg_w_nv), 4), "\n")



# ============================================================
# Professor Feedback Additions
# Views Only, Market Benchmark, Transaction Costs,
# Rolling 5-Year Sharpe, and Weight Stability
# ============================================================

# Set manually depending on the file:
# For BL_TRAILING.R:
view_type <- "trailing"

# For BL_MOMENTUM.R, use:
# view_type <- "momentum"

# For BL_MEANREVERSION.R, use:
# view_type <- "meanrev"


# -----------------------------
# Graph export folder
# -----------------------------

graph_dir <- "/Users/ibro/Downloads/GRAPHS"

if (!dir.exists(graph_dir)) {
  dir.create(graph_dir, recursive = TRUE)
}


# -----------------------------
# Make sure original BL/MV objects exist
# -----------------------------

if (!exists("wBL") | !exists("wNV")) {
  stop("wBL or wNV not found. Paste this section AFTER the original rolling backtest.")
}

if (!exists("retBL_raw")) {
  if (exists("retBL_excess")) {
    retBL_raw <- retBL_excess + rf[(window + 1):length(rf)]
  } else {
    stop("Neither retBL_raw nor retBL_excess exists.")
  }
}

if (!exists("retNV_raw")) {
  if (exists("retNV_excess")) {
    retNV_raw <- retNV_excess + rf[(window + 1):length(rf)]
  } else {
    stop("Neither retNV_raw nor retNV_excess exists.")
  }
}

if (!exists("backtest_dates")) {
  backtest_dates <- dates[(window + 1):length(dates)]
}

nroll <- length(retBL_raw)


# -----------------------------
# Storage for additional strategies
# -----------------------------

wVIEW <- matrix(NA, nrow = nroll, ncol = nasset)
wMKT  <- matrix(NA, nrow = nroll, ncol = nasset)

retVIEW_raw    <- rep(NA, nroll)
retVIEW_excess <- rep(NA, nroll)

retMKT_raw    <- rep(NA, nroll)
retMKT_excess <- rep(NA, nroll)


# -----------------------------
# Build Views Only and Market Size strategies
# -----------------------------

for (t in window:(nrow(indxret) - 1)) {
  
  k <- t - window + 1
  
  Rwin  <- indxret[(t - window + 1):t, , drop = FALSE]
  xrm_w <- xrm[(t - window + 1):t]
  
  # Annualized covariance matrix
  V <- cov(Rwin) * dpy
  V <- V + diag(ridge, nasset)
  
  # Sample mean expected returns
  muN <- colMeans(Rwin) * dpy
  
  # CAPM fallback for assets without direct views
  beta_i <- rep(NA, nasset)
  
  for (j in 1:nasset) {
    beta_i[j] <- cov(Rwin[, j], xrm_w) / var(xrm_w)
  }
  
  muCAPM <- beta_i * mean(xrm_w) * dpy
  names(muCAPM) <- ind_names
  
  # Recent returns used to form views
  recent_mu <- colMeans(Rwin[(window - view_window + 1):window, , drop = FALSE]) * dpy
  
  if (view_type == "trailing") {
    
    Q <- matrix(c(
      recent_mu["Softw"] - recent_mu["Telcm"],
      recent_mu["Chips"] - recent_mu["Util"],
      recent_mu["Hlth"]  - recent_mu["Steel"]
    ), ncol = 1)
    
  } else if (view_type == "momentum") {
    
    spread1 <- recent_mu["Softw"] - recent_mu["Telcm"]
    spread2 <- recent_mu["Chips"] - recent_mu["Util"]
    spread3 <- recent_mu["Hlth"]  - recent_mu["Steel"]
    
    view_size <- 0.03
    
    Q <- matrix(c(
      ifelse(spread1 > 0, view_size, -view_size),
      ifelse(spread2 > 0, view_size, -view_size),
      ifelse(spread3 > 0, view_size, -view_size)
    ), ncol = 1)
    
  } else if (view_type == "meanrev") {
    
    spread_hist1 <- Rwin[, "Softw"] - Rwin[, "Telcm"]
    spread_hist2 <- Rwin[, "Chips"] - Rwin[, "Util"]
    spread_hist3 <- Rwin[, "Hlth"]  - Rwin[, "Steel"]
    
    z1 <- mean(tail(spread_hist1, view_window)) / sd(spread_hist1)
    z2 <- mean(tail(spread_hist2, view_window)) / sd(spread_hist2)
    z3 <- mean(tail(spread_hist3, view_window)) / sd(spread_hist3)
    
    scale_view <- 0.02
    cap_view   <- 0.05
    
    q1 <- max(min(-z1 * scale_view, cap_view), -cap_view)
    q2 <- max(min(-z2 * scale_view, cap_view), -cap_view)
    q3 <- max(min(-z3 * scale_view, cap_view), -cap_view)
    
    Q <- matrix(c(q1, q2, q3), ncol = 1)
  }
  
  # Views Only expected return vector:
  # - assets with views use sample means plus view tilts
  # - assets without views use CAPM expected returns
  muVIEW <- muCAPM
  
  viewed_assets <- colSums(abs(P)) > 0
  muVIEW[viewed_assets] <- muN[viewed_assets]
  
  muVIEW <- muVIEW + as.vector(t(P) %*% Q)
  
  # Market-size benchmark weights
  wM_t <- as.numeric(wM_daily[t, ind_names])
  wM_t <- wM_t / sum(wM_t)
  
  # Optimize views-only portfolio
  w_view <- mv_weights_longonly(muVIEW, V, gamma = gamma, ridge = ridge)
  w_mkt  <- wM_t
  
  wVIEW[k, ] <- w_view
  wMKT[k, ]  <- w_mkt
  
  # Next-day returns
  next_excess <- indxret[t + 1, ]
  next_raw    <- ret[t + 1, ]
  
  retVIEW_excess[k] <- sum(w_view * next_excess)
  retMKT_excess[k]  <- sum(w_mkt  * next_excess)
  
  retVIEW_raw[k] <- sum(w_view * next_raw)
  retMKT_raw[k]  <- sum(w_mkt  * next_raw)
}

colnames(wVIEW) <- ind_names
colnames(wMKT)  <- ind_names


# -----------------------------
# Helper functions
# -----------------------------

perf_stats <- function(r, w) {
  
  ann_mean <- mean(r, na.rm = TRUE) * dpy
  ann_vol  <- sd(r, na.rm = TRUE) * sqrt(dpy)
  sharpe   <- ann_mean / ann_vol
  
  turnover <- mean(rowSums(abs(w[-1, ] - w[-nrow(w), ])), na.rm = TRUE)
  conc     <- mean(rowSums(w^2), na.rm = TRUE)
  maxw     <- mean(apply(w, 1, max), na.rm = TRUE)
  
  c(ann_return = ann_mean,
    ann_vol = ann_vol,
    sharpe = sharpe,
    turnover = turnover,
    concentration = conc,
    avg_max_weight = maxw)
}

turnover_series <- function(w) {
  c(0, rowSums(abs(w[-1, ] - w[-nrow(w), ])))
}

rolling_sharpe <- function(r, roll_window) {
  
  out <- rep(NA, length(r))
  
  for (i in roll_window:length(r)) {
    rr <- r[(i - roll_window + 1):i]
    out[i] <- mean(rr, na.rm = TRUE) / sd(rr, na.rm = TRUE) * sqrt(dpy)
  }
  
  out
}

roll_mean <- function(x, k) {
  
  out <- rep(NA, length(x))
  
  for (i in k:length(x)) {
    out[i] <- mean(x[(i - k + 1):i], na.rm = TRUE)
  }
  
  out
}


# -----------------------------
# Main performance table
# -----------------------------

feedback_stats <- rbind(
  Black_Litterman = perf_stats(retBL_raw, wBL),
  Naive_MV        = perf_stats(retNV_raw, wNV),
  Views_Only      = perf_stats(retVIEW_raw, wVIEW),
  Market_Size     = perf_stats(retMKT_raw, wMKT)
)

cat("\n==============================\n")
cat("PERFORMANCE COMPARISON\n")
cat("==============================\n")
print(round(feedback_stats, 4))


# -----------------------------
# Transaction costs
# -----------------------------

tcost <- 0.001   # 10 basis points per unit of turnover

retBL_net   <- retBL_raw   - tcost * turnover_series(wBL)
retNV_net   <- retNV_raw   - tcost * turnover_series(wNV)
retVIEW_net <- retVIEW_raw - tcost * turnover_series(wVIEW)
retMKT_net  <- retMKT_raw  - tcost * turnover_series(wMKT)

net_stats <- rbind(
  Black_Litterman = perf_stats(retBL_net, wBL),
  Naive_MV        = perf_stats(retNV_net, wNV),
  Views_Only      = perf_stats(retVIEW_net, wVIEW),
  Market_Size     = perf_stats(retMKT_net, wMKT)
)

cat("\n==============================\n")
cat("PERFORMANCE AFTER TRANSACTION COSTS\n")
cat("==============================\n")
print(round(net_stats, 4))


# -----------------------------
# CAPM attribution
# -----------------------------

capm_feedback <- data.frame(
  BL_excess       = retBL_raw - rf[(window + 1):length(rf)],
  NV_excess       = retNV_raw - rf[(window + 1):length(rf)],
  VIEW_excess     = retVIEW_raw - rf[(window + 1):length(rf)],
  MKT_Port_excess = retMKT_raw - rf[(window + 1):length(rf)],
  MKT             = xrm[(window + 1):length(xrm)]
)

fit_bl_feedback   <- lm(BL_excess ~ MKT, data = capm_feedback)
fit_nv_feedback   <- lm(NV_excess ~ MKT, data = capm_feedback)
fit_view_feedback <- lm(VIEW_excess ~ MKT, data = capm_feedback)
fit_mkt_feedback  <- lm(MKT_Port_excess ~ MKT, data = capm_feedback)

capm_table <- rbind(
  Black_Litterman = coef(fit_bl_feedback),
  Naive_MV        = coef(fit_nv_feedback),
  Views_Only      = coef(fit_view_feedback),
  Market_Size     = coef(fit_mkt_feedback)
)

colnames(capm_table) <- c("alpha_daily", "beta")

cat("\n==============================\n")
cat("CAPM ATTRIBUTION\n")
cat("==============================\n")
print(round(capm_table, 6))


# -----------------------------
# Rolling 5-year Sharpe ratios by decade
# Export each decade to PNG
# -----------------------------

roll_window <- round(5 * dpy)

roll_bl   <- rolling_sharpe(retBL_raw, roll_window)
roll_nv   <- rolling_sharpe(retNV_raw, roll_window)
roll_view <- rolling_sharpe(retVIEW_raw, roll_window)
roll_mkt  <- rolling_sharpe(retMKT_raw, roll_window)

plot_dates <- as.Date(as.character(backtest_dates), format = "%Y%m%d")
plot_years <- as.numeric(format(plot_dates, "%Y"))
plot_decades <- floor(plot_years / 10) * 10
decades <- sort(unique(plot_decades))

for (dec in decades) {
  
  idx <- which(plot_decades == dec)
  
  # Skip early decade if not enough rolling observations
  if (sum(!is.na(roll_bl[idx])) < 10) {
    next
  }
  
  y_all <- c(roll_bl[idx], roll_nv[idx], roll_view[idx], roll_mkt[idx])
  y_all <- y_all[!is.na(y_all)]
  
  png(
    filename = file.path(graph_dir, paste0("rolling_sharpe_", view_type, "_", dec, "s.png")),
    width = 1200,
    height = 800
  )
  
  plot(plot_dates[idx], roll_bl[idx],
       type = "l", lwd = 2, col = "red",
       ylim = range(y_all, na.rm = TRUE),
       xlab = "Date", ylab = "Rolling 5-Year Sharpe",
       main = paste0("Rolling 5-Year Sharpe Ratios: ", dec, "s"))
  
  lines(plot_dates[idx], roll_nv[idx],   lwd = 2, lty = 2, col = "purple")
  lines(plot_dates[idx], roll_view[idx], lwd = 2, lty = 3, col = "blue")
  lines(plot_dates[idx], roll_mkt[idx],  lwd = 2, lty = 4, col = "darkgreen")
  
  legend("topright",
         legend = c("Black-Litterman", "Naive MV", "Views Only", "Market Size"),
         col = c("red", "purple", "blue", "darkgreen"),
         lwd = 2,
         lty = c(1, 2, 3, 4),
         bty = "n")
  
  dev.off()
}


# -----------------------------
# Rolling 5-year return by decade
# Export each decade to PNG
# -----------------------------

rolling_return <- function(r, roll_window) {
  
  out <- rep(NA, length(r))
  
  for (i in roll_window:length(r)) {
    rr <- r[(i - roll_window + 1):i]
    out[i] <- prod(1 + rr, na.rm = TRUE) - 1
  }
  
  out
}

roll_ret_bl   <- rolling_return(retBL_raw, roll_window)
roll_ret_nv   <- rolling_return(retNV_raw, roll_window)
roll_ret_view <- rolling_return(retVIEW_raw, roll_window)
roll_ret_mkt  <- rolling_return(retMKT_raw, roll_window)

for (dec in decades) {
  
  idx <- which(plot_decades == dec)
  
  if (sum(!is.na(roll_ret_bl[idx])) < 10) {
    next
  }
  
  y_all <- c(roll_ret_bl[idx], roll_ret_nv[idx], roll_ret_view[idx], roll_ret_mkt[idx])
  y_all <- y_all[!is.na(y_all)]
  
  png(
    filename = file.path(graph_dir, paste0("rolling_5yr_return_", view_type, "_", dec, "s.png")),
    width = 1200,
    height = 800
  )
  
  plot(plot_dates[idx], roll_ret_bl[idx],
       type = "l", lwd = 2, col = "red",
       ylim = range(y_all, na.rm = TRUE),
       xlab = "Date", ylab = "Rolling 5-Year Return",
       main = paste0("Rolling 5-Year Return: ", dec, "s"))
  
  lines(plot_dates[idx], roll_ret_nv[idx],   lwd = 2, lty = 2, col = "purple")
  lines(plot_dates[idx], roll_ret_view[idx], lwd = 2, lty = 3, col = "blue")
  lines(plot_dates[idx], roll_ret_mkt[idx],  lwd = 2, lty = 4, col = "darkgreen")
  
  legend("topright",
         legend = c("Black-Litterman", "Naive MV", "Views Only", "Market Size"),
         col = c("red", "purple", "blue", "darkgreen"),
         lwd = 2,
         lty = c(1, 2, 3, 4),
         bty = "n")
  
  dev.off()
}


# -----------------------------
# Rolling 5-year volatility by decade
# Export each decade to PNG
# -----------------------------

rolling_vol <- function(r, roll_window) {
  
  out <- rep(NA, length(r))
  
  for (i in roll_window:length(r)) {
    rr <- r[(i - roll_window + 1):i]
    out[i] <- sd(rr, na.rm = TRUE) * sqrt(dpy)
  }
  
  out
}

roll_vol_bl   <- rolling_vol(retBL_raw, roll_window)
roll_vol_nv   <- rolling_vol(retNV_raw, roll_window)
roll_vol_view <- rolling_vol(retVIEW_raw, roll_window)
roll_vol_mkt  <- rolling_vol(retMKT_raw, roll_window)

for (dec in decades) {
  
  idx <- which(plot_decades == dec)
  
  if (sum(!is.na(roll_vol_bl[idx])) < 10) {
    next
  }
  
  y_all <- c(roll_vol_bl[idx], roll_vol_nv[idx], roll_vol_view[idx], roll_vol_mkt[idx])
  y_all <- y_all[!is.na(y_all)]
  
  png(
    filename = file.path(graph_dir, paste0("rolling_5yr_vol_", view_type, "_", dec, "s.png")),
    width = 1200,
    height = 800
  )
  
  plot(plot_dates[idx], roll_vol_bl[idx],
       type = "l", lwd = 2, col = "red",
       ylim = range(y_all, na.rm = TRUE),
       xlab = "Date", ylab = "Rolling 5-Year Annualized Volatility",
       main = paste0("Rolling 5-Year Volatility: ", dec, "s"))
  
  lines(plot_dates[idx], roll_vol_nv[idx],   lwd = 2, lty = 2, col = "purple")
  lines(plot_dates[idx], roll_vol_view[idx], lwd = 2, lty = 3, col = "blue")
  lines(plot_dates[idx], roll_vol_mkt[idx],  lwd = 2, lty = 4, col = "darkgreen")
  
  legend("topright",
         legend = c("Black-Litterman", "Naive MV", "Views Only", "Market Size"),
         col = c("red", "purple", "blue", "darkgreen"),
         lwd = 2,
         lty = c(1, 2, 3, 4),
         bty = "n")
  
  dev.off()
}


# -----------------------------
# Log cumulative wealth plot
# Export to PNG
# -----------------------------

cum_bl   <- cumprod(1 + retBL_raw)
cum_nv   <- cumprod(1 + retNV_raw)
cum_view <- cumprod(1 + retVIEW_raw)
cum_mkt  <- cumprod(1 + retMKT_raw)

png(
  filename = file.path(graph_dir, paste0("log_cumulative_wealth_", view_type, ".png")),
  width = 1200,
  height = 800
)

plot(plot_dates, log(cum_bl),
     type = "l", lwd = 2, col = "red",
     xlab = "Date", ylab = "Log Cumulative Wealth",
     main = "Log Cumulative Wealth Comparison")

lines(plot_dates, log(cum_nv),   lwd = 2, lty = 2, col = "purple")
lines(plot_dates, log(cum_view), lwd = 2, lty = 3, col = "blue")
lines(plot_dates, log(cum_mkt),  lwd = 2, lty = 4, col = "darkgreen")

legend("topleft",
       legend = c("Black-Litterman", "Naive MV", "Views Only", "Market Size"),
       col = c("red", "purple", "blue", "darkgreen"),
       lwd = 2,
       lty = c(1, 2, 3, 4),
       bty = "n")

dev.off()


# -----------------------------
# Weight stability plot by decade
# All strategies on one graph
# Export each decade to PNG
# -----------------------------

max_w_bl   <- apply(wBL, 1, max, na.rm = TRUE)
max_w_nv   <- apply(wNV, 1, max, na.rm = TRUE)
max_w_view <- apply(wVIEW, 1, max, na.rm = TRUE)
max_w_mkt  <- apply(wMKT, 1, max, na.rm = TRUE)

smooth_window <- round(dpy)

max_w_bl_s   <- roll_mean(max_w_bl, smooth_window)
max_w_nv_s   <- roll_mean(max_w_nv, smooth_window)
max_w_view_s <- roll_mean(max_w_view, smooth_window)
max_w_mkt_s  <- roll_mean(max_w_mkt, smooth_window)

for (dec in decades) {
  
  idx <- which(plot_decades == dec)
  
  if (sum(!is.na(max_w_bl_s[idx])) < 10) {
    next
  }
  
  y_all <- c(max_w_bl_s[idx], max_w_nv_s[idx],
             max_w_view_s[idx], max_w_mkt_s[idx])
  y_all <- y_all[!is.na(y_all)]
  
  png(
    filename = file.path(graph_dir, paste0("weight_stability_all_", view_type, "_", dec, "s.png")),
    width = 1200,
    height = 800
  )
  
  plot(plot_dates[idx], max_w_bl_s[idx],
       type = "l", lwd = 2, col = "red",
       ylim = range(y_all, na.rm = TRUE),
       xlab = "Date", ylab = "1-Year Avg Max Weight",
       main = paste0("Weight Stability: All Strategies, ", dec, "s"))
  
  lines(plot_dates[idx], max_w_nv_s[idx],
        lwd = 2, lty = 2, col = "purple")
  
  lines(plot_dates[idx], max_w_view_s[idx],
        lwd = 2, lty = 3, col = "blue")
  
  lines(plot_dates[idx], max_w_mkt_s[idx],
        lwd = 2, lty = 4, col = "darkgreen")
  
  legend("topright",
         legend = c("Black-Litterman", "Naive MV", "Views Only", "Market Size"),
         col = c("red", "purple", "blue", "darkgreen"),
         lwd = 2,
         lty = c(1, 2, 3, 4),
         bty = "n")
  
  dev.off()
}


# -----------------------------
# Average weights table
# -----------------------------

avg_weights_feedback <- data.frame(
  Industry = ind_names,
  BL_AvgWeight = colMeans(wBL, na.rm = TRUE),
  NaiveMV_AvgWeight = colMeans(wNV, na.rm = TRUE),
  ViewsOnly_AvgWeight = colMeans(wVIEW, na.rm = TRUE),
  MarketSize_AvgWeight = colMeans(wMKT, na.rm = TRUE)
)

cat("\n==============================\n")
cat("TOP 10 AVERAGE WEIGHTS\n")
cat("==============================\n")

cat("\nBlack-Litterman\n")
print(head(avg_weights_feedback[order(-avg_weights_feedback$BL_AvgWeight), ], 10))

cat("\nNaive MV\n")
print(head(avg_weights_feedback[order(-avg_weights_feedback$NaiveMV_AvgWeight), ], 10))

cat("\nViews Only\n")
print(head(avg_weights_feedback[order(-avg_weights_feedback$ViewsOnly_AvgWeight), ], 10))

cat("\nMarket Size\n")
print(head(avg_weights_feedback[order(-avg_weights_feedback$MarketSize_AvgWeight), ], 10))

cat("\nPNG files exported to:", graph_dir, "\n")
