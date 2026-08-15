rm(list = ls())

# ============================================================
# Black-Litterman vs Naive Mean-Variance MOMENTUM VIEWS
# 47 Fama-French Industry Portfolios - Daily
# Uses LN19 formula: mu_E = g V w_M
# Then mu_BL = [ (tau V)^-1 + P' D^-1 P ]^-1 [ (tau V)^-1 mu_E + P' D^-1 Q ]
# ============================================================

# -----------------------------
# Load data
# -----------------------------

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
  
  # Relative views Q from momentum signal
  recent_mu <- colMeans(Rwin[(window - view_window + 1):window, , drop = FALSE]) * dpy
  
  spread1 <- recent_mu["Softw"] - recent_mu["Telcm"]
  spread2 <- recent_mu["Chips"] - recent_mu["Util"]
  spread3 <- recent_mu["Hlth"]  - recent_mu["Steel"]
  
  view_size <- 0.03   # 3% annualized relative view
  
  Q <- matrix(c(
    ifelse(spread1 > 0,  view_size, -view_size),
    ifelse(spread2 > 0,  view_size, -view_size),
    ifelse(spread3 > 0,  view_size, -view_size)
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
