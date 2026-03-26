# =============================================================================
# 01_log_returns.R
#
# Bypassing the Truncation Problem of Truncated Levy Flights
#
# PURPOSE:
#   Compute log-returns w(Delta.t) at multiple time scales and verify that
#   their distribution is consistent with a (truncated) Levy stable process.
#
# THEORY:
#   If X_t follows a Levy-stable process with index alpha, the log-price
#   x_t = log(X_t) has increments:
#
#       w(Delta.t) = x_{t+Delta.t} - x_t
#
#   whose mean scales as  E[w] ~ omega * Delta.t^beta  (drift exponent beta)
#   and whose density at zero satisfies  P(0|Delta.t) ~ Delta.t^{-1/alpha}.
#
#   This script:
#     [1] Computes w(Delta.t) for Delta.t in {1, 60, 120, 240, 360, 720, 1440} min
#     [2] Fits the drift scaling:  log E[w] ~ beta * log(Delta.t)
#     [3] Estimates omega, the drift coefficient
#     [4] Centers the returns:  W = w - omega * Delta.t^beta
#     [5] Plots the log-density of W (should collapse across scales for stable)
# =============================================================================

source("R/00_load_data.R")

# Only use the 2015 sub-sample (Jan-Jun) as in the original paper
# (the SNB crisis of Jan 15, 2015 is the key event)
DATA <- with(data.1[data.1$date >= "2015-01-01", ], V4)  # ask price

# -- [1] Log-price and returns at multiple Delta.t ----------------------------

x       <- log(DATA)
n       <- length(x)

# Time scales in minutes
delta.t <- c(60, 120, 240, 360, 720, 1440)

# w(Delta.t)_t = x_{t+Delta.t} - x_t  (lagged difference)
log_return <- function(dt, x, n) x[(1 + dt):n] - x[1:(n - dt)]

w.0 <- log_return(1,          x, n)   #  1 min
w.1 <- log_return(delta.t[1], x, n)   #  1 h
w.2 <- log_return(delta.t[2], x, n)   #  2 h
w.3 <- log_return(delta.t[3], x, n)   #  4 h
w.4 <- log_return(delta.t[4], x, n)   #  6 h
w.5 <- log_return(delta.t[5], x, n)   # 12 h
w.6 <- log_return(delta.t[6], x, n)   #  1 d

# -- [2] KDE of log-returns ---------------------------------------------------

f.hat <- list(
  w.0 = density(w.0),
  w.1 = density(w.1),
  w.2 = density(w.2),
  w.3 = density(w.3),
  w.4 = density(w.4),
  w.5 = density(w.5),
  w.6 = density(w.6)
)

# -- [3] Moments table --------------------------------------------------------

library(moments)

skwns   <- sapply(list(w.0, w.1, w.2, w.3, w.4, w.5, w.6), skewness)
krt     <- sapply(list(w.0, w.1, w.2, w.3, w.4, w.5, w.6), kurtosis)
w.mean  <- sapply(list(w.1, w.2, w.3, w.4, w.5, w.6), mean)
w.meansd <- w.mean / sapply(list(w.1, w.2, w.3, w.4, w.5, w.6), sd)

# -- [4] Drift scaling: log E[w] ~ beta * log(Delta.t) -----------------------
# Figure 1

model.drift <- lm(log(w.mean) ~ log(delta.t))
beta.hat    <- coef(model.drift)[2]   # drift exponent

plot(log(delta.t), log(w.mean),
     ylab = "log E[w]", xlab = "log(Delta.t)",
     main = "Drift scaling (Figure 1a)")
abline(model.drift, col = "red")

model.cv <- lm(log(w.meansd) ~ log(delta.t))

plot(log(delta.t), log(w.meansd),
     ylab = "log(mean/sd)", xlab = "log(Delta.t)",
     main = "Coefficient of variation scaling (Figure 1b)")
abline(model.cv, col = "red")

# -- [5] Omega (drift coefficient) at each scale ------------------------------
# omega(Delta.t) = E[w(Delta.t)] / Delta.t^beta

get_omega <- function(w, beta, dt) mean(w) / (dt^beta)

omega.1 <- get_omega(w.1, beta.hat, delta.t[1])
omega.2 <- get_omega(w.2, beta.hat, delta.t[2])
omega.3 <- get_omega(w.3, beta.hat, delta.t[3])
omega.4 <- get_omega(w.4, beta.hat, delta.t[4])
omega.5 <- get_omega(w.5, beta.hat, delta.t[5])
omega.6 <- get_omega(w.6, beta.hat, delta.t[6])

omega <- mean(c(omega.1, omega.2, omega.3, omega.4, omega.5, omega.6))

# -- [6] Centered returns W = w - omega * Delta.t^beta -----------------------

center_return <- function(w, dt, omega, beta) w - omega * (dt^beta)

W.0 <- center_return(w.0,          1, omega, beta.hat)
W.1 <- center_return(w.1, delta.t[1], omega, beta.hat)
W.2 <- center_return(w.2, delta.t[2], omega, beta.hat)
W.3 <- center_return(w.3, delta.t[3], omega, beta.hat)
W.4 <- center_return(w.4, delta.t[4], omega, beta.hat)
W.5 <- center_return(w.5, delta.t[5], omega, beta.hat)
W.6 <- center_return(w.6, delta.t[6], omega, beta.hat)

# -- [7] Log-density of centered returns (Figure 3) --------------------------
# Under a Levy-stable process, these log-densities should collapse (be parallel)
# after rescaling by Delta.t^{1/alpha}.

fw <- lapply(list(W.0, W.1, W.2, W.3, W.4, W.5, W.6), density)

plot(fw[[1]]$x, log(fw[[1]]$y),
     type = "l", col = "red",
     xlim = c(-0.05, 0.20),
     ylab = "log P(W)", xlab = "W",
     main = "Log-density of centered returns (Figure 3)")
for (i in 2:7) lines(fw[[i]]$x, log(fw[[i]]$y))

# -- [8] Estimate alpha from P(0 | Delta.t) -----------------------------------
# P(0|Delta.t) is estimated as the KDE evaluated at 0, using each scale's bandwidth.

KDE_at_zero <- function(bw, w) {
  z <- (0 - w) / bw
  mean(dnorm(z)) / bw
}

p.0 <- sapply(seq_along(delta.t), function(i) {
  KDE_at_zero(f.hat[[i + 1]]$bw, list(w.1, w.2, w.3, w.4, w.5, w.6)[[i]])
})

model.alpha <- lm(log(p.0) ~ log(delta.t))
alpha.hat   <- -1 / coef(model.alpha)[2]   # KEY RESULT

plot(log(delta.t), log(p.0),
     ylab = "log P(0|Delta.t)", xlab = "log(Delta.t)",
     main = sprintf("Alpha estimation: alpha = %.3f (Figure 4)", alpha.hat))
abline(model.alpha, col = "red")

cat(sprintf("\nEstimated Levy index:  alpha = %.4f\n", alpha.hat))
cat(sprintf("Estimated drift exp:   beta  = %.4f\n", beta.hat))
cat(sprintf("Estimated drift coeff: omega = %.6f\n", omega))

# -- [9] Levy-rescaled returns ------------------------------------------------
# W_S(Delta.t) = W(Delta.t) / Delta.t^{1/alpha}
# Under a pure Levy stable process, all W_S should have the same distribution.

levy_rescale <- function(W, dt, alpha) W / (dt^(1/alpha))

WS.0 <- levy_rescale(W.0,          1, alpha.hat)
WS.1 <- levy_rescale(W.1, delta.t[1], alpha.hat)
WS.2 <- levy_rescale(W.2, delta.t[2], alpha.hat)
WS.3 <- levy_rescale(W.3, delta.t[3], alpha.hat)
WS.4 <- levy_rescale(W.4, delta.t[4], alpha.hat)
WS.5 <- levy_rescale(W.5, delta.t[5], alpha.hat)
WS.6 <- levy_rescale(W.6, delta.t[6], alpha.hat)

fws <- lapply(list(WS.0, WS.1, WS.2, WS.3, WS.4, WS.5, WS.6), density)

plot(fws[[1]]$x, log(fws[[1]]$y),
     type = "l", col = "red",
     ylab = "log P(W_S)", xlab = "W_S",
     main = "Log-density after Levy rescaling (Figure 5)\nShould collapse to one curve")
for (i in 2:7) lines(fws[[i]]$x, log(fws[[i]]$y))
