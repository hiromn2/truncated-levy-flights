# =============================================================================
# 02_main_analysis.R
#
# Bypassing the Truncation Problem of Truncated Levy Flights
#
# PURPOSE:
#   Main empirical pipeline for the paper "On the Boundedness of an Exchange Rate."
#   Implements:
#     [A] Intraday returns r_t = X_t / m_t - 1  (normalized by daily mean)
#     [B] Daily log-variance sigma^2_t
#     [C] Daily truncation bound c_t = max_t |r_t|
#     [D] Daily volatility nu_t = sqrt(Var(X_t)) / m_t
#     [E] Power-law fit: log c_t ~ beta * log nu_t  (identifies beta)
#     [F] Normalized bound zeta_t = c_t / nu_t^beta
#     [G] Empirical CDF of zeta
#     [H] Frechet distribution fit:  F(z) = exp(-gamma * z^{-alpha})
#     [I] Tail index via the NLS scaling method (bypassing truncation)
#
# THEORY:
#   Under the paper's model, the normalized bound Z_t = c_t / nu_t^beta
#   follows a Frechet distribution (Type II extreme value):
#
#       F(z; alpha, gamma) = exp(-gamma * z^{-alpha}),   z > 0
#
#   The Frechet arises because c_t is the maximum of a heavy-tailed intraday
#   return process. The log-log transformation linearizes the CDF:
#
#       log(-log F(z)) = log(gamma) - alpha * log(z)
#
#   so a simple OLS on the empirical log(-log CDF) vs log(z) identifies
#   alpha and gamma.
#
#   The NLS method in block [I] additionally exploits the self-similarity of
#   the tail:  P(|Z| > eta*z) = eta^{-alpha} * P(|Z| > z)  for Pareto tails,
#   fitting alpha directly via nonlinear least squares.
# =============================================================================

source("R/00_load_data.R")

# -----------------------------------------------------------------------------
# [A] Intraday return r_t = X_t / m_t - 1
#     where m_t = daily mean of the bid price X_t
# -----------------------------------------------------------------------------
# daily mean
data.mean           <- aggregate(V3 ~ date, data = data.1, mean, na.rm = TRUE)
names(data.mean)[2] <- "m.t"

# merge back and compute intraday return
data.rt             <- merge(data.1, data.mean, by = "date", all.x = TRUE)
data.rt$r.t         <- data.rt$V3 / data.rt$m.t - 1
data.rt$abs.r.t     <- abs(data.rt$r.t)
data.rt$time        <- seq_len(nrow(data.rt))

# -----------------------------------------------------------------------------
# [B] Daily log-variance sigma^2_t = Var(log X_t) within each day
# -----------------------------------------------------------------------------
data.rt$log.X.t      <- log(data.rt$V3)
data.sigma           <- aggregate(log.X.t ~ date, data = data.rt, var, na.rm = TRUE)
names(data.sigma)[2] <- "sigma2.t"
data.sigma$time      <- seq_len(nrow(data.sigma))

# -----------------------------------------------------------------------------
# [C] Daily truncation bound c_t = max_t |r_t|
# -----------------------------------------------------------------------------
data.ct           <- aggregate(abs.r.t ~ date, data = data.rt, max, na.rm = TRUE)
names(data.ct)[2] <- "c.t"

# -----------------------------------------------------------------------------
# [D] Daily volatility nu_t = sqrt(Var(X_t)) / m_t
#     (coefficient of variation of the price level within each day)
# -----------------------------------------------------------------------------
data.M2t            <- aggregate(V3 ~ date, data = data.rt, var, na.rm = TRUE)
names(data.M2t)[2]  <- "M2.t"
nu.t                <- sqrt(data.M2t$M2.t) / data.mean$m.t
data.nut            <- data.frame(nu.t, time = seq_len(length(nu.t)))

# -----------------------------------------------------------------------------
# [E] Power-law fit: log c_t ~ beta * log nu_t
#     Identifies beta (the volatility scaling exponent)
# -----------------------------------------------------------------------------
log.c.t   <- log(data.ct$c.t)
log.nu.t  <- log(data.nut$nu.t)
sigma.t.2 <- data.nut$nu.t^2

data.zeta.fit <- data.frame(log.c.t, log.nu.t, sigma.t.2)
data.zeta.fit <- data.zeta.fit[is.finite(data.zeta.fit$log.nu.t), ]

line <- lm(log.c.t ~ log.nu.t, data = data.zeta.fit, na.action = na.omit)
data.zeta.fit$fitted.line <- fitted(line)
summary(line)

beta <- coef(line)[2]
cat(sprintf("\nEstimated scaling exponent: beta = %.4f\n", beta))

# -----------------------------------------------------------------------------
# [F] Normalized bound zeta_t = c_t / nu_t^beta
# -----------------------------------------------------------------------------
zeta      <- data.ct$c.t / (data.nut$nu.t^beta)
data.zeta <- data.frame(zeta, time = seq_len(length(zeta)))

# -----------------------------------------------------------------------------
# [G] Empirical CDF of zeta
# -----------------------------------------------------------------------------
freqs        <- table(zeta)
z.values     <- as.numeric(names(freqs))
log.z.values <- log(z.values)
F.z          <- cumsum(freqs) / sum(freqs)
F.list       <- data.frame(log.z.values, z.values, F.z)

# Keep only interior values (avoid F=1 giving log(-log(0)) = -Inf)
F.list <- F.list[F.list$F.z < 1 & F.list$z.values > 1, ]

# -----------------------------------------------------------------------------
# [H] Frechet fit: log(-log F(z)) = log(gamma) - alpha * log(z)
#
#     The Frechet CDF is  F(z) = exp(-gamma * z^{-alpha}).
#     Taking log twice:   log(-log F(z)) = log(gamma) - alpha * log(z).
#     So OLS of log(-log F) on log(z) gives:
#       intercept -> log(gamma) -> gamma = exp(intercept)
#       slope     -> -alpha     -> alpha = -slope
# -----------------------------------------------------------------------------
F.list$log.z      <- log(F.list$z.values)
F.list$loglog.F.z <- log(-log(F.list$F.z))

fit     <- lm(loglog.F.z ~ log.z, data = F.list)
F.list$fitted.line <- fitted(fit)
summary(fit)

gamma.hat <- exp(coef(fit)[1])
alpha.hat <- -coef(fit)[2]

cat(sprintf("Frechet parameters:  alpha = %.4f,  gamma = %.4f\n",
            alpha.hat, gamma.hat))

# Implied survival probabilities at selected quantiles
z.grid <- c(1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
F.z.fit <- exp(-gamma.hat * (z.grid^(-alpha.hat)))
print(data.frame(z = z.grid,
                 F.z = round(F.z.fit, 6),
                 tail.prob = format(1 - F.z.fit, scientific = FALSE)))

# -----------------------------------------------------------------------------
# [I] Tail-index via NLS self-similarity: R(eta*z) = eta^{-alpha} * R(z)
#
#     For a Pareto tail, P(|Z| > z) ~ z^{-alpha}, so:
#         P(|Z| > eta*z) = eta^{-alpha} * P(|Z| > z)
#
#     We estimate R(z) = P(|Z| > z) from the empirical distribution of
#     |z_t| = |r_t| / nu_t^beta, then fit alpha by NLS over a grid of
#     (z, eta) pairs. This is the "bypassing" step: it works even when the
#     distribution is truncated, because truncation only removes the far tail
#     and the scaling still holds in the intermediate range.
# -----------------------------------------------------------------------------
data.nut$date      <- data.mean$date
data.nut$nu.t.beta <- data.nut$nu.t^beta

data.rt.2          <- merge(data.rt, data.nut, by = "date", all.x = TRUE)
data.rt.2$abs.z    <- data.rt.2$abs.r.t / data.rt.2$nu.t.beta
data.rt.2$z.t      <- data.rt.2$r.t / data.rt.2$nu.t.beta
data.rt.2          <- data.rt.2[!is.na(data.rt.2$abs.z), ]

max.z   <- round(max(data.rt.2$abs.z), 3)
zeta.0  <- 1.0
eta.seq <- seq(0.7, 1.2, by = 0.1)

# Build the (q, eta) grid and compute empirical tail probabilities
k <- 1
q.abs.z <- R.z <- R.z.eta <- eta.z <- NULL

for (eta.grid in eta.seq) {
  for (q.z in seq(zeta.0, max.z, by = 0.01)) {
    R.z[k]     <- mean(data.rt.2$abs.z > q.z)
    R.z.eta[k] <- mean(data.rt.2$abs.z > eta.grid * q.z)
    q.abs.z[k] <- q.z
    eta.z[k]   <- eta.grid
    k <- k + 1
  }
}

data.R.fit       <- data.frame(q.abs.z, eta.z, R.z.eta, R.z)
data.R.fit$eta.g <- as.factor(data.R.fit$eta.z)

# NLS fit: R(z) = eta^alpha * R(eta*z)
fit.alpha <- nls(R.z ~ (eta.z^alpha) * R.z.eta,
                 data  = data.R.fit,
                 start = list(alpha = 2))

alpha.nls <- coef(fit.alpha)[1]
cat(sprintf("NLS tail index estimate:  alpha = %.4f\n", alpha.nls))

# Drop rows where R.z.eta == 0 (log undefined) before plotting
data.R.fit.plot <- data.R.fit[data.R.fit$R.z.eta > 0, ]

p_tail_scaled <- ggplot(data.R.fit.plot,
                        aes(x = log(q.abs.z),
                            y = log((eta.z^alpha.nls) * R.z.eta),
                            group = eta.g)) +
  geom_point(aes(shape = eta.g)) +
  xlab(expression("ln[" * italic("z") * "]")) +
  ylab(expression("ln[" * hat("P") * "(|" * italic("Z")[italic("t")] *
                    "| > " * italic("z") * ")]")) +
  labs(shape = expression(eta)) +
  theme(text = element_text(family = "Times New Roman", size = 21, colour = "black"),
        axis.title.y = element_text(angle = 90, size = 22, colour = "black",
                                    margin = margin(r = 15, unit = "pt")),
        axis.title.x = element_text(vjust = -1.5,
                                    margin = margin(b = 15, unit = "pt")),
        legend.title = element_text(colour = "black"), legend.title.align = 0.5)

print(p_tail_scaled)

# Save
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("output/figures/CHF_CaudalIndex2.png", p_tail_scaled,
       width = 14, height = 12, units = "cm", dpi = 300)
