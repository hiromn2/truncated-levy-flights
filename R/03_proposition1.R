# =============================================================================
# 03_proposition1.R
#
# Bypassing the Truncation Problem of Truncated Levy Flights
#
# PURPOSE:
#   Numerical verification of Proposition 1 in the paper.
#   Studies the variance of the truncated Cauchy and truncated Normal
#   distributions as a function of the truncation parameter c.
#
# THEORY (Truncated Cauchy):
#   Let X ~ Cauchy(0,1) truncated to [-c, c]. The normalizing constant is:
#
#       eta(c) = 2 * arctan(c)
#
#   The (unnormalized) second moment integral is:
#
#       V(c) = integral_{-c}^{c} x^2 / (pi*(1+x^2)) dx = c - arctan(c)
#
#   So the variance of the truncated Cauchy is:
#
#       Var(c) = 2*V(c)/eta(c) = 2*(c - arctan(c)) / (2*arctan(c))
#              = (c - arctan(c)) / arctan(c)
#              = c / arctan(c) - 1
#
#   Proposition 1 establishes that the second derivative d^2 Var / d(eta^2)
#   has a specific sign, which is key to the identification argument in the paper.
#
# THEORY (Truncated Normal):
#   Let X ~ N(0,1) truncated to [-c, c]. The variance is:
#
#       Var(c) = 1 - c*phi(c)/Phi(c)*2/(2*Phi(c)-1)... (standard formula)
#
#   where phi = standard normal PDF, Phi = standard normal CDF.
#   For symmetric truncation [-c, c]:
#
#       Var(c) = 1 + (a*phi(a) - b*phi(b))/(Phi(b)-Phi(a))
#                  - ((phi(a)-phi(b))/(Phi(b)-Phi(a)))^2
#   with a = -c, b = c  (the second term vanishes since phi(-c) = phi(c)).
# =============================================================================

# =============================================================================
# PART 1: Truncated Cauchy
# =============================================================================

c.grid <- seq(0, 5, by = 0.001)

# -- Variance of truncated Cauchy: Var(c) = c/arctan(c) - 1 -----------------
var.c <- c.grid / atan(c.grid) - 1

# Standard deviation as a function of c
plot(c.grid, sqrt(var.c),
     type = "l", xlab = "c", ylab = "SD(c)",
     main = "Standard deviation of Truncated Cauchy")

# -- Derivatives of eta(c) = 2*arctan(c) -------------------------------------
# eta  = 2*arctan(c)   [normalizing constant]
# eta' = 2/(1+c^2)
# eta''= -4c/(1+c^2)^2

eta.c    <- 2 * atan(c.grid)
d.eta.c  <- 2 / (1 + c.grid^2)
d2.eta.c <- -2 * c.grid / (1 + c.grid^2)^2

# -- Second moment integral V(c) = c - arctan(c) and its derivatives ---------
# V(c)  = c - arctan(c)
# V'(c) = 1 - 1/(1+c^2) = c^2/(1+c^2)
# V''(c)= 2c/(1+c^2)^2

V.c    <- c.grid - atan(c.grid)
d.V.c  <- 1 - 1 / (1 + c.grid^2)
d2.V.c <- 2 * c.grid / (1 + c.grid^2)^2

# -- Variance and its derivative with respect to eta --------------------------
# Var(c) = 2*V(c)/eta(c)
# d Var/d eta = (d Var/dc) / (d eta/dc)
#
# d Var/dc = (d.V.c * eta.c - V.c * d.eta.c) / eta.c^2   [quotient rule on 2V/eta]
# Multiply by 2 since Var = 2V/eta:

Var.c  <- 2 * V.c / eta.c

# dVar/dc
d.Var.dc <- 2 * (d.V.c * eta.c - V.c * d.eta.c) / eta.c^2

# dVar/deta  = (dVar/dc) / (deta/dc)
d.Var.d.eta <- d.Var.dc / d.eta.c

# Second derivative d^2 Var / d eta^2 (by chain and quotient rule)
# Components of the second derivative (labeled to match original code):
ratio.4 <- d2.V.c / eta.c                       # from V'' term
ratio.5 <- d.V.c * d.eta.c / (eta.c^2)          # cross term
ratio.6 <- d.Var.dc * d.eta.c / eta.c            # Var' * eta' term
ratio.7 <- Var.c * ((d2.eta.c / eta.c) - (d.eta.c / eta.c)^2)  # Var * (eta''/eta - (eta'/eta)^2)

# d^2 Var / d eta^2
d2.Var <- ratio.4 - ratio.5 - ratio.6 - ratio.7

plot(c.grid, d2.Var,
     type = "l", xlab = "c", ylab = "d²Var/deta²",
     main = "Second derivative of Var w.r.t. eta (Truncated Cauchy)")
abline(h = 0, lty = 2, col = "red")

# Diagnostic: ratio Var / (c^2)  -- measures how close Var is to c^2
ratio.8 <- 2 * Var.c / (c.grid^2)
plot(c.grid, ratio.8,
     type = "l", xlab = "c", ylab = "Var / c²",
     main = "Var(c) relative to c² (Truncated Cauchy)")

# =============================================================================
# PART 2: Truncated Normal
# =============================================================================

c.grid <- seq(0, 6, by = 0.1)
a      <- -c.grid   # lower truncation point
b      <-  c.grid   # upper truncation point

phi.a <- dnorm(a)   # phi(-c) = phi(c)
phi.b <- dnorm(b)
Phi.a <- pnorm(a)   # Phi(-c)
Phi.b <- pnorm(b)   # Phi(c)

# General truncated normal variance formula (Greene, 2003):
#   Var = 1 + (a*phi(a) - b*phi(b))/(Phi(b)-Phi(a))
#           - ((phi(a) - phi(b))/(Phi(b)-Phi(a)))^2
# For symmetric [-c,c]: phi(a)=phi(b) so the last term vanishes.
Var.tn <- 1 + (a * phi.a - b * phi.b) / (Phi.b - Phi.a) -
              ((phi.a - phi.b) / (Phi.b - Phi.a))^2

# Alternative closed-form for symmetric case:
#   Var = 1 - c*sqrt(2) * phi(c) / ((1 - 2*Phi(-c)) * sqrt(pi))  ... (approx)
Var.tn2 <- 1 - (c.grid * sqrt(2)) / ((1 - 2 * Phi.a) * sqrt(pi)) * exp(-(c.grid^2) / 2)

# SD of truncated normal
plot(c.grid, sqrt(Var.tn),
     type = "l", xlab = "c", ylab = "SD(c)",
     main = "Standard deviation of Truncated Normal")

# Log-log plot to inspect power-law behavior at large c
plot(log(c.grid[-1]), log(sqrt(Var.tn[-1])),
     type = "l", xlab = "log c", ylab = "log SD(c)",
     main = "Log-log: SD of Truncated Normal")

# Power-law approximation: SD ~ c^gamma
fit.loglog <- lm(log(sqrt(Var.tn[-1])) ~ log(c.grid[-1]))
cat(sprintf("\nLog-log slope for truncated Normal SD: %.4f\n",
            coef(fit.loglog)[2]))

# Linear approximation: c ~ k * SD  (used in the paper's calibration)
fit.linear <- lm(c.grid ~ sqrt(Var.tn) + 0)   # force through origin
k.hat      <- coef(fit.linear)[1]
cat(sprintf("Linear approximation: c ≈ %.4f * SD (theoretical: sqrt(3) ≈ %.4f)\n",
            k.hat, sqrt(3)))

# Overlay on SD plot
plot(sqrt(Var.tn), c.grid,
     type = "l", xlab = "SD(c)", ylab = "c",
     main = "c vs SD: data and linear fit")
lines(sqrt(Var.tn), k.hat * sqrt(Var.tn), col = "red", lty = 2)
legend("topleft", legend = c("Exact", sprintf("c ≈ %.2f*SD", k.hat)),
       col = c("black", "red"), lty = 1:2)
