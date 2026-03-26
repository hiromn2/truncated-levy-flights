# =============================================================================
# 04_figures.R
#
# Bypassing the Truncation Problem of Truncated Levy Flights
#
# PURPOSE:
#   Generate all publication-quality figures (ggplot2).
#   Requires that 02_main_analysis.R has been run first (objects in memory).
#
# FIGURES:
#   Figure 1(a): Full USD/CHF time series X_t
#   Figure 1(b): Intraday series on the SNB crisis day (2015-01-15)
#   Figure 2(a): Intraday return r_t = X_t/m_t - 1
#   Figure 2(b): Normalized return z_t = r_t / nu_t^beta
#   Figure 3:    log c_t vs log nu_t (power-law scaling)
#   Figure 4:    Frechet fit: log(-log F(z)) vs log(z)
# =============================================================================

source("R/00_load_data.R")
# NOTE: Run 02_main_analysis.R first to populate data.rt, data.zeta.fit,
# F.list, alpha.hat, gamma.hat, beta, data.rt.2 in the global environment.
# If running this file standalone, uncomment the next line:
# source("R/02_main_analysis.R")

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# -- Shared theme and utilities ------------------------------------------------

# Scientific notation axis labels: 4.69e6 -> 4.69 x 10^6
scientific_10 <- function(x) {
  parse(text = gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

# Publication theme (Times New Roman, 21pt; falls back if font not installed)
pub_theme <- theme(
  text         = element_text(family = "Times New Roman", size = 21, colour = "black"),
  axis.title.y = element_text(family = "Times New Roman", angle = 90, size = 22,
                               colour = "black",
                               margin = margin(t = 0, r = 15, b = 0, l = 0, unit = "pt")),
  axis.title.x = element_text(family = "Times New Roman", vjust = -1.5,
                               margin = margin(t = 0, r = 0, b = 15, l = 0, unit = "pt"))
)

save_fig <- function(fig, filename, width = 14, height = 12, units = "cm", dpi = 300) {
  ggsave(file.path("output/figures", filename), fig,
         width = width, height = height, units = units, dpi = dpi)
  message("Saved: output/figures/", filename)
}

# =============================================================================
# Figure 1(a): Full USD/CHF bid price series X_t
# =============================================================================

n           <- nrow(data.1)
data.1$time <- seq_len(n)
bks         <- c(1, 1500000, 3000000, 4500000)

fig1a <- ggplot(data.1, aes(x = time, y = V3)) +
  geom_line() +
  xlab(expression(italic("t"))) +
  ylab(expression(italic("X")[italic("t")])) +
  scale_x_continuous(labels = scientific_10, breaks = bks) +
  pub_theme

print(fig1a)
save_fig(fig1a, "CHF_Xt.png")

# =============================================================================
# Figure 1(b): Intraday series on the SNB crisis day (2015-01-15)
#   On this day, the Swiss National Bank abandoned its EUR/CHF floor,
#   causing one of the largest single-day FX moves in history.
# =============================================================================

D.crisis.15 <- data.1[data.1$date == "2015-01-15", ]
bks.crisis  <- c(4689800, 4690300, 4690800)

fig1b <- ggplot(D.crisis.15, aes(x = time, y = V3)) +
  geom_line() +
  xlab(expression(italic("t"))) +
  ylab(expression(italic("X")[italic("t")])) +
  scale_x_continuous(labels = scientific_10, breaks = bks.crisis) +
  pub_theme

print(fig1b)
save_fig(fig1b, "CHF_Xt_150115.png")

# =============================================================================
# Figure 2(a): Intraday return r_t = X_t / m_t - 1
# =============================================================================

bks <- c(1, 1500000, 3000000, 4500000)

fig2a <- ggplot(data.rt, aes(x = time, y = r.t)) +
  geom_line() +
  xlab(expression(italic("t"))) +
  ylab(expression(italic("r")[italic("t")])) +
  scale_x_continuous(labels = scientific_10, breaks = bks) +
  pub_theme

print(fig2a)
save_fig(fig2a, "CHF_rt.png")

# =============================================================================
# Figure 2(b): Normalized return z_t = r_t / nu_t^beta
# =============================================================================

fig2b <- ggplot(data.rt.2, aes(x = time.x, y = z.t)) +
  geom_line() +
  xlab(expression(italic("t"))) +
  ylab(expression(italic("Z")[italic("t")])) +
  scale_x_continuous(labels = scientific_10, breaks = bks) +
  pub_theme

print(fig2b)
save_fig(fig2b, "CHF_Zt.png")

# =============================================================================
# Figure 3: log c_t vs log nu_t  (power-law scaling of the truncation bound)
# =============================================================================

# Sort for clean line
data.plot3 <- data.zeta.fit[order(data.zeta.fit$log.nu.t), ]
data.plot3 <- na.omit(data.plot3)

fig3 <- ggplot(data.plot3, aes(x = log.nu.t, y = log.c.t)) +
  geom_point(colour = "black", size = 0.8) +
  geom_line(aes(y = fitted.line), colour = "red", linewidth = 1) +
  xlab(expression(paste("log ") * hat(sigma)[italic("t")])) +
  ylab(expression(paste("log ") * hat(italic("l"))[italic("t")])) +
  pub_theme

print(fig3)
save_fig(fig3, "CHF_log_ct_vs_nu.png")

# =============================================================================
# Figure 4: Frechet fit — log(-log F(z)) vs log(z)
# =============================================================================

fig4 <- ggplot(F.list, aes(x = log.z, y = loglog.F.z)) +
  geom_point() +
  geom_line(aes(y = fitted.line), colour = "red", linetype = "solid", linewidth = 1) +
  xlab(expression(paste("ln[") * italic("z") * paste("]"))) +
  ylab(expression(paste("ln[-ln ") * hat(italic("F")) *
                    paste("(") * italic("z") * paste("]"))) +
  pub_theme

print(fig4)
save_fig(fig4, "CHF_Frechet_fit.png")

message("\nAll figures saved to output/figures/")
