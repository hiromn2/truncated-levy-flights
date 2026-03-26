# =============================================================================
# 05_animations.R
#
# Bypassing the Truncation Problem of Truncated Levy Flights
#
# PURPOSE:
#   Animated versions of the key figures using gganimate + magick.
#   Requires that 02_main_analysis.R has been run first.
#
# OUTPUT:
#   output/figures/animation_crisis.gif  -- intraday X_t on 2015-01-15
#   output/figures/animation_scaling.gif -- animated log c_t vs log nu_t
# =============================================================================

library(gganimate)
library(magick)

# NOTE: Run 00_load_data.R and 02_main_analysis.R first.
# source("R/00_load_data.R")
# source("R/02_main_analysis.R")

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# -- Shared utilities ----------------------------------------------------------

scientific_10 <- function(x) {
  parse(text = gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

pub_theme <- theme(
  text         = element_text(family = "Times New Roman", size = 21, colour = "black"),
  axis.title.y = element_text(angle = 90, size = 22, colour = "black",
                               margin = margin(r = 15, unit = "pt")),
  axis.title.x = element_text(vjust = -1.5,
                               margin = margin(b = 15, unit = "pt"))
)

# =============================================================================
# Animation 1: USD/CHF intraday bid price on the SNB crisis day (2015-01-15)
# =============================================================================

D.crisis.15 <- data.1[data.1$date == "2015-01-15", ]
bks.crisis  <- c(4689800, 4690300, 4690800)

fig_crisis <- ggplot(D.crisis.15, aes(x = time, y = V3)) +
  geom_line() +
  xlab(expression(italic("t"))) +
  ylab(expression(italic("X")[italic("t")])) +
  scale_x_continuous(labels = scientific_10, breaks = bks.crisis) +
  pub_theme

fig_crisis_anim <- fig_crisis +
  geom_point() +
  transition_reveal(time)   # reveal data along the time axis

anim1 <- animate(fig_crisis_anim, fps = 10, duration = 8,
                 width = 600, height = 500)

magick::image_write(anim1, path = "output/figures/animation_crisis.gif")
message("Saved: output/figures/animation_crisis.gif")

# =============================================================================
# Animation 2: Animated scatter of log c_t vs log nu_t
#   Points are revealed in order of increasing log nu_t,
#   showing how the power-law scaling builds up.
# =============================================================================

data.plot <- data.zeta.fit[order(data.zeta.fit$log.nu.t), ]
data.plot <- na.omit(data.plot)

# Scatter layer
data.scatter <- data.frame(
  x     = data.plot$log.nu.t,
  y     = data.plot$log.c.t,
  group = "data"
)

# Fitted line layer
data.line <- data.frame(
  x     = data.plot$log.nu.t,
  y     = data.plot$fitted.line,
  group = "fit"
)

fig_scaling <- ggplot(data.scatter, aes(x = x, y = y)) +
  geom_point(colour = "black", size = 0.8) +
  geom_line(data = data.line, aes(x = x, y = y), colour = "red", linewidth = 1) +
  xlab(expression(paste("log ") * hat(sigma)[italic("t")])) +
  ylab(expression(paste("log ") * hat(italic("l"))[italic("t")])) +
  pub_theme

fig_scaling_anim <- fig_scaling +
  transition_time(data.scatter$x) +
  shadow_mark()

anim2 <- animate(fig_scaling_anim, fps = 10, duration = 8,
                 width = 600, height = 500)

magick::image_write(anim2, path = "output/figures/animation_scaling.gif")
message("Saved: output/figures/animation_scaling.gif")
