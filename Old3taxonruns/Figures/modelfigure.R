# Load packages
library(ggplot2)
library(gridExtra)
library(ggpattern)

# Reversed x-axis limits
x_limits <- c(4.25, 0)
y_max <- 0.5  # consistent max height for all plots

# Blue shades
blue1 <- "#7ec8e3"  # light
blue2 <- "#2a6f97"  # dark

# -------- Graph 1: One striped uniform from 0.1 to 3.5 --------
df1 <- data.frame(
  xmin = 0.1,
  xmax = 3.5,
  ymin = 0,
  ymax = y_max
)

p1 <- ggplot(df1) +
  geom_rect_pattern(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = blue1,           # Light blue background
    color = NA,             # No rectangle outline
    pattern = 'stripe',
    pattern_fill = blue2,   # Dark blue stripes
    pattern_colour = NA,    # No stripe outline
    pattern_angle = 45,
    pattern_density = 0.2,
    pattern_spacing = 0.25,
    alpha = 0.8
  ) +
  scale_x_reverse(limits = x_limits) +
  coord_cartesian(ylim = c(0, y_max)) +
  ggtitle("Unknown") +
  labs(y = "Density", x = "") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),  # removes minor gridlines
    panel.grid.major = element_line(color = "gray80", size = 0.2)  # keep major gridlines
  ) +
  geom_vline(xintercept = 3, color = "darkred", linetype = "dotted", size = 1) +
  geom_vline(xintercept = 1, color = "darkgreen", linetype = "dotted", size = 1) +
  guides(fill = "none")  # Remove the group legend

# -------- Graph 2: Two uniforms from 3.5–2.5 and 1.5–0.5 --------
df2 <- data.frame(
  xmin = c(2.5, 0.5),
  xmax = c(3.5, 1.5),
  ymin = 0,
  ymax = y_max,
  fill = c(blue1, blue2)
)

p2 <- ggplot(df2) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    alpha = 0.7
  ) +
  scale_fill_identity() +
  scale_x_reverse(limits = x_limits) +
  coord_cartesian(ylim = c(0, y_max)) +
  ggtitle("Informed (Uniform)") +
  labs(y = "Density", x = "") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),  # removes minor gridlines
    panel.grid.major = element_line(color = "gray80", size = 0.2)  # keep major gridlines
  ) +
  geom_vline(xintercept = 3, color = "darkred", linetype = "dotted", size = 1) +
  geom_vline(xintercept = 1, color = "darkgreen", linetype = "dotted", size = 1) +
  guides(fill = "none")  # Remove the group legend

# -------- Graph 3: Normals at 3 and 1, SD = 0.1 --------
x_vals <- seq(0, 4.25, length.out = 1000)
df3 <- data.frame(
  x = rep(x_vals, 2),
  y = c(dnorm(x_vals, mean = 3, sd = 0.1), dnorm(x_vals, mean = 1, sd = 0.1)),
  group = rep(c("mean3", "mean1"), each = length(x_vals)),
  fill = rep(c(blue1, blue2), each = length(x_vals))
)

p3 <- ggplot(df3, aes(x = x, y = y, fill = fill)) +
  geom_area(data = subset(df3, group == "mean3"), fill = blue1, alpha = 0.6) +
  geom_area(data = subset(df3, group == "mean1"), fill = blue2, alpha = 0.6) +
  scale_x_reverse(limits = x_limits) +
  scale_y_continuous(limits = c(0, max(df3$y) * 1.1)) +
  ggtitle("Informed (Normal)") +
  labs(y = "Density", x = "") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),  # removes minor gridlines
    panel.grid.major = element_line(color = "gray80", size = 0.2)  # keep major gridlines
  ) +
  geom_vline(xintercept = 3, color = "darkred", linetype = "dotted", size = 1) +
  geom_vline(xintercept = 1, color = "darkgreen", linetype = "dotted", size = 1) +
  guides(fill = "none")  # Remove the group legend

# -------- Graph 4: Two uniforms from 3.5–2.5 and 2.5–1.5 --------
df4 <- data.frame(
  xmin = c(2.5, 1.5),
  xmax = c(3.5, 2.5),
  ymin = 0,
  ymax = y_max,
  fill = c(blue1, blue2)
)

p4 <- ggplot(df4) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    alpha = 0.7
  ) +
  scale_fill_identity() +
  scale_x_reverse(limits = x_limits) +
  coord_cartesian(ylim = c(0, y_max)) +
  ggtitle("Incorrect (Uniform)") +
  labs(y = "Density", x = "") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),  # removes minor gridlines
    panel.grid.major = element_line(color = "gray80", size = 0.2)  # keep major gridlines
  ) +
  geom_vline(xintercept = 3, color = "darkred", linetype = "dotted", size = 1) +
  geom_vline(xintercept = 1, color = "darkgreen", linetype = "dotted", size = 1) +
  guides(fill = "none")  # Remove the group legend

# -------- Graph 5: Normals at 3 and 2, SD = 0.1 --------
x_vals_2 <- seq(0, 4.25, length.out = 1000)
df5 <- data.frame(
  x = rep(x_vals_2, 2),
  y = c(dnorm(x_vals_2, mean = 3, sd = 0.1), dnorm(x_vals_2, mean = 2, sd = 0.1)),
  group = rep(c("mean3", "mean2"), each = length(x_vals_2)),
  fill = rep(c(blue1, blue2), each = length(x_vals_2))
)

p5 <- ggplot(df5, aes(x = x, y = y, fill = fill)) +
  geom_area(data = subset(df5, group == "mean3"), fill = blue1, alpha = 0.6) +
  geom_area(data = subset(df5, group == "mean2"), fill = blue2, alpha = 0.6) +
  scale_x_reverse(limits = x_limits) +
  scale_y_continuous(limits = c(0, max(df5$y) * 1.1)) +
  ggtitle("Incorrect (Normal)") +
  labs(y = "Density", x = "") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),  # removes minor gridlines
    panel.grid.major = element_line(color = "gray80", size = 0.2)  # keep major gridlines
  ) +
  geom_vline(xintercept = 3, color = "darkred", linetype = "dotted", size = 1) +
  geom_vline(xintercept = 1, color = "darkgreen", linetype = "dotted", size = 1) +
  guides(fill = "none")  # Remove the group legend

# -------- Arrange Plots --------
grid.arrange(p1, p2, p3, p4, p5, ncol = 1, bottom = "Years before present (millions)")


