# Loading libraries

library(tidyverse)
library(patchwork)
library(extrafont)

font_import(pattern = "Arial")
loadfonts(device = "pdf")

source("Desktop/git/functions/custom_theme.R")

# mean monthly data

spieke <- read.csv("Desktop/git/Chapter 2/resulting data/weather spiekeroog.csv")
brach <- read.csv("Desktop/git/Chapter 2/resulting data/weather brachwitz.csv")

# Subsetting for relevant columns

spieke_df <- spieke[, c(1:4, 7, 9:12)]
brach_df <- brach[, c(1:4, 7, 9:12)]

# Computing means across years (baseline mean)

spieke_means <- spieke_df %>%
  group_by(month) %>%
  summarise(across(
    precipitation_sum:sunshine_duration..s., 
    mean, na.rm = TRUE), .groups = "drop")
brach_means <- brach_df %>%
  group_by(month) %>%
  summarise(across(
    precipitation_sum:sunshine_duration..s., 
    mean, na.rm = TRUE), .groups = "drop")

# Computing anomalies

spieke_df2 <- spieke_df %>%
  left_join(spieke_means, by = "month", suffix = c("", "_baseline"))
spieke_df3 <- spieke_df2 %>%
  mutate(across(precipitation_sum:sunshine_duration..s.,
                ~ .x - get(paste0(cur_column(), "_baseline")),
                .names = "{.col}_anomaly"))

brach_df2 <- brach_df %>%
  left_join(brach_means, by = "month", suffix = c("", "_baseline"))
brach_df3 <- brach_df2 %>%
  mutate(across(precipitation_sum:sunshine_duration..s.,
                ~ .x - get(paste0(cur_column(), "_baseline")),
                .names = "{.col}_anomaly"))

# Plotting the anomalies

month_vec <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
cols <- c("2021" = "skyblue", "2022" = "blue", 
          "2023" = "green", "2024" = "red2", "2025" = "orange")  # defining colors
cols_main <- c("2021" = "skyblue", "2022" = "grey80", 
               "2023" = "grey40", "2024" = "red2", "2025" = "grey60")  # defining colors

# Brachwitz

b1 <- ggplot(brach_df3, aes(x = month, y = temperature_2m_mean...C._anomaly, color = factor(year))) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = cols, name = "Year") +
  scale_x_continuous(breaks = c(1:12), labels = month_vec) +
  scale_y_continuous(n.breaks = 6) +
  labs(x = "Month", y = "Air temperature anomaly (°C)") +
  theme_pub3()
b2 <- ggplot(brach_df3, aes(x = month, y = soil_temperature_0_to_7cm_mean...C._anomaly, color = factor(year))) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = cols, name = "Year") +
  scale_x_continuous(breaks = c(1:12), labels = month_vec) +
  scale_y_continuous(n.breaks = 6) +  
  labs(x = "Month", y = "Soil temperature anomaly (°C)") +
  theme_pub3()
b3 <- ggplot(brach_df3, aes(x = month, y = relative_humidity_2m_mean...._anomaly, color = factor(year))) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = cols, name = "Year") +
  scale_x_continuous(breaks = c(1:12), labels = month_vec) +
  scale_y_continuous(n.breaks = 6) +  
  labs(x = "Month", y = "Relative humidity anomaly (%)") +
  theme_pub3()
b4 <- ggplot(brach_df3, aes(x = month, y = precipitation_sum_anomaly, color = factor(year))) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = cols, name = "Year") +
  scale_x_continuous(breaks = c(1:12), labels = month_vec) +
  scale_y_continuous(n.breaks = 6) +  
  labs(x = "Month", y = "Precipitation anomaly (mm)") +
  theme_pub3()
b5 <- ggplot(brach_df3, aes(x = month, y = dew_point_2m_mean...C._anomaly, color = factor(year))) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = cols, name = "Year") +
  scale_x_continuous(breaks = c(1:12), labels = month_vec) +
  #scale_y_continuous(n.breaks = 6) +  
  labs(x = "Month", y = "Dew point temperature anomaly (°C)") +
  theme_pub3()
b6 <- ggplot(brach_df3, aes(x = month, y = sunshine_duration..s._anomaly, color = factor(year))) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = cols, name = "Year") +
  scale_x_continuous(breaks = c(1:12), labels = month_vec) +
  #scale_y_continuous(n.breaks = 6) +  
  labs(x = "Month", y = "Sunshine duration anomaly (s)") +
  theme_pub4()
b7 <- ggplot(brach_df3, aes(x = month, y = wind_speed_10m_mean..km.h._anomaly, color = factor(year))) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = cols, name = "Year") +
  scale_x_continuous(breaks = c(1:12), labels = month_vec) +
  #scale_y_continuous(n.breaks = 6) +  
  labs(x = "Month", y = "Wind speed anomaly (km/h)") +
  theme_pub3()

# Spiekeroog
#line_types <- c("2021" = "solid", "2022" = "dotted", 
#               "2023" = "dashed", "2024" = "solid", 
#                "2025" = "dotdash") # defining line type

s1 <- ggplot(spieke_df3, aes(x = month,
                             y = temperature_2m_mean...C._anomaly,
                             color = factor(year))) +
  annotate(
    "rect", xmin = 1.5, xmax = 5.5, ymin = -Inf, 
    ymax = Inf, alpha = 0.2, fill = "darkgrey") +   # Highlight region between February (2) and May (5)
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = cols_main, name = "Year") +
  #scale_linetype_manual(values = line_types, name = "Year") +
  scale_x_continuous(breaks = 1:12, labels = month_vec) +
  scale_y_continuous(n.breaks = 6) +
  labs(x = "Month", y = "Air temperature anomaly (°C)") +
  theme_pub4()
s2 <- ggplot(spieke_df3, aes(x = month, y = soil_temperature_0_to_7cm_mean...C._anomaly, 
                             color = factor(year))) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = cols, name = "Year") +
  scale_x_continuous(breaks = c(1:12), labels = month_vec) +
  scale_y_continuous(n.breaks = 6) +  
  labs(x = "Month", y = "Soil temperature anomaly (°C)") +
  theme_pub3()
s3 <- ggplot(spieke_df3, aes(x = month, y = relative_humidity_2m_mean...._anomaly, color = factor(year))) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = cols, name = "Year") +
  scale_x_continuous(breaks = c(1:12), labels = month_vec) +
  scale_y_continuous(n.breaks = 6) +  
  labs(x = "Month", y = "Relative humidity anomaly (%)") +
  theme_pub3()
s4 <- ggplot(spieke_df3, aes(x = month, y = precipitation_sum_anomaly, color = factor(year))) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = cols, name = "Year") +
  scale_x_continuous(breaks = c(1:12), labels = month_vec) +
  scale_y_continuous(n.breaks = 6) +  
  labs(x = "Month", y = "Precipitation anomaly (mm)") +
  theme_pub3()
s5 <- ggplot(spieke_df3, aes(x = month, y = dew_point_2m_mean...C._anomaly, color = factor(year))) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = cols, name = "Year") +
  scale_x_continuous(breaks = c(1:12), labels = month_vec) +
  #scale_y_continuous(n.breaks = 6) +  
  labs(x = "Month", y = "Dew point temperature anomaly (°C)") +
  theme_pub3()
s6 <- ggplot(spieke_df3, aes(x = month, y = sunshine_duration..s._anomaly, color = factor(year))) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = cols, name = "Year") +
  scale_x_continuous(breaks = c(1:12), labels = month_vec) +
  #scale_y_continuous(n.breaks = 6) +  
  labs(x = "Month", y = "Sunshine duration anomaly (s)") +
  theme_pub3()
s7 <- ggplot(spieke_df3, aes(x = month, y = wind_speed_10m_mean..km.h._anomaly, color = factor(year))) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = cols, name = "Year") +
  scale_x_continuous(breaks = c(1:12), labels = month_vec) +
  #scale_y_continuous(n.breaks = 6) +  
  labs(x = "Month", y = "Wind speed anomaly (km/h)") +
  theme_pub4()

# joinging plots

anomalies_brach <- b1 + b2 + b3 + b4 + b5 + b6 + b7 +
  patchwork::plot_layout(nrow = 3) +
  plot_annotation(tag_levels = "a", title = "Anomalies (2022 - 2025)")

anomalies_spieke <- s2 + s3 + s5 + s4 + s5 + s7 + s6 +
  patchwork::plot_layout(nrow = 3) +
  plot_annotation(tag_levels = "a", title = "Anomalies (2021 - 2025)")

# saving plots

ggsave(filename = "Desktop/git/Chapter 2/plots/Figure-2a.pdf", 
       plot = s1, width = 7, height = 4, units = "in", dpi = 450)
ggsave(filename = "Desktop/git/Chapter 2/plots/anomaly - brachwitz.pdf", 
       plot = anomalies_brach, width = 12, height = 8, units = "in", dpi = 450)
ggsave(filename = "Desktop/git/Chapter 2/plots/anomaly - spiekeroog.pdf", 
       plot = anomalies_spieke, width = 12, height = 8, units = "in", dpi = 450)


# Saving files
write.csv(spieke_df3, "Desktop/git/Chapter 2/resulting data/Anomalies spiekeroog.csv", row.names = FALSE)
write.csv(brach_df3, "Desktop/git/Chapter 2/resulting data/Anomalies brachwitz.csv", row.names = FALSE)
