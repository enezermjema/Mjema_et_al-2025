# Required packages

library(tidyverse)
library(patchwork)
library(readxl)
library(extrafont)

font_import(pattern = "Arial")
loadfonts(device = "pdf")

source("Desktop/git/functions/custom_theme.R")

# Climate information 
# These are daily measurements for different meteorological paramaters

spieke <- read_xlsx("Desktop/git/Chapter 2/data/open-meteo-53.81N7.80E0m.xlsx", 
                    skip = 3, col_names = TRUE)
brach <- read_xlsx("Desktop/git/Chapter 2/data/open-meteo-51.56N11.92E88m.xlsx", 
                   skip = 3, col_names = TRUE)

# Monthly means

spieke_month <- spieke %>%
  mutate(date_col = as.Date(time), 
         year = year(date_col),
         month = month(date_col),
         day = day(date_col)) %>%
  group_by(year, month) %>%
  summarise(
    precipitation_sum = sum(`precipitation_sum (mm)`, na.rm = TRUE),
    across(c(3:12), ~mean(.x, na.rm = TRUE)))

brach_month <- brach %>%
  mutate(date_col = as.Date(time), 
         year = year(date_col),
         month = month(date_col),
         day = day(date_col)) %>%
  group_by(year, month) %>%
  summarise(
    precipitation_sum = sum(`precipitation_sum (mm)`, na.rm = TRUE),
    across(c(3:12), ~mean(.x, na.rm = TRUE)))

# Adding seasons

spieke_month <- spieke_month %>%
  mutate(season = case_when(
    month %in% c(12, 1, 2) ~ "Winter",
    month %in% c(3, 4, 5) ~ "Spring",
    month %in% c(6, 7, 8) ~ "Summer",
    month %in% c(9, 10, 11) ~ "Autumn"
  ))

brach_month <- brach_month %>%
  mutate(season = case_when(
    month %in% c(12, 1, 2) ~ "Winter",
    month %in% c(3, 4, 5) ~ "Spring",
    month %in% c(6, 7, 8) ~ "Summer",
    month %in% c(9, 10, 11) ~ "Autumn"
  ))

# Simple line plot for both location across years
month_vec <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

cols <- c("2020" = "black", "2021" = "skyblue", "2022" = "blue", 
          "2023" = "green", "2024" = "orange", "2025" = "red")

## For spiekeroog

airtemp_spieke <- ggplot(spieke_month) +
  aes(x = month, y = `temperature_2m_mean (°C)`, group = year, color = as.factor(year)) +
  geom_line(size = .5) +
  geom_point() +
  scale_colour_manual(values = cols, name = "Year") +
  labs(title = NULL,
       x = NULL, #"Month",
       y = "Air temperature (°C)",
       color = "Year") +
  scale_x_continuous(breaks = c(1:12), labels = month_vec) +
  scale_y_continuous(limits = c(0, 20)) +
  theme_pub3()
soiltemp_spieke <- ggplot(spieke_month) +
  aes(x = month, y = `soil_temperature_0_to_7cm_mean (°C)`, group = year, color = as.factor(year)) +
  geom_line(size = .5) +
  geom_point() +
  scale_colour_manual(values = cols, name = "Year") +
  labs(title = NULL,
       x = NULL, #"Month",
       y = "Soil temperature (°C)",
       color = "Year") +
  scale_x_continuous(breaks = c(1:12), labels = month_vec) +
  scale_y_continuous(limits = c(0, 20)) +
  theme_pub3()
dewpoint_spieke <- ggplot(spieke_month) +
  aes(x = month, y = `dew_point_2m_mean (°C)`, group = year, color = as.factor(year)) +
  geom_line(size = .5) +
  geom_point() +
  scale_colour_manual(values = cols, name = "Year") +
  labs(title = NULL,
       x = "Month", 
       y = "Dew point temperature (°C)",
       color = "Year") +
  scale_x_continuous(breaks = c(1:12), labels = month_vec) +
  theme_pub3()
hum_spieke <- ggplot(spieke_month) +
  aes(x = month, y = `relative_humidity_2m_mean (%)`, group = year, color = as.factor(year)) +
  geom_line(size = .5) +
  geom_point() +
  scale_colour_manual(values = cols, name = "Year") +
  labs(title = NULL,
       x = "Month",
       y = "Relative humidity (%)",
       color = "Year") +
  scale_x_continuous(breaks = c(1:12), labels = month_vec) +
  scale_y_continuous(limits = c(74, 89)) +
  theme_pub3()
prec_spieke <- ggplot(spieke_month) +
  aes(x = month, y = precipitation_sum, group = year, color = as.factor(year)) +
  geom_line(size = .5) +
  geom_point() +
  scale_colour_manual(values = cols, name = "Year") +
  labs(title = NULL,
       x = "Month",
       y = "Precipitation (mm)",
       color = "Year") +
  scale_x_continuous(breaks = c(1:12), labels = month_vec) +
  #scale_y_continuous(limits = c(0, 0.33)) +
  theme_pub4()

weather_spieke <- soiltemp_spieke + airtemp_spieke + prec_spieke + 
  dewpoint_spieke + hum_spieke  +
  patchwork::plot_layout(nrow = 2) +
  plot_annotation(tag_levels = "a")

weather_spieke

## For Brachwitz

airtemp_brach <- ggplot(brach_month) +
  aes(x = month, y = `temperature_2m_mean (°C)`, group = year, color = as.factor(year)) +
  geom_line(size = .5) +
  geom_point() +
  scale_colour_manual(values = cols, name = "Year") +
  labs(title = NULL,
       x = NULL, #"Month",
       y = "Air temperature (°C)",
       color = "Year") +
  scale_x_continuous(breaks = c(1:12), labels = month_vec) +
  scale_y_continuous(limits = c(0, 21.8)) +
  theme_pub3()
soiltemp_brach <- ggplot(brach_month) +
  aes(x = month, y = `soil_temperature_0_to_7cm_mean (°C)`, group = year, color = as.factor(year)) +
  geom_line(size = .5) +
  geom_point() +
  scale_colour_manual(values = cols, name = "Year") +
  labs(title = NULL,
       x = "Month", #"Month",
       y = "Soil temperature (°C)",
       color = "Year") +
  scale_x_continuous(breaks = c(1:12), labels = month_vec) +
  scale_y_continuous(limits = c(0, 22.5)) +
  theme_pub3()
dewpoint_brach <- ggplot(brach_month) +
  aes(x = month, y = `dew_point_2m_mean (°C)`, group = year, color = as.factor(year)) +
  geom_line(size = .5) +
  geom_point() +
  scale_colour_manual(values = cols, name = "Year") +
  labs(title = NULL,
       x = NULL, 
       y = "Dew point temperature (°C)",
       color = "Year") +
  scale_x_continuous(breaks = c(1:12), labels = month_vec) +
  theme_pub3()
hum_brach <- ggplot(brach_month) +
  aes(x = month, y = `relative_humidity_2m_mean (%)`, group = year, color = as.factor(year)) +
  geom_line(size = .5) +
  geom_point() +
  scale_colour_manual(values = cols, name = "Year") +
  labs(title = NULL,
       x = "Month",
       y = "Relative humidity (%)",
       color = "Year") +
  scale_x_continuous(breaks = c(1:12), labels = month_vec) +
  scale_y_continuous(limits = c(57, 87)) +
  theme_pub3()
prec_brach <- ggplot(brach_month) +
  aes(x = month, y = precipitation_sum, group = year, color = as.factor(year)) +
  geom_line(size = .5) +
  geom_point() +
  scale_colour_manual(values = cols, name = "Year") +
  labs(title = NULL,
       x = "Month",
       y = "Precipitation (mm)",
       color = "Year") +
  scale_x_continuous(breaks = c(1:12), labels = month_vec) +
  #scale_y_continuous(limits = c(0, 0.16)) +
  theme_pub4()

weather_brach <- soiltemp_brach + airtemp_brach + prec_brach + 
  dewpoint_brach + hum_brach + 
  plot_layout(nrow = 2) +
  plot_annotation(tag_levels = "a")

weather_brach

# saving

ggsave(filename = "Desktop/git/Chapter 2/plots/weather - brachwitz.pdf", 
       plot = weather_brach, width = 12, height = 8, units = "in", dpi = 450)
ggsave(filename = "Desktop/git/Chapter 2/plots/weather - spiekeroog.pdf", 
       plot = weather_spieke, width = 12, height = 8, units = "in", dpi = 450)

# For RDA analysis 

## spring

spieke_spring <- spieke_month %>%
  filter(month %in% c(11, 12, 1, 2, 3, 4, 5)) %>%
  filter(!(year == 2020 & month == 1)) %>%
  filter(!(year == 2020 & month == 2)) %>%
  filter(!(year == 2020 & month == 3)) %>%
  filter(!(year == 2020 & month == 4)) %>%
  filter(!(year == 2020 & month == 5))

brach_spring <- brach_month %>%
  filter(month %in% c(11, 12, 1, 2, 3, 4, 5)) %>%
  filter(!(year == 2022 & month == 1)) %>%
  filter(!(year == 2022 & month == 2)) %>%
  filter(!(year == 2022 & month == 3)) %>%
  filter(!(year == 2022 & month == 4)) %>%
  filter(!(year == 2022 & month == 5))

# Winter

spieke_winter <- spieke_month %>%
  filter(year == 2021) %>%
  filter(month %in% c(7, 8, 9, 10, 11))

write.csv(spieke_winter, "Desktop/git/Chapter 2/resulting data/spieke-winter RDA.csv", row.names = FALSE)
write.csv(spieke_spring, "Desktop/git/Chapter 2/resulting data/spieke-spring RDA.csv", row.names = FALSE)
write.csv(brach_spring, "Desktop/git/Chapter 2/resulting data/brach-spring RDA.csv", row.names = FALSE)

write.csv(spieke_month, "Desktop/git/Chapter 2/resulting data/weather spiekeroog.csv", row.names = FALSE)
write.csv(brach_month, "Desktop/git/Chapter 2/resulting data/weather brachwitz.csv", row.names = FALSE)
