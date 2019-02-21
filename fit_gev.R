# Fitting the GEV to historical data
library(extRemes) #provides GEV functions
# read in the data (data from NOAA's HURDAT 2 database v03.10)
hurr_data <- read.csv("storm_list.csv", header = TRUE)
colnames(hurr_data) <- c("serial_number", "count", "name", "max_sustained_wind_speed_kt", "year", "start_month", "end_month", "start_day", "end_day")
hurr_data <- hurr_data[hurr_data$year >= 1899, ]

# fit GEV to historic data
fit <- fevd(hurr_data$max_sustained_wind_speed_kt, type = "GEV", units = "kt")

# examine fit with qqplot
plot(fit, type = "primary")

# GEV parameter estimates
mu <- fit$results$par[1] #82.24106
sigma <- fit$results$par[2] #35.87302
xi <- fit$results$par[3] #-0.3032049

# historical hurricane frequency
p <- 0.039

# check (1+xi*y_i/sigma) > 0 for all y_i (Coles 2001)
1 + xi * ((hurr_data$max_sustained_wind_speed_kt - mu) / sigma) > 0

# Appendix S1 Figure S1 : make return level vs return period plot
tiff('return_level_plot.tiff', units = "in", width = 6, height = 4, res = 600)
plot(fit, type = "rl", main = "Return Level Plot")
legend("topleft", c("Historical Data", "GEV Model", "95% CI"), lty = c(0, 1, 2), pch = c(1, NA, NA), col = c("black", "black", "grey"), lwd = c(NA, 1, 2))
dev.off()

#Appendix S1 Figure S2: empirical cdf plot, checking the GEV fit
wind_breaks <- seq(32.5, 172.5, by = 5)
curr_hurr <- hurr_data[hurr_data$max_sustained_wind_speed_kt > 64, ]
pdf_hist <- hist(hurr_data$max_sustained_wind_speed_kt, plot = FALSE, breaks = wind_breaks)$counts
x <- seq(35, 265, by = 5)
cdf_x <- exp(-((1 + xi * ((x - mu) / sigma)) ^ (-1 / xi)))
tiff('cdf.tiff', units = "in", width = 6, height = 4, res = 600)
plot(cumsum(pdf_hist) / nrow(hurr_data) ~ seq(35, 170, by = 5), xlab = "Wind Speed (kt)", ylab = "Probability", ylim = c(0, 1), main = "Empirical vs model CDF")
lines(cdf_x ~ x)
dev.off()

# find category matrix
cat_list <- c("TS", "1", "2", "3", "4", "5")
category_threshold <- c(35, 64, 83, 96, 113, 137, max(hurr_data$max_sustained_wind_speed_kt))
category_mat <- hist(hurr_data$max_sustained_wind_speed_kt, breaks = category_threshold, freq = TRUE, xlab = "Wind Speed (kt)", ylab = "Number of Storms", main = "Historic Storm Counts 1899-2016", col = "grey")
u <- 64
category_mat$density/sum(category_mat$density)

hurr_only <- hurr_data[hurr_data$max_sustained_wind_speed_kt > 64, ]
hurr_only_hist <- hist(hurr_only$max_sustained_wind_speed_kt, breaks = category_threshold, plot = FALSE)
hurr_only_hist$density/sum(hurr_only_hist$density)
