#### ESTIMATE SENSITIVITIES####
#setwd("YOUR WORKING DIRECTORY")

# read in growth rate values
file_names <- list.files()
p <- 0.039
pset_of_interest <- c(p*0.98, p*0.99, p, p*1.01, p*1.02)
pset_of_interest <- sort(pset_of_interest)
percent_change_vec <- c(-0.02, -0.01, 0, 0.01, 0.02, seq(0.1, 2, by=0.1))
S_levels <- c("down02", "down01", "original", "01", "02")
S_labels <- c("-0.02", "-0.01", "0", "0.01", "0.02")
d_levels <- c("original", "01", "02")

all_data <- c()
for(ifile in 1:length(file_names)) {
  curr_file_name <- file_names[ifile]
  curr_data <- read.csv(curr_file_name, header=TRUE)
  curr_data <- curr_data[,-1]
  curr_param <- strsplit(curr_file_name, "_")
  curr_freq <- curr_param[[1]][2]
  curr_freq <- as.numeric(as.character(curr_freq))
  curr_S <- curr_param[[1]][length(curr_param[[1]])-1]
  curr_d <- curr_param[[1]][length(curr_param[[1]])-3]
  curr_freq <- curr_freq[[1]][length(curr_freq[[1]])]
  curr_freq <- as.numeric(paste0(curr_freq, collapse=""))
  curr_data <- cbind(curr_data, rep(curr_freq, nrow(curr_data)), rep(curr_S, nrow(curr_data)), rep(curr_d, nrow(curr_data)))
  all_data <- rbind(all_data, curr_data)
  print(ifile)
}
colnames(all_data) <- c("aest", "vest", "GEV_percent_change", "frequency", "S", "d")


#create numerical levels for shifts in d and S
all_data$S <- ordered(all_data$S, levels=S_levels, labels=c(-0.02, -0.01, 0, 0.01, 0.02))
all_data <- all_data[is.na(all_data$S) == FALSE,]
all_data$d <- ordered(all_data$d, levels=d_levels, labels=c(0, 0.01, 0.02))
all_data <- all_data[is.na(all_data$d) == FALSE,]

#select only frequencies of interest
all_data <- all_data[which(all_data$frequency %in% pset_of_interest),]
#select only data with historical frequency
curr_data_hist <- subset(all_data, frequency == p)

#select GEV, frequency, and d = historical, S changing
curr_data_hist_S <- subset(curr_data_hist, GEV_percent_change == 0)
curr_data_hist_S <- subset(curr_data_hist_S, d == 0)
curr_data_hist_S$S <- as.numeric(as.character(curr_data_hist_S$S))
#find the slope
S_change.lm <- summary(lm(aest ~ S, data=curr_data_hist_S)); S_change.lm

#select GEV, frequency, S = historical, d changing
curr_data_hist_d <- subset(curr_data_hist, GEV_percent_change == 0)
curr_data_hist_d <- subset(curr_data_hist_d, S == 0)
curr_data_hist_d$d <- as.numeric(as.character(curr_data_hist_d$d))
#find the slope
d_change.lm <- summary(lm(aest ~ d, data=curr_data_hist_d)); d_change.lm

#select frequency, S, d = historical, GEV changing
curr_data_hist_GEV <- subset(curr_data_hist, S == 0)
curr_data_hist_GEV <- subset(curr_data_hist_GEV, d == 0)
curr_data_hist_GEV <- subset(curr_data_hist_GEV, GEV_percent_change %in% seq(-0.02, 0.02, by=0.01))
#find the slope
GEV_change.lm <- summary(lm(aest ~ GEV_percent_change, data=curr_data_hist_GEV)); GEV_change.lm

#select GEV, S, d = historical, frequency changing
curr_data_hist_freq <- subset(all_data, GEV_percent_change == 0)
curr_data_hist_freq <- subset(curr_data_hist_freq, d == 0)
curr_data_hist_freq <- subset(curr_data_hist_freq, S == 0)
#find the slope
freq_change.lm <- summary(lm(aest ~ seq(-0.02, 0.02, by=0.01), data=curr_data_hist_freq)); freq_change.lm

#### slope additivity comparison ####
# get stochastic growth rate estimate for historical parameters
hist_aest <- all_data$aest[all_data$S == 0 & all_data$frequency == p & all_data$GEV_percent_change == 0 & all_data$d == 0]
S_slope <- S_change.lm$coefficients[[2]]
d_slope <- d_change.lm$coefficients[[2]]
gev_slope <- GEV_change.lm$coefficients[[2]]
freq_slope <- freq_change.lm$coefficients[[2]]

#check to see if slopes are additive

#S and d
#increase 0.01
est <- hist_aest + (S_slope*0.01) + (d_slope*0.01)
sim_aest <- all_data$aest[all_data$S == 0.01 & all_data$frequency==p & all_data$GEV_percent_change == 0 & all_data$d == 0.01]
100*((est-sim_aest)/(est+sim_aest))

#increase 0.02
est <- hist_aest + (S_slope*0.02) + (d_slope*0.02)
sim_aest <- all_data$aest[all_data$S == 0.02 & all_data$frequency==p & all_data$GEV_percent_change == 0 & all_data$d == 0.02]
100*((est-sim_aest)/(est+sim_aest))


#S and frequency
#increase 0.01
est <- hist_aest + (S_slope*0.01) + (freq_slope*0.01)
sim_aest <- all_data$aest[all_data$S == 0.01 & all_data$frequency==p*1.01 & all_data$GEV_percent_change == 0 & all_data$d == 0]
100*((est-sim_aest)/(est+sim_aest))

#increase 0.02
est <- hist_aest + (S_slope*0.02) + (freq_slope*0.02)
sim_aest <- all_data$aest[all_data$S == 0.02 & all_data$frequency==p*1.02 & all_data$GEV_percent_change == 0 & all_data$d == 0]
100*((est-sim_aest)/(est+sim_aest))

#decrease 0.01
est <- hist_aest - (S_slope*0.01) - (freq_slope*0.01)
sim_aest <- all_data$aest[all_data$S == -0.01 & all_data$frequency==p*0.99 & all_data$GEV_percent_change == 0 & all_data$d == 0]
100*((est-sim_aest)/(est+sim_aest))

#decrease 0.02
est <- hist_aest - (S_slope*0.02) - (freq_slope*0.02)
sim_aest <- all_data$aest[all_data$S == -0.02 & all_data$frequency==p*0.98 & all_data$GEV_percent_change == 0 & all_data$d == 0]
100*((est-sim_aest)/(est+sim_aest))


#S and intensity
#increase 0.01
est <- hist_aest + (S_slope*0.01) + (gev_slope*0.01)
sim_aest <- all_data$aest[all_data$S == 0.01 & all_data$frequency==p & all_data$GEV_percent_change == 0.01 & all_data$d == 0]
100*((est-sim_aest)/(est+sim_aest))

#increase 0.02
est <- hist_aest + (S_slope*0.02) + (gev_slope*0.02)
sim_aest <- all_data$aest[all_data$S == 0.02 & all_data$frequency==p & all_data$GEV_percent_change == 0.02 & all_data$d == 0]
100*((est-sim_aest)/(est+sim_aest))

#decrease 0.01
est <- hist_aest - (S_slope*0.01) - (gev_slope*0.01)
sim_aest <- all_data$aest[all_data$S == -0.01 & all_data$frequency==p & all_data$GEV_percent_change == -0.01 & all_data$d == 0]
100*((est-sim_aest)/(est+sim_aest))

#decrease 0.02
est <- hist_aest - (S_slope*0.02) - (gev_slope*0.02)
sim_aest <- all_data$aest[all_data$S == -0.02 & all_data$frequency==p & all_data$GEV_percent_change == -0.02 & all_data$d == 0]
100*((est-sim_aest)/(est+sim_aest))


#d and frequency
#increase 0.01
est <- hist_aest + (d_slope*0.01) + (freq_slope*0.01)
sim_aest <- all_data$aest[all_data$S == 0 & all_data$frequency==p*1.01 & all_data$GEV_percent_change == 0 & all_data$d == 0.01]
100*((est-sim_aest)/(est+sim_aest))

#increase 0.02
est <- hist_aest + (d_slope*0.02) + (freq_slope*0.02)
sim_aest <- all_data$aest[all_data$S == 0 & all_data$frequency==p*1.02 & all_data$GEV_percent_change == 0 & all_data$d == 0.02]
100*((est-sim_aest)/(est+sim_aest))



#d and intensity
#increase 0.01
est <- hist_aest + (d_slope*0.01) + (gev_slope*0.01)
sim_aest <- all_data$aest[all_data$S == 0 & all_data$frequency==p & all_data$GEV_percent_change == 0.01 & all_data$d == 0.01]
100*((est-sim_aest)/(est+sim_aest))

#increase 0.02
est <- hist_aest + (d_slope*0.02) + (gev_slope*0.02)
sim_aest <- all_data$aest[all_data$S == 0 & all_data$frequency==p & all_data$GEV_percent_change == 0.02 & all_data$d == 0.02]
100*((est-sim_aest)/(est+sim_aest))


#frequency and intensity
#increase 0.01
est <- hist_aest + (freq_slope*0.01) + (gev_slope*0.01)
sim_aest <- all_data$aest[all_data$S == 0 & all_data$frequency==p*1.01 & all_data$GEV_percent_change == 0.01 & all_data$d == 0]
100*((est-sim_aest)/(est+sim_aest))

#increase 0.02
est <- hist_aest + (freq_slope*0.02) + (gev_slope*0.02)
sim_aest <- all_data$aest[all_data$S == 0 & all_data$frequency==p*1.02 & all_data$GEV_percent_change == 0.02 & all_data$d == 0]
100*((est-sim_aest)/(est+sim_aest))

#decrease 0.01
est <- hist_aest - (freq_slope*0.01) - (gev_slope*0.01)
sim_aest <- all_data$aest[all_data$S == 0 & all_data$frequency==p*0.99 & all_data$GEV_percent_change == -0.01 & all_data$d == 0]
100*((est-sim_aest)/(est+sim_aest))

#decrease 0.02
est <- hist_aest - (freq_slope*0.02) - (gev_slope*0.02)
sim_aest <- all_data$aest[all_data$S == 0 & all_data$frequency==p*0.98 & all_data$GEV_percent_change == -0.02 & all_data$d == 0]
100*((est-sim_aest)/(est+sim_aest))

#S, d and frequency
#increase 0.01
est <- hist_aest + (S_slope*0.01) + (d_slope*0.01) + (freq_slope*0.01)
sim_aest <- all_data$aest[all_data$S == 0.01 & all_data$frequency==p*1.01 & all_data$GEV_percent_change == 0 & all_data$d == 0.01]
100*((est-sim_aest)/(est+sim_aest))

#increase 0.02
est <- hist_aest + (S_slope*0.02) + (d_slope*0.02) + (freq_slope*0.02)
sim_aest <- all_data$aest[all_data$S == 0.02 & all_data$frequency==p*1.02 & all_data$GEV_percent_change == 0 & all_data$d == 0.02]
100*((est-sim_aest)/(est+sim_aest))


#S, d and intensity
#increase 0.01
est <- hist_aest + (S_slope*0.01) + (d_slope*0.01) + (gev_slope*0.01)
sim_aest <- all_data$aest[all_data$S == 0.01 & all_data$frequency==p & all_data$GEV_percent_change == 0.01 & all_data$d == 0.01]
100*((est-sim_aest)/(est+sim_aest))

#increase 0.02
est <- hist_aest + (S_slope*0.02) + (d_slope*0.02) + (gev_slope*0.02)
sim_aest <- all_data$aest[all_data$S == 0.02 & all_data$frequency==p & all_data$GEV_percent_change == 0.02 & all_data$d == 0.02]
100*((est-sim_aest)/(est+sim_aest))


#d, frequency and intensity
#increase 0.01
est <- hist_aest + (d_slope*0.01) + (freq_slope*0.01) + (gev_slope*0.01)
sim_aest <- all_data$aest[all_data$S == 0 & all_data$frequency==p*1.01 & all_data$GEV_percent_change == 0.01 & all_data$d == 0.01]
100*((est-sim_aest)/(est+sim_aest))

#increase 0.02
est <- hist_aest + (d_slope*0.02) + (freq_slope*0.02) + (gev_slope*0.02)
sim_aest <- all_data$aest[all_data$S == 0 & all_data$frequency==p*1.02 & all_data$GEV_percent_change == 0.02 & all_data$d == 0.02]
100*((est-sim_aest)/(est+sim_aest))


#S, frequency and intensity
#increase 0.01
est <- hist_aest + (S_slope*0.01) + (freq_slope*0.01) + (gev_slope*0.01)
sim_aest <- all_data$aest[all_data$S == 0.01 & all_data$frequency==p*1.01 & all_data$GEV_percent_change == 0.01 & all_data$d == 0]
100*((est-sim_aest)/(est+sim_aest))

#increase 0.02
est <- hist_aest + (S_slope*0.02) + (freq_slope*0.02) + (gev_slope*0.02)
sim_aest <- all_data$aest[all_data$S == 0.02 & all_data$frequency==p*1.02 & all_data$GEV_percent_change == 0.02 & all_data$d == 0]
100*((est-sim_aest)/(est+sim_aest))

#decrease 0.01
est <- hist_aest - (S_slope*0.01) - (freq_slope*0.01) - (gev_slope*0.01)
sim_aest <- all_data$aest[all_data$S == -0.01 & all_data$frequency==p*0.99 & all_data$GEV_percent_change == -0.01 & all_data$d == 0]
100*((est-sim_aest)/(est+sim_aest))

#decrease 0.02
est <- hist_aest - (S_slope*0.02) - (freq_slope*0.02) - (gev_slope*0.02)
sim_aest <- all_data$aest[all_data$S == -0.02 & all_data$frequency==p*0.98 & all_data$GEV_percent_change == -0.02 & all_data$d == 0]
100*((est-sim_aest)/(est+sim_aest))


#S, d, frequency and intensity
est <- hist_aest + (S_slope*0.01) + (d_slope*0.01) + (freq_slope*0.01) + (gev_slope*0.01)
sim_aest <- all_data$aest[all_data$S == 0.01 & all_data$frequency==p*1.01 & all_data$GEV_percent_change == 0.01 & all_data$d == 0.01]
100*((est-sim_aest)/(est+sim_aest))

#increase 0.02
est <- hist_aest + (S_slope*0.02) + (d_slope*0.02) + (freq_slope*0.02) + (gev_slope*0.02)
sim_aest <- all_data$aest[all_data$S == 0.02 & all_data$frequency==p*1.02 & all_data$GEV_percent_change == 0.02 & all_data$d == 0.02]
100*((est-sim_aest)/(est+sim_aest))
