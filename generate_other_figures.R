#generate figures
library(ggplot2)

# read in the vital rate matrices from Pascarella and Horvitz 1998
source("vital_rate_matrices.R")
# source canopy damage matrices to test
source("damage_matrices.R")
# source canopy recovery matrices to test
source("recovery_matrices.R")
# source functions to create Markov chain for environmental states
source("mchain.R")
source("fit_gev.R") #fits GEV to historical data; provides GEV parameter estimates

# Main text Figure 2
# Plotting intensity vs frequency
# Multiple plot function from RCookbook http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
                     ncol = cols,
                     nrow = ceiling(numPlots / cols))
  }

  if (numPlots == 1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]],
            vp = viewport(
              layout.pos.row = matchidx$row,
              layout.pos.col = matchidx$col
            ))
    }
  }
}



pset <- c(2*p, 1.9*p, 1.8*p, 1.7*p, 1.6*p, 1.5*p, 1.4*p, 1.3*p, 1.2*p, 1.1*p, 1.05*p, 1.02*p, 1.01*p, p, 0.99*p, 0.98*p, 0.95*p, 0.9*p, 0.8*p, 0.7*p, 0.6*p, 0.5*p, 0.4*p, 0.3*p, 0.2*p, 0.1*p)

percent_change_vec <- c(seq(-0.02, 0.02, by = 0.01), seq(0.1, 2, by = 0.1))

#read in the files
#setwd(YOUR WORKING DIRECTORY)

mu <- 82.24106
sigma <- 35.87302
xi <- -0.3032049

all_data_both <- c()
for (ip in 1:length(pset)) {
  currp <- pset[ip]
  curr_data <-
    read.csv(
      paste0(
        "p_",
        currp,
        "_d_original_S_original_both.csv"
      )
    )
  curr_data[,1] <- rep("both", nrow(curr_data))
  colnames(curr_data) <- c("variables", "aest", "vest", "percent_change")
  new_mu <- mu * (1 + percent_change_vec)
  new_sigma <- sigma * (1 + percent_change_vec)
  curr_data <-
    cbind(curr_data, new_mu, new_sigma, p = rep(currp, nrow(curr_data)))
  all_data_both <- rbind(all_data_both, curr_data)
}

freq_of_interest <- c(0.5*p, 0.8*p, p, 1.2*p, 1.5*p)
data_of_interest_freq <- subset(all_data_both, p %in% freq_of_interest)
data_of_interest_freq$p <- as.factor(data_of_interest_freq$p)

f_plot <- ggplot(data_of_interest_freq, aes(100 * percent_change, aest, group = p)) +
  geom_line(aes(linetype = p), size = 0.5) +
  ggtitle("Change in Intensity") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5),
    legend.key = element_blank()) +
  xlab("Percent Change in GEV Parameters") +
  ylab(expression(paste("log" ~ lambda[S]))) +
  scale_linetype_discrete(name = "Frequency     ")

percent_of_interest <- c(0, 2)
data_of_interest_percent <- subset(all_data_both, percent_change %in% percent_of_interest)
data_of_interest_percent$percent_change <- as.factor(data_of_interest_percent$percent_change)

I_plot <- ggplot(data_of_interest_percent, aes(p, aest, group = percent_change)) +
  geom_line(aes(linetype = percent_change), size = 0.5) +
  ggtitle("Change in Frequency") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5),
    legend.key = element_blank()) +
  xlab("Frequency") +
  ylab(expression(paste("log" ~ lambda[S]))) +
  scale_linetype_discrete(name = "GEV % Change", labels = c("0", "200"))

tiff('growth_rate_plots.tiff', units = "in", width = 6, height = 7, res = 600)
multiplot(f_plot, I_plot, cols = 1)
dev.off()


#### Appendix S3 Figure S1 ####
N_list <- 1:500
upper_limit <- mu + (sigma / abs(xi))
x <- seq(u, upper_limit)
t_of_x <- (1 + ((x - mu) / sigma) * xi) ^ (-1 / xi)
pdf_x <- (1 / sigma) * (t_of_x ^ (xi + 1)) * exp(-t_of_x)
normalized_pdf_x <- pdf_x / sum(pdf_x)
plot(normalized_pdf_x ~ x, type = "l", lwd = 2, ylab = "Probability", xlab = "Wind Speed (kt)", main = "Probability Density Function for Hurricane Wind Speed", xlim = c(60, 250))

curr_percent <- 0.15
new_mu <- mu * (1 + curr_percent)
new_sigma <- sigma * (1 + curr_percent)
new_upper_limit <- new_mu + (new_sigma / abs(xi))
new_x <- seq(u, new_upper_limit)
new_t_of_x <- (1 + ((new_x - new_mu) / new_sigma) * xi) ^ (-1 / xi)
new_pdf_x <- (1 / new_sigma) * (new_t_of_x ^ (xi + 1)) * exp(-new_t_of_x)
new_normalized_pdf_x <- new_pdf_x / sum(new_pdf_x)
lines(new_normalized_pdf_x ~ new_x, lty = 2)
legend("topright", legend = c("Historical", "15% Increase"), lty = c(1, 2), bty = "n")

tiff('pdf_plot_with_shifted_GEV.tiff', units = "in", width = 6, height = 4, res = 600)
plot(normalized_pdf_x ~ x, type = "l", lwd = 2, ylab = "Probability", xlab = "Wind Speed (kt)", main = "Probability Density Function for Hurricane Wind Speed", xlim = c(60, 250))
lines(new_normalized_pdf_x ~ new_x, lty = 2)
legend("topright", legend = c("Historical", "15% Increase"), lty = c(1, 2), bty = "n")
dev.off()



#### Appendix S3 Figure S2 canopy states plot ####
d <- d_original

S <- S_original
percent_change_vec <- c(0, 0.5, 1, 1.5, 2)
nstate <- 7
bigP <- array(0, c(length(percent_change_vec), nstate, nstate))

xxvals <- read.csv("xx_values_for_mchain.csv", header = TRUE)
xxvals <- xxvals[, 2]
state_summary_stats <-
  matrix(nrow = length(percent_change_vec) * length(pset),
         ncol = 4)

for (ip in 1:length(pset)) {
  p <- pset[ip]
  for (j in 1:length(percent_change_vec)) {
    curr_percent <- percent_change_vec[j]
    new_mu <- mu * (1 + curr_percent)
    new_sigma <- sigma * (1 + curr_percent)
    upper_limit <- new_mu + (new_sigma / abs(xi))
    x <- seq(u, upper_limit)
    t_of_x <- (1 + xi * ((x - new_mu) / new_sigma)) ^ (-1 / xi)
    pdf_x <- (1 / new_sigma) * (t_of_x ^ (xi + 1)) * exp(-t_of_x)
    cdf_x <- exp(-(t_of_x))
    cdf_data <- as.data.frame(cbind(x, cdf_x))

    #discretize random draws into hurricane categories on the Saffir-Simpson scale
    category_threshold <-  c(64, 83, 96, 113, 137)
    w <- rep(0, length(category_threshold))
    for (icat in 1:length(category_threshold)) {
      if (icat == length(category_threshold)) {
        curr_lower <- category_threshold[icat]
        prob1 <- cdf_data$cdf_x[which(cdf_data$x == curr_lower)]
        prob2 <- 1
        w[icat] <- prob2 - prob1
      } else{
        curr_lower <- category_threshold[icat]
        curr_upper <- category_threshold[icat + 1]
        prob1 <- cdf_data$cdf_x[which(cdf_data$x == curr_lower)]
        prob2 <- cdf_data$cdf_x[which(cdf_data$x == curr_upper)]
        w[icat] <- prob2 - prob1
      }
    }
    w <- w / sum(w)

    h <- d %*% w * p

    pl <- h[1]
    pm <- h[2]
    ps <- h[3]

    nh <- 1 - sum(p)

    D <- matrix(
      c(p, p, ps+pm+.5*pl, ps+pm, ps+.5*pm, ps,    .5*ps,
        0,    0,    .5*pl,       .5*pl, .5*pm,    .5*pm, .5*ps,
        0,    0,    0,           .5*pl, .5*pl,    .5*pm, .5*pm,
        0,    0,    0,           0,     .5*pl,    .5*pl, .5*pm,
        0,    0,    0,           0,     0,        .5*pl, .5*pl,
        0,    0,    0,           0,     0,        0,     .5*pl,
        0,    0,    0,           0,     0,        0,     0),
      nrow = nstate,
      ncol = nstate,
      byrow = TRUE
    )
    D <- D / p

    bigP[j, , ] <- p * D + (1 - p) * S
  }

  nstate_vec <- dim(bigP)

  nstate <- nstate_vec[2]
  eq_vecs <- matrix(0, nrow = length(percent_change_vec), ncol = nstate)

  rho <- rep(0, length(percent_change_vec))

  dc <- matrix(0, nrow = length(percent_change_vec), ncol = nstate)

  for (i in 1:length(percent_change_vec)) {
    P <- drop(bigP[i, , ])
    vr <- eigen(P)$vectors
    dc1 <- diag(eigen(P)$values)
    d_sort <- sort(diag(dc1))
    #this is a sort in ascending order
    v_sort <- apply(vr, 1, rev) #this puts cols in right places
    v_eq <- v_sort[, nstate]

    v_eq <- v_eq / sum(v_eq)
    # in case the eigenvec is negative
    eq_vecs[i, ] <- t(v_eq) #store eqm vecs
    dc[i, ] <- t(d_sort)
  }
  rho <- abs(dc[, nstate - 1])


  #REQUIRES MCHAIN.R AND MAKEHURSET.R AND POPMATS2.R

  #input is time interval;
  #and initial pop vector; here set to an even vector
  ts <- 100000 #the length we use, can be changed

  #FOR each hurricane intensity
  aest <- matrix(0, nrow = 1, ncol = length(percent_change_vec))

  vest <- aest

  for (i_percent in 1:length(percent_change_vec)) {
    P_subset <- matrix(nrow = nstate, ncol = nstate, 1000)
    for (istate in 1:nstate) {
      P_subset[, istate] <-
        as.numeric(as.character(bigP[i_percent, , istate]))
    }
    P <- P_subset
    states <- mchain(t(P), ts, xxvals)

    state_percent <- matrix(nrow = length(states), ncol = 1)
    for (istate1 in 1:length(states)) {
      if (states[istate1] == 1) {
        state_percent[istate1] = 65
      }
      if (states[istate1] == 2) {
        state_percent[istate1] = 55
      }
      if (states[istate1] == 3) {
        state_percent[istate1] = 45
      }
      if (states[istate1] == 4) {
        state_percent[istate1] = 35
      }
      if (states[istate1] == 5) {
        state_percent[istate1] = 25
      }
      if (states[istate1] == 6) {
        state_percent[istate1] =  15
      }
      if (states[istate1] == 7) {
        state_percent[istate1] = 5
      }
    }
    state_summary_stats[(ip - 1) * length(percent_change_vec) + i_percent, ] <-
      c(mean(state_percent),
        var(state_percent),
        percent_change_vec[i_percent],
        p)
  }
}


state_summary_stats <- as.data.frame(state_summary_stats)
colnames(state_summary_stats) <-
  (c("state_mean", "state_var", "percent_change", "frequency"))
state_summary_stats$frequency <-
  as.factor(state_summary_stats$frequency)

data_of_interest_freq <-
  subset(state_summary_stats, frequency %in% freq_of_interest)
data_of_interest_freq$frequency <-
  as.factor(data_of_interest_freq$frequency)


mean_plot <- ggplot(data_of_interest_freq, aes(100 * percent_change, state_mean, group = frequency)) +
  geom_line(aes(lty = frequency), size = 0.5) +
  ggtitle("Mean Canopy State") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5),
    legend.key = element_blank()) +
  xlab("Percent Change in GEV Parameters") +
  ylab("Mean Canopy State") +
  scale_linetype_discrete(name = "Frequency     ")

var_plot <- ggplot(data_of_interest_freq, aes(100 * percent_change, sqrt(state_var), group = frequency)) +
  geom_line(aes(lty = frequency), size = 0.5) +
  ggtitle("Standard Deviation in Canopy State") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5),
    legend.key = element_blank()) +
  xlab("Percent Change in GEV Parameters") +
  ylab("Standard Deviation of Canopy State") +
  scale_linetype_discrete(name = "Frequency     ")


tiff('canopy_state_plots.tiff', units = "in", width = 6, height = 7, res = 600)
multiplot(mean_plot, var_plot, cols = 1)
dev.off()

#### S4 FIGURES 2 AND 3 ####

#### mu only ####
all_data_mu_only <- c()
for (ip in 1:length(pset)) {
  currp <- pset[ip]
  curr_data <-
    read.csv(
      paste0(
        "p_",
        currp,
        "_d_original_S_original_mu_only.csv"
      )
    )
  curr_data[,1] <- rep("mu_only", nrow(curr_data))
  colnames(curr_data) <- c("variables", "aest", "vest", "percent_change")
  new_mu <- mu * (1 - percent_change_vec)
  #new_sigma <- sigma * (1 + percent_change_vec)
  new_sigma <- rep(sigma, nrow(curr_data))
  curr_data <-
    cbind(curr_data, new_mu, new_sigma, p = rep(currp, nrow(curr_data)))
  all_data_mu_only <- rbind(all_data_mu_only, curr_data)
}
p <- 0.039
freq_of_interest <- c(0.5*p, 0.8*p, p, 1.2*p, 1.5*p)
data_of_interest_freq <- subset(all_data_mu_only, p %in% freq_of_interest)
data_of_interest_freq$p <- as.factor(data_of_interest_freq$p)

f_plot <- ggplot(data_of_interest_freq, aes(100 * percent_change, aest, group = p)) +
  geom_line(aes(linetype = p), size = 0.5) +
  ggtitle("Change in Intensity") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5),
    legend.key = element_blank()) +
  xlab("Percent Change in GEV Parameters") +
  ylab(expression(paste("log" ~ lambda[S]))) +
  scale_linetype_discrete(name = "Frequency     ")

percent_of_interest <- c(0, 2)
data_of_interest_percent <- subset(all_data_mu_only, percent_change %in% percent_of_interest)
data_of_interest_percent$percent_change <- as.factor(data_of_interest_percent$percent_change)

I_plot <- ggplot(data_of_interest_percent, aes(p, aest, group = percent_change)) +
  geom_line(aes(linetype = percent_change), size = 0.5) +
  ggtitle("Change in Frequency") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5),
    legend.key = element_blank()) +
  xlab("Frequency") +
  ylab(expression(paste("log" ~ lambda[S]))) +
  scale_linetype_discrete(name = "GEV % Change", labels = c("0", "200"))

tiff('growth_rate_plots_mu_only.tiff', units = "in", width = 6, height = 7, res = 600)
multiplot(f_plot, I_plot, cols = 1)
dev.off()


#### sigma only ####
percent_change_vec <- c(-0.02, -0.01, 0, 0.01, 0.02, seq(0.1, 2, by=0.1))

all_data_sigma_only <- c()
for (ip in 1:length(pset)) {
  currp <- pset[ip]
  curr_data <-
    read.csv(
      paste0(
        "p_",
        currp,
        "_d_original_S_original_sigma_only.csv"
      )
    )
  curr_data[,1] <- rep("sigma_only", nrow(curr_data))
  colnames(curr_data) <- c("variables", "aest", "vest", "percent_change")
  #new_mu <- mu * (1 - percent_change_vec)
  new_mu <- rep(mu, nrow(curr_data))
  new_sigma <- sigma * (1 + percent_change_vec)
  curr_data <-
    cbind(curr_data, new_mu, new_sigma, p = rep(currp, nrow(curr_data)))
  all_data_sigma_only <- rbind(all_data_sigma_only, curr_data)
}

freq_of_interest <- c(0.5*p, 0.8*p, p, 1.2*p, 1.5*p)
data_of_interest_freq <- subset(all_data_sigma_only, p %in% freq_of_interest)
data_of_interest_freq$p <- as.factor(data_of_interest_freq$p)

f_plot <- ggplot(data_of_interest_freq, aes(100 * percent_change, aest, group = p)) +
  geom_line(aes(linetype = p), size = 0.5) +
  ggtitle("Change in Intensity") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5),
    legend.key = element_blank()) +
  xlab("Percent Change in GEV Parameters") +
  ylab(expression(paste("log" ~ lambda[S]))) +
  scale_linetype_discrete(name = "Frequency     ")

percent_of_interest <- c(0, 2)
data_of_interest_percent <- subset(all_data_sigma_only, percent_change %in% percent_of_interest)
data_of_interest_percent$percent_change <- as.factor(data_of_interest_percent$percent_change)

I_plot <- ggplot(data_of_interest_percent, aes(p, aest, group = percent_change)) +
  geom_line(aes(linetype = percent_change), size = 0.5) +
  ggtitle("Change in Frequency") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5),
    legend.key = element_blank()) +
  xlab("Frequency") +
  ylab(expression(paste("log" ~ lambda[S]))) +
  scale_linetype_discrete(name = "GEV % Change", labels = c("0", "200"))

tiff('growth_rate_plots_sigma_only.tiff', units = "in", width = 6, height = 7, res = 600)
multiplot(f_plot, I_plot, cols = 1)
dev.off()
