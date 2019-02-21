#add the extremes package
library(fExtremes)

#set the working directory
setwd("/Users/Diana/Box/DIANA/HURRICANES/ms/revised_ms/revised_code/")

#read in the stochasticity file -- we use the same stochastic perturbations each combination of parameters so we are measuring the effect of the parameter choices, not the influence of environmental stochasticity on the results
xxvals <- read.csv("xx_values_for_mchain.csv", header=TRUE)
xxvals <- xxvals[,2]

# read in the vital rate matrices (B) from Pascarella and Horvitz 1998
source("vital_rate_matrices.R")
# source canopy damage matrices (d) to test
source("damage_matrices.R")
# source canopy recovery matrices (S) to test
source("recovery_matrices.R")
# source functions to create Markov chain for environmental states
source("mchain.R")
source("fit_gev.R") #fits GEV to historical data; provides GEV parameter estimates

# set some constants
u <- 64 # minimum maximum sustained wind speed to be classified as a hurricane
nstage <- 8 # define number of stages in vital rate matrices
nstate <- 7 # number of environmental states

####SHIFTING THE ENVIRONMENTAL TRANSITION MATRIX P TO REFLECT POTENTIAL CLIMATE CHANGE SCENARIOS####
# matrix convention P(a,b) is prob a <-- b
# this means that sum_a [P(a,b)] = sum of each column = 1
# so e is a left eigenvector with eval 1

# set of hurricane frequencies to test
p <- 0.039
pset <- c(2*p, 1.9*p, 1.8*p, 1.7*p, 1.6*p, 1.5*p, 1.4*p, 1.3*p, 1.2*p, 1.1*p, 1.05*p, 1.02*p, 1.01*p, p, 0.99*p, 0.98*p, 0.95*p, 0.9*p, 0.8*p, 0.7*p, 0.6*p, 0.5*p, 0.4*p, 0.3*p, 0.2*p, 0.1*p)

# set of percent change in GEV parameters to test
percent_change_vec <- c(-0.02, -0.01, 0, 0.01, 0.02, seq(0.1, 2, by=0.1))


# # create megamatrix to hold all C matrices
# bigP <- array(0, c(length(percent_change_vec), nstate, nstate))
#
# # run simulations for combinations of d, S, p, and percent change in GEV mu and sigma
# for(id in 1:length(d_list)) { #go through all damage scenarios
#   d <- get(d_list[id])
#   for(iS in 1:length(S_list)) { #go through all recovery scenarios
#     S <- get(S_list[iS])
#     for(ip in 1:length(pset)) {
#       p <- pset[ip]
#       for(j in 1:length(percent_change_vec)){
#         curr_percent <- percent_change_vec[j]
#         #new_mu <- mu*(1 + curr_percent) #shift mu by percent of interest
#         new_sigma <- sigma*(1 + curr_percent) #shift sigma by percent of interest
#         new_mu <- mu
#         upper_limit <- new_mu + (new_sigma/abs(xi)) #calculate upper limit of shifted GEV
#         x <- seq(u, upper_limit) #x values above hurricane threshold to maximum wind speed
#         cdf_x <- exp(-(1 + xi*((x - new_mu)/new_sigma))^(-1/xi)) #calculate cumulative distribution function for corresponding x values
#         cdf_data <- as.data.frame(cbind(x, cdf_x)) #store in dataframe
#
#         #discretize random draws into hurricane categories on the Saffir-Simpson scale
#         category_threshold <-  c(64, 83, 96, 113, 137) #Saffir Simpson hurricane category thresholds (minimum)
#         w <- rep(0, length(category_threshold)) #create vector to store probabilities by category
#         for(icat in 1:length(category_threshold)) {
#           if(icat == length(category_threshold)) {
#             curr_lower <- category_threshold[icat]
#             prob1 <- cdf_data$cdf_x[which(cdf_data$x == curr_lower)]
#             prob2 <- 1
#             w[icat] <- prob2 - prob1
#           }else{
#             curr_lower <- category_threshold[icat]
#             curr_upper <- category_threshold[icat+1]
#             prob1 <- cdf_data$cdf_x[which(cdf_data$x == curr_lower)]
#             prob2 <- cdf_data$cdf_x[which(cdf_data$x == curr_upper)]
#             w[icat] <- prob2 - prob1
#           }}
#         w <- w/sum(w) #normalize to sum of all categories is 1; CDF for GEV includes values below hurricane threshold 64kt, but we are using intensity GIVEN a hurricane happens, so we remove all values below 64kt, thus causing the need for normalization
#         h <- d %*% w*p
#         ps <- h[3,1] #severe damage
#         pm <- h[2,1] #medium damage
#         pl <- h[1,1] #low damage
#
#         #environmental transition matrix if there IS a hurricane
#         D <- matrix(c(p, p, ps+pm+0.5*pl, ps+pm, ps+0.5*pm, ps, 0.5*ps,
#                       0, 0, 0.5*pl, 0.5*pl, 0.5*pm, 0.5*pm, 0.5*ps,
#                       0, 0, 0, 0.5*pl, 0.5*pl, 0.5*pm, 0.5*pm,
#                       0, 0, 0, 0, 0.5*pl, 0.5*pl, 0.5*pm,
#                       0, 0, 0, 0, 0, 0.5*pl, 0.5*pl,
#                       0, 0, 0, 0, 0, 0, 0.5*pl,
#                       0, 0, 0, 0, 0, 0, 0), nrow=nstate, ncol=nstate, byrow=TRUE)
#         D <- D/p #normalize so columns sum to 1
#
#         #combine D and S to form P; matrix models overall patch transition probabilities over time
#         bigP[j,,] <- p*D + (1 - p)*S
#       }
#
#       nstate_vec <- dim(bigP)
#       nstate <- nstate_vec[2]
#       eq_vecs <- matrix(0, nrow=length(percent_change_vec), ncol=nstate)
#       rho <- rep(0, length(percent_change_vec))
#       dc <- matrix(0, nrow=length(percent_change_vec), ncol=nstate)
#       for(i in 1:length(percent_change_vec)) { #run populations simulations for a given intensity
#         P <- drop(bigP[i,,])
#         vr <- eigen(P)$vectors #get right eigenvector
#         dc1 <- diag(eigen(P)$values)
#         d_sort <- sort(diag(dc1)) #sorts in ascending order
#         v_sort <- apply(vr, 1, rev) #puts columns in right places
#         v_eq <- v_sort[,nstate]
#         v_eq <- v_eq/sum(v_eq) # in case the eigenvector is negative
#         eq_vecs[i,] <- t(v_eq) #store equilibrium vectors
#         dc[i,] <- t(d_sort)
#       }
#       rho <- abs(dc[, nstate-1])
#
#       #input is time interval;
#       #and initial pop vector; here set to an even vector
#       ts <- 100000 #the length of time for the simulation, can be changed
#
#       #for each hurricane intensity
#       aest <- matrix(0, nrow=1, ncol=length(percent_change_vec)) #aest will hold the stochastic growth rate estimates
#       stdest <- aest #stdest will hold the standard deviation of the stochastic growth rate estimates
#       for (iPercent in 1:length(percent_change_vec)) {
#         P_subset <- matrix(nrow=nstate, ncol=nstate, 1000)
#         for(istate in 1:nstate) {
#           P_subset[,istate] <-  as.numeric(as.character(bigP[iPercent,,istate]))
#         }
#         P <- P_subset
#         states <- mchain(t(P), ts, xxvals) #list of environmental states based on Markov chain
#         #returns a vector of length (t2+1) states
#         #set up vectors that have the number of rows defined by the number of stages, set to one
#         vec1<-rep(1, nstage)
#         vec1<-vec1/sum(vec1)
#         vec2 <- vec1
#         #set up a new cumulative growth
#         mgrowth <- 0
#         mgrowth2 <- 0 #stores square of log growths
#         for(i in 1:ts) { #loop through by year
#           mat1 <- get(paste0("B", states[i]))
#           vec1 <- mat1 %*% vec1
#           growth1 <- sum(vec1)
#           vec1<- vec1/growth1
#           mgrowth <- mgrowth + log(growth1)
#           mgrowth2 <- mgrowth2 + log(growth1)*log(growth1)
#         }
#         aest[iPercent] <- mgrowth/ts
#         ve <- (mgrowth2 - aest[iPercent]*aest[iPercent])/ts
#         stdest[iPercent] <- 1.96*sqrt(ve/ts) #convert from variance to standard deviation
#         print(iPercent)
#         print(p)
#       }
#       #write results to file
#       write.csv(cbind(aest[1,], stdest[1,], percent_change_vec), paste0("p_", p, "_", d_list[id], "_", S_list[iS], "_sigma_only.csv"))
#     }
#   }
# }


# set of percent change in GEV parameters to test
percent_change_vec <- c(-0.02, -0.01, 0, 0.01, 0.02, seq(0.1, 2, by=0.1))

setwd("/Users/Diana/Desktop/test_June_28_2018/")
# create megamatrix to hold all C matrices
bigP <- array(0, c(length(percent_change_vec), nstate, nstate))

# run simulations for combinations of d, S, p, and percent change in GEV mu and sigma
for(id in 1:length(d_list)) { #go through all damage scenarios
  d <- get(d_list[id])
  for(iS in 1:length(S_list)) { #go through all recovery scenarios
    S <- get(S_list[iS])
    for(ip in 1:length(pset)) {
      p <- pset[ip]
      for(j in 1:length(percent_change_vec)){
        curr_percent <- percent_change_vec[j]
        new_mu <- mu*(1 + curr_percent) #shift mu by percent of interest
        new_sigma <- sigma*(1 + curr_percent) #shift sigma by percent of interest
        upper_limit <- new_mu + (new_sigma/abs(xi)) #calculate upper limit of shifted GEV
        x <- seq(u, upper_limit) #x values above hurricane threshold to maximum wind speed
        cdf_x <- exp(-(1 + xi*((x - new_mu)/new_sigma))^(-1/xi)) #calculate cumulative distribution function for corresponding x values
        cdf_data <- as.data.frame(cbind(x, cdf_x)) #store in dataframe

        #discretize random draws into hurricane categories on the Saffir-Simpson scale
        category_threshold <-  c(64, 83, 96, 113, 137) #Saffir Simpson hurricane category thresholds (minimum)
        w <- rep(0, length(category_threshold)) #create vector to store probabilities by category
        for(icat in 1:length(category_threshold)) {
          if(icat == length(category_threshold)) {
            curr_lower <- category_threshold[icat]
            prob1 <- cdf_data$cdf_x[which(cdf_data$x == curr_lower)]
            prob2 <- 1
            w[icat] <- prob2 - prob1
          }else{
            curr_lower <- category_threshold[icat]
            curr_upper <- category_threshold[icat+1]
            prob1 <- cdf_data$cdf_x[which(cdf_data$x == curr_lower)]
            prob2 <- cdf_data$cdf_x[which(cdf_data$x == curr_upper)]
            w[icat] <- prob2 - prob1
          }}
        w <- w/sum(w) #normalize to sum of all categories is 1; CDF for GEV includes values below hurricane threshold 64kt, but we are using intensity GIVEN a hurricane happens, so we remove all values below 64kt, thus causing the need for normalization
        h <- d %*% w*p
        ps <- h[3,1] #severe damage
        pm <- h[2,1] #medium damage
        pl <- h[1,1] #low damage

        #environmental transition matrix if there IS a hurricane
        D <- matrix(c(p, p, ps+pm+0.5*pl, ps+pm, ps+0.5*pm, ps, 0.5*ps,
                      0, 0, 0.5*pl, 0.5*pl, 0.5*pm, 0.5*pm, 0.5*ps,
                      0, 0, 0, 0.5*pl, 0.5*pl, 0.5*pm, 0.5*pm,
                      0, 0, 0, 0, 0.5*pl, 0.5*pl, 0.5*pm,
                      0, 0, 0, 0, 0, 0.5*pl, 0.5*pl,
                      0, 0, 0, 0, 0, 0, 0.5*pl,
                      0, 0, 0, 0, 0, 0, 0), nrow=nstate, ncol=nstate, byrow=TRUE)
        D <- D/p #normalize so columns sum to 1

        #combine D and S to form P; matrix models overall patch transition probabilities over time
        bigP[j,,] <- p*D + (1 - p)*S
      }

      nstate_vec <- dim(bigP)
      nstate <- nstate_vec[2]
      eq_vecs <- matrix(0, nrow=length(percent_change_vec), ncol=nstate)
      rho <- rep(0, length(percent_change_vec))
      dc <- matrix(0, nrow=length(percent_change_vec), ncol=nstate)
      for(i in 1:length(percent_change_vec)) { #run populations simulations for a given intensity
        P <- drop(bigP[i,,])
        vr <- eigen(P)$vectors #get right eigenvector
        dc1 <- diag(eigen(P)$values)
        d_sort <- sort(diag(dc1)) #sorts in ascending order
        v_sort <- apply(vr, 1, rev) #puts columns in right places
        v_eq <- v_sort[,nstate]
        v_eq <- v_eq/sum(v_eq) # in case the eigenvector is negative
        eq_vecs[i,] <- t(v_eq) #store equilibrium vectors
        dc[i,] <- t(d_sort)
      }
      rho <- abs(dc[, nstate-1])

      #input is time interval;
      #and initial pop vector; here set to an even vector
      ts <- 100000 #the length of time for the simulation, can be changed

      #for each hurricane intensity
      aest <- matrix(0, nrow=1, ncol=length(percent_change_vec)) #aest will hold the stochastic growth rate estimates
      stdest <- aest #stdest will hold the standard deviation of the stochastic growth rate estimates
      for (iPercent in 1:length(percent_change_vec)) {
        P_subset <- matrix(nrow=nstate, ncol=nstate, 1000)
        for(istate in 1:nstate) {
          P_subset[,istate] <-  as.numeric(as.character(bigP[iPercent,,istate]))
        }
        P <- P_subset
        states <- mchain(t(P), ts, xxvals) #list of environmental states based on Markov chain
        #returns a vector of length (t2+1) states
        #set up vectors that have the number of rows defined by the number of stages, set to one
        vec1<-rep(1, nstage)
        vec1<-vec1/sum(vec1)
        vec2 <- vec1
        #set up a new cumulative growth
        mgrowth <- 0
        mgrowth2 <- 0 #stores square of log growths
        for(i in 1:ts) { #loop through by year
          mat1 <- get(paste0("B", states[i]))
          vec1 <- mat1 %*% vec1
          growth1 <- sum(vec1)
          vec1<- vec1/growth1
          mgrowth <- mgrowth + log(growth1)
          mgrowth2 <- mgrowth2 + log(growth1)*log(growth1)
        }
        aest[iPercent] <- mgrowth/ts
        ve <- (mgrowth2 - aest[iPercent]*aest[iPercent])/ts
        stdest[iPercent] <- 1.96*sqrt(ve/ts) #convert from variance to standard deviation
        print(iPercent)
        print(p)
      }
      #write results to file
      write.csv(cbind(aest[1,], stdest[1,], percent_change_vec), paste0("p_", p, "_", d_list[id], "_", S_list[iS], "_both.csv"))
    }
  }
}



#add the extremes package
library(fExtremes)

#set the working directory
#setwd("YOUR WORKING DIRECTORY")

#read in the stochasticity file -- we use the same stochastic perturbations each combination of parameters so we are measuring the effect of the parameter choices, not the influence of environmental stochasticity on the results
xxvals <- read.csv("xx_values_for_mchain.csv", header=TRUE)
xxvals <- xxvals[,2]

# read in the vital rate matrices (B) from Pascarella and Horvitz 1998
source("vital_rate_matrices.R")
# source canopy damage matrices (d) to test
source("damage_matrices.R")
# source canopy recovery matrices (S) to test
source("recovery_matrices.R")
# source functions to create Markov chain for environmental states
source("mchain.R")
source("fit_gev.R") #fits GEV to historical data; provides GEV parameter estimates

# set some constants
u <- 64 # minimum maximum sustained wind speed to be classified as a hurricane
nstage <- 8 # define number of stages in vital rate matrices
nstate <- 7 # number of environmental states

####SHIFTING THE ENVIRONMENTAL TRANSITION MATRIX P TO REFLECT POTENTIAL CLIMATE CHANGE SCENARIOS####
# matrix convention P(a,b) is prob a <-- b
# this means that sum_a [P(a,b)] = sum of each column = 1
# so e is a left eigenvector with eval 1

# set of hurricane frequencies to test
p <- 0.039
pset <- c(2*p, 1.9*p, 1.8*p, 1.7*p, 1.6*p, 1.5*p, 1.4*p, 1.3*p, 1.2*p, 1.1*p, 1.05*p, 1.02*p, 1.01*p, p, 0.99*p, 0.98*p, 0.95*p, 0.9*p, 0.8*p, 0.7*p, 0.6*p, 0.5*p, 0.4*p, 0.3*p, 0.2*p, 0.1*p)

# set of percent change in GEV parameters to test
#percent_change_vec <- c(-0.01, 0, 0.01, 0.02, seq(0.1, 0.7, by=0.1))


# create megamatrix to hold all C matrices
bigP <- array(0, c(length(percent_change_vec), nstate, nstate))

# run simulations for combinations of d, S, p, and percent change in GEV mu and sigma
for(id in 1:length(d_list)) { #go through all damage scenarios
  d <- get(d_list[id])
  for(iS in 1:length(S_list)) { #go through all recovery scenarios
    S <- get(S_list[iS])
    for(ip in 1:length(pset)) {
      p <- pset[ip]
      for(j in 1:length(percent_change_vec)){
        curr_percent <- percent_change_vec[j]
        new_mu <- mu*(1 + curr_percent) #shift mu by percent of interest
        #new_sigma <- sigma*(1 + curr_percent) #shift sigma by percent of interest
        new_sigma <- sigma
        upper_limit <- new_mu + (new_sigma/abs(xi)) #calculate upper limit of shifted GEV
        x <- seq(u, upper_limit) #x values above hurricane threshold to maximum wind speed
        cdf_x <- exp(-(1 + xi*((x - new_mu)/new_sigma))^(-1/xi)) #calculate cumulative distribution function for corresponding x values
        cdf_data <- as.data.frame(cbind(x, cdf_x)) #store in dataframe

        #discretize random draws into hurricane categories on the Saffir-Simpson scale
        category_threshold <-  c(64, 83, 96, 113, 137) #Saffir Simpson hurricane category thresholds (minimum)
        w <- rep(0, length(category_threshold)) #create vector to store probabilities by category
        for(icat in 1:length(category_threshold)) {
          if(icat == length(category_threshold)) {
            curr_lower <- category_threshold[icat]
            prob1 <- cdf_data$cdf_x[which(cdf_data$x == curr_lower)]
            prob2 <- 1
            w[icat] <- prob2 - prob1
          }else{
            curr_lower <- category_threshold[icat]
            curr_upper <- category_threshold[icat+1]
            prob1 <- cdf_data$cdf_x[which(cdf_data$x == curr_lower)]
            prob2 <- cdf_data$cdf_x[which(cdf_data$x == curr_upper)]
            w[icat] <- prob2 - prob1
          }}
        w <- w/sum(w) #normalize to sum of all categories is 1; CDF for GEV includes values below hurricane threshold 64kt, but we are using intensity GIVEN a hurricane happens, so we remove all values below 64kt, thus causing the need for normalization
        h <- d %*% w*p
        ps <- h[3,1] #severe damage
        pm <- h[2,1] #medium damage
        pl <- h[1,1] #low damage

        #environmental transition matrix if there IS a hurricane
        D <- matrix(c(p, p, ps+pm+0.5*pl, ps+pm, ps+0.5*pm, ps, 0.5*ps,
                      0, 0, 0.5*pl, 0.5*pl, 0.5*pm, 0.5*pm, 0.5*ps,
                      0, 0, 0, 0.5*pl, 0.5*pl, 0.5*pm, 0.5*pm,
                      0, 0, 0, 0, 0.5*pl, 0.5*pl, 0.5*pm,
                      0, 0, 0, 0, 0, 0.5*pl, 0.5*pl,
                      0, 0, 0, 0, 0, 0, 0.5*pl,
                      0, 0, 0, 0, 0, 0, 0), nrow=nstate, ncol=nstate, byrow=TRUE)
        D <- D/p #normalize so columns sum to 1

        #combine D and S to form P; matrix models overall patch transition probabilities over time
        bigP[j,,] <- p*D + (1 - p)*S
      }

      nstate_vec <- dim(bigP)
      nstate <- nstate_vec[2]
      eq_vecs <- matrix(0, nrow=length(percent_change_vec), ncol=nstate)
      rho <- rep(0, length(percent_change_vec))
      dc <- matrix(0, nrow=length(percent_change_vec), ncol=nstate)
      for(i in 1:length(percent_change_vec)) { #run populations simulations for a given intensity
        P <- drop(bigP[i,,])
        vr <- eigen(P)$vectors #get right eigenvector
        dc1 <- diag(eigen(P)$values)
        d_sort <- sort(diag(dc1)) #sorts in ascending order
        v_sort <- apply(vr, 1, rev) #puts columns in right places
        v_eq <- v_sort[,nstate]
        v_eq <- v_eq/sum(v_eq) # in case the eigenvector is negative
        eq_vecs[i,] <- t(v_eq) #store equilibrium vectors
        dc[i,] <- t(d_sort)
      }
      rho <- abs(dc[, nstate-1])

      #input is time interval;
      #and initial pop vector; here set to an even vector
      ts <- 100000 #the length of time for the simulation, can be changed

      #for each hurricane intensity
      aest <- matrix(0, nrow=1, ncol=length(percent_change_vec)) #aest will hold the stochastic growth rate estimates
      stdest <- aest #stdest will hold the standard deviation of the stochastic growth rate estimates
      for (iPercent in 1:length(percent_change_vec)) {
        P_subset <- matrix(nrow=nstate, ncol=nstate, 1000)
        for(istate in 1:nstate) {
          P_subset[,istate] <-  as.numeric(as.character(bigP[iPercent,,istate]))
        }
        P <- P_subset
        states <- mchain(t(P), ts, xxvals) #list of environmental states based on Markov chain
        #returns a vector of length (t2+1) states
        #set up vectors that have the number of rows defined by the number of stages, set to one
        vec1<-rep(1, nstage)
        vec1<-vec1/sum(vec1)
        vec2 <- vec1
        #set up a new cumulative growth
        mgrowth <- 0
        mgrowth2 <- 0 #stores square of log growths
        for(i in 1:ts) { #loop through by year
          mat1 <- get(paste0("B", states[i]))
          vec1 <- mat1 %*% vec1
          growth1 <- sum(vec1)
          vec1<- vec1/growth1
          mgrowth <- mgrowth + log(growth1)
          mgrowth2 <- mgrowth2 + log(growth1)*log(growth1)
        }
        aest[iPercent] <- mgrowth/ts
        ve <- (mgrowth2 - aest[iPercent]*aest[iPercent])/ts
        stdest[iPercent] <- 1.96*sqrt(ve/ts) #convert from variance to standard deviation
        print(iPercent)
        print(p)
      }
      #write results to file
      write.csv(cbind(aest[1,], stdest[1,], percent_change_vec), paste0("p_", p, "_", d_list[id], "_", S_list[iS], "_mu_only.csv"))
    }
  }
}


