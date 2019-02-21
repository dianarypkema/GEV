#### CANOPY DAMAGE LEVEL MATRICES ####
#historic damage level matrix from Pascarella and Horvitz 1998
#columns correspond to hurricane categories (1 through 5)
#rows correspond to low, medium, severe damage to the canopy
d_original <- matrix(c(0.8, 0.6, 0.4, 0.333, 0.1,
                           0.2, 0.3, 0.4, 0.333, 0.4,
                           0,   0.1, 0.2, 0.334, 0.5), nrow = 3, ncol = 5, byrow = TRUE)

#increase damage level by 0.01 to severe from low (0.01 from row 1 to row 3)
d_01 <- matrix(c(0.79, 0.59, 0.39, 0.323, 0.09,
                     0.2,  0.3,  0.4,  0.333, 0.4,
                     0.01, 0.11, 0.21, 0.344, 0.51), nrow = 3, ncol = 5, byrow = TRUE)

#increase damage level by 0.02 to severe from low
d_02=matrix(c(0.78, 0.58, 0.38, 0.313, 0.08,
                  0.20, 0.30, 0.40, 0.333, 0.40,
                  0.02, 0.12, 0.22, 0.354, 0.52), nrow = 3, ncol = 5, byrow=TRUE)

#create list of perturbed damage matrices to test
d_list <- c("d_original", "d_01", "d_02")