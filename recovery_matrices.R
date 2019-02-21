#### RECOVERY DYNAMICS MATRICES ####
#historic canopy closure (succession) S matrix (called B matrix in Pascarella and Horvitz 1998)
#environmental transition matrix if there is NOT a hurricane

nstate <- 7 # number of environmental states of canopy openness

S_original=matrix(c(0.75, 0, 0,	0, 0, 0, 0,
                      0.25, 0.75,	0, 0,	0, 0,	0,
                      0, 0.25, 0.5, 0, 0,	0, 0,
                      0, 0, 0.5, 0.25, 0,	0, 0,
                      0, 0, 0, 0.75, 0.25, 0, 0,
                      0, 0, 0, 0, 0.75, 0.25, 0,
                      0, 0, 0, 0, 0, 0.75, 1), nrow=nstate, ncol=nstate, byrow=TRUE)

#lengthen recovery time by increasing probability of remaining in a given state by 0.01 (and decreasing probability of transitioning to a more closed canopy state by 0.01)
S_01=matrix(c(0.76, 0, 0,	0, 0,	0, 0,
                0.24, 0.76,	0, 0,	0, 0, 0,
                0, 0.24, 0.51, 0,	0, 0,	0,
                0, 0, 0.49, 0.26, 0, 0,	0,
                0, 0, 0, 0.74, 0.26, 0, 0,
                0, 0, 0, 0, 0.74, 0.26, 0,
                0, 0, 0, 0, 0, 0.74, 1), nrow=nstate, ncol=nstate, byrow=TRUE)

#lengthen recovery time by increasing probability of remaining in a given state by 0.02
S_02=matrix(c(0.77, 0, 0,	0, 0,	0, 0,
                0.23, 0.77,	0, 0,	0, 0,	0,
                0, 0.23, 0.52, 0,	0, 0,	0,
                0, 0, 0.48, 0.27, 0, 0,	0,
                0, 0, 0, 0.73, 0.27, 0, 0,
                0, 0, 0, 0, 0.73, 0.27, 0,
                0, 0, 0, 0, 0, 0.73, 1), nrow=nstate, ncol=nstate, byrow=TRUE)

#decrease recovery time by decreasing probability of remaining in a given state by 0.01
S_down01=matrix(c(0.74, 0, 0,	0, 0,	0, 0,
                    0.26, 0.74,	0, 0,	0, 0,	0,
                    0, 0.26, 0.49, 0,	0, 0,	0,
                    0, 0, 0.51, 0.24, 0, 0,	0,
                    0, 0, 0, 0.76, 0.24, 0, 0,
                    0, 0, 0, 0, 0.76, 0.24, 0,
                    0, 0, 0, 0, 0, 0.76, 1), nrow=nstate, ncol=nstate, byrow=TRUE)

#decrease recovery time by decreasing probability of remaining in a given state by 0.02
S_down02=matrix(c(0.73, 0, 0,	0, 0, 0, 0,
                    0.27, 0.73,	0, 0,	0, 0,	0,
                    0, 0.27, 0.48, 0,	0, 0,	0,
                    0, 0, 0.52, 0.23, 0, 0,	0,
                    0, 0, 0, 0.77, 0.23, 0, 0,
                    0, 0, 0, 0, 0.77, 0.23, 0,
                    0, 0, 0, 0, 0, 0.77, 1), nrow=nstate, ncol=nstate, byrow=TRUE)

#create list of perturbed B matrices to test
S_list <- c("S_original", "S_01", "S_02", "S_down01", "S_down02")