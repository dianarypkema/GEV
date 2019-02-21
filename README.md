# GEV

Rypkema, D., Horvitz, C., and Tuljapurkar, S. 2019. How climate affects extreme events and hence ecological population models. Ecology.
________________________________________
DESCRIPTION

R code and data files for fitting the generalized extreme value distribution to historical hurricane data in our SE Florida study area, conducting Ardisia escallonioides population dynamics simulations, and generating figures.
________________________________________
AUTHORS

Diana Rypkema, Carol Horvitz, and Shripad Tuljapurkar

________________________________________
FILES

damage_matrices.R - historical canopy damage matrix from Pascarella and Horvitz (1998) and perturbed canopy damage matrices used for sensitivity analysis

fit_gev.R - fits GEV model to historic data (contained in storm_list.csv) providing GEV location (mu), scale (sigma), and shape (xi) parameters; generates Appendix S1 Figures S1 and S2 and category matrix

generate_other_figures.R - generates Figure 2 and Appendix S3 Figures S1 and S2

mchain.R - generates Markov chain of environmental states

population_simulations.R - simulates Ardisia escallonioides population dynamics for different combinations of hurricane frequency, intensity, canopy damage, and canopy recovery (bmat) scenarios; for increasing intensity simulations, changes both mu and sigma simultaneously (as discussed in main text)

population_simulations_with_mu_and_sigma_change.R - simulates Ardisia escallonioides population dynamics as in population_simulations.R but changes intensity by shifting only mu or only sigma (as in Appendix S4)

recovery_matrices.R - historical canopy recovery matrix (bmat) from Pascarella and Horvitz (1998) and perturbed recovery matrices for sensitivity analysis

sensitivity_calculations.R - calculates sensitivity of the stochastic growth rate with respect to perturbations in hurricane frequency, intensity and canopy damage and canopy recovery rate (as in Table 2)

vital_rate_matrices.R - vital rate matrices for Ardisia escallonioides in different canopy states from Pascarella and Horvitz (1998)

storm_list.csv - list of historical storms used to parameterize the GEV

xx_values_for_mchain.csv - randomly generated values used in DataS4.R
