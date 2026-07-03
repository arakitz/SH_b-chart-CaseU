########################################################################
# OOC (Out-of-Control) Performance with Estimated Limits - Scenario 1
# 
# This script evaluates the performance of control charts for Beta 
# distributed processes when limits must be estimated from Phase I data.
# It compares three adjustment methods (A, B with epsilon=0, B with epsilon=0.2)
# against the plug-in approach and known-limits case.
########################################################################

# Load required packages
library(MASS)          # For statistical functions
library(gamlss.dist)   # For Beta distribution functions
library(gamlss)        # For GAMLSS modeling (estimation of distribution parameters)

########################################################################
# SECTION 1: Simulation Parameters
########################################################################

# True in-control parameter values for the Beta distribution
# Beta distribution is parameterized by mu (mean) and sigma (scale)
v01 <- 0.35  # In-control mu parameter
v02 <- 0.30  # In-control sigma parameter

# Simulation settings
sims <- 25000   # Number of simulation runs for estimating control limits
m <- 200        # Size of the Phase I sample (used for parameter estimation)

# Nominal False Alarm Rate (FAR) - corresponds to 3-sigma limits for normal
a <- 0.0027

# Adjusted FAR values (from tables in the referenced paper)
aa1 <- 0.00276  # Adjustment A (a')
aa2 <- 0.00101  # Adjustment B, epsilon = 0 (a'')
aa3 <- 0.00155  # Adjustment B, epsilon = 0.2 (a'')

########################################################################
# SECTION 2: Case K - Known Parameters (Benchmark)
########################################################################

# Calculate control limits when parameters are known (no estimation error)
UCL0 <- qBE(1 - a/2, mu = v01, sigma = v02)  # Upper Control Limit
LCL0 <- qBE(a/2, mu = v01, sigma = v02)      # Lower Control Limit
cat('The limits in Case K are ', 'LCL:', LCL0, ' and UCL:', UCL0, '\n')

# Define shifts in process parameters to evaluate out-of-control performance
d1s <- seq(0.8, 1.2, by = 0.1)  # Multiplicative shifts in mu
d2s <- d1s                       # Multiplicative shifts in sigma

########################################################################
# SECTION 3: Calculate ARL for Known Parameters (Case K)
########################################################################

listARL <- c()  # Empty vector to store results

# Loop over all combinations of shifts in mu and sigma
for(d1 in d1s) {
  for(d2 in d2s) {
    v11 <- d1 * v01  # Out-of-control mu
    v12 <- d2 * v02  # Out-of-control sigma
    
    # Probability of an observation falling outside control limits
    pout <- pBE(LCL0, mu = v11, sigma = v12) + 1 - pBE(UCL0, mu = v11, sigma = v12)
    
    # Average Run Length = 1 / probability of out-of-control signal
    arl <- 1 / pout
    
    # Store the shift parameters and resulting ARL
    listARL <- c(listARL, d1, d2, arl)
  }
}

# Convert to matrix for easier handling
M1 <- matrix(listARL, ncol = 3, byrow = TRUE)
# Columns: [d1 shift, d2 shift, ARL]

########################################################################
# SECTION 4: Adjustment Method A
########################################################################

# This method uses the adjusted FAR (aa1) to compensate for estimation error
limits1 <- c()  # Store estimated limits from each simulation

# Generate Phase I samples and estimate parameters for each run
for(j1 in 1:sims) {
  # Generate Phase I data from Beta distribution with known parameters
  X0 <- rBE(m, mu = v01, sigma = v02)
  
  # Fit Beta distribution using GAMLSS (Generalized Additive Models)
  fitX <- gamlss(X0 ~ 1, family = BE, trace = FALSE)
  
  # Extract estimated parameters (link function is logit)
  
  v1est <- 1 - 1/(1 + exp(fitX$mu.coefficients[[1]][1])) # mu hat
  v2est <- 1 - 1/(1 + exp(fitX$sigma.coefficients[[1]][1])) # sigma hat
  
  # Calculate control limits using adjusted FAR
  UCL1 <- qBE(1 - aa1/2, mu = v1est, sigma = v2est)
  LCL1 <- qBE(aa1/2, mu = v1est, sigma = v2est)
  
  # Store the limits
  limits1 <- c(limits1, LCL1, UCL1)
}

# Reshape limits into matrix format
conlim1 <- matrix(limits1, ncol = 2, byrow = TRUE)
# Columns: [LCL, UCL]
head(conlim1)  # Display first few rows for checking

# Pre-calculate shift combinations for efficiency
listd1 <- rep(d1s, each = length(d2s))
listd2 <- rep(d2s, times = length(d1s))

# Calculate ARL performance for each set of estimated limits
# Matrix to store ARL values: rows = shifts, columns = simulation runs
AARLa <- matrix(NA, nrow = length(listd1), ncol = sims)

for(jj1 in 1:sims) {
  LCL1 <- conlim1[jj1, 1]
  UCL1 <- conlim1[jj1, 2]
  
  for(kk1 in 1:length(listd2)) {
    v11 <- listd1[kk1] * v01
    v12 <- listd2[kk1] * v02
    
    pout <- pBE(LCL1, mu = v11, sigma = v12) + 1 - pBE(UCL1, mu = v11, sigma = v12)
    arlA <- 1 / pout
    
    AARLa[kk1, jj1] <- arlA
  }
}

# Calculate Average ARL (AARL) and Standard Deviation of ARL (SDARL)
# across all simulation runs for each shift combination
listARLadjA <- apply(AARLa, 1, mean)
listSDARLadjA <- apply(AARLa, 1, sd)

# Combine results into a matrix
M2 <- cbind(listd1, listd2, listARLadjA, listSDARLadjA)
# Columns: [d1, d2, AARL, SDARL]

########################################################################
# SECTION 5: Adjustment Method B with epsilon = 0
########################################################################

# This method uses adjusted FAR aa2 (epsilon = 0)
limits2 <- c()

for(j2 in 1:sims) {
  X0 <- rBE(m, mu = v01, sigma = v02)
  fitX <- gamlss(X0 ~ 1, family = BE, trace = FALSE)
  v1est <- 1 - 1/(1 + exp(fitX$mu.coefficients[[1]][1]))
  v2est <- 1 - 1/(1 + exp(fitX$sigma.coefficients[[1]][1]))
  
  UCL2 <- qBE(1 - aa2/2, mu = v1est, sigma = v2est)
  LCL2 <- qBE(aa2/2, mu = v1est, sigma = v2est)
  
  limits2 <- c(limits2, LCL2, UCL2)
}

conlim2 <- matrix(limits2, ncol = 2, byrow = TRUE)
head(conlim2)

# ARL performance evaluation for Method B (epsilon = 0)
AARLb <- matrix(NA, nrow = length(listd1), ncol = sims)

for(jj2 in 1:sims) {
  LCL2 <- conlim2[jj2, 1]
  UCL2 <- conlim2[jj2, 2]
  
  for(ii0 in 1:length(listd2)) {
    v11 <- listd1[ii0] * v01
    v12 <- listd2[ii0] * v02
    
    pout <- pBE(LCL2, mu = v11, sigma = v12) + 1 - pBE(UCL2, mu = v11, sigma = v12)
    arlB <- 1 / pout
    
    AARLb[ii0, jj2] <- arlB
  }
}

listARLadjB <- apply(AARLb, 1, mean)
listSDARLadjB <- apply(AARLb, 1, sd)

M3 <- cbind(listd1, listd2, listARLadjB, listSDARLadjB)

########################################################################
# SECTION 6: Adjustment Method B with epsilon = 0.20
########################################################################

# This method uses adjusted FAR aa3 (epsilon = 0.20)
limits3 <- c()

for(j2 in 1:sims) {
  X0 <- rBE(m, mu = v01, sigma = v02)
  fitX <- gamlss(X0 ~ 1, family = BE, trace = FALSE)
  v1est <- 1 - 1/(1 + exp(fitX$mu.coefficients[[1]][1]))
  v2est <- 1 - 1/(1 + exp(fitX$sigma.coefficients[[1]][1]))
  
  UCL3 <- qBE(1 - aa3/2, mu = v1est, sigma = v2est)
  LCL3 <- qBE(aa3/2, mu = v1est, sigma = v2est)
  
  limits3 <- c(limits3, LCL3, UCL3)
}

conlim3 <- matrix(limits3, ncol = 2, byrow = TRUE)
head(conlim3)

# ARL performance evaluation for Method B (epsilon = 0.20)
AARL3 <- matrix(NA, nrow = length(listd1), ncol = sims)

for(jj2 in 1:sims) {
  LCL3 <- conlim3[jj2, 1]
  UCL3 <- conlim3[jj2, 2]
  
  for(ii0 in 1:length(listd1)) {
    v11 <- listd1[ii0] * v01
    v12 <- listd2[ii0] * v02
    
    pout <- pBE(LCL3, mu = v11, sigma = v12) + 1 - pBE(UCL3, mu = v11, sigma = v12)
    arl3 <- 1 / pout
    
    AARL3[ii0, jj2] <- arl3
  }
}

listARL3 <- apply(AARL3, 1, mean)
listSDARL3 <- apply(AARL3, 1, sd)

M4 <- cbind(listd1, listd2, listARL3, listSDARL3)

########################################################################
# SECTION 7: Plug-in Limits (Naive Approach)
########################################################################

# This method uses the nominal FAR directly on estimated parameters
# (no adjustment for estimation error)
limitsPI <- c()

for(j1 in 1:sims) {
  X0 <- rBE(m, mu = v01, sigma = v02)
  fitX <- gamlss(X0 ~ 1, family = BE, trace = FALSE)
  v1est <- 1 - 1/(1 + exp(fitX$mu.coefficients[[1]][1]))
  v2est <- 1 - 1/(1 + exp(fitX$sigma.coefficients[[1]][1]))
  
  UCL <- qBE(1 - a/2, mu = v1est, sigma = v2est)
  LCL <- qBE(a/2, mu = v1est, sigma = v2est)
  
  limitsPI <- c(limitsPI, LCL, UCL)
}

conlimPI <- matrix(limitsPI, ncol = 2, byrow = TRUE)
head(conlimPI)

# ARL performance evaluation for plug-in limits
AARLpi <- matrix(NA, nrow = length(listd1), ncol = sims)

for(jj1 in 1:sims) {
  LCL <- conlimPI[jj1, 1]
  UCL <- conlimPI[jj1, 2]
  
  for(ii0 in 1:length(listd2)) {
    v11 <- listd1[ii0] * v01
    v12 <- listd2[ii0] * v02
    
    pout <- pBE(LCL, mu = v11, sigma = v12) + 1 - pBE(UCL, mu = v11, sigma = v12)
    arlPI <- 1 / pout
    
    AARLpi[ii0, jj1] <- arlPI
  }
}

listARLpi <- apply(AARLpi, 1, mean)
listSDARLpi <- apply(AARLpi, 1, sd)

M5 <- cbind(listd1, listd2, listARLpi, listSDARLpi)

########################################################################
# SECTION 8: Compile Final Results
########################################################################

# Combine all results into a single data frame for easy export and analysis
ARL.performance.OOC <- data.frame(
  d1_shift = listd1,          # Shift in mu (multiplicative factor)
  d2_shift = listd2,          # Shift in sigma (multiplicative factor)
  ARL_Known = M1[, 3],        # ARL with known parameters (Case K)
  AARL_MethodA = M2[, 3],     # AARL for Adjustment Method A
  SDARL_MethodA = M2[, 4],    # SDARL for Adjustment Method A
  AARL_MethodB_eps0 = M3[, 3],# AARL for Method B with epsilon=0
  SDARL_MethodB_eps0 = M3[, 4],# SDARL for Method B with epsilon=0
  AARL_MethodB_eps02 = M4[, 3],# AARL for Method B with epsilon=0.2
  SDARL_MethodB_eps02 = M4[, 4],# SDARL for Method B with epsilon=0.2
  AARL_PlugIn = M5[, 3],      # AARL for plug-in approach
  SDARL_PlugIn = M5[, 4]      # SDARL for plug-in approach
)

# Display the results
print(ARL.performance.OOC)

########################################################################
# END OF SCRIPT
########################################################################