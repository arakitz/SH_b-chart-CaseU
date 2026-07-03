########################################################################
# SH[B] Control Chart for Sydney Relative Humidity Data
# 
# This script analyzes monthly relative humidity data from Sydney,
# estimates control limits using Beta distribution, and compares
# different methods for constructing control limits when parameters
# must be estimated from Phase I data.
# 
# Methods compared:
# 1. Plug-in limits (naive approach)
# 2. Adjustment Method A
# 3. Adjustment Method B (epsilon = 0)
# 4. Adjustment Method B (epsilon = 0.2)
# 5. Bootstrap-based limits
########################################################################

########################################################################
# SECTION 1: Load Required Packages and Data
########################################################################

# Required packages for analysis
library(dplyr)          # For data manipulation
library(readxl)         # For reading Excel files
library(gamlss.dist)    # For Beta distribution functions
library(gamlss)         # For GAMLSS modeling

# Note: The commented code below was used for initial exploration
# with the rattle package's weather dataset
# library(rattle)
# data1 <- as.data.frame(weatherAUS)
# data1[data1$Location == "Sydney", ]

########################################################################
# SECTION 2: Data Preparation and Processing
########################################################################

# Read the relative humidity data from Excel file
# The file should be in the current working directory
RHsydney <- read_excel("RHsydney.xls")

# Convert to data frame for easier manipulation
RHsydney1 <- as.data.frame(RHsydney)

# Attach for direct column access (use with caution)
attach(RHsydney1)

# Calculate minimum relative humidity for each month
# Group by Year and Month, then find minimum RH value
min_humidity <- RHsydney1 %>%
  group_by(Year, Month) %>%
  summarise(
    min_relative_humidity = min(RH, na.rm = TRUE),
    .groups = "drop"
  )

# Display the monthly minimum humidity values
print(min_humidity, n = 61)

# Extract the minimum humidity values as the main data series
# X represents monthly minimum relative humidity
x <- min_humidity$min_relative_humidity

########################################################################
# SECTION 3: Create Time Series Plot (Figure 1)
########################################################################

# Plot the humidity time series with customized axes
plot(x[1:61], type = 'l', axes = FALSE, ylim = c(0, 60),
     ylab = 'Relative Humidity (%)', xlab = 'time')

# Custom x-axis with meaningful date labels
axis(1, c(1, 10, 20, 30, 40, 50, 60),
     labels = c("10-2010", '08-2011', '06-2012', '06-2013', 
                '04-2014', '02-2015', '12-2015'))
axis(2)  # Default y-axis

# Add vertical line to separate Phase I and Phase II periods
abline(v = 50, lty = 2, col = 'grey', lwd = 2)

# Add points to the line plot
points(x[1:61], pch = 16)

########################################################################
# SECTION 4: Descriptive Statistics and Parameter Estimation
########################################################################

# Summary statistics of the humidity data
summary(x)

# Use the first 50 observations as Phase I (calibration) sample
# Scale to [0,1] for Beta distribution (requires values between 0 and 1)
x0 <- x[1:50] / 100

# Fit Beta distribution to the Phase I sample using GAMLSS
# BE = Beta distribution with mu and sigma parameterization
A <- gamlss(x0 ~ 1, family = BE, trace = FALSE)

# Extract estimated parameters
# Note: GAMLSS uses logit link for mu and log link for sigma
# Back-transform to original scale
v1est <- 1 - 1 / (1 + exp(A$mu.coefficients[[1]][1]))    # Estimated mu
v2est <- 1 - 1 / (1 + exp(A$sigma.coefficients[[1]][1])) # Estimated sigma

# Round for cleaner output
mu.hat <- round(v1est, digits = 5)
sigma.hat <- round(v2est, digits = 5)

# Display estimated parameters
paste('The estimated mu is:', mu.hat, ' and the estimated sigma is', sigma.hat)

########################################################################
# SECTION 5: Calculate Control Limits Using Different Methods
########################################################################

# Nominal False Alarm Rate (FAR) for control chart
alpha <- 0.0027

#------------------------------------------------------------------------
# 5.1: Plug-in Limits (naive approach)
#------------------------------------------------------------------------
# Uses estimated parameters directly with nominal FAR
UCLest <- qBE(1 - alpha/2, mu = mu.hat, sigma = sigma.hat)
LCLest <- qBE(alpha/2, mu = mu.hat, sigma = sigma.hat)

#------------------------------------------------------------------------
# 5.2: Adjustment Method A
#------------------------------------------------------------------------
# Uses adjusted FAR (a') from the paper's tables
aa1 <- 0.00340
UCL1 <- qBE(1 - aa1/2, mu = mu.hat, sigma = sigma.hat)
LCL1 <- qBE(aa1/2, mu = mu.hat, sigma = sigma.hat)

#------------------------------------------------------------------------
# 5.3: Adjustment Method B with epsilon = 0
#------------------------------------------------------------------------
# Uses adjusted FAR (a'') for epsilon = 0 from the paper's tables
aa2 <- 0.00022
UCL2 <- qBE(1 - aa2/2, mu = mu.hat, sigma = sigma.hat)
LCL2 <- qBE(aa2/2, mu = mu.hat, sigma = sigma.hat)

#------------------------------------------------------------------------
# 5.4: Adjustment Method B with epsilon = 0.2
#------------------------------------------------------------------------
# Uses adjusted FAR (a'') for epsilon = 0.2 from the paper's tables
aa3 <- 0.00052
UCL3 <- qBE(1 - aa3/2, mu = mu.hat, sigma = sigma.hat)
LCL3 <- qBE(aa3/2, mu = mu.hat, sigma = sigma.hat)

#------------------------------------------------------------------------
# 5.5: Bootstrap-based Limits
#------------------------------------------------------------------------
# Use bootstrap resampling to estimate limits that account for 
# estimation uncertainty in Phase I parameters

simsB <- 25000  # Number of bootstrap samples

# Initialize matrix to store bootstrap results
# Columns: [estimated mu, estimated sigma, LCL, UCL]
B <- matrix(0, ncol = 4, nrow = simsB)

# Bootstrap loop: resample from Phase I data with replacement
for(j in 1:simsB) {
  # Resample from Phase I data
  xB <- sample(x0, size = length(x0), replace = TRUE)
  
  # Fit Beta distribution to bootstrap sample
  AB <- gamlss(xB ~ 1, family = BE, trace = FALSE)
  
  # Extract estimated parameters
  v1estB <- 1 - 1 / (1 + exp(AB$mu.coefficients[[1]][1]))
  v2estB <- 1 - 1 / (1 + exp(AB$sigma.coefficients[[1]][1]))
  
  # Store estimated parameters
  B[j, 1] <- v1estB
  B[j, 2] <- v2estB
  
  # Calculate control limits using bootstrap estimates
  UCLestB <- qBE(1 - alpha/2, mu = B[j, 1], sigma = B[j, 2])
  LCLestB <- qBE(alpha/2, mu = B[j, 1], sigma = B[j, 2])
  
  # Store limits
  B[j, 3] <- LCLestB
  B[j, 4] <- UCLestB
}

# Use percentiles of bootstrap distribution for control limits
# Lower limit: 2.5th percentile of bootstrap LCLs
LCLstarB <- quantile(B[, 3], c(0.025))
# Upper limit: 97.5th percentile of bootstrap UCLs
UCLstarB <- quantile(B[, 4], c(0.975))

########################################################################
# SECTION 6: Create Phase I/II Control Chart with All Limits
########################################################################

# Plot the full time series with all control limit methods overlaid
plot(x[1:61], type = 'l', axes = FALSE, ylim = c(0, 90),
     ylab = 'Relative Humidity (%)', xlab = 'time',
     main = expression(SH[B] * "-chart for Sydney RH data"))

# Custom x-axis labels
axis(1, c(1, 10, 20, 30, 40, 50, 60),
     labels = c("01-2010", '08-2011', '06-2012', '06-2013', 
                '04-2014', '02-2015', '12-2015'))
axis(2)  # Default y-axis

# Vertical line separating Phase I (first 50 observations) from Phase II
abline(v = 50, lty = 1, lwd = 1)

# Add data points
points(x[1:61], pch = 16, cex = 1.5)

# Add control limits for each method (scaled back to percentage)
# Plug-in limits (solid line)
abline(h = (100 * UCLest), lty = 1, lwd = 1)
abline(h = (100 * LCLest), lty = 1, lwd = 1)

# Method A (dashed line)
abline(h = (100 * UCL1), lty = 2, lwd = 1)
abline(h = (100 * LCL1), lty = 2, lwd = 1)

# Method B, eps=0 (dotted line)
abline(h = (100 * UCL2), lty = 3, lwd = 1)
abline(h = (100 * LCL2), lty = 3, lwd = 1)

# Method B, eps=0.2 (dot-dash line)
abline(h = (100 * UCL3), lty = 4, lwd = 1)
abline(h = (100 * LCL3), lty = 4, lwd = 1)

# Bootstrap limits (long-dash line)
abline(h = (100 * UCLstarB), lty = 5, lwd = 1)
abline(h = (100 * LCLstarB), lty = 5, lwd = 1)

# Add legend identifying all methods
legend('topright', 
       c('Plug-in', 'AdjA', 'AdjB, eps=0', 'AdjB, eps=0.2', 'Bootstrap'),
       lty = 1:5, lwd = c(1, 1, 1, 1, 1), cex = 0.75)

# Add text labels for Phase I and Phase II periods
text(42, 70, "Phase I")
text(53, 70, "Phase II")

########################################################################
# SECTION 7: Out-of-Control Scenario - Shift in Mu
########################################################################

# Define OOC data representing a shift in the mu parameter
# (These values correspond to a process where the mean has shifted)
y <- c(0.30609619, 0.25720683, 0.29109522,
       0.37668141, 0.50084142, 0.35030771,
       0.38610248, 0.31420693, 0.52388937,
       0.21405942, 0.17209139, 0.14343708,
       0.31683782, 0.11597719, 0.15894399,
       0.31830692, 0.33706144, 0.16367507,
       0.27481577, 0.15476869, 0.53893121,
       0.31816739, 0.18837781, 0.40100539,
       0.38775144, 0.36651343, 0.19352635,
       0.20936857, 0.33370933, 0.22000786,
       0.18772538, 0.56640794, 0.16146448,
       0.33519866, 0.44298871, 0.12165616,
       0.47767826, 0.35177966, 0.22916423,
       0.06598294, 0.10970764, 0.42896833,
       0.28774871, 0.34526811, 0.09874269,
       0.17692411, 0.35719114, 0.30169661,
       0.24562682, 0.57983994, 0.32418274,
       0.11011865, 0.37507631, 0.44071340,
       0.16974362, 0.21863353, 0.15890170,
       0.16431547, 0.32858313, 0.30379717,
       0.22211175, 0.59030055, 0.24904023,
       0.35161064, 0.34471037, 0.17930245,
       0.39895906, 0.17159297, 0.18355666,
       0.21480203)

# Create control chart for OOC scenario with shift in mu
plot(y, type = 'l', axes = FALSE, ylim = c(0, 0.9),
     xlim = c(0, 70), ylab = '', xlab = 'time',
     main = expression(SH[B] * "-chart for OOC data, shift in " * mu))

# Custom axes
axis(1, c(0, 10, 20, 30, 40, 50, 60, 70))
axis(2, c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))

# Add data points
points(y, pch = 16, cex = 1.5)

# Overlay all control limit methods (using original scale, not percentages)
# Plug-in limits (solid)
abline(h = (UCLest), lty = 1, lwd = 1)
abline(h = (LCLest), lty = 1, lwd = 1)

# Method A (dashed)
abline(h = (UCL1), lty = 2, lwd = 1)
abline(h = (LCL1), lty = 2, lwd = 1)

# Method B, eps=0 (dotted)
abline(h = (UCL2), lty = 3, lwd = 1)
abline(h = (LCL2), lty = 3, lwd = 1)

# Method B, eps=0.2 (dot-dash)
abline(h = (UCL3), lty = 4, lwd = 1)
abline(h = (LCL3), lty = 4, lwd = 1)

# Bootstrap limits (long-dash)
abline(h = (UCLstarB), lty = 5, lwd = 1)
abline(h = (LCLstarB), lty = 5, lwd = 1)

# Add legend
legend('topright', 
       c('Plug-in', 'AdjA', 'AdjB, eps=0', 'AdjB, eps=0.2', 'Bootstrap'),
       lty = 1:5, lwd = c(1, 1, 1, 1, 1), cex = 0.75)

# Add vertical line separating IC and OOC periods
abline(v = 20, lty = 2, col = 'darkgrey')

# Add period labels
text(12, 0.75, "IC period")
text(28, 0.75, "OOC period")

########################################################################
# SECTION 8: Out-of-Control Scenario - Shift in Sigma
########################################################################

# Define OOC data representing a shift in the sigma parameter
# (These values correspond to a process where the variability has changed)
y <- c(0.40890923, 0.20711377, 0.48565290,
       0.41342518, 0.14514276, 0.28553799,
       0.36154274, 0.17148259, 0.25833854,
       0.44677121, 0.20605114, 0.20945256,
       0.33630733, 0.33168441, 0.35256853,
       0.47304573, 0.17692752, 0.19831857,
       0.23240015, 0.23455703, 0.13494582,
       0.21935793, 0.25749392, 0.26100636,
       0.09398065, 0.22532295, 0.35740527,
       0.37906132, 0.17381027, 0.10006924,
       0.18711972, 0.24825863, 0.31336222,
       0.16107138, 0.24012080, 0.40356519,
       0.16573125, 0.48085574, 0.13661168,
       0.17103081, 0.22010640, 0.50253456,
       0.16494011, 0.60033655, 0.30104390,
       0.22931254, 0.59511069, 0.18845682,
       0.34120781, 0.56904273, 0.20563767,
       0.64621361, 0.24440051, 0.20950493,
       0.12702669, 0.18782709, 0.40409009,
       0.15381170, 0.33829964, 0.26553242,
       0.11436794, 0.32380358, 0.23113463,
       0.20722643, 0.47252334, 0.26407311,
       0.19531615, 0.43436949, 0.04664434,
       0.16630466)

# Create control chart for OOC scenario with shift in sigma
plot(y, type = 'l', axes = FALSE, ylim = c(0, 0.9),
     xlim = c(0, 70), ylab = '', xlab = 'time',
     main = expression(SH[B] * "-chart for OOC data, shift in " * sigma))

# Custom axes
axis(1, c(0, 10, 20, 30, 40, 50, 60, 70))
axis(2, c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))

# Add data points
points(y, pch = 16, cex = 1.5)

# Overlay all control limit methods
# Plug-in limits (solid)
abline(h = (UCLest), lty = 1, lwd = 1)
abline(h = (LCLest), lty = 1, lwd = 1)

# Method A (dashed)
abline(h = (UCL1), lty = 2, lwd = 1)
abline(h = (LCL1), lty = 2, lwd = 1)

# Method B, eps=0 (dotted)
abline(h = (UCL2), lty = 3, lwd = 1)
abline(h = (LCL2), lty = 3, lwd = 1)

# Method B, eps=0.2 (dot-dash)
abline(h = (UCL3), lty = 4, lwd = 1)
abline(h = (LCL3), lty = 4, lwd = 1)

# Bootstrap limits (long-dash)
abline(h = (UCLstarB), lty = 5, lwd = 1)
abline(h = (LCLstarB), lty = 5, lwd = 1)

# Add legend
legend('topright', 
       c('Plug-in', 'AdjA', 'AdjB, eps=0', 'AdjB, eps=0.2', 'Bootstrap'),
       lty = 1:5, lwd = c(1, 1, 1, 1, 1), cex = 0.75)

# Add vertical line separating IC and OOC periods
abline(v = 20, lty = 2, col = 'darkgrey')

# Add period labels
text(12, 0.75, "IC period")
text(28, 0.75, "OOC period")

########################################################################
# END OF SCRIPT
########################################################################
