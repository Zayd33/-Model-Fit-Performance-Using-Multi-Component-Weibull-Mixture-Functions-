# Load required libraries
library(MASS)
library(survival)
library(fitdistrplus)
library(stats4)
library(MASS)
library(ggplot2)

# Load data
data <- read.csv("lag.vtc.2dp.csv")
wind_speed <- na.omit(data$Windspeed)

# Define modified bimodal Weibull with shift (delta)
bimodal_weibull_shifted <- function(x, alpha1, beta1, alpha2, beta2, p, delta) {
  x_shifted <- x - delta
  x_shifted[x_shifted <= 0] <- 1e-6  # to avoid log(0) or negative powers
  p * dweibull(x_shifted, shape=beta1, scale=alpha1) +
    (1 - p) * dweibull(x_shifted, shape=beta2, scale=alpha2)
}

# Negative log-likelihood with delta
nll_shifted <- function(alpha1, beta1, alpha2, beta2, p, delta) {
  if (any(c(alpha1, beta1, alpha2, beta2, p, delta) <= 0) || p >= 1) return(1e6)
  -sum(log(bimodal_weibull_shifted(wind_speed, alpha1, beta1, alpha2, beta2, p, delta)))
}

# Fit model with shift parameter
start_vals <- list(alpha1=2, beta1=2, alpha2=5, beta2=2, p=0.5, delta=0.5)
fit_shift <- mle(nll_shifted, start=start_vals, method="L-BFGS-B",
                 lower=c(0.1, 0.1, 0.1, 0.1, 0.01, 0.01), 
                 upper=c(20, 10, 20, 10, 0.99, 3))

# Extract parameters
params <- coef(fit_shift)
alpha1 <- params["alpha1"]
beta1 <- params["beta1"]
alpha2 <- params["alpha2"]
beta2 <- params["beta2"]
p <- params["p"]
delta <- params["delta"]

# Histogram + fitted curve
hist(wind_speed, probability=TRUE, breaks=20, col="lightgreen",
     main="Modified Mixture Weibull Fit - Lag2dp", xlab="Wind Speed (m/s)")
x_seq <- seq(min(wind_speed), max(wind_speed), length.out=100)
lines(x_seq, bimodal_weibull_shifted(x_seq, alpha1, beta1, alpha2, beta2, p, delta),
      col="red", lwd=2)

# K-S Test
ks_stat <- ks.test(wind_speed, function(x) {
  x_shifted <- x - delta
  x_shifted[x_shifted <= 0] <- 1e-6
  p * pweibull(x_shifted, shape=beta1, scale=alpha1) +
    (1 - p) * pweibull(x_shifted, shape=beta2, scale=alpha2)
})
print(ks_stat)

# Chi-square and RMSE
obs_counts <- hist(wind_speed, breaks=10, plot=FALSE)$counts
bin_edges <- hist(wind_speed, breaks=10, plot=FALSE)$breaks
expected_probs <- diff(sapply(bin_edges, function(b) {
  b_shifted <- b - delta
  b_shifted[b_shifted <= 0] <- 1e-6
  p * pweibull(b_shifted, shape=beta1, scale=alpha1) +
    (1 - p) * pweibull(b_shifted, shape=beta2, scale=alpha2)
}))
expected_counts <- expected_probs * length(wind_speed)
chisq_stat <- sum((obs_counts - expected_counts)^2 / expected_counts)
p_value_chisq <- pchisq(chisq_stat, df=length(obs_counts) - 6, lower.tail=FALSE)
cat("Chi-square =", chisq_stat, "p =", p_value_chisq, "\n")

# RMSE
fitted_vals <- bimodal_weibull_shifted(wind_speed, alpha1, beta1, alpha2, beta2, p, delta)
rmse <- sqrt(mean((fitted_vals - density(wind_speed)$y[1:length(wind_speed)])^2))
cat("RMSE =", rmse, "\n")

# Monte Carlo Simulation
set.seed(123)
sim_data <- c(
  rweibull(round(p * 10000), shape=beta1, scale=alpha1),
  rweibull(round((1 - p) * 10000), shape=beta2, scale=alpha2)
) + delta

# Validation Plot
plot(density(wind_speed), col="blue", lwd=2, main="Validation - Modified Weibull_Lag2",
     xlab="Wind Speed (m/s)")
lines(density(sim_data), col="darkred", lwd=2, lty=2)
legend("topright", legend=c("Observed", "Simulated"), col=c("blue", "darkred"), lwd=2, lty=1:2)