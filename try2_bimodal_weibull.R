# Load necessary libraries
library(MASS)
library(survival)
library(fitdistrplus)
library(VGAM)
library(stats4)
library(splines)
library(ggplot2)
library(stats4)

# Read the data
data <- read.csv("abk.vtc.2dp.csv")
wind_speed <- na.omit(data$Windspeed)

# Plot histogram
hist(wind_speed, probability=TRUE, breaks=20, col="skyblue",
     main="Ogbomoso", xlab="WindSpeed (m/s)")

# Define Bimodal Weibull PDF
bimodal_weibull <- function(x, alpha1, beta1, alpha2, beta2, p) {
  p * dweibull(x, shape=beta1, scale=alpha1) +
    (1 - p) * dweibull(x, shape=beta2, scale=alpha2)
}

# Negative Log-Likelihood Function
nll <- function(alpha1, beta1, alpha2, beta2, p) {
  if (any(c(alpha1, beta1, alpha2, beta2, p) <= 0) || p >= 1) return(1e6)
  -sum(log(bimodal_weibull(wind_speed, alpha1, beta1, alpha2, beta2, p)))
}

# Initial parameter guesses
start_vals <- list(alpha1=2, beta1=2, alpha2=6, beta2=2, p=0.5)

# Maximum Likelihood Estimation
fit <- mle(nll, start=start_vals, method="L-BFGS-B",
           lower=c(0.1, 0.1, 0.1, 0.1, 0.01), upper=c(20, 10, 20, 10, 0.99))

# Extract fitted parameters
params <- coef(fit)
alpha1 <- params["alpha1"]
beta1 <- params["beta1"]
alpha2 <- params["alpha2"]
beta2 <- params["beta2"]
p <- params["p"]

# Overlay fitted PDF
x_seq <- seq(min(wind_speed), max(wind_speed), length.out=100)
y_fit <- bimodal_weibull(x_seq, alpha1, beta1, alpha2, beta2, p)
lines(x_seq, y_fit, col="red", lwd=2)

# Goodness-of-fit
# Kolmogorov-Smirnov Test
ks_stat <- ks.test(wind_speed, function(x) {
  p * pweibull(x, shape=beta1, scale=alpha1) +
    (1 - p) * pweibull(x, shape=beta2, scale=alpha2)
})
print(ks_stat)

# Chi-Square Test (using histogram bins)
obs_counts <- hist(wind_speed, breaks=10, plot=FALSE)$counts
bin_edges <- hist(wind_speed, breaks=10, plot=FALSE)$breaks
expected_probs <- diff(sapply(bin_edges, function(b) {
  p * pweibull(b, shape=beta1, scale=alpha1) +
    (1 - p) * pweibull(b, shape=beta2, scale=alpha2)
}))
expected_counts <- expected_probs * length(wind_speed)
chisq_stat <- sum((obs_counts - expected_counts)^2 / expected_counts)
df_chisq <- length(obs_counts) - 5  # parameters estimated
p_value_chisq <- pchisq(chisq_stat, df=df_chisq, lower.tail=FALSE)
cat("Chi-square statistic =", chisq_stat, "p-value =", p_value_chisq, "\n")

# RMSE
fitted_values <- bimodal_weibull(wind_speed, alpha1, beta1, alpha2, beta2, p)
rmse <- sqrt(mean((fitted_values - density(wind_speed)$y[1:length(wind_speed)])^2))
cat("RMSE =", rmse, "\n")

# Monte Carlo Simulation for Validation
set.seed(123)
sim_data <- c(
  rweibull(round(p * 10000), shape=beta1, scale=alpha1),
  rweibull(round((1 - p) * 10000), shape=beta2, scale=alpha2)
)

# Compare density plots
plot(density(wind_speed), col="blue", lwd=2, main="Validation-Ogbomoso",
     xlab="Wind Speed (m/s)")
lines(density(sim_data), col="darkgreen", lwd=2, lty=2)
legend("topright", legend=c("observed", "simulated"), col=c("blue", "darkgreen"), lwd=2, lty=1:2)