data <- read.table("gauge.txt", header = TRUE)

# Question 1
linear_model <- lm(gain ~ density, data = data)
summary(linear_model)
par(mfrow = c(1, 2))
plot(data$density, data$gain, xlab = "Density", ylab = "Gain", main = "Linear Fit of Gain vs. Density")
abline(linear_model, col = "red", lwd = 2)
plot(linear_model$fitted.values, resid(linear_model), main = "Residuals vs Fitted Values", xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)
par(mfrow = c(1, 1))

# Question 2
exp_model <- nls(gain ~ A * exp(beta * density), data = data, start = list(A = max(data$gain), beta = -5))
summary(exp_model)
A <- coef(exp_model)["A"]
beta <- coef(exp_model)["beta"]
par(mfrow = c(1, 2))
plot(data$density, data$gain, main = "Exponential Fit of Gain vs. Density", xlab = "Density", ylab = "Gain")
curve(A * exp(beta * x), from = min(data$density), to = max(data$density), add = TRUE, col = "red", lwd = 2)
residuals_exp <- resid(exp_model)
plot(fitted(exp_model), residuals_exp, main = "Residuals vs Fitted Values (Exponential Model)", xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)
par(mfrow = c(1, 1))

# Question 3
set.seed(123)
data$noisy_density <- data$density + rnorm(n = nrow(data), mean = 0, sd = 0.01)
noisy_exp_model <- nls(gain ~ A * exp(beta * noisy_density), data = data, start = list(A = 430, beta = -5))
summary(noisy_exp_model)
par(mfrow = c(1, 2))
plot(fitted(noisy_exp_model), resid(noisy_exp_model), main = "Residuals (Noisy Densities)", xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)
plot(fitted(exp_model), resid(exp_model), main = "Residuals (Original Densities)", xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)
par(mfrow = c(1, 2))
plot(data$density, data$gain, main = "Original Fitted Model", xlab = "Density", ylab = "Gain", col = "blue")
curve(coef(exp_model)["A"] * exp(coef(exp_model)["beta"] * x), add = TRUE, col = "red", lwd = 2, lty = 1)
legend("topright", legend = c("Original Fit"), col = "red", lty = 1, lwd = 2)
plot(data$density, data$gain, main = "Noisy Fitted Model", xlab = "Density", ylab = "Gain", col = "blue")
curve(coef(noisy_exp_model)["A"] * exp(coef(noisy_exp_model)["beta"] * x), add = TRUE, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Noisy Fit"), col = "red", lty = 2, lwd = 2)

# Question 4
predict_gain <- function(density, model, conf_level = 0.95) {
    A <- coef(model)["A"]
    beta <- coef(model)["beta"]
    predicted_gain <- A * exp(beta * density)
    param_var <- vcov(model)
    var_A <- param_var["A", "A"]
    var_beta <- param_var["beta", "beta"]
    cov_A_beta <- param_var["A", "beta"]
    se <- sqrt((exp(beta * density))^2 * var_A + (A * density * exp(beta * density))^2 * var_beta + 2 * (exp(beta * density)) * (A * density * exp(beta * density)) * cov_A_beta)
    z <- qnorm(1 - (1 - conf_level) / 2)
    lower <- predicted_gain - z * se
    upper <- predicted_gain + z * se
    return(list(predicted_gain = predicted_gain, lower = lower, upper = upper))
}
densities <- c(0.508, 0.001)
predictions <- lapply(densities, function(d) predict_gain(density = d, model = exp_model))
for (i in seq_along(densities)) {
    cat(sprintf("Density: %.3f\n", densities[i]))
    cat(sprintf("Predicted Gain: %.2f\n", predictions[[i]]$predicted_gain))
    cat(sprintf("95%% Prediction Interval: [%.2f, %.2f]\n\n", predictions[[i]]$lower, predictions[[i]]$upper))
}
plot(data$density, data$gain, main = "Forward Prediction with Confidence Intervals", xlab = "Density", ylab = "Gain", pch = 19, col = "blue")
curve(coef(exp_model)["A"] * exp(coef(exp_model)["beta"] * x), from = min(data$density), to = max(data$density), add = TRUE, col = "red", lwd = 2)
for (i in seq_along(densities)) {
    points(densities[i], predictions[[i]]$predicted_gain, col = "pink", pch = 19)
    arrows(densities[i], predictions[[i]]$lower, densities[i], predictions[[i]]$upper, angle = 90, code = 3, length = 0.05, col = "pink")
}
legend("topright", legend = c("Fitted Model", "Prediction"), col = c("red", "pink"), lty = 1, pch = c(NA, 19), bty = "n")

# Question 5
predict_density <- function(gain, model, conf_level = 0.95) {
    A <- coef(model)["A"]
    beta <- coef(model)["beta"]
    predicted_density <- log(gain / A) / beta
    param_var <- vcov(model)
    var_A <- param_var["A", "A"]
    var_beta <- param_var["beta", "beta"]
    cov_A_beta <- param_var["A", "beta"]
    se <- sqrt((1 / (gain * beta))^2 * var_A + ((log(gain / A) / beta^2))^2 * var_beta + 2 * (1 / (gain * beta)) * (log(gain / A) / beta^2) * cov_A_beta)
    z <- qnorm(1 - (1 - conf_level) / 2)
    lower <- predicted_density - z * se
    upper <- predicted_density + z * se
    return(list(predicted_density = predicted_density, lower = lower, upper = upper))
}
gains <- c(38.6, 426.7)
reverse_predictions <- lapply(gains, function(g) predict_density(gain = g, model = exp_model))
for (i in seq_along(gains)) {
    cat(sprintf("Gain: %.1f\n", gains[i]))
    cat(sprintf("Predicted Density: %.3f\n", reverse_predictions[[i]]$predicted_density))
    cat(sprintf("95%% Prediction Interval: [%.3f, %.3f]\n\n", reverse_predictions[[i]]$lower, reverse_predictions[[i]]$upper))
}

# Question 6
cross_validate <- function(omit_density, data, conf_level = 0.95) {
    cv_data <- subset(data, density != omit_density)
    cv_model <- nls(gain ~ A * exp(beta * density), data = cv_data, start = list(A = max(cv_data$gain), beta = -5))
    avg_gain <- mean(subset(data, density == omit_density)$gain)
    prediction <- predict_density(gain = avg_gain, model = cv_model, conf_level = conf_level)
    list(omitted_density = omit_density, avg_gain = avg_gain, prediction = prediction)
}
cv_results <- lapply(c(0.508, 0.001), function(d) cross_validate(d, data))
for (result in cv_results) {
    cat(sprintf("Omitted Density: %.3f\n", result$omitted_density))
    cat(sprintf("Average Gain: %.2f\n", result$avg_gain))
    cat(sprintf("Predicted Density: %.3f\n", result$prediction$predicted_density))
    cat(sprintf("95%% Prediction Interval: [%.3f, %.3f]\n\n", result$prediction$lower, result$prediction$upper))
}

# Advanced Analysis
density_q1 <- quantile(data$density, 0.25)
density_q3 <- quantile(data$density, 0.75)
intermediate_mask <- data$density > density_q1 & data$density < density_q3
intermediate_data <- data[intermediate_mask, ]
if (nrow(intermediate_data) == 0) {
    stop("test")
}
intermediate_data$predicted_gain <- fitted(exp_model)[intermediate_mask]
intermediate_data$bias <- intermediate_data$predicted_gain - intermediate_data$gain
mean_bias <- mean(intermediate_data$bias, na.rm = TRUE)
bias_summary <- summary(intermediate_data$bias)
list(mean_bias = mean_bias, bias_summary = bias_summary)
