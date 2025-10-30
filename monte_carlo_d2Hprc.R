# Monte Carlo Simulation for δ²H Precipitation Reconstruction
# Purpose: Reconstruct δ²H precipitation (δ²Hprc) from plant wax δ²H (n-C29 alkane)
#          using vegetation-corrected fractionation (R.A index) for the Rotsee sediment
#          core, as described in Santos et al. (2025, GRL submission).
# Author: R. N. Santos
# Created: 2025-09-10
# Last Modified: 2025-09-19
# Dependencies: R (>= 4.4.1), ggplot2, readxl, writexl, dplyr, tidyr

# Load required libraries
library(ggplot2)
library(readxl)
library(writexl)
library(dplyr)
library(tidyr)

# Note: Set your working directory to the folder containing the input data
# setwd("path/to/your/directory") # Uncomment and modify as needed

# Read input data
#You may adjust the fGR variable. Here, we use R.A_fGR, which is available on the "R.A_fGR" spreadsheet of the ROT21_Inputs_MonteCarlo_simulation_d2Hprc.xlsx file. 
input_file <- "ROT21_Inputs_MonteCarlo_simulation_d2Hprc.xlsx"
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}
data <- read_excel(input_file)

# Verify required columns
required_cols <- c("Age_BP_ka", "d2H_C29", "R.A_f_GR")
if (!all(required_cols %in% colnames(data))) {
  stop("Missing required columns in input data: ", paste(required_cols, collapse = ", "))
}

# Validate f_GR values (must be between 0 and 1)
if (any(data$R.A_f_GR < 0 | data$R.A_f_GR > 1, na.rm = TRUE)) {
  stop("R.A_f_GR values must be between 0 and 1")
}

# Define fractionation values and uncertainties (See Santos et al., 025 and references therein)
epsilon_C29_WP <- -110  # Mean εapp for C29 in Woody Plants (‰)
sigma_C29_WP <- 21      # Uncertainty for C29 in Woody Plants (‰)
epsilon_C29_GR <- -165  # Mean εapp for C29 in Grasses (‰)
sigma_C29_GR <- 25      # Uncertainty for C29 in Grasses (‰)

# Monte Carlo simulation parameters
n_simulations <- 10000
set.seed(12345)  # For reproducibility

# Initialize list to store simulation results
simulations_list <- list()
skipped_rows <- 0  # Track skipped rows

# Perform Monte Carlo simulations for each sample
for (i in 1:nrow(data)) {
  # Extract sample data
  age <- data$Age_BP_ka[i]
  d2H_C29 <- data$d2H_C29[i]
  f_grasses_RA <- data$R.A_f_GR[i]
  
  # Skip if any required value is NA
  if (any(is.na(c(d2H_C29, f_grasses_RA)))) {
    skipped_rows <- skipped_rows + 1
    message("Skipping row ", i, " due to NA values")
    next
  }
  
  # Generate Monte Carlo samples for fractionation factors
  epsilon_C29_WP_samples <- rnorm(n_simulations, mean = epsilon_C29_WP, sd = sigma_C29_WP)
  epsilon_C29_GR_samples <- rnorm(n_simulations, mean = epsilon_C29_GR, sd = sigma_C29_GR)
  
  # Compute vegetation-weighted fractionation factor (ε_mix)
  epsilon_mix_RA <- (f_grasses_RA * epsilon_C29_GR_samples) + ((1 - f_grasses_RA) * epsilon_C29_WP_samples)
  
  # Compute δ²H precipitation values
  d2H_prec_RA_sim <- ((d2H_C29 + 1000) / ((epsilon_mix_RA / 1000) + 1)) - 1000
  
  # Store results
  simulations_df <- data.frame(
    Age = rep(age, n_simulations),
    vegetation_correction_RA = epsilon_mix_RA,
    d2H_prec_RA = d2H_prec_RA_sim
  )
  
  simulations_list[[i]] <- simulations_df
}

# Report skipped rows
if (skipped_rows > 0) {
  message("Total rows skipped due to NA values: ", skipped_rows)
}

# Combine simulation results
final_results <- do.call(rbind, simulations_list)

# Export raw simulation results
output_sim <- "d2H_precipitation_simulations.csv"
if (file.exists(output_sim)) {
  warning("Output file already exists: ", output_sim, ". Skipping write.")
} else {
  write.csv(final_results, output_sim, row.names = FALSE)
}

# Compute confidence intervals
results_ci <- final_results %>%
  group_by(Age) %>%
  summarize(
    Mean_veg_corr_RA = mean(vegetation_correction_RA, na.rm = TRUE),
    CI_95_low_veg_RA = quantile(vegetation_correction_RA, 0.025, na.rm = TRUE),
    CI_95_high_veg_RA = quantile(vegetation_correction_RA, 0.975, na.rm = TRUE),
    Mean_RA = mean(d2H_prec_RA, na.rm = TRUE),
    CI_95_low_RA = quantile(d2H_prec_RA, 0.025, na.rm = TRUE),
    CI_95_high_RA = quantile(d2H_prec_RA, 0.975, na.rm = TRUE),
    CI_68_low_RA = quantile(d2H_prec_RA, 0.16, na.rm = TRUE),
    CI_68_high_RA = quantile(d2H_prec_RA, 0.84, na.rm = TRUE)
  ) %>%
  ungroup()

# Add LOESS-smoothed values
loess_fit <- loess(Mean_RA ~ Age, data = results_ci, span = 0.2)
results_ci$LOESS_Mean_RA <- predict(loess_fit, newdata = data.frame(Age = results_ci$Age))

# Export summarized results with LOESS values
output_ci <- "d2H_prec_reconstructions_with_CI.xlsx"
if (file.exists(output_ci)) {
  warning("Output file already exists: ", output_ci, ". Skipping write.")
} else {
  write_xlsx(results_ci, output_ci)
}

# Plot δ²H precipitation with confidence intervals and LOESS smoothing
p <- ggplot(results_ci, aes(x = Age)) +
  geom_ribbon(aes(ymin = CI_95_low_RA, ymax = CI_95_high_RA), fill = "grey80", alpha = 1) +
  geom_ribbon(aes(ymin = CI_68_low_RA, ymax = CI_68_high_RA), fill = "grey60", alpha = 1) +
  geom_point(aes(y = Mean_RA), size = 1.5, color = "black") +
  geom_line(aes(y = LOESS_Mean_RA), color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Age (ka cal BP)",
    y = expression(delta^2*H[prc] ~ "(\u2030)"),
    title = expression("δ"^2*"H"[prc] ~ "Reconstruction (R.A-Based)"),
    caption = "Shaded: 68% CI (dark gray), 95% CI (light gray); Red dashed = LOESS (span = 0.2)"
  ) +
  theme_minimal(base_size = 14)

# Display plot at the end of the script
print(p)

# Save plots
ggsave("d2H_precip_plot_with_CI_and_LOESS.png", plot = p, width = 10, height = 6, dpi = 300)
ggsave("d2H_precip_plot_with_CI_and_LOESS.pdf", plot = p, width = 10, height = 6)
