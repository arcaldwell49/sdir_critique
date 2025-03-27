# Simulation of a parallel group RCT with randomized assignment
# Testing boundary conditions where SDs of change are equal despite correlation < 1

# Set seed for reproducibility
set.seed(12345)

# Function to generate data with proper randomization to treatment or control
# Key parameters:
#   cor_potential_outcomes: correlation between Y(0) and Y(1)
#   control_error_sd: SD of within-subject variation in control group
#   treatment_error_sd: SD of within-subject variation in treatment group
simulate_parallel_rct <- function(n_total = 1000, 
                                  mean_effect = 5, 
                                  cor_potential_outcomes = 0.7,
                                  sd_potential_outcomes = 8,
                                  control_error_sd = 4,
                                  treatment_error_sd = NULL) {
  
  # If treatment_error_sd is not specified, use the same as control
  if (is.null(treatment_error_sd)) {
    treatment_error_sd <- control_error_sd
  }
  
  # Generate correlated potential outcomes for ALL subjects
  z1 <- rnorm(n_total, mean = 0, sd = 1)
  z2 <- rnorm(n_total, mean = 0, sd = 1)
  
  # Create correlation between potential outcomes
  y0 <- z1 * sd_potential_outcomes  # Potential outcome under control
  y1 <- cor_potential_outcomes * y0 + 
    sqrt(1 - cor_potential_outcomes^2) * z2 * sd_potential_outcomes + 
    mean_effect  # Potential outcome under treatment
  
  # Calculate the theoretical SD of individual treatment effects
  theoretical_sd_effects <- sqrt(2 * (1 - cor_potential_outcomes)) * sd_potential_outcomes
  
  # Calculate individual treatment effects for all subjects
  individual_effects <- y1 - y0
  
  # Randomize subjects to treatment or control
  # 1 = treatment, 0 = control
  treatment_assignment <- sample(c(0, 1), n_total, replace = TRUE)
  
  # Generate baseline measurements for all subjects
  # Based on potential outcome y0 (everyone starts at control state)
  baseline_error <- rnorm(n_total, 0, ifelse(treatment_assignment == 1, 
                                             treatment_error_sd, 
                                             control_error_sd))
  baseline <- y0 + baseline_error
  
  # Generate follow-up measurements based on assignment
  # Control subjects: measured y0 at follow-up
  # Treatment subjects: measured y1 at follow-up
  followup_error <- rnorm(n_total, 0, ifelse(treatment_assignment == 1, 
                                             treatment_error_sd, 
                                             control_error_sd))
  
  # Initialize follow-up vector
  followup <- numeric(n_total)
  
  # Assign follow-up values based on treatment assignment
  for (i in 1:n_total) {
    if (treatment_assignment[i] == 1) {
      # Subject received treatment
      followup[i] <- y1[i] + followup_error[i]
    } else {
      # Subject received control
      followup[i] <- y0[i] + followup_error[i]
    }
  }
  
  # Calculate change for each subject
  change <- followup - baseline
  
  # Extract data for each group
  control_indices <- which(treatment_assignment == 0)
  treatment_indices <- which(treatment_assignment == 1)
  
  # Calculate changes for each group
  change_control <- change[control_indices]
  change_treatment <- change[treatment_indices]
  
  # Calculate SDs of changes
  sd_change_control <- sd(change_control)
  sd_change_treatment <- sd(change_treatment)
  
  # Calculate SDir
  sdir <- sign(sd_change_treatment^2 - sd_change_control^2) * sqrt(abs(sd_change_treatment^2 - sd_change_control^2))
  
  # Return results
  list(
    cor_potential_outcomes = cor_potential_outcomes,
    empirical_cor = cor(y0, y1),
    treatment_error_sd = treatment_error_sd,
    control_error_sd = control_error_sd,
    ratio_error_sd = treatment_error_sd / control_error_sd,
    n_control = length(control_indices),
    n_treatment = length(treatment_indices),
    true_effects = individual_effects,
    sd_true_effects = sd(individual_effects),
    theoretical_sd_effects = theoretical_sd_effects,
    sdir = sdir,
    sd_control = sd_change_control,
    sd_treatment = sd_change_treatment,
    sd_ratio = sd_change_treatment / sd_change_control,
    change_control = change_control,
    change_treatment = change_treatment,
    treatment_assignment = treatment_assignment,
    baseline = baseline,
    followup = followup
  )
}

# Function to find the treatment error SD that produces equal SDs of change
find_equal_sd_factor <- function(cor_potential_outcomes = 0.7,
                                 sd_potential_outcomes = 8,
                                 control_error_sd = 4,
                                 n_total = 5000,
                                 tolerance = 0.01) {
  
  # Start with a treatment error SD that's likely too small
  lower_bound <- 0.1 * control_error_sd
  upper_bound <- control_error_sd
  
  max_iterations <- 20
  current_iteration <- 0
  
  while ((upper_bound - lower_bound) > tolerance * control_error_sd && 
         current_iteration < max_iterations) {
    
    current_iteration <- current_iteration + 1
    mid_point <- (lower_bound + upper_bound) / 2
    
    # Try this error SD
    sim <- simulate_parallel_rct(
      n_total = n_total,
      cor_potential_outcomes = cor_potential_outcomes,
      sd_potential_outcomes = sd_potential_outcomes,
      control_error_sd = control_error_sd,
      treatment_error_sd = mid_point
    )
    
    sd_ratio <- sim$sd_treatment / sim$sd_control
    
    if (sd_ratio > 1.0) {
      # Treatment SD still too large
      upper_bound <- mid_point
    } else {
      # Treatment SD too small or equal
      lower_bound <- mid_point
    }
    
    cat("Iteration:", current_iteration, 
        "Treatment Error SD:", mid_point, 
        "SD Ratio:", sd_ratio, 
        "Control SD:", sim$sd_control,
        "Treatment SD:", sim$sd_treatment, "\n")
    
    # If we're close enough to equal SDs, break
    if (abs(sd_ratio - 1.0) < tolerance) {
      break
    }
  }
  
  return(mid_point)
}

# Let's test various correlations
correlations <- c(0.9, 0.8, 0.7, 0.6, 0.5)
results <- data.frame(
  correlation = numeric(),
  treatment_error_sd = numeric(),
  control_error_sd = numeric(),
  error_ratio = numeric(),
  sd_control = numeric(),
  sd_treatment = numeric(),
  sd_ratio = numeric(),
  sdir = numeric(),
  theoretical_sd_effects = numeric(),
  n_simulations = numeric()
)

# Find the treatment error SD for each correlation
# Use multiple simulations for more stable results
n_simulations <- 5

for (cor_val in correlations) {
  cat("\n===== Finding measurement error SD for correlation =", cor_val, "=====\n")
  
  # First find the approximate treatment error SD that gives equal SDs
  optimal_treatment_error_sd <- find_equal_sd_factor(
    cor_potential_outcomes = cor_val,
    n_total = 5000  # Larger sample for more stable estimation
  )
  
  # Run multiple simulations with this error SD
  sd_controls <- numeric(n_simulations)
  sd_treatments <- numeric(n_simulations)
  sdirs <- numeric(n_simulations)
  
  for (i in 1:n_simulations) {
    sim <- simulate_parallel_rct(
      n_total = 2000,
      cor_potential_outcomes = cor_val,
      control_error_sd = 4,
      treatment_error_sd = optimal_treatment_error_sd
    )
    
    sd_controls[i] <- sim$sd_control
    sd_treatments[i] <- sim$sd_treatment
    sdirs[i] <- sim$sdir
  }
  
  # Add to results - use means across simulations
  results <- rbind(results, data.frame(
    correlation = cor_val,
    treatment_error_sd = optimal_treatment_error_sd,
    control_error_sd = 4,
    error_ratio = optimal_treatment_error_sd / 4,
    sd_control = mean(sd_controls),
    sd_treatment = mean(sd_treatments),
    sd_ratio = mean(sd_treatments) / mean(sd_controls),
    sdir = mean(sdirs),
    theoretical_sd_effects = sqrt(2 * (1 - cor_val)) * 8,
    n_simulations = n_simulations
  ))
}

# Print the results
print(results)

# Let's analyze one specific case in detail
detailed_cor <- 0.7
detailed_treatment_error_sd <- results$treatment_error_sd[results$correlation == detailed_cor]

# Run a detailed simulation
detailed_sim <- simulate_parallel_rct(
  n_total = 10000,  # Larger sample for more stable results
  cor_potential_outcomes = detailed_cor,
  control_error_sd = 4,
  treatment_error_sd = detailed_treatment_error_sd
)

# Print detailed results
cat("\n===== Detailed Simulation with Correlation =", detailed_cor, "=====\n")
cat("Control Error SD:", detailed_sim$control_error_sd, "\n")
cat("Treatment Error SD:", detailed_sim$treatment_error_sd, "\n")
cat("Ratio of Error SDs (treatment/control):", detailed_sim$ratio_error_sd, "\n")
cat("Number of control subjects:", detailed_sim$n_control, "\n")
cat("Number of treatment subjects:", detailed_sim$n_treatment, "\n")
cat("SD of change in control group:", detailed_sim$sd_control, "\n")
cat("SD of change in treatment group:", detailed_sim$sd_treatment, "\n")
cat("Ratio of SDs (treatment/control):", detailed_sim$sd_ratio, "\n")
cat("Estimated SDir:", detailed_sim$sdir, "\n")
cat("True SD of individual effects:", detailed_sim$sd_true_effects, "\n")
cat("Theoretical SD of individual effects:", detailed_sim$theoretical_sd_effects, "\n")

# Demonstrate the problem with SDir in this scenario
cat("\nProblem with SDir: When SDs of change are equal, SDir ≈ 0\n")
cat("But true individual response heterogeneity exists (SD =", 
    detailed_sim$sd_true_effects, ")\n")

# Also create a scenario where control has larger SD than treatment
reverse_error_sd <- detailed_treatment_error_sd * 0.8  # Make treatment error even smaller

reverse_sim <- simulate_parallel_rct(
  n_total = 10000,
  cor_potential_outcomes = detailed_cor,
  control_error_sd = 4,
  treatment_error_sd = reverse_error_sd
)

cat("\n===== Simulation with Larger Control SD =====\n")
cat("Control Error SD:", reverse_sim$control_error_sd, "\n")
cat("Treatment Error SD:", reverse_sim$treatment_error_sd, "\n")
cat("Ratio of Error SDs (treatment/control):", reverse_sim$ratio_error_sd, "\n")
cat("SD of change in control group:", reverse_sim$sd_control, "\n")
cat("SD of change in treatment group:", reverse_sim$sd_treatment, "\n")
cat("Ratio of SDs (treatment/control):", reverse_sim$sd_ratio, "\n")
cat("Estimated SDir:", reverse_sim$sdir, "\n") 
cat("True SD of individual effects:", reverse_sim$sd_true_effects, "\n")
cat("Theoretical SD of individual effects:", reverse_sim$theoretical_sd_effects, "\n")
cat("Problem: SDir would be reported as 0 (or undefined for negative values)\n")
cat("despite true individual response heterogeneity (SD =", 
    reverse_sim$sd_true_effects, ")\n")

# Create histograms to visualize the distribution of changes
par(mfrow=c(2, 1))
hist(detailed_sim$change_control, 
     main="Distribution of Changes in Control Group",
     xlab="Change", 
     col="lightblue",
     xlim=c(-20, 20),
     breaks=30)
abline(v=mean(detailed_sim$change_control), col="red", lwd=2)
sd_text <- paste("SD =", round(detailed_sim$sd_control, 2))
text(10, 200, sd_text, pos=4)

hist(detailed_sim$change_treatment, 
     main="Distribution of Changes in Treatment Group",
     xlab="Change", 
     col="lightgreen",
     xlim=c(-20, 20),
     breaks=30)
abline(v=mean(detailed_sim$change_treatment), col="red", lwd=2)
sd_text <- paste("SD =", round(detailed_sim$sd_treatment, 2))
text(10, 200, sd_text, pos=4)

# Reset the plotting layout
par(mfrow=c(1, 1))

# Plot baseline vs. change for both groups
plot(detailed_sim$baseline[detailed_sim$treatment_assignment == 0], 
     detailed_sim$followup[detailed_sim$treatment_assignment == 0] - 
       detailed_sim$baseline[detailed_sim$treatment_assignment == 0],
     main="Baseline vs. Change by Group",
     xlab="Baseline",
     ylab="Change",
     pch=19,
     col=rgb(0, 0, 1, 0.3),  # Blue for control
     xlim=range(detailed_sim$baseline),
     ylim=range(c(detailed_sim$change_control, detailed_sim$change_treatment)))

points(detailed_sim$baseline[detailed_sim$treatment_assignment == 1], 
       detailed_sim$followup[detailed_sim$treatment_assignment == 1] - 
         detailed_sim$baseline[detailed_sim$treatment_assignment == 1],
       pch=19,
       col=rgb(1, 0, 0, 0.3))  # Red for treatment

abline(h=mean(detailed_sim$change_control), col="blue", lwd=2)
abline(h=mean(detailed_sim$change_treatment), col="red", lwd=2)
legend("topright", 
       legend=c("Control", "Treatment"), 
       col=c("blue", "red"), 
       pch=19)

# Calculate and display the variance components
var_control_change <- var(detailed_sim$change_control)
var_treatment_change <- var(detailed_sim$change_treatment)
var_individual_effects <- var(detailed_sim$true_effects)

cat("\n===== Variance Components =====\n")
cat("Variance of changes in control group:", var_control_change, "\n")
cat("Variance of changes in treatment group:", var_treatment_change, "\n")
cat("Variance of true individual effects:", var_individual_effects, "\n")

# Run a standard scenario for comparison (equal measurement error)
standard_sim <- simulate_parallel_rct(
  n_total = 10000,
  cor_potential_outcomes = detailed_cor,
  control_error_sd = 4,
  treatment_error_sd = 4  # Equal to control error
)

cat("\n===== Standard Scenario (Equal Measurement Error) =====\n")
cat("Control Error SD:", standard_sim$control_error_sd, "\n")
cat("Treatment Error SD:", standard_sim$treatment_error_sd, "\n")
cat("SD of change in control group:", standard_sim$sd_control, "\n")
cat("SD of change in treatment group:", standard_sim$sd_treatment, "\n")
cat("Ratio of SDs (treatment/control):", standard_sim$sd_ratio, "\n")
cat("Estimated SDir:", standard_sim$sdir, "\n")
cat("True SD of individual effects:", standard_sim$sd_true_effects, "\n")
cat("In this case, SDir correctly estimates the SD of true effects\n")