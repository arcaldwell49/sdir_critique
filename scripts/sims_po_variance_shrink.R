# Simulation of a parallel group RCT with variance shrinkage effect
# Using the potential outcomes framework to demonstrate how treatments
# can have a "normalizing" effect that reduces variance

# Set seed for reproducibility
set.seed(12345)

# Function to generate data with a normalizing treatment effect
# Key parameters:
#   baseline_dependency: controls how much treatment effect depends on baseline (0-1)
#     - 0: constant treatment effect regardless of baseline
#     - 1: complete normalization toward a target value
simulate_normalizing_treatment <- function(n_total = 1000, 
                                           mean_effect = 5, 
                                           baseline_dependency = 0.5,
                                           target_value = 50,
                                           sd_baseline = 10,
                                           within_subject_sd = 4) {
  
  # Generate baseline values (true values before any measurement error)
  true_baseline <- rnorm(n_total, mean = 40, sd = sd_baseline)
  
  # Calculate treatment effects that depend on baseline values
  # For high baseline_dependency:
  #   - Individuals with high baseline values get larger negative effects
  #   - Individuals with low baseline values get larger positive effects
  #   - Result is "normalizing" toward the target value
  
  # Pure constant effect (when baseline_dependency = 0)
  constant_effect_component <- mean_effect
  
  # Normalization effect (when baseline_dependency > 0)
  # This moves values toward the target value
  normalization_component <- (target_value - true_baseline)
  
  # Combine the components based on the baseline_dependency parameter
  individual_effects <- (1 - baseline_dependency) * constant_effect_component + 
    baseline_dependency * normalization_component
  
  # Calculate potential outcomes
  potential_outcome_control <- true_baseline  # Y(0)
  potential_outcome_treatment <- true_baseline + individual_effects  # Y(1)
  
  # Calculate the correlation between potential outcomes
  cor_potential_outcomes <- cor(potential_outcome_control, potential_outcome_treatment)
  
  # Calculate the theoretical SD of individual treatment effects
  sd_y0 <- sd(potential_outcome_control)
  sd_y1 <- sd(potential_outcome_treatment)
  cov_y0_y1 <- cov(potential_outcome_control, potential_outcome_treatment)
  theoretical_var_D <- sd_y0^2 + sd_y1^2 - 2 * cov_y0_y1
  theoretical_sd_D <- sqrt(theoretical_var_D)
  
  # Randomize subjects to treatment or control
  # 1 = treatment, 0 = control
  treatment_assignment <- sample(c(0, 1), n_total, replace = TRUE)
  
  # Generate baseline measurements for all subjects (with measurement error)
  baseline_error <- rnorm(n_total, 0, within_subject_sd)
  baseline <- true_baseline + baseline_error
  
  # Generate follow-up measurements based on assignment
  followup_error <- rnorm(n_total, 0, within_subject_sd)
  
  # Initialize follow-up vector
  followup <- numeric(n_total)
  
  # Assign follow-up values based on treatment assignment
  for (i in 1:n_total) {
    if (treatment_assignment[i] == 1) {
      # Subject received treatment
      followup[i] <- potential_outcome_treatment[i] + followup_error[i]
    } else {
      # Subject received control
      followup[i] <- potential_outcome_control[i] + followup_error[i]
    }
  }
  
  # Calculate change for each subject
  change <- followup - baseline
  
  # Extract data for each group
  control_indices <- which(treatment_assignment == 0)
  treatment_indices <- which(treatment_assignment == 1)
  
  # Calculate SDs for each group
  sd_baseline_control <- sd(baseline[control_indices])
  sd_followup_control <- sd(followup[control_indices])
  sd_baseline_treatment <- sd(baseline[treatment_indices])
  sd_followup_treatment <- sd(followup[treatment_indices])
  
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
    baseline_dependency = baseline_dependency,
    cor_potential_outcomes = cor_potential_outcomes,
    n_control = length(control_indices),
    n_treatment = length(treatment_indices),
    true_effects = individual_effects,
    sd_true_effects = sd(individual_effects),
    theoretical_sd_effects = theoretical_sd_D,
    sdir = sdir,
    sd_baseline_control = sd_baseline_control,
    sd_followup_control = sd_followup_control,
    sd_baseline_treatment = sd_baseline_treatment,
    sd_followup_treatment = sd_followup_treatment,
    sd_change_control = sd_change_control,
    sd_change_treatment = sd_change_treatment,
    baseline_control = baseline[control_indices],
    followup_control = followup[control_indices],
    baseline_treatment = baseline[treatment_indices],
    followup_treatment = followup[treatment_indices],
    change_control = change_control,
    change_treatment = change_treatment,
    treatment_assignment = treatment_assignment,
    true_baseline = true_baseline,
    potential_outcome_control = potential_outcome_control,
    potential_outcome_treatment = potential_outcome_treatment
  )
}

# Let's test various levels of baseline dependency
dependency_levels <- seq(0, 1, by = 0.2)
n_sims <- 5  # Number of simulations per dependency level

results <- data.frame(
  baseline_dependency = numeric(),
  cor_potential_outcomes = numeric(),
  sd_y0 = numeric(),
  sd_y1 = numeric(),
  sd_true_effects = numeric(),
  theoretical_sd_effects = numeric(),
  sd_baseline_control = numeric(),
  sd_followup_control = numeric(),
  variance_ratio_control = numeric(),
  sd_baseline_treatment = numeric(),
  sd_followup_treatment = numeric(),
  variance_ratio_treatment = numeric(),
  sd_change_control = numeric(),
  sd_change_treatment = numeric(),
  sdir = numeric(),
  sim = numeric()
)

for (dep in dependency_levels) {
  for (sim in 1:n_sims) {
    # Run the simulation
    sim_result <- simulate_normalizing_treatment(
      n_total = 2000,
      mean_effect = 5,
      baseline_dependency = dep,
      target_value = 50,
      sd_baseline = 10,
      within_subject_sd = 4
    )
    
    # Add to results
    results <- rbind(results, data.frame(
      baseline_dependency = dep,
      cor_potential_outcomes = sim_result$cor_potential_outcomes,
      sd_y0 = sd(sim_result$potential_outcome_control),
      sd_y1 = sd(sim_result$potential_outcome_treatment),
      sd_true_effects = sim_result$sd_true_effects,
      theoretical_sd_effects = sim_result$theoretical_sd_effects,
      sd_baseline_control = sim_result$sd_baseline_control,
      sd_followup_control = sim_result$sd_followup_control,
      variance_ratio_control = (sim_result$sd_followup_control^2) / (sim_result$sd_baseline_control^2),
      sd_baseline_treatment = sim_result$sd_baseline_treatment,
      sd_followup_treatment = sim_result$sd_followup_treatment,
      variance_ratio_treatment = (sim_result$sd_followup_treatment^2) / (sim_result$sd_baseline_treatment^2),
      sd_change_control = sim_result$sd_change_control,
      sd_change_treatment = sim_result$sd_change_treatment,
      sdir = sim_result$sdir,
      sim = sim
    ))
  }
}

# Calculate average results for each dependency level
avg_results <- aggregate(
  . ~ baseline_dependency, 
  data = results[, !names(results) %in% "sim"], 
  FUN = mean
)

# Print the results
print(avg_results)

# Analyze one specific level of baseline dependency in detail
detailed_dep <- 0.6  # Higher dependency = stronger normalizing effect
detailed_sim <- simulate_normalizing_treatment(
  n_total = 5000,
  mean_effect = 5,
  baseline_dependency = detailed_dep,
  target_value = 50,
  sd_baseline = 10,
  within_subject_sd = 4
)

# Print detailed results
cat("\n===== Detailed Simulation with Baseline Dependency =", detailed_dep, "=====\n")
cat("Correlation between potential outcomes:", detailed_sim$cor_potential_outcomes, "\n")
cat("SD of potential outcome under control:", sd(detailed_sim$potential_outcome_control), "\n")
cat("SD of potential outcome under treatment:", sd(detailed_sim$potential_outcome_treatment), "\n")
cat("Variance ratio of potential outcomes (treatment/control):", 
    (sd(detailed_sim$potential_outcome_treatment)^2) / (sd(detailed_sim$potential_outcome_control)^2), "\n")
cat("\nSD of baseline in control group:", detailed_sim$sd_baseline_control, "\n")
cat("SD of follow-up in control group:", detailed_sim$sd_followup_control, "\n")
cat("Variance ratio in control group (follow-up/baseline):", 
    (detailed_sim$sd_followup_control^2) / (detailed_sim$sd_baseline_control^2), "\n")
cat("\nSD of baseline in treatment group:", detailed_sim$sd_baseline_treatment, "\n")
cat("SD of follow-up in treatment group:", detailed_sim$sd_followup_treatment, "\n")
cat("Variance ratio in treatment group (follow-up/baseline):", 
    (detailed_sim$sd_followup_treatment^2) / (detailed_sim$sd_baseline_treatment^2), "\n")
cat("\nSD of changes in control group:", detailed_sim$sd_change_control, "\n")
cat("SD of changes in treatment group:", detailed_sim$sd_change_treatment, "\n")
cat("Ratio of SDs of changes (treatment/control):", 
    detailed_sim$sd_change_treatment / detailed_sim$sd_change_control, "\n")
cat("\nSD of true individual effects:", detailed_sim$sd_true_effects, "\n")
cat("Theoretical SD of individual effects:", detailed_sim$theoretical_sd_effects, "\n")
cat("Estimated SDir:", detailed_sim$sdir, "\n")

cat("\nKey finding: Despite true heterogeneity in treatment effects,")
cat("\nthe variance shrinkage effect can make the SD of changes in")
cat("\nthe treatment group smaller than in the control group.\n")

# Visualize the results

# 1. Plot the relationship between baseline values and individual effects
par(mfrow=c(2, 2))

plot(detailed_sim$true_baseline, detailed_sim$true_effects,
     main="Treatment Effect vs. Baseline",
     xlab="True Baseline Value",
     ylab="Individual Treatment Effect",
     pch=19,
     col=rgb(0, 0, 1, 0.3))  # Blue with transparency
abline(h=mean(detailed_sim$true_effects), col="red", lwd=2)
abline(lm(detailed_sim$true_effects ~ detailed_sim$true_baseline), col="green", lwd=2)

# 2. Plot the distribution of baseline and follow-up values in each group
hist(detailed_sim$baseline_control,
     main="Distribution of Baseline Values\nin Control Group",
     xlab="Baseline Value",
     col=rgb(0, 0, 1, 0.5),
     xlim=range(c(detailed_sim$baseline_control, detailed_sim$baseline_treatment,
                  detailed_sim$followup_control, detailed_sim$followup_treatment)),
     breaks=30)
abline(v=mean(detailed_sim$baseline_control), col="red", lwd=2)

hist(detailed_sim$followup_control,
     main="Distribution of Follow-up Values\nin Control Group",
     xlab="Follow-up Value",
     col=rgb(0, 0, 1, 0.5),
     xlim=range(c(detailed_sim$baseline_control, detailed_sim$baseline_treatment,
                  detailed_sim$followup_control, detailed_sim$followup_treatment)),
     breaks=30)
abline(v=mean(detailed_sim$followup_control), col="red", lwd=2)

hist(detailed_sim$baseline_treatment,
     main="Distribution of Baseline Values\nin Treatment Group",
     xlab="Baseline Value",
     col=rgb(1, 0, 0, 0.5),
     xlim=range(c(detailed_sim$baseline_control, detailed_sim$baseline_treatment,
                  detailed_sim$followup_control, detailed_sim$followup_treatment)),
     breaks=30)
abline(v=mean(detailed_sim$baseline_treatment), col="red", lwd=2)

hist(detailed_sim$followup_treatment,
     main="Distribution of Follow-up Values\nin Treatment Group",
     xlab="Follow-up Value",
     col=rgb(1, 0, 0, 0.5),
     xlim=range(c(detailed_sim$baseline_control, detailed_sim$baseline_treatment,
                  detailed_sim$followup_control, detailed_sim$followup_treatment)),
     breaks=30)
abline(v=mean(detailed_sim$followup_treatment), col="red", lwd=2)

# Reset the plotting layout
par(mfrow=c(1, 1))

# 3. Plot baseline vs. follow-up for both groups
plot(detailed_sim$baseline_control, detailed_sim$followup_control,
     main="Baseline vs. Follow-up by Group",
     xlab="Baseline Value",
     ylab="Follow-up Value",
     pch=19,
     col=rgb(0, 0, 1, 0.3),  # Blue for control
     xlim=range(c(detailed_sim$baseline_control, detailed_sim$baseline_treatment)),
     ylim=range(c(detailed_sim$followup_control, detailed_sim$followup_treatment)))

points(detailed_sim$baseline_treatment, detailed_sim$followup_treatment,
       pch=19,
       col=rgb(1, 0, 0, 0.3))  # Red for treatment

# Add identity line (y = x)
abline(0, 1, col="gray", lwd=2, lty=2)

# Add regression lines
control_lm <- lm(detailed_sim$followup_control ~ detailed_sim$baseline_control)
treatment_lm <- lm(detailed_sim$followup_treatment ~ detailed_sim$baseline_treatment)

abline(control_lm, col="blue", lwd=2)
abline(treatment_lm, col="red", lwd=2)

# Add horizontal lines at the means
abline(h=mean(detailed_sim$followup_control), col="blue", lwd=1, lty=3)
abline(h=mean(detailed_sim$followup_treatment), col="red", lwd=1, lty=3)

# Add target value line
abline(h=50, col="green", lwd=2, lty=4)

legend("topleft", 
       legend=c("Control", "Treatment", "Identity Line", "Target Value"), 
       col=c("blue", "red", "gray", "green"), 
       pch=c(19, 19, NA, NA),
       lty=c(1, 1, 2, 4),
       lwd=c(2, 2, 2, 2))

# 4. Plot the relationship between baseline dependency and variance ratio
dependency_plot_data <- data.frame(
  dependency = dependency_levels,
  variance_ratio_control = numeric(length(dependency_levels)),
  variance_ratio_treatment = numeric(length(dependency_levels))
)

for (i in 1:length(dependency_levels)) {
  dep <- dependency_levels[i]
  dep_data <- subset(results, baseline_dependency == dep)
  dependency_plot_data$variance_ratio_control[i] <- mean(dep_data$variance_ratio_control)
  dependency_plot_data$variance_ratio_treatment[i] <- mean(dep_data$variance_ratio_treatment)
}

plot(dependency_plot_data$dependency, dependency_plot_data$variance_ratio_treatment,
     type="b", col="red", lwd=2, pch=19,
     main="Effect of Baseline Dependency on Variance Ratio",
     xlab="Baseline Dependency",
     ylab="Variance Ratio (Follow-up / Baseline)",
     ylim=range(c(dependency_plot_data$variance_ratio_control, 
                  dependency_plot_data$variance_ratio_treatment)))

lines(dependency_plot_data$dependency, dependency_plot_data$variance_ratio_control,
      type="b", col="blue", lwd=2, pch=19)

abline(h=1, lty=2, col="gray")

legend("topright", 
       legend=c("Treatment Group", "Control Group", "No Change in Variance"), 
       col=c("red", "blue", "gray"), 
       lty=c(1, 1, 2),
       lwd=c(2, 2, 1),
       pch=c(19, 19, NA))

# 5. Plot the relationship between baseline dependency and SDir vs true effect SD
sdir_comparison <- data.frame(
  dependency = dependency_levels,
  sdir = numeric(length(dependency_levels)),
  true_sd = numeric(length(dependency_levels))
)

for (i in 1:length(dependency_levels)) {
  dep <- dependency_levels[i]
  dep_data <- subset(results, baseline_dependency == dep)
  sdir_comparison$sdir[i] <- mean(dep_data$sdir)
  sdir_comparison$true_sd[i] <- mean(dep_data$sd_true_effects)
}

plot(sdir_comparison$dependency, sdir_comparison$true_sd,
     type="b", col="green", lwd=2, pch=19,
     main="SDir vs True SD of Effects",
     xlab="Baseline Dependency",
     ylab="Standard Deviation",
     ylim=range(c(0, sdir_comparison$true_sd)))

lines(sdir_comparison$dependency, sdir_comparison$sdir,
      type="b", col="purple", lwd=2, pch=19)

legend("topright", 
       legend=c("True SD of Effects", "Estimated SDir"), 
       col=c("green", "purple"), 
       lty=1,
       lwd=2,
       pch=19)