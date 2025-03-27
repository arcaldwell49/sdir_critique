# Simulation to test how well SDir reflects Gadbury's heterogeneity of treatment effect (HTE)
# This simulation creates data with known individual treatment response heterogeneity
# and compares how well the SDir method captures the true heterogeneity

# Set seed for reproducibility
set.seed(997741)

# Function to generate data with known treatment response heterogeneity
simulate_data <- function(n_subjects = 200, mean_effect = 5, 
                          true_sd_response = 2, 
                          within_subject_sd = 3) {
  
  # Generate true baseline values
  true_baseline <- rnorm(n_subjects, mean = 35, sd = 8)
  
  # Generate individual treatment effects
  if(true_sd_response == 0) {
    # No heterogeneity - everyone gets the same effect
    individual_effects <- rep(mean_effect, n_subjects)
  } else {
    # With heterogeneity - effects vary according to specified SD
    individual_effects <- rnorm(n_subjects, mean = mean_effect, sd = true_sd_response)
  }
  
  # Calculate true follow-up values
  true_followup_control <- true_baseline  # No effect in control group
  true_followup_treatment <- true_baseline + individual_effects
  
  # Add random within-subject variation to create observed values
  # Each measurement has its own random error
  obs_baseline_control <- true_baseline + rnorm(n_subjects, 0, within_subject_sd)
  obs_followup_control <- true_followup_control + rnorm(n_subjects, 0, within_subject_sd)
  
  obs_baseline_treatment <- true_baseline + rnorm(n_subjects, 0, within_subject_sd)
  obs_followup_treatment <- true_followup_treatment + rnorm(n_subjects, 0, within_subject_sd)
  
  # Calculate observed changes
  change_control <- obs_followup_control - obs_baseline_control
  change_treatment <- obs_followup_treatment - obs_baseline_treatment
  
  # Calculate SDir (as per Atkinson & Batterham)
  sd_change_control <- sd(change_control)
  sd_change_treatment <- sd(change_treatment)
  sdir <- sign(sd_change_treatment^2 - sd_change_control^2) * sqrt(abs(sd_change_treatment^2 - sd_change_control^2))
  #sqrt(max(0, sd_change_treatment^2 - sd_change_control^2))
  # Return the results
  list(
    true_sd = true_sd_response,
    sdir = sdir,
    sd_control = sd_change_control,
    sd_treatment = sd_change_treatment,
    sd_true_effects = sd(individual_effects),
    change_control = change_control,
    change_treatment = change_treatment,
    true_effects = individual_effects
  )
}

# Run simulations with varying levels of true heterogeneity
sd_values <- c(0, 1, 2, 3, 4, 5)
n_reps <- 100
results <- data.frame(
  true_sd = numeric(),
  sdir = numeric(),
  sd_control = numeric(),
  sd_treatment = numeric(),
  rep = numeric()
)

# Run multiple replications for each true SD value
for(sd_val in sd_values) {
  for(rep in 1:n_reps) {
    sim <- simulate_data(n_subjects = 200, true_sd_response = sd_val)
    results <- rbind(results, data.frame(
      true_sd = sd_val,
      sdir = sim$sdir,
      sd_control = sim$sd_control,
      sd_treatment = sim$sd_treatment,
      rep = rep
    ))
  }
}

# Calculate average SDir for each true SD value
avg_results <- aggregate(
  cbind(sdir, sd_control, sd_treatment) ~ true_sd, 
  data = results, 
  FUN = mean
)

print("Average results for each true SD value:")
print(avg_results)

# Run a regression of SDir on true SD to check accuracy
model <- lm(sdir ~ true_sd, data = results)
summary(model)

ggplot(aes(y=sdir, x=true_sd, group = true_sd), data = results) + geom_boxplot()
# Analyze a single simulation in detail
detailed_sim <- simulate_data(n_subjects = 500, 
                              true_sd_response = 3,
                              within_subject_sd = 3)

# Print key metrics
cat("\nDetailed simulation with true SD =", detailed_sim$true_sd, ":\n")
cat("True SD of individual effects:", detailed_sim$sd_true_effects, "\n")
cat("Estimated SDir:", detailed_sim$sdir, "\n")
cat("SD of change in control group:", detailed_sim$sd_control, "\n")
cat("SD of change in treatment group:", detailed_sim$sd_treatment, "\n")

# Gadbury's framework focuses on the "unit-treatment interaction"
# This represents the individual deviation from the average treatment effect
# Let's compute this for our detailed simulation
mean_effect <- mean(detailed_sim$true_effects)
unit_treatment_interaction <- detailed_sim$true_effects - mean_effect

cat("\nUnit-Treatment Interaction Statistics:\n")
cat("Mean:", mean(unit_treatment_interaction), "\n")  # Should be close to 0
cat("SD:", sd(unit_treatment_interaction), "\n")      # Should be close to true_sd

# Compare with variance-component approach from Gadbury
# Gadbury defined D = Y(T) - Y(C), where Y(T) and Y(C) are potential outcomes
# The variance of D has components for the average treatment effect and individual variation
# Let's calculate the theoretical variance components
var_Y_T <- var(detailed_sim$change_treatment)  # Variance of potential outcomes under treatment
var_Y_C <- var(detailed_sim$change_control)    # Variance of potential outcomes under control
cor_potential_outcomes <- cor(detailed_sim$true_followup_treatment, detailed_sim$true_followup_control)

# Gadbury's decomposition: var(D) = var(Y(T)) + var(Y(C)) - 2*cov(Y(T),Y(C))
theoretical_var_D <- var_Y_T + var_Y_C - 2*cor_potential_outcomes*sqrt(var_Y_T*var_Y_C)
cat("\nTheoretical variance of treatment effect (var(D)):", theoretical_var_D, "\n")
cat("Square root (SD of individual effects):", sqrt(theoretical_var_D), "\n")

# In our simulation, we know the true individual effects directly
empirical_var_D <- var(detailed_sim$true_effects)
cat("Empirical variance of individual effects:", empirical_var_D, "\n")
cat("Square root (SD of individual effects):", sqrt(empirical_var_D), "\n")

# Compare with SDir
cat("\nEstimated SDir:", detailed_sim$sdir, "\n")
cat("True SD of individual effects:", detailed_sim$sd_true_effects, "\n")