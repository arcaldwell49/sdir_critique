# Simulation to test SDir using the potential outcomes framework
# This allows specifying correlation between potential outcomes directly
# rather than specifying the SD of individual responses

# Set seed for reproducibility
set.seed(12345)
#library(simstudy)
# Function to generate data using the potential outcomes framework
# Parameters:
#   n_subjects: number of subjects
#   mean_effect: mean treatment effect 
#   cor_potential_outcomes: correlation between Y(1) and Y(0)
#   sd_potential_outcomes: standard deviation of potential outcomes
#   within_subject_sd: additional within-subject random error
simulate_potential_outcomes <- function(n_subjects = 200, 
                                        mean_effect = 5, 
                                        cor_potential_outcomes = 0.9,
                                        sd_potential_outcomes = 8,
                                        within_subject_sd = 3) {
  
  # Generate correlated potential outcomes Y(0) and Y(1)
  # Here Y(0) is the potential outcome under control
  # Y(1) is the potential outcome under treatment
  #n_subjects2 = n_subjects*2
  # First, generate uncorrelated normal random variables
  z1 <- rnorm(n_subjects, mean = 0, sd = 1)
  z2 <- rnorm(n_subjects, mean = 0, sd = 1)
  
  # Create correlation by combining the uncorrelated variables
  # Y(0) = standard normal
  # Y(1) = cor*Y(0) + sqrt(1-cor^2)*z2
  y0 <- z1 * sd_potential_outcomes
  y1 <- cor_potential_outcomes * y0 + 
    sqrt(1 - cor_potential_outcomes^2) * z2 * sd_potential_outcomes + 
    mean_effect
  # These are the true potential outcomes
  # The individual treatment effects are Y(1) - Y(0)
  individual_effects <- y1 - y0
  # Calculate the theoretical SD of individual effects
  # Based on Gadbury's formula: var(D) = var(Y(1)) + var(Y(0)) - 2*cov(Y(1),Y(0))
  var_y0 <- var(y0)
  var_y1 <- var(y1)
  cov_y0_y1 <- cov(y0, y1)
  cor_y0_y1 <- cor(y0, y1)
  theoretical_var_D <- var_y1 + var_y0 - 2 * cov_y0_y1
  sample_theoretical_sd_D <- sqrt(theoretical_var_D)
  theoretical_sd_D = sqrt(sd_potential_outcomes^2 + sd_potential_outcomes^2 - 
    2 * sd_potential_outcomes * sd_potential_outcomes * cor_potential_outcomes)
  # Now add within-subject variability to create observed values
  # Each person has two observations in each condition (baseline and follow-up)
  # assume the true value doesn't change in control condition
  # obs_baseline_control <- y0 + rnorm(n_subjects, 0, within_subject_sd)
  # obs_followup_control <- y0 + rnorm(n_subjects, 0, within_subject_sd)
  # 
  # obs_baseline_treatment <- y0 + rnorm(n_subjects, 0, within_subject_sd)
  # obs_followup_treatment <- y1 + rnorm(n_subjects, 0, within_subject_sd)
  # Randomly assign subjects to treatment or control
  treatment_assignment <- sample(c(0,1), n_subjects, replace=TRUE)
  
  # Observe baseline for all subjects (assuming baseline is based on y0)
  obs_baseline <- y0 + rnorm(n_subjects, 0, within_subject_sd)
  
  # Observe follow-up based on assignment
  obs_followup <- numeric(n_subjects)
  for(i in 1:n_subjects) {
    if(treatment_assignment[i] == 1) {
      # Subject received treatment
      obs_followup[i] <- y1[i] + rnorm(1, 0, within_subject_sd)
    } else {
      # Subject received control
      obs_followup[i] <- y0[i] + rnorm(1, 0, within_subject_sd)
    }
  }
  
  # Calculate changes for each group
  change_control <- obs_followup[treatment_assignment == 0] - 
    obs_baseline[treatment_assignment == 0]
  change_treatment <- obs_followup[treatment_assignment == 1] - 
    obs_baseline[treatment_assignment == 1]
  # # Calculate observed changes
  # change_control <- obs_followup_control - obs_baseline_control
  # change_treatment <- obs_followup_treatment - obs_baseline_treatment 
  #cor(obs_followup_control, obs_baseline_control)
  # Calculate SDir (as per Atkinson & Batterham)
  sd_change_control <- sd(change_control)
  sd_change_treatment <- sd(change_treatment)
  sdir <- sign(sd_change_treatment^2 - sd_change_control^2) * sqrt(abs(sd_change_treatment^2 - sd_change_control^2))
  
  # Return the results
  list(
    cor_potential_outcomes = cor_potential_outcomes,
    empirical_cor = cor(y0, y1),
    true_effects = individual_effects,
    sd_true_effects = sd(individual_effects),
    theoretical_sd_effects = theoretical_sd_D,
    sdir = sdir,
    sd_control = sd_change_control,
    sd_treatment = sd_change_treatment,
    change_control = change_control,
    change_treatment = change_treatment,
    potential_outcome_control = y0,
    potential_outcome_treatment = y1
  )
}

# Run simulations with varying correlations between potential outcomes
# Note: Higher correlation = less heterogeneity in treatment effects
correlation_values <- seq(from = 1.0, to = -1.0, by = -0.1)
n_reps <- 100
results <- data.frame(
  correlation = numeric(),
  sdir = numeric(),
  sd_true_effects = numeric(),
  theoretical_sd_effects = numeric(),
  sd_control = numeric(),
  sd_treatment = numeric(),
  rep = numeric()
)

# Run multiple replications for each correlation
for(cor_val in correlation_values) {
  for(rep in 1:n_reps) {
    sim <- simulate_potential_outcomes(
      n_subjects = 20, 
      cor_potential_outcomes = cor_val,
      sd_potential_outcomes = 8,
      within_subject_sd = 3
    )
    
    results <- rbind(results, data.frame(
      correlation = cor_val,
      sdir = sim$sdir,
      sd_true_effects = sim$sd_true_effects,
      theoretical_sd_effects = sim$theoretical_sd_effects,
      sd_control = sim$sd_control,
      sd_treatment = sim$sd_treatment,
      rep = rep
    ))
  }
}

# Calculate average results for each correlation value
avg_results <- aggregate(
  cbind(sdir, sd_true_effects, theoretical_sd_effects, sd_control, sd_treatment) ~ correlation, 
  data = results, 
  FUN = mean
)

print("Average results for each correlation between potential outcomes:")
print(avg_results)

# Let's visualize the relationship between correlation and heterogeneity
# Run a specific example for detailed analysis
cor_examples <- c(1.0, 0.9, 0.7, 0.5, 0.3)
detailed_examples <- list()

for(cor_val in cor_examples) {
  detailed_examples[[as.character(cor_val)]] <- simulate_potential_outcomes(
    n_subjects = 500, 
    cor_potential_outcomes = cor_val,
    sd_potential_outcomes = 8,
    within_subject_sd = 3
  )
}

# Print key metrics for each correlation example
for(cor_val in cor_examples) {
  sim <- detailed_examples[[as.character(cor_val)]]
  cat("\n===== Correlation between potential outcomes =", cor_val, "=====\n")
  cat("Empirical correlation:", sim$empirical_cor, "\n")
  cat("Theoretical SD of individual effects:", sim$theoretical_sd_effects, "\n")
  cat("Empirical SD of individual effects:", sim$sd_true_effects, "\n")
  cat("Estimated SDir:", sim$sdir, "\n")
  cat("SD of change in control group:", sim$sd_control, "\n")
  cat("SD of change in treatment group:", sim$sd_treatment, "\n\n")
}

# Demonstrate the relationship between correlation and the SD of individual effects
# When correlation = 1, there should be no heterogeneity (all individual effects equal)
# When correlation < 1, heterogeneity emerges
# When correlation = 0, maximum heterogeneity given the SD of potential outcomes

# The relationship between correlation and the SD of individual effects is:
# SD(D) = sqrt(var(Y1) + var(Y0) - 2*cov(Y1,Y0))
# For standardized variables with variances = 1:
# SD(D) = sqrt(2*(1-cor(Y1,Y0)))

# Display this relationship
correlations <- seq(from = 0, to = 1, by = 0.01)
theoretical_sds <- sqrt(2*(1-correlations)) * 8  # Scale by SD of potential outcomes

print("Theoretical relationship between correlation and SD of individual effects:")
head_correlations <- head(data.frame(correlation = correlations, sd_effects = theoretical_sds), 10)
tail_correlations <- tail(data.frame(correlation = correlations, sd_effects = theoretical_sds), 10)
print(head_correlations)
print(tail_correlations)