# Load required packages
library(metafor)  # Required for log_vr_test
library(tidyverse)
library(here)
# Import other functions -----
# # Source the functions from var_functions.R
source(
  here("functions",
       "var_functions.R")
)




# Define simulation parameters
set.seed(1914)  # For reproducibility
n_sims <- 10000  # Number of simulations for each condition
sample_sizes <- seq(10, 100, by = 10)  # Sample sizes to test
variance_ratios <- seq(1, 3, by = 0.25)  # Variance ratios to test
alpha <- 0.05  # Significance level

# Initialize results matrix
results <- expand.grid(
  sample_size = sample_sizes,
  variance_ratio = variance_ratios,
  test = c("diff_sdir_test", "log_vr_test", "var.test"),
  power = NA
)

# Run simulations
pb <- txtProgressBar(min = 0, max = nrow(results), style = 3)
counter <- 0

for (i in 1:nrow(results)) {
  n <- results$sample_size[i]
  var_ratio <- results$variance_ratio[i]
  test_type <- results$test[i]
  
  # Track rejections
  rejections <- 0
  
  for (sim in 1:n_sims) {
    # Generate data
    # Group 1 with variance 1
    sd1 <- 1
    # Group 2 with variance based on ratio
    sd2 <- sqrt(1/var_ratio)
    
    # Generate samples
    group1 <- rnorm(n, mean = 0, sd = sd1)
    group2 <- rnorm(n, mean = 0, sd = sd2)
    
    # Compute test statistics
    p_value <- switch(
      test_type,
      "diff_sdir_test" = diff_sdir_test(sd(group1), n, sd(group2), n, alternative = "two.sided")$p.value,
      "log_vr_test" = log_vr_test(sd(group1), n, sd(group2), n, alternative = "two.sided")$p.value,
      "var.test" = var.test(group1, group2)$p.value
    )
    
    # Count rejections
    if (p_value < alpha) {
      rejections <- rejections + 1
    }
  }
  
  # Calculate power
  results$power[i] <- rejections / n_sims
  
  # Update progress bar
  counter <- counter + 1
  setTxtProgressBar(pb, counter)
}

close(pb)

# Reshape results for easier plotting
results_wide <- reshape(results, 
                        idvar = c("sample_size", "variance_ratio"),
                        timevar = "test", 
                        direction = "wide") %>%
  mutate(difference_f_sdir = power.var.test - power.diff_sdir_test)

# Save results
saveRDS(results, "variance_test_power_results.rds")

# Generate plots
library(ggplot2)

# Plot 1: Power by sample size for each test at different variance ratios
ggplot(results, aes(x = sample_size, y = power, color = test)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ variance_ratio, labeller = label_both) +
  labs(
    title = "Power of Variance Tests by Sample Size",
    x = "Sample Size per Group",
    y = "Power (1-β)",
    color = "Test"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(here::here("figures","power_by_sample_size.png"), width = 10, height = 8)

# Plot 2: Power by variance ratio for each test at different sample sizes
ggplot(results, aes(x = variance_ratio, y = power, color = test)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ sample_size, labeller = label_both) +
  labs(
    title = "Power of Variance Tests by Variance Ratio",
    x = "Variance Ratio",
    y = "Power (1-β)",
    color = "Test"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(here::here("figures","power_by_variance_ratio.png"), width = 10, height = 8)


# Plot 3: Difference in power (F - SDir) at different sample sizes
ggplot(results_wide, aes(x = sample_size, y = difference_f_sdir)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ variance_ratio, labeller = label_both) +
  labs(
    title = "Difference in Power between F-test and SDir",
    x = "Sample Size",
    y = "Difference in Power (Variance Ratio F-test - SDir",
    color = "Test"
  ) +
  geom_hline(yintercept = 0)+
  theme_bw() +
  scale_y_continuous(labels = scales::percent)+
  theme(legend.position = "bottom")

ggsave(here::here("figures","power_difference.png"), width = 10, height = 8)

# Create heatmaps for each test
for (test_name in unique(results$test)) {
  test_results <- subset(results, test == test_name)
  
  ggplot(test_results, aes(x = sample_size, y = variance_ratio, fill = power)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(
      title = paste("Power of", test_name),
      x = "Sample Size per Group",
      y = "Variance Ratio",
      fill = "Power"
    ) +
    theme_minimal()
  
  ggsave(here::here("figures",paste0("power_heatmap_", test_name, ".png")), width = 7, height = 6)
}

