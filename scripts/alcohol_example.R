# Reproduction from: Gary L. Gadbury, Hari K. Iyer & David B. Allison (2001) EVALUATING
# SUBJECT-TREATMENT INTERACTION WHEN COMPARING TWO TREATMENTS, Journal of
# Biopharmaceutical Statistics, 11:4, 313-333, DOI: 10.1081/BIP-120008851
# link to this article: https://doi.org/10.1081/BIP-120008851
# Example usage with the alcohol intake data from Example 1
alcohol_example <- function() {
  # Data from Table 1
  sst <- c(874, 389, 612, 798, 1152, 893, 541, 741, 1064, 862, 213)
  control <- c(1042, 1617, 1180, 973, 1552, 1251, 1151, 1511, 728, 1079, 951, 1319)
  
  # Run sensitivity analysis
  sens_mle <- sensitivity_analysis(sst, control, method = "mle")
  sens_boot <- sensitivity_analysis(sst, control, method = "bootstrap")
  
  return(list(mle_results = sens_mle, bootstrap_results = sens_boot))
}

test = alcohol_example()

library(tidyverse)
# plot MLE results

test$mle_results %>%
  ggplot(aes(x=rho,y=sigma_d,
             ymin = sigma_d_lower,
             ymax = sigma_d_upper)) +
  geom_hline(yintercept = 0, color="darkred", alpha = .8) +
  geom_ribbon(fill = "grey",
              alpha = .2) +
  geom_line(linetype=4, #color = "red",
            size = 1.1) +
  ggprism::theme_prism()


test$mle_results %>%
  ggplot(aes(x=rho,y=p_minus,
             ymin = p_minus_lower,
             ymax = p_minus_upper)) +
  geom_hline(yintercept = 0, color="darkred", alpha = .8) +
  geom_ribbon(fill = "grey",
              alpha = .2) +
  geom_line(linetype=4, #color = "red",
            size = 1.1) +
  ggprism::theme_prism()
