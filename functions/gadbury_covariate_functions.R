#' Prepare data for covariate analysis
#' @param data Data frame containing the response, treatment indicator, and covariate
#' @param response Name of response variable
#' @param treatment Name of treatment indicator (1=treatment, 0=control)
#' @param covariate Name of covariate
#' @return List containing reformatted data and model fits
prepare_covariate_data <- function(data, response, treatment, covariate) {
  # Ensure inputs are character strings
  response <- as.character(substitute(response))
  treatment <- as.character(substitute(treatment))
  covariate <- as.character(substitute(covariate))
  
  # Fit main effect and interaction models
  formula_main <- as.formula(paste(response, "~", treatment, "+", covariate))
  formula_int <- as.formula(paste(response, "~", treatment, "*", covariate))
  
  fit_main <- lm(formula_main, data = data)
  fit_int <- lm(formula_int, data = data)
  
  # Split data into treatment and control groups
  treat_data <- subset(data, treatment == 1)
  ctrl_data <- subset(data, treatment == 0)
  
  # Fit separate models for each group
  formula_sep <- as.formula(paste(response, "~", covariate))
  fit_treat <- lm(formula_sep, data = treat_data)
  fit_ctrl <- lm(formula_sep, data = ctrl_data)
  
  list(
    data = data,
    treat_data = treat_data,
    ctrl_data = ctrl_data,
    fit_main = fit_main,
    fit_int = fit_int,
    fit_treat = fit_treat,
    fit_ctrl = fit_ctrl,
    vars = list(
      response = response,
      treatment = treatment,
      covariate = covariate
    )
  )
}

#' Estimate conditional treatment effects using linear models
#' @param model_data Output from prepare_covariate_data
#' @param z0 Covariate value at which to evaluate conditional effects
#' @param rho_xy_z Partial correlation between potential outcomes given Z
#' @param conf.level Confidence level for intervals
#' @return List containing conditional estimates and confidence intervals
estimate_conditional_effects_lm <- function(model_data, z0, rho_xy_z, conf.level = 0.95) {
  # Extract relevant components
  fit_treat <- model_data$fit_treat
  fit_ctrl <- model_data$fit_ctrl
  treat_data <- model_data$treat_data
  ctrl_data <- model_data$ctrl_data
  vars <- model_data$vars
  
  # Get residual standard deviations
  s_x_z <- summary(fit_treat)$sigma
  s_y_z <- summary(fit_ctrl)$sigma
  
  # Get regression coefficients
  beta_x <- coef(fit_treat)[2]  # slope for treatment group
  beta_y <- coef(fit_ctrl)[2]   # slope for control group
  
  # Calculate predicted values at z0
  new_data <- data.frame(
    temp = z0
  )
  names(new_data) <- vars$covariate
  
  pred_x <- predict(fit_treat, newdata = new_data, se.fit = TRUE)
  pred_y <- predict(fit_ctrl, newdata = new_data, se.fit = TRUE)
  
  # Calculate conditional mean treatment effect
  mu_d_z <- pred_x$fit - pred_y$fit
  
  # Calculate standard error of conditional mean effect
  se_mu_d_z <- sqrt(pred_x$se.fit^2 + pred_y$se.fit^2)
  
  # Calculate conditional variance of treatment effects
  sigma_d_z_sq <- s_x_z^2 + s_y_z^2 - 2 * s_x_z * s_y_z * rho_xy_z
  
  # Calculate variance of sigma_d_z estimate
  n1 <- nrow(treat_data)
  n2 <- nrow(ctrl_data)
  
  var_sigma_d_z_sq <- 2 * (
    (s_x_z^2/n1) * (s_x_z - rho_xy_z * s_y_z)^2 +
      (s_y_z^2/n2) * (s_y_z - rho_xy_z * s_x_z)^2
  )
  
  # Calculate confidence intervals
  z_crit <- qnorm(1 - (1 - conf.level)/2)
  
  mu_d_z_ci <- c(
    mu_d_z - z_crit * se_mu_d_z,
    mu_d_z + z_crit * se_mu_d_z
  )
  
  sigma_2_d_z_ci <- (c(
    sigma_d_z_sq - z_crit * sqrt(var_sigma_d_z_sq),
    sigma_d_z_sq + z_crit * sqrt(var_sigma_d_z_sq)
  ))
  
  sigma_d_z_ci <- sqrt(abs(sigma_2_d_z_ci)) * sign(sigma_2_d_z_ci)
  
  # Test for interaction
  int_test <- anova(model_data$fit_main, model_data$fit_int)
  
  return(list(
    mu_d_z = mu_d_z,
    mu_d_z_ci = mu_d_z_ci,
    mu_d_z_se = se_mu_d_z,
    sigma_d_z = sqrt(sigma_d_z_sq),
    sigma_d_z_ci = sigma_d_z_ci,
    beta_x = beta_x,
    beta_y = beta_y,
    s_x_z = s_x_z,
    s_y_z = s_y_z,
    interaction_test = int_test,
    model_summary = list(
      treatment = summary(fit_treat),
      control = summary(fit_ctrl),
      interaction = summary(model_data$fit_int)
    )
  ))
}

#' Calculate conditional probability of unfavorable effect using linear models
#' @param model_data Output from prepare_covariate_data
#' @param z0 Covariate value at which to evaluate
#' @param rho_xy_z Partial correlation between potential outcomes given Z
#' @param conf.level Confidence level for intervals
#' @return List containing P-(Z=z0) estimate and confidence interval
estimate_conditional_p_minus_lm <- function(model_data, z0, rho_xy_z, conf.level = 0.95) {
  # Get conditional estimates
  cond_est <- estimate_conditional_effects_lm(model_data, z0, rho_xy_z, conf.level)
  
  # Calculate P-(Z=z0)
  p_minus_z <- pnorm(cond_est$mu_d_z / cond_est$sigma_d_z)
  
  # Calculate variance components
  phi_term <- dnorm(cond_est$mu_d_z / cond_est$sigma_d_z)
  var_p_minus_z <- (phi_term^2 / cond_est$sigma_d_z^2) * (
    cond_est$mu_d_z_se^2 + 
      (cond_est$mu_d_z^2 * var(cond_est$sigma_d_z_ci^2)) / (4 * cond_est$sigma_d_z^4)
  )
  
  # Calculate confidence intervals
  z_crit <- qnorm(1 - (1 - conf.level)/2)
  ci_lower <- p_minus_z - z_crit * sqrt(var_p_minus_z)
  ci_upper <- p_minus_z + z_crit * sqrt(var_p_minus_z)
  
  return(list(
    p_minus_z = p_minus_z,
    ci = c(ci_lower, ci_upper)
  ))
}

#' Covariate-based sensitivity analysis using linear models
#' @param model_data Output from prepare_covariate_data
#' @param z0 Covariate value at which to evaluate
#' @param rho_seq Sequence of partial correlation values
#' @param conf.level Confidence level for intervals
#' @return Data frame of conditional estimates across rho values
covariate_sensitivity_analysis_lm <- function(model_data, z0, 
                                              rho_seq = seq(-1, 1, by = 0.1),
                                              conf.level = 0.95) {
  
  results <- lapply(rho_seq, function(rho) {
    cond_effects <- estimate_conditional_effects_lm(model_data, z0, rho, conf.level)
    p_minus <- estimate_conditional_p_minus_lm(model_data, z0, rho, conf.level)
    
    data.frame(
      rho = rho,
      mu_d_z = cond_effects$mu_d_z,
      mu_d_z_lower = cond_effects$mu_d_z_ci[1],
      mu_d_z_upper = cond_effects$mu_d_z_ci[2],
      sigma_d_z = cond_effects$sigma_d_z,
      sigma_d_z_lower = cond_effects$sigma_d_z_ci[1],
      sigma_d_z_upper = cond_effects$sigma_d_z_ci[2],
      p_minus_z = p_minus$p_minus_z,
      p_minus_z_lower = p_minus$ci[1],
      p_minus_z_upper = p_minus$ci[2]
    )
  })
  
  do.call(rbind, results)
}

# Example usage with blood pressure data
bp_example_lm <- function() 
  {
  # Create data frame from Example 2
  bp_data <- data.frame(
    change = c(-7, 4, -18, -17, 3, 5, -1, -10, -11, 2,  # Treatment
               1, -12, 1, 3, -3, 5, -5, 2, 11, 1, 3),   # Control
    baseline = c(107, 110, 123, 129, 112, 111, 107, 112, 136, 102,  # Treatment
                 123, 109, 112, 102, 98, 114, 119, 112, 110, 117, 130), # Control
    treatment = c(rep(1, 10), rep(0, 11))
  )
  
  # Prepare data and models
  model_data <- prepare_covariate_data(bp_data, "change", "treatment", "baseline")
  
  # Run sensitivity analysis at mean baseline (z=114)
  sens_results <- covariate_sensitivity_analysis_lm(model_data, 114)
  
  # Run sensitivity analysis at high baseline (z=130)
  sens_results_high <- covariate_sensitivity_analysis_lm(model_data, 130)
  
  return(list(
    model_data = model_data,
    z114_results = sens_results,
    z130_results = sens_results_high
  ))
}

test = bp_example_lm()

# Check plots... not sure this is working yet...
test$z114_results %>%
  ggplot(aes(x=rho,y=sigma_d_z,
             ymin = sigma_d_z_lower,
             ymax = sigma_d_z_upper)) +
  geom_hline(yintercept = 0, color="darkred", alpha = .8) +
  geom_ribbon(fill = "grey",
              alpha = .2) +
  geom_line(linetype=4, #color = "red",
            size = 1.1) +
  ggprism::theme_prism()


test$z130_results %>%
  ggplot(aes(x=rho,y=sigma_d_z,
             ymin = sigma_d_z_lower,
             ymax = sigma_d_z_upper)) +
  geom_hline(yintercept = 0, color="darkred", alpha = .8) +
  geom_ribbon(fill = "grey",
              alpha = .2) +
  geom_line(linetype=4, #color = "red",
            size = 1.1) +
  ggprism::theme_prism()

