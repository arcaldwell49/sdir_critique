# Functions for evaluating subject-treatment interaction between two treatments

#' Calculate MLEs and confidence intervals for treatment effect variance
#' @param treatment Numeric vector of observations from treatment group
#' @param control Numeric vector of observations from control group
#' @param rho_xy Correlation between potential outcomes (-1 to 1)
#' @param conf.level Confidence level for intervals (default 0.95)
#' @return List containing estimates and confidence intervals
estimate_sigma_d <- function(treatment, 
                             control, 
                             rho_xy, 
                             conf.level = 0.95) {
  # Sample sizes
  n1 <- length(treatment)
  n2 <- length(control)
  
  # Calculate sample statistics
  sigma_x_hat <- sd(treatment)
  sigma_y_hat <- sd(control)
  
  # Estimate sigma_d^2
  sigma_d_sq_hat <- sigma_x_hat^2 + sigma_y_hat^2 - 
    2 * sigma_x_hat * sigma_y_hat * rho_xy
  
  # Calculate variance of sigma_d^2 estimate (equation 4 from paper)
  var_sigma_d_sq <- 2 * (
    (sigma_x_hat^2/n1) * (sigma_x_hat - rho_xy * sigma_y_hat)^2 +
      (sigma_y_hat^2/n2) * (sigma_y_hat - rho_xy * sigma_x_hat)^2
  )
  
  # Calculate confidence intervals
  z_crit <- qnorm(1 - (1 - conf.level)/2)
  ci_lower <- sigma_d_sq_hat - z_crit * sqrt(var_sigma_d_sq)
  ci_upper <- sigma_d_sq_hat + z_crit * sqrt(var_sigma_d_sq)
  
  # Take square root for sigma_d
  sigma_d_hat <- sqrt(sigma_d_sq_hat)
  # remove sign errors
  sigma_d_ci <- sqrt(abs(c(ci_lower, ci_upper))) * sign(c(ci_lower,ci_upper))
  
  return(list(
    sigma_d = sigma_d_hat,
    sigma_d_ci = sigma_d_ci,
    sigma_d_sq = sigma_d_sq_hat,
    sigma_d_sq_ci = c(ci_lower, ci_upper)
  ))
}

#' Calculate probability of unfavorable treatment effect
#' @param treatment Numeric vector of observations from treatment group
#' @param control Numeric vector of observations from control group
#' @param ATE Average treatment effect. A mean effect, like those from a ANCOVA, can be supplied in lieu of using the raw data
#' @param rho_xy Correlation between potential outcomes (-1 to 1)
#' @param conf.level Confidence level for intervals (default 0.95)
#' @param lower.tail determines whether the area under the normal distribution curve is returned to the left or right of a specified value.
#'  Default is FALSE which infers higher scores in treatment group are positive (i.e., positivie ATE is good)
#' @return List containing P- estimate and confidence interval
estimate_p_minus <- function(treatment, 
                             control, 
                             ATE = NULL,
                             rho_xy, 
                             conf.level = 0.95,
                             lower.tail = FALSE) {
  # Sample sizes
  n1 <- length(treatment)
  n2 <- length(control)
  lt = lower.tail
  # Calculate sample statistics
  sigma_x_hat <- sd(treatment)
  sigma_y_hat <- sd(control)
  # Set the average effect of treatment
  if(is.null(ATE)){
    mu_d_hat <- mean(treatment) - mean(control)
  } else{
    mu_d_hat = ATE
    }
  
  
  # Get sigma_d estimates
  sigma_est <- estimate_sigma_d(treatment = treatment, 
                                control = control, 
                                rho_xy = rho_xy, 
                                conf.level = conf.level)
  sigma_d_hat <- sigma_est$sigma_d
  
  # Calculate P- (probability of unfavorable effect)
  p_minus_hat <- pnorm(mu_d_hat/sigma_d_hat, lower.tail = lt)
  
  
  # Calculate variance components (equation 5 from paper)
  var_mu_d <- sigma_x_hat^2/n1 + sigma_y_hat^2/n2
  phi_term <- dnorm(mu_d_hat/sigma_d_hat)
  
  var_p_minus <- (phi_term^2/sigma_d_hat^2) * (
    var_mu_d + 
      (mu_d_hat^2 * sigma_est$sigma_d_sq_ci[2])/(4 * sigma_d_hat^4)
  )
  
  # Calculate confidence intervals
  z_crit <- qnorm(1 - (1 - conf.level)/2)
  ci_lower <- p_minus_hat - z_crit * sqrt(var_p_minus)
  ci_upper <- p_minus_hat + z_crit * sqrt(var_p_minus)
  
  return(list(
    p_minus = p_minus_hat,
    ci = c(ci_lower, ci_upper)
  ))
}

#' Bootstrap version for small samples
#' @param treatment Numeric vector of observations from treatment group
#' @param control Numeric vector of observations from control group
#' @param rho_xy Correlation between potential outcomes (-1 to 1)
#' @param conf.level Confidence level for intervals (default 0.95)
#' @param B Number of bootstrap samples (default 1999)
#' @return List containing bootstrap estimates and confidence intervals
bootstrap_sigma_d <- function(treatment, 
                              control,
                              rho_xy, 
                              conf.level = 0.95, 
                              B = 1999) {
  n1 <- length(treatment)
  n2 <- length(control)
  
  # Function to compute sigma_d for a bootstrap sample
  boot_sigma_d <- function(treat_sample, ctrl_sample) {
    sigma_x <- sd(treat_sample) * sqrt(n1/(n1-1))
    sigma_y <- sd(ctrl_sample) * sqrt(n2/(n2-1))
    sigma_d_sq <- sigma_x^2 + sigma_y^2 - 2*sigma_x*sigma_y*rho_xy
    return(sqrt(sigma_d_sq))
  }
  
  # Generate bootstrap samples and compute sigma_d
  boot_estimates <- replicate(B, {
    treat_sample <- sample(treatment, size = n1, replace = TRUE)
    ctrl_sample <- sample(control, size = n2, replace = TRUE)
    boot_sigma_d(treat_sample, ctrl_sample)
  })
  
  # Calculate percentile confidence intervals
  ci <- quantile(boot_estimates, probs = c((1-conf.level)/2, 1-(1-conf.level)/2))
  
  return(list(
    sigma_d = mean(boot_estimates),
    ci = ci,
    boot_samples = boot_estimates
  ))
}

#' Sensitivity analysis across range of rho values
#' @param treatment Numeric vector of observations from treatment group
#' @param control Numeric vector of observations from control group
#' @param ATE Average treatment effect. A mean effect, like those from a ANCOVA, can be supplied in lieu of using the raw data
#' @param rho_seq Sequence of rho values to evaluate
#' @param conf.level Confidence level for intervals
#' @param method Either "mle" or "bootstrap"
#' @param B Number of bootstrap samples if method="bootstrap", default is 1999.
#' @param lower.tail determines whether the area under the normal distribution curve is returned to the left or right of a specified value.
#'  Default is FALSE which infers higher scores in treatment group are positive (i.e., positivie ATE is good)
#' @return Data frame of estimates across rho values
sensitivity_analysis <- function(treatment, control, 
                                 ATE = NULL,
                                 rho_seq = seq(-1, 1, by = 0.1),
                                 conf.level = 0.95,
                                 method = "mle",
                                 B = 1999,
                                 lower.tail = FALSE) {
  
  results <- lapply(rho_seq, function(rho) {
    if(method == "mle") {
      sigma_d <- estimate_sigma_d(treatment = treatment, 
                                  control = control, 
                                  rho_xy = rho, 
                                  conf.level = conf.level)
      p_minus <- estimate_p_minus(treatment = treatment, 
                                  control = control, ATE = ATE,
                                  rho_xy = rho, 
                                  conf.level = conf.level,
                                  lower.tail = lower.tail)
      
      data.frame(
        rho = rho,
        sigma_d = sigma_d$sigma_d,
        sigma_d_lower = sigma_d$sigma_d_ci[1],
        sigma_d_upper = sigma_d$sigma_d_ci[2],
        p_minus = p_minus$p_minus,
        p_minus_lower = p_minus$ci[1],
        p_minus_upper = p_minus$ci[2]
      )
    } else {
      boot_results <- bootstrap_sigma_d(treatment, control, rho, conf.level, B)
      data.frame(
        rho = rho,
        sigma_d = boot_results$sigma_d,
        sigma_d_lower = boot_results$ci[1],
        sigma_d_upper = boot_results$ci[2]
      )
    }
  })
  
  do.call(rbind, results)
}
