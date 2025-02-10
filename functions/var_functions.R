p_from_z = function(x,
                    alternative = "two.sided"){
  if(alternative == "two.sided"){
    2*pnorm(-abs(unlist(x)))
  } else  if (alternative  == "greater"){
    pnorm(x, lower.tail = FALSE)
  } else if (alternative  == "less"){
    pnorm(x, lower.tail = TRUE)
  } else{
    stop("alternative must be two.sided, greater, or less")
  }
  
}

# SDir (normal based test from Hopkins variances estimate)
diff_sdir_test = function(sd1, n1,
                          sd2, n2,
                          alternative = c("two.sided",
                                          "greater",
                                          "less")){
  
  alternative = match.arg(alternative)
  var1 = sd1^2
  var1_se <- var1*(sqrt(2/(n1 - 1)))
  var2 = sd2^2
  var2_se <- var2*(sqrt(2/(n2 - 1)))
  df1 = n1 -1
  df2 = n2 - 1
  if((var1-var2)> 0){
    sd_ir = sqrt(var1-var2)
    
  } else{
    sd_ir = -1*sqrt(var2-var1)
    
  }
  est_diff <- var1-var2
  
  est_diff_SE <- sqrt(2 * (sd1^4/df1 + sd2^4/df2)) 
  teststat <- est_diff/est_diff_SE
  pval = p_from_z(teststat, alternative = alternative)
  
  statistic = teststat
  names(statistic) = "z"
  null1 = 0
  names(null1) = "Standard Deviation of the Individual Response"
  estimate = sd_ir
  names(estimate)= "Standard Deviation of the Individual Response"
  std_err = sqrt(est_diff_SE)
  sum_stats = paste0("SD1 = ", sd1, ", SD2 = ", sd2)
  
  rval <- list(statistic = statistic,
               p.value = as.numeric(pval),
               estimate = estimate,
               stderr = std_err,
               null.value = null1,
               alternative = alternative,
               method = "Difference in Variances",
               data.name = sum_stats)
  class(rval) <- "htest"
  return(rval)
}
# ratio of variances (normal based test?) logVR

log_vr_test = function(sd1, n1,
                       sd2, n2,
                       alternative = c("two.sided",
                                       "greater",
                                       "less")){
  
  # old code
  # var1 = sd1^2
  # var1_se <- var1*(sqrt(2/(n1 - 1)))
  # var2 = sd2^2
  # var2_se <- var2*(sqrt(2/(n2 - 1)))
  # yi <- log(sqrt(var1)/sqrt(var2)) + 1/(2 * (n1 - 1)) - 1/(2 * 
  #                                                                     (n2 - 1))
  # 
  # vi <- 1/(2 * (n1 - 1)) + 1/(2 * (n2 - 1))
  # 
  # teststat <- yi/sqrt(vi)
  alternative = match.arg(alternative)
  # metafor code
  es_est = metafor::escalc(
    measure = "VR",
    n1i = n1,
    n2i = n2,
    sd1i = sd1,
    sd2i = sd2
  )
  
  teststat <- es_est$yi/sqrt(es_est$vi)
  pval = p_from_z(teststat, 
                  alternative = alternative)
  
  statistic = teststat
  names(statistic) = "z"
  null1 = 0
  names(null1) = "log ratio of SD ratio"
  estimate = es_est$yi
  names(estimate)= "log VR"
  std_err = sqrt(es_est$vi)
  sum_stats = paste0("SD1 = ", sd1, ", SD2 = ", sd2)
  
  rval <- list(statistic = statistic,
               p.value = as.numeric(pval),
               estimate = estimate,
               stderr = std_err,
               null.value = null1,
               alternative = alternative,
               method = "log transformed variability ratio",
               data.name = sum_stats)
  class(rval) <- "htest"
  return(rval)
}
