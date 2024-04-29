rm(list=ls());gc()
graphics.off()
# loads package
library(sensemakr)
library(ggplot2)
library(dplyr)
library(tidyr)

invasion_growth <- function(t, r, i0) {
  i0+r*(t-t0)
}

sigmoid_function <- function(t, yoi, beta, rate, p0) {
  1-exp(-rate * (t - yoi)^beta - p0)
}

calculate_Fb <- function(a, I, pd_beta, r_pd, pr_beta, r_pr, pd0, pr0, t0, b) {
  Fb <- 0
  for (i in t0:a) {
    pd_product <- if (i < a) {
      prod(1 - sapply(i:(a-1), sigmoid_function, i, pd_beta, r_pd, pd0))
    } else {
      1
    }
    pd_current <- sigmoid_function(a, i, pd_beta, r_pd, pd0)
    pr_product <- 1-prod(1 - sapply(a:b, sigmoid_function, a, pr_beta, r_pr, pr0))
    
    Fb <- Fb + I[i - t0 + 1] * pd_product * pd_current * pr_product
  }
  return(Fb)
}

predict_F <- function(i0, r, pd_beta, r_pd, pr_beta, r_pr, pd0, pr0, t0, b){
  cat("Start predicting--------------\n")
  years <- t0:b
  I <- sapply(years, invasion_growth, r, i0)
  # Calculate Fb_values 
  Fb_values <- sapply(years, function(a) {
    Fb <- calculate_Fb(a, I, pd_beta, r_pd, pr_beta, r_pr, pd0, pr0, t0, b)
    # cat("Year:", a, "Fb:", Fb, "\n")
    return(Fb)
  })
  cat("Finish predicting-------------\n")
  result <- data.frame(year = years, pre = Fb_values, invasion = I)
  return(result)
}


s_analysis = function(param_ranges, diff_obs, t0, i){
  # Generate random parameters
  # Use these parameters to predict F values for b1 and b2
  params = c(i0 = param_ranges$i0[i], r = param_ranges$r[i],
             pd_beta = pd_beta, r_pd = r_pd, 
             pr_beta = pr_beta, r_pr = r_pr, 
             pd0 = pd0, pr0 = pr0)
  results_b1 <- predict_F(i0 = param_ranges$i0[i], r = param_ranges$r[i],
                          pd_beta = pd_beta, r_pd = r_pd, 
                          pr_beta = pr_beta, r_pr = r_pr, 
                          pd0 = pd0, pr0 = pr0, 
                          t0 = t0, b = 2016)
  
  results_b2 <- predict_F(i0 = param_ranges$i0[i], r = param_ranges$r[i],
                          pd_beta = pd_beta, r_pd = r_pd, 
                          pr_beta = pr_beta, r_pr = r_pr, 
                          pd0 = pd0, pr0 = pr0, 
                          t0 = t0, b = 2021)
  
  # Compute the difference in F values and the correlation with actual invasion
  common_years <- intersect(results_b1$year, results_b2$year)
  diff_F <- results_b2$pre[common_years - t0 + 1] - results_b1$pre[common_years - t0 + 1]
  y = diff_F[396:417]
  x = common_years[396:417]
  m = lm(y ~ x)
  slp = coef(m)["x"]
  list(parameters = params, slope = slp, result_b1 = results_b1, results_b2 = results_b2, diff_F = diff_F)
}


param_ranges <- list(
  i0 = c(417, 417, 417, 417, 417),     
  r = c(-1, -0.5, 0, 0.5, 1)    # Replace with actual ranges
  # Replace with actual ranges
)


t0 = 1600
b1 = 2016
b2 = 2021

pd_beta <- 1.599110e+00
r_pd <- 0.000000e+00 
pd0 = 5.407369e-03
pr_beta <- 1.368035e+00
r_pr <- 1.438073e-02
pr0 = 8.712529e-03


results = list()
for (i in 1:length(param_ranges$i0)){
  res = s_analysis(param_ranges, diff_obs, t0, i)
  results[[i]] = res
}

I_slopes = param_ranges$r
dF_slopes = c(results[[1]]$slope, results[[2]]$slope, results[[3]]$slope, results[[4]]$slope, results[[5]]$slope)
