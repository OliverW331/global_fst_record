rm(list=ls());gc()
graphics.off()
# loads package
library(sensemakr)
library(ggplot2)
library(dplyr)
library(tidyr)
v1.2 = read.csv("../data/GlobalAlienSpeciesFirstRecordDatabase_v1.2.csv")
v2 = read.csv("../data/GlobalAlienSpeciesFirstRecordDatabase_v2.csv")
v3.1 = read.csv("../data/GlobalAlienSpeciesFirstRecordDatabase_v3.1_freedata.csv")

tr_df = function(v){
  df = v %>%
    filter(FirstRecord >= t0)
  
  species_count_df <- df %>%
    group_by(FirstRecord) %>%
    summarise(Count = n(), .groups = 'drop')
  
  b <- max(species_count_df$FirstRecord)
  all_years_df <- data.frame(FirstRecord = t0:b)
  
  species_count_df <- merge(all_years_df, species_count_df, by = "FirstRecord", all.x = TRUE)
  
  species_count_df$Count[is.na(species_count_df$Count)] <- 0
  return(species_count_df)
}

invasion_growth <- function(t, K, r, t_m, r_d, t_d) {
  if (t<t_d) {
    K / (1 + exp(-r * (t - t_m)))
  }else{
    result = K / (1 + exp(-r * (t - t_m))) - exp(r_d * (t-t_d))
    if (result < 0){
      0
    }else{
      result
    }
  }
}

# Define the sigmoid functions for Pd and Pr
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

predict_F <- function(K, r, t_m, r_d, t_d, pd_beta, r_pd, pr_beta, r_pr, pd0, pr0, t0, b){
  cat("Start predicting--------------\n")
  years <- t0:b
  I <- sapply(years, invasion_growth, K, r, t_m, r_d, t_d)
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

generate_random_params <- function(param_ranges) {
  params <- sapply(param_ranges, function(range) {
    runif(1, min = range[1], max = range[2])
  })
  
  return(params)
}

s_analysis = function(param_ranges, diff_obs, t0){
  # Generate random parameters
  params <- generate_random_params(param_ranges)
  # Use these parameters to predict F values for b1 and b2
  results_b1 <- predict_F(K = params['K'], r = params['r'], t_m = params['t_m'], 
                          r_d = params['r_d'], t_d = params['t_d'], 
                          pd_beta = params['pd_beta'], r_pd = params['r_pd'], 
                          pr_beta = params['pr_beta'], r_pr = params['r_pr'], 
                          pd0 = params['pd0'], pr0 = params['pr0'], 
                          t0 = t0, b = 2016)
  
  results_b2 <- predict_F(K = params['K'], r = params['r'], t_m = params['t_m'], 
                          r_d = params['r_d'], t_d = params['t_d'], 
                          pd_beta = params['pd_beta'], r_pd = params['r_pd'], 
                          pr_beta = params['pr_beta'], r_pr = params['r_pr'], 
                          pd0 = params['pd0'], pr0 = params['pr0'], 
                          t0 = t0, b = 2021)
  
  # Compute the difference in F values and the correlation with actual invasion
  common_years <- intersect(results_b1$year, results_b2$year)
  diff_F <- results_b2$pre[common_years - t0 + 1] - results_b1$pre[common_years - t0 + 1]
  
  cr = cor(diff_F[396:417], common_years[396:417])
  list(parameters = params, correlation = cr, result_b1 = results_b1, results_b2 = results_b2, diff_F = diff_F)
}

t0 = 1600
b1 = 2016
b2 = 2021

years <- t0:b2

#Sensitivity analysis
#change the parameters of invasions, keeping the pd and pr
# Define the range for each parameter
param_ranges <- list(
  K = c(100, 1000),     
  r = c(0.01, 0.5),     # Replace with actual ranges
  t_m = c(1800, 2024),       # Replace with actual ranges
  r_d = c(0, 1), # Replace with actual ranges
  t_d = c(1970, 2024),
  pd_beta = c(0, 1),
  r_pd = c(0.01, 0.2),  
  pd0 = c(0.001, 0.01),
  pr_beta = c(0, 1),  
  r_pr = c(0.001, 0.2),
  pr0 = c(0.001, 0.01)
)


df_v1 = tr_df(v1.2)
df_v2 = tr_df(v2)
diff_obs = df_v2$Count[1:length(df_v1$Count)]-df_v1$Count

# correlations <- replicate(10, s_analysis(param_ranges, diff_obs, t0))

# saveRDS(correlations, "../model_result/sensitivity_analysis.rds")

increase_correlations = readRDS("../model_result/sensitivity_analysis_corrected_300_I_increase.rds")
decrease_correlations = readRDS("../model_result/sensitivity_analysis_corrected_300_I_decline.rds")
hist(unlist(increase_correlations[2,]), breaks = 50, main = "Histogram of Correlations when no decline occurs in Invasion", xlab = "Correlation")
hist(unlist(decrease_correlations[2,]), breaks = 50, main = "Histogram of Correlations when a decline occurs in Invasion", xlab = "Correlation")

d_slopes = numeric(length(decrease_correlations[5,]))
dI_slopes = numeric(length(decrease_correlations[5,]))
i_slopes = numeric(length(increase_correlations[5,]))
iI_slopes = numeric(length(increase_correlations[5,]))
start_y = 396
end_y = 417
time = years[start_y:end_y]
for (i in 1:length(d_slopes)){
  y1 = decrease_correlations[5,i]$diff_F[start_y:end_y]
  y2 = increase_correlations[5,i]$diff_F[start_y:end_y]
  yI1 = decrease_correlations[3,i]$results_b1$invasion[start_y:end_y]
  yI2 = increase_correlations[3,i]$results_b1$invasion[start_y:end_y]
  m1 = lm(y1 ~ time)
  m2 = lm(y2 ~ time)
  I1 = lm(yI1 ~ time)
  I2 = lm(yI2 ~ time)
  d_slopes[i] = coef(m1)["time"]
  i_slopes[i] = coef(m2)["time"]
  dI_slopes[i] = coef(I1)["time"]
  iI_slopes[i] = coef(I2)["time"]
}
hist(i_slopes, breaks = 30, main = "Histogram of Slopes when no decline occurs in Invasion", xlab = "Slopes")
hist(d_slopes, breaks = 30, main = "Histogram of Slopes when a decline occurs in Invasion", xlab = "Slopes")
hist(iI_slopes, breaks = 30, main = "Histogram of Slopes of Invasion when no decline occurs in Invasion", xlab = "Slopes")
hist(dI_slopes, breaks = 30, main = "Histogram of Slopes of Invasion when a decline occurs in Invasion", xlab = "Slopes")
dF_slopes = c(d_slopes, i_slopes)
I_slopes = c(dI_slopes, iI_slopes)
dF_I_lm = lm(dF_slopes ~ I_slopes)
dF_I_slope = coef(dF_I_lm)["I_slopes"]
obs_dF_lm = lm(diff_obs[start_y:end_y]~ time)
obs_df_slope = coef(obs_dF_lm)["time"]
plot(I_slopes, dF_slopes)

