rm(list=ls());gc()
#graphics.off()
library(dplyr)
library(ggplot2)
library(tidyr)
# Define the starting year t0
t0 <- 1600
# These are determined by the data we have
save_par_v1 = "../model_result/pars_v1-2_logis_sig_sig.rds"
save_par_v2 = "../model_result/pars_v2_logis_sig_sig.rds"
save_par_v3 = "../model_result/pars_v3-1_logis_sig_sig.rds"
fit_v1 = FALSE
fit_v2 = FALSE
fit_v3 = FALSE
save_model_v1 = "../model/model_v1-2_logis_sig_sig.rds"
save_model_v2 = "../model/model_v2_logis_sig_sig.rds"
save_model_v3 = "../model/model_v3-1_logis_sig_sig.rds"
# v1.1 = read.csv("../data/GlobalAlienSpeciesFirstRecordDatabase_v1.1.csv")
v1.2 = read.csv("../data/GlobalAlienSpeciesFirstRecordDatabase_v1.2.csv")
v2 = read.csv("../data/GlobalAlienSpeciesFirstRecordDatabase_v2.csv")

v3.1 = read.csv("../data/GlobalAlienSpeciesFirstRecordDatabase_v3.1_freedata.csv")

#create the training data
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

#set the final year

# Constants for the logistic function
# These are the parameters to fit

### Here we define the how actual invasions vary through years
#define a logistic growth function for Ii
invasion_growth <- function(t, K, r, t_m, r_d, t_d) {
  if (t<t_d) {
    K / (1 + exp(-r * (t - t_m)))
  }else{
    K / (1 + exp(-r * (t - t_m))) - exp(r_d * (t-t_d))
  }
}

#Simulate that the probability of detection increases exponentially due to abundance
sigmoid_function <- function(t, yoi, beta, rate, p0) {
  1-exp(-rate * (t - yoi)^beta - p0)
}

#formula

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


# Generate I values using the logistic growth function

# Loss function definition
loss <- function(params) {
  K <- params[1]
  r <- params[2]
  t_m <- params[3]
  r_d <- params[4]
  t_d <- params[5]
  pd_beta <- params[6]
  r_pd <- params[7]
  pd0 <- params[8]
  pr_beta <- params[9]
  r_pr <- params[10]
  pr0 <- params[11]
  
  years <- t0:b
  I <- sapply(years, invasion_growth, K, r, t_m, r_d, t_d)
  
  Fb_values <- sapply(years, function(a) calculate_Fb(a, I, pd_beta, r_pd, pr_beta, r_pr, pd0, pr0, t0, b))
  
  obs <- species_count_df$Count
  ls = sum(abs(Fb_values - obs))
  print(params)
  print(ls)
  if (is.na(ls)) {
    ls = 1e10
  }
  return(ls)
}

predict = function(params, b){
  K <- params[1]
  r <- params[2]
  t_m <- params[3]
  r_d <- params[4]
  t_d <- params[5]
  pd_beta <- params[6]
  r_pd <- params[7]
  pd0 <- params[8]
  pr_beta <- params[9]
  r_pr <- params[10]
  pr0 <- params[11]
  
  years <- t0:b
  I <- sapply(years, invasion_growth, K, r, t_m, r_d, t_d)
  # Calculate Fb_values
  Fb_values <- sapply(years, function(a) {
    Fb <- calculate_Fb(a, I, pd_beta, r_pd, pr_beta, r_pr, pd0, pr0, t0, b)
    cat("Year:", a, "Fb:", Fb, "\n")
    return(Fb)
  })
  result = data.frame(year = years, pre = Fb_values, invasion = I)
  return(result)
}
### DF1 ###
# Initial parameter values for optimisation
#(K, r, t_m, r_d, t_d, pd_beta, r_pd, pd0, pr_beta, r_pr, pr0)
lower_bounds <- c(1, 0.0001, t0, 0.0001, t0,  -Inf, 0, 0, -Inf, 0, 0)
upper_bounds <- c(1e6, 1, 2050, 1, 2050,  Inf, Inf, Inf, Inf, Inf, Inf)

if(fit_v1){
  cat("Start Fitting V1 \n")
  init_params <- c(K = 1e+03, r = 2e-02, t_m = 1.872854e+03, r_d=0.2, t_d = 1995, pd_beta = 1, pd_rate = 1e-03, pd0 = 0.005, pr_beta = 1, pr_rate = 1e-03, pr0 = 0.005)
  species_count_df = tr_df(v1.2)
  b = max(species_count_df$FirstRecord)
  # Call optim to minimize the loss function
  opt_results <- optim(init_params, loss, method = "L-BFGS-B", lower = lower_bounds, upper = upper_bounds)
  #opt_results <- optim(init_params, loss, method = "L-BFGS-B", lower = c(1, optimized_params[2], optimized_params[3], optimized_params[4], optimized_params[5]), upper = c(1e6, optimized_params[2], optimized_params[3], optimized_params[4], optimized_params[5]))
  # Optimized parameter values
  optimized_params <- opt_results$par
  # Print the optimized parameters
  print(optimized_params)
  saveRDS(optimized_params, save_par_v1)
  saveRDS(opt_results, save_model_v1)
}


### DF2 ###
if(fit_v2){
  cat("Start Fitting V2 \n")
  init_params <- readRDS(save_par_v1)
  species_count_df = tr_df(v2)
  b = max(species_count_df$FirstRecord)
  # Call optim to minimize the loss function
  opt_results <- optim(init_params, loss, method = "L-BFGS-B", lower = lower_bounds, upper = upper_bounds)
  #opt_results <- optim(init_params, loss, method = "L-BFGS-B", lower = c(1, optimized_params[2], optimized_params[3], optimized_params[4], optimized_params[5]), upper = c(1e6, optimized_params[2], optimized_params[3], optimized_params[4], optimized_params[5]))
  # Optimized parameter values
  optimized_params <- opt_results$par
  # Print the optimized parameters
  print(optimized_params)
  saveRDS(optimized_params, save_par_v2)
  saveRDS(opt_results, save_model_v2)
}


### DF3 ###
if(fit_v3){
  cat("Start Fitting V3 \n")
  init_params <- readRDS(save_par_v2)
  species_count_df = tr_df(v3.1)
  b = max(species_count_df$FirstRecord)
  # Call optim to minimize the loss function
  opt_results <- optim(init_params, loss, method = "L-BFGS-B", lower = lower_bounds, upper = upper_bounds)
  #opt_results <- optim(init_params, loss, method = "L-BFGS-B", lower = c(1, optimized_params[2], optimized_params[3], optimized_params[4], optimized_params[5]), upper = c(1e6, optimized_params[2], optimized_params[3], optimized_params[4], optimized_params[5]))
  # Optimized parameter values
  optimized_params <- opt_results$par
  # Print the optimized parameters
  print(optimized_params)
  saveRDS(optimized_params, save_par_v3)
  saveRDS(opt_results, save_model_v3)
}
df_v1 = tr_df(v1.2)
df_v2 = tr_df(v2)
df_v3 = tr_df(v3.1)
pars_v1 = readRDS(save_par_v1)
pars_v2 = readRDS(save_par_v2)
pars_v3 = readRDS(save_par_v3)
pvo = predict(pars_v1, max(v1.2$FirstRecord))
pvo_v2 = predict(pars_v2, max(v2$FirstRecord))
pvo_v3 = predict(pars_v3, max(v3.1$FirstRecord))
pvo$pre_v2 = pvo_v2$pre[1:length(pvo$pre)]
pvo$pre_v3 = pvo_v3$pre[1:length(pvo$pre)]
pvo$pd_v2 = pvo_v2$pd[1:length(pvo$pre)]
pvo$pd_v3 = pvo_v3$pd[1:length(pvo$pre)]
pvo$cum_Fb1 <- cumsum(pvo$pre)
pvo$cum_Fb2 <- cumsum(pvo$pre_v2)
pvo$cum_Fb3 <- cumsum(pvo$pre_v3)
pvo$obs_v1 = df_v1$Count
pvo$obs_v2 = df_v2$Count[1:length(pvo$pre)]
pvo$obs_v3 = df_v3$Count[1:length(pvo$pre)]
# pvo$obs = species_count_df$Count
ggplot(pvo) +
  geom_line(aes(x = year, y = pre, color = "pre")) +
  geom_line(aes(x = year, y = obs_v1, color = "obs")) +
  geom_line(aes(x = year, y = pre_v2, color = "pre2")) +
  geom_line(aes(x = year, y = obs_v2, color = "obs2")) +
  geom_line(aes(x = year, y = pre_v3, color = "pre3")) +
  geom_line(aes(x = year, y = obs_v3, color = "obs3")) +
  # geom_line(aes(x = year, y = invasion, color = "I")) +
  scale_color_manual(values = c("pre" = "blue", "pre3" = "purple", "pre2" = "green", "obs" = "red3", "obs2" = "orange3", "obs3" = "yellow3")) +
  ggtitle("Prediction vs. Observation with Logistic Growth Invasions (I) and Sigmoid Pd") +
  xlab("Year") +
  ylab("Number of invasions") +
  theme_bw() +
  guides(color = guide_legend(title = "Series"))
ggsave("../model_graph/pre_invasion_logis_sig.png")


ggplot(pvo) +
  # geom_line(aes(x = Year, y = Fb1, color = "F(1)")) +
  # geom_line(aes(x = Year, y = Fb2, color = "F(2)")) +
  # geom_line(aes(x = Year, y = Fb3, color = "F(3)")) +
  # geom_line(aes(x = Year, y = I, color = "I")) +
  geom_line(aes(x = year, y = cum_Fb1, color = "Cumulative F(1)")) +
  geom_line(aes(x = year, y = cum_Fb2, color = "Cumulative F(2)")) +
  geom_line(aes(x = year, y = cum_Fb3, color = "Cumulative F(3)")) +
  scale_color_manual(values = c("F(1)" = "blue", "F(2)" = "green", "F(3)" = "purple", 
                                "I" = "red", "Cumulative F(1)" = "blue4", 
                                "Cumulative F(2)" = "green4", "Cumulative F(3)" = "purple4")) +
  ggtitle("Change in Cumulative Recorded Invasions (F(b)) and Logistic Growth Invasions (I) and Sigmoid Pd") +
  xlab("Year 'a'") +
  ylab("Value") +
  theme_bw() +
  guides(color = guide_legend(title = "Series"))
#ggsave("../model_graph/cumu_pre_invasion_logis_sig.png")


years = t0:1750
#the curves of pd and pr
pd_v1 <- sigmoid_function(years, 1600, pars_v1['pd_beta'], pars_v1['pd_rate'], pars_v1['pd0'])
pr_v1 <- sigmoid_function(years, 1600, pars_v1['pr_beta'], pars_v1['pr_rate'], pars_v1['pr0'])
pd_v2 <- sigmoid_function(years, 1600, pars_v2['pd_beta'], pars_v2['pd_rate'], pars_v2['pd0'])
pr_v2 <- sigmoid_function(years, 1600, pars_v2['pr_beta'], pars_v2['pr_rate'], pars_v2['pr0'])
pd_v3 <- sigmoid_function(years, 1600, pars_v3['pd_beta'], pars_v3['pd_rate'], pars_v3['pd0'])
pr_v3 <- sigmoid_function(years, 1600, pars_v3['pr_beta'], pars_v3['pr_rate'], pars_v3['pr0'])
# Repeat the calculation for versions 2 and 3 with their respective parameters

ggplot() +
  # geom_line(aes(x = Year, y = Fb1, color = "F(1)")) +
  # geom_line(aes(x = Year, y = Fb2, color = "F(2)")) +
  # geom_line(aes(x = Year, y = Fb3, color = "F(3)")) +
  # geom_line(aes(x = Year, y = I, color = "I")) +
  geom_line(aes(x = years, y = pd_v1, color = "p1")) +
  geom_line(aes(x = years, y = pd_v2, color = "p2")) +
  geom_line(aes(x = years, y = pd_v3, color = "p3")) +
  scale_color_manual(values = c("p1" = "blue", "p2" = "red2", "p3" = "purple")) +
  ggtitle("Probability of Detection (pd) Over Time") +
  xlab("Year") +
  ylab("Probability of Detection (pd)") +
  theme_bw() +
  guides(color = guide_legend(title = "Series"))

ggplot() +
  # geom_line(aes(x = Year, y = Fb1, color = "F(1)")) +
  # geom_line(aes(x = Year, y = Fb2, color = "F(2)")) +
  # geom_line(aes(x = Year, y = Fb3, color = "F(3)")) +
  # geom_line(aes(x = Year, y = I, color = "I")) +
  geom_line(aes(x = years, y = pr_v1, color = "p1")) +
  geom_line(aes(x = years, y = pr_v2, color = "p2")) +
  geom_line(aes(x = years, y = pr_v3, color = "p3")) +
  scale_color_manual(values = c("p1" = "blue", "p2" = "red2", "p3" = "purple")) +
  ggtitle("Probability of Reporting (pr) Over Time") +
  xlab("Year") +
  ylab("Probability of Reporting (pr)") +
  theme_bw() +
  guides(color = guide_legend(title = "Series"))


