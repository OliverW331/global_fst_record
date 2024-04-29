rm(list=ls());gc()
graphics.off()
library(dplyr)
library(ggplot2)
library(tidyr)
# Define the starting year t0
t0 <- 1600
# These are determined by the data we have
save_par_v1 = "../model_result/pars_v1-2_logis.rds"
save_par_v2 = "../model_result/pars_v2_logis.rds"
save_par_v3 = "../model_result/pars_v3-1_logis.rds"
fit_v1 = FALSE
fit_v2 = FALSE
fit_v3 = FALSE
save_model_v1 = "../model/model_v1-2_logis.rds"
save_model_v2 = "../model/model_v2_logis.rds"
save_model_v3 = "../model/model_v3-1_logis.rds"
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
logistic_growth <- function(t, K, r, t_m) {
  K / (1 + exp(-r * (t - t_m)))
}


#formula
calculate_Fb <- function(a, I, pd, pr, t0, b) {
  #Fb <- Fb
  Fb = 0
  for (i in t0:a) {
    index <- i - t0 + 1  # Adjust index for I vector
    Fb <- Fb + I[index] * (1 - pd)^(a - i) * pd * (1 - (1 - pr)^(b - a))
  }
  return(Fb)
}

# Generate I values using the logistic growth function

# Loss function definition
loss <- function(params) {
  # Unpack parameters
  K <- params[1]
  r <- params[2]
  t_m <- params[3]
  pd <- params[4]
  pr <- params[5]
  #Fb <- params[6]
  # Calculate I using logistic growth function
  years <- t0:b
  I <- sapply(years, logistic_growth, K, r, t_m)
  
  # Calculate Fb_values
  Fb_values <- sapply(years, function(a) calculate_Fb(a, I, pd, pr, t0, b))
  
  # Calculate loss as the sum of absolute differences
  obs <- species_count_df$Count
  ls = sum(abs(Fb_values - obs))
  print(params)
  print(ls)
  return (ls)
}

predict = function(params, b){
  K <- params[1]
  r <- params[2]
  t_m <- params[3]
  pd <- params[4]
  pr <- params[5]
  #Fb <- params[6]
  years <- t0:b
  I <- sapply(years, logistic_growth, K, r, t_m)
  
  # Calculate Fb_values
  Fb_values <- sapply(years, function(a) calculate_Fb(a, I, pd, pr, t0, b))
  result = data.frame(year = years, pre = Fb_values, invasion = I)
  return(result)
}
### DF1 ###
# Initial parameter values for optimisation

if(fit_v1){
  init_params <- c(K = 1800, r = 0.01, t_m = 1950, pd = 0.1, pr = 0.2)
  species_count_df = tr_df(v1.2)
  b = max(species_count_df$FirstRecord)
  # Call optim to minimize the loss function
  opt_results <- optim(init_params, loss, method = "L-BFGS-B", lower = c(1, 0.0001, t0, 0.000001, 0.000001), upper = c(1e6, 1, 2021, 1, 1))
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
  init_params <- readRDS(save_par_v1)
  species_count_df = tr_df(v2)
  b = max(species_count_df$FirstRecord)
  # Call optim to minimize the loss function
  opt_results <- optim(init_params, loss, method = "L-BFGS-B", lower = c(1, 0.0001, t0, 0.000001, 0.000001), upper = c(1e6, 1, 2021, 1, 1))
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
  init_params <- readRDS(save_par_v2)
  species_count_df = tr_df(v3.1)
  b = max(species_count_df$FirstRecord)
  # Call optim to minimize the loss function
  opt_results <- optim(init_params, loss, method = "L-BFGS-B", lower = c(1, 0.0001, t0, 0.000001, 0.000001), upper = c(1e6, 1, 2021, 1, 1))
  #opt_results <- optim(init_params, loss, method = "L-BFGS-B", lower = c(1, optimized_params[2], optimized_params[3], optimized_params[4], optimized_params[5]), upper = c(1e6, optimized_params[2], optimized_params[3], optimized_params[4], optimized_params[5]))
  # Optimized parameter values
  optimized_params <- opt_results$par
  # Print the optimized parameters
  print(optimized_params)
  saveRDS(optimized_params, save_par_v3)
  saveRDS(opt_results, save_model_v3)
}

pars_v1 = readRDS(save_par_v1)
pars_v2 = readRDS(save_par_v2)
pars_v3 = readRDS(save_par_v3)
pvo = predict(pars_v1, max(v1.2$FirstRecord))
pvo_v2 = predict(pars_v2, max(v2$FirstRecord))
pvo_v3 = predict(pars_v3, max(v3.1$FirstRecord))
pvo$pre_v2 = pvo_v2$pre[1:length(pvo$pre)]
pvo$pre_v3 = pvo_v3$pre[1:length(pvo$pre)]
df_v1 = tr_df(v1.2)
df_v2 = tr_df(v2)
df_v3 = tr_df(v3.1)
pvo$cum_Fb1 <- cumsum(pvo$pre)
pvo$cum_Fb2 <- cumsum(pvo$pre_v2)
pvo$cum_Fb3 <- cumsum(pvo$pre_v3)
pvo$obs_v1 = df_v1$Count
pvo$obs_v2 = df_v2$Count[1:length(pvo$pre)]
pvo$obs_v3 = df_v3$Count[1:length(pvo$pre)]

# pvo$obs = species_count_df$Count
ggplot(pvo) +
  geom_line(aes(x = year, y = pre, color = "pre")) +
  geom_point(aes(x = year, y = obs_v1, color = "obs")) +
  geom_line(aes(x = year, y = pre_v2, color = "pre2")) +
  geom_point(aes(x = year, y = obs_v2, color = "obs2")) +
  geom_line(aes(x = year, y = pre_v3, color = "pre3")) +
  geom_point(aes(x = year, y = obs_v3, color = "obs3")) +
  #geom_line(aes(x = year, y = invasion, color = "I")) +
  scale_color_manual(values = c("pre" = "blue", "obs" = "blue4", "pre3" = "purple","obs3" = "purple4", "pre2" = "green", "obs2" = "green4")) +
  ggtitle("Prediction vs. Observation with Logistic Growth Invasions (I)") +
  xlab("Year") +
  ylab("Number of invasions") +
  theme_bw() +
  guides(color = guide_legend(title = "Series"))
ggsave("../model_graph/pre_invasion_logis.png")


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
  ggtitle("Change in Cumulative Recorded Invasions (F(b)) and Logistic Growth Invasions (I)") +
  xlab("Year 'a'") +
  ylab("Value") +
  theme_bw() +
  guides(color = guide_legend(title = "Series"))
ggsave("../model_graph/cumu_pre_invasion_logis.png")

### Find potential pd change over time
pd_df = data.frame(pd = c(pars_v1[4], pars_v2[4], pars_v3[4]), b = c(2016, 2020, 2021))

