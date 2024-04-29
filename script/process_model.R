rm(list=ls());gc()
graphics.off()
library(dplyr)
library(ggplot2)
library(tidyr)
# Define the starting year t0 and end year b
# These are determined by the data we have
save_par = "../model_result/pars_v1.rds"
save_graph = "../model_graph/pre_obs_invasion_v1.png"


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
# Initial parameter values for optimisation
init_params <- c(K = 1800, r = 0.01, t_m = 1950, pd = 0.1, pr = 0.2)
t0 <- 1600
species_count_df = tr_df(v1.2)

b = max(species_count_df$FirstRecord)


# Call optim to minimize the loss function
opt_results <- optim(init_params, loss, method = "L-BFGS-B", lower = c(1, 0.0001, t0, 0.000001, 0.000001), upper = c(1e6, 1, 2021, 1, 1))
#opt_results <- optim(init_params, loss, method = "L-BFGS-B", lower = c(1, optimized_params[2], optimized_params[3], optimized_params[4], optimized_params[5]), upper = c(1e6, optimized_params[2], optimized_params[3], optimized_params[4], optimized_params[5]))
# Optimized parameter values
optimized_params <- opt_results$par

# Print the optimized parameters
print(optimized_params)
saveRDS(optimized_params, save_par)

pvo = predict(optimized_params, b)
pvo_v2 = predict(optimized_params, max(v2$FirstRecord))
pvo_v3 = predict(optimized_params, max(v3.1$FirstRecord))
pvo$pre_v2 = pvo_v2$pre[1:length(pvo$pre)]
pvo$pre_v3 = pvo_v3$pre[1:length(pvo$pre)]
pvo$obs = species_count_df$Count
ggplot(pvo) +
  geom_line(aes(x = year, y = pre, color = "pre")) +
  #geom_line(aes(x = year, y = obs, color = "obs")) +
  geom_line(aes(x = year, y = pre_V2, color = "pre2")) +
  geom_line(aes(x = year, y = pre_v3, color = "pre3")) +
  #geom_line(aes(x = year, y = invasion, color = "I")) +
  scale_color_manual(values = c("pre" = "blue", "obs" = "red", "pre3" = "purple", "pre2" = "green")) +
  ggtitle("Prediction vs. Observation with Logistic Growth Invasions (I)") +
  xlab("Year") +
  ylab("Number of invasions") +
  theme_bw() +
  guides(color = guide_legend(title = "Series"))
ggsave(save_graph)

pvo <- pvo %>%
  mutate(cum_pre = cumsum(pre),
         cum_pre2 = cumsum(pre_v2),
         cum_pre3 = cumsum(pre_v3))

# Reshape the dataframe to long format
pvo_long <- pvo[400:417,] %>%
  select(year, cum_pre, cum_pre2, cum_pre3) %>%
  pivot_longer(cols = -year, names_to = "predictor", values_to = "cumulative_invasions")

# Plot the data
ggplot(pvo_long, aes(x = year, y = cumulative_invasions, color = predictor)) +
  geom_line() +
  ggtitle("Cumulative Number of Invasions Across Years") +
  xlab("Year") +
  ylab("Cumulative Number of Invasions") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

### Round2
#now let's train the model on the second version
# init_params <- optimized_params
# t0 <- 1600
# species_count_df = tr_df(v2)
# b = max(species_count_df$FirstRecord)
# opt_results_rnd2 <- optim(init_params, loss, method = "L-BFGS-B", lower = c(1, 0.0001, t0, 0.000001, 0.000001), upper = c(1e6, 1, 2021, 1, 1))
# optimized_params_rnd2 = opt_results_rnd2$par
# pvo_v2 = predict(optimized_params_rnd2, species_count_df$Count)
# pvo$prev2 = pvo_v2$pre[1:length(pvo$pre)]
# pvo$obsv2 = pvo_v2$obs[1:length(pvo$obs)]
# 
# 
# #plot
# ggplot(pvo) +
#   geom_line(aes(x = year, y = pre, color = "pre")) +
#   geom_line(aes(x = year, y = obs, color = "obs")) +
#   geom_line(aes(x = year, y = prev2, color = "prev2")) +
#   geom_line(aes(x = year, y = obsv2, color = "obsv2")) +
#   #geom_line(aes(x = year, y = invasion, color = "I")) +
#   scale_color_manual(values = c("pre" = "blue", "obs" = "red", "I" = "purple", "prev2" = "green", "obsv2" = "orange")) +
#   ggtitle("Prediction vs. Observation with Logistic Growth Invasions (I)") +
#   xlab("Year") +
#   ylab("Number of invasions") +
#   theme_bw() +
#   guides(color = guide_legend(title = "Series"))
# 
# #Use version 3 as a testing df
# species_count_df = tr_df(v3.1)
# species_count_df = species_count_df[1:(b-t0+1),]
# pvo_v3.1 = predict(optimized_params_rnd2, species_count_df$Count)
# pvo$prev3.1 = pvo_v3.1$pre[1:length(pvo$pre)]
# pvo$obsv3.1 = pvo_v3.1$obs[1:length(pvo$obs)]
# ggplot(pvo) +
#   geom_line(aes(x = year, y = pre, color = "pre")) +
#   #geom_line(aes(x = year, y = obs, color = "obs")) +
#   geom_line(aes(x = year, y = prev2, color = "prev2")) +
#   #geom_line(aes(x = year, y = obsv2, color = "obsv2")) +
#   #geom_line(aes(x = year, y = prev3.1, color = "prev3")) +
#   geom_line(aes(x = year, y = obsv3.1, color = "obsv3")) +
#   #geom_line(aes(x = year, y = invasion, color = "I")) +
#   scale_color_manual(values = c("pre" = "blue", "obs" = "red", "prev2" = "green", "obsv2" = "orange", "prev3" = "purple", "obsv3" = "yellow")) +
#   ggtitle("Prediction vs. Observation with Logistic Growth Invasions (I)") +
#   xlab("Year") +
#   ylab("Number of invasions") +
#   theme_bw() +
#   guides(color = guide_legend(title = "Series"))
