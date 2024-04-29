rm(list=ls());gc()
#graphics.off()
library(ggplot2)
library(dplyr)
library(tidyr)
v1.2 = read.csv("../data/GlobalAlienSpeciesFirstRecordDatabase_v1.2.csv")
v2 = read.csv("../data/GlobalAlienSpeciesFirstRecordDatabase_v2.csv")
v3.1 = read.csv("../data/GlobalAlienSpeciesFirstRecordDatabase_v3.1_freedata.csv")

t0 = 1990
b1 = 2016
b2 = 2020

#parameters for I
# K <- 0.7e+03     # Carrying capacity
# r <- 2e-2 # Growth rate

#if I is constant
K = 989.7734
r = 0.01902461
t_m <- 1823.122    # Midpoint of logistic growth
r_d = 0.1975846
t_d = 1916.837
#if I is declining

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
  K / (1 + exp(-r * (t - t_m)))
  # if (t<t_d) {
  #   K / (1 + exp(-r * (t - t_m)))
  # }else{
  #   K / (1 + exp(-r * (t - t_m))) - exp(r_d * (t-t_d))
  # }
}

# Define the sigmoid functions for Pd and Pr
sigmoid_function <- function(t, yoi, beta, rate, p0) {
  1-exp(-rate * (t - yoi)^beta - p0)
}

# Parameters for the sigmoid functions
pd_beta <- 1.544489e+00
r_pd <- 9.615813e-07 
pd0 = 0.005175256
pr_beta <- 1.327955
r_pr <- 0.01395847
pr0 = 0.008437069
years <- t0:b2

# Function to calculate Delta F
# calculate_delta_F <- function(a, I, pd_beta, r_pd, pr_beta, r_pr, pd0, pr0, b1, b2) {
#   Fb <- 0
#   #i is the year of invasion
#   #a is the year of detection
#   for (i in t0:a) {
#     pd_product <- if (i < a) {
#       prod(1 - sapply(i:(a-1), sigmoid_function, i, pd_beta, r_pd, pd0))
#     } else {
#       1
#     }
#     pd_current <- sigmoid_function(a, i, pd_beta, r_pd, pd0)
#     pr_product <- prod(1 - sapply(a:b1, sigmoid_function, a, pr_beta, r_pr, pr0))
#     pr_diff = 1-prod(1-sapply((b1+1):b2, sigmoid_function, a, pr_beta, r_pr, pr0))
#     Fb <- Fb + I[i - t0 + 1] * pd_product * pd_current * pr_product * pr_diff
#   }
#   return(Fb)
# }

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
  years <- t0:b
  I <- sapply(years, invasion_growth, K, r, t_m, r_d, t_d)
  # Calculate Fb_values 
  Fb_values <- sapply(years, function(a) {
    Fb <- calculate_Fb(a, I, pd_beta, r_pd, pr_beta, r_pr, pd0, pr0, t0, b)
    cat("Year:", a, "Fb:", Fb, "\n")
    return(Fb)
  })
  
  result <- data.frame(year = years, pre = Fb_values, invasion = I)
  return(result)
}

df_v1 = tr_df(v1.2)
df_v2 = tr_df(v2)
df_v3 = tr_df(v3.1)
# Assuming I is a vector of invasion counts
I <- sapply(years, invasion_growth, K, r, t_m, r_d, t_d)
# pb_1600_2016 <- sapply(years, sigmoid_function, 1600, pd_beta, r_pd,  pd0)
# pr_1600_2016 <- sapply(years, sigmoid_function, 1600, pr_beta, r_pr,  pr0)
invasion = data.frame(year = years, invasion = I, pb = pb_1600_2016, pr = pr_1600_2016)
ggplot(invasion) +
  # geom_line(aes(x = year, y = invasion, color = "I")) +
  geom_line(aes(x = year, y = pb, color = "pb")) +
  geom_line(aes(x = year, y = pr, color = "pr")) +
  scale_color_manual(values = c("I" = "green4", "pb" = "blue3", "pr" = "red3")) +
  ggtitle("Invasions over Years") +
  xlab("Year") +
  ylab("Number of invasions") +
  theme_bw() +
  guides(color = guide_legend(title = "Series"))

pvo <- predict_F(K, r, t_m, r_d, t_d, pd_beta, r_pd, pr_beta, r_pr, pd0, pr0, t0, b1)
pvo2 <- predict_F(K, r, t_m, r_d, t_d, pd_beta, r_pd, pr_beta, r_pr, pd0, pr0, t0, b2)
pvo$pre_v2 = pvo2$pre[1:length(pvo$pre)]
pvo$diff = pvo$pre_v2-pvo$pre
pvo$obs_v1 = df_v1$Count
pvo$obs_v2 = df_v2$Count[1:length(pvo$pre)]
pvo$obs_v3 = df_v3$Count[1:length(pvo$pre)]
pvo$obs_diff = pvo$obs_v2-pvo$obs_v1
pvo$obs_diff2 = pvo$obs_v3-pvo$obs_v1
pvo$obs_diff3 = pvo$obs_v3-pvo$obs_v2
pvo$cum_Fb1 <- cumsum(pvo$pre)
pvo$cum_Fb2 <- cumsum(pvo$pre_v2)
pvo$cum_obs1 <- cumsum(pvo$obs_v1)
pvo$cum_obs2 <- cumsum(pvo$obs_v2)

# Calculate Delta F for each year a
# delta_F_values <- sapply(t0:b1, function(a) calculate_delta_F(a, I, pd_values, pr_values, b1, b2))
# 
# # Create a dataframe for plotting
# delta_F_df <- data.frame(Year = t0:b1, Delta_F = delta_F_values)
# delta_F_df$delf_act = pvo$diff
# Plot Delta F over a

# ggplot(delta_F_df) +
#   geom_line(aes(x = Year, y = Delta_F, color = "pre")) +
#   scale_color_manual(values = c("pre" = "blue1"))+
#   ggtitle("Change in Cumulative Recorded Invasions (Delta F) Over Time") +
#   xlab("Year") +
#   ylab("Number of invasions") +
#   theme_bw() +
#   guides(color = guide_legend(title = "Series"))

ggplot(pvo) +
  geom_line(aes(x = year, y = pre, color = "pre")) +
  geom_line(aes(x = year, y = obs_v1, color = "obs")) +
  geom_line(aes(x = year, y = obs_v3, color = "obs")) +
  geom_line(aes(x = year, y = pre_v2, color = "pre2")) +
  #geom_line(aes(x = year, y = invasion, color = "I")) +
  scale_color_manual(values = c("pre" = "blue", "obs" = "red", "pre2" = "purple", "I" = "green4")) +
  ggtitle("Prediction vs. Observation with Logistic Growth Invasions (I) and sigmoid Pd") +
  xlab("Year") +
  ylab("Number of invasions") +
  theme_bw() +
  guides(color = guide_legend(title = "Series"))
#ggsave("../model_graph/pre_invasion_logis_expo.png")
#[386:417, ]
#[200:300, ]
ggplot(pvo) +
  geom_line(aes(x = year, y = diff, color = "pre")) +
  geom_line(aes(x = year, y = obs_diff, color = "obs")) +
  # geom_line(aes(x = year, y = obs_diff2, color = "obs2")) +
  # geom_line(aes(x = year, y = obs_diff3, color = "obs3")) +
  scale_color_manual(values = c("obs3" = "blue1", "obs" = "red1", "obs2" = "red3"))+
  ggtitle("Prediction diff with Logistic Growth Invasions (I)") +
  xlab("Year") +
  ylab("Number of invasions") +
  theme_bw() +
  guides(color = guide_legend(title = "Series"))

ggplot(pvo) +
  geom_line(aes(x = year, y = cum_Fb2-cum_Fb1, color = "pre")) +
  geom_line(aes(x = year, y = cum_obs2-cum_obs1, color = "obs")) +
  scale_color_manual(values = c("pre" = "blue1", "obs" = "red1"))+
  ggtitle("Cumulative difference of Predictio and Observation with Logistic Growth Invasions (I)") +
  xlab("Year") +
  ylab("Number of invasions") +
  theme_bw() +
  guides(color = guide_legend(title = "Series"))


###test how f changed according to parameters
baseline_params <- list(K = 1.000001e+03,
                        r = 2.863829e-02,
                        t_m = 1.872851e+03,   # Midpoint of logistic growth
                        r_d = 2.000035e-01,
                        t_d = 1.995000e+03,
                        pd_beta = 1.599110e+00,
                        r_pd = 0.000000e+00,
                        pd0 = 5.407369e-03,
                        pr_beta = 1.368035e+00,
                        r_pr = 1.438073e-02,
                        pr0 = 8.712529e-03,
                        t0 = 1600, 
                        b = 2016)
results_list<- list()
K_values <- seq(800, 1000, by = 25)  # Example sequence
results_K <- list()

for (K in K_values) {
  result <- predict_F(K, baseline_params$r, baseline_params$t_m, baseline_params$r_d, baseline_params$t_d, 
                      baseline_params$pd_beta, baseline_params$r_pd, baseline_params$pr_beta, baseline_params$r_pr, 
                      baseline_params$pd0, baseline_params$pr0, baseline_params$t0, baseline_params$b)
  results_K[[paste("K", K, sep = "_")]] <- result
}

results_list$K <- results_K

plots <- lapply(results_list$K, function(df) {
  ggplot(df, aes(x = year, y = pre)) + 
    geom_line() +
    ggtitle(paste("Impact of varying K on Fb over years"))
})

# Display plots (You might need to adjust this part based on your actual data structure)
print(plots[[1]])  # Example to print the first plot


r_values <- seq(0.01, 0.1, by = 0.01)  # Example sequence
results_r <- list()

for (r in r_values) {
  result <- predict_F(baseline_params$K, r, baseline_params$t_m, baseline_params$r_d, baseline_params$t_d, 
                      baseline_params$pd_beta, baseline_params$r_pd, baseline_params$pr_beta, baseline_params$r_pr, 
                      baseline_params$pd0, baseline_params$pr0, baseline_params$t0, baseline_params$b)
  results_r[[paste("r", r, sep = "_")]] <- result
}

results_list$r <- results_r

pd_beta_values <- seq(0.5, 2, by = 0.02)  # Example sequence
results_pd_beta <- list()

for (pd_beta in pd_beta_values) {
  result <- predict_F(baseline_params$K, baseline_params$r, baseline_params$t_m, baseline_params$r_d, baseline_params$t_d, 
                      pd_beta, baseline_params$r_pd, baseline_params$pr_beta, baseline_params$r_pr, 
                      baseline_params$pd0, baseline_params$pr0, baseline_params$t0, baseline_params$b)
  results_r[[paste("pd_beta", pd_beta, sep = "_")]] <- result
}

results_list$pd_beta <- results_pd_beta

r_pd_values <- seq(0.01, 0.5, by = 0.01)  # Example sequence
results_r_pd <- list()

for (r_pd in r_pd_values) {
  result <- predict_F(baseline_params$K, baseline_params$r, baseline_params$t_m, baseline_params$r_d, baseline_params$t_d, 
                      baseline_params$pd_beta, r_pd, baseline_params$pr_beta, baseline_params$r_pr, 
                      baseline_params$pd0, baseline_params$pr0, baseline_params$t0, baseline_params$b)
  results_r[[paste("r_pd", r_pd, sep = "_")]] <- result
}

results_list$r_pd <- results_pd_beta