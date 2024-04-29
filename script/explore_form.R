rm(list=ls());gc()
graphics.off()
library(dplyr)
library(ggplot2)

# Define the starting and ending years
v1.2 = read.csv("../data/GlobalAlienSpeciesFirstRecordDatabase_v1.2.csv")
v2 = read.csv("../data/GlobalAlienSpeciesFirstRecordDatabase_v2.csv")

v3.1 = read.csv("../data/GlobalAlienSpeciesFirstRecordDatabase_v3.1_freedata.csv")

t0 <- 1600
b <- 2016
b2 <- 2020
b3 <- 2023
# Define parameters for the logistic invasion model
K <- 1.904281e+03     # Carrying capacity
r <- 2.585832e-02      # Growth rate
t_m <- 1.872854e+03    # Midpoint of logistic growth

# Define parameters for the exponential probability of detection
#pd0_v1 <- 0.002308069
pd0 <- 0.468309e-03# Initial probability of detection at t0
r_pd <- 4.343891e-03   # Growth rate of probability of detection

pr0 <- 1.667608e-01   # Initial probability of reporting at t0
r_pr <- 0 
# Define a constant probability of reporting
#pr <- 1.667608e-01

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

logistic_growth <- function(t, K, r, t_m) {
  K / (1 + exp(-r * (t - t_m)))
}

exponential_pd <- function(t, pd0, r_pd) {
  pd0 * exp(r_pd * (t - t0))
}

exponential_pr <- function(t, pr0, r_pr) {
  pr0 * exp(r_pr * (t - t0))
}


calculate_Fb <- function(a, I, pd0, r_pd, pr0, r_pr, t0, b) {
  Fb <- 0
  for (i in t0:a) {
    pd_product <- if (i < a) {
      prod(1 - sapply(i:(a-1), exponential_pd, pd0, r_pd))
    } else {
      1
    }
    pd_current <- exponential_pd(a, pd0, r_pd)
    pr_product <- 1-prod(1 - sapply(a:b, exponential_pr, pr0, r_pr))
    
    Fb <- Fb + I[i - t0 + 1] * pd_product * pd_current * pr_product
  }
  return(Fb)
}

predict_F <- function(K, r, t_m, pd0, r_pd, pr0, r_pr, b){
  years <- t0:b
  I <- sapply(years, logistic_growth, K, r, t_m)
  pd_values <- sapply(years - t0, exponential_pd, pd0, r_pd)
  pr_values <- sapply(years - t0, exponential_pr, pr0, r_pr)
  # Calculate Fb_values and log calculations
  Fb_values <- sapply(years, function(a) {
    Fb <- calculate_Fb(a, I, pd0, r_pd, pr0, r_pr, t0, b)
    cat("Year:", a, "Fb:", Fb, "\n")
    return(Fb)
  })
  
  # Log the range of pd_values
  cat("pd_values range:", range(pd_values), "\n")
  
  # Calculate the effect of the reporting term for different 'a' and 'b' values
  cat("Reporting effect range:", range(pr_values), "\n")
  
  result <- data.frame(year = years, pre = Fb_values, invasion = I, pd = pd_values, pr = pr_values)
  return(result)
}

# Call the predict_F function for different 'b' values and store the results
pvo <- predict_F(K, r, t_m, pd0, r_pd, pr0, r_pr, b)
pvo2 <- predict_F(K, r, t_m, pd0, r_pd, pr0, r_pr, b2)
pvo3 <- predict_F(K, r, t_m, pd0, r_pd, pr0, r_pr, b3)

df_v1 = tr_df(v1.2)
df_v2 = tr_df(v2)
df_v3 = tr_df(v3.1)

pvo$pre_v2 = pvo2$pre[1:length(pvo$pre)]
pvo$pre_v3 = pvo3$pre[1:length(pvo$pre)]
pvo$cum_Fb1 <- cumsum(pvo$pre)
pvo$cum_Fb2 <- cumsum(pvo$pre_v2)
pvo$cum_Fb3 <- cumsum(pvo$pre_v3)
pvo$obs_v1 = df_v1$Count
ggplot(pvo) +
  geom_line(aes(x = year, y = pre, color = "pre")) +
  geom_line(aes(x = year, y = pre_v2, color = "pre2")) +
  geom_line(aes(x = year, y = pre_v3, color = "pre3")) +
  #geom_line(aes(x = year, y = invasion, color = "I")) +
  geom_line(aes(x = year, y = obs_v1, color = "pd")) +
  #geom_line(aes(x = year, y = pd, color = "pd"))
  scale_color_manual(values = c("pre" = "blue1", "pre2" = "blue3", "pre3" = "blue4", "obs" = "red", "pd" = "purple")) +
  ggtitle("Prediction vs. Observation with Logistic Growth Invasions (I)") +
  xlab("Year") +
  ylab("Number of invasions") +
  theme_bw() +
  guides(color = guide_legend(title = "Series"))

ggplot(pvo) +
  geom_line(aes(x = year, y = cum_Fb1, color = "Cumulative F(1)")) +
  geom_line(aes(x = year, y = cum_Fb2, color = "Cumulative F(2)")) +
  geom_line(aes(x = year, y = cum_Fb3, color = "Cumulative F(3)")) +
  scale_color_manual(values = c("F(1)" = "blue", "F(2)" = "green", "F(3)" = "purple", 
                                "I" = "red", "Cumulative F(1)" = "blue4", 
                                "Cumulative F(2)" = "green4", "Cumulative F(3)" = "purple4")) +
  ggtitle("Change in Cumulative Recorded Invasions (F(b)) and Logistic Growth Invasions (I) and Exponential Pd") +
  xlab("Year 'a'") +
  ylab("Value") +
  theme_bw() +
  guides(color = guide_legend(title = "Series"))

ggplot(pvo) +
  geom_line(aes(x = year, y = pre_v2-pre, color = "pre")) +
  # geom_line(aes(x = year, y = pre_v2, color = "pre2")) +
  # geom_line(aes(x = year, y = pre_v3, color = "pre3")) +
  #geom_line(aes(x = year, y = invasion, color = "I")) +
  #geom_line(aes(x = year, y = obs_v1, color = "pd")) +
  #geom_line(aes(x = year, y = pd, color = "pd"))
  scale_color_manual(values = c("pre" = "blue1", "pre2" = "blue3", "pre3" = "blue4", "obs" = "red", "pd" = "purple")) +
  ggtitle("Prediction diff with Logistic Growth Invasions (I)") +
  xlab("Year") +
  ylab("Number of invasions") +
  theme_bw() +
  guides(color = guide_legend(title = "Series"))
