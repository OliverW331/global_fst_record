rm(list=ls());gc()
#graphics.off()
library(ggplot2)
library(dplyr)
library(tidyr)

invasion_growth <- function(t, K, r, t_m, r_d, t_d) {
  if (t<t_d) {
    K / (1 + exp(-r * (t - t_m)))
  }else{
    K / (1 + exp(-r * (t - t_m))) - exp(r_d * (t-t_d))
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

baseline_params <- list(K = 1.000001e+03,
                        r = 2.863829e-02,
                        t_m = 1.872851e+03,   # Midpoint of logistic growth
                        r_d = 2.000035e-01,
                        t_d = 1.995000e+03,
                        pd_beta = 1.599110e+00,
                        r_pd = 0.001e+00,
                        pd0 = 5.407369e-03,
                        pr_beta = 1.368035e+00,
                        r_pr = 1.438073e-02,
                        pr0 = 8.712529e-03,
                        t0 = 1995, 
                        b = 2016)

results_list<- list()

### K ###
K_values <- seq(800, 1000, by = 25)  # Example sequence
results_K <- list()

for (K in K_values) {
  result <- predict_F(K, baseline_params$r, baseline_params$t_m, baseline_params$r_d, baseline_params$t_d, 
                      baseline_params$pd_beta, baseline_params$r_pd, baseline_params$pr_beta, baseline_params$r_pr, 
                      baseline_params$pd0, baseline_params$pr0, baseline_params$t0, baseline_params$b)
  results_K[[paste("K", K, sep = "_")]] <- result
}

results_list$K <- results_K
gc()
### r ###
r_values <- seq(0.01, 0.1, by = 0.01)  # Example sequence
results_r <- list()

for (r in r_values) {
  result <- predict_F(baseline_params$K, r, baseline_params$t_m, baseline_params$r_d, baseline_params$t_d, 
                      baseline_params$pd_beta, baseline_params$r_pd, baseline_params$pr_beta, baseline_params$r_pr, 
                      baseline_params$pd0, baseline_params$pr0, baseline_params$t0, baseline_params$b)
  results_r[[paste("r", r, sep = "_")]] <- result
}

results_list$r <- results_r
gc()
### pd_beta ###
pd_beta_values <- seq(0.5, 2, by = 0.1)  # Example sequence
results_pd_beta <- list()

for (pd_beta in pd_beta_values) {
  result <- predict_F(baseline_params$K, baseline_params$r, baseline_params$t_m, baseline_params$r_d, baseline_params$t_d, 
                      pd_beta, baseline_params$r_pd, baseline_params$pr_beta, baseline_params$r_pr, 
                      baseline_params$pd0, baseline_params$pr0, baseline_params$t0, baseline_params$b)
  results_pd_beta[[paste("pd_beta", pd_beta, sep = "_")]] <- result
}

results_list$pd_beta <- results_pd_beta
gc()
### r_pd ###
r_pd_values <- seq(0.01, 0.1, by = 0.01)  # Example sequence
results_r_pd <- list()

for (r_pd in r_pd_values) {
  result <- predict_F(baseline_params$K, baseline_params$r, baseline_params$t_m, baseline_params$r_d, baseline_params$t_d, 
                      baseline_params$pd_beta, r_pd, baseline_params$pr_beta, baseline_params$r_pr, 
                      baseline_params$pd0, baseline_params$pr0, baseline_params$t0, baseline_params$b)
  results_r_pd[[paste("r_pd", r_pd, sep = "_")]] <- result
}

results_list$r_pd <- results_r_pd
gc()
### pd0 ###
pd0_values <- seq(0.001, 0.01, by = 0.001)  # Example sequence
results_pd0 <- list()

for (pd0 in pd0_values) {
  result <- predict_F(baseline_params$K, baseline_params$r, baseline_params$t_m, baseline_params$r_d, baseline_params$t_d, 
                      baseline_params$pd_beta, baseline_params$r_pd, baseline_params$pr_beta, baseline_params$r_pr, 
                      pd0, baseline_params$pr0, baseline_params$t0, baseline_params$b)
  results_pd0[[paste("pd0", pd0, sep = "_")]] <- result
}

results_list$pd0 <- results_pd0
gc()
### pd_beta ###
pr_beta_values <- seq(0.5, 2, by = 0.02)  # Example sequence
results_pr_beta <- list()

for (pr_beta in pr_beta_values) {
  result <- predict_F(baseline_params$K, baseline_params$r, baseline_params$t_m, baseline_params$r_d, baseline_params$t_d, 
                      baseline_params$pd_beta, baseline_params$r_pd, pr_beta, baseline_params$r_pr, 
                      baseline_params$pd0, baseline_params$pr0, baseline_params$t0, baseline_params$b)
  results_pr_beta[[paste("pr_beta", pr_beta, sep = "_")]] <- result
}

results_list$pr_beta <- results_pr_beta
gc()
### r_pr ###
r_pr_values <- seq(0.01, 0.1, by = 0.01)  # Example sequence
results_r_pr <- list()

for (r_pr in r_pr_values) {
  result <- predict_F(baseline_params$K, baseline_params$r, baseline_params$t_m, baseline_params$r_d, baseline_params$t_d, 
                      baseline_params$pd_beta, baseline_params$r_pd, baseline_params$pr_beta, r_pr, 
                      baseline_params$pd0, baseline_params$pr0, baseline_params$t0, baseline_params$b)
  results_r_pr[[paste("r_pr", r_pr, sep = "_")]] <- result
}

results_list$r_pr <- results_r_pr
gc()
### pr0 ###
pr0_values <- seq(0.001, 0.01, by = 0.001)  # Example sequence
results_pr0 <- list()

for (pr0 in pr0_values) {
  result <- predict_F(baseline_params$K, baseline_params$r, baseline_params$t_m, baseline_params$r_d, baseline_params$t_d, 
                      baseline_params$pd_beta, baseline_params$r_pd, baseline_params$pr_beta, baseline_params$r_pr, 
                      baseline_params$pd0, pr0, baseline_params$t0, baseline_params$b)
  results_pr0[[paste("pr0", pr0, sep = "_")]] <- result
}

results_list$pr0 <- results_pr0

saveRDS(results_list, "../model_result/sensitivity_analysis_results_list.rds")


if (1) {
  ###Plot How F Varies with different K###
  results_list = readRDS("../model_result/sensitivity_analysis_results_list.rds")
  combined_df <- do.call(rbind, lapply(names(results_list$K), function(k_name) {
    df <- results_list$K[[k_name]]
    df$K_value <- k_name  # Add a column to identify the K value for each row
    return(df)
  }))
  
  combined_df$K_value <- as.factor(combined_df$K_value)
  p_K <-ggplot(combined_df, aes(x = year, y = pre, color = K_value)) + 
    geom_line() + 
    labs(title = "Impact of varying K on Fb over years",
         x = "Year",
         y = "Predicted Fb values",
         color = "K Value") +
    theme_bw()
  print(p_K)
  ggsave("../model_graph/sensitivity_analysis/K.png", plot = p_K, width = 10, height = 6)
  
  ###Plot How F Varies with different r###
  # Assuming each data frame in results_list$r has 'year', 'pre', and 'invasion' columns
  combined_df_r <- do.call(rbind, lapply(names(results_list$r), function(r_name) {
    df <- results_list$r[[r_name]]
    df$r_value <- r_name  # Add a column to identify the r value for each row
    return(df)
  }))
  
  # Make sure 'r_value' is a factor for plotting purposes
  combined_df_r$r_value <- as.factor(combined_df_r$r_value)
  p_r <- ggplot(combined_df_r, aes(x = year, y = pre, color = r_value)) + 
    geom_line() + 
    labs(title = "Impact of varying r on Fb over years",
         x = "Year",
         y = "Predicted Fb values",
         color = "r Value") +
    theme_bw()
  print(p_r)
  ggsave("../model_graph/sensitivity_analysis/r.png", plot = p_r, width = 10, height = 6)
  
  
  ###Plot How F Varies with different pd_beta###
  # Assuming each data frame in results_list$pd_beta has 'year', 'pre', and 'invasion' columns
  combined_df_pd_beta <- do.call(rbind, lapply(names(results_list$pd_beta), function(pd_beta_name) {
    df <- results_list$pd_beta[[pd_beta_name]]
    df$pd_beta_value <- pd_beta_name  # Add a column to identify the pd_beta value for each row
    return(df)
  }))
  
  # Convert pd_beta_value to a factor for plotting
  combined_df_pd_beta$pd_beta_value <- as.factor(combined_df_pd_beta$pd_beta_value)
  p_pd_beta <- ggplot(combined_df_pd_beta, aes(x = year, y = pre, color = pd_beta_value)) + 
    geom_line() + 
    labs(title = "Impact of varying pd_beta on Fb over years",
         x = "Year",
         y = "Predicted Fb values",
         color = "pd_beta Value") +
    theme_bw()
  print(p_pd_beta)
  ggsave("../model_graph/sensitivity_analysis/pd_beta.png", plot = p_pd_beta, width = 10, height = 6)
  
  # Assuming each data frame in results_list$r_pd has 'year', 'pre', and 'invasion' columns
  combined_df_r_pd <- do.call(rbind, lapply(names(results_list$r_pd), function(r_pd_name) {
    df <- results_list$r_pd[[r_pd_name]]
    df$r_pd_value <- r_pd_name  # Add a column to identify the r_pd value for each row
    return(df)
  }))
  
  # Convert r_pd_value to a factor for plotting
  combined_df_r_pd$r_pd_value <- as.factor(combined_df_r_pd$r_pd_value)
  p_r_pd <- ggplot(combined_df_r_pd, aes(x = year, y = pre, color = r_pd_value)) + 
    geom_line() + 
    labs(title = "Impact of varying r_pd on Fb over years",
         x = "Year",
         y = "Predicted Fb values",
         color = "r_pd Value") +
    theme_bw()
  
  print(p_r_pd)
  ggsave("../model_graph/sensitivity_analysis/r_pd.png", plot = p_r_pd, width = 10, height = 6)
  
  ###plot pd0###
  combined_df_pd0 <- do.call(rbind, lapply(names(results_list$pd0), function(pd0_name) {
    df <- results_list$pd0[[pd0_name]]
    df$pd0_value <- pd0_name  # Add a column to identify the pd0 value for each row
    return(df)
  }))
  
  # Convert pd0_value to a factor for plotting
  combined_df_pd0$pd0_value <- as.factor(combined_df_pd0$pd0_value)
  p_pd0 <- ggplot(combined_df_pd0, aes(x = year, y = pre, color = pd0_value)) + 
    geom_line() + 
    labs(title = "Impact of varying pd0 on Fb over years",
         x = "Year",
         y = "Predicted Fb values",
         color = "pd0 Value") +
    theme_bw()
  
  # Display the plot
  print(p_pd0)
  ggsave("../model_graph/sensitivity_analysis/pd0.png", plot = p_pd0, width = 10, height = 6)
  
  ###plot pr_beta###
  combined_df_pr_beta <- do.call(rbind, lapply(names(results_list$pr_beta), function(pr_beta_name) {
    df <- results_list$pr_beta[[pr_beta_name]]
    df$pr_beta_value <- pr_beta_name  # Add a column to identify the pr_beta value for each row
    return(df)
  }))
  
  # Convert pr_beta_value to a factor for plotting
  combined_df_pr_beta$pr_beta_value <- as.factor(combined_df_pr_beta$pr_beta_value)
  p_pr_beta <- ggplot(combined_df_pr_beta, aes(x = year, y = pre, color = pr_beta_value)) + 
    geom_line() + 
    labs(title = "Impact of varying pr_beta on Fb over years",
         x = "Year",
         y = "Predicted Fb values",
         color = "pr_beta Value") +
    theme_bw()
  
  # Display the plot
  print(p_pr_beta)
  ggsave("../model_graph/sensitivity_analysis/pr_beta.png", plot = p_pr_beta, width = 10, height = 6)
  
  combined_df_r_pr <- do.call(rbind, lapply(names(results_list$r_pr), function(r_pr_name) {
    df <- results_list$r_pr[[r_pr_name]]
    df$r_pr_value <- r_pr_name  # Add a column to identify the r_pr value for each row
    return(df)
  }))
  
  # Convert r_pr_value to a factor for plotting
  combined_df_r_pr$r_pr_value <- as.factor(combined_df_r_pr$r_pr_value)
  p_r_pr <- ggplot(combined_df_r_pr, aes(x = year, y = pre, color = r_pr_value)) + 
    geom_line() + 
    labs(title = "Impact of varying r_pr on Fb over years",
         x = "Year",
         y = "Predicted Fb values",
         color = "r_pr Value") +
    theme_bw()
  
  # Display the plot
  print(p_r_pr)
  ggsave("../model_graph/sensitivity_analysis/r_pr.png", plot = p_r_pr, width = 10, height = 6)
  
  combined_df_pr0 <- do.call(rbind, lapply(names(results_list$pr0), function(pr0_name) {
    df <- results_list$pr0[[pr0_name]]
    df$pr0_value <- pr0_name  # Add a column to identify the pr0 value for each row
    return(df)
  }))
  
  # Convert pr0_value to a factor for plotting
  combined_df_pr0$pr0_value <- as.factor(combined_df_pr0$pr0_value)
  p_pr0 <- ggplot(combined_df_pr0, aes(x = year, y = pre, color = pr0_value)) + 
    geom_line() + 
    labs(title = "Impact of varying pr0 on Fb over years",
         x = "Year",
         y = "Predicted Fb values",
         color = "pr0 Value") +
    theme_bw()
  
  # Display the plot
  print(p_pr0)
  ggsave("../model_graph/sensitivity_analysis/pr0.png", plot = p_pr0, width = 10, height = 6)
}