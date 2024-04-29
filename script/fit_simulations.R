rm(list=ls());gc()
#graphics.off()
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(foreach)
library(doParallel)

params = FALSE
simulate = FALSE
fit = FALSE
prediction = TRUE

invasion_growth <- function(t, K, r, t_m, r_d, t_d) {
  K / (1 + exp((-r) * (t - t_m)))
  # if (t<t_d) {
  #   K / (1 + exp((-r) * (t - t_m)))
  # }else{
  #   K / (1 + exp((-r) * (t - t_m))) - exp(r_d * (t-t_d))
  #   #K / (1 + exp((-r) * (t - t_m))-) - exp(r_d * (t-t_d))
  # }
}

#Simulate that the probability of detection increases logistically due to abundance
sigmoid_function <- function(t, yoi, beta, rate, p0) {
  if ( beta == 0 && rate == 0 && p0 == 0){
    return(1)
  }else{
    return(1-exp(-rate * (t - yoi)^beta - p0))
  }
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
loss <- function(params, species_count_df) {
  params1=null_params
  params1[names(params)]=params
  
  K <- params1[1]
  r <- params1[2]
  t_m <- params1[3]
  r_d <- 0
  t_d <- 0
  pd_beta <- params1[4]
  r_pd <- params1[5]
  pd0 <- params1[6]
  pr_beta <- params1[7]
  r_pr <- params1[8]
  pr0 <- params1[9]
  
  years <- t0:b
  I <- sapply(years, invasion_growth, K, r, t_m, r_d, t_d)
  
  Fb_values <- sapply(years, function(a) calculate_Fb(a, I, pd_beta, r_pd, pr_beta, r_pr, pd0, pr0, t0, b))
  
  obs <- species_count_df$pre
  ls = sum((Fb_values - obs)^2)
  # print(params)
  # print(ls)
  if (is.na(ls) || !is.finite(ls)) {
    ls = 1e12
  }
  return(ls)
}

predict = function(params, b){
  K <- params[1]
  r <- params[2]
  t_m <- params[3]
  r_d <- 0
  t_d <- 0
  pd_beta <- params[4]
  r_pd <- params[5]
  pd0 <- params[6]
  pr_beta <- params[7]
  r_pr <- params[8]
  pr0 <- params[9]
  
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

plot_pre_over_year <- function(prediction_df, results_df, index) {
  # Combine the two data frames by year
  combined_df <- merge(prediction_df, results_df, by = "year", suffixes = c("_pred", "_res"))
  
  # Plotting
  p <- ggplot(combined_df, aes(x = year)) +
    geom_line(aes(y = pre_pred, color = "Predicted")) +
    geom_line(aes(y = pre_res, color = "Results")) +
    labs(title = paste("Parameter Set", index), y = "Pre", x = "Year") +
    scale_color_manual(values = c("Predicted" = "blue", "Results" = "red")) +
    theme_bw()
  
  return(p)
}

identify_outliers <- function(x) {
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  iqr <- q3 - q1
  lower_bound <- q1 - 1.5 * iqr
  upper_bound <- q3 + 1.5 * iqr
  x < lower_bound | x > upper_bound
}


set.seed(123) # Ensure reproducibility
n <- 15 # Number of random selections
t0 = 1990
b = 2020
lp = 0.75
hp = 1.25
if(params){
  baseline_params <- list(K = 1.000e+03,
                          r = 2e-02,
                          t_m = 1.872e+03,   # Midpoint of logistic growth
                          # r_d = 2.000e-01,
                          # t_d = 1.995e+03,
                          pd_beta = 1.599e+00,
                          r_pd = 0.000001e+00,
                          pd0 = 5.407e-03,
                          pr_beta = 0,
                          r_pr = 0,
                          pr0 = 0)
  
  #generate the different cases of simulated real world
  
  
  random_params <- replicate(n, {
    list(
      K = runif(1, 500, 1500),
      r = baseline_params$r * runif(1, lp, hp),
      t_m = runif(1, 1800, 2000),   # Midpoint of logistic growth
      # r_d = baseline_params$r_d * runif(1, lp, hp),
      # t_d = baseline_params$t_d * runif(1, lp, hp),
      pd_beta = baseline_params$pd_beta * runif(1, lp, hp),
      r_pd = baseline_params$r_pd * runif(1, lp, hp),
      pd0 = baseline_params$pd0 * runif(1, lp, hp),
      pr_beta = baseline_params$pr_beta * runif(1, lp, hp),
      r_pr = baseline_params$r_pr * runif(1, lp, hp),
      pr0 = baseline_params$pr0 * runif(1, lp, hp)
    )
  }, simplify = FALSE)
  saveRDS(random_params, paste0("../model_result/fit_simulation_params_", n, "_", t0, "_", lp*100, "-", hp*100, ".rds"))
}
random_params = readRDS(paste0("../model_result/fit_simulation_params_", n, "_", t0, "_", lp*100, "-", hp*100, ".rds")) 


if (simulate){
  results_n <- mapply(function(params, i, b) {
    # Construct the name using the index i
    name <- paste("n", i, sep = "_")
    
    # Convert the list of parameters to a numeric vector in the correct order
    params_vector <- unlist(params)
    
    # Call the predict function with the parameters vector and b
    result <- predict(params_vector, b)
    
    # Return a named list with the result
    setNames(list(result), name)
  }, random_params, seq_along(random_params), MoreArgs = list(b = b), SIMPLIFY = FALSE)
  
  saveRDS(results_n, paste0("../model_result/fit_simulation_results_", n, "_", t0, "_", lp*100, "-", hp*100, ".rds"))
}

results_n = readRDS(paste0("../model_result/fit_simulation_results_", n, "_", t0, "_", lp*100, "-", hp*100, ".rds"))

if (fit){
  init_params = c(K = 800, r = 2e-02, t_m = 1.800e+03, pd_beta = 1, r_pd = 1e-03, pd0 = 0.0050)
  null_params = c(K = 800, r = 2e-02, t_m = 1.800e+03, pd_beta = 0, r_pd = 0, pd0 = 0, pr_beta = 0, r_pr = 0, pr0 = 0)
  # lower_bounds <- c(1, 0.000001, 0, 0.000001, 0,  -Inf, 0, 0, -Inf, 0, 0)
  # upper_bounds <- c(1e5, 1, 2050, 1, 2050,  Inf, Inf, Inf, Inf, Inf, Inf)
  # numCores <- detectCores()
  # cl <- makeCluster(numCores - 1)
  # registerDoParallel(cl)
  registerDoParallel()
  options(cores=15) 
  fitted_n <- foreach(i = seq_along(results_n), .combine = 'c') %dopar% {
    result <- results_n[[i]]
    name <- paste("fitted", i, sep = "_")
    species_count_df <- result[[1]]  # Ensure this is correctly extracting species_count_df
    l <- 1e10
    not_best_fit <- TRUE
    # init_params <- random_params[[i]]  # Assuming this needs to be specific for each iteration
    # opt_results <- optim(init_params, function(params) loss(params, species_count_df), control = list(maxit = 10))
    while (l > 1e-5 && not_best_fit) {
      old_init_params <- init_params
      opt_results <- optim(init_params, function(params) loss(params, species_count_df), control = list(maxit = 2000))
      # opt_results <- optim(init_params, function(params) loss(params, species_count_df), method = "CG", control = list(maxit = 20000))
      # opt_results <- optim(init_params, function(params) loss(params, species_count_df), method = "BFGS", control = list(maxit = 20000))
      # opt_results <- optim(init_params, function(params) loss(params, species_count_df), method = "SANN", control = list(maxit = 20000, temp = 15, tmax = 50))
      # opt_results <- optim(init_params, function(params) loss(params, species_count_df), control = list(maxit = 20000))
      # opt_results <- optim(init_params, function(params) loss(params, species_count_df), control = list(maxit = 2000))
      # opt_results <- optim(init_params, function(params) loss(params, species_count_df), method = "L-BFGS-B", lower = lower_bounds, upper = upper_bounds, control = list(maxit = 20000))
      # cat("new l\n")
      l <- opt_results$value
      init_params <- opt_results$par
      if (identical(old_init_params, init_params)){
        not_best_fit <- FALSE
      }
    }

    # Return named list of optimization results
    setNames(list(opt_results), name)
  }
  # stopCluster(cl)
  # fitted_n <- mapply(function(result, i, b){
  #   name <- paste("fitted", i, sep = "_")
  #   init_params = unlist(random_params[[i]])
  #   species_count_df = result[[1]]  # Ensure this is correctly extracting species_count_df
  #   l = 1e10
  #   not_best_fit = TRUE
  #   opt_results <- optim(init_params, function(params) loss(params, species_count_df), control = list(maxit = 10))
  #   # while (l > 1e-5 && not_best_fit) {
  #   #   old_init_params = init_params
  #   #   opt_results <- optim(init_params, function(params) loss(params, species_count_df), control = list(maxit = 1000))
  #   #   print("new l")
  #   #   l = opt_results$value
  #   #   init_params = opt_results$par
  #   #   if (identical(old_init_params, init_params)){
  #   #     not_best_fit = FALSE
  #   #   }
  #   # }
  
  #   # mcmc
  #   #print(result)
  #   # print(init_params), parscale = parscale_values
  #   setNames(list(opt_results), name)
  # }, results_n, seq_along(results_n), MoreArgs = list(b = b), SIMPLIFY = FALSE)
  
  saveRDS(fitted_n, paste0("../model_result/fit_simulation_fitted_results_", n, "_", t0, "_", lp*100, "-", hp*100, ".rds"))
}

fitted_n = readRDS(paste0("../model_result/fit_simulation_fitted_results_", n, "_", t0, "_", lp*100, "-", hp*100, ".rds"))

# Assuming random_params is a list of lists with each inner list containing parameter values
# Convert random_params to a dataframe
sim_params_df <- do.call(rbind, lapply(random_params, function(x) unlist(x)))

# For fitted_n, extract optimized parameters and convert to a dataframe
fit_params_matrix <- t(sapply(fitted_n, function(x) x$par))
fit_params_df <- as.data.frame(fit_params_matrix)

# Ensure column names for sim_params_df and fit_params_df match and indicate simulated vs. fitted
# Assuming the order and number of parameters are the same for sim and fit

colnames(sim_params_df) <- paste0("sim-", colnames(sim_params_df))
colnames(fit_params_df) <- paste0("fit-", colnames(fit_params_df))


# Now both dataframes should have matching column names for corresponding parameters

# Combine the simulated and fitted parameters dataframes
combined_params_df <- cbind(sim_params_df, fit_params_df)

# combined_params_df <- combined_params_df %>%
#   filter(!if_any(where(is.numeric), remove_outliers))

# Convert the dataframe from wide to long format
long_params_df <- pivot_longer(combined_params_df, 
                               cols = everything(), 
                               names_to = c(".value", "parameter"), 
                               names_sep = "-")

long_params_df <- long_params_df %>%
  group_by(parameter) %>%
  mutate(outlier = identify_outliers(sim) | identify_outliers(fit)) %>%
  ungroup()

long_params_df <- long_params_df %>%
  filter(!outlier)

# 'long_params_df' now has a column for 'parameter' indicating the parameter name,
# and two columns for 'sim' and 'fit' representing the simulated and fitted values.
# Plot simulated vs. fitted values for each parameter
ggplot(long_params_df, aes(x = fit, y = sim, color = parameter)) +
  geom_point() +  # Add points
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +  # Add y=x line
  facet_wrap(~ parameter, scales = "free") +  # Create a separate plot for each parameter
  theme_bw() +
  labs(x = "Fitted Value", y = "Simulated Value", title = "Simulated vs. Fitted Values for Each Parameter") +
  theme(legend.position = "none")  # Hide legend because the facet titles already indicate the parameter
#ggsave("../model_graph/Simulated vs. Fitted Values for Each Parameter_95-10.png")

if (prediction) {
  prediction_n <- mapply(function(params, i, b) {
    # Construct the name using the index i
    name <- paste("n", i, sep = "_")
    
    # Convert the list of parameters to a numeric vector in the correct order
    params_vector <- unlist(params[[1]]$par)
    
    # Call the predict function with the parameters vector and b
    result <- predict(params_vector, b)
    
    # Return a named list with the result
    setNames(list(result), name)
  }, fitted_n, seq_along(fitted_n), MoreArgs = list(b = b), SIMPLIFY = FALSE)
  saveRDS(prediction_n, paste0("../model_result/fit_simulation_predicted_results_", n, "_", t0, "_", lp*100, "-", hp*100, ".rds"))
}

prediction_n = readRDS(paste0("../model_result/fit_simulation_predicted_results_", n, "_", t0, "_", lp*100, "-", hp*100, ".rds"))
plots_list <- lapply(1:length(prediction_n), function(i) {
  prediction_df <- prediction_n[[i]][[1]]
  results_df <- results_n[[i]][[1]]
  plot_pre_over_year(prediction_df, results_df, i)
})

# If you want to display all plots in one grid (consider this might be impractical with 50 plots)
# grid.arrange(grobs = plots_list, ncol = 5)

# Alternatively, save each plot to a file
for (i in seq_along(plots_list)) {
  # ggsave(paste0("../model_graph/simulation_", n, "_", t0, "_", lp*100, "-", hp*100, "/plot_", i, ".png", sep = ""), plots_list[[i]], width = 10, height = 6)
}

