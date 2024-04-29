rm(list=ls());gc()
#graphics.off()
library(ggplot2)
library(dplyr)
library(tidyr)

invasion_probability <- function(t, pi) {
  pi
}

# Define the sigmoid functions for Pd and Pr
sigmoid_function <- function(t, yoi, beta, rate, p0) {
  1-exp(-rate * (t - yoi)^beta - p0)
}

lh_record <- function(a, IP, pd_beta, r_pd, pr_beta, r_pr, pd0, pr0, t0, b) {
  Fb <- 0
  for (i in t0:a) {
    pi_product <- if (i == t0){
      1
    }else{
      prod(1 - IP[1:(i-t0)])
    }
    pi_current <- IP[i - t0 + 1]
    pd_product <- if (i < a) {
      prod(1 - sapply(i:(a-1), sigmoid_function, i, pd_beta, r_pd, pd0))
    } else {
      1
    }
    pd_current <- sigmoid_function(a, i, pd_beta, r_pd, pd0)
    pr_product <- 1-prod(1 - sapply(a:b, sigmoid_function, a, pr_beta, r_pr, pr0))
    Fb <- Fb + pi_current * pi_product * pd_product * pd_current * pr_product
  }
  return(Fb)
}

lh_no_record <- function(a, IP, pd_beta, r_pd, pr_beta, r_pr, pd0, pr0, t0, b){
  Fb <- prod(1 - IP[1:(a-t0+1)])
  for (i in t0:a) {
    #case 1: no invasion
    # no_invade = prod(1 - IP[1:(i-t0+1)])
    #case 2: no detection
    pi_product <- if (i == t0){
      1
    }else{
      prod(1 - IP[1:(i-t0)])
    }
    pi_current <- IP[i - t0 + 1]
    # no_pd = prod(1 - sapply(i:a, sigmoid_function, i, pd_beta, r_pd, pd0))
    # no_dect = pi_product * pi_current * no_pd
    #case 3: no report
    pd_product <- if (i < a) {
      prod(1 - sapply(i:(a-1), sigmoid_function, i, pd_beta, r_pd, pd0))
    } else {
      1
    }
    pd_current <- sigmoid_function(a, i, pd_beta, r_pd, pd0)
    no_dect = pi_product * pi_current * (1-pd_product*pd_current)
    no_pr <- prod(1 - sapply(a:b, sigmoid_function, a, pr_beta, r_pr, pr0))
    no_report = pi_product*pi_current*pd_product*pd_current*no_pr
    Fb <- Fb + no_dect + no_report
  }
  return(Fb)
}

pi = 0.0001
t0 = 1600
b = 2016
years <- t0:b
IP <- sapply(years, invasion_probability, pi)



baseline_params <- list(K = 1.000e+03,
                        r = 2e-02,
                        t_m = 1.872e+03,   # Midpoint of logistic growth
                        r_d = 2.000e-01,
                        t_d = 1.995e+03,
                        pd_beta = 1.599e+00,
                        r_pd = 0.100e+00,
                        pd0 = 5.407e-03,
                        pr_beta = 1.368e+00,
                        r_pr = 1.438e-02,
                        pr0 = 8.712e-03,
                        t0 = 1600, 
                        b = 2016)

set.seed(123) # Ensure reproducibility
n <- 100 # Number of random selections

random_params <- replicate(n, {
  list(
    K = baseline_params$K * runif(1, 0.9, 1.1),
    r = baseline_params$r * runif(1, 0.9, 1.1),
    t_m = baseline_params$t_m * runif(1, 0.9, 1.1),   # Midpoint of logistic growth
    r_d = baseline_params$r_d * runif(1, 0.9, 1.1),
    t_d = baseline_params$t_d * runif(1, 0.9, 1.1),
    pd_beta = baseline_params$pd_beta * runif(1, 0.9, 1.1),
    r_pd = baseline_params$r_pd * runif(1, 0.9, 1.1),
    pd0 = baseline_params$pd0 * runif(1, 0.9, 1.1),
    pr_beta = baseline_params$pr_beta * runif(1, 0.9, 1.1),
    r_pr = baseline_params$r_pr * runif(1, 0.9, 1.1),
    pr0 = baseline_params$pr0 * runif(1, 0.9, 1.1),
    a = sample(1600:2016, 1) # Randomly choosing a year for 'a'
  )
}, simplify = FALSE)


# Adjusted part of the calculation loop
results <- lapply(random_params, function(params) {
  list(
    lh_no_record = lh_no_record(params$a, IP, params$pd_beta, params$r_pd, params$pr_beta, params$r_pr, params$pd0, params$pr0, t0, b),
    # lh_record = lh_record(params$a, IP, params$pd_beta, params$r_pd, params$pr_beta, params$r_pr, params$pd0, params$pr0, t0, b),
    one_minus_lh_record = 1 - lh_record(params$a, IP, params$pd_beta, params$r_pd, params$pr_beta, params$r_pr, params$pd0, params$pr0, t0, b)
  )
})

results_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    iteration = i,
    lh_no_record = results[[i]]["lh_no_record"],
    # lh_record = results[[i]]["lh_record"],
    one_minus_lh_record = results[[i]]["one_minus_lh_record"]
  )
}))


results_df <- pivot_longer(results_df, cols = c("lh_no_record", "one_minus_lh_record"), names_to = "type", values_to = "value")

ggplot(results_df, aes(x = iteration, y = value, color = type)) +
  geom_line() +
  labs(title = "Likelihoods of No Record vs. 1 - Record", x = "Iteration", y = "Likelihood") +
  theme_bw()

# Assuming 'results' is your list
results_df <- do.call(rbind, lapply(results, function(x) {
  data.frame(lh_no_record = x$lh_no_record, one_minus_lh_record = x$one_minus_lh_record)
}))


ggplot(results_df, aes(x = lh_no_record, y = one_minus_lh_record)) +
  geom_point() +  # This adds the scatter plot points
  geom_line() +   # This connects the points with lines
  labs(title = "Relationship between lh_no_record and 1 - lh_record",
       x = "lh_no_record",
       y = "1 - lh_record") +
  theme_minimal()


