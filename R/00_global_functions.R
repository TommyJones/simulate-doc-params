
library(tidyverse)
library(tmsamples)
library(SigOptR)

# library(poweRlaw)

# draws a DTM from lda(ish) parameters
draw <- function(
  Nk,
  Nv,
  alpha_shape,
  alpha_sum,
  eta_shape,
  eta_sum,
  doc_lengths, # Nd implicit here
  threads = parallel::detectCores() - 1,
  verbose = TRUE
) {
  
  alpha <- rbeta(Nk, alpha_shape[1], alpha_shape[2])
  
  alpha <- alpha / sum(alpha) * alpha_sum
  
  eta <- tmsamples::generate_zipf(
    vocab_size = Nv, 
    magnitude = eta_sum, 
    zipf_par = abs(eta_shape)
  )
  
  param <- 
    tmsamples::sample_parameters(
      alpha = alpha,
      beta = eta,
      num_documents = length(doc_lengths)
    )
  
  dtm <- 
    tmsamples::sample_documents(
      theta = param$theta,
      phi = param$phi,
      doc_lengths = doc_lengths,
      threads = threads,
      verbose = verbose
    )
  
  dtm
}


# estimate parameters from zipfs, heaps, and taylors laws
# uses log-log regression for now. MLE would be better...
calc_law_params <- function(
  dtm
) {
  
  # get Zipf
  zipf_data <- tibble(
    rank = 1:ncol(dtm),
    freq = sort(Matrix::colSums(dtm), decreasing = TRUE)
  )
  
  zipf_model <- lm(
    I(log10(freq)) ~ I(log10(rank)),
    data = zipf_data %>% filter(freq > 0)
  )
  
  # get Heaps
  heaps_data <- tibble(
    n = Matrix::rowSums(dtm),
    v = Matrix::rowSums(dtm > 0)
  )
  
  heaps_model <- lm(
    I(log10(v)) ~ I(log10(n)),
    data = heaps_data
  )
  
  # get Taylor
  taylor_data <- tibble(
    mu = Matrix::rowMeans(dtm)
  ) %>%
    mutate(
      sigma = Matrix::rowSums((dtm - mu) ^ 2) / ncol(dtm)
    )
  
  taylor_model <- lm(
    I(log10(sigma)) ~ I(log10(mu)),
    data = taylor_data
  )
  
  # format output
  tibble(
    zipf_intercept = zipf_model$coefficients[1],
    zipf_coef = zipf_model$coefficients[2],
    heaps_intercept = heaps_model$coefficients[1],
    heaps_coef = heaps_model$coefficients[2],
    taylor_intercept = taylor_model$coefficients[1],
    taylor_coef = taylor_model$coefficients[2]
  )
  
}


# run the sigopt optimization loop and record results
# note this requires a sigopt api key
run_sigopt <- function(
  dtm_true,
  exp_name,
  optimization_var = c("zipf", "heaps", "taylor", "Nv"),
  observations = 100,
  threads = parallel::detectCores() - 1,
  verbose = TRUE
) {
  
  # declare global (w/i function scope) variables
  Nv <- ncol(dtm_true)
  
  Nd <- nrow(dtm_true)
  
  doc_lengths <- Matrix::rowSums(dtm_true)
  
  # format optimization metric varable
  if (optimization_var == "zipf") {
    metric_list <- list(
      list(name = "zipf_intercept_err", objective = "minimize", strategy = "optimize"),
      list(name = "zipf_coef_err", objective = "minimize", strategy = "optimize")
    )
  } else if (optimization_var == "heaps") {
    metric_list <- list(
      list(name = "heaps_intercept_err", objective = "minimize", strategy = "optimize"),
      list(name = "heaps_coef_err", objective = "minimize", strategy = "optimize")
    )
  } else if (optimization_var == "taylor") {
    metric_list <- list(
      list(name = "taylor_intercept_err", objective = "minimize", strategy = "optimize"),
      list(name = "taylor_coef_err", objective = "minimize", strategy = "optimize")
    )
  } else if (optimization_var == "Nv") {
    metric_list <- list(
      list(name = "Nv_err", objective = "minimize", strategy = "optimize")
    )
  } else {
    stop("optimization var must be one of 'zipf', 'heaps', 'taylor', or 'Nv'")
  }
  
  # create sigopt experiment
  exp <- list(
    name = exp_name,
    parameters = list(
      list(name = "alpha_shape1", type = "double", bounds = list(min = 0.05, max = 100.0)),
      list(name = "alpha_shape2", type = "double", bounds = list(min = 0.05, max = 100.0)),
      list(name = "alpha_sum", type = "double", bounds = list(min = 1.0, max = 100)),
      list(name = "eta_sum", type = "double", bounds = list(min = 10.0, max = Nv)),
      list(name = "Nk", type = "int", bounds = list(min = 2, max = 200))
    ),
    observation_budget = observations,
    metrics = metric_list,
    project = "dissertation"
  ) %>%
    SigOptR::create_experiment()
  
  # calculate power law pars
  true_pars <- calc_law_params(dtm_true) %>%
    mutate(
      Nd = Nd,
      Nv = Nv,
      mean_doc_length = mean(doc_lengths),
      sd_doc_length = sd(doc_lengths)
    ) %>%
    mutate(
      exp_name = exp_name
    )
  
  # run the optimization loop
  simulated_pars <- 
    seq_len(exp$observation_budget) %>%
    map(function(j){
      
      # get parameters from sigopt
      suggestion <- SigOptR::create_suggestion(exp$id)
      
      # record parameters in output data frame
      output <- tibble(
        alpha_shape1 = suggestion$assignments$alpha_shape1,
        alpha_shape2 = suggestion$assignments$alpha_shape2,
        alpha_sum = suggestion$assignments$alpha_sum,
        eta_shape = true_pars$zipf_coef,
        eta_sum = suggestion$assignments$eta_sum,
        Nk = suggestion$assignments$Nk
      )
      
      # simulate corpus
      dtm_sim <- draw(
        Nv = Nv,
        Nk = output$Nk,
        alpha_shape = c(output$alpha_shape1, output$alpha_shape2),
        alpha_sum = output$alpha_sum,
        eta_shape = output$eta_shape, 
        eta_sum = output$eta_sum,
        doc_lengths = doc_lengths, 
        threads = threads,
        verbose = verbose
      )
      
      # calculate simulation power law pars
      output <- 
        output %>%
        cbind(
          calc_law_params(dtm_sim)
        ) %>%
        as_tibble() %>%
        mutate(
          exp_name = exp_name,
          sim_Nv = sum(colSums(dtm_sim) > 0)
        )
      
      # record errors
      errors <- 
        tibble(
          zipf_intercept_err = abs(output$zipf_intercept - true_pars$zipf_intercept),
          zipf_coef_err = abs(output$zipf_coef - true_pars$zipf_coef),
          heaps_intercept_err = abs(output$heaps_intercept - true_pars$heaps_intercept),
          heaps_coef_err = abs(output$heaps_coef - true_pars$heaps_coef),
          taylor_intercept_err = abs(output$taylor_intercept - true_pars$taylor_intercept),
          taylor_coef_err = abs(output$taylor_coef - true_pars$taylor_coef),
          Nv_err = abs(output$sim_Nv - true_pars$Nv)
        )
      
      # get values to report
      if (optimization_var == "zipf") {
        values = list(
          list(name = "zipf_intercept_err", value = errors$zipf_intercept_err),
          list(name = "zipf_coef_err", value = errors$zipf_coef_err)
        )
      } else if (optimization_var == "heaps") {
        values = list(
          list(name = "heaps_intercept_err", value = errors$heaps_intercept_err),
          list(name = "heaps_coef_err", value = errors$heaps_coef_err)
        )
      } else if (optimization_var == "taylor") {
        values = list(
          list(name = "taylor_intercept_err", value = errors$taylor_intercept_err),
          list(name = "taylor_coef_err", value = errors$taylor_coef_err)
        )
      } else { # must be Nv
        values = list(
          list(name = "Nv_err", value = errors$Nv_err)
        )
      }
      
      
      # report errors to sigopt
      SigOptR::create_observation(exp$id, list(
        suggestion = suggestion$id,
        values = values
      ))
      
      # return output data frame
      output
    }) %>%
    bind_rows()
  
  # return results tibble
  results <-
    simulated_pars %>%
    full_join(
      true_pars %>%
        select(
          exp_name,
          Nd,
          Nv,
          mean_doc_length,
          sd_doc_length,
          true_zipf_intercept = zipf_intercept,
          true_zipf_coef = zipf_coef,
          true_heaps_intercept = heaps_intercept,
          true_heaps_coef = heaps_coef,
          true_taylor_intercept = taylor_intercept,
          true_taylor_coef = taylor_coef
        ),
      by = c("exp_name" = "exp_name")
    )
  
  results
}
