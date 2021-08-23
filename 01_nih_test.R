

source("R/00_global_functions.R")

# SigOpt API token
Sys.setenv(
  SIGOPT_API_TOKEN =
    scan("sigopt-api-token", what = "character", sep = "\n", quiet = TRUE)
)

nih_dtm <- 
  textmineR::CreateDtm(
    doc_vec = textmineR::nih_sample$ABSTRACT_TEXT,
    stopword_vec = c()
  )

sample_zipf <- 
  run_sigopt(
    dtm_true = nih_dtm,
    exp_name = "nih-test-zipf",
    observations = 20,
    optimization_var = "zipf"
  )

sample_heaps <- 
  run_sigopt(
    dtm_true = nih_dtm,
    exp_name = "nih-test-heaps",
    observations = 20,
    optimization_var = "heaps"
  )

sample_taylor <- 
  run_sigopt(
    dtm_true = nih_dtm,
    exp_name = "nih-test-tyalor",
    observations = 20,
    optimization_var = "taylor"
  )

sample_nv <- 
  run_sigopt(
    dtm_true = nih_dtm,
    exp_name = "nih-test-nv",
    observations = 20,
    optimization_var = "Nv"
  )


nih_test <- bind_rows(
  sample_zipf,
  sample_heaps,
  sample_taylor,
  sample_nv
)

# experiment with distance
v1 <- nih_test %>%
  select(
    Nv,
    true_zipf_intercept,
    true_zipf_coef,
    true_heaps_intercept,
    true_heaps_coef,
    true_taylor_intercept,
    true_taylor_coef
  ) %>%
  .[1, ] %>%
  as.numeric()

nih_dist <- nih_test %>%
  select(
    alpha_shape1,
    alpha_shape2,
    alpha_sum,
    eta_shape,
    eta_sum,
    Nk,
    Nk,
    Nv
  ) %>% 
  mutate(
    exp_dist = apply(
      nih_test %>%
        select(
          sim_Nv,
          zipf_intercept,
          zipf_coef,
          heaps_intercept,
          heaps_coef,
          taylor_intercept,
          taylor_coef
          ),1, function(x) {
              vec <- (x - v1) / v1
              
              d <- dist(rbind(vec, rep(0, length(vec)))) %>%
                as.matrix() %>%
                .[1, 2]
              
              d
            })
  )


save(
  nih_test,
  nih_dist,
  file = "data-derived/nih-test-results.RData"
)
