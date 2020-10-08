test_that("parallel vs sequential Kt give same results", {
  data(grp_cor_data)
  
  forkoption = getOption("future.fork.enable")
  options(future.fork.enable = TRUE)
  
  library(furrr)
  plan(multiprocess(workers = 2))
  
  grp1 <- grp_cor_data$data
  
  seq_cor = visualizationQualityControl:::visqc_ici_kendallt_ref(t(grp1), exclude_0 = FALSE, scale_max = FALSE, diag_good = FALSE)
  par_cor = visqc_ici_kendallt(t(grp1), exclude_0 = FALSE, scale_max = FALSE, diag_good = FALSE)
  
  options(future.fork.enable = forkoption)
  
  expect_equal(seq_cor, par_cor)
})
