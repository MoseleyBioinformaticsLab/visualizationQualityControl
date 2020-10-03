test_that("parallel vs sequential give same results", {
  data(grp_cor_data)
  
  library(furrr)
  
  grp1 <- grp_cor_data$data
  
  seq_cor = visqc_ici_kendallt(t(grp1), exclude_0 = FALSE, scale_max = FALSE, diag_good = FALSE)
  par_cor = visqc_ici_kendallt_splitup(t(grp1), exclude_0 = FALSE, scale_max = FALSE, diag_good = FALSE)
  
  expect_equal(seq_cor, par_cor)
})
