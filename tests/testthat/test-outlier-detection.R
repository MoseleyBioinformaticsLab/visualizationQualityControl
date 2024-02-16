data(grp_cor_data)

grp1 <- grp_cor_data$data[, 1:10]
grp2 <- grp_cor_data$data[, 11:20]

grp1_class <- grp_cor_data$class[1:10]
grp_class <- grp_cor_data$class

all_samples_no_out <- cbind(grp1, grp2)
all_samples_out = all_samples_no_out
all_samples_out[, 3] <- grp2[, 3]
all_samples_out[, 15] <- grp1[, 3]

cor_no_out <- cor(all_samples_no_out)
cor_out = cor(all_samples_out)

med_no_out = median_correlations(cor_no_out, grp_class)
med_out = median_correlations(cor_out, grp_class)
set.seed(1234)
no_outlier <- (vapply(seq(1, 20), function(x){rnorm(200, 0, 1)}, numeric(200)))
has_outlier = (vapply(seq(1, 2), function(x){rnorm(200, 5, 0.5)}, numeric(200)))

noout_frac = outlier_fraction(no_outlier, grp_class)

has_out2 = no_outlier
has_out2[, 1] = has_outlier[, 1]
has_out2[, 12] = has_outlier[, 2]
hasout_frac = outlier_fraction(has_out2, grp_class)

test_that("no outliers works", {
  
  noout_samples_cor = determine_outliers(median_correlations = med_no_out)
  expect_equal(sum(noout_samples_cor$outlier), 0)

  noout_samples_frac = determine_outliers(outlier_fraction = noout_frac)
  expect_equal(sum(noout_samples_frac$outlier), 0)
  
  noout_samples_both = determine_outliers(med_no_out, noout_frac)
  expect_equal(sum(noout_samples_both$outlier), 0)
  
})

test_that("two outliers correlation works", {
  out2_samples_cor = determine_outliers(median_correlations = med_out)
  
  expect_equal(which(out2_samples_cor$outlier), c(3, 15))
  
  out2_samples_cor0 = determine_outliers(median_correlations = med_out, cor_weight = 0)
  expect_equal(sum(out2_samples_cor0$outlier), 0)
})

test_that("two outliers fraction works", {
  out2_samples_frac = determine_outliers(outlier_fraction = hasout_frac)
  expect_equal(which(out2_samples_frac$outlier), c(1, 12))
  
  out2_samples_frac0 = determine_outliers(outlier_fraction = hasout_frac, frac_weight = 0)
  expect_equal(sum(out2_samples_frac0$outlier), 0)
})

test_that("combined works", {
  out4_samples = determine_outliers(med_out, hasout_frac)
  expect_equal(which(out4_samples$outlier), c(1, 3, 12, 15))
  
  out4_samples_cor = determine_outliers(med_out, hasout_frac, frac_weight = 0)
  expect_equal(which(out4_samples_cor$outlier), c(3, 15))
  
  out4_samples_frac = determine_outliers(med_out, hasout_frac, cor_weight = 0)
  expect_equal(which(out4_samples_frac$outlier), c(1, 12))
})
