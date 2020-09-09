test_that("basic kendall-tau matches base R", {
  x = seq(1, 10)
  y = seq(1, 10)
  expect_equal(kendallt(x, y), cor(x, y, method = "kendall"))
  
  y[2] = 15
  expect_equal(kendallt(x, y), cor(x, y, method = "kendall"))
  
  y = seq(10, 1) # should give -1
  expect_equal(kendallt(x, y), cor(x, y, method = "kendall"))
  
  y[2] = 15
  expect_equal(kendallt(x, y), cor(x, y, method = "kendall"))
})
