test_that("multiplication works", {
  expect_equal(test_fun(3), 4)
  expect_identical(test_fun(2), 3)
})
