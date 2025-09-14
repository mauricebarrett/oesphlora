library(testthat)

test_that("basic math works", {
  expect_equal(2 + 2, 4)
})

test_that("string uppercasing works", {
  expect_equal(toupper("pixi"), "PIXI")
})
