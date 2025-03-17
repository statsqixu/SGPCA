test_that("simulator1 works", {
  data <- simulator1(n = 100, G = 100, C = 5, seed = 123)
  expect_equal(nrow(data$X), 100)
  expect_equal(ncol(data$X), 500)
  expect_equal(length(unique(data$group_label)), 100)
  expect_equal(max(table(data$group_label)), 5)
  expect_equal(norm(data$pc1, type = "2"), 1)
})

test_that("simulator2 works", {
  data <- simulator2(n = 100, G = 100, cond = 1, seed = 123)
  expect_equal(nrow(data$X), 100)
  expect_equal(ncol(data$X), 1000)
  expect_equal(length(unique(data$group_label)), 100)
  expect_equal(max(table(data$group_label)), 10)
  expect_equal(norm(data$pc1, type = "2"), 1)
})

test_that("simulator3 works", {
  data <- simulator3(n = 100, G = 100, seed = 123)
  expect_equal(nrow(data$X), 100)
  expect_equal(ncol(data$X), 1000)
  expect_equal(length(unique(data$group_label)), 100)
  expect_equal(max(table(data$group_label)), 10)
  expect_equal(norm(data$pc1, type = "2"), 1)
  expect_equal(norm(data$pc2, type = "2"), 1)
  expect_equal(norm(data$pc3, type = "2"), 1)
})
