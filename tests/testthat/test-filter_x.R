context("filter_x")

#### Setup ---------------------------------------------------------------------
library(dplyr)
library(tidyr)

# Simulate data with missing values
# See https://stackoverflow.com/questions/50528719
set.seed(2)

N <- 20
prc_missing = 0.20

data <- data.frame(
  protein = 1:N,
  sample1 = rnorm(N, -5:5),
  sample2 = rnorm(N, -5:5),
  sample3 = rnorm(N, -5:5),
  sample4 = rnorm(N, -5:5)
)

data_na <- data %>%
  pivot_longer(cols = -protein, names_to = "sample", values_to = "logfc") %>%
  mutate(r = runif(nrow(.)),
         logfc = ifelse(r <= prc_missing, NA, logfc)) %>%
  select(-r) %>%
  pivot_wider(names_from = "sample", values_from = "logfc")

data_zero <- data_na %>%
  mutate_all(~ replace(., is.na(.), 0))

#### Tests ---------------------------------------------------------------------

test_that("filter_na works", {
  data_na_filt <- filter_na(
    data_na,
    op = "==",
    setop = union,
    regex = "sample",
    value = 2)

  expect_equal(nrow(data_na_filt), 2)
})

test_that("filter_na with multiple regex works", {
  data_na_filt <- filter_na(
    data_na,
    op = "==",
    setop = intersect,
    regex = c("sample[1-2]", "sample[3-4]"),
    value = c(1, 0)
    )

  expect_equal(nrow(data_na_filt), 5)
})

test_that("filter_zero works", {
  data_zero_filt <- filter_zero(
    data_zero,
    op = "==",
    setop = union,
    regex = "sample",
    value = 2)

  expect_equal(nrow(data_zero_filt), 2)
})

test_that("filter_zero with multiple regex works", {
  data_zero_filt <- filter_zero(
    data_zero,
    op = "==",
    setop = intersect,
    regex = c("sample[1-2]", "sample[3-4]"),
    value = c(1, 0)
  )

  expect_equal(nrow(data_zero_filt), 5)
})

test_that("filter_val works", {
  data_filt <- filter_val(
    data,
    op = ">=",
    setop = union,
    regex = "sample",
    value = 1)

  expect_equal(nrow(data_filt), 8)
})

test_that("filter_val with multiple regex works", {
  data_filt <- filter_val(
    data,
    op = ">=",
    setop = intersect,
    regex = c("sample[1-2]", "sample[3-4]"),
    value = c(1, 2)
  )

  expect_equal(nrow(data_filt), 7)
})

#### Sanity checks -------------------------------------------------------------
