test_that("base2int works", {
  expect_equal(sum(base2int("ACGT")), 6)
})
