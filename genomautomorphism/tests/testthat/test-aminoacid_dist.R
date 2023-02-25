test_that("aminoacid_dist works", {
    x <- "ALRWTC"
    y <- "AMWMT-"
    res <- aminoacid_dist(aa1 = x, aa2 = y)
    expect_equal(as.numeric(res[6]), 2)
})
