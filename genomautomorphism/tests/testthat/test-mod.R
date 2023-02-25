test_that("mode works", {
    n <- diag(x=1, nrow = 4, ncol = 4) * c(43,125,2,112)
    m <- c(64,4,4,64)
    m <- mod(n = n, m = m)
    expect_true(all(diag(m) == c(43,1,2,48)))
})
