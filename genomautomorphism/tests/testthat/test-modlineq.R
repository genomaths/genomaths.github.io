test_that("modlineq works", {
    x <- c(9,32,24,56,60,27,28,5)
    y <-  c(8,1,0,56,60,0,28,2)
    modulo <- c(64,125,64,64,64,64,64,64)
    mt <- modlineq(a = y, b = x, n = modulo, no.sol = 1L)
    expect_true(all((y %*% diag(mt$diag) + mt$translation) %% modulo == x))
})
