test_that("autZ125 function test", {
    data(aln)
    autms <- autZ125(
        seq = aln,
        start = 28,
        end = 30,
        verbose = FALSE
    )
    test1 <- autms$autm == 106
    expect_true(test1)
})
