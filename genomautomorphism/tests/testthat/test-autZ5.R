test_that("autZ5 function test", {
    data(aln)
    autms <- autZ5(
        seq = aln,
        start = 28,
        end = 30,
        verbose = FALSE
    )
    test1 <- autms$autm[1] == 2
    expect_true(test1)
})
