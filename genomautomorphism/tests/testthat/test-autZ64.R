test_that("autZ64 function test", {
    data(aln)
    autms <- autZ64(
        seq = aln,
        start = 28,
        end = 30,
        verbose = FALSE
    )
    test1 <- autms$autm == 41
    expect_true(test1)
})
