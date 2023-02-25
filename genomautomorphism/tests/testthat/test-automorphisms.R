test_that("automorphisms function test", {
    data(aln)
    autms <- automorphisms(
        seq = aln,
        group = "Z5^3",
        start = 25,
        end = 30,
        verbose = FALSE
    )
    test1 <- autms$autm[2] == "2,1,1"
    autms <- automorphisms(
        seq = aln,
        group = "Z125",
        start = 28,
        end = 30,
        verbose = FALSE
    )
    test2 <- autms$autm == 106
    autms <- automorphisms(
        seq = aln,
        group = "Z64",
        start = 28,
        end = 30,
        verbose = FALSE
    )
    test3 <- autms$autm == 41
    expect_true(all(test1, test2, test3))
})
