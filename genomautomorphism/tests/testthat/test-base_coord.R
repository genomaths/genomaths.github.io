test_that("base_coord function test", {
    data(aln)
    test1 <- base_coord(
        base = aln,
        cube = "ACGT"
    )
    test1 <- test1$coord1[2] == 1 && test1$coord2[2] == 3
    expect_true(test1)
})
