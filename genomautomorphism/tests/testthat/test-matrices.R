test_that("matrices function test", {
    data(aln, package = "GenomAutomorphism")
    test1 <- matrices(
        x = aln,
        base_seq = TRUE,
        filepath = NULL,
        cube = "ACGT",
        group = "Z4",
    )

    test2 <- matrices(
        x = aln,
        base_seq = FALSE,
        filepath = NULL,
        cube = "ACGT",
        group = "Z5^3",
    )

    test1 <- sum(test1$coord1[1:10]) == 18
    test2 <- sum(test2$coord1[1:10]) == 21

    expect_true(all(test1, test2))
})
