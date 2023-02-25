test_that("get_coord function test", {
    data(aln, package = "GenomAutomorphism")
    test1 <- get_coord(
        x = aln,
        cube = "ACGT",
        group = "Z5"
    )
    test1 <- sum(test1@CoordList$coord1) == 127
    test2 <- get_coord(
        x = aln,
        base_seq = FALSE,
        cube = "ACGT",
        group = "Z64"
    )
    test2 <- sum(test2@CoordList$coord1[2:5]) == 168
    test3 <- get_coord(
        x = aln,
        base_seq = FALSE,
        cube = "ACGT",
        group = "Z125"
    )
    test3 <- sum(test3@CoordList$coord1[1:5]) == 428
    expect_true(all(test1, test2, test3))
})
