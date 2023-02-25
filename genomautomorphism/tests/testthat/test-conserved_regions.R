test_that("conserved_regions function test", {
    data(autm, package = "GenomAutomorphism")
    test1 <- conserved_regions(autm[24:30])
    expect_true(test1$aa1[1] == "D" && test1$seq2[4] == "GGC")
})
