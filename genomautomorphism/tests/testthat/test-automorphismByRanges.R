test_that("automorphismByRanges function test", {
    data(autm, package = "GenomAutomorphism")
    test1 <- automorphismByRanges(x = autm[c(24:35)])
    expect_true(test1$cube[1] == "TGCA")
})
