test_that("automorphism_bycoef function test", {
    data(autm, package = "GenomAutomorphism")
    test1 <- automorphism_bycoef(x = autm[1:10])
    expect_true(test1$autm[9] == 1 && test1$aa2[10] == "E")
})
