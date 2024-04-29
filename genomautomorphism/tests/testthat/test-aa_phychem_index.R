test_that("aa_phychem_index works", {
    data("aaindex1", package = "GenomAutomorphism" )
    m <- mean(aa_phychem_index(acc = "ARGP820101", aaindex = "aaindex1"))
    expect_equal(round(m,4), 0.9975)
})
