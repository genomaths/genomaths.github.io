test_that("multiplication works", {
    data("aaindex2", package = "GenomAutomorphism" )
    s <- get_mutscore("A", "C", acc = "MIYS930101", aaindex = "aaindex2")
    expect_equal(s, -0.51)
})
