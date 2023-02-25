test_that("aut3D function test", {
    data(aln, package = "GenomAutomorphism")
    autms <- aut3D(
        seq = aln,
        start = 25,
        end = 30,
        verbose = FALSE
    )
    expect_true(autms$autm[2] == "2,1,1" && autms$aa1[1] == "Q")
})
