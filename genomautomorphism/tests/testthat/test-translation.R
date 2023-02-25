test_that("multiplication works", {
    data("aln", package = "GenomAutomorphism")
    aa <- translation(aln)
    expect_true(grepl("IYVGITTLQFLGS", as.character(aa[2])))
})
