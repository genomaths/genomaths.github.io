test_that("peptide_phychem_index works", {
    aa <- peptide_phychem_index('ACGTCATCAAGT', acc = "EISD840101")
    expect_equal(sum(aa), 0.9)
})
