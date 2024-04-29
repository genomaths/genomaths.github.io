test_that("multiplication works", {
     data("aaindex2", package = "GenomAutomorphism" )
     mat <- aa_mutmat(aaindex = "aaindex2", acc_list = TRUE)
     mat <- grepl("MIYS930101", mat[37])
     aa <- aa_mutmat(acc = "MIYS930101", aaindex = "aaindex2")
     aa <- sum(aa[1,]) == -2.96
     expect_true(mat && aa)
})
