test_that("base2codon function test", {
    seq <- c("ACCTCA")
    test1 <- base2codon(x = seq)
    test1 <- all(test1 == c("ACC", "TCA"))
    expect_true(test1)
})
