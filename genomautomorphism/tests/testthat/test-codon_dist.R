test_that("codon_dist works", {
    x <- c("ACG","CGT","GTA","CCG","TGA","CTG","ACG")
    y <- c("TGC","GCC","CGT","GAC","---","TGA","A-G")
    x <- codon_dist(x, y, group = "Z4")
    test1 <- is.na(x[5])
    test2 <- (x[1] == 0.8750)
    expect_true(test1 && test2)
})
