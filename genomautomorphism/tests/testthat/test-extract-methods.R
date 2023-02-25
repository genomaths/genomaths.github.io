test_that("extract-methods works", {
    data(brca1_autm)
    test1 <- (names(brca1_autm[1]) == "human_1.human_2")
    test2 <- (length(brca1_autm[[3]]) == 761)
    expect_true(test1 && test2)
})
