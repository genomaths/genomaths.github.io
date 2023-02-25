test_that("AutomorphismList-methods works", {
    data(brca1_autm)
    test1 <- names(brca1_autm[1]) == "human_1.human_2"
    test2 <- inherits(brca1_autm[[1]], "Automorphism")
    test3 <- inherits(brca1_autm$human_1.gorilla_1, "Automorphism")
    expect_true(all(test1,test2,test3))
})
