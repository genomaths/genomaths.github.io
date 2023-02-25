test_that("as.AutomorphismList works", {
    data(autm)
    lista <- list(a1 = autm, a2 = autm, a3 = autm, a4 = autm)
    test1 <- inherits(as.AutomorphismList(lista), "AutomorphismList")
    lista <- as(lista, "GRangesList")
    test2 <- inherits(as.AutomorphismList(lista), "AutomorphismList")
    expect_true(test1,test2)
})
