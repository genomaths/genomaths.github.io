test_that("slapply works", {
    x <- list(a = 1:10, beta = exp(-3:3), logic = c(TRUE, FALSE, FALSE, TRUE))
    class(x) <- "nice"
    expect_equal(
        class(slapply(
            x,
            mean,
            keep.attr = TRUE,
            simplify = TRUE
        )),
        "nice"
    )
})
