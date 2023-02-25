test_that("codon_dist_matrix works", {
    x <- codon_dist_matrix(genetic_code = "5", cube = "TGCA", group = "Z4",
                        output = "vector")
    expect_equal(x[[63]], 0.375)
})
