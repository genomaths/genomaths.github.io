test_that("seqranges works", {
    data(aln)
    sq <- seqranges(
        x = aln,
        base_seq = TRUE,
        filepath = NULL,
    )
    sq <- paste(sq$seq1[seq(5)], collapse = "")
    expect_equal(sq, "ACCTA")
})
