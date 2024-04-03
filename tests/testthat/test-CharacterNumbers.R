testthat::test_that("Test CharacterNumbers", {
    ##-----------------------------------------------------------------------##
    testthat::expect_error(CharacterNumbers(Vect.number=c("hello", "test")),
                           "'Vect.number' must be a numeric vector",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    res1_ChNum <- CharacterNumbers(Vect.number=c(0, 1, 2, 10))
    res2_ChNum <- CharacterNumbers(Vect.number=c(0, 1, 2))

    testthat::expect_equal(res1_ChNum, c("00", "01", "02", "10"))
    testthat::expect_equal(res2_ChNum, c("0", "1", "2"))
})
