test_that("Test CharacterNumbers", {
    #
    expect_error(CharacterNumbers(Vect.number=c("hello", "test")),
                 "'Vect.number' must be a numeric vector",
                 fixed=TRUE)
    #
    expect_equal(CharacterNumbers(Vect.number=c(0, 1, 2, 10)),
                 c("00", "01", "02", "10"))
})
#
