
context("Test methods")

test_that("set setsf", {
    n = 5
    m = 10
    gsd = makeExample(n = n, m = m)
    sf = as.factor(sample(1:5,m, replace = TRUE))
    expect_error((setsf(gsd) = sf))
    sf = sample(1:5,m, replace = TRUE)
    sf[5] = NA
    expect_error((setsf(gsd) = sf))
    sf = sample(1:5,m-1, replace = TRUE)
    expect_error((setsf(gsd) = sf),
                 "Size factor length should be the same as number of samples.")
    sf = sample(1:5,m, replace = TRUE)
    sf[5] = -2
    expect_error((setsf(gsd) = sf), "Size factor cannot be 0 or negative.")

})



