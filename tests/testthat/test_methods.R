
context("Test methods")

library(DESeq2)

test_that("Test getsf", {
    n = 100
    m = 50
    gsd = makeExample(n = n, m = m)
    sf = getsf(gsd)
    expect_true(is.null(sf))
    gsd = NBAMSeq(gsd, parallel = TRUE)
    sf = getsf(gsd)
    expect_true(is.numeric(sf))
    expect_true(length(sf) == m)
    expect_true(all(names(sf) == paste0("sample", 1:m)))

    dds = DESeqDataSetFromMatrix(countData = assay(gsd),
                                 colData = colData(gsd), design = ~pheno)
    dds = estimateSizeFactors(dds)
    expect_true(all(sizeFactors(dds)==sf))

})


test_that("set getsf", {
    n = 100
    m = 50
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


    sf = sample(1:5,m, replace = TRUE)
    setsf(gsd) = sf
    expect_true(all(getsf(gsd) == sf))
    gsd = NBAMSeq(gsd, parallel = TRUE)
    expect_true(all(getsf(gsd) == sf))

})



