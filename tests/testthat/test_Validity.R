context("Test validity")

test_that("Test invalid input", {

    n = 200
    m = 30
    counts = matrix(rnbinom(n * m, 10, 0.3), n)

    var1 = rnorm(m)
    var2 = rnorm(m)
    var3 = rnorm(m)
    var4 = as.factor(sample(c(0,1,2), m, replace = TRUE))
    colData = data.frame(var1 = var1, var2 = var2,
                         var3 = var3, var4 = var4)

    counts1 = counts
    counts1[1,1] = "1"
    expect_error(NBAMSeqDataSet(counts1, colData, design = ~ s(var1)))
    counts2 = counts
    counts2[1,1] = NA
    expect_error(NBAMSeqDataSet(counts2, colData, design = ~ s(var1)))
    counts3 = counts
    counts3[1,1] = 3.5
    expect_message(NBAMSeqDataSet(counts3, colData, design = ~ s(var1)),
                   "countData contains non-integers, rounded to the nearest
                integer automatically.")
    counts4 = counts
    counts4[1,1] = -5
    expect_error(NBAMSeqDataSet(counts4, colData, design = ~ s(var1)))

    design1 = matrix(sample(m*4), m, 4)
    expect_error(NBAMSeqDataSet(counts, colData, design = design1))
    expect_error(NBAMSeqDataSet(countData = counts, colData = colData,
        design = ~ var4 + s(var1) + var2 + var3 + var5))

    colData2 = colData
    colData2$var1 = sample(c("normal", "tumor"), m, replace = TRUE)
    expect_error(NBAMSeqDataSet(countData = counts, colData = colData2,
                                design = ~ var4 + var1 + s(var2) + var3))
    expect_error(NBAMSeqDataSet(countData = counts, colData = colData,
        design = ~ var4 + var1 + var2 + var3))
    expect_error(NBAMSeqDataSet(countData = counts, colData = colData,
                                design = ~ s(var4) + var1 + var2 + var3),
                 "Non linear term can not be a factor.")

    counts5 = counts
    counts5[1:2,] = 0
    expect_message(NBAMSeqDataSet(counts5, colData, design = ~ s(var1)),
"2 genes have all 0 counts for all samples. Consider filtering out these genes
before differential expression analysis.")
    counts6 = matrix(0,n,m)
    expect_error(NBAMSeqDataSet(counts6, colData, design = ~ s(var1)))
})

