require(robustlmm)

testData <- generateAnovaDatasets(1, c(2, 5), c(3, 5), 10)$generateData(1)
rlmer(y ~ (1 + Var1|Var2) + (1 + Var4|Var3), testData, max.iter = 1)
## test is that this runs without errors
