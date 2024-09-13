set.seed(0)
X = matrix(rnorm(10000), ncol = 100)
control_indices = 1
train_indices = 1:50
test_indices = 51:100
baseline_predictor = colMeans(X[train_indices, ])
correlations = c()
for(i in test_indices){
  correlations[i-50] = cor(baseline_predictor - X[control_indices, ], X[i, ] - X[control_indices, ])
}
mean(correlations)
# 0.649