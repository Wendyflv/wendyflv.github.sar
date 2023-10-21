#带迭代算法的空间模型LASSO对Boston房价估计
library(spdep)

# 加载Boston数据集
data(Boston)
summary(Boston)
# Create a 506x506 matrix filled with zeros
matrix <- matrix(0, nrow = 506, ncol = 506)

# Assign values based on the condition
for (i in 1:506) {
  for (j in 1:506) {
    if (abs(i - j) == 1) {
      matrix[i, j] <- 1
    }
  }
}
mat_normalized <- lapply(1:nrow(matrix), function(i) matrix[i, ] / sum(matrix[i, ]))
mat_normalized <- do.call(rbind, mat_normalized)
ww <- mat_normalized
trans_y0 <- log(yy0)
trans_x0 <- cbind(xx0[,1],xx0[,2],xx0[,3],xx0[,4],xx0[,5]^2,xx0[,6]^2,xx0[,7],log(xx0[,8]),log(xx0[,9]),xx0[,10],xx0[,11],xx0[,12],log(xx0[,13]))
trans_x0[is.infinite(trans_x0)] <- 0



calculateBIC2 <- function(lambda_values, xx0, yy0,k,b0) {
  n <- nrow(xx0) # 数据的观测数
  BIC_values <- numeric(length(lambda_values)) # 存储不同lambda对应的BIC值
  b_new_values <- list() # 存储不同lambda对应的b_new
  
  for (i in seq_along(lambda_values)) {
    result <- LASSO.LQA(b0,xx0,yy0,lambda = lambda_values[i],1e-4,10000,1e-4)
    b_new <- result
    #fit <- flasso.fit.cvxr1(xx0, yy0, lambda = lambda_values[i], rho = rhoo, w = ww, Cmat1 = CC1, dvec1 = dd1,Cmat2 = CC2, dvec2 = dd2,Cmat3 = CC3, dvec3 = dd3,Cmat4 = CC4, dvec4 = dd4,Cmat5 = CC5, dvec5 = dd5, Emat = EE, fvec = ff)  # 进行回归拟合
    #b_new <- unlist(fit[12])  # 提取回归系数
    #y_hat <- solve(diag(n) - rhoo * ww) %*% xx0 %*% b_new  # 预测值
    y_hat <- xx0%*%b_new
    SSE <- sum((yy0 - y_hat)^2)  # 计算误差平方和
    BIC_values[i] <- -2*(-(n/2)*log(SSE/n)-(n/2)*log(2*pi))+k*log(n)  # 计算BIC值
    b_new_values[[i]] <- b_new  # 存储b_new
  }
  
  min_BIC_index <- which.min(BIC_values)  # 找到最小的BIC值的索引
  min_BIC <- BIC_values[min_BIC_index]  # 最小的BIC值
  min_b_new <- b_new_values[[min_BIC_index]]  # 对应最小BIC值的b_new
  best_lambda <- lambda_values[min_BIC_index]  # 对应最小BIC值的lambda
  
  return(list(b_new = min_b_new, lambda = best_lambda, BIC = min_BIC))
}

# 定义lambda的范围
lambda_values <- seq(0.01, 0.1, by = 0.01)

# 调用函数计算最小BIC值和对应的b_new、lambda
result <- calculateBIC2(lambda_values, trans_x0, trans_y0,  13, b0)

# 提取结果
min_b_new <- result$b_new
best_lambda <- result$lambda
min_BIC <- result$BIC