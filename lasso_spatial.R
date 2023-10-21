#sigma_square目标函数
objective2 <- function(sigma_square, rhoo_old, b_old,xx0,yy0,ww) {
  En <- yy0 - rhoo_old * ww %*% yy0 - xx0 %*% b_old;
  Sn <- In - rhoo_old * ww;
  log_likelihood <- -n/2 * log(2*pi*sigma_square) - (t(En) %*% En) / (2*sigma_square) + try(log(det(Sn)))
  return(-log_likelihood)
}

objective3 <- function(sigma_square_new, rhoo, b_old,xx0,yy0,ww) {
  Sn1 <- In - rhoo * ww
  En1 <- yy0 - rhoo * ww %*% yy0 - xx0 %*% b_old
  # 检查det(Sn1)的值是否为正，并且非NA和Inf
  det_Sn1 <- det(Sn1)
  if (!is.na(det_Sn1) && is.finite(det_Sn1) && det_Sn1 > 0) {
    log_likelihood <- -n/2 * log(2 * pi * sigma_square_new) - (t(En1) %*% En1) / (2 * sigma_square_new) + log(det_Sn1)
    return(-log_likelihood)
  } else {
    return(Inf)
  }
}



calculateBIC0_lasso_spatial <- function(lambda_values, xx0, yy0,k,b0,rhoo,ww) {
  n <- nrow(xx0) # 数据的观测数
  BIC_values <- numeric(length(lambda_values)) # 存储不同lambda对应的BIC值
  b_new_values <- list() # 存储不同lambda对应的b_new
  
  for (i in seq_along(lambda_values)) {
    yy <- yy0 - rhoo * ww %*% yy0
    result <- LASSO.LQA(b0,xx0,yy,lambda = lambda_values[i],1e-4,10000,1e-4)
    b_new <- result
    y_hat <- solve(diag(n) - rhoo * ww) %*% xx0 %*% b_new  # 预测值
    #y_hat <- xx0%*%b_new
    SSE <- sum((yy0- y_hat)^2)  # 计算误差平方和
    BIC_values[i] <- -2*(-(n/2)*log(SSE/n)-(n/2)*log(2*pi))+k*log(n)  # 计算BIC值
    b_new_values[[i]] <- b_new  # 存储b_new
  }
  
  min_BIC_index <- which.min(BIC_values)  # 找到最小的BIC值的索引
  min_BIC <- BIC_values[min_BIC_index]  # 最小的BIC值
  min_b_new <- b_new_values[[min_BIC_index]]  # 对应最小BIC值的b_new
  best_lambda <- lambda_values[min_BIC_index]  # 对应最小BIC值的lambda
  
  return(list(b_new = min_b_new, lambda = best_lambda, BIC = min_BIC))
}

#step2~step5
fun_stepall0_lasso_spatial <-  function(sigma_square,rhoo,b,xx0,yy0,ww){
  
  
  iterations_count <- 0
  # 初始化旧变量
  sigma_square_old <- sigma_square
  rhoo_old <- rhoo
  b_old <- b
  sigma_square_new <- sigma_square
  rhoo_new <- rhoo
  b_new<- b
  # 迭代更新步骤2到步骤4，直到收敛
  for (iteration in 1:max_iterations){
    iterations_count <- iterations_count + 1
    
    
    # Step 2: 更新 sigma_square
    lower2 <- 0   # 最小值约束
    upper2 <- 1000  # 最大值约束
    
    result2 <- optimize(function(sigma_square){
      objective2(sigma_square, rhoo_old, b_old, xx0, yy0, ww)
    }, interval = c(lower2, upper2))
    sigma_square_new <- result2$minimum
    
    # Step 3: 更新 rhoo
    library(Matrix)
    lower3 <- -1   # 最小值约束
    upper3 <- 0.99  # 最大值约束
    result3 <- optimize(function(rhoo) {
      objective3(sigma_square_new, rhoo, b_old, xx0, yy0, ww)
    }, interval = c(lower3, upper3))
    rhoo_new <- result3$minimum
    
    
    
    # Step 4: 更新 b
    yy <- yy0- rhoo_new * ww %*% yy0
    #求解带约束的b
    # 定义lambda的范围
    lambda_values <- seq(0.5, 0.8, by = 0.1)
    
    # 调用函数计算最小BIC值和对应的b_new、lambda
    #result <- calculateBIC0_lasso_spatial(lambda_values, xx0, yy0, 3+q, b_old,rhoo_new,ww)
    
    # 提取结果
    
    #best_lambda <- result$lambda
    best_lambda <- 0.08
    
    
    b_new <-LASSO.LQA(b_old, xx0, yy, lambda = best_lambda,1e-4,10000,1e-4)
    
    
    
    
    # 结束条件，判断是否收敛
    if (max(abs(b_new - b_old)) < tolerance &&
        abs(rhoo_new - rhoo_old) < tolerance &&
        abs(sigma_square_new - sigma_square_old) < tolerance) {
      break  # 如果满足收敛条件，跳出迭代循环
    }
    
    # 更新旧变量
    b_old <- b_new
    rhoo_old <- rhoo_new
    sigma_square_old <- sigma_square_new
    #print(paste(iterations_count,":   b:",b_old,", rhoo: ",rhoo_old,", sigma_square: ",sigma_square_old))
    
  }
  
  final_list <- list(sigma_square_new, rhoo_new, b_new,iterations_count)
  return (final_list)
  
}



result0 <- fun_stepall0_lasso_spatial(sigma_square,rhoo,b,xx,yy,Wn)
result0 <- fun_stepall0(sigma_square,rhoo,b0,trans_x0,trans_y0,ww)
b_final <- result0[[3]]
rhoo_final <- result0[[2]]
sigma_square_final <- result0[[1]]
iterations <- result0[[4]]