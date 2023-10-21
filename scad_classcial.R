# First order derivate of penalized function
dp <- function(beta, lambda, hyper){
  ifelse(abs(beta) <= lambda, lambda*sign(beta), lambda*sign(beta)* max(hyper*lambda - beta, 0)/(hyper-1)*lambda)
}

LASSO.LQA3 <- function(beta0, x, y, lambda, thd1, M, thd2){
  # 用LQA算法计算lasso的数值解
  # beta0为初始系数向量（p维），n为样本量，p为变量维数，lambda为调节参数的待选范围，thd1为将β分量压缩为0的阈值，thd2为停止迭代的阈值
  # 在迭代过程中将β小于阈值的分量压缩为零
  n <- dim(x)[1] 
  p <- dim(x)[2]
  # 重写sigma的计算函数
  sigma.func <- function(beta, lambda, a) {
    p <- length(beta)
    sigma <- matrix(0, nrow = p, ncol = p)
    
    di <- c()
    for(i in 1:p){
      di[i] <- dp(abs(beta[i]), lambda, a)/abs(beta[i])
    }
    diag(sigma) <- di
    return(sigma)
  }
  
  beta <- beta0                                   # 初值为最小二乘估计
  k = 0
  repeat{
    sigma <- sigma.func(beta,lambda, 3.7) 
    newbeta <- ginv( t(x)%*%x + n*sigma ) %*% t(x) %*% y # 计算beta_k+1
    for(beta.i in newbeta){                                # 当其分量足够小时将其压缩为0
      if(abs(beta.i) < thd1) beta.i <- 0
    }
    if(t(newbeta-beta)%*%(newbeta-beta) < thd2 ||  (k >= 5000)) break      # 当两次迭代的向量差距(差的二范数)足够小时停止迭代
    beta <- newbeta
    k = k + 1
    print(paste("循环：",k))
    #cat("\r", k)                                           # 输出迭代次数
  }
  return(beta) 
}




calculateBIC0_scad_classical <- function(lambda_values, xx0, yy0,k,b0,rhoo,ww) {
  n <- nrow(xx0) # 数据的观测数
  BIC_values <- numeric(length(lambda_values)) # 存储不同lambda对应的BIC值
  b_new_values <- list() # 存储不同lambda对应的b_new
  
  for (i in seq_along(lambda_values)) {
    #yy <- yy0 - rhoo * ww %*% yy0
    result <- LASSO.LQA3(b0,xx0,yy0,lambda = lambda_values[i],1e-4,10000,1e-4)
    b_new <- result
    #y_hat <- solve(diag(n) - rhoo * ww) %*% xx0 %*% b_new  # 预测值
    y_hat <- xx0%*%b_new
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


b0 <- ginv(t(xx)%*%xx)%*%t(xx)%*%yy
result <- LASSO.LQA3(b0,xx,yy,0.1,1e-4,10000,1e-4)



# 开始蒙特卡洛模拟循环
num_simulations <- 100;
t <- 0;
b_zero_count_correct <- 0;
b_zero_count_incorrect <- 0;
mse_sum <- 0;
b1_final_total <- rep(0,100);
b2_final_total <- rep(0,100);
b3_final_total <- rep(0,100);
sigma_square_total <- rep(0,100);

for (i in 1:num_simulations) {
  set.seed(i)
  q <- 5  # 用于表示维度减去3的值
  n <- 60
  rhoo <- 0.8
  R <- 20
  # 计算协方差矩阵
  cov_mat <- matrix(NA, nrow = q + 3, ncol = q + 3)
  for (i in 1:(q + 3)) {
    for (j in 1:(q + 3)) {
      cov_mat[i, j] <- 0.5^abs(i - j)
    }
  }
  # 生成随机样本，其中协变量考虑了协方差矩阵
  sample <- MASS::mvrnorm(n , mu = rep(0, q + 3), Sigma = cov_mat)
  
  xx <- sample
  #生成空间矩阵
  vec1m <- c(1, 1, 1)
  In3 <- diag(3)
  Bm <- (1/(3-1))*(vec1m %*% t(vec1m)-In3)
  InR <- diag(R)
  Wn <- InR %x% Bm
  #b取值
  zeroq <- rep(0,q)
  bb <- c(3,2,1.6,zeroq)
  #误差项
  e0 = rnorm(n, 0, 1.5)
  e <- e0;
  
  #生成y
  yy <- solve(diag(n) - rhoo * Wn) %*% (xx %*% bb +e);
  #约束项
  EE <- matrix(c(0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0), nrow = 2, ncol = 8, byrow = TRUE)
  ff <- c(1.6,3)
  CC <- matrix(c(1,0,1,0,0,0,0,0,0,-1,0,0,0,-1,0,0),nrow = 2, ncol = 8, byrow = TRUE)
  dd <- c(4,-2.5)
  #sigma_square
  sigma_square_0 = var(e);
  # 设置收敛的阈值和最大迭代次数
  tolerance <- 1e-6  # 收敛阈值
  max_iterations <- 100  # 最大迭代次数
  #获取BIC对应值
  bic_values <- c()
  # 定义lambda的范围
  lambda_values <- seq(0.5, 1.0, by = 0.1)
  bb0 <-  ginv(t(xx)%*%xx)%*%t(xx)%*%yy
  # 调用函数计算最小BIC值和对应的b_new、lambda
 # result <- calculateBIC0_scad_classical(lambda_values, xx, yy, 3+q, bb0,rhoo,ww)
  
  # 提取结果
  
  #best_lambda <- result$lambda
  best_lambda <- 0.15
  
  
  b_new <-LASSO.LQA3(bb0, xx, yy, lambda = best_lambda,1e-4,10000,1e-4)
  b_final <- b_new
  
  
  t=t+1;
  b1_final_total[t] <- b_final[1]
  b2_final_total[t] <- b_final[2]
  b3_final_total[t] <- b_final[3]
  
  b0_count_correct <- sum(b_final[4:length(b_final)] < 0.05)
  b0_count_incorrect <- sum(b_final[1:3]<0.005)
  b_zero_count_correct <- b_zero_count_correct+b0_count_correct
  b_zero_count_incorrect <- b_zero_count_incorrect+b0_count_incorrect
  mse <- sum(abs(b_final - bb))
  mse_sum <- mse_sum +mse
}

correct <-  b_zero_count_correct/100
incorrect <- b_zero_count_incorrect/100
mse <- mse_sum/100
#b1
b1_mean <- mean(b1_final_total)
b1_mad <- median( abs(b1_final_total - median((b1_final_total))))
b1_sd <- sd(b1_final_total)
#b2
b2_mean <- mean(b2_final_total)
b2_mad <- median( abs(b2_final_total - median((b2_final_total))))
b2_sd <- sd(b2_final_total)
#b3
b3_mean <- mean(b3_final_total)
b3_mad <- median( abs(b3_final_total - median((b3_final_total))))
b3_sd <- sd(b3_final_total)
#sigma
#sigma_mad <- median( abs(sqrt(sigma_square_total) - median((sqrt(sigma_square_total)))))

