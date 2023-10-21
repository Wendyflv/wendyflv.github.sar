# rhoo <- 0.2 sigma1  SCAD  spatial
#n <- 60 q <- 5


LASSO.LQA <- function(beta0, x, y, lambda, thd1, M, thd2){
  # 用LQA算法计算lasso的数值解
  # beta0为初始系数向量（p维），n为样本量，p为变量维数，lambda为调节参数的待选范围，thd1为将β分量压缩为0的阈值，thd2为停止迭代的阈值
  # 在迭代过程中将β小于阈值的分量压缩为零
  n <- dim(x)[1] 
  p <- dim(x)[2]
  sigma.func <- function(beta){                   # Σ矩阵
    beta.copy = beta
    for(i in 1:p){
      if(beta.copy[i] == 0) beta.copy[i] <- M     # 为使矩阵逆存在，将1/0用一个很大的数M代替
      else beta.copy[i] <- 1/abs(beta.copy[i])
    }
    sigma <- diag(as.vector(lambda*beta.copy))    # 对角线为lambda_i / |beta|
    return(sigma)
  }
  
  beta <- beta0                                   # 初值为最小二乘估计
  k = 0
  repeat{
    sigma <- sigma.func(beta) 
    newbeta <- ginv( t(x)%*%x + n*sigma ) %*% t(x) %*% y # 计算beta_k+1
    for(beta.i in newbeta){                                # 当其分量足够小时将其压缩为0
      if(abs(beta.i) < thd1) beta.i <- 0
    }
    if(t(newbeta-beta)%*%(newbeta-beta) < thd2  ||  (k >= 10000)) break      # 当两次迭代的向量差距(差的二范数)足够小时停止迭代
    beta <- newbeta
    k = k + 1
    print(paste("循环：",k))
    #cat("\r", k)                                           # 输出迭代次数
  }
  return(beta) 
}





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
  q <- 50 # 用于表示维度减去3的值
  n <- 120
  rhoo <- 0.0
  R <- 40
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
  e0 = rnorm(n, 0, 1.0)
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
  
  #写似然函数
  In <- diag(n);
  Sn <- In - rhoo * Wn;
  En <- yy - rhoo * Wn %*% yy - xx %*% bb;
  log_likelihood <- -n/2 * log(2*pi*sigma_square_0) - (t(En) %*% En) / (2*sigma_square_0)+log(det(Sn));
  
  
  #step 1 初始化
  #对给定rhoo
  b = solve(t(xx)%*%xx) %*% t(xx) %*% Sn %*% yy
  sigma_square = (t(En) %*% En) / n
  # 初始化变量
  b_old <- b
  rhoo_old <- rhoo
  sigma_square_old <- sigma_square
  
  #获取BIC对应值
  bic_values <- c()
  result <- fun_stepall0_lasso_spatial(sigma_square_old,rhoo_old,b_old, xx, yy, Wn)
  
  #b,rhoo,sigma_square最终估计
  b_final <- result[[3]]
  t=t+1;
  b1_final_total[t] <- b_final[1]
  b2_final_total[t] <- b_final[2]
  b3_final_total[t] <- b_final[3]
  rhoo_final <- result[[2]]
  sigma_square_final <- result[[1]]
  sigma_square_total[t] <- sigma_square_final
  b0_count_correct <- sum(b_final[4:length(b_final)] < 0.01)
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
sigma_mad <- median( abs(sqrt(sigma_square_total) - median((sqrt(sigma_square_total)))))

