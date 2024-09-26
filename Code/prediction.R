library(pfica)
library(fda)
library(plot.matrix)
library(forecast)
library(wavelets)
library(caret)
library(rlist)
library(pls)
library(foreach)
library(MFPCA)
library(doParallel)
library(diagram)
library(MLmetrics)
library(data.table) # for copy
library(ggplot2)
library(ips)
library(RcppAlgos)
library(writexl)
library(Matrix)
library(expm)
library(parallel) # faster computing
library(assist)   # find smoother matrix
library(wavethresh) # wavelet basis
library(fda)
library(far) # for Gram-schmidt, but too slow!
library(quadprog)
library(Matrix)
library(fda.usc)
library(splines)
library(grplasso)
library(glmnet)
library(ggpmisc)
library(refund)

source("1.function_define.R")


#################################
#### Scenario I (stationary) ####
#################################
load("../Data/stationary_true.RData")
result.df.total.stationary.list <- list()
for(ttt in 1:20){
  print(paste("TOTAL iteration:", ttt))
  stationary_simuldata.list.tmp <- list()
  x_err_std <- 3
  set.seed(10)
  for(i in 1:300){
    # stationary_simuldata.list[[i]] <- t(t(stationary_true)+rnorm(ncol(stationary_true), 0, x_err_std))
    stationary_simuldata.list.tmp[[i]] <- stationary_true + matrix(rnorm(5*1024, 0, x_err_std), nrow=5, ncol=1024)
  }
  num_datasets <- 300
  X_total_stationary.tmp <- array(NA, dim = c(num_datasets, 1024, 5))
  
  # Generate 300 different datasets
  for (i in 1:num_datasets) {
    X_total_stationary.tmp[i, , ] <- t(stationary_simuldata.list.tmp[[i]])
  }
  
  X_total_stationary_time_wav.tmp <- array(NA, dim = c(num_datasets, 5, 1024))
  for (i in 1:num_datasets) {
    X_total_stationary_time_wav.tmp[i, , ] <- stationary_simuldata.list.tmp[[i]]
  }
  
  num_datasets <- dim(X_total_stationary.tmp)[1]
  Y_stationary <- numeric(num_datasets)
  Z_stationary.tmp <- numeric(num_datasets)
  # Create time and wavelength axes with a small random noise added
  time_axis <- sin(2*pi/5*c(c(1:5)-1/3))
  wavelength_axis <- sin(2*pi/256*c(c(1:1024)-1/3))
  # Compute the outer product for the i-th dataset
  beta <- outer(time_axis, wavelength_axis, FUN = "*")
  
  for (i in 1:num_datasets) {
    Z_stationary.tmp[i] <- sum((X_total_stationary.tmp[i,,]-mean(X_total_stationary.tmp[i,,])) * t(beta))
  }
  
  SNR_list <- c(2,5,10,20)
  y_err_std_list <- sqrt(var(Z_stationary.tmp)/SNR_list)
  result.df.total.stationary.tmp <- data.frame(row.names = paste(c("RMSE", "MAE"), "(SNR=", rep(SNR_list, each=2), ")", sep=""))
  
  ################################################################################
  # Two-way model averaging
  ################################################################################
  ncomp = 5 # SVD (u,v) pair 2개로 가정
  ## PenFSVD ##
  penfsvd = apply(X_total_stationary_time_wav.tmp, 1, tw_svd_smu, ncomp=ncomp)
  
  ## basisFSVD ##
  basfsvd = apply(X_total_stationary_time_wav.tmp, 1, basis_svd, nb_by = 1, ncmp=ncomp)
  
  nt <- dim(X_total_stationary_time_wav.tmp)[2]
  nlambda <- dim(X_total_stationary_time_wav.tmp)[3]
  
  ## sign adjustment (flip) ##
  pen_u1 = matrix(0, nrow=num_datasets, ncol=nt)
  for (i in 1:num_datasets) {
    pen_u1[i,] = penfsvd[[i]]$pen_u[1:nt,1]
  }
  for(i in which(pen_u1[,1]<0)){
    penfsvd[[i]]$pen_u[1:nt,1] <- - penfsvd[[i]]$pen_u[1:nt,1]
  }
  # matplot(t(pen_u1), type="l") # sign이 다르게 나오는 wafer가 있으면 조정해주기
  pen_u2 = matrix(0, nrow=num_datasets, ncol=nt)
  for (i in 1:num_datasets) {
    pen_u2[i,] = penfsvd[[i]]$pen_u[1:nt,2]
  }
  # matplot(t(pen_u2), type="l")
  for(i in which(pen_u2[,1]<0)){
    penfsvd[[i]]$pen_u[1:nt,2] <- - penfsvd[[i]]$pen_u[1:nt,2]
  }
  
  pen_u3 = matrix(0, nrow=num_datasets, ncol=nt)
  for (i in 1:num_datasets) {
    pen_u3[i,] = penfsvd[[i]]$pen_u[1:nt,3]
  }
  # matplot(t(pen_u3), type="l")
  for(i in which(pen_u3[,1]<0)){
    penfsvd[[i]]$pen_u[1:nt,3] <- - penfsvd[[i]]$pen_u[1:nt,3]
  }
  
  pen_u4 = matrix(0, nrow=num_datasets, ncol=nt)
  for (i in 1:num_datasets) {
    pen_u4[i,] = penfsvd[[i]]$pen_u[1:nt,4]
  }
  # matplot(t(pen_u4), type="l")
  for(i in which(pen_u4[,2]<0)){
    penfsvd[[i]]$pen_u[1:nt,4] <- - penfsvd[[i]]$pen_u[1:nt,4]
  }
  
  pen_u5 = matrix(0, nrow=num_datasets, ncol=nt)
  for (i in 1:num_datasets) {
    pen_u5[i,] = penfsvd[[i]]$pen_u[1:nt,5]
  }
  # matplot(t(pen_u5), type="l")
  for(i in which(pen_u5[,1]<0)){
    penfsvd[[i]]$pen_u[1:nt,5] <- - penfsvd[[i]]$pen_u[1:nt,5]
  }
  
  pen_v1 = matrix(0, nrow=num_datasets, ncol=nlambda)
  for (i in 1:num_datasets) {
    pen_v1[i,] = penfsvd[[i]]$pen_v[,1]
  }
  for(i in which(pen_v1[,1]<0)){
    penfsvd[[i]]$pen_v[1:nlambda,1] <- - penfsvd[[i]]$pen_v[1:nlambda,1]
  }
  # matplot(t(pen_v1), type="l") # sign이 다르게 나오는 wafer가 있으면 조정해주기
  pen_v2 = matrix(0, nrow=num_datasets, ncol=nlambda)
  for (i in 1:num_datasets) {
    pen_v2[i,] = penfsvd[[i]]$pen_v[,2]
  }
  # matplot(t(pen_v2), type="l")
  for(i in which(pen_v2[,1]<0)){
    penfsvd[[i]]$pen_v[1:nlambda,2] <- - penfsvd[[i]]$pen_v[1:nlambda,2]
  }
  
  pen_v3 = matrix(0, nrow=num_datasets, ncol=nlambda)
  for (i in 1:num_datasets) {
    pen_v3[i,] = penfsvd[[i]]$pen_v[,3]
  }
  # matplot(t(pen_v3), type="l")
  for(i in which(pen_v3[,1]<0)){
    penfsvd[[i]]$pen_v[1:nlambda,3] <- - penfsvd[[i]]$pen_v[1:nlambda,3]
  }
  # for(i in which(pen_v3[,15]<0)){
  #   penfsvd[[i]]$pen_v[,3] <- - penfsvd[[i]]$pen_v[,3]
  # }
  
  pen_v4 = matrix(0, nrow=num_datasets, ncol=nlambda)
  for (i in 1:num_datasets) {
    pen_v4[i,] = penfsvd[[i]]$pen_v[,4]
  }
  # matplot(t(pen_v4), type="l")
  for(i in which(pen_v4[,1]<0)){
    penfsvd[[i]]$pen_v[1:nlambda,4] <- - penfsvd[[i]]$pen_v[1:nlambda,4]
  }
  
  pen_v5 = matrix(0, nrow=num_datasets, ncol=nlambda)
  for (i in 1:num_datasets) {
    pen_v5[i,] = penfsvd[[i]]$pen_v[,5]
  }
  # matplot(t(pen_v5), type="l")
  for(i in which(pen_v5[,1]<0)){
    penfsvd[[i]]$pen_v[1:nlambda,5] <- - penfsvd[[i]]$pen_v[1:nlambda,5]
  }
  
  # for(i in which(pen_v5[,1]<0)){
  #   penfsvd[[i]]$pen_v[,5] <- - penfsvd[[i]]$pen_v[,5]
  # }
  
  # basis fsvd 결과도 필요하면 sign 조정
  
  ### join the result ###
  pen_u1 = do.call(rbind, lapply(penfsvd, function(x) x$pen_u[,1]))
  pen_u2 = do.call(rbind, lapply(penfsvd, function(x) x$pen_u[,2]))
  pen_v1 = do.call(rbind, lapply(penfsvd, function(x) x$pen_v[,1]))
  pen_v2 = do.call(rbind, lapply(penfsvd, function(x) x$pen_v[,2]))
  
  pen_u3 = do.call(rbind, lapply(penfsvd, function(x) x$pen_u[,3]))
  pen_u4 = do.call(rbind, lapply(penfsvd, function(x) x$pen_u[,4]))
  pen_v3 = do.call(rbind, lapply(penfsvd, function(x) x$pen_v[,3]))
  pen_v4 = do.call(rbind, lapply(penfsvd, function(x) x$pen_v[,4]))
  
  pen_u5 = do.call(rbind, lapply(penfsvd, function(x) x$pen_u[,5]))
  pen_v5 = do.call(rbind, lapply(penfsvd, function(x) x$pen_v[,5]))
  
  
  oes_penfsvd <- list(u = vector("list", length = ncomp),
                      v = vector("list", length = ncomp))
  oes_penfsvd$u[[1]] = pen_u1
  oes_penfsvd$u[[2]] = pen_u2
  oes_penfsvd$v[[1]] = pen_v1
  oes_penfsvd$v[[2]] = pen_v2
  
  oes_penfsvd$u[[3]] = pen_u3
  oes_penfsvd$u[[4]] = pen_u4
  oes_penfsvd$v[[3]] = pen_v3
  oes_penfsvd$v[[4]] = pen_v4
  
  oes_penfsvd$u[[5]] = pen_u5
  oes_penfsvd$v[[5]] = pen_v5
  
  bas_u1 = do.call(rbind, lapply(basfsvd, function(x) x$u[,1]))
  bas_u2 = do.call(rbind, lapply(basfsvd, function(x) x$u[,2]))
  bas_v1 = do.call(rbind, lapply(basfsvd, function(x) x$v[,1]))
  bas_v2 = do.call(rbind, lapply(basfsvd, function(x) x$v[,2]))
  
  bas_u3 = do.call(rbind, lapply(basfsvd, function(x) x$u[,3]))
  bas_u4 = do.call(rbind, lapply(basfsvd, function(x) x$u[,4]))
  bas_v3 = do.call(rbind, lapply(basfsvd, function(x) x$v[,3]))
  bas_v4 = do.call(rbind, lapply(basfsvd, function(x) x$v[,4]))
  
  bas_u5 = do.call(rbind, lapply(basfsvd, function(x) x$u[,5]))
  bas_v5 = do.call(rbind, lapply(basfsvd, function(x) x$v[,5]))
  
  oes_basisfsvd <- list(u = vector("list", length = ncomp),
                        v = vector("list", length = ncomp))
  oes_basisfsvd$u[[1]] = bas_u1
  oes_basisfsvd$u[[2]] = bas_u2
  oes_basisfsvd$v[[1]] = bas_v1
  oes_basisfsvd$v[[2]] = bas_v2
  
  oes_basisfsvd$u[[3]] = bas_u3
  oes_basisfsvd$u[[4]] = bas_u4
  oes_basisfsvd$v[[3]] = bas_v3
  oes_basisfsvd$v[[4]] = bas_v4
  
  oes_basisfsvd$u[[5]] = bas_u5
  oes_basisfsvd$v[[5]] = bas_v5
  
  for(ii in 1:length(y_err_std_list)){
    print(paste("iteration:", ii))
    y_err_std <- y_err_std_list[ii]
    set.seed(10)
    for (i in 1:num_datasets) {
      Y_stationary[i] <- sum((X_total_stationary.tmp[i,,]-mean(X_total_stationary.tmp[i,,])) * t(beta)) + rnorm(1, mean = 0, sd = y_err_std)  # adding some noise
    }
    
    ################################################################################
    # MFPLS
    ################################################################################
    set.seed(ttt)
    train_index = sample(1:num_datasets, num_datasets*0.8) 
    test_index = (1:num_datasets)[-train_index]
    
    y_train <- Y_stationary[train_index]
    y_test <- Y_stationary[test_index]
    
    X_train <- X_total_stationary.tmp[train_index,,]
    X_test <- X_total_stationary.tmp[test_index,,]
    
    
    basis_num <- seq(30,60,15) ; basis_num
    Basis_name <- c('Basis 30','Basis 45','Basis 60')
    parameter.change.result.df <- data.frame(row.names = c('pls_num','mean_RMSE'))
    num_folds <- 10 # 10은 시간이 너무 많이 걸림 2로했을 때 약 28분 소요
    n <- dim(X_train)[1] 
    fold_indices <- createFolds(1:n, k = num_folds) 
    
    parameter.change.result.df_stationary <- parameter_selection_MFPLS(X_train,y_train,basis_num,Basis_name,num_folds,fold_indices)
    i_th <- which.min(parameter.change.result.df_stationary[2,]) # select argmin RMSE
    basis <- basis_num[i_th]
    pls_component <- parameter.change.result.df_stationary[1,i_th]
    
    # see test errors
    MFPLS_model <- MFPLS(conv2fda(X_train),y_train, ncomp=pls_component,K=basis)
    prd <- predict.MFPLS(MFPLS_model, conv2fda(X_test))
    y_pred <- prd$y_pred
    # calculate accuracy
    result.df.total.stationary.tmp[2*(ii-1)+1,'proposed'] <- rmse_score(y_pred,y_test)
    result.df.total.stationary.tmp[2*(ii-1)+2,'proposed'] <- mae_score(y_pred,y_test)
    
    ################################################################################
    # PCR
    ################################################################################
    oes_peak_stationary = get_peak_wv(X_total_stationary_time_wav.tmp, span=3, npeaks=20)
    
    ### oes_peak에 PCR, PLS, LASSO 적용하면 됨 ###
    peak_pca_train <- princomp(oes_peak_stationary[train_index,])
    
    peak_compdata <- peak_pca_train$scores[, 1:5]
    
    colnames(peak_compdata) <- c(paste(rep("PC", 5), 1:5, sep = ""))
    # peak_compdata %>% str() #dim: 250 * 45
    
    lm_fit <- lm(y~., data = data.frame(y = y_train, peak_compdata))
    
    peak_compdata_test <- predict(peak_pca_train, newdata = oes_peak_stationary[test_index,])[, 1:5]
    colnames(peak_compdata_test) <- c(paste(rep("PC", 5), 1:5, sep = ""))
    # str(peak_compdata_test)
    # View(peak_compdata_test)
    
    a_pca <- predict(lm_fit, data.frame(peak_compdata_test))
    result.df.total.stationary.tmp[2*(ii-1)+1,'PCR'] <- sqrt(mean((y_test-a_pca)^2))
    result.df.total.stationary.tmp[2*(ii-1)+2,'PCR'] <- mean(abs(y_test - a_pca))
    
    ################################################################################
    # PC+LASSO
    ################################################################################
    pc_lasso = cv.glmnet(x=peak_compdata,
                         y=y_train)
    a_pc_lasso = predict(pc_lasso, newx=peak_compdata_test, s=pc_lasso$lambda.min) 
    result.df.total.stationary.tmp[2*(ii-1)+1,'PCLASSO'] <- sqrt(mean((y_test-a_pc_lasso)^2))
    result.df.total.stationary.tmp[2*(ii-1)+2,'PCLASSO'] <- mean(abs(y_test - a_pc_lasso))
    
    ################################################################################
    # PLS
    ################################################################################
    pls_fit = pls::plsr(y ~ ., ncomp=4,  data = data.frame(y = y_train, peak_compdata),
                        validation="LOO")
    a_pls = predict(pls_fit, ncomp=3, data.frame(peak_compdata_test))
    a_pls = a_pls[,1,1]
    result.df.total.stationary.tmp[2*(ii-1)+1,'PLS'] <- sqrt(mean((y_test-a_pls)^2))
    result.df.total.stationary.tmp[2*(ii-1)+2,'PLS'] <- mean(abs(y_test - a_pls))
    
    ################################################################################
    # Two-way model averaging
    ################################################################################
    
    # str(oes_penfsvd) # u1, u2 (each: 307 x time point(=9)) 
    # str(oes_basisfsvd) # v1, v2 (307 x 1024)
    
    ### PenFSVD + Bspline  ###
    # time point가 너무 적어서 basis 개수 조정이 어려움 -> error발생
    penfsvd_result = jk_bsp(
      X = c(oes_penfsvd$u, oes_penfsvd$v), 
      y_ = Y_stationary, 
      x_nb = c(rep(5,ncomp), rep(50,ncomp)),
      b_nb = c(rep(5,ncomp), rep(50,ncomp)),
      tr = train_index, tst = test_index
    )
    
    result.df.total.stationary.tmp[2*(ii-1)+1,'PenFSVD_Bsp'] <- penfsvd_result$rmse
    result.df.total.stationary.tmp[2*(ii-1)+2,'PenFSVD_Bsp'] <- penfsvd_result$mae
    
    ### BasFSVD + Bspline  ###
    # time point가 너무 적어서 basis 개수 조정이 어려움 -> error발생
    basisfsvd_result = jk_bsp(
      X = c(oes_basisfsvd$u, oes_basisfsvd$v), 
      y_ = Y_stationary, 
      x_nb = c(rep(5,ncomp), rep(50,ncomp)),
      b_nb = c(rep(5,ncomp), rep(50,ncomp)),
      tr = train_index, tst = test_index
    )
    
    result.df.total.stationary.tmp[2*(ii-1)+1,'BasFSVD_Bsp'] <- basisfsvd_result$rmse
    result.df.total.stationary.tmp[2*(ii-1)+2,'BasFSVD_Bsp'] <- basisfsvd_result$mae
    
    jma_oes_res <- list()
    jma_oes_res$PenBsp <- penfsvd_result
    jma_oes_res$BasBsp <- basisfsvd_result
    
    ### PenFSVD + Wv ###
    jma_oes_res$PenWv150 <- jk_wv(
      X = oes_penfsvd$v, y_ = Y_stationary, q = 0.15,
      tr = train_index, tst = test_index
    )
    
    jma_oes_res$PenMerged150 <- jk_merge(
      error_mat1 = jma_oes_res$PenBsp$error_mat[, 1:ncomp],
      error_mat2 = jma_oes_res$PenWv150$error_mat,
      y_m1 = jma_oes_res$PenBsp$y_m[, 1:ncomp],
      y_m2 = jma_oes_res$PenWv150$y_m,
      y_ = Y_stationary, tr = train_index, tst = test_index )
    
    result.df.total.stationary.tmp[2*(ii-1)+1,'PenFSVD_Wv'] <- jma_oes_res$PenMerged150$rmse
    result.df.total.stationary.tmp[2*(ii-1)+2,'PenFSVD_Wv'] <- jma_oes_res$PenMerged150$mae
    
    ### Bas + Wv ###
    jma_oes_res$BasWv150 <- jk_wv(
      X = oes_basisfsvd$v, y_ = Y_stationary, q = 0.15,
      tr = train_index, tst = test_index
    )
    
    jma_oes_res$BasMerged150 <- jk_merge(
      error_mat1 = jma_oes_res$BasBsp$error_mat[, 1:ncomp], 
      error_mat2 = jma_oes_res$BasWv150$error_mat,
      y_m1 = jma_oes_res$BasBsp$y_m[, 1:ncomp], 
      y_m2 = jma_oes_res$BasWv150$y_m,
      y_ = Y_stationary, tr = train_index, tst = test_index)
    
    result.df.total.stationary.tmp[2*(ii-1)+1,'BasFSVD_Wv'] <- jma_oes_res$BasMerged150$rmse
    result.df.total.stationary.tmp[2*(ii-1)+2,'BasFSVD_Wv'] <- jma_oes_res$BasMerged150$mae
    
    ################################################################################
    # GFLM
    ################################################################################
    #### Setup ####
    ncomp = 5 # number of component function/2
    nX = 1 # number of matrix X
    nt = 5  # number of time grid points on [0,1]
    ns = 1024 # number of wavelength grid points [0,1]
    
    ### grid for GFLM ###
    p = ncomp*nX
    Tps = vector("list", length=p)
    for (j in 1:p) {
      Tps[[j]] = (0:(nt-1))/(nt-1)
    }
    for (j in 1:p+p) {
      Tps[[j]] = (0:(ns-1))/(ns-1)
    }
    lambda = 10^seq(2,-3,by=-0.5)
    phi = 10^seq(3,-3,by=-1)
    nphi = length(phi)
    
    X_grpl <- list()
    for(i in 1:length(c(oes_basisfsvd$u, oes_basisfsvd$v))){
      X_grpl[[i]] <- (c(oes_basisfsvd$u, oes_basisfsvd$v)[[i]])[train_index,]
    }
    
    X_grpl_test <- list()
    for(i in 1:length(c(oes_basisfsvd$u, oes_basisfsvd$v))){
      X_grpl_test[[i]] <- (c(oes_basisfsvd$u, oes_basisfsvd$v)[[i]])[test_index,]
    }
    
    cv_grpl = cv.grplFlinear(k=5, Y=y_train, X=X_grpl, Tps=Tps, lambda=lambda, phi=phi,dfs = 40)
    cv_grpl_error = apply(cv_grpl, c(2,3), sum)
    minidx = which.min(cv_grpl_error)
    minphi_id = ifelse(minidx %% nphi !=0, minidx %% nphi, minidx %% nphi + nphi)
    minlam_id = ifelse(minidx %% nphi !=0, minidx %/% nphi + 1, minidx %/% nphi)
    grpl_fit = grplFlinear(Y=y_train, X=X_grpl, Tps=Tps, lambda=lambda[minlam_id], phi=phi[minphi_id], dfs=40)
    pred = rep(grpl_fit$intercept, length(test_index))
    for (jj in 1:p) {
      pred <- pred + (1/(nt-1))*X_grpl_test[[jj]] %*% grpl_fit$Coef[[jj]][,1]
    }
    for (jj in 1:p+p) {
      pred <- pred + (1/(ns-1))*X_grpl_test[[jj]] %*% grpl_fit$Coef[[jj]][,1]
    }
    
    pred_GFLM <- as.numeric(pred)
    result.df.total.stationary.tmp[2*(ii-1)+1,'GFLM'] <- sqrt(mean((y_test-pred_GFLM)^2))
    result.df.total.stationary.tmp[2*(ii-1)+2,'GFLM'] <- mean(abs(y_test-pred_GFLM))
    
    ################################################################################
    # PFR
    ################################################################################
    pfr_df = data.frame(y=y_train, 
                        u1=I(X_grpl[[1]]), u2=I(X_grpl[[2]]), u3=I(X_grpl[[3]]), u4=I(X_grpl[[4]]), u5=I(X_grpl[[5]]),
                        v1=I(X_grpl[[6]]), v2=I(X_grpl[[7]]), v3=I(X_grpl[[8]]), v4=I(X_grpl[[9]]), v5=I(X_grpl[[10]])) 
    
    pfr_fit = pfr(y ~ lf(u1, k=10, bs="ps")+lf(u2, k=10, bs="ps")+lf(u3, k=10, bs="ps") + lf(u4, k=10, bs="ps") + lf(u5, k=10, bs="ps") + 
                    lf(v1, k=10, bs="ps") + lf(v2, k=10, bs="ps")+lf(v3, k=10, bs="ps") + lf(v4, k=10, bs="ps") + lf(v5, k=10, bs="ps"), data=pfr_df)
    
    pfr_df_test = data.frame(y=y_test, 
                             u1=I(X_grpl_test[[1]]), u2=I(X_grpl_test[[2]]), u3=I(X_grpl_test[[3]]), u4=I(X_grpl_test[[4]]), u5=I(X_grpl_test[[5]]),
                             v1=I(X_grpl_test[[6]]), v2=I(X_grpl_test[[7]]), v3=I(X_grpl_test[[8]]), v4=I(X_grpl_test[[9]]), v5=I(X_grpl_test[[10]])) 
    
    pred_pfr = as.numeric(predict(pfr_fit, newdata=pfr_df_test))
    result.df.total.stationary.tmp[2*(ii-1)+1,'PFR'] <- sqrt(mean((y_test-pred_pfr)^2))
    result.df.total.stationary.tmp[2*(ii-1)+2,'PFR'] <- mean(abs(y_test-pred_pfr))
  }
  result.df.total.stationary.list[[ttt]] <- result.df.total.stationary.tmp
}


######################################
#### Scenario II (non-stationary) ####
######################################
load("../Data/ns_stationary_true.RData")
result.df.total.ns_stationary.list <- list()

for(ttt in 1:20){
  print(paste("TOTAL iteration:", ttt))
  ns_stationary_simuldata.list.tmp <- list()
  x_err_std <- 3
  set.seed(20)
  for(i in 1:300){
    # ns_stationary_simuldata.list[[i]] <- t(t(ns_stationary_true)+rnorm(ncol(ns_stationary_true), 0, x_err_std))
    ns_stationary_simuldata.list.tmp[[i]] <- ns_stationary_true + matrix(rnorm(5*1024, 0, x_err_std), nrow=5, ncol=1024)
  }
  num_datasets <- 300
  X_total_ns_stationary.tmp <- array(NA, dim = c(num_datasets, 1024, 5))
  
  # Generate 300 different datasets
  for (i in 1:num_datasets) {
    X_total_ns_stationary.tmp[i, , ] <- t(ns_stationary_simuldata.list.tmp[[i]])
  }
  
  X_total_ns_stationary_time_wav.tmp <- array(NA, dim = c(num_datasets, 5, 1024))
  for (i in 1:num_datasets) {
    X_total_ns_stationary_time_wav.tmp[i, , ] <- ns_stationary_simuldata.list.tmp[[i]]
  }
  
  num_datasets <- dim(X_total_ns_stationary.tmp)[1]
  Y_ns_stationary <- numeric(num_datasets)
  Z_ns_stationary.tmp <- numeric(num_datasets)
  # Create time and wavelength axes with a small random noise added
  time_axis <- sin(2*pi/5*c(c(1:5)-1/3))
  wavelength_axis <- sin(2*pi/256*c(c(1:1024)-1/3))
  # Compute the outer product for the i-th dataset
  beta <- outer(time_axis, wavelength_axis, FUN = "*")
  
  for (i in 1:num_datasets) {
    Z_ns_stationary.tmp[i] <- sum((X_total_ns_stationary.tmp[i,,]-mean(X_total_ns_stationary.tmp[i,,])) * t(beta))
  }
  
  SNR_list <- c(2,5,10,20)
  y_err_std_list <- sqrt(var(Z_ns_stationary.tmp)/SNR_list)
  result.df.total.ns_stationary.tmp <- data.frame(row.names = paste(c("RMSE", "MAE"), "(SNR=", rep(SNR_list, each=2), ")", sep=""))
  
  ################################################################################
  # Two-way model averaging
  ################################################################################
  ncomp = 5 # SVD (u,v) pair 2개로 가정
  ## PenFSVD ##
  penfsvd = apply(X_total_ns_stationary_time_wav.tmp, 1, tw_svd_smu, ncomp=ncomp)
  
  ## basisFSVD ##
  basfsvd = apply(X_total_ns_stationary_time_wav.tmp, 1, basis_svd, nb_by = 1, ncmp=ncomp)
  
  nt <- dim(X_total_ns_stationary_time_wav.tmp)[2]
  nlambda <- dim(X_total_ns_stationary_time_wav.tmp)[3]
  
  ## sign adjustment (flip) ##
  pen_u1 = matrix(0, nrow=num_datasets, ncol=nt)
  for (i in 1:num_datasets) {
    pen_u1[i,] = penfsvd[[i]]$pen_u[1:nt,1]
  }
  for(i in which(pen_u1[,1]<0)){
    penfsvd[[i]]$pen_u[1:nt,1] <- - penfsvd[[i]]$pen_u[1:nt,1]
  }
  # matplot(t(pen_u1), type="l") # sign이 다르게 나오는 wafer가 있으면 조정해주기
  pen_u2 = matrix(0, nrow=num_datasets, ncol=nt)
  for (i in 1:num_datasets) {
    pen_u2[i,] = penfsvd[[i]]$pen_u[1:nt,2]
  }
  # matplot(t(pen_u2), type="l")
  for(i in which(pen_u2[,1]<0)){
    penfsvd[[i]]$pen_u[1:nt,2] <- - penfsvd[[i]]$pen_u[1:nt,2]
  }
  
  pen_u3 = matrix(0, nrow=num_datasets, ncol=nt)
  for (i in 1:num_datasets) {
    pen_u3[i,] = penfsvd[[i]]$pen_u[1:nt,3]
  }
  # matplot(t(pen_u3), type="l")
  for(i in which(pen_u3[,1]<0)){
    penfsvd[[i]]$pen_u[1:nt,3] <- - penfsvd[[i]]$pen_u[1:nt,3]
  }
  
  pen_u4 = matrix(0, nrow=num_datasets, ncol=nt)
  for (i in 1:num_datasets) {
    pen_u4[i,] = penfsvd[[i]]$pen_u[1:nt,4]
  }
  # matplot(t(pen_u4), type="l")
  for(i in which(pen_u4[,2]<0)){
    penfsvd[[i]]$pen_u[1:nt,4] <- - penfsvd[[i]]$pen_u[1:nt,4]
  }
  
  pen_u5 = matrix(0, nrow=num_datasets, ncol=nt)
  for (i in 1:num_datasets) {
    pen_u5[i,] = penfsvd[[i]]$pen_u[1:nt,5]
  }
  # matplot(t(pen_u5), type="l")
  for(i in which(pen_u5[,1]<0)){
    penfsvd[[i]]$pen_u[1:nt,5] <- - penfsvd[[i]]$pen_u[1:nt,5]
  }
  
  pen_v1 = matrix(0, nrow=num_datasets, ncol=nlambda)
  for (i in 1:num_datasets) {
    pen_v1[i,] = penfsvd[[i]]$pen_v[,1]
  }
  for(i in which(pen_v1[,1]<0)){
    penfsvd[[i]]$pen_v[1:nlambda,1] <- - penfsvd[[i]]$pen_v[1:nlambda,1]
  }
  # matplot(t(pen_v1), type="l") # sign이 다르게 나오는 wafer가 있으면 조정해주기
  pen_v2 = matrix(0, nrow=num_datasets, ncol=nlambda)
  for (i in 1:num_datasets) {
    pen_v2[i,] = penfsvd[[i]]$pen_v[,2]
  }
  # matplot(t(pen_v2), type="l")
  for(i in which(pen_v2[,1]<0)){
    penfsvd[[i]]$pen_v[1:nlambda,2] <- - penfsvd[[i]]$pen_v[1:nlambda,2]
  }
  
  pen_v3 = matrix(0, nrow=num_datasets, ncol=nlambda)
  for (i in 1:num_datasets) {
    pen_v3[i,] = penfsvd[[i]]$pen_v[,3]
  }
  # matplot(t(pen_v3), type="l")
  for(i in which(pen_v3[,1]<0)){
    penfsvd[[i]]$pen_v[1:nlambda,3] <- - penfsvd[[i]]$pen_v[1:nlambda,3]
  }
  # for(i in which(pen_v3[,15]<0)){
  #   penfsvd[[i]]$pen_v[,3] <- - penfsvd[[i]]$pen_v[,3]
  # }
  
  pen_v4 = matrix(0, nrow=num_datasets, ncol=nlambda)
  for (i in 1:num_datasets) {
    pen_v4[i,] = penfsvd[[i]]$pen_v[,4]
  }
  # matplot(t(pen_v4), type="l")
  for(i in which(pen_v4[,1]<0)){
    penfsvd[[i]]$pen_v[1:nlambda,4] <- - penfsvd[[i]]$pen_v[1:nlambda,4]
  }
  
  pen_v5 = matrix(0, nrow=num_datasets, ncol=nlambda)
  for (i in 1:num_datasets) {
    pen_v5[i,] = penfsvd[[i]]$pen_v[,5]
  }
  # matplot(t(pen_v5), type="l")
  for(i in which(pen_v5[,1]<0)){
    penfsvd[[i]]$pen_v[1:nlambda,5] <- - penfsvd[[i]]$pen_v[1:nlambda,5]
  }
  
  # for(i in which(pen_v5[,1]<0)){
  #   penfsvd[[i]]$pen_v[,5] <- - penfsvd[[i]]$pen_v[,5]
  # }
  
  # basis fsvd 결과도 필요하면 sign 조정
  
  ### join the result ###
  pen_u1 = do.call(rbind, lapply(penfsvd, function(x) x$pen_u[,1]))
  pen_u2 = do.call(rbind, lapply(penfsvd, function(x) x$pen_u[,2]))
  pen_v1 = do.call(rbind, lapply(penfsvd, function(x) x$pen_v[,1]))
  pen_v2 = do.call(rbind, lapply(penfsvd, function(x) x$pen_v[,2]))
  
  pen_u3 = do.call(rbind, lapply(penfsvd, function(x) x$pen_u[,3]))
  pen_u4 = do.call(rbind, lapply(penfsvd, function(x) x$pen_u[,4]))
  pen_v3 = do.call(rbind, lapply(penfsvd, function(x) x$pen_v[,3]))
  pen_v4 = do.call(rbind, lapply(penfsvd, function(x) x$pen_v[,4]))
  
  pen_u5 = do.call(rbind, lapply(penfsvd, function(x) x$pen_u[,5]))
  pen_v5 = do.call(rbind, lapply(penfsvd, function(x) x$pen_v[,5]))
  
  
  oes_penfsvd <- list(u = vector("list", length = ncomp),
                      v = vector("list", length = ncomp))
  oes_penfsvd$u[[1]] = pen_u1
  oes_penfsvd$u[[2]] = pen_u2
  oes_penfsvd$v[[1]] = pen_v1
  oes_penfsvd$v[[2]] = pen_v2
  
  oes_penfsvd$u[[3]] = pen_u3
  oes_penfsvd$u[[4]] = pen_u4
  oes_penfsvd$v[[3]] = pen_v3
  oes_penfsvd$v[[4]] = pen_v4
  
  oes_penfsvd$u[[5]] = pen_u5
  oes_penfsvd$v[[5]] = pen_v5
  
  bas_u1 = do.call(rbind, lapply(basfsvd, function(x) x$u[,1]))
  bas_u2 = do.call(rbind, lapply(basfsvd, function(x) x$u[,2]))
  bas_v1 = do.call(rbind, lapply(basfsvd, function(x) x$v[,1]))
  bas_v2 = do.call(rbind, lapply(basfsvd, function(x) x$v[,2]))
  
  bas_u3 = do.call(rbind, lapply(basfsvd, function(x) x$u[,3]))
  bas_u4 = do.call(rbind, lapply(basfsvd, function(x) x$u[,4]))
  bas_v3 = do.call(rbind, lapply(basfsvd, function(x) x$v[,3]))
  bas_v4 = do.call(rbind, lapply(basfsvd, function(x) x$v[,4]))
  
  bas_u5 = do.call(rbind, lapply(basfsvd, function(x) x$u[,5]))
  bas_v5 = do.call(rbind, lapply(basfsvd, function(x) x$v[,5]))
  
  oes_basisfsvd <- list(u = vector("list", length = ncomp),
                        v = vector("list", length = ncomp))
  oes_basisfsvd$u[[1]] = bas_u1
  oes_basisfsvd$u[[2]] = bas_u2
  oes_basisfsvd$v[[1]] = bas_v1
  oes_basisfsvd$v[[2]] = bas_v2
  
  oes_basisfsvd$u[[3]] = bas_u3
  oes_basisfsvd$u[[4]] = bas_u4
  oes_basisfsvd$v[[3]] = bas_v3
  oes_basisfsvd$v[[4]] = bas_v4
  
  oes_basisfsvd$u[[5]] = bas_u5
  oes_basisfsvd$v[[5]] = bas_v5
  
  for(ii in 1:length(y_err_std_list)){
    print(paste("iteration:", ii))
    y_err_std <- y_err_std_list[ii]
    set.seed(20)
    for (i in 1:num_datasets) {
      Y_ns_stationary[i] <- sum((X_total_ns_stationary.tmp[i,,]-mean(X_total_ns_stationary.tmp[i,,])) * t(beta)) + rnorm(1, mean = 0, sd = y_err_std)  # adding some noise
    }
    
    ################################################################################
    # MFPLS
    ################################################################################
    set.seed(ttt)
    train_index = sample(1:num_datasets, num_datasets*0.8) 
    test_index = (1:num_datasets)[-train_index]
    
    y_train <- Y_ns_stationary[train_index]
    y_test <- Y_ns_stationary[test_index]
    
    X_train <- X_total_ns_stationary.tmp[train_index,,]
    X_test <- X_total_ns_stationary.tmp[test_index,,]
    
    
    basis_num <- seq(30,60,15) ; basis_num
    Basis_name <- c('Basis 30','Basis 45','Basis 60')
    parameter.change.result.df <- data.frame(row.names = c('pls_num','mean_RMSE'))
    num_folds <- 10 # 10은 시간이 너무 많이 걸림 2로했을 때 약 28분 소요
    n <- dim(X_train)[1] 
    fold_indices <- createFolds(1:n, k = num_folds) 
    
    parameter.change.result.df_ns_stationary <- parameter_selection_MFPLS(X_train,y_train,basis_num,Basis_name,num_folds,fold_indices)
    i_th <- which.min(parameter.change.result.df_ns_stationary[2,]) # select argmin RMSE
    basis <- basis_num[i_th]
    pls_component <- parameter.change.result.df_ns_stationary[1,i_th]
    
    # see test errors
    MFPLS_model <- MFPLS(conv2fda(X_train),y_train, ncomp=pls_component,K=basis)
    prd <- predict.MFPLS(MFPLS_model, conv2fda(X_test))
    y_pred <- prd$y_pred
    # calculate accuracy
    result.df.total.ns_stationary.tmp[2*(ii-1)+1,'proposed'] <- rmse_score(y_pred,y_test)
    result.df.total.ns_stationary.tmp[2*(ii-1)+2,'proposed'] <- mae_score(y_pred,y_test)
    
    ################################################################################
    # PCR
    ################################################################################
    oes_peak_ns_stationary = get_peak_wv(X_total_ns_stationary_time_wav.tmp, span=3, npeaks=20)
    
    ### oes_peak에 PCR, PLS, LASSO 적용하면 됨 ###
    peak_pca_train <- princomp(oes_peak_ns_stationary[train_index,])
    
    peak_compdata <- peak_pca_train$scores[, 1:5]
    
    colnames(peak_compdata) <- c(paste(rep("PC", 5), 1:5, sep = ""))
    # peak_compdata %>% str() #dim: 250 * 45
    
    lm_fit <- lm(y~., data = data.frame(y = y_train, peak_compdata))
    
    peak_compdata_test <- predict(peak_pca_train, newdata = oes_peak_ns_stationary[test_index,])[, 1:5]
    colnames(peak_compdata_test) <- c(paste(rep("PC", 5), 1:5, sep = ""))
    # str(peak_compdata_test)
    # View(peak_compdata_test)
    
    a_pca <- predict(lm_fit, data.frame(peak_compdata_test))
    result.df.total.ns_stationary.tmp[2*(ii-1)+1,'PCR'] <- sqrt(mean((y_test-a_pca)^2))
    result.df.total.ns_stationary.tmp[2*(ii-1)+2,'PCR'] <- mean(abs(y_test - a_pca))
    
    ################################################################################
    # PC+LASSO
    ################################################################################
    pc_lasso = cv.glmnet(x=peak_compdata,
                         y=y_train)
    a_pc_lasso = predict(pc_lasso, newx=peak_compdata_test, s=pc_lasso$lambda.min) 
    result.df.total.ns_stationary.tmp[2*(ii-1)+1,'PCLASSO'] <- sqrt(mean((y_test-a_pc_lasso)^2))
    result.df.total.ns_stationary.tmp[2*(ii-1)+2,'PCLASSO'] <- mean(abs(y_test - a_pc_lasso))
    
    ################################################################################
    # PLS
    ################################################################################
    pls_fit = pls::plsr(y ~ ., ncomp=4,  data = data.frame(y = y_train, peak_compdata),
                        validation="LOO")
    a_pls = predict(pls_fit, ncomp=3, data.frame(peak_compdata_test))
    a_pls = a_pls[,1,1]
    result.df.total.ns_stationary.tmp[2*(ii-1)+1,'PLS'] <- sqrt(mean((y_test-a_pls)^2))
    result.df.total.ns_stationary.tmp[2*(ii-1)+2,'PLS'] <- mean(abs(y_test - a_pls))
    
    ################################################################################
    # Two-way model averaging
    ################################################################################
    
    # str(oes_penfsvd) # u1, u2 (each: 307 x time point(=9)) 
    # str(oes_basisfsvd) # v1, v2 (307 x 1024)
    
    ### PenFSVD + Bspline  ###
    # time point가 너무 적어서 basis 개수 조정이 어려움 -> error발생
    penfsvd_result = jk_bsp(
      X = c(oes_penfsvd$u, oes_penfsvd$v), 
      y_ = Y_ns_stationary, 
      x_nb = c(rep(5,ncomp), rep(50,ncomp)),
      b_nb = c(rep(5,ncomp), rep(50,ncomp)),
      tr = train_index, tst = test_index
    )
    
    result.df.total.ns_stationary.tmp[2*(ii-1)+1,'PenFSVD_Bsp'] <- penfsvd_result$rmse
    result.df.total.ns_stationary.tmp[2*(ii-1)+2,'PenFSVD_Bsp'] <- penfsvd_result$mae
    
    ### BasFSVD + Bspline  ###
    # time point가 너무 적어서 basis 개수 조정이 어려움 -> error발생
    basisfsvd_result = jk_bsp(
      X = c(oes_basisfsvd$u, oes_basisfsvd$v), 
      y_ = Y_ns_stationary, 
      x_nb = c(rep(5,ncomp), rep(50,ncomp)),
      b_nb = c(rep(5,ncomp), rep(50,ncomp)),
      tr = train_index, tst = test_index
    )
    
    result.df.total.ns_stationary.tmp[2*(ii-1)+1,'BasFSVD_Bsp'] <- basisfsvd_result$rmse
    result.df.total.ns_stationary.tmp[2*(ii-1)+2,'BasFSVD_Bsp'] <- basisfsvd_result$mae
    
    jma_oes_res <- list()
    jma_oes_res$PenBsp <- penfsvd_result
    jma_oes_res$BasBsp <- basisfsvd_result
    
    ### PenFSVD + Wv ###
    jma_oes_res$PenWv150 <- jk_wv(
      X = oes_penfsvd$v, y_ = Y_ns_stationary, q = 0.15,
      tr = train_index, tst = test_index
    )
    
    jma_oes_res$PenMerged150 <- jk_merge(
      error_mat1 = jma_oes_res$PenBsp$error_mat[, 1:ncomp],
      error_mat2 = jma_oes_res$PenWv150$error_mat,
      y_m1 = jma_oes_res$PenBsp$y_m[, 1:ncomp],
      y_m2 = jma_oes_res$PenWv150$y_m,
      y_ = Y_ns_stationary, tr = train_index, tst = test_index )
    
    result.df.total.ns_stationary.tmp[2*(ii-1)+1,'PenFSVD_Wv'] <- jma_oes_res$PenMerged150$rmse
    result.df.total.ns_stationary.tmp[2*(ii-1)+2,'PenFSVD_Wv'] <- jma_oes_res$PenMerged150$mae
    
    ### Bas + Wv ###
    jma_oes_res$BasWv150 <- jk_wv(
      X = oes_basisfsvd$v, y_ = Y_ns_stationary, q = 0.15,
      tr = train_index, tst = test_index
    )
    
    jma_oes_res$BasMerged150 <- jk_merge(
      error_mat1 = jma_oes_res$BasBsp$error_mat[, 1:ncomp], 
      error_mat2 = jma_oes_res$BasWv150$error_mat,
      y_m1 = jma_oes_res$BasBsp$y_m[, 1:ncomp], 
      y_m2 = jma_oes_res$BasWv150$y_m,
      y_ = Y_ns_stationary, tr = train_index, tst = test_index)
    
    result.df.total.ns_stationary.tmp[2*(ii-1)+1,'BasFSVD_Wv'] <- jma_oes_res$BasMerged150$rmse
    result.df.total.ns_stationary.tmp[2*(ii-1)+2,'BasFSVD_Wv'] <- jma_oes_res$BasMerged150$mae
    
    ################################################################################
    # GFLM
    ################################################################################
    #### Setup ####
    ncomp = 5 # number of component function/2
    nX = 1 # number of matrix X
    nt = 5  # number of time grid points on [0,1]
    ns = 1024 # number of wavelength grid points [0,1]
    
    ### grid for GFLM ###
    p = ncomp*nX
    Tps = vector("list", length=p)
    for (j in 1:p) {
      Tps[[j]] = (0:(nt-1))/(nt-1)
    }
    for (j in 1:p+p) {
      Tps[[j]] = (0:(ns-1))/(ns-1)
    }
    lambda = 10^seq(2,-3,by=-0.5)
    phi = 10^seq(3,-3,by=-1)
    nphi = length(phi)
    
    X_grpl <- list()
    for(i in 1:length(c(oes_basisfsvd$u, oes_basisfsvd$v))){
      X_grpl[[i]] <- (c(oes_basisfsvd$u, oes_basisfsvd$v)[[i]])[train_index,]
    }
    
    X_grpl_test <- list()
    for(i in 1:length(c(oes_basisfsvd$u, oes_basisfsvd$v))){
      X_grpl_test[[i]] <- (c(oes_basisfsvd$u, oes_basisfsvd$v)[[i]])[test_index,]
    }
    
    cv_grpl = cv.grplFlinear(k=5, Y=y_train, X=X_grpl, Tps=Tps, lambda=lambda, phi=phi,dfs = 40)
    cv_grpl_error = apply(cv_grpl, c(2,3), sum)
    minidx = which.min(cv_grpl_error)
    minphi_id = ifelse(minidx %% nphi !=0, minidx %% nphi, minidx %% nphi + nphi)
    minlam_id = ifelse(minidx %% nphi !=0, minidx %/% nphi + 1, minidx %/% nphi)
    grpl_fit = grplFlinear(Y=y_train, X=X_grpl, Tps=Tps, lambda=lambda[minlam_id], phi=phi[minphi_id], dfs=40)
    pred = rep(grpl_fit$intercept, length(test_index))
    for (jj in 1:p) {
      pred <- pred + (1/(nt-1))*X_grpl_test[[jj]] %*% grpl_fit$Coef[[jj]][,1]
    }
    for (jj in 1:p+p) {
      pred <- pred + (1/(ns-1))*X_grpl_test[[jj]] %*% grpl_fit$Coef[[jj]][,1]
    }
    
    pred_GFLM <- as.numeric(pred)
    result.df.total.ns_stationary.tmp[2*(ii-1)+1,'GFLM'] <- sqrt(mean((y_test-pred_GFLM)^2))
    result.df.total.ns_stationary.tmp[2*(ii-1)+2,'GFLM'] <- mean(abs(y_test-pred_GFLM))
    
    ################################################################################
    # PFR
    ################################################################################
    pfr_df = data.frame(y=y_train, 
                        u1=I(X_grpl[[1]]), u2=I(X_grpl[[2]]), u3=I(X_grpl[[3]]), u4=I(X_grpl[[4]]), u5=I(X_grpl[[5]]),
                        v1=I(X_grpl[[6]]), v2=I(X_grpl[[7]]), v3=I(X_grpl[[8]]), v4=I(X_grpl[[9]]), v5=I(X_grpl[[10]])) 
    
    pfr_fit = pfr(y ~ lf(u1, k=10, bs="ps")+lf(u2, k=10, bs="ps")+lf(u3, k=10, bs="ps") + lf(u4, k=10, bs="ps") + lf(u5, k=10, bs="ps") + 
                    lf(v1, k=10, bs="ps") + lf(v2, k=10, bs="ps")+lf(v3, k=10, bs="ps") + lf(v4, k=10, bs="ps") + lf(v5, k=10, bs="ps"), data=pfr_df)
    
    pfr_df_test = data.frame(y=y_test, 
                             u1=I(X_grpl_test[[1]]), u2=I(X_grpl_test[[2]]), u3=I(X_grpl_test[[3]]), u4=I(X_grpl_test[[4]]), u5=I(X_grpl_test[[5]]),
                             v1=I(X_grpl_test[[6]]), v2=I(X_grpl_test[[7]]), v3=I(X_grpl_test[[8]]), v4=I(X_grpl_test[[9]]), v5=I(X_grpl_test[[10]])) 
    
    pred_pfr = as.numeric(predict(pfr_fit, newdata=pfr_df_test))
    result.df.total.ns_stationary.tmp[2*(ii-1)+1,'PFR'] <- sqrt(mean((y_test-pred_pfr)^2))
    result.df.total.ns_stationary.tmp[2*(ii-1)+2,'PFR'] <- mean(abs(y_test-pred_pfr))
  }
  result.df.total.ns_stationary.list[[ttt]] <- result.df.total.ns_stationary.tmp
}


######################################
#### Scenario III (Fourier basis) ####
######################################
load("../Data/sc_true.RData")
result.df.total.sc.list <- list()

for(ttt in 1:20){
  print(paste("TOTAL iteration:", ttt))
  sc_simuldata.list.tmp <- list()
  x_err_std <- 3
  set.seed(30)
  for(i in 1:300){
    # sc_simuldata.list[[i]] <- t(t(sc_true)+rnorm(ncol(sc_true), 0, x_err_std))
    sc_simuldata.list.tmp[[i]] <- sc_true + matrix(rnorm(5*1024, 0, x_err_std), nrow=5, ncol=1024)
  }
  num_datasets <- 300
  X_total_sc.tmp <- array(NA, dim = c(num_datasets, 1024, 5))
  
  # Generate 300 different datasets
  for (i in 1:num_datasets) {
    X_total_sc.tmp[i, , ] <- t(sc_simuldata.list.tmp[[i]])
  }
  
  X_total_sc_time_wav.tmp <- array(NA, dim = c(num_datasets, 5, 1024))
  for (i in 1:num_datasets) {
    X_total_sc_time_wav.tmp[i, , ] <- sc_simuldata.list.tmp[[i]]
  }
  
  num_datasets <- dim(X_total_sc.tmp)[1]
  Y_sc <- numeric(num_datasets)
  Z_sc.tmp <- numeric(num_datasets)
  # Create time and wavelength axes with a small random noise added
  time_axis <- sin(2*pi/5*c(c(1:5)-1/3))
  wavelength_axis <- sin(2*pi/256*c(c(1:1024)-1/3))
  # Compute the outer product for the i-th dataset
  beta <- outer(time_axis, wavelength_axis, FUN = "*")
  
  for (i in 1:num_datasets) {
    Z_sc.tmp[i] <- sum((X_total_sc.tmp[i,,]-mean(X_total_sc.tmp[i,,])) * t(beta))
  }
  
  SNR_list <- c(2,5,10,20)
  y_err_std_list <- sqrt(var(Z_sc.tmp)/SNR_list)
  result.df.total.sc.tmp <- data.frame(row.names = paste(c("RMSE", "MAE"), "(SNR=", rep(SNR_list, each=2), ")", sep=""))
  
  ################################################################################
  # Two-way model averaging
  ################################################################################
  ncomp = 5 # SVD (u,v) pair 2개로 가정
  ## PenFSVD ##
  penfsvd = apply(X_total_sc_time_wav.tmp, 1, tw_svd_smu, ncomp=ncomp)
  
  ## basisFSVD ##
  basfsvd = apply(X_total_sc_time_wav.tmp, 1, basis_svd, nb_by = 1, ncmp=ncomp)
  
  nt <- dim(X_total_sc_time_wav.tmp)[2]
  nlambda <- dim(X_total_sc_time_wav.tmp)[3]
  
  ## sign adjustment (flip) ##
  pen_u1 = matrix(0, nrow=num_datasets, ncol=nt)
  for (i in 1:num_datasets) {
    pen_u1[i,] = penfsvd[[i]]$pen_u[1:nt,1]
  }
  for(i in which(pen_u1[,1]<0)){
    penfsvd[[i]]$pen_u[1:nt,1] <- - penfsvd[[i]]$pen_u[1:nt,1]
  }
  # matplot(t(pen_u1), type="l") # sign이 다르게 나오는 wafer가 있으면 조정해주기
  pen_u2 = matrix(0, nrow=num_datasets, ncol=nt)
  for (i in 1:num_datasets) {
    pen_u2[i,] = penfsvd[[i]]$pen_u[1:nt,2]
  }
  # matplot(t(pen_u2), type="l")
  for(i in which(pen_u2[,1]<0)){
    penfsvd[[i]]$pen_u[1:nt,2] <- - penfsvd[[i]]$pen_u[1:nt,2]
  }
  
  pen_u3 = matrix(0, nrow=num_datasets, ncol=nt)
  for (i in 1:num_datasets) {
    pen_u3[i,] = penfsvd[[i]]$pen_u[1:nt,3]
  }
  # matplot(t(pen_u3), type="l")
  for(i in which(pen_u3[,1]<0)){
    penfsvd[[i]]$pen_u[1:nt,3] <- - penfsvd[[i]]$pen_u[1:nt,3]
  }
  
  pen_u4 = matrix(0, nrow=num_datasets, ncol=nt)
  for (i in 1:num_datasets) {
    pen_u4[i,] = penfsvd[[i]]$pen_u[1:nt,4]
  }
  # matplot(t(pen_u4), type="l")
  for(i in which(pen_u4[,2]<0)){
    penfsvd[[i]]$pen_u[1:nt,4] <- - penfsvd[[i]]$pen_u[1:nt,4]
  }
  
  pen_u5 = matrix(0, nrow=num_datasets, ncol=nt)
  for (i in 1:num_datasets) {
    pen_u5[i,] = penfsvd[[i]]$pen_u[1:nt,5]
  }
  # matplot(t(pen_u5), type="l")
  for(i in which(pen_u5[,1]<0)){
    penfsvd[[i]]$pen_u[1:nt,5] <- - penfsvd[[i]]$pen_u[1:nt,5]
  }
  
  pen_v1 = matrix(0, nrow=num_datasets, ncol=nlambda)
  for (i in 1:num_datasets) {
    pen_v1[i,] = penfsvd[[i]]$pen_v[,1]
  }
  for(i in which(pen_v1[,1]<0)){
    penfsvd[[i]]$pen_v[1:nlambda,1] <- - penfsvd[[i]]$pen_v[1:nlambda,1]
  }
  # matplot(t(pen_v1), type="l") # sign이 다르게 나오는 wafer가 있으면 조정해주기
  pen_v2 = matrix(0, nrow=num_datasets, ncol=nlambda)
  for (i in 1:num_datasets) {
    pen_v2[i,] = penfsvd[[i]]$pen_v[,2]
  }
  # matplot(t(pen_v2), type="l")
  for(i in which(pen_v2[,1]<0)){
    penfsvd[[i]]$pen_v[1:nlambda,2] <- - penfsvd[[i]]$pen_v[1:nlambda,2]
  }
  
  pen_v3 = matrix(0, nrow=num_datasets, ncol=nlambda)
  for (i in 1:num_datasets) {
    pen_v3[i,] = penfsvd[[i]]$pen_v[,3]
  }
  # matplot(t(pen_v3), type="l")
  for(i in which(pen_v3[,1]<0)){
    penfsvd[[i]]$pen_v[1:nlambda,3] <- - penfsvd[[i]]$pen_v[1:nlambda,3]
  }
  # for(i in which(pen_v3[,15]<0)){
  #   penfsvd[[i]]$pen_v[,3] <- - penfsvd[[i]]$pen_v[,3]
  # }
  
  pen_v4 = matrix(0, nrow=num_datasets, ncol=nlambda)
  for (i in 1:num_datasets) {
    pen_v4[i,] = penfsvd[[i]]$pen_v[,4]
  }
  # matplot(t(pen_v4), type="l")
  for(i in which(pen_v4[,1]<0)){
    penfsvd[[i]]$pen_v[1:nlambda,4] <- - penfsvd[[i]]$pen_v[1:nlambda,4]
  }
  
  pen_v5 = matrix(0, nrow=num_datasets, ncol=nlambda)
  for (i in 1:num_datasets) {
    pen_v5[i,] = penfsvd[[i]]$pen_v[,5]
  }
  # matplot(t(pen_v5), type="l")
  for(i in which(pen_v5[,1]<0)){
    penfsvd[[i]]$pen_v[1:nlambda,5] <- - penfsvd[[i]]$pen_v[1:nlambda,5]
  }
  
  # for(i in which(pen_v5[,1]<0)){
  #   penfsvd[[i]]$pen_v[,5] <- - penfsvd[[i]]$pen_v[,5]
  # }
  
  # basis fsvd 결과도 필요하면 sign 조정
  
  ### join the result ###
  pen_u1 = do.call(rbind, lapply(penfsvd, function(x) x$pen_u[,1]))
  pen_u2 = do.call(rbind, lapply(penfsvd, function(x) x$pen_u[,2]))
  pen_v1 = do.call(rbind, lapply(penfsvd, function(x) x$pen_v[,1]))
  pen_v2 = do.call(rbind, lapply(penfsvd, function(x) x$pen_v[,2]))
  
  pen_u3 = do.call(rbind, lapply(penfsvd, function(x) x$pen_u[,3]))
  pen_u4 = do.call(rbind, lapply(penfsvd, function(x) x$pen_u[,4]))
  pen_v3 = do.call(rbind, lapply(penfsvd, function(x) x$pen_v[,3]))
  pen_v4 = do.call(rbind, lapply(penfsvd, function(x) x$pen_v[,4]))
  
  pen_u5 = do.call(rbind, lapply(penfsvd, function(x) x$pen_u[,5]))
  pen_v5 = do.call(rbind, lapply(penfsvd, function(x) x$pen_v[,5]))
  
  
  oes_penfsvd <- list(u = vector("list", length = ncomp),
                      v = vector("list", length = ncomp))
  oes_penfsvd$u[[1]] = pen_u1
  oes_penfsvd$u[[2]] = pen_u2
  oes_penfsvd$v[[1]] = pen_v1
  oes_penfsvd$v[[2]] = pen_v2
  
  oes_penfsvd$u[[3]] = pen_u3
  oes_penfsvd$u[[4]] = pen_u4
  oes_penfsvd$v[[3]] = pen_v3
  oes_penfsvd$v[[4]] = pen_v4
  
  oes_penfsvd$u[[5]] = pen_u5
  oes_penfsvd$v[[5]] = pen_v5
  
  bas_u1 = do.call(rbind, lapply(basfsvd, function(x) x$u[,1]))
  bas_u2 = do.call(rbind, lapply(basfsvd, function(x) x$u[,2]))
  bas_v1 = do.call(rbind, lapply(basfsvd, function(x) x$v[,1]))
  bas_v2 = do.call(rbind, lapply(basfsvd, function(x) x$v[,2]))
  
  bas_u3 = do.call(rbind, lapply(basfsvd, function(x) x$u[,3]))
  bas_u4 = do.call(rbind, lapply(basfsvd, function(x) x$u[,4]))
  bas_v3 = do.call(rbind, lapply(basfsvd, function(x) x$v[,3]))
  bas_v4 = do.call(rbind, lapply(basfsvd, function(x) x$v[,4]))
  
  bas_u5 = do.call(rbind, lapply(basfsvd, function(x) x$u[,5]))
  bas_v5 = do.call(rbind, lapply(basfsvd, function(x) x$v[,5]))
  
  oes_basisfsvd <- list(u = vector("list", length = ncomp),
                        v = vector("list", length = ncomp))
  oes_basisfsvd$u[[1]] = bas_u1
  oes_basisfsvd$u[[2]] = bas_u2
  oes_basisfsvd$v[[1]] = bas_v1
  oes_basisfsvd$v[[2]] = bas_v2
  
  oes_basisfsvd$u[[3]] = bas_u3
  oes_basisfsvd$u[[4]] = bas_u4
  oes_basisfsvd$v[[3]] = bas_v3
  oes_basisfsvd$v[[4]] = bas_v4
  
  oes_basisfsvd$u[[5]] = bas_u5
  oes_basisfsvd$v[[5]] = bas_v5
  
  for(ii in 1:length(y_err_std_list)){
    print(paste("iteration:", ii))
    y_err_std <- y_err_std_list[ii]
    set.seed(30)
    for (i in 1:num_datasets) {
      Y_sc[i] <- sum((X_total_sc.tmp[i,,]-mean(X_total_sc.tmp[i,,])) * t(beta)) + rnorm(1, mean = 0, sd = y_err_std)  # adding some noise
    }
    
    ################################################################################
    # MFPLS
    ################################################################################
    set.seed(ttt)
    train_index = sample(1:num_datasets, num_datasets*0.8) 
    test_index = (1:num_datasets)[-train_index]
    
    y_train <- Y_sc[train_index]
    y_test <- Y_sc[test_index]
    
    X_train <- X_total_sc.tmp[train_index,,]
    X_test <- X_total_sc.tmp[test_index,,]
    
    
    basis_num <- seq(30,60,15) ; basis_num
    Basis_name <- c('Basis 30','Basis 45','Basis 60')
    parameter.change.result.df <- data.frame(row.names = c('pls_num','mean_RMSE'))
    num_folds <- 10 # 10은 시간이 너무 많이 걸림 2로했을 때 약 28분 소요
    n <- dim(X_train)[1] 
    fold_indices <- createFolds(1:n, k = num_folds) 
    
    parameter.change.result.df_sc <- parameter_selection_MFPLS(X_train,y_train,basis_num,Basis_name,num_folds,fold_indices)
    i_th <- which.min(parameter.change.result.df_sc[2,]) # select argmin RMSE
    basis <- basis_num[i_th]
    pls_component <- parameter.change.result.df_sc[1,i_th]
    
    # see test errors
    MFPLS_model <- MFPLS(conv2fda(X_train),y_train, ncomp=pls_component,K=basis)
    prd <- predict.MFPLS(MFPLS_model, conv2fda(X_test))
    y_pred <- prd$y_pred
    # calculate accuracy
    result.df.total.sc.tmp[2*(ii-1)+1,'proposed'] <- rmse_score(y_pred,y_test)
    result.df.total.sc.tmp[2*(ii-1)+2,'proposed'] <- mae_score(y_pred,y_test)
    
    ################################################################################
    # PCR
    ################################################################################
    oes_peak_sc = get_peak_wv(X_total_sc_time_wav.tmp, span=3, npeaks=20)
    
    ### oes_peak에 PCR, PLS, LASSO 적용하면 됨 ###
    peak_pca_train <- princomp(oes_peak_sc[train_index,])
    
    peak_compdata <- peak_pca_train$scores[, 1:5]
    
    colnames(peak_compdata) <- c(paste(rep("PC", 5), 1:5, sep = ""))
    # peak_compdata %>% str() #dim: 250 * 45
    
    lm_fit <- lm(y~., data = data.frame(y = y_train, peak_compdata))
    
    peak_compdata_test <- predict(peak_pca_train, newdata = oes_peak_sc[test_index,])[, 1:5]
    colnames(peak_compdata_test) <- c(paste(rep("PC", 5), 1:5, sep = ""))
    # str(peak_compdata_test)
    # View(peak_compdata_test)
    
    a_pca <- predict(lm_fit, data.frame(peak_compdata_test))
    result.df.total.sc.tmp[2*(ii-1)+1,'PCR'] <- sqrt(mean((y_test-a_pca)^2))
    result.df.total.sc.tmp[2*(ii-1)+2,'PCR'] <- mean(abs(y_test - a_pca))
    
    ################################################################################
    # PC+LASSO
    ################################################################################
    pc_lasso = cv.glmnet(x=peak_compdata,
                         y=y_train)
    a_pc_lasso = predict(pc_lasso, newx=peak_compdata_test, s=pc_lasso$lambda.min) 
    result.df.total.sc.tmp[2*(ii-1)+1,'PCLASSO'] <- sqrt(mean((y_test-a_pc_lasso)^2))
    result.df.total.sc.tmp[2*(ii-1)+2,'PCLASSO'] <- mean(abs(y_test - a_pc_lasso))
    
    ################################################################################
    # PLS
    ################################################################################
    pls_fit = pls::plsr(y ~ ., ncomp=4,  data = data.frame(y = y_train, peak_compdata),
                        validation="LOO")
    a_pls = predict(pls_fit, ncomp=3, data.frame(peak_compdata_test))
    a_pls = a_pls[,1,1]
    result.df.total.sc.tmp[2*(ii-1)+1,'PLS'] <- sqrt(mean((y_test-a_pls)^2))
    result.df.total.sc.tmp[2*(ii-1)+2,'PLS'] <- mean(abs(y_test - a_pls))
    
    ################################################################################
    # Two-way model averaging
    ################################################################################
    
    # str(oes_penfsvd) # u1, u2 (each: 307 x time point(=9)) 
    # str(oes_basisfsvd) # v1, v2 (307 x 1024)
    
    ### PenFSVD + Bspline  ###
    # time point가 너무 적어서 basis 개수 조정이 어려움 -> error발생
    penfsvd_result = jk_bsp(
      X = c(oes_penfsvd$u, oes_penfsvd$v), 
      y_ = Y_sc, 
      x_nb = c(rep(5,ncomp), rep(50,ncomp)),
      b_nb = c(rep(5,ncomp), rep(50,ncomp)),
      tr = train_index, tst = test_index
    )
    
    result.df.total.sc.tmp[2*(ii-1)+1,'PenFSVD_Bsp'] <- penfsvd_result$rmse
    result.df.total.sc.tmp[2*(ii-1)+2,'PenFSVD_Bsp'] <- penfsvd_result$mae
    
    ### BasFSVD + Bspline  ###
    # time point가 너무 적어서 basis 개수 조정이 어려움 -> error발생
    basisfsvd_result = jk_bsp(
      X = c(oes_basisfsvd$u, oes_basisfsvd$v), 
      y_ = Y_sc, 
      x_nb = c(rep(5,ncomp), rep(50,ncomp)),
      b_nb = c(rep(5,ncomp), rep(50,ncomp)),
      tr = train_index, tst = test_index
    )
    
    result.df.total.sc.tmp[2*(ii-1)+1,'BasFSVD_Bsp'] <- basisfsvd_result$rmse
    result.df.total.sc.tmp[2*(ii-1)+2,'BasFSVD_Bsp'] <- basisfsvd_result$mae
    
    jma_oes_res <- list()
    jma_oes_res$PenBsp <- penfsvd_result
    jma_oes_res$BasBsp <- basisfsvd_result
    
    ### PenFSVD + Wv ###
    jma_oes_res$PenWv150 <- jk_wv(
      X = oes_penfsvd$v, y_ = Y_sc, q = 0.15,
      tr = train_index, tst = test_index
    )
    
    jma_oes_res$PenMerged150 <- jk_merge(
      error_mat1 = jma_oes_res$PenBsp$error_mat[, 1:ncomp],
      error_mat2 = jma_oes_res$PenWv150$error_mat,
      y_m1 = jma_oes_res$PenBsp$y_m[, 1:ncomp],
      y_m2 = jma_oes_res$PenWv150$y_m,
      y_ = Y_sc, tr = train_index, tst = test_index )
    
    result.df.total.sc.tmp[2*(ii-1)+1,'PenFSVD_Wv'] <- jma_oes_res$PenMerged150$rmse
    result.df.total.sc.tmp[2*(ii-1)+2,'PenFSVD_Wv'] <- jma_oes_res$PenMerged150$mae
    
    ### Bas + Wv ###
    jma_oes_res$BasWv150 <- jk_wv(
      X = oes_basisfsvd$v, y_ = Y_sc, q = 0.15,
      tr = train_index, tst = test_index
    )
    
    jma_oes_res$BasMerged150 <- jk_merge(
      error_mat1 = jma_oes_res$BasBsp$error_mat[, 1:ncomp], 
      error_mat2 = jma_oes_res$BasWv150$error_mat,
      y_m1 = jma_oes_res$BasBsp$y_m[, 1:ncomp], 
      y_m2 = jma_oes_res$BasWv150$y_m,
      y_ = Y_sc, tr = train_index, tst = test_index)
    
    result.df.total.sc.tmp[2*(ii-1)+1,'BasFSVD_Wv'] <- jma_oes_res$BasMerged150$rmse
    result.df.total.sc.tmp[2*(ii-1)+2,'BasFSVD_Wv'] <- jma_oes_res$BasMerged150$mae
    
    ################################################################################
    # GFLM
    ################################################################################
    #### Setup ####
    ncomp = 5 # number of component function/2
    nX = 1 # number of matrix X
    nt = 5  # number of time grid points on [0,1]
    ns = 1024 # number of wavelength grid points [0,1]
    
    ### grid for GFLM ###
    p = ncomp*nX
    Tps = vector("list", length=p)
    for (j in 1:p) {
      Tps[[j]] = (0:(nt-1))/(nt-1)
    }
    for (j in 1:p+p) {
      Tps[[j]] = (0:(ns-1))/(ns-1)
    }
    lambda = 10^seq(2,-3,by=-0.5)
    phi = 10^seq(3,-3,by=-1)
    nphi = length(phi)
    
    X_grpl <- list()
    for(i in 1:length(c(oes_basisfsvd$u, oes_basisfsvd$v))){
      X_grpl[[i]] <- (c(oes_basisfsvd$u, oes_basisfsvd$v)[[i]])[train_index,]
    }
    
    X_grpl_test <- list()
    for(i in 1:length(c(oes_basisfsvd$u, oes_basisfsvd$v))){
      X_grpl_test[[i]] <- (c(oes_basisfsvd$u, oes_basisfsvd$v)[[i]])[test_index,]
    }
    
    cv_grpl = cv.grplFlinear(k=5, Y=y_train, X=X_grpl, Tps=Tps, lambda=lambda, phi=phi,dfs = 40)
    cv_grpl_error = apply(cv_grpl, c(2,3), sum)
    minidx = which.min(cv_grpl_error)
    minphi_id = ifelse(minidx %% nphi !=0, minidx %% nphi, minidx %% nphi + nphi)
    minlam_id = ifelse(minidx %% nphi !=0, minidx %/% nphi + 1, minidx %/% nphi)
    grpl_fit = grplFlinear(Y=y_train, X=X_grpl, Tps=Tps, lambda=lambda[minlam_id], phi=phi[minphi_id], dfs=40)
    pred = rep(grpl_fit$intercept, length(test_index))
    for (jj in 1:p) {
      pred <- pred + (1/(nt-1))*X_grpl_test[[jj]] %*% grpl_fit$Coef[[jj]][,1]
    }
    for (jj in 1:p+p) {
      pred <- pred + (1/(ns-1))*X_grpl_test[[jj]] %*% grpl_fit$Coef[[jj]][,1]
    }
    
    pred_GFLM <- as.numeric(pred)
    result.df.total.sc.tmp[2*(ii-1)+1,'GFLM'] <- sqrt(mean((y_test-pred_GFLM)^2))
    result.df.total.sc.tmp[2*(ii-1)+2,'GFLM'] <- mean(abs(y_test-pred_GFLM))
    
    ################################################################################
    # PFR
    ################################################################################
    pfr_df = data.frame(y=y_train, 
                        u1=I(X_grpl[[1]]), u2=I(X_grpl[[2]]), u3=I(X_grpl[[3]]), u4=I(X_grpl[[4]]), u5=I(X_grpl[[5]]),
                        v1=I(X_grpl[[6]]), v2=I(X_grpl[[7]]), v3=I(X_grpl[[8]]), v4=I(X_grpl[[9]]), v5=I(X_grpl[[10]])) 
    
    pfr_fit = pfr(y ~ lf(u1, k=10, bs="ps")+lf(u2, k=10, bs="ps")+lf(u3, k=10, bs="ps") + lf(u4, k=10, bs="ps") + lf(u5, k=10, bs="ps") + 
                    lf(v1, k=10, bs="ps") + lf(v2, k=10, bs="ps")+lf(v3, k=10, bs="ps") + lf(v4, k=10, bs="ps") + lf(v5, k=10, bs="ps"), data=pfr_df)
    
    pfr_df_test = data.frame(y=y_test, 
                             u1=I(X_grpl_test[[1]]), u2=I(X_grpl_test[[2]]), u3=I(X_grpl_test[[3]]), u4=I(X_grpl_test[[4]]), u5=I(X_grpl_test[[5]]),
                             v1=I(X_grpl_test[[6]]), v2=I(X_grpl_test[[7]]), v3=I(X_grpl_test[[8]]), v4=I(X_grpl_test[[9]]), v5=I(X_grpl_test[[10]])) 
    
    pred_pfr = as.numeric(predict(pfr_fit, newdata=pfr_df_test))
    result.df.total.sc.tmp[2*(ii-1)+1,'PFR'] <- sqrt(mean((y_test-pred_pfr)^2))
    result.df.total.sc.tmp[2*(ii-1)+2,'PFR'] <- mean(abs(y_test-pred_pfr))
  }
  result.df.total.sc.list[[ttt]] <- result.df.total.sc.tmp
}


#####################################
#### Scenario IV (Wavelet basis) ####
#####################################
load("../Data/wavelet_true.RData")
result.df.total.wavelet.list <- list()

for(ttt in 1:20){
  print(paste("TOTAL iteration:", ttt))
  wavelet_simuldata.list.tmp <- list()
  x_err_std <- 3
  set.seed(40)
  for(i in 1:300){
    # wavelet_simuldata.list[[i]] <- t(t(wavelet_true)+rnorm(ncol(wavelet_true), 0, x_err_std))
    wavelet_simuldata.list.tmp[[i]] <- wavelet_true + matrix(rnorm(5*1024, 0, x_err_std), nrow=5, ncol=1024)
  }
  num_datasets <- 300
  X_total_wavelet.tmp <- array(NA, dim = c(num_datasets, 1024, 5))
  
  # Generate 300 different datasets
  for (i in 1:num_datasets) {
    X_total_wavelet.tmp[i, , ] <- t(wavelet_simuldata.list.tmp[[i]])
  }
  
  X_total_wavelet_time_wav.tmp <- array(NA, dim = c(num_datasets, 5, 1024))
  for (i in 1:num_datasets) {
    X_total_wavelet_time_wav.tmp[i, , ] <- wavelet_simuldata.list.tmp[[i]]
  }
  
  num_datasets <- dim(X_total_wavelet.tmp)[1]
  Y_wavelet <- numeric(num_datasets)
  Z_wavelet.tmp <- numeric(num_datasets)
  # Create time and wavelength axes with a small random noise added
  time_axis <- sin(2*pi/5*c(c(1:5)-1/3))
  wavelength_axis <- sin(2*pi/256*c(c(1:1024)-1/3))
  # Compute the outer product for the i-th dataset
  beta <- outer(time_axis, wavelength_axis, FUN = "*")
  
  for (i in 1:num_datasets) {
    Z_wavelet.tmp[i] <- sum((X_total_wavelet.tmp[i,,]-mean(X_total_wavelet.tmp[i,,])) * t(beta))
  }
  
  SNR_list <- c(2,5,10,20)
  y_err_std_list <- sqrt(var(Z_wavelet.tmp)/SNR_list)
  result.df.total.wavelet.tmp <- data.frame(row.names = paste(c("RMSE", "MAE"), "(SNR=", rep(SNR_list, each=2), ")", sep=""))
  
  ################################################################################
  # Two-way model averaging
  ################################################################################
  ncomp = 5 # SVD (u,v) pair 2개로 가정
  ## PenFSVD ##
  penfsvd = apply(X_total_wavelet_time_wav.tmp, 1, tw_svd_smu, ncomp=ncomp)
  
  ## basisFSVD ##
  basfsvd = apply(X_total_wavelet_time_wav.tmp, 1, basis_svd, nb_by = 1, ncmp=ncomp)
  
  nt <- dim(X_total_wavelet_time_wav.tmp)[2]
  nlambda <- dim(X_total_wavelet_time_wav.tmp)[3]
  
  ## sign adjustment (flip) ##
  pen_u1 = matrix(0, nrow=num_datasets, ncol=nt)
  for (i in 1:num_datasets) {
    pen_u1[i,] = penfsvd[[i]]$pen_u[1:nt,1]
  }
  for(i in which(pen_u1[,1]<0)){
    penfsvd[[i]]$pen_u[1:nt,1] <- - penfsvd[[i]]$pen_u[1:nt,1]
  }
  # matplot(t(pen_u1), type="l") # sign이 다르게 나오는 wafer가 있으면 조정해주기
  pen_u2 = matrix(0, nrow=num_datasets, ncol=nt)
  for (i in 1:num_datasets) {
    pen_u2[i,] = penfsvd[[i]]$pen_u[1:nt,2]
  }
  # matplot(t(pen_u2), type="l")
  for(i in which(pen_u2[,1]<0)){
    penfsvd[[i]]$pen_u[1:nt,2] <- - penfsvd[[i]]$pen_u[1:nt,2]
  }
  
  pen_u3 = matrix(0, nrow=num_datasets, ncol=nt)
  for (i in 1:num_datasets) {
    pen_u3[i,] = penfsvd[[i]]$pen_u[1:nt,3]
  }
  # matplot(t(pen_u3), type="l")
  for(i in which(pen_u3[,1]<0)){
    penfsvd[[i]]$pen_u[1:nt,3] <- - penfsvd[[i]]$pen_u[1:nt,3]
  }
  
  pen_u4 = matrix(0, nrow=num_datasets, ncol=nt)
  for (i in 1:num_datasets) {
    pen_u4[i,] = penfsvd[[i]]$pen_u[1:nt,4]
  }
  # matplot(t(pen_u4), type="l")
  for(i in which(pen_u4[,2]<0)){
    penfsvd[[i]]$pen_u[1:nt,4] <- - penfsvd[[i]]$pen_u[1:nt,4]
  }
  
  pen_u5 = matrix(0, nrow=num_datasets, ncol=nt)
  for (i in 1:num_datasets) {
    pen_u5[i,] = penfsvd[[i]]$pen_u[1:nt,5]
  }
  # matplot(t(pen_u5), type="l")
  for(i in which(pen_u5[,1]<0)){
    penfsvd[[i]]$pen_u[1:nt,5] <- - penfsvd[[i]]$pen_u[1:nt,5]
  }
  
  pen_v1 = matrix(0, nrow=num_datasets, ncol=nlambda)
  for (i in 1:num_datasets) {
    pen_v1[i,] = penfsvd[[i]]$pen_v[,1]
  }
  for(i in which(pen_v1[,1]<0)){
    penfsvd[[i]]$pen_v[1:nlambda,1] <- - penfsvd[[i]]$pen_v[1:nlambda,1]
  }
  # matplot(t(pen_v1), type="l") # sign이 다르게 나오는 wafer가 있으면 조정해주기
  pen_v2 = matrix(0, nrow=num_datasets, ncol=nlambda)
  for (i in 1:num_datasets) {
    pen_v2[i,] = penfsvd[[i]]$pen_v[,2]
  }
  # matplot(t(pen_v2), type="l")
  for(i in which(pen_v2[,1]<0)){
    penfsvd[[i]]$pen_v[1:nlambda,2] <- - penfsvd[[i]]$pen_v[1:nlambda,2]
  }
  
  pen_v3 = matrix(0, nrow=num_datasets, ncol=nlambda)
  for (i in 1:num_datasets) {
    pen_v3[i,] = penfsvd[[i]]$pen_v[,3]
  }
  # matplot(t(pen_v3), type="l")
  for(i in which(pen_v3[,1]<0)){
    penfsvd[[i]]$pen_v[1:nlambda,3] <- - penfsvd[[i]]$pen_v[1:nlambda,3]
  }
  # for(i in which(pen_v3[,15]<0)){
  #   penfsvd[[i]]$pen_v[,3] <- - penfsvd[[i]]$pen_v[,3]
  # }
  
  pen_v4 = matrix(0, nrow=num_datasets, ncol=nlambda)
  for (i in 1:num_datasets) {
    pen_v4[i,] = penfsvd[[i]]$pen_v[,4]
  }
  # matplot(t(pen_v4), type="l")
  for(i in which(pen_v4[,1]<0)){
    penfsvd[[i]]$pen_v[1:nlambda,4] <- - penfsvd[[i]]$pen_v[1:nlambda,4]
  }
  
  pen_v5 = matrix(0, nrow=num_datasets, ncol=nlambda)
  for (i in 1:num_datasets) {
    pen_v5[i,] = penfsvd[[i]]$pen_v[,5]
  }
  # matplot(t(pen_v5), type="l")
  for(i in which(pen_v5[,1]<0)){
    penfsvd[[i]]$pen_v[1:nlambda,5] <- - penfsvd[[i]]$pen_v[1:nlambda,5]
  }
  
  # for(i in which(pen_v5[,1]<0)){
  #   penfsvd[[i]]$pen_v[,5] <- - penfsvd[[i]]$pen_v[,5]
  # }
  
  # basis fsvd 결과도 필요하면 sign 조정
  
  ### join the result ###
  pen_u1 = do.call(rbind, lapply(penfsvd, function(x) x$pen_u[,1]))
  pen_u2 = do.call(rbind, lapply(penfsvd, function(x) x$pen_u[,2]))
  pen_v1 = do.call(rbind, lapply(penfsvd, function(x) x$pen_v[,1]))
  pen_v2 = do.call(rbind, lapply(penfsvd, function(x) x$pen_v[,2]))
  
  pen_u3 = do.call(rbind, lapply(penfsvd, function(x) x$pen_u[,3]))
  pen_u4 = do.call(rbind, lapply(penfsvd, function(x) x$pen_u[,4]))
  pen_v3 = do.call(rbind, lapply(penfsvd, function(x) x$pen_v[,3]))
  pen_v4 = do.call(rbind, lapply(penfsvd, function(x) x$pen_v[,4]))
  
  pen_u5 = do.call(rbind, lapply(penfsvd, function(x) x$pen_u[,5]))
  pen_v5 = do.call(rbind, lapply(penfsvd, function(x) x$pen_v[,5]))
  
  
  oes_penfsvd <- list(u = vector("list", length = ncomp),
                      v = vector("list", length = ncomp))
  oes_penfsvd$u[[1]] = pen_u1
  oes_penfsvd$u[[2]] = pen_u2
  oes_penfsvd$v[[1]] = pen_v1
  oes_penfsvd$v[[2]] = pen_v2
  
  oes_penfsvd$u[[3]] = pen_u3
  oes_penfsvd$u[[4]] = pen_u4
  oes_penfsvd$v[[3]] = pen_v3
  oes_penfsvd$v[[4]] = pen_v4
  
  oes_penfsvd$u[[5]] = pen_u5
  oes_penfsvd$v[[5]] = pen_v5
  
  bas_u1 = do.call(rbind, lapply(basfsvd, function(x) x$u[,1]))
  bas_u2 = do.call(rbind, lapply(basfsvd, function(x) x$u[,2]))
  bas_v1 = do.call(rbind, lapply(basfsvd, function(x) x$v[,1]))
  bas_v2 = do.call(rbind, lapply(basfsvd, function(x) x$v[,2]))
  
  bas_u3 = do.call(rbind, lapply(basfsvd, function(x) x$u[,3]))
  bas_u4 = do.call(rbind, lapply(basfsvd, function(x) x$u[,4]))
  bas_v3 = do.call(rbind, lapply(basfsvd, function(x) x$v[,3]))
  bas_v4 = do.call(rbind, lapply(basfsvd, function(x) x$v[,4]))
  
  bas_u5 = do.call(rbind, lapply(basfsvd, function(x) x$u[,5]))
  bas_v5 = do.call(rbind, lapply(basfsvd, function(x) x$v[,5]))
  
  oes_basisfsvd <- list(u = vector("list", length = ncomp),
                        v = vector("list", length = ncomp))
  oes_basisfsvd$u[[1]] = bas_u1
  oes_basisfsvd$u[[2]] = bas_u2
  oes_basisfsvd$v[[1]] = bas_v1
  oes_basisfsvd$v[[2]] = bas_v2
  
  oes_basisfsvd$u[[3]] = bas_u3
  oes_basisfsvd$u[[4]] = bas_u4
  oes_basisfsvd$v[[3]] = bas_v3
  oes_basisfsvd$v[[4]] = bas_v4
  
  oes_basisfsvd$u[[5]] = bas_u5
  oes_basisfsvd$v[[5]] = bas_v5
  
  for(ii in 1:length(y_err_std_list)){
    print(paste("iteration:", ii))
    y_err_std <- y_err_std_list[ii]
    set.seed(40)
    for (i in 1:num_datasets) {
      Y_wavelet[i] <- sum((X_total_wavelet.tmp[i,,]-mean(X_total_wavelet.tmp[i,,])) * t(beta)) + rnorm(1, mean = 0, sd = y_err_std)  # adding some noise
    }
    
    ################################################################################
    # MFPLS
    ################################################################################
    set.seed(ttt)
    train_index = sample(1:num_datasets, num_datasets*0.8) 
    test_index = (1:num_datasets)[-train_index]
    
    y_train <- Y_wavelet[train_index]
    y_test <- Y_wavelet[test_index]
    
    X_train <- X_total_wavelet.tmp[train_index,,]
    X_test <- X_total_wavelet.tmp[test_index,,]
    
    
    basis_num <- seq(30,60,15) ; basis_num
    Basis_name <- c('Basis 30','Basis 45','Basis 60')
    parameter.change.result.df <- data.frame(row.names = c('pls_num','mean_RMSE'))
    num_folds <- 10 # 10은 시간이 너무 많이 걸림 2로했을 때 약 28분 소요
    n <- dim(X_train)[1] 
    fold_indices <- createFolds(1:n, k = num_folds) 
    
    parameter.change.result.df_wavelet <- parameter_selection_MFPLS(X_train,y_train,basis_num,Basis_name,num_folds,fold_indices)
    i_th <- which.min(parameter.change.result.df_wavelet[2,]) # select argmin RMSE
    basis <- basis_num[i_th]
    pls_component <- parameter.change.result.df_wavelet[1,i_th]
    
    # see test errors
    MFPLS_model <- MFPLS(conv2fda(X_train),y_train, ncomp=pls_component,K=basis)
    prd <- predict.MFPLS(MFPLS_model, conv2fda(X_test))
    y_pred <- prd$y_pred
    # calculate accuracy
    result.df.total.wavelet.tmp[2*(ii-1)+1,'proposed'] <- rmse_score(y_pred,y_test)
    result.df.total.wavelet.tmp[2*(ii-1)+2,'proposed'] <- mae_score(y_pred,y_test)
    
    ################################################################################
    # PCR
    ################################################################################
    oes_peak_wavelet = get_peak_wv(X_total_wavelet_time_wav.tmp, span=3, npeaks=20)
    
    ### oes_peak에 PCR, PLS, LASSO 적용하면 됨 ###
    peak_pca_train <- princomp(oes_peak_wavelet[train_index,])
    
    peak_compdata <- peak_pca_train$scores[, 1:5]
    
    colnames(peak_compdata) <- c(paste(rep("PC", 5), 1:5, sep = ""))
    # peak_compdata %>% str() #dim: 250 * 45
    
    lm_fit <- lm(y~., data = data.frame(y = y_train, peak_compdata))
    
    peak_compdata_test <- predict(peak_pca_train, newdata = oes_peak_wavelet[test_index,])[, 1:5]
    colnames(peak_compdata_test) <- c(paste(rep("PC", 5), 1:5, sep = ""))
    # str(peak_compdata_test)
    # View(peak_compdata_test)
    
    a_pca <- predict(lm_fit, data.frame(peak_compdata_test))
    result.df.total.wavelet.tmp[2*(ii-1)+1,'PCR'] <- sqrt(mean((y_test-a_pca)^2))
    result.df.total.wavelet.tmp[2*(ii-1)+2,'PCR'] <- mean(abs(y_test - a_pca))
    
    ################################################################################
    # PC+LASSO
    ################################################################################
    pc_lasso = cv.glmnet(x=peak_compdata,
                         y=y_train)
    a_pc_lasso = predict(pc_lasso, newx=peak_compdata_test, s=pc_lasso$lambda.min) 
    result.df.total.wavelet.tmp[2*(ii-1)+1,'PCLASSO'] <- sqrt(mean((y_test-a_pc_lasso)^2))
    result.df.total.wavelet.tmp[2*(ii-1)+2,'PCLASSO'] <- mean(abs(y_test - a_pc_lasso))
    
    ################################################################################
    # PLS
    ################################################################################
    pls_fit = pls::plsr(y ~ ., ncomp=4,  data = data.frame(y = y_train, peak_compdata),
                        validation="LOO")
    a_pls = predict(pls_fit, ncomp=3, data.frame(peak_compdata_test))
    a_pls = a_pls[,1,1]
    result.df.total.wavelet.tmp[2*(ii-1)+1,'PLS'] <- sqrt(mean((y_test-a_pls)^2))
    result.df.total.wavelet.tmp[2*(ii-1)+2,'PLS'] <- mean(abs(y_test - a_pls))
    
    ################################################################################
    # Two-way model averaging
    ################################################################################
    
    # str(oes_penfsvd) # u1, u2 (each: 307 x time point(=9)) 
    # str(oes_basisfsvd) # v1, v2 (307 x 1024)
    
    ### PenFSVD + Bspline  ###
    # time point가 너무 적어서 basis 개수 조정이 어려움 -> error발생
    penfsvd_result = jk_bsp(
      X = c(oes_penfsvd$u, oes_penfsvd$v), 
      y_ = Y_wavelet, 
      x_nb = c(rep(5,ncomp), rep(50,ncomp)),
      b_nb = c(rep(5,ncomp), rep(50,ncomp)),
      tr = train_index, tst = test_index
    )
    
    result.df.total.wavelet.tmp[2*(ii-1)+1,'PenFSVD_Bsp'] <- penfsvd_result$rmse
    result.df.total.wavelet.tmp[2*(ii-1)+2,'PenFSVD_Bsp'] <- penfsvd_result$mae
    
    ### BasFSVD + Bspline  ###
    # time point가 너무 적어서 basis 개수 조정이 어려움 -> error발생
    basisfsvd_result = jk_bsp(
      X = c(oes_basisfsvd$u, oes_basisfsvd$v), 
      y_ = Y_wavelet, 
      x_nb = c(rep(5,ncomp), rep(50,ncomp)),
      b_nb = c(rep(5,ncomp), rep(50,ncomp)),
      tr = train_index, tst = test_index
    )
    
    result.df.total.wavelet.tmp[2*(ii-1)+1,'BasFSVD_Bsp'] <- basisfsvd_result$rmse
    result.df.total.wavelet.tmp[2*(ii-1)+2,'BasFSVD_Bsp'] <- basisfsvd_result$mae
    
    jma_oes_res <- list()
    jma_oes_res$PenBsp <- penfsvd_result
    jma_oes_res$BasBsp <- basisfsvd_result
    
    ### PenFSVD + Wv ###
    jma_oes_res$PenWv150 <- jk_wv(
      X = oes_penfsvd$v, y_ = Y_wavelet, q = 0.15,
      tr = train_index, tst = test_index
    )
    
    jma_oes_res$PenMerged150 <- jk_merge(
      error_mat1 = jma_oes_res$PenBsp$error_mat[, 1:ncomp],
      error_mat2 = jma_oes_res$PenWv150$error_mat,
      y_m1 = jma_oes_res$PenBsp$y_m[, 1:ncomp],
      y_m2 = jma_oes_res$PenWv150$y_m,
      y_ = Y_wavelet, tr = train_index, tst = test_index )
    
    result.df.total.wavelet.tmp[2*(ii-1)+1,'PenFSVD_Wv'] <- jma_oes_res$PenMerged150$rmse
    result.df.total.wavelet.tmp[2*(ii-1)+2,'PenFSVD_Wv'] <- jma_oes_res$PenMerged150$mae
    
    ### Bas + Wv ###
    jma_oes_res$BasWv150 <- jk_wv(
      X = oes_basisfsvd$v, y_ = Y_wavelet, q = 0.15,
      tr = train_index, tst = test_index
    )
    
    jma_oes_res$BasMerged150 <- jk_merge(
      error_mat1 = jma_oes_res$BasBsp$error_mat[, 1:ncomp], 
      error_mat2 = jma_oes_res$BasWv150$error_mat,
      y_m1 = jma_oes_res$BasBsp$y_m[, 1:ncomp], 
      y_m2 = jma_oes_res$BasWv150$y_m,
      y_ = Y_wavelet, tr = train_index, tst = test_index)
    
    result.df.total.wavelet.tmp[2*(ii-1)+1,'BasFSVD_Wv'] <- jma_oes_res$BasMerged150$rmse
    result.df.total.wavelet.tmp[2*(ii-1)+2,'BasFSVD_Wv'] <- jma_oes_res$BasMerged150$mae
    
    ################################################################################
    # GFLM
    ################################################################################
    #### Setup ####
    ncomp = 5 # number of component function/2
    nX = 1 # number of matrix X
    nt = 5  # number of time grid points on [0,1]
    ns = 1024 # number of wavelength grid points [0,1]
    
    ### grid for GFLM ###
    p = ncomp*nX
    Tps = vector("list", length=p)
    for (j in 1:p) {
      Tps[[j]] = (0:(nt-1))/(nt-1)
    }
    for (j in 1:p+p) {
      Tps[[j]] = (0:(ns-1))/(ns-1)
    }
    lambda = 10^seq(2,-3,by=-0.5)
    phi = 10^seq(3,-3,by=-1)
    nphi = length(phi)
    
    X_grpl <- list()
    for(i in 1:length(c(oes_basisfsvd$u, oes_basisfsvd$v))){
      X_grpl[[i]] <- (c(oes_basisfsvd$u, oes_basisfsvd$v)[[i]])[train_index,]
    }
    
    X_grpl_test <- list()
    for(i in 1:length(c(oes_basisfsvd$u, oes_basisfsvd$v))){
      X_grpl_test[[i]] <- (c(oes_basisfsvd$u, oes_basisfsvd$v)[[i]])[test_index,]
    }
    
    cv_grpl = cv.grplFlinear(k=5, Y=y_train, X=X_grpl, Tps=Tps, lambda=lambda, phi=phi,dfs = 40)
    cv_grpl_error = apply(cv_grpl, c(2,3), sum)
    minidx = which.min(cv_grpl_error)
    minphi_id = ifelse(minidx %% nphi !=0, minidx %% nphi, minidx %% nphi + nphi)
    minlam_id = ifelse(minidx %% nphi !=0, minidx %/% nphi + 1, minidx %/% nphi)
    grpl_fit = grplFlinear(Y=y_train, X=X_grpl, Tps=Tps, lambda=lambda[minlam_id], phi=phi[minphi_id], dfs=40)
    pred = rep(grpl_fit$intercept, length(test_index))
    for (jj in 1:p) {
      pred <- pred + (1/(nt-1))*X_grpl_test[[jj]] %*% grpl_fit$Coef[[jj]][,1]
    }
    for (jj in 1:p+p) {
      pred <- pred + (1/(ns-1))*X_grpl_test[[jj]] %*% grpl_fit$Coef[[jj]][,1]
    }
    
    pred_GFLM <- as.numeric(pred)
    result.df.total.wavelet.tmp[2*(ii-1)+1,'GFLM'] <- sqrt(mean((y_test-pred_GFLM)^2))
    result.df.total.wavelet.tmp[2*(ii-1)+2,'GFLM'] <- mean(abs(y_test-pred_GFLM))
    
    ################################################################################
    # PFR
    ################################################################################
    pfr_df = data.frame(y=y_train, 
                        u1=I(X_grpl[[1]]), u2=I(X_grpl[[2]]), u3=I(X_grpl[[3]]), u4=I(X_grpl[[4]]), u5=I(X_grpl[[5]]),
                        v1=I(X_grpl[[6]]), v2=I(X_grpl[[7]]), v3=I(X_grpl[[8]]), v4=I(X_grpl[[9]]), v5=I(X_grpl[[10]])) 
    
    pfr_fit = pfr(y ~ lf(u1, k=10, bs="ps")+lf(u2, k=10, bs="ps")+lf(u3, k=10, bs="ps") + lf(u4, k=10, bs="ps") + lf(u5, k=10, bs="ps") + 
                    lf(v1, k=10, bs="ps") + lf(v2, k=10, bs="ps")+lf(v3, k=10, bs="ps") + lf(v4, k=10, bs="ps") + lf(v5, k=10, bs="ps"), data=pfr_df)
    
    pfr_df_test = data.frame(y=y_test, 
                             u1=I(X_grpl_test[[1]]), u2=I(X_grpl_test[[2]]), u3=I(X_grpl_test[[3]]), u4=I(X_grpl_test[[4]]), u5=I(X_grpl_test[[5]]),
                             v1=I(X_grpl_test[[6]]), v2=I(X_grpl_test[[7]]), v3=I(X_grpl_test[[8]]), v4=I(X_grpl_test[[9]]), v5=I(X_grpl_test[[10]])) 
    
    pred_pfr = as.numeric(predict(pfr_fit, newdata=pfr_df_test))
    result.df.total.wavelet.tmp[2*(ii-1)+1,'PFR'] <- sqrt(mean((y_test-pred_pfr)^2))
    result.df.total.wavelet.tmp[2*(ii-1)+2,'PFR'] <- mean(abs(y_test-pred_pfr))
  }
  result.df.total.wavelet.list[[ttt]] <- result.df.total.wavelet.tmp
}

