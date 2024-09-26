# 1. pre defined function
#################################################################################
# (1) for EDA
#################################################################################

log.transformation<-function(data, finalrnum){
  log.data<-list()
  for (i in 1:finalrnum) {
    log.data[[i]] <- log(data[[i]])
  }
  return(log.data)
}
# time averaging function
timeaveraging<-function(data,finalrownum){
  cnum<-length(colMeans(data[[1]])); #1024
  result<-matrix(ncol=cnum)
  for (i in 1:finalrownum){
    result<-rbind(result,colMeans(data[[i]]))
  }
  result<-result[-1,]
  result
}
# change array data into funData object
conv2fda=function(data, arg=NULL){ 
  fdata=NULL
  if (length(dim(data))==2){
    fdata=funData( 1:dim(data)[2] ,data)
  }
  else {
    d<-dim(data)[3]
    for( j in 1:d ){
      fdata=c(fdata, funData( 1:dim(data)[2] , rbind(data[, ,j] ) ) )
    }
    return(multiFunData(fdata))
  }
}
# for burn-in removing , trim 25% of data from start time point
time25per.removing<- function(data,finalrownum){ #finalrownum = 307
  for (i in 1:finalrownum){
    data[[i]]<-data[[i]][(round(dim(data[[i]])[1]*0.25, digits = 0)+1):dim(data[[i]])[1],]
  }
  data
}
# for finding the location having important information, make split
time.half.spliting<-function(data,finalrownum){
  data1<-list();data2<-list();
  for (i in 1:finalrownum){
    data1[[i]]<-copy(data[[i]][1: round(dim(data[[i]])[1]*0.5, digits = 0)  ,])
    data2[[i]]<-copy(data[[i]][(round(dim(data[[i]])[1]*0.5, digits = 0)+1):dim(data[[i]])[1], ])
  }
  return(list(data1,data2))
}
time.quarter.spliting<-function(data,finalrownum){
  data1<-list();data2<-list();data3<-list();data4<-list();
  for (i in 1:finalrownum){
    data1[[i]]<-copy(data[[i]][1: round(dim(data[[i]])[1]*0.25, digits = 0)  ,])
    data2[[i]]<-copy(data[[i]][(round(dim(data[[i]])[1]*0.25, digits = 0)+1):(round(dim(data[[i]])[1]*0.5, digits = 0)), ])
    data3[[i]]<-copy(data[[i]][(round(dim(data[[i]])[1]*0.5, digits = 0)+1):(round(dim(data[[i]])[1]*0.75, digits = 0)),])
    data4[[i]]<-copy(data[[i]][(round(dim(data[[i]])[1]*0.75, digits = 0)+1):dim(data[[i]])[1], ])
  }
  return(list(data1,data2,data3,data4))
}

time.triple.spliting<-function(data,finalrownum){
  data1<-list();data2<-list();data3<-list();
  for (i in 1:finalrownum){
    data1[[i]]<-copy(data[[i]][1: round(dim(data[[i]])[1]*0.3333, digits = 0)  ,])
    data2[[i]]<-copy(data[[i]][(round(dim(data[[i]])[1]*0.3333, digits = 0)+1):(round(dim(data[[i]])[1]*0.6666, digits = 0)), ])
    data3[[i]]<-copy(data[[i]][(round(dim(data[[i]])[1]*0.6666, digits = 0)+1):dim(data[[i]])[1], ])
  }
  return(list(data1,data2,data3))
}

# for matching time length difference in some wafer. included time25per.removing 
time.25per.removing.endpoint.matching<-function(data,finalrownum,bwlength){
  for (i in 1:finalrownum){
    start<-(round(dim(data[[i]])[1]*0.25, digits = 0)+1)
    data[[i]]<-data[[i]][start:(start+bwlength-1),]
  }
  return(data)
}

#################################################################################
# (2) for Modeling (MFPLS)
#################################################################################
demean_fun<-function(func_X){
  m_X=meanFunction(func_X) #for each j=1...d, calculate colMean
  func_X=func_X-m_X
  func_X
}


prod_=function(x,y){
  ## Produit terme par terme! 
  if(dim(y)[1]==1){ 
    matrix(rep(x, ncol(y) ), byrow=F, ncol=dim(y)[2] )  *matrix(rep(y, nrow(x) ) , byrow = T, ncol=dim(y)[2])
  }else{ 
    x*y
  }
  
}
get_score=function(Mf_data, X_mean=NULL, uniExpansions, verbose=F ){ 
  if( is.null(X_mean) ){ 
    X_=Mf_data
  }else{ 
    X_=Mf_data-X_mean
    
  }
  d_f=length(Mf_data)
  if(verbose)message('base Projection ...')
  uniBasis <- mapply(function(expansion, data){do.call(univDecomp, c(list(funDataObject = data), expansion))},
                     expansion = uniExpansions, data = X_, SIMPLIFY = FALSE)
  
  allScores <- foreach::foreach(j = seq_len(d_f), .combine = "cbind")%do%{uniBasis[[j]]$scores}
  allScores
}



predict.MFPLS=function(out, X_test){ 
  u_E=list()
  d<-length(X_test)
  if(d==1 & attr(X_test, 'class')=='funData' )X_test<-multiFunData(X_test)
  
  
  for( u in 1:d){ 
    if(dimSupp(X_test)[u]==1){
      u_E= append(u_E, list(list(type = "splines1D", k = out$K[u]) ) )
    }else{ 
      u_E= append(u_E, list(list(type = "splines2D", k = c(out$K[u], out$K[u]), m=c(out$m_image, out$m_image) ) ) )
    }
  }
  
  X_s=score2fd(get_score(X_test ,out$Mean_X,u_E), u_E, arg_=out$arg_, K_j=get_K_j(out))
  list("y_pred"=scalarProduct(X_s,  get_b_PLS(out) )+ out$Mean_Y, "pls"=X_s)
}


score2fd=function(Scores, uniExpansions, arg_, K_j){ 
  d=length(uniExpansions)
  u_E=uniExpansions
  par=list()
  for( j in 1:d){
    ## Seulement pour B-spline 
    
    if(uniExpansions[[j]]$type=='splines1D') par=append(par, list(list(bs='ps', k = uniExpansions[[j]]$k, m= c(2, 2) )) )
    if(uniExpansions[[j]]$type=='splines2D') par=append(par, list(list(bs='ps', k = uniExpansions[[j]]$k, m= uniExpansions[[j]]$m )) )
    
    
  }
  
  
  if(dim(Scores)[1]>1){
    out_funct<- foreach::foreach(j = seq_len(length(u_E))) %do% {
      univExpansion(type = uniExpansions[[j]]$type,
                    scores = as.matrix(Scores[,K_j[[j]] ]),
                    argvals = arg_[[j]],
                    functions = uniExpansions[[j]]$functions,
                    params = par[[j]])
    }
  }else{ 
    out_funct<- foreach::foreach(j = seq_len(length(u_E))) %do% {
      univExpansion(type = uniExpansions[[j]]$type,
                    scores = matrix(Scores[ (uniExpansions[[j]]$k*(j-1)+1):(j*uniExpansions[[j]]$k)  ], nrow=1),
                    argvals = arg_[[j]],
                    functions = uniExpansions[[j]]$functions,
                    params = par[[j]])
    }
    
  }
  
  
  out_funct=multiFunData(out_funct)
  out_funct
}

get_K_j=function(obj_MFPLS){
  G=cumsum(lapply(lapply(obj_MFPLS$Scores, dim), function(x)x[2]))
  K_j=list()
  for( j in 1:length(G)  ){ 
    if(j>1){ 
      K_j=append(K_j, list(seq(G[j-1]+1, G[j])) )  
    }else{ 
      K_j=append(K_j, list(seq(G[j])) )  
    }
    
  }
  K_j
}

num_matrix=function(M){ 
  
  matrix(as.numeric(M), ncol=ncol(M))
}

## 1-4 steps :  Procedure!
unit_mfpls=function(alpha, Y, F_r, arg_){
  
  N=nrow(alpha[[1]]) ## Sample size
  d=length(alpha) ## Number of dimensions 
  
  ## Step _ 1 
  T_=list() ## Partial components 
  W_=list() ## Partial weight
  
  for( j in 1:d){ 
    
    out_j=oscorespls.fit(Y=Y,X=tcrossprod(alpha[[j]], F_r[[j]]) , ncomp=1, center=FALSE)
    t_j=out_j$scores
    T_=append(T_, list(t_j))
    w_j=out_j$projection
    W_=append(W_, list(w_j))
    
  }
  
  
  T_=matrix(unlist(T_), ncol=d) ## Concatenation of partial components : N*d 
  W_=matrix(unlist(W_), ncol=d)
  out=oscorespls.fit(Y=Y, X=T_, ncomp=1, center=FALSE) ## Ordinary PLS on Partial components 
  
  U=out$projection ## U Vector 
  
  #W_F=matrix(unlist(W_F), nrow=1)## True Weight 
  ## Step 2 
  
  
  W_F=W_*0
  for( j in 1:d){ 
    W_F[, j]= U[j]*W_[, j]#fpls to MFPLS using relationship
  }
  
  T_F=T_%*%U ## First component 
  #return(W_, U)
  
  
  # 
  # #Functional w^{(j)} for display purpose 
  # w_=list()
  # for( j in 1:d){ 
  #   if(length(arg_[[j]])==1){ 
  #   w_j=univExpansion(type='splines1D', matrix(W_F[, j], nrow =1), argvals = arg_[[j]], functions=NULL, params =  list(bs = "ps", m = 2, k = ncol(alpha[[j]]) ))
  #   }else{ 
  #     w_j=univExpansion(type='splines1D', matrix(W_F[, j], nrow =1), argvals = arg_[[j]], functions=NULL, params =  list(bs = "ps", m = 2, k = c(K,K)) )
  #     }
  #   
  #   w_=append(w_, w_j)
  #   
  # }
  # w_=multiFunData(w_ ) 
  # 
  
  ## Step 3-4  
  
  sd_1=sd(T_F)*sqrt(N-1) ## Var of First Components *(N-1), (N-1) is for having T^tT 
  
  alpha_1=list()
  # Compute P-Coefficients (Scores on Psi(t) is enough ) and residuals of X (alpha_1)    
  P=list()
  for( j in 1:d ){
    p_j=(1/(sd_1^2))*t(T_F)%*%alpha[[j]]
    a_1=alpha[[j]]-prod_(T_F, p_j)
    alpha_1=append(alpha_1, list(a_1) )
    P=append(P, list(p_j) )
  }
  
  # Compute C coefficient and residuals of Y ( Y_1)
  C=(1/(sd_1^2))*t(T_F)%*%Y 
  Y_1=Y-prod_(T_F, C)
  
  
  ## Return results 
  
  res=list()
  res$Y_1=Y_1
  res$alpha_1=alpha_1
  #res$w=w_
  res$W=matrix(W_F, ncol=1)
  res$P=matrix(unlist(P), ncol=1)  
  res$C=C
  res$T_F=matrix(T_F, ncol=1)
  res
}

## One MFPLS  Iteration 

Int_MFPLS<-function(X_alpha, Y_0, H, F_f, F_r, arg_, F_r_inv){ 
  Xi=list()
  P=list()
  W=list()
  C=list()
  for( i in 1:H){ 
    
    d_i=unit_mfpls(X_alpha, Y_0, F_r, arg_)
    
    Xi=append(Xi, list(d_i$T_F))
    C=append(C, list(d_i$C))
    P=append(P, list(d_i$P))
    W=append(W, list(d_i$W))
    
    ## Mise ? jour
    
    Y_0=d_i$Y_1
    X_alpha=d_i$alpha_1
    
    
  }
  
  
  Xi=list.cbind(Xi)
  P=list.cbind(P)
  W=F_r_inv%*%list.cbind(W) ## Correction des scores 
  
  P_2<-P; P_2[, H]<-0## The last coef doesn't contribute.
  P_W=t(W)%*%F_f%*%P_2
  rm(P_2)
  P_W[upper.tri(P_W, diag = T)]=0
  
  I_P=P_W+diag(x=1, nrow = nrow(P_W), ncol=ncol(P_W))
  
  Inv_P=solve(I_P)
  
  V=W%*%t(Inv_P)
  res=list()
  
  res$Xi=Xi ## Components
  res$C=C ## Coef_Y
  res$P=P ## Coef_P
  res$W=W ## Weight 
  res$V=V ## V
  res
}




MFPLS=function(X_t, Y, ncomps, K, m_image=1){ 
  
  d=length(X_t)
  if(d==1 & attr(X_t, 'class')=='funData' )X_t<-multiFunData(X_t)
  m_X=meanFunction(X_t)
  
  if(length(K)==1) K=rep(K, d) ## M?me nombre de base
  m_Y=mean(Y)
  
  X_t=X_t-m_X
  Y=Y-m_Y
  
  beta=list()
  F_r=list()
  F_=list()
  F_r_i=list()
  for( j in 1:d){
    
    decomp_j=univDecomp('splines1D', X_t[[j]], bs = "ps",m=2, k = K[j]) ## splines cubic 
    
    beta_j=decomp_j$scores 
    F_j=decomp_j$B
    F_r_j=expm::sqrtm(decomp_j$B)
    F_r_i=append(F_r_i, list(solve(F_r_j)))
    beta=append(beta, list(beta_j) )
    F_r=append(F_r,list(F_r_j) )
    F_=append(F_, list(F_j))
  }
  
  N=nrow(beta[[1]]) ## Sample size
  
  #d=sum(dimSupp(X_t)) ## Number of dimensions 
  #return(F_)
  ## Build F_f using F_$
  F_f=as.matrix(Matrix::bdiag(F_))
  
  #for(j in 1:d){
  #  K_j=dim(F_[[j]])[1]
  #  dim_j=(K_j*(j-1)+1):(K_j*j) 
  #  F_f[dim_j,dim_j ]=F_[[j]]
  
  #}
  
  #F_r_inv=matrix(0, ncol=K*d, nrow=K*d)
  F_r_inv=as.matrix(Matrix::bdiag(F_r_i))
  #for(j in 1:d){
  #  dim_j=(K*(j-1)+1):(K*j) 
  #  F_r_inv[dim_j,dim_j ]=solve(F_r[[j]])
  
  #}
  
  #return(beta)
  arg_=funData::argvals(X_t)
  res=Int_MFPLS(beta, Y,H=ncomps, F_f, F_r, arg_, F_r_inv)
  
  res$K=K
  res$ncomps=ncomps
  res$d=d
  res$arg_=arg_
  res$Mean_X=m_X
  res$Mean_Y=m_Y
  res$Scores=beta
  #res$F_f=F_f
  #res$F_r_inv=F_r_inv
  ## Corriger les scores 
  #res$W=res$F_r_inv%*%res$W
  #res$V=res$F_r_inv%*%res$V
  res$m_image=m_image
  class(res)<-'MFPLS'
  res
}



fun=function(Scores, arg_, d, K, K_j, m_image){
  
  Scores=as.matrix(Scores)
  f_t=list()
  for( j in 1:d){ 
    if(length(arg_[[j]])==1 ){ 
      f_t_j=univExpansion(type='splines1D', t(Scores[K_j[[j]], ]), argvals = arg_[[j]], functions=NULL, params =  list(bs = "ps", m = 2, k = K[j]))
    }else{ 
      f_t_j=univExpansion(type='splines2D', t(Scores[K_j[[j]], ]), argvals = arg_[[j]], functions=NULL, params =  list(bs = "ps", m = list(rep(m_image, m_image),rep(m_image, m_image)), k = c(K[j],K[j])  ))
    }
    
    f_t=append(f_t, f_t_j)
  }
  f_t=multiFunData(f_t)
  f_t
}



get_functions=function(obj_MFPLS, fun_name='w'){ 
  ## 
  
  d=obj_MFPLS$d  
  
  K=obj_MFPLS$K
  ds=dimSupp(obj_MFPLS$Mean_X) ## Format des des donn?es 
  arg_=obj_MFPLS$arg_
  K_j=get_K_j(obj_MFPLS)
  Scores=eval(parse(text=paste0('obj_MFPLS', '$', toupper(fun_name))))
  #return(Scores)
  # Correction par F_r 
  #Scores=obj_MFPLS$F_r_inv%*%Scores
  
  f_t=fun(Scores, arg_, d, K, K_j, obj_MFPLS$m_image)
  f_t
  
}


##  https://github.com/ClaraHapp/MFPCA/blob/master/R/univariateExpansions.R
expandBasisFunction <- function(scores, argvals = functions@argvals, functions){
  if(dim(scores)[2] != nObs(functions))
    stop("expandBasisFunction: number of scores for each observation and number of eigenfunctions does not match.")
  
  # collapse higher-dimensional functions, multiply with scores and resize the result
  d <- dim(functions@X)
  nd <- length(d)
  
  if(nd == 2)
    resX <- scores %*% functions@X
  
  if(nd == 3)
  {
    resX <- array(NA, dim = c(dim(scores)[1], d[-1]))
    
    for(i in seq_len(d[2]))
      resX[,i,] <- scores %*% functions@X[,i,]
  }
  
  if(nd == 4)
  {
    resX <- array(NA, dim = c(dim(scores)[1], d[-1]))
    
    for(i in seq_len(d[2]))
      for(j in seq_len(d[3]))
        resX[,i,j,] <- scores %*% functions@X[,i,j,]
  }
  
  if(nd > 4) # slow solution due to aperm
  {
    resX <- aperm(plyr::aaply(.data = functions@X, .margins = 3:nd, 
                              .fun = function(x,y){y %*% x}, y = scores), 
                  c(nd-1,nd, seq_len((nd-2))))
    dimnames(resX) <- NULL
  }  
  
  return( funData(argvals, resX) )
}


get_X_PLS=function(obj_MFPLS){ 
  
  p_t=get_functions(obj_MFPLS, fun_name = 'p')
  d=obj_MFPLS$d
  arg_=obj_MFPLS$arg_
  Scores=obj_MFPLS$Xi
  X_ncomp=list()
  for(j in 1:d){ 
    X_ncomp=append(X_ncomp, expandBasisFunction(Scores, argvals=arg_[[j]] , p_t[[j]]) )
    
  }
  
  multiFunData(X_ncomp)
}

get_b_PLS=function(obj_MFPLS){ 
  d=obj_MFPLS$d
  arg_=obj_MFPLS$arg_
  C=as.matrix(as.numeric(obj_MFPLS$C) ) 
  
  K_j=get_K_j(obj_MFPLS)
  
  Scores=obj_MFPLS$V%*%C
  
  b=fun(Scores, obj_MFPLS$arg_, obj_MFPLS$d, obj_MFPLS$K, K_j, obj_MFPLS$m_image)
  b
}




get_X_PLS=function(obj_MFPLS){ 
  
  v_t=get_functions(obj_MFPLS, fun_name = 'v')
  d=obj_MFPLS$d
  arg_=obj_MFPLS$arg_
  Scores=obj_MFPLS$Xi
  X_ncomp=list()
  for(j in 1:d){ 
    X_ncomp=append(expandBasisFunction(Scores, arg_[[j]], v_t[[j]]), X_ncomp )
    
  }
  
  multiFunData(X_ncomp)
  
}


cv_MFPLS=function(Y, X_t, Ncomps=1:(K_s-1), ncv=10, class=F, K_s=40, approx=T, par=T){ 
  
  N=nObs(X_t)#This functions returns the number of observations in a funData,
  k_folds<-caret::createFolds(1:N, ncv)
  error<-MLmetrics::RMSE
  pred_func<-function(x, y)predict(x, y)
  
  pred=c()
  for( ncomp in Ncomps){
    error_cv=c()
    for( k in 1:ncv){ 
      
      ind_train=-k_folds[[k]]
      out_train=MFPLS(X_t[ind_train], Y[ind_train], K=K_s, ncomps = ncomp)
      
      Y_hat_=pred_func(out_train, X_t[-ind_train])
      
      error_cv=c(error_cv, error( Y_hat_, Y[-ind_train]) )
      
      
    }
    pred=c(pred, mean(error_cv) )
  }
  
  
  if(approx) pred<-round(pred, 5)
  
  res=list(details=data.frame(Ncomps, pred), n_opt=Ncomps[which.min(pred)] )
  res
}


plot.MFPLS<-function(obj_mfpls, fun='b', ...){
  fun<-tolower(fun)
  if(!any(fun %in% c('b', 'x', 'w', 'p', 'v') ) ) stop( paste('function  "', fun, '" is not recognized') )
  
  if(fun=='b'){ 
    plot(get_b_PLS(obj_mfpls), ...)
  }else{ 
    if(fun=='x'){ 
      plot(get_X_PLS(obj_mfpls), ...)
    }else{ 
      plot(get_functions(obj_mfpls, fun), ...)
    }
    
  }
  
}

## accuracy score
mae_score<-MLmetrics::MAE
rmse_score<-MLmetrics::RMSE
r2_score<-MLmetrics::R2_Score

# function

Forward_selection_MFPLS <- function(p,n,basis_num,num_folds,fold_indices,X_train,y_train){
  
  # Step 3-5: Forward selection with 10-fold cross-validation
  selected_vars <- c()
  cum_vars<-c()
  added_var<-c()
  rmse_values <- vector()
  pls_component_values <- vector() 
  variable_results <- list()
  
  start_time<-Sys.time()
  for (i in 1:p) {
    rmse_cv <- rep(NA, num_folds)
    pls_component_cv <- rep(NA, num_folds)
    print('==================================')
    print('i = variables select') ; print(i)
    print('==================================')
    best_rmse_for_j<-Inf
    rmse_for_j<-vector()
    pls_for_j<-vector()
    for (j in 1:p) {
      if (!(j %in% cum_vars)) {
        selected_vars <- c(cum_vars, j)
        print(selected_vars)
        for (fold in 1:num_folds) {
          train_indices <- unlist(fold_indices[-fold])
          test_indices <- unlist(fold_indices[fold])
          
          train_X_fold <- conv2fda(X_train[train_indices,,selected_vars])
          test_X_fold <- conv2fda(X_train[test_indices,,selected_vars])
          
          train_y_fold <- y_train[train_indices]
          test_y_fold <- y_train[test_indices]
          
          # find best pls components
          error_1<-c();
          for (ncomp in 1:basis_num){
            # for computational aspects
            MFPLS_model<-MFPLS(train_X_fold, ncomp=ncomp, train_y_fold, K=basis_num)
            prd<-predict.MFPLS(MFPLS_model, test_X_fold)
            y_pred_1<-prd$y_pred;
            rmse_1 <-rmse_score(test_y_fold, y_pred_1)
            error_1<-c(error_1,rmse_1)
          }
          
          print('fold')
          print(fold)
          print(paste('best pls components : ' , which.min(error_1)))
          print(paste(' rmse : ' , error_1[which.min(error_1)]))
          
          rmse_cv[fold] <-error_1[which.min(error_1)]
          pls_component_cv[fold] <-which.min(error_1)
        }# 10-fold 끝남
        #10-fold 평균내었을 때 좋은 오차값,pls값 저장
        mean_rmse_cv<-mean(rmse_cv)
        if (mean_rmse_cv < best_rmse_for_j) { # j가 돎에따라 가장 좋은 변수를 added_var로 기억하는 과정
          print('=============renewed! =================')
          added_var<-j
          print(added_var)
          print('=======================================')
          best_rmse_for_j <- mean_rmse_cv
          pls_for_j<-mean(pls_component_cv)
          print('pls')
          print(pls_for_j)}
        
      }# if   문
      
    }# j가 다 돎.
    cum_vars<-c(cum_vars, added_var)
    mean_rmse <- best_rmse_for_j
    mean_pls_components <- pls_for_j
    print('mean_pls_compo')
    print(mean_pls_components)
    # Store the results
    rmse_values <- c(rmse_values, mean_rmse)
    pls_component_values <- c(pls_component_values, mean_pls_components)
    
    variable_results <-c(variable_results, list(var=cum_vars, rmse=mean_rmse, pls=mean_pls_components))
    print(variable_results)
    
    small_p<-length(variable_results)/3
    rmse_<-variable_results[seq(3,3*small_p,by=3)-1]
    if(small_p>1){
      if(rmse_[length(rmse_)-1]$rmse<rmse_[length(rmse_)]$rmse){break;break;break}
    }
    
  } 
  end_time<-Sys.time()
  consumed_time<-end_time-start_time
  print('consumed time for model selection is ... ')
  print(consumed_time)
  
  return(variable_results)
}

parameter_selection_MFPLS <- function(X_train,y_train,basis_num,Basis_name,num_folds,fold_indices){
  start_time<-Sys.time()
  for (a in 1:length(basis_num)){
    print('-----------------------------------------------')
    error_basis<-c();
    num_basis<-basis_num[a];
    basis_col<-Basis_name[a] ;print(basis_col)
    error_pls<-c();
    rmse_cv <- c();
    pls_component_cv <-c();
    for (fold in 1:num_folds){
      train_indices <- unlist(fold_indices[-fold])
      valid_indices <- unlist(fold_indices[fold])
      train_X_fold <- conv2fda(X_train[train_indices,,])
      valid_X_fold <- conv2fda(X_train[valid_indices,,])
      train_y_fold <- y_train[train_indices]
      valid_y_fold <- y_train[valid_indices]
      
      error_fold<-c()
      for (ncomp in 1:num_basis) {
        MFPLS_model<-MFPLS(train_X_fold, ncomp=ncomp, train_y_fold, K=num_basis)
        prd<-predict.MFPLS(MFPLS_model, valid_X_fold)
        y_pred_pls<-prd$y_pred;
        rmse_1 <-rmse_score(y_pred_pls, valid_y_fold)
        error_fold<-c(error_fold,rmse_1)
      }
      rmse_cv <-c(rmse_cv,error_fold[which.min(error_fold)])
      pls_component_cv <- c(pls_component_cv ,which.min(error_fold))
    }
    pls_for_basis<-mean(pls_component_cv)
    error_for_basis<- mean(rmse_cv)
    parameter.change.result.df[1,Basis_name[a]]<- round(pls_for_basis)
    parameter.change.result.df[2,Basis_name[a]]<- error_for_basis
    print('pls components')
    print(round(pls_for_basis))
  }
  end_time<-Sys.time()
  consumed_time<-end_time-start_time
  print('consumed time is ... ')
  print(consumed_time)
  return(parameter.change.result.df)
}


## Find peaks and calculate corresponding average and sd
get_peak_wv <- function(data, span = 5, npeaks = 20){
  nobs = dim(data)[1]
  nt = dim(data)[2]
  nwv = dim(data)[3] #1024
  compressed_data <- data.frame(matrix(NA, nrow = nobs, ncol = 2*npeaks))
  data_tmp <- array(dim = c(nobs, nt, npeaks))
  for(i in 1:nobs){
    for(t in 1:nt){
      peaks_allind <- ggpmisc:::find_peaks(data[i, t, ], span = span)
      data_tmp[i,t,] <- sort(data[i, t, peaks_allind], decreasing = T)[1:npeaks]
    }
    compressed_data[i,1:npeaks] <- apply(data_tmp[i,,], 2, mean, na.rm=TRUE)
    compressed_data[i,npeaks+1:npeaks] <- apply(data_tmp[i,,], 2, sd, na.rm=TRUE)
  }
  names(compressed_data) <- c(paste(rep("m", npeaks), 1:npeaks, sep = ""), 
                              paste(rep("sd", npeaks), 1:npeaks, sep = ""))
  return(compressed_data)
}


### two-way model avg
# PenFSVD
svd_reg_u2 <- function(data, s_u, i){
  X_tilde <- s_u %*% data   # smooth data once. no big difference empirically..
  tmp <- svd(X_tilde)
  return(list(u = tmp$u[,i], v = tmp$v[,i]))
}

tw_svd_smu <- function(data, ncomp = 5){
  plain_svd <- svd(data)
  res <- list(
    pen_u = matrix(0, nrow = nrow(data), ncol = ncomp),
    pen_v = matrix(0, nrow = ncol(data), ncol = ncomp),
    d = plain_svd$d)
  for(i in 1:ncomp){
    ssr_u <- assist::ssr(u ~ times,
                         data = data.frame(u = plain_svd$u[, i], 
                                           times = 1:nrow(data)),
                         scale = T, rk = cubic(times))
    s_u <- assist::hat.ssr(ssr_u) 
    tmp <- svd_reg_u2(data, s_u, i)
    res$pen_u[, i] <- tmp$u
    res$pen_v[, i] <- tmp$v
  }
  return(res)
}

orth_bu <- function(evalarg, rangeval = NULL, nb){
  if(is.null(rangeval)){
    rangeval <- c(evalarg[1], evalarg[length(evalarg)])
  }
  B <- getbasismatrix(
    evalarg, 
    create.bspline.basis(rangeval, nbasis = nb)
  )
  OB <- orthonormalization(B)[,1:nb]
  return(OB)
}

basis_svd <- function(data, nb_by = 5, ncmp = 5){
  #### data of size ntime * nwv(=1024)
  ntime = dim(data)[1]
  nwv = dim(data)[2]
  
  nb = round(ntime / nb_by)
  Bu = orth_bu(1:ntime, nb = nb)
  Bv = t(GenW(n = nwv)) #vanishing mmt = 10 by default
  
  X_tilde = t(Bu) %*% data %*% Bv
  
  ## rank one svd on X_tilde
  tmp <- svd(X_tilde)
  u = Bu %*% tmp$u[,1:ncmp]
  v = Bv %*% tmp$v[,1:ncmp]
  u_coef = tmp$u[,1:ncmp]
  v_coef = tmp$v[,1:ncmp]
  return(list(u = u, v = v, d = tmp$d,
              u_coef = u_coef, v_coef = v_coef))
}

#### Decomposition of X ####
gen_decompose = function(X, compfun, N, nt, ns, ncomp) {
  
  true_comps1 = compfun
  
  ##### PenFSVD on X_{i}(t_j, s_k) #####
  sim_tw_res1 = apply(X, 1, tw_svd_smu, ncomp=ncomp)
  #sim_tw_res1 = apply(X, 1, tw_svd2_smu, ncomp=ncomp)
  sim_penfsvd1 = list(u = vector("list", length = ncomp),
                      v = vector("list", length = ncomp))
  for (j in 1:ncomp) {
    for (i in 1:N) {
      sim_penfsvd1$u[[j]] = rbind(sim_penfsvd1$u[[j]], sim_tw_res1[[i]]$pen_u[,j])
      sim_penfsvd1$v[[j]] = rbind(sim_penfsvd1$v[[j]], sim_tw_res1[[i]]$pen_v[,j])
    }
  }  
  ### flip ###
  for (j in 1:ncomp) {
    for (i in 1:N) {
      sgn = sum(true_comps1$u[[j]][i,] %*% sim_penfsvd1$u[[j]][i,])
      if (sgn < 0) {
        sim_penfsvd1$u[[j]][i,] = (-1)*sim_penfsvd1$u[[j]][i,]
        sim_penfsvd1$v[[j]][i,] = (-1)*sim_penfsvd1$v[[j]][i,]
      } 
    }
  }
  
  ##### BasFSVD on X_{i}(t_j, s_k) #####
  sim_basis_res1 = apply(X, 1, basis_svd, nb_by=ncomp, ncmp=ncomp)
  sim_basfsvd1 = list(u = vector("list", length = ncomp),
                      v = vector("list", length = ncomp))
  
  for (j in 1:ncomp) {
    for (i in 1:N) {
      sim_basfsvd1$u[[j]] = rbind(sim_basfsvd1$u[[j]], sim_basis_res1[[i]]$u[,j])
      sim_basfsvd1$v[[j]] = rbind(sim_basfsvd1$v[[j]], sim_basis_res1[[i]]$v[,j])
    }
  }
  ### flip ###
  for (j in 1:ncomp) {
    for (i in 1:N) {
      sgn = sum(true_comps1$u[[j]][i,] %*% sim_basfsvd1$u[[j]][i,])
      if (sgn < 0) {
        sim_basfsvd1$u[[j]][i,] = (-1)*sim_basfsvd1$u[[j]][i,]
        sim_basfsvd1$v[[j]][i,] = (-1)*sim_basfsvd1$v[[j]][i,]
      } 
    }
  }
  
  ##### plain SVD on X_{i}(t_j, s_k) #####
  sim_plain_res1 = apply(X, 1, plain_svd, ncomp=ncomp)
  sim_plain1 = list(u = vector("list", length = ncomp),
                    v = vector("list", length = ncomp))
  
  for (j in 1:ncomp) {
    for (i in 1:N) {
      sim_plain1$u[[j]] = rbind(sim_plain1$u[[j]], sim_plain_res1[[i]]$u[,j])
      sim_plain1$v[[j]] = rbind(sim_plain1$v[[j]], sim_plain_res1[[i]]$v[,j])
    }
  }
  ### flip ###
  for (j in 1:ncomp) {
    for (i in 1:N) {
      sgn = sum(true_comps1$u[[j]][i,] %*% sim_plain1$u[[j]][i,])
      if (sgn < 0) {
        sim_plain1$u[[j]][i,] = (-1)*sim_plain1$u[[j]][i,]
        sim_plain1$v[[j]][i,] = (-1)*sim_plain1$v[[j]][i,]
      } 
    }
  }
  
  #### Estimation error for component function ####
  est_pen_u = est_bas_u = est_plain_u = c()
  est_pen_v = est_bas_v = est_plain_v = c()
  
  for (j in 1:ncomp) {
    est_pen_u[j] = mean(apply((true_comps1$u[[j]]-sim_penfsvd1$u[[j]]), 1, function(x) sum(x^2)*(1/(nt-1))))
    est_bas_u[j] = mean(apply((true_comps1$u[[j]]-sim_basfsvd1$u[[j]]), 1, function(x) sum(x^2)*(1/(nt-1))))
    est_plain_u[j] = mean(apply((true_comps1$u[[j]]-sim_plain1$u[[j]]), 1, function(x) sum(x^2)*(1/(nt-1))))
    
    est_pen_v[j] = mean(apply((true_comps1$v[[j]]-sim_penfsvd1$v[[j]]), 1, function(x) sum(x^2)*(1/(ns-1))))
    est_bas_v[j] = mean(apply((true_comps1$v[[j]]-sim_basfsvd1$v[[j]]), 1, function(x) sum(x^2)*(1/(ns-1))))
    est_plain_v[j] = mean(apply((true_comps1$v[[j]]-sim_plain1$v[[j]]), 1, function(x) sum(x^2)*(1/(ns-1))))
  }
  est_error=list(pen_u=est_pen_u, bas_u=est_bas_u, plain_u=est_plain_u,
                 pen_v=est_pen_v, bas_v=est_bas_v, plain_v=est_plain_v)
  
  res = list(true_cp=true_comps1, pen_cp=sim_penfsvd1, bas_cp=sim_basfsvd1, plain_cp=sim_plain1,
             ise=est_error)
  res
}

#### Bspline Basis ####
jk_bsp_errmat <- function(X, y_, x_nb, b_nb,
                          tr = train_ind, tst = test_ind, verbose = T){
  #### X : list of length p(=number of variables), each is nobs*argvals matrix
  #### y_: response variable
  #### x_nb, b_nb : nbasis for x and beta, respectively. each is vector of length p
  nobs <- length(y_)
  p <- length(X)
  
  # Find error matrix
  error_mat <- matrix(0, nrow = length(tr), ncol = p)
  for(j in 1:p){
    X_fd <- fdata(mdata = X[[j]])
    for(i in 1:length(tr)){
      jn_ind <- tr[setdiff(1:length(tr), i)] ### jackknife index
      X_fd_tmp <- X_fd[jn_ind]
      tmp_fit <- fregre.basis(
        fdataobj = X_fd_tmp, y = y_[jn_ind],
        basis.x = create.bspline.basis(rangeval = X_fd$rangeval, nbasis = x_nb[j]),
        basis.b = create.bspline.basis(rangeval = X_fd$rangeval, nbasis = b_nb[j]),
        lambda = 0)
      #print(tmp_fit)
      tmp_pred <- predict(tmp_fit, X_fd[tr[i]])
      #print(tmp_pred)
      # tmp_pred <- fit01(tmp_pred)
      error_mat[i, j] <- y_[tr[i]] - tmp_pred
    }
    if(verbose){
      cat("Model", j, "done\n")
    }
  }
  return(error_mat)
}

jk_bsp <- function(X, y_, x_nb, b_nb,
                   tr = train_ind, tst = test_ind, verbose = T){
  nobs <- length(y_)
  p <- length(X)
  
  res <- list()
  # Find error matrix
  res$error_mat <- jk_bsp_errmat(X = X, y_ = y_, x_nb = x_nb, b_nb = b_nb,
                                 tr = tr, tst = tst, verbose = verbose)
  # Find weight (QP)
  pd.err <- nearPD(t(res$error_mat) %*% res$error_mat / length(tr), doSym=TRUE)
  
  res$weight <- solve.QP(
    Dmat = as.matrix(pd.err$mat),
    dvec = rep(0, p),
    Amat = t(rbind(rep(1, p), diag(1, nrow = p, ncol = p))),
    bvec = c(1, rep(0, p)),
    meq = 1
  )
  
  res$y_m <- matrix(0, nrow = length(tst), ncol = p)
  for(j in 1:p){
    X_fd <- fdata(mdata = X[[j]])
    X_fd_tr <- X_fd[tr]
    tmp_fit <- fregre.basis(
      fdataobj = X_fd_tr, y = y_[tr],
      basis.x = create.bspline.basis(rangeval = X_fd$rangeval, nbasis = x_nb[j]),
      basis.b = create.bspline.basis(rangeval = X_fd$rangeval, nbasis = b_nb[j]),
      lambda = 0)
    tmp_pred <- predict(tmp_fit, X_fd[tst])
    # res$y_m[, j] <- fit01(tmp_pred)
    res$y_m[, j] <- tmp_pred
  }
  res$y_jpa <- res$y_m %*% res$weight$solution
  res$rmse <- sqrt(mean((res$y_jpa - y_[tst])^2))
  res$mae <- mean(abs(res$y_jpa - y_[tst]))
  res$r2 <- 1 - (sum((y_[tst] - res$y_jpa)^2)/sum((y_[tst] - mean(y_[tst]))^2))
  return(res)
}

fit01 <- function(x){
  #### Make y fit into [0,1]
  res <- x
  #print(sum(is.na(x)))
  if(any(x<0)){
    res[which(x<0)] <- 0
  }
  if(any(x>1)){
    res[which(x>1)] <- 1
  }
  return(res)
}

jk_wv_errmat <- function(X, y_, q = 0.1,
                         tr = train_ind, tst = test_ind, verbose = T){
  #### X : list of length p(=number of variables), each is nobs*argvals matrix
  #### y_: response variable
  #### x_nb, b_nb : nbasis for x and beta, respectively. each is vector of length p
  #### q: screening threshold
  nobs <- length(y_)
  p <- length(X)
  
  # Find error matrix
  error_mat <- matrix(0, nrow = length(tr), ncol = p)
  for(j in 1:p){
    tmp_x <- decomp1d(x = X[[j]], min.scale = 2)
    for(i in 1:length(tr)){
      jn_ind <- tr[setdiff(1:length(tr), i)] ### jackknife index
      tmp_scx <- scr_mag(tmp_x[jn_ind,], thres_q = q)
      tmp_fit <- tmp_fit <- lm(
        y~., data = data.frame(y = y_[jn_ind], tmp_scx$screened_x))
      tmp_pred <- predict(
        tmp_fit, newdata = data.frame(matrix(tmp_x[tr[i], tmp_scx$ind], nrow=1)))
      # tmp_pred <- fit01(tmp_pred)
      error_mat[i, j] <- y_[tr[i]] - tmp_pred
    }
    if(verbose){
      cat("Model", j, "done\n")
    }
  }
  return(error_mat)
}


jk_wv <- function(X, y_, q = 0.1,
                  tr = train_ind, tst = test_ind, verbose = T){
  #### X : list of length p(=number of variables), each is nobs*argvals matrix
  #### y_: response variable
  #### x_nb, b_nb : nbasis for x and beta, respectively. each is vector of length p
  #### q: screening threshold
  nobs <- length(y_)
  p <- length(X)
  
  res <- list()
  # Find error matrix
  res$error_mat <- jk_wv_errmat(X = X, y_ = y_, q = q,
                                tr = tr, tst = tst, verbose = verbose)
  pd.err <- nearPD(t(res$error_mat) %*% res$error_mat / length(tr), doSym=TRUE)
  # Find weight (QP)
  res$weight <- solve.QP(
    Dmat = as.matrix(pd.err$mat),
    dvec = rep(0, p),
    Amat = t(rbind(rep(1, p), diag(1, nrow = p, ncol = p))),
    bvec = c(1, rep(0, p)),
    meq = 1
  )
  
  res$y_m <- matrix(0, nrow = length(tst), ncol = p)
  for(j in 1:p){
    tmp_x <- decomp1d(x = X[[j]], min.scale = 2)
    tmp_scx <- scr_mag(tmp_x[tr,], thres_q = q)
    tmp_fit <- lm(y~., data = data.frame(y = y_[tr], tmp_scx$screened_x))
    tmp_pred <- predict(tmp_fit, newdata = data.frame(matrix(tmp_x[tst, tmp_scx$ind], nrow=length(tst))))
    # res$y_m[, j] <- fit01(tmp_pred)
    res$y_m[, j] <- tmp_pred
  }
  res$y_jpa <- res$y_m %*% res$weight$solution
  res$rmse <- sqrt(mean((res$y_jpa - y_[tst])^2))
  res$mae <- mean(abs(res$y_jpa - y_[tst]))
  res$r2 <- 1 - (sum((y_[tst] - res$y_jpa)^2)/sum((y_[tst] - mean(y_[tst]))^2))
  return(res)
}

decomp1d<-function(x, family="DaubLeAsymm", filter.number=4, 
                   bc="periodic", min.scale=2){
  nsample<-length(x[,1])
  N<-length(x[1,])
  wdt<-matrix(0, nrow=nsample, ncol=N)
  
  wds<-apply(x,1,wd, filter.number=filter.number, min.scale=min.scale, family=family, bc=bc)
  
  for(i in 1:nsample){
    base<-0
    for(j in ((log2(N)-1):min.scale)){
      temp<-accessD(wds[[i]], level=j, boundary=F)
      wdt[i,((base+1):(base+2^j))]<-temp
      base<-base+2^j
    }
    wdt[i,((N-2^min.scale+1):N)]<-accessC(wds[[i]], level=min.scale, boundary=F)
  }
  return(wdt)
}

###### Screening by mean Magnitude
scr_mag <- function(x, thres_q = 0.75, thres_var = 0.95, use_var = F){
  abs_mag <- abs(apply(x, 2, mean))
  x_var <- apply(x, 2, var)
  mag_order <- order(abs_mag, decreasing = T)
  var_cumsum <- cumsum(x_var[mag_order]) / sum(x_var)
  if(use_var){
    end_ord_ind <- which(var_cumsum >= thres_var)[1]
  } else {
    end_ord_ind <- round(dim(x)[2] * thres_q)
  }
  ind <- mag_order[1:end_ord_ind]
  var_explained <- var_cumsum[end_ord_ind]
  return(list(screened_x = x[, ind], ind = ind, var_explained = var_explained))
}


######### Merging ##########
#### Merging models ####
jk_merge <- function(error_mat1, error_mat2, y_m1, y_m2, 
                     y_, tr, tst){
  res <- list()
  res$error_mat <- cbind(error_mat1, error_mat2)
  p <- ncol(res$error_mat)
  
  pd.err <- nearPD(t(res$error_mat) %*% res$error_mat / length(tr), doSym=TRUE)
  # Find weight (QP)
  res$weight <- solve.QP(
    Dmat = as.matrix(pd.err$mat),
    dvec = rep(0, p),
    Amat = t(rbind(rep(1, p), diag(1, nrow = p, ncol = p))),
    bvec = c(1, rep(0, p)),
    meq = 1
  )
  
  res$y_m <- cbind(y_m1, y_m2)
  res$y_jpa <- res$y_m %*% res$weight$solution
  res$rmse <- sqrt(mean((res$y_jpa - y_[tst])^2))
  res$mae <- mean(abs(res$y_jpa - y_[tst]))
  res$r2 <- 1 - (sum((y_[tst] - res$y_jpa)^2)/sum((y_[tst] - mean(y_[tst]))^2))
  return(res)
}

# Multiple Functional Linear Model

grplFlinear <- function(Y, X, Tps, lambda, phi, dfs = 10,
                        adapt1 = NULL, adapt2 = NULL, ...){
  
  #### Observed functions x_{ij}(t) at possibly different support grids t
  ## Y = nx1 vector of responses
  ## X = list of matrices, each element corresponds to observations of one functional predictor.
  ##		Thus, [[jj]] would provide the jj-th observed function with subjects as rows and the support grid as columns
  ##		Note that each function (i.e. component of X) can have different support
  ## Tps = list of vectors, each element is a suppport grid of the corresponding observed function.
  ##		Thus, Tps[[jj]] would give the support grid of the jj-th observed function. This grid is
  ##		assumed to be the same for all subjects ii, and equidistant
  ## lambda = vector of penalty parameters
  ## phi = penalty parameter for smoothing
  ## dfs = dfs used for basis expansions of coefficient functions (can be a vector)
  ## adapt1, adapt2 = adaptive weights for the coefficient functions and second derivatives, respectively.
  ##     (if given, vectors with length matching the number of predictors)
  
  
  nsub = length(Y) ## number of subjects
  nfunc = length(Tps) ## number of functions per subject
  
  
  
  #### We use bsplines as basis functions for the corrsponding beta functions
  if (length(dfs) == 1)
    dfs = rep(dfs, nfunc) ## vector of intended df of each spline basis
  if (length(dfs) != nfunc)
    stop("length of dfs does not match number of predictors")
  
  B <- Psi <- Omega <- K <- iR <- eK <- list()
  delt <- rep(NA, nfunc)
  
  
  for (jj in 1:nfunc){
    spj = diff(range(Tps[[jj]]))#/(dfs[jj]-2)
    bknj = c(min(Tps[[jj]]) - spj, max(Tps[[jj]]) + spj) ## boundary knots
    B[[jj]] = bs(Tps[[jj]], df=dfs[jj], Boundary.knots=bknj) ## basis spline set up
    delt[jj] = Tps[[jj]][2] - Tps[[jj]][1] ## differences in Tps
    
    Psi[[jj]] = delt[jj] * t(B[[jj]]) %*% B[[jj]] ## approximate norm of bsplines assuming dense design
    if (length(adapt1) == nfunc)
      Psi[[jj]] = adapt1[jj]*Psi[[jj]]
    
    dBj <- matrix(NA,nrow(B[[jj]]),ncol(B[[jj]]))
    for (k in 1:ncol(B[[jj]])) ## computation of second derivatives
    {
      iS <- interpSpline(Tps[[jj]],B[[jj]][,k])
      dBj[,k] <- predict(iS, Tps[[jj]], deriv = 2)$y
    }
    Omega[[jj]] = delt[jj] * t(dBj) %*% dBj ## approximate norm of 2nd deriv of bsplines assuming dense design
    if (length(adapt2) == nfunc)
      Omega[[jj]] = adapt2[jj]*Omega[[jj]]
    
    K[[jj]] = Psi[[jj]] + phi * Omega[[jj]] ## K matrix
    K[[jj]] = as.matrix(nearPD(K[[jj]], doSym=TRUE)$mat)
    eK[[jj]] <- eigen(K[[jj]])
    #iR[[jj]] <- t((1/sqrt(eK[[jj]]$values))*t(eK[[jj]]$vectors))
    iR[[jj]] = backsolve(chol(K[[jj]]), x = diag(ncol(K[[jj]])))  ## inverse of cholesky of K
  }
  
  
  ## covariates for the linear model
  Z = 1
  for (jj in 1:nfunc)
  {
    tmp = delt[jj]*(X[[jj]]%*%B[[jj]])
    Z = cbind(Z, tmp%*%iR[[jj]])
  }
  
  
  ## group lasso
  index = c(NA,rep(1:nfunc,dfs))
  grpl = grplasso(x = Z, y = Y, index = index, model = LinReg(), lambda = lambda,
                  standardize = F, control=grpl.control(trace=0), ...)
  
  
  ## output: intercept and fitted coefficient functions
  intercept = grpl$coef[1,]
  Coef <- list()
  index[1] = 0
  for (jj in 1:nfunc)
  {
    Coef[[jj]] <- B[[jj]]%*%iR[[jj]]%*%grpl$coef[index == jj,]
  }
  
  out = list("intercept" = intercept, "Coef" = Coef)
  return(out)
}



cv.grplFlinear <- function(k, Y, X, Tps, lambda, phi, dfs = 10, ...)
{
  # K-fold cross-validation function
  n <- length(Y)
  p <- length(X)
  os <- sample(n,n)
  
  cvError <- array(NA, c(k,length(phi),length(lambda)))
  nk <- floor(n/k)
  for (wk in 1:k)
  {
    cat("fold",wk,"\n")
    if (wk < k)
      inds <- os[(wk-1)*nk+(1:nk)]
    else
      inds <- os[-(1:((wk-1)*nk))]
    
    # training/test obs
    Xk <- X
    Xkt <- X
    delt = rep(NA, p)
    for (jj in 1:p)
    {
      Xk[[jj]] <- X[[jj]][-inds,]
      Xkt[[jj]] <- X[[jj]][inds,]
      delt[jj] = Tps[[jj]][2] - Tps[[jj]][1]
    }
    Yk <- Y[-inds]
    Ykt <- Y[inds]
    
    # phi walues
    for (wp in 1:length(phi))
    {
      # model fitting
      grplk <- grplFlinear(Y = Yk, X = Xk, Tps = Tps, lambda = lambda, 
                           phi = phi[wp], dfs = dfs, ...)
      
      # test set prediction
      predk <- matrix(NA,length(inds),length(lambda))
      for (ll in 1:length(lambda))
      {
        predll <- rep(grplk$intercept[ll],nrow(predk))
        for (jj in 1:p)
        {
          predll <- predll + delt[jj]*Xkt[[jj]]%*%grplk$Coef[[jj]][,ll]
        }
        # error
        cvError[wk,wp,ll] <- sum((Ykt - predll)^2)
      }
    }
  }
  return(cvError)
}

mother.wavelet.haar <- function(x){
  res <- vector(length=length(x))
  for(i in 1:length(x)){
    if((0 <= x[i]) & (x[i]<1/2)){
      res[i] <- 1
    }
    else if((1/2 <= x[i]) & (x[i]<1)){
      res[i] <- -1
    }
    else res[i] <- 0
  }
  return(res)
}

haar.wavelet <- function(x, j, k){
  return(2^(j/2)*mother.wavelet.haar(2^j*x-k))
}

father.wavelet.haar <- function(x){
  res <- vector(length=length(x))
  for(i in 1:length(x)){
    if((0 <= x[i]) & (x[i]<1)){
      res[i] <- 1
    }
    else res[i] <- 0
  }
  return(res)
}

haar.scaling <- function(x, j, k){
  return(2^(j/2)*father.wavelet.haar(2^j*x-k))
}
