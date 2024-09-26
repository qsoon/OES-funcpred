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
# stationary_true : 5 x 1024 matrix (underlying true X matrix)
stationary_simuldata.list <- list()
x_err_std <- 3
set.seed(10)
for(i in 1:303){
  # stationary_simuldata.list[[i]] <- t(t(stationary_true)+rnorm(ncol(stationary_true), 0, x_err_std))
  stationary_simuldata.list[[i]] <- stationary_true + matrix(rnorm(5*1024, 0, x_err_std), nrow=5, ncol=1024)
}
par(mfrow=c(1,2))
persp(1:5, 1:1024, stationary_simuldata.list[[1]], theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(stationary_simuldata.list[[1]]),
      cex.lab=3)
persp(1:5, 1:1024, stationary_simuldata.list[[2]], theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(stationary_simuldata.list[[2]]),
      cex.lab=3)


# outlier generate
# stationary
stationary_simuldata.list[[301]][3,] <- stationary_simuldata.list[[301]][3,] + max(stationary_simuldata.list[[301]]) # +로 주면 fica 도메인에서 outlier 위로 올라감
stationary_simuldata.list[[302]][,500] <- stationary_simuldata.list[[302]][,500] - 3*max(stationary_simuldata.list[[302]]) # +로 주면 fica 도메인에서 outlier 위로 올라감
stationary_simuldata.list[[303]][2:4,450:550] <- stationary_simuldata.list[[303]][2:4,450:550] + max(stationary_simuldata.list[[303]]) # +로 주면 fica 도메인에서 outlier 위로 올라감

par(mfrow=c(2,2))
persp(1:5, 1:1024, stationary_simuldata.list[[1]], theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(stationary_simuldata.list[[1]]),
      cex.lab=3)
persp(1:5, 1:1024, stationary_simuldata.list[[301]], theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(stationary_simuldata.list[[301]]),
      cex.lab=3)
persp(1:5, 1:1024, stationary_simuldata.list[[302]], theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(stationary_simuldata.list[[302]]),
      cex.lab=3)
persp(1:5, 1:1024, stationary_simuldata.list[[303]], theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(stationary_simuldata.list[[303]]),
      cex.lab=3)


# outlier detection using FICA (multivariate를 univariate으로 만들어야 돌릴 수 있음)
stationary.datamat.simul <- NULL
for(i in 1:303){
  stationary.datamat.simul <- cbind(stationary.datamat.simul, as.numeric(stationary_simuldata.list[[i]]))
}
arg <- 1:nrow(stationary.datamat.simul)

q <- 14
B <- create.bspline.basis(rangeval=c(min(arg),max(arg)), nbasis=q)
x.simul.stationary <- Data2fd(stationary.datamat.simul, argvals = arg, B)

Lfdobj <- int2Lfd(max(0, norder(B)-2))
penf.simul.stationary <- fdPar(B, Lfdobj, lambda=200) #penalty parameter (lambda)
ica.fd.simul.stationary <- ffobi(x.simul.stationary,  eigenfPar = penf.simul.stationary, w = "Cholesky")
## Plot by minimmum smoothed kurtosis
sc.simul.stationary <- ica.fd.simul.stationary$ICA.scores


# coloring
col.simul.stationary <- rep("black", 303)
col.simul.stationary[301] <- "red"
col.simul.stationary[302] <- "blue"
col.simul.stationary[303] <- "green"
# pairs(sc.simul1[,1:q], col=col.simul1, pch=20, cex=2, labels = paste("IC", 1:q), lower.panel=NULL)

# FICA domain plot (outlier identification)
par(mfrow=c(1,1))
plot(sc.simul.stationary[,2], sc.simul.stationary[,1], pch = 20, col=col.simul.stationary, cex=2,
     ylab = "IC 1", xlab = "IC 2")



######################################
#### Scenario II (non-stationary) ####
######################################

load("../Data/ns_stationary_true.RData")
# ns_stationary_true : 5 x 1024 matrix (underlying true X matrix)

ns_stationary_simuldata.list <- list()
x_err_std <- 3
set.seed(20)
for(i in 1:303){
  # ns_stationary_simuldata.list[[i]] <- t(t(ns_stationary_true)+rnorm(ncol(ns_stationary_true), 0, x_err_std))
  ns_stationary_simuldata.list[[i]] <- ns_stationary_true + matrix(rnorm(5*1024, 0, x_err_std), nrow=5, ncol=1024)
}
par(mfrow=c(1,2))
persp(1:5, 1:1024, ns_stationary_simuldata.list[[1]], theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(ns_stationary_simuldata.list[[1]]),
      cex.lab=3)
persp(1:5, 1:1024, ns_stationary_simuldata.list[[2]], theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(ns_stationary_simuldata.list[[2]]),
      cex.lab=3)


# outlier generate
ns_stationary_simuldata.list[[301]][3,] <- ns_stationary_simuldata.list[[301]][3,] + max(ns_stationary_simuldata.list[[301]]) # +로 주면 fica 도메인에서 outlier 위로 올라감
ns_stationary_simuldata.list[[302]][,500] <- ns_stationary_simuldata.list[[302]][,500] - 3*max(ns_stationary_simuldata.list[[302]]) # +로 주면 fica 도메인에서 outlier 위로 올라감
ns_stationary_simuldata.list[[303]][2:4,450:550] <- ns_stationary_simuldata.list[[303]][2:4,450:550] + max(ns_stationary_simuldata.list[[303]]) # +로 주면 fica 도메인에서 outlier 위로 올라감

par(mfrow=c(2,2))
persp(1:5, 1:1024, ns_stationary_simuldata.list[[1]], theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(ns_stationary_simuldata.list[[1]]),
      cex.lab=3)
persp(1:5, 1:1024, ns_stationary_simuldata.list[[301]], theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(ns_stationary_simuldata.list[[301]]),
      cex.lab=3)
persp(1:5, 1:1024, ns_stationary_simuldata.list[[302]], theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(ns_stationary_simuldata.list[[302]]),
      cex.lab=3)
persp(1:5, 1:1024, ns_stationary_simuldata.list[[303]], theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(ns_stationary_simuldata.list[[303]]),
      cex.lab=3)


ns_stationary.datamat.simul <- NULL
for(i in 1:303){
  ns_stationary.datamat.simul <- cbind(ns_stationary.datamat.simul, as.numeric(ns_stationary_simuldata.list[[i]]))
}
arg <- 1:nrow(ns_stationary.datamat.simul)

q <- 14
B <- create.bspline.basis(rangeval=c(min(arg),max(arg)), nbasis=q)
x.simul.ns_stationary <- Data2fd(ns_stationary.datamat.simul, argvals = arg, B)
#plot(x) #plot data
Lfdobj <- int2Lfd(max(0, norder(B)-2))
penf.simul.ns_stationary <- fdPar(B, Lfdobj, lambda=200) #penalty parameter (lambda)
ica.fd.simul.ns_stationary <- ffobi(x.simul.ns_stationary,  eigenfPar = penf.simul.ns_stationary, w = "Cholesky")
## Plot by minimmum smoothed kurtosis
sc.simul.ns_stationary <- ica.fd.simul.ns_stationary$ICA.scores



col.simul.ns_stationary <- rep("black", 303)
col.simul.ns_stationary[301] <- "red"
col.simul.ns_stationary[302] <- "blue"
col.simul.ns_stationary[303] <- "green"
# pairs(sc.simul1[,1:q], col=col.simul1, pch=20, cex=2, labels = paste("IC", 1:q), lower.panel=NULL)

# FICA domain plot (outlier identification)
par(mfrow=c(1,1))
plot(sc.simul.ns_stationary[,2], sc.simul.ns_stationary[,1], pch = 20, col=col.simul.ns_stationary, cex=2,
     ylab = "IC 1", xlab = "IC 2")


######################################
#### Scenario III (Fourier basis) ####
######################################

load("../Data/sc_true.RData")
# sc_true : 5 x 1024 matrix (underlying true X matrix)
sc_simuldata.list <- list()
x_err_std <- 3
set.seed(30)
for(i in 1:303){
  # stationary_simuldata.list[[i]] <- t(t(stationary_true)+rnorm(ncol(stationary_true), 0, x_err_std))
  sc_simuldata.list[[i]] <- sc_true + matrix(rnorm(5*1024, 0, x_err_std), nrow=5, ncol=1024)
}
par(mfrow=c(1,2))
persp(1:5, 1:1024, sc_simuldata.list[[1]], theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(sc_simuldata.list[[1]]),
      cex.lab=3)
persp(1:5, 1:1024, sc_simuldata.list[[2]], theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(sc_simuldata.list[[2]]),
      cex.lab=3)

# outlier generate
# stationary
sc_simuldata.list[[301]][3,] <- sc_simuldata.list[[301]][3,] + max(sc_simuldata.list[[301]]) # +로 주면 fica 도메인에서 outlier 위로 올라감
sc_simuldata.list[[302]][,500] <- sc_simuldata.list[[302]][,500] - 5*max(sc_simuldata.list[[302]]) # +로 주면 fica 도메인에서 outlier 위로 올라감
sc_simuldata.list[[303]][2:4,450:550] <- sc_simuldata.list[[303]][2:4,450:550] + max(sc_simuldata.list[[303]]) # +로 주면 fica 도메인에서 outlier 위로 올라감

par(mfrow=c(2,2))
persp(1:5, 1:1024, sc_simuldata.list[[1]], theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(sc_simuldata.list[[1]]),
      cex.lab=3)
persp(1:5, 1:1024, sc_simuldata.list[[301]], theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(sc_simuldata.list[[301]]),
      cex.lab=3)
persp(1:5, 1:1024, sc_simuldata.list[[302]], theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(sc_simuldata.list[[302]]),
      cex.lab=3)
persp(1:5, 1:1024, sc_simuldata.list[[303]], theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(sc_simuldata.list[[303]]),
      cex.lab=3)


# outlier detection using FICA (multivariate를 univariate으로 만들어야 돌릴 수 있음)
sc.datamat.simul <- NULL
for(i in 1:303){
  sc.datamat.simul <- cbind(sc.datamat.simul, as.numeric(sc_simuldata.list[[i]]))
}
arg <- 1:nrow(sc.datamat.simul)

q <- 14
B <- create.bspline.basis(rangeval=c(min(arg),max(arg)), nbasis=q)
x.simul.sc <- Data2fd(sc.datamat.simul, argvals = arg, B)

Lfdobj <- int2Lfd(max(0, norder(B)-2))
penf.simul.sc <- fdPar(B, Lfdobj, lambda=200) #penalty parameter (lambda)
ica.fd.simul.sc <- ffobi(x.simul.sc,  eigenfPar = penf.simul.sc, w = "Cholesky")
## Plot by minimmum smoothed kurtosis
sc.simul.sc <- ica.fd.simul.sc$ICA.scores


# coloring
col.simul.sc <- rep("black", 303)
col.simul.sc[301] <- "red"
col.simul.sc[302] <- "blue"
col.simul.sc[303] <- "green"
# pairs(sc.simul1[,1:q], col=col.simul1, pch=20, cex=2, labels = paste("IC", 1:q), lower.panel=NULL)

# FICA domain plot (outlier identification)
par(mfrow=c(1,1))
plot(sc.simul.sc[,2], sc.simul.sc[,1], pch = 20, col=col.simul.sc, cex=2,
     ylab = "IC 1", xlab = "IC 2")


#####################################
#### Scenario IV (Wavelet basis) ####
#####################################

load("../Data/wavelet_true.RData")
# wavelet_true : 5 x 1024 matrix (underlying true X matrix)
wavelet_simuldata.list <- list()
x_err_std <- 3
set.seed(40)
for(i in 1:303){
  # stationary_simuldata.list[[i]] <- t(t(stationary_true)+rnorm(ncol(stationary_true), 0, x_err_std))
  wavelet_simuldata.list[[i]] <- wavelet_true + matrix(rnorm(5*1024, 0, x_err_std), nrow=5, ncol=1024)
}
par(mfrow=c(1,2))
persp(1:5, 1:1024, wavelet_simuldata.list[[1]], theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(wavelet_simuldata.list[[1]]),
      cex.lab=3)
persp(1:5, 1:1024, wavelet_simuldata.list[[2]], theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(wavelet_simuldata.list[[2]]),
      cex.lab=3)

# outlier generate
# stationary
wavelet_simuldata.list[[301]][3,] <- wavelet_simuldata.list[[301]][3,] + max(wavelet_simuldata.list[[301]]) # +로 주면 fica 도메인에서 outlier 위로 올라감
wavelet_simuldata.list[[302]][,500] <- wavelet_simuldata.list[[302]][,500] - 3*max(wavelet_simuldata.list[[302]]) # +로 주면 fica 도메인에서 outlier 위로 올라감
wavelet_simuldata.list[[303]][2:4,450:550] <- wavelet_simuldata.list[[303]][2:4,450:550] + max(wavelet_simuldata.list[[303]]) # +로 주면 fica 도메인에서 outlier 위로 올라감

par(mfrow=c(2,2))
persp(1:5, 1:1024, wavelet_simuldata.list[[1]], theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(wavelet_simuldata.list[[1]]),
      cex.lab=3)
persp(1:5, 1:1024, wavelet_simuldata.list[[301]], theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(wavelet_simuldata.list[[301]]),
      cex.lab=3)
persp(1:5, 1:1024, wavelet_simuldata.list[[302]], theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(wavelet_simuldata.list[[302]]),
      cex.lab=3)
persp(1:5, 1:1024, wavelet_simuldata.list[[303]], theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(wavelet_simuldata.list[[303]]),
      cex.lab=3)


# outlier detection using FICA (multivariate를 univariate으로 만들어야 돌릴 수 있음)
wavelet.datamat.simul <- NULL
for(i in 1:303){
  wavelet.datamat.simul <- cbind(wavelet.datamat.simul, as.numeric(wavelet_simuldata.list[[i]]))
}
arg <- 1:nrow(wavelet.datamat.simul)

q <- 14
B <- create.bspline.basis(rangeval=c(min(arg),max(arg)), nbasis=q)
x.simul.wavelet <- Data2fd(wavelet.datamat.simul, argvals = arg, B)

Lfdobj <- int2Lfd(max(0, norder(B)-2))
penf.simul.wavelet <- fdPar(B, Lfdobj, lambda=200) #penalty parameter (lambda)
ica.fd.simul.wavelet <- ffobi(x.simul.wavelet,  eigenfPar = penf.simul.wavelet, w = "Cholesky")
## Plot by minimmum smoothed kurtosis
sc.simul.wavelet <- ica.fd.simul.wavelet$ICA.scores


# coloring
col.simul.wavelet <- rep("black", 303)
col.simul.wavelet[301] <- "red"
col.simul.wavelet[302] <- "blue"
col.simul.wavelet[303] <- "green"
# pairs(sc.simul1[,1:q], col=col.simul1, pch=20, cex=2, labels = paste("IC", 1:q), lower.panel=NULL)

# FICA domain plot (outlier identification)
par(mfrow=c(1,1))
plot(sc.simul.wavelet[,2], sc.simul.wavelet[,1], pch = 20, col=col.simul.wavelet, cex=2,
     ylab = "IC 1", xlab = "IC 2")
