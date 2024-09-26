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
# stationary
fitarima1 <- auto.arima(X_307_1024_13[1,,1], stationary=TRUE)

plot(X_307_1024_13[1,,1], type="l")
lines(fitted(fitarima1), col="red")

set.seed(10)
arimasim1 <- arima.sim(model = list(order = c(3,0,3), ar = c(2.2480033,-1.7549856,0.4987829),
                                    ma=c(0.3682579,-0.4684280,-0.3069292) ), n = 1024, innov=rnorm(1024,0,0.3))

set.seed(100)
arimasim2 <- arima.sim(model = list(order = c(5,0,3), ar = c(0.4167164,0.7269218,-0.1134622,-0.410328,0.2882262),
                                    ma=c(2.2354155,1.9175431,0.5824619)), n = 1024, innov=rnorm(1024,0,0.3))

set.seed(100)
arimasim3 <- arima.sim(model = list(order = c(4,0,3), ar = c(0.63477711,0.47947458,-0.06916927,-0.05711675),
                                    ma=c(2.1863652,1.8133879,0.5736847)), n = 1024, innov=rnorm(1024,0,0.3))

set.seed(1000)
arimasim4 <- arima.sim(model = list(order = c(3,0,1), ar = c(1.9940181,-1.4385340,0.4301553),
                                    ma=0.6203851), n = 1024, innov=rnorm(1024,0,0.3))

set.seed(1)
arimasim5 <- arima.sim(model = list(order = c(5,0,4), ar = c(1.3793139,-0.2306818,-0.1118268,-0.2145350,0.1670442),
                                    ma=c(1.6057971,0.7638062,-0.2105626,-0.2448167)), n = 1024, innov=rnorm(1024,0,0.3))

par(mfrow=c(2,3))
plot(-min(arimasim1)+arimasim1, type="l", main="ARIMA(3,0,3)")
plot(-min(arimasim2)+arimasim2, type="l", main="ARIMA(5,0,3)")
plot(-min(arimasim3)+arimasim3, type="l", main="ARIMA(4,0,3)")
plot(-min(arimasim4)+arimasim4, type="l", main="ARIMA(3,0,1)")
plot(-min(arimasim5)+arimasim5, type="l", main="ARIMA(5,0,4)")


# underlying X0 for Scenario I
stationary_true <- rbind(-min(arimasim1)+arimasim1,
                         -min(arimasim2)+arimasim2,
                         -min(arimasim3)+arimasim3,
                         -min(arimasim4)+arimasim4,
                         -min(arimasim5)+arimasim5)

persp(1:5, 1:1024, stationary_true, theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(stationary_true),
      cex.lab=3)

save(stationary_true, file="../Data/stationary_true.RData")


######################################
#### Scenario II (non-stationary) ####
######################################
# non-stationary

set.seed(1000)
ns_arimasim1 <- arima.sim(model = list(order = c(3,1,3), ar = c(1.6469007,-0.9761981,0.2360489),
                                       ma=c(-0.01821266,-0.63557982,-0.23570489)), n = 1024, innov=rnorm(1024,0,0.3))

set.seed(10000)
ns_arimasim2 <- arima.sim(model = list(order = c(3,1,2), ar = c(1.7593077,-1.1293089,0.2983178),
                                       ma=c(-0.1939220,-0.7359101)), n = 1024, innov=rnorm(1024,0,0.3))

set.seed(10)
ns_arimasim3 <- arima.sim(model = list(order = c(5,1,1), ar = c(1.31880572,-1.13581723,0.78435894,-0.46190033,0.09723474),
                                       ma=0.4810081), n = 1024, innov=rnorm(1024,0,0.3))

set.seed(100)
ns_arimasim4 <- arima.sim(model = list(order = c(2,1,1), ar = c(1.0214451,-0.4591397),
                                       ma=0.6122778), n = 1024, innov=rnorm(1024,0,0.3))

set.seed(1)
ns_arimasim5 <- arima.sim(model = list(order = c(4,1,1), ar = c(1.3686410,-1.0466829,0.5709685,-0.2390273),
                                       ma=0.6070021), n = 1024, innov=rnorm(1024,0,0.3))

par(mfrow=c(2,3))
plot(-min(ns_arimasim1[-1])+ns_arimasim1[-1], type="l", main="ARIMA(3,1,3)")
plot(-min(ns_arimasim2[-1])+ns_arimasim2[-1], type="l", main="ARIMA(3,1,2)")
plot(-min(ns_arimasim3[-1])+ns_arimasim3[-1], type="l", main="ARIMA(5,1,1)")
plot(-min(ns_arimasim4[-1])+ns_arimasim4[-1], type="l", main="ARIMA(2,1,1)")
plot(-min(ns_arimasim5[-1])+ns_arimasim5[-1], type="l", main="ARIMA(4,1,1)")

# non-stationary
ns_stationary_true <- rbind(-min(ns_arimasim1[-1])+ns_arimasim1[-1],
                            -min(ns_arimasim2[-1])+ns_arimasim2[-1],
                            -min(ns_arimasim3[-1])+ns_arimasim3[-1],
                            -min(ns_arimasim4[-1])+ns_arimasim4[-1],
                            -min(ns_arimasim5[-1])+ns_arimasim5[-1])

persp(1:5, 1:1024, ns_stationary_true, theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(ns_stationary_true),
      cex.lab=3)

save(ns_stationary_true, file="../Data/ns_stationary_true.RData")


######################################
#### Scenario III (Fourier basis) ####
######################################
# Fourier basis

x1 <- (15*sin(2*pi/400*c(1:1024-100))+7*sin(2*pi/100*c(1:1024-400)))
x2 <- (15*sin(2*pi/400*c(1:1024-100))+5*cos(2*pi/300*c(1:1024-500)))
x3 <- (5*sin(2*pi/100*c(1:1024-300))+7*cos(2*pi/150*c(1:1024-200)))
x4 <- (10*cos(2*pi/500*c(1:1024-100))+5*cos(2*pi/700*c(1:1024-200))) 
x5 <- (6*cos(2*pi/500*c(1:1024-100))+5*sin(2*pi/100*c(1:1024-400)))

par(mfrow=c(2,3))
plot(x1) ; plot(x2) ; plot(x3)
plot(x4) ; plot(x5)

sc_true <- rbind(x1,x2,x3,x4,x5)

par(mfrow=c(1,1))
persp(1:5, 1:1024, sc_true, theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(sc_true),
      cex.lab=3)

save(sc_true, file="sc_true.RData")


#####################################
#### Scenario IV (Wavelet basis) ####
#####################################
# wavelet basis

load("../Data/wvcoefs.RData")

x1_fit <- vector(length=1024)
for(i in 1:1024){
  x1_fit[i] <- sum(haar.scaling((i-1)/1024, 9, 0:511)* wvcoefs[[1]])
}
x1_fit <- x1_fit - min(x1_fit)

x2_fit <- vector(length=1024)
for(i in 1:1024){
  x2_fit[i] <- sum(haar.scaling((i-1)/1024, 9, 0:511)*  wvcoefs[[2]])
}
x2_fit <- x2_fit - min(x2_fit)

x3_fit <- vector(length=1024)
for(i in 1:1024){
  x3_fit[i] <- sum(haar.scaling((i-1)/1024, 9, 0:511)*  wvcoefs[[3]])
}
x3_fit <- x3_fit - min(x3_fit)

x4_fit <- vector(length=1024)
for(i in 1:1024){
  x4_fit[i] <- sum(haar.scaling((i-1)/1024, 9, 0:511)*  wvcoefs[[4]])
}
x4_fit <- x4_fit - min(x4_fit)

x5_fit <- vector(length=1024)
for(i in 1:1024){
  x5_fit[i] <- sum(haar.scaling((i-1)/1024, 9, 0:511)*  wvcoefs[[5]])
}
x5_fit <- x5_fit - min(x5_fit)

wavelet_true <- rbind(x1_fit, x2_fit, x3_fit, x4_fit, x5_fit)

par(mfrow=c(1,1))
persp(1:5, 1:1024, wavelet_true, theta = 40, phi = 20, expand = 0.5, col = "gray",
      xlab="time", ylab="wavelength", zlab="value", zlim=range(wavelet_true),
      cex.lab=3)

save(wavelet_true, file="wavelet_true.RData")
