# Compare models

library(BGLR)
library(reshape2)
library(dplyr)
library(R2jags)
library(ggplot2)

# Check the code for the other models here:
# https://github.com/gdlc/BGLR-R/blob/master/inst/md/GxE_usingInteractions.md

remove(list = ls())

# gmBAMMIT model

gmBammit <- function(data, G, hParam){

  N <- nrow(data)
  Y <- data$y
  B1 <- length(levels(data$var1))
  B2 <- length(levels(data$var2))
  var1 <- data$var1
  var2 <- data$var2

  mb1 <- hParam$mb1
  mmu <- hParam$mmu
  smu <- hParam$smu
  a <- hParam$a
  b <- hParam$b
  Q <- hParam$Q
  stheta <-  hParam$stheta

  n.thin <- hParam$nThin
  n.burnin <- hParam$burnIn


  modelCode <- "
    model{
      # Likelihood
      for (n in 1:N) {
          Y[n] ~ dnorm(mu[n], sy)
          mu[n] = muall + b1[var1[n]] + b2[var2[n]] + int[n]
          int[n] = sum(lambda[1:Q] * beta1[var1[n],1:Q] * beta2[var2[n],1:Q])
      }

      # Priors

      muall ~ dnorm(mmu, smu^-2) # grand mean

      #for(i in 1:B1) {
          b1 ~ dmnorm(mb1[1:B1], sb1[1:B1, 1:B1]) # Prior on effect 1
      #}

      for(i in 1:B2) {
          b2[i] ~ dnorm(0, sb2^-2) # Prior on effect 2
      }

      for(q in 1:Q){
          for(i in 1:B1){
            thetaBeta1[i,q] ~ dnorm(0,sb1[i,i])
          }
          mBeta1[q] = sum(thetaBeta1[1:B1,q])/B1

          for(i in 1:B1){
            thetaBeta1New[i,q] = thetaBeta1[i,q] - mBeta1[q]
          }

          sqrtThetaBeta1[q] = sqrt(1/(sum(thetaBeta1New[1:B1,q]^2 + 0.000001)))
          for(i in 1:B1){
            beta1[i,q] = thetaBeta1New[i,q]*sqrtThetaBeta1[q]
          }
      }


       for(q in 1:Q){
          for(i in 1:B2){
            thetaBeta2[i,q] ~ dnorm(0,stheta)
          }
          mBeta2[q] = sum(thetaBeta2[1:B2,q])/B2

          for(i in 1:B2){
            thetaBeta2New[i,q] = thetaBeta2[i,q] - mBeta2[q]
          }

          sqrtThetaBeta2[q] = sqrt(1/(sum(thetaBeta2New[1:B2,q]^2 + 0.000001)))
          for(i in 1:B2){
            beta2[i,q] = thetaBeta2New[i,q]*sqrtThetaBeta2[q]
          }
      }

      # Prior on eigenvalues
      for(q in 1:Q) {
        lambda_raw[q] ~ dnorm(0, slambda^-2)T(0,)
      }
      lambda = sort(lambda_raw)

      #Priors
      sb2 ~ dt(0, 1, 1)T(0,)
      slambda ~ dt(0, 1, 1)T(0,)

      sy ~ dgamma(a, b) # Prior on residual standard deviation - inverse of sy

    }
    "
  # Set up the data
  modelData <- list(
    N = N,
    Y = Y,
    B1 = B1,
    B2 = B2,
    Q = Q,
    var1 = var1,
    var2 = var2,
    sb1 = sb1,
    mmu = mmu,
    smu = smu,
    stheta = stheta,
    a = a,
    b = b,
    mb1 = mb1
  )

  # Choose the parameters to watch
  modelParameters <- c(
    "b1", "b2",  "lambda", "beta1", "beta2", "sy", "muall", "int", "mu"
  )

  # Run the model
  modelRun <- jags(
    data = modelData,
    parameters.to.save = modelParameters,
    model.file = textConnection(modelCode),
    #progress.bar = "none",
    n.thin = n.thin,
    n.burnin = n.burnin,
    n.iter = n.thin * n.burnin + 2000
  )

  return(modelRun)
}

# Model 1 (M1): Environment and main additive genomic effects

m1_jags <- function(data, G, hParam) {

  N <- nrow(data)
  Y <- data$y
  B1 <- length(levels(data$var1))
  B2 <- length(levels(data$var2))
  var1 <- data$var1
  var2 <- data$var2

  mb1 <- hParam$mb1
  mmu <- hParam$mmu
  smu <- hParam$smu
  a <- hParam$a
  b <- hParam$b

  n.thin <- hParam$nThin
  n.burnin <- hParam$burnIn

  modelCode <- "
    model{
      # Likelihood
      for (n in 1:N) {
          Y[n] ~ dnorm(mu[n], sy)
          mu[n] = muall + b1[var1[n]] + b2[var2[n]]
      }

      # Priors

      muall ~ dnorm(mmu, smu^-2) # grand mean
      b1 ~ dmnorm(mb1[1:B1], G[1:B1, 1:B1]*(sb1^-2)) # Prior on genotype

      for(i in 1:B2) {
          b2[i] ~ dnorm(0, sb2^-2) # Prior on effect 2
      }

      #Priors
      sb1 ~ dt(0, 1, 1)T(0,)
      sb2 ~ dt(0, 1, 1)T(0,)
      slambda ~ dt(0, 1, 1)T(0,)

      sy ~ dgamma(a, b) # Prior on residual standard deviation - inverse of sy
    }
    "

  # Set up the data
  modelData <- list(
    N = N,
    Y = Y,
    B1 = B1,
    B2 = B2,
    G = G,
    var1 = var1,
    var2 = var2,
    mb1 = mb1,
    mmu = mmu,
    smu = smu,
    a = a,
    b = b
  )

  # Choose the parameters to watch
  modelParameters <- c(
    "b1", "b2", "sy", "muall", "mu"
  )

  # Run the model
  modelRun <- jags(
    data = modelData,
    parameters.to.save = modelParameters,
    model.file = textConnection(modelCode),
    # progress.bar = "none",
    n.thin = n.thin,
    n.burnin = n.burnin,
    n.iter = n.thin * n.burnin + 2000
  )

  return(modelRun)
}

# Model 2(M2): As a random regression on markers

m2 <- function(Y = wheat.Y,
               M = wheat.X,
               nIter = 3000,
               burnIn = 1000) {

  y2 <- Y[, 2]
  y3 <- Y[, 3]
  y <- c(y2, y3)

  X <- scale(M) / sqrt(ncol(M))
  X0 <- matrix(nrow = nrow(X), ncol = ncol(X), 0) # a matrix full of zeros
  X_main <- rbind(X, X)
  X_1 <- rbind(X, X0)
  X_2 <- rbind(X0, X)
  # X_3 <- rbind(X0, X0, X, X0)
  # X_4 <- rbind(X0, X0, X0, X)

  fm <- BGLR(y = y, ETA = list(
    main = list(X = X_main, model = "BRR"),
    int1 = list(X = X_1, model = "BRR"),
    int2 = list(X = X_2, model = "BRR")),
    # int3 = list(X = X_3, model = "BRR"),
    # int4 = list(X = X_4, model = "BRR")),
    nIter = nIter, burnIn = burnIn, saveAt = "GxE_", groups = rep(1:2, each = nrow(X)))

  return(fm)
}


# Model 3(M3): Using genomic relationships
m3 <- function(Y = wheat.Y,
               M = wheat.X,
               G = BGData::getG(wheat.X, center = TRUE, scaleG = TRUE, scale = TRUE),
               nIter = 3000,
               burnIn = 1000) {
  y2 <- Y[, 2]
  y3 <- Y[, 3]
  y <- c(y2, y3)

  EVD <- eigen(G)
  PC <- EVD$vectors[, EVD$values > 1e-5]
  for (i in 1:ncol(PC)) {
    PC[, i] <- EVD$vectors[, i] * sqrt(EVD$values[i])
  }

  X <- scale(M) / sqrt(ncol(M))
  XMain <- rbind(PC, PC)
  X0 <- matrix(nrow = nrow(X), ncol = ncol(PC), 0) # a matrix full of zeros
  X1 <- rbind(PC, X0)
  X2 <- rbind(X0, PC)

  LP <- list(main = list(X = XMain, model = "BRR"), int1 = list(X = X1, model = "BRR"), int2 = list(X = X2, model = "BRR"))
  fmGRM <- BGLR(y = y, ETA = LP, nIter = nIter, burnIn = burnIn, saveAt = "GRM_", groups = rep(1:2, each = nrow(X)))

  return(fmGRM)
}


# Data --------------------------------------------------------------------

data(wheat)
# wheat.X |> dim()
# M <- wheat.X
#
# N <- nrow(M)
# m <- ncol(M)
# p <- colMeans(M) / 2

WWG <- function(M, p) {
  w <- scale(x = M, center = T, scale = F)

  S <- ((M == 2) * 1) * -rep(2 * (1 - p)^2, each = N) + ((M == 1) * 1) * rep(2 * p * (1 - p), each = N) + ((M == 0) * 1) * (-rep(2 * p^2, each = N))

  WWl <- w %*% t(w)
  Ga <- WWl / (sum(diag(WWl)) / N) + diag(1e-6, nrow(WWl))

  SSl <- S %*% t(S)
  Gd <- SSl / (sum(diag(SSl)) / N)

  return(list(Ga = Ga, Gd = Gd))
}
#G <- WWG(M = M, p = p)

# Getting a small sample of the dataset
Y <- wheat.Y[1:10, 1:2]
M <- wheat.X[1:10, ]
N <- nrow(M)
p <- colMeans(M) / 2
G <- WWG(M = M, p = p)

dat <- Y |>
  as.data.frame() |>
  melt() |>
  rename(c(
    var2 = variable,
    y = value
  )) |>
  mutate(var1 = as.factor(rep(c(1:nrow(Y)), 2)))

# hyper parameters

mmu <- mean(dat$y)
smu <- sd(dat$y)
stheta <- 1
a <- 0.1
b <- 0.1
Q <- 1
mb1 <- rep(0, length(levels(dat$var1)))
sb1 <- G$Ga

nIter <- 100
burnIn <- 100
nThin <- 2

hParam <- list(
  mmu = mmu,
  smu = smu,
  sb1 = sb1,
  a = a,
  b = b,
  mb1 = mb1,
  Q = Q,
  stheta = stheta,
  nIter = nIter,
  burnIn = burnIn,
  nThin = nThin
)


# Fitting models ----------------------------------------------------------

# gmBammit

gmBammit_fit <- gmBammit(data = dat, G = G$Ga, hParam = hParam)
plot(gmBammit_fit)

caret::RMSE(dat$y, gmBammit_fit$BUGSoutput$median$mu)
qplot(dat$y, gmBammit_fit$BUGSoutput$median$mu) +
  geom_abline() + theme_bw() + labs(x = "Y", y = expression(hat(Y)))


# M1
m1_fit <- m1_jags(data = dat, G = G$Ga, hParam = hParam)
plot(m1_fit)

caret::RMSE(dat$y, m1_fit$BUGSoutput$median$mu)
qplot(dat$y, m1_fit$BUGSoutput$median$mu) +
  geom_abline() + theme_bw() + labs(x = "Y", y = expression(hat(Y)))


# M2
m2_fit <- m2(Y = wheat.Y, M = wheat.X, nIter = nIter, burnIn = burnIn)

caret::RMSE(m2_fit$y, m2_fit$yHat)
qplot(m2_fit$y, m2_fit$yHat) +
  geom_abline() + theme_bw() + labs(x = "Y", y = expression(hat(Y)))

# M3

m3_fit <- m3(Y = wheat.Y,
             M = wheat.X,
             G = BGData::getG(wheat.X, center = TRUE, scaleG = TRUE, scale = TRUE),
             nIter = nIter,
             burnIn = burnIn)

caret::RMSE(m3_fit$y, m3_fit$yHat)
qplot(m3_fit$y, m3_fit$yHat) +
  geom_abline() + theme_bw() + labs(x = "Y", y = expression(hat(Y)))


