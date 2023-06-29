library(BGLR)
library(reshape2)
library(dplyr)
library(R2jags)
library(ggplot2)

remove(list = c())

# getting the genomic relationship matrix (2 methods)

data(wheat)
wheat.X |> dim()
M <- wheat.X

# Method 1: using the package AGHmatrix
#G <- AGHmatrix::Gmatrix(SNPmatrix = M, method = 'VanRaden')

# Second method
N <- nrow(M)
m <- ncol(M)
p <- colMeans(M)/2

WWG <- function(M, p){
  w <- scale(x = M, center = T, scale = F)

  S <- ((M==2)*1) * - rep(2*(1-p)^2, each=N) + ((M==1)*1) * rep(2*p*(1-p), each=N) + ((M==0)*1) * (-rep(2*p^2, each=N))

  WWl <- w %*% t(w)
  Ga <- WWl/(sum(diag(WWl))/N) + diag(1e-6, nrow(WWl))

  SSl <- S %*% t(S)
  Gd <- SSl/(sum(diag(SSl))/N)

  return(list(Ga=Ga,Gd=Gd))
}

G2 <- WWG(M = M, p = p)



# Normalize the genotypic data
# genotype_norm <- scale(M, center = TRUE, scale = TRUE)
# # Compute the VanRaden genomic relationship matrix
# G3 <- tcrossprod(genotype_norm) / ncol(genotype_norm)


# Taking G as the Cov matrix
# bammit with just g + e + int


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

      for(i in 1:B1) {
          b1[i] ~ dnorm(0, sb1[i,i]) # Prior on effect 1
      }

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


# Test with a small data set first:
Y <- wheat.Y[1:20, ]
M <- wheat.X[1:20, ]
G <- AGHmatrix::Gmatrix(SNPmatrix = M, method = 'VanRaden')

dat <- Y |>
  as.data.frame() |>
  melt() |>
  rename(c(var2 = variable,
           y = value)) |>
  mutate(var1 = as.factor(rep(c(1:nrow(Y)), 4)))

N <- nrow(dat)
Y <- dat$y
B1 <- length(levels(dat$var1))
B2 <- length(levels(dat$var2))
var1 <- dat$var1
var2 <- dat$var2
mmu <- mean(dat$y)
smu <- sd(dat$y)
stheta <- 1
a <- 0.1
b <- 0.1
Q <- 1
sb1 <- solve(G2$Ga)


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
  b = b
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
  n.thin = 2,
  n.burnin = 2000,
  n.iter = 2000 * 2 + 2000
)


# see convergence
plot(modelRun)

caret::RMSE(Y,modelRun$BUGSoutput$median$mu)
qplot(Y,modelRun$BUGSoutput$median$mu) +
  geom_abline() + theme_bw() + labs(x = 'Y', y = expression(hat(Y)))


# -------------------------------------------------------------------------

# considering fill cov matrix

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

Y <- wheat.Y[1:100, ]
M <- wheat.X[1:100, ]
#G <- AGHmatrix::Gmatrix(SNPmatrix = M, method = 'VanRaden')

dat <- Y |>
  as.data.frame() |>
  melt() |>
  rename(c(var2 = variable,
           y = value)) |>
  mutate(var1 = as.factor(rep(c(1:nrow(Y)), 4)))

N <- nrow(dat)
Y <- dat$y
B1 <- length(levels(dat$var1))
B2 <- length(levels(dat$var2))
var1 <- dat$var1
var2 <- dat$var2
mmu <- mean(dat$y)
smu <- sd(dat$y)
stheta <- 1
a <- 0.1
b <- 0.1
Q <- 1
mb1 <- rep(0, B1)
#sb1 <- solve(G2$Ga)
sb1 <- G2$Ga


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
  n.thin = 2,
  n.burnin = 2000,
  n.iter = 2000 * 2 + 2000
)

plot(modelRun)
caret::RMSE(Y,modelRun$BUGSoutput$median$mu)

qplot(Y,modelRun$BUGSoutput$median$mu) +
  geom_abline() + theme_bw() + labs(x = 'Y', y = expression(hat(Y)))

# test data


Y_test <- wheat.Y[101:nrow(wheat.Y), ]

dat_test <- Y_test |>
  as.data.frame() |>
  melt() |>
  rename(c(var2 = variable,
           y = value)) |>
  mutate(var1 = as.factor(rep(c(1:nrow(Y_test)), 4)))


# Prediction function


