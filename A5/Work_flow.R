# 6017 ASSIGNMENT 5
# DUE DATE: 3 NOVEMBER 09:00

###### Packages Loading ######
library(mvnfast)
library(RMTstat)
###### Packages Loading ######


# Question 1
## (a)
c = 3/4
y = 1/4
h = sqrt(c + y - c*y)
a = (1 - h)^2 / (1 - y)^2
b = (1 + h)^2 / (1 - y)^2

## (b)
# c = 3/4
# y = 1/4
# h = sqrt(c + y - c*y)
gamma = 1 / (1 - y)
phi <- function(x ,gamma, c){
  if (x == gamma) {
    print("The x is equal to gamma so the value is not defined")
    return(NA)  # Return NA if x is equal to gamma
  }
  numerator = gamma * x * (x - 1 + c)
  denominator = x - gamma
  result = numerator / denominator
  return(result)
}

crit_left = gamma * (1 - h)
crit_right = gamma* (1 + h)

delta_1 = -1/10
delta_2 = 1/10

a_1 = 1 / (1 + delta_1)
a_2 = 1 / (1 + delta_2)

1 < a_1 && a < crit_right
crit_left < a_1 && a < 1

# Thus, this 2 values will converge to the a and b in the limiting.

## (c)
n = 400
# y_n = 1/4
# c_n = 3/4
sim.F.max.eigenvalues = function(delta_1, delta_2, n = 400, c, y, n_sims = 1) {
# 
#   Return the maximum eigenvalues of Fishier random matrices.
#   Assume the covariance structure: Sigma1 = diag(p)
#   Sigma2 = diag(p) + diag(x = c(delta_1, delta_2, rep(0, p - 2)), p)
# 
  p = n * y
  m = p / c
  evalues = c()
  for (sim in 1:n_sims) {
    Sigma1 = diag(p) 
    Sigma2 = diag(p) + diag(x = c(delta_1, delta_2, rep(0, p - 2)), p)
    mu = rep(0, p)
    X = rmvn(m+1, mu, Sigma1)
    Z = rmvn(n+1, mu, Sigma2)
    S1 = cov(X)
    S2 = cov(Z)
    FF = S1 %*% solve(S2)
    evalues = c(evalues, max(eigen(FF, only.values = TRUE)$values)) # only get the maximum eigenvalues
    }
  evalues 
}
p = n * y
m = p / c
s_p = (1/m * (sqrt(m) + sqrt(p)) * (1/sqrt(m) + 1/sqrt(p)))^{1/3}
lambdas = sim.F.max.eigenvalues(delta_1 = -1/10, delta_2 = 1/10, c = c, y = y, n_sims = 1000)

F1 = (lambdas - b) / s_p
range(F1)
mean(F1)
sd(F1)
hist(F1, breaks=80, xlim=c(-4,4)) -> h
plot(h$mids, h$density, type="s", xlab="", ylab="", ylim=c(0,1))
curve(dtw, -4, 4, lwd=2, lty=3, col="darkmagenta", add=TRUE)
title(main="Distribution of largest eigenvalue vs. TW density")


## (d)
# c = 3/4
# y = 1/4
h = sqrt(c + y - c*y)
# gamma = 1 / (1 - y) 
gamma > 1
# crit_left = gamma * (1 - h)
# crit_right = gamma* (1 + h)


kapa = 1 / (gamma * (1 + h)) - 1
ell = 1 / (gamma * (1 - h)) - 1


## (e)

delta1 = ell - 1/10
delta2 = kapa + 1/10

a_1 = 1 / (delta1 + 1)
a_2 = 1 / (delta2 + 1)

lambda_1 <- phi(x = a_1, gamma, c)
lambda_2 <- phi(x = a_2, gamma, c)
a
b
lambda_1
lambda_2
lambda_1 < a
lambda_2 > b


n = 400
# y_n = 1/4
# c_n = 3/4

sim.F.eigenvalues = function(delta_1, delta_2, n = 400, c = c_n, y = y_n, n_sims = 1) {
  # 
  #   Return the eigenvalues, maximum and minimum eigenvalues of Fishier random matrices.
  # 
  p = n * y
  m = p / c
  evalues = c()
  eval_max = c()
  eval_min = c()
  for (sim in 1:n_sims) {
    # Sigma1 = diag(p)
    # Sigma2 = diag(p) + diag(x = c(delta_1, delta_2, rep(0, p - 2)), p)
    Sigma1 = diag(x = 1, p)
    Sigma2 = diag(p) + diag(x = c(delta_1, delta_2, rep(0, p - 2)), p)
    # Sigma2 = diag(x = c(delta_1+1, delta_2+1, rep(1, p - 2)), p)
    mu = rep(0, p)
    X = rmvn(m+1, mu, Sigma1)
    Z = rmvn(n+1, mu, Sigma2)
    S1 = cov(X)
    S2 = cov(Z)
    FF = S1 %*% solve(S2)
    eval = eigen(FF, only.values = TRUE)$values
    eval_max = c(eval_max, max(eval)) 
    eval_min = c(eval_min, min(eval))
    evalues = c(evalues, eval)
  }
  return(list(evalues, eval_max, eval_min))
}

dfisher = Vectorize(function(x, c, y) { 
  h = sqrt(c + y - c*y)
  a = (1 - h)^2 / (1 - y)^2
  b = (1 + h)^2 / (1 - y)^2
  ifelse(x <= a | x >= b, 0, suppressWarnings(sqrt((x - a) * (b - x))*(1-y)/(2 * pi * x *(c+y*x))))
},"x")

eval.list = sim.F.eigenvalues(delta1, delta2, n_sims = 1000)
evalues = unlist(eval.list[1])
eval_max = unlist(eval.list[2])
eval_min = unlist(eval.list[3])
hh = hist(evalues, breaks=100, xlim=c(-1,1.2*max(evalues)), freq=FALSE, col=8, main='')
x = seq(hh$breaks[1], hh$breaks[length(hh$breaks)], length.out=200)
lines(x, dfisher(x, c, y), type='l', lwd=2, col="blue", ylab="")
abline(v = lambda_1, lty=2, lwd=2, col='red')
abline(v = lambda_2, lty=2, lwd=2, col='red')

# abline(v = a, lty=2, lwd=2, col='orange')
# abline(v = b, lty=2, lwd=2, col='orange')
abline(v = mean(eval_max), lty=2, lwd=2, col='darkgreen')
abline(v = mean(eval_min), lty=2, lwd=2, col='darkgreen')

hist(eval_max, breaks=100, xlim=c(0,1.2*max(eval_max)), freq=FALSE, col=8, main='')
hist(eval_min, breaks=100, xlim=c(0,1.2*max(eval_min)), freq=FALSE, col=8, main='')

eval.list2 = sim.F.eigenvalues(0, 0, n_sims = 1000)
evalues2 = unlist(eval.list2[1])
eval_max2 = unlist(eval.list2[2])
eval_min2 = unlist(eval.list2[3])
hist(evalues2, breaks=100, xlim=c(-1,1.2*max(evalues)), freq=FALSE, col=8, main='')

## (f)


