sim.F.eigenvalues = function(a1, a2, n = 400, c, y, n_sims = 1) {
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
    Sigma1 = diag(x = 1, p) + diag(x = c(a1, a2, rep(0, p - 2)), p)
    Sigma2 = diag(p) 
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

c = 3/4
y = 1/4
h = sqrt(c + y - c*y)
a = (1 - h)^2 / (1 - y)^2
b = (1 + h)^2 / (1 - y)^2
gamma = 1 / (1 - y)

crit_left = gamma * (1 - h)
crit_right = gamma* (1 + h)
crit_left
crit_right

a1 = crit_left - 1
a2 = crit_right + 1

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

phi(a1, gamma, c) < a
phi(a2, gamma, c) > b

eval.list = sim.F.eigenvalues(a1, a2, c = c, y = y, n_sims = 1000)
evalues = unlist(eval.list[1])
eval_max = unlist(eval.list[2])
eval_min = unlist(eval.list[3])
hh = hist(evalues, breaks=100, xlim=c(-1,1.2*max(evalues)), freq=FALSE, col=8, main='')
x = seq(hh$breaks[1], hh$breaks[length(hh$breaks)], length.out=200)
lines(x, dfisher(x, c, y), type='l', lwd=2, col="blue", ylab="")
abline(v = phi(a1, gamma, c), lty=2, lwd=2, col='red')
abline(v = phi(a2, gamma, c), lty=2, lwd=2, col='red')

abline(v = a, lty=2, lwd=2, col='orange')
abline(v = b, lty=2, lwd=2, col='orange')
abline(v = mean(eval_max), lty=2, lwd=2, col='darkgreen')
abline(v = mean(eval_min), lty=2, lwd=2, col='darkgreen')

hist(eval_max, breaks=100, xlim=c(0,1.2*max(eval_max)), freq=FALSE, col=8, main='')
hist(eval_min, breaks=100, xlim=c(0,1.2*max(eval_min)), freq=FALSE, col=8, main='')

eval.list2 = sim.F.eigenvalues(0, 0, n_sims = 1000)
evalues2 = unlist(eval.list2[1])
eval_max2 = unlist(eval.list2[2])
eval_min2 = unlist(eval.list2[3])
hist(evalues2, breaks=100, xlim=c(-1,1.2*max(evalues)), freq=FALSE, col=8, main='')
