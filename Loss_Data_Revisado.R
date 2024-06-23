library(actuar)
library(mixtools)
library(fBasics)

######################################################### Loading Data ############################################################

loss <- na.omit(loss)
loss <- loss[which(loss > 0)]

################################################### Descriptive Statistics ########################################################

## Summary Statistics
 
options(scipen=999)
basicStats(loss)
options(scipen=0)

## Limits to be applied in the graphs

max_lim <- 20
min_lim <- 0

N <- length(loss)

## Histogram

hist(loss[which(loss < max_lim)],plot=TRUE,prob=TRUE, breaks = 250, main = "", xlab = "Loss Severity", ylab = "Frequency")

## Exponential

func_ML_exp <- function(coef){
  
  lambda <- coef[1]
  
  ln_L = -(N*log(lambda) - lambda*sum(loss))
  
  return(ln_L)
  
}

lambda.c <- 0.5
exp <- optim(lambda.c, func_ML_exp,control = list(reltol = 1e-15))
hess <- numDeriv:::hessian(func_ML_exp, exp$par)
ep_exp <- sqrt(diag(solve(hess)))

# Goodness-of-Fit Measures

exp_loglik <- -as.numeric(exp$value)
exp_aic <- 2*(length(exp$par) - exp_loglik)
exp_bic <- length(exp$par) * log(N) - 2 * exp_loglik

# Histogram

dexp <- function(x,lambda){
  
  dexp <- lambda*exp(-(lambda*x))
  
}

hist(loss[which(loss < max_lim)],plot=TRUE,prob=TRUE, breaks = 250, main = "", xlab = "Loss Severity", ylab = "Frequency")
curve(dexp(x,exp$par),add=TRUE,lwd=2,lty=4)
legend("topright","Exponential",lty=4,lwd=2)

# Distribution Function

pexp <- function(x,lambda){
  
 pexp <- 1 - exp(-(lambda*x)) 
  
}

plot(ecdf(loss[which(loss < max_lim | loss > min_lim)]),do.points=FALSE,verticals=TRUE,xlim=c(min_lim,max_lim),main="Cumulative Distribution Function - Loss")
curve(pexp(x,exp$par),add=TRUE,lwd=2,lty=4)
legend("bottomright","Exponential",lty=2,lwd=2)

exp_ks <- as.numeric(ks.test(loss,"pexp",exp$par)$statistic)

# Quantile Function

qexp <- function(x,lambda){
  
  qexp <- -(log(1-x)/lambda)
  
}

plot(qexp((1:N-0.5)/N, exp$par), sort(loss), main="", xlab="Exponential Quantile", ylab= "Sample Quantile")
abline(0,1)

## Boxplot

boxplot(loss,main="",ylab="Loss Severity")

## Empirical Quantiles

quantile(loss,probs = c(0.9, 0.95, 0.975, 0.99, 0.995, 0.999, 0.9995, 0.9999),type=8)

############################################################## Fitting #############################################################

### Parametric Models

# Log-Normal

# Likelihood Function

func_ML_lnorm <- function(coef){
  
  mu <- coef[1]
  sigma <- coef[2]
  
  ln_L = -(-(N/2)*log(2*pi) - (N/2)*log(sigma^2) - (1/(2*sigma^2))* sum((log(loss)-mu)^2) - sum(log(loss)))
  
  return(ln_L)
  
}

mu.c <- 1
sigma.c <- 1
lnorm <- optim(c(mu.c,sigma.c),func_ML_lnorm,control = list(reltol = 1e-15))
hess <- numDeriv:::hessian(func_ML_lnorm, lnorm$par)
ep_lnorm <- sqrt(diag(solve(hess)))

# Goodness-of-Fit Measures

lnorm_loglik <- -lnorm$value
lnorm_aic <- 2*(length(lnorm$par) - lnorm_loglik)
lnorm_bic <- length(lnorm$par) * log(N) - 2 * lnorm_loglik

# Histogram

dlnorm <- function(x,mu,sigma){
  
  dlnorm <- (1/(x*sigma*sqrt(2*pi)))*exp(-((log(x) - mu)^2)/(2*(sigma^2)))
  
  return(dlnorm)
  
}

hist(loss[which(loss < max_lim)],plot=TRUE,prob=TRUE,breaks =150, main = "Histogram of Losses", xlab = "Severity", ylab = "Frequency")
curve(dlnorm(x,lnorm$par[1],lnorm$par[2]),add=TRUE,lwd=2,lty=4)
legend("topright","Log-Normal",lty=2,lwd=2)

# Distribution Function

plnorm <- function(x,mu,sigma){
  
  plnorm <- pnorm((log(x) - mu)/sigma)
  
  return(plnorm)
  
}

plot(ecdf(loss[which(loss < max_lim | loss > min_lim)]),do.points=FALSE,verticals=TRUE,xlim=c(min_lim,max_lim),main="Cumulative Distribution Function - Loss")
curve(plnorm(x,lnorm$par[1],lnorm$par[2]),add=TRUE,lwd=2,lty=4)
legend("bottomright","Log-Normal",lty=2,lwd=2)

lnorm_ks <- as.numeric(ks.test(loss,"plnorm",lnorm$par[1],lnorm$par[2])$statistic)

# Quantile Function

qlnorm <- function(x,mu,sigma){
  
  qlnorm <- exp(mu + sigma*qnorm(x))
  
  return(qlnorm)
  
}

qlnorm(c(0.9,0.95,0.975,0.99,0.995,0.999,0.9995,0.9999), lnorm$par[1],lnorm$par[2])
plot(qlnorm((1:N-0.5)/N, lnorm$par[1],lnorm$par[2]), sort(loss), main = "Lognormal", xlab="Theoretical Quantiles", ylab="Sample Quantiles", ylim = c(0,max(loss)))
abline(0,1)

# Gamma

# Likelihood Function

func_ML_gamma <- function(coef){
  
  alpha <- coef[1]
  beta <- coef[2]
  
  ln_L = -((alpha-1)*sum(log(loss)) -(1/beta)*sum(loss) - N*(alpha*log(beta) + log(gamma(alpha))))
  
  return(ln_L)
  
}

alpha.c <- 0.5
beta.c <- 0.5
gamma <- optim(c(alpha.c,beta.c),func_ML_gamma,control = list(reltol = 1e-15))
hess <- numDeriv:::hessian(func_ML_gamma, gamma$par)
ep_gamma <- sqrt(diag(solve(hess)))

# Goodness-of-Fit Measures

gamma_loglik <- -gamma$value
gamma_aic <- 2*(length(gamma$par) - gamma_loglik)
gamma_bic <- length(gamma$par) * log(N) - 2 * gamma_loglik

# Histogram

dgamma <- function(x,alpha,beta){

  dgamma <- (1/((beta^alpha)*gamma(alpha)))*(x^(alpha-1))*exp(-(x/beta))  
  
}

hist(loss[which(loss < max_lim)],plot=TRUE,prob=TRUE,breaks = 150, main = "Histogram of Losses", xlab = "Severity", ylab = "Frequency")
curve(dgamma(x,gamma$par[1],gamma$par[2]),add=TRUE,lwd=2,lty=4)
legend("topright","Gamma",lty=2,lwd=2)

# Distribution Function

plot(ecdf(loss[which(loss < max_lim | loss > min_lim)]),do.points=FALSE,verticals=TRUE,xlim=c(min_lim,max_lim),main="Cumulative Distribution Function - Loss")
curve(pgamma(x,gamma$par[1],1/gamma$par[2]),add=T,lwd=2,lty=2)
legend("bottomright","Gamma",lty=2,lwd=2)

gamma_ks <- as.numeric(ks.test(loss,"pgamma",gamma$par[1],1/gamma$par[2])$statistic)

# Quantile Function

qgamma(c(0.90,0.95,0.975,0.99,0.995,0.999,0.9995,0.9999), gamma$par[1], 1/gamma$par[2])
plot(qgamma((1:N-0.5)/N, gamma$par[1], 1/gamma$par[2]), sort(loss), main = "Gamma" ,xlab = "Theoretical Quantiles", ylab="Sample Quantiles", ylim = c(0,max(loss)))
abline(0,1)

# Weibull

func_ML_weibull <- function(coef){
  
  alpha <- coef[1]
  beta <- coef[2]
  
  ln_L = -((N*(log(alpha)-alpha*log(beta))) + ((alpha - 1)*sum(log(loss))) - sum((loss/beta)^alpha))
  
  return(ln_L)
  
}

alpha.c <- 1
beta.c <- 3
weibull <- optim(c(alpha.c,beta.c),func_ML_weibull,control = list(reltol = 1e-35))
hess <- numDeriv:::hessian(func_ML_weibull, weibull$par)
ep_weibull <- sqrt(diag(solve(hess)))

# Goodness-of-Fit Measures

weibull_loglik <- -weibull$value
weibull_aic <- 2*(length(weibull$par) - weibull_loglik)
weibull_bic <- length(weibull$par) * log(N) - 2 * weibull_loglik

# Histogram

dweibull <- function(x,alpha,beta){
  
  dweibull <- (alpha/beta^alpha)*(x^(alpha-1))*exp(-(x/beta)^alpha)
  
  
}

hist(loss[which(loss < max_lim)],plot=TRUE,prob=TRUE,breaks = 150, main = "Histogram of Losses", xlab = "Severity", ylab = "Frequency")
curve(dweibull(x,weibull$par[1],weibull$par[2]),add=TRUE,lwd=2,lty=4)
legend("topright","Weibull",lty=2,lwd=2)

# Distribution Function

pweibull <- function(x,alpha,beta){
  
  pweibull <- 1 - exp(-(x/beta)^alpha)
  
}

plot(ecdf(loss[which(loss < max_lim | loss > min_lim)]),do.points=FALSE,verticals=TRUE,xlim=c(min_lim,max_lim),main="Cumulative Distribution Function - Loss")
curve(pweibull(x,weibull$par[1],weibull$par[2]),add=TRUE,lwd=2,lty=4)
legend("bottomright","Weibull",lty=2,lwd=2)

weibull_ks <- as.numeric(ks.test(loss, "pweibull", weibull$par[1], weibull$par[2])$statistic)

# Quantile Function

qweibull <- function(x,alpha,beta){
  
  qweibull <- beta*((-log(1-x))^(1/alpha))
  
  return(qweibull)
  
}

qweibull(c(0.9,0.95,0.975,0.99,0.995,0.999,0.9995,0.9999), weibull$par[1], weibull$par[2])
plot(qweibull((1:N-0.5)/N, weibull$par[1], weibull$par[2]), sort(loss), main = "Weibull" ,xlab = "Theoretical Quantiles", ylab="Sample Quantiles", ylim = c(0,max(loss)))
abline(0,1)

# Pareto

# Likelihood Function

func_ML_pareto <- function(coef){
  
  a <- coef[1]
  c <- min(loss)
  
  ln_L = -((N * log(a)) + (N * a * log(c)) - (a + 1) * sum(log(loss)))
    
  return(ln_L)
  
}

a.c <- 0.1

pareto <- nlm(func_ML_pareto, a.c, steptol = 1e-10, gradtol = 1e-30, iterlim = 100000, ndigit = 12, fscale = 0, print.level = 1)
pareto <- optim(a.c,func_ML_pareto)
a.pareto <- pareto$par
c.pareto <- min(loss)
hess <- numDeriv:::hessian(func_ML_pareto, pareto$par)
ep_pareto <- sqrt(diag(solve(hess)))

# Goodness-of-Fit Measures

pareto_loglik <- -pareto$value
pareto_aic <- 2*((length(pareto$par) + 1) - pareto_loglik)
pareto_bic <- (length(pareto$par) + 1) * log(length(loss)) - 2 * pareto_loglik

# Histogram

dpareto <- function(x,a,c){
  
  dpareto <- (a*c^a)/(x^(a+1))
  
}

hist(loss[which(loss < max_lim)],plot=TRUE,prob=TRUE,breaks = 150, main = "Histogram of Losses", xlab = "Severity", ylab = "Frequency")
curve(dpareto(x,a.pareto,c.pareto),xlim=c(c.pareto, max_lim), add=TRUE,lwd=2,lty=4)
legend("topright","Pareto",lty=2,lwd=2)
abline(v=c.pareto,lty=2)

# Distribution Function

ppareto <- function(x,a,c){
  
  ppareto <- 1 - ((c/x)^a)
  
}

plot(ecdf(loss[which(loss < max_lim | loss > min_lim)]),do.points=FALSE,verticals=TRUE,xlim=c(min_lim,max_lim),main="Cumulative Distribution Function - Loss")
curve(ppareto(x,a.pareto,c.pareto),xlim=c(c.pareto, max_lim),add=TRUE,lwd=2,lty=4)
abline(v=c.pareto,lty=2)
legend("bottomright","Pareto",lty=2,lwd=2)

pareto_ks <- as.numeric(ks.test(loss, "ppareto", a.pareto, c.pareto)$statistic)

# Quantile Function

qpareto <- function(x,a,c){
  
  qpareto <- (c/((1-x)^(1/a)))
  
  return(qpareto)
  
}

qpareto(c(0.9,0.95,0.975,0.99,0.995,0.999,0.9995,0.9999),a.pareto,c.pareto)
plot(qpareto((1:N-0.5)/N, a.pareto, c.pareto), sort(loss), main = "Pareto" ,xlab = "Theoretical Quantiles", ylab="Sample Quantiles")
abline(0,1)

# Inverse Gaussian

func_ML_invnorm <- function(coef){
  
  mu <- coef[1]
  lambda <- coef[2]
  
  ln_L = -((N/2)*log(lambda/(2*pi)) - (3/2)*sum(log(loss)) - lambda*sum(((loss-mu)^2)/(2*(mu^2)*loss)))
  
  return(ln_L)
  
}

invnorm <- optim(c(mean(loss),1/((sum(1/loss)/N) - (1/mean(loss)))),func_ML_invnorm)
hess <- numDeriv:::hessian(func_ML_invnorm,c(invnorm$par[1],invnorm$par[2]))
ep_invnorm <- sqrt(diag(solve(hess)))

# Goodness-of-Fit Measures

invnorm_loglik <- -invnorm$value
invnorm_aic <- 2*(length(invnorm$par) - invnorm_loglik)
invnorm_bic <- length(invnorm$par) * log(N) - 2 * invnorm_loglik

# Histogram

dinvgauss <- function(x,mu,lambda){
  
  dinvgauss <- sqrt(lambda/(2*pi*x^3))*exp((-lambda*(x-mu)^2)/(2*mu^2*x))
  
}

hist(loss[which(loss < max_lim)],plot=TRUE,prob=TRUE,breaks = 150, main = "Histogram of Losses", xlab = "Severity", ylab = "Frequency")
curve(dinvgauss(x,invnorm$par[1],invnorm$par[2]),add=TRUE,lwd=2,lty=4)
legend("topright","Inverse Gaussian",lty=2,lwd=2)

# Distribution Function

pinvgauss <- function(x,mu,lambda){
  
  pinvgauss <- pnorm(sqrt(lambda/x)*((x/mu) - 1)) + exp((2*lambda)/mu)*pnorm(-sqrt(lambda/x)*((x/mu) + 1))
  
}

plot(ecdf(loss[which(loss < max_lim | loss > min_lim)]),do.points=FALSE,verticals=TRUE,xlim=c(min_lim,max_lim),main="Cumulative Distribution Function - Loss")
curve(pinvgauss(x,invnorm$par[1],invnorm$par[2]),add=TRUE,lwd=2,lty=4)
legend("bottomright","Inverse Gaussian",lty=2,lwd=2)

invnorm_ks <- as.numeric(ks.test(loss,"pinvgauss",invnorm$par[1],invnorm$par[2])$statistic)

# Quantile Function

actuar:::qinvgauss(c(0.9,0.95,0.975,0.99,0.995,0.999,0.9995,0.9999),invnorm$par[1],invnorm$par[2])
plot(actuar:::qinvgauss((1:N-0.5)/N, invnorm$par[1], invnorm$par[2]), sort(loss), main = "Inverse Gaussian" ,xlab = "Theoretical Quantiles", ylab="Sample Quantiles", ylim = c(0,max(loss)))
abline(0,1)

# Inverse Gamma

func_ML_invgamma <- function(coef){
  
  alpha <- coef[1]
  beta <- coef[2]
  
  ln_L = -(-N*(alpha+1)*mean(log(loss)) - N*log(gamma(alpha)) + N*alpha*log(beta) - sum(beta/loss))
  
  return(ln_L)
  
}

invgamma <- optim(c(1,5),func_ML_invgamma)
hess <- numDeriv:::hessian(func_ML_invgamma,c(invgamma$par[1],invgamma$par[2]))
ep_invgamma <- sqrt(diag(solve(hess)))

# Goodness-of-Fit Measures

invgamma_loglik <- -invgamma$value
invgamma_aic <- 2*(length(invgamma$par) - invgamma_loglik)
invgamma_bic <- (length(invgamma$par) * log(length(loss))) - 2 * invgamma_loglik

# Histogram

dinvgamma <- function(x,alpha,beta){
  
  dinvgamma <- ((beta^alpha)/gamma(alpha))*(x^(-alpha - 1))*exp(-(beta/x))
  
}

hist(loss[which(loss < max_lim)],plot=TRUE,prob=TRUE,breaks = 150, main = "Histogram of Losses", xlab = "Severity", ylab = "Frequency")
curve(dinvgamma(x,invgamma$par[1], invgamma$par[2]),add=TRUE,lwd=2,lty=4)
legend("topright","Inverse Gamma",lty=2,lwd=2)

# Distribution Function

pinvgamma <- function(x,alpha,beta){
  
  incgamma <- function(u){
    
    incgamma <- (u^(alpha-1)) * exp(-u)
    
    return(incgamma)
    
  }
  
  pinvgamma <- tryCatch(
    {
      integrate(incgamma,as.numeric(beta)/x,Inf)$value/gamma(as.numeric(alpha)) 
    },
    error=function(e) {
      0
    })
  
  return(pinvgamma)
  
}

pinvgamma <- Vectorize(pinvgamma)

plot(ecdf(loss[which(loss < max_lim | loss > min_lim)]),do.points=FALSE,verticals=TRUE,xlim=c(min_lim,max_lim),main="Cumulative Distribution Function - Loss")
curve(pinvgamma(x,invgamma$par[1], invgamma$par[2]),add=TRUE,lwd=2,lty=4)
legend("bottomright","Inverse Gamma",lty=2,lwd=2)

invgamma_ks <- as.numeric(ks.test(loss, "pinvgamma", invgamma$par[1],invgamma$par[2])$statistic)

# Quantile Function

actuar:::qinvgamma(c(0.9,0.95,0.975,0.99,0.995,0.999,0.9995,0.9999),invgamma$par[1], 1/invgamma$par[2])
plot(actuar:::qinvgamma((1:N-0.5)/N, invgamma$par[1], 1/invgamma$par[2]), sort(loss), main = "Inverse Gamma" , xlab = "Theoretical Quantiles", ylab="Sample Quantiles", ylim = c(0,max(loss)))
abline(0,1)

# Burr Type XII

func_ML_burr <- function(coef){
  
  lambda <- coef[1]
  k <- coef[2]
  alpha <- coef[3]
  
  ln_L = -(N*log(k) + N*log(lambda) + (lambda-1)*sum(log(loss)) - N*lambda*log(alpha) - (1+k)*sum(log(1+(loss/alpha)^lambda)))
  
  return(ln_L)
  
}

burr <- optim(c(1,1,20),func_ML_burr,control = list(reltol = 1e-25))
hess <- numDeriv::hessian(func_ML_burr,c(burr$par[1],burr$par[2],burr$par[3]))
ep_burr <- sqrt(diag(solve(hess)))

# Goodness-of-Fit Measures

burr_loglik <- -burr$value
burr_aic <- 2*(length(burr$par) - burr_loglik)
burr_bic <- (length(burr$par) * log(N)) - 2 * burr_loglik

# Histogram

dburr <- function(x,lambda,k,alpha){
  
  dburr <- (((k*lambda)/alpha)*((x/alpha)^(lambda-1)))/((1+(x/alpha)^lambda)^(k+1))
  
}

hist(loss[which(loss < max_lim)],plot=TRUE,prob=TRUE,breaks = 150, main = "Histogram of Losses", xlab = "Severity", ylab = "Frequency")
curve(dburr(x,burr$par[1], burr$par[2],burr$par[3]),add=TRUE,lwd=2,lty=4)
legend("topright","Burr Type XII",lty=2,lwd=2)

# Distribution Function

pburr <- function(x,lambda,k,alpha){
  
  pburr <- 1 - (1/((1+(x/alpha)^lambda)^k))
  
}

plot(ecdf(loss[which(loss < max_lim | loss > min_lim)]),do.points=FALSE,verticals=TRUE,xlim=c(min_lim,max_lim),main="Cumulative Distribution Function - Loss")
curve(pburr(x,burr$par[1], burr$par[2], burr$par[3]),add=TRUE,lwd=2,lty=4)
legend("bottomright","Burr Type XII",lty=2,lwd=2)

burr_ks <- as.numeric(ks.test(loss, "pburr", burr$par[1], burr$par[2], burr$par[3])$statistic)

# Quantile Function

qburr <- function(x,lambda,k,alpha){
  
  qburr <- alpha*(((1-x)^(-1/k) - 1)^(1/lambda))
  
  return(qburr)
  
}

qburr(c(0.9,0.95,0.975,0.99,0.995,0.999,0.9995,0.9999),burr$par[1], burr$par[2], burr$par[3])
plot(qburr((1:N-0.5)/N, burr$par[1], burr$par[2], burr$par[3]), sort(loss), main = "Burr Type XII" , xlab = "Theoretical Quantiles", ylab="Sample Quantiles", ylim = c(0,max(loss)))
abline(0,1)

# Generalized Lambda

glambda <- fBasics:::gldFit(loss, lambda1 = 1, lambda2 = -0.4, lambda3 = -1, lambda4 = -5, method = "mle", doplot = FALSE)

func_ML_gld <- function(coef){
  
  lambda1 <- coef[1]
  lambda2 <- coef[2]
  lambda3 <- coef[3]
  lambda4 <- coef[4]
  
  ln_L <- -(sum(log(fBasics:::dgld(loss,lambda1,lambda2,lambda3,lambda4))))
  
}

# Goodness-of-Fit Measures

glambda_loglik <- glambda@fit$minimum
glambda_aic <- 2*(length(glambda@fit$estimate) - glambda@fit$minimum)
glambda_bic <- length(glambda@fit$estimate) * log(length(loss)) - 2 * glambda@fit$minimum
glambda_ks <- as.numeric(ks.test(loss,"pgld",lambda1 = glambda@fit$estimate[1], lambda2 = glambda@fit$estimate[2], lambda3 = glambda@fit$estimate[3], lambda4 = glambda@fit$estimate[4])$statistic)

# Histogram

hist(loss[which(loss < max_lim)],plot=TRUE,prob=TRUE,breaks = 150, main = "Histogram of Losses", xlab = "Severity", ylab = "Frequency")
curve(dgld(x,lambda1 = glambda@fit$estimate[1], lambda2 = glambda@fit$estimate[2], lambda3 = glambda@fit$estimate[3], lambda4 = glambda@fit$estimate[4]),add=TRUE,lwd=2,lty=4)
legend("topright","Generalized Lambda",lty=2,lwd=2)

# Distribtuion Function

plot(ecdf(loss[which(loss < max_lim | loss > min_lim)]),do.points=FALSE,verticals=TRUE,xlim=c(min_lim,max_lim),main="Cumulative Distribution Function - Loss")
curve(fBasics:::pgld(x,lambda1 = glambda@fit$estimate[1], lambda2 = glambda@fit$estimate[2], lambda3 = glambda@fit$estimate[3], lambda4 = glambda@fit$estimate[4]),add=TRUE,lwd=2,lty=4)
legend("bottomright","Generalized Lambda",lty=2,lwd=2)

# Quantile Function

qgld <- function(x,lambda1,lambda2,lambda3,lambda4){
  
  qgld <- lambda1 + (x^lambda3 - (1-x)^lambda4)/lambda2
  
  return(qgld)
  
}

qgld(c(0.9,0.95,0.975,0.99,0.995,0.999,0.9995,0.9999), lambda1 = glambda@fit$estimate[1], lambda2 = glambda@fit$estimate[2], lambda3 = glambda@fit$estimate[3], lambda4 = glambda@fit$estimate[4])
plot(qgld((1:N-0.5)/N, lambda1 = glambda@fit$estimate[1], lambda2 = glambda@fit$estimate[2], lambda3 = glambda@fit$estimate[3], lambda4 = glambda@fit$estimate[4]), sort(loss), main = "Generalized Lambda" ,xlab = "Theoretical Quantiles", ylab="Sample Quantiles", ylim = c(0,max(loss)))
abline(0,1)

# Mixture Gamma

mixgamma <- mixtools:::gammamixEM(loss,k=6,maxit=100000)

# Goodness-of-Fit Measures

mixgamma_loglik <- mixgamma$loglik
mixgamma_aic <- 2*((length(mixgamma$gamma.pars) + (length(mixgamma$lambda) - 1)) - mixgamma_loglik)
mixgamma_bic <- (length(mixgamma$gamma.pars) + (length(mixgamma$lambda) - 1)) * log(length(loss)) - 2 * mixgamma_loglik
mixgamma_ks <- as.numeric(ks.test(loss,"pmgamma",mgshape=c(mixgamma$gamma.pars[1,]),mgscale=c(mixgamma$gamma.pars[2,]),mgweight=c(mixgamma$lambda))$statistic)

# Histogram

hist(loss[which(loss < max_lim)],plot=TRUE,prob=TRUE, 250, main = "Histogram of Losses", xlab = "Severity", ylab = "Frequency")
curve(evmix:::dmgamma(x, mgshape = c(mixgamma$gamma.pars[1,]), mgscale = c(mixgamma$gamma.pars[2,]), mgweight = c(mixgamma$lambda)),add=TRUE,lwd=2,lty=4)
legend("topright","Mixture Gamma",lty=2,lwd=2)

# Distribution Function

plot(ecdf(loss[which(loss < max_lim | loss > min_lim)]),do.points=FALSE,verticals=TRUE,xlim=c(min_lim,max_lim),main="Cumulative Distribution Function - Loss")
curve(evmix:::pmgamma(x, mgshape = c(mixgamma$gamma.pars[1,]), mgscale = c(mixgamma$gamma.pars[2,]), mgweight = c(mixgamma$lambda)),add=TRUE,lwd=2,lty=4)
legend("bottomright","Mixture Gamma",lty=2,lwd=2)

# Quantile Function

evmix:::qmgamma(c(0.9,0.95,0.975,0.99,0.995,0.999,0.9995,0.9999),mgshape = c(mixgamma$gamma.pars[1,]), mgscale = c(mixgamma$gamma.pars[2,]), mgweight = c(mixgamma$lambda))
plot(evmix:::qmgamma((1:N-0.5)/N, mgshape = c(mixgamma$gamma.pars[1,]), mgscale = c(mixgamma$gamma.pars[2,]), mgweight = c(mixgamma$lambda)), sort(loss), main = "Mixture Gamma" ,xlab = "Theoretical Quantiles", ylab="Sample Quantiles", ylim = c(0,max(loss)))
abline(0,1)

### Semi-parametric Models

## Corrected Generalized Champernowne Distribuition (Buch-Kromann, 2006)

# Step 1

## Maximum-Likelihood

M <- median(loss)

func_ML_champer <- function(coef){
  
  a <- coef[1]
  c <- coef[2]
  
  ln_L = -(N*log(a) + N*log((M + c)^(a) - c^a) + (a - 1)*sum(log(loss + c)) - 2* sum(log((loss +c)^(a) + (M + c)^(a) - 2*c^a)))
  
  return(ln_L)
  
}

a.c <- 1
c.c <- 50

champernowne <- nlm(func_ML_champer,c(a.c,c.c), steptol = 1e-10, gradtol = 1e-10, iterlim = 100000, ndigit = 12, fscale = 0, print.level = 1)

a.hat <- champernowne$estimate[1]
c.hat <- champernowne$estimate[2]
M.hat <- median(loss)

# Distribution Function

T_hat <- function(x){
  
  T_hat = ((x + c.hat)^(a.hat) - c.hat^a.hat)/((x + c.hat)^(a.hat) + (M.hat + c.hat)^a.hat - 2*c.hat^a.hat)
  
  return(T_hat)
  
}

# Density Function

t_hat <- function(x){
  
  t_hat = (a.hat*(x + c.hat)^(a.hat - 1)*((M.hat + c.hat)^(a.hat) - c.hat^a.hat))/(((x + c.hat)^(a.hat) + (M.hat + c.hat)^a.hat - 2*c.hat^a.hat)^2)
  
  return(t_hat)
  
}

# Quantile Function

T_hat_inv <- function(x){
  
  T_hat_inv = (((x*((M.hat + c.hat)^a.hat) - (2*x - 1)*(c.hat^a.hat))/(1 - x))^(1/a.hat)) - c.hat
  
  return(T_hat_inv)
  
}

# Histogram

hist(loss[which(loss < max_lim)],plot=TRUE,prob=TRUE,breaks = 150, main = "", xlab = "Loss Severity", ylab = "Frequency")
curve(t_hat, add = TRUE,lwd=2,lty = 2,col=1)
legend("topright", c("Modified Champernowne"), lty=2,lwd=2,col=1,bty="n")

# Distribtution Function

plot(ecdf(loss[which(loss < max_lim | loss > min_lim)]),xlim=c(min_lim,max_lim),do.points=FALSE,verticals=TRUE)
curve(T_hat,add = TRUE,lty=2,lwd=2,col=1)
legend("bottomright", c("Modified Champernowne"), lty=2,lwd=2,col=1)

## Quantile-Mean Method

# Empirical Quantile

quantile(loss,0.95)

# Theoretical Quantile

T_hat_inv(0.95)

# Shape Parameter

a.hat <- 1.36

# Empirical Mean

mean(loss)

# Theoretical Mean

integrand <- function(x){
  
  integrand = x * t_hat(x)
  
  return(integrand)
  
}

E_hat <- integrate(integrand,lower=0,upper=Inf)$value
E_hat

# Location Parameter

c.hat <- 1

# Step 2

transf_loss <- T_hat(loss)

# Step 3

b_transf <- (((40*pi^0.5)/N)^(0.2))*sd(transf_loss)
x <- seq(0, 1, 1/length(transf_loss))
g.hat.corrected <- evmix:::dbckden(x,transf_loss,kernel = "epanechnikov",lambda=b_transf)
g.hat.corrected <- cbind(x,g.hat.corrected)

# Histogram Transformed Losses

hist(transf_loss, prob = TRUE, breaks = round(length(transf_loss)^0.5), main = "",xlab="Transformed Loss Severity",ylab="Frequency")
lines(g.hat.corrected, lwd = 1)
abline(h=1,lwd=1,lty=2)
legend("topright",c("Correction Estimator"), lty=1, lwd=1, bty="n")

# Step 4

grid <- 1/100
x_min <- 0
x_max <- max(loss)
x <- seq(x_min,x_max,grid)

sx <- sd(loss)
sy <- sd(transf_loss)

scal <- sx/sy

y <- T_hat(loss)
yscal <- y * scal

tgrid <- T_hat(x)
tgscal <- tgrid *scal

dtg <- t_hat(x)
dtgscal <- dtg * scal

# Bandwidth Selection

b <- (((40*pi^(1/2))/N)^(1/5)) * (sd(loss))

# Kernel Function

ker <- function(x){kader:::epanechnikov(x)}

# Transformed Kernel Density Function

f_hat <- c()

for (i in 1:length(x)){
  
  diff <- (tgscal[i] - yscal)/b
  
  f_hat[i] <- dtgscal[i] * (1/(N*b)) * sum(ker(diff))
  
}

# Transformed Kernel Distribution Function

F_hat <- c()

F_hat[1] <- f_hat[1]

for (i in 2:length(x)){
  
  F_hat[i] <- F_hat[i-1] + f_hat[i]
   
}

F_hat <- F_hat/sum(f_hat)
F_hat <- append(F_hat,0,after=0)

# Transformed Kernel Quantile Function

options(scipen = 999)
Q_hat <- cbind(x,F_hat[-length(x)])
options(scipen = 0)

# Histogram

hist(loss[which(loss < max_lim)],plot=TRUE,prob=TRUE,breaks = 150, main = "", xlab = "Loss Severity", ylab = "Frequency")
lines(f_hat~x, lwd=2, lty=1, col=1)
curve(t_hat, add = TRUE, lwd=2, lty=2, col=1)
legend("topright", c("Transformed Kernel Density","Modified Champernowne"), lwd = c(2,2), lty=c(1,2), col=c(1,1), bty = "n")

# Quantile Function

plot(ecdf(loss[which(loss < max_lim | loss > min_lim)]),xlim=c(min_lim,max_lim),do.points=FALSE,verticals=TRUE, main="")
lines(F_hat[-length(x)]~x,lwd=2,lty=1,col=1)
legend("bottomright", c("Transformed Kernel Density"), lty=1,col=1,lwd=2)

# Kolmogorov-Smirnov Goodness-of-Fit Test

F_emp <- function(sample,x){
  count = 0
  for(n in sample){
    if(n <= x){
      count = count + 1
    }
  }
  return(count/length(sample))
}

F_hat_emp <- c()
D.aux <- c()

for(i in 1:length(x)){
  
  F_hat_emp[i] <- F_emp(loss,x[i])
  D.aux[i] <- abs(F_hat_emp[i] - F_hat[-length(x)][i])
   
}

transfkernel_ks <- max(D.aux)

## Threshold Selection Method

# Hill Estimator

hill <- function(k){
  
  order_loss <- sort(loss, decreasing = T)
  
  H <- mean((log(order_loss[1:k]) - log(order_loss[k + 1])))
  
  return(H)
  
}

H <- sapply(1:(N-1),hill)

# Hill Plot

plot(1/H[15:length(H)], type = "l", xlim = c(1,length(H)), main="", xlab="Order Statistics", ylab="Hill Estimate")
axis(side = 3, at = c(1,0.1*length(loss),0.2*length(loss),0.3*length(loss),0.4*length(loss),0.5*length(loss),0.6*length(loss),0.7*length(loss),0.8*length(loss),0.9*length(loss),length(loss)), labels = round(sort(loss, decreasing = T)[c(15,0.1*length(loss),0.2*length(loss),0.3*length(loss),0.4*length(loss),0.5*length(loss),0.6*length(loss),0.7*length(loss),0.8*length(loss),0.9*length(loss),length(loss))],2))

## Minimum AMSE Method (Caeiro and Gomes, 2016)

K <- c(floor(N^0.99),floor(N^0.995))

M_hill <- function(k,j){
  
  order_loss <- sort(loss, decreasing = T)
  M <- mean((log(order_loss[1:k]) - log(order_loss[k + 1]))^j)
  
  return(M)
  
}

M11 <- M_hill(K[1],1)
M12 <- M_hill(K[1],2)
M13 <- M_hill(K[1],3)

M21 <- M_hill(K[2],1)
M22 <- M_hill(K[2],2)
M23 <- M_hill(K[2],3)

# tau = 1

W11 <- (M11-sqrt(M12/2))/(sqrt(M12/2)-(M13/6)^(1/3))
W21 <- (M21-sqrt(M22/2))/(sqrt(M22/2)-(M23/6)^(1/3))
rho11 <- -abs(3*(W11-1)/(W11-3))
rho21 <- -abs(3*(W21-1)/(W21-3))

# tau = 0

W10 <- (log(M11)-(1/2)*log(M12/2))/((1/2)*log(M12/2)-(1/3)*log(M13/6))
W20 <- (log(M21)-(1/2)*log(M22/2))/((1/2)*log(M22/2)-(1/3)*log(M23/6))
rho10 <- -abs(3*(W10-1)/(W10-3))
rho20 <- -abs(3*(W20-1)/(W20-3))

chi1 <- median(c(rho11,rho21))
chi0 <- median(c(rho10,rho20))

I1 <- c((rho11-chi1)^2,(rho21-chi1)^2)
I0 <- c((rho10-chi0)^2,(rho20-chi0)^2)

if (sum(I0) <= sum(I1)){tau <- 0; rho <- rho20} else {tau <- 1; rho <- rho21}

U <- c()

for (i in 1:K[2]){
  order_loss <- sort(loss, decreasing = TRUE)
  U[i] <- i*(log(order_loss[i]) - log(order_loss[i+1]))
}

i <- 1:K[2]
d <- mean((i/K[2])^(-rho))

D <- function(a){
  
  D <- mean((i/K[2])^(-a)*U[i])
  
  return(D)
}

beta <- (K[2]/N)^rho*(d*D(0)-D(rho))/(d*D(rho)-D(2*rho))

k_0 <- floor((((1-rho)^2*N^(-2*rho))/(-2*rho*beta^2))^(1/(1-2*rho)))

# Threshold

threshold <- sort(loss,decreasing = T)[k_0]

## Kernel-EVT Mixture Method (MacDonald et al., 2011)

# GPD

phi_u <- length(loss[which(loss > threshold)])/length(loss)

A <- length(loss[which(loss < threshold)])

B <- length(loss[which(loss > threshold)])

### GPD - Likelihood Function

func_ML_gpd <- function(coef){
  
  sigma <- coef[1]
  xi <- coef[2]
  
  ln_L = -(B*log(phi_u) + sum(log((1/sigma)*(1 + (xi*(loss[which(loss > threshold)] - threshold)/sigma))^(-1-1/xi))))
  
  return(ln_L)
  
}

sigma.c <- 1

xi.c <- 1

gpd <- nlm(func_ML_gpd,c(sigma.c,xi.c), steptol = 1e-15, gradtol = 1e-15, iterlim = 100000, ndigit = 12, fscale = 0, print.level = 1)
gpd <- optim(c(sigma.c,xi.c),func_ML_gpd)

l_GPD <- gpd$value
GPD_parameters <- gpd$par

hess <- numDeriv:::hessian(func_ML_gpd,GPD_parameters)
ep_gpd <- sqrt(diag(solve(hess)))

### Kernel Density - Likelihood Cross Validation Function

func_ML_kernel <- function(coef){
  
  lambda <- coef[1]
  
  ker_cross <- function(x){
    
    ker_cross <- c()
    ker_cross.aux <- c()
    
    A.aux <- which(x < threshold)
    
    for (i in A.aux){
      
      seq <- seq(1,length(x),1)[-i]
      
      for (j in seq){
        
        ker_cross.aux[j] <- dnorm((x[i] - x[j])/lambda)
        
      }
      
      ker_cross[i] <- (1/(length(x) - 1)) * (1/lambda) * sum(ker_cross.aux[-i])
      
    }
    
    return(ker_cross[A.aux])
    
  }
  
  ker_int <- function(x){
    
    ker_int.aux <- c()
    
    for (i in 1:length(x)){
      
      ker_int.aux[i] <- pnorm((threshold - x[i])/lambda)
      
    }
    
    ker_int <- sum(ker_int.aux) * (1/length(x))
    
    return(ker_int)
    
  }

  ln_L <- -(A * (log(1 - phi_u) - log(ker_int(loss))) + sum(log(ker_cross(loss)))) 
  
  return(ln_L)
  
}

lambda.c <- 1

kden <- optim(lambda.c,func_ML_kernel)
l_kernel <- kden$value
lambda <- kden$par
hess <- numDeriv:::hessian(func_ML_kernel,kden$par)
ep_kden <- sqrt(diag(solve(hess)))

## Extreme Value Kernel Mixture Likelihood - Goodness-of-Fit Measures

KGPD_loglik <- -(l_GPD + l_kernel)
KGPD_aic <- 2*((length(gpd$par)+length(kden$par)) - KGPD_loglik)
KGPD_bic <- ((length(gpd$par)+length(kden$par))* log(length(loss))) - 2*KGPD_loglik

# Histogram

dgpd <- function(x,xi,mu,sigma){
  
  dgpd <- (1/sigma)*((1 + xi*((x-mu)/sigma))^(-1-1/xi))
  
  return(dgpd)
  
}

dkden <- function(x,lambda){
  
  diff <- c()
  
  for (i in 1:length(loss)){diff[i] <- dnorm((x - loss[i])/lambda,0,1)}
  
  dkden <- (1/N)*(1/lambda)*(sum(diff))
  
  return(dkden)
  
}

dkden <- Vectorize(dkden)

pkden <- function(x,lambda){
  
  pkden.aux <- c()
  
  for (i in 1:length(loss)){
    
    pkden.aux[i] <- pnorm((x - loss[i])/lambda)
    
  }
  
  pkden <- sum(pkden.aux) * (1/length(loss))
  
  return(pkden)
  
}

pkden <- Vectorize(pkden)

hist(loss[which(loss < max_lim)], plot=T, prob=T, breaks = 250, main = "", xlab = "Loss Severity", ylab = "Frequency")
curve((phi_u)*dgpd(x, GPD_parameters[2],threshold,GPD_parameters[1]),add=T,lty=1,lwd=2,col=1,xlim=c(threshold,max_lim))
curve((1-phi_u)*(dkden(x,lambda)/pkden(threshold,lambda)),add=T,xlim=c(0,threshold),lty=1,lwd=2,col=1)
segments(threshold,(1-phi_u)*(dkden(threshold,lambda)/ker_int(loss)),threshold,(phi_u)*dgpd(threshold,GPD_parameters[2],threshold,GPD_parameters[1]),lty=1,lwd=2,col=1)
abline(v = c(threshold), lty = 2, lwd = 2, col = "black")
legend("topright", c("Extreme Value Kernel Mixture"), lty= 1,col="black", lwd = 2, bty = "n")

# Distribution Function

pgpd <- function(x,xi,mu,sigma){
  
  pgpd <- 1 - ((1 + xi*((x - mu)/sigma))^(-1/xi))
  
  return(pgpd)
  
}

plot(ecdf(loss[which(loss < max_lim | loss > min_lim)]), do.points=FALSE, verticals=TRUE, xlim=c(min_lim,max_lim), main="", ylab="Emprirical d.f.")
curve((1-phi_u) + phi_u*pgpd(x,GPD_parameters[2],threshold,GPD_parameters[1]),add=T,lty=2,lwd=2,col=1,xlim=c(threshold,max_lim))
curve((1-phi_u)*(pkden(x,lambda)/pkden(threshold,lambda)),add=T,lty=2,lwd=2,col=1,xlim=c(0,threshold))
legend("bottomright", c("Extreme Value Kernel Mixture"), lty=2, col=1, lwd=2, bty = "n")

# Kolmogorov-Smirnov Goodness-of-Fit Test

D.aux <- c()

for(i in 1:length(x)){
  
  if(x < threshold){
  D.aux[i] <- abs(F_hat_emp[i] - (1-phi_u)*(pkden(x[i],lambda)/pkden(threshold,lambda)))
  }
  else{
  D.aux[i] <- abs(F_hat_emp[i] - ((1-phi_u)+phi_u*pgpd(x[i],GPD_parameters[2],threshold,GPD_parameters[1])))  
  }

}

KGPD_ks <- max(D.aux)

# Quantile Function

qgpd <- function(x,xi,mu,sigma){
  
  qgpd <- mu + (sigma/xi)*((((length(loss)/length(loss[which(loss > mu)]))*(1-x))^(-xi))-1)
  
  return(qgpd)
}

qgpd(c(0.9,0.95,0.975,0.99,0.995,0.999,0.9995,0.9999),GPD_parameters[2],threshold,GPD_parameters[1])
seq <- (1:N-0.5)/N
plot(qgpd(seq[which(seq > (1 - phi_u))],GPD_parameters[2],threshold,GPD_parameters[1]), sort(loss)[which(seq > (1 - phi_u))], xlab="Theoretical Quantiles", ylab="Sample Quantiles", ylim = c(0,max(loss)))
abline(0,1)
abline(v=threshold,lty=3)

###################################################### K-S Parametric Bootstrap ####################################################

B <- 10000

## Log-Normal

func_ML_lnorm <- function(coef){
  
  mu <- coef[1]
  sigma <- coef[2]
  
  ln_L = -(-(N/2)*log(2*pi) - (N/2)*log(sigma^2) - (1/(2*sigma^2))* sum((log(sample)-mu)^2) - sum(log(sample)))
  
  return(ln_L)
  
}

mu.c <- 1
sigma.c <- 1

rlnorm <- function(n,mu,sigma){
  
  sample <- runif(n,0,1)
  rlnorm <- sapply(sample,function(x) qlnorm(x,mu,sigma))
  
  return(rlnorm)
  
}

bs_lnorm <- c()

for (i in 1:B) {
  sample <- rlnorm(N,lnorm$par[1],lnorm$par[2])
  m <- optim(c(mu.c,sigma.c),func_ML_lnorm,control = list(reltol = 1e-15))
  bs_lnorm[i] <- as.numeric(ks.test(sample,"plnorm",m$par[1],m$par[2])$statistic)
}

pcks_lnorm <- length(bs_lnorm[which(bs_lnorm >= lnorm_ks)])/B
hist(bs_lnorm,breaks = 250,main="Sampling Distribtuion - Lognormal",prob=T,xlab="K-S Statistics",xlim=c(0,0.2))
abline(v=lnorm_ks,col=2,lwd=2)

## Gamma

func_ML_gamma <- function(coef){
  
  alpha <- coef[1]
  beta <- coef[2]
  
  ln_L = -((alpha-1)*sum(log(sample)) -(1/beta)*sum(sample) - N*(alpha*log(beta) + log(gamma(alpha))))
  
  return(ln_L)
  
}

alpha.c <- 0.5
beta.c <- 0.5

rgamma <- function(n,alpha,beta){
  
  sample <- runif(n,0,1)
  rgamma <- sapply(sample,function(x) qgamma(x,alpha,1/beta))
  
  return(rgamma)
  
}

bs_gamma <- c()

for (i in 1:B) {
  sample <- rgamma(N,gamma$par[1],gamma$par[2])
  m <- optim(c(alpha.c,beta.c),func_ML_gamma,control = list(reltol = 1e-15))
  bs_gamma[i] <- as.numeric(ks.test(sample,"pgamma",m$par[1],1/m$par[2])$statistic)
}

pcks_gamma <- length(bs_gamma[which(bs_gamma >= gamma_ks)])/B
hist(bs_gamma,breaks = 250,main="Sampling Distribtuion - Gamma",prob=T,xlab="K-S Statistics",xlim=c(0,0.25))
abline(v=gamma_ks,col=2,lwd=2)

## Weibull

func_ML_weibull <- function(coef){
  
  alpha <- coef[1]
  beta <- coef[2]
  
  ln_L = -((N*(log(alpha)-alpha*log(beta))) + ((alpha - 1)*sum(log(sample))) - sum((sample/beta)^alpha))
  
  return(ln_L)
  
}

alpha.c <- 1
beta.c <- 30

rweibull <- function(n,alpha,beta){
  
  sample <- runif(n,0,1)
  rweibull <- sapply(sample,function(x) qweibull(x,alpha,beta))
  
  return(rweibull)
  
}

bs_weibull <- c()

for (i in 1:B) {
  
  sample <- rweibull(N,weibull$par[1],weibull$par[2])
  m <- optim(c(alpha.c,beta.c),func_ML_weibull,control = list(reltol = 1e-25))
  bs_weibull[i] <- as.numeric(ks.test(sample,"pweibull",m$par[1],m$par[2])$statistic)
  
}

pcks_weibull <- length(bs_weibull[which(bs_weibull >= weibull_ks)])/B
hist(bs_weibull,breaks = 250,main="Sampling Distribtuion - Weibull",prob=T,xlab="K-S Statistics",xlim=c(0,0.3))
abline(v=weibull_ks,col=2,lwd=2)

## Pareto

func_ML_pareto <- function(coef){
  
  a <- coef[1]
  c <- min(sample)
  
  ln_L = -((N * log(a)) + (N * a * log(c)) - (a + 1) * sum(log(sample)))
  
  return(ln_L)
  
}

a.c <- 0.1

rpareto <- function(n,a,c){
  
  sample <- runif(n,0,1)
  rpareto <- sapply(sample,function(x) qpareto(x,a,c))
  
  return(rpareto)
  
}

bs_pareto <- c()

options(warn=-1)

for (i in 1:B){
  sample <- rpareto(N,pareto$par,min(loss))
  m <- optim(a.c,func_ML_pareto,control = list(reltol = 1e-15))
  bs_pareto[i] <- as.numeric(ks.test(sample,"ppareto",m$par[1],min(sample))$statistic)
}

pcks_pareto <- length(bs_pareto[which(bs_pareto >= pareto_ks)])/B
hist(bs_pareto,breaks = 250,main="Sampling Distribtuion - Pareto",prob=T,xlab="K-S Statistics",xlim=c(0,0.5))
abline(v=pareto_ks,col=2,lwd=2)

## Inverse Gaussian

func_ML_invnorm <- function(coef){
  
  mu <- coef[1]
  lambda <- coef[2]
  
  ln_L = -((N/2)*log(lambda/(2*pi)) - (3/2)*sum(log(sample)) - lambda*sum(((sample-mu)^2)/(2*(mu^2)*sample)))
  
  return(ln_L)
  
}

rinvgauss <- function(n,mu,lambda){
  
  sample <- runif(n,0,1)
  rinvgauss <- sapply(sample,function(x) actuar:::qinvgauss(x,mu,lambda))
  
  return(rinvgauss)
  
}

bs_invnorm <- c()

for (i in 1:B){
  
  sample <- rinvgauss(N,invnorm$par[1],invnorm$par[2])
  m <- optim(c(mean(loss),1/((sum(1/loss)/N) - (1/mean(loss)))),func_ML_invnorm)
  bs_invnorm[i] <- as.numeric(ks.test(sample,"pinvgauss",m$par[1],m$par[2])$statistic)
  
}

pcks_invnorm <- length(bs_invnorm[which(bs_invnorm >= invnorm_ks)])/B
hist(bs_invnorm,breaks = 250,main="Sampling Distribtuion - Inverse Gaussian",prob=T,xlab="K-S Statistics",xlim=c(0,0.3))
abline(v=invnorm_ks,col=2,lwd=2)

## Inverse Gamma

func_ML_invgamma <- function(coef){
  
  alpha <- coef[1]
  beta <- coef[2]
  
  ln_L = -(-N*(alpha+1)*mean(log(sample)) - N*log(gamma(alpha)) + N*alpha*log(beta) - sum(beta/sample))
  
  return(ln_L)
  
}

rinvgamma <- function(n,alpha,beta){
  
  sample <- runif(n,0,1)
  rinvgamma <- sapply(sample,function(x) actuar:::qinvgamma(x,alpha,1/beta))
  
  return(rinvgamma)
  
}

bs_invgamma <- c()

for (i in 1:B){
  sample <- rinvgamma(N,invgamma$par[1],invgamma$par[2])
  m <- optim(c(1,5),func_ML_invgamma)
  bs_invgamma[i] <- as.numeric(ks.test(sample,"pinvgamma",m$par[1],m$par[2])$statistic)
}

pcks_invgamma <- length(bs_invgamma[which(bs_invgamma >= invgamma_ks)])/B
hist(bs_invgamma,breaks=250,main="Sampling Distribtuion - Inverse Gamma",prob=T,xlab="K-S Statistics",xlim=c(0,0.25))
abline(v=invgamma_ks,col=2,lwd=2)

## Burr Type XII

func_ML_burr <- function(coef){
  
  lambda <- coef[1]
  k <- coef[2]
  alpha <- coef[3]
  
  ln_L = -(N*log(k) + N*log(lambda) + (lambda-1)*sum(log(sample)) - N*lambda*log(alpha) - (1+k)*sum(log(1+(sample/alpha)^lambda)))
  
  return(ln_L)
  
}

rburr <- function(n,lambda,k,alpha){
  
  sample <- runif(n,0,1)
  rburr <- sapply(sample,function(x) qburr(x,lambda,k,alpha))
  
  return(rburr)
  
}

bs_burr <- c()

for (i in 1:B) {
  sample <- rburr(N,burr$par[1],burr$par[2],burr$par[3])
  m <- optim(c(20,1,1),func_ML_burr,control = list(reltol = 1e-25))
  bs_burr[i] <- as.numeric(ks.test(sample,"pburr",m$par[1],m$par[2],m$par[3])$statistic)
}

pcks_burr <- length(bs_burr[which(bs_burr >= burr_ks)])/B
hist(bs_burr,breaks=250,main="Sampling Distribtuion - Burr Type XII",prob=T,col=1,xlab="K-S Statistics",xlim=c(0,0.1))
abline(v=burr_ks,col=2,lwd=2)

## Generalized Lambda

bs_glambda <- c()

for (i in i:B){
  sample <- fBasics:::rgld(N,as.numeric(glambda@fit$estimate[1]),as.numeric(glambda@fit$estimate[2]),as.numeric(glambda@fit$estimate[3]),as.numeric(glambda@fit$estimate[4]))
  m <- fBasics:::gldFit(sample,lambda1=1,lambda2=-0.4,lambda3=-5,lambda4=-1,method="mle",doplot=F,trace=F)
  bs_glambda[i] <- as.numeric(ks.test(sample,"pgld",lambda1=m@fit$estimate[1],lambda2=m@fit$estimate[2],lambda3=m@fit$estimate[3],lambda4=m@fit$estimate[4])$statistic)
}

pcks_glambda <- length(bs_glambda[which(bs_glambda >= glambda_ks)])/i
hist(bs_glambda,breaks=550,main="Sampling Distribtuion - Generalized Lambda",prob=T,xlab="K-S Statistics",col=1,xlim=c(0,0.5))
abline(v=glambda_ks,col=2,lwd=2)

## Gamma Mixture

library(evmix)

bs_mixgamma <- c()

for (i in 1:B){
  sample <- evmix:::rmgamma(N,mgshape=c(mixgamma$gamma.pars[1,]),mgscale=c(mixgamma$gamma.pars[2,]),mgweight=c(mixgamma$lambda))
  m <- mixtools:::gammamixEM(loss,k=4,maxit=10000)
  bs_mixgamma[i] <- as.numeric(ks.test(sample,"pmgamma",mgshape=c(mixgamma$gamma.pars[1,]),mgscale=c(mixgamma$gamma.pars[2,]),mgweight=c(mixgamma$lambda))$statistic)
}

pcks_mixgamma <- length(bs_mixgamma[which(bs_mixgamma >= mixgamma_ks)])/B
hist(bs_mixgamma,breaks=250,main="Sampling Distribtuion - Gamma Mixture",prob=T,col=1,xlab="K-S Statistics",xlim=c(0,0.1))
abline(v=mixgamma_ks,col=2,lwd=2)

### Results

# Table

distribution <- c("Lognormal","Gamma","Weibull","Pareto","Inverse Gaussian","Inverse Gamma","Burr Type XII","Generalized Lambda","Gamma Mixture","Extreme Value Kernel Mixture")
NLL <- -as.numeric(rbind(lnorm_loglik,gamma_loglik,weibull_loglik,pareto_loglik,invnorm_loglik,invgamma_loglik,burr_loglik,glambda_loglik,mixgamma_loglik,KGPD_loglik))
aic <- as.numeric(rbind(lnorm_aic,gamma_aic,weibull_aic,pareto_aic,invnorm_aic,invgamma_aic,burr_aic,glambda_aic,mixgamma_aic,KGPD_aic))
bic <- as.numeric(rbind(lnorm_bic,gamma_bic,weibull_bic,pareto_bic,invnorm_bic,invgamma_bic,burr_bic,glambda_bic,mixgamma_bic,KGPD_bic))
delta_aic <- aic - min(aic)
delta_bic <- bic - min(bic)
K_S <- as.numeric(rbind(lnorm_ks,gamma_ks,weibull_ks,pareto_ks,invnorm_ks,invgamma_ks,burr_ks,glambda_ks,mixgamma_ks,KGPD_ks))
table <- cbind(distribution,NLL,aic,bic,delta_aic,delta_bic,K_S)

# Quantile Plots

par(mfrow=c(3,3))

plot(qlnorm((1:N-0.5)/N, lnorm$par[1],lnorm$par[2]), sort(loss), main = "Lognormal", xlab="Theoretical Quantiles", ylab="Sample Quantiles", ylim = c(0,max(loss)))
abline(0,1)

plot(qgamma((1:N-0.5)/N, gamma$par[1], 1/gamma$par[2]), sort(loss), main = "Gamma" ,xlab = "Theoretical Quantiles", ylab="Sample Quantiles", ylim = c(0,max(loss)))
abline(0,1)

plot(qweibull((1:N-0.5)/N, weibull$par[1], weibull$par[2]), sort(loss), main = "Weibull" ,xlab = "Theoretical Quantiles", ylab="Sample Quantiles", ylim = c(0,max(loss)))
abline(0,1)

plot(qpareto((1:N-0.5)/N, a.pareto, c.pareto), sort(loss), main = "Pareto" ,xlab = "Theoretical Quantiles", ylab="Sample Quantiles", ylim = c(0,max(loss)))
abline(0,1)

plot(actuar:::qinvgauss((1:N-0.5)/N, invnorm$par[1], invnorm$par[2]), sort(loss), main = "Inverse Gaussian" ,xlab = "Theoretical Quantiles", ylab="Sample Quantiles", ylim = c(0,max(loss)))
abline(0,1)

plot(actuar:::qinvgamma((1:N-0.5)/N, invgamma$par[1], 1/invgamma$par[2]), sort(loss), main = "Inverse Gamma" , xlab = "Theoretical Quantiles", ylab="Sample Quantiles", ylim = c(0,max(loss)))
abline(0,1)

plot(qburr((1:N-0.5)/N, burr$par[1], burr$par[2], burr$par[3]), sort(loss), main = "Burr Type XII" ,xlab = "Theoretical Quantiles", ylab="Sample Quantiles", ylim = c(0,max(loss)))
abline(0,1)

plot(qgld((1:N-0.5)/N, lambda1 = glambda@fit$estimate[1], lambda2 = glambda@fit$estimate[2], lambda3 = glambda@fit$estimate[3], lambda4 = glambda@fit$estimate[4]), sort(loss), main = "Generalized Lambda" ,xlab = "Theoretical Quantiles", ylab="Sample Quantiles", ylim = c(0,max(loss)))
abline(0,1)

plot(evmix:::qmgamma((1:N-0.5)/N, mgshape = c(mixgamma$gamma.pars[1,]), mgscale = c(mixgamma$gamma.pars[2,]), mgweight = c(mixgamma$lambda)), sort(loss), main = "Mixture Gamma" ,xlab = "Theoretical Quantiles", ylab="Sample Quantiles", ylim = c(0,max(loss)))
abline(0,1)

########################################################### Robust Inference #########################################################

## Method of Probability-Weighted Moments (Hosking and Wallis, 1987)

B.aux <- which(loss > threshold)

B <- length(B.aux)
  
a_r <- function(r){

  order_loss <- sort(loss[B.aux])
  
  gamma <- -0.35
  delta <- 0
  
  a_r.aux <- c()
  
  for(j in (1:length(order_loss))){
    
    a_r.aux[j] <- ((1-((j+gamma)/(B+delta)))^r)*(order_loss[j]-threshold)
    
  }
  
  a_r <- (1/length(order_loss))*sum(a_r.aux)
  
  return(a_r)
  
}

GPD_parameters_PWM <- c()
GPD_parameters_PWM[1] <- (2*a_r(0)*a_r(1))/(a_r(0)-2*a_r(1))
GPD_parameters_PWM[2] <- 2 - (a_r(0)/(a_r(0)-2*a_r(1)))

if(GPD_parameters_PWM[2] > 0.5){denom <- NA 
warning("Asymptotic standard errors not available for"," PWM Method when xi > 0.5")} else{denom <- B*(1 - 2*GPD_parameters_PWM[2])*(3 - 2*GPD_parameters_PWM[2])

var1 <- (1 - GPD_parameters_PWM[2])*(1 - GPD_parameters_PWM[2] + 2*GPD_parameters_PWM[2]^2)*(2 - GPD_parameters_PWM[2])^2
var2 <- (7 - 18*GPD_parameters_PWM[2] + 11*GPD_parameters_PWM[2]^2 - 2*GPD_parameters_PWM[2]^3)*GPD_parameters_PWM[1]^2
cov <- GPD_parameters_PWM[1]*(2 - GPD_parameters_PWM[2])*(2 - 6*GPD_parameters_PWM[2] + 7*GPD_parameters_PWM[2]^2 - 2*GPD_parameters_PWM[2]^3)
varcov_PWM <- matrix(c(var2,cov,cov,var1),2)/denom
ep_gpd_PWM <- sqrt(diag(varcov_PWM))
}

## Method of Trimmed Moments (Brazauskas and Kleefeld, 2009)

order_loss <- sort(loss[B.aux])

a1 <- 0.10
b1 <- 0.55

a2 <- 0.70
b2 <- 0.05

LBP <- min(a1,a2)
UBP <- min(b1,b2)

order_loss[floor(B*LBP)]
order_loss[floor(B*(1-UBP))]

mu_hat1 <- (1/(B-floor(a1*B)-floor(b1*B)))*sum(order_loss[(floor(a1*B)+1):(B-floor(b1*B))]-threshold)
mu_hat2 <- (1/(B-floor(a2*B)-floor(b2*B)))*sum(order_loss[(floor(a2*B)+1):(B-floor(b2*B))]-threshold)

func_MTM_gpd <- function(coef){

xi <- coef[1]
sigma <- mu_hat1*xi*(1/(1-((((1-a1)^(xi+1))-(b1^(xi+1)))/((xi+1)*(1-a1-b1)))))
  
mu1 <- sigma*(1/xi)*(1-((((1-a1)^(xi+1))-(b1^(xi+1)))/((xi+1)*(1-a1-b1))))
mu2 <- sigma*(1/xi)*(1-((((1-a2)^(xi+1))-(b2^(xi+1)))/((xi+1)*(1-a2-b2))))

MTM <- (mu_hat1/mu_hat2) - (mu1/mu2)

return(MTM)

}

MTM_gpd <- uniroot(func_MTM_gpd,interval=c(-2,2))

GPD_parameters_MTM <- c()
GPD_parameters_MTM[1] <- (mu_hat1*MTM_gpd$root*(1/(1-((((1-a1)^(MTM_gpd$root+1))-(b1^(MTM_gpd$root+1)))/((MTM_gpd$root+1)*(1-a1-b1))))))
GPD_parameters_MTM[2] <- MTM_gpd$root

mu1 <- GPD_parameters_MTM[1]*(1/GPD_parameters_MTM[2])*(1-((((1-a1)^(GPD_parameters_MTM[2]+1))-(b1^(GPD_parameters_MTM[2]+1)))/((GPD_parameters_MTM[2]+1)*(1-a1-b1))))
mu2 <- GPD_parameters_MTM[1]*(1/GPD_parameters_MTM[2])*(1-((((1-a2)^(GPD_parameters_MTM[2]+1))-(b2^(GPD_parameters_MTM[2]+1)))/((GPD_parameters_MTM[2]+1)*(1-a2-b2))))

dmu1 <- (-GPD_parameters_MTM[1])*(1/(GPD_parameters_MTM[2]*(GPD_parameters_MTM[2]+1)))*((((2*GPD_parameters_MTM[2])+1)*(mu1/GPD_parameters_MTM[1]))-1+(((((1-a1)^(GPD_parameters_MTM[2]+1))*(log(1-a1)))-((b1^(GPD_parameters_MTM[2]+1))*log(b1)))/(1-a1-b1)))
dmu2 <- (-GPD_parameters_MTM[1])*(1/(GPD_parameters_MTM[2]*(GPD_parameters_MTM[2]+1)))*((((2*GPD_parameters_MTM[2])+1)*(mu2/GPD_parameters_MTM[1]))-1+(((((1-a2)^(GPD_parameters_MTM[2]+1))*(log(1-a2)))-((b2^(GPD_parameters_MTM[2]+1))*log(b2)))/(1-a2-b2)))
d <- ((dmu1*mu2)-(dmu2*mu1))

D11 <- GPD_parameters_MTM[1]*((d-(dmu1*mu2))/(mu1*d))
D12 <- GPD_parameters_MTM[1]*(dmu1/d)
D21 <- mu2/d
D22 <- -mu1/d

D <- matrix(c(D11,D21,D12,D22),2)

composite.trapezoid <- function(a,b,c,d,k){

  hx <- ((1-b)-a)/k
  hy <- ((1-d)-c)/k
  i <- 0:k
  j <- 0:k
  xi <- a+(i*hx)
  yj <- c+(j*hy)
  
  g <- function(x,y){min(x,y)-(x*y)}
  f2 <- function(x){evir:::qgpd(x,-GPD_parameters_MTM[2],threshold,GPD_parameters_MTM[1])}

  I_k <- matrix(0,k+1,k+1)
  
  for(l in 1:(k+1)){
    
    for(m in 1:(k+1)){
      
      wx <- if(l==1 || l==(k+1)){wx = hx*(1/2)} else {wx = hx}
      wy <- if(m==1 || m==(k+1)){wy = hy*(1/2)} else {wy = hy}
      
      I_k[l,m] <- g(xi[l],yj[m])*numDeriv:::grad(f2,xi[l])*numDeriv:::grad(f2,yj[m])*wx*wy
      
    }
  }
  return(sum(I_k))
}

S11 <- (1/((1-a1-b1)*(1-a1-b1)))*composite.trapezoid(a1,b1,a1,b1,500)
S12 <- (1/((1-a1-b1)*(1-a2-b2)))*composite.trapezoid(a1,b1,a2,b2,500)
S21 <- (1/((1-a2-b2)*(1-a1-b1)))*composite.trapezoid(a2,b2,a1,b1,500)
S22 <- (1/((1-a2-b2)*(1-a2-b2)))*composite.trapezoid(a2,b2,a2,b2,500)

S <- matrix(c(S11,S21,S12,S22),2)

varcov_MTM <- (D%*%S%*%t(D))/B
ep_gpd_MTM <- sqrt(diag(varcov_MTM))

GPD_parameters_MTM[2] <- -GPD_parameters_MTM[2]

## Method of Medians (Peng and Welsh, 2001)

func_MM_gpd <- function(coef){

xi <- coef[1]
sigma <- (xi/((2^xi)-1))*median(order_loss-threshold)

func_MM_gpd.aux <- function(coef){

y1 <- coef[1]
  
z <- ((-log(y1)/xi)-(((1+xi)/(xi^2))*(1-(y1^xi)))) - ((-log(y1+0.5)/xi)-(((1+xi)/(xi^2))*(1-((y1+0.5)^xi))))

return(z)

}

y1 <- uniroot(func_MM_gpd.aux,interval=c(0.1,2))$root

MM <- median((log(1+(xi*(order_loss-threshold)/sigma))/(xi^2))-(((1+xi)*(order_loss-threshold))/((sigma*xi)+((xi^2)*(order_loss-threshold))))) - ((-log(y1)/xi)-(((1+xi)/(xi^2))*(1-(y1^xi))))

return(MM)
  
}

MM_gpd <- uniroot(func_MM_gpd,interval=c(0.1,2))

GPD_parameters_MM <- c()
GPD_parameters_MM[1] <- (MM_gpd$root/((2^MM_gpd$root)-1))*median(order_loss-threshold)
GPD_parameters_MM[2] <- MM_gpd$root

func_MM_gpd.aux <- function(coef){
  
  y1 <- coef[1]
  
  z <- ((-log(y1)/GPD_parameters_MM[2])-(((1+GPD_parameters_MM[2])/(GPD_parameters_MM[2]^2))*(1-(y1^GPD_parameters_MM[2])))) - ((-log(y1+0.5)/GPD_parameters_MM[2])-(((1+GPD_parameters_MM[2])/(GPD_parameters_MM[2]^2))*(1-((y1+0.5)^GPD_parameters_MM[2]))))
  
  return(z)
  
}

y1 <- uniroot(func_MM_gpd.aux,interval=c(0.1,2))$root

v22 <- -((2*y1*log(y1))/GPD_parameters_MM[2])+((2*(y1+0.5)*log(y1+0.5))/GPD_parameters_MM[2])+((2*(y1^(GPD_parameters_MM[2]+1)))/(GPD_parameters_MM[2]^2))-((2*((y1+0.5)^(GPD_parameters_MM[2]+1)))/(GPD_parameters_MM[2]^2))+(1/(GPD_parameters_MM[2]^2))
v21 <- ((2/(GPD_parameters_MM[1]*GPD_parameters_MM[2]))*(((y1+0.5)^(GPD_parameters_MM[2]+1))-(y1^(GPD_parameters_MM[2]+1))))-(1/(GPD_parameters_MM[1]*GPD_parameters_MM[2]))
v12 <- (((2^(-GPD_parameters_MM[2]))+((GPD_parameters_MM[2]*log(2))-1))/(GPD_parameters_MM[2]^2))*sign(GPD_parameters_MM[2]+1)
v11 <- (((2^(-GPD_parameters_MM[2]))-1)/(GPD_parameters_MM[2]*GPD_parameters_MM[1]))*sign(GPD_parameters_MM[2]+1)
V <- matrix(c(v11,v21,v12,v22),2)

b11 <- 1
b12 <- 0
b21 <- b12
b22 <- 1
C <- matrix(c(b11,b21,b12,b22),2)

varcov_MM <- (solve(V)%*%C%*%solve(t(V)))/B
ep_gpd_MM <- sqrt(diag(varcov_MM))

## Method of Minimum Density Power Divergence Estimator (Juaréz and Schucany, 2004)

a <- 0.10

H_a <- function(coef){
  
  sigma <- coef[1]
  xi <- coef[2]
  
  H_a <- (1/((sigma^a)*(1+a-a*xi))) - (1+(1/a))*(1/B)*sum((1/(sigma^a))*((1-(xi*(order_loss-threshold)/sigma))^(((1/xi)-1)*a)))
  
}

MDPDE_gpd <- optim(c(150,-0.3),H_a)

GPD_parameters_MDPDE <- c()
GPD_parameters_MDPDE[1] <- MDPDE_gpd$par[1]
GPD_parameters_MDPDE[2] <- MDPDE_gpd$par[2]

## Density Function

g <- function(x){
  
  g <- (1/GPD_parameters_MDPDE[1])*((1-(GPD_parameters_MDPDE[2]*x/GPD_parameters_MDPDE[1]))^((1/GPD_parameters_MDPDE[2])-1))
  
  return(g)
  
}

## Score Function Vector

S <- function(x){
  
 S <- matrix(0,2,1)

 S[1] <- ((-1/(GPD_parameters_MDPDE[2]^2))*log(1-(GPD_parameters_MDPDE[2]*x/GPD_parameters_MDPDE[1])))+
  ((1/GPD_parameters_MDPDE[2])*((1/GPD_parameters_MDPDE[2])-1)*
  (1-((1-(GPD_parameters_MDPDE[2]*x/GPD_parameters_MDPDE[1]))^(-1))))
 
 S[2] <- (-1/(GPD_parameters_MDPDE[1]*GPD_parameters_MDPDE[2]))+((1/GPD_parameters_MDPDE[1])*
  ((1/GPD_parameters_MDPDE[2])-1)*((1-(GPD_parameters_MDPDE[2]*x/GPD_parameters_MDPDE[1]))^(-1))) 
 
return(S)

}

## Information Matrix

I <- function(x){
  
  I <- matrix(0,2,2)
  
  I[1,1] <- ((-2/(GPD_parameters_MDPDE[2]^3))*log(1-(GPD_parameters_MDPDE[2]*x/GPD_parameters_MDPDE[1])))+
    ((3-GPD_parameters_MDPDE[2])/(GPD_parameters_MDPDE[2]^3))-(((2*(2-GPD_parameters_MDPDE[2]))/(GPD_parameters_MDPDE[2]^3))*
    ((1-(GPD_parameters_MDPDE[2]*x/GPD_parameters_MDPDE[1]))^(-1)))+(((1-GPD_parameters_MDPDE[2])/(GPD_parameters_MDPDE[2]^3))*
    ((1-(GPD_parameters_MDPDE[2]*x/GPD_parameters_MDPDE[1]))^(-2)))
  
  I[1,2] <- (-1/(GPD_parameters_MDPDE[1]*(GPD_parameters_MDPDE[2]^2)))+(((2-GPD_parameters_MDPDE[2])/(GPD_parameters_MDPDE[1]*(GPD_parameters_MDPDE[2]^2)))*
    ((1-(GPD_parameters_MDPDE[2]*x/GPD_parameters_MDPDE[1]))^(-1)))-(((1-GPD_parameters_MDPDE[2])/(GPD_parameters_MDPDE[1]*(GPD_parameters_MDPDE[2]^2)))*
    ((1-(GPD_parameters_MDPDE[2]*x/GPD_parameters_MDPDE[1]))^(-2)))
  
  I[2,1] <- I[1,2]
  
  I[2,2] <- (-1/(GPD_parameters_MDPDE[2]*(GPD_parameters_MDPDE[1]^2)))+((1/(GPD_parameters_MDPDE[1]^2))*
    ((1/GPD_parameters_MDPDE[2])-1)*((1-(GPD_parameters_MDPDE[2]*x/GPD_parameters_MDPDE[1]))^(-2)))
  
  return(I)
  
}

## Matrix K

S_1 <- matrix(0,2,2)
S_2 <- matrix(0,2,1)

for (i in (1:B)){

S.aux_1 <- (S(order_loss[i]-threshold)%*%t(S(order_loss[i]-threshold)))*(g(order_loss[i]-threshold)^(2*a))
S.aux_2 <- S(order_loss[i]-threshold)*(g(order_loss[i]-threshold)^a)

S_1 <- S_1 + S.aux_1
S_2 <- S_2 + S.aux_2

}

K_1 <- matrix(0,2,2)
K_2 <- matrix(0,2,2)

K_1 <- (1/B)*S_1
K_2 <- (1/(B^2))*(S_2%*%t(S_2))

K <- K_1-K_2

## Matrix J

#### J_1

J_1 <- matrix(0,2,2)

S_gamma <- function(x){((-1/(GPD_parameters_MDPDE[2]^2))*log(1-(GPD_parameters_MDPDE[2]*x/GPD_parameters_MDPDE[1])))+
    ((1/GPD_parameters_MDPDE[2])*((1/GPD_parameters_MDPDE[2])-1)*
    (1-((1-(GPD_parameters_MDPDE[2]*x/GPD_parameters_MDPDE[1]))^(-1))))}

J_1.aux.1 <- function(x){(S_gamma(x)^2)*(g(x)^(a+1))}

J_1[1,1] <- pracma:::integral(J_1.aux.1,xmin=0,xmax=Inf,reltol = 1e-50000)

S_sigma <- function(x){(-1/(GPD_parameters_MDPDE[1]*GPD_parameters_MDPDE[2]))+((1/GPD_parameters_MDPDE[1])*
    ((1/GPD_parameters_MDPDE[2])-1)*((1-(GPD_parameters_MDPDE[2]*x/GPD_parameters_MDPDE[1]))^(-1)))}

J_1.aux.2 <- function(x){S_gamma(x)*S_sigma(x)*(g(x)^(a+1))}

J_1[1,2] <- pracma:::integral(J_1.aux.2,xmin=0,xmax=Inf,reltol = 1e-50000)  

J_1[2,1] <- J_1[1,2]

J_1.aux.3 <- function(x){(S_sigma(x)^2)*(g(x)^(a+1))}

J_1[2,2] <- pracma:::integral(J_1.aux.3,xmin=0,xmax=Inf,reltol = 1e-50000)

J_1 <- (1+a)*J_1

#### J_2

J_2 <- matrix(0,2,2)

I_gamma <- function(x){((-2/(GPD_parameters_MDPDE[2]^3))*log(1-(GPD_parameters_MDPDE[2]*x/GPD_parameters_MDPDE[1])))+
    ((3-GPD_parameters_MDPDE[2])/(GPD_parameters_MDPDE[2]^3))-(((2*(2-GPD_parameters_MDPDE[2]))/(GPD_parameters_MDPDE[2]^3))*
    ((1-(GPD_parameters_MDPDE[2]*x/GPD_parameters_MDPDE[1]))^(-1)))+(((1-GPD_parameters_MDPDE[2])/(GPD_parameters_MDPDE[2]^3))*
    ((1-(GPD_parameters_MDPDE[2]*x/GPD_parameters_MDPDE[1]))^(-2)))}

J_2.aux.1 <- function(x){I_gamma(x)*(g(x)^(a+1))}

J_2[1,1] <- pracma:::integral(J_2.aux.1,xmin=0,xmax=Inf,reltol = 1e-50000)

I_gamma_sigma <- function(x){(-1/(GPD_parameters_MDPDE[1]*(GPD_parameters_MDPDE[2]^2)))+(((2-GPD_parameters_MDPDE[2])/(GPD_parameters_MDPDE[1]*(GPD_parameters_MDPDE[2]^2)))*
                 ((1-(GPD_parameters_MDPDE[2]*x/GPD_parameters_MDPDE[1]))^(-1)))-(((1-GPD_parameters_MDPDE[2])/(GPD_parameters_MDPDE[1]*(GPD_parameters_MDPDE[2]^2)))*
                 ((1-(GPD_parameters_MDPDE[2]*x/GPD_parameters_MDPDE[1]))^(-2)))}

J_2.aux.2 <- function(x){I_gamma_sigma(x)*(g(x)^(a+1))}

J_2[1,2] <- pracma:::integral(J_2.aux.2,xmin=0,xmax=Inf,reltol = 1e-50000)

J_2[2,1] <- J_2[1,2]

I_sigma <- function(x){(-1/(GPD_parameters_MDPDE[2]*(GPD_parameters_MDPDE[1]^2)))+((1/(GPD_parameters_MDPDE[1]^2))*
           ((1/GPD_parameters_MDPDE[2])-1)*((1-(GPD_parameters_MDPDE[2]*x/GPD_parameters_MDPDE[1]))^(-2)))}

J_2.aux.3 <- function(x){I_sigma(x)*(g(x)^(a+1))}

J_2[2,2] <- pracma:::integral(J_2.aux.3,xmin=0,xmax=Inf,reltol = 1e-50000)

#### J_3

J_3.aux <- matrix(0,2,2)
J_3 <- matrix(0,2,2)

for (i in 1:B){

J_3.aux <- (I(order_loss[i]-threshold)-a*(S(order_loss[i]-threshold)%*%t(S(order_loss[i]-threshold))))*(g(order_loss[i]-threshold)^a)
J_3 <- J_3 + J_3.aux

}

J_3 <- (1/B)*J_3

J <- J_1-J_2+J_3

varcov_MDPDE <- (1/B)*(solve(J)%*%K%*%solve(J))
ep_gpd_MDPDE <- c(sqrt(diag(varcov_MDPDE))[2],sqrt(diag(varcov_MDPDE))[1])

GPD_parameters_MDPDE[2] <- -GPD_parameters_MDPDE[2]
varcov_MDPDE.aux <- varcov_MDPDE[2,2]
varcov_MDPDE[2,2] <- varcov_MDPDE[1,1]
varcov_MDPDE[1,1] <- varcov_MDPDE.aux

## Asymptotic Relative Efficient (ARE)

varcov_MLE <- solve(numDeriv:::hessian(func_ML_gpd,GPD_parameters))

ARE_PWM_MLE <- (det(varcov_MLE)/det(varcov_PWM))^(1/2)
ARE_MM_MLE <- (det(varcov_MLE)/det(varcov_MM))^(1/2)
ARE_MTM_MLE <- (det(varcov_MLE)/det(varcov_MTM))^(1/2)
ARE_MDPDE_MLE <- (det(varcov_MLE)/det(varcov_MDPDE))^(1/2)

## Trimmed Mean Absolute Deviation (tMAD)

b_n <- c()
  
for (j in 1:B){
  
  parameters <- GPD_parameters_MDPDE
  
  b_n[j] <- abs(sort(order_loss)[j] - evir:::qgpd((j-0.5)/B,parameters[2],threshold,parameters[1]))
  
}

b_n <- sort(b_n)

delta <- c(0.5,0.75,0.9,0.95,1)
tMAD <- function(delta){sum(b_n[1:floor(B*delta)])/(floor(B*delta))}
sapply(delta,tMAD)

## Percentile-Residual Plot (PR Plot)

seq <- seq(1:B)
s_n.aux <- c()
s_n <- c()
R_n <- c()
percentile <- c()

for (j in 1:B){

  varcov <- varcov_MTM
  parameters <- GPD_parameters_MTM
  p <- (j-0.5)/B
  
  s_n.aux[j] <- sort(order_loss)[j] - evir:::qgpd((j-0.5)/B,parameters[2],threshold,parameters[1])
  
  D <- c((((1-p)^(-parameters[2]))-1)/parameters[2],(parameters[1]*((1-p)^(-parameters[2]))/(parameters[2]^2))*(((1-p)^parameters[2])-(parameters[2]*log(1-p))-1))
    
  s_n[j] <- sqrt(t(D)%*%varcov%*%D)
  
  percentile[j] <- (j-0.5)/B
  
  R_n[j] <- s_n.aux[j]/s_n[j]
  
}

plot(R_n~percentile,type="p",ylab="Standardized Residuals",xlab="Empirical Percentile Levels",main=expression(bold('MDPDE'[2])),ylim=c(-5,5))
abline(h=2.5,lty=2,col=2,lwd=2)
abline(h=-2.5,lty=2,col=2,lwd=2)
abline(h=0,lty=3,lwd=2)
