# Moment Estimator

M_hill <- function(k){
  
  order_loss <- sort(loss, decreasing = T)
  
  M <- mean((log(order_loss[1:k]) - log(order_loss[k + 1]))^2)
  
  return(M)
  
}

M <- sapply(1:(N-1),M_hill)

## Guillou and Hall (2001)

Q <- function(k){
  
  U <- c()
  aux <- c()
  
  T_hill <- function(k){
    
    order_loss <- sort(loss,decreasing = TRUE)
    
    for (i in 1:k){
      
      U[i] <- i*(log(order_loss[i]) - log(order_loss[i+1]))
      aux[i] <- (k-2*i+1)*U[i]
      
    }
    
    T_hill <- sqrt(3/k^3)*(sum(aux)/mean(U))
    T_hill
  }
  
  start <- k-floor(k/2)
  end <- k+floor(k/2)
  y <- start:end
  x <- sapply(y,T_hill)
  erg <- 1/(2*floor(k/2)+1)*sum(x^2)
  erg2 <- sqrt(erg)
  crit <- 1.25
  as.numeric(erg2>=crit)
}

kmax <- floor(N/1.5)

k_0 <- 1

while (Q(k_0)==0){
  
  k_0 <- k_0+1
  
  if (k_0==kmax) break
}

## Two-step bootstrap procedure (Danielsson et al., 2001)

grid_boot <- seq(100,2000,100)

k_1 <- c()
k_2 <- c()
R <- c()

for (i in grid_boot){
  
  n1 <- i
  n2 <- round((n1^2)/N)
  B <- 1000
  
  Q_boot_1 <- function (k) { 
    
    order_loss <- sort(boot_sample_1, decreasing = T)
    M <- mean((log(order_loss[1:k]) - log(order_loss[k + 1]))^2)
    H <- mean((log(order_loss[1:k]) - log(order_loss[k + 1])))
    Q <- M - (2*H^2)
    
    return(Q)
    
  }
  
  Q_boot_2 <- function (k) { 
    
    order_loss <- sort(boot_sample_2, decreasing = T)
    M <- mean((log(order_loss[1:k]) - log(order_loss[k + 1]))^2)
    H <- mean((log(order_loss[1:k]) - log(order_loss[k + 1])))
    Q <- M - (2*H^2)
    
    return(Q)
    
  }
  
  Q_1 <- matrix(nrow = B, ncol = (n1-1))
  Q_2 <- matrix(nrow = B, ncol = (n2-1))
  
  for (l in 1:B){
    
    boot_sample_1 <- sample(loss, n1, replace=T)
    boot_sample_2 <- sample(loss, n2, replace=T)
    
    Q_1[l,] <- sapply(1:(n1-1),Q_boot_1)
    Q_2[l,] <- sapply(1:(n2-1),Q_boot_2)
    
  }
  
  Q_1 <- Q_1^2
  Q_2 <- Q_2^2
  
  Q_1 <- colMeans(Q_1)
  Q_2 <- colMeans(Q_2)
  
  k_1_min <- which.min(Q_1)
  k_2_min <- which.min(Q_2)
  
}

k_0 <- floor(((k_1_min^2)/k_2_min)*(((log(k_1_min)^2)/((2*log(n1) - log(k_1_min))^2))^((log(n1) - log(k_1_min))/log(n1))))
rho <- log(k_1_min)/(-2*log(n1) + 2*log(k_1_min))