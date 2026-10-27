
geraRegGamma<-function(n, epsilon.alavanca=0, epsilon.resposta=0, beta, sigma) # apenas um diferente de zero ou os dois epsilons são iguais.
{
  X<-matrix(1,nrow=n,ncol=2)
  X[,2]<-runif(n,min=1, max=10)
  if(epsilon.alavanca == 0 | epsilon.resposta == 0)
  {
    n_outliers <- floor(n * epsilon.alavanca)
    # Ponto de alavanca C2
    if (n_outliers > 0) {
      idx_outliers <- sample(1:n, size = n_outliers, replace = FALSE)
      X[idx_outliers,2]<- runif(n_outliers, min = 16, max = 17)
    }
    
    mu<-exp(X%*%beta) # ligacao log
    y <- rGA(n = n, mu = mu, sigma = sigma)
    
    # Outilers na resposta C1
    n_outliers <- floor(n * epsilon.resposta)
    if (n_outliers > 0) {
      idx_outliers <- sample(1:n, size = n_outliers, replace = FALSE)
      mu_outlier <- mu[idx_outliers] * 5
      y[idx_outliers] <- rGA(n = n_outliers, mu = mu_outlier, sigma = sigma)
    }
  }else if(epsilon.alavanca == epsilon.resposta & epsilon.resposta != 0){
    n_outliers <- floor(n * epsilon.alavanca)
    # Ponto de alavanca e outlier C3
    if (n_outliers > 0) {
      idx_outliers <- sample(1:n, size = n_outliers, replace = FALSE)
      X[idx_outliers,2]<- runif(n_outliers, min = 16, max = 17)
      mu<-exp(X%*%beta) # ligacao log
      y <- rGA(n = n, mu = mu, sigma = sigma)
      mu_outlier <- mu[idx_outliers] * 5
      y[idx_outliers] <- rGA(n = n_outliers, mu = mu_outlier, sigma = sigma)
    }
    
  }
  
  
  return(dados=data.frame(resp=y,x=X[,1],x2=X[,2]))
}