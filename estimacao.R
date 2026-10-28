gera.amostras <- function(n, alpha, epsilon.alavanca = 0, epsilon.resposta = 0.1,  beta.verd = c(0.3, 0.1), sigma_fixo = 0.5, nome.amostra) {
  
  amostras.estimacao <- foreach(
    i = 1:M, 
    .packages = "gamlss",
    .export = c("M", "nome.amostra") 
  ) %dopar% {
    source("funcoes.R")
    
    dados <- geraRegGamma(n, epsilon.alavanca, epsilon.resposta, beta=beta.verd,sigma=sigma_fixo)
    
    # EMV
    ajuste <- gamlss(resp ~ x2, family = GA, data = dados, trace = FALSE,
                     sigma.start = sigma_fixo, sigma.fix = TRUE)
    
    theta.emv <- ajuste$mu.coefficients
    
    # MDPDE
    beta0<-as.matrix(theta.emv)
    phi<-1/sigma_fixo^2
    X<-as.matrix(dados[,-1])
    Y<-dados$resp

    criterio<-criterio<-1E5

    while(criterio>=1E-5)
    {
     
      eta0 <- X %*% beta0
      mu0  <- as.numeric(exp(eta0))
      
      w <- as.numeric(dGA(Y, mu = mu0, sigma = sigma_fixo)^alpha)
      W <- diag(w)
      
      Z_adj <- eta0 + (Y - mu0) / mu0
      
      b.soma <- calcular_b_soma(X, mu0, alpha, phi = phi)
     
      Hessiana <- phi * t(X) %*% W %*% X
      Vetor_Z  <- phi * t(X) %*% W %*% Z_adj
      
      beta <- solve(Hessiana) %*% (Vetor_Z - b.soma)
      
      criterio <- sum((beta - beta0)^2)
      beta0 <- beta

    }
    theta.emdpde<-beta
    
    amostras.estimacao<-list(theta.emv = theta.emv,theta.emdpde=as.numeric(theta.emdpde))
  }
  
  saveRDS(amostras.estimacao, file = paste0("Amostras/", nome.amostra, "_", n,"(",alpha,")", ".RDS"))
}
