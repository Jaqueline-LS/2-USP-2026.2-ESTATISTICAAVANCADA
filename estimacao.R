gera.amostras <- function(n, epsilon.alavanca = 0, epsilon.resposta = 0,  beta = c(0.3, 0.1), sigma_fixo = 0.5) {
  
  amostras.estimacao <- foreach(
    i = 1:M, 
    .packages = "gamlss",
    .export = c("geraRegGamma", "M", "nome.amostra") 
  ) %dopar% {
    source("funcoes.R")
    
    dados <- geraRegGamma(n, epsilon.alavanca, epsilon.resposta, beta=beta,sigma=sigma_fixo)
    
    ajuste <- gamlss(resp ~ x2, family = GA, data = dados, trace = FALSE,
                     sigma.start = sigma_fixo, sigma.fix = TRUE)
    
    theta.emv <- ajuste$mu.coefficients
    amostras.estimacao<-list(theta.emv = theta.emv)
  }
  
  saveRDS(amostras.estimacao, file = paste0("Amostras/", nome.amostra, "_", n, ".RDS"))
}
