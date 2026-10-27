suppressMessages(source("estimacao.R"))
suppressMessages(library("doParallel"))

M<-3
tam.n<-c(50,200,1000)
beta<-c(0.3,0.1)
phi<- 4  # phi = 1/sigma2
sigma<-1/sqrt(phi)
alfas<-c(0,0.1,0.25,0.5,0.75,1)

## C0
estima.alpha<-function(alpha, nome.amostra,ea = 0, er = 0)
{
  num_cores <- detectCores() - 3
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  sapply(tam.n, gera.amostras, alpha=alpha, epsilon.alavanca = ea, epsilon.resposta = er,  beta = c(0.3, 0.1), sigma_fixo = 0.5, nome.amostra)
  stopCluster(cl)
}
lapply(alfas, estima.alpha, nome.amostra="C0")


## C1
# 5%

lapply(alfas, estima.alpha,nome.amostra="C1.5",ea = 0, er = 0.05)

# 10%

lapply(alfas, estima.alpha, nome.amostra="C1.10", ea = 0, er = 0.10)




