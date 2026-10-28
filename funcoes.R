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

calcular_b_soma <- function(X, mu, alpha, phi) {

  termo_num <- alpha * phi^(1 - alpha) * gamma((1 + alpha) * phi - alpha)
  termo_den <- (1 + alpha)^((1 + alpha) * phi - alpha + 2) * (gamma(phi))^(1 + alpha)
  constante <- - (termo_num / termo_den)
  
  pesos_b <- constante * (mu^(-alpha))
  
  b_soma <- t(X) %*% pesos_b
  
  return(b_soma)
}


# Tabelas



recupera.intercepto.emv<-function(n,alpha, nome.amostra)
{
  
  arquivo<-paste0("Amostras/", nome.amostra, "_", n,"(",alpha,")", ".RDS")
  amostras<-readRDS(arquivo)
  interceptos.emv<-sapply(1:M,FUN=function(x){amostras[[x]]$theta.emv[1]})
  return(interceptos.emv)
}

recupera.intercepto.emdpde<-function(n,alpha, nome.amostra)
{
  
  arquivo<-paste0("Amostras/", nome.amostra, "_", n,"(",alpha,")", ".RDS")
  amostras<-readRDS(arquivo)
  interceptos.emdpde<-sapply(1:M,FUN=function(x){amostras[[x]]$theta.emdpde[1]})
  return(interceptos.emdpde)
}

recupera.beta1.emv<-function(n,alpha, nome.amostra)
{
  
  arquivo<-paste0("Amostras/", nome.amostra, "_", n,"(",alpha,")", ".RDS")
  amostras<-readRDS(arquivo)
  betas.emv<-sapply(1:M,FUN=function(x){amostras[[x]]$theta.emv[2]})
  return(betas.emv)
}

recupera.beta1.emdpde<-function(n,alpha, nome.amostra)
{
  arquivo<-paste0("Amostras/", nome.amostra, "_", n,"(",alpha,")", ".RDS")
  amostras<-readRDS(arquivo)
  betas.emdpde<-sapply(1:M,FUN=function(x){amostras[[x]]$theta.emdpde[2]})
  return(betas.emdpde)
}

calcula_desempenho_completo.intercepto <- function(alfas, tam.n, valor_real, nome.amostra) {
  
  # Itera sobre todos os valores de alpha
  resultados_lista <- lapply(alfas, function(alpha) {
    
    # 1. Recupera as estimativas para EMV e MDPDE repassando nome.amostra
    interceptos.emv <- lapply(tam.n, function(n) {
      recupera.intercepto.emv(n = n, alpha = alpha, nome.amostra = nome.amostra)
    })
    
    interceptos.emdpde <- lapply(tam.n, function(n) {
      recupera.intercepto.emdpde(n = n, alpha = alpha, nome.amostra = nome.amostra)
    })
    
    # 2. Calcula as métricas para o EMV
    metricas_emv <- sapply(1:length(tam.n), function(i) {
      est <- interceptos.emv[[i]]
      vies <- mean(est - valor_real)
      eqm  <- mean((est - valor_real)^2)
      var  <- var(est)
      return(c(Vies = vies, EQM = eqm, Variancia = var))
    })
    
    # 3. Calcula as métricas para o MDPDE
    metricas_mdpde <- sapply(1:length(tam.n), function(i) {
      est <- interceptos.emdpde[[i]]
      vies <- mean(est - valor_real)
      eqm  <- mean((est - valor_real)^2)
      var  <- var(est)
      return(c(Vies = vies, EQM = eqm, Variancia = var))
    })
    
    # 4. Formata os resultados do alpha atual em data.frames
    df_emv <- as.data.frame(t(metricas_emv))
    colnames(df_emv) <- c("Vies_EMV", "EQM_EMV", "Var_EMV")
    
    df_mdpde <- as.data.frame(t(metricas_mdpde))
    colnames(df_mdpde) <- c("Vies_MDPDE", "EQM_MDPDE", "Var_MDPDE")
    
    # Junta as colunas do alpha, n, EMV e MDPDE
    res_alpha <- cbind(
      Alpha = alpha,
      n = tam.n,
      df_emv,
      df_mdpde
    )
    
    return(res_alpha)
  })
  
  # Empilha os resultados de todos os alfas numa única tabela
  tabla_final <- do.call(rbind, resultados_lista)
  
  alpha_counts <- table(tabla_final$Alpha)
  alpha_labels <- paste0("$\\alpha = ", names(alpha_counts), "$")
  row_indices <- setNames(as.numeric(alpha_counts), alpha_labels)
  
  # Formatação da tabela em PDF
  tabla_final %>%
    dplyr::select(-Alpha) %>% 
    kbl(
      format = "latex",
      booktabs = TRUE,
      longtable = FALSE,
      digits = 4,
      escape = FALSE,
      caption = paste0("(Bias, MSE, and Variance) between MLE and MDPDE across sample sizes ($n$) and tuning parameters ($\\alpha$) - Intercepto"),
      col.names = c("n", "Bias", "MSE", "Variance", "Bias", "MSE", "Variance"),
      align = c("c", "r", "r", "r", "r", "r", "r")
    ) %>%
    add_header_above(c(" " = 1, "MLE" = 3, "MDPDE" = 3)) %>%
    pack_rows(index = row_indices, escape = FALSE) %>%
    kable_styling(
      latex_options = c("striped", "hold_position", "scale_down"),
      font_size = 7
    )
}


calcula_desempenho_completo.beta1 <- function(alfas, tam.n, valor_real, nome.amostra) {
  
  # Itera sobre todos os valores de alpha
  resultados_lista <- lapply(alfas, function(alpha) {
    
    # 1. Recupera as estimativas para EMV e MDPDE repassando nome.amostra
    beta1.emv <- lapply(tam.n, function(n) {
      recupera.beta1.emv(n = n, alpha = alpha, nome.amostra = nome.amostra)
    })
    
    beta1.emdpde <- lapply(tam.n, function(n) {
      recupera.beta1.emdpde(n = n, alpha = alpha, nome.amostra = nome.amostra)
    })
    
    # 2. Calcula as métricas para o EMV
    metricas_emv <- sapply(1:length(tam.n), function(i) {
      est <- beta1.emv[[i]]
      vies <- mean(est - valor_real)
      eqm  <- mean((est - valor_real)^2)
      var  <- var(est)
      return(c(Vies = vies, EQM = eqm, Variancia = var))
    })
    
    # 3. Calcula as métricas para o MDPDE
    metricas_mdpde <- sapply(1:length(tam.n), function(i) {
      est <- beta1.emdpde[[i]]
      vies <- mean(est - valor_real)
      eqm  <- mean((est - valor_real)^2)
      var  <- var(est)
      return(c(Vies = vies, EQM = eqm, Variancia = var))
    })
    
    # 4. Formata os resultados do alpha atual em data.frames
    df_emv <- as.data.frame(t(metricas_emv))
    colnames(df_emv) <- c("Vies_EMV", "EQM_EMV", "Var_EMV")
    
    df_mdpde <- as.data.frame(t(metricas_mdpde))
    colnames(df_mdpde) <- c("Vies_MDPDE", "EQM_MDPDE", "Var_MDPDE")
    
    # Junta as colunas do alpha, n, EMV e MDPDE
    res_alpha <- cbind(
      Alpha = alpha,
      n = tam.n,
      df_emv,
      df_mdpde
    )
    
    return(res_alpha)
  })
  
  # Empilha os resultados de todos os alfas numa única tabela
  tabla_final <- do.call(rbind, resultados_lista)
  
  alpha_counts <- table(tabla_final$Alpha)
  alpha_labels <- paste0("$\\alpha = ", names(alpha_counts), "$")
  row_indices <- setNames(as.numeric(alpha_counts), alpha_labels)
  
  # Formatação da tabela em PDF
  tabla_final %>%
    dplyr::select(-Alpha) %>% 
    kbl(
      format = "latex",
      booktabs = TRUE,
      longtable = FALSE,
      digits = 4,
      escape = FALSE,
      caption = paste0("(Bias, MSE, and Variance) between MLE and MDPDE across sample sizes ($n$) and tuning parameters ($\\alpha$) - Beta1"),
      col.names = c("n", "Bias", "MSE", "Variance", "Bias", "MSE", "Variance"),
      align = c("c", "r", "r", "r", "r", "r", "r")
    ) %>%
    add_header_above(c(" " = 1, "MLE" = 3, "MDPDE" = 3)) %>%
    pack_rows(index = row_indices, escape = FALSE) %>%
    kable_styling(
      latex_options = c("striped", "hold_position", "scale_down"),
      font_size = 7
    )
}





geraBoxplot<-function(nome.amostra)
{
  jpeg(file=paste0("Resultados/",paste0(nome.amostra,"_Intercepto" ,".png")), width = 1000, 
       height = 500, quality = 100, pointsize = 20)
  par(mfrow=c(2,3), mar=c(3,3,3, 2))
  
  boxplot.intercepto<-function(alpha)
  {
    
    interceptos.emv<-lapply(tam.n, recupera.intercepto.emv, alpha, nome.amostra)
    
    
    interceptos.emdpde<-lapply(tam.n, recupera.intercepto.emdpde, alpha, nome.amostra)
    
    paleta_completa <- c("#EADDF0", "#C2A5CF", "#8B629C", "#E8F5D6", "#C9E69E", "#8CB858")
    boxplot(c(interceptos.emv,interceptos.emdpde),boxwex=0.2,main=paste0(nome.amostra," Intercepto" ," (", alpha,")"), 
            pch=1,at=c(1,1.25,1.5, 1.75, 2, 2.25), xaxt = "n",col=paleta_completa)
    abline(h=beta.verd[1], lty=2, lwd=2)
    
    if(alpha==0)legend("topright",legend = paste0(rep(c("EMV","MDPDE"),each=3)," - n = ",tam.n), 
                       col=paleta_completa, bty="n", fill=paleta_completa, cex=0.5)
  }
  
  lapply(alfas, boxplot.intercepto)
  dev.off()
  
  jpeg(file=paste0("Resultados/",paste0(nome.amostra,"_Beta1" ,".png")), width = 1000, height = 500, quality = 100, pointsize = 20)
  par(mfrow=c(2,3), mar=c(3,3,3, 2))
  boxplot.beta1<-function(alfa)
  {
    recupera.beta1.emv<-function(n)
    {
      
      arquivo<-paste0("Amostras/", nome.amostra, "_", n,"(",alfa,")", ".RDS")
      amostras<-readRDS(arquivo)
      betas.emv<-sapply(1:M,FUN=function(x){amostras[[x]]$theta.emv[2]})
      return(betas.emv)
    }
    betas.emv<-lapply(tam.n, recupera.beta1.emv)
    
    recupera.beta1.emdpde<-function(n)
    {
      
      arquivo<-paste0("Amostras/", nome.amostra, "_", n,"(",alfa,")", ".RDS")
      amostras<-readRDS(arquivo)
      betas.emdpde<-sapply(1:M,FUN=function(x){amostras[[x]]$theta.emdpde[2]})
      return(betas.emdpde)
    }
    betas.emdpde<-lapply(tam.n, recupera.beta1.emdpde)
    
    paleta_completa <- c("#EADDF0", "#C2A5CF", "#8B629C", "#E8F5D6", "#C9E69E", "#8CB858")
    boxplot(c(betas.emv,betas.emdpde),boxwex=0.2,main=paste0(nome.amostra," Beta1" ," (", alfa,")"), pch=1,at=c(1,1.25,1.5, 1.75, 2, 2.25), xaxt = "n",col=paleta_completa)
    abline(h=beta.verd[2], lty=2, lwd=2)
    
    if(alfa==0)legend("topright",legend = paste0(rep(c("EMV","MDPDE"),each=3)," - n = ",tam.n), col=paleta_completa, bty="n", fill=paleta_completa, cex=0.5)
  }
  lapply(alfas, boxplot.beta1)
  dev.off()
  
}

