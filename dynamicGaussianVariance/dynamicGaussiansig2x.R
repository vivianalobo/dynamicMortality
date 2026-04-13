#######################################################################
#### author: Viviana G R Lobo
#### univariate model 
#### mortality dlm function w/ variance for each age x - sig2x
#### t-student
#### seguindo capitulo 10.8 do west and harisson
#### stochastic changes in variance 
#######################################################################

library(dplyr)
library(ggplot2)
library(BayesMortalityPlus)
library(patchwork)
#library(tidyverse)
require(stringr)
library(prodlim)
library(purrr)
library(tidyr)
library(scales)
require(gridExtra)


## filtering
ff <- function(y, Ft, Gt, m0, C0, delta, delta.sig2, n0, d0) {
  
  N <- length(y)
  p <- length(m0)
  
  resultado.m <- matrix(NA_real_, N, p)
  resultado.C <- array(NA_real_, c(N, p, p))
  resultado.W <- array(NA_real_, c(N, p, p))
  resultado.a <- matrix(NA_real_, N, p)
  resultado.R <- array(NA_real_, c(N, p, p))
  resultado.f <- numeric(N)
  resultado.Q <- numeric(N)
  resultado.S <- numeric(N)
  resultado.n <- numeric(N)
  resultado.d <- numeric(N)
  
  if (length(delta) == 1) delta <- rep(delta, N)
  if (length(delta.sig2) == 1) delta.sig2 <- rep(delta.sig2, N)
  
  ## inicialização
  S0 <- 0.01
  
  Wt <- C0 * (1 - delta[1]) / delta[1]
  at <- Gt %*% m0
  Rt <- Gt %*% C0 %*% t(Gt) + Wt
  
  ft <- as.numeric(Ft %*% at)
  Qt <- as.numeric(Ft %*% Rt %*% t(Ft) + S0)
  et <- y[1] - ft
  
  #nt tamanho efetivo da amostra (força da informação)
  #dt escala, ligada à soma de quadrados dos erros
  
  nt <- delta.sig2[1] * n0 + 1
  dt <- delta.sig2[1] * d0 + S0 * et^2 / Qt
  St <- dt / nt
  
  At <- Rt %*% t(Ft) / Qt
  mt <- at + At * et
  Ct <- (St / S0) * (Rt - At %*% Ft %*% Rt)
  
  resultado.m[1, ] <- mt
  resultado.C[1, , ] <- Ct
  resultado.W[1, , ] <- Wt
  resultado.a[1, ] <- at
  resultado.R[1, , ] <- Rt
  resultado.f[1] <- ft
  resultado.Q[1] <- Qt
  resultado.S[1] <- St
  resultado.n[1] <- nt
  resultado.d[1] <- dt
  
  ## filtro recursivo
  for (j in 2:N) {
    
    Wt <- Ct * (1 - delta[j]) / delta[j]
    at <- Gt %*% mt
    Rt <- Gt %*% Ct %*% t(Gt) + Wt
    
    ft <- as.numeric(Ft %*% at)
    Qt <- as.numeric(Ft %*% Rt %*% t(Ft) + St)
    et <- y[j] - ft
    
    nt <- delta.sig2[j] * nt + 1
    dt <- delta.sig2[j] * dt + St * et^2 / Qt
    St.new <- dt / nt
    
    At <- Rt %*% t(Ft) / Qt
    mt <- at + At * et
    Ct <- (St.new / St) * (Rt - At %*% Ft %*% Rt)
    
    
    resultado.m[j, ] <- mt
    resultado.C[j, , ] <- Ct
    resultado.W[j, , ] <- Wt
    resultado.a[j, ] <- at
    resultado.R[j, , ] <- Rt
    resultado.f[j] <- ft
    resultado.Q[j] <- Qt
    resultado.S[j] <- St.new
    resultado.n[j] <- nt
    resultado.d[j] <- dt
    
    St <- St.new
  }
  
  return(list(
    m = resultado.m,
    C = resultado.C,
    a = resultado.a,
    R = resultado.R,
    W = resultado.W,
    f = resultado.f,
    Q = resultado.Q,
    S = resultado.S,
    n = resultado.n,
    d = resultado.d
  ))
}


bs <- function(m, C, a, R, Gt) {
  
  t <- nrow(m)
  p <- ncol(m)
  
  as <- matrix(NA, t, p)
  Rs <- array(NA, c(t, p, p))
  
  as[t,] <- m[t,]
  Rs[t,,] <- C[t,,]
  
  for (i in (t-1):1) {
    Bt <- C[i,,] %*% t(Gt) %*% solve(R[i+1,,])
    as[i,] <- m[i,] + Bt %*% (as[i+1,] - a[i+1,])
    Rs[i,,] <- C[i,,] + Bt %*% (Rs[i+1,,] - R[i+1,,]) %*% t(Bt)
  }
  
  return(list(as = as, Rs = Rs))
}

dlm.sig2x <- function(y,
                      Ft = matrix(c(1,0), 1),
                      Gt = matrix(c(1,0,1,1), 2),
                      delta = 0.85,
                      delta.sig2 = c(rep(0.90, 10), rep(0.99, length(y)-10)),
                      prior = list(m0 = c(0,0), C0 = diag(c(100,100))),
                      prior.sig2 = list(n0 =20, d0 =20*var(y[1:20])),
                      M = 5000,
                      ages = 0:(length(y)-1)) {
  
  ## Filtering
  filter <- ff(
    y = y,
    Ft = Ft,
    Gt = Gt,
    m0 = prior$m0,
    C0 = prior$C0,
    delta = delta,
    delta.sig2 = delta.sig2,
    n0 = prior.sig2$n0,
    d0 = prior.sig2$d0
  )
  
  ## Smoothing
  smooth <- bs(
    m = filter$m,
    C = filter$C,
    a = filter$a,
    R = filter$R,
    Gt = Gt
  )
  
  ## Sampling
  t <- length(y)
  p <- length(prior$m0)
  fit <- list()
  
  
  theta <- array(NA, c(M, t, p))
  mu <- matrix(NA, M, t)
  sig2 <- matrix(NA, M, t)
  
  for (i in 1:t) {
    ## sigma^2_t | D_t
    sig2[,i] <- 1 / rgamma(M,shape = filter$n[i]/2, rate  = filter$d[i]/2 )
    ## theta_t | sigma^2_t
    
    ## Rt vem do modelo de estado
    ##  dt/nt e a escala media da variancia E(sig2t) = dt/nt === E(sig2) = beta/alpha caso variancia constante
    ## variancia media em t
    ## dt/nt Rt variancia media do erro depois de integrar o sig2 fora
    theta[, i, ] <- mvtnorm::rmvt(M,
                                  delta = smooth$as[i, ],
                                  sigma = smooth$Rs[i,,] * (filter$d[i] / filter$n[i]),
                                  df = filter$n[i], type = "shifted")
    
    mu[,i] <- theta[,i,] %*% t(Ft)
    
  }
  
  
  fit$mu = mu
  fit$theta = theta
  fit$sig2 = sig2
  # fit$Wt = filter$Wt
  fit$param = list(mt = filter$m, Ct = filter$C,
                   ft = filter$f, Qt = filter$Q,
                   as = smooth$as, Rs = filter$Rs,nt = filter$n,
                   dt = filter$d,
                   St = filter$S
                   )
  fit$info = list(y = y,
                  ages = ages,
                  Ft = Ft,
                  Gt = Gt,
                  delta = delta,
                  delta.sig2 = delta.sig2,
                  prior = prior,
                  prior.sig2 = prior.sig2)
  
  
  
  
  return(structure(fit, class = "DLM"))
}


predict.DLM.sig2x <- function(object, h, prob = 0.95){
  
  fit <- object
  
  ## dimensões
  N <- length(fit$info$y)
  p = length(fit$info$prior$m0)
  y = fit$info$y
  Ft <- fit$info$Ft
  Gt <- fit$info$Gt
  delta <- fit$info$delta[length(fit$info$delta)]
  delta.sig2 <- fit$info$delta.sig2[length(fit$info$delta.sig2)]
  
  aux = fit$param
  n<- dim(fit$sig2)[[1]]
  
  ## simulações
  sim <- matrix(NA_real_, nrow = n, ncol = h)
  
  
  Wt = aux$Ct[N,,] * (1 - delta) / delta
  at = Gt %*% aux$mt[N,]
  Rt = Gt %*% aux$Ct[N,,] %*% t(Gt) + Wt
  ft = as.numeric(Ft %*% at)
  Qt = as.numeric(Ft %*% Rt %*% t(Ft) + 1)
  
  nt <-   delta.sig2 * aux$nt[N] 
  dt <-  delta.sig2 * aux$dt[N] 
  St <- dt / nt

  ## preditiva t-Student
  sim[,1] <- ft + sqrt(Qt*St) * rt(n, df = nt)
  

  if(h>1) for(k in 2:h){
    #Wt = Ct
    
    at = Gt %*% at
    Rt = Gt %*% Rt %*% t(Gt) + Wt
    ft = as.numeric(Ft %*% at)
    Qt = as.numeric(Ft %*% Rt %*% t(Ft) + 1)
    
    nt <- delta.sig2 * nt 
    dt <-  delta.sig2 *dt 
    St <- dt / nt
    
    ## preditiva t-Student
    sim[,k] <- ft + sqrt(Qt*St) * rt(n, df = nt)
  }
  
  
  qx_sim = exp(sim)
  qx_fitted = apply(qx_sim, 2, median, na.rm = T)
  qx_lim = apply(qx_sim, 2, quantile, probs = c((1-prob)/2, (1+prob)/2), na.rm = T)
  
  qx_fitted = data.frame(Ages = (fit$info$ages[N]+1):(fit$info$ages[N]+h), qx_fitted = qx_fitted)
  ret = data.frame(age = qx_fitted$Ages, qx.fitted = 1 - exp(-qx_fitted$qx_fitted),
                   qx.lower = 1 - exp(-qx_lim[1,]), qx.upper = 1 - exp(-qx_lim[2,]))
  ret[ret[,2:4] < 0,2:4] = 0
  ret[ret[,2:4] > 1,2:4] = 1
  
  return(ret)
}


