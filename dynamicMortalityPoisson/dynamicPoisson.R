#################################################################
##### Author: Viviana G R Lobo
##### Modelo Poisson Dinamico + Extrapolacao h passos a frente
##### Linear Bayes 
##### Proposta do artigo: West, Harrison e Migon (1985)
### data: England + Wales, male and female, 2010-2012
#################################################################

setwd("~/Dropbox/Semestre 2025.1 UFRJ/Poisson Dinamico")

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

grad <- seq(0,110,20)
grad[(grad/5)%%2 != 0] <- ""
point <- format_format(big.mark = " ", decimal.mark = ".", scientific = FALSE)


#####==========#####
##### dataset  #####
#####==========#####
dx <- read.csv("deaths_MYB2.csv")
Ex <- read.csv("MYB2_mais_OLDAGE.csv")

dx.all = dx %>%
  dplyr::select(!X) %>% 
  pivot_longer(cols = starts_with("X"), values_to = "dx", names_to = "year") %>%
  mutate(year = substring(year, 2)) %>%
  filter(Age != '105+' & Age != '0') %>%
  mutate(age = as.numeric(Age)) %>%
  dplyr::select(age, sex, year, dx)

Ex.all = Ex %>%
  dplyr::select(!X) %>% 
  pivot_longer(cols = starts_with("population_"), values_to = "Ex", names_to = "year") %>%
  mutate(year = str_remove(year, "population_")) %>%
  filter(age != 105 & age !=0 )

mx = dx.all %>% 
  left_join(Ex.all, by = c("sex", "age", "year")) %>%
  mutate(mx = dx/Ex) %>%
  mutate(sex = as.character(sex))


### masculino, 2010
dt.male <- mx %>%
  filter(sex == 1, year== 2010) %>%
  dplyr::select(age, dx, Ex)

Dx.male<- dt.male$dx
Ex.male<- dt.male$Ex
ages<- dt.male$age 

########============================================================
#######  Filtering + Smoothing
p<-2
mt.s <- matrix(0, p, length(ages))      
Ct.s <- array(0, dim = c(p, p, length(ages)))
at.s<- matrix(NA, nrow = p, ncol = length(ages))         
Rt.s <- array(NA, dim = c(p, p, length(ages)))            
Wt.s <- array(NA, dim = c(p, p, length(ages)))  

# media.post = var.post = ICupper = IClower <- NULL

Gt <- matrix(c(1, 0, 1, 1), p, p)
Ft <- matrix(c(1, 0), p, 1)
delta<- rep(c(0.99, 0.8, 0.85, 0.99), c(5, 30, 50, 19))
m0 <- rep(0, nrow(Gt))
C0 <- diag(100, nrow(Gt))

## Filtro de Kalman
### passo 1 (inicializacao)
Wt <- C0 * ((1 - delta[1]) / delta[1])
### prior eta_t | D_t-1 ; theta_t | D_t-1
at <- Gt %*% m0            
Rt <- Gt %*% C0 %*% t(Gt) + Wt
### E([lambda_t | D_t-1]; Var[lambda_t | D_t-1]; Cov[lambda_t, theta_t, | D_t-1]
ft <- as.numeric(t(Ft) %*% at) +log(Ex.male[1])
st<- Rt%*%Ft
### compatibilizando priori
qt <- as.numeric(t(Ft) %*%st)
alpha<- 1/qt 
beta<-  alpha*exp(-ft)

### posteriori
### parte do eta_t - caso poisson-gamma (ver west, harrison, migon 1985)
alpha_post <- alpha + Dx.male[1]
beta_post <- beta + 1#Ex.male[1]

# media.post[1] <- alpha_post/beta_post
# var.post[1] <- alpha_post/(beta_post^2)
# ICupper[1] <- media.post[1] + 2*sqrt(var.post[1])
# IClower[1] <- media.post[1] - 2*sqrt(var.post[1])

# gt = E[g(eta_t) | Dt] and pt = V[g(eta_t) | Dt]  - posterior CP
gt <- digamma(alpha_post) - log(beta_post)
pt <- trigamma(alpha_post)

mt <- at + st %*% ((gt - ft)/qt)
Ct <- Rt - st%*%t(st)*((1-(pt/qt))/qt)

mt.s[,1] <- mt
Ct.s[,,1]<- Ct
Wt.s[,,1]<- Ct * (1 - delta[1]) / delta[1]
at.s[,1]<- at
Rt.s[,,1]<- Rt

### passo 2 (updating)
for (t in 2:length(ages)) {
  
  Wt.s[,,t] =  Ct.s[,,t-1]* ((1 - delta[t]) / delta[t])  
  ### state transition
  at.s[,t] <- Gt %*% mt.s[,t-1]
  Rt.s[,,t] <- Gt %*% Ct.s[,,t-1] %*% t(Gt) + Wt.s[,,t]
  ft <- as.numeric(t(Ft) %*% at.s[,t]) + log(Ex.male[t])
  st<- Rt.s[,,t]%*%Ft
  qt <- as.numeric(t(Ft) %*%st)
  alpha<- 1/qt 
  beta<-  alpha*exp(-ft )
  
  #### posterior eta_t | Dt according to Migon1985
  # where gt = E[g(eta_t) | Dt] and pt = V[g(eta_t) | Dt] 
  alpha_post <- alpha + Dx.male[t]
  beta_post <- beta +  1#Ex.male[t]
  gt <- digamma(alpha_post) - log(beta_post)
  pt <- trigamma(alpha_post)
  
  ### state update theta_x | Dx 
  mt.s[,t] <- at.s[,t] + st %*% ((gt - ft)/qt)
  Ct.s[,,t] <- Rt.s[,,t] - st%*%t(st)* ((1-(pt/qt))/qt)
  
  # media.post[t] <- alpha_post/beta_post
  # var.post[t] <- alpha_post/(beta_post^2)
  # ICupper[t] <- media.post[t] + 2*sqrt(var.post[t])
  # IClower[t] <- media.post[t] - 2*sqrt(var.post[t])
  
}

#### backward sampling - linear Bayes
N <- length(ages)
theta.smooth <- matrix(NA, nrow = p, ncol = N)
C.smooth <- array(NA, dim = c(p, p, N))
ICupper.s<- IClower.s<- NULL

theta.smooth[, N] <- mt.s[, N]
C.smooth[, , N] <- Ct.s[, , N]

ICupper.s[N] <-theta.smooth[1,N] + 2*sqrt(C.smooth[1,1,N])
IClower.s[N] <- theta.smooth[1,N] - 2*sqrt(C.smooth[1,1,N])

for (t in (N - 1):1) {
  Bt <- Ct.s[, , t] %*% t(Gt) %*% chol2inv(chol(Rt.s[, , t + 1]))
  theta.smooth[, t] <- mt.s[, t] + Bt %*% (theta.smooth[, t + 1] - at.s[, t + 1])
  C.smooth[, , t] <- Ct.s[, , t] -Bt %*% ( C.smooth[, , t + 1] - Rt.s[, , t + 1]) %*% t(Bt)
  
  ICupper.s[t] <-theta.smooth[1,t] + 2*sqrt(C.smooth[1,1,t])
  IClower.s[t] <- theta.smooth[1,t] - 2*sqrt(C.smooth[1,1,t])
}


# plot(log(Dx.male/Ex.male),pch=19, col="gray", ylab = "Nível", xlab = "Idade", main = "", ylim=c(-10,1))
# #lines(log(ICupper), lty=2, col="blue")
# #lines(log(IClower), lty=2, col="blue")
# lines(ages, mt.s[1, ], col = "red", lty = 2)
# lines(theta.smooth[1,],col="darkgreen")
# lines(ICupper.s,col="darkgreen", lty=2)
# lines(IClower.s,col="darkgreen", lty=2)
# legend("topleft", c("Suavizado", "Filtrado"), col = c("darkgreen", "red"),
#        lty = c(1, 2), lwd = 2)

#----------------------------------------------------------------------------
#### Previsao h passos a frente
M <- 1000  # número de amostras
N <- length(ages)
#theta.s <- array(NA, dim = c(2, N, M))
# considerando previsao ate os 120 anos - caso lognormal/forster paper
h <- 16
# Vetores para armazenar forecast
a.f <- matrix(NA, nrow = p, ncol = h)
m.f <- matrix(0, p, h)        # estados: nível e tendência
C.f <- array(0, dim = c(p, p, h))  # matriz de covariância dos estados
W.f <- array(NA, dim = c(p, p, h))
R.f <- array(NA, dim = c(p, p, h))
f.f <- numeric(h)
q.f <- numeric(h)
alpha.f <- numeric(h)
beta.f <- numeric(h)
lambda.forecast <- matrix(NA, nrow = M, ncol = h)
Dx.forecast <- matrix(NA, nrow = M, ncol = h)
#### projetando exposicao para idades avancadas - extapolacao
# Estima a taxa de decaimento da população (exposição) a partir de uma regressão log-linear da própria exposição,
# nas idades altas (init = w-35 até o máximo da série, w).
# Usa essa taxa estimada para projetar a exposição futuramente com:
# Se decay_rate < 1, sua população naquelas idades está decaindo exponencialmente com a idade.
# Se decay_rate ≈ 1, ela está mais ou menos constante.
# Se decay_rate > 1 (raro nas idades avançadas), estaria crescendo.
w<- length(ages)
init<- w-35 ## como escolher o melhor ponto de corte aqui? Nao sei. mas estou pensando em pegar sempre 35 anos pra tras para entender
### melhor o comportamento da exposicao.
age.cut <- init:w
log.Ex <- log(Ex.male[age.cut])
idade <- age.cut
fitE <- lm(log.Ex ~ idade)
decay_rate<- exp(fitE$coefficients[2]) #taxa de decaimento exponencial.
Ex.forecast<- c(Ex.male, rep(NA, h))
idade_proj <- (w+1):(w+h)
Ex.forecast[(w+1):(w+h)] <- Ex.forecast[w] * decay_rate^(1:h)


for (m in 1:M) {
  
  # theta.s[, N, m] <- MASS::mvrnorm(1, mt.s[,N], Ct.s[,,N])
  # 
  # for (t in (N-1):1) {
  #   Bt <- Ct.s[,,t] %*% t(Gt) %*% chol2inv(chol(Rt.s[,,t+1])) #(solve(Rt.s[,,t+1]))
  #   ht <- mt.s[,t] + Bt %*% (theta.s[, t+1, m] - at.s[,t+1])
  #   Ht <- Ct.s[,,t] - Bt %*% Rt.s[,,t+1] %*% t(Bt)
  #   theta.s[, t, m] <- MASS::mvrnorm(1, ht, Ht)
  # }
  
  #### forecast - passo x+k
  if(length(delta) > 1){ delta_pred = delta[N] } else{ delta_pred = delta }
  # inicializando
  W.f[, , 1] <- Ct.s[, , N] * ((1 - delta_pred) / delta_pred)  
  a.f[, 1] <- Gt %*% mt.s[, N]
  R.f[, , 1] <- Gt %*% Ct.s[, , N] %*% t(Gt) + W.f[, , 1]
  f.f[1] <- as.numeric(t(Ft) %*% a.f[, 1])
  s.sf <- R.f[, , 1] %*% Ft
  q.f[1] <- as.numeric(t(Ft) %*% s.sf)
  
  alpha.f[1] <- 1 / q.f[1]
  beta.f[1] <- alpha.f[1] * exp(-f.f[1])
  
  ## lambda_x+k and D_x+k
  lambda.forecast[m, 1] <- rgamma(1, shape = alpha.f[1], rate = beta.f[1])
  Dx.forecast[m, 1] <- rpois(1, lambda = Ex.forecast[N+1] * lambda.forecast[m, 1])
  
  ## posterior
  alpha_postf <- alpha.f[1] + Dx.forecast[m, 1]
  beta_postf <- beta.f[1] + Ex.forecast[N+1]
  gt.f <- digamma(alpha_postf) - log(beta_postf)
  pt.f <- trigamma(alpha_postf)
  
  ### state update theta_x+k | Dx+k~ (mx, Cx)
  m.f[,1] <- a.f[,1] + s.sf %*% ((gt.f - f.f[1])/q.f[1])
  C.f[,,1] <- R.f[, , 1] - s.sf%*%t(s.sf)* ((1-(pt.f/q.f[1]))/q.f[1])
  
  if(h > 1) for(k in 2:h){
    ### theta_t+k | D_t-1(k) ~(at(k), Rt(k))
    W.f[, , k] <- C.f[, , k-1] * ((1 - delta_pred) / delta_pred) 
    a.f[, k] <- Gt %*% a.f[, k-1]
    R.f[, , k] <- Gt %*% R.f[, , k-1] %*% t(Gt) + W.f[, , k]
    
    ### eta_t+k | D_t-1(k) ~CP(alphat(k), betat(k))
    f.f[k] <- as.numeric(t(Ft) %*% a.f[, k])
    s.sf <- R.f[, , k] %*% Ft
    q.f[k] <- as.numeric(t(Ft) %*% s.sf)
    
    ### p(Y)_t+k | D_t-1(K) parameters alphat(k) betat(k) - CP conjugate prior
    alpha.f[k] <- 1 / q.f[k]
    beta.f[k] <- alpha.f[k] * exp(-f.f[k])
    
    lambda.forecast[m, k] <- rgamma(1, shape = alpha.f[k], rate = beta.f[k])
    Dx.forecast[m, k] <- rpois(1, lambda = Ex.forecast[N+k] * lambda.forecast[m, k])
    
    alpha_postf <- alpha.f[k] + Dx.forecast[m,k]
    beta_postf <- beta.f[k] +  Ex.forecast[N+k]
    
    gt.f <- digamma(alpha_postf) - log(beta_postf)
    pt.f <- trigamma(alpha_postf)
    
    ### state update theta_x | Dx ~ (mx, Cx)
    m.f[,k] <- a.f[,k] + s.sf %*% ((gt.f - f.f[k])/q.f[k])
    C.f[,,k] <- R.f[, , k] - s.sf%*%t(s.sf)* ((1-(pt.f/q.f[k]))/q.f[k])
  }
}

### --- plot
# Taxa observada (óbitos / exposição)
df.obs <- data.frame( age = ages, mx.obs = Dx.male / Ex.male, rate.obs = log(Dx.male / Ex.male))

# rate.smooth <- exp(theta.s[1, , ])  # exp(mu.s)
# df.fitted<- data.frame(age=ages, rate.m = apply(rate.smooth, 1, mean),
#                        rate.l=apply(rate.smooth, 1, quantile, probs = 0.025),
#                        rate.u = apply(rate.smooth, 1, quantile, probs = 0.975), model = "poisson")


#rate.smooth <- theta.smooth[1,]
df.fitted<- data.frame(age=ages, rate.m = exp(theta.smooth[1,]),
                       rate.l= exp(IClower.s),
                       rate.u = exp(ICupper.s), model = "poisson")

# forecast
h <- ncol(Dx.forecast)  
M <- nrow(Dx.forecast)  
mx.forecast <- matrix(NA, nrow = M, ncol = h)
for (k in 1:h) {mx.forecast[, k] <- Dx.forecast[, k] / Ex.forecast[w + k]}
mean.f <- apply(mx.forecast, 2, mean)
q025.f <- apply(mx.forecast, 2, quantile, probs = 0.025)
q975.f <- apply(mx.forecast, 2, quantile, probs = 0.975)
df.forecast <- data.frame(age = max(ages) + 1:h, rate.m = mean.f, rate.l = q025.f, rate.u = q975.f,model = "poisson")

# plot joint : fit + forecast - na escala mx
df.poisson <- bind_rows( df.fitted, df.forecast)
ggplot() + geom_point(data = df.obs, aes(x = age, y = mx.obs), color = "gray") +
  geom_line(data = df.poisson, aes(x = age, y = (rate.m)), color = "steelblue", linewidth = 0.8) +
  geom_ribbon(data = df.poisson,aes(x = age, ymin = (rate.l), ymax = (rate.u)), fill = "steelblue", alpha = 0.3) +
  scale_y_continuous(expression(m[x]), trans = 'log10', labels = scales::comma) +
  scale_x_continuous("Age", breaks = seq(0, 120, by = 20)) +
  theme_bw() + theme(legend.position = "none")


## na escala do qx
df.fitted2<- data.frame(age=ages, rate.m =  1- exp(-exp(theta.smooth[1,])),
                       rate.l= 1- exp(- exp(IClower.s)),
                       rate.u = 1- exp(- exp(ICupper.s)), model = "poisson")
df.forecast2 <- data.frame(age = max(ages) + 1:h, rate.m = 1- exp(-mean.f), rate.l =  1- exp(-q025.f), 
                           rate.u = 1- exp(-q975.f),model = "poisson")
df.poisson2 <- bind_rows( df.fitted2, df.forecast2)
ggplot() + geom_point(data = df.obs, aes(x = age, y = mx.obs), color = "gray") +
  geom_line(data = df.poisson2, aes(x = age, y = rate.m), color = "steelblue", linewidth = 0.8) +
  geom_ribbon(data = df.poisson2,aes(x = age, ymin = (rate.l), ymax = (rate.u)), fill = "steelblue", alpha = 0.3) +
  scale_y_continuous(expression(q[x]), trans = 'log10', labels = scales::comma) +
  scale_x_continuous("Age", breaks = seq(0, 120, by = 20)) +
  theme_bw() + theme(legend.position = "none")


### na escala do log
df.fitted3<- data.frame(age=ages, rate.m = theta.smooth[1,],
                        rate.l= IClower.s,
                        rate.u = ICupper.s, model = "poisson")
df.forecast3 <- data.frame(age = max(ages) + 1:h, rate.m = log(mean.f), rate.l =  log(q025.f), 
                           rate.u = log(q975.f),model = "poisson")
df.poisson3 <- bind_rows( df.fitted3, df.forecast3)
ggplot() + geom_point(data = df.obs, aes(x = age, y = log(mx.obs)), color = "gray") +
  geom_line(data = df.poisson3, aes(x = age, y = rate.m), color = "steelblue", linewidth = 0.8) +
  geom_ribbon(data = df.poisson3,aes(x = age, ymin = (rate.l), ymax = (rate.u)), fill = "steelblue", alpha = 0.3) +
  scale_y_continuous(expression(log(m[x])),  labels = scales::comma) +
  scale_x_continuous("Age", breaks = seq(0, 120, by = 20)) +
  theme_bw() + theme(legend.position = "none")
####=====================================================
### intervalo preditivo - poisson

### --- gerar intervalo preditivo sobre idades observadas ---
# lambda suavizado nas idades observadas
lambda.smooth <- exp(theta.smooth[1, ])  # passa para a escala exponencial
# lambda forecast já tem dimensão (M x h)
# combinando a suavizacao + forecast
lambda.s <- cbind(matrix(rep(lambda.smooth, each = M), nrow = M), lambda.forecast)
dim(lambda.s)  # (M, N + h)
# --- gerar D simulados
Dx.pred <- matrix(NA, nrow = M, ncol = w + h)
for (i in 1:(w + h)) {
  Dx.pred[, i] <- rpois(M, lambda = Ex.forecast[i] * lambda.s[,i ])
}
# Dx.pred
# plot(apply(Dx.pred,2,mean), t='l')

# --- calcular taxas m[x]
mx.pred <- sweep(Dx.pred, 2, Ex.forecast, FUN = "/")

# --- organizar dataframe
df.pred.poisson <- data.frame(
  age = c(ages, idade_proj),
  rate.m = apply(mx.pred, 2, mean),
  rate.l = apply(mx.pred, 2, quantile, probs = 0.025),
  rate.u = apply(mx.pred, 2, quantile, probs = 0.975),
  model = "poisson pred"
)

df.pred.poisson2 <- data.frame(
  age = c(ages, idade_proj),
  rate.m = apply(1- exp(-mx.pred), 2, mean), #apply(mx.pred, 2, mean),
  rate.l = apply(1- exp(-mx.pred), 2, quantile, probs = 0.025),
  rate.u = apply(1- exp(-mx.pred), 2, quantile, probs = 0.975),
  model = "poisson pred"
)

ggplot() + 
  geom_point(data = df.obs, aes(x = age, y = mx.obs), color = "gray", size = 1.2) +
  geom_ribbon(data = df.pred.poisson,  aes(x = age, ymin = rate.l, ymax = rate.u, fill = "preditivo"), 
              alpha = 0.3) +
  geom_line(data = df.pred.poisson, aes(x = age, y = rate.m), color = "steelblue", linewidth = 0.8) +
  geom_ribbon(data = df.poisson,  aes(x = age, ymin = rate.l, ymax = rate.u, fill = "ajuste"),  
              alpha = 0.3) +
  scale_fill_manual(name = "", values = c("preditivo" = "darkorange", "ajuste" = "steelblue")) +
  scale_y_continuous(expression(m[x]),  trans = 'log10', labels = scales::comma) +
  scale_x_continuous("Age", breaks = seq(0, 120, by = 20)) +
  theme_classic(base_size = 13) +
  theme( legend.position = c(0.8, 0.20),
         strip.background = element_rect(colour = "black", fill = "gray87"),
         panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
         legend.title = element_blank(),
         legend.key.height = unit(0.6, "cm"),
         legend.text = element_text(color = "black", size = 16),
         axis.title = element_text(color = "black", size = 14),
         axis.text = element_text(color = "black", size = 14))

#O modelo ajustado (incerteza do estado latente, via distribuição do log(m_x) suavizado — linhas e faixas azuis);
#A incerteza preditiva total, considerando a variabilidade estocástica dos óbitos D_x via Poisson (faixa laranja).

# 
# #==============================================================
# ###################### Comparando com o modelo log-normal dlm
# #### dlm - lognormal BayesMortalityPlus
# source("predict_dlm.R")
# log.mx <- log(Dx.male / Ex.male)
# fit.dlm <- BayesMortalityPlus::dlm(log.mx, delta = delta,  prior = list(m0 = m0, C0 = C0), ages = 1:w)
# 
# qx.m <- fitted(fit.dlm)
# qx.pred <- predict.DLM(fit.dlm, h = h)
# qx.all <- rbind(qx.m, qx.pred)
# ages.all <- c(ages, max(ages) + 1:h)
# 
# df.dlm <- data.frame(
#   age = ages.all,
#   rate.m = -log(1 - (qx.all$qx.fitted)),
#   rate.l = -log(1 - (qx.all$qx.lower)),
#   rate.u = -log(1 - (qx.all$qx.upper)),
#   model = "log-normal"
# )
# 
# df.plot <- bind_rows(df.poisson,df.pred.poisson, df.dlm)
# df.plot$model <- factor(df.plot$model, levels = c("poisson", "poisson pred", "log-normal"))
# 
# #### acho que o dlm nao esta fornecendo o intervalo de credibilidade preditivo na parte que os dados foram ajustados
# ggplot() +
#   geom_point(data = df.obs, aes(x = age, y = mx.obs), color = "gray", size = 1) +
#   geom_line(data = df.plot, aes(x = age, y = rate.m, color = model), linewidth = 0.8) +
#   geom_ribbon(data = df.plot, aes(x = age, ymin = rate.l, ymax = rate.u, fill = model),   alpha = 0.25) +
#   scale_color_manual(values = c("steelblue", "darkorange", "darkgreen")) +
#   scale_fill_manual(values = c("steelblue", "darkorange", "darkgreen")) +
#   scale_y_continuous(expression(m[x]),                    trans = 'log10',labels = scales::comma) +
#   scale_x_continuous("Age", breaks = seq(0, 120, by = 20)) +
#   theme_classic(base_size = 13) +
#   theme(
#     legend.position = c(0.8, 0.20),
#     strip.background = element_rect(colour = "black", fill = "gray87"),
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
#     legend.title = element_blank(),
#     legend.key.height = unit(0.6, "cm"),
#     legend.text = element_text(color = "black", size = 16),
#     axis.title = element_text(color = "black", size = 14),
#     axis.text = element_text(color = "black", size = 14))
# 
# 

