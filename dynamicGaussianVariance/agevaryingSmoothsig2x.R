#################################################################
##### Author: Viviana G R Lobo
#### section 3.3 Age-Varying smoothness - Gaussian
#### with discounting variance
### data: England + Wales, male and female, 2010-2012
#################################################################


setwd("~/Dropbox/Semestre 2024.1 UFRJ/Projeto Multivariado (viviana)")

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


source("dynamicGaussiansig2x.R")

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


### Usando ano de 2010 para construir cenarios
ytime=2010

### Dados do cenario 1
df.mx<- mx %>%
  mutate(log.mx=log(mx))

y.2010m<- df.mx %>%
  filter(year=='2010', sex=='1') %>%
  dplyr:: select(age,log.mx)

y.2010m<-y.2010m$log.mx
mx.m2010<- mx %>%
  filter(year=='2010', sex=='1') %>%
  dplyr::select(age,mx)

p=2 ; w=length(y.2010m)
Gt <- matrix(c(1, 0, 1, 1), p, p)
Ft <- matrix(c(1, 0), p, 1)
delta <- rep(c(0.99, 0.8, 0.85, 0.99), c(5, 30, 50, 19))
m0 <- rep(0, p)
C0 <- diag(100, p)

fit1= dlm.sig2x(y.2010m, delta = 0.75,  prior.sig2 = list(n0=40, d0=40*var(y.2010m[1:20])))
fit4= dlm.sig2x(y.2010m, delta = 0.9999,  prior.sig2 = list(n0=40, d0=40*var(y.2010m[1:20])))
fit5= dlm.sig2x(y.2010m, delta = delta, ages = 1:w, prior = list(m0 = m0, C0 = C0),
                                    prior.sig2 = list(n0=40, d0=40*var(y.2010m[1:20])))


d1<- fitted(fit1)
d4<- fitted(fit4)
d5<- fitted(fit5)

qxall<- bind_rows(d1,d4,d5, .id="id")
qxall



mx.m2010. <- mx.m2010 %>%
  filter(age<=60)

qxall. = qxall %>% 
  filter(age <=60)



p1<-ggplot(NULL,aes(x = 0:60)) + 
  geom_point(data =mx.m2010., aes(x = age, y = mx), col = "steelblue") +
  geom_line(data=qxall., aes(x=age, y =  -log(1-qx.fitted), color = id),col = "steelblue") +
  theme_classic(base_size = 20) + 
  geom_ribbon(data = qxall., aes(x = age, ymin =  -log(1-qx.lower), ymax =  -log(1-qx.upper), fill="id"), alpha = 0.25) + 
  scale_y_continuous(expression(m[x]), limits = c(0.00003, 0.01),
                     trans = 'log10', labels = scales::comma) +
  scale_x_continuous("Age", breaks = seq(0, 60, by = 10)) +  
  scale_fill_manual(values=c("steelblue","steelblue","steelblue","steelblue","steelblue"), labels=c("1", "2", "3", "4","5"), guide="none") + 
  
  theme(legend.position.inside = c(0.05, 0.86),
        strip.background=element_rect(colour="black",
                                      fill="gray87"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title=element_blank(),
        legend.key.height = unit(.6, "cm") ,
        legend.text=element_text(color="black", size=16),
        axis.title = element_text(color = "black", size = 16),
        axis.text = element_text(color="black",size=16))+
  facet_wrap(~id, ncol=5, labeller = labeller(id = 
                                                c("1" = "0.75",
                                                  
                                                  "2" = "0.999",
                                                  "3" = "different per age"
                                                )))

p1

pdf("FigAgeVarying1sig2x.pdf", width=12, height=5)
p1
graphics.off()
