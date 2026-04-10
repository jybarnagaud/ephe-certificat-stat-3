#------------------------------------------------------------------------------#
# Certificat en analyse de données - niveau perfectionnement
# séquence: inférence bayésienne
# auteur : Jean-Yves Barnagaud (jean-yves.barnagaud@ephe.psl.eu)
# révisé : 04/2026

# réplication des analyses de la séquence de cours
#------------------------------------------------------------------------------#


## installer les packages nécessaires ------------------------------------------

# install.packages("V8", dependencies = TRUE)
# install.packages("rstan", dependencies = TRUE)
# install.packages("brms", dependencies = TRUE)

## charger les packages --------------------------------------------------------

library(ggplot2)
library(brms)
library(rjags)
library(bayesplot)
library(tidyverse)

## tutoriels -------------------------------------------------------------------

# BRMS :  https://bookdown.org/ajkurz/DBDA_recoded/dichotomous-predicted-variable.html
# JAGS : https://bookdown.org/kevin_davisross/bayesian-reasoning-and-methods/introduction-to-jags.html

## simulation - survie des ours ------------------------------------------------

n <- 60
y <- 19

# dénominateur du théorème de Bayes

num <- function(theta) {dbinom(y,n,theta) * dbeta(theta,1,1)}
den <- integrate(num,0,1)$value
den

# distribution a priori
prior <- seq(0, 1, 0.01)
post <- data.frame(survie = prior, ratio <- num(prior) / den)

ggplot(post)+
  aes(x = survie,y = ratio)+
  geom_line(linewidth = 1.5)+
  labs(x = "Probabilité de survie",y = "Densité a posteriori") +
  theme_minimal()


## données ---------------------------------------------------------------------

herb <- read.csv2("donnees/Diversite_vegetale_placette.csv", header = T)
herb07 <- subset(herb, Annee == 2007)
fagus07 <- subset(herb07, Nom_Latin == "Fagus sylvatica" &
                    Presence == 1)

# centrage-réduction du prédicteur (facultatif mais conseillé pour les modèles
# complexes ou pour comparer entre différents algorithmes. Si vous voulez vérifier
# la cohérence entre les résultats de STAN et JAGS dans le script ci-dessous, 
# privilégiez la variable centrée réduite. Si vous voulez des résultats interprétables
# sur l'échelle naturelle, la variable brute suffit). 

fagus07$sdivspe <- as.vector(scale(fagus07$Div_spe_veg))

## exploration -----------------------------------------------------------------

ggplot(fagus07) +
  aes(y = Div_spe_veg, x = factor(Consommation)) +
  geom_boxplot(fill = "gray70")+
  theme_classic()+
  labs(y = "Diversité du peuplement",x = "Abroutissement")+
  coord_flip()

## modèle fréquentiste ---------------------------------------------------------

frq.fagus07 <- glm(Consommation ~ sdivspe, family = "binomial", data = fagus07)
summary(frq.fagus07)

## modèle bayésien : choix des priors ------------------------------------------

# choix des priors

prior1 <- rnorm(1000, 0, 10) # densité de proba d'un prior normal à sd = 10
logit.prior <- plogis(prior1)
hist(logit.prior)

prior2 <- rnorm(1000, 0, 100) # densité de proba d'un prior normal à sd = 100
logit.prior2 <- plogis(prior2)
hist(logit.prior2)

prior3 <- rnorm(1000, 0, 2) # densité de proba d'un prior normal à sd = 2
logit.prior3 <- plogis(prior3)
hist(logit.prior3)

par(mfrow = c(1,3))
hist(logit.prior, main = "N(0,10)")
hist(logit.prior2, main = "N(0,100)")
hist(logit.prior3, main = "N(0,2)")

## modèle bayésien avec brms ---------------------------------------------------

brms.fagus07 <- brm(
  Consommation |
    trials(1) ~ 1 + sdivspe,
  family = "binomial",
  data = fagus07,
  prior = c(prior(normal(0, 2), class = Intercept), prior(normal(0, 2), class = b)),
  iter = 2500,
  warmup = 500,
  chains = 3,
  cores = 3,
  seed = 21
)

summary(brms.fagus07)
plot(brms.fagus07)

# autre solution, équivalente (supposément plus efficace, peu sensible ici)

brms.fagus07.2 <- brm(
  Consommation  ~ 1 + Div_spe_veg,
  family = "bernoulli",
  data = fagus07,
  prior = c(prior(normal(0, 2), class = Intercept), prior(normal(0, 2), class = b)),
  iter = 2500,
  warmup = 500,
  chains = 3,
  cores = 3,
  seed = 21
)

summary(brms.fagus07.2)
plot(brms.fagus07.2)

## modèle bayésien avec jags ---------------------------------------------------

# données 

dat.jags <- list(
  consommation = fagus07$Consommation,
  div_spe_veg = fagus07$Div_spe_veg,
  n.parcelles = nrow(fagus07)
)

# paramètres à conserver

param <- c("alpha","beta")

# script du modèle (attention au codage des priors : sous JAGS, la loi normale
# est exprimée en moyenne et précision, la précision étant 1/variance. sous
# Stan, elle est exprimée en moyenne et écart-type)

mod <- "scripts/ephe-niveau3-mardi-cours-jags.R"

# run du modèle

jags.fagus07 <- jags.model(file = mod, 
                    data = dat.jags,n.chains = 3, n.adapt = 500)

update(jags.fagus07,n.iter = 500)

post.jags.fagus07 <- coda.samples(jags.fagus07,
                                 variable.names = param,
                                 n.iter = 2500)

# vérification des chaines rapide (package coda)

plot(post.jags.fagus07)

mcmc_pairs(post.jags.fagus07, pars = c("alpha", "beta"))

# inférence graphique (package bayesplot)

mcmc_areas(
  post.jags.fagus07,            
  pars = c("beta"),
  prob = 0.95) 

mcmc_trace(post.jags.fagus07, pars = "beta") 

mcmc_combo(post.jags.fagus07)

# inférence quantitative

summary(post.jags.fagus07)
