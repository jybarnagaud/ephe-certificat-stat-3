#------------------------------------------------------------------------------#
# Certificat en analyse de données - niveau perfectionnement
# séquence: n-mixture models / inférence bayésienne
# auteur : Jean-Yves Barnagaud (jean-yves.barnagaud@ephe.psl.eu)
# révisé : 04/2026

# objectif de la séquence : ajuster un n mixture model en cadre bayésien
#------------------------------------------------------------------------------#

rm(list = ls())

library(ggplot2)
library(rjags) # interfaçage de JAGS avec R
library(coda) # utilitaires pour lire des chaînes MCMC
library(bayesplot) # d'autres utilitaires (pour les amateurs de ggplot)
library(MCMCvis) # encore d'autres utilitaires

## modèle écologique sur le Pinson----------------------------------------------

# Préparation des données

# jeu de données: abondance de pinsons dans les forêts de Nouvelle Zélande: 3 réplicats par plot, variables année, jour de l'année, altitude, proportion de forêt native

pinson <- read.csv2("donnees/nouvelle_zelande_pinson.csv")

# exploration des variables explicatives : altitude et proportion de forêt native

ggplot(pinson) +
  aes(x = altitude, y = PROPNATFOR) +
  geom_point()

# préparer la variable de réponse : comptage maxi sur 3*5 minutes

obs <- as.matrix(pinson[, 3:5])
obs.naif <- apply(pinson[, 3:5], 1, "max")

# préparer les variables explicatives

PROPNATFOR <- as.vector(scale(pinson$PROPNATFOR)) # on centre-réduit pour limiter les risques de mauvaises convergences (également applicable aux modèles fréquentistes)
altitude <- as.vector(scale(pinson$altitude))
year <- as.vector(as.numeric(factor(pinson$year))) # on a passé les années en facteur. 1=2005, 2=2006
julian_date <- as.vector(scale(pinson$julian_date))

# les indices dont on aura besoin

npoints <- nrow(obs)
nrepl <- ncol(obs)
npoints.naif <- length(obs.naif)
nyears <- length(unique(year))

# liste des données

nzdata.naif <- list(
  N = obs.naif,
  PROPNATFOR = PROPNATFOR,
  altitude = altitude,
  year = year,
  npoints = npoints,
  nyears = nyears
)

# choisir les paramètres à suivre

monitor.naif <- c("alpha", "beta", "delta", "esp", "r.naif", "E", "E.rep")

# lancer le modèle (run d'adaptation)

jnz.naif <- jags.model(
  file = "donnees/nouvelle_zelande_modele_ecol.R",
  data = nzdata.naif,
  n.chains = 3,
  n.adapt = 1000
)

# burn in

update(jnz.naif, n.iter = 5000)

# codas (= chaines de Markov)

jnz.naif.coda <- coda.samples(jnz.naif,
                             variable.names = monitor.naif,
                             n.iter = 5000,
                             thin = 50)

# explorer les chaines

color_scheme_set("viridis")
mcmc_trace(jnz.naif.coda, pars = c("alpha[1]", "alpha[2]", "beta", "delta"))

# explorer les paramètres de plusieurs manières

mcmc_areas(jnz.naif.coda, pars = c("alpha[1]", "alpha[2]", "beta", "delta"))
mcmc_intervals(jnz.naif.coda, pars = c("alpha[1]", "alpha[2]", "beta", "delta"))
mcmc_dens_overlay(jnz.naif.coda, pars = c("alpha[1]", "alpha[2]"))

# explorer les corrélations entre paramètres (! on peut s'inquiéter de la relation
# entre altitude et forêt native!)

mcmc_pairs(
  jnz.naif.coda,
  pars = c("alpha[1]", "alpha[2]", "beta", "delta"),
  off_diag_args = list(size = 1.5)
)

# explorer les résidus

resid.nz <- MCMCchains(jnz.naif.coda, params = c("r.naif"))
med.resid.nz <- apply(resid.nz, 2, median)
esp.nz <- MCMCchains(jnz.naif.coda, params = c("esp"))
med.esp.nz <- apply(esp.nz, 2, median)

# calcul des résidus de pearson

r.pearson <- med.resid.nz / sqrt(med.esp.nz)
plot(med.esp.nz, r.pearson) # on retrouve le diagnostic classique de résidus

# ajustement

E.rep.nz <- MCMCchains(jnz.naif.coda, params = c("E.rep"))
E.nz <- MCMCchains(jnz.naif.coda, params = c("E"))
E.med <- apply(E.nz, 2, median)
E.rep.med <- apply(E.rep.nz, 2, median)
plot(E.med, E.rep.med) # cet ajustement n'est vraiment pas bon.

## modèle n-mixture sur le pinson-----------------------------------------------

# liste des données

nzdata <- list(
  obs = obs,
  PROPNATFOR = PROPNATFOR,
  altitude = altitude,
  year = year,
  julian_date = julian_date,
  npoints = npoints,
  nrepl = nrepl,
  nyears = nyears
)

# choisir les paramètres à suivre

monitor <- c(
  "pmean",
  "fit",
  "fit.new",
  "N",
  "beta.det",
  "beta",
  "delta",
  "mualpha",
  "muyear",
  "det.alpha"
)

# lancer le modèle

jnz <- jags.model(
  file = "scripts/nouvelle_zelande_modele.R",
  data = nzdata,
  n.chains = 3,
  n.adapt = 1000
)

# ce message d'erreur suggère qu'il va falloir initialiser quelque chose - en l'occurrence, N
# on choisit d'initialiser N avec la valeur maximum observée par point

N <- apply(obs, 1, max)

inis <- list(N = N)

jnz <- jags.model(
  file = "donnees/nouvelle_zelande_modele.R",
  data = nzdata,
  inits = inis,
  n.chains = 3,
  n.adapt = 1000
)

update(jnz, n.iter = 5000) # ça va prendre un peu de temps!! bien réfléchir à son modèle...

jnz2 <- coda.samples(jnz,
                    variable.names = monitor,
                    n.iter = 5000,
                    thin = 50)

# mesure de convergence

gelman.plot(jnz2)

# on va aussi regarder les chaines --> on pourrait faire un peu plus stable

color_scheme_set("viridis")
mcmc_trace(
  jnz2,
  pars = c(
    "beta",
    "delta",
    "beta.det[1]",
    "beta.det[2]",
    "mualpha[1]",
    "mualpha[2]",
    "muyear",
    "pmean"
  )
)

# corrélations entre les paramètres

mcmc_pairs(
  jnz2,
  pars = c(
    "beta",
    "delta",
    "beta.det[1]",
    "beta.det[2]",
    "mualpha[1]",
    "mualpha[2]",
    "muyear",
    "pmean"
  ),
  off_diag_args = list(size = 1.5)
)

# mesure d'ajustement

fit.new.nz <- MCMCchains(jnz2, params = c("fit.new"))

fit.nz <- MCMCchains(jnz2, params = c("fit"))

plot(fit.nz, fit.new.nz)
abline(0, 1, col = "blue") # bon ajustement, qu'on peut quantifier:
pval = sum(fit.nz > fit.new.nz) / length(fit.nz) # posterior predictive check pour quantifier l'ajustement - on veut être autour de 0.5
pval

# on veut quand même voir si les estimations sont réalistes

Npred.nz <- MCMCchains(jnz2, params = "N")

med.N <- apply(Npred.nz, 2, median)
hist(med.N) # on a de 0 à 12 pinsons avec surtout des 0 - 6, réaliste

plot(obs.naif, med.N) # bonne relation entre le max sur 3 comptages et le N estimé

# quelle est la probabilité de détection

pdet.nz <- MCMCchains(jnz2, params = "pmean")

hist(pdet.nz)  # probabilité de détection moyenne sur les points
quantile(pdet.nz, p = c(0.025, 0.5, 0.975)) # intervalle de crédibilité

# à chaque passage, on détecte un peu plus d'un pinson sur deux.

# étude des résultats

summary(jnz2)

# et graphiquement

# on commence par regarder les variables de la couche de process (effets de la forêt et de l'altitude)

mcmc_intervals(jnz2, pars = c("beta", "delta"))

# ensuite, on va voir celles de la couche d'obs (intercept et interaction date*année)

mcmc_intervals(jnz2, pars = c("beta.det[1]", "beta.det[2]", "det.alpha"))

# on va maintenant vérifier les effets aléatoires année sur la couche de process

mcmc_intervals(jnz2, pars = c("mualpha[1]", "mualpha[2]"))
