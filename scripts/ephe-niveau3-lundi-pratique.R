#-----------------------------------------------------------------#
### Formation "Certificat en analyse de données pour les écologues"
### Jean-Yves Barnagaud : jean-yves.barnagaud@ephe.psl.eu ###
# version 2023
# Cas pratique : modèles mixtes
#-----------------------------------------------------------------#

rm(list = ls())

# packages utiles
library(ggplot2)
library(ggeffects)
library(sjPlot)
library(cowplot)
library(lme4)
library(MuMIn)
library(visreg)

#-------------------#
#### Tétras lyre ####
#-------------------#

# données
tetras = read.table("donnees/comptages_tetras_lyre.txt",
                    header = T,
                    sep = "\t")

# création des indices de sites UN, RN, RG pour les sorties graphiques
tetras$UN = factor(paste("ID", tetras$UN, sep = ""))
tetras$RN = factor(paste("ID", tetras$RN, sep = ""))
tetras$RG = factor(paste("ID", tetras$RG, sep = ""))

# on vérifie la disponibilité des données par niveau de facteur
table(tetras$UN)
table(tetras$RN)
table(tetras$RG) # attention au niveau régions, qui n'a que deux niveaux

table(unique(tetras[, c("UN", "RN")])$RN)

table(unique(tetras[, c("RN", "RG")])$RG)

table(unique(tetras[, c("UN", "RG")])$RG)

# modèle 1 : effets aléatoires UN sur la pente et sur l'intercept
logpoules = log(tetras$Poules)
mod.tetras = glmer(
  Jeunes ~ NAOdjfm + I(NAOdjfm^2) + (1 + NAOdjfm + I(NAOdjfm^2) |
                                       RN / UN),
  family = "poisson",
  data = tetras,
  offset = logpoules
)

# NB : si un message "boundary (singular) fit sort, ce n'est pas inquiétant : c'est classique dans un modèle
# complexe comme celui ci, et est lié à la tendance de lme4 à sur-reporter des erreurs. Ici, la raison de ce message est que vu
# l'échelle de variation de certaines variables, on se retrouve en bordure de surface de vraisemblance. Cela n'affecte pas les
# paramètres

# NB2 : si vous avez bien fait le travail, vous avez commencé par ajuster un modèle de Poisson sans l'effet UN, puis avec juste un
# effet aléatoire UN sur l'intercept, puis avec les pentes aléatoires. Nous n'avons pas mis ces étapes dans le script afin de ne pas
# le surcharger, mais elles sont indispensables. Ne cherchez pas à comparer ces différents modèles - il existe des outils de comparaison
# mais ce n'est pas l'objectif : ces premiers modèles vous servent surtout à vérifier la stabilité des paramètres.

# résidus
plot(mod.tetras)

# on va chercher à comprendre le point tout en haut
tetras[which(resid(mod.tetras, type = "pearson") > 5), ]
# il s'agit manifestement d'une UN où la reproduction est particulièrement élevée

# R²
r.squaredGLMM(mod.tetras)
# on constate la forte contribution de l'effet aléatoire

# effets
res.p = ggemmeans(mod.tetras, terms = "NAOdjfm [all]")
plot(res.p, show_residuals = T) + ylim(0, 100)

# pour représenter les effets aléatoires
plot_model(mod.tetras, type = "eff", terms = "NAOdjfm")
plot_model(
  mod.tetras,
  type = "pred",
  terms = c("NAOdjfm [all]", "UN"),
  pred.type = "re"
)

# on rescale
plot(res.p, residuals = T) + ylim(0, 15)

# modèle 2
mod.tetras2 = glmer(
  Jeunes ~ NAOdjfm + I(NAOdjfm^2) + (1 + NAOdjfm + I(NAOdjfm^2) |
                                       RG / RN / UN),
  family = "poisson",
  data = tetras,
  offset = logpoules
)

r.squaredGLMM(mod.tetras2)
# on ne gagne pas grand chose avec l'emboîtement d'effets aléatoires : ce n'est pas très étonnant vu le peu d'information qu'ils contiennent

res.p = ggemmeans(mod.tetras2, terms = "NAOdjfm [all]")
plot(res.p, residuals = T) + ylim(0, 100)

plot_model(mod.tetras2, type = "eff", terms = "NAOdjfm")
plot_model(
  mod.tetras2,
  type = "pred",
  terms = c("NAOdjfm [all]", "UN"),
  pred.type = "re"
)
# le modèle commence à être vraiment complexe, on n'arrive plus à afficher les effets aléatoires emboîtés

# autre vue de l'effet fixe NAO avec les graphiques R de base
library(effects)
ef = effect("NAOdjfm", mod.tetras)
plot(ef)

### comparaison avec le modèle sans effets aléatoires ###
mod.tetras3 = glm(
  Jeunes ~ NAOdjfm + I(NAOdjfm^2),
  family = "poisson",
  data = tetras,
  offset = logpoules
)
plot_model(mod.tetras3, type = "eff", terms = "NAOdjfm")

p1 = plot_model(mod.tetras2, type = "eff", terms = "NAOdjfm")
p2 = plot_model(mod.tetras3, type = "eff", terms = "NAOdjfm")
cowplot::plot_grid(p1, p2)

### optimum climatique ###

# modèle mixte
y = function(x) {
  0.38795 + 0.08736 * x - 0.08166 * x * x
}
xmax.mix = optimize(y, interval = range(tetras$NAOdjfm), maximum = TRUE)
xmax.mix

# glm
y = function(x) {
  0.46064 + 0.04435 * x - 0.08280 * x * x
}
xmax.glm = optimize(y, interval = range(tetras$NAOdjfm), maximum = TRUE)
xmax.glm
