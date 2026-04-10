#------------------------------------------------------------------------------#
# Certificat en analyse de données - niveau perfectionnement
# séquence: données opportunistes et indices de diversité
# auteur : Jean-Yves Barnagaud (jean-yves.barnagaud@ephe.psl.eu)
# révisé : 04/2026

# objectif de la séquence : réfléchir à une manière de prendre en compte les
# effets spatiaux et les hétérogénéités d'échantillonnage dans une étude basée
# sur des données opportunistes

# voir : https://doi.org/10.1111/ecog.04122
#------------------------------------------------------------------------------#

library(vegan)
library(mgcv)
library(plyr)

## Données ---------------------------------------------------------------------

# reptiles

msp <- read.table("donnees/reptiles_sites_especes.txt",
                  header = T,
                  sep = "\t")

# Environnement

luse <- read.table(
  "donnees/reptiles_covariables_environnement.txt",
  header = T,
  sep = "\t"
)

# grille de référence

grille <- read.table(
  "donnees/reptiles_grille.txt",
  header = T,
  sep = "\t"
)

## Diversité taxonomique -------------------------------------------------------

# richesse spécifique

sr <- specnumber(msp)

# cartographie

sr.carte <- merge(grille,
                  sr,
                  by.x = "ID",
                  by.y = 0,
                  all = F)

plot(sr.carte$X_centro, sr.carte$Y_centro, cex = sr.carte$y / 10)

## Effort d'échantillonnage ----------------------------------------------------

# exploration

effort.total <- sum(msp)
hist(as.matrix(msp))

# cartographie

effort <- apply(msp, 1, "sum")
effort.carte <- merge(grille,
                      effort,
                      by.x = "ID",
                      by.y = 0,
                      all = F)
plot(effort.carte$X_centro, effort.carte$Y_centro, cex = effort.carte$y /
       1000)

plot(effort.carte$X_centro,
     effort.carte$Y_centro,
     cex = log(effort.carte$y / 100))
plot(effort, sr)

## Courbes de raréfaction ------------------------------------------------------

# courbe de raréfaction sur les points d'observation

crar <- rarecurve(msp,
                  label = F,
                  xlab = "Nombre de données",
                  ylab = "Richesse spécifique")

crar <- rarecurve(
  msp,
  xlim = c(0, 100),
  label = F,
  xlab = "Nombre de données",
  ylab = "Richesse spécifique"
)

# variation du nombre de mailles selon l'effort consenti

barplot(
  c(
    sum(effort > 20),
    sum(effort > 40),
    sum(effort > 60),
    sum(effort > 80),
    sum(effort > 100)
  ),
  names.arg = c(20, 40, 60, 80, 100),
  main = "nombre de sites à xx données",
  ylab = "nombre de sites"
)

# carte du nombre de mailles à un seuil d'effort donné

effort.carte100 <- subset(effort.carte, y >= 100)
plot(effort.carte$X_centro, effort.carte$Y_centro, col = "gray70")
points(
  effort.carte100$X_centro,
  effort.carte100$Y_centro,
  pch = 21,
  bg = "darkred",
  col = "darkred"
)

# échantillonnage des covariables à seuil d'effort donné

effort.cov <-
  merge(effort.carte, luse, by = "ID", all = F)
hist(effort.cov$alti_median)
effort.cov100 =
  subset(effort.cov, y >= 100)
hist(effort.cov100$alti_median,
     add = T,
     col = "gray80")

hist(effort.cov$urbanisation_CORINE)
hist(effort.cov100$urbanisation_CORINE,
     add = T,
     col = "gray80")

# deyxu-le essai

effort.carte50 <- subset(effort.carte, y >= 50)
plot(effort.carte$X_centro, effort.carte$Y_centro, col = "gray70")
points(
  effort.carte50$X_centro,
  effort.carte50$Y_centro,
  pch = 21,
  bg = "darkred",
  col = "darkred"
)

hist(effort.cov$alti_median)
effort.cov50 = subset(effort.cov, y >= 50)
hist(effort.cov50$alti_median,
     add = T,
     col = "gray80")

hist(effort.cov$urbanisation_CORINE)
hist(effort.cov50$urbanisation_CORINE,
     add = T,
     col = "gray80")

# Richesse raréfiée à 50 données
                                                                                                                                                                                                                                                                                                                                                                                                                                                }
crar <- rarecurve(
  msp,
  xlim = c(0, 100),
  label = F,
  xlab = "Nombre de données",
  ylab = "Richesse spécifique"
)
abline(v = 50, lty = "dashed", col = "blue")

msp50 <- msp[which(apply(msp, 1, "sum") >= 50), ]
sr50 <- rarefy(msp50, sample = 50)

# comparaison avec la richesse brute

sr.brute <- specnumber(msp50)
plot(sr.brute, sr50)

# carte des mailles à effort = 50 données

effort.carte50$richesse.rarefiee <- sr50

plot(effort.carte50$y, sr50, main = "richesse raréfiée à 50 points vs nombre de points initial")

## Effets des variables environnementales --------------------------------------

data.50 <-
  merge(effort.cov50[, c("ID",
                         "X_centro",
                         "Y_centro",
                         "alti_median",
                         "urbanisation_CORINE")],
        sr50,
        by.x = "ID",
        by.y = 0,
        all = F)
colnames(data.50)[ncol(data.50)] =
  "richesse.rarefiee"


# relations bivariées

par(mfrow =
      c(1, 2))
plot(data.50$alti_median, data.50$richesse.rarefiee)
plot(data.50$urbanisation_CORINE, data.50$richesse.rarefiee)

# modèle sur la richesse brute

data.50b <- merge(data.50,
                 sr,
                 by.x = "ID",
                 by.y = 0,
                 all = F)
colnames(data.50b)[ncol(data.50b)] = "richesse.brute"

mod.brut <- lm(log(richesse.brute + 1) ~ alti_median + urbanisation_CORINE,
              data = data.50b)

par(mfrow = c(2, 2))
plot(mod.brut)

summary(mod.brut)

# modèle sur la richesse raréfiée

mod.rar <- lm(log(richesse.rarefiee + 1) ~ alti_median + urbanisation_CORINE,
             data = data.50b)

par(mfrow = c(2, 2))
plot(mod.rar)

summary(mod.rar)

## Contrôle de l'autocorrélation spatiale résiduelle

# contrôle visuel (cartographique) 

plot(data.50b$X_centro, data.50b$Y_centro, cex = residuals(mod.rar) * 10)

# correction par smoother

mod.sp.rar <- gam(
  log(richesse.rarefiee + 1) ~ alti_median + urbanisation_CORINE + s(X_centro, Y_centro),
  data = data.50b
)

# vérification des résidus 

  plot(data.50b$X_centro, data.50b$Y_centro, cex = residuals(mod.sp.rar) *
         10)

# inférence

summary(mod.sp.rar)

## Indice de Shannon -----------------------------------------------------------

# communautés raréfiées à 50 données

k <-
    10
  d <- 
    50
  rarlist <-
    replicate(rrarefy(msp50, sample = d), n = k, simplify = F)

# indice de Shannon sur les communautés raréfiées (moyenne et SD)
  
rarH <- ldply(rarlist, .fun = diversity)
rarHm <- apply(rarH, 2, "mean")
rarHs <- apply(rarH, 2, "sd")
  
# indice de Shannon sur matrice observée (non raréfiée)

obsH <- diversity(msp50)

plot(obsH,rarHm)
abline(0,1,col="red")