#------------------------------------------------------------------------------#
## simulation : influence de l'autocorrélation spatiale sur une régression ##
# cet exemple a été généré à partir de ChatGPT 
# Jean-Yves Barnagaud (jean-yves.barnagaud@ephe.psl.eu)
# 08 avril 2026
#------------------------------------------------------------------------------#

# paramètres de simulation -----------------------------------------------------

set.seed(42)
n <- 200
s <- seq(0, 10, length.out = n)

beta0 <- 1
beta1 <- 2
epsilon <- rnorm(n, 0, 0.2)

## cas 1 : pas d'effet spatial -------------------------------------------------

# ici, les résidus ne sont pas structurés spatialement

X <- rnorm(n)
Y <- beta0 + beta1 * X + epsilon

model1 <- lm(Y ~ X)
coef(model1)

# pente bonne et résidus non autocorrélés

## cas 2 : ACS résiduelle sans lien avec le prédicteur -------------------------

spatial_effect <- sin(s)

X <- rnorm(n) 
Y <- beta0 + beta1 * X + spatial_effect + epsilon

model2 <- lm(Y ~ X)
coef(model2)

# pente bonne mais les résidus sont spatialement autocorrélés

## cas 3 : le prédicteur est structuré spatialement ----------------------------

X <- scale(s) + rnorm(n, 0, 0.05)

spatial_effect <- 5 * sin(s)

Y <- beta0 + beta1 * X + spatial_effect + epsilon

model3 <- lm(Y ~ X)
coef(model3)

# la pente n'est plus bonne 

## cas 4: effet spatial fort et prédicteur structuré ---------------------------

base_spatial <- scale(sin(s))

X <- 1.5 * base_spatial + rnorm(n, 0, 0.2)

spatial_effect <- 5 * base_spatial + rnorm(n, 0, 0.2)

Y <- beta0 + beta1 * X + spatial_effect + epsilon

model4 <- lm(Y ~ X)
coef(model4)

# la pente n'est plus bonne du tout

## visualisation ---------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(42)

# ----------------------------
# 1. Grille spatiale 2D
# ----------------------------
n_side <- 20
grid <- expand.grid(
  x = seq(0, 10, length.out = n_side),
  y = seq(0, 10, length.out = n_side)
)

n <- nrow(grid)

# Structure spatiale de base
base_spatial <- with(grid, scale(sin(x/2) + cos(y/3)))

# Paramètres
beta0 <- 1
beta1 <- 2

# ----------------------------
# 2. Fonction pour générer un cas
# ----------------------------
generate_case <- function(case_name) {
  
  epsilon <- rnorm(n, 0, 0.2)
  
  if (case_name == "Cas 1") {
    X <- rnorm(n)
    spatial_effect <- 0
    
  } else if (case_name == "Cas 2") {
    X <- rnorm(n)
    spatial_effect <- 3 * base_spatial
    
  } else if (case_name == "Cas 3") {
    X <- scale(grid$x) + rnorm(n, 0, 0.05)
    spatial_effect <- 5 * base_spatial
    
  } else if (case_name == "Cas 4") {
    X <- 1.5 * base_spatial + rnorm(n, 0, 0.2)
    spatial_effect <- 5 * base_spatial + rnorm(n, 0, 0.2)
  }
  
  Y <- beta0 + beta1 * X + spatial_effect + epsilon
  
  df <- grid %>%
    mutate(
      X = as.numeric(X),
      Y = as.numeric(Y)
    )
  
  model <- lm(Y ~ X, data = df)
  
  df$residuals <- resid(model)
  df$case <- case_name
  
  return(df)
}

# ----------------------------
# 3. Générer les 4 cas
# ----------------------------
df_all <- bind_rows(
  generate_case("Cas 1"),
  generate_case("Cas 2"),
  generate_case("Cas 3"),
  generate_case("Cas 4")
)

# ----------------------------
# 4. Format long pour ggplot
# ----------------------------
df_long <- df_all %>%
  select(x, y, case, Y, X, residuals) %>%
  pivot_longer(cols = c(Y, X, residuals),
               names_to = "variable",
               values_to = "value")

# ----------------------------
# 5. Figure finale
# ----------------------------
ggplot(df_long, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  coord_fixed() +
  scale_fill_viridis_c() +
  facet_grid(case ~ variable) +
  labs(
    title = "Structures spatiales : Y, X et résidus selon les cas",
    fill = "Valeur"
  ) +
  theme_minimal(base_size = 12)

#+
 # scale_fill_viridis_c(option = "C", limits = range(df_long$value))

## graphiques attendu vs observé -----------------------------------------------

# ----------------------------
# 1. Grille et structure
# ----------------------------
n_side <- 20
grid <- expand.grid(
  x = seq(0, 10, length.out = n_side),
  y = seq(0, 10, length.out = n_side)
)

n <- nrow(grid)
base_spatial <- with(grid, scale(sin(x/2) + cos(y/3)))

beta0 <- 1
beta1 <- 2

# ----------------------------
# 2. Fonction génération
# ----------------------------
generate_case <- function(case_name) {
  
  epsilon <- rnorm(n, 0, 0.2)
  
  if (case_name == "Cas 1") {
    X <- rnorm(n)
    spatial_effect <- 0
    
  } else if (case_name == "Cas 2") {
    X <- rnorm(n)
    spatial_effect <- 3 * base_spatial
    
  } else if (case_name == "Cas 3") {
    X <- scale(grid$x) + rnorm(n, 0, 0.05)
    spatial_effect <- 5 * base_spatial
    
  } else if (case_name == "Cas 4") {
    X <- 1.5 * base_spatial + rnorm(n, 0, 0.2)
    spatial_effect <- 5 * base_spatial + rnorm(n, 0, 0.2)
  }
  
  Y <- beta0 + beta1 * X + spatial_effect + epsilon
  
  df <- data.frame(X = as.numeric(X), Y = as.numeric(Y))
  
  model <- lm(Y ~ X, data = df)
  
  df$case <- case_name
  df$Y_true <- beta0 + beta1 * df$X  # vraie relation
  
  # prédictions modèle observé + IC
  pred <- predict(model, interval = "confidence")
  df$Y_hat <- pred[, "fit"]
  df$Y_hat_low <- pred[, "lwr"]
  df$Y_hat_high <- pred[, "upr"]
  
  return(df)
}

# ----------------------------
# 3. Générer les données
# ----------------------------
df_all <- bind_rows(
  generate_case("Cas 1"),
  generate_case("Cas 2"),
  generate_case("Cas 3"),
  generate_case("Cas 4")
)

# ----------------------------
# 4. Graphique
# ----------------------------
ggplot(df_all, aes(x = X, y = Y)) +
  
  # Points
  geom_point(alpha = 0.4) +
  
  # IC modèle estimé
  geom_ribbon(aes(ymin = Y_hat_low, ymax = Y_hat_high),
              fill = "grey70", alpha = 0.5) +
  
  # Droite estimée (biaisée ou non)
  geom_line(aes(y = Y_hat), color = "blue", linewidth = 1.2) +
  
  # Vraie droite
  geom_line(aes(y = Y_true), color = "red", linewidth = 1.2, linetype = "dashed") +
  
  facet_wrap(~ case, scales = "free") +
  
  labs(
    title = "Effet de l'autocorrélation spatiale sur la pente estimée",
    subtitle = "Rouge = vraie relation | Bleu = régression estimée",
    x = "X",
    y = "Y"
  ) +
  
  theme_minimal(base_size = 13)
