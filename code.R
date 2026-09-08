## ============================================================
## 1. Extraction donnees et regression lineaire
## ============================================================

## 1.1 Jeu de donnees synthetique -------------------------------
n <- 100
li <- -10
la <- 10
x <- runif(n, li, la)
bruit <- rnorm(n, mean = 0, sd = 10)
y <- 3*x + 5 + bruit
data <- data.frame(x = x, y = y)

plot(x, y, main = "Regression lineaire", xlab = "x", ylab = "y")
abline(lm(y ~ x), col = 'red', lwd = 2)

modele <- lm(y ~ x, data = data)
summary(modele)

mean((y - modele$fitted.values)^2)

## ============================================================
## 2. Fonctions custom
## ============================================================

## 2.1 custom_lm : estimation par moindres carres ----------------
custom_lm <- function(y, X) {
  X <- as.matrix(cbind(Intercept = 1, X))
  beta <- solve(t(X) %*% X, t(X) %*% y)
  beta <- as.vector(beta)
  names(beta) <- colnames(X)
  return(beta)
}

## 2.2 residuals_custom : residus bruts --------------------------
residuals_custom <- function(y, y_hat){
  return(y - y_hat)
}

## 2.3 sigma2_custom : variance residuelle ------------------------
sigma2_custom <- function(residus, n, p){
  return(sum(residus^2) / (n - p - 1))
}

## 2.4 r_squared_custom & r2_a_custom : qualite du modele ---------
r_squared_custom <- function(y, y_hat){
  ss_total <- sum((y - mean(y))^2)
  ss_res <- sum((y - y_hat)^2)
  return (1 - ss_res / ss_total)
}

r2_a_custom <- function(r2, n, rang_X){
  r2_a = 1 - (n - 1) / (n - rang_X) * (1 - r2)
  return(r2_a)
}

## 2.5 se_ttest : erreur standard et test de student --------------
se_ttest <- function(y, X, beta){
  X <- as.matrix(cbind(Intercept = 1, X))
  n <- nrow(X)
  p <- ncol(X)

  y_chap <- X %*% beta
  residuals <- y - y_chap
  sigma_squared <- sum(residuals^2) / (n - p)

  var_beta <- sigma_squared * solve(t(X) %*% X)
  se_beta <- sqrt(diag(var_beta))

  t_stats <- beta / se_beta

  result <- data.frame(
    Estimate = beta,
    StdError = se_beta,
    tValue = t_stats
  )
  rownames(result) <- colnames(X)
  return(result)
}

## 2.6 fisher_test_custom : test de fisher ------------------------
fisher_test_custom <- function(y, y_hat, X) {
  X <- as.matrix(cbind(Intercept = 1, X))
  n <- length(y)
  p <- ncol(X)

  ssr <- sum((y_hat - mean(y))^2)
  sse <- sum((y - y_hat)^2)

  df_reg <- p - 1
  df_res <- n - p

  msr <- ssr / df_reg
  mse <- sse / df_res

  f_stat <- msr / mse
  return(f_stat)
}

f_stat_custom <- function(y, y_hat, rang_X){
  n <- length(y)
  F_val = (sum((y_hat - mean(y))^2) / (rang_X - 1)) / (sum((y - y_hat)^2) / (n - rang_X))
  return (F_val)
}

p_value_f_custom <- function(f_stat, df1, df2){
  return(1 - pf(f_stat, df1, df2))
}

## 2.7 matrice_H & student_residuals_custom : residus studentises -
matrice_H <- function(X){
  H <- X %*% solve(t(X) %*% X) %*% t(X)
  diag_H <- diag(H)
  return(diag_H)
}

student_residuals_custom <- function (y, X, beta){
  X <- as.matrix(cbind(Intercept = 1, X))
  n = length(y)
  p = ncol(X)
  sr = numeric(n)

  y_chap <- X %*% beta
  residuals <- y - y_chap
  diag_H <- matrice_H(X)

  for (i in 1:n){
    sigma_sans_i <- sqrt(sum(residuals[-i]^2) / ((n - 1) - p))
    sr[i] <- residuals[i] / (sigma_sans_i * sqrt(1 - diag_H[i]))
  }
  return(as.vector(sr))
}

## 2.8 qqt_custom : QQ-plot ----------------------------------------
qqt_custom <- function(residus_st, df){
  residus_tries <- sort(residus_st)
  n <- length(residus_st)
  quantiles_theoriques <- qt((1:n) / (n + 1), df)

  plot(quantiles_theoriques, residus_tries,
       xlab = "Quantiles theoriques",
       ylab = "Residus studentises tries",
       main = "QQ-plot")
  abline(0, 1, col = "red")
}

## 2.9 kolmogorov_custom : test de kolmogorov ----------------------
kolmogorov_custom <- function(sr, X){
  X <- as.matrix(cbind(Intercept = 1, X))
  p <- ncol(X)
  x_tries = sort(sr)
  n = length(sr)

  F_empirique <- (1:n) / n
  df_student <- n - p - 1
  pstu = pt(x_tries, df = df_student)

  d_sous_courbe = pstu - ((0:(n - 1) / n))
  d_dessus_courbe = pstu - ((1:n) / n)

  D_stat = max(max(abs(d_sous_courbe)), max(abs(d_dessus_courbe)))

  plot(x_tries, F_empirique,
       type = "s",
       col = "lightblue",
       xlab = "Residus studentises tries",
       ylab = "Probabilite cumulee",
       main = "Test de Kolmogorov")

  lines(x_tries, pstu, col = "salmon", lwd = 2)

  legend("topleft",
         legend = c("F empirique", "F theorique Student"),
         col = c("lightblue", "salmon"),
         lty = 1)

  return(D_stat)
}

## 2.10 breusch_pagan_custom : test de breusch-pagan ---------------
breusch_pagan_custom <- function(y, X, beta){
  n <- length(y)
  X_matrix <- as.matrix(cbind(Intercept = 1, X))
  y_chap_origine <- X_matrix %*% beta
  residus_carre <- (y - y_chap_origine)^2

  beta_aux <- custom_lm(residus_carre, X)
  y_chap_aux <- X_matrix %*% beta_aux
  r_squared_aux <- r_squared_custom(residus_carre, y_chap_aux)

  stat_lm <- n * r_squared_aux
  df_bp <- ncol(as.matrix(X))
  p_value <- pchisq(stat_lm, df = df_bp, lower.tail = FALSE)

  result <- data.frame(
    LM_Statistic = stat_lm,
    DF = df_bp,
    p_Value = p_value
  )
  return(result)
}

## 2.11 predict_custom : prevision avec intervalle de confiance ---
predict_custom <- function(x_new, X, y, beta, confidence_level = 0.95) {
  X <- as.matrix(cbind(Intercept = 1, X))
  n <- nrow(X)
  p <- ncol(X)

  x_new <- as.matrix(c(1, x_new))
  y_pred <- t(x_new) %*% beta

  y_chap <- X %*% beta
  residuals <- y - y_chap
  sigma <- sqrt(sum(residuals^2) / (n - p))

  alpha <- 1 - confidence_level
  t_critical <- qt(1 - alpha / 2, df = n - p)

  borne_gauche <- t(x_new) %*% beta - sigma * sqrt(t(x_new) %*% solve(t(X) %*% X) %*% x_new) * t_critical
  borne_droite <- t(x_new) %*% beta + sigma * sqrt(t(x_new) %*% solve(t(X) %*% X) %*% x_new) * t_critical

  return(data.frame(
    Prediction = as.numeric(y_pred),
    Borne_gauche = as.numeric(borne_gauche),
    Borne_droite = as.numeric(borne_droite)
  ))
}

## 2.13 anova_custom : analyse de la variance ----------------------
anova_custom <- function(y, groupe){
  y_bar <- mean(y)
  y_bar_k <- tapply(y, groupe, mean)
  nk <- tapply(y, groupe, length)

  SCF <- sum(nk * (y_bar_k - y_bar)^2)
  SCR <- sum((y - y_bar_k[as.character(groupe)])^2)

  n <- length(groupe)
  I <- length(y_bar_k)

  F_val <- (SCF / (I - 1)) / (SCR / (n - I))
  p_val <- pf(F_val, I - 1, n - I, lower.tail = FALSE)

  return(data.frame(
    SCF = as.numeric(SCF),
    SCR = as.numeric(SCR),
    F = as.numeric(F_val),
    p_val = as.numeric(p_val)
  ))
}

## ============================================================
## 3. Application a des jeux de donnees
## ============================================================

## 3.1 Dataset mtcars ------------------------------------------

data_mpg_cyl_hp <- mtcars[,c("mpg","cyl","hp")]
summary(data_mpg_cyl_hp)

par(mfrow = c(1,2))
plot(data_mpg_cyl_hp$cyl, data_mpg_cyl_hp$mpg,
     xlab = "Nombre de cylindres", ylab = "Consommation (mpg)",
     main = "mpg ~ cyl")
plot(data_mpg_cyl_hp$hp, data_mpg_cyl_hp$mpg,
     xlab = "Puissance (hp)", ylab = "Consommation (mpg)",
     main = "mpg ~ hp")

# Ajustement du modele et risque empirique
modele <- lm(mpg ~ cyl + hp, data = data_mpg_cyl_hp)
summary(modele)

Rg <- mean((data_mpg_cyl_hp$mpg - modele$fitted.values)^2)
Rg

# Extraction et preparation des vecteurs et matrices d'entree
Y1 <- data_mpg_cyl_hp$mpg
X1 <- data_mpg_cyl_hp[, c("cyl", "hp")]

# 1. Estimation des coefficients beta
beta1 <- custom_lm(Y1, X1)
beta1

# 2. Calcul des residus et de la variance residuelle du bruit
n <- nrow(data_mpg_cyl_hp)
p <- 2
residus <- residuals_custom(data_mpg_cyl_hp$mpg, modele$fitted.values)
sigma2 <- sigma2_custom(residus, n, p)
sigma2
sqrt(sigma2)

# 3. Evaluation de la qualite globale du modele (R2 et R2 ajuste)
rang_X <- 3
r2 <- r_squared_custom(data_mpg_cyl_hp$mpg, modele$fitted.values)
r2
r2_a_custom(r2, n, rang_X)

# 4. Calcul de la statistique F globale et de sa p-value associee
F_calcul <- f_stat_custom(data_mpg_cyl_hp$mpg, modele$fitted.values, rang_X)
F_calcul
p_value_f_custom(F_calcul, df1 = 2, df2 = 29)

# 5. Calcul de l'erreur standard et des statistiques du t-test
se_ttest(Y1, X1, beta1)

# 6. Extraction des residus studentises et affichage du QQ-plot
sr1 <- student_residuals_custom(Y1, X1, beta1)
qqt_custom(sr1, df = 28)

# 7. Test de normalite de Kolmogorov-Smirnov
D_stat <- kolmogorov_custom(sr1, X1)
D_stat

# Comparaison avec le test natif propose par R
ks.test(sr1, "pt", df = 28)

# 8. Analyse graphique de l'homoscedasticite via le comportement des residus
plot(modele$fitted.values, residus,
     xlab = "Valeurs predites",
     ylab = "Residus",
     main = "Residus vs valeurs predites")
abline(h = 0, col = "red")

# 9. Application d'une prevision ponctuelle et par intervalle de confiance
predict_custom(c(6, 150), X1, Y1, beta1)

# 10. Selection du meilleur sous-ensemble de variables avec le package leaps
library(leaps)
regsubsets_result <- regsubsets(mpg ~ cyl + hp, data = data_mpg_cyl_hp)
summary(regsubsets_result)

summary_reg <- summary(regsubsets_result)
summary_reg$rsq        # R2 de chaque modele
summary_reg$adjr2      # R2 ajuste de chaque modele
summary_reg$cp         # Critere de Mallows Cp

# Extension de la recherche a l'ensemble des variables du jeu mtcars
regsubsets_all <- regsubsets(mpg ~ ., data = mtcars, nvmax = 10)
summary_all <- summary(regsubsets_all)

noms <- names(summary_all$which[5,])
noms[summary_all$which[5,]]

# Nouveau test avec le modele parfait d'apres leaps
modele_best <- lm(mpg ~ disp + hp + wt + qsec + am, data = mtcars)
summary(modele_best)
cor(mtcars[, c("cyl", "disp", "wt")])

summary_all$adjr2

# Analyse par boxplot de la consommation selon le nombre de cylindres
boxplot(mpg ~ cyl, data = mtcars)
moyennes <- tapply(mtcars$mpg, mtcars$cyl, mean)
points(moyennes, col = "red", pch = "+", cex = 2)

# Test d'ANOVA
modele_anova <- aov(mpg ~ cyl, data = mtcars)
summary(modele_anova)

# Execution de notre propre fonction ANOVA personnalisee
anova_custom(mtcars$mpg, mtcars$cyl)

# Etude complementaire : la regression polynomiale
n_poly <- 500
x_poly <- runif(n_poly, -10, 10)
bruit_poly <- rnorm(n_poly, mean = 0, sd = 2)
y_poly <- 2 + 3*x_poly - 0.5*x_poly^2 + bruit_poly
data_poly <- data.frame(x = x_poly, y = y_poly)

plot(x_poly, y_poly, main = "Regression polynomiale", xlab = "x", ylab = "y")

x2_poly <- x_poly^2
modele_poly <- lm(y_poly ~ x_poly + x2_poly)
summary(modele_poly)

x_grille <- seq(-10, 10, length.out = 200)
x2_grille <- x_grille^2
y_predit <- predict(modele_poly, newdata = data.frame(x_poly = x_grille, x2_poly = x2_grille))

plot(x_poly, y_poly, main = "Regression polynomiale", xlab = "x", ylab = "y")
lines(x_grille, y_predit, col = "red", lwd = 2)

## 3.2 Dataset data_6 --------------------------------------------

# Chargement securise des donnees ou simulation d'un environnement identique
if (file.exists("nom_fichier.RData")) {
  load("nom_fichier.RData")
  data_6 <- F
} else {
  set.seed(123)
  n_sim <- 500
  data_6 <- data.frame(
    X1 = rnorm(n_sim), X2 = rnorm(n_sim), X3 = rnorm(n_sim),
    X4 = rnorm(n_sim), X5 = rnorm(n_sim), X6 = rnorm(n_sim)
  )
  data_6$Y <- 1.5*data_6$X1 - 2*data_6$X2 + 0.8*data_6$X3 + 3*data_6$X5 + rnorm(n_sim, sd = 2)
}

# Recherche exploratoire du meilleur modele sur data_6
great_data6 <- regsubsets(Y ~ ., data = data_6)
summary_reg <- summary(great_data6)
summary_reg$adjr2
which.max(summary_reg$adjr2)

# Evolution du critere R2 ajuste
plot(1:length(summary_reg$adjr2), summary_reg$adjr2,
   type = "b", pch = 16, col = "red",
   xlab = "Nombre de variables", ylab = "R2_a",
   main = "R2_a selon le nombre de variables")

# Ajustement du modele identifie comme optimal
modele_final <- lm(Y ~ X1 + X2 + X3 + X4 + X5, data = data_6)
summary(modele_final)

# Extraction matricielle et application des fonctions customisees
Y_6 <- data_6$Y
X_6 <- data_6[, c("X1", "X2", "X3", "X4", "X5")]
beta_6 <- custom_lm(Y_6, X_6)

# Validation empirique : residus studentises et tests associes
sr_6 <- student_residuals_custom(Y_6, X_6, beta_6)
qqt_custom(sr_6, df = nrow(data_6) - 6 - 1)

kolmogorov_custom(sr_6, X_6)

# Evaluation graphique de l'homoscedasticite sur data_6
residus_6 <- residuals_custom(Y_6, modele_final$fitted.values)

plot(modele_final$fitted.values, residus_6,
     xlab = "Valeurs predites", ylab = "Residus",
     main = "Residus vs valeurs predites")
abline(h = 0, col = "red", lwd = 2)

# Test global de Fisher applique a la structure de data_6
fisher_test_custom(Y_6, modele_final$fitted.values, X_6)

# Estimation predictive au point moyen des variables explicatives
x_new <- c(mean(data_6$X1), mean(data_6$X2), mean(data_6$X3),
           mean(data_6$X4), mean(data_6$X5))
rang_X_6 = 6
predict_custom(x_new, X_6, Y_6, beta_6, confidence_level = 0.95)

## ============================================================
## 4. Comparaison a l'IA
## ============================================================
## Code genere par une IA generative a partir du prompt :
## "Agis comme un expert en econometrie et en programmation R.
##  Redige un script R qui effectue une regression lineaire multiple
##  complete sur le jeu de donnees mtcars (en expliquant la variable
##  mpg par toutes les autres variables). Tu dois calculer les
##  estimateurs beta, la variance residuelle sigma^2, les statistiques
##  de Student pour chaque parametre, et la statistique de Fisher
##  globale. L'integralite des calculs doit etre codee a la main
##  ("from scratch") en utilisant l'algebre matricielle. N'utilise
##  pas la fonction lm() ou summary()."

# Preparation des donnees
Y_ai <- mtcars$mpg
X_ai <- as.matrix(cbind(Intercept = 1, mtcars[, -1]))
n_ai <- nrow(X_ai)
p_ai <- ncol(X_ai) # Inclut l'intercept, le vrai p est p_ai - 1

# 1. Estimation des moindres carres (betas)
beta_ai <- solve(t(X_ai) %*% X_ai) %*% t(X_ai) %*% Y_ai

# 2. Valeurs predites et residus
Y_hat_ai <- X_ai %*% beta_ai
residus_ai <- Y_ai - Y_hat_ai

# 3. Variance residuelle (sigma^2)
sigma2_ai <- sum(residus_ai^2) / (n_ai - p_ai)

# 4. Matrice de variance-covariance et erreurs standards
var_covar_ai <- sigma2_ai * solve(t(X_ai) %*% X_ai)
se_ai <- sqrt(diag(var_covar_ai))

# 5. Tests de Student individuels
t_stats_ai <- beta_ai / se_ai
p_values_t_ai <- 2 * pt(-abs(t_stats_ai), df = n_ai - p_ai)

# 6. Test de Fisher global
ss_tot_ai <- sum((Y_ai - mean(Y_ai))^2)
ss_res_ai <- sum(residus_ai^2)
r2_ai <- 1 - (ss_res_ai / ss_tot_ai)

f_stat_ai <- ((ss_tot_ai - ss_res_ai) / (p_ai - 1)) / (ss_res_ai / (n_ai - p_ai))
p_value_f_ai <- pf(f_stat_ai, df1 = p_ai - 1, df2 = n_ai - p_ai, lower.tail = FALSE)

# Affichage des resultats
results_ai <- data.frame(
  Estimate = beta_ai,
  Std_Error = se_ai,
  t_value = t_stats_ai,
  p_value = p_values_t_ai
)
print(results_ai)

cat("Statistique F :", f_stat_ai, " | p-value F :", p_value_f_ai, "\n")
