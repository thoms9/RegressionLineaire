# Régression Linéaire — Implémentation From Scratch en R

Projet portant sur la régression linéaire, l'un des outils statistiques les plus fondamentaux en pratique. Plutôt que d'utiliser directement les fonctions R existantes (`lm()`, `summary()`...), l'ensemble des estimateurs et tests statistiques est ré-implémenté à partir des formules matricielles vues en cours, afin de comprendre en profondeur chaque étape de la régression linéaire.

**Auteur :** Thomas Begotti

## Démarche

1. Génération d'un jeu de données synthétique respectant le cadre du modèle linéaire (`Y = 3x + 5 + ε`, bruit gaussien centré) pour valider la compréhension du modèle.
2. Application sur le jeu de données réel `mtcars` (consommation des véhicules expliquée par le nombre de cylindres et la puissance).
3. Application sur un second jeu de données `data_6` pour consolider la démarche.
4. Comparaison du code produit avec celui généré par une IA générative à partir du même cahier des charges.

## Fonctions implémentées (`code.R`)

Toutes les fonctions sont codées à partir des formules matricielles (`(X'X)⁻¹X'Y`, etc.), sans appel à `lm()` :

| Fonction | Rôle |
|---|---|
| `custom_lm` | Estimation des coefficients β par moindres carrés |
| `residuals_custom` | Résidus bruts |
| `sigma2_custom` | Estimateur sans biais de la variance résiduelle σ² |
| `r_squared_custom` / `r2_a_custom` | R² et R² ajusté (qualité du modèle) |
| `se_ttest` | Erreur standard des coefficients et statistiques de Student |
| `fisher_test_custom` / `f_stat_custom` / `p_value_f_custom` | Test de Fisher global de significativité |
| `matrice_H` / `student_residuals_custom` | Matrice de projection (leviers) et résidus studentisés |
| `qqt_custom` | QQ-plot des résidus studentisés vs loi de Student |
| `kolmogorov_custom` | Test de normalité de Kolmogorov-Smirnov |
| `breusch_pagan_custom` | Test d'homoscédasticité de Breusch-Pagan |
| `predict_custom` | Prédiction ponctuelle avec intervalle de confiance à 95% |
| `anova_custom` | Analyse de la variance (décomposition SCT = SCF + SCR) |

## Résultats principaux

- **mtcars** (`mpg ~ cyl + hp`) : les coefficients retrouvés par `custom_lm` correspondent exactement à ceux de `lm()` natif (σ ≈ 3.17, risque empirique R(g) ≈ 9.12). Le QQ-plot et le test de Kolmogorov (p-valeur = 0.77) valident l'hypothèse de gaussianité du bruit.
- Une sélection de sous-ensemble de variables (`leaps::regsubsets`) et une analyse ANOVA/boxplot par nombre de cylindres complètent l'étude, ainsi qu'un exemple de régression polynomiale.
- **data_6** (jeu simulé à 6 variables) : le meilleur sous-modèle est identifié via le R² ajusté, puis validé par les mêmes diagnostics (résidus studentisés, Kolmogorov, homoscédasticité, test de Fisher, prédiction au point moyen).
- **Comparaison à l'IA** : à partir d'un prompt identique, une IA générative produit un script condensé sans fonctions intermédiaires. L'approche humaine privilégie la modularité (fonctions réutilisables), ce qui facilite l'application à un nouveau jeu de données sans réécrire les formules matricielles.

## Contenu du dépôt

- `code.R` — ensemble des fonctions custom et des applications aux jeux de données `mtcars` et `data_6`, ainsi que le script de comparaison avec le code généré par IA.
