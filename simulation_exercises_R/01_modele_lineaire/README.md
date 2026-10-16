# Exercice 01 — Régression linéaire (OLS)

**Cible :** `true_prs` (continu) · **Prédicteurs :** 100 SNPs codés 0/1/2 (nombre d'allèles)

## Objectifs
- Ajuster une régression linéaire avec `lm()` sur un puis sur tous les SNPs
- Interpréter un coefficient sous le modèle additif
- Calculer RMSE et R² sur train et sur test, repérer le surajustement
- Lire les diagnostics des résidus

## Déroulé
1. Exploration du PRS
2. Modèle à un SNP (`rs3135388`)
3. Modèle à 100 SNPs, évaluation train/test
4. Coefficients significatifs : SNPs `candidate` ou `background` ?
5. Diagnostics (`plot(lm)`)
6. À vous : modèle restreint aux SNPs candidats, interaction

## Lancer
Ouvrir `script_etudiant.R` dans RStudio, faire `setwd()` sur ce dossier, puis exécuter section par section.

## Notions clés
OLS · moindres carrés · R² · RMSE · découpage train/test · correction de Bonferroni
