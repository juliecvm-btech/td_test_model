# Exercice 03 — Régression pénalisée : Ridge, Lasso, ElasticNet

**Cible :** `phenotype` · **Package :** `glmnet`

## Objectifs
- Comprendre les pénalités L2 (Ridge), L1 (Lasso) et mixte (ElasticNet)
- Lire un chemin de régularisation
- Choisir λ par validation croisée (`lambda.min` ou `lambda.1se`)
- Utiliser le Lasso pour **sélectionner des SNPs** et évaluer la stabilité de cette sélection

## Déroulé
1. Chemins de régularisation Ridge et Lasso
2. `cv.glmnet` avec l'AUC comme critère
3. Comparaison des AUC test avec le GLM non pénalisé
4. SNPs retenus par le Lasso : combien de `candidate` et de `background` ?
5. Stabilité de la sélection par bootstrap
6. À vous : grille sur α, version gaussienne sur `true_prs`

## Notions clés
Biais–variance · pénalité L1/L2 · λ · validation croisée k-fold · sélection de variables
