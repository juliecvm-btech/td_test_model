# Exercice 04 — Random Forest

**Cibles :** `phenotype` (forêt de probabilités) puis `true_prs` · **Package :** `ranger` (installé automatiquement si absent)

## Objectifs
- Comprendre le bagging et le tirage aléatoire des variables (`mtry`)
- Utiliser l'erreur **OOB** pour régler les hyperparamètres sans jeu de validation
- Voir l'effet du nombre d'arbres et de `min.node.size`
- Interpréter l'importance par permutation

## Déroulé
1. Forêt par défaut, erreur OOB et AUC test
2. AUC en fonction de `num.trees`
3. Grille `mtry` × `min.node.size` évaluée par l'OOB
4. Importance des SNPs (candidats en bleu)
5. Forêt en régression sur le PRS, comparaison avec l'OLS
6. À vous : importance par impureté et par permutation, forêt réduite

## Notions clés
Arbre de décision · bootstrap · bagging · décorrélation des arbres · OOB · importance des variables
