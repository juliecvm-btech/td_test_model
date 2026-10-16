# Exercice 05 — XGBoost (gradient boosting)

**Cible :** `phenotype` · **Package :** `xgboost` (API ≥ 3.0 : argument `evals`)

## Objectifs
- Comprendre le boosting : des arbres ajoutés un par un pour corriger les erreurs des précédents
- Observer le surajustement au fil des itérations
- Utiliser l'**early stopping** sur un jeu de validation interne
- Régler `eta` et `max_depth` avec `xgb.cv`
- Interpréter l'importance **Gain**, comparée à Cover et Frequency

## Déroulé
1. 1000 itérations sans arrêt : courbes AUC fit et validation
2. Early stopping, puis réentraînement sur tout le train
3. Grille `eta` × `max_depth` en validation croisée
4. Importance des SNPs
5. À vous : régularisation L1/L2, régression sur `true_prs`

## Notions clés
Boosting · learning rate · early stopping · validation et test · souches (max_depth = 1) · Gain
