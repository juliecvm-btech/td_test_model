# Exercice 06 — Deep Learning (MLP) et comparaison finale

**Cible :** `phenotype` · **Packages :** `nnet` (obligatoire), `keras` (optionnel, nécessite Python/TensorFlow)

## Objectifs
- Construire un perceptron multicouche et compter ses paramètres
- Standardiser les entrées sans fuite d'information (statistiques calculées sur le train)
- Régulariser avec le weight decay, puis avec dropout et early stopping (keras)
- **Comparer toutes les méthodes du TD** sur le même jeu de test

## Déroulé
1. MLP à une couche cachée avec `nnet`
2. Grille `size` × `decay` sur une validation interne
3. (Optionnel) MLP profond avec keras ; cette partie est ignorée si TensorFlow n'est pas installé
4. Tableau et graphique comparatifs : GLM, Ridge, Lasso, RF, XGBoost, MLP, et le vrai PRS comme référence

## Notions clés
Neurone · fonction d'activation · rétropropagation · weight decay et Ridge · dropout · rapport performance / interprétabilité
