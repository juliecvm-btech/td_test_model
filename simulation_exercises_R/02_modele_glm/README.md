# Exercice 02 — GLM gaussien et logistique

**Cibles :** `true_prs` (gaussien) puis `phenotype` (binaire, 25 % de cas) · **Package :** `pROC`

## Objectifs
- Voir qu'un GLM gaussien est identique à un `lm()`
- Comprendre pourquoi une régression linéaire ne convient pas à une variable 0/1
- Ajuster une régression logistique et interpréter des odds ratios
- Évaluer un classifieur : matrice de confusion, sensibilité/spécificité, **courbe ROC et AUC**

## Déroulé
1. GLM gaussien et LM
2. Modèle de probabilité linéaire : prédictions hors de [0, 1]
3. Régression logistique et odds ratios
4. Choix du seuil de décision (0.5 ou la prévalence)
5. ROC / AUC sur train et sur test
6. Modèle réduit à 10 SNPs (et piège de la fuite d'information)
7. À vous : AUC de référence avec le vrai PRS

## Notions clés
Fonction de lien · logit · odds ratio · accuracy et classes déséquilibrées · AUC · AIC
