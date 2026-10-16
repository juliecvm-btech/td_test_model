# =============================================================================
# Exercice 04 — Random Forest (ranger)
# Objectif : découvrir les méthodes d'ensemble (bagging + sous-échantillonnage
#            des variables) et l'importance des variables
# =============================================================================

# ---- 0. Chargement des données ---------------------------------------------
chemins <- c("../shared_data", "shared_data", "simulation_exercises_R/shared_data")
dossier_data <- chemins[dir.exists(chemins)][1]
if (is.na(dossier_data)) stop("shared_data introuvable : faites setwd() dans le dossier de l'exercice")

for (p in c("ranger", "pROC")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(ranger)
library(pROC)

data <- read.csv(file.path(dossier_data, "data_complet.csv"))
snp_info <- read.csv(file.path(dossier_data, "snp_info_top100.csv"))  # résumé GWAS : SNP, CHR, BP, A1, A2, BETA, SE, P
snp_info$snp_id <- snp_info$SNP
snp_info$type <- ifelse(grepl("^rs_bg", snp_info$SNP), "background", "candidate")
snps <- setdiff(names(data), c("ID", "true_prs", "phenotype"))

set.seed(42)
idx_train <- sample(nrow(data), 0.7 * nrow(data))
train <- data[idx_train, ]
test  <- data[-idx_train, ]

# Pour la classification, ranger attend un facteur
train_c <- train[, c(snps, "phenotype")]; train_c$phenotype <- factor(train_c$phenotype)
auc_test <- function(proba) as.numeric(auc(roc(test$phenotype, proba, quiet = TRUE)))

# ---- 1. Une forêt par défaut ------------------------------------------------
rf <- ranger(phenotype ~ ., data = train_c, num.trees = 500,
             probability = TRUE, importance = "permutation", seed = 42)
rf
rf$prediction.error   # erreur OOB (Brier score en mode probability)
proba_test <- predict(rf, test)$predictions[, "1"]
auc_test(proba_test)
# QUESTION 1 : qu'est-ce que l'erreur OOB (out-of-bag) ? Pourquoi est-ce
#              une estimation "gratuite" de l'erreur de généralisation ?

# ---- 2. Nombre d'arbres -----------------------------------------------------
n_arbres <- c(10, 50, 100, 250, 500, 1000)
auc_arbres <- sapply(n_arbres, function(nt) {
  m <- ranger(phenotype ~ ., data = train_c, num.trees = nt, probability = TRUE, seed = 42)
  auc_test(predict(m, test)$predictions[, "1"])
})
plot(n_arbres, auc_arbres, type = "b", log = "x", xlab = "num.trees", ylab = "AUC test")
# QUESTION 2 : augmenter le nombre d'arbres provoque-t-il du surajustement ?

# ---- 3. Réglage de mtry et min.node.size (via l'OOB) ------------------------
grille <- expand.grid(mtry = c(3, 10, 20, 50), min.node.size = c(1, 10, 50))
grille$oob <- apply(grille, 1, function(g) {
  ranger(phenotype ~ ., data = train_c, num.trees = 500, probability = TRUE,
         mtry = g["mtry"], min.node.size = g["min.node.size"], seed = 42)$prediction.error
})
grille[order(grille$oob), ]
meilleur <- grille[which.min(grille$oob), ]
# QUESTION 3 : que contrôle mtry ? Pourquoi le défaut est-il sqrt(p) ≈ 10 ?

rf_best <- ranger(phenotype ~ ., data = train_c, num.trees = 500, probability = TRUE,
                  mtry = meilleur$mtry, min.node.size = meilleur$min.node.size,
                  importance = "permutation", seed = 42)
auc_test(predict(rf_best, test)$predictions[, "1"])

# ---- 4. Importance des variables --------------------------------------------
imp <- data.frame(snp_id = names(rf_best$variable.importance),
                  importance = rf_best$variable.importance)
imp <- merge(imp, snp_info[, c("snp_id", "type")])
imp <- imp[order(-imp$importance), ]
head(imp, 20)
top <- head(imp, 20)
barplot(rev(top$importance), names.arg = rev(top$snp_id), horiz = TRUE, las = 1,
        cex.names = 0.6, col = ifelse(rev(top$type) == "candidate", "#3366AA", "grey70"),
        main = "Importance par permutation (top 20)")
# QUESTION 4 : comparez ce classement avec les SNPs sélectionnés par le Lasso (ex. 03).

# ---- 5. Régression sur le PRS ------------------------------------------------
rf_reg <- ranger(true_prs ~ ., data = train[, c(snps, "true_prs")], num.trees = 500, seed = 42)
pred <- predict(rf_reg, test)$predictions
1 - sum((test$true_prs - pred)^2) / sum((test$true_prs - mean(test$true_prs))^2)
# QUESTION 5 : le R² test de la forêt est-il meilleur que celui de l'OLS (ex. 01) ?
#              Que vous apprend ce résultat sur la nature (additive ?) du PRS simulé ?

# ---- 6. À vous ---------------------------------------------------------------
# a) Comparez importance = "impurity" et "permutation". Mêmes SNPs en tête ?
# b) Entraînez une forêt sur les 20 SNPs les plus importants seulement.
