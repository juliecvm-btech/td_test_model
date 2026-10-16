# =============================================================================
# Exercice 05 — XGBoost (gradient boosting)
# Objectif : boosting d'arbres, early stopping, importance "Gain"
# =============================================================================

# ---- 0. Chargement des données ---------------------------------------------
chemins <- c("../shared_data", "shared_data", "simulation_exercises_R/shared_data")
dossier_data <- chemins[dir.exists(chemins)][1]
if (is.na(dossier_data)) stop("shared_data introuvable : faites setwd() dans le dossier de l'exercice")

for (p in c("xgboost", "pROC")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(xgboost)
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

# Jeu de validation interne (pris dans train) pour l'early stopping :
# le jeu de test ne doit JAMAIS servir à régler le modèle.
idx_val <- sample(nrow(train), 0.2 * nrow(train))
d_fit <- xgb.DMatrix(as.matrix(train[-idx_val, snps]), label = train$phenotype[-idx_val])
d_val <- xgb.DMatrix(as.matrix(train[idx_val, snps]),  label = train$phenotype[idx_val])
d_train <- xgb.DMatrix(as.matrix(train[, snps]), label = train$phenotype)
d_test  <- xgb.DMatrix(as.matrix(test[, snps]),  label = test$phenotype)

auc_test <- function(proba) as.numeric(auc(roc(test$phenotype, proba, quiet = TRUE)))

params <- list(objective = "binary:logistic", eval_metric = "auc",
               eta = 0.05, max_depth = 3, subsample = 0.8, colsample_bytree = 0.8)

# ---- 1. Boosting sans early stopping ----------------------------------------
set.seed(42)
xgb_long <- xgb.train(params, d_fit, nrounds = 1000,
                      evals = list(fit = d_fit, val = d_val), verbose = 0)
hist_eval <- attributes(xgb_long)$evaluation_log
plot(hist_eval$iter, hist_eval$fit_auc, type = "l", col = "red", ylim = c(0.5, 1),
     xlab = "Itération (nombre d'arbres)", ylab = "AUC")
lines(hist_eval$iter, hist_eval$val_auc, col = "blue")
legend("bottomright", c("fit", "validation"), col = c("red", "blue"), lty = 1)
# QUESTION 1 : décrivez les deux courbes. À partir de quand le modèle surajuste-t-il ?

# ---- 2. Early stopping ------------------------------------------------------
set.seed(42)
xgb_es <- xgb.train(params, d_fit, nrounds = 1000,
                    evals = list(fit = d_fit, val = d_val),
                    early_stopping_rounds = 50, verbose = 0)
n_opt <- xgb.attr(xgb_es, "best_iteration")
n_opt <- as.integer(n_opt) + 1   # best_iteration est indexé à partir de 0
n_opt
# Réentraînement sur tout `train` avec le nombre d'arbres retenu
set.seed(42)
xgb_final <- xgb.train(params, d_train, nrounds = n_opt, verbose = 0)
auc_test(predict(xgb_final, d_test))
# QUESTION 2 : pourquoi utiliser un jeu de validation plutôt que le jeu de test ?

# ---- 3. Effet des hyperparamètres -------------------------------------------
grille <- expand.grid(eta = c(0.01, 0.05, 0.3), max_depth = c(1, 3, 6))
grille$auc_cv <- NA; grille$n_arbres <- NA
for (i in seq_len(nrow(grille))) {
  p <- modifyList(params, as.list(grille[i, c("eta", "max_depth")]))
  set.seed(42)
  cv <- xgb.cv(p, d_train, nrounds = 2000, nfold = 5, early_stopping_rounds = 50, verbose = 0)
  log_cv <- cv$evaluation_log
  grille$auc_cv[i]   <- max(log_cv$test_auc_mean)
  grille$n_arbres[i] <- which.max(log_cv$test_auc_mean)
}
grille[order(-grille$auc_cv), ]
# QUESTION 3 : quel lien entre eta (learning rate) et le nombre d'arbres optimal ?
# QUESTION 4 : max_depth = 1 correspond à des "souches" (stumps) : modèle additif.
#              Est-il compétitif ici ? Qu'est-ce que cela suggère ?

# ---- 4. Importance des variables --------------------------------------------
imp <- xgb.importance(model = xgb_final)
imp <- merge(as.data.frame(imp), snp_info[, c("snp_id", "type")],
             by.x = "Feature", by.y = "snp_id")
imp <- imp[order(-imp$Gain), ]
head(imp[, c("Feature", "Gain", "Cover", "Frequency", "type")], 20)
xgb.plot.importance(xgb.importance(model = xgb_final), top_n = 20, measure = "Gain")
# QUESTION 5 : différence entre Gain, Cover et Frequency ?
#              Pourquoi Gain est-il en général la mesure la plus informative ?

# ---- 5. À vous ---------------------------------------------------------------
# a) Ajoutez une pénalité (lambda = 10, alpha = 1) dans params : effet sur l'AUC ?
# b) Entraînez un XGBoost en régression (objective = "reg:squarederror") sur
#    true_prs et comparez le R² test aux exercices 01, 03 et 04.
