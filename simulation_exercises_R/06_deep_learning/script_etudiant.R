# =============================================================================
# Exercice 06 — Deep Learning : perceptron multicouche (nnet / keras)
# Objectif : réseaux de neurones, puis comparaison de TOUTES les méthodes
# =============================================================================

# ---- 0. Chargement des données ---------------------------------------------
chemins <- c("../shared_data", "shared_data", "simulation_exercises_R/shared_data")
dossier_data <- chemins[dir.exists(chemins)][1]
if (is.na(dossier_data)) stop("shared_data introuvable : faites setwd() dans le dossier de l'exercice")

for (p in c("nnet", "pROC", "glmnet", "ranger", "xgboost"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(nnet)
library(pROC)

data <- read.csv(file.path(dossier_data, "data_complet.csv"))
snps <- setdiff(names(data), c("ID", "true_prs", "phenotype"))

set.seed(42)
idx_train <- sample(nrow(data), 0.7 * nrow(data))
train <- data[idx_train, ]
test  <- data[-idx_train, ]

auc_test <- function(proba) as.numeric(auc(roc(test$phenotype, as.vector(proba), quiet = TRUE)))

# Standardisation : indispensable pour les réseaux de neurones.
# Moyenne et écart-type calculés sur train UNIQUEMENT, puis appliqués au test.
mu <- colMeans(train[, snps]); sigma <- apply(train[, snps], 2, sd)
X_train <- scale(as.matrix(train[, snps]), center = mu, scale = sigma)
X_test  <- scale(as.matrix(test[, snps]),  center = mu, scale = sigma)
y_train <- train$phenotype

# ---- 1. Un MLP simple (1 couche cachée) -------------------------------------
set.seed(42)
mlp <- nnet(X_train, y_train, size = 5, decay = 0, entropy = TRUE,
            maxit = 500, trace = FALSE)
c(AUC_train = as.numeric(auc(roc(y_train, as.vector(predict(mlp, X_train)), quiet = TRUE))),
  AUC_test  = auc_test(predict(mlp, X_test)))
# QUESTION 1 : combien de paramètres compte ce réseau ? (100 entrées, 5 neurones, 1 sortie)
#              Comparez à 2000 observations. Que constatez-vous entre train et test ?

# ---- 2. Régularisation : weight decay & taille ------------------------------
# Validation interne pour choisir les hyperparamètres
idx_val <- sample(nrow(X_train), 0.2 * nrow(X_train))
grille <- expand.grid(size = c(1, 3, 5, 10), decay = c(0, 0.01, 0.1, 1, 5))
grille$auc_val <- apply(grille, 1, function(g) {
  set.seed(42)
  m <- nnet(X_train[-idx_val, ], y_train[-idx_val], size = g["size"], decay = g["decay"],
            entropy = TRUE, maxit = 500, trace = FALSE, MaxNWts = 5000)
  as.numeric(auc(roc(y_train[idx_val], as.vector(predict(m, X_train[idx_val, ])), quiet = TRUE)))
})
grille[order(-grille$auc_val), ][1:10, ]
best <- grille[which.max(grille$auc_val), ]
# QUESTION 2 : quel est l'effet du weight decay ? À quelle méthode de l'ex. 03
#              correspond-il ?

set.seed(42)
mlp_best <- nnet(X_train, y_train, size = best$size, decay = best$decay,
                 entropy = TRUE, maxit = 500, trace = FALSE, MaxNWts = 5000)
auc_mlp <- auc_test(predict(mlp_best, X_test))
auc_mlp

# ---- 3. (Optionnel) MLP profond avec keras ----------------------------------
# Nécessite une installation Python/TensorFlow fonctionnelle (keras::install_keras()).
if (requireNamespace("keras", quietly = TRUE) && keras::is_keras_available()) {
  library(keras)
  tensorflow::set_random_seed(42)
  model <- keras_model_sequential() |>
    layer_dense(units = 32, activation = "relu", input_shape = ncol(X_train),
                kernel_regularizer = regularizer_l2(0.01)) |>
    layer_dropout(rate = 0.3) |>
    layer_dense(units = 16, activation = "relu") |>
    layer_dropout(rate = 0.3) |>
    layer_dense(units = 1, activation = "sigmoid")
  model |> compile(optimizer = optimizer_adam(1e-3), loss = "binary_crossentropy",
                   metrics = list(metric_auc(name = "auc")))
  historique <- model |> fit(
    X_train, y_train, epochs = 200, batch_size = 64, validation_split = 0.2, verbose = 0,
    callbacks = list(callback_early_stopping(monitor = "val_auc", mode = "max",
                                             patience = 20, restore_best_weights = TRUE)))
  plot(historique)
  auc_keras <- auc_test(predict(model, X_test))
  auc_keras
  # QUESTION 3 : rôle du dropout et de l'early stopping ?
} else {
  message("keras/TensorFlow non disponible : partie 3 ignorée.")
  auc_keras <- NA
}

# ---- 4. Comparaison de toutes les méthodes ----------------------------------
library(glmnet); library(ranger); library(xgboost)
Xr_train <- as.matrix(train[, snps]); Xr_test <- as.matrix(test[, snps])
train_c <- train[, c(snps, "phenotype")]; train_c$phenotype <- factor(train_c$phenotype)

m_glm <- glm(phenotype ~ ., data = train[, c(snps, "phenotype")], family = binomial)
set.seed(42); cv_lasso <- cv.glmnet(Xr_train, y_train, family = "binomial", alpha = 1)
set.seed(42); cv_ridge <- cv.glmnet(Xr_train, y_train, family = "binomial", alpha = 0)
rf <- ranger(phenotype ~ ., data = train_c, num.trees = 500, probability = TRUE, seed = 42)
set.seed(42)
xgb <- xgb.train(list(objective = "binary:logistic", eta = 0.05, max_depth = 2,
                      subsample = 0.8, colsample_bytree = 0.8),
                 xgb.DMatrix(Xr_train, label = y_train), nrounds = 300, verbose = 0)
m_ref <- glm(phenotype ~ true_prs, data = train, family = binomial)

comparaison <- data.frame(
  methode = c("GLM logistique", "Ridge", "Lasso", "Random Forest", "XGBoost",
              "MLP (nnet)", "MLP (keras)", "Référence : vrai PRS"),
  AUC_test = c(
    auc_test(predict(m_glm, test, type = "response")),
    auc_test(predict(cv_ridge, Xr_test, s = "lambda.1se", type = "response")),
    auc_test(predict(cv_lasso, Xr_test, s = "lambda.1se", type = "response")),
    auc_test(predict(rf, test)$predictions[, "1"]),
    auc_test(predict(xgb, xgb.DMatrix(Xr_test))),
    auc_mlp, auc_keras,
    auc_test(predict(m_ref, test, type = "response"))
  )
)
comparaison <- comparaison[order(-comparaison$AUC_test), ]
print(comparaison, row.names = FALSE)
ok <- !is.na(comparaison$AUC_test)
par(mar = c(4, 11, 2, 1))
barplot(rev(comparaison$AUC_test[ok]), names.arg = rev(comparaison$methode[ok]),
        horiz = TRUE, las = 1, xlim = c(0.5, 1), xpd = FALSE, xlab = "AUC test",
        main = "Comparaison des méthodes")
par(mar = c(5, 4, 4, 2) + 0.1)

# QUESTION 4 : quelle méthode est la meilleure ? L'écart est-il significatif ?
#              (indice : pROC::roc.test pour comparer deux AUC)
# QUESTION 5 : pourquoi les modèles complexes (RF, XGBoost, MLP) ne battent-ils
#              pas nettement les modèles linéaires pénalisés sur ces données ?
# QUESTION 6 : quel modèle choisiriez-vous en pratique, en tenant compte de la
#              performance, de l'interprétabilité et du coût de calcul ?
