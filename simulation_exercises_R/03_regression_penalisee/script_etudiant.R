# =============================================================================
# Exercice 03 — Régression pénalisée : Ridge, Lasso, ElasticNet (glmnet)
# Objectif : limiter le surajustement et sélectionner des SNPs
# =============================================================================

# ---- 0. Chargement des données ---------------------------------------------
chemins <- c("../shared_data", "shared_data", "simulation_exercises_R/shared_data")
dossier_data <- chemins[dir.exists(chemins)][1]
if (is.na(dossier_data)) stop("shared_data introuvable : faites setwd() dans le dossier de l'exercice")

for (p in c("glmnet", "pROC")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(glmnet)
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

# glmnet attend une matrice X et un vecteur y
X_train <- as.matrix(train[, snps]); y_train <- train$phenotype
X_test  <- as.matrix(test[, snps]);  y_test  <- test$phenotype

auc_test <- function(proba) as.numeric(auc(roc(y_test, as.vector(proba), quiet = TRUE)))

# ---- 1. Chemin de régularisation --------------------------------------------
fit_ridge <- glmnet(X_train, y_train, family = "binomial", alpha = 0)
fit_lasso <- glmnet(X_train, y_train, family = "binomial", alpha = 1)
par(mfrow = c(1, 2))
plot(fit_ridge, xvar = "lambda", main = "Ridge (alpha = 0)")
plot(fit_lasso, xvar = "lambda", main = "Lasso (alpha = 1)")
par(mfrow = c(1, 1))
# QUESTION 1 : comment les coefficients évoluent-ils quand lambda augmente ?
#              Quelle différence fondamentale entre Ridge et Lasso ?

# ---- 2. Choix de lambda par validation croisée ------------------------------
set.seed(42)
cv_ridge <- cv.glmnet(X_train, y_train, family = "binomial", alpha = 0,   type.measure = "auc", nfolds = 10)
set.seed(42)
cv_lasso <- cv.glmnet(X_train, y_train, family = "binomial", alpha = 1,   type.measure = "auc", nfolds = 10)
set.seed(42)
cv_enet  <- cv.glmnet(X_train, y_train, family = "binomial", alpha = 0.5, type.measure = "auc", nfolds = 10)
plot(cv_lasso)
c(lambda_min = cv_lasso$lambda.min, lambda_1se = cv_lasso$lambda.1se)
# QUESTION 2 : différence entre lambda.min et lambda.1se ? Lequel préférer ?

# ---- 3. Comparaison sur le jeu de test --------------------------------------
m_glm <- glm(phenotype ~ ., data = train[, c(snps, "phenotype")], family = binomial)
resultats <- data.frame(
  modele = c("GLM (non pénalisé)", "Ridge", "Lasso", "ElasticNet (0.5)"),
  AUC_test = c(
    auc_test(predict(m_glm, test, type = "response")),
    auc_test(predict(cv_ridge, X_test, s = "lambda.1se", type = "response")),
    auc_test(predict(cv_lasso, X_test, s = "lambda.1se", type = "response")),
    auc_test(predict(cv_enet,  X_test, s = "lambda.1se", type = "response"))
  ),
  n_SNPs = c(length(snps),
             sum(coef(cv_ridge, s = "lambda.1se")[-1] != 0),
             sum(coef(cv_lasso, s = "lambda.1se")[-1] != 0),
             sum(coef(cv_enet,  s = "lambda.1se")[-1] != 0))
)
print(resultats)
# QUESTION 3 : la pénalisation améliore-t-elle l'AUC test ? Pourquoi ?

# ---- 4. SNPs sélectionnés par le Lasso --------------------------------------
b <- coef(cv_lasso, s = "lambda.1se")[-1, 1]
selection <- data.frame(snp_id = names(b), beta = b)
selection <- merge(selection[selection$beta != 0, ], snp_info[, c("snp_id", "type")])
selection[order(-abs(selection$beta)), ]
table(selection$type)
# QUESTION 4 : le Lasso retrouve-t-il les SNPs "candidate" ? Combien de
#              "background" (faux positifs) garde-t-il ?

# ---- 5. Stabilité de la sélection -------------------------------------------
n_rep <- 20
compte <- setNames(numeric(length(snps)), snps)
for (r in seq_len(n_rep)) {
  i <- sample(nrow(X_train), replace = TRUE)
  cv_r <- cv.glmnet(X_train[i, ], y_train[i], family = "binomial", alpha = 1, nfolds = 5)
  compte <- compte + (coef(cv_r, s = "lambda.1se")[-1, 1] != 0)
}
stab <- data.frame(snp_id = snps, freq = compte / n_rep)
stab <- merge(stab, snp_info[, c("snp_id", "type")])
head(stab[order(-stab$freq), ], 20)
# QUESTION 5 : quels SNPs sont sélectionnés de façon stable (> 80 %) ?

# ---- 6. À vous ---------------------------------------------------------------
# a) Testez alpha dans c(0, 0.25, 0.5, 0.75, 1) : lequel donne la meilleure AUC CV ?
# b) Refaites l'exercice en régression (family = "gaussian") sur true_prs
#    et comparez le R² test avec celui de l'exercice 01.
