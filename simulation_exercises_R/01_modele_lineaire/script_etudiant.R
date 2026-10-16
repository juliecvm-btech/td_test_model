# =============================================================================
# Exercice 01 — Régression linéaire (OLS)
# Objectif : prédire le PRS continu (true_prs) à partir des 100 SNPs
# Ouvrir ce script depuis son dossier (ou setwd() dans 01_modele_lineaire/)
# =============================================================================

# ---- 0. Chargement des données ---------------------------------------------
chemins <- c("../shared_data", "shared_data", "simulation_exercises_R/shared_data")
dossier_data <- chemins[dir.exists(chemins)][1]
if (is.na(dossier_data)) stop("shared_data introuvable : faites setwd() dans le dossier de l'exercice")

data <- read.csv(file.path(dossier_data, "data_complet.csv"))
snp_info <- read.csv(file.path(dossier_data, "snp_info_top100.csv"))  # résumé GWAS : SNP, CHR, BP, A1, A2, BETA, SE, P
snp_info$snp_id <- snp_info$SNP
snp_info$type <- ifelse(grepl("^rs_bg", snp_info$SNP), "background", "candidate")
snps <- setdiff(names(data), c("ID", "true_prs", "phenotype"))

dim(data)
summary(data$true_prs)
hist(data$true_prs, breaks = 40, main = "Distribution du PRS", xlab = "true_prs")

# Découpage train / test (identique dans tous les exercices)
set.seed(42)
idx_train <- sample(nrow(data), 0.7 * nrow(data))
train <- data[idx_train, ]
test  <- data[-idx_train, ]

# Métriques
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2))
r2   <- function(y, yhat) 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)

# ---- 1. Un seul SNP ---------------------------------------------------------
m1 <- lm(true_prs ~ rs3135388, data = train)
summary(m1)
boxplot(true_prs ~ rs3135388, data = train,
        xlab = "Nombre d'allèles (rs3135388)", ylab = "true_prs")
# QUESTION 1 : interprétez le coefficient de rs3135388 (modèle additif).
# QUESTION 2 : quelle part de variance du PRS ce SNP explique-t-il (R²) ?

# ---- 2. Tous les SNPs -------------------------------------------------------
formule <- reformulate(snps, response = "true_prs")
m_all <- lm(formule, data = train)
summary(m_all)$r.squared

pred_train <- predict(m_all, train)
pred_test  <- predict(m_all, test)

resultats <- data.frame(
  jeu  = c("train", "test"),
  RMSE = c(rmse(train$true_prs, pred_train), rmse(test$true_prs, pred_test)),
  R2   = c(r2(train$true_prs, pred_train),   r2(test$true_prs, pred_test))
)
print(resultats)
# QUESTION 3 : comparez R² train et test. Y a-t-il du surajustement ?

plot(test$true_prs, pred_test, pch = 16, col = "#3366AA80",
     xlab = "PRS observé", ylab = "PRS prédit", main = "LM — jeu de test")
abline(0, 1, col = "red", lwd = 2)

# ---- 3. Quels SNPs comptent ? -----------------------------------------------
coefs <- summary(m_all)$coefficients[-1, ]
coefs <- data.frame(snp_id = rownames(coefs), beta = coefs[, 1], pval = coefs[, 4])
coefs <- merge(coefs, snp_info[, c("snp_id", "type")], by = "snp_id")
coefs <- coefs[order(coefs$pval), ]
head(coefs, 15)
table(significatif = coefs$pval < 0.05 / nrow(coefs), type = coefs$type)
# QUESTION 4 : les SNPs significatifs (Bonferroni) sont-ils des "candidate"
#              ou des "background" ? Qu'en concluez-vous ?

# ---- 4. Diagnostics ----------------------------------------------------------
par(mfrow = c(2, 2)); plot(m_all); par(mfrow = c(1, 1))
# QUESTION 5 : les hypothèses de l'OLS (linéarité, normalité, homoscédasticité)
#              semblent-elles respectées ?

# ---- 5. À vous ---------------------------------------------------------------
# a) Ajustez un modèle avec uniquement les SNPs de type "candidate".
#    Son R² test est-il meilleur ou moins bon que le modèle complet ?
# b) Ajoutez un terme d'interaction entre les deux premiers SNPs. Utile ?
