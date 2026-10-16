# =============================================================================
# Exercice 02 — GLM : Gaussien et Logistique
# Objectif : généraliser la régression à la classification (phenotype 0/1)
# =============================================================================

# ---- 0. Chargement des données ---------------------------------------------
chemins <- c("../shared_data", "shared_data", "simulation_exercises_R/shared_data")
dossier_data <- chemins[dir.exists(chemins)][1]
if (is.na(dossier_data)) stop("shared_data introuvable : faites setwd() dans le dossier de l'exercice")

if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
library(pROC)

data <- read.csv(file.path(dossier_data, "data_complet.csv"))
snp_info <- read.csv(file.path(dossier_data, "snp_info_top100.csv"))  # résumé GWAS : SNP, CHR, BP, A1, A2, BETA, SE, P
snp_info$snp_id <- snp_info$SNP
snp_info$type <- ifelse(grepl("^rs_bg", snp_info$SNP), "background", "candidate")
snps <- setdiff(names(data), c("ID", "true_prs", "phenotype"))

table(data$phenotype)   # 1500 témoins / 500 cas
boxplot(true_prs ~ phenotype, data = data, ylab = "true_prs")

set.seed(42)
idx_train <- sample(nrow(data), 0.7 * nrow(data))
train <- data[idx_train, ]
test  <- data[-idx_train, ]

# ---- 1. GLM Gaussien = LM ---------------------------------------------------
m_gauss <- glm(reformulate(snps, "true_prs"), data = train, family = gaussian)
m_lm    <- lm(reformulate(snps, "true_prs"), data = train)
all.equal(coef(m_gauss), coef(m_lm))
# QUESTION 1 : pourquoi les coefficients sont-ils identiques ?
#              (fonction de lien, loi de Y, maximum de vraisemblance vs moindres carrés)

# ---- 2. Pourquoi pas un LM sur une variable 0/1 ? ---------------------------
m_lpm <- lm(reformulate(snps, "phenotype"), data = train)
range(predict(m_lpm, test))
# QUESTION 2 : que remarquez-vous sur les valeurs prédites ?

# ---- 3. Régression logistique -----------------------------------------------
m_logit <- glm(reformulate(snps, "phenotype"), data = train, family = binomial)
summary(m_logit)

# Odds ratios
or <- exp(cbind(OR = coef(m_logit), confint.default(m_logit)))
head(or[order(-abs(log(or[, "OR"]))), ], 10)
# QUESTION 3 : interprétez l'OR du SNP le plus associé.

# ---- 4. Prédiction & seuil de décision --------------------------------------
proba_test <- predict(m_logit, test, type = "response")
hist(proba_test, breaks = 30, main = "Probabilités prédites (test)")

classe_05 <- as.integer(proba_test > 0.5)
table(predit = classe_05, observe = test$phenotype)
mean(classe_05 == test$phenotype)   # accuracy
# QUESTION 4 : une accuracy de 75 % s'obtient en prédisant toujours "0".
#              Pourquoi l'accuracy est-elle trompeuse ici ?

seuil <- mean(train$phenotype)      # seuil = prévalence (0.25)
classe_prev <- as.integer(proba_test > seuil)
tab <- table(predit = classe_prev, observe = test$phenotype)
tab
c(sensibilite = tab["1", "1"] / sum(tab[, "1"]),
  specificite = tab["0", "0"] / sum(tab[, "0"]))

# ---- 5. Courbe ROC & AUC ----------------------------------------------------
roc_test  <- roc(test$phenotype, proba_test, quiet = TRUE)
roc_train <- roc(train$phenotype, predict(m_logit, train, type = "response"), quiet = TRUE)
c(AUC_train = auc(roc_train), AUC_test = auc(roc_test))
plot(roc_test, main = sprintf("GLM logistique — AUC test = %.3f", auc(roc_test)))
# QUESTION 5 : que mesure l'AUC ? Pourquoi ne dépend-elle pas du seuil ?

# ---- 6. Modèle réduit -------------------------------------------------------
top10 <- head(snp_info$SNP[order(snp_info$P)], 10)   # 10 meilleurs SNPs du GWAS externe
m_top10 <- glm(reformulate(top10, "phenotype"), data = train, family = binomial)
auc(roc(test$phenotype, predict(m_top10, test, type = "response"), quiet = TRUE))
AIC(m_logit, m_top10)
# QUESTION 6 : le modèle à 10 SNPs fait-il mieux ou moins bien que celui à 100 ?
#              (top10 vient du GWAS externe : pas de fuite d'information ici)

# ---- 7. À vous ---------------------------------------------------------------
# a) Refaites la sélection des 10 SNPs à partir de tests marginaux sur `train`
#    (et non sur tout le jeu, sinon fuite d'information). Mêmes SNPs que le GWAS ?
# b) Ajustez phenotype ~ true_prs : quelle AUC obtient le "vrai" PRS ?
#    C'est une borne de référence pour tous les exercices suivants.
