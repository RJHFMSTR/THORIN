library(data.table)
library(parallel)

N_CORE <- 32

## Input files
p_file <- "Phenotype_manifest.txt"    # should contain at least columns: PGS, PHENO (or similar)
m_file <- "Mate_pairs.txt"            # two columns: IDs of partners
c_file <- "Covariates.txt.gz"         # columns: ID, PC1–PC10, etc.

## Read phenotype manifest
pheno <- fread(p_file)
## adapt these if your column names are different
## e.g. pheno$PGS and pheno$PHENO
pgs_list <- pheno$PGS

## Read mate pairs
couples <- fread(m_file)
setnames(couples, 1:2, c("ID1", "ID2"))

## Read covariates
cov <- fread(c_file)  # assumes a column "ID" and PCs PC1–PC10

pc_vars <- paste0("PC", 1:10)

run <- function(pgs) {

  ## Read mate-pair PGS file for this trait
  infile <- paste0("Out/PGS_calc/Scores/", pgs, ".mate_pairs.txt") # one score per individual for the 22 autosomes.
  scores <- fread(infile)
  setnames(scores, 1:2, c("ID", "SCORE"))

  ## Merge PGS for both partners
  d <- merge(couples, scores, by.x = "ID1", by.y = "ID", all.x = TRUE)
  setnames(d, "SCORE", "score_id1")
  d <- merge(d, scores, by.x = "ID2", by.y = "ID", all.x = TRUE)
  setnames(d, "SCORE", "score_id2")

  ## Add PCs for partner 1
  d <- merge(d, cov[, c("ID", pc_vars), with = FALSE],
             by.x = "ID1", by.y = "ID", all.x = TRUE)

  ## Keep complete cases for scores + PCs
  d_cc <- d[complete.cases(score_id1, score_id2, d[, ..pc_vars])]

  N_pairs <- nrow(d_cc)

  ## Raw correlation between partners' PGS
  cor_res <- cor.test(d_cc$score_id1, d_cc$score_id2)

  ## PC-adjusted model: score_id2 ~ score_id1 + PC1–PC10
  formula_pc <- as.formula(
    paste("score_id2 ~ score_id1 +", paste(pc_vars, collapse = " + "))
  )
  fit_pc <- lm(formula_pc, data = d_cc)
  cf_pc  <- coef(summary(fit_pc))
  ss_pc  <- cf_pc["score_id1", ]

  out <- data.frame(
    PGS   = pgs,
    PHENO = pheno$PHENO[match(pgs, pheno$PGS)],   # adapt PHENO column name if needed
    N     = N_pairs,
    R     = unname(cor_res$estimate),
    RP    = cor_res$p.value,
    beta_pc = unname(ss_pc["Estimate"]),
    se_pc   = unname(ss_pc["Std. Error"]),
    p_pc    = unname(ss_pc["Pr(>|t|)"])
  )

  return(out)
}

## Parallel run over all PGS
results_list <- mclapply(pgs_list, run, mc.cores = N_CORE)
results <- rbindlist(results_list)
fwrite(results, "Mate_pairs_GAM.txt", sep = "\t")



