## compute_gam_haplotypes.R
library(data.table)
library(parallel)

N_CORE <- 32

## Input files
p_file <- "Phenotype_manifest.txt"    # should contain at least columns: PGS, PHENO (or similar)
c_file <- "Covariates.txt.gz"         # columns: ID, PC1–PC10, etc.

## Read phenotype manifest
pheno <- fread(p_file)
## adapt these if your column names are different
pgs_list <- pheno$PGS

## Read covariates
cov <- fread(c_file)  # assumes a column "ID" and PCs PC1–PC10
pc_vars <- paste0("PC", 1:10)

run <- function(pgs) {

  ## Read haplotype-based PGS file for this trait
  ## Adapt the path/name if needed (e.g. include chr or other components)
  infile <- paste0( pgs, "_haploid_based_PGS.txt")
  scores <- fread(infile)
  setnames(scores, 1:3, c("ID", "PGS_h1", "PGS_h2"))  # sample, hap1_score, hap2_score

  ## hap1_score and hap2_score must be genome-wide scores. If you run the PGS per chromosome, make sure to sum up the scores for the 22 autosomes into a single genome-wide score.
  ## if you have missing chromosome for some individuals, replace the score for the missing chromosome using the average score for that chromosome computed from the reste of the data.

  ## Merge with covariates
  d <- merge(scores, cov[, c("ID", pc_vars), with = FALSE],
             by = "ID", all.x = TRUE)

  ## Keep complete cases for hap1, hap2, and PCs
  d_cc <- d[complete.cases(PGS_h1, PGS_h2, d[, ..pc_vars])]

  N_ind <- nrow(d_cc)

  ## Raw correlation between haplotypes
  cor_res <- cor.test(d_cc$PGS_h1, d_cc$PGS_h2)

  ## PC-adjusted model: PGS_h2 ~ PGS_h1 + PC1–PC10
  formula_pc <- as.formula(
    paste("PGS_h2 ~ PGS_h1 +", paste(pc_vars, collapse = " + "))
  )
  fit_pc <- lm(formula_pc, data = d_cc)
  cf_pc  <- coef(summary(fit_pc))
  ss_pc  <- cf_pc["PGS_h1", ]

  out <- data.frame(
    PGS    = pgs,
    PHENO  = pheno$PHENO[match(pgs, pheno$PGS)],   # adapt PHENO column name if needed
    N      = N_ind,
    R      = unname(cor_res$estimate),
    RP     = cor_res$p.value,
    beta_pc = unname(ss_pc["Estimate"]),
    se_pc   = unname(ss_pc["Std. Error"]),
    p_pc    = unname(ss_pc["Pr(>|t|)"])
  )

  return(out)
}

## Parallel run over all PGS
results_list <- mclapply(pgs_list, run, mc.cores = N_CORE)
results <- rbindlist(results_list)
fwrite(results, "Haplotype_GAM.txt", sep = "\t")

