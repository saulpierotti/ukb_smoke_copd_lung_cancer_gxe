#!/usr/bin/env Rscript

library("pgenlibr")

args <- commandArgs(trailingOnly = TRUE)
pgen_basename <- args[[1]]
pheno_cov_file <- args[[2]]
var_range_low <- as.numeric(args[[3]])
var_range_high <- as.numeric(args[[4]])

pheno_cov <- read.table(
    pheno_cov_file,
    sep = ",",
    header = TRUE,
    check.names = FALSE
)
pvar <- pgenlibr::NewPvar(sprintf("%s.pvar.zst", pgen_basename))
pgen <- pgenlibr::NewPgen(sprintf("%s.pgen", pgen_basename), pvar = pvar)
psam <- fread(sprintf("%s.psam", pgen_basename))
pheno_cov[["g"]] <- pgenlibr::Buf(pgen)
pc_names <- sprintf("pc%s", 1:10)
PC <- pheno_cov[, pc_names] |> as.matrix()

# output file
outname <- sprintf("%s.tsv.gwas.gz", pgen_basename)
out_con <- gzfile(outname, "w")
header <- "var_id\tmaf\tlrt_chisq_gxe\tlrt_df_gxe\tpval_lrt_gxe\tpval_wald_gxe\tlrt_chisq_g\tlrt_df_g\tpval_lrt_g\tpval_wald_g"
writeLines(header, out_con)

# the E fit does not depend on G  
fit_e <- glm(lung_cancer ~ age + bmi + PC + previous_or_current_smoker, data = pheno_cov, family = "binomial")
ll_e <- logLik(fit_e)

for (i in var_range_low:(var_range_high - 1)) {
  if ((i %% 1000) == 0){
    message(sprintf("%s of %s", i, var_range_high - 1))
  }

  # read genetic variant
  pgenlibr::Read(pgen, pheno_cov[["g"]], i)
  var_id <- pgenlibr::GetVariantId(pvar, i)
  alt_af <- mean(pheno_cov[["g"]], na.rm = TRUE) / 2
  maf <- min(alt_af, 1 - alt_af)

  # fit models
  fit_g_e_gxe <- glm(lung_cancer ~ age + bmi + PC + previous_or_current_smoker + g + g:previous_or_current_smoker, data = pheno_cov, family = "binomial")
  fit_g_e <- glm(lung_cancer ~ age + bmi + PC + previous_or_current_smoker + g, data = pheno_cov, family = "binomial")

  # log likelihoods
  ll_g_e_gxe <- logLik(fit_g_e_gxe)
  ll_g_e <- logLik(fit_g_e)

  # likelihood ratio test - degrees of freedom
  lrt_df_gxe <- attributes(ll_g_e_gxe)[["df"]] - attributes(ll_g_e)[["df"]]
  lrt_df_g <- attributes(ll_g_e)[["df"]] - attributes(ll_e)[["df"]]
  
  # Chi-squared statistics
  lrt_chisq_gxe <- 2 * as.numeric(ll_g_e_gxe - ll_g_e)
  lrt_chisq_g <- 2 * as.numeric(ll_g_e - ll_e)

  # p-values
  pval_lrt_gxe <- (
      if (lrt_df_gxe > 0) pchisq(lrt_chisq_gxe, df = lrt_df_gxe, lower.tail = FALSE) else NA
  )
  pval_lrt_g <- (
      if (lrt_df_g > 0) pchisq(lrt_chisq_g, df = lrt_df_g, lower.tail = FALSE) else NA
  )

  # wald tests
  coef_g_e <- fit_g_e |> summary() |> coef()
  pval_wald_g <- coef_g_e["g", "Pr(>|z|)"]

  coef_g_e_gxe <- fit_g_e_gxe |> summary() |> coef()
  pval_wald_gxe <- coef_g_e_gxe["previous_or_current_smokerTRUE:g", "Pr(>|z|)"]

  # write line to output
  lineout <- sprintf(
    "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
    var_id,
    maf,
    lrt_chisq_gxe,
    lrt_df_gxe,
    pval_lrt_gxe,
    pval_wald_gxe,
    lrt_chisq_g,
    lrt_df_g,
    pval_lrt_g,
    pval_wald_g
  )
  writeLines(lineout, out_con)
}

pgenlibr::ClosePgen(pgen)
close(out_con)