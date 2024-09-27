packages <- c("data.table", "pgenlibr")
install.packages(setdiff(packages, rownames(installed.packages()))) 

library("data.table")
library("pgenlibr")

args <- commandArgs(trailingOnly = TRUE)
pgen_basename <- args[[1]]
pheno_cov_file <- args[[2]]

pheno_cov <- fread(pheno_cov_file)
pvar <- pgenlibr::NewPvar(sprintf("%s.pvar.zst", pgen_basename))
pgen <- pgenlibr::NewPgen(sprintf("%s.pgen", pgen_basename), pvar = pvar)
psam <- fread(sprintf("%s.psam", pgen_basename))
nvars <- pgenlibr::GetVariantCt(pgen)
pheno_cov[["g"]] <- pgenlibr::Buf(pgen)
pc_names <- sprintf("pc%s", 1:10)
PC <- pheno_cov[, ..pc_names] |> as.matrix()

# output file
outname <- sprintf("%s.tsv.gwas.gz", pgen_basename)
out_con <- gzfile(outname, "w")
header <- "var_id\tlrt_chisq_gxe\tlrt_df_gxe\tpval_gxe\tlrt_chisq_g\tlrt_df_g\tpval_g"
writeLines(header, out_con)

nvars <- 10
pb <- txtProgressBar(1, nvars, style = 3)
for (i in 1:nvars) {
  setTxtProgressBar(pb, i)

  # read genetic variant
  pgenlibr::Read(pgen, pheno_cov[["g"]], i)
  var_id <- pgenlibr::GetVariantId(pvar, i)

  # fit models
  fit_g_e_gxe <- glm(lung_cancer ~ age + bmi + PC + previous_or_current_smoker + g + g:previous_or_current_smoker, data = pheno_cov, family = "binomial")
  fit_g_e <- glm(lung_cancer ~ age + bmi + PC + previous_or_current_smoker + g, data = pheno_cov, family = "binomial")
  fit_e <- glm(lung_cancer ~ age + bmi + PC + previous_or_current_smoker, data = pheno_cov, family = "binomial")

  # log likelihoods
  ll_g_e_gxe <- logLik(fit_g_e_gxe)
  ll_g_e <- logLik(fit_g_e)
  ll_e <- logLik(fit_e)

  # likelihood ratio test - degrees of freedom
  lrt_df_gxe <- attributes(ll_g_e_gxe)[["df"]] - attributes(ll_g_e)[["df"]]
  lrt_df_g <- attributes(ll_g_e)[["df"]] - attributes(ll_e)[["df"]]
  
  # Chi-squared statistics
  lrt_chisq_gxe <- 2 * as.numeric(ll_g_e_gxe - ll_g_e)
  lrt_chisq_g <- 2 * as.numeric(ll_g_e - ll_e)

  # p-values
  pval_gxe <- (
      if (lrt_df_gxe > 0) pchisq(lrt_chisq_gxe, df = lrt_df_gxe, lower.tail = FALSE) else NA
  )
  pval_g <- (
      if (lrt_df_g > 0) pchisq(lrt_chisq_g, df = lrt_df_g, lower.tail = FALSE) else NA
  )

  # write line to output
  lineout <- sprintf(
    "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s",
    var_id,
    lrt_chisq_gxe,
    lrt_df_gxe,
    pval_gxe,
    lrt_chisq_g,
    lrt_df_g,
    pval_g
  )
  writeLines(lineout, out_con)
}
