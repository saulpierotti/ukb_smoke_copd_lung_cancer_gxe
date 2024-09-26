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
df[["g"]] <- pgenlibr::Buf(pgen)
pc_names <- colnames(df)[grepl("^22009-0\\.", colnames(df))]
PC <- df[, ..pc_names] |> as.matrix()

# output file
outname <- sprintf("%s.tsv.gwas.gz", pgen_basename)
out_con <- gzfile(outname, "w")

pb <- txtProgressBar(1, nvars, style = 3)
for (i in 1:nvars) {
  setTxtProgressBar(pb, i)

  # read genetic variant
  pgenlibr::Read(pgen, df[["g"]], i)
  var_id <- pgenlibr::GetVariantId(pvar, i)
  var_info <- strsplit(var_id, '_')[[1]]
  chr <- var_info[[1]]
  pos <- var_info[[2]]
  ref <- var_info[[3]]
  alt <- var_info[[4]]

  # fit models
  fit_g_e_gxe <- glm(lung_cancer ~ age + bmi + PC + previous_or_current_smoker + g + g:previous_or_current_smoker, data = df, family = "binomial")
  fit_g_e <- glm(lung_cancer ~ age + bmi + PC + previous_or_current_smoker + g, data = df, family = "binomial")
  fit_e <- glm(lung_cancer ~ age + bmi + PC + previous_or_current_smoker, data = df, family = "binomial")

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
    "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s",
    chr,
    pos,
    var_id,
    ref,
    alt,
    lrt_chisq_gxe,
    lrt_df_gxe,
    pval_gxe,
    lrt_chisq_g,
    lrt_df_g,
    pval_g
  )
  writeLines(lineout, out_con)
}
