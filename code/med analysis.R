setwd("xx/slaughter")

slaughter <- read.table(
  "slaughter.txt",
  header = TRUE,
  sep = "\t",
  na.strings = c("NA")
)
traits <- c(
  "UW100",
  "EN85.100",
  "ESCL100",
  "EST100",
  "ESS100",
  "ESW100",
  "EW100",
  "YW100",
  "YP100"
)
INT <- function(x){
  n <- sum(!is.na(x))
  r <- rank(x, na.last = "keep", ties.method = "average")
  qnorm((r - 0.5) / n)
}
pheno_INT <- slaughter[, c("FID", "IID")]

for(tr in traits){
  pheno_INT[[paste0(tr, "_INT")]] <- INT(slaughter[[tr]])
}
write.table(
  pheno_INT,
  file = "Phenotypes_INT_with_ID.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE



library(dplyr)
library(mediation)


med_results <- list()  # save result

for(snp in sig_snp_cols){
  for(pheno in phenotypes){
    
    # select col
    cols_to_keep <- c("FID","IID", snp, "UW100_INT", pheno, covariates)
    df <- data_all[ , cols_to_keep, drop=FALSE]
    df <- na.omit(df)
    
    # skip few samples
    if(nrow(df) < 20) next
    
    df$SNP_var <- df[[snp]]
    
    med_model <- lm(UW100_INT ~ SNP_var + BW100 + PC1 + PC2 + PC3, data=df)
    
    out_formula <- as.formula(paste0(pheno, " ~ SNP_var + UW100_INT + BW100 + PC1 + PC2 + PC3"))
    out_model <- lm(out_formula, data=df)
    
    # med analysis
    med_fit <- mediate(med_model, out_model, treat="SNP_var", mediator="UW100_INT",
                       boot=TRUE, sims=1000)
    
    # save 
    med_summary <- summary(med_fit)
    med_results[[paste0(snp,"_",pheno)]] <- med_summary
    
    # ADE、ACME、Proportion
    cat("SNP:", snp, "Phenotype:", pheno, "\n")
    cat("ACME (mediation):", med_summary$d0, "p =", med_summary$d0.p, "\n")
    cat("ADE (direct):", med_summary$z0, "p =", med_summary$z0.p, "\n")
    cat("Proportion mediated:", med_summary$n0, "\n\n")
  }
}

med_table <- do.call(rbind, lapply(names(med_results), function(x){
  res <- med_results[[x]]
  data.frame(
    SNP_PHENO = x,
    ACME = res$d0,
    ACME_p = res$d0.p,
    ADE = res$z0,
    ADE_p = res$z0.p,
    Prop_Mediated = res$n0
  )
}))

# output txt
write.table(med_table, "mediation_results_UW100.txt", sep="\t", quote=FALSE, row.names=FALSE)







)
