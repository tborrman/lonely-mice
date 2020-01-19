calc_pval <- function(r) {
  # Perform two-sample t-test
  # for jSI vs GH samples
  # Args:
  #   r: row of oe data for all samples
  # Return:
  #   p: p-value from t-test
  
  jSI <- as.numeric(c(r["X282.PV.PFC_inter_KR_oe"], 
           r["X311.PV.PFC_inter_KR_oe"],
           r["X712.PV.PFC_inter_KR_oe"]))
  GH <- as.numeric(c(r["X071.PV.PFC_inter_KR_oe"],
          r["X073.PV.PFC_inter_KR_oe"],
          r["X931.PV.PFC_inter_KR_oe"]))
  res <- t.test(jSI, GH, alternative = "two.sided")
  p <- res$p.value
  return(p)
}


# Isolated samples (jSI)
df_282 <- read.table("clean_master_requested_loops_282-PV-PFC_inter_oe_intxns.txt", sep="\t", header=TRUE)
df_311 <- read.table("clean_master_requested_loops_311-PV-PFC_inter_oe_intxns.txt", sep="\t", header=TRUE)
df_712 <- read.table("clean_master_requested_loops_712-PV-PFC_inter_oe_intxns.txt", sep="\t", header=TRUE)

# Group Housed samples (GH)
df_71 <- read.table("clean_master_requested_loops_071-PV-PFC_inter_oe_intxns.txt", sep="\t", header=TRUE)
df_73 <- read.table("clean_master_requested_loops_073-PV-PFC_inter_oe_intxns.txt", sep="\t", header=TRUE)
df_931 <- read.table("clean_master_requested_loops_931-PV-PFC_inter_oe_intxns.txt", sep="\t", header=TRUE)

df_oe <- cbind(df_282, df_311["X311.PV.PFC_inter_KR_oe"], df_712["X712.PV.PFC_inter_KR_oe"],
               df_71["X071.PV.PFC_inter_KR_oe"], df_73["X073.PV.PFC_inter_KR_oe"],
               df_931["X931.PV.PFC_inter_KR_oe"])
pvals <- c()
# Parse loop list
for (i in 1:nrow(df_oe)) {
  p <- calc_pval(df_oe[i,])
  pvals <- c(pvals, p)
}

qvals <- p.adjust(pvals, method="fdr")
df_qvals <- cbind(df_oe, pvals, qvals)
df_order_p <- df_qvals[order(df_qvals$pvals),] 

write.table(df_order_p, "clean_master_requested_loops_oe_intxns.txt", sep="\t",
            col.names=TRUE, row.names=FALSE, quote=FALSE)
