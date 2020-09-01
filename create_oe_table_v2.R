calc_pval <- function(r, w) {
  # Perform two-sample t-test
  # for jSI vs GH samples
  # Args:
  #   r: row of oe data for all samples
  #   w: boolean if TRUE perform Wilcoxon test
  #      otherwise T-test
  # Return:
  #   p: p-value from t-test
  
  jSI <- as.numeric(c(r["X714_inter_KR_oe"], 
           r["X712.PV.PFC_inter_KR_oe"],
           r["X421_inter_KR_oe"]))
  GH <- as.numeric(c(r["X173_inter_KR_oe"],
          r["X171_inter_KR_oe"],
          r["X73_combined_KR_oe"]))
  if (w) {
    res <- wilcox.test(jSI, GH, alternative = "two.sided")
    p <- res$p.value
  }
  else {
    res <- t.test(jSI, GH, alternative = "two.sided")
    p <- res$p.value
    
  }
  return(p)
}


# Isolated samples (jSI)
df_714 <- read.table("clean_master_requested_loops_714_inter_oe_intxns.txt", sep="\t", header=TRUE)
df_712 <- read.table("clean_master_requested_loops_712-PV-PFC_inter_oe_intxns.txt", sep="\t", header=TRUE)
df_421 <- read.table("clean_master_requested_loops_421_inter_oe_intxns.txt", sep="\t", header=TRUE)

# Group Housed samples (GH)
df_173 <- read.table("clean_master_requested_loops_173_inter_oe_intxns.txt", sep="\t", header=TRUE)
df_171 <- read.table("clean_master_requested_loops_171_inter_oe_intxns.txt", sep="\t", header=TRUE)
df_73 <- read.table("clean_master_requested_loops_73_combined_oe_intxns.txt", sep="\t", header=TRUE)

df_oe <- cbind(df_714, df_712["X712.PV.PFC_inter_KR_oe"], df_421["X421_inter_KR_oe"],
               df_173["X173_inter_KR_oe"], df_171["X171_inter_KR_oe"],
               df_73["X73_combined_KR_oe"])
t_pvals <- c()
w_pvals <- c()

do_wilcoxon <- TRUE

# Parse loop list
for (i in 1:nrow(df_oe)) {
  p <- calc_pval(df_oe[i,], !do_wilcoxon)
  t_pvals <- c(t_pvals, p)
  p <- calc_pval(df_oe[i,], do_wilcoxon)
  w_pvals <- c(w_pvals, p)
}

t_qvals <- p.adjust(t_pvals, method="fdr")
w_qvals <- p.adjust(w_pvals, method="fdr")
df_qvals <- cbind(df_oe, t_pvals, t_qvals, w_pvals, w_qvals)
df_order_p <- df_qvals[order(df_qvals$t_pvals),] 

write.table(df_order_p, "clean_master_requested_loops_oe_intxns_qvals.txt", sep="\t",
            col.names=TRUE, row.names=FALSE, quote=FALSE)

# # Example boxplot
# library(ggplot2)
# ex_df <- read.table("clean_master_requested_loops_oe_intxns_ttest.txt", sep="\t", header=TRUE)
# oe <- as.numeric(ex_df[1,24:29])
# group <- as.factor(c(rep('jSI', 3), rep('GH', 3)))
# 
# box_df <- data.frame(oe, group)
# png("oe_loop_example_boxplot.png", height=1800, width=1100, res=300)
# ggplot(box_df, aes(x=group, y=oe, fill=group)) + 
#   geom_boxplot() + 
#   geom_jitter() +
#   ylab("observed / expected") +
#   ggtitle("Loop\nchr6:53770000-53780000\nchr6:53950000-53960000") +
#   scale_fill_manual(values=c("darkorange", "dodgerblue"))
# dev.off()
#   
# # Histogram of all oe values
# oe_values <- c(as.numeric(ex_df$X282.PV.PFC_inter_KR_oe), 
#                as.numeric(ex_df$X311.PV.PFC_inter_KR_oe), 
#                as.numeric(ex_df$X712.PV.PFC_inter_KR_oe),
#                as.numeric(ex_df$X071.PV.PFC_inter_KR_oe),
#                as.numeric(ex_df$X073.PV.PFC_inter_KR_oe), 
#                as.numeric(ex_df$X931.PV.PFC_inter_KR_oe))
# png("loop_oe_hist.png", width=1500, height=1100, res=300)
# hist(oe_values, col="grey", breaks=10000, main="", xlab="observed/expected values",
#      xlim=c(0,60))
# dev.off()
# 
# 
# 
