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
  res <- wilcox.test(jSI, GH, alternative = "two.sided")
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

write.table(df_order_p, "clean_master_requested_loops_oe_intxns_wilcoxon.txt", sep="\t",
            col.names=TRUE, row.names=FALSE, quote=FALSE)

# Example boxplot
library(ggplot2)
ex_df <- read.table("clean_master_requested_loops_oe_intxns_ttest.txt", sep="\t", header=TRUE)
oe <- as.numeric(ex_df[1,24:29])
group <- as.factor(c(rep('jSI', 3), rep('GH', 3)))

box_df <- data.frame(oe, group)
png("oe_loop_example_boxplot.png", height=1800, width=1100, res=300)
ggplot(box_df, aes(x=group, y=oe, fill=group)) + 
  geom_boxplot() + 
  geom_jitter() +
  ylab("observed / expected") +
  ggtitle("Loop\nchr6:53770000-53780000\nchr6:53950000-53960000") +
  scale_fill_manual(values=c("darkorange", "dodgerblue"))
dev.off()
  
# Histogram of all oe values
oe_values <- c(as.numeric(ex_df$X282.PV.PFC_inter_KR_oe), 
               as.numeric(ex_df$X311.PV.PFC_inter_KR_oe), 
               as.numeric(ex_df$X712.PV.PFC_inter_KR_oe),
               as.numeric(ex_df$X071.PV.PFC_inter_KR_oe),
               as.numeric(ex_df$X073.PV.PFC_inter_KR_oe), 
               as.numeric(ex_df$X931.PV.PFC_inter_KR_oe))
png("loop_oe_hist.png", width=1500, height=1100, res=300)
hist(oe_values, col="grey", breaks=10000, main="", xlab="observed/expected values",
     xlim=c(0,60))
dev.off()



