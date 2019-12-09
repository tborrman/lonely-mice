library(pheatmap)
library(RColorBrewer)

df <- read.table('../loops/clean_master_requested_loops', header=TRUE, sep="\t")
qvals <- df[c("X712.PV.PFC_fdrDonut", "X931.PV.PFC_fdrDonut")]
colnames(qvals) <- c("712-PV-PFC", "931-PV-PFC")

png("../loops/qvalue_heatmap.png", width=1200, height=3500, res=300)
pheatmap(qvals, color = rev(brewer.pal(9,"Reds")), show_rownames=FALSE)
dev.off()
