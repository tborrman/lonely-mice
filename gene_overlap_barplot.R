# Manually entering in gene numbers overlapping loops in sample specific loops (jSI_712, GH_931)
# See line counts for files: 
# X712.PV.PFC_fdrDonut_overlapping_genes_anchor_sample_specific.txt
# X712.PV.PFC_fdrDonut_overlapping_genes_full_loop_sample_specific.txt
# X712.PV.PFC_fdrDonut_overlapping_genes_tss_sample_specific.txt
# X931.PV.PFC_fdrDonut_overlapping_genes_anchor_sample_specific.txt
# X931.PV.PFC_fdrDonut_overlapping_genes_full_loop_sample_specific.txt
# X931.PV.PFC_fdrDonut_overlapping_genes_tss_sample_specific.txt

full_loop <- c(214, 298)
anchor_loop <- c(72, 94)
tss_loop <- c(10, 21)

png("gene_overlap_barplot.png", width=1300, height=1800, res=300)
par(mar=c(5,6,4,1) + 0.1)
df <- data.frame(full_loop, anchor_loop, tss_loop, row.names = c("jSI_712", "GH_931"))
tab <- as.table(t(as.matrix(df)))
barplot(tab, beside=TRUE, legend.text=rownames(tab), ylab="Number of diff expressed genes\n overlapping sample specific loops",
        args.legend=list(x="topleft"))
dev.off()

jSI <- read.table("X712.PV.PFC_fdrDonut_overlapping_genes_anchor_sample_specific.txt", sep="\t")
GH <- read.table("X931.PV.PFC_fdrDonut_overlapping_genes_anchor_sample_specific.txt", sep="\t")

jSI_712 <- jSI$V6
GH_931 <- GH$V6
ttest <- t.test(jSI_712, GH_931)

png("anchor_loops_FC_boxplots.png", width=1200, height=1800, res=300)
par(cex.axis=0.75)
boxplot(jSI_712, GH_931,  col=c("dodgerblue", "darkorange"), ylab="log2( jSI / GH )", 
        names=c("Genes in jSI 712\nspecific_loops", "Genes in GH 931\nspecific loops"))
text(1.5,2, paste("p = ", sprintf('%.2f', round(ttest$p.value, 2)), sep=""))
dev.off()

jSI <- read.table("X712.PV.PFC_fdrDonut_overlapping_genes_full_loop_sample_specific.txt", sep="\t")
GH <- read.table("X931.PV.PFC_fdrDonut_overlapping_genes_full_loop_sample_specific.txt", sep="\t")

jSI_712 <- jSI$V6
GH_931 <- GH$V6
ttest <- t.test(jSI_712, GH_931)

png("full_loops_FC_boxplots.png", width=1200, height=1800, res=300)
par(cex.axis=0.75)
boxplot(jSI_712, GH_931,  col=c("dodgerblue", "darkorange"), ylab="log2( jSI / GH )", 
        names=c("Genes in jSI 712\nspecific_loops", "Genes in GH 931\nspecific loops"))
text(1.5,2, paste("p = ", sprintf('%.2f', round(ttest$p.value, 2)), sep=""))
dev.off()

jSI <- read.table("X712.PV.PFC_fdrDonut_overlapping_genes_tss_sample_specific.txt", sep="\t")
GH <- read.table("X931.PV.PFC_fdrDonut_overlapping_genes_tss_sample_specific.txt", sep="\t")

jSI_712 <- jSI$V6
GH_931 <- GH$V6
ttest <- t.test(jSI_712, GH_931)

png("tss_loops_FC_boxplots.png", width=1200, height=1800, res=300)
par(cex.axis=0.75)
boxplot(jSI_712, GH_931,  col=c("dodgerblue", "darkorange"), ylab="log2( jSI / GH )", 
        names=c("Genes in jSI 712\nspecific_loops", "Genes in GH 931\nspecific loops"))
text(1.5,2, paste("p = ", sprintf('%.2f', round(ttest$p.value, 2)), sep=""))
dev.off()
