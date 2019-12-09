
df <- read.table('../loops/clean_master_requested_loops', header=TRUE, sep="\t")

png("../loops/loop_distance.png", height=1500, width=2500, res=300)
par(mar=c(5, 5, 4, 2) + 0.1,lwd=2)

GH_loops <- df[df$X931.PV.PFC_fdrDonut < -1,]
GH_distances <- GH_loops$y2 - GH_loops$x1
GH_ecdf <- ecdf(log10(GH_distances))
plot(GH_ecdf, verticals=TRUE, do.points=FALSE, col="darkorange", 
     xlim=c(4.5, 6.7), 
     main="", xlab=expression(log[10]* '(bp distance)'),
     col.01line = NULL, cex.lab=1.5, ylab=expression(F[n]*'(x)'))

jSI_loops <- df[df$X712.PV.PFC_fdrDonut < -1,]
jSI_distances <- jSI_loops$y2 - jSI_loops$x1
jSI_ecdf <- ecdf(log10(jSI_distances))
plot(jSI_ecdf, verticals=TRUE, do.points=FALSE, col="dodgerblue", add=TRUE,
     col.01line = NULL)


legend("bottomright", c("931-PV-PFC", "712-PV-PFC"),
       col=c("darkorange", "dodgerblue"),
       lty=rep(1,2), bty="n")

dev.off()


#KS
print("K-S test results:")
print(paste("Two_samples", "p-value", sep="    "))
ksr <- ks.test(GH_distances, jSI_distances)
print(paste("GH_jSI", ksr$p.value, sep="    "))



