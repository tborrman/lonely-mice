#!/usr/bin/env Rscript
library(argparse)
library(VennDiagram)

parser <- ArgumentParser(description= "Count overlap of called loops for 2 merged_loops files output by HICCUPS")
parser$add_argument("-f1", help= "merged_loops file 1", type="character", required=TRUE)
parser$add_argument("-f2", help= "merged_loops file 2", type="character", required=TRUE)
parser$add_argument("-r", help= "radius for merging in kb", type="integer", default=0)
args <- parser$parse_args()


count_overlaps <- function(F1_row, F2, radius) {
	# Get data
	chr_1 <- F1_row$chr1
	chr_2 <- F1_row$chr2
	x1 <- F1_row$x1
	y1 <- F1_row$y1
	# Count overlaps
	overlaps <- sum(any((chr_1 == F2$chr1) & (chr_2 == F2$chr2) & ((abs(x1 - F2$x1) <= radius)) & (abs(y1 - F2$y1) <= radius)))
	if (overlaps > 1) {
		print("WARNING : Multiple overlaps found for loop")
		stop()
	}
	return(overlaps)
}

F1 <- read.table(args$f1, header=TRUE, sep="\t")
F2 <- read.table(args$f2, header=TRUE, sep="\t")

total_F1 <- length(F1$x1)
total_F2 <- length(F2$x1)

count <- 0
for (row_idx in 1:nrow(F1)) {
		count <- count + count_overlaps(F1[row_idx,], F2, args$r * 1000)
	}

png("overlap_loops_venn.png", width=2500, height=2500, res=300)
draw.pairwise.venn(area1 = total_F1, area2 = total_F2, cross.area = count, 
	category = c(args$f1, args$f2), fill = c("dodgerblue", "darkorange"),
	cat.pos = c(350,170))
dev.off()

print("Done with Venn")


