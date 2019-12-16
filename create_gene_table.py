#!/usr/bin/env python
import argparse
import sys

parser=argparse.ArgumentParser(description='add gene position to gene expression table')
parser.add_argument('-i', help='gene expression table (ex. nodupsjSI_vs_GH_all_5_2.csv)', type=str, required=True)
args = parser.parse_args()

def main():
	OUT = open('../genes/gene_table.txt', 'w')
	with open(args.i, 'r') as f:
		f.readline()
		for line in f:
			gene = line.split(',')[0]
			found = False
			with open('../genes/mm10_GENCODE_VM11_comp_pseudo_genes.txt', 'r') as e:
				for eline in e:
					esplitline = eline.split('\t')
					egene = esplitline[12]
					if gene == egene:
						found =True
						chrom = esplitline[2]
						start = esplitline[4]
						end = esplitline[5]
						splitline = line.split(',')
						tabline = '\t'.join(splitline)
						OUT.write('\t'.join([chrom, start, end])+'\t'+tabline)
						break
			if not found:
				print gene

if __name__ == '__main__':
	main()

