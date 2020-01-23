#!/usr/bin/env python
import argparse
import sys
import os

# Flush STOUT continuously
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

parser=argparse.ArgumentParser(description='Get genes overlapping loop loci for significantly different loops')
parser.add_argument('-g', help= 'Gene table (ex. gene_table)', type=str, required=True)
parser.add_argument('-l', help='Loop table (ex. clean_master_requested_loops)', type=str, required=True)
parser.add_argument('-c', help='loop significant diff column (ex: pvals)', type=str, required=True)
parser.add_argument('-o', help='overlap protocols: full (full gene : full loop), '+
	'anchor (full gene : anchor loop), tss (tss gene : anchor loop)', type=str, default='tss')
args = parser.parse_args()

def counter(i):
	if i%100 == 0:
		print 'On line: '+str(i)

def overlap(a, b):
        '''
        Return boolean for whether interval b
        overlaps interval a
        '''
        if (a[0] <= b[0] <= a[1]) or (b[0] <= a[0] <= b[1]):
                return True
        else:
                return False

def check_start_below_stop(start, stop):
	if start >= stop:
		print 'ERROR: start >= stop'
		sys.exit()

def is_loop_significant(line, p):
	col_id = {'pvals': 29,
				'qvals': 30}
	splitline = line.rstrip().split('\t')
	sig_val = float(splitline[col_id[p]])
	if sig_val < 0.05:
		return True
	else:
		return False

def write_overlap_genes(geneFile, loop_chrom, loop_interval, outfile, gene_set, p):
	'''
	Write genes overlapping loop calls
	to file
	'''
	loop_chrom = 'chr' + loop_chrom
	with open(geneFile, 'r') as g:
		g.readline()
		for line in g:
			splitline = line.split('\t')
			ENSG = splitline[3]
			if ENSG not in gene_set:
				gene_chrom = splitline[0]
				if p == 'tss':
					gene_start = int(splitline[1])
					gene_stop = int(splitline[1])
				else:
					gene_start = int(splitline[1])
					gene_stop = int(splitline[2])
					check_start_below_stop(gene_start, gene_stop)
				gene_interval = [gene_start, gene_stop]
				
				# Check if gene overlaps loop region 
				if (gene_chrom == loop_chrom) and overlap(loop_interval, gene_interval):
					gene_set.add(ENSG)
					outfile.write(line)
	return gene_set


				

def main():
	# Initialize empty gene list
	genes = set()
	# Significant  loop counter
	num_loops_sig = 0
	# OUTPUT file
	if args.o == 'anchor':
		OUT = open('sig_loops_oe_ttest_overlapping_genes_anchor.txt', 'w' )
	elif args.o == 'full':
			OUT = open('sig_loops_oe_ttest_overlapping_genes_full_loop.txt', 'w')
	elif args.o == 'tss':
			OUT = open('sig_loops_oe_ttest_overlapping_genes_tss.txt', 'w')
	else:
		print 'ERROR: incorrect overlap protocol'
		sys.exit()

	with open(args.l, 'r') as l:
		# Skip header
		l.readline()
		# Parse loop call table
		for i, loop_line in enumerate(l):
			# counter(i)
			# Check if loop exists in cell type
			if is_loop_significant(loop_line, args.c):
				num_loops_sig += 1
				if args.o == 'anchor':
					split_loop_line = loop_line.split('\t')
					loop_chrom = split_loop_line[0]
					loop_start1 = int(split_loop_line[1])
					loop_stop1 = int(split_loop_line[2])
					loop_start2 = int(split_loop_line[4])
					loop_stop2 = int(split_loop_line[5])
					check_start_below_stop(loop_start1, loop_stop1)
					check_start_below_stop(loop_start2, loop_stop2)
					loop1_interval = [loop_start1, loop_stop1]
					loop2_interval = [loop_start2, loop_stop2]
					# Parse gene table for overlapping genes
					genes = write_overlap_genes(args.g, loop_chrom, loop1_interval, OUT, genes, args.o)
					genes = write_overlap_genes(args.g, loop_chrom, loop2_interval, OUT, genes, args.o)

				elif args.o == 'full':				
					split_loop_line = loop_line.split('\t')
					loop_chrom = split_loop_line[0]
					loop_start = int(split_loop_line[1])
					loop_stop = int(split_loop_line[5])
					check_start_below_stop(loop_start, loop_stop)
					loop_interval = [loop_start, loop_stop]
					# Parse gene table for overlapping genes
					genes = write_overlap_genes(args.g, loop_chrom, loop_interval, OUT, genes, args.o)

				elif args.o == 'tss':
					split_loop_line = loop_line.split('\t')
					loop_chrom = split_loop_line[0]
					loop_start1 = int(split_loop_line[1])
					loop_stop1 = int(split_loop_line[2])
					loop_start2 = int(split_loop_line[4])
					loop_stop2 = int(split_loop_line[5])
					check_start_below_stop(loop_start1, loop_stop1)
					check_start_below_stop(loop_start2, loop_stop2)
					loop1_interval = [loop_start1, loop_stop1]
					loop2_interval = [loop_start2, loop_stop2]
					# Parse gene table for overlapping genes
					genes = write_overlap_genes(args.g, loop_chrom, loop1_interval, OUT, genes, args.o)
					genes = write_overlap_genes(args.g, loop_chrom, loop2_interval, OUT, genes, args.o)


				
	OUT.close()
	print 'Number of signifcant loops: ' + str(num_loops_sig) 
	print 'Number of genes overlapping significant loops: '+str(len(genes))

if __name__ == '__main__':
	main()
