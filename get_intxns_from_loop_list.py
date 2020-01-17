#!/usr/bin/env python2

import argparse
import subprocess
import sys
import os

parser = argparse.ArgumentParser(description = 'Print intxn count for pixels of .hic map specified by HiCCUPS loop list')
parser.add_argument('-i', help = 'input .hic file', type = str, required = True)
parser.add_argument('-l', help = 'HiCCUPS loop list', type = str, required = True)
parser.add_argument('-o', help = 'intxn measure <observed/oe>', type = str, default = 'oe')
parser.add_argument('-b', help = 'balancing method <NONE/VC/VC_SQRT/KR>', type = str, default = 'KR')
parser.add_argument('-p', help = 'path to juicer_tools.jar', type=str, 
	default='/home/borrmant/zusers/arima_juicer/opt/juicer/scripts/juicer_tools.jar')
args = parser.parse_args()


def format_coords(l):
	'''
	Extract loop coordinates from file line
	and format for input into juicer_tools.jar 
	dump query

	Args:
		l: str loop line

	Returns:
		fc: str formatted loop coordinates
	'''
	c = l.split('\t')[:6]
	if 'chr' in c[0] or 'chr' in c[3]:
		print 'chromosome naming must be changed'
		sys.exit()
	# Note here the start coordinates are used twice for juicer dump formatting
	fc = c[0] + ':' + c[1] + ':' + c[1] + ' ' + c[3] + ':' + c[4] + ':' + c[4]
	return fc

def loop_resolution(l):
	'''
	Extract loop resolution from file line

	Args:
		l: str loop line

	Returns:
		r: loop resolution in bp
	'''
	c = map(int, l.split('\t')[1:3])
	r = c[1] - c[0]
	if r != 10000:
		print 'WARNING loop is not 10kb res'
		sys.exit()
	if r % 5000 != 0:
		print 'WARNING loop resolution not divisible by 5kb'
		sys.exit()
	return r

def get_intxns(p, o, b, i, c, r):
	'''
	Call juicer_tools.jar dump to get
	intxn counts from .hic file for 
	given pixel loop coordinate
	
	Args:
		p: path to juicer_tools.jar
		o: intxn measure <observed/oe>
		b: balancing method <NONE/VC/VC_SQRT/KR>
		i: input .hic file
		c: pixel loop coordinate
		r: loop resolution
	
	Returns:
		itx: float intxn counts from .hic file
	'''
	command = 'java -jar ' + p + ' dump ' + o + ' ' + b + \
		' ' + i + ' ' + c + ' BP ' + str(r) 
	p = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE)
	output = p.communicate()[0]
	# If output is empty the pixel had zero reads mapped
	if output == '':
		itx = 0.0
	# Make sure only single pixel is found
	elif output.count('\n') != 1:
		print 'ERROR: more than single loop pixel extracted'
		sys.exit()
	else:
		itx = float(output.strip().split('\t')[2])
	return itx


def main():
	
	
	hicfile = os.path.basename(args.i)[:-4]
	loopfile = os.path.basename(args.l)
	# Open output file
	OUT = open(loopfile + '_' + hicfile + '_' + args.o + '_intxns.txt', 'w')
	# Read loop list
	IN = open(args.l, 'r')
	# Format header
	OUT.write(IN.readline().strip() + '\t' + hicfile + '_' + args.b + 
		'_' + args.o + '\n')
	# Parse loop list
	for i, line in enumerate(IN):
		if i % 50 == 0:
			print 'On row ' + str(i)
		coords = format_coords(line)
		res = loop_resolution(line)
		intxns = get_intxns(args.p, args.o, args.b, args.i, coords, res)
		# Write results to file
		OUT.write(line.strip() + '\t' + str(intxns) + '\n')
	IN.close()
	OUT.close()

if __name__ == '__main__':
	main()




