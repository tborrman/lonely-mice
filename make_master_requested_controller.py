#!/usr/bin/env python

import subprocess

for i in range(9):
	subprocess.call('bsub -q long -R "rusage[mem=5000]" -n 1 -W 10:00 -o ' + str(i) + '.out ' +
		'-e ' + str(i) + '.err ./make_master_requested_loops.py -i master_loops_' + 
		'{:02}'.format(i), shell = True)



