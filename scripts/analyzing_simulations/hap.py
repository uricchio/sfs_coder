#!/usr/bin/python
#
#

import sys
from sfscoder import sfs

# an sfs_code output file that we will analyze
f =sys.argv[1]

# initialize a data object and set the file path
data = sfs.SFSData(file=f)

# get all the data from the simulations in the file
data.get_sims()


# make a haplotype and print it
for sim in data.sims:
    sim.haplotype()
    sim.print_hap(pop=0)
