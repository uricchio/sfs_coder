#!/usr/bin/python
#
#

import sys
import sfs

# an sfs_code output file that we will analyze
f =sys.argv[1]

# initialize a data object and set the file path
data = sfs.SFSData(file=f)

# get all the data from the simulations in the file
data.get_sims()

for sim in data.sims:
    S= sim.calc_S(multi_skip=True)
    print S
    
