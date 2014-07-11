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

# The simulations in the file are stored in the data.sims attribute.
# For each Simulation object in data.sims, we can calculate pi for
# a set of loci

for sim in data.sims:
    t =  sim.calc_fit()
    print t
    #pis = sim.calc_pi(loci=[0],multi_skip=False)
    #print pis

  
    
