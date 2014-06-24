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

# The simulations in the file are stored in the data.sims attribute.
# For each Simulation object in data.sims, we can calculate pi for
# a set of loci

for sim in data.sims:
    sfs =  sim.get_sfs(input_log='/netapp/home/lawrence.uricchio/sfs_coder/scripts/running_simulations/sims/snm/err/log.build_input.txt',start=4194000,stop=4195373)
    print sfs
    #pis = sim.calc_pi(loci=[0],multi_skip=False)
    #print pis

  
    
