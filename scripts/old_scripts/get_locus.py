#!/usr/bin/python
#
#

import sys
import re
import readsfs

f = sys.argv[1]

data = readsfs.FileData()
data.set_file(f)
data.get_sims()

# calculating total distance simulated
#n = 0
#for sim in data.sims:
#    n += sim.command.n_sites
#    for thing in sim.command.gap:
#        n+= thing
#    print n
#    exit()

for sim in data.sims:
    for mut in sim.muts:
        if mut.locus == 212 and mut.fit > 0.:
            i = 0
            for chr in mut.chrs[1]:
                if chr == -1:
                    i = 200
                    break
                i+=1
            i/=200.
            print mut.t_fix[1] -mut.t_init[1], i 

#for locus in sim.loci:
#        for pos in sim.loci[locus]:
#            for mut in sim.loci[locus][pos]:
#                
#                print mut.chrs[0]
