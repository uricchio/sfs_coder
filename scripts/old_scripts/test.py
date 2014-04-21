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
    pis = sim.calc_pi_by_locus()[1]
    bit = 0
    for pi in pis:
        if pi > 0.:
            bit += 1
    if bit == 0:
        continue
    for pi in pis:    
        print pi,
    print

exit()
for sim in data.sims:
    for thing in sim.calc_pi():
        print thing ,
    print       

#for locus in sim.loci:
#        for pos in sim.loci[locus]:
#            for mut in sim.loci[locus][pos]:
#                
#                print mut.chrs[0]
