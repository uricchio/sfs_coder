#!/usr/bin/python
#
#

import sys
import re
from sfs import *
from readsfs import *

file = sys.argv[1]
prefix = sys.argv[2]

n_chr = 500
n_pop = 3
locus_len = 1000
rho = 0.001
# ancestral population size in scaled or non-scaled population? Going with non-scaled for now
N_init = 7300

data = FileData()
data.set_file(file)
data.get_sims()

for sim in data.sims:
    sim.make_map_and_hap(n_chr,n_pop,prefix,locus_len,rho,N_init)
