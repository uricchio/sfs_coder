#!/usr/bin/python
#
#

import sys
import re
from readsfs import *

f =sys.argv[1]

data = FileData()
data.set_file(f)
data.get_sims()

for sim in data.sims:
    num_sub = sim.count_pos_sub()
    if num_sub[0] != 'NA':
        for thing in num_sub:
            print thing,
        print
    
