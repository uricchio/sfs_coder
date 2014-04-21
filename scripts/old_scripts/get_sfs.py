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
    for thing in sim.get_sfs(pop=2):
        print thing,
    print
    
