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
    sim.calc_fit(pop=1)
    
