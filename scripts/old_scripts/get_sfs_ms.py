#!/usr/bin/python
#
#

import ms
import sys
from readsfs import *

f =sys.argv[1]

data = msData(file=f)
data.get_sims()

s = int(sys.argv[2])
e = int(sys.argv[3])

for sim in data.sims:
    sim.get_sfs(start = s, end=e)
    for thing in sim.sfs:
        print thing,
    print
