#!/usr/bin/python
#$ -e sim.div.log
#$ -o sim.div.log
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -l h_rt=240:00:00
#$ -t 1-100
#$ -l arch=linux-x64
#$ -l mem_free=1G
#$ -l netapp=1G

import command
import os
from random import randint

# initialize a new SFS_CODE command
com = command.SFSCommand(prefix='tennessen.lactase')

# build the command line for the gutenkunst model with lactase
com.genomic(N=100,model='tennessen',sel=False)

# set the location of sfs_code, set the prefix of the out files
com.sfs_code_loc = os.path.join(os.path.expanduser('~'),'rotations/hernandez/software/sfs_code/bin/sfs_code')

# execute the command 
com.execute(rand=randint(1,100000))
