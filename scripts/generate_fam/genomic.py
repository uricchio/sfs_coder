#!/usr/bin/python
#$ -e sim.div.log
#$ -o sim.div.log
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -l h_rt=240:00:00
#$ -t 1-1
#$ -l arch=linux-x64
#$ -l mem_free=100M
#$ -l netapp=1G

import os
from sfscoder import command
from random import randint

# initialize a new SFS_CODE command
com = command.SFSCommand(prefix='snm') #outdir=os.path.join(os.path.expanduser('~'),'test'))

# build the command line for the gutenkunst model with lactase
com.genomic(N=500,model='snm',sel=True,dense_dist=0)

# set the location of sfs_code, set the prefix of the out files
com.sfs_code_loc = os.path.join(os.path.expanduser('~'),'rotations/hernandez/software/sfs_code/bin/sfs_code')

# execute the command 
com.execute(rand=randint(1,100000))
