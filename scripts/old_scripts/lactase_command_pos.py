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

import sys
sys.path.append('/netapp/home/lawrence.uricchio/software_working')
import command
import random
import os

alpha = 1000

# initialize a new SFS_CODE command
com = command.SFSCommand(prefix='guten.lact.pos.N2000')

# set the location of sfs_code, set the prefix of the out files
com.sfs_code_loc = os.path.join(os.path.expanduser('~'),'rotations/hernandez/software/sfs_code/bin/sfs_code')

# add a mutation event
mutation = ['--mutation','0.618','P','1','L','212','G',str(alpha)]

# build the command line for the gutenkunst model
com.genomic(N=2000,mutation=mutation)

# execute the command 
com.execute_command(random.randint(0,1000000))
