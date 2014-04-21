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
import random
sys.path.append('/netapp/home/lawrence.uricchio/software_working')
from sfs import *
from command import *


if(len(sys.argv) < 3):
    print 'usage: python',sys.argv[0],'<alpha> <delta>'  
    exit()

sge_task_id = 1
if os.environ.has_key('SGE_TASK_ID'):
    sge_task_id = os.environ['SGE_TASK_ID']

alpha = float(sys.argv[1])
delta = float(sys.argv[2])

Nsim = 10

command = SFSCommand()

home_dir = os.path.expanduser('~')

command.set_sfs_code_loc(os.path.join(home_dir, 'rotations/hernandez/software/sfs_code/bin/sfs_code'))

command.build_RHH(alpha=alpha,delta=delta,L0=10**6,theta_neut=0.002,neg_sel_rate = 0.000001)

for i in range(1,Nsim+1):

    command.set_task_id(Nsim*(int(sge_task_id)-1)+i)

    rand = random.randint(0,10**8)

    command.set_prefix('alg2.al'+str(int(alpha))+'.del'+str(delta))

    command.execute_command(rand)
