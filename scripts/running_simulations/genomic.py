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
com = command.SFSCommand(prefix='snm') 

# build the command line for snm for a locus on chromosome 1
# uncomment the end of the line to start the simulation with the ancestral sequence
# of course you will need to put the path to your sequence file in place of the 
# path currently in the script
com.genomic(N=100,model='snm',sel=True,dense_dist=0,chr=1,begpos=4190291,endpos=4195373) # withseq=1,seqfile='../../req/hg19/chr1.fa.gz')

# set the location of sfs_code, set the prefix of the out files
com.sfs_code_loc = os.path.join(os.path.expanduser('~'),'rotations/hernandez/software/sfs_code/bin/sfs_code')

# execute the command with a random number
com.execute(rand=randint(1,100000))
