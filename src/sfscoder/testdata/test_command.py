# 

from sfscoder import command

com = command.SFSCommand()

com.sfs_code_loc='/netapp/home/lawrence.uricchio/rotations/hernandez/software/sfs_code/bin/sfs_code'

com.build_RHH(delta=0.02,r_within=True)

com.execute()


