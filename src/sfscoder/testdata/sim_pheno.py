from sfscoder import sfs

data = sfs.SFSData(file='/netapp/home/lawrence.uricchio/SKATpower/sims/gutenkunst.cod/gutenkunst.cod.11.txt')

data.get_sims()

for sim in data.sims:
    phenos = sim.sim_pheno(pops=[0,1,2],h_sq=0.01,method='SIMONS',rho=0.0,causal_pop=0)

