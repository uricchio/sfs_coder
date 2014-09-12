from sfscoder import sfs
from collections import defaultdict

data = sfs.SFSData(file='testdata/test.sfs_code.txt')

data.get_sims()

pops = [0,1]

stats = defaultdict(dict)

def print_stat(stats,name):

    print name,
    for pop in stats[name]:
        print stats[name][pop],
    print
    

for sim in data.sims:

    stats['D'] = sim.calc_tajD(pops=[0,1],start=0,stop=1001)
    stats['pi'] = sim.calc_pi(pops=[0,1],loci=[0])
    stats['S'] = sim.calc_S(pops=[0,1])
    stats['H'] = sim.calc_theta_H(pops=[0,1])
    stats['watt'] = sim.calc_watt(pops=[0,1])
    stats['ZnS'] = sim.calc_ZnS(pops=[0,1])
    stats['f'] = sim.calc_fit(pop=0)
    sim.haplotype()
    stats['h'] = sim.haplo
    
    print_stat(stats,'D')
    print_stat(stats,'H')
    print_stat(stats,'pi')
    print_stat(stats,'S')
    print_stat(stats,'ZnS')
    print_stat(stats,'watt')

    print 'f', 
    for thing in stats['f']:
        print thing,
    print

    print 'h', 
    for haplo in sim.haplo[0]:
        for all in haplo:
            print all,
        print
        break

