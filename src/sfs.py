#!/usr/bin/python

import subprocess
import os
import math
import command
import ms
import re
import matplotlib
from collections import defaultdict
from os.path import exists

class Simulation:
 
    """
    A class to store the data from SFS_CODE simulations.

    * Attributes:
    
       * *self.command = command.SFSCommand()*

       * *self.data = ''*
     
       * *self.loci = defaultdict(lambda: defaultdict(list))*

          A dictionary of dictionaries indexed by locus and position.
          Each dictionary of dictionaries is a list of Mutation
          objects that occur at the corresponding locus and position.
   
          e.g., self.loci[0][0] is a list of Mutation objects that
          occur at the first position in the first locus.

          Note that these keys will only exist in self.loci if there
          were mutations at this particular point in the sequence
          in the sample from the simulation. 

       * *self.muts = []*

          A list of all the Mutation objects in the simulation.
 
    """

    def __init__(self):
   
        self.command = command.SFSCommand()
        self.data = ''
        self.loci = defaultdict(lambda: defaultdict(list))
        self.muts = []
        self.haplo = {}
        self.causal = {}

    def set_data(self,data):
        self.data = data

    def set_command(self, command):
        self.command = command

    def make_muts(self):

        arr = self.data.split(';')
        for snp in arr:
            
            # test if the mutation has no data
            if len(snp) < 1:
                continue
            
            mut_data = snp.split(',')
            
            newmut = Mutation()
            newmut.set_all(mut_data,self.command.P,self.command.n)

            self.loci[int(newmut.locus)][int(newmut.pos)].append(newmut)

        self.merge_muts()

        return

    def merge_muts(self):

        def merge(m,multi):
            if len(m) < 2:
                newmut = m[0]
                newmut.multiallelic = multi
                return newmut
        
            chrs = defaultdict(dict)

            t_fix = {}
            t_init = {}
            fixed_pop = {}
            pops_numchr = {}
            for mut in m:
                for pop in mut.chrs:
                   t_fix[pop] = mut.t_fix[pop]
                   t_init[pop] = mut.t_init[pop]
                   fixed_pop[pop] = mut.fixed_pop[pop]
                   pops_numchr[pop] = mut.pops_numchr[pop]
                   for chr in mut.chrs[pop]:
                        chrs[pop][chr] = True

            newmut = m[0]
            newmut.chrs = chrs
            newmut.t_fix = t_fix
            newmut.t_init = t_init
            newmut.fixed_pop = fixed_pop
            newmut.pops_numchr = pops_numchr
            newmut.multiallelic = multi
            return newmut

        newloci = defaultdict(lambda: defaultdict(list))

        for locus in sorted(self.loci):
            for pos in sorted(self.loci[locus]):
                muts =defaultdict(list)
                for mut in self.loci[locus][pos]:
                    pop = -1
                    for p in mut.chrs:
                        pop =p
                        break
                    # all the muts with same t_init at same location
                    # are pushed into a single list, then merged into a single
                    # entity    
                    muts[mut.t_init[pop]].append(mut)
                
                # check whether site is multiallelic

                num_alleles = 0
                for t in muts:
                    for mut in muts[t]:
                        for pop in mut.t_fix:
                            if float(mut.t_fix[pop]) == self.command.TE:
                                num_alleles += 1
                                break        
                
                multi = False
                if num_alleles > 1:
                    multi = True
                
                for t in muts:
                    newmut =  merge(muts[t],multi)
                    self.muts.append(newmut)
                    newloci[locus][pos].append(newmut)
        self.loci = newloci  
           
        #for mut in self.muts:
        #    print mut.multiallelic
        
    def calc_S(self,multi_skip=True,loci=[]):
     
        """
        calculate the number of segretating sites within all
        populations.

        * Parameters:
          
          * *multi_skip = True*
             skip sites that are more than biallelic if true
          * *loci = []*
             A list of loci over which to calculate the number of segregating
             sites.  Uses all loci if this is left blank.
        """

        if len(loci) == 0:
            loci = self.loci
        nums = []

        for pop in range(0, self.command.n_pops):
            num = 0
            for locus in loci:
               for pos in self.loci[locus]:
                   for mut in self.loci[locus][pos]:
                       if pop not in mut.chrs:
                           continue
                       if float(mut.t_fix[pop]) < self.command.TE:
                           continue
                       if multi_skip:
                           if mut.multiallelic:
                               continue
                       num += 1
            nums.append(num)    
        
        return nums  

    def calc_tajD(self, loci=[]):
        return

    def calc_ZnS(self, loci=[]):
        return 
 
    def calc_pi(self,multi_skip=True,loci=[]):
        
        """ 
        calculate the mean pairwise diversity per site
        bewteen pairs of sequences across a set of loci.
        If the loci parameter is left undefined by the
        user, then this method calculates :math:`\pi` 
        over all loci in the simulation.

        * Parameters:
              
          * *multi_skip=True*
            
             If True, skip loci that are multiallelic.
             Otherwise lump all the derived alleles together.

          * *loci=[]*
            
             An array of loci over which to calculate 
             :math:`\pi`.  If left blank, all loci
             are used in the calculation.
         
        * Return value: An array of values of :math:`\pi`
          values indexed by population number.
  
        """

        for locus in loci:
            if locus not in self.loci:
                print "Locus",locus,"is not in self.loci!"
                print "restrict input to calc_pi (loci=[...]) to loci that",
                print "were actually simulated"
                exit()

        pis = []
        tot = 0
        if len(loci) == 0:
            loci=self.loci  
            tot = self.command.n_sites
        else:
            locidict = {}
            for locus in loci:
                tot += self.command.L[locus]
                locidict[locus] = self.loci[locus]

            loci = locidict
                
        for pop in range(0,self.command.n_pops):
            pi = 0.0
            for locus in loci:
                for pos in loci[locus]:
                    sum=0
                    for mut in loci[locus][pos]:
                        if pop not in mut.chrs:
                           continue
                        if mut.fixed_pop[pop] == True:
                           continue
                        if mut.multiallelic == True and multi_skip == True:
                           continue
                        sum += mut.pops_numchr[pop]
                    pi_var = (self.command.n[pop]-sum)*sum
                    pi += pi_var
            pi = 2. * pi / (tot*self.command.n[pop]*(self.command.n[pop]-1))
            pis.append(pi)
        return pis
   
    def calc_pi_by_locus(self):

        """  
        calculate the value of :math:`\pi` independently for each locus.

        
        * Return value: An array of arrays of pi values, indexed by 
          population and then locus number.

          e.g., if the return value is stored in the variable pi, 
          pi[0][1] is the value of :math:`\pi` in population 0 at locus 1.

        """


        pis = [[0.0 for i in range(0, len(self.command.L))] for j in range(0,self.command.n_pops)]

        for pop in range(0, self.command.n_pops):
            #pi = [0.0 for i in range(0,len(self.command.L))]
            for locus in self.loci:
                for pos in self.loci[locus]:
                    for mut in self.loci[locus][pos]:
                        if pop not in mut.chrs:
                           continue
                        if mut.fixed_pop[pop] == True:
                           continue
                        sum = mut.pops_numchr[pop]
                        pis[pop][locus] += (self.command.n[pop]-sum)*sum
                pis[pop][locus] = 2.*pis[pop][locus] / (self.command.L[locus]*self.command.n[pop]*(self.command.n[pop]-1))

        return pis

    def calc_fit(self,pop=0):

        """
        calculate the fitness of the sampled chromosomes
        within a population.

        * Parameters:
  
          * *pop=0* 
             population number
             
        """


        fit = {}             

        if pop not in self.command.n:
            return

        for chr in range(0,self.command.n[pop]):
            fit[chr] = 1
        
        for i in range(0,len(self.muts)):
           
            #print i, self.muts[i]

            mut = self.muts[i]

            if mut.fit == 0.:
                continue

            if pop not in mut.t_fix:
                continue

            if mut.t_fix[pop]<self.command.TE:
                continue

            #print pop, mut.t_fix[pop]

            #print pop, mut.locus, mut.pos, mut.fit, mut.t_fix[pop]

            if -1 in mut.chrs[pop]:
                if self.command.Z == 0:
                    for chr in range(0,self.command.n[pop]): 
                        fit[chr] *= (1+mut.fit)
                else:
                    for chr in range(0,self.command.n[pop]): 
                        fit[chr] += mut.fit
                continue

            if self.command.Z == 0:
                for chr in mut.chrs[pop]:
                    fit[chr] *= (1+mut.fit)
            else:
                for chr in mut.chrs[pop]:
                    fit[chr] += mut.fit

        for chr in range(0, self.command.n[pop]):
            print fit[chr],
        print

        return

    def get_sfs(self, pop=0, NS=True, SYN=True):

        """
        compute the site frequency specutrum for a population

        * Parameters:
          
          * *pop=0*
             population number
        """

        sfs = [0 for i in range(0,self.command.n[pop])]

        for mut in self.muts:
            if pop in mut.chrs and mut.fixed_pop[pop] == 0 and -1 not in mut.chrs[pop]:
                if NS == True and mut.non_or_syn==1:
                    sfs[mut.pops_numchr[pop]-1]+=1
                elif SYN == True and mut.non_or_syn==0:
                    sfs[mut.pops_numchr[pop]-1]+=1
        return sfs

    def haplotype(self,pop=0):
  
        self.haplo[pop] = [ [0 for i in range(0,len(self.muts)+1)] for j in range(0, self.command.n[pop])]
 
        i = 0
        for mut in self.muts:
            if pop not in mut.fixed_pop:
                continue
            if mut.fixed_pop[pop] == True:
                continue
            elif -1 in mut.chrs[pop]:
                continue
            for chr in mut.chrs[pop]:
                self.haplo[pop][chr][i] = 1
            i += 1
     
    def print_hap(self,pop=0):

        if pop not in self.haplo:
            print "Warning: must define the haplotype first!"
            print "use self.haplotype(pop=0)"
            print "exiting the function"
            return

        for arr in self.haplo[pop]:
            for var in arr:
                print var,
            print 

    def sim_pheno(self,div=0.1,thresh=2.,c=0.,effect_small=0.01,effect_big=0.1):
             
        if int(self.command.W[0]) != 2:
            print "only gamma distributed selection coefficients are",
            print "currently supported"
            return

        from scipy.stats import gamma
        from math import floor
        import random

        # for now assuming only negative selection coefficients

        #shape = float(self.command.W[4])
        #scale = 1./float(self.command.W[5])
 
        #p_div = gamma.cdf(div, shape, scale=scale)
    
        thresh /= (2.*self.command.N)

        sel_eff = []
        for mut in self.muts:
            if abs(mut.fit) > thresh:
                sel_eff.append(float(mut.fit))

        med = 0
        if len(sel_eff) > 0:

            sel_eff.sort()
            med =sel_eff[int(floor(len(sel_eff)/2.))]
        else:
            return        

        for mut in self.muts:
            if mut.fit < med: 
                if (random.random() > 0.5*(1-c)):
                    mut.effect = effect_big
                else:
                    mut.effect = effect_small
            elif abs(mut.fit) > thresh:
                if(random.random() > 0.5*(1+c)):
                    mut.effect = effect_big
                else:
                    mut.effect = effect_small

        for mut in self.muts:
            print mut.effect

class Mutation:
    
    """
        a class to store the data associated with a variant in an 
        SFS_CODE output file.  Note, both mutations and substitutions
        are stored as instances of this class.  

        * Attributes:

          * *self.locus=-1*
             The locus number of the variant

          * *self.AXY='?'* 
             'A' for autosomal, 'X' or 'Y' for the corresponding sex
             chromosomes

          * *self.pos=-1*
             The position within the locus.  Note that the positions
             within each locus start from 0.

          * *self.t_init={}* 
             A dictionary, keyed by population number, and storing 
             the time that the variant arose  

          * *self.t_fix={}*
             A dictionary, keyed by population number, and storing 
             the time that the variant fixed within the population.
             If the variant is segregating, the time stored is the
             time of sampling.

          * *self.tri_nuc='NNN'*
             The ancestral trinucleotide (the middle base is the 
             mutated base, so this is not necessarily a codon!)
 
          * *self.deriv_n='N'*
             The derived nucleotide
             
          * *self.non_or_syn='?'*
             Is the mutation synonymous (0) or nonsynonymous (1).
             0 also is used to indicate non-coding.

          * *self.ancest='?'*
             Ancestral amino acid
 
          * *self.deriv_aa='?'*
             Derived amino acid

          * *self.fit='?'*
             fitness effect of the mutation (0 for neutral)

          * *self.chrs = defaultdict(dict)*
             A dictionary of dictionaries that is keyed by population
             and chromosome number.  

             E.g., if the derived allele is present on chromosome 11 in
             population 2, then 
 
             self.chrs[2][11] = True

          * *self.pops_numchr = {}*
 
             A dictionary that stores the number of chromosomes that
             carry the derived allele in each population.

    """
   
    def __init__(self):


        self.locus = -1
        self.AXY = '?'
        self.pos = -1
        self.t_init = {}
        self.t_fix = {}
        self.tri_nuc = 'NNN'
        self.deriv_n = 'N'
        self.non_or_syn = '?'
        self.ancest = '?'
        self.deriv_aa = '?'
        self.fit = '?'
        self.chrs = defaultdict(dict)
        self.fixed_pop = {}
        self.pops_numchr = {}
        self.multiallelic = False
        self.effect = 0.

    # setters
  
    def set_locus(self,locus):
        self.locus=int(locus)
    def set_AXY(self,AXY):
        self.AXY=AXY
    def set_pos(self,pos):
        self.pos=int(pos)    
    def set_t_init(self,t_init,pop):
        self.t_init[pop]=int(t_init)    
    def set_t_fix(self,t_fix,pop):
        self.t_fix[pop]=int(t_fix)
    def set_tri_nuc(self,tri_nuc):
        self.tri_nuc=tri_nuc
    def set_deriv_n(self,deriv_n):
        self.deriv_n=deriv_n
    def set_non_or_syn(self,non_or_syn):
        self.non_or_syn=int(non_or_syn)
    def set_ancest(self,ancest):
        self.ancest=ancest
    def set_deriv_aa(self,deriv_aa):
        self.deriv_aa=deriv_aa
    def set_fit(self,fit):
        self.fit=float(fit)
    def set_num(self,num):
        self.num=float(num)
    def set_chrs(self,chrs):
        for chr in chrs:
           self.chrs[int(chr.split('.')[0])][int(chr.split('.')[1])] = True
     
    # master setter

    def set_all(self,all,P,n):
        self.set_locus(all.pop(0))        
        self.set_AXY(all.pop(0))        
        self.set_pos(all.pop(0))        
        t_init = all.pop(0)        
        t_fix =all.pop(0)        
        self.set_tri_nuc(all.pop(0))        
        self.set_deriv_n(all.pop(0))        
        self.set_non_or_syn(all.pop(0))        
        self.set_ancest(all.pop(0))        
        self.set_deriv_aa(all.pop(0))        
        self.set_fit(all.pop(0))        
        self.set_num(all.pop(0))        
        self.set_chrs(all)           

        for popu in self.chrs:
            self.set_t_init(t_init,int(popu))
            self.set_t_fix(t_fix,int(popu))
            self.pops_numchr[popu] = 0
            if -1 in self.chrs[popu]:
                self.fixed_pop[popu] = True
                self.pops_numchr[popu]=P[0]*n[popu]
                continue
            
            self.fixed_pop[popu] = False
            for chr in self.chrs[popu]:  
                self.pops_numchr[popu]+=1


class SFSData:

    """
    A class that handles the basic parsing of sfs_code output file
    data.

    * Parameters:
  
      * *file= ''*
         the path to the file that is to be read.
    
    * Attributes:

      * *self.file = file*
         the path to the file that is to be read.
          
      * *sims = []*
         an array of sfs.Simulation objects.

    """

    def __init__(self,file=''):
        self.file = file
        self.sims = []

    def set_file(self,file):
        self.file=file

    def get_sims(self):
     
        """
        A method that reads sfs_code output files and stores all the data in
        sfs.Simulation objects.

        """

        sim = Simulation()
       
        iter_num = 0
        line_num = 0
        seq_next = -1

        try: 
            f = open(self.file, 'r')
        except:
            print 'Error: cannot open file', self.file, 'for reading!' 
            exit(-1)

        command_file = command.SFSCommand()

        for line in f:

            if seq_next == 0:
               seq_next-=1
               continue

            if line_num == 0:
                line_num += 1
                sim.command.line=line.rstrip()
                sim.command.com_string = line.rstrip()
                sim.command.parse_string()
                command_file = sim.command
                continue

            if(re.search('iteration',line)):
                if iter_num == 0:
                    iter_num+=1
                    continue
                else:
                    sim.set_command(command_file)
                    sim.make_muts()
                    self.sims.append(sim)
                    sim = Simulation()
                    continue                            
            line = line.strip()
            if re.search('locus',line):
                continue
                seq_next = 0
            if re.search('Nc',line):
                continue
            if re.search('MALES',line):
                continue
            if re.search(',',line):
                sim.data+=line
        
        # now just add the last one

        sim.set_command(command_file)
        sim.make_muts()
        self.sims.append(sim)

        f.close()


class msData:

    # this class is a bit underdeveloped at the moment 
    # and hence is not present in the documentation

    def __init__(self,file=''):

        self.file = file
        self.sims = []

    def get_sims(self):

        try:
            fh = open(self.file,'r')

        except:
            print "Cannot open file ", 
            print self.file,
            exit(-1)

        start = 0
        first = 0
        sim = ms.Simulation()
        for line in fh:
            if re.search('//',line):
                if first != 0:
                    self.sims.append(sim)
                first += 1 
                sim = ms.Simulation()
                start = 0
                continue
            if (start < 2):
                start += 1
                continue
            chr = line.strip()
            sim.chrs.append(chr)
        self.sims.append(sim) 

