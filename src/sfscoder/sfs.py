#!/usr/bin/python

import subprocess
import os
import math
import command
import ms
import re
import sys
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
            newmut.set_all(mut_data,self.command.n)

            self.loci[int(newmut.locus)][int(newmut.pos)].append(newmut)

        self.merge_muts()

        return

    def merge_muts(self):

        # SFS_CODE sometimes stores mutations
        # separately in different populations.
        # this is to handle the case whene mutations
        # fix in one population but still segregate in
        # other populations.  This function merges 
        # these mutations such that they are internally
        # stored together

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
                    # are pushed into a single list, muts[t_init]
                    muts[mut.t_init[pop]].append(mut)
                
                # check whether site is multiallelic

                alleles = {}
                # assuming here that at most one mutation occurs at
                # each position in each generation
                for t_init in muts:
                    for mut in muts[t_init]:
                        for pop in mut.chrs:
                            if -1 in mut.chrs[pop]:
                                alleles[mut.deriv_n] = 0
                            else:
                                alleles[mut.deriv_n] = 0
                                alleles[mut.tri_nuc[1]] = 0
                       
                
                multi = False
                if len(alleles) > 2:
                    multi = True
                
                for t in muts:
                    newmut =  merge(muts[t],multi)
                    self.muts.append(newmut)
                    newloci[locus][pos].append(newmut)
        self.loci = newloci  
           
        #for mut in self.muts:
        #    print mut.multiallelic
        
    def calc_S(self,multi_skip=True,loci=[],pop=0,input_log='',start=-1,stop=-1):
     
        """
        calculate the number of segretating sites within all
        populations.

        * Parameters:
          
          * *multi_skip = True*
             skip sites that are more than biallelic if True
          * *loci = []*
             Array of loci over which to calculate the number of segregating
             sites.  Uses all loci if this is left blank.
        """
 
        if start != -1:
            loci = []
            regstart = self.get_region_start(input_log=input_log)
            loci = self.get_loci_within_coordinates(regstart=regstart,start=start,stop=stop)
            if len(loci) == 0:
                return None

        if len(loci) == 0:
            loci = range(0,len(self.command.L))

        if max(loci) >= len(self.command.L):
                print "No such locus, ", max(loci), "!"
                exit()

        nums = []

        for pop in range(0,self.command.n_pops):
            num = 0
            for locus in loci:
                if locus not in self.loci:
                    continue
                for pos in self.loci[locus]:
                     for mut in self.loci[locus][pos]:
                        if pop not in mut.chrs:
                            continue
                        if float(mut.t_fix[pop]) < self.command.TE:
                            continue
                        if multi_skip:
                            if mut.multiallelic:
                                continue
                        if -1 in mut.chrs[pop]:
                            continue
                        num += 1
            nums.append(num)    
        return nums[pop]  

    def calc_watt(self,pop=0,loci=[],input_log='',start=-1,stop=-1):
        
        """
        Calculate Watterson's estimator of pi for a set of loci

        * Parameters:
    
          * *pop=0*
             Population of interest

          * *loci=[]*
             Array of loci to use in the calculation. If empty,
             uses all loci.
        """
        
        if start != -1:
            loci = []
            regstart = self.get_region_start(input_log=input_log)
            loci = self.get_loci_within_coordinates(regstart=regstart,start=start,stop=stop)
            if len(loci == 0):
                return None

        if len(loci) == 0:
            loci = range(0,len(self.command.L))
        
        if max(loci) >= len(self.command.L):
            print "No such locus,", max(loci), "!"
            exit()

        tot_sites = 0
        for index in loci:
            tot_sites+=self.command.L[index]

        S = self.calc_S(loci=loci,pop=pop)
  
        watt = 0
 
        for i in range(1,self.command.n[pop]):
            watt += 1./(i+0.)
        
        watt = (S+0.)/watt
        watt /= tot_sites
        
        return watt        

    def calc_theta_H(self, pop=0, multi_skip=False,loci=[],input_log='',start=-1,stop=-1):

        """
        Calculate theta_H.

        * Parameters:

          * *pop=0*
             Population of interest

          * *loci=[]*
             Array of loci to use in the calculation.  If empty,
             calculates over all loci.

          * *mutliskip=False*
             Skip multiallelic sites
 
        """

        if start != -1:
            loci = []
            regstart = self.get_region_start(input_log=input_log)
            loci = self.get_loci_within_coordinates(regstart=regstart,start=start,stop=stop)
            if len(loci) == 0:
                return None

        if len(loci) == 0:
            loci = range(0, len(self.command.L))

        if max(loci) >= len(self.command.L):
            print "No such locus,", max(loci), "!"

              
        tot_sites = 0
        for index in loci:
            tot_sites+=self.command.L[index]

        tot = 0
        for locus in loci:
            if locus not in self.loci:
                continue
            for pos in self.loci[locus]:
                for mut in self.loci[locus][pos]: 
                    if pop not in mut.chrs:
                       continue   
                    if (mut.multiallelic == True and multi_skip==True):
                        continue
                    if -1 in mut.chrs[pop]:
                        continue
                    tot += (mut.pops_numchr[pop])**2            
 
        tot /= (tot_sites*self.command.n[pop]*(self.command.n[pop]-1.)+0.)
 
        tot *= 2

        return tot

    def calc_tajD(self,pop=0,loci=[],input_log='',start=-1,stop=-1):
  
        """
        Calculate Tajima's D

        * Parameters:

          * *pop=0*
             Population of interest

          * *loci=[]*
             Array of loci to use in the calculation.  If empty,
             calculates over all loci.

        """

        if start != -1:
            loci = []
            regstart = self.get_region_start(input_log=input_log)
            loci = self.get_loci_within_coordinates(regstart=regstart,start=start,stop=stop)
            if(len(loci) == 0):
               return None

        if len(loci) == 0:
            loci = range(0,len(self.command.L)) 
        
        if max(loci) >= len(self.command.L):
            print "No such locus,", max(loci), "!"
            exit()

        tot_sites = 0
        for l in loci:
            tot_sites+=self.command.L[l]

        S = self.calc_S(loci=loci,pop=pop)
        watt = self.calc_watt(pop=pop,loci=loci)*tot_sites
        pi = self.calc_pi(loci=loci)[pop]*tot_sites

        if (S == 0):
            return None

        numerator = pi - watt

        # now compute normalizing stats 

        a1 = 0
        a2 = 0

        nsam = self.command.n[pop]

        for i in range(1,nsam):
            a1 += 1./(i+0.)
            a2 += 1./(i**2 + 0.)
        b1 = (nsam+1.)/(3.*(nsam-1))

        b2 = (2.*(nsam**2 + nsam + 3.))/(9.*nsam*(nsam-1.))
   
        c1 = b1 - 1./a1

        c2 = b2 - (nsam+2.)/(a1*nsam+0.) + a2/(a1**2)

        e1 = c1/a1
  
        e2 = c2 /(a1**2 + a2 + 0.) 

        denom = math.sqrt(e1*(S+0.) + e2*S*(S-1.))

        return numerator/denom

    def calc_ZnS(self, pop = 0, loci=[],input_log='',start=-1,stop=-1):

        """
        Calculate ZnS.

        * Parameters:

          * *pop=0*
             Population of interest

          * *loci=[]*
             Array of loci to use in the calculation.  If empty,
             calculates over all loci.

        """

        if start != -1:
            loci = []
            regstart = self.get_region_start(input_log=input_log)
            loci = self.get_loci_within_coordinates(regstart=regstart,start=start,stop=stop)
            if len(loci) == 0:
                return None

        if len(loci) == 0:
            loci = range(0,len(self.command.L))
        
        if max(loci) >= len(self.command.L):
            print "No such locus, ", max(loci), "!"
            exit() 

        S = self.calc_S(loci=loci,pop=pop)
   
        if S < 2:
            return

   
        tot = 0

        for i in range(0, len(self.muts)):
            for j in range(i+1, len(self.muts)):
                mut = self.muts[i]
                mut2 = self.muts[j]
                if mut.locus not in loci:
                    continue
                if mut2.locus not in loci:
                    continue
                delij = mut.calc_delij(mut2,pop,n=self.command.n[pop])
                if delij is not None:
                    tot += delij

        ZnS = tot * 2. / (S*(S-1.))

        return ZnS
 
    def calc_pi(self,multi_skip=True,loci=[],input_log='',start=-1,stop=-1):
        
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

        if start != -1:
            loci = []
            regstart = self.get_region_start(input_log=input_log)
            loci = self.get_loci_within_coordinates(regstart=regstart,start=start,stop=stop)
            if(len(loci) == 0):
                return []

        pis = []
        tot = 0
        
        if len(loci) == 0:
            loci= self.loci  
            tot = self.command.n_sites
        else:
            if max(loci) >= len(self.command.L):
                print "No such locus,", max(loci), "!"
            locidict = {}
            for locus in loci:
                tot += self.command.L[locus]
                if locus not in self.loci:
                    continue
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
                        if -1 in mut.chrs:
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
        
        for locus in self.loci:
            for pos in self.loci[locus]:
                for mut in self.loci[locus][pos]:
                    if pop not in mut.t_fix:
                        continue
                    num_fixed = 0

                    for j in range(0,self.command.n_pops):
                        if (j in mut.t_fix):
                            if (mut.t_fix[j] < self.command.TE):
                                num_fixed +=1
                    if num_fixed == self.command.n_pops:
                        continue
  
                    if mut.fit == 0.:
                        continue

                    if -1 in mut.chrs[pop]:
                        if self.command.Z == 0:
                            for chr in range(0,self.command.n[pop]):
                                if mut.fit < -1:
                                    fit[chr] = 0
                                    continue
                                fit[chr] *= (1+mut.fit)
                        else:
                            for chr in range(0,self.command.n[pop]): 
                                fit[chr] += mut.fit
                        continue

                    else: 
                        if self.command.Z == 0:
                            for chr in mut.chrs[pop]:
                                if mut.fit < -1:
                                    fit[chr] = 0
                                    continue
                                fit[chr] *= (1+mut.fit)
                   
                        else:
                            for chr in mut.chrs[pop]:
                                fit[chr] += mut.fit
        f = 1
        i = 0
 
        finalfits = []
  
        for chr in range(0, self.command.n[pop]):
            if self.command.Z != 0:
                f += fit[chr]
                i += 1
                if i % int(self.command.P[0]) == 0:
                    newfit = f - int(self.command.P[0])
                    if newfit < 0:
                        newfit = 0
                    finalfits.append(newfit)
                    f = 1
            else:
                i+=1               
                f *= fit[chr]
                if i % int(self.command.P[0]) == 0: 
                    finalfits.append(f)
                    f = 1

        return finalfits

    def get_region_start(self,input_log=''):
        
        if(input_log == ''):
            print 'you must specify an input_log file!'
            print 'This file is typically created by the genomic method,'
            print 'and is stored in the \'err\' directory inside the sims directory'
        
        regstart = 0 
        r= open(input_log, 'r')
        next = 0
        
        for line in r:
            if re.search('simulating bases:', line):
                next += 1
                continue
            if next > 0:
                data = line.strip().split(' ')
                regstart = int(data[0])
                break
        return regstart

    def get_genomic_coordinates(self,regstart=-1):
        
        posdict = {}
        tot_sites = regstart-1

        f = open(self.command.a[1],'r')
        locus = 0
        for line in f:
            if tot_sites == regstart-1:
                tot_sites += 1
                continue
            data = line.strip().split(';')
            fields = data[0].split(',')
            if(re.search(',',data[0])):
                posdict[locus] = [tot_sites,tot_sites+int(fields[0])]
                locus+=1
            tot_sites += int(fields[0])
        
        return posdict
    
    def get_loci_within_coordinates(self,regstart=-1,start=-1,stop=-1):
        
        # get the loci that are within a set of
        # genomic cooridinates
     
        def overlaps(s1,e1,s2,e2):
        
            if(s2 >= s1 and s2 <= e1):
                 return True
            if(s1 >= s2 and s1 <= e2):
                 return True
            return False

        locdict = {}
        tot_sites = regstart-1

        f = open(self.command.a[1],'r')
        locus = 0
        for line in f:
            if tot_sites == regstart-1:
                tot_sites += 1
                continue
            data = line.strip().split(';')
            fields = data[0].split(',')
            if(re.search(',',data[0])):
                if overlaps(tot_sites,tot_sites+int(fields[0]),start,stop):
                    locdict[locus] = 0
                locus+=1
            
            tot_sites += int(fields[0])
            if tot_sites > stop:
                break         
        loci = [] 
        for locus in sorted(locdict):
            loci.append(locus)

        return loci
  
    def get_sfs(self, pop=0, NS=True, SYN=True, NC=True,input_log='',start=-1,stop=-1):

        """
        compute the site frequency specutrum for a population

        * Parameters:
          
          * *pop=0*
             population number
        """         
        
        import re
        
        regstart = -1
        posdict = {}

        # fisrt, get the actual start point of the simulation
        
        if start != -1:
            regstart = self.get_region_start(input_log=input_log)
            posdict = self.get_genomic_coordinates(regstart=regstart)
        
        sfs = [0 for i in range(0,self.command.n[pop]-1)]

        for mut in self.muts:
            if start !=  -1:
                if not (mut.pos+posdict[mut.locus][0] >= start and mut.pos+posdict[mut.locus][0] <= stop):
                    continue
            if pop in mut.chrs and mut.fixed_pop[pop] == 0 and -1 not in mut.chrs[pop] and mut.pops_numchr[pop] != self.command.n[pop]:
                if NS == True and mut.ancest != mut.deriv_aa:
                    sfs[mut.pops_numchr[pop]-1]+=1
                elif SYN == True and mut.ancest ==  mut.deriv_aa and mut.ancest != 'X':
                    sfs[mut.pops_numchr[pop]-1]+=1
        return sfs

    def haplotype(self,pop=0, input_log ='',start=-1,stop=-1):
  
        """
        Build haplotypes for each sampled chromosome.

        * Parameters:

          * *pop=0*
             Population of interest

        """
        loci  = []

        if start != -1:
            regstart = self.get_region_start(input_log=input_log)
            loci = self.get_loci_within_coordinates(regstart=regstart,start=start,stop=stop)
            if len(loci) == 0:
                return None

        num_muts = 0        
        muts = []

        if(len(loci) == 0):
            num_muts = len(self.muts)
            muts = self.muts
        else:
            for mut in self.muts:
                if mut.locus in loci:
                    muts.append(mut)
                    
            num_muts = len(muts)
                
        self.haplo[pop] = [ [0 for i in range(0,num_muts)] for j in range(0, self.command.n[pop])]
 
        i = 0
        for mut in muts:
            if pop not in mut.fixed_pop:
                i+=1
                continue
            elif mut.fixed_pop[pop] == True:
                for chr in range(0,self.command.n[pop]):
                    self.haplo[pop][chr][i] = 1
                i+=1
                continue
            elif -1 in mut.chrs[pop]:
                for chr in range(0,self.command.n[pop]):
                    self.haplo[pop][chr][i] = 1
                i +=1 
                continue
            for chr in mut.chrs[pop]:
                self.haplo[pop][chr][i] = 1
            i += 1
    
    def haplotype_SKAT(self,pop=0,input_log ='',start=-1,stop=-1):
  
        """
        Build haplotypes for each sampled chromosome in the 
        SKAT fashion (number of minor alleles, not necessarily
        number of derived alleles).

        * Parameters:

          * *pop=0*
             Population of interest

        """
        loci  = []

        if start != -1:
            regstart = self.get_region_start(input_log=input_log)
            loci = self.get_loci_within_coordinates(regstart=regstart,start=start,stop=stop)
            if len(loci) == 0:
                return None

        num_muts = 0
        muts = []

        if(len(loci) == 0):
            num_muts = len(self.muts)
            muts = self.muts
        else:
            for mut in self.muts:
                if mut.locus in loci:
                    muts.append(mut)
                    
            num_muts = len(muts)

        self.haplo[pop] = [[0 for i in range(0,num_muts)] for j in range(0, self.command.n[pop])]
 
        i = 0
        for mut in muts:
            tot = len(mut.chrs[pop])
            if pop not in mut.fixed_pop:
                i+=1
                continue
            # SKAT wants minor alleles, so that's why we use 0 here
            elif mut.fixed_pop[pop] == True:
                i+=1
                continue
            # same reasoning here
            elif -1 in mut.chrs[pop]:
                i +=1 
                continue
            #and same here
            elif tot < self.command.n[pop]/2:
                for chr in mut.chrs[pop]:
                    self.haplo[pop][chr][i] = 1
                i+=1 
                continue
            else:
                for chr in range(0,self.command.n[pop]):
                    if chr not in mut.chrs[pop]:
                        self.haplo[pop][chr][i] = 1        
            i += 1
            
    def write_fam(self,pop=0,file=''):

        if file == '':
            print "must specify file!: self.write_fam(pop=0,file='path/to/file')"
            exit()
        
        f = open(file,'w')

        if pop not in self.haplo:
            print "Warning: must define the haplotype first!"
            print "use self.haplotype(pop=0)"
            print "OR use self.haplotype_SKAT(pop=0) (minor alleles treated in SKAT fashion)"
            print "exiting the function"
            return

        i = 0
        while i < len(self.haplo[pop]):
            #f.write(str(i /2))
            #f.write('\t1\t')
            for k in range(0, len(self.haplo[pop][0])):
                var = (self.haplo[pop][i][k] + self.haplo[pop][i+1][k])
                f.write(str(var)+' ')
            f.write('\n')
            i+=2        

        f.close()

    def sim_pheno_simons(self,thresh=2.,c=0.,effect_small=0.01,effect_big=0.1):
             
        from scipy.stats import gamma
        from math import floor
        import random

        # for now assuming only negative selection coefficients
    
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

    def sim_pheno_EW_alt(self,pop=0,rho=1.,var_prop=1):

        # like the EW model, but allowing us to preselect rho
        # and also having same variance on every variant 
        # (bigfer variants do not have noisier reltationship
        # with effect size)

        from scipy.stats import norm,pearsonr
        from numpy import mean,var

        effects = []
        fits = []

        for mut in self.muts:
            fits.append(mut.fit)

        mean_fit = mean(fits)
        sd_fits = var(fits)**0.5
        
        sd_fit_effect = (sd_fits)*(1./(rho**2)-1)**0.5
  

        for mut in self.muts:
            mut.effect = mut.fit + norm.rvs(scale=sd_fit_effect)
            effects.append(mut.effect)

        if pop not in self.haplo:
            self.haplotype(pop)

        num = 0
        phenos = [0 for i in range(0, self.command.n[pop]/2)]
        pheno =0

        for chr in self.haplo[pop]:
            if(len(chr) != len(effects)):
                print >> sys.stderr, "Error: not every site has an effect size"
                exit()
            for j in range(0, len(chr)):
                pheno += chr[j]*effects[j]
            if num % 2 != 0:
                phenos[num/2-1] = pheno
                pheno = 0
            num += 1

        gen_var = var(phenos)

        if gen_var > 0.:
            env_var = gen_var*(1-var_prop)/(var_prop)
        else:
            env_var = 1.

        for i in range(0, len(phenos)):
            phenos[i] += norm.rvs(scale=env_var**0.5)

        for pheno in phenos:
            print pheno,
        print
    
    def sim_pheno_EW(self,pop=0,sigma=1.,tau=1.,var_prop=0.5):

        # like the EW model, but allowing us to preselect rho
        # and also having same variance on every variant 
        # (bigfer variants do not have noisier reltationship
        # with effect size)

        from scipy.stats import norm,pearsonr
        from numpy import mean,var

        effects = []

        for mut in self.muts:
            effects.append((mut.fit**tau)*(1+norm.rvs(scale=sigma)))

        if pop not in self.haplo:
            self.haplotype(pop)

        num = 0
        phenos = [0 for i in range(0, self.command.n[pop]/2)]
        pheno =0

        for chr in self.haplo[pop]:
            if(len(chr) != len(effects)):
                print >> sys.stderr, "Error: not every site has an effect size"
                exit()
            for j in range(0, len(chr)):
                pheno += chr[j]*effects[j]
            if num % 2 != 0:
                phenos[num/2-1] = pheno
                pheno = 0
            num += 1

        gen_var = var(phenos)

        if gen_var > 0.:
            env_var = gen_var*(1-var_prop)/(var_prop)
        else:
            env_var = 1.

        for i in range(0, len(phenos)):
            phenos[i] += norm.rvs(scale=env_var**0.5)

        for pheno in phenos:
            print pheno,
        print

    def sim_pheno(self,c=1.,spread=0.01,pop=0,var_prop=0.01):

        from scipy.stats import norm,gamma
        from numpy import mean, var

        if var_prop <0. or var_prop >1.:
            print "var_prop is the proportion of variance in the",
            print "phenotype due to sequences being simulated;",
            print "must be in [0,1]"
            exit()

        if 'W' not in self.command.__dict__ or len(self.command.W) == 0:
            sys.stderr.write("Is this simulation not under selection?")
            sys.stderr.write("setting var_prop to 0!\n")
            var_prop = 0.
            mean = 0.

        else :
            a = float(self.command.W[4])
            scale = 1./float(self.command.W[5])
            mean = gamma.mean(a,scale=scale)/(2.*self.command.N)

        effects = []

        for mut in self.muts:
            mut.effect = mut.fit*c + norm.rvs(scale=spread*(c*mean))
            effects.append(mut.effect)

        if pop not in self.haplo:
            self.haplotype(pop)  

        num = 0 
        phenos = [0 for i in range(0, self.command.n[pop]/2)]
        pheno =0
        
        for chr in self.haplo[pop]:
            if(len(chr) != len(effects)):
                print >> sys.stderr, "Error: not every site has an effect size"
                exit()
            for j in range(0, len(chr)):
                pheno += chr[j]*effects[j]
            if num % 2 != 0:
                phenos[num/2-1] = pheno
                pheno = 0
            num += 1     
        
        gen_var = var(phenos)
         
        if gen_var > 0.:
            env_var = gen_var*(1-var_prop)/(var_prop)
        else:
            env_var = 1.

        for i in range(0, len(phenos)):
            phenos[i] += norm.rvs(scale=env_var**0.5)

        #phen_var = var(phenos)

        for pheno in phenos:
            print pheno,
        print
 
    def print_freq_sel(self,pop=0):
    
        for mut in self.muts:
            
            if pop not in mut.fixed_pop:
                continue
            if mut.fixed_pop[pop]:
                continue
            if -1 in mut.chrs[pop]:
                continue
            if mut.fit == 0.:
                continue
            freq = len(mut.chrs[pop])/(self.command.n[pop]+0.)
            print freq, mut.fit
           
    def get_sel(self,file=file):

        fh = open(file, 'w')

        for mut in self.muts:
            fh.write(str(mut.fit)+'\n')
 
        fh.close()
 
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

    def set_all(self,all,n):
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
                self.pops_numchr[popu]=n[popu]
                continue
            
            self.fixed_pop[popu] = False
            for chr in self.chrs[popu]:  
                self.pops_numchr[popu]+=1

    def calc_delij(self,mut,pop,n):
              
       if pop not in self.chrs or pop not in mut.chrs:
           return None

       pi = (self.pops_numchr[pop]+0.)/n
       pj = (mut.pops_numchr[pop]+0.)/n
       
       if pi == 0 or pj ==0 or pi ==1 or pj == 1:
           return None

       pij = 0 
       
       for chr in self.chrs[pop]:
           if chr in mut.chrs[pop]:
               pij +=1
       pij /= (n+0.)

       Dij = pij -pi*pj

       delij = (Dij**2)/(pi*(1.-pi)*pj*(1.-pj))

       return delij

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
            print >> sys.stderr, 'Error: cannot open file', self.file, 'for reading!' 
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

    def p_fix(self, s,alpha):
        pfix = 0
        if (s >= 0.02):
            pfix = math.exp(-(1+s))
            lim = 0
            while(lim < 2000):
                pfix = math.exp((1+s)*(pfix-1))
                lim +=1
            pfix = 1-pfix
        else:
            pfix = (1-math.exp(-2*s))/(1-math.exp(-2*alpha))
        return pfix

class msData:

    # A class for ms style ourput

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

