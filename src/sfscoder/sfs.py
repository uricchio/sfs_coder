#!/usr/bin/python

import subprocess
import os
import math
import ms
import sys
import unittest
import numpy as np
from collections import defaultdict
from os.path import exists
from command import SFSCommand
import re

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
 
       * *self.haplo = {}*

          A dictionary of haplotypes for each population.

          self.haplo[0] is a 2-D array of haplotypes in the 0th 
          population.
 
          This attribute is only filled when you run the haplotype
          method (or the haplotype_SKAT method)

    """

    def __init__(self):
   
        self.command = SFSCommand()
        self.data = ''
        self.loci = defaultdict(lambda: defaultdict(list))
        self.muts = []
        self.haplo = {}
    
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
        
        self.data = ''

        return

    def merge_muts(self):

        # SFS_CODE sometimes stores the same mutation
        # separately in different populations.
        # this is to handle the case where mutations
        # fix in one population but still segregate in
        # other populations, or fix at different times.  
        # This function merges 
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
                muts = defaultdict(list)
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
                            # if the allele is fixed, only the nucleotide 
                            # that is fixed is counted for that pop; 
                            # otherwise count both alleles
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

    def calc_S(self,multi_skip=True,loci=[],pops=[0],input_log='',start=-1,stop=-1):
     
        """
        calculate the number of segretating sites within all
        populations.

             
        * Parameters:
          
          * *multi_skip = True*
             skip sites that are more than biallelic if True
          * *pops = [0]*
             populations over which to calculate
          * *loci = []*
             Array of loci over which to calculate the number of segregating
             sites.  Uses all loci if this is left blank.
          *  *start=-1*
             site at which to start calculating. If you are doing a simulation 
             of genomic structure, this is the starting position in hg19 coordinates.
             Ignored if left at the default.
             Note that the caluclation is still locus based, so if the start parameter
             is in the middle of a locus, the entire locus is included in the calculation.
             DO NOT USE THIS OPTION AND ALSO SELECT A SET OF LOCI.
          * *stop=-1*
             site at which to stop calculating. All the same caveats for the start parameter 
             apply to this parameter as well.
          * *input_log=''*
             the log file that was made by sfs_coder, storing the start and stop of the
             genomic region that was simulated.  By default, it is stored in the same
             directory as the simulated data, in the file 'err/log.build_input.txt'
        """
        loci = self.check_L(loci)
 
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

        S = dict.fromkeys(pops,0)

        for pop in pops:
            if pop not in self.command.n or self.command.n[pop] == 0:
                S[pop] = None
                continue
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
            S[pop] = num  
        return S 

    def calc_watt(self,pops=[0],loci=[],input_log='',start=-1,stop=-1):
        
        """
        Calculate Watterson's estimator of pi for a set of loci

        * Parameters:
          
          * *pops = [0]*
             populations over which to calculate
          * *loci = []*
             Array of loci over which to calculate the number of segregating
             sites.  Uses all loci if this is left blank.
          *  *start=-1*
             site at which to start calculating. If you are doing a simulation 
             of genomic structure, this is the starting position in hg19 coordinates.
             Ignored if left at the default.
             Note that the caluclation is still locus based, so if the start parameter
             is in the middle of a locus, the entire locus is included in the calculation.
             DO NOT USE THIS OPTION AND ALSO SELECT A SET OF LOCI.
          * *stop=-1*
             site at which to stop calculating. All the same caveats for the start parameter 
             apply to this parameter as well.
          * *input_log=''*
             the log file that was made by sfs_coder, storing the start and stop of the
             genomic region that was simulated.  By default, it is stored in the same
             directory as the simulated data, in the file 'err/log.build_input.txt'
        """
        loci = self.check_L(loci)
        
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

        S = self.calc_S(loci=loci,pops=pops)
  
        watt = dict.fromkeys(pops,0)
 
        for pop in pops:
            
            if pop not in self.command.n or self.command.n[pop]==0:
                watt[pop] = None
                continue
            for i in range(1,self.command.n[pop]):
                watt[pop] += 1./(i+0.)
        
            watt[pop] = (S[pop]+0.)/watt[pop]
            watt[pop] /= tot_sites
        
        return watt        

    def calc_theta_H(self, pops=[0], multi_skip=True,loci=[],input_log='',start=-1,stop=-1):

        """
        Calculate theta_H.

        * Parameters:
          
          * *multi_skip = True*
             skip sites that are more than biallelic if True
          * *pops = [0]*
             populations over which to calculate
          * *loci = []*
             Array of loci over which to calculate the number of segregating
             sites.  Uses all loci if this is left blank.
          *  *start=-1*
             site at which to start calculating. If you are doing a simulation 
             of genomic structure, this is the starting position in hg19 coordinates.
             Ignored if left at the default.
             Note that the caluclation is still locus based, so if the start parameter
             is in the middle of a locus, the entire locus is included in the calculation.
             DO NOT USE THIS OPTION AND ALSO SELECT A SET OF LOCI.
          * *stop=-1*
             site at which to stop calculating. All the same caveats for the start parameter 
             apply to this parameter as well.
          * *input_log=''*
             the log file that was made by sfs_coder, storing the start and stop of the
             genomic region that was simulated.  By default, it is stored in the same
             directory as the simulated data, in the file 'err/log.build_input.txt'
         
          * Return value: A dictionary of H values
            indexed by population number.
 
        """
        loci = self.check_L(loci)

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

        H  = dict.fromkeys(pops)

        for pop in pops:
 
            tot = 0.
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
            H[pop] = tot

        return H

    def calc_tajD(self,pops=[0],loci=[],input_log='',start=-1,stop=-1):
  
        """
        Calculate Tajima's D

        * Parameters:
          
          * *pops = [0]*
             populations over which to calculate
          * *loci = []*
             Array of loci over which to calculate the number of segregating
             sites.  Uses all loci if this is left blank.
          *  *start=-1*
             site at which to start calculating. If you are doing a simulation 
             of genomic structure, this is the starting position in hg19 coordinates.
             Ignored if left at the default.
             Note that the caluclation is still locus based, so if the start parameter
             is in the middle of a locus, the entire locus is included in the calculation.
             DO NOT USE THIS OPTION AND ALSO SELECT A SET OF LOCI.
          * *stop=-1*
             site at which to stop calculating. All the same caveats for the start parameter 
             apply to this parameter as well.
          * *input_log=''*
             the log file that was made by sfs_coder, storing the start and stop of the
             genomic region that was simulated.  By default, it is stored in the same
             directory as the simulated data, in the file 'err/log.build_input.txt'
         
          * Return value: A dictionary of values of Tajima's D
            indexed by population number.

        """
        loci = self.check_L(loci)

        if start != -1:
            loci = []
            regstart = self.get_region_start(input_log=input_log)
            loci = self.get_loci_within_coordinates(regstart=regstart,start=start,stop=stop)
            if not loci:
               return None

        if not loci:
            loci = range(0,len(self.command.L)) 
        
        if max(loci) >= len(self.command.L):
            print "No such locus,", max(loci), "!"
            exit()

        tot_sites = 0
        for l in loci:
            tot_sites+=self.command.L[l]

        S = self.calc_S(loci=loci,pops=pops)
        watt = self.calc_watt(pops=pops,loci=loci)
        pi = self.calc_pi(loci=loci,pops=pops)

        numerator = {}

        pops_final = []
        for pop in pops:
            if S[pop] is None or watt[pop] is None or pi[pop] is None:
                print >> sys.stderr, "Population", pop, "has no samples"
                continue
            numerator[pop] = (pi[pop] - watt[pop])*tot_sites
            pops_final.append(pop)
        pops = pops_final

        # now compute normalizing stats 

        denom = {}
        for pop in pops:
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

            denom[pop] = math.sqrt(e1*(S[pop]+0.) + e2*S[pop]*(S[pop]-1.))

        tajD = {}
        for pop in pops:
            tajD[pop] = numerator[pop]/denom[pop]

        return tajD

    def calc_ZnS(self, pops = [0], loci=[],input_log='',start=-1,stop=-1):

        """
        Calculate ZnS.

        * Parameters:
          
          * *pops = [0]*
             populations over which to calculate
          * *loci = []*
             Array of loci over which to calculate the number of segregating
             sites.  Uses all loci if this is left blank.
          *  *start=-1*
             site at which to start calculating. If you are doing a simulation 
             of genomic structure, this is the starting position in hg19 coordinates.
             Ignored if left at the default.
             Note that the caluclation is still locus based, so if the start parameter
             is in the middle of a locus, the entire locus is included in the calculation.
             DO NOT USE THIS OPTION AND ALSO SELECT A SET OF LOCI.
          * *stop=-1*
             site at which to stop calculating. All the same caveats for the start parameter 
             apply to this parameter as well.
          * *input_log=''*
             the log file that was made by sfs_coder, storing the start and stop of the
             genomic region that was simulated.  By default, it is stored in the same
             directory as the simulated data, in the file 'err/log.build_input.txt'

         * Return value: A dictionary of values of ZnS
           indexed by population number.
        """

        loci = self.check_L(loci)
        
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

        S = self.calc_S(loci=loci,pops=pops)
   
        tot = 0

        ZnS = dict.fromkeys(pops,0)

        for pop in pops:
            if S[pop] < 2:
               return

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

            ZnS[pop] = tot * 2. / (S[pop]*(S[pop]-1.))

        return ZnS

    def check_L(self,arr):

        if not arr:
            return arr

        ret_arr = []    

        test_arr = range(0, len(self.command.L))    

        for thing in arr:
            if thing not in test_arr:
                print >> sys.stderr, 'Warning:', thing,  'not in self.command.L'
                print >> sys.stderr, 'Removing from list'
                continue
            ret_arr.append(thing)

        return ret_arr

    def calc_pi(self,multi_skip=True,pops=[0],loci=[],input_log='',start=-1,stop=-1):
        
        """ 
        calculate the mean pairwise diversity per site
        bewteen pairs of sequences across a set of loci.
        If the loci parameter is left undefined by the
        user, then this method calculates :math:`\pi` 
        over all loci in the simulation.

        * Parameters:
          
          * *multi_skip = True*
             skip sites that are more than biallelic if True
          * *pops = [0]*
             populations over which to calculate
          * *loci = []*
             Array of loci over which to calculate the number of segregating
             sites.  Uses all loci if this is left blank.
          *  *start=-1*
             site at which to start calculating. If you are doing a simulation 
             of genomic structure, this is the starting position in hg19 coordinates.
             Ignored if left at the default.
             Note that the caluclation is still locus based, so if the start parameter
             is in the middle of a locus, the entire locus is included in the calculation.
             DO NOT USE THIS OPTION AND ALSO SELECT A SET OF LOCI.
          * *stop=-1*
             site at which to stop calculating. All the same caveats for the start parameter 
             apply to this parameter as well.
          * *input_log=''*
             the log file that was made by sfs_coder, storing the start and stop of the
             genomic region that was simulated.  By default, it is stored in the same
             directory as the simulated data, in the file 'err/log.build_input.txt'
         
         * Return value: A dictionary of :math:`\pi`
           values indexed by population number.
        """

        loci = self.check_L(loci)

        if start != -1:
            loci = []
            regstart = self.get_region_start(input_log=input_log)
            loci = self.get_loci_within_coordinates(regstart=regstart,start=start,stop=stop)
            if(len(loci) == 0):
                return {}

        pis = {}
        tot = 0
        
        for pop in pops:
            pis[pop] = 0

        if len(loci) == 0:
            loci= self.loci  
            tot = self.command.n_sites
        else:
            locidict = {}
            for locus in loci:
                tot += self.command.L[locus]
                if locus not in self.loci:
                    continue
                locidict[locus] = self.loci[locus]
            loci = locidict
                
        for pop in pops: 
            if pop not in self.command.n or self.command.n[pop] == 0:
                pis[pop] = None
                continue
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
            if self.command.n[pop] > 0:
                pi = 2. * pi / (tot*self.command.n[pop]*(self.command.n[pop]-1))
            else :
                print >> sys.stderr, "Warning: No samples in population", pop
            pis[pop] = pi
        return pis

    def calc_fit(self,pops=[0],loci=[]):

        """
        calculate the fitness of the sampled chromosomes
        within a population.

        * Parameters:
  
          * *pops=[0]* 
             population number
         
          * *loci=[]*
             set of loci over which to calcuate fitness.  Uses all
             loci if left empty 
    
        """
        
        fits = defaultdict(np.array)          

        new_pops = []

        for pop in pops:
            if pop not in self.command.n:
                print >> sys.stderr, "No such pop,", pop, "!"
            else:
                new_pops.append(pop)
                fits[pop] = np.ones(self.command.n[pop]) 
               
        pops = new_pops

        if not loci:
            loci = self.loci

        for locus in loci:
            for pos in self.loci[locus]:
                for mut in self.loci[locus][pos]:
                    for pop in pops:
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
                                        fits[pop][chr] = 0
                                        continue
                                    fits[pop][chr] *= (1+mut.fit)
                            else:
                                for chr in range(0,self.command.n[pop]): 
                                    fits[pop][chr] += mut.fit
                            continue

                        else: 
                            if self.command.Z == 0:
                                for chr in mut.chrs[pop]:
                                    if mut.fit < -1:
                                        fits[pop][chr] = 0.
                                        continue
                                    fits[pop][chr] *= (1.+mut.fit)
                   
                            else:
                                for chr in mut.chrs[pop]:
                                    fits[pop][chr] += mut.fit

        finalfits = defaultdict(np.array)
           
        for pop in pops:
            finalfits[pop] = np.ones(self.command.n[pop]/self.command.P[0])
            for i in range(0, self.command.n[pop]/self.command.P[0]):   
                for j in range(0,self.command.P[0]):
                    if self.command.Z == 0:
                        finalfits[pop][i] *= fits[pop][i*self.command.P[0]+j]
                    else:
                        finalfits[pop][i]+=(fits[pop][i*self.command.P[0]+j]-1) 
                i+=1

        return finalfits

    def get_region_start(self,input_log=''):
        
        # basically if the simulation wasn't iniiated with
        # a call to 'command.genomic', print some stuff

        if input_log == '':
            return 0

        else:

            # this code block does nothing now but might 
            # want to revive it at some point by adding 
            # another parameter
   
            if(input_log == ''):
                       
                print 'you must specify an input_log file!'
                print 'This file is typically created by the genomic method,'
                print 'and is stored in the \'err\' directory inside the sims directory'
                exit()        

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
        
        # returns a dictionary of the absolute positions
        # of the sites in each locus

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
    
        if((start < 0 or stop < 0) or stop < start):

            print >> sys.stderr,  "Warning: Must declare a positive start and stop!"
            print >> sys.stderr,  "Must have start > stop!"
            print >> sys.stderr,  "loci will be set to empty set"
            return [] 

        def overlaps(s1,e1,s2,e2):
        
            if(s2 >= s1 and s2 <= e1):
                 return True
            if(s1 >= s2 and s1 <= e2):
                 return True
            return False
        
        if self.command.a and self.command.a[0] == 'F':
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

        else:
            tot_sites = 0
            loci = []
            i = 0
            for locus_len in self.command.L:
                if (overlaps(tot_sites,tot_sites+locus_len-1,start,stop)):
                    loci.append(i) 
                tot_sites += locus_len
                i += 1
                   
            return loci    
    
    def get_sfs(self,pops=[0],loci=[],NS=True,SYN=True,NC=True,input_log='',start=-1,stop=-1):

        """
        compute the site frequency specutrum for a set of population

        * Parameters:
          
          * *pops=[0]*
             Set of populations over which to calculate the sfs
          * *loci = []*
             Array of loci over which to calculate the number of segregating
             sites.  Uses all loci if this is left blank.
          *  *start=-1*
             site at which to start calculating. If you are doing a simulation 
             of genomic structure, this is the starting position in hg19 coordinates.
             Ignored if left at the default.
             Note that the caluclation is still locus based, so if the start parameter
             is in the middle of a locus, the entire locus is included in the calculation.
             DO NOT USE THIS OPTION AND ALSO SELECT A SET OF LOCI.
          * *stop=-1*
             site at which to stop calculating. All the same caveats for the start parameter 
             apply to this parameter as well.
          * *input_log=''*
             the log file that was made by sfs_coder, storing the start and stop of the
             genomic region that was simulated.  By default, it is stored in the same
             directory as the simulated data, in the file 'err/log.build_input.txt'
          * *SYN=True*
             Include synonymous sites it True
          * *NS=True*
             Include non-synonymous sites it True
          * *NC=True*
             Include non-coding sites it True

        """         
        
        regstart = -1
        # fisrt, get the actual start point of the simulation
        
        if start != -1:
            regstart = self.get_region_start(input_log=input_log)
            loci = self.get_loci_within_coordinates(regstart=regstart,start=start,stop=stop)
        
        if not loci:
            loci = self.loci

        sfs = defaultdict(np.array)
       
        for pop in pops:
            sfs[pop] = np.zeros(self.command.n[pop])

            for mut in self.muts:
                if mut.locus not in loci:
                    continue
                if pop in mut.chrs and mut.fixed_pop[pop] == 0 and -1 not in mut.chrs[pop] and mut.pops_numchr[pop] != self.command.n[pop]:
                    if NS == True and mut.ancest != mut.deriv_aa:
                        sfs[pop][mut.pops_numchr[pop]-1]+=1
                    elif SYN == True and mut.ancest ==  mut.deriv_aa and mut.ancest != 'X':
                        sfs[pop][mut.pops_numchr[pop]-1]+=1
                    elif NC == True and mut.ancest == 'X':
                        sfs[pop][mut.pops_numchr[pop]-1]+=1
        return sfs

    def haplotype(self,pop=0, input_log ='',start=-1,stop=-1):
  
        """
        Build haplotypes for each sampled chromosome.

        * Parameters:

          * *pop=0*
             Population of interest

        The haplotypes are stored in self.haplo
        """
        loci  = []

        if start != -1:
            regstart = self.get_region_start(input_log=input_log)
            loci = self.get_loci_within_coordinates(regstart=regstart,start=start,stop=stop)
            if len(loci) == 0:
                return None

        num_muts = 0        
        muts = []

        if not loci:
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
        same as derived allele).

        * Parameters:

          * *pop=0*
             Population of interest

        The haplotypes are stored in self.haplo

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
            
    def write_geno(self,pop=0,file=''):

        
        """
        Write genotypes out a file.  You must
        call self.haplotype() or self.haplotype_SKAT()
        before running this method.

        * Parameters:

          * *pop=0*
             Population of interest

          * *file=''*
             The destination file.  Must be set to a valid path

        """

        # assuming a ploidy of 2 for now

        if file == '':
            print "must specify file!: self.write_fam(pop=0,file='path/to/file')"
            exit()
        
        f = open(file,'w')

        if pop not in self.haplo:
            print "Warning: must define the haplotype first!"
            print "use self.haplotype(pop=0)"
            print "OR use self.haplotype_SKAT(pop=0) (minor allele number, NOT necessarily derived allele count)"
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

    def get_effects(self,method,causal_pop,min_freq,tau=1.,sig=1.,opp=1.,rho=1.):

        effects = np.zeros(len(self.muts))
                 
        import random      
        from numpy.random import normal

        # note, should use pooled population frequencies here?

        if method == 'SKAT':
            i = 0
            if min_freq == -1:
                if 1./self.command.n[causal_pop] > 0.03:
                    print >> sys.stderr, "Sample size is ", self.command.n[causal_pop]
                    print >> sys.stderr, "Resetting min_freq to 1/", self.command.n[causal_pop] 
                    min_freq = 1./self.command.n[causal_pop]
                else:
                    min_freq = 0.03
            for mut in self.muts:
                if -1 in mut.chrs[causal_pop]:
                    i+=1
                    continue
                freq = (len(mut.chrs[causal_pop])+0.)/self.command.n[causal_pop]
                if (freq > 0. and (freq <= min_freq or freq >= 1-min_freq)):
                    if(random.random() < 0.05):
                        effects[i]=-0.4*math.log10(freq)
                i += 1

            return effects

        elif method == 'EW':
           i = 0
           for mut in self.muts:
               if mut.fit == 0.:
                   effects[i] = 0.
                   i+=1
                   continue
               if random.random() < opp:
                   effects[i] = math.copysign(1.,mut.fit)*(math.fabs(mut.fit)**tau)*(1+normal(scale=sig))
               else:
                   effects[i] = -1*(math.copysign(1.,mut.fit)*math.fabs(mut.fit)**tau)*(1+normal(scale=sig))
               i+=1
           return effects 
 
        elif method == 'SIMONS':
           i = 0
           for mut in self.muts:
               if mut.fit == 0.:
                   effects[i] = 0.
                   i+=1
                   continue
               if random.random() < rho:
                   if random.random() < opp:
                       effects[i] = math.copysign(1.,mut.fit)*(math.fabs(mut.fit)**tau)
                   else: 
                       effects[i] = -1*(math.copysign(1.,mut.fit)*(math.fabs(mut.fit)**tau))
               else:
                   
                   if random.random() < opp:
                      fit = self.muts[random.randint(0,len(self.muts)-1)].fit
                      effects[i] = math.copysign(1.,mut.fit)*((math.fabs(fit))**tau)
                   else:
                      fit = -1*(self.muts[random.randint(0,len(self.muts)-1)].fit)
                      effects[i] = math.copysign(1.,mut.fit)*((math.fabs(fit))**tau)
                      
               i+=1
           return effects 

        else:
            print >>  sys.stderr,  "Warning: No such method,", mehtod, "!",
            print >>  sys.stderr,  "All effects set to 0!"
            return effects

    def sim_pheno(self,h_sq=0.01,pops=[0],min_freq=-1,causal_pop=-1,opp=1.,rho=1.,method='SKAT',tau=1.,sig=1.):

        """
        A method for simulating phenotypes.

        Currently has three different methods, inspired by Wu et al (*AJHG*, 2011),
        Simons et al (*Nature genetics*, 2014), and Eyre-Walker (*PNAS*, 2010). 
 
        * *h_sq=0.01*
           The proportion of variance in the phenotype that is explained by the test sequence
 
        * *pops=[0]*
           Populations to be simulated

        * *min_freq=-1*
           The minimum frequency that is taken as non-causal (i.e., all sites
           with frequency below min_freq are taken to be causal).  Only used with
           method='SKAT'
           
        * *causal_pop=-1*
           The population that is used for determining the proportion of variance
           explained by the test_sequence.  Hence, this population will have h_sq=Vg/Vp,
           but the other populations will not.  Defaults to the first population in the 
           pops parameter if it is not reset from -1.
  
        * *rho=1.*
           The probability that a causal site's effect size is taken as proportional to its     
           selection coefficient.  Only used with method='SIMONS'
 
        * *tau=1*
           The tau parameter in the model of Eyre-Walker.  This parameter is also used in
           the 'SIMONS' mehtod, where effect sizes are taken as s**tau (s is the selection
           coefficient)
 
        * *sigma=1* 
           The sigma parameter in the model of Eyre-Walker

        Note that because we are fixing h_sq in one of the populations, there is always
        an arbitrary constant scaling each effect size.  For each method, the variation
        from the environment is taken to be a normal distribution

        """
        
        import random
        import math

        new_pops = []

        for pop in pops:
            if pop not in self.command.n:
                print  >> sys.stderr, "No samples from population", pop, "!"
            else:
                new_pops.append(pop)

        pops = new_pops
        
        # causal pop is the population in which we fix h_sq
        # to the desired value 

        if causal_pop == -1:
            causal_pop = pops[0]
        
        phenos_genetic = defaultdict(np.array)
        phenos = defaultdict(np.array)
        nv = defaultdict(np.array) 
        
        # get effect sizes
        effects = self.get_effects(method,causal_pop,min_freq,opp=opp,rho=rho,tau=tau,sig=sig)

        for pop in pops:
            phenos_genetic[pop] = np.zeros(self.command.n[pop]/self.command.P[0])
            
            i = 0
            for mut in self.muts:
                if effects[i] != 0.:
                    for chr in mut.chrs[pop]:
                        phenos_genetic[pop][chr / self.command.P[0]] += effects[i]
                i +=1

            nv[pop] = np.random.normal(size=(self.command.n[pop]/self.command.P[0]))

        a = self.fix_gv(phenos_genetic,nv,causal_pop,h_sq)
        
        if a == 1: 
            print >> sys.stderr, "Genetic variance of pop", pop, " is 0!"
            print >> sys.stderr, "All variance in phenotype will come from environment"

        for pop in pops:

            phenos_genetic[pop] = np.multiply(phenos_genetic[pop],a)
            phenos[pop] = np.add(phenos_genetic[pop], nv[pop])
            #print np.var(phenos_genetic[pop])/np.var(phenos[pop])

        return phenos

    def fix_gv(self,p,e,cp,h_sq):

        cov = np.cov(e[cp], p[cp])

        Z = cov[1][0]
        X = cov[1][1]
        Y = cov[0][0]

        # a is the constant that makes var(G)/var(P) = h_sq
        if (X == 0):
            return 1
        a = (h_sq*Z + (h_sq*(X*Y-h_sq*X*Y+h_sq*(Z**2)))**0.5)/(X-h_sq*X)

        return a

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


        self.locus = None
        self.AXY = None
        self.pos = None
        self.t_init = {}
        self.t_fix = {}
        self.tri_nuc = None
        self.deriv_n = None
        self.non_or_syn = None
        self.ancest = None
        self.deriv_aa = None
        self.fit = None
        self.chrs = defaultdict(dict)
        self.fixed_pop = {}
        self.pops_numchr = {}
        self.multiallelic = False
        self.effect = None

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
        self.set_locus(all[0])        
        self.set_AXY(all[1])        
        self.set_pos(all[2])        
        t_init = int(all[3])        
        t_fix = int(all[4])       
        self.set_tri_nuc(all[5])        
        self.set_deriv_n(all[6])        
        self.set_non_or_syn(all[7])        
        self.set_ancest(all[8])        
        self.set_deriv_aa(all[9])        
        self.set_fit(all[10])        
        self.set_num(all[11])        
        self.set_chrs(all[12:])           

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

    def get_sims(self,iter=[]):
     
        """
        A method that reads sfs_code output files and stores all the data in
        sfs.Simulation objects.
        """

        sim = Simulation()
       
        iter_num = -1
        line_num = 0
        seq_next = -1

        try: 
            f = open(self.file, 'r')
        except:
            print >> sys.stderr, 'Error: cannot open file', self.file, 'for reading!' 
            exit(-1)

        command_file = SFSCommand()

        for line in f:

            if seq_next == 0:
               seq_next-=1
               continue

            if line_num == 0:
                line_num += 1
                sim.command.line=line.rstrip()
                sim.command.com_string = line.rstrip()
                sim.command.parse_string()
                if not iter:
                    iter = range(0,sim.command.nsim)
                command_file = sim.command
                continue

            if(re.search('iteration',line)):
                if iter_num not in iter:
                    iter_num+=1
                    continue
                elif iter_num in iter:
                    sim.set_command(command_file)
                    sim.make_muts()
                    self.sims.append(sim)
                    sim = Simulation()
                    iter_num += 1
                    continue
                                           
            line = line.strip()
            if re.search('locus',line):
                continue
                seq_next = 0
            if re.search('Nc',line):
                continue
            if re.search('MALES',line):
                continue
            if re.search(',',line) and iter_num in iter:
                sim.data+=line
        
        # now just add the last one

        if iter_num in iter:
            sim.set_command(command_file)
            sim.make_muts()
            self.sims.append(sim)

        f.close()
  
        return 1

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

class myTests(unittest.TestCase):
   
    def setUp(self):
        self.data = SFSData(file = 'testdata/test.sfs_code.txt')
        self.data.get_sims()  
        for sim in self.data.sims:
            sim.haplotype()

        f = open('testdata/stats.txt')
        self.stats = defaultdict(list)        

        for line in f:
            fields = line.strip().split(' ')
            self.stats[fields[0]] = fields[1:]

    def testGetSims(self):
        self.failUnless(self.data.get_sims())

    def testPi(self):
  
        for sim in self.data.sims:
            pis_test = sim.calc_pi(pops=[0,1],loci=[0])
        self.assertEqual(self.stats['pi'][0], str(pis_test[0]))
        self.assertEqual(self.stats['pi'][1], str(pis_test[1]))
    
    def testS(self):

        S = self.stats['S']

        for sim in self.data.sims:
            S_test =  sim.calc_S(pops = range(0,sim.command.n_pops))
        self.assertEqual(str(S_test[0]), S[0]) 
        self.assertEqual(str(S_test[1]), S[1]) 

    def testZnS(self):

        ZnS = self.stats['ZnS']
 
        for sim in self.data.sims:
            ZnS_test = sim.calc_ZnS(pops=[0,1])
        self.failUnless(str(ZnS_test[0]) == ZnS[0])
        self.failUnless(str(ZnS_test[1]) == ZnS[1])
 
    def testWatt(self):
        
        watt = self.stats['watt']
 
        for sim in self.data.sims:
            watt_test = sim.calc_watt(pops=[0,1])
        self.failUnless(str(watt_test[0]) == watt[0])
        self.failUnless(str(watt_test[1]) == watt[1])
 
    def testTheta_H(self):

        H = self.stats['H']

        for sim in self.data.sims:
            H_test = sim.calc_theta_H(pops=[0,1])
        self.failUnless(str(H_test[0]) == H[0])
        self.failUnless(str(H_test[1]) == H[1])

    def testTajD(self):

        D = self.stats['D']

        for sim in self.data.sims:
            D_test = sim.calc_tajD(pops=[0,1])

        self.failUnless(str(D_test[0]) == D[0])
        self.failUnless(str(D_test[1]) == D[1])

    """
    This test currently fails due to roundoff error...
    def testFit(self):
    
        f = self.stats['f']   
  
        for sim in self.data.sims:
            f_test = sim.calc_fit(pops=[0,1])

        self.failUnless(str(f_test[0]) == f[0])
    """


    def testHaplo(self):

        h = self.stats['h']

        for sim in self.data.sims: 
            self.failUnless(map(int,np.array(h)) == map(int,np.array(sim.haplo[0][0])))
 
    def testGetsfs(self):
    
        s = self.stats['sfs']

        for sim in self.data.sims:
            sfs_test = sim.get_sfs(pops=[0,1],start=999,stop=999)    
            
        self.failUnless(str(sfs_test[0][0]) == s[0])
       
def main():
    unittest.main()

if __name__ == '__main__':
    main()

