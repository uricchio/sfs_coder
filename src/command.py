import math
import os
import sys
import subprocess
from collections import defaultdict

class Command:
    

    """
    This class stores information from parsed command lines
    and provides the mechanics to call SFS_CODE commands.
    """ 

    def add_out(self):
        
        try:
            os.stat(os.path.join(self.outdir,self.prefix))
        except:
            os.makedirs(os.path.join(self.outdir,self.prefix))

        out_path = os.path.join(self.outdir,self.prefix,self.prefix+'.'+
                                str(self.task_id)+'.txt')
      
        self.out_path = out_path

        if isinstance(self, SFSCommand):
        
            self.line.extend(['-o',out_path])

    def execute(self,rand=1): 

        """
        execute a simulation command
        
        * Parameters:
            
            * *rand=1* 
               a random integer. If the value is not reset by the user
               then a new random number is rolled within self.execute.
               This value is used as the random seed for SFS_CODE.
        """

        # this is admittedly implemented in such a way as to 
        # obviate the use of inheriting anything and
        # should be cleaned up to improve its generaility 

        from random import randint

        if rand == 1:
            rand = randint(0,1000000)

        if(isinstance(self,SFSCommand)):
            if self.sfs_code_loc == '':
                print 'sfs_code location not defined!'
                print 'set location with: <SFSCommand>.sfs_code_loc =',
                print '/path/to/sfs_code'
                exit()
 
            else:
                self.line.insert(0,self.sfs_code_loc)
 
            self.add_rand(rand)

            self.add_out()

        out_h = ''

        if(isinstance(self, MSCommand)):
            if self.ms_loc == '':
                print 'ms location not defined!'
                print 'set location with: <MSCommand>.ms_loc =',
                print '/path/to/ms'
                exit()
             
            else:
                self.line.insert(0,self.ms_loc)

            self.add_out()

            try:
                out_h = open(self.out_path,'w')
            except:
                print "Cannot open file ", self.out_path
                exit()


            
        err = os.path.join(self.outdir,self.prefix,self.err)

        try:
            os.stat(err)
        except:
            os.makedirs(err)

        err_file= os.path.join(err,'log.'+str(self.task_id)+'.e')    
            
        err_h = open(err_file, 'w')
        
        if out_h == '':
            out_h = err_h

        for word in self.line:
            print word,
        print

        process = subprocess.Popen(self.line, stdout=out_h, 
                                   stderr=err_h, close_fds = True)
        returncode = process.wait()
    
        # pop off the stuff we added so that the line doesn't grow next time 
        # we run the command

        if(isinstance(self, MSCommand)):
            return returncode
        
        self.line.pop()
        self.line.pop()
        self.line.pop()
        self.line.pop()

        return returncode

class SFSCommand(Command):

    """
    This class is used to store, parse, and convert SFS_CODE command lines.
     
    Upon initialization, an object of the SFSCommand class sets the values
    of many of its attributes to the SFS_CODE defaults.
         
    * Parameters:
 
       * *outdir=os.path.join(os.getcwd(), 'sims')* 
          A directory containing
          subdirectories with sfs_code simulations. 
       * *prefix='out'* 
          The prefix of the out directory and the data files 
          within the out directory.
       * *err='err'*
          The name of the directory that contains all the stderr ouput from
          calling sfs_code.  

    * Attributes:
     
       * *self.com_string=''*
          the entire command stored as a single string.
       * *self.outdir= outdir*   
          the parent directory of output directories for sets of SFS_CODE 
          simulations
       * *self.sfs_code_loc = ''*   
          the location of the SFS_CODE binary.
       * *self.N = 500*  
          the number of individuals in the ancestral population.   
       * *self.P = [2]*  
          the ploidy of the individuals in each population.
       * *self.t = 0.001*
          :math:`\\theta = 4Nu = 0.001`.  This is the value of 
          :math:`\\theta` in the ancestral population
       * *self.L = [5000]* 
          an array containing the length of each simulate locus.
       * *self.B = 5 self.p[0] self.N* 
          the length of the burn in (generations).
       * *self.prefix= prefix* 
          the prefix for the output file directory and each simulation file.
       * *self.r=0.0.* 
          :math:`\\rho = 4Nr = 0.0`.  The value of :math:`\\rho` in the 
          ancestral population.
       * *self.n_pops=1* 
          number of populations.
       * *self.nsim=1* 
          number of simulations.
       * *self.line=[]* 
          an array of strings, each of which is an argument
          to SFS_CODE.  This is the attribute that is used to execute SFS_CODE
          commands.

    """

    def __init__(self,outdir=os.path.join(os.getcwd(), 'sims'),prefix='out',
                 err = 'err'):
       
        """
        
        """

        self.com_string = ''
        self.outdir = outdir
        self.sfs_code_loc = ''
        self.N = 500
        self.P = [2]
        self.f = 0.5
        self.n = {0:12}    # this will break if ploidy is not 2
        self.m = {}        # need to worry about this one (migration)
        self.y = 0.5       # proportion of migrants that are male
        self.t = 0.001       
        self.r = 0.0
        self.M = 0         # substitution model
        self.V = 1
        self.v = 1
        self.L = [5000]   
        self.l = []        # complicated    
        self.a = []        # another complicated one, low priority though
        self.W = []        # and again, low priority
        self.x = 0
        self.i = 0.
        self.B = 5*self.P[0]*self.N
        self.prefix = prefix   
        self.recomb_file=''
        self.line = []
        self.err = err
        if os.environ.has_key('SGE_TASK_ID'):
            self.task_id = os.environ['SGE_TASK_ID']
        else:
            self.task_id = 1
        self.sfs_code_loc = ''
        self.TE = 0
        self.n_pops = 1
        self.nsim = 1
        self.Z = 0
        self.n_sites = 5000.
        self.gap = []

    def parse_string(self):
 
        """
        A method to parse SFS_CODE command lines.  By default, every switch
        is stored as an array with the exception of certain special cases that
        are stored as dictionaries.  

        Note, **only the short form of SFS_CODE options are currently fully 
        supported!**  For example, *-t 0.002* is supported but *--theta 0.002*
        is not. Hence, if you run SFS_CODE with the long forms or wish to 
        analyze code that used the long form to run, you may have issues.
        """

 
        temp= self.com_string.split(" -")

        # get location of sfs_code version called and number of SNPs

        loc = temp[0].split(' ')
        self.sfs_code_loc=loc[0]
        self.n_pops=int(loc[1])
        self.nsim=int(loc[2])
  
        self.attrs = {}

        for i in range(0,len(temp)):

            arr = temp[i].split()
            
            if arr[0] not in self.attrs:
                self.attrs[arr[0]] = 0
                setattr(self, arr[0], arr[1:])
            else:
                curval = getattr(self, arr[0]) 
                if isinstance(curval[0],str):
                    curval = [curval]
                curval.append(arr[1:])
                setattr(self, arr[0], curval)

        # treat the following ones as a special cases because we may often 
        # need these values for calculating stuff

        # We are making self.attrs instead of accessing self.__dict__ because 
        # we want to check whether these have been accessed yet or not

        if self.n_pops > 1 and 'n' not in self.attrs:
            self.n = {}
            for i in range(0,self.n_pops):
                self.n[i] = 12

        if 'r' in self.attrs:
            if len(self.r) == 1:
                self.r = float(self.r[0])

        if 'n' in self.attrs:
            temp = self.n
            if len(temp) == 1 and self.n_pops > 1:
               while(len(temp) < self.n_pops):
                   temp.append(temp[0])
            self.n= {}
            for i in range(0,len(temp)):
                self.n[i] = int(self.P[0])*int(temp[i])

        if 'N' in self.attrs:
            self.N = int(self.N[0])
     
        if 'TE' in self.attrs:
            if 'N' in self.__dict__:
                self.TE = math.floor(float(self.TE[0])*2*self.N)
         
        if 'L' in self.attrs:
            temp = self.L
            self.L = []
            self.n_sites = 0
            self.num_loci = int(temp.pop(0))
            if temp[-1] == 'R':
                temp.pop()
                for i in range(0,self.num_loci/(len(temp))):
                    for thing in temp:
                        self.L.append(float(thing))
                        self.n_sites += float(thing)
            elif len(temp[1:]) != self.num_loci:
                for j in range(0, self.num_loci/(len(temp))):
                    for thing in temp:
                        self.L.append(float(thing))
                        self.n_sites += float(thing)
            else:
                for j in range(0, self.num_loci):
                    self.L.append(float(temp[j]))
                    self.n_sites += float(temp[j])

        if 'Z' in self.attrs:
            self.Z = 1

        if 'a' in self.attrs:
            if self.a[0] == 'F':
                   
                fh = open(self.a[1],'r')
                linenum = 0
                self.L= []
                self.n_sites = 0.
                for line in fh:
                    if linenum == 0:
                        linenum+=1
                        continue
                    else:
                        stuff = line.strip().split(',')
                        if stuff[0][-1] != ';':
                            self.L.append(int(stuff[0]))
                            self.n_sites += float(stuff[0])
                        else:
                            self.gap.append(int(stuff[0][0:len(stuff)-2]))
        return

    def build_RHH(self,alpha=1000.,N0=5000.,rho0=0.001,lam0=(10**-10),
                  delta=0.01, L0 = -1,L1=10.**5,loop_max=10,L_neut=1000.,
                  theta_neut=0.001,minpop=100,recomb_dir='recombfiles',
                  outdir='sims',TE=2,r_within=True, neg_sel_rate = 0., 
                  alpha_neg = 5, additive=1, Lextend=1,mutation=[],bottle=[],
                  expansion=[],nreps=10, Boyko=False,Lmult=False,lmultnum=50,
                  Torg=False,non_coding=True):

        """

        A method to build a recurrent hitchhiking command line using the 
        method of Uricchio & Hernandez (2014, *Genetics*). 

        * Dependencies:
          
          * scipy
          * mpmath

        * Parameters

          * *alpha = 1000* 
             :math:`\\alpha = 2Ns`, the ancestral population 
             scaled strength of selection.  Note that demographic events 
             can change N, and hence they also changle alpha.
          * *N0 = 5000* 
             the ancestral population size
          * *rho0 = 0.001* 
             the population scaled recombination coefficient
             in the ancestral population.
          * *lam0 = 10**-10 
             the rate of positive substitutions per generation
             per site in the population.
          * *delta = 0.01* 
             a single parameter that encapsulates both delta 
             parameters from Uricchio & Hernandez (Genetics, 2014).
             Smaller values of delta result in dynamics that are a better
             match for the original population of size N0, but are more 
             computationally expensive.  We do not recommend using values
             of delta greater than 0.1.  For more information please see the
             paper referenced above.
          * *L0 = -1* 
             the length of the flanking sequence on each side of the
             neutral locus. If L0 is not reset from it's default value, 
             it is automatically set to L0 = s0/r0, where s0 and r0 are 
             alpha/2N0 and rho/4N0, respectively.
          * *theta_neut = 0.001* 
             the neutral value of theta.
          * *TE=2*
             the ending time of the simulation in units of 
             2\*N0\*self.P[0] generations.
          * *r_within=False*
             Currently only works with this option set to False, 
             but in the future will allow for recombination within the 
             neutral locus.
         """      

        from scipy.integrate import quad
        from scipy.optimize import brentq
        from mpmath import gammainc

        recomb_dir = os.path.join(os.getcwd(),recomb_dir)
        outdir = os.path.join(os.getcwd(), outdir)

        self.outdir = outdir

        #internal parameters
        
        r0 = rho0/(4.*N0)
        s0 = (alpha/(2.*N0))

        if L0==-1:
            L0 = s0/(2*r0)

        a =  s0/(r0*L0)

        def integ_diff(s):

            u_max = 2/a
    
            SWL = lambda u: math.fabs((math.exp(-2*(1-math.exp(-s0*u)))*
                                      (1 - ((1-math.exp(-s0*u))/s0)*
                                      (alpha**(-(1-math.exp(-s0*u))/s0))*
                                      gammainc(-(1-math.exp(-s0*u))/s0,
                                               1/(alpha+0.))))-
                                      (math.exp(-2*(1-math.exp(-s*u)))*
                                      (1 - ((1-math.exp(-s*u))/s)*
                                      (alpha**(-(1-math.exp(-s*u))/s))*
                                      gammainc(-(1-math.exp(-s*u))/s,
                                               1/(alpha+0.)))))

            return quad(SWL,0,u_max)[0]

        def integ(s):

            u_max = 2/a
    
            SWL = lambda u: (math.exp(-2*(1-math.exp(-s0*u)))*
                             (1 - ((1-math.exp(-s0*u))/s0)*
                             (alpha**(-(1-math.exp(-s0*u))/s0))*
                             gammainc(-(1-math.exp(-s0*u))/s0,1/(alpha+0.))))
                                       
            return quad(SWL,0,u_max)[0]
        
        init = integ(s0)

        def optim(s):
            return math.fabs((integ_diff(s)/init))-delta

        s = alpha/(2.*minpop)
        
        s = brentq(optim,s0,s)
       
        # now extend the sequence length, after computing s

        a = a/float(Lextend)
        L0 = L0*Lextend

        # compute all remaining perameters

        r = s/(a*L1)
        lam = r*lam0/r0
        N = int(alpha/(2*s))
        rho = 4*N*r

        pfix = 0
        if (s >= 0.1):
            pfix = math.exp(-(1+s))
            lim = 0
            while(lim < 200):
                pfix = math.exp((1+s)*(pfix-1))
                lim +=1
            pfix = 1-pfix
        else:
            pfix = (1-math.exp(-2*s))/(1-math.exp(-2*alpha))

        # nu is the mutation rate of selected alleles
        nu = 2*lam/pfix

        # for negative selection, we are interested in all the sites that 
        # enter, not just those that fix.  The mutation rate is therefore 
        # just scaled by the contraction in sequence length
        nu_neg_unscaled = neg_sel_rate*theta_neut/2
        nu_neg_scaled = nu_neg_unscaled*L0/L1

        p_pos = nu/(nu+nu_neg_scaled)

        nu = nu + nu_neg_scaled

        # need to solve for simulation parameters

        theta_sim = (theta_neut*L_neut + nu*2*L1)/(L_neut + 2*L1)
        theta_sim_ratio = (theta_neut*L_neut)/(nu*L1)

        #paths of interest

        home_dir = os.path.expanduser('~')
        out_path = os.path.join(home_dir, self.outdir)

        try:
            os.stat(out_path)
        except:
            os.makedirs(out_path)

        middist = math.ceil(L1+L_neut)
        finaldist = math.ceil(2*L1+L_neut)
        
        if len(self.recomb_file) == 0:
            self.recomb_file = os.path.join(recomb_dir,'recomb.L0'+str(int(L0))+'.al'+str(int(round(alpha)))+'.lam'+str(lam0)+'.txt')

        try:
            os.stat(recomb_dir)
        except:
            os.makedirs(recomb_dir)

        # proportion of recombination events that fall in 
        # neutral region
        prop_neutral = L_neut/(L_neut+2.*L0)

        first_seg = 0.5-0.5*prop_neutral
        second_seg = 0.5+0.5*prop_neutral

        try:    
            os.stat(self.recomb_file)
        except:
            rfile = open(self.recomb_file, 'w')
            rfile.write('3\n')
            rfile.write(str(int(math.ceil(L1))))
            
            if r_within == True:
                rfile.write(' '+str(first_seg)+'\n')
                rfile.write(str(int(middist)))
                rfile.write(' '+str(second_seg)+'\n')
                rfile.write(str(int(finaldist)))
                rfile.write(' 1.0\n')
            else:
                rfile.write('0.5\n')
                rfile.write(str(int(middist)))
                rfile.write('0.5000000000000001\n')
                rfile.write(str(int(finaldist)))
                rfile.write(' 1.0\n')
  
            rfile.close()

        #self.line = ['/usr/bin/time', '-f', '%e %M', self.sfs_code_loc]
        self.line = []        
  
        # one population, 10 replicates
        self.line.extend(['1',str(nreps)])
  
        # Don't print ancestral sequence
        self.line.append('-A')

        # set theta to the simulated value
        self.line.extend(['-t', str(theta_sim)])

        # Make the simulations additive
        if (additive == 1):
            self.line.append('-Z')


        #set up recombination file
        self.line.extend(['-r','F',self.recomb_file,str(rho)])
        
        #N simulated individuals
        self.line.extend(['-N',str(N)])

        #n smapled individuals at the end of the sim
        self.line.extend(['-n','10'])
        
        #3 Loci of lengths L1, L_neut, L1
        if Lmult == False:
            self.line.extend(['-L','3',str(L1),str(L_neut),str(L1)])
        else:        
            self.line.extend(['-L',str(lmultnum+2),str(L1)])
            for j in range(0,lmultnum):
                self.line.append(str(L_neut/lmultnum))       
            self.line.append(str(L1))

        
        # set noncoding for all sites
        if(non_coding == True):
            self.line.extend(['-a','N','R'])
        
        # set all mutation rates
        self.line.extend(['-v','L','A','1'])

        # set mutation rate at Locus 1 to be different
        if Lmult == False:
             self.line.extend(['-v','L','1',str(theta_sim_ratio)])
        else: 
             for j in range(0,lmultnum):
                 self.line.extend(['-v','L',str(j+1),str(theta_sim_ratio/lmultnum)])
        
        if(neg_sel_rate == 0.):

            # selection on the flanking sequences
            if Lmult == False:
                self.line.extend(['-W','L','0','1',str(alpha),'1','0'])
                self.line.extend(['-W','L','2','1',str(alpha),'1','0'])
            else:
                self.line.extend(['-W','L','0','1',str(alpha),'1','0'])
                self.line.extend(['-W','L',str(lmultnum+1),'1',str(alpha),'1','0'])
                 
  
        elif (Boyko == False and Torg == False):
            
            # selection on the flanking sequences, neutral theta
            l_pos = 1000000000
            a_pos = alpha*l_pos

            l_neg = 1000000000
            a_neg = alpha_neg*l_neg

            if Lmult == False:

                self.line.extend(['-W','L','0','2', str(p_pos), str(a_pos),
                              str(l_pos), str(a_neg), str(l_neg)])
                self.line.extend(['-W','L','2','2', str(p_pos), str(a_pos),
                              str(l_pos), str(a_neg), str(l_neg)])
            else:    
                self.line.extend(['-W','2', 'L','0',str(p_pos), str(a_pos),
                              str(l_pos), str(a_neg), str(l_neg)])
                self.line.extend(['-W','2', 'L',str(lmultnum+1),str(p_pos), str(a_pos),
                              str(l_pos), str(a_neg), str(l_neg)])
     
        elif (Torg == True):

            l_pos = 1000000000
            a_pos = alpha*l_pos

            l_neg = 0.0015625
            a_neg = 0.0415

            if Lmult == False:

                self.line.extend(['-W','L','0','2', str(p_pos), str(a_pos),
                              str(l_pos), str(a_neg), str(l_neg)])
                self.line.extend(['-W','L','2','2', str(p_pos), str(a_pos),
                              str(l_pos), str(a_neg), str(l_neg)])
            else:
                self.line.extend(['-W','L','0','2',str(p_pos), str(a_pos),
                              str(l_pos), str(a_neg), str(l_neg)])
                self.line.extend(['-W','L',str(lmultnum+1),'2',str(p_pos), str(a_pos),
                              str(l_pos), str(a_neg), str(l_neg)])

        elif (Boyko == True):
            l_pos = 1000000000
            a_pos = alpha*l_pos

            l_neg = 0.00040244
            a_neg = 0.184

            if Lmult == False:

                self.line.extend(['-W','L','0','2', str(p_pos), str(a_pos),
                              str(l_pos), str(a_neg), str(l_neg)])
                self.line.extend(['-W','L','2','2', str(p_pos), str(a_pos),
                              str(l_pos), str(a_neg), str(l_neg)])
            else:       
                # selection on all loci
                self.line.extend(['-W','L','0','2',str(p_pos), str(a_pos),
                              str(l_pos), str(a_neg), str(l_neg)])
                self.line.extend(['-W','L',str(lmultnum+1),'2',str(p_pos), str(a_pos),
                              str(l_pos), str(a_neg), str(l_neg)])
                

        # add any demographic events
        for thing in bottle:
            self.line.extend([thing])

        for thing in expansion:
            self.line.extend([thing])

        # add any mutation events
        for thing in mutation:
            self.line.extend([thing])

        # simulation end time
        self.line.extend(['-TE',str(TE)])
 
    def convert_sfs_ms(self):

        from random import randint

        #pseudocode:
        # 1: get population size at each time point
        # for each event, forwards in time
        # 2 : rescale times
        # 3 : build ms commands

        ##########
     
        def get_sizes(s,ti,N,N_b,mig_events,N_m,G_b):

            events = defaultdict(list)
            ordered_events = []
 
            N0 = self.N+0.

            for key,value in s.__dict__.items():
                if not key.startswith('__'):
                    if  key in ti:
                        if isinstance(value, list):
                            if not isinstance(value[0], list):
                                value.append(key)
                                events[key].append(value)
                                ordered_events.append(value)
                            else:
                                for thing in value:    
                                    thing.append(key)
                                    events[key].append(thing)
                                    ordered_events.append(thing)
                        else:
                            if key == 'TE':
                                value/=(self.P[0]*self.N+0.)
                            value = [value]
                            value.append(key)
                            ordered_events.append(value)

            ordered_events = sorted(ordered_events, key=lambda x:(float(x[0]),str(x[-1])))
            t0 = float(ordered_events[0][0])
            N[0][t0] = N0  

            t_dem_b = defaultdict(dict)
            t_cur_pop = {}
            t_dem_list = defaultdict(list)

            for pop in range(0,self.n_pops):
                t_cur_pop[pop] = t0
            
            # get times of previous events for each population
            for event in ordered_events:
                if event[-1] == 'TS':
                    t_dem_list[int(event[1])].append(event)
                    t_dem_list[int(event[2])].append(event)
                    t_dem_b[int(event[1])][float(event[0])] = t_cur_pop[int(event[1])]
                    t_dem_b[int(event[2])][float(event[0])] = float(event[0])
                    t_cur_pop[int(event[2])] = float(event[0])
                    t_cur_pop[int(event[1])] = float(event[0])
                    G_b[int(event[2])][float(event[0])]['TS'] = 0
                if event[-1] == 'Td':
                    if event[1] == 'P':
                        t_dem_list[int(event[2])].append(event)
                        t_dem_b[int(event[2])][float(event[0])] = t_cur_pop[int(event[2])]
                        t_cur_pop[int(event[2])] = float(event[0])
                        G_b[int(event[2])][float(event[0])]['Td']=0
                    else:
                        for pop in range(0,self.n_pops):
                            t_dem_list[pop].append(event)
                            t_dem_b[pop][float(event[0])] = t_cur_pop[pop]
                            t_cur_pop[pop] = float(event[0])
                            G_b[pop][float(event[0])]['Td']=0
                if event[-1] == 'Tg':  
                    if event[1] == 'P':
                        t_dem_list[int(event[2])].append(event)
                        t_dem_b[int(event[2])][float(event[0])] = t_cur_pop[int(event[2])]
                        t_cur_pop[int(event[2])] = float(event[0])
                    else:
                        for pop in range(0,self.n_pops):
                            t_dem_list[pop].append(event)
                            t_dem_b[pop][float(event[0])] = t_cur_pop[pop]
                            t_cur_pop[pop] = float(event[0])
                if event[-1] == 'TE':
                    for pop in range(0,self.n_pops):
                        t_dem_list[pop].append(event)
                        t_dem_b[pop][float(event[0])] = t_cur_pop[pop]
                        t_cur_pop[pop] = float(event[0])
                if event[-1] == 'Tm':
                    mig_events.append(event) 
                    if event[1] == 'L':
                        for pop in range(0,self.n_pops):
                          t_dem_list[pop].append(event)
                    else:
                        t_dem_list[int(event[2])].append(event)

            # get population sizes for every event
            N_cur = {}
            N_cur[0] = self.N+0.
            N_init = defaultdict(dict)

            for pop in t_dem_list:
              for event in t_dem_list[pop]:
                if event[-1] == 'TS':
                    if(int(event[1]) == pop):
                        N_init[int(event[1])][int(event[2])] = N_cur[int(event[1])]
                    N[int(event[2])][float(event[0])] = N_init[int(event[1])][int(event[2])]
                    N_cur[int(event[2])] = N[int(event[2])][float(event[0])]
                if event[-1] == 'Td':
                    if event[1] == 'P':
                        N_b[int(event[2])][float(event[0])] = N_cur[int(event[2])]
                        N[int(event[2])][float(event[0])] =  N_cur[int(event[2])]*float(event[3]) 
                        N_cur[int(event[2])] = N[int(event[2])][float(event[0])]                       
                    else:
                        for pop1 in range(0,self.n_pops):
                            N_b[pop][float(event[0])] = N_cur[pop]
                            N[pop1][float(event[0])] = N_cur[pop1]*float(event[3])
                            N_cur[pop1] = N[pop1][float(event[0])]                       
                if event[-1] == 'Tg':
                    if event[1] == 'P':
                        t = float(event[0])
                        N_m[pop][t] = N_cur[pop]
                        i = 0
                        for thing in t_dem_list[int(event[2])]:
                            if float(thing[0]) == t:
                                while(float(t_dem_list[int(event[2])][i][0]) <= t):
                                    i+=1
                                break
                            i+=1
                        if i < len(t_dem_list[int(event[2])]):
                            tf = float(t_dem_list[int(event[2])][i][0])
                        else:
                            tf = self.TE/(self.P[0]*self.N*1.)
                        t = round(t,len(str(int(2*self.N)))-1)
                        event[0] = tf
                        N[int(event[2])][tf] =  N_cur[int(event[2])]*math.exp(float(event[3])*(tf-t))
                        N_cur[int(event[2])] = N[int(event[2])][tf]
                if event[-1] == 'Tm':
                    if event[1] == 'L' and float(event[0]) not in N_m[pop]:                        
                        N_m[pop][float(event[0])] = N_cur[pop]
                    elif float(event[0]) not in N_m[pop]:
                        N_m[int(event[2])][float(event[0])] = N_cur[int(event[2])]

            tf = (self.TE)/(self.P[0]*self.N+0.)
            for pop in N:
                N[pop][tf] = N_cur[pop]

            return t_dem_list,mig_events,N_m

        ##########

        if self.com_string == '':
            for word in self.line:
                self.com_string+= word+' '
        if self.com_string == '':
            print >> sys.stderr, 'No command string to parse'
            exit(-1)

        self.parse_string()
         
        items = ['N','n','n_sites','t','r','n_pops']
        time_items = ['TS','Td','TE','TJ','Tk','Tg','Tm']
        items = set(items)
        time_items = set(time_items)
        N = defaultdict(dict)
        N_b = defaultdict(dict)
        N_m = defaultdict(dict)
        G_b = defaultdict(lambda: defaultdict(dict))
        mig_events = []    
            
        events,mig_events,N_m =get_sizes(self,time_items,N,N_b,mig_events,N_m,G_b)
        #for pop in N:
        #    for time in N[pop]:
        #        N[pop][time] = int(math.floor(N[pop][time]))
        #        print pop, time, N[pop][time]
         
        #for pop in N_b:
        #    for time in N_b[pop]:
        #        N_b[pop][time] = int(math.floor(N_b[pop][time]))
        #        print pop, time, N_b[pop][time]

        # Now, build ms commands

        tf = self.TE/(2.*self.N)

        ms_com = MSCommand(prefix=self.prefix)
 
        # add MSCommand location
        ms_com.line = [ms_com.ms_loc]

        # nsam and nsim
        nsam = 0
        for pop in self.n:
            nsam += self.n[pop]

        ms_com.line.extend([str(nsam),str(self.nsim)])

        # theta
        theta = self.t*N[0][tf]/self.N
        ms_com.line.extend(['-t',str(self.n_sites*theta)])

        # rho
        rho = self.r*N[0][tf]/self.N
        ms_com.line.extend(['-r',str(self.n_sites*rho),str(self.n_sites)])

        # -I
        ms_com.line.extend(['-I', str(self.n_pops)])
        for pop in self.n:
            ms_com.line.append(str(self.n[pop]))

        #ms_com.line.append('1.0')

        #set initial population sizes
        for pop in N:
            if pop > 0:
                size = N[pop][tf]
                size /= (N[0][tf]+0.)
                ms_com.line.extend(['-n',str(pop+1), str(size)])

        # TO DO
        # don't keep '-en' events that happen when there is a merger; seems to screw things up with migration
        
        for pop in events:
            events[pop].reverse()
            for event in events[pop]:
                #print pop,event
                if event[-1] == 'TS':
                    if (int(event[2]) > pop): 
                        t_e = (tf-float(event[0]))*(self.N)/(2.*N[0][tf])    
                        ms_com.line.extend(['-ej',str(t_e), str(int(event[2])+1), str(int(event[1])+1)]) 

                elif event[-1] == 'Td':
                    if 'TS' in G_b[pop][float(event[0])]:
                        continue
                    if event[1] == 'P':
                        change = N_b[pop][float(event[0])]
                        change /= (N[0][tf]+0.)
                        t_e = (tf-float(event[0]))*(self.N)/(2.*N[0][tf])    
                        ms_com.line.extend(['-en',str(t_e), str(int(event[2])+1), str(change)]) 
                elif event[-1] == 'Tg':
                    if event[1] == 'P':
                        t_e = (tf-float(event[0]))*(self.N)/(N[0][tf]*2.)    
                        rate = float(event[3])
                        rate *=( N[0][tf]*2.)/(self.N)
                        ms_com.line.extend(['-eg',str(t_e), str(int(event[2])+1),str(rate)])                  
        for event in mig_events:           
            if event[1] == 'P':
                t_e = -1
                t_for = float(event[0])
                # search for the next event that changes the migration rate for these two pops
                for event2 in mig_events:
                    if(float(event2[0]) > t_for and
                       ((event2[1] == 'P' and str(event2[2]) == str(event[2]) and str(event2[3]) == str(event[3])) 
                         or (event2[1] == 'L'))): 
                        t_e = float(event2[0])
                        break                
                   
                if t_e == -1:
                    t_e = self.TE/(self.P[0]*self.N+0.)
                t_e = (tf-t_e)*(self.N)/(N[0][tf]*2.)

                N_sfs = N_m[int(event[2])][t_for]
                rate = (float(event[4])/N_sfs)*2*N[0][tf]
                #print (float(event[4])/(2.*N_m[int(event[2])][t_for]))*(500./(7300)) 
                ms_com.line.extend(['-em', str(t_e), str(int(event[2])+1),str(int(event[3])+1),str(rate)])
            elif event[1] == 'L':
                t_e = -1
                t_for = float(event[0])
                # search for the next event that changes the migration rate for any two pops ; this is not guaranteed to work in general!
                # Since the rates are set forwards in time in SFS_CODE, we would need to search forwards in time to get all instances 
                # of all changes and decide which ones make for general changes general solution will be a super pain in the ass to code.
                for event2 in mig_events:
                    if(float(event2[0]) > t_for):
                        t_e = float(event[0])
                if t_e == -1:
                    t_e = self.TE/(self.P[0]*self.N+0.)
                t_e = (tf-t_e)*(self.N)/(N[0][tf]*2.)  
                com_array = []
                i = 0
                j= 0
                counter = -1
                for rate in event[2:(self.n_pops**2-self.n_pops+2)]:
                    if i == j: 
                        j += 1
                    #print N_m[i][t_for]
                    #print i+1, j+1, (float(rate)/(2.*N_m[i][t_for]))*(500./(7300))# => seems I am getting the right underlying rates accoring 
                    # to Gutenkunst et al (2009); so why is SFS so far off in MS?
                    rate = (float(rate)/(N_m[i][t_for]))*2*N[0][tf]
                    com_array.extend(['-em', str(t_e), str(i+1), str(j+1), str(rate)])
                    j += 1
                    if (counter % 2 == 0):
                        i += 1
                        j = 0
                    counter += 1
                ms_com.line.extend(com_array)
                    
        # print ms_com.line
        return ms_com

    def add_rand(self,rand):
 
        self.line.extend(['-s',str(int(math.floor(rand)))])

    def add_event(self, event):
 
        """
        add an event to a command line.
        
        * Parameters:
          
          * *event* 
             an array of strings corresponding to the event.

             e.g., adding a mutation at locus 0 at time 0 with
             
             .. code-block:: python
             
                event = ['--mutation', '0','L','0']

        """
        new_event = []

        for thing in event:
            new_event.append(str(thing))

        self.line.extend(new_event)

    def build_BGS(self, n_sim=10,theta=0.0001,recomb=True,rho=0.001,N=250,
                  nsam=[10],alpha=5,L=10**5,Lmid=10**5):
        
        self.line = [];
        
        self.line.extend(['1',str(n_sim),'-A'])

        self.line.extend(['-t', str(theta)])

        if recomb == True:
            self.line.extend(['-r',str(rho)])

        self.line.extend(['-N', str(N)])

        self.line.append('-n')
        
        for thing in nsam:
             self.line.append(str(thing))

        self.line.extend(['-L', '3', str(L),str(Lmid),str(L)])
        
        self.line.extend(['-a','N'])
        
        self.line.extend(['-W','L','0','1',str(alpha),'0','1'])
        self.line.extend(['-W','L','1','1',str(alpha),'0','1'])
        self.line.extend(['-W','L','2','1',str(alpha),'0','1'])
     
 
    def genomic(self,basedir=os.path.join(os.path.dirname(__file__),'../src/req'),
                      indir=os.path.join(os.getcwd(),'input_files'),
                      datafile='hg19_gencode.v14.gtf.gz',
                      phast_file ='hg19_phastCons_mammal.wig.gz',dense_dist=5000,
                      begpos=134545415,endpos=138594750,db=-1,chr=2,
                      de=-1,withseq=0,fafile='',
                      N=2000,mutation=[],sel=True,model='',t=0.001,rho=0.001,
                      nsim=10,nsam=[20]):

        """

        A method for running simulations of genomic elements using realistic
        genome structure and demographic models.  The demographic models of 
        Gutenkunst (2009, *PLoS Genetics*), Gravle (2011, *PNAS*), Tennessenn
        (2012,*Science*), and the standard neutral model are included.

        The default data sources and options are all
        human-centric, but in principle these methods could be used to simulate
        sequences from any population for which the relevant data sources are 
        available (recombination map, conserved elements, exon positions).   

        This function calls a number of perl scripts, originally 
        implemented by Ryan Hernandez, to build the input to SFS_CODE.  These
        perl scripts are bundled with sfs_coder in the directory src/req

        * Parameters:

          * *basedir=os.path.join(os.path.dirname(__file__),'../src/req')*
             the directory where all the datasources and perl for this method
             are located.  You shouldn't have to change this unless you're
             moving around the source files relative to each other.
          
          * *outdir=os.path.join(os.getcwd(),'input_files')*
             the directory where the sfs_code annotation files will be 
             written

          * *chr=2*
             chromosome number to query

          * *begpos=134545415*
             genomic coordinate of the beginning of the simulated sequence
 
          * *endpos=138594750*
             genomic coordinate of the end of the simulated sequence

          * *db=begpos*
             genomic coordinate of the region in which to include dense 
             neutral sequences.  Neutral sites within this region will be 
             simulated if they are within a specified distance of one
             of the simulated genomic elements (given by *dense_dist*)

          * *de=begpos*
             genomic coordinate of the end of the dense neutral region

          * *dense_dist=5000*
             The amount of neutral sequence to pad onto the end of each 
             conserved element or exon.  For example, if two adjaacent 
             conserved elemets are 20,000 base pairs apart, and the 
             dense_dist=5,000, then 5,000 base pairs are padded onto the 
             end of each of the conserved elements and the middle 10,000
             base pairs are not simulated.
 
             To include no neutral sequence, use dense_dist = 0.
      
             To include every neutral base pair in the region, use 
             dense_dist=-1 (potentially very computationally expensive
             for large regions!) 

          * *sel=True*
             If true, draw selection coefficients from a gamma distribution 
             of selection coefficients for conserved elements and coding
             regions.  These distributions are taken from Boyko et al 
             (coding, 2008, *PLoS Genetics*) and Torgerson et al (conserved 
             non-coding, 2009, *PLoS Genetics*).

        
        """

        if de < 0:
            de = endpos
        if db < 0:
            db = begpos

        files = 0

        if fafile == '':

            fafile = 'chr'+str(chr)+'.fa.gz'

        regname = 'chr'+str(chr)+'_'+str(begpos)+'-'+str(endpos)

        if(sel==True):
            sel = 'SEL'
        else:
            sel='NEUT'

        annotfile = os.path.join(indir,
                            'sfscode_annotation_'+regname+'_CDSbuff_'+str(sel)+'.txt')
         
        recombout = os.path.join(indir, 'sfscode_rec_'+regname+'.txt')

        try:
            os.stat(annotfile)
        except:
            files += 1
 
        try:
            os.stat(recombout)
        except:
            files+=1
 
        if files != 0:
            self.build_genomic(basedir=basedir,indir=indir,datafile=datafile,                                
                               chr=chr,begpos=begpos,endpos=endpos,db=db,de=de,
                               withseq=withseq,fafile=fafile,
                               phast_file=phast_file,dense_dist=dense_dist,
                               regname=regname,annotfile=annotfile,
                               recombout=recombout,sel=sel)
        if model != 'snm':

            self.three_pop(N=N,nsim=nsim,recombfile=recombout,nsam=nsam,
                                t=t,rho=rho,model=model)

        elif model == 'snm':
            self.snm(N=N,nsim=nsim,recombfile=recombout,nsam=nsam,t=t,
                                rho=rho,no_L=True)

        else:
            print 'model', model, 'not recognized!'
            print 'supported models include snm (standard neutral)',
            print 'gutenkunst, gravel, and tennessen'
            exit(-1)

        self.line.extend(['-a','F',annotfile])
  
        self.line.extend(mutation)

        # Boyko et al selection coefficients
        # Apparently this is already included
        # self.line.extend(['-W', '2', '0', '0', '0', '0.184', '0.00040244'])         


    def snm(self,N=500,nsim=10,t=0.001,rho=0.001,nsam=[10],recombfile='',
            mutation=[],sel=[],event=[],TE=0,L=[5000],non_coding=False,
            no_L=False):

        self.line = []

        self.line.extend(['1',str(nsim)])

        self.line.extend(['-t', str(t)])

        self.line.extend(['-N', str(N)])

        self.line.append('-n')

        for thing in nsam:
            self.line.append(str(thing))
        
        if no_L == False:
            self.line.append('-L')
            if str(L[0]) == 'R':
                for thing in L:
                    self.line.append(str(thing))
            else:
                self.line.append(str(len(L)))
                for thing in L:
                    self.line.append(str(thing))

        self.line.extend(['-A'])
 
        if non_coding ==True:
            self.line.extend(['-a', 'N'])
            if int(L[1]) > 1:
                self.line.append('R')


        for thing in mutation:
            self.line.append(str(thing))

        for thing in sel:
            self.line.append(str(thing))
        
        for thing in event:
            self.line.append(str(thing))

        if recombfile =='':
            self.line.extend(['-r',str(rho)])

        else:
            self.line.extend(['-r', 'F', recombfile, str(rho)])

        self.line.extend(['-TE', str(TE)])


    def build_genomic(self,basedir=os.path.join(os.path.dirname(__file__),'../src/req'),
                      indir=os.path.join(os.getcwd(),'input_files'),
                      datafile='hg19_gencode.v14.gtf.gz',chr=2,
                      begpos=134545415,endpos=138594750,db=-1,
                      de=-1,withseq=0,fafile='chr2.fa.gz',
                      phast_file ='hg19_phastCons_mammal.wig.gz',dense_dist=5000,
                      regname='',annotfile='',recombout='',sel='SEL'):
       

        import signal

        if regname=='':
            regname = 'chr'+str(chr)+'_'+str(begpos)+'-'+str(endpos)
     
        datadir=os.path.join(basedir,'hg19')

        if db < 0:
           db = begpos

        if de < 0:
           de = endpos

        # first, generate any directories that need to be made/check exist

        try:
           os.stat(os.path.join(basedir,'ExtractExons.pl'))
        except:
           print 'No such file', 
           print os.path.join(basedir,'ExtractExons.pl')
           exit(-1)

        try: 
            os.stat(indir)
        except:
            os.makedirs(indir)

        try:
            os.stat(datadir)
        except:
            print 'datadir', datadir, 'does not exist!'
            exit(-1)
      
        try:
            os.stat(os.path.join(datadir,datafile))
        except:
            print 'No such file', os.path.join(datadir,datafile), '!'
            exit(-1)
        
        err = os.path.join(self.outdir,self.prefix,self.err)

        try:
            os.stat(err)
        except:
            os.makedirs(err)
        
        log_exons = open(os.path.join(err,'log.exons'),'w')
        eh_exons = open(os.path.join(err,'err.exons'),'w')

        # now build first perl command, for generating the exons regions
  
        exons = ['perl', os.path.join(basedir,'ExtractExons.pl'), 
                 os.path.join(datadir,datafile), 
                 os.path.join(indir,regname+'.UTRs'), 
                 os.path.join(indir,regname+'.exons'), str(chr), str(begpos), 
                 str(endpos)] 

        ex_exons = subprocess.Popen(exons, stdout=log_exons, stderr=eh_exons, 
                                    preexec_fn=lambda: signal.signal(
                                        signal.SIGPIPE, signal.SIG_DFL), 
                                    close_fds =False) 
        returncode = ex_exons.wait()

        eh_exons.close()
        log_exons.close()
        
        eh_ncs = open(os.path.join(err,'log.cncs'),'w')

        # next perl command, conserved non-coding

        ncs = ['perl', os.path.join(basedir,'ExtractConservedNCs.pl'), 
                os.path.join(datadir,phast_file), 
                os.path.join(indir,regname+'.exons'), 
                os.path.join(indir,regname+'.UTRs'), 
                os.path.join(indir,regname+'.CNCs'), str(chr), str(begpos), 
                str(endpos)]


        ex_ncs = subprocess.Popen(ncs, stdout=eh_ncs, stderr=subprocess.STDOUT, 
                                  preexec_fn=lambda: signal.signal(
                                        signal.SIGPIPE, signal.SIG_DFL),
                                  close_fds = True)  
 
        returncode = ex_ncs.wait()

        eh_ncs.close()

        # next perl command, make annotation file for sfs_code

        if annotfile == '': 
            annotfile =os.path.join(indir, 'sfscode_annotation_'+regname+'_CDSbuff_'+str(sel)+'.txt')

        annot = ['perl', os.path.join(basedir,'makeSFSCODEannotationFile.pl'), 
                 str(begpos), str(endpos), str(dense_dist), 
                 os.path.join(annotfile), 
                 os.path.join(indir,regname+'.exons'), 
                 os.path.join(indir,regname+'.UTRs'), 
                 os.path.join(indir,regname+'.CNCs'), str(db), str(de), 
                 str(sel), '0']

        logfile = os.path.join(self.outdir,self.prefix,self.err,'log.build_input.txt')

        log = open(logfile,'w+')
         
        ex_annot = subprocess.Popen(annot, stdout=log, 
                                    preexec_fn=lambda: signal.signal(
                                        signal.SIGPIPE, signal.SIG_DFL),
                                    stderr=subprocess.STDOUT, close_fds = True)
        returncode = ex_annot.wait()

        # last perl command, recombination maps

        # first, get the min and max for the dense regions
        
        cmd = "cat "+logfile
        result = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, 
                                  stderr=subprocess.PIPE)
       
        out,err  = result.communicate()

        out_arr = out.split('\n')

        db_de = out_arr[len(out_arr)-2].split(' ')

        sim_min= db_de[0]
        sim_max= db_de[1]
        
        recombin = os.path.join(datadir,
                          'geneticMap/genetic_map_GRCh37_chr'+str(chr)+'.txt.gz')
   
        if recombout == '':
            recombout = os.path.join(indir, 'sfscode_rec_'+regname+'.txt')

        map = ['perl', os.path.join(basedir,'makeMapFile4sfscode.pl'), 
               recombin,str(sim_min), str(sim_max),recombout]

        ex_map = subprocess.Popen(map, stdout=log, stderr=subprocess.STDOUT, 
                                  preexec_fn=lambda: signal.signal(
                                        signal.SIGPIPE, signal.SIG_DFL),
                                  close_fds = True)
        returncode = ex_map.wait()

        return returncode

    # What I don't yet understand is how I am going to build in phenotype
    # If phenotype has no realtionship with selection, then I need a 
    # linear model for the outcome based on genes + environment
    # How will this be built into the simulation? Need an is_causal bit?
  
    def gutenkunst(self,add_on=False,nsim=1,N=10000,non_coding=False,
                   recombfile='',nsam=[100],mutation=[],t=0.001,rho=0.001,
                   loci=[], sel=[], L=[]):

        """
        A method that adds the Gutenkunst (2009, *PLoS Genetics*) model to
        an SFS_CODE command line. 
        """

        if not(add_on):

            #if self.sfs_code_loc == '':
            #    print "you must specify the location of the sfs_code binary",
            #    print "before calling this method!"
            #    print "self.sfs_code_loc = '/path/to/sfs_code'"
            #    exit(-1)

            self.line = []
            #self.line.append(self.sfs_code_loc)

        # 3 populations, nim
        self.line.extend(['3', str(nsim)])

        # N individuals
        self.line.extend(['-N', str(N)])

        # samples
        self.line.append('-n') 
        for thing in nsam:
            self.line.append(str(thing))

        # non_coding?
        if (non_coding):
            self.line.extend(['-a'])
    
        # no seq
        self.line.append('-A')

        if len(L) > 0:
            for thing in L:
                self.line.append(thing)

        # recomb

        if(recombfile == ''):
            self.line.extend(['-r',str(rho)])
    
        else:
            self.line.extend(['-r','F',recombfile,str(rho)])

        # Demographic events inferred under the model
        self.line.extend(['-TS', '0.359116', '0', '1'])
        self.line.extend(['-TS', '0.566244', '1', '2'])
        self.line.extend(['-TE', '0.629098'])
        self.line.extend(['-Td', '0', 'P', '0', '1.815229'])
        self.line.extend(['-Td', '0.359116', 'P', '1', '0.158014'])
        self.line.extend(['-Td', '0.566244', 'P', '1', '0.523478'])
        self.line.extend(['-Tg', '0.566244', 'P', '1', '49.49258'])
        self.line.extend(['-Td', '0.566244', 'P', '2', '0.280234'])
        self.line.extend(['-Tg', '0.566244', 'P', '2', '68.063324'])
        self.line.extend(['-Tm', '0.359116', 'P', '0', '1', '6.340369'])
        self.line.extend(['-Tm', '0.359116', 'P', '1', '0', '1.0018685'])
        self.line.extend(['-Tm', '0.566244', 'L', '0.759772', '0.546351', 
                          '0.062846', '0.204877', '0.024193', '0.10967725'])
          
        #loci
        if len(loci) > 0:
            for thing in mutation:
                self.line.append(thing)

        #selection model
        if len(sel) > 0:
            for thing in sel:
                self.line.append(thing)

        # mutation events
        if len(mutation) > 0:
            for thing in mutation:
                self.line.append(thing)
 
    def three_pop(self,add_on=False,nsim=1,N=10000,non_coding=False,
                   recombfile='',nsam=[100,100,100],mutation=[],t=0.001,rho=0.001,
                   loci=[], sel=[], L=[],t_end=0.60274,t_expand_p0=0,
                   expand_p0=1.68493, t_split_p0_p1=0.219178,
                   t_split_p1_p2=0.544658, bottle_p1_0=0.170732,
                   bottle_p1_1=0.47619,bottle_p2=0.242857,
                   growth_p1=58.4,growth_p2=80.3, mig_p0_p1=6.15,
                   mig_p1_p0=0.5, mig_all = [0.738,0.4674,0.06,0.192,
                   0.01938,0.09792], t_super=0.391465,model='gutenkunst'):
        
        """
        A general three population model with growth events.
        Inspired by Gutenkunst et al (2009, *PLoS Genetics*)
        and Gravel et al (2011, *PNAS*).  For a pictorial 
        representation see Figure 3A of the Gutenkunst paper.

        The default parameters are set to Gutenkunst et al.
        The maximum likelihood estimates of Gravel et al
        and Tennessen et al (2012, *Science*) are also 
        included.

        Use model='gravel', model='tennessen' or model='gutenkunst'
        to explicitly choose one of the models.

        The user can specify the model parameters as desired.
        To use a user specified model, simply use model = ''.
        Any unspecified parameters are set by default to the 
        parameters of the Gutenkunst model.       
        """
        
        if model == 'gutenkunst':
            t_end = 0.60274
            t_expand_p0=0
            expand_p0 = 1.68493
            t_split_p0_p1 = 0.219178
            t_split_p1_p2=0.544658
            bottle_p1_0=0.170732
            bottle_p1_1=0.47619
            bottle_p2=0.242857
            growth_p1=58.4
            growth_p2=80.3
            mig_p0_p1 = 6.15
            mig_p1_p0 = 0.5
            mig_all = [0.738, 0.4674,
                   0.06, 0.192, 0.01938, 0.09792]

        # this uses the esimates from the combined exons and low coverage
        elif model == 'gravel':
            t_end = 0.405479
            t_expand_p0=0
            expand_p0 = 1.982738
            t_split_p0_p1 = 0.265753
            t_split_p1_p2=0.342466
            bottle_p1_0=0.128575
            bottle_p1_1=0.554541
            bottle_p2=0.29554
            growth_p1=55.48
            growth_p2=70.08
            mig_p0_p1 = 4.3422
            mig_p1_p0 = 0.5583
            mig_all = [0.7237, 0.225794,
                   0.09305, 0.115754, 0.00858, 0.03421]            

        
        # these are the prepublication parameters from Dadi (provided by Ryan Hernandez)
        elif model == 'guten_old':
            t_end=0.629098
            t_expand_p0=0
            expand_p0 = 1.815229
            t_split_p0_p1 = 0.359116
            t_split_p1_p2=0.566244
            bottle_p1_0=0.158014
            bottle_p1_1=0.523478
            bottle_p2=0.280234
            growth_p1=49.49258
            growth_p2=68.0633424
            mig_p0_p1 = 6.340369
            mig_p1_p0 = 1.0018685
            mig_all = [0.759772, 0.546351,
                   0.062846, 0.204877, 0.024193, 0.10967725]


        elif model == 'tennessen':
            t_end = 0.405479
            t_expand_p0=0
            expand_p0 = 1.982738
            t_split_p0_p1 = 0.265753
            t_split_p1_p2=0.342466
            bottle_p1_0=0.128575
            bottle_p1_1=0.554541
            bottle_p2=0.29554
            growth_p1=44.822
            growth_p2=70.08
            growth_p0 = 242.36
            growth_p1_1= 284.7 
            mig_p0_p1 = 4.3422
            mig_p1_p0 = 0.5583
            mig_all = [0.7237, 0.225794,
                   0.09305, 0.115754, 0.00858, 0.03421]

        elif model != '':
            print "model", model, "not recongnized"
            print "please use one of the supported models"
            print "supported models: 'gutenkunst', 'gravel', 'tennessen', and ''",
            print "(user specified parameters)"
            exit()  
        
        if not(add_on):

            #if self.sfs_code_loc == '':
            #    print "you must specify the location of the sfs_code binary",
            #    print "before calling this method!"
            #    print "self.sfs_code_loc = '/path/to/sfs_code'"
            #    exit(-1)

            self.line = []
            #self.line.append(self.sfs_code_loc)

        # 3 populations, nim
        self.line.extend(['3', str(nsim)])

        # N individuals
        self.line.extend(['-N', str(N)])

        # samples
        self.line.extend(['-n'])

        if not isinstance(nsam, list) :
            print "Error: nsam must be a list"
            exit()

        for n in nsam:
            self.line.append(str(n))

        # non_coding?
        if (non_coding):
            self.line.extend(['-a','N'])
            if len(self.command.L) > 1:
                self.line.append('R')

        # no seq
        self.line.append('-A')

        if len(L) > 0:
            for thing in L:
                self.line.append(thing)

        self.line.extend(['-t',str(t)])

        # recomb

        if(recombfile == ''):
            self.line.extend(['-r',str(rho)])

        else:
            self.line.extend(['-r','F',recombfile,str(rho)])

        
        # Demographic events inferred under the model
        self.line.extend(['-TS', str(t_split_p0_p1), '0', '1'])
        self.line.extend(['-TS', str(t_split_p1_p2), '1', '2'])
        self.line.extend(['-TE', str(t_end)])
        self.line.extend(['-Td', str(t_expand_p0), 'P', '0', str(expand_p0)])
        self.line.extend(['-Td', str(t_split_p0_p1), 'P', '1', str(bottle_p1_0)])
        self.line.extend(['-Td', str(t_split_p1_p2), 'P', '1', str(bottle_p1_1)])
        self.line.extend(['-Tg', str(t_split_p1_p2), 'P', '1', str(growth_p1)])
        self.line.extend(['-Td', str(t_split_p1_p2), 'P', '2', str(bottle_p2)])
        self.line.extend(['-Tg', str(t_split_p1_p2), 'P', '2', str(growth_p2)])
        self.line.extend(['-Tm', str(t_split_p0_p1), 'P', '0', '1', str(mig_p0_p1)])
        self.line.extend(['-Tm', str(t_split_p0_p1), 'P', '1', '0', str(mig_p1_p0)])
        self.line.extend(['-Tm', str(t_split_p1_p2), 'L', str(mig_all[0]), str(mig_all[1]),
                          str(mig_all[2]), str(mig_all[3]), str(mig_all[4]), str(mig_all[5])])


        # tennessen
        if model == 'tennessen':
            self.line.extend(['-Tg', str(t_super), 'P', '0', str(growth_p0)])
            self.line.extend(['-Tg', str(t_super), 'P', '1', str(growth_p1_1)])

        #loci
        if len(loci) > 0:
            for thing in mutation:
                self.line.append(thing)

        #selection model
        if len(sel) > 0:
            for thing in sel:
                self.line.append(thing)

        # mutation events
        if len(mutation) > 0:
            for thing in mutation:
                self.line.append(thing)


class MSCommand(Command):


    def __init__(self,prefix='out',ms_loc=os.path.join('ms'),
                 err='err',outdir=os.path.join(os.getcwd(), 'sims')):

        self.nsim = 1
        self.ms_loc = ms_loc
        self.line = []
        self.prefix = prefix
        self.outdir=outdir
        self.err='err'
        if os.environ.has_key('SGE_TASK_ID'):
            self.task_id = os.environ['SGE_TASK_ID']
        else:
            self.task_id = 1


    def set_ms_loc(self,ms_loc,time=True):     

        from re import compile,search

        self.ms_loc = ms_loc

        if time:
            
            if len(self.line) > 0: 
                pattern = compile("ms") 
                for i in range(0,len(self.line)):
                    if(pattern.search(self.line[i])):
                        self.line[i] = self.ms_loc
                        break

class Reg():

    import subprocess

    def __init__(self,reg_loc='',
                 dir = os.path.join(os.getcwd(), 'ABC')):

        self.n_params = 2
        self.n_sums = 10
        self.reg_loc = reg_loc
        self.dir = dir
        try:
            os.stat(self.dir)
        except:
            os.mkdir(self.dir)
        self.out_file_base = os.path.join(self.dir, 'ABC')
        self.prior_file =  os.path.join(self.dir,'prior_test')
        self.data_file =  os.path.join(self.dir, 'data_test')
        self.tol = 0.02
        self.line = []

    def execute(self):

        self.line = []
  
        self.line.append(self.reg_loc)
        
        self.line.extend(['-P',str(self.n_params)])
        
        self.line.extend(['-S',str(self.n_sums)])
        
        self.line.extend(['-p',self.prior_file])
        
        self.line.extend(['-d',self.data_file])
        
        self.line.extend(['-b',self.out_file_base])
        
        self.line.extend(['-t',str(self.tol)])

        self.line.append('-L') # method of Beaumont et al

        # extend with other options

        out_h = open(os.path.join(self.dir,'log.txt'),'w')
 
        err_h = open(os.path.join(self.dir,'err.txt'),'w')

        process = subprocess.Popen(self.line, stdout=out_h,
                                   stderr=err_h, close_fds = True)

        code = process.wait()
