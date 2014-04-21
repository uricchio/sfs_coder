import matplotlib
import re
from collections import defaultdict

class PlotData:
 
    def __init__(self):

        self.annot_file = ''
        self.pi_file = ''
        self.sfs_file = ''
        self.regions = defaultdict(dict)
        self.bounds = []
        self.pi = defaultdict(dict)
        self.data = []
        self.sfs = []
        self.min = 0.
        self.max = 0.

    def get_annotations(self):

        """
        get the annotations for the simulated regions from
        an SFS_CODE annotation file.

        """

        if self.annot_file == '':
            print "Must define an annotation file: self.annot_file",
            print  " = <annot_file>"
            exit()
       
        try:
            fh = open(self.annot_file, 'r')
        except:
            print "cannot open file",
            print self.annot_file+'!',
            print "Is it a regular file?"
            exit()

   
        # get all the annotations for the simulated regions

        cur_pos = 0
        for line in fh:
            ann = line.strip(';\n').split(',')
            new_pos = int(ann[0])+cur_pos 
            if re.search(',',line):
                if ann[2] == 'C':
                    self.regions[cur_pos][new_pos] = ann[2]
                elif ann[2] == 'N':
                    if ann[3] == '0':
                        self.regions[cur_pos][new_pos] = ann[2]
                    else:
                        self.regions[cur_pos][new_pos] = 'CNC'
                #print cur_pos, new_pos, self.regions[cur_pos][new_pos]
            cur_pos = new_pos
        self.max = cur_pos
       
        for reg0 in sorted(self.regions):
            for reg1 in self.regions[reg0]:
                self.bounds.append([reg0,reg1])

        fh.close()

    def get_pi(self,pi_0=1.):

        if self.pi_file == '':
            print "Must define a datafile file: self.pi_file",
            print  " = <pi_file>"
            exit()

        try:
            fh = open(self.pi_file, 'r')
        except:
            print "cannot open file",
            print self.pi_file+'!',
            print "Is it a regular file?"
            exit()

        numlines = 0
        for line in fh:
            pis = line.strip().split(' ')
            i = 0
            for pi in pis:
                if self.bounds[i][0] in self.pi:
                    self.pi[self.bounds[i][0]][self.bounds[i][1]]+=float(pi)/pi_0
                else: 
                    self.pi[self.bounds[i][0]][self.bounds[i][1]]=float(pi)/pi_0
                i+=1
            numlines += 1

        for reg0 in self.pi:
            for reg1 in self.pi[reg0]:
                self.pi[reg0][reg1]/=numlines

        fh.close()

    def get_sfs(self):
    
        if self.sfs_file == '':
            print "Must define a data file: self.sfs_file",
            print  " = <sfs_file>"
            exit()

        try:
            fh = open(self.sfs_file, 'r')
        except:
            print "cannot open file",
            print self.sfs_file+'!',
            print "Is it a regular file?"
            exit()

        numlines=0
        sfs = []

        for line in fh:
            if len(sfs) == 0:
                sfs = line.split(' ')
                for i in range(0,len(sfs)):
                    sfs[i] = int(sfs[i])
                
            else:
                i = 0
                for thing in line.split(' '):
                    sfs[i] += int(thing)
                    i += 1                    
    
        self.sfs = sfs

    def pi_plot_regions(self,xmin=-1,xmax=-1):

        import matplotlib.pyplot as plt
        from scipy.interpolate import UnivariateSpline

        if xmin < 0:
            xmin = self.min
        if xmax < 0:
            xmax = self.max

        plt.plot()
        plt.xlim([xmin,xmax])
        plt.xlabel('Genomic position (bp)',fontsize=20)
        plt.ylabel(r'$\frac{\pi}{\pi_a}$',fontsize=20)
        plt.ylim([0,1.])

        x_cnc = []
        x_n = []
        x_c = []
        y_cnc = []
        y_n = []
        y_c = []

        for reg0 in self.regions:
            for reg1 in self.regions[reg0]:
                 yvals = []
                 yvals.extend([reg0,reg1])
                 if(self.regions[reg0][reg1]=='C'):
                     #plt.plot([reg0,reg1], [self.pi[reg0][reg1],self.pi[reg0][reg1]],color='r',linewidth=4)
                     y_c.append(self.pi[reg0][reg1])
                     x_c.append(reg0+(reg1-reg0)/2.)
                 if(self.regions[reg0][reg1]=='N'):
                     #plt.plot([reg0,reg1], [self.pi[reg0][reg1],self.pi[reg0][reg1]],color='g',linewidth=4)
                     y_n.append(self.pi[reg0][reg1])
                     x_n.append(reg0+(reg1-reg0)/2.)
                 if(self.regions[reg0][reg1]=='CNC'):
                     #plt.plot([reg0,reg1], [self.pi[reg0][reg1],self.pi[reg0][reg1]],color='b',linewidth=4)
                     y_cnc.append(self.pi[reg0][reg1])
                     x_cnc.append(reg0+(reg1-reg0)/2.)
        
                 
        #print x_c
        #print y_c

        x_c_gt = []
        y_c_gt = []

        num = 0
        for thing in y_c:
            if thing > 0.:
                y_c_gt.append(thing)
                x_c_gt.append(x_c[num])
            num+=1
  
        func = UnivariateSpline(x_c_gt,y_c_gt)
        
        xs = [xmin + i*(xmax-xmin)/1000. for i in range(1,1000)]
        ys = func(xs)

        plt.scatter(x_c,y_c,color='r',label='Coding')
        plt.scatter(x_cnc,y_cnc,color='b',label='Conserved non-coding')
        plt.scatter(x_n,y_n,color='g',label='Neutral non-coding')
        plt.plot(xs,ys,color='r')
        plt.legend(loc=1,scatterpoints=1)
        plt.show()


    def plot_sfs(self):

        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
 
        ax.plot([i for i in range(0,len(self.sfs))], self.sfs, 'bs')
        ax.set_yscale('log')

        plt.show()

