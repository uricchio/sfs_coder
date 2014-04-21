class Simulation:

    def __init__(self):

        self.chrs = []
        self.sfs = []

    def get_sfs(self,start=0,end=-1):
   
        if end == -1:
            end = len(self.chrs)
        
        freqs = [0 for i in range(0,len(self.chrs[0]))]
        self.sfs = [0 for i in range(start,end)]

        for chrnum in range(start,end):
            i = 0
            for var in self.chrs[chrnum]:
                freqs[i] += int(var)
                i+=1
         
        for freq in freqs:
            if freq < end-start and freq > 0:
                self.sfs[freq-1]+=1

        return self.sfs
