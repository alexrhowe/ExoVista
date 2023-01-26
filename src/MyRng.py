import numpy as np

# Use the first prime number >10,000,000 = 10,000,019

class MyRng():
    nglobal=0
    nlist = np.zeros(10000019)
    def __init__(self,seed=0):
        fin = open('rng10M.dat','r')
        for i in range(0,10000019): self.nlist[i] = float(fin.readline())
        self.nglobal = seed
        
    def random(self,n=1):
        self.nglobal += n
        if self.nglobal > len(self.nlist): self.nglobal = n
        
        return self.nlist[self.nglobal-n:self.nglobal]
    
    def integers(self,nmin=0,nmax=1,n=1):
        self.nglobal += n
        if self.nglobal > len(self.nlist): self.nglobal = n
        
        nfloat = self.nlist[self.nglobal-n:self.nglobal]
        nint = np.array(nfloat)*(nmax-nmin+1) + nmin
        intlist = np.full(n,0)
        for i in range(0,n): intlist[i] = int(nint[i])

        return intlist

    def shuffle(self,inlist):
        n = len(inlist)
        ilist = self.random(n)
        shuffled = inlist[np.argsort(ilist)]
        
        return shuffled
