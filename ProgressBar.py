#!/usr/bin/env python
import sys

class progressbar:
    def __init__(self,toolbar_width=10):
        self.toolbar_width=toolbar_width        
        sys.stdout.write("[%s]" % (" " * toolbar_width))
        sys.stdout.flush()
        sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line, after '['
        pass
    
    def go(self,i,n):
        #print 'n,self.toolbar_width',n, self.toolbar_width
        x=int(n/self.toolbar_width)
        if x==0:
            sys.stdout.write("*")
            sys.stdout.flush()
        elif (i % x)==0:
            sys.stdout.write("*")
            sys.stdout.flush()
            pass
        if i==(n-1): sys.stdout.write("\n")
        pass
    pass
if __name__=='__main__':
    pb=progressbar(10)
    try:
        n=int(sys.argv[1])
    except:
        n=100000
        pass
    
    for i in range(n):
        pb.go(i,n)
        pass
    
            
