#!/usr/bin/env python
from Analysis import *
from matplotlib import pyplot as plt

import sys
import glob

myfolder   = sys.argv[1]
result_file=glob.glob('%s/results_*.txt' % myfolder)[0]
myName=result_file.split('_')[-1].replace('.txt','')
print myfolder,result_file,myName
l=LAT(myfolder,myName)
l.draw()
l.addFGL(print_sources=True)
l.save()
#plt.show()
