#!/usr/bin/env python
import argparse
import matplotlib
import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
from matplotlib.ticker import NullFormatter
import matplotlib.colors as colors
import matplotlib.cm as cmx

from matplotlib import rc
from check_file_exists import check_file_exists
import glob

rc('text', usetex=True)

rc('text', usetex=True)

def plotOneFile(ax,text_file_name,names=False,draw_ul=False):
    MeV2erg=1.6027e-6
    data = np.recfromtxt(text_file_name, names=True)
    x   = data['median']
    dx0 = x-data['start']
    dx1 = data['end']-x
    y   = data['EnergyFlux']*MeV2erg
    dy  = data['ErrorEF']*MeV2erg
    f=np.logical_and(dy>0,dy<y)
    ul=np.logical_not(f)
    # detections:
    name=''
    if names:
        name='GRB%s' % text_file_name.split('_')[-1].replace('.txt','')        
    if x[f].min()<1000:
        ax.errorbar(x[f],y[f],xerr=(dx0[f],dx1[f]),yerr=dy[f],ls='None',capsize=0,linewidth=3,alpha=1,label=name)
        if draw_ul: ax.errorbar(x[ul],y[ul],xerr=(dx0[ul],dx1[ul]),yerr=0.5*y[ul], uplims=True,ls='None',capsize=0,linewidth=3,alpha=1)
    
    #print '---->',text_file_name
    #print x[f].min(),x[f].max(),y[f]


if __name__=="__main__":
    desc = '''Plot series of LAT Lighlycurces'''
    parser = argparse.ArgumentParser(description=desc)
        
    parser.add_argument('--dir',help='directory containing the LC files',
                        required=True, default=None)
        
    parser.add_argument('--out_plot',help='Name for the output file (the format will be guessed from the extension)',
                        required=False, type=str)
    
    parser.add_argument('--tmin',help='minimum time for the x-axis',
                        required=False, default=None)

    parser.add_argument('--tmax',help='maximum time for the x-axis',
                        required=False, default=None)

    parser.add_argument('--names',help='Add a legend with GRB names',action='store_const',
                        const=True,default=False)
    
    parser.add_argument('--ul',help='Add also the UL',action='store_const',
                        const=True, default=False)
    
    args = parser.parse_args()
    
    fig = plt.figure(figsize=(15, 8),facecolor='w')
    
    ax = plt.axes()
    
    file_list=glob.glob('%s/*' %args.dir)
    for f in file_list:
        plotOneFile(ax,f,args.names,args.ul)
        pass
    if args.names:
        plt.legend()
    plt.yscale('log')
    plt.xscale('log')

    ax.set_xlabel('Time since t$_{GW}$ [s]',size=15)
    ax.set_ylabel('Flux [erg/s/cm$^{2}$]',size=15)
    plt.xticks(size=20)
    plt.yticks(size=20)
    
    plt.show()
