#!/usr/bin/env python
from matplotlib import pyplot as plt
from check_file_exists import check_file_exists
from matplotlib import rc
import glob
rc('text', usetex=True)

import scipy as sp
import numpy as np
import os
import argparse
from scipy.stats import chi2

def readFilie(file_in):
    data = np.recfromtxt(file_in, names=True)    
    name = data['name']
    pos  = data['pos']
    ra   = data['ra']
    dec  = data['dec']
    flux = data['flux']
    flux_err = data['flux_err']
    ts_src   = data['tssrc']
    nevt     =data['nprob']
    ts_map   =data['tsmap']
    print 'in %s found %d entries' %(file_in,len(name))
    #for i in range(len(name)): print name[i],pos[i],ts_src[i],nevt[i],ts_map[i]
    return name,pos,ra,dec,flux,flux_err,ts_src,nevt,ts_map


if __name__=="__main__":
        
    desc = '''Compare distribution with simulated and real data'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--file_data',help='Input File contianing Data results', type=check_file_exists, required=True)
    parser.add_argument('--file_simu',help='Input File contianing Simulation results', type=check_file_exists, required=False)
    parser.add_argument('--plot_tsmap', help='plot the tsmap file (if available)',
                        required=False, type=bool, default=True)
    parser.add_argument('--downloadfromSLAC', help='Download the data from SLAC (need ssh)',
                        required=False, type=bool, default=False)
    
    parser.add_argument('--cumulative', help='Draw the cumpulative TS distribution',
                        required=False, type=bool, default=False)

    
    args = parser.parse_args()
    evt_min = -1
    ts_min = 20
    tsmax  = 1
    nmax=15
    
    # data:
    data_slac_base='/afs/slac/g/glast/groups/grb/DATA-PASS8/GW151226_N8/GW151226_P8_TRANSIENT010E_-864000_864000_10.0_POS'
    data_local='DATA'
    name_d,pos_d,ra_d,dec_d,flux_d,flux_err_d,ts_src_d,nevt_d,ts_map_d=readFilie(args.file_data)
    tsmax=max(tsmax,ts_src_d.max())
    _myilter_d=sp.logical_and(nevt_d>=evt_min,ts_src_d>ts_min)
    
    print '       DATA: Total =%d Selected : %s' %(len(name_d),len(name_d[_myilter_d]))
    for i in range(len(name_d[_myilter_d])):
        #print name_d[_myilter_d][i],pos_d[_myilter_d][i],nevt_d[_myilter_d][i],ts_src_d[_myilter_d][i]
        if args.downloadfromSLAC:
            local_destination = '%s/%s_POS%03d' %(data_local,name_d[_myilter_d][i],pos_d[_myilter_d][i])
            print local_destination
            os.system('mkdir -pv %s ' %local_destination)
            cmd='scp -r rhel6-64.slac.stanford.edu:%s%03d/%s %s/%s_POS%03d' %(data_slac_base,pos_d[_myilter_d][i],name_d[_myilter_d][i],data_local,name_d[_myilter_d][i],pos_d[_myilter_d][i])
            print cmd
            os.system(cmd)
            pass
        if args.plot_tsmap:
            cmd='draw_tsmap.py %s/%s_POS%03d' %(data_local,name_d[_myilter_d][i],pos_d[_myilter_d][i])
            print cmd
            os.system(cmd)
            pass
        pass
    print args.file_simu
    if args.file_simu is not None:
        # simulation:
        simu_slac_base='/afs/slac/g/glast/groups/grb/DATA-PASS8/GW151226/GW151226_P8_TRANSIENT010E_-864000_864000_SIM_POS'
        simu_local='SIMU'    
        name_s,pos_s,ra_s,dec_s,flux_s,flux_err_s,ts_src_s,nevt_s,ts_map_s=readFilie(args.file_simu)
        tsmax=max(tsmax,ts_src_s.max())
        _myilter_s=sp.logical_and(nevt_s>=evt_min,ts_src_s>ts_min)
        print '       SIMULATIONS: Total =%d Selected : %s' %(len(name_s),len(name_s[_myilter_s]))

        for i in range(len(name_s[_myilter_s])):
            #print name_s[_myilter_s][i],pos_s[_myilter_s][i],nevt_s[_myilter_s][i],ts_src_s[_myilter_s][i]
            if downloadfromSLAC:
                cmd='scp -r rhel6-64.slac.stanford.edu:%s%03d/%s %s/%s_POS%03d' %(simu_slac_base,pos_s[_myilter_s][i],name_s[_myilter_s][i],simu_local,name_s[_myilter_s][i],pos_s[_myilter_s][i])
                print cmd
                os.system(cmd)
                pass
            if plot_tsmap:
                cmd='./draw_tsmap.py %s/%s_POS%03d' %(simu_local,name_s[_myilter_s][i],pos_s[_myilter_s][i])
                print cmd
                os.system(cmd)
                pass
            pass
        pass
    
    fig=plt.figure(figsize=(10,7),facecolor='w')
    ax=fig.add_subplot(1,1,1)
    cumulative=args.cumulative
    if cumulative:
        xbins=sp.linspace(0,tsmax,20*tsmax+1)
        plt.hist(ts_src_d,bins=xbins,label='Data',histtype='step',color='b',lw=2,cumulative=-1,normed=1)
        if args.file_simu is not None: plt.hist(ts_src_s,bins=xbins,label='Simulation',histtype='step',color='r',lw=2,linestyle='--',cumulative=-1,normed=1)
        ndf=2
        Y=chi2.sf(xbins,ndf)/2.
        plt.plot(xbins,Y,label=r'1/2 $\chi^{2}_{%d}$' % ndf, color='g')
    else:
        xbins=sp.linspace(0,tsmax,2*tsmax+1)
        plt.hist(ts_src_d,bins=xbins,label='Data',histtype='step',color='b',lw=2)#,cumulative=-1)#,normed=1
        if args.file_simu is not None: plt.hist(ts_src_s,bins=xbins,label='Simulation',histtype='step',color='r',lw=2,linestyle='--')#,cumulative=-1)#,normed=1sp.logspace(-1,3,10))
        pass
    #ndf=2
    #Y=chi2.sf(xbins,ndf)/2.
    #plt.plot(xbins,Y,label=r'1/2 $\chi^{2}_{%d}$' % ndf)
    #    ndf=1
    #    Y=chi2.sf(xbins,ndf)/1.
    #    plt.plot(xbins,Y,label=r'$\chi^{2}_{%d}$' % ndf)
    plt.xlabel("TS",size=20)
    if cumulative:
        plt.ylabel("Fraction of RoIs with $>$ TS",size=20)
    else:
        plt.ylabel("Number of ROI",size=20)
        pass
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.yscale('log')
    plt.legend(frameon=0)
    if cumulative:
        ax.set_ylim([1e-4,1.1])
    else:
        ax.set_ylim([0.1,2000])
        pass
    fig=plt.figure(figsize=(15,15),facecolor='w')

    #xbins=sp.logspace(-1,2,50)
    xbins=sp.linspace(0,tsmax,2*tsmax+1)
    
    ax=fig.add_subplot(3,2,1)
    if args.file_simu is not None: plt.hist(ts_map_s,bins=xbins,label='simulation',histtype='step')#sp.logspace(-1,3,10))
    plt.hist(ts_map_d,bins=xbins,label='data',histtype='step')#sp.logspace(-1,3,10))
    ax.set_xlabel("TS map")
    plt.yscale('log')
    plt.legend()
    
    ax=fig.add_subplot(3,2,2)
    if args.file_simu is not None: plt.hist(ts_src_s,bins=xbins,label='simulation',histtype='step',normed=1)#,cumulative=-1)#sp.logspace(-1,3,10))
    plt.hist(ts_src_d,bins=xbins,label='data',histtype='step',normed=1)#
    #ndf=2
    #Y=chi2.sf(xbins,ndf)/2.
    #plt.plot(xbins,Y,label=r'1/2 $\chi^{2}_{%d}$' % ndf)
    #    ndf=1
    #    Y=chi2.sf(xbins,ndf)/1.
    #    plt.plot(xbins,Y,label=r'$\chi^{2}_{%d}$' % ndf)
    ax.set_xlabel("TS src")
    plt.yscale('log')
    plt.legend()
    ax=fig.add_subplot(3,2,3)
    xbins=sp.linspace(0,nmax,4*nmax+1)
    print xbins
    if args.file_simu is not None: plt.hist(nevt_s,bins=xbins,label='simulation',histtype='step',normed=1,cumulative=-1)#sp.logspace(-1,3,10))
    plt.hist(nevt_d,bins=xbins,label='data',histtype='step')#sp.logspace(-1,3,10))
    ax.set_xlabel("Nevt")
    plt.yscale('log')
    plt.legend()
    
    ax=fig.add_subplot(3,2,4)
    xbins=sp.logspace(-10,-6, 20)
    if args.file_simu is not None: plt.hist(flux_s,bins=xbins,label='simulation',histtype='step')#sp.logspace(-1,3,10))
    plt.hist(flux_d,bins=xbins,label='data',histtype='step')#sp.logspace(-1,3,10))
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    
    ax=fig.add_subplot(3,2,5)
    if args.file_simu is not None: plt.plot(ts_src_s,nevt_s,label='simulation',ls='None',marker='o')#sp.logspace(-1,3,10))
    plt.plot(ts_src_d,nevt_d,label='data',ls='None',marker='o')#sp.logspace(-1,3,10))
    ax.set_xlabel("TS src")
    ax.set_ylabel("Nevt")
    plt.axhspan(evt_min, nmax, color='y', alpha=0.5, lw=0)
    plt.axvspan(ts_min, tsmax, color='y', alpha=0.5, lw=0)
    
    #plt.xscale('log')
    #plt.yscale('log')
    plt.legend()
    
    ax=fig.add_subplot(3,2,6)
    if args.file_simu is not None: plt.plot(ts_map_s,nevt_s,label='simulation',ls='None',marker='o')#sp.logspace(-1,3,10))
    plt.plot(ts_map_d,nevt_d,label='data',ls='None',marker='o')#sp.logspace(-1,3,10))
    #plt.xscale('log')
    #plt.yscale('log')
    plt.legend()
    ax.set_xlabel("TS map")
    ax.set_ylabel("Nevt")
    
    
    #plt.errorbar(flux_s,flux_d,xerr=flux_err_s,yerr=flux_err_d,label='flux')#sp.logspace(-1,3,10))
    #ax.set_xlabel("data")
    #ax.set_ylabel("simulations")
    
    plt.show()
    
