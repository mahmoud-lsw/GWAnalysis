#!/usr/bin/env python
import ROOT
import scipy as sp
from matplotlib import pyplot as plt
from scipy.stats import chi,chi2
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma', 'Bitstream Vera Sans',
                                 'Lucida Grande', 'Verdana']
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

triggertime=463917049.000
NMINEVENTS=2

tsmax=40
nbins=1
trials=1
normed=0
cumulative=0

X=sp.linspace(0,tsmax,nbins*tsmax)
#X=sp.logspace(-1,1.5,20)
print X
def p_value(x,dof=2,trials=1):
    C2=chi2(dof)
    P=C2.sf(x)
    y=1.0-sp.power(1.0-P,trials)
    return y


def cumulativeDistribution(myarray,ts):
    myarray=sp.array(myarray)
    N=len(myarray)
    Nabove=1.0*len(myarray[myarray>ts])
    return Nabove/N

def cumulativeDistribution_H(histogram,normed=0):
    nbin=len(histogram)
    cumulative=histogram
    for i in range(1,nbin):
        cumulative[nbin-i-1]+=float(cumulative[nbin-i])
        #print cumulative[nbin-i-1],cumulative[nbin-i]
    if normed==0 and not cumulative[0]==0: return cumulative
    else: return 1.0*cumulative/cumulative[0]

def fromROOTFile(myFilename):
    myFile=ROOT.TFile(myFilename,"READ")
    myTree=myFile.Get("GRBCatalog")
    N=myTree.GetEntries()
    print myFilename
    #print '%s, %d entries' % (myFilename,N)
    TSTART   = []
    TSTOP    = []
    TS_GRB   = []
    TS_TSMAP = []
    NEVT_GRB = []
    PRINT_FIRST=1
    for i in range(N):
        myTree.GetEntry(i)
        ts_grb  = myTree.LIKE_UP_TS_GRB
        evt_grb = myTree.LIKE_UP_gtsrcprob_Nthr
        tstart  = myTree.GRBMET+myTree.LIKE_UP_T0
        tstop  = myTree.GRBMET+myTree.LIKE_UP_T1
        ts_tsmap=myTree.TSMAP_MAX
        #ts_grb=ts_tsmap
        TSTART.append(tstart)
        TSTOP.append(tstop)
        TS_GRB.append(ts_grb)
        TS_TSMAP.append(ts_tsmap)
        NEVT_GRB.append(evt_grb)
        LAST=(tstart-triggertime)
        if ts_grb>20 and evt_grb>2:
            print '-----------------------------'
            print myTree.GRBNAME, tstart,myTree.LIKE_UP_T1,myTree.EMIN
            print '-----------------------------'
        if LAST>0 and PRINT_FIRST:
            PRINT_FIRST=0
            print 'FIRST:==>',myFilename,(tstart-triggertime)/86400,(tstop-triggertime)/86400,ts_grb,evt_grb
            pass
        if ts_grb>20 and evt_grb>NMINEVENTS: print myFilename,(tstart-triggertime)/86400,(tstop-triggertime)/86400,ts_grb,evt_grb
        pass
    return sp.array(TSTART),sp.array(TSTOP),sp.array(TS_GRB),sp.array(TS_TSMAP),sp.array(NEVT_GRB)




TSTART_TOT=[]
TSTOP_TOT=[]
TS_GRB_TOT    = []
TS_GRB        = []
TS_TSMAP_TOT  = []
NEVT_GRB_TOT  = []
mydirectories=['GW150914/LATBA-2-GW150914-POS00.root',
               'GW150914/LATBA-2-GW150914-POS01.root',
               'GW150914/LATBA-2-GW150914-POS02.root',
               'GW150914/LATBA-2-GW150914-POS03.root',
               'GW150914/LATBA-2-GW150914-POS04.root',
               'GW150914/LATBA-2-GW150914-POS05.root',
               'GW150914/LATBA-2-GW150914-POS06.root',
               'GW150914/LATBA-2-GW150914-POS07.root',
               'GW150914/LATBA-2-GW150914-POS08.root'
               ]


#mydirectories=['GW150914/GW150914-POS00/','GW150914/GW150914-POS01/','GW150914/GW150914-POS02/','GW150914/GW150914-POS03/','GW150914/GW150914-POS04/','GW150914/GW150914-POS05/','GW150914/GW150914-POS06/','GW150914/GW150914-POS07/','GW150914/GW150914-POS08/']

#mydirectories=['GW150914/LATBA-2-GW150914-POS01.root']

plt.figure(figsize=(10,7),facecolor='w')
colors=['k','r','g','b','m','y','c','0.75','0.5']
for i,myfile in enumerate(mydirectories):
    tstart,tstop,ts_grb,ts_tsmap,nevt_grb=fromROOTFile(myfile)

    TSTART_TOT=sp.append(TSTART_TOT,tstart)
    TSTOP_TOT=sp.append(TSTOP_TOT,tstop)
    TS_GRB_TOT=sp.append(TS_GRB_TOT,ts_grb)
    TS_TSMAP_TOT=sp.append(TS_TSMAP_TOT,ts_tsmap)
    NEVT_GRB_TOT=sp.append(NEVT_GRB_TOT,nevt_grb)

    filter_envt = (nevt_grb>NMINEVENTS)
    
    ts_grb_before=ts_grb[tstop<triggertime]
    ts_grb_after =ts_grb[tstop>triggertime]
    ts_grb_before_selected=ts_grb[(tstop<triggertime)*filter_envt]
    ts_grb_after_selected=ts_grb[(tstop>triggertime)*filter_envt]
    
    print myfile, len(ts_grb),len(ts_grb_before),len(ts_grb_after),len(ts_grb_before_selected),len(ts_grb_after_selected)
    label1='ROI (%d) Before' % i
    label2='ROI (%d) After' % i
    try:
        plt.hist((ts_grb*filter_envt)[(tstop<triggertime)],bins=X,normed=normed, histtype='step', cumulative=cumulative,label=label1,color=colors[i])
    except:
        pass
    try:
        plt.hist((ts_grb*filter_envt)[(tstop>triggertime)],bins=X,normed=normed, histtype='step', cumulative=cumulative,label=label2,color=colors[i],linestyle='--')
    except:
        pass
    pass
TS_GRB = TS_GRB_TOT*(NEVT_GRB_TOT>NMINEVENTS)
TS_GRB_BEFORE = TS_GRB[(TSTOP_TOT<triggertime)]
TS_GRB_AFTER = TS_GRB[(TSTOP_TOT>triggertime)]

plt.yscale('log')
plt.legend(fontsize=20,frameon=False)     
plt.xlabel('Test Statistic (TS$_{SRC}$)',size=20)
plt.ylabel('Number of Passages',size=20)
plt.ylim([0.5,2*len(TS_GRB)])
if normed: 
    plt.ylim([5e-5,2])
    plt.ylabel('Post-trial probability',size=20)
    pass
plt.xticks(size=20)
plt.yticks(size=20)
##################################################
print 'N TS distribution (all)...:',len(TS_GRB_TOT)
print 'N TS distribution (sel)...:',len(TS_GRB[TS_GRB>0])
print 'N TS before............:',len(TS_GRB_BEFORE)
print 'N TS after............:',len(TS_GRB_AFTER)
if len(TS_GRB_BEFORE)>0:  print 'Max TS before...........: %.1f (%.1e)' %(TS_GRB_BEFORE.max(),p_value(TS_GRB_BEFORE.max()))
if len(TS_GRB_AFTER)>0:   print 'Max TS after...........: %.1f (%.1e)' %(TS_GRB_AFTER.max(),p_value(TS_GRB_AFTER.max()))

plt.figure(figsize=(10,7),facecolor='w')

try:
    plt.hist(TS_GRB_BEFORE,bins=X,normed=normed, histtype='step', cumulative=cumulative,label='All Before',color='b')
except:
    pass
try:
    plt.hist(TS_GRB_AFTER,bins=X,normed=normed, histtype='step', cumulative=cumulative,label='All After',color='r')
except:
    pass
plt.yscale('log')
plt.legend(fontsize=20,frameon=False)     
plt.xlabel('Test Statistic (TS$_{SRC}$)',size=20)
plt.ylabel('Number of Passages',size=20)
plt.ylim([0.5,2*len(TS_GRB)])

if normed: 
    plt.ylim([5e-5,2])
    plt.ylabel('Post-trial probability',size=20)
plt.xticks(size=20)
plt.yticks(size=20)

plt.figure(figsize=(10,7),facecolor='w')

plt.hist(TS_GRB_TOT,bins=X,normed=normed, histtype='step', cumulative=cumulative,label='Data: All Events',color='b',lw=2)
plt.hist(TS_GRB,bins=X,normed=normed, histtype='step', cumulative=cumulative,label='Data: Selected Events (N$_{\gamma}>$%d)' % NMINEVENTS,color='r',lw=2)
########################################################################
## ADD SIMULATIONS:
########################################################################
TSTART_TOT    = []
TSTOP_TOT     = []
TS_GRB_TOT    = []
TS_GRB        = []
TS_TSMAP_TOT  = []
NEVT_GRB_TOT  = []
mydirectories=['GW150914/LATBA-2-GW150914-SIM-POS00.root',
               'GW150914/LATBA-2-GW150914-SIM-POS01.root',
               'GW150914/LATBA-2-GW150914-SIM-POS02.root',
               'GW150914/LATBA-2-GW150914-SIM-POS03.root',
               'GW150914/LATBA-2-GW150914-SIM-POS04.root',
               'GW150914/LATBA-2-GW150914-SIM-POS05.root',
               'GW150914/LATBA-2-GW150914-SIM-POS06.root',
               'GW150914/LATBA-2-GW150914-SIM-POS07.root',
               'GW150914/LATBA-2-GW150914-SIM-POS08.root'
               ]
for i,myfile in enumerate(mydirectories):
    tstart,tstop,ts_grb,ts_tsmap,nevt_grb=fromROOTFile(myfile)

    TSTART_TOT=sp.append(TSTART_TOT,tstart)
    TSTOP_TOT=sp.append(TSTOP_TOT,tstop)
    TS_GRB_TOT=sp.append(TS_GRB_TOT,ts_grb)
    TS_TSMAP_TOT=sp.append(TS_TSMAP_TOT,ts_tsmap)
    NEVT_GRB_TOT=sp.append(NEVT_GRB_TOT,nevt_grb)

    filter_envt = (nevt_grb>NMINEVENTS)
    
    ts_grb_before=ts_grb[tstop<triggertime]
    ts_grb_after =ts_grb[tstop>triggertime]
    ts_grb_before_selected=ts_grb[(tstop<triggertime)*filter_envt]
    ts_grb_after_selected=ts_grb[(tstop>triggertime)*filter_envt]    
    print myfile, len(ts_grb),len(ts_grb_before),len(ts_grb_after),len(ts_grb_before_selected),len(ts_grb_after_selected)
    pass
TS_GRB = TS_GRB_TOT*(NEVT_GRB_TOT>NMINEVENTS)
TS_GRB_BEFORE = TS_GRB[(TSTOP_TOT<triggertime)]
TS_GRB_AFTER = TS_GRB[(TSTOP_TOT>triggertime)]

#tstart,tstop,ts_grb,ts_tsmap,nevt_grb=fromROOTFile('Analysis_GI/LATBA-2-GRB081024-SIMULATED-2.root')
#ts_grb_sel=ts_grb*(nevt_grb>NMINEVENTS)
#print 'MC:', len(ts_grb),len(ts_grb[nevt_grb>NMINEVENTS])

plt.hist(TS_GRB_TOT,bins=X,normed=normed, histtype='step', cumulative=cumulative,label='MC: All Events',color='b',linestyle='--',lw=2)

#hy,hx=sp.histogram(TS_GRB_TOT,bins=X)
#HX=0.5*(hx[1:]+hx[:-1])
#Y1=hy+sp.sqrt(hy)
#Y1=cumulativeDistribution_H(Y,1)
#Y2=hy-sp.sqrt(hy)
#Y2=cumulativeDistribution_H(Y,1)
#plt.fill_between(HX, Y1, Y2,facecolor='yellow',alpha=0.5,label='MC')

plt.hist(TS_GRB,bins=X,normed=normed, histtype='step', cumulative=cumulative,label=r'MC: Selected Events (N$_{\gamma}>$%d)' % NMINEVENTS,color='r',linestyle='--',lw=2)
#hy,hx=sp.histogram(TS_GRB,bins=X)
#HX=0.5*(hx[1:]+hx[:-1])
#Y1=hy+sp.sqrt(hy)
#Y1=cumulativeDistribution_H(Y,1)
#Y2=hy-sp.sqrt(hy)
#Y2=cumulativeDistribution_H(Y,1)
#plt.fill_between(HX, Y1, Y2,facecolor='yellow',alpha=0.5,label='')

#hy,hx=sp.histogram(TS_GRB,bins=X)#,normed=0, histtype='step', cumulative=-1,label='Slected Events (MC)',color='k',linestyle='-',lw=2)
#plt.plot(X[:-1],hy,'b')
#plt.plot(X[:-1],hy+sp.sqrt(hy),'r')
#plt.plot(X[:-1],hy-sp.sqrt(hy),'g')

#Y=cumulativeDistribution_H(hy,1)
#plt.plot(X[:-1],Y,'b')
#Y=cumulativeDistribution_H(hy+sp.sqrt(hy),1)
#plt.plot(X[:-1],Y,'0.5')
#Y=cumulativeDistribution_H(hy-sp.sqrt(hy),1)
#plt.plot(X[:-1],Y,'0.5')

#plt.hist(TS_GRB_TOT+sp.sqrt(TS_GRB_TOT),bins=X,normed=normed, histtype='step', cumulative=-1,label='',color='0.5',linestyle='-')
#plt.hist(TS_GRB+sp.sqrt(TS_GRB),bins=X,normed=normed, histtype='step', cumulative=-1,label='',color='0.5',linestyle='-')
#plt.hist(sp.maximum(0,TS_GRB_TOT-sp.sqrt(TS_GRB_TOT)),bins=X,normed=normed, histtype='step', cumulative=-1,label='',color='0.5',linestyle='-')
#plt.hist(sp.maximum(0,TS_GRB-sp.sqrt(TS_GRB)),bins=X,normed=normed, histtype='step', cumulative=-1,label='',color='0.5',linestyle='-',lw=2)

#############################
# DRAW CHI 2
#############################

if 0:
    p_coincident=p_value(TS_coincident,dof=2,trials=trials)
    print 'TS=%.1f, p-value=%.1e' %(TS_coincident,p_coincident)
    plt.annotate('GRB (TS=%.1f)' % (TS_coincident), \
                 xy=(TS_coincident,p_coincident), \
                 xycoords='data',xytext=(0, 200),textcoords='offset points',\
                 fontsize=20,
                 horizontalalignment='center',
                 arrowprops=dict(facecolor='green', shrink=0.05))
    #             arrowprops=dict(arrowstyle="->",fc=(1.0, 0.7, 0.7)))

#for trials in sp.linspace(100,200,10):
#    dof=2
#    if 1: plt.plot(X,p_value(X,dof=dof,trials=trials),'--',label='chi(dof=%d,trials=%.d)'% (dof,trials))
#    pass
#if 0: plt.plot(X,p_value(X,dof=4,trials=2),'-',label='Distribution from MC')
#plt.plot(X,p_value(X,2,1),'r--',label='$\\chi^{2}_{dof=2}$, no trial factors')
#plt.plot(X,p_value(X,2,trials),'b--',label='$\\chi^{2}_{dof=2}$, N=%s' % trials)
#plt.plot(X,p_value(X,8,1),'b--',label='$\\chi^{2}_{dof=8}$, N=%s' % trials)

plt.yscale('log')
plt.legend(fontsize=20,frameon=False)     
plt.xlabel('Test Statistic (TS$_{SRC}$)',size=20)
plt.ylabel('Number of Passages',size=20)
plt.ylim([0.5,2*len(TS_GRB_TOT)])
if normed: 
    plt.ylim([5e-5,2])
    plt.ylabel('Post-trial probability',size=20)
plt.xticks(size=20)
plt.yticks(size=20)

plt.show();


