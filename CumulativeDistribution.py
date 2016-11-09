#!/usr/bin/env python
from Analysis import *
from scipy.stats import chi,chi2
triggertime=463917049.000
NMINEVENTS=2

def p_value(x,dof=2,trials=1):
    C2=chi2(dof)
    P=C2.sf(x)
    y=1.0-sp.power(1.0-P,trials)
    return y

tsmax=35
nbins=100
trials=1
normed=1

X=sp.linspace(0,tsmax,nbins*tsmax)

TSTART_TOT=[]
TSTOP_TOT=[]
TS_GRB_TOT    = []
TS_GRB        = []
TS_TSMAP_TOT  = []
NEVT_GRB_TOT  = []


mydirectories=['GW150914/GW150914-POS00/','GW150914/GW150914-POS01/','GW150914/GW150914-POS02/','GW150914/GW150914-POS03/','GW150914/GW150914-POS04/','GW150914/GW150914-POS05/','GW150914/GW150914-POS06/','GW150914/GW150914-POS07/','GW150914/GW150914-POS08/']

#mydirectories=['GW150914/GW150914-POS05/']
plt.figure(figsize=(10,7),facecolor='w')
colors=['k','r','g','b','m','y','c','0.75','0.5']
for i,mydir in enumerate(mydirectories):
    
    tstart   = sp.load('%s/tstart.npy' % mydir)
    tstop       = sp.load('%s/tstop.npy' % mydir)
    ts_grb     = sp.load('%s/ts_grb.npy' % mydir)
    ts_tsmap      = sp.load('%s/ts_tsmap.npy' % mydir)
    nevt_grb = sp.load('%s/nevt_grb.npy' % mydir)

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

    print mydir, len(ts_grb),len(ts_grb_before),len(ts_grb_after),len(ts_grb_before_selected),len(ts_grb_after_selected)
    label1='ROI (%d) Before' % i
    label2='ROI (%d) After' % i
#    if i==0:
#        label1='Before'
#        label2='After'
#    else:
#        label1=''
#        label2=''
#        pass
    plt.hist((ts_grb*filter_envt)[(tstop<triggertime)],bins=X,normed=normed, histtype='step', cumulative=-1,label=label1,color=colors[i])
    plt.hist((ts_grb*filter_envt)[(tstop>triggertime)],bins=X,normed=normed, histtype='step', cumulative=-1,label=label2,color=colors[i],linestyle='--')

    pass
TS_GRB = TS_GRB_TOT*(NEVT_GRB_TOT>NMINEVENTS)
TS_GRB_BEFORE = TS_GRB[(TSTOP_TOT<triggertime)]
TS_GRB_AFTER = TS_GRB[(TSTOP_TOT>triggertime)]

plt.yscale('log')
plt.legend(fontsize=20,frameon=False)     
plt.xlabel('Test Statistic (TS)',size=20)
plt.ylabel('Number of events',size=20)
plt.ylim([0.5,2*len(TS_GRB)])
if normed: 
    plt.ylim([5e-5,2])
    plt.ylabel('Post-trial probability',size=20)
    pass
plt.xticks(size=20)
plt.yticks(size=20)

print 'N TS distribution (all)...:',len(TS_GRB_TOT)
print 'N TS distribution (sel)...:',len(TS_GRB_TOT)
print 'N TS before............:',len(TS_GRB_BEFORE)
print 'N TS after............:',len(TS_GRB_AFTER)
if len(TS_GRB_BEFORE)>0:  print 'Max TS before...........: %.1f (%.1e)' %(TS_GRB_BEFORE.max(),p_value(TS_GRB_BEFORE.max()))
if len(TS_GRB_AFTER)>0:   print 'Max TS after...........: %.1f (%.1e)' %(TS_GRB_AFTER.max(),p_value(TS_GRB_AFTER.max()))

plt.figure(figsize=(10,7),facecolor='w')

plt.hist(TS_GRB_BEFORE,bins=X,normed=normed, histtype='step', cumulative=-1,label='All Before',color='b')
plt.hist(TS_GRB_AFTER,bins=X,normed=normed, histtype='step', cumulative=-1,label='All After',color='r')

plt.yscale('log')
plt.legend(fontsize=20,frameon=False)     
plt.xlabel('Test Statistic (TS)',size=20)
plt.ylabel('Number of events',size=20)
plt.ylim([0.5,2*len(TS_GRB)])

if normed: 
    plt.ylim([5e-5,2])
    plt.ylabel('Post-trial probability',size=20)
plt.xticks(size=20)
plt.yticks(size=20)

plt.figure(figsize=(10,7),facecolor='w')

plt.hist(TS_GRB_TOT,bins=X,normed=normed, histtype='step', cumulative=-1,label='All Events',color='b')
plt.hist(TS_GRB,bins=X,normed=normed, histtype='step', cumulative=-1,label='Slected Events',color='r',lw=2)

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
plt.xlabel('Test Statistic (TS)',size=20)
plt.ylabel('Number of events',size=20)
plt.ylim([0.5,2*len(TS_GRB_TOT)])
if normed: 
    plt.ylim([5e-5,2])
    plt.ylabel('Post-trial probability',size=20)
plt.xticks(size=20)
plt.yticks(size=20)

plt.show();
