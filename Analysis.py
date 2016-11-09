#!/usr/bin/env python
import sys
import json
import scipy as sp
import numpy as np
import healpy as hp
from matplotlib import pyplot as plt
import matplotlib.colors as colors
from astropy.io import fits
from astropy import coordinates as coord
from astropy import units as u
from glob import glob
import aplpy
from math import fabs
import ROOT
from ProgressBar import progressbar
import matplotlib.cm as cmx
import matplotlib as mpl
import numpy as np
from FT2 import FT2

from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma', 'Bitstream Vera Sans',
                                 'Lucida Grande', 'Verdana']
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


###
NMINEVENTS=2
TSMIN=20
###
def angsep(x1,y1,x2,y2):
    #    """
    #    spherical angle separation, aka haversine formula  
    #    input and output are in degrees               
    #    """
    dlat = np.math.radians(y2 - y1)
    dlon = np.math.radians(x2 - x1)
    y1 = np.math.radians(y1)
    y2 = np.math.radians(y2)
    a = np.sin(dlat/2.)*np.sin(dlat/2.) + np.cos(y1)*np.cos(y2)*np.sin(dlon/2.)*np.sin(dlon/2.)
    c  = 2*np.arctan2(np.sqrt(a), np.sqrt(1.-a))
    return np.math.degrees(c)

def EntryExit(ft2file,ra,dec,t0,theta_max,zenith_max,roi):
    #array of points:
    points = coord.SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    #array of reference points:
    referenceTimes=np.ones_like(points)*t0
    #get the FT2 file:
    hdulist=fits.open(ft2file)    
    SC_data=hdulist['SC_DATA'].data
    #array of times:
    time=SC_data.field('START')-t0
    time_selected=np.logical_and(time>-100, time<100)
    # SPACECRAFT:                                                                                                                                                                                                                                                           
    ra_zenith     = SC_data.field('RA_ZENITH')[time_selected]
    dec_zenith    = SC_data.field('DEC_ZENITH')[time_selected]
    #zenith        = coord.SkyCoord(ra=ra_zenith, dec=dec_zenith, unit="deg", frame='icrs')
    
    zenith_angles = angsep(x1,y1,x2,y2)
    
    print ra.shape, dec.shape, ra_zenith.shape, dec_zenith.shape,zenith.shape,zenith_angles.shape
    #zenith      = angsep(ra,dec,ra_zenith,dec_zenith)
    print zenith.shape
    # array of zeniths:
    #zenith      = coord.SkyCoord(ra=ra_zenith*u.degree, dec=dec_zenith*u.degree, frame='icrs')
    #print time_selected.shape,ra_zenith.shape,dec_zenith.shape
    
    # BORESIGHT:                                                                                                                                                                                                                                                            
    ra_scz  = SC_data.field('RA_SCZ')[time_selected]
    dec_scz = SC_data.field('DEC_SCZ')[time_selected]
    #array of z-axis:
    #scz     = np.array(coord.SkyCoord(ra=ra_scz*u.degree, dec=dec_scz*u.degree, frame='icrs'))
    
    #array of entry times:
    hdulist.close()
    in_fov = np.logical_and(scz<theta_max,zenith<zenith_max-roi)
    time_in_fov=time[in_fov]-referenceTimes
    print 'scz:',scz.shape
    print 'zenith:',zenith.shape
    print 'time:',time.shape
    print (scz<theta_max)
    print time.shape, in_fov, time_in_fov.shape
    #entry_times = sp.maximum(time[not in_fov]-referenceTimes<0))
    #exit_times  = sp.maximum(time[in_fov]-referenceTimes<0))
    #print time[time[in_fov]
    #entry_times=np.logical_and(in_fov,time>(scz<theta_max,zenith<zenith_max-roi)
    #return time, inSAA, scz.separation(point).degree,  zenith.separation(point).degree


    


def AngSeparation(ft2file,ra,dec):
    point = coord.SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    hdulist=fits.open(ft2file)    
    SC_data=hdulist['SC_DATA'].data
    # SPACECRAFT:                                                                                                                                                                                                                                                           
    ra_zenith   = SC_data.field('RA_ZENITH')
    dec_zenith  = SC_data.field('DEC_ZENITH')
    zenith      = coord.SkyCoord(ra=ra_zenith*u.degree, dec=dec_zenith*u.degree, frame='icrs')
    # BORESIGHT:                                                                                                                                                                                                                                                            
    ra_scz  = SC_data.field('RA_SCZ')
    dec_scz = SC_data.field('DEC_SCZ')
    scz     = coord.SkyCoord(ra=ra_scz*u.degree, dec=dec_scz*u.degree, frame='icrs')
    
    time=SC_data.field('START')    
    inSAA = SC_data.field('IN_SAA')
    
    hdulist.close()
    return time, inSAA, scz.separation(point).degree,  zenith.separation(point).degree

def AngSeparation_fast(point,scz,zenith):
    return scz.separation(point).degree,  zenith.separation(point).degree

def AngSeparation_ufast(ra,dec,scz_ra,scz_dec,zenith_ra,zenith_dec):
    return angsep(ra,dec,scz_ra,scz_dec),angsep(ra,dec,zenith_ra,zenith_dec)

def entry_exit(ra,dec,scz_ra,scz_dec,zenith_ra,zenith_dec,time,theta_max,zenith_max):

    def infov(theta_0,zenith_0):
        return theta_0<theta_max and zenith_0<zenith_max
    
    idx      =  sp.argmax(time[time<0])
    theta_0, zenith_0 =  AngSeparation_ufast(ra,dec,scz_ra[idx],scz_dec[idx],zenith_ra[idx],zenith_dec[idx])
    
    # IF IN THE FOV:
    if infov(theta_0, zenith_0):
        while infov(theta_0, zenith_0):
            idx-=1
            theta_0, zenith_0 =  AngSeparation_ufast(ra,dec,scz_ra[idx],scz_dec[idx],zenith_ra[idx],zenith_dec[idx])
            pass
        t0=time[idx]
        idx+=1            
    else:
        while not infov(theta_0, zenith_0):
            idx+=1
            theta_0, zenith_0 =  AngSeparation_ufast(ra,dec,scz_ra[idx],scz_dec[idx],zenith_ra[idx],zenith_dec[idx])
            pass
        t0=time[idx]
        pass
    # Now it is in the fov
    while infov(theta_0, zenith_0):
        idx+=1
        theta_0, zenith_0 =  AngSeparation_ufast(ra,dec,scz_ra[idx],scz_dec[idx],zenith_ra[idx],zenith_dec[idx])
        pass
    t1=time[idx]
    return t0,t1

def pix_to_sky(idx, nside):
    """Convert the pixels corresponding to the input indexes to sky coordinates (RA, Dec)"""
    
    theta, phi = hp.pix2ang(nside, idx)
    
    ra = np.rad2deg(phi)
    dec = np.rad2deg(0.5 * np.pi - theta)
    
    return ra, dec

def sky_to_pix(ra,dec, nside):
    """Convert the pixels corresponding to the input indexes to sky coordinates (RA, Dec)"""
    phi   = np.deg2rad(ra)
    theta = 0.5 * np.pi-np.deg2rad(dec)
    #dec = np.rad2deg(0.5 * np.pi - theta)    
    return hp.ang2pix(nside, theta, phi)

def readResult(resultFile):
    lines=file(resultFile,'r').readlines()
    myDictionary={}
    for l in lines:
        if '#' in l: continue
        p,v=l.split('=')
        #print p,v
        try:    myDictionary[p.strip()]=float(v.strip())
        except: myDictionary[p.strip()]=v.strip()
        pass
    return myDictionary

def addFGLtomap(map_fits,fgl_fits='/Users/omodei/Documents/DATA/CATALOGS/3FGL/gll_psc_v16.fit', radius=None):
    map_fitsfile=fits.open(map_fits)
    CRVAL1 = map_fitsfile[0].header['CRVAL1']
    CDELT1 = fabs(map_fitsfile[0].header['CDELT1'])
    NAXIS1 = map_fitsfile[0].header['NAXIS1']
        
    CRVAL2 = map_fitsfile[0].header['CRVAL2']
    CDELT2 = fabs(map_fitsfile[0].header['CDELT2'])
    NAXIS2 = map_fitsfile[0].header['NAXIS2']
    if radius is None: radius=(CDELT1*(NAXIS1/2.))
    print 'Adding 3FGL... RADIUS: ',radius
    fglfitsfile=fits.open(fgl_fits)
    fglfits=fglfitsfile[1]
    fgl_signifAvg=fglfits.data.field('Signif_Avg')
    fgl_flux=fglfits.data.field('Flux100_300')+fglfits.data.field('Flux300_1000')
    #for k in fglfits.header.keys(): print k,fglfits.header[k]
    fgl_ra=fglfits.data.field('RAJ2000')
    fgl_dec=fglfits.data.field('DEJ2000')
    fgl_name=fglfits.data.field('Source_Name')
    fgl_assoc=fglfits.data.field('ASSOC1')
    fgl_class=fglfits.data.field('CLASS1')
    fgl_skycoord = coord.SkyCoord(ra=fgl_ra*u.degree,dec=fgl_dec*u.degree)
    mapcenter    = coord.SkyCoord(ra=CRVAL1*u.degree,dec=CRVAL2*u.degree)
    mask_fgl     = sp.logical_and((fgl_signifAvg>20),fgl_skycoord.separation(mapcenter).degree < radius)
    #factor=1.0
    #mask_fgl=(fgl_ra>CRVAL1-factor*DX)*(fgl_ra<CRVAL1+factor*DX)*(fgl_dec>CRVAL2-factor*DY)*(fgl_dec<CRVAL2+factor*DY)
    fglfitsfile.close()
    X=fgl_ra[mask_fgl]
    Y=fgl_dec[mask_fgl]
    T=fgl_name[mask_fgl]             
    A=fgl_assoc[mask_fgl] 
    C=fgl_class[mask_fgl]
    #if print_sources:
    #    for i in range(len(X)):
    #        print X[i],Y[i],T[i],A[i],C[i]
    #        pass
    #    pass
    
    return X,Y,T,A,C

class LAT():
    def __init__(self,myDir,grbn,T0=0):
        self.gc=None
        self.T0=T0
        #print myDir,grbn
        #rootfile=glob('%s/%s_gtsrcprob_LIKE_UP_*root' % (myDir,grbn))[0]
        #print rootfile
        #self.eventsROOTfile = ROOT.TFile(rootfile)
        #EVENTS=self.eventsROOTfile.Get('Events')
        #NEVENTS    = EVENTS.GetEntries()
        #NEVENTS_90 = EVENTS.GetEntries('SRC_PROBABILITY>0.9')
        #print 'Nevents=',NEVENTS, 'P>0.9=', NEVENTS_90
        #self.eventsROOTfile.Close()
        print '=====>',myDir,grbn
        self.results=readResult(glob('%s/results_%s.txt' % (myDir,grbn))[0])
        self.eventsfile = glob('%s/%s_events_LIKE_MY_*fits' % (myDir,grbn))[0]
        self.tsmapfile  = glob('%s/%s_LAT_tsmap_LIKE_MY.fits' % (myDir,grbn))[0]
        
        tsmp_fitsfile=fits.open(self.tsmapfile)
        self.TSMAX   =tsmp_fitsfile[0].data.max()
        
        self.TS_GRB  = self.results['LIKE_UP_TS_GRB']
        self.NEVT_GRB  = self.results['LIKE_UP_gtsrcprob_Nthr']

        
        self.date_obs=tsmp_fitsfile[0].header['DATE-OBS'].split('.')[0]
        self.date_end=tsmp_fitsfile[0].header['DATE-END'].split('.')[0]
        
        self.TSTART=tsmp_fitsfile[0].header['TSTART']
        self.TSTOP=tsmp_fitsfile[0].header['TSTOP']
        
        self.CRVAL1 = tsmp_fitsfile[0].header['CRVAL1']
        self.CDELT1 = fabs(tsmp_fitsfile[0].header['CDELT1'])
        self.NAXIS1 = tsmp_fitsfile[0].header['NAXIS1']

        self.CRVAL2 = tsmp_fitsfile[0].header['CRVAL2']
        self.CDELT2 = fabs(tsmp_fitsfile[0].header['CDELT2'])
        self.NAXIS2 = tsmp_fitsfile[0].header['NAXIS2']

        tsmp_fitsfile.close()

        #self.mystring='R.A.=%.2f, Dec.=%.2f, %s -- %s (%d), TSmax=%.1f, EMAX=%.1f MeV\n' %(self.CRVAL1,self.CRVAL2, self.date_obs,self.date_end,self.TSTART-T0,self.TSMAX,self.ene_max)
        self.mystring='|%6.1f,%6.1f | %s -- %s (%8d, %8.1f, [%8.1f]) |  %.0f/%.0f/%d  |\n' %(self.CRVAL1,self.CRVAL2, self.date_obs,self.date_end,(self.TSTART-T0)/86400,(self.TSTOP-T0)/86400,
                                                                                                       self.TSTOP-self.TSTART,self.TSMAX,self.TS_GRB,self.NEVT_GRB)
        if self.TS_GRB>TSMIN and self.NEVT_GRB>NMINEVENTS:
            self.mystring='|%6.1f,%6.1f | %s -- %s (%8d, %8.1f, [%8.1f]) | *%.0f/%.0f/%d* |\n' %(self.CRVAL1,self.CRVAL2, self.date_obs,self.date_end,(self.TSTART-T0)/86400,(self.TSTOP-T0)/86400,
                                                                                                           self.TSTOP-self.TSTART,self.TSMAX,self.TS_GRB,self.NEVT_GRB)


    def __del__(self):
        if self.gc is not None: self.gc.close()
        pass
    def is_inside(self,ra,dec):
        Px,Py=self.gc.world2pixel(ra, dec)
        return sp.logical_and(sp.logical_and(Px>0,Px<=self.NAXIS1),sp.logical_and(Py>0,Py<=self.NAXIS2))
        #W0X,W0Y=self.gc.pixel2world(1, 1)
        #self.gc.show_markers(W0X,W0Y, edgecolor='m', linewidths=3, facecolor='none',marker='o', s=10, alpha=1.0)
        #W1X,W1Y=self.gc.pixel2world(self.NAXIS1, self.NAXIS2)
        #self.gc.show_markers(W1X,W1Y, edgecolor='y', linewidths=3, facecolor='none',marker='o', s=10, alpha=1.0)
        #return False
    
    def draw(self):
        map_x=10
        map_y=7
        print '***** Drawing TS map...',self.tsmapfile
        self.gc = aplpy.FITSFigure(self.tsmapfile,figsize=(map_x,map_y),facecolor='w')
        self.gc.show_colorscale(vmin=0, vmax=25)
        self.gc.add_colorbar()
        #self.gc.colorbar(cax, ticks=[-1, 0, 1])
        #cbar.ax.set_yticklabels(['< -1', '0', '> 1']) 
        # DRAW EVENTS:
        events_file=fits.open(self.eventsfile)
        self.fits_file=events_file[1].data    
        self.ra=self.fits_file.field('RA')
        self.dec=self.fits_file.field('DEC')
        self.energy=self.fits_file.field('ENERGY')
        self.ene_max  = 0
        self.ra_max   = 0
        self.dec_max  = 02
        self.mask=self.is_inside(self.ra,self.dec)
        #for i,e in enumerate self.energy:
        if len(self.energy[self.mask])>0:
            ene_id   = sp.argmax(self.energy[self.mask])
            self.ene_max  = self.energy[self.mask][ene_id]
            self.ra_max   = self.ra[self.mask][ene_id]
            self.dec_max  = self.dec[self.mask][ene_id]
            pass
        
        self.gc.show_markers(self.ra,self.dec, edgecolor='w', linewidths=3, facecolor='none',marker='o', s=self.energy/10.0, alpha=1.0)
        #self.gc.add_label(0.5,0.9,'%s - %s (TSmax=%.1f)' %(self.date_obs,self.date_end,self.TSMAX),relative=True,size=15,color='w')
        self.gc.add_label(0.5,0.9,'%s - %s' %(self.date_obs,self.date_end),relative=True,size=15,color='w')
        #self.gc.add_label(0.3,0.95,'R.A.=%.2f Dec.=%.2f' %(self.CRVAL1,self.CRVAL2),relative=True,size=15,color='w')

        #self.gc.add_label(0.5,1.02,'R.A.=%.2f Dec.=%.2f, $\\Delta$T (d)=%.1f (Nevt=%d, TS$_{max}$=%.1f, TS$_{src}$=%.1f)' %(self.CRVAL1,
        #                                                                                                                    self.CRVAL2,
        #                                                                                                                    (self.TSTART-self.T0)/86400.,
        #                                                                                                                    self.NEVT_GRB,self.TSMAX,self.TS_GRB),
        #                                                                                                                    relative=True,size=15,color='k')

        if self.T0:
            self.gc.add_label(0.5,1.02,'R.A.,Dec.=%.1f,%.1f, [%.1f days] (N$_{\gamma}$=%d, TS$_{max}$=%.0f, TS$_{src}$=%.0f)' %(self.CRVAL1,
                                                                                                                                self.CRVAL2,
                                                                                                                                (self.TSTART-self.T0)/86400.,
                                                                                                                                self.NEVT_GRB,self.TSMAX,self.TS_GRB),
                                                                                                                                relative=True,size=15,color='k')
        else:
            self.gc.add_label(0.5,1.02,'R.A.,Dec.=%.1f,%.1f (N$_{\gamma}$=%d, TS$_{max}$=%.0f, TS$_{src}$=%.0f)' %(self.CRVAL1,self.CRVAL2,self.NEVT_GRB,self.TSMAX,self.TS_GRB),relative=True,size=15,color='k')
        
        
        #if self.ene_max>0:  self.gc.add_label(self.ra_max,self.dec_max,'%d MeV' %(self.ene_max),size=20,color='w')
        self.gc.add_grid()
        self.gc.grid.set_color('white')
        self.gc.grid.set_alpha(0.5)
        self.gc.grid.set_xspacing(5.0) 
        self.gc.grid.set_yspacing(5.0)
        self.gc.tick_labels.set_xformat('dd')
        self.gc.tick_labels.set_yformat('dd')
        self.gc.tick_labels.set_font(size=20)
        self.gc.axis_labels.set_font(size=20)
        self.gc.colorbar.set_axis_label_text('TS')
        self.gc.colorbar.set_axis_label_font(size=20)
        self.gc.colorbar.set_font(size=20)
        self.gc.colorbar.set_ticks([0,5,10,15,20,25])
                             #print dir(self.gc.tick_labels)
        #self.gc.tick_labels.set_xfontsize(20)
        
        #for xw,yw in zip(self.ra,self.dec):
        #    print '----->',xw,yw,self.gc.world2pixel(xw, yw),
        events_file.close()
        return self.gc
    
    def addContour(self,D,perc):
        X,Y=plotContour(D,perc)
        self.gc.show_markers(X,Y, edgecolor='yellow', facecolor='none',marker='o', s=2, alpha=1.0)

    def addFakeContour(self):
        #plt.plot(arc_ra1,arc_dec1,'-',label='Simulated Contour 90%')
        #plt.plot(arc_ra2,arc_dec2,'-',label='Simulated Contour 50%')

        self.gc.show_markers(arc_ra1,arc_dec1, edgecolor='yellow', facecolor='none',marker='o', linewidths=1,s=2, alpha=1.0)
        self.gc.show_markers(arc_ra2,arc_dec2, edgecolor='yellow', facecolor='none',marker='o', linewidths=1,s=2, alpha=1.0)



    def addFGL(self,radius=None,print_sources=False):
        fgl_fits='/Users/omodei/Documents/DATA/CATALOGS/3FGL/gll_psc_v16.fit'
        X,Y,T,A,C = addFGLtomap(self.tsmapfile,fgl_fits, radius)
        if len(X)==0: return 
        FGL_mask=self.is_inside(X,Y)
        self.gc.show_markers(X[FGL_mask],Y[FGL_mask], edgecolor='yellow', facecolor='none',marker='*', s=100, linewidth=2,alpha=1.0)
        for i in range(len(A[FGL_mask])):
            fgl_ra=X[FGL_mask][i]
            fgl_dec=Y[FGL_mask][i]
            fgl_name=T[FGL_mask][i]
            fgl_ass=A[FGL_mask][i]
            mylabel=fgl_ass
            if print_sources:  print '===>',fgl_ra,fgl_dec,fgl_name,fgl_ass
            #self.gc.add_label(fgl_ra,fgl_dec,mylabel,size=10,color='yellow')
            pass
        
        #for i in range(len(X)):
        #    if len(A[i])>1: mylabel='%s' % (A[i])
        #    else: mylabel=T[i]
        #    self.gc.add_label(X[i],Y[i]+0.2,'%s' %(mylabel),size=10,color='yellow')
    

    def save(self):
        filename='tsmap_%d_%d_%d.png' %(self.CRVAL1,self.CRVAL2,self.TSTART-self.T0)
        print filename
        self.gc.save(filename)
        
def plotContour(D,percentile):
    X=[]
    Y=[]
    for contourId in range(len(D['contours'])):
        if ('%d' %percentile) in D['contours'][contourId]['name']:
            npoints=len(D['contours'][contourId]['coords'])
            for i in range(npoints):
                X.append(D['contours'][contourId]['coords'][i][0])
                Y.append(D['contours'][contourId]['coords'][i][1])
                #print X[i],Y[i]
                pass
            pass
        pass
    X=sp.array(X)
    Y=sp.array(Y)
    return X,Y

    
def plotCircle(RA,DEC,radius,npoints=100):
    theta=sp.linspace(0,2*sp.pi,npoints)
    x = np.remainder(RA+360+radius*sp.cos(theta),360) # shift RA values
    
    y = DEC+radius*sp.sin(theta)
    ind=y>90
    y[ind]-=180
    return x,y

def plotGrid():
    ras=sp.linspace(0,360,12)
    decs=sp.linspace(-90,90,12)
    gx=[]
    gy=[]
    for ra in ras:
        for y in sp.linspace(-90,90,180):
            gx.append(ra)
            gy.append(y)
            pass
        pass
    for dec in decs:
        for x in sp.linspace(0,360,360):
            gx.append(x)
            gy.append(dec)
            pass
    return gx,gy

def fakeContour(l0,b0,dl,db,npoints=100):
    
    theta=sp.linspace(0,sp.pi,npoints)
    arc_l=dl*sp.cos(theta)
    arc_b=db*sp.sin(theta)
    #plt.plot(l0+arc_l,b0+arc_b)
    #plt.plot(l0+arc_l,b0-arc_b)
    arc_ra=[]
    arc_dec=[] 
    for l,b in zip(l0+arc_l,b0+arc_b):
        sc=coord.SkyCoord(l=l*u.degree, b=b*u.degree, frame='galactic')
        arc_ra.append(sc.fk5.ra.degree)
        arc_dec.append(sc.fk5.dec.degree)
        pass
    if db>0:
        for l,b in zip(l0+arc_l,b0-arc_b):
            sc=coord.SkyCoord(l=l*u.degree, b=b*u.degree, frame='galactic')
            arc_ra.append(sc.fk5.ra.degree)
            arc_dec.append(sc.fk5.dec.degree)
            pass
        pass
    return arc_ra,arc_dec

def plot_mwd(RA,Dec,ax,orgx=0,orgy=0,color='k',area=1,alpha=1.0,cmap='default',label=''):
    ''' RA, Dec are arrays of the same length.
    RA takes values in [0,360), Dec in [-90,90],
    which represent angles in degrees.
    org is the origin of the plot, 0 or a multiple of 30 degrees in [0,360).
    title is the title of the figure.
    projection is the kind of projection: 'mollweide', 'aitoff', 'hammer', 'lambert'
    '''
    x = np.remainder(RA+360-orgx,360) # shift RA values
    ind = x>180
    x[ind] -=360    # scale conversion to [-180, 180]
    x=-x    # reverse the scale: East to the left

    y = sp.remainder(Dec+270-orgy,180)-90

    x_tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    x_tick_labels = np.remainder(x_tick_labels+360+orgx,360)

    y_tick_labels = np.array([-75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75])
    y_tick_labels = np.remainder(y_tick_labels+270+orgy,180)-90
    #tick_labels = np.array(['10$^{h}$', '8$^{h}$', '6$^{h}$', '4$^{h}$', '2$^{h}$', '0$^{h}$', '-2$^{h}$', '-4$^{h}$', '-6$^{h}$', '-8$^{h}$', '-10$^{h}$'])
    
    sc=ax.scatter(np.radians(x),np.radians(y), marker='.',c=color,s=area,edgecolor='face',label=label,alpha=alpha,cmap=cmap)  # convert degrees to radians
    ax.set_xticklabels(x_tick_labels)     # we add the scale on the x axis
    ax.set_yticklabels(y_tick_labels)     # we add the scale on the x axis
    ax.set_xlabel("RA")
    ax.xaxis.label.set_fontsize(22)
    ax.set_ylabel("Dec")
    ax.yaxis.label.set_fontsize(22)
    ax.grid(True)
    return sc 
##################
def MakeDistribution(mydir,T0,Contour,draw=1):
    mystring='| RA, DEC | DATE (second from T0) | TSMAX | EMAX |\n'

    TS_GRB   = []
    TS_TSMAP = []
    TSTART   = []
    TSTOP    = []
    NEVT_GRB = []
    
    print '============================================================'
    GRBs={}
    files=glob('%s/*_LAT_tsmap_LIKE_MY*'% mydir)
    for f in files:
        grb=f.replace(mydir,'').split('_LAT_tsmap_LIKE_MY')[0]
        GRBs[grb]=grb
        pass
    k=GRBs.keys()
    k.sort()
    #k=k[:10]
    pb=progressbar(20)
    N=len(k)
    for i,g in enumerate(k):
        pb.go(i,N)
        try:
            lat=LAT(mydir,g,T0=T0)
        except:
            continue
        ts_grb   = lat.TS_GRB
        ts_tsmap = lat.TSMAX
        tstart   = lat.TSTART
        tstop    = lat.TSTOP
        nevt_grb = lat.NEVT_GRB

        TS_GRB.append(lat.TS_GRB)
        TS_TSMAP.append(ts_tsmap)
        TSTART.append(tstart)
        TSTOP.append(tstop)
        NEVT_GRB.append(nevt_grb)
        
        if (ts_grb>TSMIN and nevt_grb>NMINEVENTS and draw) or (N<5 and draw):
            lat.draw()
            #lat.addFakeContour()
            lat.addContour(Contour,50)
            lat.addContour(Contour,90)
            lat.addFGL()
            lat.save()
            pass        
        mystring+=lat.mystring
        del(lat)
        pass
    print mystring
    return sp.array(TSTART),sp.array(TSTOP),sp.array(TS_GRB),sp.array(TS_TSMAP),sp.array(NEVT_GRB)
##################################################
##################################################
def analyze(json_skymap,
            healpix_map,
            ft2file,
            triggertime,
            name='GW',
            orgx=0,
            orgy=0,
            projection= 'hammer',
            map_x       = 10,
            map_y       = 7,
            perc        = 0.1,
            tmin        = -864000,
            tmax        = 864000,
            theta_max   = 70,
            zenith_max  = 100,
            degrade     = 32,
            roi=None,
            #degrade     = None
            irfs='P8_TRANSIENT010E',
            split_orbits=False,
            adaptive_time=True,tsmap=False):
    
    ligo_map   = hp.read_map(healpix_map)
    if degrade is not None: ligo_map = hp.ud_grade(ligo_map,degrade)

    NPIX=hp.get_map_size(ligo_map)
    NSIDE=hp.npix2nside(NPIX)
    RESOLUTION=hp.nside2resol(NSIDE,arcmin=True)/60.
    print 'N PIXEL      = ',NPIX
    print 'NSIDE        = ',NSIDE
    print 'RESOLUTION   = ',RESOLUTION
    
    if roi is None: masked_radius = RESOLUTION
    else: masked_radius=roi

    print len(ligo_map)
    pixels = np.arange(NPIX)
    #print ligo_map*([)ligo_map>1e-5)
    p_max=ligo_map.max()
    p_min=ligo_map.min()

    p_selected=(perc)*p_max
    ones      = np.ones(len(ligo_map))
    mask      = (ligo_map>p_selected)
    ligo_selected  = ones*(mask)
    ligo_entry=ones*(mask)
    ligo_expo=ones*(mask)
    pixels_selected = pixels[mask]
    
    #ra,dec = pix_to_sky(pixels,NSIDE)

    masked_ra,masked_dec = pix_to_sky(pixels_selected,NSIDE)
    #masked_ra_dec= coord.SkyCoord(ra=masked_ra*u.degree, dec=masked_dec*u.degree, frame='icrs')
    
    Nselected=len(pixels_selected)
    print 'Number of selected pixels:',Nselected
    myFT2=FT2(ft2file,triggertime-100000,triggertime+100000)
    myFT2.fov(theta_max,zenith_max-masked_radius)
    #EntryExit(ft2file,masked_ra,masked_dec,triggertime,theta_max,zenith_max,masked_radius)
    #hdulist=fits.open(ft2file)    
    #SC_data=hdulist['SC_DATA'].data

    #time_all = SC_data.field('START')-triggertime
    #time_sel = sp.logical_and(time_all>-10000,time_all<10000)
    #time = time_all[time_sel]
    #inSAA = SC_data.field('IN_SAA')[time_sel]

    # SPACECRAFT:                                                                                                                                                                                                                                                           
    #ra_zenith      = SC_data.field('RA_ZENITH')[time_sel]
    #dec_zenith     = SC_data.field('DEC_ZENITH')[time_sel]
    #zenith_ra_dec  = coord.SkyCoord(ra=ra_zenith*u.degree, dec=dec_zenith*u.degree, frame='icrs')
    # BORESIGHT:                                                                                                                                                                                                                                                            
    #ra_scz         = SC_data.field('RA_SCZ')[time_sel]
    #dec_scz        = SC_data.field('DEC_SCZ')[time_sel]
    #scz_ra_dec     = coord.SkyCoord(ra=ra_scz*u.degree, dec=dec_scz*u.degree, frame='icrs')
    #print 'zenith...:',zenith_ra_dec.shape
    #print 'scz......:',scz_ra_dec.shape    
    #hdulist.close()
    ######
    times_t0=[]
    times_t1=[]
    if adaptive_time or split_orbits:
        #######################################
        for i in range(Nselected):
            #if i<100: continue
            #was_inFov=0
            #time_in=[]
            #time_out=[]
            entry_time,exit_time=myFT2.getEntryExitTime(masked_ra[i],masked_dec[i],triggertime) #entry_exit(masked_ra[i],masked_dec[i],ra_scz,dec_scz,ra_zenith,dec_zenith,time,theta_max,zenith_max-masked_radius)
            #for j,t in enumerate(time):
            #    if scz[j]<theta_max and zenith[j] < zenith_max-masked_radius:
            #        inFov=1
            #    else:
            #        inFov=0
            #        pass
            #    if was_inFov==0 and inFov==1:
            #        #print 'Enter the FoV at time.......:',t
            #        was_inFov=1
            #        time_in.append(t)
            #        pass
            #    if was_inFov==1 and inFov==0:
            #        #print 'Exit from the FoV at time...:',t 
            #        was_inFov=0
            #        if len(time_in)>0 and t>time_in[-1]:
            #            time_out.append(t)                    
            #            pass
            #        if time_out[-1]>0: break
            #        pass
            #    pass
                
            pix=pixels_selected[i]

            #entry_time=time_in[-1]
            #exit_time =time_out[-1]
            expo_time=exit_time-entry_time

            ligo_entry[pix]*=entry_time
            ligo_expo[pix]*=expo_time
            times_t0.append(entry_time)
            times_t1.append(exit_time)
            print 'i/N:%7d/%7d pix:%7d ra:%10.3f dec:%10.3f entry:%10.3f exit:%10.3f exposure:%10.3f' %(i,Nselected,pix,masked_ra[i],masked_dec[i],entry_time,exit_time,expo_time)
            #a=raw_input('')
            pass
        times_t0=np.array(times_t0)
        times_t1=np.array(times_t1)
        
        fig=plt.figure(figsize=(20,10),facecolor='w')
        idx = ligo_entry == 0
        ligo_entry[idx] = np.nan

        from matplotlib import cm
        jet = plt.get_cmap('jet') 
        #jet = colors.Colormap('jet')
        cNorm  = colors.Normalize(vmin=times_t0.min(), vmax=times_t0.max())
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

        cool_cmap = cm.jet
        cool_cmap.set_under("w") # sets background to white
        cool_cmap.set_bad("w") # sets background to white
        #cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=jet,norm=cNorm,orientation='horizontal')
        #cb1.set_label ('Time since LIGO trigger [s]')
        rot=(orgx,orgy)        
        hp.mollview(ligo_map/p_max, sub=221, title='Ligo Map',cmap=cool_cmap,rot=rot)
        hp.graticule()
        hp.mollview(ligo_selected, sub=222, title='Ligo Map ($>$%.1f)' % perc,cmap=cool_cmap,norm='lin',min=0,max=1.0,rot=rot)
        hp.graticule()
        hp.mollview(ligo_entry, sub=223, title='Entry Time',cmap=jet,norm='lin',min=times_t0.min(),max=times_t0.max(),rot=rot)
        hp.graticule()
        hp.mollview(ligo_expo, sub=224, title='Exposure (s)',cmap=cool_cmap,norm='log',min=0.1,max=max(ligo_expo),rot=rot)
        hp.graticule()
        ax = plt.gca()
    
        #fig = plt.figure(figsize=(15, 8))
        #hp.mollview(ligo_entry, title='Entry Time',cmap=jet,norm='lin',min=times_t0.min(),max=times_t0.max(),rot=rot)
        #hp.graticule()
        #ax = plt.gca()
        #cb1 = mpl.colorbar.ColorbarBase(ax, cmap=jet,norm=cNorm,orientation='horizontal',ticklocation='top')
        #cb1.set_label('Time since t$_{GW}$ [s]',size=15)
        #cb1.ax.tick_params(labelsize=20)


    else:
        times_t0=tmin*np.ones(Nselected)
        times_t1=tmax*np.ones(Nselected)
        pass
    
    file_out_name='pixels_%d.txt' % NSIDE
    print 'writing the output file...',file_out_name
    print len(pixels_selected)
    output_pixels='#pixel NSIDE'
    for i,p in enumerate(pixels_selected):
        output_pixels+='%d %d\n' %(p,NSIDE)
        pass
    fout = file(file_out_name,'w')
    fout.writelines(output_pixels)
    fout.close()

    
    exposure=times_t1-times_t0
    print 'Entry Time: MINIMUM: %.1f, MAXIMUM: %.1f ' %(times_t0.min(),times_t0.max())
    print 'Exit  Time: MINIMUM: %.1f, MAXIMUM: %.1f ' %(times_t1.min(),times_t1.max())
    print 'Exposure    MINIMUM: %.1f, MAXIMUM: %.1f ' %(exposure.min(),exposure.max())

    #matplotlib.colors.Colormap(name, N=256)

    jet = plt.get_cmap('jet') 
    #jet = colors.Colormap('jet')
    cNorm  = colors.Normalize(vmin=times_t0.min(), vmax=times_t0.max())
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

    #1 PLOT the JSON countours:
    JSON_D = json.load(file(json_skymap))
    fig = plt.figure(figsize=(1.5*map_x,map_x),facecolor='w')
    ax  = fig.add_subplot(111, projection=projection)
    
    RA,Dec = plotContour(JSON_D,int((1.-perc)*100))
    plot_mwd(RA,Dec,ax,color='k',orgx=orgx,orgy=orgy,label=name)
    txt='# ra dec flux flux_err tssrc\n'
    for i in range(Nselected):
        ra0  = masked_ra[i]
        dec0 = masked_dec[i]
        tmin=times_t0[i]
        tmax=times_t1[i]
        txt+='%f %f %f %f %f\n' % (ra0,dec0,tmin,0,tmax)
        #print times[i]
        colorVal = scalarMap.to_rgba(times_t0[i])
        if NSIDE<32:
            x0,y0=plotCircle(ra0,dec0,masked_radius,100)
            plot_mwd(x0,y0,ax,color=colorVal,orgx=orgx,orgy=orgy,label='')
            pass
        else:
            plot_mwd(np.array([ra0]),np.array([dec0]),ax,color=colorVal,orgx=orgx,orgy=orgy,label='')
            pass
        if adaptive_time:            
            outname='%s_N%d/%s_%s_ORB_ROI_%.1f_POS%03d' %(name,NSIDE,name,irfs,masked_radius,i)
        else:
            outname='%s_N%d/%s_%s_%d_%d_%.1f_POS%03d' %(name,NSIDE,name,irfs,tmin,tmax,masked_radius,i)
            pass
        options='-nolle'
        if split_orbits:  options='-nolle -orbit'
        if tsmap: options+=' -tsmap'
        cmd='./GWAnalysis.py -ra %.1f -dec %.1f -t0 %.3f -t90 %.1f -irfs %s -emin 100 -zmax %.1f -roi %.1f %s -out %s' %(ra0,dec0,triggertime+tmin,tmax-tmin,irfs,zenith_max, masked_radius,options,outname)
        print cmd
        #cm.show_markers(plotCircle(masked_ra[i],masked_dec[i],masked_radius), edgecolor=circle_color, facecolor='none',marker='o', s=1, alpha=1.0)
        #print i,masked_ra[i],masked_dec[i],masked_radius
        pass
    file_out_name='times_%d.txt' % NSIDE
    fout = file(file_out_name,'w')
    fout.writelines(txt)
    fout.close()
    
    ax2 = fig.add_axes([0.15, 0.05, 0.7, 0.05])

    cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=jet,norm=cNorm,orientation='horizontal')
    cb1.set_label('Time since LIGO trigger [s]')


    #cbar = fig.colorbar(ax, ticks=[times.min(),times.max()], orientation='horizontal')
    #cbar.ax.set_xticklabels(['Low', 'Medium', 'High'])  # horizontal colorbar



    plt.show()
    exit()
    

    
