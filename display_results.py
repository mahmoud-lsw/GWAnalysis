#!/usr/bin/env python
#from mpl_toolkits.basemap import Basemap
from GWAnalysis.Analysis import *
from matplotlib import rcParams
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl
import numpy as np
import healpy as hp
import sys
myFile=sys.argv[1]
pixel_file=sys.argv[2]
pixels=[]
for l in file(pixel_file,'r'):
    v = l.split()
    pixels.append(int(v[0]))
    NSIDE = int(v[1])
    pass
print NSIDE

#'GW151226_results.txt'
ra=[]
dec=[]
flux=[]
ts_tsmap = []
ts_src   = []
nprob    = []

for l in file(myFile,'r').readlines():
    v=l.split()
    if len(v)==0: continue
    ra.append(float(v[2]))
    dec.append(float(v[3]))
    flux.append(float(v[4]))
    ts_src.append(float(v[6]))
    nprob.append(int(v[7]))
    ts_tsmap.append(float(v[8]))
    print v[0],v[1],ra[-1],dec[-1],flux[-1],ts_src[-1],nprob[-1],ts_tsmap[-1]    
    pass
N    = len(ra) 
ra   = np.array(ra)
dec  = np.array(dec)
flux = np.array(flux)
ts_tsmap = np.array(ts_tsmap)
ts_src   = np.array(ts_src)
nprob    = np.array(nprob)

#idx=sky_to_pix(ra,dec, NSIDE)
NPIXELS=hp.nside2npix(NSIDE)
flux_map     = np.zeros(NPIXELS)
ts_tsmap_map = np.zeros(NPIXELS)
ts_src_map   = np.zeros(NPIXELS)
nprob_map    = np.zeros(NPIXELS)

for i,p in enumerate(pixels):
    flux_map[p]=flux[i]
    ts_tsmap_map[p]=ts_tsmap[i]
    ts_src_map[p]=ts_src[i]
    nprob_map[p]=nprob[i]
    pass
fig=plt.figure(figsize=(15,7),facecolor='w')
#fig.add_subplot(2,2,1)
from matplotlib import cm
cool_cmap = cm.jet
cool_cmap.set_under("w") # sets background to white
fmin=flux_map[flux_map>0].min()
fmax=flux_map[flux_map>0].max()
hp.mollview(flux_map, sub=111, title="Flux Upper Limit", unit="ph/cm$^{2}$/s", format='%.2g', norm='log',min=fmin,max=fmax,cmap=cool_cmap,notext=False)
hp.graticule()
#lat=90
##for lon in range(-180,180,30):
#    txt=lon
#    if lon<0: txt=lon+360
#    hp.projtext(sp.radians(lat),sp.radians(lon),'%d' %(txt))
#    pass
#lon=180
#for lat in range(0,180,30):
#    txt=lat
#    if lat>90: txt=lat-180
#    hp.projtext(sp.radians(lat),sp.radians(lon),'%d' %(txt))
#    pass

lat=0
for lon in range(-150,180,30):
    hp.projtext(lon,lat,'%d' %(lon),lonlat=True,size=15,va='bottom')
    pass
#lon=180
#for lat in range(-75,90,15):
#    hp.projtext(lon,lat,'%d' %(lat),lonlat=True,size=20)
#    horizontalalignment='right'
#    pass
lon=179.9
for lat in range(-60,90,30):
    if lat==0:va='center'
    elif lat>0:va='bottom'
    else: va='top'

    hp.projtext(lon,lat,'%d' %(lat),lonlat=True,size=15,horizontalalignment='right',va=va)
    pass
plt.text(-2.2,0,'Dec',rotation=90,size=20)
plt.text(0,-1.1,'RA',size=20,ha='center')

fig=plt.figure(figsize=(20,10),facecolor='w')
#fig.add_subplot(2,2,2)
hp.mollview(ts_tsmap_map, sub=221, title="TS max",min=0,max=25)
hp.graticule()
#fig.add_subplot(2,2,3)
hp.mollview(ts_src_map, sub=223, title="TS src",min=0,max=25)
hp.graticule()
#fig.add_subplot(2,2,4)
hp.mollview(nprob_map, sub=224, title="Nevt")
hp.graticule()

plt.show()
