import sys, os
import skill_metrics as sm
import glob
import datetime
import itertools
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.colorbar
import netCDF4
import numpy
import cartopy
import cmocean
import imageio
from pathlib import Path
from pygifsicle import optimize
import matplotlib
matplotlib.use('Agg')
from netCDF4 import Dataset, num2date, date2num
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import cartopy.crs as ccrs
import cartopy
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from cmocean import cm as cmo
import os
import yaml
defaultcm = cmo.thermal
from matplotlib import cm
import re


def stats(ax, area, anomaly, label):
    mn = (anomaly*area).sum()/area.sum()
    ame = (abs(anomaly)*area).sum()/area.sum()
    sd = numpy.sqrt( ((anomaly-mn)**2*area).sum()/area.sum() )
    rms = numpy.sqrt( (anomaly**2*area).sum()/area.sum() )
    qmn, qmx = anomaly.min(), anomaly.max()
    #print(label, 'mean =', mn, 'sd =', sd, 'rms =', rms, 'min =', qmn, 'max =', qmx, 'ame =', ame )
    bb = ax.get_position()
    plt.gcf().text(bb.x0,bb.y1+.01,'ame=%.3f$^\circ$C'%ame, horizontalalignment='left')
    plt.gcf().text(bb.x1,bb.y1+.01,'rms=%.3f$^\circ$C'%rms, horizontalalignment='right')

#G5 = netCDF4.Dataset('/work/noaa/da/ycteng/runs/convert_hat10/soca_gridspec_hat10.nc')
G5 = netCDF4.Dataset('../../../soca-diagnostics/soca_gridspec.nc')
M5 = netCDF4.Dataset('/work/noaa/marine/Cameron.Book/ostia_check/ocean_mask_hat10.nc')

xq5 = G5.variables['lon'][0,:]
yq5 = G5.variables['lat'][0,:]
a5 = G5.variables['area'][0,:]
m5 = M5.variables['mask'][:]; a5 = a5*m5

# HAT10 Region
extent = [-98, -10, 0, 45]
extent = [-80, -50, 10, 40]

central_lon = numpy.mean(extent[:2])
central_lat = numpy.mean(extent[2:])

wdir = '/work/noaa/stmp/lliu/rtfos_ostia/fcst/'
delta=datetime.timedelta(days=1)

start_date = datetime.date(2020, 7, 28)
end_date = datetime.date(2020, 8, 3) #datetime.date(2003,  3, 31)

bias_fcst=np.zeros([3,5])
rmsd_fcst=np.zeros([3,5])
sdev_fcst=np.zeros([3,5])
print('bias_fcst',bias_fcst)
w=0
x=numpy.arange(0,5) #number of fcst days
for expname in ['all','ctrl','noda']: 
 print('expname',expname)
 expdir=wdir+expname
 ini_date=start_date
 print('ini_date',ini_date)
 n=5
 m=2
 bias=np.zeros([m,n])
 rmsd=np.zeros([m,n])
 sdev=np.zeros([m,n])
 ccoef=np.zeros([m,n])
 i=0
 while ini_date<=datetime.date(2020,7,29):
  tmp_date=ini_date+delta
  c=0
  print('date check',tmp_date,ini_date+5*delta)
  while tmp_date <= ini_date+5*delta: #end_date:
    thisdate=datetime.datetime.strptime(str(tmp_date),"%Y-%m-%d").strftime("%Y%m%d")
    ncfile_ost=wdir+'/ostia/'+thisdate+'120000-UKMO-L4_GHRSST-SSTfnd-OSTIA-GLOB-v02.0-fv02.0-sst-hat10.nc'
    OM1 = netCDF4.Dataset(ncfile_ost,'r')
    ncfile=expdir+'/ufs.hafs.cpl.hi.'+str(thisdate)+'-'+str(ini_date)+'.nc'
    W5 = netCDF4.Dataset(ncfile,'r')

    fig = plt.figure(figsize=(12,8))
    vmin,vmax,ci,cmap = -2.25,2.25,.5,plt.cm.RdYlBu_r
    axes = []
    cilev = numpy.arange(vmin-ci,vmax+ci*2,ci)
    norm = matplotlib.colors.BoundaryNorm(boundaries=cilev, ncolors=cmap.N)

    ost=OM1.variables['sst_o'][0,:,:]*m5
    ocn=W5.variables['ocnImp_So_t'][0,:,:]*m5
    ocn = ocn-273.15
    ost=ost.flatten()
    ocn=np.squeeze(ocn)
    ocn=ocn.flatten()
# first filter out ost nan
    ocn = ocn[~np.isnan(ost)]
    ost = ost[~np.isnan(ost)]

# then filter large difference
    ost0 = ost[abs(ocn-ost)<50]
    ocn = ocn[abs(ocn-ost)<50]
     
    bias[i,c]=sm.bias(ocn,ost0) 
    rmsd[i,c]=sm.rmsd(ocn,ost0)
    sdev[i,c]=np.std(ocn)
    ccoef0 = np.corrcoef(ocn,ost0) 
    ccoef[i,c] = ccoef0[0,1]      
#    bias=abs(bias)
    #print('checks',i,bias)
    c=c+1
    tmp_date += delta
  i=i+1
  ini_date += delta
  print('i,ini_date',i,ini_date)

 #experiment loop
 print('w',w)
 #print('rmsd',rmsd) # different initial date 
 bias_fcst[w,:]=np.mean(ccoef,axis=0)
 rmsd_fcst[w,:]=np.mean(rmsd,axis=0)
 sdev_fcst[w,:]=np.mean(sdev,axis=0)
 w=w+1

plt.close('all')
print('rmsd_fcst,bias_fcst,sdev_fcst',rmsd_fcst,bias_fcst,sdev_fcst) # different experiment
fig, ax=plt.subplots(figsize=(10,3.5))
plt.title('RMSD (CDEPS fcst-OSTIA) (deg)',fontsize=13)
plt.plot(x,rmsd_fcst[0,:],'-*',color='g',markerfacecolor='w',markersize=4,label='all',linewidth=2,alpha=0.7) #exp 1
plt.plot(x,rmsd_fcst[1,:],'-*',color='r',markerfacecolor='w',markersize=4,label='ctrl',linewidth=2,alpha=0.7) #exp 2
plt.plot(x,rmsd_fcst[2,:],'-*',color='k',markerfacecolor='w',markersize=4,label='noda',linewidth=2,alpha=0.7)
plt.grid()
plt.legend()
plt.savefig('rmsd_fcst_ostia.png',bbox_inches=0)
plt.close(fig)

fig, ax=plt.subplots(figsize=(10,3.5))
plt.title('Correlation (CDEPS FCST-OSTIA) (deg)',fontsize=13)
plt.plot(x,bias_fcst[0,:],'-*',color='g',markerfacecolor='w',markersize=4,label='all',linewidth=2,alpha=0.7)
plt.plot(x,bias_fcst[1,:],'-*',color='r',markerfacecolor='w',markersize=4,label='ctrl',linewidth=2,alpha=0.7)
plt.plot(x,bias_fcst[2,:],'-*',color='k',markerfacecolor='w',markersize=4,label='noda',linewidth=2,alpha=0.7)
plt.grid()
plt.legend()
plt.savefig('corr_fcst_ostia.png',bbox_inches=0)

fig, ax=plt.subplots(figsize=(10,3.5))
plt.title('Standard deviation (CDEPS FCST) (deg)',fontsize=13)
plt.plot(x,sdev_fcst[0,:],'-*',color='g',markerfacecolor='w',markersize=4,label='all',linewidth=2,alpha=0.7)
plt.plot(x,sdev_fcst[1,:],'-*',color='r',markerfacecolor='w',markersize=4,label='ctrl',linewidth=2,alpha=0.7)
plt.plot(x,sdev_fcst[2,:],'-*',color='k',markerfacecolor='w',markersize=4,label='noda',linewidth=2,alpha=0.7)
plt.grid()
plt.legend()
plt.savefig('sdev_fcst_ostia.png',bbox_inches=0)
plt.close('all')

label = {'ALL-OSTIA': 'r', 'CTRL-OSTIA': 'g', 'NODA-OSTIA': 'b'}
#### Taylor diag ####
####OSTIA 1
sdev= sdev_fcst[0,:]
crmsd= rmsd_fcst[0,:]
ccoef= bias_fcst[0,:]
sm.taylor_diagram(sdev,
                  crmsd,
                  ccoef, markercolor ='r', alpha = 0.0)
###### ostia 2
sdev1= sdev_fcst[1,:]
crmsd1= rmsd_fcst[1,:]
ccoef1= bias_fcst[1,:]
sm.taylor_diagram(sdev1,
                  crmsd1,
                  ccoef1, markercolor ='g', alpha = 0.0,
                  overlay = 'on', markerLabel = label)
sdev2= sdev_fcst[2,:]
crmsd2= rmsd_fcst[2,:]
ccoef2= bias_fcst[2,:]
sm.taylor_diagram(sdev2,
                  crmsd2,
                  ccoef2, markercolor ='b', alpha = 0.0,
                  overlay = 'on', markerLabel = label)


    # Write plot to file
plt.savefig('taylor13.png',dpi=150,facecolor='w')

    # Show plot
plt.show()


#animate
#image_path = Path('/work/noaa/marine/Cameron.Book/ostia_check/sst_diff_figs')
#images = list(image_path.glob('*.png'))
#images.sort()
#gif_path = "ostia_movie.gif"

