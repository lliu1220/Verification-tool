import mpu
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

from netCDF4 import Dataset

defaultcm = cmo.thermal
import sys, os, glob
sys.path.insert(0,'/home/PyUils/')
import matplotlib
matplotlib.use("TkAgg")
import numpy as np
import matplotlib.pyplot as plt
from utils import coast180
import matplotlib.dates as mdates
import datetime
from utils4HWRF import readTrack6hrly, plot_Saffir_Simpson_WPscale

from pathlib import Path

start_date = datetime.date(2020, 7, 28)
end_date = datetime.date(2020, 7, 30) #datetime.date(2003,  3, 31)

wdir='/work/noaa/hwrf/lliu/scrub/'
btk='bal092020.dat'
bdtb,blonb,blatb,bpmnb,bvmxb=readTrack6hrly(btk)
blonb=blonb[19:40] ###to match best track with storm track
blatb=blatb[19:40]
bdtb=bdtb[19:40]
bpmnb=bpmnb[19:40]
blonb=[-1*a for a in blonb]
#print('check',blonb)
std=np.zeros([5,21])
bias=np.zeros([5,21])
j=0
for expname in ['hafs_couplehycom_cdeps-202007-mom6_NODA_spinup','hafs_couplehycom_cdeps-202007-mom6_3dvar_spinup','hafs_couplehycom_cdeps-202007-mom6_GLD_1012','hafs_couplehycom_cdeps-202007-mom6_ALL_update_1004','hafs_develop_202102_rt_regional_C192_1n4_static_new']:
 print('expname',expname)
 expdir=wdir+expname
 ini_date=start_date
 tmp_date=ini_date
# print('ini_date',ini_date)

# bdt=np.zeros([4,21])
 blon=np.zeros([4,21])
 blat=np.zeros([4,21])
 bpmn=np.zeros([4,21])
 bvmx=np.zeros([4,21])
 dist=np.zeros([4,21])
 delta=datetime.timedelta(days=1)
 i=0
 while tmp_date<=datetime.date(2020,7,28):
  #tmp_date=ini_date
  print('date check',tmp_date,ini_date+5*delta)

 # while tmp_date <= ini_date+5*delta: #end_date:
  print('i',i)
  thisdate=datetime.datetime.strptime(str(tmp_date),"%Y-%m-%d").strftime("%Y%m%d")

  btkfile=expdir+'/com/'+str(thisdate)+'12/00L/'+'natl00l.'+str(thisdate)+'12.trak.hafs.atcfunix.all'
  print('btkfile',btkfile)
  bdt1,blon1,blat1,bpmn1,bvmx1=readTrack6hrly(btkfile)
 # print('check',np.shape(blon1))

  blon1=blon1[0:21]
  blat1=blat1[0:21]
  bpmn1=bpmn1[0:21]
  bvmx1=bvmx1[0:21]

  blon[i,:]=blon1
  blat[i,:]=blat1
  bpmn[i,:]=bpmn1
  bvmx[i,:]=bvmx1
  #print('shapes',np.shape(blon))
   #bdt[i,:],blon[i,:],blat[i,:],bpmn[i,:],bvmx[i,:]=readTrack6hrly(btkfile)

  bpmn[i,bpmn1<100]=np.nan

  blon[i,:]=[-1*a for a in blon[i,:]]
  #print('blon',blon)
  for x in range(0,21):
    dist[i,x] = mpu.haversine_distance((blatb[x],blonb[x]),(blat[i,x],blon[i,x]))
  dist=dist/1.852
#  print('shapes',np.shape(dist))
  i=i+1
  tmp_date=tmp_date+delta
 std[j,:]=np.std(dist,axis=0)
 bias[j,:]=np.mean(dist,axis=0)
 j=j+1
print('checks',(std),bias)

fig=plt.figure(figsize=(10,6.))


fig, ax=plt.subplots(figsize=(8,5.5))
plt.title('Track error Isaias 5-day',fontsize=13)
plt.title('Isaias Track 5-day forecast, 3 day cycles')
#llabs=['BT',EXPT]

plt.plot(bdtb,std[0,:],'-o',color='green',markerfacecolor='w',markersize=2,label='NODA',linewidth=3,alpha=0.7)
plt.plot(bdtb,std[1,:],'-o',color='blue',markerfacecolor='w',markersize=2,label='CTRL',linewidth=3,alpha=0.7)
plt.plot(bdtb,std[2,:],'-o',color='magenta',markerfacecolor='w',markersize=2,label='GLD',linewidth=3,alpha=0.7)
plt.plot(bdtb,std[3,:],'-o',color='red',markerfacecolor='w',markersize=2,label='ALL',linewidth=3,alpha=0.7)
plt.plot(bdtb,std[4,:],'-o',color='grey',markerfacecolor='w',markersize=2,label='HYCOM',linewidth=3,alpha=0.7)
ax.xaxis.set_major_locator(mdates.DayLocator())
ax.legend(["NODA","CTRL","GLD","ALL","HYCOM"],loc="lower center")
ax.grid()
#plt.ylim([960,1020])
#plt.ylim([1000,1020])
plt.ylabel('Absolute Track error',fontsize=12)
plt.xlabel('Year 2020',fontsize=12)
plt.savefig('error.png',bbox_inches=0)


