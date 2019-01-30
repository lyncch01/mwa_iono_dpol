from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
import RMextract.getRM as RM
import datetime
import numpy as np
from math import *
import rm_library as rml
import matplotlib.pyplot as plt



x, y, z = rml.WGS84ToITRF(-26.703319 * np.pi / 180.0, 116.670815 * np.pi / 180.0, 0.0)
MWA_antennas = [[x, y, z]]


OBJECT="EoR0"
date = datetime.datetime(2015,12,8)
sc = SkyCoord("00h00m0.0s", "-27d00m00.0s", frame=FK5)
stime = [12., 0., 0.]
etime = [14., 45., 0.]

radec = (sc.ra.radian, sc.dec.radian)
START_TIME = "%4d-%02d-%02d %02d:%02d:%02d" %(date.year, date.month, date.day, stime[0], stime[1], stime[2])
END_TIME =  "%4d-%02d-%02d %02d:%02d:%02d" %(date.year, date.month, date.day, etime[0], etime[1], etime[2])
TIME_STEP = 120.0
print("Processing for %s; %s - %s" %(OBJECT, START_TIME, END_TIME))


result = RM.getRM(use_azel=False, start_time=START_TIME, end_time=END_TIME, object=OBJECT,radec=radec, timestep=TIME_STEP, stat_positions=MWA_antennas, useEMM=True, EMMpath='/home/corvus/Research/RMextract/EMM/')

# set time relative to zero
timegrid = result['times']
timegrid = (timegrid - timegrid[0]) / 3600.0
RMs = result['RM'] 
TECs = result['STEC']
times = result['times']
selected_key = False
for key in TECs.keys():
	if not selected_key:
		use_key = key
		selected_key = True
rmvals = RMs[use_key]
times = []
for t in timegrid:
	aa0 = stime[0] + (stime[1]/60.) + (stime[2]/3600.)
	aa1 = aa0 + t
	times.append(aa1)

plt.figure()
plt.scatter(times, rmvals)
#plt.ylim(-3.5, -0.5)
plt.xlabel('UT (hr)')
plt.ylabel(r'Ionospheric RM (rad m$^{-2}$)')
plt.show()


fout = open("%s_%s_.csv" %(OBJECT, START_TIME.split()[0]), "wt")
fout.write("#UTC,RM,STEC\n")
for i in range(len(timegrid)):
	if np.abs(RMs[use_key][i]) > 50.0:
		continue
	fout.write("%f,%f,%f\n" %(timegrid[i], RMs[use_key][i], TECs[use_key][i]))
fout.close()
