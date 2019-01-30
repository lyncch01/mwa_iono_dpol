import RMextract.PosTools as PosTools
import RMextract.getIONEX as ionex
import ephem
import os
import datetime
import pytz
import numpy as np
from math import *
import math

#====================USEFUL CONSTANTS==================================================
C   =2.99792458e8
r2h = 12.0/np.pi
r2d = 180.0/np.pi
h2d = 360.0/24.0
d2h = 24.0/360.0
ARCMIN2RAD = np.pi / (60.0 * 180.0)


# the following parameter is used to extend the observation range by 120 sec
# before and after the actual specified time. If we want to correct an
# actual data set, this is required for the scipy 1d interpolation routine to
# ensure that the calculated range of data exceeds the time range actually
# observed - I have used the same value in ALBUS - Tony
TIME_OFFSET = 120.0

# height of thin layer ionosphere
ION_HEIGHT=450.e3
sm_a = 6378137.0
invf = 298.257223563
ff = 1.0 / invf
HAS_PYRAP = False
HAS_EPHEM = True

#====================SOME CRAZINESS ABOUT LEAPSECONDS==================================


Leapseconds=[datetime.datetime(1972, 6, 30,0,0,0,0,pytz.utc),
             datetime.datetime(1972, 12, 31,0,0,0,0,pytz.utc),
             datetime.datetime(1973, 12, 31,0,0,0,0,pytz.utc),
             datetime.datetime(1974, 12, 31,0,0,0,0,pytz.utc),
             datetime.datetime(1975, 12, 31,0,0,0,0,pytz.utc),
             datetime.datetime(1976, 12, 31,0,0,0,0,pytz.utc),
             datetime.datetime(1977, 12, 31,0,0,0,0,pytz.utc),
             datetime.datetime(1978, 12, 31,0,0,0,0,pytz.utc),
             datetime.datetime(1979, 12, 31,0,0,0,0,pytz.utc),
             datetime.datetime(1981, 6, 30,0,0,0,0,pytz.utc),
             datetime.datetime(1982, 6, 30,0,0,0,0,pytz.utc),
             datetime.datetime(1983, 6, 30,0,0,0,0,pytz.utc),
             datetime.datetime(1985, 6, 30,0,0,0,0,pytz.utc),
             datetime.datetime(1987, 12, 31,0,0,0,0,pytz.utc),
             datetime.datetime(1989, 12, 31,0,0,0,0,pytz.utc),
             datetime.datetime(1990, 12, 31,0,0,0,0,pytz.utc),
             datetime.datetime(1992, 6, 30,0,0,0,0,pytz.utc),
             datetime.datetime(1993, 6, 30,0,0,0,0,pytz.utc),
             datetime.datetime(1994, 6, 30,0,0,0,0,pytz.utc),
             datetime.datetime(1995, 12, 31,0,0,0,0,pytz.utc),
             datetime.datetime(1997, 6, 30,0,0,0,0,pytz.utc),
             datetime.datetime(1998, 12, 31,0,0,0,0,pytz.utc),
             datetime.datetime(2005, 12, 31,0,0,0,0,pytz.utc),
             datetime.datetime(2008, 12, 31,0,0,0,0,pytz.utc),
             datetime.datetime(2012, 6, 30,0,0,0,0,pytz.utc),
	     datetime.datetime(2015, 1,  1,0,0,0,0,pytz.utc)]

# this is when GPStime=0 is defined
GPSzero=datetime.datetime(1980,1,6,0,0,0,0,pytz.utc)
Leapinterval_Start=[datetime.datetime(1960,1,1,0,0,0,0,pytz.utc)]
Leapinterval_End=[]
SumLeapseconds=[0]
for ls in Leapseconds:
	Leapinterval_End.append(ls)
	SumLeapseconds.append(SumLeapseconds[-1]+1)
	Leapinterval_Start.append(ls)

Leapinterval_End.append(datetime.datetime(2016,1,1,0,0,0,0,pytz.utc)) ##### DO I NEED TO CHANGE?

# number of leap seconds that had already elapsed at GPStime=0
GPSzero_leapseconds=SumLeapseconds[np.where(np.array(Leapinterval_Start)<=GPSzero)[0][-1]]
Leapinterval_Start_GPS=[]
Leapinterval_End_GPS=[]
ls=0
for i in range(len(Leapinterval_Start)):
	dt=Leapinterval_Start[i]-GPSzero+datetime.timedelta(seconds=ls-GPSzero_leapseconds)
	Leapinterval_Start_GPS.append(dt.days*86400 + dt.seconds)
	dt=Leapinterval_End[i]-GPSzero+datetime.timedelta(seconds=ls-GPSzero_leapseconds)
	Leapinterval_End_GPS.append(dt.days*86400 + dt.seconds)
	ls=SumLeapseconds[i]    


#====================FUNCTIONS==========================================================


################################################################################
## datetime=gps_datetime(gps); converts from gps seconds to datetime.datetime

def gps_datetime(gps):
	ls=SumLeapseconds[np.where(np.array(Leapinterval_Start_GPS)<=gps)[0][-1]]
	dt=datetime.timedelta(seconds=gps-(ls-GPSzero_leapseconds))
	return GPSzero + dt

################################################################################
### Convert latitude (radians), longitude (radians) and elevation (metres) to ITRF XYZ and vice versa

def WGS84ToITRF(lat, lon, h): # WGS-84 to ITRF (input in radians)
	SINK = math.sin(lat)
	COSK = math.cos(lat)
	e2 = 2.0 * ff - ff * ff
	v = sm_a / math.sqrt(1.0 - e2 * SINK * SINK)
	x = (v + h) * COSK * math.cos(lon)
	y = (v + h) * COSK * math.sin(lon)
	z = ((1 - e2) * v + h) * SINK
	return x, y, z

def ITRFToWGS84(x, y, z): #ITRF to WGS-84
	e2 = 2.0 * ff - ff * ff
	E = e2 / (1.0 - e2)
	b = sm_a * (1.0 - ff)
	p = math.sqrt(x * x + y * y)
	q = math.atan2(z * sm_a, (p * b))
	lat = math.atan2((z + E * b * math.sin(q) * math.sin(q) * math.sin(q)), (p - e2 * sm_a * math.cos(q) * math.cos(q) * math.cos(q)))
	v = sm_a / math.sqrt(1.0 - e2 * math.sin(lat) * math.sin(lat))
	lon = math.atan2(y, x)
	h = (p / math.cos(lat)) - v
	lat = math.degrees(lat)
	lon = math.degrees(lon)
	return lat, lon, h       # output in degrees here


################################################################################
### Convert geodetic latitude to geocentric latitude input and output in radians

def GeodeticToGeocentricLat(geodetic_lat, height):
	l_sin = math.sin(geodetic_lat)
	e2 = 2.0 * ff - ff * ff
	div = math.sqrt(1 - e2 * l_sin**2)
	rn =  sm_a / div 
	rn_div = rn + height
	ratio = 1 - e2 * rn / rn_div
	tan_geocentric_lat = ratio * math.tan(geodetic_lat) 
	geocentric_lat = math.atan(tan_geocentric_lat)
	return geocentric_lat


################################################################################
### See http://stackoverflow.com/questions/15954978/ecef-from-azimuth-elevation-range-and-observer-lat-lon-alt

def aer2ecef(azimuthDeg, elevationDeg, slantRange, obs_lat, obs_lon, obs_alt): # input expected in degrees here
	obs_lat_r = math.radians(obs_lat)
	obs_lon_r = math.radians(obs_lon)
	sitex, sitey, sitez = WGS84ToITRF(obs_lat_r,obs_lon_r,obs_alt)

	#some needed calculations
	slat = math.sin(math.radians(obs_lat))
	slon = math.sin(math.radians(obs_lon))
	clat = math.cos(math.radians(obs_lat))
	clon = math.cos(math.radians(obs_lon))

	azRad = math.radians(azimuthDeg)
	elRad = math.radians(elevationDeg)

	# az,el,range to sez convertion
	south  = -slantRange * math.cos(elRad) * math.cos(azRad)
	east   =  slantRange * math.cos(elRad) * math.sin(azRad)
	zenith =  slantRange * math.sin(elRad)


	x = ( slat * clon * south) + (-slon * east) + (clat * clon * zenith) + sitex
	y = ( slat * slon * south) + ( clon * east) + (clat * slon * zenith) + sitey
	z = (-clat *        south) + ( slat * zenith) + sitez
	lat,lon,h = ITRFToWGS84(x, y, z)
    	return lat, lon, h  # data are in units of degrees


################################################################################
### Get the year, month, day, day_fraction from an MJD; Taken from _Astronomical Algorithms_, Meeus, 1991
def get_ymdf_from_JD(JD):   
	JD2 = JD + 0.5
    	Z = int(JD2)
    	F = JD2 - Z
    	A = Z
    	if(Z>= 2299161):
    	    alpha = int((Z-1867216.25)/36524.25)
    	    A = Z + 1 + alpha - int(alpha/4)
    	B = A + 1524
    	C = int((B-122.1)/365.25)
    	D = int(365.25*C)
    	E = int((B-D)/30.6001)
    	day_total = B - D - int(30.6001*E) + F
    	month = E - 1
    	if(E >= 14): month = E - 13
	year = C - 4716
    	if(month <= 2): year += 1
   	day = int(day_total)
   	day_fraction = day_total - day
    	return year, month, day, day_fraction

################################################################################
### Get hours, minues, seconds from a fractional day; does not worry about leap seconds.
def get_hms_from_frac(day_fraction):
    	h = day_fraction * 24.0
    	hour = int(h+2E-13)
    	m = (h - hour) * 60.0
    	minute = int(m+1E-11)
    	second = (m - minute) * 60.0
    	return hour, minute, second

################################################################################
def get_ymdh_from_JD(JD):
    """get hours, minues, seconds from a fractional day.
"""
    year, month, day, day_fraction = get_ymdf_from_JD(JD)
    hour, minute, second = get_hms_from_frac(day_fraction)
    return year, month, day, hour, minute, second

################################################################################
def obtain_observation_year_month_day_fraction(start_time):
    	julian_day = (start_time / 86400.0) + 2400000.5
    	result = get_ymdf_from_JD(julian_day)
    	return result

################################################################################
# Get the day of year from the Year, month, day for the start of observations
def obtain_observation_year_month_day_hms(start_time):
	julian_day = (start_time / 86400.0) + 2400000.5
	year, month, day, hour, minute, second  = get_ymdh_from_JD(julian_day)
	return (year, month, day,hour, minute, second)

