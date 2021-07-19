#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 16:00:02 2021

@author: mccallum
"""

import healpy as hp
import ephem
import numpy as np
import concurrent.futures
from functools import partial
import os
import sys
from ephem import newton

sys.path.append("../")
import IMUtils

#Location of Observatory
long_in = '21:24.40'
lat_in = '-30:43.16'
elevation_in = 1062.5


#Scan Inputs  Based on MeerKAT paper
#Slew in azimuth at rate degrees/second
az_freq = 0.083333333

delta_az= 0.013333333

#Sampling frequency in Hz
#Make this at least enough to sample within beam FWHM size ~10 times
FWHM=60./60.
fsamp=0.5

#Survey time
surveytime=60.*60.*1.5#23.9345
time_array = np.arange(0,surveytime,1./fsamp)

#Setup the pyephem observer
pb = ephem.Observer()
pb.elevation = elevation_in
pb.long = long_in
pb.lat = lat_in




RA=[]
DEC=[]
AZ=[]
PSI=[]


NSIDE=128         



#Run the MeerKAT Paper Scan Strategy (2011.13789)
MEERKATelevations=True
if MEERKATelevations==True:
    #Scan 1
    elev_array=[42.0]
    pb.date=ephem.Date('2019/02/24 19:48:28.00')
    azstart=41.6
    azend=59.6
    ra,dec,az,psi=IMUtils.generate_scan(elev_array[0],azend,azstart,pb,surveytime,fsamp,az_freq,delta_az)
    RA.append(ra)
    DEC.append(dec)
    AZ.append(az)
    PSI.append(psi)
    
    #Scan 2
    elev_array=[41.0]
    pb.date=ephem.Date('2019/02/25 00:40:11.00')
    azstart=-61
    azend=-43
    ra,dec,az,psi=IMUtils.generate_scan(elev_array[0],azend,azstart,pb,surveytime,fsamp,az_freq,delta_az)
    RA.append(ra)
    DEC.append(dec)
    AZ.append(az)
    PSI.append(psi)
    
    #Scan 3
    elev_array=[40.5]
    pb.date=ephem.Date('2019/03/30 17:19:02.00')
    azstart=43.7
    azend=61.7
    ra,dec,az,psi=IMUtils.generate_scan(elev_array[0],azend,azstart,pb,surveytime,fsamp,az_freq,delta_az)
    RA.append(ra)
    DEC.append(dec)
    AZ.append(az)
    PSI.append(psi)
    
    #Scan 4
    elev_array=[41.0]
    pb.date=ephem.Date('2019/04/01 22:06:17.00')
    azstart=-61
    azend=-43
    ra,dec,az,psi=IMUtils.generate_scan(elev_array[0],azend,azstart,pb,surveytime,fsamp,az_freq,delta_az)
    RA.append(ra)
    DEC.append(dec)
    AZ.append(az)
    PSI.append(psi)
    
    #Scan 5
    elev_array=[41.5]
    pb.date=ephem.Date('2019/04/23 20:41:56.00')
    azstart=-60.3
    azend=-42.3
    ra,dec,az,psi=IMUtils.generate_scan(elev_array[0],azend,azstart,pb,surveytime,fsamp,az_freq,delta_az)
    RA.append(ra)
    DEC.append(dec)
    AZ.append(az)
    PSI.append(psi)
    
    #Scan 6
    elev_array=[41.5]
    pb.date=ephem.Date('2019/04/24 20:39:57.00')
    azstart=-60.3
    azend=-42.3
    ra,dec,az,psi=IMUtils.generate_scan(elev_array[0],azend,azstart,pb,surveytime,fsamp,az_freq,delta_az)
    RA.append(ra)
    DEC.append(dec)
    AZ.append(az)
    PSI.append(psi)
    
    #Scan 7
    elev_array=[43.4]
    pb.date=ephem.Date('2019/07/11 15:09:53.00')
    azstart=-55.3
    azend=-37.3
    ra,dec,az,psi=IMUtils.generate_scan(elev_array[0],azend,azstart,pb,surveytime,fsamp,az_freq,delta_az)
    RA.append(ra)
    DEC.append(dec)
    AZ.append(az)
    PSI.append(psi)









extrahighelevation=False
if extrahighelevation==True:
    '''Target Field Limits'''
    LowerRA = 140+11
    UpperRA = 175+11
    LowerDec = -1
    UpperDec= 8
    
    #Rising Scan
    elev_array=[49.0]
    #Get time visible and az limits
    pb,azstart,azend = IMUtils.find_field_timeandaz(LowerRA,UpperRA,LowerDec,UpperDec,pb,elev_array[0],True,lat_in,long_in,elevation_in)
    while azstart>180:
        azstart-=360
    while azend>180:
        azend-=360
    ra,dec,az,psi=IMUtils.generate_scan(elev_array[0],azend,azstart,pb,surveytime,fsamp,az_freq,delta_az)
    RA.append(ra)
    DEC.append(dec)
    AZ.append(az)
    PSI.append(psi)
        

    #Setting Scan
    elev_array=[49.0]
    #Get time visible and az limits
    pb,azstart,azend = IMUtils.find_field_timeandaz(LowerRA,UpperRA,LowerDec,UpperDec,pb,elev_array[0],False,lat_in,long_in,elevation_in)
    azend = -np.abs(azend)#quick fix for now
    while azstart>180:
        azstart-=360
    while azend>180:
        azend-=360
    ra,dec,az,psi=IMUtils.generate_scan(elev_array[0],azend,azstart,pb,surveytime,fsamp,az_freq,delta_az)
    RA.append(ra)
    DEC.append(dec)
    AZ.append(az)
    PSI.append(psi)



extralowelevation=False
if extralowelevation==True:
    '''Target Field Limits'''
    LowerRA = 140+11
    UpperRA = 175+11
    LowerDec = -1
    UpperDec= 8

    #Rising Scan
    elev_array=[39.0]
    #Get time visible and az limits
    pb,azstart,azend = IMUtils.find_field_timeandaz(LowerRA,UpperRA,LowerDec,UpperDec,pb,elev_array[0],True,lat_in,long_in,elevation_in)
    while azstart>180:
        azstart-=360
    while azend>180:
        azend-=360
    
    ra,dec,az,psi=IMUtils.generate_scan(elev_array[0],azend,azstart,pb,surveytime,fsamp,az_freq,delta_az)
    RA.append(ra)
    DEC.append(dec)
    AZ.append(az)
    PSI.append(psi)


    #Setting Scan
    elev_array=[39.0]
    #Get time visible and az limits
    pb,azstart,azend = IMUtils.find_field_timeandaz(LowerRA,UpperRA,LowerDec,UpperDec,pb,elev_array[0],False,lat_in,long_in,elevation_in)
    azend = -np.abs(azend)#quick fix for now
    while azstart>180:
        azstart-=360
    while azend>180:
        azend-=360
    
    ra,dec,az,psi=IMUtils.generate_scan(elev_array[0],azend,azstart,pb,surveytime,fsamp,az_freq,delta_az)
    RA.append(ra)
    DEC.append(dec)
    AZ.append(az)
    PSI.append(psi)    




#Save the data
RA = np.ravel(RA)
DEC = np.ravel(DEC)
AZ = np.ravel(AZ)
PSI=np.ravel(PSI)

TOD = np.array([RA,DEC,AZ,PSI])

if extralowelevation == True and MEERKATelevations==True and extrahighelevation==False:
    np.savetxt('TOD_extralowelevs.txt',TOD)
elif MEERKATelevations==True and extralowelevation == False and extrahighelevation==False:
    np.savetxt('TOD.txt',TOD)
elif MEERKATelevations==False and extralowelevation == True and extrahighelevation==False:
    np.savetxt('TOD_justlowelevs.txt',TOD)
elif extrahighelevation==True and MEERKATelevations==True and extralowelevation==False:
    np.savetxt('TOD_extrahighelevs.txt',TOD)
elif extrahighelevation==True and MEERKATelevations==False and extralowelevation==False:
    np.savetxt('TOD_justextrahighelevs.txt',TOD)