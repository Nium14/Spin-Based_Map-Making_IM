#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 14:29:35 2021

@author: mccallum
"""
import healpy as hp
import numpy as np
import sys
import matplotlib.pyplot as plt


sys.path.append("../")
import IMUtils

#Load in the relevant standard MeerKAT TOD information
ra,dec,az,psi = np.loadtxt('./TOD.txt')
az[az==360] = 0.

NSIDE = 128
N = 1
I=np.ones((N,hp.nside2npix(NSIDE)))
pixel_index = hp.ang2pix(NSIDE,dec,ra)



#Calculate Hitmap and hmaps
hitmap = np.zeros((N,len(I[0])))
for i in range(0,len(I)):
    np.add.at(hitmap[i],pixel_index,1)

hmaps = IMUtils.get_hmaps(psi, pixel_index, np.array([0,1,2,3,4]), NSIDE)
hmaps/=hmaps[0]


#Instantiate arrays for the summed TOD to be put in
tod_A = np.ones((N,len(az)))
tod_B = np.ones((N,len(az)))
#Sum the detector pairs and add them to the appropriate TOD array
sum_AB = 0.5*(tod_A + tod_B)




#Perform spin-0 map-making
T_map = np.zeros((N,1,len(I[0])))
pixcond_T = np.zeros((N,len(I[0])))
for i in range(0,len(I)):
    T_map[i],pixcond_T[i] = IMUtils.mapmake_simplebin_spink(T_map[i], pixcond_T[i], sum_AB[i], pixel_index, psi, NSIDE, spins=np.array([0]),mask=None)


#Perform spin-0 and 2 map-making
IQU_Out = np.zeros((N,3,len(I[0])))
pixcond_IQU = np.zeros((N,len(I[0])))
for i in range(0,len(I)):
    IQU_Out[i],pixcond_IQU[i] = IMUtils.mapmake_simplebin_spink(IQU_Out[i],pixcond_IQU[i], sum_AB[i], pixel_index, psi, NSIDE, spins=np.array([0,2]),mask=None)



#Perform spin-0 and 1 and 3 map-making
IQU_Out013 = np.zeros((N,5,len(I[0])))
pixcond_IQU013 = np.zeros((N,len(I[0])))
for i in range(0,len(I)):
    IQU_Out013[i],pixcond_IQU013[i] = IMUtils.mapmake_simplebin_spink(IQU_Out013[i],pixcond_IQU013[i], sum_AB[i], pixel_index, psi, NSIDE, spins=np.array([0,1,3]),mask=None)


#Find masks based on how well conditioned the map-making is
connum = 1e11
mask0 = 1.*pixcond_T[0]/pixcond_T[0]
mask0[pixcond_T[0]>connum]=np.nan

mask02 = 1.*pixcond_IQU[0]/pixcond_IQU[0]
mask02[pixcond_IQU[0]>connum]=np.nan

mask013 = 1.*pixcond_IQU013[0]/pixcond_IQU013[0]
mask013[pixcond_IQU013[0]>connum]=np.nan



#####As above but for MeerKAT plus low elevations#####

#Load in the relevant standard MeerKAT TOD information
ra,dec,az,psi = np.loadtxt('./TOD_extralowelevs.txt')
az[az==360] = 0.

I=np.ones((N,hp.nside2npix(NSIDE)))
pixel_index = hp.ang2pix(NSIDE,dec,ra)



#Calculate Hitmap and hmaps
hitmap_Low = np.zeros((N,len(I[0])))
for i in range(0,len(I)):
    np.add.at(hitmap_Low[i],pixel_index,1)

hmaps_Low = IMUtils.get_hmaps(psi, pixel_index, np.array([0,1,2,3,4]), NSIDE)
hmaps_Low/=hmaps_Low[0]

#Instantiate arrays for the summed TOD to be put in
tod_A = np.ones((N,len(az)))
tod_B = np.ones((N,len(az)))
#Sum the detector pairs and add them to the appropriate TOD array
sum_AB = 0.5*(tod_A + tod_B)


#Instantiate arrays for the summed TOD to be put in
tod_A = np.ones((N,len(az)))
tod_B = np.ones((N,len(az)))
#Sum the detector pairs and add them to the appropriate TOD array
sum_AB = 0.5*(tod_A + tod_B)




#Perform spin-0 map-making
T_map = np.zeros((N,1,len(I[0])))
pixcond_T = np.zeros((N,len(I[0])))
for i in range(0,len(I)):
    T_map[i],pixcond_T[i] = IMUtils.mapmake_simplebin_spink(T_map[i], pixcond_T[i], sum_AB[i], pixel_index, psi, NSIDE, spins=np.array([0]),mask=None)


#Perform spin-0 and 2 map-making
IQU_Out = np.zeros((N,3,len(I[0])))
pixcond_IQU = np.zeros((N,len(I[0])))
for i in range(0,len(I)):
    IQU_Out[i],pixcond_IQU[i] = IMUtils.mapmake_simplebin_spink(IQU_Out[i],pixcond_IQU[i], sum_AB[i], pixel_index, psi, NSIDE, spins=np.array([0,2]),mask=None)



#Perform spin-0 and 1 and 3 map-making
IQU_Out013 = np.zeros((N,5,len(I[0])))
pixcond_IQU013 = np.zeros((N,len(I[0])))
for i in range(0,len(I)):
    IQU_Out013[i],pixcond_IQU013[i] = IMUtils.mapmake_simplebin_spink(IQU_Out013[i],pixcond_IQU013[i], sum_AB[i], pixel_index, psi, NSIDE, spins=np.array([0,1,3]),mask=None)


#Find masks based on how well conditioned the map-making is
connum = 1e11
mask0_Low = 1.*pixcond_T[0]/pixcond_T[0]
mask0_Low[pixcond_T[0]>connum]=np.nan

mask02_Low = 1.*pixcond_IQU[0]/pixcond_IQU[0]
mask02_Low[pixcond_IQU[0]>connum]=np.nan

mask013_Low = 1.*pixcond_IQU013[0]/pixcond_IQU013[0]
mask013_Low[pixcond_IQU013[0]>connum]=np.nan






#####As above but for MeerKAT plus high elevations#####

#Load in the relevant standard MeerKAT TOD information
ra,dec,az,psi = np.loadtxt('./TOD_extrahighelevs.txt')
az[az==360] = 0.
I=np.ones((N,hp.nside2npix(NSIDE)))
pixel_index = hp.ang2pix(NSIDE,dec,ra)


#Calculate Hitmap and hmaps
hitmap_High = np.zeros((N,len(I[0])))
for i in range(0,len(I)):
    np.add.at(hitmap_High[i],pixel_index,1)

hmaps_High = IMUtils.get_hmaps(psi, pixel_index, np.array([0,1,2,3,4]), NSIDE)
hmaps_High/=hmaps_High[0]

#Instantiate arrays for the summed TOD to be put in
tod_A = np.ones((N,len(az)))
tod_B = np.ones((N,len(az)))
#Sum the detector pairs and add them to the appropriate TOD array
sum_AB = 0.5*(tod_A + tod_B)




#Perform spin-0 map-making
T_map = np.zeros((N,1,len(I[0])))
pixcond_T = np.zeros((N,len(I[0])))
for i in range(0,len(I)):
    T_map[i],pixcond_T[i] = IMUtils.mapmake_simplebin_spink(T_map[i], pixcond_T[i], sum_AB[i], pixel_index, psi, NSIDE, spins=np.array([0]),mask=None)


#Perform spin-0 and 2 map-making
IQU_Out = np.zeros((N,3,len(I[0])))
pixcond_IQU = np.zeros((N,len(I[0])))
for i in range(0,len(I)):
    IQU_Out[i],pixcond_IQU[i] = IMUtils.mapmake_simplebin_spink(IQU_Out[i],pixcond_IQU[i], sum_AB[i], pixel_index, psi, NSIDE, spins=np.array([0,2]),mask=None)



#Perform spin-0 and 1 and 3 map-making
IQU_Out013 = np.zeros((N,5,len(I[0])))
pixcond_IQU013 = np.zeros((N,len(I[0])))
for i in range(0,len(I)):
    IQU_Out013[i],pixcond_IQU013[i] = IMUtils.mapmake_simplebin_spink(IQU_Out013[i],pixcond_IQU013[i], sum_AB[i], pixel_index, psi, NSIDE, spins=np.array([0,1,3]),mask=None)


#Find masks based on how well conditioned the map-making is
connum = 1e11
mask0_High = 1.*pixcond_T[0]/pixcond_T[0]
mask0_High[pixcond_T[0]>connum]=np.nan

mask02_High = 1.*pixcond_IQU[0]/pixcond_IQU[0]
mask02_High[pixcond_IQU[0]>connum]=np.nan

mask013_High = 1.*pixcond_IQU013[0]/pixcond_IQU013[0]
mask013_High[pixcond_IQU013[0]>connum]=np.nan



#Set extent of plot
lonra=[-22.5,17.5]
latra=[-6,6]

#Generate hitmaps*mask plots
fig = plt.figure(figsize=(20, 6))
ax1 = fig.add_subplot(3,2,1)
plt.axes(ax1)
hit_02=hp.cartview(hitmap[0]*mask02,rot=(162.5,3,0),lonra=lonra,latra=latra,hold=True,title='',cbar=None,return_projected_map=True)
ax2 = fig.add_subplot(3,2,2)
plt.axes(ax2)
hit_02_Low=hp.cartview(hitmap_Low[0]*mask02_Low,rot=(162.5,3,0),lonra=lonra,latra=latra,hold=True,title='',cbar=None,return_projected_map=True)
ax3 = fig.add_subplot(3,2,3)
plt.axes(ax3)
hit_013=hp.cartview(hitmap[0]*mask013,rot=(162.5,3,0),lonra=lonra,latra=latra,hold=True,title='',cbar=None,return_projected_map=True)
ax4 = fig.add_subplot(3,2,4)
plt.axes(ax4)
hit_013_High=hp.cartview(hitmap_High[0]*mask013_High,rot=(162.5,3,0),lonra=lonra,latra=latra,hold=True,title='',cbar=None,return_projected_map=True)
ax5 = fig.add_subplot(3,2,5)
plt.axes(ax5)
hit_013=hp.cartview(hitmap[0]*mask013,rot=(162.5,3,0),lonra=lonra,latra=latra,hold=True,title='',cbar=None,return_projected_map=True)
ax6 = fig.add_subplot(3,2,6)
plt.axes(ax6)
hit_013_Low=hp.cartview(hitmap_Low[0]*mask013_Low,rot=(162.5,3,0),lonra=lonra,latra=latra,hold=True,title='',cbar=None,return_projected_map=True)
plt.close()






plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 24})


#Spin 0, 2 MeerKAT comparison to MeerKAT plus Low elevation
fig = plt.figure(figsize=(8, 3))
ax1 = fig.add_subplot(1,2,1)
plt.axes(ax1)
im1 = ax1.imshow(hit_02, origin='lower',extent=(160+lonra[1],160+lonra[0],3+latra[0],3+latra[1]))
plt.xlabel('R.A.')
plt.ylabel('Dec.')
plt.yticks(np.arange(-2,10,2))
clb=fig.colorbar(im1, orientation="horizontal")
clb.set_label('Hits')
ax2 = fig.add_subplot(1,2,2)
plt.axes(ax2)
im2 = ax2.imshow(hit_02_Low, origin='lower',extent=(160+lonra[1],160+lonra[0],3+latra[0],3+latra[1]))
plt.xlabel('R.A.')
plt.ylabel('Dec.')
plt.yticks(np.arange(-2,10,2))
clb=fig.colorbar(im2, orientation="horizontal")
clb.set_label('Hits')


#Spin 0, 1, 3 MeerKAT comparison to MeerKAT plus High elevation
fig = plt.figure(figsize=(8, 3))
ax3 = fig.add_subplot(1,2,1)
plt.axes(ax3)
im3 = ax3.imshow(hit_013, origin='lower',extent=(160+lonra[1],160+lonra[0],3+latra[0],3+latra[1]))
plt.xlabel('R.A.')
plt.ylabel('Dec.')
plt.yticks(np.arange(-2,10,2))
clb=fig.colorbar(im3, orientation="horizontal")
clb.set_label('Hits')
ax4 = fig.add_subplot(1,2,2)
plt.axes(ax4)
im4 = ax4.imshow(hit_013_High, origin='lower',extent=(160+lonra[1],160+lonra[0],3+latra[0],3+latra[1]))
plt.xlabel('R.A.')
plt.ylabel('Dec.')
plt.yticks(np.arange(-2,10,2))
clb=fig.colorbar(im4, orientation="horizontal")
clb.set_label('Hits')


#Spin 0, 1, 3 MeerKAT comparison to MeerKAT plus Low elevation
fig = plt.figure(figsize=(8, 3))
ax3 = fig.add_subplot(1,2,1)
plt.axes(ax3)
im3 = ax3.imshow(hit_013, origin='lower',extent=(160+lonra[1],160+lonra[0],3+latra[0],3+latra[1]))
plt.xlabel('R.A.')
plt.ylabel('Dec.')
plt.yticks(np.arange(-2,10,2))
clb=fig.colorbar(im3, orientation="horizontal")
clb.set_label('Hits')
ax4 = fig.add_subplot(1,2,2)
plt.axes(ax4)
im4 = ax4.imshow(hit_013_Low, origin='lower',extent=(160+lonra[1],160+lonra[0],3+latra[0],3+latra[1]))
plt.xlabel('R.A.')
plt.ylabel('Dec.')
plt.yticks(np.arange(-2,10,2))
clb=fig.colorbar(im4, orientation="horizontal")
clb.set_label('Hits')






#Print |h_k|^2
print('|h_1|^2 = '+np.str(np.nanmean(np.abs(hmaps[1])**2)))
print('|h_2|^2 = '+np.str(np.nanmean(np.abs(hmaps[2])**2)))
print('|h_3|^2 = '+np.str(np.nanmean(np.abs(hmaps[3])**2)))
print('|h_4|^2 = '+np.str(np.nanmean(np.abs(hmaps[4])**2)))


print('Low |h_1|^2 = '+np.str(np.nanmean(np.abs(hmaps_Low[1])**2)))
print('Low |h_2|^2 = '+np.str(np.nanmean(np.abs(hmaps_Low[2])**2)))
print('Low |h_3|^2 = '+np.str(np.nanmean(np.abs(hmaps_Low[3])**2)))
print('Low |h_4|^2 = '+np.str(np.nanmean(np.abs(hmaps_Low[4])**2)))
