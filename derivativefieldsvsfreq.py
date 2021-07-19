#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 11:13:49 2021

@author: mccallum
"""

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys
sys.path.append("../")
import IMUtils

#Load the maps
nside_crime = 512
start = 1
end = 150
N = end-start+1
gsyncI = np.zeros((N,hp.nside2npix(nside_crime)))
gsyncQ = np.zeros((N,hp.nside2npix(nside_crime)))
gsyncU = np.zeros((N,hp.nside2npix(nside_crime)))
psource = np.zeros((N,hp.nside2npix(nside_crime)))
egfreefree = np.zeros((N,hp.nside2npix(nside_crime)))
gfreefree = np.zeros((N,hp.nside2npix(nside_crime)))
h1cosm = np.zeros((N,hp.nside2npix(nside_crime)))

#Make sure to update paths to where data from CRIME is stored
for i in range(start,end+1):
    if i < 10:
        val = '00'+np.str(i)
    elif i<100:
        val = '0'+np.str(i)
    else:
        val = np.str(i)
    gsyncI[i-start],gsyncQ[i-start],gsyncU[i-start] = hp.read_map('./IM_CRIME_Data/gsync/gsync_'+val+'_b.fits',field=(0,1,2))
    psource[i-start] = hp.read_map('./IM_CRIME_Data/psources/psources_'+val+'.fits')
    egfreefree[i-start] = hp.read_map('./IM_CRIME_Data/egfree/egfree_'+val+'.fits')
    gfreefree[i-start] = hp.read_map('./IM_CRIME_Data/gfree/gfree_'+val+'.fits')
    h1cosm[i-start] = hp.read_map('./IM_CRIME_Data/cosmo/sim_2048_'+val+'.fits')


#Degrade to correct NSIDE
NSIDE = 128
gsyncI = hp.ud_grade(gsyncI, NSIDE)
gsyncQ = hp.ud_grade(gsyncQ, NSIDE)
gsyncU = hp.ud_grade(gsyncU, NSIDE)
psource = hp.ud_grade(psource, NSIDE)
egfreefree = hp.ud_grade(egfreefree, NSIDE)
gfreefree = hp.ud_grade(gfreefree, NSIDE)
h1cosm = hp.ud_grade(h1cosm, NSIDE)



#Beam smooth the map
beaminarcmin=60

for i in range(0,len(gsyncI)):
    #Smooth by beam
    gsyncI[i] = hp.smoothing(gsyncI[i],fwhm=np.deg2rad(beaminarcmin/60.))
    gsyncQ[i] = hp.smoothing(gsyncQ[i],fwhm=np.deg2rad(beaminarcmin/60.))
    gsyncU[i] = hp.smoothing(gsyncU[i],fwhm=np.deg2rad(beaminarcmin/60.))
    psource[i] = hp.smoothing(psource[i],fwhm=np.deg2rad(beaminarcmin/60.))
    egfreefree[i] = hp.smoothing(egfreefree[i],fwhm=np.deg2rad(beaminarcmin/60.))
    gfreefree[i] = hp.smoothing(gfreefree[i],fwhm=np.deg2rad(beaminarcmin/60.))
    h1cosm[i] = hp.smoothing(h1cosm[i],fwhm=np.deg2rad(beaminarcmin/60.))

    print(i)

I = h1cosm + gsyncI + psource + egfreefree + gfreefree
Q = gsyncQ
U = gsyncU



#Calculate the derivative of the intensity and polarization fields
ellmax=300

dI=np.zeros(shape=I.shape)
dU=np.zeros(shape=I.shape)
dQ=np.zeros(shape=I.shape)
for i in range(0,len(I)):
    #Spin raising/lowering operators
    dI[i] = IMUtils.get_derivativeoffield(I[i],ellmax)#(dIdth + 1j*dIdph)
    dQ[i] = IMUtils.get_derivativeoffield(Q[i],ellmax)#(dQdth + 1j*dQdph)
    dU[i] = IMUtils.get_derivativeoffield(U[i],ellmax)#(dUdth + 1j*dUdph)

dP = dQ + 1j *dU


#Load Crime frequency band info
n_freq, nu_0, nu_f, z_0, z_f = np.loadtxt('./IM_CRIME_Data/nuTable.txt').T



#Set the beam squint for the simulation
randomzeta=False
linearzeta=True
if randomzeta:
    seed=1
    np.random.seed(seed)
    zeta = np.random.random(N)*0.6 -0.3
    zeta = np.deg2rad(zeta/60)

elif linearzeta:
    zeta = np.linspace(-0.2,0.6,N)
    zeta = np.deg2rad(zeta/60)



#Plot the derivative fields as a function of frequency at RA=155degrees and DEC=3degrees
RA=155
DEC=3

plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 20})


fig,(ax1,ax2,ax3)=plt.subplots(nrows=3,sharex=True)
plt.subplot(311)
plt.plot(((nu_f+nu_0)/2)[start-1:end], zeta*np.real(dI.T[hp.ang2pix(NSIDE,np.deg2rad(DEC),np.deg2rad(RA),lonlat=True)]),label='$\Re \{\\zeta \, \eth I\}$')# (RA,DEC)=(155,3)')
plt.plot(((nu_f+nu_0)/2)[start-1:end], zeta*np.imag(dI.T[hp.ang2pix(NSIDE,np.deg2rad(DEC),np.deg2rad(RA),lonlat=True)]),label='$\Im \{\\zeta \, \eth I\}$')# (RA,DEC)=(155,3)')
plt.legend(loc='best',ncol=3)
plt.ylabel('mK')
plt.yticks(np.array([-200,0,200]))

plt.subplot(312)
plt.plot(((nu_f+nu_0)/2)[start-1:end], zeta*np.real(dP.T[hp.ang2pix(NSIDE,np.deg2rad(DEC),np.deg2rad(RA),lonlat=True)]),label='$\Re \{\\zeta \, \eth P\}$')# (RA,DEC)=(155,3)')
plt.plot(((nu_f+nu_0)/2)[start-1:end], zeta*np.imag(dP.T[hp.ang2pix(NSIDE,np.deg2rad(DEC),np.deg2rad(RA),lonlat=True)]),label='$\Im \{\\zeta \, \eth P\}$')# (RA,DEC)=(155,3)')
plt.legend(loc='best',ncol=3)
plt.ylabel('mK')
plt.subplots_adjust(hspace=.0)
plt.yticks(np.arange(-5,5,2.5))

plt.subplot(313)
plt.plot(((nu_f+nu_0)/2)[start-1:end],np.rad2deg(zeta)*60)
plt.xlabel('Frequency (MHz)')
plt.ylabel('$\\zeta$ (arcmin)')
if randomzeta:
    plt.yticks(np.arange(-0.2,0.3,0.2))
if linearzeta:
    plt.yticks(np.arange(-0.2,0.8,0.2))

