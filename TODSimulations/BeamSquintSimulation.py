#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 11:13:11 2021

@author: mccallum
"""

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys
sys.path.append("../")
import IMUtils


#Set ScanType to 'MeerKAT' for standard MeerKAT survey
#Set ScanType to 'MeerKATPlusLowScan' for standard MeerKAT survey with additional low elevation
#Set ScanType to 'MeerKATPlusHighScan' for standard MeerKAT survey with additional high elevation
ScanType='MeerKATPlusLowScan'


#Load in the scan information and change az of 360 to be 0 since they are the same
if ScanType == 'MeerKAT':
    ra,dec,az,psi = np.loadtxt('../ScanDataGeneration/TOD.txt')
elif ScanType == 'MeerKATPlusHighScan':
    ra,dec,az,psi = np.loadtxt('../ScanDataGeneration/TOD_extrahighelevs.txt')
elif ScanType == 'MeerKATPlusLowScan':
    ra,dec,az,psi = np.loadtxt('../ScanDataGeneration/TOD_extralowelevs.txt')
az[az==360] = 0.



#Read in input maps from CRIME
crime_nside = 512
start = 1
end = 150
N = end-start+1
gsyncI = np.zeros((N,hp.nside2npix(crime_nside)))
gsyncQ = np.zeros((N,hp.nside2npix(crime_nside)))
gsyncU = np.zeros((N,hp.nside2npix(crime_nside)))
psource = np.zeros((N,hp.nside2npix(crime_nside)))
egfreefree = np.zeros((N,hp.nside2npix(crime_nside)))
gfreefree = np.zeros((N,hp.nside2npix(crime_nside)))
h1cosm = np.zeros((N,hp.nside2npix(crime_nside)))

for i in range(start,end+1):
    if i < 10:
        val = '00'+np.str(i)
    elif i<100:
        val = '0'+np.str(i)
    else:
        val = np.str(i)
    gsyncI[i-start],gsyncQ[i-start],gsyncU[i-start] = hp.read_map('../IM_CRIME_Data/gsync/gsync_'+val+'_b.fits',field=(0,1,2))
    psource[i-start] = hp.read_map('../IM_CRIME_Data/psources/psources_'+val+'.fits')
    egfreefree[i-start] = hp.read_map('../IM_CRIME_Data/egfree/egfree_'+val+'.fits')
    gfreefree[i-start] = hp.read_map('../IM_CRIME_Data/gfree/gfree_'+val+'.fits')
    h1cosm[i-start] = hp.read_map('../IM_CRIME_Data/cosmo/sim_2048_'+val+'.fits')

#Degrade to NSIDE
NSIDE = 128
gsyncI = hp.ud_grade(gsyncI, NSIDE)
gsyncQ = hp.ud_grade(gsyncQ, NSIDE)
gsyncU = hp.ud_grade(gsyncU, NSIDE)
psource = hp.ud_grade(psource, NSIDE)
egfreefree = hp.ud_grade(egfreefree, NSIDE)
gfreefree = hp.ud_grade(gfreefree, NSIDE)
h1cosm = hp.ud_grade(h1cosm, NSIDE)


#Beam smooth map
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



#Set lmax for simulation - make sure it is not too far past the beam scale
ellmax=300
#Calculate the derivative of the intensity and polarization fields
dI=np.zeros(shape=I.shape)
dU=np.zeros(shape=I.shape)
dQ=np.zeros(shape=I.shape)
for i in range(0,len(I)):
    #Spin raising/lowering operators
    #Coordinate convention choice, this can either conjugated or not
    dI[i] = np.conjugate(IMUtils.get_derivativeoffield(I[i],ellmax))
    dQ[i] = np.conjugate(IMUtils.get_derivativeoffield(Q[i],ellmax))
    dU[i] = np.conjugate(IMUtils.get_derivativeoffield(U[i],ellmax))

dP = dQ + 1j *dU
dPconj = dQ - 1j *dU

dbarQ = np.conjugate(dQ)
dbarU = np.conjugate(dU)
dbarP = dbarQ + 1j *dbarU
dbarPconj = dbarQ - 1j *dbarU


#Instantiate arrays for the summed TOD to be put in
tod_A = np.zeros((N,len(ra)))
tod_B = np.zeros((N,len(ra)))


#Array for hitmap
hitmap = np.zeros((N,len(I[0])))


#Set the beam squint level
constantbs = False
linearbs = False
randombs = False
linearbslarge = False

if constantbs:
    rhoa = np.ones(N)*0.4
    rhob = np.zeros(N)
elif randombs:
    seed=1
    np.random.seed(seed)
    rhoa = np.random.random(N)*0.6 -0.3
    rhob = np.zeros(N)
elif linearbs:
    rhoa = np.linspace(-0.2,0.6,N)
    rhob = np.zeros(N)
elif linearbslarge:
    rhoa = 10.*np.linspace(-0.2,0.6,N)
    rhob = np.zeros(N)
else:
    rhoa=np.zeros(N)
    rhob=np.zeros(N)

chi_A = 0.0
chi_B = 0.0
rho_A=np.deg2rad(rhoa/60.)
rho_B=np.deg2rad(rhob/60.)


#Loop over scan data and perform TOD sim
#Find pixel in map based on input ra and dec from scan
pixel_index = hp.ang2pix(NSIDE,dec,ra)

for i in range(0,len(I)):
    #Define the sky signal for the two feeds +
    Sky_theta_A = (1.*I[i][pixel_index] + Q[i][pixel_index] *np.cos(2*psi) + U[i][pixel_index] *np.sin(2*psi)) +0J
    Sky_theta_A += 0.5*rho_A[i]*dI[i][pixel_index]*(np.cos(psi+chi_A) -1j*np.sin(psi+chi_A)) + 0.5*rho_A[i]*np.conjugate(dI[i][pixel_index])*(np.cos(psi+chi_A) + 1j*np.sin(psi+chi_A)) +0J
    Sky_theta_A += -0.25*rho_A[i]*dP[i][pixel_index]*(np.cos(3.*psi+chi_A) - 1j*np.sin(3.*psi+chi_A)) - 0.25*rho_A[i]*dbarP[i][pixel_index]*(np.cos(psi-chi_A) -1j*np.sin(psi-chi_A)) - 0.25*rho_A[i]*dPconj[i][pixel_index]*(np.cos(psi+chi_A) +1j*np.sin(psi+chi_A)) - 0.25*rho_A[i]*dbarPconj[i][pixel_index]*(np.cos(3.*psi+chi_A) + 1j*np.sin(3.*psi+chi_A))
    
    Sky_theta_B = (1.*I[i][pixel_index] + Q[i][pixel_index] *np.cos(2*(psi+np.pi/2.)) + U[i][pixel_index] *np.sin(2*(psi+np.pi/2.))) +0J
    Sky_theta_B += 0.5*rho_B[i]*dI[i][pixel_index]*(np.cos(psi+chi_B) +1j*np.sin(psi+chi_B)) + 0.5*rho_B[i]*np.conjugate(dI[i][pixel_index])*(np.cos(psi+chi_B) -1j*np.sin(psi+chi_B)) +0J
    Sky_theta_B += -0.25*rho_B[i]*dP[i][pixel_index]*(np.cos(3.*psi+chi_B) +1j*np.sin(3.*psi+chi_B)) - 0.25*rho_B[i]*dbarP[i][pixel_index]*(np.cos(psi-chi_B) +1j*np.sin(psi-chi_B)) - 0.25*rho_B[i]*dPconj[i][pixel_index]*(np.cos(-psi+chi_B) +1j*np.sin(-psi+chi_B)) - 0.25*rho_B[i]*dbarPconj[i][pixel_index]*(np.cos(3.*psi+chi_B) -1j*np.sin(3.*psi+chi_B))
    
    
    tod_A[i] = 1.*Sky_theta_A
    tod_B[i] = 1.*Sky_theta_B
    
    np.add.at(hitmap[i],pixel_index,1)
    print(i)


#Sum the detector pairs and add them to the appropriate TOD array
sum_AB = 0.5*(tod_A + tod_B)



mask = hitmap/hitmap


#Calculate Spin-0 Map-making
T_map = np.zeros((N,1,len(I[0])))
pixcond_T = np.zeros((N,len(I[0])))
for i in range(0,len(I)):
    T_map[i],pixcond_T[i] = IMUtils.mapmake_simplebin_spink(T_map[i], pixcond_T[i], sum_AB[i], pixel_index.astype(int),psi,NSIDE,spins=np.array([0,]),mask=None)



#Calculate Spin-0 and 1 Map-making
IQU_Out = np.zeros((N,3,len(I[0])))
pixcond_IQU = np.zeros((N,len(I[0])))
for i in range(0,len(I)):
    IQU_Out[i],pixcond_IQU[i] = IMUtils.mapmake_simplebin_spink(IQU_Out[i], pixcond_IQU[i], sum_AB[i], pixel_index,psi,NSIDE,spins=np.array([0,1]),mask=None)



#Calculate spin-0 and 3
IQU_Out03 = np.zeros((N,3,len(I[0])))
pixcond_IQU03 = np.zeros((N,len(I[0])))
for i in range(0,len(I)):
    IQU_Out03[i],pixcond_IQU03[i] = IMUtils.mapmake_simplebin_spink(IQU_Out03[i], pixcond_IQU03[i], sum_AB[i], pixel_index,psi,NSIDE,spins=np.array([0,3]),mask=None)


#Calculate spin-0 and 1 and 3
IQU_Out013 = np.zeros((N,5,len(I[0])))
pixcond_IQU013 = np.zeros((N,len(I[0])))
for i in range(0,len(I)):
    IQU_Out013[i],pixcond_IQU013[i] = IMUtils.mapmake_simplebin_spink(IQU_Out013[i], pixcond_IQU013[i], sum_AB[i], pixel_index,psi,NSIDE,spins=np.array([0,1,3]),mask=None)



#Set the condition number limit for map-making to remove ill-conditoned pixels.
connum = 1e11
T_map[:,0,:][pixcond_T>connum] = np.nan
for i in range(0,len(IQU_Out)):
    for j in range(0,3):
        IQU_Out[i][j][pixcond_IQU[i]>connum] = np.nan
        IQU_Out013[i][j][pixcond_IQU013[i]>connum] = np.nan
        IQU_Out03[i][j][pixcond_IQU03[i]>connum] = np.nan

Ispin0mapmaking = T_map[:,0]
Ispin01mapmaking = IQU_Out[:,0]
Ispin03mapmaking = IQU_Out03[:,0]
Ispin013mapmaking = IQU_Out013[:,0]


'''
plt.figure()
hp.mollview(Ispin0mapmaking[142]-mask[142]*I[142],title='Spin-0 Map-Making 915.5MHz Residual Out-In')
hp.mollview(Ispin01mapmaking[142]-mask[142]*I[142],title='Spin-0,1 Map-Making 915.5MHz Residual Out-In')
hp.mollview(Ispin013mapmaking[142]-mask[142]*I[142],title='Spin-0,1,3 Map-Making 915.5MHz Residual Out-In')

hp.mollview(Ispin0mapmaking[142]-mask[142]*I[142],title='Spin-0 Map-Making 915.5MHz Residual Out-In',sub=(311))
hp.mollview(Ispin01mapmaking[142]-mask[142]*I[142],title='Spin-0,1 Map-Making 915.5MHz Residual Out-In',sub=(312))
hp.mollview(Ispin013mapmaking[142]-mask[142]*I[142],title='Spin-0,1,3 Map-Making 915.5MHz Residual Out-In',sub=(313))

hp.mollview((Ispin0mapmaking[142]-mask[142]*I[142])/(Ispin0mapmaking[142]-mask[142]*I[142]),title='Spin-0 Map-Making 915.5MHz Mask',sub=(311))
hp.mollview((Ispin01mapmaking[142]-mask[142]*I[142])/(Ispin01mapmaking[142]-mask[142]*I[142]),title='Spin-0,1 Map-Making 915.5MHz Mask',sub=(312))
hp.mollview((Ispin013mapmaking[142]-mask[142]*I[142])/(Ispin013mapmaking[142]-mask[142]*I[142]),title='Spin-0,1,3 Map-Making 915.5MHz Mask',sub=(313))
'''



if constantbs:
    hp.write_map('Ispin0mapmaking_BeamSquint_rhoConstantCondNumCut1e11_'+ScanType+'.fits', Ispin0mapmaking)
    hp.write_map('Ispin0and1mapmaking_BeamSquint_rhoConstantCondNumCut1e11_'+ScanType+'.fits', Ispin01mapmaking)
    hp.write_map('Ispin0and3mapmaking_BeamSquint_rhoConstantCondNumCut1e11_'+ScanType+'.fits', Ispin03mapmaking)
    hp.write_map('Ispin0and1and3mapmaking_BeamSquint_rhoConstantCondNumCut1e11_'+ScanType+'.fits', Ispin013mapmaking)
elif randombs:
    hp.write_map('Ispin0mapmaking_BeamSquint_rhoRandomCondNumCut1e11_'+ScanType+'.fits', Ispin0mapmaking)
    hp.write_map('Ispin0and1mapmaking_BeamSquint_rhoRandomCondNumCut1e11_'+ScanType+'.fits', Ispin01mapmaking)
    hp.write_map('Ispin0and3mapmaking_BeamSquint_rhoRandomCondNumCut1e11_'+ScanType+'.fits', Ispin03mapmaking)
    hp.write_map('Ispin0and1and3mapmaking_BeamSquint_rhoRandomCondNumCut1e11_'+ScanType+'.fits', Ispin013mapmaking)
elif linearbs:
    hp.write_map('Ispin0mapmaking_BeamSquint_rhoLinearCondNumCut1e11_'+ScanType+'.fits', Ispin0mapmaking)
    hp.write_map('Ispin0and1mapmaking_BeamSquint_rhoLinearCondNumCut1e11_'+ScanType+'.fits', Ispin01mapmaking)
    hp.write_map('Ispin0and3mapmaking_BeamSquint_rhoLinearCondNumCut1e11_'+ScanType+'.fits', Ispin03mapmaking)
    hp.write_map('Ispin0and1and3mapmaking_BeamSquint_rhoLinearCondNumCut1e11_'+ScanType+'.fits', Ispin013mapmaking)
elif linearbslarge:
    hp.write_map('Ispin0mapmaking_BeamSquint_rhoLinearLargeCondNumCut1e11_'+ScanType+'.fits', Ispin0mapmaking)
    hp.write_map('Ispin0and1mapmaking_BeamSquint_rhoLinearLargeCondNumCut1e11_'+ScanType+'.fits', Ispin01mapmaking)
    hp.write_map('Ispin0and3mapmaking_BeamSquint_rhoLinearLargeCondNumCut1e11_'+ScanType+'.fits', Ispin03mapmaking)
    hp.write_map('Ispin0and1and3mapmaking_BeamSquint_rhoLinearLargeCondNumCut1e11_'+ScanType+'.fits', Ispin013mapmaking)
else:
    hp.write_map('Ispin0mapmaking_BeamSquint_rho0CondNumCut1e11_'+ScanType+'.fits', Ispin0mapmaking)
    hp.write_map('Ispin0and1mapmaking_BeamSquint_rho0CondNumCut1e11_'+ScanType+'.fits', Ispin01mapmaking)
    hp.write_map('Ispin0and3mapmaking_BeamSquint_rho0CondNumCut1e11_'+ScanType+'.fits', Ispin03mapmaking)
    hp.write_map('Ispin0and1and3mapmaking_BeamSquint_rho0CondNumCut1e11_'+ScanType+'.fits', Ispin01mapmaking)
