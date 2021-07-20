#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 10:27:04 2021

@author: mccallum
"""

import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

import sys
#Set path to where gmca4im is installed
sys.path.insert(1, '/mirror/scratch/mccallum/IntensityMapping/ForegroundCleaning/gmca4im/scripts')
import gmca4im_lib2 as g4i
sys.path.append("../")
import IMUtils


#Set the frequency channel range and load relevant data
channel1=2
channel2=148
n_freq, nu_0, nu_f, z_0, z_f = np.loadtxt('../IM_CRIME_Data/nuTable.txt').T
freqs = (nu_f+nu_0)/2

freqs=freqs[channel1:channel2]


#Set the gain levels used in the TOD sims
g_A=0.01
g_B=0.0
T_calibration = 1. + (g_A+g_B)/2.

#Choose which scan type
ScanType='MeerKAT'

#Load the TOD maps
obs_maps_spin0and2mapmake = hp.read_map('../TODSimulations/Ispin0and2mapmaking_Gain_g1a'+np.str(g_A)+'g1b'+np.str(g_B)+'CondNumCut1e11_'+ScanType+'.fits',field=np.arange(channel1,channel2))
obs_maps_spin0mapmake = hp.read_map('../TODSimulations/Ispin0mapmaking_Gain_g1a'+np.str(g_A)+'g1b'+np.str(g_B)+'CondNumCut1e11_'+ScanType+'.fits',field=np.arange(channel1,channel2))

obs_maps_spin0and2mapmake[np.isnan(obs_maps_spin0and2mapmake)]=0.
obs_maps_spin0mapmake[np.isnan(obs_maps_spin0mapmake)]=0.

#Remove the Temperature amplification - we assumed perfect calibration
obs_maps_spin0and2mapmake/=T_calibration
obs_maps_spin0mapmake/=T_calibration

#Load input HI signal
#(Make sure to load the one which is the same NSIDE and beam smeared as the TOD)
h1cosm = hp.read_map('../TODSimulations/h1cosm.fits',field=np.arange(channel1,channel2))

#Load no systematics TOD maps
obs_maps_spin0and2nosyst = hp.read_map('../TODSimulations/Ispin0and2mapmaking_Gain_g1a0.0g1b0.0CondNumCut1e11_'+ScanType+'.fits',field=np.arange(channel1,channel2))
obs_maps_spin0nosyst = hp.read_map('../TODSimulations/Ispin0mapmaking_Gain_g1a0.0g1b0.0CondNumCut1e11_'+ScanType+'.fits',field=np.arange(channel1,channel2))
obs_maps_spin0and2nosyst[np.isnan(obs_maps_spin0and2nosyst)]=0.
obs_maps_spin0nosyst[np.isnan(obs_maps_spin0nosyst)]=0.


#Make mask for data - this will be limited by the least well conditioned map-making
#Makes for easy comparison as we will use same mask on all data
mask = (obs_maps_spin0and2mapmake/obs_maps_spin0and2mapmake)*(obs_maps_spin0mapmake/obs_maps_spin0mapmake)*(obs_maps_spin0nosyst/obs_maps_spin0nosyst)*(obs_maps_spin0nosyst/obs_maps_spin0nosyst)
mask[np.isnan(mask)]=0.


#Mask the data
obs_maps_spin0and2nosyst*=mask
obs_maps_spin0nosyst*=mask
obs_maps_spin0and2mapmake*=mask
obs_maps_spin0mapmake*=mask


#Remove the mean of the masked field from the data (see gmca4im for details why)
for i in range(0,len(obs_maps_spin0and2mapmake)):
    obs_maps_spin0and2mapmake[i][obs_maps_spin0and2mapmake[i]!=0] = obs_maps_spin0and2mapmake[i][obs_maps_spin0and2mapmake[i]!=0] - np.mean(obs_maps_spin0and2mapmake[i][obs_maps_spin0and2mapmake[i]!=0])
    obs_maps_spin0mapmake[i][obs_maps_spin0mapmake[i]!=0] = obs_maps_spin0mapmake[i][obs_maps_spin0mapmake[i]!=0] - np.mean(obs_maps_spin0mapmake[i][obs_maps_spin0mapmake[i]!=0])
    
    obs_maps_spin0and2nosyst[i][obs_maps_spin0and2nosyst[i]!=0] = obs_maps_spin0and2nosyst[i][obs_maps_spin0and2nosyst[i]!=0] - np.mean(obs_maps_spin0and2nosyst[i][obs_maps_spin0and2nosyst[i]!=0])
    obs_maps_spin0nosyst[i][obs_maps_spin0nosyst[i]!=0] = obs_maps_spin0nosyst[i][obs_maps_spin0nosyst[i]!=0] - np.mean(obs_maps_spin0nosyst[i][obs_maps_spin0nosyst[i]!=0])



#Generate or load the wavelets
genwavelets=False
if genwavelets:
    X_wt_spin0and2mapmake = g4i.wavelet_transform(obs_maps_spin0and2mapmake)
    X_wt_spin0mapmake = g4i.wavelet_transform(obs_maps_spin0mapmake)
    
    X_wt_spin0and2nosyst = g4i.wavelet_transform(obs_maps_spin0and2nosyst)
    X_wt_spin0nosyst = g4i.wavelet_transform(obs_maps_spin0nosyst)
else:
    X_wt_spin0and2mapmake = np.loadtxt('./WaveletsIspin0and2mapmaking_Gain_g1a0.01g1b0.0CondNumCut1e11.txt')
    X_wt_spin0mapmake = np.loadtxt('./WaveletsIspin0mapmaking_Gain_g1a0.01g1b0.0CondNumCut1e11.txt')

    X_wt_spin0and2nosyst = np.loadtxt('./WaveletsIspin0and2mapmaking_Gain_g1a0.0g1b0.0CondNumCut1e11.txt')
    X_wt_spin0nosyst = np.loadtxt('./WaveletsIspin0and2mapmaking_Gain_g1a0.0g1b0.0CondNumCut1e11.txt')



#Save the wavelets
savewavelets = False
if savewavelets:
    np.savetxt('WaveletsIspin0and2mapmaking_Gain_g1a0.01g1b0.0CondNumCut1e11.txt',X_wt_spin0and2mapmake)
    np.savetxt('WaveletsIspin0mapmaking_Gain_g1a0.01g1b0.0CondNumCut1e11.txt',X_wt_spin0mapmake)
    np.savetxt('WaveletsIspin0and2mapmaking_Gain_g1a0.0g1b0.0CondNumCut1e11.txt',X_wt_spin0and2nosyst)
    np.savetxt('WaveletsIspin0mapmaking_Gain_g1a0.0g1b0.0CondNumCut1e11.txt',X_wt_spin0nosyst)



#Perform the GMCA
#####GMCA PARAMETERS#####
n_s   = 5  # number of sources to estimate
mints = 0.1 # minimum threshold
nmax  = 100 # number of iterations (usually 100 is safe)
L0    = 0   # L0 norm (1) or L1 norm (0)


#Initial guess for the mixing matrix
AInit = None
ColFixed = None

#Choose whether to whiten the data
whitening = False; epsi = 1e-3

#Estimate mixing matrix
Ae_spin0and2mapmake = g4i.run_GMCA(X_wt_spin0and2mapmake,AInit,n_s,mints,nmax,L0,ColFixed,whitening,epsi)
Ae_spin0mapmake = g4i.run_GMCA(X_wt_spin0mapmake,AInit,n_s,mints,nmax,L0,ColFixed,whitening,epsi)

Ae_spin0and2nosyst = g4i.run_GMCA(X_wt_spin0and2nosyst,AInit,n_s,mints,nmax,L0,ColFixed,whitening,epsi)
Ae_spin0nosyst = g4i.run_GMCA(X_wt_spin0nosyst,AInit,n_s,mints,nmax,L0,ColFixed,whitening,epsi)

#Reconstructed maps by GMCA
piA_spin0and2mapmake = np.linalg.inv(Ae_spin0and2mapmake.T@Ae_spin0and2mapmake)@Ae_spin0and2mapmake.T
#Reproject onto original map coordinates
Se_sph_spin0and2mapmake = piA_spin0and2mapmake@obs_maps_spin0and2mapmake # LS estimate of the sources in the pixel domain
X_gmca_spin0and2mapmake = Ae_spin0and2mapmake@Se_sph_spin0and2mapmake; del Se_sph_spin0and2mapmake, piA_spin0and2mapmake

#Reconstructed maps by GMCA
piA_spin0mapmake = np.linalg.inv(Ae_spin0mapmake.T@Ae_spin0mapmake)@Ae_spin0mapmake.T
#Reproject onto original map coordinates
Se_sph_spin0mapmake = piA_spin0mapmake@obs_maps_spin0mapmake # LS estimate of the sources in the pixel domain
X_gmca_spin0mapmake = Ae_spin0mapmake@Se_sph_spin0mapmake; del Se_sph_spin0mapmake, piA_spin0mapmake

#Reconstructed maps by GMCA
piA_spin0and2nosyst = np.linalg.inv(Ae_spin0and2nosyst.T@Ae_spin0and2nosyst)@Ae_spin0and2nosyst.T
#Reproject onto original map coordinates
Se_sph_spin0and2nosyst = piA_spin0and2nosyst@obs_maps_spin0and2nosyst # LS estimate of the sources in the pixel domain
X_gmca_spin0and2nosyst = Ae_spin0and2nosyst@Se_sph_spin0and2nosyst; del Se_sph_spin0and2nosyst, piA_spin0and2nosyst

#Reconstructed maps by GMCA
piA_spin0nosyst = np.linalg.inv(Ae_spin0nosyst.T@Ae_spin0nosyst)@Ae_spin0nosyst.T
#Reproject onto original map coordinates
Se_sph_spin0nosyst = piA_spin0nosyst@obs_maps_spin0nosyst # LS estimate of the sources in the pixel domain
X_gmca_spin0nosyst = Ae_spin0nosyst@Se_sph_spin0nosyst; del Se_sph_spin0nosyst, piA_spin0nosyst





mask = mask[0]

#Pick a channel number to examine 140 is ~915.5MHz
ich = 140



#Map Space
map_input = 1.*h1cosm[ich]
mask[np.isnan(mask)]=0.
map_input*=mask
map_input[np.isnan(map_input)]=0.
map_input[map_input!=0] -= np.mean(map_input[map_input!=0])



#Choose whether to apodize. Don't really need to in this case.
Apodize=False
if Apodize:
    mask = IMUtils.ApodizeMask(mask,2.0)


#Mask the input map
map_input*=mask


#Subtract the GMCA calculated foregrounds from input maps giving residual maps
residuals_spin0and2mapmake = obs_maps_spin0and2mapmake[ich]-X_gmca_spin0and2mapmake[ich]
residuals_spin0mapmake = obs_maps_spin0mapmake[ich]-X_gmca_spin0mapmake[ich]

residuals_spin0and2nosyst = obs_maps_spin0and2nosyst[ich]-X_gmca_spin0and2nosyst[ich]
residuals_spin0nosyst = obs_maps_spin0nosyst[ich]-X_gmca_spin0nosyst[ich]


residuals_spin0and2mapmake*=mask
residuals_spin0mapmake*=mask
residuals_spin0and2nosyst*=mask
residuals_spin0nosyst*=mask


#Settings for cartview size of field to look at
lonra=[-22.5,17.5]
latra=[-6,6]

fig = plt.figure(figsize=(20, 6))
fig.suptitle('Channel = '+np.str(np.round(freqs[ich],2))+'MHz Spin 0 (and 2) Map-Making Lower (Upper)',fontsize=20)
ax1 = fig.add_subplot(2,3,1)
plt.axes(ax1)
input1 = hp.cartview(map_input*(mask/mask),title='Input Cosmological Signal',unit=r'$T$ [mK]',hold=True,rot=(162.5,3,0),lonra=lonra,latra=latra,return_projected_map=True)
ax2 = fig.add_subplot(2,3,2)
plt.axes(ax2)
output02 = hp.cartview(residuals_spin0and2mapmake*(mask/mask),title='GMCA residuals, Num Components = 5',unit=r'$T$ [mK]',hold=True,rot=(162.5,3,0),lonra=lonra,latra=latra,return_projected_map=True)
ax3 = fig.add_subplot(2,3,3)
plt.axes(ax3)
output02_ns = hp.cartview(residuals_spin0and2nosyst*(mask/mask),title='GMCA residuals, Num Components = 5 No Systematic',unit=r'$T$ [mK]',hold=True,rot=(162.5,3,0),lonra=lonra,latra=latra,return_projected_map=True)


#fig = plt.figure(figsize=(20, 6))
#fig.suptitle('Channel = '+np.str(np.round(freqs[ich],2))+'MHz Spin 0 Map-Making',fontsize=20)
ax4 = fig.add_subplot(2,3,4)
plt.axes(ax4)
#hp.cartview(map_input,title='Cosmological Signal',unit=r'$T$ [mK]',hold=True,rot=(162.5,3,0),lonra=lonra,latra=latra)
ax5 = fig.add_subplot(2,3,5)
plt.axes(ax5)
output0 = hp.cartview(residuals_spin0mapmake*(mask/mask),title='GMCA residuals, Num Components = 5',unit=r'$T$ [mK]',hold=True,rot=(162.5,3,0),lonra=lonra,latra=latra,return_projected_map=True)
ax6 = fig.add_subplot(2,3,6)
plt.axes(ax6)
output0_ns = hp.cartview(residuals_spin0nosyst*(mask/mask),title='GMCA residuals, Num Components = 5 No Systematic',unit=r'$T$ [mK]',hold=True,rot=(162.5,3,0),lonra=lonra,latra=latra,return_projected_map=True)




###Paper map plots
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 20})
left=lonra[1]
right=lonra[0]

fig, (ax1) = plt.subplots(ncols=1,nrows=1,figsize=(4, 3))
#fig, (ax1) = plt.subplots(ncols=1,nrows=1)
#fig.suptitle('Channel = '+np.str(np.round(freqs[ich],2))+'MHz Spin 0 (and 2) Map-Making Lower (Upper)',fontsize=20)
plt.axes(ax1)
#plt.title('Input Cosmological Signal')
im1 = ax1.imshow(input1, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none',vmin=-0.1,vmax=0.28)
plt.xlabel('R.A.')
plt.ylabel('Dec.')
plt.xticks(np.arange(175,135,-5))
clb=fig.colorbar(im1, orientation="horizontal", pad=0.2)
clb.set_label('mK')

plt.rcParams.update({'font.size': 20})
fig, (ax1,ax2) = plt.subplots(ncols=2,nrows=1)
#fig.suptitle('Channel = '+np.str(np.round(freqs[ich],2))+'MHz Spin 0 (and 2) Map-Making Lower (Upper)',fontsize=20)
plt.axes(ax1)
#plt.title('GMCA residuals, Num Components = 5, Spin-0 and 2')
im2 = ax1.imshow(output02, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none',vmin=-0.1,vmax=0.28)
plt.xlabel('R.A.')
plt.ylabel('Dec.')
plt.xticks(np.arange(175,135,-5))
clb=fig.colorbar(im2, orientation="horizontal", pad=0.2)
clb.set_label('mK')

plt.axes(ax2)
#plt.title('GMCA residuals, Num Components = 5, Spin-0')
im5 = ax2.imshow(output0, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none')
plt.xlabel('R.A.')
plt.ylabel('Dec.')
plt.xticks(np.arange(175,135,-5))
clb=fig.colorbar(im5, orientation="horizontal", pad=0.2)
clb.set_label('mK')


fig, (ax1,ax2) = plt.subplots(ncols=2,nrows=1)
plt.axes(ax1)
#plt.title('GMCA residuals, Num Components = 5, Spin-0 and 2, No Systematic')
im3 = ax1.imshow(output02_ns, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none',vmin=-0.1,vmax=0.28)
plt.xlabel('R.A.')
plt.ylabel('Dec.')
clb=fig.colorbar(im3, orientation="horizontal", pad=0.2)
clb.set_label('mK')

plt.axes(ax2)
#plt.title('GMCA residuals, Num Components = 5, Spin-0, No Systematic')
im6 = ax2.imshow(output0_ns, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none',vmin=-0.1,vmax=0.28)
plt.xlabel('R.A.')
plt.ylabel('Dec.')
clb=fig.colorbar(im6, orientation="horizontal", pad=0.2)
clb.set_label('mK')


del residuals_spin0and2mapmake, residuals_spin0mapmake






#C_l space
residuals_spin0and2mapmake = mask*(obs_maps_spin0and2mapmake[ich]-X_gmca_spin0and2mapmake[ich])
residuals_spin0mapmake = mask*(obs_maps_spin0mapmake[ich]-X_gmca_spin0mapmake[ich])

residuals_spin0and2nosyst = mask*(obs_maps_spin0and2nosyst[ich]-X_gmca_spin0and2nosyst[ich])
residuals_spin0nosyst = mask*(obs_maps_spin0nosyst[ich]-X_gmca_spin0nosyst[ich])


plt.figure()
#plt.rcParams["figure.figsize"] = (8,6)
plt.rcParams["axes.labelsize"] = 16
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 16})
ell, y = g4i.plot_cl(map_input)
plt.semilogy(ell, y,'-',c='k',label='Cosmological Signal');


ellnosyst, ynosyst = g4i.plot_cl(residuals_spin0and2nosyst)
plt.semilogy(ellnosyst,ynosyst,'o',c='r',label='Spin-0 and 2 map-making, No Systematic')

ell_spin0nosyst, y_spin0nosyst = g4i.plot_cl(residuals_spin0nosyst)
plt.semilogy(ell_spin0nosyst,y_spin0nosyst,'.',c='m',label='Spin-0 map-making, No Systematic')


ell, y = g4i.plot_cl(residuals_spin0and2mapmake)
plt.semilogy(ell,y,'+',c='b',label='Spin-0 and 2 map-making, Effective Gain Mismatch Included')

ell_spin0mapmake, y_spin0mapmake = g4i.plot_cl(residuals_spin0mapmake)
plt.semilogy(ell_spin0mapmake,y_spin0mapmake,'x',c='g',label='Spin-0 map-making, Effective Gain Mismatch Included')


#plt.title('Channel = '+np.str(np.round(freqs[ich],2))+'MHz')
plt.legend()
ax = plt.gca()
ax.set(xlim=[15,200],xlabel="$\\ell$",ylabel="$\\ell(\\ell+1)C_{\\ell}/(2\\pi)$ [mK$^2$]");





#K Space
map_input = 1.*h1cosm
mask[np.isnan(mask)]=0.
map_input*=(mask/mask)
map_input[map_input==0] = np.nan
map_input -= np.array([np.nanmean(map_input,axis=1)]).T
map_input[np.isnan(map_input)]=0.
map_input*=mask


NSIDE=128
pix = np.arange(0,hp.nside2npix(NSIDE))
indexes_los = pix[mask!=0]



residuals_syst02 = obs_maps_spin0and2mapmake-X_gmca_spin0and2mapmake
residuals_syst0 = obs_maps_spin0mapmake-X_gmca_spin0mapmake
residuals_nosyst02 = obs_maps_spin0and2nosyst-X_gmca_spin0and2nosyst
residuals_nosyst0 = obs_maps_spin0nosyst-X_gmca_spin0nosyst

residuals_syst02*=mask
residuals_syst0*=mask
residuals_nosyst02*=mask
residuals_nosyst0*=mask

plt.figure()
plt.rcParams["figure.figsize"] = (8,6)
plt.rcParams["axes.labelsize"] = 16
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 16})

k, P = g4i.plot_nuPk(map_input,indexes_los,freqs)
plt.semilogy(k, P,'-',c='k',label='Cosmological Signal')

k, P = g4i.plot_nuPk(residuals_nosyst02,indexes_los,freqs)
plt.semilogy(k, P,'o',c='r',label='Spin-0 and 2 map-making, No Systematic')
         
k, P = g4i.plot_nuPk(residuals_nosyst0,indexes_los,freqs)
plt.semilogy(k, P,'.',c='m',label='Spin-0 map-making, No Systematic')

k, P = g4i.plot_nuPk(residuals_syst02,indexes_los,freqs)
plt.semilogy(k, P,'+',c='b',label='Spin-0 and 2 map-making, Effective Gain Mismatch Included')

k, P = g4i.plot_nuPk(residuals_syst0,indexes_los,freqs)
plt.semilogy(k, P,'x',c='g',label='Spin-0 map-making, Effective Gain Mismatch Included')
         

plt.legend(fontsize=16)
ax = plt.gca()
ax.set(xlabel="$k_{\\nu}$ [MHz$^{-1}$]",ylabel="$P$ [mK$^2$ MHz$^2$]");





