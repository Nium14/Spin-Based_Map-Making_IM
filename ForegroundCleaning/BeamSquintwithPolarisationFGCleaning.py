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



#Choose which scan type
ScanType='MeerKATPlusLowScan'


#Load the TOD maps
constantbs=False
randombs=False
linearbs=True
linearbslarge=False

if constantbs:
    obs_maps_spin0mapmake = hp.read_map('../TODSimulations/Ispin0mapmaking_BeamSquint_rhoConstantCondNumCut1e11_'+ScanType+'.fits',field=np.arange(channel1,channel2))
    obs_maps_spin0and1and3mapmake = hp.read_map('../TODSimulations/Ispin0and1and3mapmaking_BeamSquint_rhoConstantCondNumCut1e11_'+ScanType+'.fits',field=np.arange(channel1,channel2))
if randombs:
    obs_maps_spin0mapmake = hp.read_map('../TODSimulations/Ispin0mapmaking_BeamSquint_rhoRandomCondNumCut1e11_'+ScanType+'.fits',field=np.arange(channel1,channel2))
    obs_maps_spin0and1and3mapmake = hp.read_map('../TODSimulations/Ispin0and1and3mapmaking_BeamSquint_rhoRandomCondNumCut1e11_'+ScanType+'.fits',field=np.arange(channel1,channel2))
if linearbs:
    obs_maps_spin0mapmake = hp.read_map('../TODSimulations/Ispin0mapmaking_BeamSquint_rhoLinearCondNumCut1e11_'+ScanType+'.fits',field=np.arange(channel1,channel2))
    obs_maps_spin0and1and3mapmake = hp.read_map('../TODSimulations/Ispin0and1and3mapmaking_BeamSquint_rhoLinearCondNumCut1e11_'+ScanType+'.fits',field=np.arange(channel1,channel2))
if linearbslarge:
    obs_maps_spin0mapmake = hp.read_map('../TODSimulations/Ispin0mapmaking_BeamSquint_rhoLinearLargeCondNumCut1e11_'+ScanType+'.fits',field=np.arange(channel1,channel2))
    obs_maps_spin0and1and3mapmake = hp.read_map('../TODSimulations/Ispin0and1and3mapmaking_BeamSquint_rhoLinearLargeCondNumCut1e11_'+ScanType+'.fits',field=np.arange(channel1,channel2))


obs_maps_spin0and1and3mapmake[np.isnan(obs_maps_spin0and1and3mapmake)]=0.
obs_maps_spin0mapmake[np.isnan(obs_maps_spin0mapmake)]=0.


#Load input HI signal
#(Make sure to load the one which is the same NSIDE and beam smeared as the TOD)
h1cosm = hp.read_map('../TODSimulations/h1cosm.fits',field=np.arange(channel1,channel2))


#Load no systematics TOD maps
obs_maps_spin0and1and3nosyst = hp.read_map('../TODSimulations/Ispin0and1mapmaking_BeamSquint_rho0CondNumCut1e11_'+ScanType+'.fits',field=np.arange(channel1,channel2))
obs_maps_spin0nosyst = hp.read_map('../TODSimulations/Ispin0mapmaking_BeamSquint_rho0CondNumCut1e11_'+ScanType+'.fits',field=np.arange(channel1,channel2))
obs_maps_spin0and1and3nosyst[np.isnan(obs_maps_spin0and1and3nosyst)]=0.
obs_maps_spin0nosyst[np.isnan(obs_maps_spin0nosyst)]=0.


#Make mask for data - this will be limited by the least well conditioned map-making
#Makes for easy comparison as we will use same mask on all data
mask = (obs_maps_spin0and1and3mapmake/obs_maps_spin0and1and3mapmake)*(obs_maps_spin0mapmake/obs_maps_spin0mapmake)*(obs_maps_spin0nosyst/obs_maps_spin0nosyst)*(obs_maps_spin0and1and3nosyst/obs_maps_spin0and1and3nosyst)
mask[np.isnan(mask)]=0.


#Mask the data
obs_maps_spin0and1and3nosyst*=mask
obs_maps_spin0nosyst*=mask
obs_maps_spin0and1and3mapmake*=mask
obs_maps_spin0mapmake*=mask


#Remove the mean of the masked field from the data (see gmca4im for details why)
for i in range(0,len(obs_maps_spin0and1and3mapmake)):
    obs_maps_spin0and1and3mapmake[i][obs_maps_spin0and1and3mapmake[i]!=0] = obs_maps_spin0and1and3mapmake[i][obs_maps_spin0and1and3mapmake[i]!=0] - np.mean(obs_maps_spin0and1and3mapmake[i][obs_maps_spin0and1and3mapmake[i]!=0])
    obs_maps_spin0mapmake[i][obs_maps_spin0mapmake[i]!=0] = obs_maps_spin0mapmake[i][obs_maps_spin0mapmake[i]!=0] - np.mean(obs_maps_spin0mapmake[i][obs_maps_spin0mapmake[i]!=0])
    
    obs_maps_spin0and1and3nosyst[i][obs_maps_spin0and1and3nosyst[i]!=0] = obs_maps_spin0and1and3nosyst[i][obs_maps_spin0and1and3nosyst[i]!=0] - np.mean(obs_maps_spin0and1and3nosyst[i][obs_maps_spin0and1and3nosyst[i]!=0])
    obs_maps_spin0nosyst[i][obs_maps_spin0nosyst[i]!=0] = obs_maps_spin0nosyst[i][obs_maps_spin0nosyst[i]!=0] - np.mean(obs_maps_spin0nosyst[i][obs_maps_spin0nosyst[i]!=0])



#Generate or load the wavelets
genwavelets=True
if genwavelets:
    X_wt_spin0and1and3mapmake = g4i.wavelet_transform(obs_maps_spin0and1and3mapmake)
    #X_wt_spin0and1mapmake = g4i.wavelet_transform(obs_maps_spin0and1mapmake)
    X_wt_spin0mapmake = g4i.wavelet_transform(obs_maps_spin0mapmake)

    X_wt_spin0and1and3nosyst = g4i.wavelet_transform(obs_maps_spin0and1and3nosyst)
    #X_wt_spin0and1nosyst = g4i.wavelet_transform(obs_maps_spin0and1nosyst)
    X_wt_spin0nosyst = g4i.wavelet_transform(obs_maps_spin0nosyst)
else:
    X_wt_spin0and1and3mapmake = np.loadtxt('./WaveletsIspin0and1mapmaking_BeamSquint_rhoLinearCondNumCut1e11.txt')
    X_wt_spin0mapmake = np.loadtxt('./WaveletsIspin0mapmaking_BeamSquint_rhoLinearCondNumCut1e11.txt')

    X_wt_spin0and1and3nosyst = np.loadtxt('./WaveletsIspin0and1mapmaking_BeamSquint_rho0CondNumCut1e11.txt')
    X_wt_spin0nosyst = np.loadtxt('./WaveletsIspin0mapmaking_BeamSquint_rho0CondNumCut1e11.txt')


#Save the wavelets
savewavelets = True
if savewavelets:
    np.savetxt('WaveletsIspin0and1mapmaking_BeamSquint_rhoRandomCondNumCut1e11.txt',X_wt_spin0and1and3mapmake)
    np.savetxt('WaveletsIspin0mapmaking_BeamSquint_rhoRandomCondNumCut1e11.txt',X_wt_spin0mapmake)
    np.savetxt('WaveletsIspin0and1mapmaking_BeamSquint_rho0CondNumCut1e11.txt',X_wt_spin0and1and3nosyst)
    np.savetxt('WaveletsIspin0mapmaking_BeamSquint_rho0CondNumCut1e11.txt',X_wt_spin0nosyst)


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
Ae_spin0and1and3mapmake = g4i.run_GMCA(X_wt_spin0and1and3mapmake,AInit,n_s,mints,nmax,L0,ColFixed,whitening,epsi)
Ae_spin0mapmake = g4i.run_GMCA(X_wt_spin0mapmake,AInit,n_s,mints,nmax,L0,ColFixed,whitening,epsi)

Ae_spin0and1and3nosyst = g4i.run_GMCA(X_wt_spin0and1and3nosyst,AInit,n_s,mints,nmax,L0,ColFixed,whitening,epsi)
Ae_spin0nosyst = g4i.run_GMCA(X_wt_spin0nosyst,AInit,n_s,mints,nmax,L0,ColFixed,whitening,epsi)


#Reconstructed maps by GMCA
piA_spin0and1and3mapmake = np.linalg.inv(Ae_spin0and1and3mapmake.T@Ae_spin0and1and3mapmake)@Ae_spin0and1and3mapmake.T
#Reproject onto original map coordinates
Se_sph_spin0and1and3mapmake = piA_spin0and1and3mapmake@obs_maps_spin0and1and3mapmake # LS estimate of the sources in the pixel domain
X_gmca_spin0and1and3mapmake = Ae_spin0and1and3mapmake@Se_sph_spin0and1and3mapmake; del Se_sph_spin0and1and3mapmake, piA_spin0and1and3mapmake

#Reconstructed maps by GMCA
piA_spin0mapmake = np.linalg.inv(Ae_spin0mapmake.T@Ae_spin0mapmake)@Ae_spin0mapmake.T
#Reproject onto original map coordinates
Se_sph_spin0mapmake = piA_spin0mapmake@obs_maps_spin0mapmake # LS estimate of the sources in the pixel domain
X_gmca_spin0mapmake = Ae_spin0mapmake@Se_sph_spin0mapmake; del Se_sph_spin0mapmake, piA_spin0mapmake

#Reconstructed maps by GMCA
piA_spin0and1and3nosyst = np.linalg.inv(Ae_spin0and1and3nosyst.T@Ae_spin0and1and3nosyst)@Ae_spin0and1and3nosyst.T
#Reproject onto original map coordinates
Se_sph_spin0and1and3nosyst = piA_spin0and1and3nosyst@obs_maps_spin0and1and3nosyst # LS estimate of the sources in the pixel domain
X_gmca_spin0and1and3nosyst = Ae_spin0and1and3nosyst@Se_sph_spin0and1and3nosyst; del Se_sph_spin0and1and3nosyst, piA_spin0and1and3nosyst

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



#Mask the input map
map_input*=mask


#Subtract the GMCA calculated foregrounds from input maps giving residual maps
residuals_spin0and1and3mapmake = obs_maps_spin0and1and3mapmake[ich]-X_gmca_spin0and1and3mapmake[ich]
residuals_spin0mapmake = obs_maps_spin0mapmake[ich]-X_gmca_spin0mapmake[ich]

residuals_spin0and1and3nosyst = obs_maps_spin0and1and3nosyst[ich]-X_gmca_spin0and1and3nosyst[ich]
residuals_spin0nosyst = obs_maps_spin0nosyst[ich]-X_gmca_spin0nosyst[ich]

residuals_spin0and1and3mapmake*=mask
residuals_spin0mapmake*=mask
residuals_spin0and1and3nosyst*=mask
residuals_spin0nosyst*=mask


#Settings for cartview size of field to look at
lonra=[-22.5,17.5]
latra=[-6,6]


fig = plt.figure(figsize=(20, 6))
fig.suptitle('Channel = '+np.str(np.round(freqs[ich],2))+'MHz Spin 0 (and 1 and 3) Map-Making Lower (Upper)',fontsize=20)
ax1 = fig.add_subplot(2,3,1)
plt.axes(ax1)
input1 = hp.cartview(map_input*(mask/mask),title='Input Cosmological Signal',unit=r'$T$ [mK]',hold=True,rot=(162.5,3,0),lonra=lonra,latra=latra,return_projected_map=True)
ax2 = fig.add_subplot(2,3,2)
plt.axes(ax2)
output013 = hp.cartview(residuals_spin0and1and3mapmake*(mask/mask),title='GMCA residuals, Num Components = 5',unit=r'$T$ [mK]',hold=True,rot=(162.5,3,0),lonra=lonra,latra=latra,return_projected_map=True)
ax3 = fig.add_subplot(2,3,3)
plt.axes(ax3)
output013_ns = hp.cartview(residuals_spin0and1and3nosyst*(mask/mask),title='GMCA residuals, Num Components = 5 No Systematic',unit=r'$T$ [mK]',hold=True,rot=(162.5,3,0),lonra=lonra,latra=latra,return_projected_map=True)


#fig = plt.figure(figsize=(20, 6))
#fig.suptitle('Channel = '+np.str(np.round(freqs[ich],2))+'MHz Spin 0 Map-Making',fontsize=20)
ax4 = fig.add_subplot(2,3,4)
plt.axes(ax4)
#hp.cartview(map_input*(mask/mask),title='Cosmological Signal',unit=r'$T$ [mK]',hold=True,rot=(162.5,3,0),lonra=lonra,latra=latra)
ax5 = fig.add_subplot(2,3,5)
plt.axes(ax5)
output0 = hp.cartview(residuals_spin0mapmake*(mask/mask),title='GMCA residuals, Num Components = 5',unit=r'$T$ [mK]',hold=True,rot=(162.5,3,0),lonra=lonra,latra=latra,return_projected_map=True)
ax6 = fig.add_subplot(2,3,6)
plt.axes(ax6)
output0_ns = hp.cartview(residuals_spin0nosyst*(mask/mask),title='GMCA residuals, Num Components = 5 No Systematic',unit=r'$T$ [mK]',hold=True,rot=(162.5,3,0),lonra=lonra,latra=latra,return_projected_map=True)



###Paper map plots
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 24})
plt.rcParams["axes.labelsize"] = 24

left=lonra[1]
right=lonra[0]

fig, (ax1) = plt.subplots(ncols=1,nrows=1)
plt.axes(ax1)
im1 = ax1.imshow(input1, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none',vmin=-0.1,vmax=0.28)
plt.xlabel('R.A.')
plt.ylabel('Dec.')
plt.xticks(np.arange(175,135,-5))
#plt.yticks(np.arange(-2,10,2))
clb=fig.colorbar(im1, orientation="horizontal", pad=0.2)
clb.set_label('mK')

plt.rcParams.update({'font.size': 24})
plt.rcParams["axes.labelsize"] = 24
fig, (ax1,ax2) = plt.subplots(ncols=2,nrows=1)
plt.axes(ax1)
im2 = ax1.imshow(output013, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none',vmin=-0.1,vmax=0.28)
plt.xlabel('R.A.')
plt.ylabel('Dec.')
clb=fig.colorbar(im2, orientation="horizontal", pad=0.2)
clb.set_label('mK')
plt.xticks(np.arange(175,135,-5))
plt.axes(ax2)
im5 = ax2.imshow(output0, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none')
plt.xlabel('R.A.')
plt.ylabel('Dec.')
clb=fig.colorbar(im5, orientation="horizontal", pad=0.2)
clb.set_label('mK')
plt.xticks(np.arange(175,135,-5))

fig, (ax1,ax2) = plt.subplots(ncols=2,nrows=1)
plt.axes(ax1)
im3 = ax1.imshow(output013_ns, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none',vmin=-0.1,vmax=0.28)
plt.xlabel('R.A.')
plt.ylabel('Dec.')
clb=fig.colorbar(im3, orientation="horizontal", pad=0.2)
clb.set_label('mK')

plt.axes(ax2)
im6 = ax2.imshow(output0_ns, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none',vmin=-0.1,vmax=0.28)
plt.xlabel('R.A.')
plt.ylabel('Dec.')
clb=fig.colorbar(im6, orientation="horizontal", pad=0.2)
clb.set_label('mK')

del residuals_spin0and1and3mapmake, residuals_spin0mapmake, residuals_spin0nosyst, residuals_spin0and1and3nosyst








#C_l space

residuals_spin0and1and3mapmake = mask*(obs_maps_spin0and1and3mapmake[ich]-X_gmca_spin0and1and3mapmake[ich])
residuals_spin0mapmake = mask*(obs_maps_spin0mapmake[ich]-X_gmca_spin0mapmake[ich])

residuals_spin0and1and3nosyst = mask*(obs_maps_spin0and1and3nosyst[ich]-X_gmca_spin0and1and3nosyst[ich])
residuals_spin0nosyst = mask*(obs_maps_spin0nosyst[ich]-X_gmca_spin0nosyst[ich])


plt.figure()
#plt.rcParams["figure.figsize"] = (8,6)


plt.rcParams["axes.labelsize"] = 16
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 16})
ell, y = g4i.plot_cl(map_input)
plt.semilogy(ell, y,'-',c='#1f77b4',label='Cosmological Signal');


ellnosyst, ynosyst = g4i.plot_cl(residuals_spin0and1and3nosyst)
plt.semilogy(ellnosyst,ynosyst,'o',c='r',label='Spin-0, 1 and 3 map-making, No Systematic')

ell_spin0nosyst, y_spin0nosyst = g4i.plot_cl(residuals_spin0nosyst)
plt.semilogy(ell_spin0nosyst,y_spin0nosyst,'.',c='m',label='Spin-0 map-making, No Systematic')#+np.str(n_s))


ell, y = g4i.plot_cl(residuals_spin0and1and3mapmake)
plt.semilogy(ell,y,'+',c='b',label='Spin-0, 1 and 3 map-making, Beam Squint Included')

ell_spin0mapmake, y_spin0mapmake = g4i.plot_cl(residuals_spin0mapmake)
plt.semilogy(ell_spin0mapmake,y_spin0mapmake,'x',c='g',label='Spin-0 map-making, Beam Squint Included')#+np.str(n_s))


#plt.title('Channel = '+np.str(np.round(freqs[ich],2))+'MHz')
plt.legend()
ax = plt.gca()
ax.set(xlim=[15,200],xlabel="$\\ell$",ylabel="$\\ell(\\ell+1)C_{\\ell}/(2\\pi)$ [mK$^2$]");










#K Space
map_input = 1.*h1cosm
mask[np.isnan(mask)]=0.
map_input*=(mask/mask)
map_input[map_input==0] = np.nan
map_input = (map_input.T - np.nanmean(map_input,axis=1)).T
map_input[np.isnan(map_input)]=0.
map_input*=mask


NSIDE=128
pix = np.arange(0,hp.nside2npix(NSIDE))
indexes_los = pix[mask!=0]


residuals_syst013 = obs_maps_spin0and1and3mapmake-X_gmca_spin0and1and3mapmake
residuals_syst0 = obs_maps_spin0mapmake-X_gmca_spin0mapmake
residuals_nosyst013 = obs_maps_spin0and1and3nosyst-X_gmca_spin0and1and3nosyst
residuals_nosyst0 = obs_maps_spin0nosyst-X_gmca_spin0nosyst

residuals_syst013*=mask
residuals_syst0*=mask
residuals_nosyst013*=mask
residuals_nosyst0*=mask

plt.figure()
plt.rcParams["figure.figsize"] = (8,6)
plt.rcParams["axes.labelsize"] = 16
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 16})

k, P = g4i.plot_nuPk(map_input,indexes_los,freqs)
plt.semilogy(k, P,'-',c='#1f77b4',label='Cosmological Signal')

k, P = g4i.plot_nuPk(residuals_nosyst013,indexes_los,freqs)
plt.semilogy(k, P,'o',c='r',label='Spin-0, 1 and 3 map-making, No Systematic')
         
k, P = g4i.plot_nuPk(residuals_nosyst0,indexes_los,freqs)
plt.semilogy(k, P,'.',c='m',label='Spin-0 map-making, No Systematic')

k, P = g4i.plot_nuPk(residuals_syst013,indexes_los,freqs)
plt.semilogy(k, P,'+',c='b',label='Spin-0, 1 and 3 map-making, Beam Squint Included')

k, P = g4i.plot_nuPk(residuals_syst0,indexes_los,freqs)
plt.semilogy(k, P,'x',c='g',label='Spin-0 map-making, Beam Squint Included')

plt.legend(loc='lower right')
ax = plt.gca()
ax.set(xlabel="$k_{\\nu}$ [MHz$^{-1}$]",ylabel="$P$ [mK$^2$ MHz$^2$]");



#Plot Mask
plt.figure()
hp.cartview(mask*(mask/mask),title='Apodized Mask',hold=True,rot=(162.5,3,0),lonra=lonra,latra=latra)















sys.exit()
####Plots in other orientation####
mask = (obs_maps_spin0and1and3mapmake/obs_maps_spin0and1and3mapmake)*(obs_maps_spin0mapmake/obs_maps_spin0mapmake)*(obs_maps_spin0nosyst/obs_maps_spin0nosyst)*(obs_maps_spin0and1and3nosyst/obs_maps_spin0and1and3nosyst)
mask[np.isnan(mask)]=0.
mask = mask[0]
# pick a channel number
ich = 140

#Map Space
map_input = 1.*h1cosm[ich]
mask[np.isnan(mask)]=0.
map_input*=mask
map_input[np.isnan(map_input)]=0.
map_input[map_input!=0] -= np.mean(map_input[map_input!=0])



Apodize=False
if Apodize:
    '''
    hp.write_map('Mask.fits',mask[0])
    
    mask_file = "mask_file=Mask.fits\n"
    hole_min_size = "hole_min_size=0\n"
    distance_file= "distance_file=Distance.fits"
    
    File_Object = open(r"./apodizefile.txt","w+")
    File_Object.writelines([mask_file,hole_min_size,distance_file])
    File_Object.close()
        
    command = '/mirror/scratch/mccallum/Healpix_3.60/bin/process_mask ' + './apodizefile.txt'
    exit_process = subprocess.call(command, shell=True)
    '''
    dist = hp.read_map('Distance.fits')
    
    aposcale=2.*np.mean(dist[dist>0])
    mask_apo = 1.*dist
    mask_apo[dist>aposcale]=1.
    mask_apo[dist<=aposcale]/=aposcale

    mask*=mask_apo
    mask[np.isnan(mask)]=0.


map_input*=mask


residuals_spin0and1and3mapmake = obs_maps_spin0and1and3mapmake[ich]-X_gmca_spin0and1and3mapmake[ich]
residuals_spin0mapmake = obs_maps_spin0mapmake[ich]-X_gmca_spin0mapmake[ich]

residuals_spin0and1and3nosyst = obs_maps_spin0and1and3nosyst[ich]-X_gmca_spin0and1and3nosyst[ich]
residuals_spin0nosyst = obs_maps_spin0nosyst[ich]-X_gmca_spin0nosyst[ich]


residuals_spin0and1and3mapmake*=mask
residuals_spin0mapmake*=mask
residuals_spin0and1and3nosyst*=mask
residuals_spin0nosyst*=mask


lonra=[-22.5,17.5]
latra=[-6,6]


fig, (ax1,ax2,ax3) = plt.subplots(ncols=3)
plt.axes(ax1)
input1 = hp.cartview(map_input*(mask/mask),rot=(162.5,3,0),return_projected_map=True,lonra=lonra,latra=latra,hold=True,title='Input HI')
plt.axes(ax2)
output0 = hp.cartview(residuals_spin0mapmake*(mask/mask),rot=(162.5,3,0),return_projected_map=True,lonra=lonra,latra=latra,hold=True,title='GMCA residuals 0')
plt.axes(ax3)
output02 = hp.cartview(residuals_spin0and1and3mapmake*(mask/mask),rot=(162.5,3,0),return_projected_map=True,lonra=lonra,latra=latra,hold=True,title='GMCA residuals 0 and 1')

fig, (ax1,ax2,ax3) = plt.subplots(ncols=3)
plt.axes(ax1)
input1 = hp.cartview(map_input*(mask/mask),rot=(162.5,3,0),return_projected_map=True,lonra=lonra,latra=latra,hold=True,title='Input HI')
plt.axes(ax2)
output0_ns = hp.cartview(residuals_spin0nosyst*(mask/mask),rot=(162.5,3,0),return_projected_map=True,lonra=lonra,latra=latra,hold=True,title='GMCA residuals 0 No Syst')
plt.axes(ax3)
output02_ns = hp.cartview(residuals_spin0and1and3nosyst*(mask/mask),rot=(162.5,3,0),return_projected_map=True,lonra=lonra,latra=latra,hold=True,title='GMCA residuals 0 and 1 No Syst')




input1flip = input1#np.fliplr(input1)
output0flip = output0#np.fliplr(output0)
output02flip = output02#np.fliplr(output02)
output0_nsflip = output0_ns#np.fliplr(output0_ns)
output02_nsflip = output02_ns#np.fliplr(output02_ns)




left=lonra[1]
right=lonra[0]

fig = plt.figure(figsize=(20, 6))
fig.suptitle('Channel = '+np.str(np.round(freqs[ich],2))+'MHz Spin 0 (and 2) Map-Making Lower (Upper)',fontsize=20)
ax1 = fig.add_subplot(2,3,1)
plt.axes(ax1)
plt.title('Input Cosmological Signal')
im1 = ax1.imshow(input1flip, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none')
plt.xlabel('R.A.')
plt.ylabel('Dec.')
clb=fig.colorbar(im1, orientation="horizontal", pad=0.2)
clb.set_label('mK')
ax2 = fig.add_subplot(2,3,2)
plt.axes(ax2)
plt.title('GMCA residuals, Num Components = 5, Spin-0 and 1')
im2 = ax2.imshow(output02flip, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none')
plt.xlabel('R.A.')
plt.ylabel('Dec.')
clb=fig.colorbar(im2, orientation="horizontal", pad=0.2)
clb.set_label('mK')
ax3 = fig.add_subplot(2,3,3)
plt.axes(ax3)
plt.title('GMCA residuals, Num Components = 5, Spin-0 and 1, No Systematic')
im3 = ax3.imshow(output02_nsflip, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none')
plt.xlabel('R.A.')
plt.ylabel('Dec.')
clb=fig.colorbar(im3, orientation="horizontal", pad=0.2)
clb.set_label('mK')



ax5 = fig.add_subplot(2,3,5)
plt.axes(ax5)
plt.title('GMCA residuals, Num Components = 5, Spin-0')
im5 = ax5.imshow(output0flip, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none')
plt.xlabel('R.A.')
plt.ylabel('Dec.')
clb=fig.colorbar(im5, orientation="horizontal", pad=0.2)
clb.set_label('mK')
ax6 = fig.add_subplot(2,3,6)
plt.axes(ax6)
plt.title('GMCA residuals, Num Components = 5, Spin-0, No Systematic')
im6 = ax6.imshow(output0_nsflip, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none')
plt.xlabel('R.A.')
plt.ylabel('Dec.')
clb=fig.colorbar(im6, orientation="horizontal", pad=0.2)
clb.set_label('mK')







###Paper map plots
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 24})
plt.rcParams["axes.labelsize"] = 24

fig, (ax1) = plt.subplots(ncols=1,nrows=1)
#fig.suptitle('Channel = '+np.str(np.round(freqs[ich],2))+'MHz Spin 0 (and 2) Map-Making Lower (Upper)',fontsize=20)
plt.axes(ax1)
#plt.title('Input Cosmological Signal')
im1 = ax1.imshow(input1flip, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none',vmin=-0.1,vmax=0.28)
plt.xlabel('R.A.')
plt.ylabel('Dec.')
plt.xticks(np.arange(175,135,-5))
#plt.yticks(np.arange(-2,10,2))
clb=fig.colorbar(im1, orientation="horizontal", pad=0.2)
clb.set_label('mK')



plt.rcParams.update({'font.size': 24})
plt.rcParams["axes.labelsize"] = 24
fig, (ax1,ax2) = plt.subplots(ncols=2,nrows=1)
#fig.suptitle('Channel = '+np.str(np.round(freqs[ich],2))+'MHz Spin 0 (and 1) Map-Making Lower (Upper)',fontsize=20)
plt.axes(ax1)
#plt.title('GMCA residuals, Num Components = 5, Spin-0 and 1')
im2 = ax1.imshow(output02flip, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none',vmin=-0.1,vmax=0.28)
plt.xlabel('R.A.')
plt.ylabel('Dec.')
clb=fig.colorbar(im2, orientation="horizontal", pad=0.2)
clb.set_label('mK')
plt.xticks(np.arange(175,135,-5))

plt.axes(ax2)
#plt.title('GMCA residuals, Num Components = 5, Spin-0')
im5 = ax2.imshow(output0flip, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none')
plt.xlabel('R.A.')
plt.ylabel('Dec.')
clb=fig.colorbar(im5, orientation="horizontal", pad=0.2)
clb.set_label('mK')
plt.xticks(np.arange(175,135,-5))


fig, (ax1,ax2) = plt.subplots(ncols=2,nrows=1)
plt.axes(ax1)
#plt.title('GMCA residuals, Num Components = 5, Spin-0 and 1, No Systematic')
im3 = ax1.imshow(output02_nsflip, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none',vmin=-0.1,vmax=0.28)
plt.xlabel('R.A.')
plt.ylabel('Dec.')
clb=fig.colorbar(im3, orientation="horizontal", pad=0.2)
clb.set_label('mK')

plt.axes(ax2)
#plt.title('GMCA residuals, Num Components = 5, Spin-0, No Systematic')
im6 = ax2.imshow(output0_nsflip, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none',vmin=-0.1,vmax=0.28)
plt.xlabel('R.A.')
plt.ylabel('Dec.')
clb=fig.colorbar(im6, orientation="horizontal", pad=0.2)
clb.set_label('mK')
