#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 26 18:23:03 2021

@author: mccallum
"""
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt



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
    gsyncI[i-start],gsyncQ[i-start],gsyncU[i-start] = hp.read_map('../IM_CRIME_Data/gsync/gsync_'+val+'_b.fits',field=(0,1,2))
    psource[i-start] = hp.read_map('../IM_CRIME_Data/psources/psources_'+val+'.fits')
    egfreefree[i-start] = hp.read_map('../IM_CRIME_Data/egfree/egfree_'+val+'.fits')
    gfreefree[i-start] = hp.read_map('../IM_CRIME_Data/gfree/gfree_'+val+'.fits')
    h1cosm[i-start] = hp.read_map('../IM_CRIME_Data/cosmo/sim_2048_'+val+'.fits')


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

    
n_freq, nu_0, nu_f, z_0, z_f = np.loadtxt('../IM_CRIME_Data/nuTable.txt').T



#Plot the fields as a function of frequency at RA=155degrees and DEC=3degrees
RA=155
DEC=3
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 24})

plt.semilogy(((nu_f+nu_0)/2)[start-1:end], h1cosm.T[hp.ang2pix(NSIDE,np.deg2rad(DEC),np.deg2rad(RA),lonlat=True)],label='HI Cosmological Signal',linewidth=1)
plt.semilogy(((nu_f+nu_0)/2)[start-1:end], gsyncI.T[hp.ang2pix(NSIDE,np.deg2rad(DEC),np.deg2rad(RA),lonlat=True)],label='Galactic Synchrotron I',linewidth=1)
plt.semilogy(((nu_f+nu_0)/2)[start-1:end], np.abs(gsyncQ.T[hp.ang2pix(NSIDE,np.deg2rad(DEC),np.deg2rad(RA),lonlat=True)]),label='Galactic Synchrotron Q',linewidth=1)
plt.semilogy(((nu_f+nu_0)/2)[start-1:end], np.abs(gsyncU.T[hp.ang2pix(NSIDE,np.deg2rad(DEC),np.deg2rad(RA),lonlat=True)]),label='Galactic Synchrotron U',linewidth=1)

plt.semilogy(((nu_f+nu_0)/2)[start-1:end], np.abs(psource.T[hp.ang2pix(NSIDE,np.deg2rad(DEC),np.deg2rad(RA),lonlat=True)]),'c-',label='Point Sources',linewidth=1)
plt.semilogy(((nu_f+nu_0)/2)[start-1:end], np.abs(egfreefree.T[hp.ang2pix(NSIDE,np.deg2rad(DEC),np.deg2rad(RA),lonlat=True)]),'k-',label='Extragalactic Free-Free',linewidth=1)
plt.semilogy(((nu_f+nu_0)/2)[start-1:end], np.abs(gfreefree.T[hp.ang2pix(NSIDE,np.deg2rad(DEC),np.deg2rad(RA),lonlat=True)]),'m-',label='Galactic Free-Free',linewidth=1)


plt.semilogy(((nu_f+nu_0)/2)[start-1:end], -(psource.T[hp.ang2pix(NSIDE,np.deg2rad(DEC),np.deg2rad(RA),lonlat=True)]),'c.',markersize=6)

plt.semilogy(((nu_f+nu_0)/2)[start-1:end], -(gsyncQ.T[hp.ang2pix(NSIDE,np.deg2rad(DEC),np.deg2rad(RA),lonlat=True)]),'g.',markersize=6)
plt.semilogy(((nu_f+nu_0)/2)[start-1:end], -(gsyncU.T[hp.ang2pix(NSIDE,np.deg2rad(DEC),np.deg2rad(RA),lonlat=True)]),'r.',markersize=6)


plt.semilogy(((nu_f+nu_0)/2)[start-1:end], -(egfreefree.T[hp.ang2pix(NSIDE,np.deg2rad(DEC),np.deg2rad(RA),lonlat=True)]),'k.',markersize=6)
plt.semilogy(((nu_f+nu_0)/2)[start-1:end], -(gfreefree.T[hp.ang2pix(NSIDE,np.deg2rad(DEC),np.deg2rad(RA),lonlat=True)]),'m.',markersize=6)


plt.legend(loc='best',ncol=2)
plt.ylabel('mK')
plt.xlabel('Frequency (MHz)')





plt.figure()
#Plot input maps at frequency of ~915.5MHz
freq_band=142
mask = hp.read_map('./Survey_Mask_NSIDE'+np.str(NSIDE)+'.fits')
lonra=[-22.5,17.5]
latra=[-6,6]
hI_in = hp.cartview(h1cosm[freq_band]/mask,rot=(162.5,3,0),return_projected_map=True,lonra=lonra,latra=latra,hold=True)
gsync_in_I = hp.cartview(gsyncI[freq_band]/mask,rot=(162.5,3,0),return_projected_map=True,lonra=lonra,latra=latra,hold=True)
gsync_in_Q = hp.cartview(gsyncQ[freq_band]/mask,rot=(162.5,3,0),return_projected_map=True,lonra=lonra,latra=latra,hold=True)
gsync_in_U = hp.cartview(gsyncU[freq_band]/mask,rot=(162.5,3,0),return_projected_map=True,lonra=lonra,latra=latra,hold=True)
psource_in = hp.cartview(psource[freq_band]/mask,rot=(162.5,3,0),return_projected_map=True,lonra=lonra,latra=latra,hold=True)
egfreefree_in = hp.cartview(egfreefree[freq_band]/mask,rot=(162.5,3,0),return_projected_map=True,lonra=lonra,latra=latra,hold=True)
gfreefree_in = hp.cartview(gfreefree[freq_band]/mask,rot=(162.5,3,0),return_projected_map=True,lonra=lonra,latra=latra,hold=True)
total_I = hI_in + gsync_in_I + psource_in + egfreefree_in + gfreefree_in
plt.close()


left=lonra[1]
right=lonra[0]



plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 15})

fig,((ax1,ax2),(ax3,ax4),(ax5,ax6),(ax7,ax8)) = plt.subplots(ncols=2,nrows=4,figsize=(4, 3),sharex=True)

sh=0.4

plt.axes(ax1)
im1 = ax1.imshow(total_I, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none')
plt.ylabel('Dec.')
clb=fig.colorbar(im1,orientation="vertical", shrink=sh)
clb.set_label('mK')
plt.axes(ax2)
im2 = ax2.imshow(hI_in, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none')
plt.ylabel('Dec.')
clb = fig.colorbar(im2, orientation="vertical",shrink=sh)
clb.set_label('mK')
plt.axes(ax3)
im3 = ax3.imshow(psource_in, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none')
plt.ylabel('Dec.')
clb = fig.colorbar(im3, orientation="vertical",shrink=sh)
clb.set_label('mK')
plt.axes(ax4)
im4 = ax4.imshow(gfreefree_in, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none')
plt.ylabel('Dec.')
clb = fig.colorbar(im4, orientation="vertical",shrink=sh)
clb.set_label('mK')
plt.axes(ax5)
im5 = ax5.imshow(egfreefree_in, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none')
plt.ylabel('Dec.')
clb = fig.colorbar(im5, orientation="vertical",shrink=sh)
clb.set_label('mK')
plt.axes(ax6)
im6 = ax6.imshow(gsync_in_I, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none')
plt.ylabel('Dec.')
clb = fig.colorbar(im6, orientation="vertical",shrink=sh)
clb.set_label('mK')
plt.axes(ax7)
im7 = ax7.imshow(gsync_in_Q, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none')
plt.xlabel('R.A.')
plt.ylabel('Dec.')
clb = fig.colorbar(im7, orientation="vertical",shrink=sh)
clb.set_label('mK')
plt.axes(ax8)
im8 = ax8.imshow(gsync_in_U, origin='lower',extent=(160+left,160+right,3+latra[0],3+latra[1]), interpolation = 'none')
plt.xlabel('R.A.')
plt.ylabel('Dec.')
clb = fig.colorbar(im8, orientation="vertical",shrink=sh)
clb.set_label('mK')


#plt.subplots_adjust(wspace=0, hspace=0)
plt.subplots_adjust(hspace = 0)
plt.subplots_adjust(hspace = -0.5)