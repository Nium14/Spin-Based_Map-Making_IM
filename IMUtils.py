#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 14:38:02 2021

@author: mccallum
"""
import healpy as hp



#####BEAM SQUINT#####
def get_derivativeoffield(A,lmax):
    '''
    Get the derivative of an input field A.

    Parameters
    ----------
    A : Array - healpix map
        Input field.
    lmax : integer
        The maximum l to go to when working out the derivative.

    Returns
    -------
    dA : Array - healpix maps
        Outputs the derivative of the input field. As d/dtheta + i*d/dphi

    '''
    NSIDE=hp.get_nside(A)
    
    almA = hp.map2alm(A,lmax=lmax)

    A,dAdth,dAdph = hp.alm2map_der1(almA, NSIDE)
    
    dA = (dAdth + 1j*dAdph)

    return dA



#####SCANNING STRATEGIES #####
import ephem
import numpy as np
from functools import partial
from ephem import newton

def find_time_of_elevation(time,RA_in,DEC_in,elev,lat_in, long_in, elevation_in):
    DEC_in_deg=ephem.degrees(np.deg2rad(DEC_in))
    RA_in_deg=ephem.degrees(np.deg2rad(RA_in))
    Target=ephem.FixedBody()
    pb = ephem.Observer()
    Target._ra=RA_in_deg
    Target._dec=DEC_in_deg
    pb.date=time
    pb.long = long_in
    pb.lat = lat_in
    pb.elevation = elevation_in
    out=Target.compute(pb)
    root=np.rad2deg(Target.alt)-np.rad2deg(elev)
    return root


def dec_endaz(Az,elev,dec_end,pb):
    ra,dec=np.rad2deg(pb.radec_of(Az,elev))
    return dec-np.rad2deg(dec_end)


def find_field_timeandaz(LowerRA,UpperRA,LowerDec,UpperDec,pb,elev,rising,lat_in,long_in,elevation_in):
    if rising:
        field = ephem.FixedBody()
        field._ra = ephem.degrees(np.str(LowerRA))
        field._dec = ephem.degrees(np.str(LowerDec))
        #if -1 setting, if +1 rising
        risingorsetting=1
        dec_end=np.deg2rad(UpperDec)
    else:
        #Setting
        field = ephem.FixedBody()
        field._ra = ephem.degrees(np.str(LowerRA))
        field._dec = ephem.degrees(np.str(UpperDec))
        #field._ra = ephem.degrees('40')
        #field._dec = ephem.degrees('-30')
    
        #if -1 setting, if +1 rising
        risingorsetting=-1
        dec_end=np.deg2rad(LowerDec)
    
    #Find Date/Time Field Visible
    guess2=pb.next_transit(field)
    guess1=ephem.Date(guess2-4*risingorsetting*ephem.hour)
    elevtime_partial = partial(find_time_of_elevation,RA_in=np.rad2deg(field.ra),DEC_in=np.rad2deg(field.dec),elev=np.deg2rad(elev),lat_in=lat_in,long_in=long_in,elevation_in=elevation_in)
    out=newton(elevtime_partial,guess1,guess2)
    pb.date=out
    field.compute(pb)
    #Find Azimuth start and end based on field decs
    azstart=np.rad2deg(field.az)
    guess1=field.az
    guess2=field.az+np.deg2rad(40)
    dec_endaz_partial = partial(dec_endaz,elev=np.deg2rad(elev),dec_end=dec_end,pb=pb)
    out=newton(dec_endaz_partial,guess1,guess2)
    azend=np.rad2deg(out)
    
    while azstart>180:
        azstart-=360
    while azend>180:
        azend-=360
    
    print(azstart)
    print(azend)
    
    return pb,azstart,azend


def crossingangle(RA_in,DEC_in,az,delta_az,elev,pb):
    new_az=az+delta_az
    ra_new, dec_new = pb.radec_of(np.deg2rad(new_az),np.deg2rad(elev))
    delta_dec=float(repr(dec_new))-DEC_in
    delta_ra= float(repr(ra_new))-RA_in
    if delta_ra > np.pi:
        delta_ra-=2.*np.pi
    elif delta_ra<-np.pi:
        delta_ra+=2.*np.pi
    psi = np.arctan2(delta_dec,delta_ra*np.cos(DEC_in))
    #if delta_az<0:
    #    psi*=-1
    return psi


def generate_scan(elevation,azend,azstart,pb,surveytime,fsamp,az_freq,delta_az):
    #Generate lists to add to for scan data
    RA=[]
    DEC=[]
    AZ=[]
    PSI=[]
    
    time_array = np.arange(0,surveytime,1./fsamp)
    
    azaim = 1*azend
    az = 1*azstart
    for time in time_array:
        if azend>azstart:
            if azaim == azend:
                if az+(1./fsamp)*az_freq < azend:
                    azaim = azend
                else:
                    azaim = azstart
                print(az)
                az += (1./fsamp)*az_freq
                print(az)
                #print(azaim)
                #print(azend)
                
            if azaim == azstart:
                if az-(1./fsamp)*az_freq > azstart:
                    azaim = azstart
                else:
                    azaim = azend
                print(az)
                az -= (1./fsamp)*az_freq
                print(az)
                #print(azaim)
                #print(azstart)
                
        else:
            if azaim == azend:
                if az-(1./fsamp)*az_freq > azend:
                    azaim = azend
                else:
                    azaim = azstart
                print(az)
                az -= (1./fsamp)*az_freq
                print(az)
                #print(azaim)
                #print(azend)
                
            if azaim == azstart:
                if az-(1./fsamp)*az_freq < azstart:
                    azaim = azstart
                else:
                    azaim = azend
                print(az)
                az += (1./fsamp)*az_freq
                print(az)
                #print(azaim)
                #print(azstart)
        
                
        ra,dec=pb.radec_of(np.deg2rad(az),np.deg2rad(elevation))
        RA.append(ra)
        DEC.append(dec)
        AZ.append(az)
        #print(ra,dec)
        #print(az)
        
        pb.date += (1./fsamp)*ephem.second
        #print((1./fsamp)*ephem.second)
        #print(pb.date.datetime())
    
        #delta_az always positive due to focal plane orientation to sky - makes calculation easier
        psi=crossingangle(ra,dec,az,delta_az,elevation,pb)
        PSI.append(psi)
            

    #Change to numpy arrays
    RA = np.array(RA)
    DEC = np.array(DEC)
    AZ = np.array(AZ)
    PSI=np.array(PSI)
    
    #Convert from pyephem to healpy dec convention
    DEC = np.pi/2. - DEC
    return RA,DEC,AZ,PSI





#####MAP-MAKING#####

#Define a function to map make simple bin style for spin-k signals
def mapmake_simplebin_spink(maps,pixcond,tod,pixel_index,psi,NSIDE,spins=np.array([0]),mask=None):
    '''
    Simple Binned Map-making function for arbitrary spin fields.

    Parameters
    ----------
    maps : Numpy array of Healpix maps
        The signals will be output here -> spin-0 needs one map, other spins need 2 maps each
    pixcond : Numpy array - Healpix Map
        Single map to which the condition of the map-making matrix for each pixel is returned
    tod : numpy array
        Input time ordered data -> single detector, summed, pair etc. (Needs sufficient psi coverage such that matrix is not singular).
    pixel_index : Numpy array - Healpix Map
        Map of the pixel indices.
    psi : numpy array
        Crossing angles of the time ordered data.
    NSIDE : Integer
        NSIDE of the maps.
    spins : Numpy array, optional
        The spin field to solve for e.g. spins=np.array([0,1,3]) would solve for fields of spin 0, 1, and 3. The default is np.array([0]).
    mask : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    maps : Numpy array of Healpix maps
        Array of maps in ascending spin order (spin-0 just one row) (spin k!=0 has Z_k^Q row followed by Z_k^U row).
    pixcond : Numpy array - Healpix Map
        Array of condition of map-making matrix in each pixel.
    '''
    
    #Calculate hitmap
    hitmap=np.zeros(hp.nside2npix(NSIDE))
    np.add.at(hitmap,pixel_index.astype(int),1)
    if mask!=None:
        hitmap*=mask
        
    #if spin-0 included as this only needs one row/column in the matrix
    if len(np.argwhere(spins==0))>0:
        cos_sin_matrix = np.zeros((hp.nside2npix(NSIDE),len(spins)*2-1,len(spins)*2-1))
        data_vector = np.zeros((len(spins)*2 -1,hp.nside2npix(NSIDE)))
        spin_temp = np.sort(np.concatenate((spins[spins>0],spins[spins>0])))

        #Data in first row of vector is just averaged tod - no psi dependence
        np.add.at(data_vector[0], pixel_index.astype(int), tod)
        data_vector[0] /= hitmap
        #If other spins included add cos and sin rows in so
        if len(np.argwhere(spins!=0))>0:
            for i in np.arange(1,len(spin_temp)+1,2):
                for j in np.arange(1,len(spin_temp)+1,2):
                    np.add.at(cos_sin_matrix.T[i,j,:], pixel_index.astype(int), np.cos(spin_temp[j]*psi)*np.cos(spin_temp[i]*psi))
                    np.add.at(cos_sin_matrix.T[i+1,j,:], pixel_index.astype(int), np.cos(spin_temp[j]*psi)*np.sin(spin_temp[i]*psi))
                    
                    np.add.at(cos_sin_matrix.T[i+1,j+1,:], pixel_index.astype(int), np.sin(spin_temp[j]*psi)*np.sin(spin_temp[i]*psi))
                    np.add.at(cos_sin_matrix.T[i,j+1,:], pixel_index.astype(int), np.sin(spin_temp[j]*psi)*np.cos(spin_temp[i]*psi))
                
                np.add.at(cos_sin_matrix.T[i,0,:], pixel_index.astype(int), np.cos(spin_temp[i]*psi))
                np.add.at(cos_sin_matrix.T[0,i,:], pixel_index.astype(int), np.cos(spin_temp[i]*psi))
                np.add.at(cos_sin_matrix.T[i+1,0,:], pixel_index.astype(int), np.sin(spin_temp[i]*psi))
                np.add.at(cos_sin_matrix.T[0,i+1,:], pixel_index.astype(int), np.sin(spin_temp[i]*psi))
                
                np.add.at(data_vector[i], pixel_index.astype(int), tod*np.cos(spin_temp[i]*psi))
                np.add.at(data_vector[i+1], pixel_index.astype(int), tod*np.sin(spin_temp[i]*psi))
        
            cos_sin_matrix=cos_sin_matrix.T
            cos_sin_matrix/=hitmap
            cos_sin_matrix=cos_sin_matrix.T
            
            data_vector[1:]/=hitmap
        cos_sin_matrix.T[0,0,:]=1.
    #If not zero then the spins all have cos and sin rows in matrix
    else:
        cos_sin_matrix = np.zeros((hp.nside2npix(NSIDE),len(spins)*2,len(spins)*2))
        data_vector = np.zeros((len(spins)*2,hp.nside2npix(NSIDE)))
        spin_temp = np.sort(np.concatenate((spins[spins>0],spins[spins>0])))

        for i in np.arange(0,len(spin_temp),2):
            for j in np.arange(0,len(spin_temp),2):
                np.add.at(cos_sin_matrix.T[i,j,:], pixel_index.astype(int), np.cos(spin_temp[j]*psi)*np.cos(spin_temp[i]*psi))
                np.add.at(cos_sin_matrix.T[i+1,j,:], pixel_index.astype(int), np.cos(spin_temp[j]*psi)*np.sin(spin_temp[i]*psi))
                
                np.add.at(cos_sin_matrix.T[i+1,j+1,:], pixel_index.astype(int), np.sin(spin_temp[j]*psi)*np.sin(spin_temp[i]*psi))
                np.add.at(cos_sin_matrix.T[i,j+1,:], pixel_index.astype(int), np.sin(spin_temp[j]*psi)*np.cos(spin_temp[i]*psi))                
            
            np.add.at(data_vector[i], pixel_index.astype(int), tod*np.cos(spin_temp[i]*psi))
            np.add.at(data_vector[i+1], pixel_index.astype(int), tod*np.sin(spin_temp[i]*psi))
        cos_sin_matrix=cos_sin_matrix.T
        cos_sin_matrix/=hitmap
        cos_sin_matrix=cos_sin_matrix.T
        data_vector/=hitmap
    
    #Loop over the pixels and map-make
    #Probably possible to speed up somehow by removing loop here??
    pixtemp=np.arange(0,hp.nside2npix(NSIDE))[hitmap==0]
    maps[:,pixtemp] = np.nan
    pixtemp=np.arange(0,hp.nside2npix(NSIDE))[hitmap!=0]
    for i in pixtemp:
        A = np.matrix(cos_sin_matrix[i])
        D = np.matrix([[j] for j in data_vector.T[i].tolist()])
        #Get the conditioning of the matrix
        pixcond[i] = np.linalg.cond(A)
        
        try:
            #Perform the map-making
            maps[:,i] = np.squeeze(np.asarray(A.I * D))
        except np.linalg.LinAlgError as err:
            if 'Singular matrix' in str(err):
                maps[:,i] = np.nan
                print ('Singular Matrix at pixel '+np.str(i)+' setting to nan.')
            else:
                raise
            
    
    return maps,pixcond




#Define a function to map make simple bin style for spin-k signals
def mapmake_h_kbin_spink(maps,pixcond,tod,pixel_index,psi,NSIDE,spins=np.array([0]),mask=None):
    """
    Simple Binned Map-making function for arbitrary spin fields using h_k approach
    
    Inputs:
    maps is an array of healpy maps to which the signals will be output -> spin-0 needs one map, other spins need 2 maps each
    pixcond is a single map to which the condition of the map-making matrix for each pixel is returned
    tod is the input time ordered data -> single detector, summed, pair etc. (Needs sufficient psi coverage such that matrix is not singular)
    pixel_index is a map of the pixel indices
    hitmap is the map of the hits in a pixel
    NSIDE is the nside of the maps
    spins is an array of the spin field to solve for e.g. spins=np.array([0,1,3]) would solve for fields of spin 0, 1, and 3
    
    returns
    -array of maps in ascending spin order (spin-0 just one row) (spin k!=0 has spin -k field (Z_k^Q-iZ_k^U) row followed by spin +k field (Z_k^Q+iZ_k^U) row)
    -array of condition of map-making matrix in each pixel
    """
    
    #Calculate hitmap
    hitmap=np.zeros(hp.nside2npix(NSIDE))
    np.add.at(hitmap,pixel_index.astype(int),1)
    if mask!=None:
        hitmap*=mask
        
    #if spin-0 included as this only needs one row/column in the matrix
    if len(np.argwhere(spins==0))>0:
        h_n_matrix = np.zeros((hp.nside2npix(NSIDE),len(spins)*2-1,len(spins)*2-1))*0J
        data_vector = np.zeros((len(spins)*2 -1,hp.nside2npix(NSIDE)))*0J
        spin_temp = np.sort(np.concatenate((spins,spins[spins>0])))
        spin_temp[::2]*=-1
    else:
        h_n_matrix = np.zeros((hp.nside2npix(NSIDE),len(spins)*2,len(spins)*2))*0J
        data_vector = np.zeros((len(spins)*2,hp.nside2npix(NSIDE)))*0J
        spin_temp = np.sort(np.concatenate((spins[spins>0],spins[spins>0])))
        spin_temp[1::2]*=-1

    for i in np.arange(0,len(spin_temp),1):
        for j in np.arange(0,len(spin_temp),1):
            np.add.at(h_n_matrix.T[i,j,:], pixel_index.astype(int), np.cos((spin_temp[i]+spin_temp[j])*psi)+1j*np.sin((spin_temp[i]+spin_temp[j])*psi))
        
        if spin_temp[i]==0:
            np.add.at(data_vector[i], pixel_index.astype(int), tod*(np.cos(spin_temp[i]*psi)+1j*np.sin(spin_temp[i]*psi)))
        else:
            np.add.at(data_vector[i], pixel_index.astype(int), 0.5*tod*(np.cos(spin_temp[i]*psi)+1j*np.sin(spin_temp[i]*psi)))

    h_n_matrix=h_n_matrix.T
    h_n_matrix/=hitmap
    h_n_matrix=h_n_matrix.T
        
    data_vector[0:]/=hitmap
    if len(np.argwhere(spins==0))>0:
        h_n_matrix.T[0,0,:]=1.
        h_n_matrix.T[1:,0,:]*=0.5
        h_n_matrix.T[0,1:,:]*=0.5
        if len(np.argwhere(spins!=0))>0:
            h_n_matrix.T[1:,1:,:]*=0.25
    else:
        h_n_matrix.T[0:,0:,:]*=0.25
    
    
    #Loop over the pixels and map-make
    #Probably possible to speed up somehow by removing loop here??
    pixtemp=np.arange(0,hp.nside2npix(NSIDE))[hitmap==0]
    maps[:,pixtemp] = np.nan
    pixtemp=np.arange(0,hp.nside2npix(NSIDE))[hitmap!=0]
    
    for i in pixtemp:
        A = np.matrix(h_n_matrix[i])
        D = np.matrix([[j] for j in data_vector.T[i].tolist()])
        #Get the conditioning of the matrix
        pixcond[i] = np.linalg.cond(A)
        
        try:
            #Perform the map-making
            maps[:,i] = np.squeeze(np.asarray(A.I * D))
        except np.linalg.LinAlgError as err:
            if 'Singular matrix' in str(err):
                maps[:,i] = np.nan
                print ('Singular Matrix at pixel '+np.str(i)+' setting to nan.')
            else:
                raise
    return maps,pixcond






#####FG Removal#####
import subprocess
def ApodizeMask(mask,apodisation_scale):
    hp.write_map('Mask.fits',mask)
    
    mask_file = "mask_file=Mask.fits\n"
    hole_min_size = "hole_min_size=0\n"
    distance_file= "distance_file=Distance.fits"
    
    File_Object = open(r"./apodizefile.txt","w+")
    File_Object.writelines([mask_file,hole_min_size,distance_file])
    File_Object.close()
        
    command = '/mirror/scratch/mccallum/Healpix_3.60/bin/process_mask ' + './apodizefile.txt'
    exit_process = subprocess.call(command, shell=True)
    dist = hp.read_map('Distance.fits')
    
    aposcale=apodisation_scale*np.mean(dist[dist>0])
    mask_apo = 1.*dist
    mask_apo[dist>aposcale]=1.
    mask_apo[dist<=aposcale]/=aposcale

    mask*=mask_apo
    mask[np.isnan(mask)]=0.
    
    return mask