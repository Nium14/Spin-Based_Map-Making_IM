#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 14:38:02 2021

@author: mccallum
"""
import healpy as hp


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