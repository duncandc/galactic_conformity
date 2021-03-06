#!/usr/bin/env python

#Author: Duncan Campbell
#Written: Feburary 20, 2014
#Yale University
#Description: construct a mock redshift survey

import numpy as np
import sys
import matplotlib.pyplot as plt
import math
from astropy import table
from astropy.io import ascii
from astropy import cosmology
from scipy.interpolate import interp1d
import custom_utilities as cu

def main():

    filepath  = cu.get_output_path()+'processed_data/hearin_mocks/custom_catalogues/'
    savepath = filepath+'tinker_radec_mocks/'
    filename  = sys.argv[1]
    catalogue = filename[0:-33]
    
    #file format: ID x y z Vx Vy Vz Mag with headers
    data = table.Table.read(filepath+filename, format='ascii')
    data = np.array(data)
    dtype = [('ID','>i8'),('ra','>f8'),('dec','>f8'),('z','>f8'),('app_mag','>f8'),('abs_mag','>f8')]
    dtype = np.dtype(dtype)
    mock  = np.recarray((len(data),), dtype=dtype)
    mock['ID'] = data['ID']
    #mock['k'] = np.arange(0,len(mock),1,dtype=int)

    app_mag_lim = 17.77
    Lbox = 250.0 #Mpc/h
    cosmo = cosmology.core.FlatLambdaCDM(H0=100.0,Om0=0.27) #h=1
    c = 299792.458 #speed of light in km/s

    #compute comoving distance from observer
    r = np.sqrt(data['x']**2+data['y']**2+data['z']**2)

    #compute radial velocity
    ct = data['z']/r
    st = np.sqrt(1.0-ct**2)
    cp = data['x']/np.sqrt(data['x']**2+data['y']**2)
    sp = data['y']/np.sqrt(data['x']**2+data['y']**2)
    vr = data['Vx']*st*cp + data['Vy']*st*sp + data['Vz']*ct

    #compute cosmological redshift and add contribution from perculiar velocity
    y = np.arange(0,2.0,0.01)
    x = cosmology.funcs.comoving_distance(y, cosmo=cosmo)
    f = interp1d(x, y, kind='cubic')
    z_cos = f(r)
    mock['z'] = z_cos+(vr/c)*(1.0+z_cos)

    #calculate spherical coordinates
    theta = np.arccos(data['z']/r)
    phi   = np.arccos(cp) #atan(y/x)
    
    #convert spherical coordinates into ra,dec in radians
    mock['ra']  = phi
    mock['dec'] = (math.pi/2.0) - theta

    #remove galaxies beyond r=L_box
    keep_1 = (r < Lbox)

    #calculate the luminosity distance
    dl = (1.0+z_cos)*r #flat cosmology
    #compute the apparant magnitude
    mock['app_mag'] = data['Mag']+ 5.0*np.log10(dl) + 25.0
    #calculate absolute magnitudes using z_obs
    mock['abs_mag'] = mock['app_mag'] - dist_mod(mock['z'],cosmo=cosmo)
    
    #apply apparant magnitude limit
    keep_2 = (mock['app_mag'] < app_mag_lim)
    keep = (keep_2&keep_1)
    mock = mock[keep]

    print max(mock['z']), f(Lbox), (max(mock['z'])-f(Lbox))*c

    #convert ra,dec into degrees
    mock['ra'] = mock['ra']*180.0/math.pi
    mock['dec'] = mock['dec']*180.0/math.pi

    #save output
    print 'saving ascii version of the catalogue...'
    filename = catalogue+'_radec_mock'
    print filename
    data_table = table.table.Table(data=mock)
    data_table.write(savepath+filename+'.dat', format='ascii')
    print data_table

def dist_mod(z, cosmo):
    ld = cosmology.funcs.luminosity_distance(z, cosmo=cosmo).value
    dist_mod = 5.0*(np.log10(ld*1000.0*1000.0)-1.0)
    return dist_mod

if __name__ == '__main__':
    main()
