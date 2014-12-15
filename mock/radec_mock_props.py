#!/usr/bin/python

#Author: Duncan Campbell
#Written: December, 2014
#Yale University
#Description: look at the properties of ra dec mock

###packages###
import numpy as np
from astropy.io import ascii
import h5py
import sys
import custom_utilities as cu
import matplotlib.pyplot as plt
from astropy import cosmology
from scipy.interpolate import interp1d

def main():

    if len(sys.argv)>1:
        catalogue = sys.argv[1]
    else: catalogue = 'Mr19_age_distribution_matching_mock'

    #open the ra,dec mock
    filepath   = cu.get_output_path()+'processed_data/hearin_mocks/custom_catalogues/'
    filename   = catalogue+'_radec_mock.dat'
    mock_radec = ascii.read(filepath+filename, delimiter='\s', Reader=ascii.Basic, data_start=1)
    mock_radec = np.array(mock_radec)
    for name in mock_radec.dtype.names: print '\t', name
    
    #open full mock
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
    print 'opening mock catalogue:', catalogue+'.hdf5'
    #open catalogue
    f1 = h5py.File(filepath_mock+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f1.get(catalogue)
    mock = np.array(mock)
    for name in mock.dtype.names: print '\t', name
    
    app_mag_lim = 17.77
    Lbox = 250.0 #Mpc/h
    cosmo = cosmology.core.FlatLambdaCDM(H0=100.0,Om0=0.27) #h=1
    c = 299792.458 #speed of light in km/s
    
    #determine abs mag limit as afunction of z to overplot
    app_mag_r_lim = 17.77 
    z_bins        = np.arange(1.0,0.001,-0.001)
    dm            = dist_mod(z_bins,cosmo)
    abs_mag_r_lim = app_mag_r_lim - dm
    f_r           = interp1d(abs_mag_r_lim, z_bins, kind='cubic', bounds_error=False, fill_value=max(z_bins))
    
    z_lim =  f_r(-19)
    print z_lim
    
    plt.figure()
    plt.plot(mock_radec['z'],mock['M_r,0.1'][mock_radec['k']],'.',ms=2)
    plt.plot(z_bins,abs_mag_r_lim,color='green')
    plt.plot([0,z_lim],[-19,-19],color='red')
    plt.plot([z_lim,z_lim],[-19,-24],color='red')
    plt.plot([0.02,0.02],[-19,-24],color='red')
    plt.plot([0.068,0.068],[-19,-24],'--',color='red')
    plt.ylim([-18.5,-23])
    plt.xlim([0,0.1])
    plt.show()

def dist_mod(z, cosmo):
    ld = cosmo.luminosity_distance(z).value
    dist_mod = 5.0*(np.log10(ld*1000.0*1000.0)-1.0)
    return dist_mod

if __name__ == '__main__':
    main()