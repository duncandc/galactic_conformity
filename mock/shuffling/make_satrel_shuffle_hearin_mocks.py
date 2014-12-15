#!/usr/bin/python

#Author: Duncan Campbell
#Written: December, 11 2013
#Yale University
#Description: Read in hdf5 verions of Andrew Hearin's age-matching mock catalogues.
# Shuffle satellites amongts ALL haloes in a mass bin, preserving the relative 
# position of the satellites to the central of the halo.

###packages###
import numpy as np
import h5py
import sys
import custom_utilities as cu
from astropy.io import ascii
from astropy import table

def main():
    ###make sure to change these when running in a new enviorment!###
    #location of data directory
    filepath = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
    savepath = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
    #################################################################

    catalogues=['Mr19_age_distribution_matching_mock','sm9.8_age_matching_mock']

    catalogue = sys.argv[1]
    f =  h5py.File(filepath+catalogue+'.hdf5', 'r')
    GC = f.get(catalogue)
    GC = np.array(GC)

    #make new catalogue and copy over values from original catalogue
    dtype = GC.dtype.descr
    dtype = np.dtype(dtype)
    GC_new = np.recarray((len(GC),), dtype=dtype)
    GC_new.fill(0.0)
    GC_new = np.array(GC, copy=True)

    #identify central galaxies, host haloes
    centrals = np.where(GC['ID_host'] == -1)[0] #indices of the central galaxies
    print 'number of centrals, host haloes:', len(centrals)
    satellites = np.where(GC['ID_host'] != -1)[0] #indices of the satellite galaxies
    print 'number of satellites, host haloes:', len(satellites)

    #define mass bins, and which central are in each mass bin
    mass_bins = np.arange(8.0,16.0,0.1) #log mass bins
    mass_hist, bins = np.histogram(GC['M_host'][centrals], bins=mass_bins) #group histogram by log(host_mass)
    mass_bin_ind = np.digitize(GC['M_host'][centrals], bins=mass_bins) #indices of groups in log(host_mass) bins

    #go through each mass bin
    for i in range(0,len(mass_bins)-1):
        print i, 'mass bin:', mass_bins[i], mass_bins[i+1]
        ind = np.where(mass_bin_ind==i+1)[0] #indices of host haloes in this mass bin
        if len(ind)>0: #if there are any haloes in the mass bin
            print 'number of groups:', len(ind)
            ids = GC['ID_halo'][centrals[ind]]
            sat_galaxy_members = np.in1d(GC['ID_host'],ids)      #satellite galaxies in the mass bin
            sat_galaxy_members = np.where(sat_galaxy_members)[0] #indicies of galaxies
            cen_galaxy_members = np.in1d(GC['ID_halo'],ids)      #central galaxies in the mass bin
            cen_galaxy_members = np.where(cen_galaxy_members)[0] #indicies of galaxies
            galaxy_members = np.hstack((sat_galaxy_members,cen_galaxy_members))
            print 'number of galaxies:', len(galaxy_members)
            satellite_members = np.where(GC['ID_host'][galaxy_members]!=-1)[0] #satellites
            satellite_members = galaxy_members[satellite_members]              #indices of satellites
            central_members = np.where(GC['ID_host'][galaxy_members]==-1)[0]   #centrals
            central_members = galaxy_members[central_members]                  #indices of centrals
            print 'number of centrals:',   len(central_members)
            print 'number of satellites:', len(satellite_members)
            print 'check:',  len(central_members) + len(satellite_members) == len(galaxy_members)
            #shuffle list of host haloes in mass bin
            shuffle = np.random.permutation(np.arange(0,len(central_members),1))
            shuffled_central_members = central_members[shuffle]
            unshuffle = np.arange(0,len(central_members),1)
            ran_index = np.random.random_integers(0,len(central_members)-1,len(satellite_members))
            #shuffle satellite systems --> leave gal props alone, change halo props
            for i in range(0,len(satellite_members)):
                print "\r",i,
                sys.stdout.flush()
                old_host_ID  = GC['ID_host'][satellite_members[i]]         #old host halo ID
                old_host_ind = np.where(GC['ID_halo']==old_host_ID)[0]    #index of old host central
                new_host_ind = ran_index[i] #location in central members list
                new_host_ind = central_members[new_host_ind]      #new host index
                #assign a new host properties
                GC_new['ID_host'][satellite_members[i]] = GC['ID_halo'][new_host_ind]
                GC_new['M_host'][satellite_members[i]]  = GC['M_host'][new_host_ind]
                GC_new['M_vir'][satellite_members[i]]   = GC['M_vir'][new_host_ind]
                GC_new['R200'][satellite_members[i]]    = GC['R200'][new_host_ind]
                #calculate satellite positions
                x_new_cen = GC['x'][new_host_ind]
                y_new_cen = GC['y'][new_host_ind]
                z_new_cen = GC['z'][new_host_ind]
                x_old_cen = GC['x'][old_host_ind]
                y_old_cen = GC['y'][old_host_ind]
                z_old_cen = GC['z'][old_host_ind]
                GC_new['x'][satellite_members[i]] = GC_new['x'][satellite_members[i]] - x_old_cen + x_new_cen
                GC_new['y'][satellite_members[i]] = GC_new['y'][satellite_members[i]] - y_old_cen + y_new_cen
                GC_new['z'][satellite_members[i]] = GC_new['z'][satellite_members[i]] - z_old_cen + z_new_cen
                #calculate satellite velocities
                Vx_new_cen = GC['Vx'][new_host_ind]
                Vy_new_cen = GC['Vy'][new_host_ind]
                Vz_new_cen = GC['Vz'][new_host_ind]
                Vx_old_cen = GC['Vx'][old_host_ind]
                Vy_old_cen = GC['Vy'][old_host_ind]
                Vz_old_cen = GC['Vz'][old_host_ind]
                GC_new['Vx'][satellite_members[i]] = GC_new['Vx'][satellite_members[i]] - Vx_old_cen + Vx_new_cen
                GC_new['Vy'][satellite_members[i]] = GC_new['Vy'][satellite_members[i]] - Vy_old_cen + Vy_new_cen
                GC_new['Vz'][satellite_members[i]] = GC_new['Vz'][satellite_members[i]] - Vz_old_cen + Vz_new_cen
                  
    
    #Fix any boundary condition issues, Lbox=250 Mpc
    fix = np.where(GC_new['x'] < 0.0)[0]
    GC_new['x'][fix] = 250.0 - np.absolute(GC_new['x'][fix])
    fix = np.where(GC_new['y'] < 0.0)[0]
    GC_new['y'][fix] = 250.0 - np.absolute(GC_new['y'][fix])
    fix = np.where(GC_new['z'] < 0.0)[0]
    GC_new['z'][fix] = 250.0 - np.absolute(GC_new['z'][fix])
    fix = np.where(GC_new['x'] > 250.0)[0]
    GC_new['x'][fix] = GC_new['x'][fix] - 250.0
    fix = np.where(GC_new['y'] > 250.0)[0]
    GC_new['y'][fix] = GC_new['y'][fix] - 250.0
    fix = np.where(GC_new['z'] > 250.0)[0]
    GC_new['z'][fix] = GC_new['z'][fix] - 250.0

      
    catalogue = catalogue+'_satrel_shuffle'
    print 'saving hdf5 version of the catalogue...'
    filename = catalogue
    print filename
    f = h5py.File(savepath+filename+'.hdf5', 'w')
    dset = f.create_dataset(catalogue, data=GC_new)
    f.close()

    print 'saving ascii version of the catalogue...'
    print filename
    data_table = table.table.Table(data=GC_new)
    ascii.write(data_table, savepath+filename+'.dat')
    print data_table

if __name__ == '__main__':
  main()
