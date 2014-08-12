#!/usr/bin/python

#Author: Duncan Campbell
#Written: August 14, 2013
#Yale University
#Description: Read in Andrew Hearin's age matching mock catalogue.
#  Add empty haloes to the catalogue if they have M_host >= min(M_host)
#  greater than the minimum included in the unaltered mock.
#  Also add R_200, R_proj, r to each entry
#  Also add an empty column for stellar mass (for uniformity with the sm
#  catalogue)

###packages###
import numpy as np
import math
import h5py
import sys
import custom_utilities as cu
from astropy.io import ascii
from astropy import table

def main():
    ###make sure to change these when running in a new enviorment!###
    #location of data directory
    filepath = cu.get_output_path() + 'processed_data/hearin_mocks/'
    savepath = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
    #################################################################

    halo_catalogue = 'hlist_1.00030.list'
    mock_catalogue = sys.argv[1]
    #mock_catalogue='Mr19_age_distribution_matching_mock'

    #read in the mock catalogue
    f =  h5py.File(filepath+'custom_catalogues/'+mock_catalogue+'.hdf5', 'r')
    GC = f.get(mock_catalogue) #mock catalogue
    GC = np.array(GC)
    for name in GC.dtype.names: print name
    #read in the halo catalogue
    f =  h5py.File(filepath+halo_catalogue+'.hdf5', 'r')
    HC = f.get(halo_catalogue) #halo catalogue
    HC = np.array(HC)
    for name in HC.dtype.names: print name

    #identify central/satellite galaxies in the mock
    centrals   = np.where(GC['ID_host'] == -1)[0] #central galaxies (primary haloes)
    satellites = np.where(GC['ID_host'] != -1)[0] #satellite galaxies (subhaloes)
    print 'number of centrals:',   len(centrals)
    print 'number of satellites:', len(satellites)

    #idenitify haloes in the halo catalogue that are present in the mock
    host_haloes = np.where(HC['upid'] == -1)[0]
    subhaloes   = np.where(HC['pid']  != -1)[0]
    print 'number of host haloes in halo catalogue:', len(host_haloes)
    print 'number of subhaloes in halo catalogue:',   len(subhaloes)
    mock_halo_inds = np.in1d(HC['id'][host_haloes],GC['ID_halo'][centrals])
    mock_halo_inds = host_haloes[mock_halo_inds] #indices of haloes in halo catalogue
    mock_halo_ids  = HC['id'][mock_halo_inds] #halo ids
    print 'number of host haloes in mock:', len(mock_halo_ids)
    print 'test:', len(mock_halo_ids)==len(centrals)
    
    min_halo_mass = np.min(HC['M200b'][mock_halo_inds]) #minimum halo mass in the mock
    print 'minimum log halo mass:', np.log10(min_halo_mass)

    #determine which haloes should be included given the minimum halo mass
    mass_lim_haloes = np.where(((HC['M200b'] >= min_halo_mass) & (HC['pid']==-1)))[0]
    print 'number of haloes >= min_mass:', len(mass_lim_haloes)
    halo_ids = HC['id'][mass_lim_haloes] #IDs of haloes >= minimum halo mass

    #idenitfy which haloes are ***not*** already in the mock
    missing_halo_inds = np.in1d(halo_ids, mock_halo_ids)
    missing_halo_inds = np.where(missing_halo_inds == False)[0]
    missing_halo_inds = mass_lim_haloes[missing_halo_inds]
    #missing_halo_ids  = halo_ids[missing_halo_inds]
    missing_halo_ids  = HC['id'][missing_halo_inds]
    print 'number of haloes to be added to mock:', len(missing_halo_ids)

    #make new catalogue
    dtype = [('ID_halo', '<i8'), ('x', '<f8'), ('y', '<f8'), ('z', '<f8'), ('Vx', '<f8'), ('Vy', '<f8'),\
             ('Vz', '<f8'), ('M_vir', '<f8'), ('V_peak', '<f8'), ('M_r,0.1', '<f8'), ('M_star', '<f8'),\
             ('g-r', '<f8'), ('M_host', '<f8'), ('ID_host', '<i8'), ('r', '<f8'), ('R_proj', '<f8'),\
             ('R200', '<f8')]
    print dtype
    dtype = np.dtype(dtype)
    N_entries = len(GC)+len(missing_halo_inds)
    GC_new = np.recarray((N_entries,), dtype=dtype)
    GC_new.fill(-99) #empty value indicator
    #fill with the mock
    ind_begin = 0
    ind_end   = len(GC)
    GC_new['ID_halo'][ind_begin:ind_end] = np.array(GC['ID_halo'], copy=True)
    GC_new['x'][ind_begin:ind_end]       = np.array(GC['x'], copy=True)
    GC_new['y'][ind_begin:ind_end]       = np.array(GC['y'], copy=True)
    GC_new['z'][ind_begin:ind_end]       = np.array(GC['z'], copy=True)
    GC_new['Vx'][ind_begin:ind_end]      = np.array(GC['Vx'], copy=True)
    GC_new['Vy'][ind_begin:ind_end]      = np.array(GC['Vy'], copy=True)
    GC_new['Vz'][ind_begin:ind_end]      = np.array(GC['Vz'], copy=True)
    GC_new['M_vir'][ind_begin:ind_end]   = np.array(GC['M_vir'], copy=True)
    GC_new['V_peak'][ind_begin:ind_end]  = np.array(GC['V_peak'], copy=True)
    GC_new['M_r,0.1'][ind_begin:ind_end] = np.array(GC['M_r,0.1'], copy=True)
    GC_new['g-r'][ind_begin:ind_end]     = np.array(GC['g-r'], copy=True)
    GC_new['M_host'][ind_begin:ind_end]  = np.array(GC['M_host'], copy=True)
    GC_new['ID_host'][ind_begin:ind_end] = np.array(GC['ID_host'], copy=True)
    GC_new['M_star'][ind_begin:ind_end] = np.array(GC['M_star'], copy=True)
    GC_new['R200'][ind_begin:ind_end] = np.array(GC['R200'], copy=True)
    GC_new['r'][ind_begin:ind_end] = np.array(GC['r'], copy=True)
    GC_new['R_proj'][ind_begin:ind_end] = np.array(GC['R_proj'], copy=True)
    #add empty haloes
    ind_begin = len(GC)
    ind_end   = len(GC)+len(missing_halo_inds)
    GC_new['ID_halo'][ind_begin:ind_end] = HC['id'][missing_halo_inds]
    GC_new['ID_host'][ind_begin:ind_end] = HC['pid'][missing_halo_inds]
    GC_new['x'][ind_begin:ind_end]       = HC['x'][missing_halo_inds]
    GC_new['y'][ind_begin:ind_end]       = HC['y'][missing_halo_inds]
    GC_new['z'][ind_begin:ind_end]       = HC['z'][missing_halo_inds]
    GC_new['Vx'][ind_begin:ind_end]      = HC['vx'][missing_halo_inds]
    GC_new['Vy'][ind_begin:ind_end]      = HC['vy'][missing_halo_inds]
    GC_new['Vz'][ind_begin:ind_end]      = HC['vz'][missing_halo_inds]
    GC_new['M_host'][ind_begin:ind_end]  = np.log10(HC['M200b'][missing_halo_inds])
    GC_new['M_vir'][ind_begin:ind_end]   = np.log10(HC['mvir'][missing_halo_inds])
    GC_new['V_peak'][ind_begin:ind_end]  = HC['Vpeak'][missing_halo_inds]
    
    
    print 'saving hdf5 version of the extended catalogue...'
    filename = mock_catalogue+'_extended'
    print savepath+filename+'.hdf5'
    f = h5py.File(savepath+filename+'.hdf5', 'w')
    dset = f.create_dataset(mock_catalogue+'_extended', data=GC_new)
    f.close()

    print 'saving ascii version of the extended catalogue...'
    filename = mock_catalogue+'_extended'
    print savepath+filename+'.dat'
    data_table = table.table.Table(data=GC_new)
    ascii.write(data_table, savepath+filename+'.dat')
    print data_table

if __name__ == '__main__':
  main()
