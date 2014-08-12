#!/usr/bin/python

#Author: Duncan Campbell
#Written: July 9, 2013
#Yale University
#Description: Read in ascii group catalogues and save as HDF5 files.

###packages###
import numpy as np
from astropy.io import ascii
from astropy import table
import h5py
import gc


def main():
    ###make sure to change these when running in a new enviorment!###
    #location of data directory
    filepath = '/scratch/dac29/data/Hearin_mocks/'
    #save data to directory...
    savepath = '/scratch/dac29/output/processed_data/hearin_mocks/'
    #################################################################

    catalogue = 'Mr19_age_distribution_matching_mock'
    name = ['ID_halo','x','y','z','Vx','Vy','Vz','M_vir','V_peak','M_r,0.1',\
            'g-r','M_host','ID_host']
    print 'reading in:', catalogue
    filename = catalogue+'.dat'
    data = ascii.read(filepath+filename, delimiter='\s', names=name, \
                      guess=False, Reader=ascii.Basic, data_start=0)
    f = h5py.File(savepath+catalogue+'.hdf5', 'w')
    dset = f.create_dataset(catalogue, data=data)
    f.close()
    filename = catalogue
    data_table = table.table.Table(data=data)
    ascii.write(data_table, savepath+filename+'.dat')
    gc.collect()

    catalogue = 'sm9.8_age_matching_mock'
    name = ['ID_halo','x','y','z','Vx','Vy','Vz','M_star','g-r','M_host','ID_host']
    print 'reading in:', catalogue
    filename = catalogue+'.dat'
    data = ascii.read(filepath+filename, delimiter='\s', names=name, \
                      guess=False, Reader=ascii.Basic, data_start=0)
    f = h5py.File(savepath+catalogue+'.hdf5', 'w')
    dset = f.create_dataset(catalogue, data=data)
    f.close()
    filename = catalogue
    data_table = table.table.Table(data=data)
    ascii.write(data_table, savepath+filename+'.dat')
    gc.collect()
    
    '''
    catalogue = 'hlist_1.00030.list'
    name = ['scale', 'id', 'desc_scale', 'desc_id', 'num_prog', 'pid', 'upid',\
            'desc_pid', 'phantom', 'sam_mvir', 'mvir', 'rvir', 'rs', 'vrms',\
            'mmp?', 'scale_of_last_MM', 'vmax', 'x', 'y', 'z', 'vx', 'vy', 'vz',\
            'Jx', 'Jy', 'Jz', 'Spin', 'Breadth_first_ID', 'Depth_first_ID',\
            'Tree_root_ID', 'Orig_halo_ID', 'Snap_num', 'Next_coprogenitor_depthfirst_ID',\
            'Last_progenitor_depthfirst_ID', 'Rs_Klypin', 'Mvir_all', 'M200b', 'M200c',\
            'M500c', 'M2500c', 'Xoff', 'Voff', 'Spin_Bullock', 'b_to_a', 'c_to_a',\
            'A[x]', 'A[y]', 'A[z]', 'T/|U|', 'Macc', 'Mpeak', 'Vacc', 'Vpeak']
    print 'reading in:', catalogue
    filename = catalogue
    data = ascii.read(filepath+filename, delimiter='\s', names=name, \
                      guess=False, Reader=ascii.Basic, data_start = 0)
    f = h5py.File(savepath+catalogue+'.hdf5', 'w')
    dset = f.create_dataset(catalogue, data=data)
    f.close()
    filename = catalogue
    data_table = table.table.Table(data=data)
    ascii.write(data_table, savepath+filename+'.dat')
    gc.collect()
    '''


if __name__ == '__main__':
  main()
