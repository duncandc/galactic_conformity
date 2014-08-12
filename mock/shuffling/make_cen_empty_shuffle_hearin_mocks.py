#!/usr/bin/python

#Author: Duncan Campbell
#Written: August 14, 2013
#Yale University
#Description: Read in Andrew hearin's Age-matching mock catalogues.
#  All centrals in a mass bin are shuffeld amongst occupied primary
#  haloes in the mass bin. 

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
    print 'running central shuffle for:', catalogue
    f =  h5py.File(filepath+catalogue+'.hdf5', 'r')
    GC = f.get(catalogue)

    #make new catalogue
    dtype = GC.dtype.descr
    dtype = np.dtype(dtype)
    GC_new = np.recarray((len(GC),), dtype=dtype)
    GC_new = np.array(GC,copy=True)


    #identify centrals
    centrals = np.where(GC['ID_host'] == -1)[0] #central galaxies
    print 'number of centrals:', len(centrals)

    #define mass bins, and which central are in each mass bin
    mass_bins = np.arange(8.0,14.0,0.1) #log mass bins
    mass_hist, bins = np.histogram(GC['M_host'][centrals], bins=mass_bins) #group histogram by log(mass)
    mass_bin_ind = np.digitize(GC['M_host'][centrals], bins=mass_bins) #indices of groups in log(mass) bins

    #go through each mass bin
    for i in range(0,len(mass_bins)-1):
        print 'mass bin:', i, mass_bins[i], mass_bins[i+1]
        ind = np.where(mass_bin_ind==i+1)[0]
        if len(ind)>0: #if there are any haloes in the mass bin
            print 'number of groups:', len(ind)
            ids = GC['ID_halo'][centrals[ind]]
            sat_galaxy_members = np.in1d(GC['ID_host'],ids) #satellite galaxies in the mass bin
            sat_galaxy_members = np.where(sat_galaxy_members)[0] #indicies of galaxies
            cen_galaxy_members = np.in1d(GC['ID_halo'],ids) #central galaxies in the mass bin
            cen_galaxy_members = np.where(cen_galaxy_members)[0] #indicies of galaxies
            galaxy_members = np.hstack((sat_galaxy_members,cen_galaxy_members))
            print 'number of galaxies:', len(galaxy_members)
            satellite_members = np.where(GC['ID_host'][galaxy_members]!=-1)[0] #satellites
            satellite_members = galaxy_members[satellite_members] #indices of satellites
            central_members = np.where(GC['ID_host'][galaxy_members]==-1)[0] #centrals
            central_members = galaxy_members[central_members] #indices of centrals
            print 'number of centrals:', len(central_members)
            print 'number of satellites:', len(satellite_members)
            print 'check:',  len(central_members) + len(satellite_members) == len(galaxy_members)
            #now shuffle centrals amongst groups in mass bin
            shuffle = np.random.permutation(np.arange(0,len(central_members),1))
            GC_new['M_r,0.1'][central_members] = GC['M_r,0.1'][central_members[shuffle]]
            GC_new['M_star'][central_members] = GC['M_star'][central_members[shuffle]]
            GC_new['g-r'][central_members] = GC['g-r'][central_members[shuffle]]

    #only use haloes that were occupied in the original mocks
    #remove empty haloes
    keep = np.where(GC_new['g-r']!=-99)[0]
    dtype = GC.dtype.descr
    dtype = np.dtype(dtype)
    GC_new_short = np.recarray((len(keep),), dtype=dtype)
    GC_new_short = np.array(GC_new[keep], copy=True)

    print 'saving hdf5 version of the catalogue...'
    filename = catalogue+'_cen_shuffle'
    f = h5py.File(savepath+filename+'.hdf5', 'w')
    dset = f.create_dataset(catalogue, data=GC_new_short)
    f.close()

    print 'saving ascii version of the catalogue...'
    filename = catalogue+'_cen_shuffle'
    data_table = table.table.Table(data=GC_new_short)
    ascii.write(data_table, savepath+filename+'.dat')
    print data_table

if __name__ == '__main__':
  main()
