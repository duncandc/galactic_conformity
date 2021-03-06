#!/usr/bin/python

#Author: Duncan Campbell
#Written: August 14, 2013
#Yale University
#Description: Read in hdf5 hearin mock catalogues and scramble satelite galaxies among groups in halo mass bins

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
    filepath = cu.get_output_path() + 'processed_data/hearin_mocks/'
    savepath = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
    #################################################################

    catalogues=['Mr19_age_distribution_matching_mock','sm9.8_age_matching_mock']

    catalogue = sys.argv[1]
    print 'running satellite and central shuffle on:', catalogue
    f =  h5py.File(filepath+catalogue+'.hdf5', 'r')
    GC = f.get(catalogue)

    #only use haloes that were occupied in the original mocks
    #remove empty haloes
    keep = np.where(GC_new['g-r']!=-99)[0]
    GC=GC[keep]

    print GC.dtype
    #make new catalogue
    dtype = GC.dtype.descr
    dtype = np.dtype(dtype)
    GC_new = np.recarray((len(GC),), dtype=dtype)
    GC_new = np.array(GC,copy=True)

    #identify centrals
    centrals = np.where(GC['ID_host'] == -1)[0] #central galaxies
    print 'number of centrals:', len(centrals)

    #define mass bins, and which central are in each mass bin
    mass_bins = np.arange(8.0,16.0,0.1) #log mass bins
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
            new_central_members = np.random.permutation(central_members)
            GC_new['M_r,0.1'][central_members] = GC['M_r,0.1'][new_central_members]
            GC_new['M_star'][central_members] = GC['M_star'][new_central_members]
            GC_new['g-r'][central_members] = GC['g-r'][new_central_members]
            #now shuffle satellites amongst groups in mass bin
            shuffle = np.random.permutation(np.arange(0,len(satellite_members),1))
            if len(satellite_members)>0:
                GC_new['M_r,0.1'][satellite_members] = GC['M_r,0.1'][satellite_members[shuffle]]
                GC_new['M_star'][satellite_members] = GC['M_star'][satellite_members[shuffle]]
                GC_new['g-r'][satellite_members] = GC['g-r'][satellite_members[shuffle]]

    print 'saving hdf5 version of the catalogue...'
    filename = catalogue+'_satcen_shuffle'
    f = h5py.File(savepath+filename+'.hdf5', 'w')
    dset = f.create_dataset(catalogue, data=GC_new)
    f.close()

    print 'saving ascii version of the catalogue...'
    filename = catalogue+'_satcen_shuffle'
    data_table = table.table.Table(data=GC_new)
    ascii.write(data_table, savepath+filename+'.dat')
    print data_table

if __name__ == '__main__':
  main()
