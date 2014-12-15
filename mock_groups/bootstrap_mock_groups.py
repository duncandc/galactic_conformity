#!/usr/bin/env python

import numpy as np
import h5py
import sys

def main():

    group = sys.argv[1]
    catalogue = sys.argv[2]
    seed = int(sys.argv[3])
    np.random.seed(seed)
     
    filepath = '/scratch/dac29/output/processed_data/'+group+'_groupcat/mock_runs/4th_run/custom_catalogues/'
    savepath = '/scratch/dac29/output/processed_data/'+group+'_groupcat/mock_runs/4th_run/custom_catalogues/bootstraps/'

    print 'opening group catalogue:', catalogue
    #open catalogue
    f  =  h5py.File(filepath+catalogue+'.hdf5', 'r') #open catalogue file
    GC = f.get(catalogue)
    GC = np.array(GC)

    #create a new array to store bootstrapped version.  Make it too big.
    data = np.recarray((len(GC)*2,), dtype=GC.dtype)
    data.fill(-99)
    
    #get central/satellite indices
    centrals = np.where(GC['RANK']==1)[0]
    satellites = np.where(GC['RANK']==0)[0]
    
    #randomly sample the centrals(groups)
    sample = np.random.random_integers(0,len(centrals)-1,len(centrals))

    #loop through the catalogue
    prev_ind = 0
    for i in range(0,len(sample)):
        group_ID  = GC['GROUP_ID'][centrals[sample[i]]]
        members   = np.where(GC['GROUP_ID']==group_ID)[0]
        N_members = len(members)
        new_group_ID = i
        data[prev_ind:prev_ind+N_members] = GC[members]
        data['GROUP_ID'][prev_ind:prev_ind+N_members] = new_group_ID
        prev_ind = prev_ind + N_members

    #remove empty entries
    keep = (data['GROUP_ID']!=-99)
    data = data[keep]

    print 'saving hdf5 version of the catalogue...'
    filename = catalogue+'_'+str(seed)
    print filename
    f = h5py.File(savepath+filename+'.hdf5', 'w')
    dset = f.create_dataset(filename, data=data)
    f.close()

   
if __name__ == '__main__':
    main() 
