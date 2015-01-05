#!/usr/bin/env python

#Duncan Campbell
#April 8, 2014
#Yale University
#examine the velocity distribution of satellite galaxies

#load packages
import numpy as np
import h5py
import matplotlib.pyplot as plt
import custom_utilities as cu
import sys

def main():
    
    if len(sys.argv)>1: catalogue = sys.argv[1]
    else: catalogue = 'Mr19_age_distribution_matching_mock'
    
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
    print 'opening mock catalogue:', catalogue+'.hdf5'
    #open catalogue
    f1 = h5py.File(filepath_mock+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f1.get(catalogue)

    host = np.where(mock['ID_host']==-1)[0]
    sub = np.where(mock['ID_host']!=-1)[0]
    
    ind1, ind2 = cu.match(mock['ID_host'][sub],mock['ID_halo'][host])
    dv = np.zeros(len(mock))
    dv[sub[ind1]] = mock['Vx'][sub[ind1]] - mock['Vx'][host[ind2]]
    
    color = mock['g-r']
    LHS   = 0.7 - 0.032*(mock['M_r,0.1']+16.5) #Weinmann 2006
    blue  = np.where((color<LHS) & (mock['M_r,0.1']!=-99))[0]
    red   = np.where((color>LHS) & (mock['M_r,0.1']!=-99))[0]
    
    red_host = np.in1d(host,red)
    red_host = host[red_host]
    blue_host = np.in1d(host,blue)
    blue_host = host[blue_host]
    
    ind1, ind2 = cu.match(mock['ID_host'][sub],mock['ID_halo'][host])
    sat_w_host = sub[ind1]
    sat_w_host_host = host[ind2]
    
    ind1_r, ind2_r = cu.match(mock['ID_host'][sub],mock['ID_halo'][red_host])
    sat_w_red_host = sub[ind1_r]
    sat_w_red_host_host = red_host[ind2_r]
    
    ind1_b, ind2_b = cu.match(mock['ID_host'][sub],mock['ID_halo'][blue_host])
    sat_w_blue_host = sub[ind1_b]
    sat_w_blue_host_host = blue_host[ind2_b]
    
    '''
    plt.figure()
    plt.plot(mock['M_host'],mock['Vx'],'.')
    plt.show()
    '''
    
    w=0.2
    bins=np.arange(11.5,15,w)
    bin_centers = (bins[:-1]+bins[1:])/2.0
    result = np.digitize(mock['M_host'],bins=bins)
    mean_red = np.zeros(len(bins)-1)
    std_red = np.zeros(len(bins)-1)
    mean_blue = np.zeros(len(bins)-1)
    std_blue = np.zeros(len(bins)-1)
    for i in range(0,len(bins)-1):
        inds = np.where(result==i+1)[0]
        red_inds = np.in1d(inds,sat_w_red_host)
        red_inds = inds[red_inds]
        blue_inds = np.in1d(inds,sat_w_blue_host)
        blue_inds = inds[blue_inds]
        mean_red[i] = np.mean(dv[red_inds])
        std_red[i] = np.std(dv[red_inds])
        mean_blue[i] = np.mean(dv[blue_inds])
        std_blue[i] = np.std(dv[blue_inds])
    
    plt.figure()
    plt.plot(mock['M_host'][sat_w_red_host],dv[sat_w_red_host],'.', ms=3, color='grey', alpha=0.1)
    plt.plot(mock['M_host'][sat_w_blue_host],dv[sat_w_blue_host],'.', ms=3, color='grey', alpha=0.1)
    plt.plot(bin_centers,mean_red+std_red,'-', color='red')
    plt.plot(bin_centers,mean_red-std_red,'-', color='red')
    plt.plot(bin_centers,mean_red,'o-', color='red')
    #plt.errorbar(bin_centers,mean_red,yerr=std_red, fmt='o-', color='red', linewidth=3)
    plt.plot(bin_centers,mean_blue-std_blue,'-', color='blue')
    plt.plot(bin_centers,mean_blue+std_blue,'-', color='blue')
    plt.plot(bin_centers,mean_blue,'o-', color='blue')
    #plt.errorbar(bin_centers+w/2.0,mean_blue,yerr=std_blue, fmt='o-', color='blue', linewidth=3)
    plt.xlim([11,15])
    plt.xlabel(r'M host')
    plt.ylabel(r'$\Delta$ v km/s')
    plt.show()
    
    w=0.05
    bins=np.arange(2,4,w)
    bin_centers = (bins[:-1]+bins[1:])/2.0
    result = np.digitize(np.log10(mock['V_peak']),bins=bins)
    mean_red = np.zeros(len(bins)-1)
    std_red = np.zeros(len(bins)-1)
    mean_blue = np.zeros(len(bins)-1)
    std_blue = np.zeros(len(bins)-1)
    for i in range(0,len(bins)-1):
        inds = np.where(result==i+1)[0]
        red_inds = np.in1d(inds,sat_w_red_host_host)
        red_inds = inds[red_inds]
        blue_inds = np.in1d(inds,sat_w_blue_host_host)
        blue_inds = inds[blue_inds]
        mean_red[i] = np.mean(dv[red_inds])
        std_red[i] = np.std(dv[red_inds])
        mean_blue[i] = np.mean(dv[blue_inds])
        std_blue[i] = np.std(dv[blue_inds])
    
    plt.figure()
    plt.plot(mock['V_peak'][sat_w_red_host_host],dv[sat_w_red_host],'.', ms=3, color='red', alpha=0.2)
    plt.plot(mock['V_peak'][sat_w_blue_host_host],dv[sat_w_blue_host],'.', ms=3, color='blue', alpha=0.2)
    '''
    plt.plot(10**(bin_centers),mean_red+std_red,'-', color='red')
    plt.plot(10**(bin_centers),mean_red-std_red,'-', color='red')
    plt.plot(10**(bin_centers),mean_red,'o-', color='red')
    #plt.errorbar(10**(bin_centers),mean_red,yerr=std_red, fmt='o-', color='red')
    plt.plot(10**(bin_centers),mean_blue+std_blue,'-', color='blue')
    plt.plot(10**(bin_centers),mean_blue-std_blue,'-', color='blue')
    plt.plot(10**(bin_centers),mean_blue,'o-', color='blue')
    #plt.errorbar(10**(bin_centers+w/2.0),mean_blue,yerr=std_blue, fmt='o-', color='blue')
    '''
    plt.xlim([100,1500])
    plt.xscale('log')
    plt.xlabel(r'V peak km/s')
    plt.ylabel(r'$\Delta$ v km/s')
    plt.show()
    
    
    
    
    
    
if __name__ == '__main__':
    main()
    
    