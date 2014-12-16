import numpy as np
import h5py
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams
import custom_utilities as cu
from astropy import cosmology
import sys

def main():

    if len(sys.argv)>1: catalogue = sys.argv[1]
    esle: catalogue = 'sample3_L_model.mr19'
    plotpath = cu.get_plot_path() + 'analysis/groupcats/'
    
    filepath_cat = cu.get_output_path() + 'analysis/yang_groupcat/'
    filename = catalogue+'_blue_central_conformity.npy'
    arr1_1 = np.load(filepath_cat+filename)
    filename = catalogue+'_red_central_conformity.npy'
    arr2_1 = np.load(filepath_cat+filename)

    filepath_cat = cu.get_output_path() + 'analysis/berlind_groupcat/'
    filename = catalogue+'_blue_central_conformity.npy'
    arr1_2 = np.load(filepath_cat+filename)
    filename = catalogue+'_red_central_conformity.npy'
    arr2_2 = np.load(filepath_cat+filename)

    filepath_cat = cu.get_output_path() + 'analysis/berlind_groupcat/'
    filename = catalogue+'_blue_central_conformity.npy'
    arr1_3 = np.load(filepath_cat+filename)
    filename = catalogue+'_red_central_conformity.npy'
    arr2_3 = np.load(filepath_cat+filename)

    catalogue_mock = 'Mr19_age_distribution_matching_mock'
    filepath_cat = cu.get_output_path() + 'analysis/hearin_mocks/'
    filename = catalogue_mock+'_blue_central_conformity.npy'
    arr1_mock = np.load(filepath_cat+filename)
    filename = catalogue_mock+'_red_central_conformity.npy'
    arr2_mock = np.load(filepath_cat+filename)


    mass_bins = np.arange(10,15.2,0.2) #log mass bins

    fig = plt.figure() 
    #Yang
    f_cen_blue = arr1_1[0]
    N_cen_blue = arr1_1[2]
    f_cen_red = arr2_1[0]
    N_cen_red = arr2_1[2]
    selection = np.where( N_cen_blue > 50)[0]
    p1a, = plt.plot(mass_bins[0:-1][selection],f_cen_blue[selection], color='blue',linewidth=2)
    selection = np.where( N_cen_red > 50)[0]
    p2a, = plt.plot(mass_bins[0:-1][selection],f_cen_red[selection], color='red',linewidth=2 )
    p3a, = plt.plot([0,1],[-1,-1], color='black',linewidth=2 )
    

    #Berlind
    f_cen_blue = arr1_2[0]
    N_cen_blue = arr1_2[2]
    f_cen_red = arr2_2[0]
    N_cen_red = arr2_2[2]
    selection = np.where( N_cen_blue > 50)[0]
    p1b, = plt.plot(mass_bins[0:-1][selection], f_cen_blue[selection], color='blue', alpha=0.5,linewidth=2)
    selection = np.where( N_cen_red > 50)[0]
    p2b, = plt.plot(mass_bins[0:-1][selection], f_cen_red[selection], color='red', alpha=0.5,linewidth=2)
    p3b, = plt.plot([0,1],[-1,-1], color='black',alpha=0.5,linewidth=2 )

    #Wetzel
    f_cen_blue = arr1_3[0]
    N_cen_blue = arr1_3[2]
    f_cen_red = arr2_3[0]
    N_cen_red = arr2_3[2]
    selection = np.where( N_cen_blue > 50)[0]
    p1c, = plt.plot(mass_bins[0:-1][selection], f_cen_blue[selection], color='blue', ls='--',linewidth=2)
    selection = np.where( N_cen_red > 50)[0]
    p2c, = plt.plot(mass_bins[0:-1][selection], f_cen_red[selection], color='red', ls='--',linewidth=2)
    p3c, = plt.plot([0,1],[-1,-1], color='black', ls='--', linewidth=2 )

    #Hearin Mock
    '''
    f_cen_blue = arr1_mock[0]
    N_cen_blue = arr1_mock[2]
    f_cen_red = arr2_mock[0]
    N_cen_red = arr2_mock[2]
    selection = np.where( N_cen_blue > 50)[0]
    p1d, = plt.plot(mass_bins[0:-1][selection], f_cen_blue[selection], color='green', lw=3)
    selection = np.where( N_cen_red > 50)[0]
    p2d, = plt.plot(mass_bins[0:-1][selection], f_cen_red[selection], color='magenta', lw=3)
    '''

    font = {'size':19}
    matplotlib.rc('font', **font)

    plt.ylim([0,0.8])
    plt.xlim([11,15])
    plt.xlabel(r'$\mathrm{log}(M)$ $[{h}^{-1}{M}_{\odot}]$')
    plt.ylabel(r'$f_{\mathrm{late}}$',fontsize=25)
    #plt.title(catalogue)
    #plt.legend([p1a, p1b, p1c, p1d], ["Yang blue centrals", "Berlind blue centrals", "Tinker blue centrals", "mock blue centrals"])
    #plt.legend([p1a, p1b, p1c], ["Yang et al. blue centrals", "Berlind et al. blue centrals", "Tinker et al. blue centrals"])
    plt.legend([p3a, p3b, p3c], ["Yang et al.", "Berlind et al.", "Tinker et al."],prop={'size':20})
    filename = catalogue+'_compare_conformity'
    fig.savefig(plotpath+filename+'.eps')
    plt.show(block=True)

if __name__ == '__main__':
    main() 
