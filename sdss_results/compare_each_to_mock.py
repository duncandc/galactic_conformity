#!/usr/bin/env python

import numpy as np
import h5py
import matplotlib.pyplot as plt
import custom_utilities as cu
import sys

def main():
    
    if len(sys.argv)>1: mock_1 = sys.argv[1]
    else: mock_1 = 'Mr19_age_distribution_matching_mock'
    if len(sys.argv)>2: group_cat = sys.argv[2]
    else: group_cat = 'yang_groupcat'
    
    #set up plotting space
    fig, axes = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True, figsize=(6.95, 2.66))
    fig.subplots_adjust(hspace=0, wspace=0.05)
    fig.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.9)
    axes = axes.flatten()
    #plot 1, left
    ax = axes[0]
    ax.set_ylim([0,1])
    ax.set_xlim([12,15])
    ax.set_ylabel(r'$f_{\rm blue}$')
   # ax.set_xlabel(r'$M~M_{\odot}h^{-1}$')
    #plot 2, middle
    ax = axes[1]
    ax.set_ylim([0,1])
    ax.set_xlim([10**12,10**15])
    ax.set_xlabel(r'$\log[M/(h^{-1}M_{\odot})]$')
    #plot 3, right
    ax = axes[2]
    ax.set_ylim([0,1])
    ax.set_xlim([12,15])
    #ax.set_xlabel(r'$M~M_{\odot}h^{-1}$')
    ax.set_xticks([12,13,14,15])
    ax.set_xticklabels([r'$12$', r'$13$', r'$14$',r' '])
    
    mass_bins = np.arange(10,15.2,0.2) #log mass bins
    
    #open age-matching mock
    filepath = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
    print 'opening mock catalogue:', mock_1+'.hdf5'
    f1 = h5py.File(filepath+mock_1+'.hdf5', 'r') #open catalogue file
    GC = f1.get(mock_1)
    
    #galaxy color
    h=1.0
    color = GC['g-r']
    #LHS = 0.7 - 0.032*(GC['M_r,0.1']-5.0*np.log10(h)+16.5) #Weinmann 2006
    LHS = 0.21-0.03*GC['M_r,0.1']
    blue = np.where(color<LHS)[0] #indices of blue galaxies
    red = np.where(color>LHS)[0] #indicies of red galaxies
    lookup = color<LHS #true if blue, false if red
    flip = (lookup == False)

    centrals = np.where(GC['ID_host'] == -1)[0] #central galaxies
    N_groups = len(centrals)
    
    mass_hist, bins = np.histogram(GC['M_host'][centrals], bins=mass_bins) #group histogram by log(mass)
    mass_bin_ind = np.digitize(GC['M_host'][centrals], bins=mass_bins) #indices of groups in log(mass) bins
    selection = np.where(lookup==True)[0]
    selection = np.intersect1d(selection,centrals)
    mass_hist_blue, bins = np.histogram(GC['M_host'][selection], bins=mass_bins) #group histogram by log(mass)
    mass_bin_ind_blue = np.digitize(GC['M_host'][selection], bins=mass_bins) #indices of groups in log(mass) bins
    selection = np.where(lookup==False)[0]
    selection = np.intersect1d(selection,centrals)
    mass_hist_red, bins = np.histogram(GC['M_host'][selection], bins=mass_bins) #group histogram by log(mass)
    mass_bin_ind_red = np.digitize(GC['M_host'][selection], bins=mass_bins) #indices of groups in log(mass) bins
    width = 1.0 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    #normalize the hists
    mass_hist_red = mass_hist_red/float(N_groups)
    mass_hist_blue = mass_hist_blue/float(N_groups)

    f_cen_blue   = np.zeros(len(mass_bins)-1) #fraction of blue satellites for blue centrals
    f_cen_red    = np.zeros(len(mass_bins)-1) #fraction of blue satellites for red centrals
    err_cen_blue = np.zeros(len(mass_bins)-1) #poisson error on blue fraction for blue centrals
    err_cen_red  = np.zeros(len(mass_bins)-1) #poisson error on blue fraction for red centrals
    N_cen_blue   = np.zeros(len(mass_bins)-1) #number of satellites for blue centrals
    N_cen_red    = np.zeros(len(mass_bins)-1) #number of satellites for red centrals
    N_groups_blue= np.zeros(len(mass_bins)-1) #number of satellites for red centrals
    N_groups_red = np.zeros(len(mass_bins)-1) #number of satellites for red centrals
    for i in range(0,len(mass_bins)-1):
        print i, mass_bins[i], mass_bins[i+1]
        ind = np.where(mass_bin_ind==i+1)[0]
        if len(ind)>0:
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
            blue_central_members = np.where(lookup[central_members]==True)[0] #blue centrals
            blue_central_members = central_members[blue_central_members] #indicies of blue centrals
            red_central_members = np.where(lookup[central_members]==False)[0] #red centrals
            red_central_members = central_members[red_central_members] #indicies of red centrals

            print 'number of blue centrals:', len(blue_central_members)
            print 'number of red centrals:', len(red_central_members) 
            print 'check:', len(blue_central_members)+len(red_central_members) == len(central_members)  
            
            blue_central_satellites = np.in1d(GC['ID_host'][satellite_members],GC['ID_halo'][blue_central_members])
            blue_central_satellites = np.where(blue_central_satellites)[0]
            blue_central_satellites = satellite_members[blue_central_satellites]
            red_central_satellites = np.in1d(GC['ID_host'][satellite_members],GC['ID_halo'][red_central_members])
            red_central_satellites = np.where(red_central_satellites)[0]
            red_central_satellites = satellite_members[red_central_satellites]
            
            print 'number of blue central satellites:', len(blue_central_satellites)
            print 'number of red central satellites:', len(red_central_satellites)
            print 'check:', len(blue_central_satellites) + len(red_central_satellites) == len(satellite_members)

            blue_central_blue_satellites = np.where(lookup[blue_central_satellites]==True)[0]
            blue_central_blue_satellites = blue_central_satellites[blue_central_blue_satellites]
            blue_central_red_satellites = np.where(lookup[blue_central_satellites]==False)[0]
            blue_central_red_satellites = blue_central_satellites[blue_central_red_satellites]
            red_central_blue_satellites = np.where(lookup[red_central_satellites]==True)[0]
            red_central_blue_satellites = red_central_satellites[red_central_blue_satellites]
            red_central_red_satellites = np.where(lookup[red_central_satellites]==False)[0]
            red_central_red_satellites = red_central_satellites[red_central_red_satellites]

            N_blue_central_blue_satellites  = float(len(blue_central_blue_satellites))
            N_blue_central_red_satellites = float(len(blue_central_red_satellites))
            N_red_central_blue_satellites  = float(len(red_central_blue_satellites))
            N_red_central_red_satellites = float(len(red_central_red_satellites))
            
            print 'check:', N_blue_central_blue_satellites+N_blue_central_red_satellites==len(blue_central_satellites)
            print 'check:', N_red_central_blue_satellites+N_red_central_red_satellites==len(red_central_satellites)

            if len(blue_central_satellites)>0: 
                f_cen_blue[i] = N_blue_central_blue_satellites / \
                             (N_blue_central_blue_satellites + N_blue_central_red_satellites)
                print 'f_cen_blue:', f_cen_blue[i]
            if len(red_central_satellites)>0:
                f_cen_red[i] = N_red_central_blue_satellites / \
                            (N_red_central_blue_satellites + N_red_central_red_satellites)
                print 'f_cen_red:', f_cen_red[i]

            #number satellites in the bin
            N_cen_blue[i] = float(len(blue_central_satellites))
            N_cen_red[i] = float(len(red_central_satellites))
            #number of centrals in the bin
            N_groups_blue[i] = float(len(blue_central_members))
            N_groups_red[i] = float(len(red_central_members))
   

    ax=axes[0]
    selection = np.where( N_cen_blue > 50)[0]
    p0a_mock, = ax.plot(mass_bins[0:-1][selection],f_cen_blue[selection], '-', color='blue', alpha=1, ms=4, mec='blue', mfc='none')
    p0a_mock_fake, = ax.plot(10*mass_bins[0:-1][selection],f_cen_blue[selection], '-', color='black', alpha=1, ms=4, mec='black', mfc='none')
    selection = np.where( N_cen_red > 50)[0]
    p0b_mock, = ax.plot(mass_bins[0:-1][selection],f_cen_red[selection], '-', color='red', alpha=1, ms=4, mec='red', mfc='none')
    p0b_mock_fake, = ax.plot(10**mass_bins[0:-1][selection],f_cen_red[selection], '-', color='black', alpha=1, ms=4, mec='black', mfc='none')
    ax=axes[1]
    selection = np.where( N_cen_blue > 50)[0]
    p0a_mock, = ax.plot(mass_bins[0:-1][selection],f_cen_blue[selection], '-', color='blue', alpha=1, ms=4, mec='blue', mfc='none')
    p0a_mock_fake, = ax.plot(10*mass_bins[0:-1][selection],f_cen_blue[selection], '-', color='black', alpha=1, ms=4, mec='black', mfc='none')
    selection = np.where( N_cen_red > 50)[0]
    p0b_mock, = ax.plot(mass_bins[0:-1][selection],f_cen_red[selection], '-', color='red', alpha=1, ms=4, mec='red', mfc='none')
    p0b_mock_fake, = ax.plot(10**mass_bins[0:-1][selection],f_cen_red[selection], '-', color='black', alpha=1, ms=4, mec='black', mfc='none')
    ax=axes[2]
    selection = np.where( N_cen_blue > 50)[0]
    p2a, = ax.plot(mass_bins[0:-1][selection],f_cen_blue[selection], '-', color='blue', alpha=1, ms=4, mec='blue', mfc='none')
    p2a_mock_fake, = ax.plot(10*mass_bins[0:-1][selection],f_cen_blue[selection], '-', color='black', alpha=1, ms=4, mec='black', mfc='none')
    selection = np.where( N_cen_red > 50)[0]
    p2b_mock, = ax.plot(mass_bins[0:-1][selection],f_cen_red[selection], '-', color='red', alpha=1, ms=4, mec='red', mfc='none')
    p2b_mock_fake, = ax.plot(10**mass_bins[0:-1][selection],f_cen_red[selection], '-', color='black', alpha=1, ms=4, mec='black', mfc='none')
    
    
    group_cat = 'yang'
    sample = 'sample3_L_model.mr19'
    filepath_cat = cu.get_output_path() + 'processed_data/yang_groupcat/custom_catalogues/'

    catalogue = sample

    N_boots = 50

    f_cen_blue   = np.zeros((N_boots,len(mass_bins)-1)) #fraction of blue satellites for blue centrals
    f_cen_red    = np.zeros((N_boots,len(mass_bins)-1)) #fraction of blue satellites for red centrals
    N_cen_blue   = np.zeros((N_boots,len(mass_bins)-1)) #number of satellites for blue centrals
    N_cen_red    = np.zeros((N_boots,len(mass_bins)-1)) #number of satellites for red centrals
    N_groups_blue= np.zeros((N_boots,len(mass_bins)-1)) #number of satellites for red centrals
    N_groups_red = np.zeros((N_boots,len(mass_bins)-1)) #number of satellites for red centrals

    for boot in range(0,N_boots):

        print 'opening group catalogue:', catalogue
        #open catalogue
        filepath = filepath_cat+'bootstraps/'
        catalogue_1 = catalogue+'_'+str(boot)
        print catalogue_1
        f =  h5py.File(filepath+catalogue_1+'.hdf5', 'r') #open catalogue file
        GC = f.get(catalogue_1)
        print len(GC)

        #do we have a measure of M_u?
        good = np.where(GC['M_g,0.1']!=-99.9)[0]
        bad  = np.where(GC['M_g,0.1']==-99.9)[0]
        print 'number of galaxies without M_u:', len(bad)    

        #galaxy color
        h=1.0
        color = GC['M_g,0.1'] - GC['M_r,0.1']
        #LHS = 0.7 - 0.032*(GC['M_r,0.1']-5.0*np.log10(h)+16.5) #Weinmann 2006
        LHS = 0.21-0.03*GC['M_r,0.1']
        blue = np.where(color<LHS)[0] #indices of blue galaxies
        red = np.where(color>LHS)[0] #indicies of red galaxies
        lookup = color<LHS #true if blue, false if red
        flip = (lookup == False)

        print 'number of blue galaxies:', len(blue)
        print 'number of red galaxies:', len(red)

        centrals = np.where(GC['RPROJ'][good] == 0)[0] #central galaxies
        centrals = good[centrals]
        bad_centrals = np.where(GC['RPROJ'][bad] == 0)[0] #central galaxies with ssfr measured
        bad_centrals = bad[bad_centrals]
        N_groups = len(centrals)+len(bad_centrals)
    
        mass_hist, bins = np.histogram(GC['MGROUP'][centrals], bins=mass_bins) #group histogram by log(mass)
        mass_bin_ind = np.digitize(GC['MGROUP'][centrals], bins=mass_bins) #indices of groups in log(mass) bins
        selection = np.where(lookup==True)[0]
        selection = np.intersect1d(selection,centrals)
        mass_hist_blue, bins = np.histogram(GC['MGROUP'][selection], bins=mass_bins) #group histogram by log(mass)
        mass_bin_ind_blue = np.digitize(GC['MGROUP'][selection], bins=mass_bins) #indices of groups in log(mass) bins
        selection = np.where(lookup==False)[0]
        selection = np.intersect1d(selection,centrals)
        mass_hist_red, bins = np.histogram(GC['MGROUP'][selection], bins=mass_bins) #group histogram by log(mass)
        mass_bin_ind_red = np.digitize(GC['MGROUP'][selection], bins=mass_bins) #indices of groups in log(mass) bins
        selection = bad_centrals
        mass_hist_grey, bins = np.histogram(GC['MGROUP'][selection], bins=mass_bins) #group histogram by log(mass)
        if len(bad_centrals)>0:
            mass_bin_ind_grey = np.digitize(GC['MGROUP'][selection], bins=mass_bins) #indices of groups in log(mass) bins
        width = 1.0 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        #normalize the hists
        mass_hist_grey = mass_hist_grey/float(N_groups)
        mass_hist_red = mass_hist_red/float(N_groups)
        mass_hist_blue = mass_hist_blue/float(N_groups)

    
        for i in range(0,len(mass_bins)-1):
            print i, mass_bins[i], mass_bins[i+1]
            ind = np.where(mass_bin_ind==i+1)[0]
            if len(ind)>0:
                print 'number of groups:', len(ind)
                ids = GC['GROUP_ID'][centrals[ind]]
                galaxy_members = np.in1d(GC['GROUP_ID'],ids) #galaxies in the mass bin
                galaxy_members = np.where(galaxy_members)[0] #indicies of galaxies
                print 'number of galaxies:', len(galaxy_members)
                satellite_members = np.where(GC['RPROJ'][galaxy_members]>0)[0] #satellites
                satellite_members = galaxy_members[satellite_members] #indices of satellites
                central_members = np.where(GC['RPROJ'][galaxy_members]==0)[0] #centrals
                central_members = galaxy_members[central_members] #indices of centrals
                print 'number of centrals:', len(central_members)
                print 'number of satellites:', len(satellite_members)
                print 'check:',  len(central_members) + len(satellite_members) == len(galaxy_members)
                blue_central_members = np.where(lookup[central_members]==True)[0] #blue centrals
                blue_central_members = central_members[blue_central_members] #indicies of blue centrals
                red_central_members = np.where(lookup[central_members]==False)[0] #red centrals
                red_central_members = central_members[red_central_members] #indicies of red centrals

                print 'number of blue centrals:', len(blue_central_members)
                print 'number of red centrals:', len(red_central_members) 
                print 'check:', len(blue_central_members)+len(red_central_members) == len(central_members)  
            
                blue_central_satellites = np.in1d(GC['GROUP_ID'][satellite_members],GC['GROUP_ID'][blue_central_members])
                blue_central_satellites = np.where(blue_central_satellites)[0]
                blue_central_satellites = satellite_members[blue_central_satellites]
                red_central_satellites = np.in1d(GC['GROUP_ID'][satellite_members],GC['GROUP_ID'][red_central_members])
                red_central_satellites = np.where(red_central_satellites)[0]
                red_central_satellites = satellite_members[red_central_satellites]
            
                print 'number of blue central satellites:', len(blue_central_satellites)
                print 'number of red central satellites:', len(red_central_satellites)
                print 'check:', len(blue_central_satellites) + len(red_central_satellites) == len(satellite_members)

                blue_central_blue_satellites = np.where(lookup[blue_central_satellites]==True)[0]
                blue_central_blue_satellites = blue_central_satellites[blue_central_blue_satellites]
                blue_central_red_satellites = np.where(lookup[blue_central_satellites]==False)[0]
                blue_central_red_satellites = blue_central_satellites[blue_central_red_satellites]
                red_central_blue_satellites = np.where(lookup[red_central_satellites]==True)[0]
                red_central_blue_satellites = red_central_satellites[red_central_blue_satellites]
                red_central_red_satellites = np.where(lookup[red_central_satellites]==False)[0]
                red_central_red_satellites = red_central_satellites[red_central_red_satellites]

                N_blue_central_blue_satellites  = float(len(blue_central_blue_satellites))
                N_blue_central_red_satellites = float(len(blue_central_red_satellites))
                N_red_central_blue_satellites  = float(len(red_central_blue_satellites))
                N_red_central_red_satellites = float(len(red_central_red_satellites))
            
                print 'check:', N_blue_central_blue_satellites+N_blue_central_red_satellites==len(blue_central_satellites)
                print 'check:', N_red_central_blue_satellites+N_red_central_red_satellites==len(red_central_satellites)

                if len(blue_central_satellites)>0: 
                    f_cen_blue[boot,i] = N_blue_central_blue_satellites / \
                             (N_blue_central_blue_satellites + N_blue_central_red_satellites)
                    print 'f_cen_blue:', f_cen_blue[boot,i]
                if len(red_central_satellites)>0:
                    f_cen_red[boot,i] = N_red_central_blue_satellites / \
                            (N_red_central_blue_satellites + N_red_central_red_satellites)
                    print 'f_cen_red:', f_cen_red[boot,i]

                #number satellites in the bin
                N_cen_blue[boot,i] = float(len(blue_central_satellites))
                N_cen_red[boot,i] = float(len(red_central_satellites))
                #number of centrals in the bin
                N_groups_blue[boot,i] = float(len(blue_central_members))
                N_groups_red[boot,i] = float(len(red_central_members))

    f_blue = np.mean(f_cen_blue, axis=0)
    f_red = np.mean(f_cen_red, axis=0)
    N_blue = np.mean(N_cen_blue, axis=0)
    N_red = np.mean(N_cen_red, axis=0)
    err_blue = np.std(f_cen_blue, axis=0)
    err_red = np.std(f_cen_red, axis=0)

    print f_blue
    print err_blue
    print N_blue
    
    ax=axes[0]
    selection = np.where( N_blue > 50)[0]
    p2a_yang = ax.errorbar(mass_bins[0:-1][selection],f_blue[selection],yerr=err_blue[selection], fmt='o', color='blue', ms=4, mec='blue', mfc='none')
    p2a_yang_fake = ax.errorbar(10**mass_bins[0:-1][selection],f_blue[selection],yerr=err_blue[selection], fmt='o', color='black', ms=4, mec='black', mfc='none')
    selection = np.where( N_red > 50)[0]
    p2b = ax.errorbar(mass_bins[0:-1][selection],f_red[selection], yerr=err_red[selection], fmt='o', color='red', ms=4, mec='red', mfc='none')
    ax.legend((p2a_mock_fake, p2a_yang_fake),('age-matching','group finder on\n SDSS'), fontsize=10, numpoints=1, frameon=False )
    
    group_cat = 'tinker'
    sample = 'sample3_L_model.mr19'
    filepath_cat = cu.get_output_path() + 'processed_data/tinker_groupcat/custom_catalogues/'

    catalogue = sample

    N_boots = 50

    f_cen_blue   = np.zeros((N_boots,len(mass_bins)-1)) #fraction of blue satellites for blue centrals
    f_cen_red    = np.zeros((N_boots,len(mass_bins)-1)) #fraction of blue satellites for red centrals
    N_cen_blue   = np.zeros((N_boots,len(mass_bins)-1)) #number of satellites for blue centrals
    N_cen_red    = np.zeros((N_boots,len(mass_bins)-1)) #number of satellites for red centrals
    N_groups_blue= np.zeros((N_boots,len(mass_bins)-1)) #number of satellites for red centrals
    N_groups_red = np.zeros((N_boots,len(mass_bins)-1)) #number of satellites for red centrals

    for boot in range(0,N_boots):

        print 'opening group catalogue:', catalogue
        #open catalogue
        filepath = filepath_cat+'bootstraps/'
        catalogue_1 = catalogue+'_'+str(boot)
        print catalogue_1
        f =  h5py.File(filepath+catalogue_1+'.hdf5', 'r') #open catalogue file
        GC = f.get(catalogue_1)
        print len(GC)

        #do we have a measure of M_u?
        good = np.where(GC['M_g,0.1']!=-99.9)[0]
        bad  = np.where(GC['M_g,0.1']==-99.9)[0]
        print 'number of galaxies without M_u:', len(bad)    

        #galaxy color
        h=1.0
        color = GC['M_g,0.1'] - GC['M_r,0.1']
        #LHS = 0.7 - 0.032*(GC['M_r,0.1']-5.0*np.log10(h)+16.5) #Weinmann 2006
        LHS = 0.21-0.03*GC['M_r,0.1']
        blue = np.where(color<LHS)[0] #indices of blue galaxies
        red = np.where(color>LHS)[0] #indicies of red galaxies
        lookup = color<LHS #true if blue, false if red
        flip = (lookup == False)

        print 'number of blue galaxies:', len(blue)
        print 'number of red galaxies:', len(red)

        centrals = np.where(GC['RPROJ'][good] == 0)[0] #central galaxies
        centrals = good[centrals]
        bad_centrals = np.where(GC['RPROJ'][bad] == 0)[0] #central galaxies with ssfr measured
        bad_centrals = bad[bad_centrals]
        N_groups = len(centrals)+len(bad_centrals)
    
        mass_hist, bins = np.histogram(GC['MGROUP'][centrals], bins=mass_bins) #group histogram by log(mass)
        mass_bin_ind = np.digitize(GC['MGROUP'][centrals], bins=mass_bins) #indices of groups in log(mass) bins
        selection = np.where(lookup==True)[0]
        selection = np.intersect1d(selection,centrals)
        mass_hist_blue, bins = np.histogram(GC['MGROUP'][selection], bins=mass_bins) #group histogram by log(mass)
        mass_bin_ind_blue = np.digitize(GC['MGROUP'][selection], bins=mass_bins) #indices of groups in log(mass) bins
        selection = np.where(lookup==False)[0]
        selection = np.intersect1d(selection,centrals)
        mass_hist_red, bins = np.histogram(GC['MGROUP'][selection], bins=mass_bins) #group histogram by log(mass)
        mass_bin_ind_red = np.digitize(GC['MGROUP'][selection], bins=mass_bins) #indices of groups in log(mass) bins
        selection = bad_centrals
        mass_hist_grey, bins = np.histogram(GC['MGROUP'][selection], bins=mass_bins) #group histogram by log(mass)
        if len(bad_centrals)>0:
            mass_bin_ind_grey = np.digitize(GC['MGROUP'][selection], bins=mass_bins) #indices of groups in log(mass) bins
        width = 1.0 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        #normalize the hists
        mass_hist_grey = mass_hist_grey/float(N_groups)
        mass_hist_red = mass_hist_red/float(N_groups)
        mass_hist_blue = mass_hist_blue/float(N_groups)

    
        for i in range(0,len(mass_bins)-1):
            print i, mass_bins[i], mass_bins[i+1]
            ind = np.where(mass_bin_ind==i+1)[0]
            if len(ind)>0:
                print 'number of groups:', len(ind)
                ids = GC['GROUP_ID'][centrals[ind]]
                galaxy_members = np.in1d(GC['GROUP_ID'],ids) #galaxies in the mass bin
                galaxy_members = np.where(galaxy_members)[0] #indicies of galaxies
                print 'number of galaxies:', len(galaxy_members)
                satellite_members = np.where(GC['RPROJ'][galaxy_members]>0)[0] #satellites
                satellite_members = galaxy_members[satellite_members] #indices of satellites
                central_members = np.where(GC['RPROJ'][galaxy_members]==0)[0] #centrals
                central_members = galaxy_members[central_members] #indices of centrals
                print 'number of centrals:', len(central_members)
                print 'number of satellites:', len(satellite_members)
                print 'check:',  len(central_members) + len(satellite_members) == len(galaxy_members)
                blue_central_members = np.where(lookup[central_members]==True)[0] #blue centrals
                blue_central_members = central_members[blue_central_members] #indicies of blue centrals
                red_central_members = np.where(lookup[central_members]==False)[0] #red centrals
                red_central_members = central_members[red_central_members] #indicies of red centrals

                print 'number of blue centrals:', len(blue_central_members)
                print 'number of red centrals:', len(red_central_members) 
                print 'check:', len(blue_central_members)+len(red_central_members) == len(central_members)  
            
                blue_central_satellites = np.in1d(GC['GROUP_ID'][satellite_members],GC['GROUP_ID'][blue_central_members])
                blue_central_satellites = np.where(blue_central_satellites)[0]
                blue_central_satellites = satellite_members[blue_central_satellites]
                red_central_satellites = np.in1d(GC['GROUP_ID'][satellite_members],GC['GROUP_ID'][red_central_members])
                red_central_satellites = np.where(red_central_satellites)[0]
                red_central_satellites = satellite_members[red_central_satellites]
            
                print 'number of blue central satellites:', len(blue_central_satellites)
                print 'number of red central satellites:', len(red_central_satellites)
                print 'check:', len(blue_central_satellites) + len(red_central_satellites) == len(satellite_members)

                blue_central_blue_satellites = np.where(lookup[blue_central_satellites]==True)[0]
                blue_central_blue_satellites = blue_central_satellites[blue_central_blue_satellites]
                blue_central_red_satellites = np.where(lookup[blue_central_satellites]==False)[0]
                blue_central_red_satellites = blue_central_satellites[blue_central_red_satellites]
                red_central_blue_satellites = np.where(lookup[red_central_satellites]==True)[0]
                red_central_blue_satellites = red_central_satellites[red_central_blue_satellites]
                red_central_red_satellites = np.where(lookup[red_central_satellites]==False)[0]
                red_central_red_satellites = red_central_satellites[red_central_red_satellites]

                N_blue_central_blue_satellites  = float(len(blue_central_blue_satellites))
                N_blue_central_red_satellites = float(len(blue_central_red_satellites))
                N_red_central_blue_satellites  = float(len(red_central_blue_satellites))
                N_red_central_red_satellites = float(len(red_central_red_satellites))
            
                print 'check:', N_blue_central_blue_satellites+N_blue_central_red_satellites==len(blue_central_satellites)
                print 'check:', N_red_central_blue_satellites+N_red_central_red_satellites==len(red_central_satellites)

                if len(blue_central_satellites)>0: 
                    f_cen_blue[boot,i] = N_blue_central_blue_satellites / \
                             (N_blue_central_blue_satellites + N_blue_central_red_satellites)
                    print 'f_cen_blue:', f_cen_blue[boot,i]
                if len(red_central_satellites)>0:
                    f_cen_red[boot,i] = N_red_central_blue_satellites / \
                            (N_red_central_blue_satellites + N_red_central_red_satellites)
                    print 'f_cen_red:', f_cen_red[boot,i]

                #number satellites in the bin
                N_cen_blue[boot,i] = float(len(blue_central_satellites))
                N_cen_red[boot,i] = float(len(red_central_satellites))
                #number of centrals in the bin
                N_groups_blue[boot,i] = float(len(blue_central_members))
                N_groups_red[boot,i] = float(len(red_central_members))

    f_blue = np.mean(f_cen_blue, axis=0)
    f_red = np.mean(f_cen_red, axis=0)
    N_blue = np.mean(N_cen_blue, axis=0)
    N_red = np.mean(N_cen_red, axis=0)
    err_blue = np.std(f_cen_blue, axis=0)
    err_red = np.std(f_cen_red, axis=0)

    print f_blue
    print err_blue
    print N_blue
    
    ax=axes[1]
    selection = np.where( N_blue > 50)[0]
    p2a_yang = ax.errorbar(mass_bins[0:-1][selection],f_blue[selection],yerr=err_blue[selection], fmt='o', color='blue', ms=4, mec='blue', mfc='none')
    p2a_yang_fake = ax.errorbar(10**mass_bins[0:-1][selection],f_blue[selection],yerr=err_blue[selection], fmt='o', color='black', ms=4, mec='black', mfc='none')
    selection = np.where( N_red > 50)[0]
    p2b = ax.errorbar(mass_bins[0:-1][selection],f_red[selection], yerr=err_red[selection], fmt='o', color='red', ms=4, mec='red', mfc='none')
    
    group_cat = 'berlind'
    sample = 'sample3_L_model.mr19'
    filepath_cat = cu.get_output_path() + 'processed_data/berlind_groupcat/custom_catalogues/'

    catalogue = sample

    N_boots = 50

    f_cen_blue   = np.zeros((N_boots,len(mass_bins)-1)) #fraction of blue satellites for blue centrals
    f_cen_red    = np.zeros((N_boots,len(mass_bins)-1)) #fraction of blue satellites for red centrals
    N_cen_blue   = np.zeros((N_boots,len(mass_bins)-1)) #number of satellites for blue centrals
    N_cen_red    = np.zeros((N_boots,len(mass_bins)-1)) #number of satellites for red centrals
    N_groups_blue= np.zeros((N_boots,len(mass_bins)-1)) #number of satellites for red centrals
    N_groups_red = np.zeros((N_boots,len(mass_bins)-1)) #number of satellites for red centrals

    for boot in range(0,N_boots):

        print 'opening group catalogue:', catalogue
        #open catalogue
        filepath = filepath_cat+'bootstraps/'
        catalogue_1 = catalogue+'_'+str(boot)
        print catalogue_1
        f =  h5py.File(filepath+catalogue_1+'.hdf5', 'r') #open catalogue file
        GC = f.get(catalogue_1)
        print len(GC)

        #do we have a measure of M_u?
        good = np.where(GC['M_g,0.1']!=-99.9)[0]
        bad  = np.where(GC['M_g,0.1']==-99.9)[0]
        print 'number of galaxies without M_u:', len(bad)    

        #galaxy color
        h=1.0
        color = GC['M_g,0.1'] - GC['M_r,0.1']
        #LHS = 0.7 - 0.032*(GC['M_r,0.1']-5.0*np.log10(h)+16.5) #Weinmann 2006
        LHS = 0.21-0.03*GC['M_r,0.1']
        blue = np.where(color<LHS)[0] #indices of blue galaxies
        red = np.where(color>LHS)[0] #indicies of red galaxies
        lookup = color<LHS #true if blue, false if red
        flip = (lookup == False)

        print 'number of blue galaxies:', len(blue)
        print 'number of red galaxies:', len(red)

        centrals = np.where(GC['RPROJ'][good] == 0)[0] #central galaxies
        centrals = good[centrals]
        bad_centrals = np.where(GC['RPROJ'][bad] == 0)[0] #central galaxies with ssfr measured
        bad_centrals = bad[bad_centrals]
        N_groups = len(centrals)+len(bad_centrals)
    
        mass_hist, bins = np.histogram(GC['MGROUP'][centrals], bins=mass_bins) #group histogram by log(mass)
        mass_bin_ind = np.digitize(GC['MGROUP'][centrals], bins=mass_bins) #indices of groups in log(mass) bins
        selection = np.where(lookup==True)[0]
        selection = np.intersect1d(selection,centrals)
        mass_hist_blue, bins = np.histogram(GC['MGROUP'][selection], bins=mass_bins) #group histogram by log(mass)
        mass_bin_ind_blue = np.digitize(GC['MGROUP'][selection], bins=mass_bins) #indices of groups in log(mass) bins
        selection = np.where(lookup==False)[0]
        selection = np.intersect1d(selection,centrals)
        mass_hist_red, bins = np.histogram(GC['MGROUP'][selection], bins=mass_bins) #group histogram by log(mass)
        mass_bin_ind_red = np.digitize(GC['MGROUP'][selection], bins=mass_bins) #indices of groups in log(mass) bins
        selection = bad_centrals
        mass_hist_grey, bins = np.histogram(GC['MGROUP'][selection], bins=mass_bins) #group histogram by log(mass)
        if len(bad_centrals)>0:
            mass_bin_ind_grey = np.digitize(GC['MGROUP'][selection], bins=mass_bins) #indices of groups in log(mass) bins
        width = 1.0 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        #normalize the hists
        mass_hist_grey = mass_hist_grey/float(N_groups)
        mass_hist_red = mass_hist_red/float(N_groups)
        mass_hist_blue = mass_hist_blue/float(N_groups)

    
        for i in range(0,len(mass_bins)-1):
            print i, mass_bins[i], mass_bins[i+1]
            ind = np.where(mass_bin_ind==i+1)[0]
            if len(ind)>0:
                print 'number of groups:', len(ind)
                ids = GC['GROUP_ID'][centrals[ind]]
                galaxy_members = np.in1d(GC['GROUP_ID'],ids) #galaxies in the mass bin
                galaxy_members = np.where(galaxy_members)[0] #indicies of galaxies
                print 'number of galaxies:', len(galaxy_members)
                satellite_members = np.where(GC['RPROJ'][galaxy_members]>0)[0] #satellites
                satellite_members = galaxy_members[satellite_members] #indices of satellites
                central_members = np.where(GC['RPROJ'][galaxy_members]==0)[0] #centrals
                central_members = galaxy_members[central_members] #indices of centrals
                print 'number of centrals:', len(central_members)
                print 'number of satellites:', len(satellite_members)
                print 'check:',  len(central_members) + len(satellite_members) == len(galaxy_members)
                blue_central_members = np.where(lookup[central_members]==True)[0] #blue centrals
                blue_central_members = central_members[blue_central_members] #indicies of blue centrals
                red_central_members = np.where(lookup[central_members]==False)[0] #red centrals
                red_central_members = central_members[red_central_members] #indicies of red centrals

                print 'number of blue centrals:', len(blue_central_members)
                print 'number of red centrals:', len(red_central_members) 
                print 'check:', len(blue_central_members)+len(red_central_members) == len(central_members)  
            
                blue_central_satellites = np.in1d(GC['GROUP_ID'][satellite_members],GC['GROUP_ID'][blue_central_members])
                blue_central_satellites = np.where(blue_central_satellites)[0]
                blue_central_satellites = satellite_members[blue_central_satellites]
                red_central_satellites = np.in1d(GC['GROUP_ID'][satellite_members],GC['GROUP_ID'][red_central_members])
                red_central_satellites = np.where(red_central_satellites)[0]
                red_central_satellites = satellite_members[red_central_satellites]
            
                print 'number of blue central satellites:', len(blue_central_satellites)
                print 'number of red central satellites:', len(red_central_satellites)
                print 'check:', len(blue_central_satellites) + len(red_central_satellites) == len(satellite_members)

                blue_central_blue_satellites = np.where(lookup[blue_central_satellites]==True)[0]
                blue_central_blue_satellites = blue_central_satellites[blue_central_blue_satellites]
                blue_central_red_satellites = np.where(lookup[blue_central_satellites]==False)[0]
                blue_central_red_satellites = blue_central_satellites[blue_central_red_satellites]
                red_central_blue_satellites = np.where(lookup[red_central_satellites]==True)[0]
                red_central_blue_satellites = red_central_satellites[red_central_blue_satellites]
                red_central_red_satellites = np.where(lookup[red_central_satellites]==False)[0]
                red_central_red_satellites = red_central_satellites[red_central_red_satellites]

                N_blue_central_blue_satellites  = float(len(blue_central_blue_satellites))
                N_blue_central_red_satellites = float(len(blue_central_red_satellites))
                N_red_central_blue_satellites  = float(len(red_central_blue_satellites))
                N_red_central_red_satellites = float(len(red_central_red_satellites))
            
                print 'check:', N_blue_central_blue_satellites+N_blue_central_red_satellites==len(blue_central_satellites)
                print 'check:', N_red_central_blue_satellites+N_red_central_red_satellites==len(red_central_satellites)

                if len(blue_central_satellites)>0: 
                    f_cen_blue[boot,i] = N_blue_central_blue_satellites / \
                             (N_blue_central_blue_satellites + N_blue_central_red_satellites)
                    print 'f_cen_blue:', f_cen_blue[boot,i]
                if len(red_central_satellites)>0:
                    f_cen_red[boot,i] = N_red_central_blue_satellites / \
                            (N_red_central_blue_satellites + N_red_central_red_satellites)
                    print 'f_cen_red:', f_cen_red[boot,i]

                #number satellites in the bin
                N_cen_blue[boot,i] = float(len(blue_central_satellites))
                N_cen_red[boot,i] = float(len(red_central_satellites))
                #number of centrals in the bin
                N_groups_blue[boot,i] = float(len(blue_central_members))
                N_groups_red[boot,i] = float(len(red_central_members))

    f_blue = np.mean(f_cen_blue, axis=0)
    f_red = np.mean(f_cen_red, axis=0)
    N_blue = np.mean(N_cen_blue, axis=0)
    N_red = np.mean(N_cen_red, axis=0)
    err_blue = np.std(f_cen_blue, axis=0)
    err_red = np.std(f_cen_red, axis=0)

    print f_blue
    print err_blue
    print N_blue
    
    ax=axes[2]
    selection = np.where( N_blue > 50)[0]
    p2a_yang = ax.errorbar(mass_bins[0:-1][selection],f_blue[selection],yerr=err_blue[selection], fmt='o', color='blue', ms=4, mec='blue', mfc='none')
    p2a_yang_fake = ax.errorbar(10**mass_bins[0:-1][selection],f_blue[selection],yerr=err_blue[selection], fmt='o', color='black', ms=4, mec='black', mfc='none')
    selection = np.where( N_red > 50)[0]
    p2b = ax.errorbar(mass_bins[0:-1][selection],f_red[selection], yerr=err_red[selection], fmt='o', color='red', ms=4, mec='red', mfc='none')
    
    plt.show()
    
    filepath = cu.get_plot_path()+'analysis/groupcats/'
    filename = 'mock_all_sdss_conformity_comparison.eps'
    fig.savefig(filepath+filename)
    
    
if __name__ == '__main__':
    main() 