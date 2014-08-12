This directory contains python scripts to shuffle Andrew Hearin's age-matching
mock catalogues.  The different shuffeling proceduers and the associated scripts
are listed in this file.

Type A shuffeling uses only primary haloes that were previously occupied by a 
galaxy in the unaltered age-matching mock catalogue. These shuffles are:
1. satellite shuffeling in sub-haloes
2. satellite shuffeling 
3. central shuffeling
4. satellite system shuffeling
5. central-satellite system shuffeling

Type B shuffeling uses all haloes with a mass greater than or equal to the 
minimum mass of all occupied haloes in the unaltered age-matching mock 
catalogue.  All subhaloes are included that have a max_mass greater than or 
equal to the minimum max_mass of all occupied sub-haloes in the unaltered 
age-matching mock.  These shuffels are:
1. satellite shuffeling in sub-haloes
2. satellite shuffeling
3. central shuffeling
4. central and satellite shuffeling in sub-haloes
5. central and satellite shuffeling
6. central shuffeling and conditional satellite shuffeling in sub-haloes
7. central shuffeling and conditional satellite shuffeling
8  satellite system shuffeling
9. central-satellite system shuffeling

Descriptions of shuffelings(All shuffeling is done in primary halo mass bins):
1A: All satellites in a mass bin are shuffeld amongst occupied sub-haloes in
    the mass bin.
2A: All centrals in a mass bin are shuffeld amongst occupied primary haloes 
    in the mass bin.
3A: All satllites systems in a mass bin are shuffeld amongst primary haloes 
    in the mass bin, preserving the relative position of the satellites to the
    center of their host.
4A: All galaxy systems as a whole are shuffeld in a mass bin amongst all
    occupied primary haloes in the mass bin, preserving the relative positions 
    of the galaxies of any individual ssystem. 
1B: All satellites in a mass bin are shuffeld amongst all sub-haloes in the
    mass bin.
2B: All centrals are shuffeld in a mass bin amongst all primary haloes in the
    mass bin.
3B: Centrals and satellites in a mass bin are shuffeld indpendantly amongst all 
    primary haloes and subhaloes in the mass bin.
4B: Centrals in a mass bin are shuffeld indpendantly amongst all primary 
    haloes in the mass bin. The satellites in the mass bin are shuffeld to a
    halo, conditional on that halo having a satellite.
5B: All satllites systems in the mass bin are shuffeld amongst primary haloes in    the mass bin preserving the relative position of the satellites to the 
    center of their host.
6B: All galaxy systems as a whole are shuffeld in a mass bin amongst all
    occupied primary haloes in the mass bin, preserving the relative positions
    of the galaxies of any individual system. 

	    
