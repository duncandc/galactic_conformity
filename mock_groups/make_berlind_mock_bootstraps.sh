#!/bin/bash

catalogues=('Mr19_age_distribution_matching_mock_groups' \
            'Mr19_age_distribution_matching_mock_sys_empty_shuffle_satrel_shuffle_groups')

for catalogue in ${catalogues[*]}
do
    wait
    echo $catalogue

    python bootstrap_mock_groups.py berlind $catalogue 0 & 
    python bootstrap_mock_groups.py berlind $catalogue 1 &
    python bootstrap_mock_groups.py berlind $catalogue 2 &
    python bootstrap_mock_groups.py berlind $catalogue 3 &
    python bootstrap_mock_groups.py berlind $catalogue 4 &
    python bootstrap_mock_groups.py berlind $catalogue 5 &
    python bootstrap_mock_groups.py berlind $catalogue 6 &
    python bootstrap_mock_groups.py berlind $catalogue 7 &
    python bootstrap_mock_groups.py berlind $catalogue 8 &
    python bootstrap_mock_groups.py berlind $catalogue 9 &
    python bootstrap_mock_groups.py berlind $catalogue 10
 
    python bootstrap_mock_groups.py berlind $catalogue 11 &
    python bootstrap_mock_groups.py berlind $catalogue 12 &
    python bootstrap_mock_groups.py berlind $catalogue 13 &
    python bootstrap_mock_groups.py berlind $catalogue 14 &
    python bootstrap_mock_groups.py berlind $catalogue 15 &
    python bootstrap_mock_groups.py berlind $catalogue 16 &
    python bootstrap_mock_groups.py berlind $catalogue 17 &
    python bootstrap_mock_groups.py berlind $catalogue 18 &
    python bootstrap_mock_groups.py berlind $catalogue 19 &
    python bootstrap_mock_groups.py berlind $catalogue 20 

    python bootstrap_mock_groups.py berlind $catalogue 21 &
    python bootstrap_mock_groups.py berlind $catalogue 22 &
    python bootstrap_mock_groups.py berlind $catalogue 23 &
    python bootstrap_mock_groups.py berlind $catalogue 24 &
    python bootstrap_mock_groups.py berlind $catalogue 25 &
    python bootstrap_mock_groups.py berlind $catalogue 26 &
    python bootstrap_mock_groups.py berlind $catalogue 27 &
    python bootstrap_mock_groups.py berlind $catalogue 28 &
    python bootstrap_mock_groups.py berlind $catalogue 29 &
    python bootstrap_mock_groups.py berlind $catalogue 30

    python bootstrap_mock_groups.py berlind $catalogue 31 &
    python bootstrap_mock_groups.py berlind $catalogue 32 &
    python bootstrap_mock_groups.py berlind $catalogue 33 &
    python bootstrap_mock_groups.py berlind $catalogue 34 &
    python bootstrap_mock_groups.py berlind $catalogue 35 &
    python bootstrap_mock_groups.py berlind $catalogue 36 &
    python bootstrap_mock_groups.py berlind $catalogue 37 &
    python bootstrap_mock_groups.py berlind $catalogue 38 &
    python bootstrap_mock_groups.py berlind $catalogue 39 &
    python bootstrap_mock_groups.py berlind $catalogue 40 

    python bootstrap_mock_groups.py berlind $catalogue 41 &
    python bootstrap_mock_groups.py berlind $catalogue 42 &
    python bootstrap_mock_groups.py berlind $catalogue 43 &
    python bootstrap_mock_groups.py berlind $catalogue 44 &
    python bootstrap_mock_groups.py berlind $catalogue 45 &
    python bootstrap_mock_groups.py berlind $catalogue 46 &
    python bootstrap_mock_groups.py berlind $catalogue 47 &
    python bootstrap_mock_groups.py berlind $catalogue 48 &
    python bootstrap_mock_groups.py berlind $catalogue 49 &
    python bootstrap_mock_groups.py berlind $catalogue 50 
done
