#!/bin/bash

catalogues=('Mr19_age_distribution_matching_mock_groups' \
            'Mr19_age_distribution_matching_mock_sys_empty_shuffle_satrel_shuffle_groups')

for catalogue in ${catalogues[*]}
do
    wait
    echo $catalogue

    python bootstrap_mock_groups.py yang $catalogue 0 & 
    python bootstrap_mock_groups.py yang $catalogue 1 &
    python bootstrap_mock_groups.py yang $catalogue 2 &
    python bootstrap_mock_groups.py yang $catalogue 3 &
    python bootstrap_mock_groups.py yang $catalogue 4 &
    python bootstrap_mock_groups.py yang $catalogue 5 &
    python bootstrap_mock_groups.py yang $catalogue 6 &
    python bootstrap_mock_groups.py yang $catalogue 7 &
    python bootstrap_mock_groups.py yang $catalogue 8 &
    python bootstrap_mock_groups.py yang $catalogue 9 &
    python bootstrap_mock_groups.py yang $catalogue 10
 
    python bootstrap_mock_groups.py yang $catalogue 11 &
    python bootstrap_mock_groups.py yang $catalogue 12 &
    python bootstrap_mock_groups.py yang $catalogue 13 &
    python bootstrap_mock_groups.py yang $catalogue 14 &
    python bootstrap_mock_groups.py yang $catalogue 15 &
    python bootstrap_mock_groups.py yang $catalogue 16 &
    python bootstrap_mock_groups.py yang $catalogue 17 &
    python bootstrap_mock_groups.py yang $catalogue 18 &
    python bootstrap_mock_groups.py yang $catalogue 19 &
    python bootstrap_mock_groups.py yang $catalogue 20 

    python bootstrap_mock_groups.py yang $catalogue 21 &
    python bootstrap_mock_groups.py yang $catalogue 22 &
    python bootstrap_mock_groups.py yang $catalogue 23 &
    python bootstrap_mock_groups.py yang $catalogue 24 &
    python bootstrap_mock_groups.py yang $catalogue 25 &
    python bootstrap_mock_groups.py yang $catalogue 26 &
    python bootstrap_mock_groups.py yang $catalogue 27 &
    python bootstrap_mock_groups.py yang $catalogue 28 &
    python bootstrap_mock_groups.py yang $catalogue 29 &
    python bootstrap_mock_groups.py yang $catalogue 30

    python bootstrap_mock_groups.py yang $catalogue 31 &
    python bootstrap_mock_groups.py yang $catalogue 32 &
    python bootstrap_mock_groups.py yang $catalogue 33 &
    python bootstrap_mock_groups.py yang $catalogue 34 &
    python bootstrap_mock_groups.py yang $catalogue 35 &
    python bootstrap_mock_groups.py yang $catalogue 36 &
    python bootstrap_mock_groups.py yang $catalogue 37 &
    python bootstrap_mock_groups.py yang $catalogue 38 &
    python bootstrap_mock_groups.py yang $catalogue 39 &
    python bootstrap_mock_groups.py yang $catalogue 40 

    python bootstrap_mock_groups.py yang $catalogue 41 &
    python bootstrap_mock_groups.py yang $catalogue 42 &
    python bootstrap_mock_groups.py yang $catalogue 43 &
    python bootstrap_mock_groups.py yang $catalogue 44 &
    python bootstrap_mock_groups.py yang $catalogue 45 &
    python bootstrap_mock_groups.py yang $catalogue 46 &
    python bootstrap_mock_groups.py yang $catalogue 47 &
    python bootstrap_mock_groups.py yang $catalogue 48 &
    python bootstrap_mock_groups.py yang $catalogue 49 &
    python bootstrap_mock_groups.py yang $catalogue 50 
done

