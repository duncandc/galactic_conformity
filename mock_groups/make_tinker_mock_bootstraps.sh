#!/bin/bash

catalogues=('Mr19_age_distribution_matching_mock_groups' \
            'Mr19_age_distribution_matching_mock_sys_empty_shuffle_satrel_shuffle_groups')

for catalogue in ${catalogues[*]}
do
    wait
    echo $catalogue

    python bootstrap_mock_groups.py tinker $catalogue 0 & 
    python bootstrap_mock_groups.py tinker $catalogue 1 &
    python bootstrap_mock_groups.py tinker $catalogue 2 &
    python bootstrap_mock_groups.py tinker $catalogue 3 &
    python bootstrap_mock_groups.py tinker $catalogue 4 &
    python bootstrap_mock_groups.py tinker $catalogue 5 &
    python bootstrap_mock_groups.py tinker $catalogue 6 &
    python bootstrap_mock_groups.py tinker $catalogue 7 &
    python bootstrap_mock_groups.py tinker $catalogue 8 &
    python bootstrap_mock_groups.py tinker $catalogue 9 &
    python bootstrap_mock_groups.py tinker $catalogue 10
 
    python bootstrap_mock_groups.py tinker $catalogue 11 &
    python bootstrap_mock_groups.py tinker $catalogue 12 &
    python bootstrap_mock_groups.py tinker $catalogue 13 &
    python bootstrap_mock_groups.py tinker $catalogue 14 &
    python bootstrap_mock_groups.py tinker $catalogue 15 &
    python bootstrap_mock_groups.py tinker $catalogue 16 &
    python bootstrap_mock_groups.py tinker $catalogue 17 &
    python bootstrap_mock_groups.py tinker $catalogue 18 &
    python bootstrap_mock_groups.py tinker $catalogue 19 &
    python bootstrap_mock_groups.py tinker $catalogue 20 

    python bootstrap_mock_groups.py tinker $catalogue 21 &
    python bootstrap_mock_groups.py tinker $catalogue 22 &
    python bootstrap_mock_groups.py tinker $catalogue 23 &
    python bootstrap_mock_groups.py tinker $catalogue 24 &
    python bootstrap_mock_groups.py tinker $catalogue 25 &
    python bootstrap_mock_groups.py tinker $catalogue 26 &
    python bootstrap_mock_groups.py tinker $catalogue 27 &
    python bootstrap_mock_groups.py tinker $catalogue 28 &
    python bootstrap_mock_groups.py tinker $catalogue 29 &
    python bootstrap_mock_groups.py tinker $catalogue 30

    python bootstrap_mock_groups.py tinker $catalogue 31 &
    python bootstrap_mock_groups.py tinker $catalogue 32 &
    python bootstrap_mock_groups.py tinker $catalogue 33 &
    python bootstrap_mock_groups.py tinker $catalogue 34 &
    python bootstrap_mock_groups.py tinker $catalogue 35 &
    python bootstrap_mock_groups.py tinker $catalogue 36 &
    python bootstrap_mock_groups.py tinker $catalogue 37 &
    python bootstrap_mock_groups.py tinker $catalogue 38 &
    python bootstrap_mock_groups.py tinker $catalogue 39 &
    python bootstrap_mock_groups.py tinker $catalogue 40 

    python bootstrap_mock_groups.py tinker $catalogue 41 &
    python bootstrap_mock_groups.py tinker $catalogue 42 &
    python bootstrap_mock_groups.py tinker $catalogue 43 &
    python bootstrap_mock_groups.py tinker $catalogue 44 &
    python bootstrap_mock_groups.py tinker $catalogue 45 &
    python bootstrap_mock_groups.py tinker $catalogue 46 &
    python bootstrap_mock_groups.py tinker $catalogue 47 &
    python bootstrap_mock_groups.py tinker $catalogue 48 &
    python bootstrap_mock_groups.py tinker $catalogue 49 &
    python bootstrap_mock_groups.py tinker $catalogue 50 
done
