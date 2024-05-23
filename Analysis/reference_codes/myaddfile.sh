#!/bin/bash

    filelist=`ls pytree_*.root`
    I=0
    J=0
    for file in $filelist
    do
        smalllist=$smalllist" "$file
        I=`expr $I + 1`
        I1=`expr $I % 10`
        if [ $I1 -eq 0 ]
        then
            echo "Adding $smalllist to pytree_MONASH13_HardQCD_noMPIISR_100M_$J.root"
            hadd pytree_MONASH13_HardQCD_noMPIISR_100M_$J.root $smalllist
            smalllist=""
            J=`expr $J + 1`
        fi
    done

    if [ -n "$smalllist" ]
    then
        echo "Adding $smalllist to pytree_MONASH13_HardQCD_noMPIISR_100M_$J.root"
        hadd pytree_MONASH13_HardQCD_noMPIISR_100M_$J.root $smalllist
        smalllist=""
        J=`expr $J + 1`
    fi
    smalllist=""

# hadd pytree_MONASH13_HardQCD_noMPIISR_10M.root pytree_MONASH13_HardQCD_noMPIISR_10M_*.root
# mv   Data_LowAndIntermediateMu_eff0_trig2_jet3_jetcut6_corrtype1.root ~/workarea/CorrelationCode/02AnaCorr/01Rootfiles/
# rm   all*.root
