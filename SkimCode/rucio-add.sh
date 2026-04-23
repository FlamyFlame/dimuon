  RSE=RAL-LCG2-ECHO_SCRATCHDISK                                                                                                                           
  for ds in \
    data23_hi:data23_hi.00463189.physics_HardProbes.merge.AOD.r16069_p6447_tid41716200_00 \
    data23_hi:data23_hi.00463308.physics_HardProbes.merge.AOD.r16069_p6447_tid41716227_00 \
    data23_hi:data23_hi.00463380.physics_HardProbes.merge.AOD.r16069_p6447_tid41716239_00 ; do                                                            
    rucio add-rule --lifetime 2592000 "$ds" 1 "$RSE"                                                                                                      
  done   
