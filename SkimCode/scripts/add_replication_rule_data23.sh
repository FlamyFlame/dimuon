# Commands to add replication rules for the data23_hi datsets in the physics_HardProbe stream
# Divided the datasets into 4 batches of < 200TB size each
# Divided each batch into 3-4 sub-batches of < 50TB size each
# Considerations: "rucio list-account-usage USERNAME" shows I have 200TB quota in total, and up to 50TB on each scratch-disk based RSE
# Choice of RSE is based on "freespace" (in unit of TB) outputted by "rucio list-rse-attributes RSE_NAME"
# Required process for each batch: 
# 		(1) staging (add replication rules)
# 		(2) data skimming (must be done quickly before file removal from scratch-disk cleaning)
# 		(3) delete replication rules to free space

rucio add-rule data23_hi:data23_hi.00462549.physics_HardProbes.merge.AOD.f1399_m2209 1 CERN-PROD_SCRATCHDISK 
rucio add-rule data23_hi:data23_hi.00462809.physics_HardProbes.merge.AOD.f1406_m2212 1 CERN-PROD_SCRATCHDISK 
rucio add-rule data23_hi:data23_hi.00462759.physics_HardProbes.merge.AOD.f1406_m2212 1 CERN-PROD_SCRATCHDISK 
rucio add-rule data23_hi:data23_hi.00462201.physics_HardProbes.merge.AOD.f1399_m2209 1 CERN-PROD_SCRATCHDISK 

rucio add-rule data23_hi:data23_hi.00463046.physics_HardProbes.merge.AOD.f1408_m2212 1 DESY-HH_SCRATCHDISK 
rucio add-rule data23_hi:data23_hi.00463308.physics_HardProbes.merge.AOD.f1408_m2212 1 DESY-HH_SCRATCHDISK 
rucio add-rule data23_hi:data23_hi.00462107.physics_HardProbes.merge.AOD.f1396_m2209 1 DESY-HH_SCRATCHDISK 

rucio add-rule data23_hi:data23_hi.00461636.physics_HardProbes.merge.AOD.f1393_m2203 1 INFN-T1_SCRATCHDISK  
rucio add-rule data23_hi:data23_hi.00463216.physics_HardProbes.merge.AOD.f1408_m2212 1 INFN-T1_SCRATCHDISK  
rucio add-rule data23_hi:data23_hi.00463364.physics_HardProbes.merge.AOD.f1408_m2212 1 INFN-T1_SCRATCHDISK  

rucio add-rule data23_hi:data23_hi.00461899.physics_HardProbes.merge.AOD.f1393_m2203 1 BNL-OSG2_SCRATCHDISK
rucio add-rule data23_hi:data23_hi.00462588.physics_HardProbes.merge.AOD.f1403_m2209 1 BNL-OSG2_SCRATCHDISK
rucio add-rule data23_hi:data23_hi.00463040.physics_HardProbes.merge.AOD.f1408_m2212 1 BNL-OSG2_SCRATCHDISK

##############################################################################

rucio add-rule data23_hi:data23_hi.00463222.physics_HardProbes.merge.AOD.f1408_m2212 1 CERN-PROD_SCRATCHDISK 
rucio add-rule data23_hi:data23_hi.00463389.physics_HardProbes.merge.AOD.f1408_m2212 1 CERN-PROD_SCRATCHDISK 
# rucio add-rule data23_hi:data23_hi.00462240.physics_HardProbes.merge.AOD.f1399_m2209 1 CERN-PROD_SCRATCHDISK 
rucio add-rule data23_hi:data23_hi.00462717.physics_HardProbes.merge.AOD.f1403_m2212 1 CERN-PROD_SCRATCHDISK 
rucio add-rule data23_hi:data23_hi.00462485.physics_HardProbes.merge.AOD.f1399_m2209 1 CERN-PROD_SCRATCHDISK 
rucio add-rule data23_hi:data23_hi.00462494.physics_HardProbes.merge.AOD.f1399_m2209 1 CERN-PROD_SCRATCHDISK 
rucio add-rule data23_hi:data23_hi.00463017.physics_HardProbes.merge.AOD.f1408_m2212 1 CERN-PROD_SCRATCHDISK 

rucio add-rule data23_hi:data23_hi.00463380.physics_HardProbes.merge.AOD.f1408_m2212 1 DESY-HH_SCRATCHDISK  
rucio add-rule data23_hi:data23_hi.00463315.physics_HardProbes.merge.AOD.f1408_m2212 1 DESY-HH_SCRATCHDISK  
rucio add-rule data23_hi:data23_hi.00462705.physics_HardProbes.merge.AOD.f1403_m2212 1 DESY-HH_SCRATCHDISK  
rucio add-rule data23_hi:data23_hi.00462016.physics_HardProbes.merge.AOD.f1393_m2206 1 DESY-HH_SCRATCHDISK  
rucio add-rule data23_hi:data23_hi.00461655.physics_HardProbes.merge.AOD.f1393_m2203 1 DESY-HH_SCRATCHDISK  

rucio add-rule data23_hi:data23_hi.00463414.physics_HardProbes.merge.AOD.f1408_m2212 1 INFN-T1_SCRATCHDISK  
rucio add-rule data23_hi:data23_hi.00462969.physics_HardProbes.merge.AOD.f1406_m2212 1 INFN-T1_SCRATCHDISK  
rucio add-rule data23_hi:data23_hi.00462022.physics_HardProbes.merge.AOD.f1396_m2206 1 INFN-T1_SCRATCHDISK  

rucio add-rule data23_hi:data23_hi.00461633.physics_HardProbes.merge.AOD.f1393_m2203 1 BNL-OSG2_SCRATCHDISK
rucio add-rule data23_hi:data23_hi.00462576.physics_HardProbes.merge.AOD.f1399_m2209 1 BNL-OSG2_SCRATCHDISK
rucio add-rule data23_hi:data23_hi.00463120.physics_HardProbes.merge.AOD.f1408_m2212 1 BNL-OSG2_SCRATCHDISK
rucio add-rule data23_hi:data23_hi.00462111.physics_HardProbes.merge.AOD.f1399_m2209 1 BNL-OSG2_SCRATCHDISK

##############################################################################


rucio add-rule data23_hi:data23_hi.00463185.physics_HardProbes.merge.AOD.f1408_m2212 1 CERN-PROD_SCRATCHDISK 
rucio add-rule data23_hi:data23_hi.00463021.physics_HardProbes.merge.AOD.f1408_m2212 1 CERN-PROD_SCRATCHDISK 
rucio add-rule data23_hi:data23_hi.00463263.physics_HardProbes.merge.AOD.f1408_m2212 1 CERN-PROD_SCRATCHDISK 

rucio add-rule data23_hi:data23_hi.00462763.physics_HardProbes.merge.AOD.f1406_m2212 1 DESY-HH_SCRATCHDISK  
rucio add-rule data23_hi:data23_hi.00462814.physics_HardProbes.merge.AOD.f1406_m2212 1 DESY-HH_SCRATCHDISK  
rucio add-rule data23_hi:data23_hi.00463155.physics_HardProbes.merge.AOD.f1408_m2212 1 DESY-HH_SCRATCHDISK  

rucio add-rule data23_hi:data23_hi.00462964.physics_HardProbes.merge.AOD.f1406_m2212 1 INFN-T1_SCRATCHDISK  
rucio add-rule data23_hi:data23_hi.00462145.physics_HardProbes.merge.AOD.f1399_m2209 1 INFN-T1_SCRATCHDISK  
rucio add-rule data23_hi:data23_hi.00463255.physics_HardProbes.merge.AOD.f1408_m2212 1 INFN-T1_SCRATCHDISK  
rucio add-rule data23_hi:data23_hi.00463189.physics_HardProbes.merge.AOD.f1408_m2212 1 INFN-T1_SCRATCHDISK  
rucio add-rule data23_hi:data23_hi.00462667.physics_HardProbes.merge.AOD.f1403_m2212 1 INFN-T1_SCRATCHDISK  
rucio add-rule data23_hi:data23_hi.00461641.physics_HardProbes.merge.AOD.f1393_m2203 1 INFN-T1_SCRATCHDISK  

rucio add-rule data23_hi:data23_hi.00461669.physics_HardProbes.merge.AOD.f1393_m2203 1 BNL-OSG2_SCRATCHDISK
rucio add-rule data23_hi:data23_hi.00461635.physics_HardProbes.merge.AOD.f1393_m2203 1 BNL-OSG2_SCRATCHDISK
rucio add-rule data23_hi:data23_hi.00462617.physics_HardProbes.merge.AOD.f1403_m2209 1 BNL-OSG2_SCRATCHDISK
rucio add-rule data23_hi:data23_hi.00462777.physics_HardProbes.merge.AOD.f1406_m2212 1 BNL-OSG2_SCRATCHDISK
rucio add-rule data23_hi:data23_hi.00462502.physics_HardProbes.merge.AOD.f1399_m2209 1 BNL-OSG2_SCRATCHDISK
rucio add-rule data23_hi:data23_hi.00462580.physics_HardProbes.merge.AOD.f1403_m2209 1 BNL-OSG2_SCRATCHDISK
rucio add-rule data23_hi:data23_hi.00462149.physics_HardProbes.merge.AOD.f1399_m2209 1 BNL-OSG2_SCRATCHDISK

##############################################################################


rucio add-rule data23_hi:data23_hi.00463427.physics_HardProbes.merge.AOD.f1408_m2212 1 CERN-PROD_SCRATCHDISK 
rucio add-rule data23_hi:data23_hi.00462441.physics_HardProbes.merge.AOD.f1399_m2209 1 CERN-PROD_SCRATCHDISK 
rucio add-rule data23_hi:data23_hi.00462995.physics_HardProbes.merge.AOD.f1408_m2212 1 CERN-PROD_SCRATCHDISK 

rucio add-rule data23_hi:data23_hi.00462533.physics_HardProbes.merge.AOD.f1399_m2209 1 DESY-HH_SCRATCHDISK 
rucio add-rule data23_hi:data23_hi.00462542.physics_HardProbes.merge.AOD.f1399_m2209 1 DESY-HH_SCRATCHDISK 
rucio add-rule data23_hi:data23_hi.00462205.physics_HardProbes.merge.AOD.f1399_m2209 1 DESY-HH_SCRATCHDISK 

rucio add-rule data23_hi:data23_hi.00461674.physics_HardProbes.merge.AOD.f1393_m2203 1 BNL-OSG2_SCRATCHDISK
rucio add-rule data23_hi:data23_hi.00462767.physics_HardProbes.merge.AOD.f1406_m2212 1 BNL-OSG2_SCRATCHDISK
rucio add-rule data23_hi:data23_hi.00461738.physics_HardProbes.merge.AOD.f1393_m2203 1 BNL-OSG2_SCRATCHDISK
rucio add-rule data23_hi:data23_hi.00462244.physics_HardProbes.merge.AOD.f1399_m2209 1 BNL-OSG2_SCRATCHDISK
rucio add-rule data23_hi:data23_hi.00462677.physics_HardProbes.merge.AOD.f1403_m2212 1 BNL-OSG2_SCRATCHDISK
rucio add-rule data23_hi:data23_hi.00463004.physics_HardProbes.merge.AOD.f1408_m2212 1 BNL-OSG2_SCRATCHDISK
rucio add-rule data23_hi:data23_hi.00463124.physics_HardProbes.merge.AOD.f1408_m2212 1 BNL-OSG2_SCRATCHDISK
