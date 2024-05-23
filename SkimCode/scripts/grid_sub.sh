#2018 Pb+Pb data
pathena TrigRates_UPC2018.py \
        --inDS data18_hi.periodAllYear.physics_UPC.PhysCont.AOD.t0pro23_v01 \
        --outDS user.soumya.TrigRates.UPC.PbPb2018.Feb2024.2. \
        --excludeFile=myfile.root \
        --excludeFile=clean.sh \
        --excludeFile=vim_backup

##2015 Pb+Pb data
#pathena TrigRates_UPC2015.py \
#        --inDS data15_hi.periodL.physics_UPC.PhysCont.AOD.pro23_v02 \
#        --outDS user.soumya.TrigRates.UPC.PbPb2015. \
#        --excludeFile=myfile.root \
#        --excludeFile=clean.sh \
#        --excludeFile=vim_backup \

