#2017 pp data
pathena TrigRates.py \
        --inDS data17_5TeV.periodM.physics_Main.PhysCont.AOD.pro23_v02 \
        --outDS user.yuhang.TrigRates.dimuon.pp2017data.May2024.2. \
        --excludeFile=myfile_286411_r9582.root \
        --excludeFile=myfile_286474_r11819.root \
        --excludeFile=clean.sh \
        --excludeFile=vim_backup

#2015 pp data
# pathena TrigRates.py \
#         --inDsTxt InDstxt_PP2015_5TeV.txt \
#         --outDS user.yuhang.TrigRates.dimuon.pp2015data.May2024.3. \
#         --excludeFile=myfile_286411_r9582.root \
#         --excludeFile=myfile_286474_r11819.root \
#         --excludeFile=clean.sh \
#         --excludeFile=vim_backup



##2015 Pb+Pb data
#pathena TrigRates_UPC2015.py \
#        --inDS data15_hi.periodL.physics_UPC.PhysCont.AOD.pro23_v02 \
#        --outDS user.soumya.TrigRates.UPC.PbPb2015. \
#        --excludeFile=myfile.root \
#        --excludeFile=clean.sh \
#        --excludeFile=vim_backup \

