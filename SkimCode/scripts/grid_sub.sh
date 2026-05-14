# 2023 Pb+Pb HardProbes data (data23_hi, sqrt(s_NN)=5.36 TeV, repro35 reprocessing)
# Inclusive container: data23_hi:data23_hi.periodAllYear2.physics_HardProbes.PhysCont.AOD.repro35_v01
# (previous single-task attempt with --inDS container + nFilesPerJob 200 failed -- too many jobs)
#
# Total: 60 datasets | 87,924 files | 660 TB | ~397.9M events
# Avg file: ~7.5 GB  | ~4,526 events/file
# --nGBPerJob MAX : lets JEDI use XRootD remote-I/O instead of local staging.
#   Without this, PanDA requires ~490 GB scratch/job to stage 60×7.5 GB files
#   locally -- no site has that. MAX activates remote-read brokerage path, which
#   only requires output+workdir disk (~44 GB), passing SARA-MATRIX (90 GB) and BNL (49 GB).
#   JEDI dynamically picks nFiles/job per site; total jobs/task may differ from
#   the 500-job estimates below, but partition boundaries remain valid.
#
# Run from: SkimCode/run_23hi/
# Setup:  source ../setup_25.sh && lsetup panda
# Full dataset list for bookkeeping: InDstxt_PbPb2023_HP.txt

# Part 1: runs 461633-462705  |  30 datasets | 479 jobs | 127.4M events
pathena --trf "TRIGRATES_RUNMODE=hi2023 athena.py TrigRates_CA.py --filesInput=%IN --evtMax=-1" \
        --inDsTxt InDstxt_PbPb2023_HP_part1.txt \
        --outDS user.yuhang.TrigRates.dimuon.PbPb2023data.May2026.v1.part1. \
        --mergeOutput \
        --nGBPerJob MAX \
        --extFile data23_hi.periodAllYear_DetStatus-v120-pro33-03_MERGED_PHYS_HeavyIon_All_Good_IgnoreBSPOT_INVALID.xml \
        --extOutFile myfile.root \
        --excludeFile=clean.sh \
        --excludeFile=vim_backup

# Part 2: runs 462717-463120  |  16 datasets | 453 jobs | 124.1M events
pathena --trf "TRIGRATES_RUNMODE=hi2023 athena.py TrigRates_CA.py --filesInput=%IN --evtMax=-1" \
        --inDsTxt InDstxt_PbPb2023_HP_part2.txt \
        --outDS user.yuhang.TrigRates.dimuon.PbPb2023data.May2026.v1.part2. \
        --mergeOutput \
        --nGBPerJob MAX \
        --extFile data23_hi.periodAllYear_DetStatus-v120-pro33-03_MERGED_PHYS_HeavyIon_All_Good_IgnoreBSPOT_INVALID.xml \
        --extOutFile myfile.root \
        --excludeFile=clean.sh \
        --excludeFile=vim_backup

# Part 3: runs 463124-463389  |  13 datasets | 451 jobs | 123.6M events
pathena --trf "TRIGRATES_RUNMODE=hi2023 athena.py TrigRates_CA.py --filesInput=%IN --evtMax=-1" \
        --inDsTxt InDstxt_PbPb2023_HP_part3.txt \
        --outDS user.yuhang.TrigRates.dimuon.PbPb2023data.May2026.v1.part3. \
        --excludedSite "*LANCS*" \
        --mergeOutput \
        --nGBPerJob MAX \
        --extFile data23_hi.periodAllYear_DetStatus-v120-pro33-03_MERGED_PHYS_HeavyIon_All_Good_IgnoreBSPOT_INVALID.xml \
        --extOutFile myfile.root \
        --excludeFile=clean.sh \
        --excludeFile=vim_backup

# Part 4: runs 463414-463427  |  2 datasets |  84 jobs |  22.8M events
pathena --trf "TRIGRATES_RUNMODE=hi2023 athena.py TrigRates_CA.py --filesInput=%IN --evtMax=-1" \
        --inDsTxt InDstxt_PbPb2023_HP_part4.txt \
        --outDS user.yuhang.TrigRates.dimuon.PbPb2023data.May2026.v1.part4. \
        --mergeOutput \
        --nGBPerJob MAX \
        --extFile data23_hi.periodAllYear_DetStatus-v120-pro33-03_MERGED_PHYS_HeavyIon_All_Good_IgnoreBSPOT_INVALID.xml \
        --extOutFile myfile.root \
        --excludeFile=clean.sh \
        --excludeFile=vim_backup
