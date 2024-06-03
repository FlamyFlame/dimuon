# For Run2 data that was processed in release 21, we need to use release 21 to skim the data
# Once logged into lxplus, do "ssh lxplus7", then "setupATLAS"

#setting up release 21.2
cp CmakeFile/CMakeLists_Rel21.2.txt source_21/HFtrigValidation/CMakeLists.txt
#edit the file source_21/HFtrigValidation/HFtrigValidation/AthenaVersion.h
rm -rf build_21
mkdir build_21
cd build_21
# acmSetup --sourcedir=../source_21 AthAnalysis,21.2,latest
acmSetup --sourcedir=../source_21 AthAnalysis,21.2.200
acm compile



#setting up release 24.2
cp CmakeFile/CMakeLists_Rel24.2.txt source_24/HFtrigValidation/CMakeLists.txt
#edit the file source_24/HFtrigValidation/HFtrigValidation/AthenaVersion.h
rm -rf build_24
mkdir build_24
cd build_24
acmSetup  --sourcedir=../source_24 AthAnalysis,24.2.40
acm compile
