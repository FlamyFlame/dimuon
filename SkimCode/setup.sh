# For Run2 data that was processed in release 21, we need to use release 21 to skim the data
# Once logged into lxplus, do "ssh lxplus7", then "setupATLAS"

#setting up release 21.2
cp CmakeFile/CMakeLists_Rel21.2.txt source/HFtrigValidation/CMakeLists.txt
#edit the file source/HFtrigValidation/HFtrigValidation/AthenaVersion.h
rm -rf build
mkdir build
cd build
acmSetup --sourcedir=../source AthAnalysis,21.2,latest
acm compile



#setting up release 24.2
cp CmakeFile/CMakeLists_Rel24.2.txt source/HFtrigValidation/CMakeLists.txt
#edit the file source/HFtrigValidation/HFtrigValidation/AthenaVersion.h
rm -rf build
mkdir build
cd build
acmSetup AthAnalysis,24.2.40
acm compile
