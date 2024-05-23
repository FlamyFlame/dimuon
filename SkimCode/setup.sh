
cd build
rm -rf .* *
acmSetup --sourcedir=../source AtlasProduction,21.0.20.1
acm clean
acm compile

cd build
rm -rf .* *
acmSetup --sourcedir=../source AthAnalysis,21.2.110
acm clean
acm compile


cd build
rm -rf .* *
acmSetup --sourcedir=../source AtlasDerivation,21.0.19.8
acm compile


cd build
rm -rf .* *
acmSetup --sourcedir=../source Athena,21.0.89
acm compile


acmSetup AthAnalysis,24.2.40
