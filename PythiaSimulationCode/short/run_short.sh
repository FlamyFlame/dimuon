cp ../RunHardQCD_short_allkin.cc .
g++ RunHardQCD_short_allkin.cc -o RunHardQCD_short_allkin -pedantic -W -Wall -Wshadow -fPIC -I$PYTHIA8/include -L$PYTHIA8/lib -lpythia8 -ldl `root-config --libs --cflags` `fastjet-config --cxxflags --libs --plugins`
rm pytreeEff.root 
./RunHardQCD_short_allkin 17 3 1