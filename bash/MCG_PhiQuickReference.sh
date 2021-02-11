mkdir -p ./exe

g++ ./src/Utilities/MCG_PhiQuickReference.C \
-o ./exe/MCG_PhiQuickReference \
-std=c++11 \
`root-config --cflags --ldflags --libs --glibs --evelibs` \
`~/Applications/pythia8303/bin/pythia8-config --cflags --libs`
