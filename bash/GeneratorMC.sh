mkdir -p ./exe

if [ `hostname` == "bownalice07.bo.infn.it" ]; then

    g++ ./root/GeneratorMC.C \
	-o ./exe/GeneratorMC \
	-std=c++11 \
	-lpythia8 \
	-L$PYTHIA_ROOT/lib/ \
	-I$PYTHIA_ROOT/include/ \
	-I$ROOTSYS/include/ \
	-L$ROOTSYS/lib/ \
	-lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -lm -ldl
    
else
  
    g++ ./src/Utilities/MCG_PhiAnalysis.C \
	-o ./exe/MCG_PhiAnalysis \
	-std=c++11 \
    `root-config --cflags --ldflags --libs --glibs --evelibs` \
    `/Applications/pythia8303/bin/pythia8-config --cflags --libs`
  
fi
