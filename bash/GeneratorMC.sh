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
  
    g++ ./src/Utilities/Anls_MonteCarloGeneratorPhiCount.C \
	-o ./exe/GeneratorMC \
	-std=c++11 \
	-lpythia8 \
	-L/Applications/pythia8303/lib/ \
	-I/Applications/pythia8303/include/ \
    -L/usr/local/Cellar/root/6.22.06_1/lib/root \
    -I/usr/local/Cellar/root/6.22.06_1/include/root \
	-lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame  -lm -ldl -stdlib=libc++
  
fi
