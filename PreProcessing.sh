g++ PreProcessing.C \
-o PreProcessing \
-std=c++11 \
-lpythia8 \
-L/Applications/pythia8243/lib/ \
-I/Applications/pythia8243/include/ \
-I/Applications/root-build/include/ \
-L/Applications/root-build/lib/ \
-lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -stdlib=libc++ -lm -ldl

