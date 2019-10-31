g++ -L $ROOTSYS/lib -lRooFit ./root/PreProcessing2D.C \
-o ./exe/PreProcessing2D \
-std=c++11 \
-I/Applications/root-build/include/ \
-L/Applications/root-build/lib/ \
-lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -stdlib=libc++ -lm -ldl 
