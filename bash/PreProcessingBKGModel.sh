g++ ./root/PreProcessingBKGModel.C \
-o ./exe/PreProcessingBKGModel \
-std=c++11 \
-I/Applications/root-build/include/ \
-L/Applications/root-build/lib/ \
-lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -stdlib=libc++ -lm -ldl 
