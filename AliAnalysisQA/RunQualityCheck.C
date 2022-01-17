#include "./inc/AliAnalysisQA.h"

void RunQualityCheck( TFile* fInput_MC, TFile* fInput_DT, TString kFolder )  {
    
    // Generate, if necessary, all the folders for the output
    fGenerateOutputDirectories();
    if ( !fSetInputHistograms() )
    
    fCheckFllTrack();
    fCheckFllTrack("Kaons");
    fCheckSelection();
    fCheckSelection("Kaons","TOF");
    fCheckPIDSigma();
    fCheckPIDSigma("Kaons","TOF");
    fCheckVertexPosition();
}
