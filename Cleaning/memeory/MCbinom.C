#include "alidock/AliAnalysisPhiCount/inc/AliAnalysisUtility.h"

void MCbinom () {
    
    Int_t   nMaximumOfStrings   = 100;
    Int_t   nEvents             = 100000;
    Int_t   nSSBreakings        = 0;
    Int_t   nPhi                = 0;
    float fBreakProbability     = 0.12;
    float fEnhancementProb      = 0.15;
    float fEnhancement          = 0.25;
    
    TH2F       *hPhi        =   new TH2F        ("hPhi",        "hPhi",         100,    -0.5, 99.5,   25.,    -0.5, 24.5);
    TProfile   *hPhiMean    =   new TProfile    ("hPhiMean",    "hPhiMean",     100,    -0.5, 99.5,   "s");
    TH1F       *hPhiGamma   =   new TH1F        ("hPhiGamma",   "hPhiGamma",    100,    0., 100.);

    
    fStartTimer("Analysis");
    for ( Int_t nStrings = 0; nStrings < nMaximumOfStrings; nStrings++ )
    {
        fPrintLoopTimer("Analysis",nStrings*nEvents,nMaximumOfStrings*nEvents,nEvents);
        for ( Int_t iEvent = 0; iEvent < nEvents; iEvent++ )
        {
            nPhi = 0;
            for ( Int_t iString = 0; iString < nStrings; iString++ )
            {
                nSSBreakings= 0;
                while ( kTRUE )
                {
                    if ( gRandom->Uniform(0.,1.) >= fEnhancementProb ) fBreakProbability *= 1+fEnhancement;
                    if ( fBreakProbability >= (0.3/2.3)*(1+fEnhancement) ) fBreakProbability = (0.3/2.3)*(1+fEnhancement);
                    if ( gRandom->Uniform(0.,1.) >= fBreakProbability ) break;
                    nSSBreakings++;
                }
                fBreakProbability     = 0.3/2.3;
                if ( nSSBreakings >= 2 )
                {
                    nPhi += nSSBreakings-1;
                }
            }
            hPhi->Fill(nStrings,nPhi);
            hPhiMean->Fill(nStrings,nPhi);
        }
        auto fMean =    hPhiMean->GetBinContent(nStrings+1);
        auto fStdv =    hPhiMean->GetBinError(nStrings+1);
        if ( fMean != 0 )        hPhiGamma->SetBinContent(nStrings+1,fStdv*fStdv/fMean -1.);
    }
    
    fStopTimer("Analysis");
    
    TFile * fout = new TFile ("MCbinom_output.root","recreate");
    hPhi->Write();
    hPhiMean->Write();
    hPhiGamma->Write();
    
    
    hPhi->Draw();
    hPhiMean->Draw("SAME");
    cout << "sigmaphi: " << (pow(hPhi->GetStdDev(),2) / hPhi->GetMean()) -1 << endl;
    
    fout->Close();
}
