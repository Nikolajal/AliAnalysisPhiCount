// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../inc/AliAnalysisPhiPair.h"
#include "RooMsgService.h"

void Analysis_PhiQuickReference ( bool fSilent = false )
{
    //---------------------//
    //  Setting up input   //
    //---------------------//
    
    //-// OPTIONS
    
    // Silencing warnings for smoother
    if ( fSilent )
    {
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
        RooMsgService::instance().setSilentMode(fSilent);
    }
    
    // Retrieving PreProcessed data histograms
    Int_t kOption = 9;
    TFile **insFile_DT_Yield            =   new TFile  *[kOption];
    insFile_DT_Yield[0]                 =   new TFile   ("./result_MCG_PhiQuickReference/MCG_PhiQuickReference_000000000.root");
    insFile_DT_Yield[1]                 =   new TFile   ("./result_MCG_PhiQuickReference/MCG_PhiQuickReference_1000000000.root");
    insFile_DT_Yield[2]                 =   new TFile   ("./result_MCG_PhiQuickReference/MCG_PhiQuickReference_2000000000.root");
    insFile_DT_Yield[3]                 =   new TFile   ("./result_MCG_PhiQuickReference/MCG_PhiQuickReference_3000000000.root");
    insFile_DT_Yield[4]                 =   new TFile   ("./result_MCG_PhiQuickReference/MCG_PhiQuickReference_4000000000.root");
    insFile_DT_Yield[5]                 =   new TFile   ("./result_MCG_PhiQuickReference/MCG_PhiQuickReference_5000000000.root");
    insFile_DT_Yield[6]                 =   new TFile   ("./result_MCG_PhiQuickReference/MCG_PhiQuickReference_6000000000.root");
    insFile_DT_Yield[7]                 =   new TFile   ("./result_MCG_PhiQuickReference/MCG_PhiQuickReference_7000000000.root");
    insFile_DT_Yield[8]                 =   new TFile   ("./result_MCG_PhiQuickReference/MCG_PhiQuickReference_000000000.root");
    
    // Recovering the histograms-------------------------------------------------------------------------------

    // >-> GAMMA ANALYSIS //

    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1D      **hPhiGamma           =   new TH1D   *[kOption];
    TH2D      **hPhiGammaMult       =   new TH2D   *[kOption];
    //
    //  Defining cumulative histogram over measurable pT
    //
    for ( Int_t iHisto = 0; iHisto < kOption; iHisto++ )
    {
        hName                   =   "hPhiGamma";
        hPhiGamma[iHisto]       =   (TH1D*)(insFile_DT_Yield[iHisto]->Get(hName));
        hName                   =   "hPhiGammaMult";
        hPhiGammaMult[iHisto]   =   (TH2D*)(insFile_DT_Yield[iHisto]->Get(hName));
    }
    //
    
    //---------------------//
    //  Setting up output  //
    //---------------------//
    
    int kMultBin    = 10;
    Double_t * fArrMult = new Double_t [kMultBin+1];
    Float_t fUtility [11] = {0,20,40,60,80,100,120,140,160,180,200};
    for ( int i = 0; i <= kMultBin; i++ )
    {
        fArrMult[i] =   fUtility[i];
    }
    
    // Creating the histograms-------------------------------------------------------------------------------

    // >-> YIELD ANALYSIS //

    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1D       *hPhiGammaOverlap;
    TH1D       *hPhiGammaMultOverlap;
    //
    hName                   =   Form("hPhiGammaOverlap");
    hTitle                  =   Form("hPhiGammaOverlap");
    hPhiGammaOverlap        =   new TH1D (hName,hTitle,kOption,-0.5,kOption-0.5);
    //
    hName                   =   Form("hPhiGammaMultOverlap");
    hTitle                  =   Form("hPhiGammaMultOverlap");
    hPhiGammaMultOverlap    =   new TH1D (hName,hTitle,kMultBin,fArrMult);
    //
    
    //-------------------------//
    //  Filling output objects //
    //-------------------------//
    
    fStartTimer("Fit_for_extrapolation");
    
    // Output File for Fit Check
    TFile*  outCheckFitYld  =   new TFile("Example.root","recreate");
    
    //fPrintLoopTimer("Fit_for_extrapolation",1,1,1);
    
    for ( Int_t iBin = 0; iBin < kOption; iBin++ )
    {
        hPhiGammaOverlap->SetBinContent (iBin+1,hPhiGamma[iBin]->GetBinContent(1));
        hPhiGammaOverlap->SetBinError   (iBin+1,hPhiGamma[iBin]->GetBinError  (1));
    }
                 
    TCanvas    *fDraw   =   new     TCanvas();
    TLegend    *fLeg    =   new     TLegend(0.1,0.5,0.3,0.9);
    
    hPhiGammaMult[0]    ->  GetXaxis()  ->  SetTitle("N_{ch} (#eta < 0.5)");
    hPhiGammaMult[0]    ->  GetYaxis()  ->  SetTitle("#gamma_{#phi}");
    
    //TString   fLabels[9]  = {"Mode "};
    
    for ( Int_t iBin = 0; iBin < kOption; iBin++ )
    {
        hPhiGammaMult[iBin]->SetFillColorAlpha(iBin+1,0.3);
        hPhiGammaMult[iBin]->SetLineColorAlpha(iBin+1,1.);
        hPhiGammaMult[iBin]->SetMarkerSize(0.);
        hPhiGammaMult[iBin]->Scale(pow(10,iBin));
        hPhiGammaMult[iBin]    ->  SetMinimum(1.e-2);
        hPhiGammaMult[iBin]    ->  SetMaximum(1.e10);
        hPhiGammaMult[iBin]->Fit("expo");
        
    }
    for ( Int_t iBin = 0; iBin < kOption; iBin++ )
    {
        hPhiGammaMult[iBin]->Draw("SAME E3");
        hPhiGammaMult[iBin]->Draw("SAME EP");
        fLeg->AddEntry(hPhiGammaMult[iBin],Form("%d",iBin),"CF");
        hPhiGammaMult[iBin]->Write();
    }
    fLeg->Draw("SAME");
    
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    
    fDraw->Write();
    
    hPhiGammaOverlap->Write();
    
    outCheckFitYld->Close();
    
    fStopTimer("Fit_for_extrapolation");
}
