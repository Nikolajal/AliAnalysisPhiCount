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
    insFile_DT_Yield[0]                 =   new TFile   ("./result_MCG_PhiQuickReference/MCG_MODE0.root");
    insFile_DT_Yield[1]                 =   new TFile   ("./result_MCG_PhiQuickReference/MCG_MODE1.root");
    insFile_DT_Yield[2]                 =   new TFile   ("./result_MCG_PhiQuickReference/MCG_MODE2.root");
    insFile_DT_Yield[3]                 =   new TFile   ("./result_MCG_PhiQuickReference/MCG_MODE3.root");
    insFile_DT_Yield[4]                 =   new TFile   ("./result_MCG_PhiQuickReference/MCG_MODE4.root");
    insFile_DT_Yield[5]                 =   new TFile   ("./result_MCG_PhiQuickReference/MCG_MODE5.root");
    insFile_DT_Yield[6]                 =   new TFile   ("./result_MCG_PhiQuickReference/MCG_MODE6.root");
    insFile_DT_Yield[7]                 =   new TFile   ("./result_MCG_PhiQuickReference/MCG_MODE7.root");
    insFile_DT_Yield[8]                 =   new TFile   ("./result_MCG_PhiQuickReference/MCG_MODE8.root");
    
    // Recovering the histograms-------------------------------------------------------------------------------

    // >-> GAMMA ANALYSIS //

    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1D      **hPhiGamma           =   new TH1D   *[kOption];
    TH2D      **hPhiGammaMult       =   new TH2D   *[kOption];
    TH1D      **hPhiYield           =   new TH1D   *[kOption];
    //
    //  Defining cumulative histogram over measurable pT
    //
    for ( Int_t iHisto = 0; iHisto < kOption; iHisto++ )
    {
        hName                   =   "hPhiGamma";
        hPhiGamma[iHisto]       =   (TH1D*)(insFile_DT_Yield[iHisto]->Get(hName));
        hName                   =   "hPhiGammaMult";
        hPhiGammaMult[iHisto]   =   (TH2D*)(insFile_DT_Yield[iHisto]->Get(hName));
        hName                   =   "hPhiYield";
        hPhiYield[iHisto]       =   (TH1D*)(insFile_DT_Yield[iHisto]->Get(hName));
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
    
    TString   fLabels[9]  = {"Monash 2013","The MPI-based original Pythia 8 scheme","The new more QCD based scheme","The new gluon-move model","The SK I e^+ e^- CR model","The SK II e^+ e^- CR model","Mode 0 from paper","Mode 2 from paper","Mode 3 from paper"};
    
    TLegend    *fLeg    =   new     TLegend(0.1,0.5,0.3,0.9);
    TLegend    *fLeg2   =   new     TLegend(0.1,0.5,0.3,0.9);
    
    for ( Int_t iBin = 0; iBin < kOption; iBin++ )
    {
        hPhiGammaMult[iBin]->SetFillColorAlpha(0,0);
        hPhiGammaMult[iBin]->SetLineColorAlpha(iBin+1,1.);
        hPhiGammaMult[iBin]->SetLineWidth(3.);
        hPhiGammaMult[iBin]->SetLineStyle(iBin+1);
        hPhiGammaMult[iBin]->Scale(1./10);
        hPhiGammaMult[iBin]->SetMarkerSize(1.);
        hPhiGammaMult[iBin]->SetMarkerColor(iBin+1);
        hPhiGammaMult[iBin]->SetMarkerStyle(iBin+39);
        
        hPhiGamma[iBin]->SetFillColorAlpha(0,0);
        hPhiGamma[iBin]->SetLineColorAlpha(iBin+1,1.);
        hPhiGamma[iBin]->SetLineWidth(3.);
        hPhiGamma[iBin]->SetLineStyle(iBin+1);
        hPhiGamma[iBin]->Scale(1./10);
        hPhiGamma[iBin]->SetMarkerSize(1.);
        hPhiGamma[iBin]->SetMarkerColor(iBin+1);
        hPhiGamma[iBin]->SetMarkerStyle(iBin+39);
        
        hPhiYield[iBin]->SetFillColorAlpha(0,0);
        hPhiYield[iBin]->SetLineColorAlpha(iBin+1,1.);
        hPhiYield[iBin]->SetLineWidth(3.);
        hPhiYield[iBin]->SetLineStyle(iBin+1);
        hPhiYield[iBin]->Scale(1./10);
        hPhiYield[iBin]->SetMarkerSize(1.);
        hPhiYield[iBin]->SetMarkerColor(iBin+1);
        hPhiYield[iBin]->SetMarkerStyle(iBin+39);
        
        fLeg->AddEntry(hPhiGamma[iBin],fLabels[iBin].Data(),"L");
        fLeg2->AddEntry(hPhiGamma[iBin],fLabels[iBin].Data(),"P");
    }
    
    
    
    
    
    
    
    TGraphErrors*hMesGamma = new TGraphErrors("hMesGamma","hMesGamma");
    hMesGamma->SetPoint(0,0.3,0.0533);
    hMesGamma->SetPointError(0,.1,0.0036);
    hMesGamma->SetLineColorAlpha(2,1.);
    hMesGamma->SetMarkerColorAlpha(2,1.);
    hMesGamma->SetMarkerStyle(21);
    
    TGraphErrors*hMesGamm2 = new TGraphErrors("hMesGamm2","hMesGamm2");
    hMesGamm2->SetPoint(0,.3,0.0533);
    hMesGamm2->SetPointError(0,.1,0.027);
    hMesGamm2->SetLineColorAlpha(2,1.);
    hMesGamm2->SetMarkerColorAlpha(2,1.);
    hMesGamm2->SetMarkerStyle(21);
    
    fLeg2->AddEntry(hMesGamma,"Stat Err","E");
    fLeg2->AddEntry(hMesGamm2,"Syst Err","F");
    fLeg2->AddEntry(hMesGamm2,"Data","P");
    
    TCanvas* fCheckGamma = new TCanvas("fCheckGamma","fCheckGamma");
    gStyle->SetOptStat(0);
    hPhiGamma[0]->SetTitle("#gamma_{#phi}");
    hPhiGamma[0]->GetYaxis()->SetTitle("#gamma_{#phi}");
    hPhiGamma[0]->GetYaxis()->SetTitleOffset(1);
    for ( Int_t iBin = 0; iBin < kOption; iBin++ )
    {
        hPhiGamma[iBin]->SetMaximum(0.10);
        hPhiGamma[iBin]->SetMinimum(0.00);
        hPhiGamma[iBin]->Draw("SAME HIST P");
    }
    hMesGamm2->Draw("SAME EP5");
    hMesGamma->Draw("SAME EP");
    fLeg2->Draw("SAME");
    fCheckGamma->SaveAs("fCheckGamma.pdf");
    
    
    
    
    
    
    
    
    TGraphErrors*hMesPhi = new TGraphErrors("hMesPhi","hMesPhi");
    hMesPhi->SetPoint(0,1.3,        0.0330);
    hMesPhi->SetPointError(0,0.1,   0.0002);
    hMesPhi->SetLineColorAlpha(2,1.);
    hMesPhi->SetMarkerColorAlpha(2,1.);
    hMesPhi->SetMarkerStyle(21);
    
    TGraphErrors*hMesPh2 = new TGraphErrors("hMesPh2","hMesPh2");
    hMesPh2->SetPoint(0,1.3,        0.033);
    hMesPh2->SetPointError(0,0.1,   0.0036);
    hMesPh2->SetLineColorAlpha(2,1.);
    hMesPh2->SetMarkerColorAlpha(2,1.);
    hMesPh2->SetMarkerStyle(21);
    
    TCanvas* fCheckPhi = new TCanvas("fCheckPhi","fCheckPhi");
    gStyle->SetOptStat(0);
    hPhiYield[0]->SetTitle("dN_{#phi}/dy");
    hPhiYield[0]->GetYaxis()->SetTitle("dN_{#phi}/dy");
    hPhiYield[0]->GetYaxis()->SetTitleOffset(1);
    for ( Int_t iBin = 0; iBin < kOption; iBin++ )
    {
        hPhiYield[iBin]->SetMaximum(0.06);
        hPhiYield[iBin]->SetMinimum(0.02);
        hPhiYield[iBin]->GetXaxis()->SetRange(2,2);
        hPhiYield[iBin]->Draw("SAME HIST P");
    }
    hMesPh2->Draw("SAME EP5");
    hMesPhi->Draw("SAME EP");
    fLeg2->Draw("SAME");
    fCheckPhi->SaveAs("fCheckPhi.pdf");
    
    
    
    
    
    
    
    TGraphErrors*hMesPhiPhi = new TGraphErrors("hMesPhiPhi","hMesPhiPhi");
    hMesPhiPhi->SetPoint(0,2.3,     0.00144);
    hMesPhiPhi->SetPointError(0,.1, 0.00005);
    hMesPhiPhi->SetLineColorAlpha(2,1.);
    hMesPhiPhi->SetMarkerColorAlpha(2,1.);
    hMesPhiPhi->SetMarkerStyle(21);
    
    TGraphErrors*hMesPhiPh2 = new TGraphErrors("hMesPhiPh2","hMesPhiPh2");
    hMesPhiPh2->SetPoint(0,2.3,     0.00144);
    hMesPhiPh2->SetPointError(0,.1, 0.00030);
    hMesPhiPh2->SetLineColorAlpha(2,1.);
    hMesPhiPh2->SetMarkerColorAlpha(2,1.);
    hMesPhiPh2->SetMarkerStyle(21);
    
    TCanvas* fCheckPhiPhi = new TCanvas("fCheckPhiPhi","fCheckPhiPhi");
    gStyle->SetOptStat(0);
    hPhiYield[0]->SetTitle("dN_{#phi#phi}/dy");
    hPhiYield[0]->GetYaxis()->SetTitle("dN_{#phi#phi}/dy");
    hPhiYield[0]->GetYaxis()->SetTitleOffset(1);
    for ( Int_t iBin = 0; iBin < kOption; iBin++ )
    {
        hPhiYield[iBin]->SetMaximum(0.004);
        hPhiYield[iBin]->SetMinimum(0.0005);
        hPhiYield[iBin]->GetXaxis()->SetRange(3,3);
        hPhiYield[iBin]->Draw("SAME HIST P");
    }
    hMesPhiPh2->Draw("SAME EP5");
    hMesPhiPhi->Draw("SAME EP");
    fLeg2->Draw("SAME");
    fCheckPhiPhi->SaveAs("fCheckPhiPhi.pdf");
    
    /*
    for ( Int_t iBin = 0; iBin < kOption; iBin++ )
    {
        hPhiGammaOverlap->SetBinContent (iBin+1,hPhiGamma[iBin]->GetBinContent(1));
        hPhiGammaOverlap->SetBinError   (iBin+1,hPhiGamma[iBin]->GetBinError  (1));
    }
                 
    TCanvas    *fDraw   =   new     TCanvas();
    
    hPhiGammaMult[0]    ->  GetXaxis()  ->  SetTitle("N_{ch} (#eta < 0.5)");
    hPhiGammaMult[0]    ->  GetYaxis()  ->  SetTitle("#gamma_{#phi}");
    
    
    for ( Int_t iBin = 0; iBin < kOption; iBin++ )
    {
        hPhiGammaMult[iBin]->SetFillColorAlpha(0,0);
        hPhiGammaMult[iBin]->SetLineColorAlpha(iBin+1,1.);
        hPhiGammaMult[iBin]->SetLineWidth(3.);
        hPhiGammaMult[iBin]->SetLineStyle(iBin+1);
        hPhiGammaMult[iBin]->SetMarkerSize(0.);
        //hPhiGammaMult[iBin]    ->  SetMinimum(1.e-2);
        //hPhiGammaMult[iBin]    ->  SetMaximum(1.e10);
        
    }
    for ( Int_t iBin = 0; iBin < kOption; iBin++ )
    {
        //hPhiGammaMult[iBin]->Draw("SAME E3");
        hPhiGammaMult[iBin]->Draw("SAME HIST L");
        fLeg->AddEntry(hPhiGammaMult[iBin],fLabels[iBin].Data(),"L");
        hPhiGammaMult[iBin]->Write();
    }
    fLeg->Draw("SAME");
    
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    
    fDraw->Write();
    
    TCanvas    *fDraw2   =   new     TCanvas();
    
    hPhiYield[0]    ->  GetXaxis()  ->  SetTitle("Yield");
    
    for ( Int_t iBin = 0; iBin < kOption; iBin++ )
    {
        hPhiYield[iBin]->SetFillColorAlpha(0,0);
        hPhiYield[iBin]->SetLineColorAlpha(iBin+1,1.);
        hPhiYield[iBin]->SetLineWidth(3.);
        hPhiYield[iBin]->SetLineStyle(iBin+1);
        hPhiYield[iBin]->SetMarkerSize(0.);
        hPhiYield[iBin]    ->  GetXaxis()  ->  SetRangeUser(0.5,5.5);
        hPhiYield[iBin]    ->  SetMinimum(1.e-6);
        hPhiYield[iBin]    ->  SetMaximum(1.);
        
    }
    for ( Int_t iBin = 0; iBin < kOption; iBin++ )
    {
        //hPhiGammaMult[iBin]->Draw("SAME E3");
        hPhiYield[iBin]->Draw("SAME HIST L");
        hPhiGammaMult[iBin]->Write();
    }
    fLeg2->Draw("SAME");
    
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    
    fDraw2->Write();
    
    
    hPhiGammaOverlap->Write();
    */
    outCheckFitYld->Close();
    
    fStopTimer("Fit_for_extrapolation");
}
