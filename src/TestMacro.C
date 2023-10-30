#include "../inc/AliAnalysisPhiPair.h"
//#include "/Users/nikolajal/alice/AliAnalysisPhiCount/AliAnalysisUtility/ReweightEfficiency.C"
#include "./MonteCarlo/MC_Plot.cxx"
#include "./MonteCarlo/MC_PlotGen.cxx"
#include "./MonteCarlo/mergetrees.C"
#include "./Analysis/AN_ShowPlots.cxx"
// !TODO: Can there really be a TODO in a Test Macro?
//!
void
kPlotd
( TGraphErrors* g1D_STD, TGraphErrors* g2D_STD, TGraphErrors* g1D_VAR, TGraphErrors* g2D_VAR, Int_t iTer, Float_t kF1, Float_t kF2, TString kL1, TString kL2 ) {
    //!
    auto k1D_STD = g1D_STD->GetY()[iTer];
    auto e1D_STD = g1D_STD->GetEY()[iTer];
    //!
    auto k1D_VAR = (g1D_VAR->GetY()[2*iTer]*kF1+g1D_VAR->GetY()[2*iTer+1]*kF2)/(kF1+kF2);
    auto e1D_VAR = sqrt((g1D_VAR->GetEY()[2*iTer]*g1D_VAR->GetEY()[2*iTer]*kF1*kF1+g1D_VAR->GetEY()[2*iTer+1]*g1D_VAR->GetEY()[2*iTer+1]*kF2*kF2)/((kF1+kF2)*(kF1+kF2)));
    //!
    auto k2D_STD = g2D_STD->GetY()[iTer];
    auto e2D_STD = g2D_STD->GetEY()[iTer];
    //!
    auto k2D_VAR = (g2D_VAR->GetY()[2*iTer]*kF1+g2D_VAR->GetY()[2*iTer+1]*kF2)/(kF1+kF2);
    auto e2D_VAR = sqrt((g2D_VAR->GetEY()[2*iTer]*g2D_VAR->GetEY()[2*iTer]*kF1*kF1+g2D_VAR->GetEY()[2*iTer+1]*g2D_VAR->GetEY()[2*iTer+1]*kF2*kF2)/((kF1+kF2)*(kF1+kF2)));
    //!
    TCanvas* c1 = new TCanvas();
    gPad->SetLogy();
    gStyle->SetOptStat(0);
    auto hSTD = uPlotDerivedQuantitiesRaw({k1D_STD,e1D_STD,0.},{k2D_STD,e2D_STD,0.});
    auto hVAR = uPlotDerivedQuantitiesRaw({k1D_VAR,e1D_VAR,0.},{k2D_VAR,e2D_VAR,0.});
    TLegend* cL = new TLegend(0.15,0.8,0.4,0.9);
    cL->SetNColumns(2);
    cL->AddEntry(hSTD,kL2.Data());
    cL->AddEntry(hVAR,kL1.Data());
    hSTD->SetLineColor(kRed);
    hVAR->SetLineColor(kBlue);
    hVAR->Draw("");
    hSTD->Draw("SAME");
    hVAR->SetName(kL1.Data());
    hSTD->SetName(kL2.Data());
    cL->Draw("SAME");
    c1->SaveAs(Form("/Users/nrubini/Analysis/ALICE/PWG-LF/PAG-RSN/_1020_Phi_Pair/result/FinalResults/%s.pdf",kL2.Data()));
}
//!
void
TestMacro
()  {
    TFile* cinput = new TFile("/Users/nrubini/Analysis/Rivet/_1020_Phi_Pair/ALICE_2023_IXXXXXX_12_5Mevs_P8_08.root");
    auto hNormalisation = uLoadHistograms<0,TH1D>(cinput,"u08_x01_y02","u08_x01_y02");
    auto hXaxis         = uLoadHistograms<0,TH1D>(cinput,"u07_x01_y02","u07_x01_y02");
    auto h1DYield       = uLoadHistograms<0,TH1D>(cinput,"x02_x02_y02","x02_x02_y02");
    auto h2DYield       = uLoadHistograms<0,TH1D>(cinput,"x03_x02_y03","x03_x02_y03");
    h1DYield->Divide(hNormalisation);
    h2DYield->Divide(hNormalisation);
    hXaxis  ->Divide(hNormalisation);
    //!
    TGraphErrors* g1DYield = new TGraphErrors();
    TGraphErrors* g2DYield = new TGraphErrors();
    TGraphErrors* gR1Yield = new TGraphErrors();
    TGraphErrors* gR2Yield = new TGraphErrors();
    TGraphErrors* gP1Yield = new TGraphErrors();
    TGraphErrors* gP2Yield = new TGraphErrors();
    g1DYield->SetName("g1DYield");
    g2DYield->SetName("g2DYield");
    gR1Yield->SetName("gR1Yield");
    gR2Yield->SetName("gR2Yield");
    gP1Yield->SetName("gP1Yield");
    gP2Yield->SetName("gP2Yield");
    for ( auto iBin = 1; iBin <= h1DYield->GetNbinsX(); iBin++ ) {
        g1DYield->SetPoint      (iBin-1, hXaxis->GetBinContent(iBin),   h1DYield->GetBinContent(iBin)   );
        g1DYield->SetPointError (iBin-1, hXaxis->GetBinError(iBin),     h1DYield->GetBinError(iBin)     );
        g2DYield->SetPoint      (iBin-1, hXaxis->GetBinContent(iBin),   h2DYield->GetBinContent(iBin)   );
        g2DYield->SetPointError (iBin-1, hXaxis->GetBinError(iBin),     h2DYield->GetBinError(iBin)     );
        auto uFinalQuantities = uCalculateDerivedQuantities( {h1DYield->GetBinContent(iBin),h1DYield->GetBinError(iBin),0}, {h2DYield->GetBinContent(iBin),h2DYield->GetBinError(iBin),0} );
        gR1Yield->SetPoint      (iBin-1, hXaxis->GetBinContent(iBin),   get<0>(uFinalQuantities["R1"])  );
        gR1Yield->SetPointError (iBin-1, hXaxis->GetBinError(iBin),     get<1>(uFinalQuantities["R1"])  );
        gR2Yield->SetPoint      (iBin-1, hXaxis->GetBinContent(iBin),   get<0>(uFinalQuantities["R2"])  );
        gR2Yield->SetPointError (iBin-1, hXaxis->GetBinError(iBin),     get<1>(uFinalQuantities["R2"])  );
        gP1Yield->SetPoint      (iBin-1, hXaxis->GetBinContent(iBin),   get<0>(uFinalQuantities["P1"])  );
        gP1Yield->SetPointError (iBin-1, hXaxis->GetBinError(iBin),     get<1>(uFinalQuantities["P1"])  );
        gP2Yield->SetPoint      (iBin-1, hXaxis->GetBinContent(iBin),   get<0>(uFinalQuantities["P2"])  );
        gP2Yield->SetPointError (iBin-1, hXaxis->GetBinError(iBin),     get<1>(uFinalQuantities["P2"])  );
    }
    //!
    TFile* fOutput = new TFile("Pythia8MonashRopes.root","RECREATE");
    g1DYield->Write();
    g2DYield->Write();
    gR1Yield->Write();
    gR2Yield->Write();
    gP1Yield->Write();
    gP2Yield->Write();
    fOutput->Close();
}

