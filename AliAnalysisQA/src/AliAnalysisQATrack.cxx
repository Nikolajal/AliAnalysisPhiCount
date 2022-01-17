#include "../inc/AliAnalysisQATrack.h"
//
void
fCheckAttribute
 ( TString fParticle = "Tracks", TString fAttribute = fQC_S_Momentum )    {
    // Speeding the process and avoid TCanvas popping up
    gROOT->SetBatch(kTRUE);
    
    auto kLogYScale =   fAttribute.Contains("Momentum");
    
    // Making the Momentum Directory
    gROOT->ProcessLine(Form(".! mkdir -p %s/%s/%s",fOutput_Directory.Data(),fOutput_Directory_Track.Data(),fParticle.Data()));
    
    // Open Data File and load List
    TFile  *fInputFile_DT           =   TFile::Open(Form("%s/%s",fInputData_Directory.Data(),fInputData_FileName.Data()));
    TList  *fInputList_DT           =   (TList*)fInputFile_DT->Get(fInputData_TList.Data());
    TFile  *fInputFile_MC           =   TFile::Open(Form("%s/%s",fInputMntC_Directory.Data(),fInputMntC_FileName.Data()));
    TList  *fInputList_MC           =   (TList*)fInputFile_MC->Get(fInputMntC_TList.Data());
    
    // Recover the Momentum Histogram
    TH2F   *fQC_Tracks_Momentum_DT  =   (TH2F*)fInputList_DT->FindObject(Form(fAttribute.Data(),fParticle.Data()));
    TH2F   *fQC_Tracks_Momentum_MC  =   (TH2F*)fInputList_MC->FindObject(Form(fAttribute.Data(),fParticle.Data()));
    
    //>->   Set Graphic Attributes
    //
    //>->->     Data Histogram
    //
    auto    fNormalisat_DT          =   fQC_Tracks_Momentum_DT->GetEntries();
    fQC_Tracks_Momentum_DT          ->  Scale(1./fNormalisat_DT);
    fQC_Tracks_Momentum_DT          ->  SetTitle("");
    fQC_Tracks_Momentum_DT          ->  SetMarkerStyle(42);
    fQC_Tracks_Momentum_DT          ->  SetMarkerColor(38);
    fQC_Tracks_Momentum_DT          ->  SetLineColor(38);
    fQC_Tracks_Momentum_DT          ->  SetFillColorAlpha(38,.33);
    //
    //>->->     Monte Carlo Histogram
    //
    auto    fNormalisat_MC          =   fQC_Tracks_Momentum_MC->GetEntries();
    fQC_Tracks_Momentum_MC          ->  Scale(1./fNormalisat_MC);
    fQC_Tracks_Momentum_MC          ->  SetMarkerStyle(27);
    fQC_Tracks_Momentum_MC          ->  SetMarkerColor(46);
    fQC_Tracks_Momentum_MC          ->  SetLineColor(46);
    fQC_Tracks_Momentum_MC          ->  SetFillColorAlpha(46,.33);
    //
    //>->->     Legend
    //
    TLegend    *fLegend             =   new TLegend(0.1,0.1,0.2,0.2);
    fLegend                         ->  AddEntry(fQC_Tracks_Momentum_DT,"Data","P");
    fLegend                         ->  AddEntry(fQC_Tracks_Momentum_MC,"MC","P");
    //
    //>->->     Latex
    //
    TLatex     *fText               =   new TLatex();
    
    //  Recovering Charged Tracks
    TH1D       *fQC_Tracks_Momentum_DT_Pos  =   fQC_Tracks_Momentum_DT->ProjectionY(Form("DT_POS_%s",fParticle.Data()),2,2);
    TH1D       *fQC_Tracks_Momentum_MC_Pos  =   fQC_Tracks_Momentum_MC->ProjectionY(Form("MC_POS_%s",fParticle.Data()),2,2);
    TH1D       *fQC_Tracks_Momentum_DT_Neg  =   fQC_Tracks_Momentum_DT->ProjectionY(Form("DT_NEG_%s",fParticle.Data()),1,1);
    TH1D       *fQC_Tracks_Momentum_MC_Neg  =   fQC_Tracks_Momentum_MC->ProjectionY(Form("MC_NEG_%s",fParticle.Data()),1,1);
    TH1D       *fQC_Tracks_Momentum_DT_Tot  =   fQC_Tracks_Momentum_DT->ProjectionY(Form("DT_TOT_%s",fParticle.Data()),1,2);
    TH1D       *fQC_Tracks_Momentum_MC_Tot  =   fQC_Tracks_Momentum_MC->ProjectionY(Form("MC_TOT_%s",fParticle.Data()),1,2);
    
    // Drawing histogram output
    TCanvas    *cCheckMomentumPos   =   new TCanvas();
    gPad                            ->  SetLogy(kLogYScale);
    gStyle                          ->  SetOptStat(0);
    fQC_Tracks_Momentum_DT_Pos      ->  Draw("E HIST L");
    fQC_Tracks_Momentum_MC_Pos      ->  Draw("E HIST L  SAME");
    fLegend                         ->  Draw("SAME");
    fText                           ->  DrawLatexNDC(0.63,0.86,"Positive Tracks");
    cCheckMomentumPos               ->  SaveAs(Form("./%s/%s/%s/%s_POS.pdf",fOutput_Directory.Data(),fOutput_Directory_Track.Data(),fParticle.Data(),Form(fAttribute.Data(),fParticle.Data())));
    
    TCanvas    *cCheckMomentumNeg   =   new TCanvas();
    gPad                            ->  SetLogy(kLogYScale);
    gStyle                          ->  SetOptStat(0);
    fQC_Tracks_Momentum_DT_Neg      ->  Draw("E HIST L");
    fQC_Tracks_Momentum_MC_Neg      ->  Draw("E HIST L  SAME");
    fLegend                         ->  Draw("SAME");
    fText                           ->  DrawLatexNDC(0.63,0.86,"Negative Tracks");
    cCheckMomentumNeg               ->  SaveAs(Form("./%s/%s/%s/%s_NEG.pdf",fOutput_Directory.Data(),fOutput_Directory_Track.Data(),fParticle.Data(),Form(fAttribute.Data(),fParticle.Data())));
    
    TCanvas    *cCheckMomentumTot   =   new TCanvas();
    gPad                            ->  SetLogy(kLogYScale);
    gStyle                          ->  SetOptStat(0);
    fQC_Tracks_Momentum_DT_Tot      ->  Draw("E HIST L");
    fQC_Tracks_Momentum_MC_Tot      ->  Draw("E HIST L  SAME");
    fLegend                         ->  Draw("SAME");
    cCheckMomentumTot               ->  SaveAs(Form("./%s/%s/%s/%s_TOT.pdf",fOutput_Directory.Data(),fOutput_Directory_Track.Data(),fParticle.Data(),Form(fAttribute.Data(),fParticle.Data())));
    
    fInputFile_DT                   ->  Close();
    fInputFile_MC                   ->  Close();
    
    // Exiting Batch Mode
    gROOT->SetBatch(kFALSE);
}
//_____________________________________________________________________________
void
fCheckDCA
 ( TString fParticle = "Tracks", TString fAttribute = fQC_S_DCAXY_PT )    {
    // Speeding the process and avoid TCanvas popping up
    gROOT->SetBatch(kTRUE);
    
    // Making the Momentum Directory
    gROOT->ProcessLine(Form(".! mkdir -p %s/%s/%s",fOutput_Directory.Data(),fOutput_Directory_Track.Data(),fParticle.Data()));
    
    // Open Data File and load List
    TFile  *fInputFile_DT           =   TFile::Open(Form("%s/%s",fInputData_Directory.Data(),fInputData_FileName.Data()));
    TList  *fInputList_DT           =   (TList*)fInputFile_DT->Get(fInputData_TList.Data());
    TFile  *fInputFile_MC           =   TFile::Open(Form("%s/%s",fInputMntC_Directory.Data(),fInputMntC_FileName.Data()));
    TList  *fInputList_MC           =   (TList*)fInputFile_MC->Get(fInputMntC_TList.Data());
    
    // Recover the Momentum Histogram
    TH3F   *fQC_Tracks_Momentum_DT  =   (TH3F*)fInputList_DT->FindObject(Form(fAttribute.Data(),fParticle.Data()));
    TH3F   *fQC_Tracks_Momentum_MC  =   (TH3F*)fInputList_MC->FindObject(Form(fAttribute.Data(),fParticle.Data()));
    
    //>->   Set Graphic Attributes
    //
    //>->->     Data Histogram
    //
    auto    fNormalisat_DT          =   fQC_Tracks_Momentum_DT->GetEntries();
    fQC_Tracks_Momentum_DT          ->  Scale(1./fNormalisat_DT);
    fQC_Tracks_Momentum_DT          ->  SetTitle("");
    fQC_Tracks_Momentum_DT          ->  SetMarkerStyle(42);
    fQC_Tracks_Momentum_DT          ->  SetMarkerColor(38);
    fQC_Tracks_Momentum_DT          ->  SetLineColor(38);
    fQC_Tracks_Momentum_DT          ->  SetFillColorAlpha(38,.33);
    //
    //>->->     Monte Carlo Histogram
    //
    auto    fNormalisat_MC          =   fQC_Tracks_Momentum_MC->GetEntries();
    fQC_Tracks_Momentum_MC          ->  Scale(1./fNormalisat_MC);
    fQC_Tracks_Momentum_MC          ->  SetTitle("");
    fQC_Tracks_Momentum_MC          ->  SetMarkerStyle(27);
    fQC_Tracks_Momentum_MC          ->  SetMarkerColor(46);
    fQC_Tracks_Momentum_MC          ->  SetLineColor(46);
    fQC_Tracks_Momentum_MC          ->  SetFillColorAlpha(46,.33);
    //
    //  Recovering Charged Tracks
    TH2D       *fQC_Tracks_Momentum_DT_Tot  =   (TH2D*)fQC_Tracks_Momentum_DT->Project3D("NOF NUF YZ")->Clone();
    fQC_Tracks_Momentum_DT_Tot->SetNameTitle("fQC_Tracks_Momentum_DT_Tot","");
    fQC_Tracks_Momentum_DT->GetXaxis()->SetRangeUser(0,1);
    TH2D       *fQC_Tracks_Momentum_DT_Pos  =   (TH2D*)fQC_Tracks_Momentum_DT->Project3D("NOF NUF YZ")->Clone();
    fQC_Tracks_Momentum_DT_Pos->SetNameTitle("fQC_Tracks_Momentum_DT_Pos","");
    fQC_Tracks_Momentum_DT->GetXaxis()->SetRangeUser(-1,0);
    TH2D       *fQC_Tracks_Momentum_DT_Neg  =   (TH2D*)fQC_Tracks_Momentum_DT->Project3D("NOF NUF YZ")->Clone();
    fQC_Tracks_Momentum_DT_Neg->SetNameTitle("fQC_Tracks_Momentum_DT_Neg","");
    //
    TH2D       *fQC_Tracks_Momentum_MC_Tot  =   (TH2D*)fQC_Tracks_Momentum_MC->Project3D("NOF NUF YZ")->Clone();
    fQC_Tracks_Momentum_MC_Tot->SetNameTitle("fQC_Tracks_Momentum_MC_Tot","");
    fQC_Tracks_Momentum_MC->GetXaxis()->SetRangeUser(0,1);
    TH2D       *fQC_Tracks_Momentum_MC_Pos  =   (TH2D*)fQC_Tracks_Momentum_MC->Project3D("NOF NUF YZ")->Clone();
    fQC_Tracks_Momentum_MC_Pos->SetNameTitle("fQC_Tracks_Momentum_MC_Pos","");
    fQC_Tracks_Momentum_MC->GetXaxis()->SetRangeUser(-1,0);
    TH2D       *fQC_Tracks_Momentum_MC_Neg  =   (TH2D*)fQC_Tracks_Momentum_MC->Project3D("NOF NUF YZ")->Clone();
    fQC_Tracks_Momentum_MC_Neg->SetNameTitle("fQC_Tracks_Momentum_MC_Neg","");
    //
    //>->->     CUT FUNCTION IN PT
    TF1* fCutTest = new TF1("fCutTest","[1]/(x-[0])",-5,5);
    fCutTest->SetParameter(0,0.0182);
    fCutTest->SetParameter(1,0.0350);
    fCutTest->SetParameter(2,1.01);
    fCutTest->SetNpx(1.e5);
    //
    //>->->     Line
    TLine      *lLine               =   new TLine();
    lLine                           ->  SetLineColor(kRed);
    lLine                           ->  SetLineWidth(2);
    //
    //>->->     Legend
    //
    TLegend    *fLegend             =   new TLegend(0.1,0.1,0.2,0.2);
    fLegend                         ->  AddEntry(fQC_Tracks_Momentum_DT,"Data","P");
    fLegend                         ->  AddEntry(fQC_Tracks_Momentum_MC,"MC","P");
    if ( fAttribute.Contains("DCAZ") ) fLegend                         ->  AddEntry(lLine,"Track Cut","L");
    
    auto kMAX   =   0.;
    auto kMIN   =   -kMAX;
    if ( fAttribute.Contains("DCAZ") ) {
        kMAX    =   2.5;
        kMIN    =   -kMAX;
    }   else if ( fAttribute.Contains("DCAXY") ) {
        kMAX    =   .4;
        kMIN    =   -kMAX;
    }
    auto kZCUT  =   2.;
    auto fLowEdge   =   0.1+0.8*(-kZCUT-kMIN)/(kMAX-kMIN);
    auto fHigEdge   =   0.1+0.8*(+kZCUT-kMIN)/(kMAX-kMIN);
    
    // Drawing histogram output
    TCanvas    *cCheckMomentumPos   =   new TCanvas();
    gPad                            ->  SetLogy(kTRUE);
    gPad                            ->  SetLogz(kTRUE);
    gStyle                          ->  SetOptStat(0);
    fQC_Tracks_Momentum_DT_Pos->GetXaxis()->SetRangeUser(kMIN,kMAX);
    fQC_Tracks_Momentum_DT_Pos->GetYaxis()->SetRangeUser(0.15,10.);
    fQC_Tracks_Momentum_DT_Pos      ->  Draw("COLZ");
    if ( fAttribute.Contains("DCAXY") ) fCutTest->Draw("same");
    if ( fAttribute.Contains("DCAZ") ) lLine    ->  DrawLineNDC(fLowEdge,0.1,fLowEdge,0.9);
    if ( fAttribute.Contains("DCAZ") ) lLine    ->  DrawLineNDC(fHigEdge,0.1,fHigEdge,0.9);
    cCheckMomentumPos               ->  SaveAs(Form("./%s/%s/%s/%s_POS.pdf",fOutput_Directory.Data(),fOutput_Directory_Track.Data(),fParticle.Data(),Form(fAttribute.Data(),fParticle.Data())));
    
    TCanvas    *cCheckMomentumNeg   =   new TCanvas();
    gPad                            ->  SetLogy(kTRUE);
    gPad                            ->  SetLogz(kTRUE);
    gStyle                          ->  SetOptStat(0);
    fQC_Tracks_Momentum_DT_Neg->GetXaxis()->SetRangeUser(kMIN,kMAX);
    fQC_Tracks_Momentum_DT_Neg->GetYaxis()->SetRangeUser(0.15,10.);
    fQC_Tracks_Momentum_DT_Neg      ->  Draw("COLZ");
    if ( fAttribute.Contains("DCAXY") ) fCutTest->Draw("same");
    if ( fAttribute.Contains("DCAZ") ) lLine    ->  DrawLineNDC(fLowEdge,0.1,fLowEdge,0.9);
    if ( fAttribute.Contains("DCAZ") ) lLine    ->  DrawLineNDC(fHigEdge,0.1,fHigEdge,0.9);
    cCheckMomentumNeg               ->  SaveAs(Form("./%s/%s/%s/%s_NEG.pdf",fOutput_Directory.Data(),fOutput_Directory_Track.Data(),fParticle.Data(),Form(fAttribute.Data(),fParticle.Data())));
    
    TCanvas    *cCheckMomentumTot   =   new TCanvas();
    gPad                            ->  SetLogy(kTRUE);
    gPad                            ->  SetLogz(kTRUE);
    gStyle                          ->  SetOptStat(0);
    fQC_Tracks_Momentum_DT_Tot->GetXaxis()->SetRangeUser(kMIN,kMAX);
    fQC_Tracks_Momentum_DT_Tot->GetYaxis()->SetRangeUser(0.15,10.);
    fQC_Tracks_Momentum_DT_Tot      ->  Draw("COLZ");
    if ( fAttribute.Contains("DCAXY") ) fCutTest->Draw("same");
    if ( fAttribute.Contains("DCAZ") ) lLine    ->  DrawLineNDC(fLowEdge,0.1,fLowEdge,0.9);
    if ( fAttribute.Contains("DCAZ") ) lLine    ->  DrawLineNDC(fHigEdge,0.1,fHigEdge,0.9);
    cCheckMomentumTot               ->  SaveAs(Form("./%s/%s/%s/%s_TOT.pdf",fOutput_Directory.Data(),fOutput_Directory_Track.Data(),fParticle.Data(),Form(fAttribute.Data(),fParticle.Data())));

    TCanvas    *cCheckFullStackPos  =   new TCanvas();
    gPad                            ->  SetLogy(kTRUE);
    gStyle                          ->  SetOptStat(0);
    fQC_Tracks_Momentum_DT_Neg->GetXaxis()->SetRangeUser(kMIN,kMAX);
    fQC_Tracks_Momentum_DT_Pos->ProjectionX("DT_POS",-1,1000)->  Draw("E HIST L  ");
    fQC_Tracks_Momentum_MC_Pos->ProjectionX("MC_POS",-1,1000)->  Draw("E HIST L  SAME");
    if ( fAttribute.Contains("DCAZ") ) lLine    ->  DrawLineNDC(fLowEdge,0.1,fLowEdge,0.9);
    if ( fAttribute.Contains("DCAZ") ) lLine    ->  DrawLineNDC(fHigEdge,0.1,fHigEdge,0.9);
    fLegend                         ->  Draw("SAME");
    cCheckFullStackPos              ->  SaveAs(Form("./%s/%s/%s/%s_FLL_POS.pdf",fOutput_Directory.Data(),fOutput_Directory_Track.Data(),fParticle.Data(),Form(fAttribute.Data(),fParticle.Data())));
    
    TCanvas    *cCheckFullStackNeg  =   new TCanvas();
    gPad                            ->  SetLogy(kTRUE);
    gStyle                          ->  SetOptStat(0);
    fQC_Tracks_Momentum_DT_Neg->GetXaxis()->SetRangeUser(kMIN,kMAX);
    fQC_Tracks_Momentum_DT_Neg->ProjectionX("DT_POS",-1,1000)->  Draw("E HIST L  ");
    fQC_Tracks_Momentum_MC_Neg->ProjectionX("MC_POS",-1,1000)->  Draw("E HIST L  SAME");
    if ( fAttribute.Contains("DCAZ") ) lLine    ->  DrawLineNDC(fLowEdge,0.1,fLowEdge,0.9);
    if ( fAttribute.Contains("DCAZ") ) lLine    ->  DrawLineNDC(fHigEdge,0.1,fHigEdge,0.9);
    fLegend                         ->  Draw("SAME");
    cCheckFullStackNeg              ->  SaveAs(Form("./%s/%s/%s/%s_FLL_NEG.pdf",fOutput_Directory.Data(),fOutput_Directory_Track.Data(),fParticle.Data(),Form(fAttribute.Data(),fParticle.Data())));
    
    TCanvas    *cCheckFullStackTot  =   new TCanvas();
    gPad                            ->  SetLogy(kTRUE);
    gStyle                          ->  SetOptStat(0);
    fQC_Tracks_Momentum_DT_Tot->GetXaxis()->SetRangeUser(kMIN,kMAX);
    fQC_Tracks_Momentum_DT_Tot->ProjectionX("DT_POS",-1,1000)->  Draw("E HIST L  ");
    fQC_Tracks_Momentum_MC_Tot->ProjectionX("MC_POS",-1,1000)->  Draw("E HIST L  SAME");
    if ( fAttribute.Contains("DCAZ") ) lLine    ->  DrawLineNDC(fLowEdge,0.1,fLowEdge,0.9);
    if ( fAttribute.Contains("DCAZ") ) lLine    ->  DrawLineNDC(fHigEdge,0.1,fHigEdge,0.9);
    fLegend                         ->  Draw("SAME");
    cCheckFullStackTot              ->  SaveAs(Form("./%s/%s/%s/%s_FLL_TOT.pdf",fOutput_Directory.Data(),fOutput_Directory_Track.Data(),fParticle.Data(),Form(fAttribute.Data(),fParticle.Data())));
    
    
    fInputFile_DT                   ->  Close();
    fInputFile_MC                   ->  Close();
    
    // Exiting Batch Mode
    gROOT->SetBatch(kFALSE);
}
//_____________________________________________________________________________
void
fCheckFllTrack
 ( TString fParticle = "Tracks" )    {
    fCheckAttribute(fParticle,fQC_S_Momentum);
    fCheckAttribute(fParticle,fQC_S_TMomentum);
    fCheckAttribute(fParticle,fQC_S_Eta);
    fCheckAttribute(fParticle,fQC_S_Phi);
    if ( fParticle.Contains("Kaons") ) return;
    fCheckDCA(fParticle,fQC_S_DCAXY_PT);
    fCheckDCA(fParticle,fQC_S_DCAZ_PT);
}
//_____________________________________________________________________________
