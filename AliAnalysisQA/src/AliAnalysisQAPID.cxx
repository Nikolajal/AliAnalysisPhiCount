#include "../inc/AliAnalysisQAPID.h"
//
//>->->->->->->->->->->->->->->->
//                              >
//  PID Checks on Selections    >
//                              >
//>->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
//
void
fCheckSelection
 ( TString fParticle = "Kaons", TString fDetector = "TPC", TString fAttribute = fQC_S_S_P )      {
    
    // Speeding the process and avoid TCanvas popping up
    gROOT->SetBatch(kTRUE);
    
    // Making the Momentum Directory
    gROOT->ProcessLine(Form(".! mkdir -p %s/%s/%s",fOutput_Directory.Data(),fOutput_Directory_PID.Data(),fDetector.Data()));
    
    // Open Data File and load List
    TFile  *fInputFile_DT           =   TFile::Open(Form("%s/%s",fInputData_Directory.Data(),fInputData_FileName.Data()));
    TList  *fInputList_DT           =   (TList*)fInputFile_DT->Get(fInputData_TList.Data());
    TFile  *fInputFile_MC           =   TFile::Open(Form("%s/%s",fInputMntC_Directory.Data(),fInputMntC_FileName.Data()));
    TList  *fInputList_MC           =   (TList*)fInputFile_MC->Get(fInputMntC_TList.Data());
    
    // Recover the Sigma Histogram
    TH3F   *fQC_PID_Attribute_DT    =   (TH3F*)fInputList_DT->FindObject(Form(fAttribute.Data(),    fParticle.Data(),   fDetector.Data()));
    TH3F   *fQC_PID_Attribute_MC    =   (TH3F*)fInputList_MC->FindObject(Form(fAttribute.Data(),    fParticle.Data(),   fDetector.Data()));
    TH3F   *fQC_PID_GeneralAt_DT    =   (TH3F*)fInputList_DT->FindObject(Form(fAttribute.Data(),    "Tracks",           fDetector.Data()));
    TH3F   *fQC_PID_GeneralAt_MC    =   (TH3F*)fInputList_MC->FindObject(Form(fAttribute.Data(),    "Tracks",           fDetector.Data()));
    
    //>->   Set Graphic Attributes
    //
    //>->->     Data Histogram
    //
    auto    fNormalisat_DT          =   fQC_PID_Attribute_DT->GetEntries();
    fQC_PID_Attribute_DT            ->  Scale(1./fNormalisat_DT);
    fQC_PID_Attribute_DT            ->  SetTitle("");
    fQC_PID_Attribute_DT            ->  SetMarkerStyle(42);
    fQC_PID_Attribute_DT            ->  SetMarkerColor(kRed);
    fQC_PID_Attribute_DT            ->  SetLineColor(kRed);
    fQC_PID_Attribute_DT            ->  SetFillColorAlpha(kRed,1.);
    //
    //>->->     Monte Carlo Histogram
    //
    auto    fNormalisat_MC          =   fQC_PID_Attribute_MC->GetEntries();
    fQC_PID_Attribute_MC            ->  Scale(1./fNormalisat_MC);
    fQC_PID_Attribute_MC            ->  SetMarkerStyle(27);
    fQC_PID_Attribute_MC            ->  SetMarkerColor(46);
    fQC_PID_Attribute_MC            ->  SetLineColor(46);
    fQC_PID_Attribute_MC            ->  SetFillColorAlpha(46,.33);
    //
    //>->->     Legend
    //
    TLegend    *fLegend             =   new TLegend(0.75,0.25,0.9,0.1);
    fLegend                         ->  AddEntry(fQC_PID_Attribute_DT,"Data","P");
    fLegend                         ->  AddEntry(fQC_PID_Attribute_MC,"MC","P");
    //
    //>->->     Line
    TLine      *lLine               =   new TLine();
    lLine                           ->  SetLineColor(kRed);
    lLine                           ->  SetLineWidth(2);
    //
    //  Recovering Charged Tracks
    TH2D       *fQC_PID_Attribute_DT_Tot  =   (TH2D*)fQC_PID_Attribute_DT->Project3D("NOF NUF ZY")->Clone();
    fQC_PID_Attribute_DT_Tot->SetNameTitle("fQC_PID_Attribute_DT_Tot","");
    fQC_PID_Attribute_DT->GetXaxis()->SetRangeUser(0,1);
    TH2D       *fQC_PID_Attribute_DT_Pos  =   (TH2D*)fQC_PID_Attribute_DT->Project3D("NOF NUF ZY")->Clone();
    fQC_PID_Attribute_DT_Tot->SetNameTitle("fQC_PID_Attribute_DT_Pos","");
    fQC_PID_Attribute_DT->GetXaxis()->SetRangeUser(-1,0);
    TH2D       *fQC_PID_Attribute_DT_Neg  =   (TH2D*)fQC_PID_Attribute_DT->Project3D("NOF NUF ZY")->Clone();
    fQC_PID_Attribute_DT_Tot->SetNameTitle("fQC_PID_Attribute_DT_Neg","");
    //
    TH2D       *fQC_PID_Attribute_MC_Tot  =   (TH2D*)fQC_PID_Attribute_MC->Project3D("NOF NUF ZY")->Clone();
    fQC_PID_Attribute_MC_Tot->SetNameTitle("fQC_PID_Attribute_MC_Tot","");
    fQC_PID_Attribute_MC->GetXaxis()->SetRangeUser(0,1);
    TH2D       *fQC_PID_Attribute_MC_Pos  =   (TH2D*)fQC_PID_Attribute_MC->Project3D("NOF NUF ZY")->Clone();
    fQC_PID_Attribute_MC_Tot->SetNameTitle("fQC_PID_Attribute_MC_Pos","");
    fQC_PID_Attribute_MC->GetXaxis()->SetRangeUser(-1,0);
    TH2D       *fQC_PID_Attribute_MC_Neg  =   (TH2D*)fQC_PID_Attribute_MC->Project3D("NOF NUF ZY")->Clone();
    fQC_PID_Attribute_MC_Tot->SetNameTitle("fQC_PID_Attribute_MC_Neg","");
    
    TH2D       *fQC_PID_GeneralAt_DT_Tot  =   (TH2D*)fQC_PID_GeneralAt_DT->Project3D("NOF NUF ZY")->Clone();
    fQC_PID_GeneralAt_DT_Tot->SetNameTitle("fQC_PID_GeneralAt_DT_Tot","");
    fQC_PID_GeneralAt_DT->GetXaxis()->SetRangeUser(0,1);
    TH2D       *fQC_PID_GeneralAt_DT_Pos  =   (TH2D*)fQC_PID_GeneralAt_DT->Project3D("NOF NUF ZY")->Clone();
    fQC_PID_GeneralAt_DT_Pos->SetNameTitle("fQC_PID_GeneralAt_DT_Pos","");
    fQC_PID_GeneralAt_DT->GetXaxis()->SetRangeUser(-1,0);
    TH2D       *fQC_PID_GeneralAt_DT_Neg  =   (TH2D*)fQC_PID_GeneralAt_DT->Project3D("NOF NUF ZY")->Clone();
    fQC_PID_GeneralAt_DT_Neg->SetNameTitle("fQC_PID_GeneralAt_DT_Neg","");
    //
    TH2D       *fQC_PID_GeneralAt_MC_Tot  =   (TH2D*)fQC_PID_GeneralAt_MC->Project3D("NOF NUF ZY")->Clone();
    fQC_PID_GeneralAt_MC_Tot->SetNameTitle("fQC_PID_GeneralAt_MC_Tot","");
    fQC_PID_GeneralAt_MC->GetXaxis()->SetRangeUser(0,1);
    TH2D       *fQC_PID_GeneralAt_MC_Pos  =   (TH2D*)fQC_PID_GeneralAt_MC->Project3D("NOF NUF ZY")->Clone();
    fQC_PID_GeneralAt_MC_Pos->SetNameTitle("fQC_PID_GeneralAt_MC_Pos","");
    fQC_PID_GeneralAt_MC->GetXaxis()->SetRangeUser(-1,0);
    TH2D       *fQC_PID_GeneralAt_MC_Neg  =   (TH2D*)fQC_PID_GeneralAt_MC->Project3D("NOF NUF ZY")->Clone();
    fQC_PID_GeneralAt_MC_Neg->SetNameTitle("fQC_PID_GeneralAt_MC_Neg","");
    
    // Drawing histogram output
    TCanvas    *cCheckMomentumTot   =   new TCanvas();
    gPad                            ->  SetLogx();
    gPad                            ->  SetLogz();
    gStyle                          ->  SetOptStat(0);
    fQC_PID_GeneralAt_DT_Tot        ->  GetXaxis()->SetRangeUser(0.15,10.);
    fQC_PID_GeneralAt_DT_Tot        ->  Draw("COLZ");
    fQC_PID_Attribute_DT_Tot        ->  Draw("CONT3 SAME");
    cCheckMomentumTot               ->  SaveAs(Form("./%s/%s/%s/%s.pdf",fOutput_Directory.Data(),fOutput_Directory_PID.Data(),fDetector.Data(),Form(fAttribute.Data(),fParticle.Data(),fDetector.Data())));
     
    fInputFile_DT                   ->  Close();
    fInputFile_MC                   ->  Close();
    
    // Exiting Batch Mode
    gROOT->SetBatch(kFALSE);
}
//_____________________________________________________________________________
void
fCheckPIDSigma
 ( TString fParticle = "Kaons", TString fDetector = "TPC" )      {
    
    // Speeding the process and avoid TCanvas popping up
    gROOT->SetBatch(kTRUE);
    
    // Making the Momentum Directory
    gROOT->ProcessLine(Form(".! mkdir -p %s/%s/%s",fOutput_Directory.Data(),fOutput_Directory_PID.Data(),fDetector.Data()));
    
    // Open Data File and load List
    TFile  *fInputFile_DT           =   TFile::Open(Form("%s/%s",fInputData_Directory.Data(),fInputData_FileName.Data()));
    TList  *fInputList_DT           =   (TList*)fInputFile_DT->Get(fInputData_TList.Data());
    TFile  *fInputFile_MC           =   TFile::Open(Form("%s/%s",fInputMntC_Directory.Data(),fInputMntC_FileName.Data()));
    TList  *fInputList_MC           =   (TList*)fInputFile_MC->Get(fInputMntC_TList.Data());
    
    // Recover the Sigma Histogram
    TH3F   *fQC_PID_Attribute_DT    =   (TH3F*)fInputList_DT->FindObject(Form(fQC_PID_S_S_P,  fDetector.Data(),    fParticle.Data()));
    TH3F   *fQC_PID_Attribute_MC    =   (TH3F*)fInputList_MC->FindObject(Form(fQC_PID_S_S_P,  fDetector.Data(),    fParticle.Data()));
    
    //>->   Set Graphic Attributes
    //
    //>->->     Data Histogram
    //
    auto    fNormalisat_DT          =   fQC_PID_Attribute_DT->GetEntries();
    fQC_PID_Attribute_DT            ->  Scale(1./fNormalisat_DT);
    fQC_PID_Attribute_DT            ->  SetTitle("");
    fQC_PID_Attribute_DT            ->  SetMarkerStyle(42);
    fQC_PID_Attribute_DT            ->  SetMarkerColor(38);
    fQC_PID_Attribute_DT            ->  SetLineColor(38);
    fQC_PID_Attribute_DT            ->  SetFillColorAlpha(38,.33);
    //
    //>->->     Monte Carlo Histogram
    //
    auto    fNormalisat_MC          =   fQC_PID_Attribute_MC->GetEntries();
    fQC_PID_Attribute_MC            ->  Scale(1./fNormalisat_MC);
    fQC_PID_Attribute_MC            ->  SetMarkerStyle(27);
    fQC_PID_Attribute_MC            ->  SetMarkerColor(46);
    fQC_PID_Attribute_MC            ->  SetLineColor(46);
    fQC_PID_Attribute_MC            ->  SetFillColorAlpha(46,.33);
    //
    //>->->     Legend
    //
    TLegend    *fLegend             =   new TLegend(0.75,0.25,0.9,0.1);
    fLegend                         ->  AddEntry(fQC_PID_Attribute_DT,"Data","P");
    fLegend                         ->  AddEntry(fQC_PID_Attribute_MC,"MC","P");
    //
    //>->->     Latex
    //
    TLatex     *fText               =   new TLatex();
    //
    //>->->     Line
    TLine      *lLine               =   new TLine();
    lLine                           ->  SetLineColor(kRed);
    lLine                           ->  SetLineWidth(4);
    auto kMAX       =   10.;
    auto kMIN       =   -kMAX;
    auto kZCUT      =   3.;
    auto fLowEdge   =   0.1+0.8*(-kZCUT-kMIN)/(kMAX-kMIN);
    auto fHigEdge   =   0.1+0.8*(+kZCUT-kMIN)/(kMAX-kMIN);
    
    //  Recovering Charged Tracks
    TH2D       *fQC_PID_Attribute_DT_Tot  =   (TH2D*)fQC_PID_Attribute_DT->Project3D("NOF NUF YZ")->Clone();
    fQC_PID_Attribute_DT_Tot->SetNameTitle("fQC_PID_Attribute_DT_Tot","");
    fQC_PID_Attribute_DT->GetXaxis()->SetRangeUser(0,1);
    TH2D       *fQC_PID_Attribute_DT_Pos  =   (TH2D*)fQC_PID_Attribute_DT->Project3D("NOF NUF YZ")->Clone();
    fQC_PID_Attribute_DT_Tot->SetNameTitle("fQC_PID_Attribute_DT_Pos","");
    fQC_PID_Attribute_DT->GetXaxis()->SetRangeUser(-1,0);
    TH2D       *fQC_PID_Attribute_DT_Neg  =   (TH2D*)fQC_PID_Attribute_DT->Project3D("NOF NUF YZ")->Clone();
    fQC_PID_Attribute_DT_Tot->SetNameTitle("fQC_PID_Attribute_DT_Neg","");
    //
    TH2D       *fQC_PID_Attribute_MC_Tot  =   (TH2D*)fQC_PID_Attribute_MC->Project3D("NOF NUF YZ")->Clone();
    fQC_PID_Attribute_MC_Tot->SetNameTitle("fQC_PID_Attribute_MC_Tot","");
    fQC_PID_Attribute_MC->GetXaxis()->SetRangeUser(0,1);
    TH2D       *fQC_PID_Attribute_MC_Pos  =   (TH2D*)fQC_PID_Attribute_MC->Project3D("NOF NUF YZ")->Clone();
    fQC_PID_Attribute_MC_Tot->SetNameTitle("fQC_PID_Attribute_MC_Pos","");
    fQC_PID_Attribute_MC->GetXaxis()->SetRangeUser(-1,0);
    TH2D       *fQC_PID_Attribute_MC_Neg  =   (TH2D*)fQC_PID_Attribute_MC->Project3D("NOF NUF YZ")->Clone();
    fQC_PID_Attribute_MC_Tot->SetNameTitle("fQC_PID_Attribute_MC_Neg","");

    // Drawing histogram output
    TCanvas    *cCheckMomentumTot   =   new TCanvas();
    gPad                            ->  SetLogy();
    gPad                            ->  SetLogz();
    gStyle                          ->  SetOptStat(0);
    fQC_PID_Attribute_DT_Tot        ->  GetYaxis()->SetRangeUser(0.15,10.);
    fQC_PID_Attribute_DT_Tot        ->  Draw("COLZ");
    lLine                           ->  DrawLineNDC(fLowEdge,0.1,fLowEdge,0.9);
    lLine                           ->  DrawLineNDC(fHigEdge,0.1,fHigEdge,0.9);
    
    if ( fDetector.Contains("TPC") )    {
        lLine                           ->  SetLineStyle(5);
        kZCUT      =   5.;
        fLowEdge   =   0.1+0.8*(-kZCUT-kMIN)/(kMAX-kMIN);
        fHigEdge   =   0.1+0.8*(+kZCUT-kMIN)/(kMAX-kMIN);
        
        kMAX       =   10.;
        kMIN       =   0.;
        kZCUT      =   1.5;
        auto fLowEdg2   =   0.1+0.8*(+kZCUT-kMIN)/(kMAX-kMIN);
        auto fHigEdg2   =   0.9;
        lLine                           ->  DrawLineNDC(fLowEdge,fLowEdg2,fLowEdge,fHigEdg2);
        lLine                           ->  DrawLineNDC(fHigEdge,fLowEdg2,fHigEdge,fHigEdg2);
        
        kMAX       =   10;
        kMIN       =   -10.;
        kZCUT      =   7.;
        fLowEdge   =   0.1+0.8*(-kZCUT-kMIN)/(kMAX-kMIN);
        kZCUT      =   5.;
        fHigEdge   =   0.1+0.8*(-kZCUT-kMIN)/(kMAX-kMIN);
        lLine                           ->  DrawLineNDC(fLowEdge,fLowEdg2,fHigEdge,fLowEdg2);
        lLine                           ->  DrawLineNDC(1-fLowEdge,fLowEdg2,1-fHigEdge,fLowEdg2);
        
        lLine                           ->  DrawLineNDC(1-fLowEdge,fLowEdg2,1-fLowEdge,0.1);
        lLine                           ->  DrawLineNDC(fLowEdge,fLowEdg2,fLowEdge,0.1);
    }
    
    cCheckMomentumTot               ->  SaveAs(Form("./%s/%s/%s/%s.pdf",fOutput_Directory.Data(),fOutput_Directory_PID.Data(),fDetector.Data(),Form(fQC_PID_S_S_P,  fDetector.Data(),    fParticle.Data())));
    delete  cCheckMomentumTot;
    
    
    auto    fQC_PID_Attribute_DT_Tot_TPC_CUT    =   fQC_PID_Attribute_DT_Tot        ->  ProjectionX("ProjDT",fQC_PID_Attribute_DT_Tot->GetYaxis()->FindBin(0.3),fQC_PID_Attribute_DT_Tot->GetYaxis()->FindBin(0.5));
    auto    fQC_PID_Attribute_MC_Tot_TPC_CUT    =   fQC_PID_Attribute_MC_Tot        ->  ProjectionX("ProjMC",fQC_PID_Attribute_MC_Tot->GetYaxis()->FindBin(0.3),fQC_PID_Attribute_MC_Tot->GetYaxis()->FindBin(0.5));
    if ( fDetector.Contains("TOF") )    {
        fQC_PID_Attribute_DT_Tot_TPC_CUT    =   fQC_PID_Attribute_DT_Tot        ->  ProjectionX("ProjDT",fQC_PID_Attribute_DT_Tot->GetYaxis()->FindBin(0.75),fQC_PID_Attribute_DT_Tot->GetYaxis()->FindBin(1.00));
        fQC_PID_Attribute_MC_Tot_TPC_CUT    =   fQC_PID_Attribute_MC_Tot        ->  ProjectionX("ProjMC",fQC_PID_Attribute_MC_Tot->GetYaxis()->FindBin(0.75),fQC_PID_Attribute_MC_Tot->GetYaxis()->FindBin(1.00));
    }
            cCheckMomentumTot   =   new TCanvas();
    gStyle                          ->  SetOptStat(0);
    fQC_PID_Attribute_DT_Tot_TPC_CUT->  GetXaxis()->    SetRangeUser(-5,5);
    fQC_PID_Attribute_DT_Tot_TPC_CUT->  Draw("E HIST L");
    fQC_PID_Attribute_MC_Tot_TPC_CUT->  Draw("E HIST L  SAME");
    auto    fDTFull                 =   fQC_PID_Attribute_DT_Tot_TPC_CUT->Integral(fQC_PID_Attribute_DT_Tot_TPC_CUT->GetXaxis()->FindBin(-5),fQC_PID_Attribute_DT_Tot_TPC_CUT->GetXaxis()->FindBin(5));
    auto    fDTPart                 =   fQC_PID_Attribute_DT_Tot_TPC_CUT->Integral(fQC_PID_Attribute_DT_Tot_TPC_CUT->GetXaxis()->FindBin(-3),fQC_PID_Attribute_DT_Tot_TPC_CUT->GetXaxis()->FindBin(3));
    auto    fDTPar2                 =   fQC_PID_Attribute_DT_Tot_TPC_CUT->Integral(fQC_PID_Attribute_DT_Tot_TPC_CUT->GetXaxis()->FindBin(-2),fQC_PID_Attribute_DT_Tot_TPC_CUT->GetXaxis()->FindBin(2));
    auto    fMCFull                 =   fQC_PID_Attribute_MC_Tot_TPC_CUT->Integral(fQC_PID_Attribute_MC_Tot_TPC_CUT->GetXaxis()->FindBin(-5),fQC_PID_Attribute_MC_Tot_TPC_CUT->GetXaxis()->FindBin(5));
    auto    fMCPart                 =   fQC_PID_Attribute_MC_Tot_TPC_CUT->Integral(fQC_PID_Attribute_MC_Tot_TPC_CUT->GetXaxis()->FindBin(-3),fQC_PID_Attribute_MC_Tot_TPC_CUT->GetXaxis()->FindBin(3));
    auto    fMCPar2                 =   fQC_PID_Attribute_MC_Tot_TPC_CUT->Integral(fQC_PID_Attribute_MC_Tot_TPC_CUT->GetXaxis()->FindBin(-2),fQC_PID_Attribute_MC_Tot_TPC_CUT->GetXaxis()->FindBin(2));
    if ( fDetector.Contains("TPC") )    fText   ->  DrawLatexNDC(0.35,0.93,Form("p_{T} [0.30,0.50] GeV/c^{2}"));
    if ( fDetector.Contains("TOF") )    fText   ->  DrawLatexNDC(0.35,0.93,Form("p_{T} [0.75,1.00] GeV/c^{2}"));
    fText                           ->  SetTextColor(38);
    fText                           ->  DrawLatexNDC(0.72,0.80,Form("#frac{N_{3#sigma}}{N_{5#sigma}} = %2.1f%s",100*fDTPart/fDTFull,"%"));
    fText                           ->  DrawLatexNDC(0.72,0.68,Form("#frac{N_{2#sigma}}{N_{5#sigma}} = %2.1f%s",100*fDTPar2/fDTFull,"%"));
    fText                           ->  SetTextColor(46);
    fText                           ->  DrawLatexNDC(0.14,0.80,Form("#frac{N_{3#sigma}}{N_{5#sigma}} = %2.1f%s",100*fMCPart/fMCFull,"%"));
    fText                           ->  DrawLatexNDC(0.14,0.68,Form("#frac{N_{2#sigma}}{N_{5#sigma}} = %2.1f%s",100*fMCPar2/fMCFull,"%"));
    fLegend                         ->  Draw("SAME");
    cCheckMomentumTot               ->  SaveAs(Form("./%s/%s/%s/%s_INT.pdf",fOutput_Directory.Data(),fOutput_Directory_PID.Data(),fDetector.Data(),Form(fQC_PID_S_S_P,  fDetector.Data(),    fParticle.Data())));
    delete  cCheckMomentumTot;
    
            fQC_PID_Attribute_DT_Tot_TPC_CUT    =   fQC_PID_Attribute_DT_Tot        ->  ProjectionX("ProjDT",fQC_PID_Attribute_DT_Tot->GetYaxis()->FindBin(0.0),fQC_PID_Attribute_DT_Tot->GetYaxis()->FindBin(0.28));
            fQC_PID_Attribute_MC_Tot_TPC_CUT    =   fQC_PID_Attribute_MC_Tot        ->  ProjectionX("ProjMC",fQC_PID_Attribute_MC_Tot->GetYaxis()->FindBin(0.0),fQC_PID_Attribute_MC_Tot->GetYaxis()->FindBin(0.28));
            cCheckMomentumTot   =   new TCanvas();
    gStyle                          ->  SetOptStat(0);
    fQC_PID_Attribute_DT_Tot_TPC_CUT->  GetXaxis()->    SetRangeUser(-10,10);
    fQC_PID_Attribute_MC_Tot_TPC_CUT->  Draw("E HIST L  SAME");
    fQC_PID_Attribute_DT_Tot_TPC_CUT->  Draw("E HIST L SAME");
    fLegend                         ->  Draw("SAME");
    cCheckMomentumTot               ->  SaveAs(Form("./%s/%s/%s/%s_CloseUp.pdf",fOutput_Directory.Data(),fOutput_Directory_PID.Data(),fDetector.Data(),Form(fQC_PID_S_S_P,  fDetector.Data(),    fParticle.Data())));
    delete  cCheckMomentumTot;
    
    fInputFile_DT                   ->  Close();
    fInputFile_MC                   ->  Close();
    
    // Exiting Batch Mode
    gROOT->SetBatch(kFALSE);
}
