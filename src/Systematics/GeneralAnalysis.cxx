// File for 1-Dimensional Analysis:
// !TODO: Set Global Counter to avoid memory loss due to overlap in names of th**
// !TODO: Make copies for TH1F/TH1D etc.
#include "../../inc/AliAnalysisPhiPair.h"
#include "RooMsgService.h"
//
void
GeneralAnalysis
( TH1F* hStandard, std::vector<TH1F*> hVariations, TString fFolder ) {
    //
    // --------- FIND THE SYSTEMATICAL RELEVANT VARIATIONS
    //
    SetStyle();
    gROOT   ->  ProcessLine ( Form(".! mkdir -p %s",(fFolder+TString("/plots/BarlowCheck/1D/")).Data()) );
    gROOT   ->  ProcessLine ( Form(".! mkdir -p %s",(fFolder+TString("/plots/BinByBinCheck/1D/")).Data()) );
    auto fRelevantVariations    =  uIsRelevantVariation(hStandard,hVariations,(fFolder+TString("/plots/BarlowCheck/1D/")).Data(),"1D");
    //
    // --------- USE THE SYSTEMATICAL RELEVANT VARIATIONS TO DETERMINE THE UNCERTAINTY
    //
    //for ( [...] ) Loop to exclude non Barlow variations
    //
    auto    uStackSystematic    =   uBuildSystematicStack(hStandard,hVariations,fRelevantVariations);
    //
    gROOT->SetBatch();
    TLegend    *lLegend =   new TLegend(0.2,0.85,0.4,0.7);
    TCanvas*c1 = new TCanvas();
    gStyle->SetOptStat(0);
    gPad->SetLogx();
    uStackSystematic->SetMinimum(0);
    uStackSystematic->SetMaximum(uStackSystematic->GetMaximum()*1.3);
    uStackSystematic->Draw("");
    lLegend->AddEntry((TH1F*)uStackSystematic->GetHists()->At(0),"#mu contr.","F");
    lLegend->AddEntry((TH1F*)uStackSystematic->GetHists()->At(1),"#sigma contr.","F");
    lLegend->Draw();
    c1->SaveAs((fFolder+TString("/plots/Full_1D_Sys.pdf")).Data());
    gROOT->SetBatch(kFALSE);
    //
    TFile      *fOutput =   new TFile   (Form("%s%s",fFolder.Data(),"/1D_Systematic.root"),"recreate");
    uBuildSystematicError(hStandard,hVariations,(fFolder+TString("/plots/BinByBinCheck/1D/")).Data(),"1D",fRelevantVariations)->Write();
    fOutput->Close();
    //
}
//
void
GeneralAnalysis
( TH2F* hStandard, std::vector<TH2F*> hVariations, TString fFolder ) {
    //
    // --------- FIND THE SYSTEMATICAL RELEVANT VARIATIONS
    //
    SetStyle();
    gROOT   ->  ProcessLine ( Form(".! mkdir -p %s",(fFolder+TString("/plots/BarlowCheck/2D/")).Data()) );
    gROOT   ->  ProcessLine ( Form(".! mkdir -p %s",(fFolder+TString("/plots/BinByBinCheck/2D/")).Data()) );
    auto fRelevantVariations    =  uIsRelevantVariation(hStandard,hVariations,(fFolder+TString("/plots/BarlowCheck/2D/")).Data(),"2D",true);
    //
    // --------- USE THE SYSTEMATICAL RELEVANT VARIATIONS TO DETERMINE THE UNCERTAINTY
    //
    //for ( [...] ) Loop to exclude non Barlow variations
    //
    auto    uStackSystematic    =   uBuildSystematicStack(hStandard,hVariations,fRelevantVariations);
    //
    gROOT->SetBatch();
    TLegend    *lLegend =   new TLegend(0.2,0.85,0.4,0.7);
    lLegend->AddEntry((TH1F*)uStackSystematic.at(0)->GetHists()->At(0),"#mu contr.","F");
    lLegend->AddEntry((TH1F*)uStackSystematic.at(0)->GetHists()->At(1),"#sigma contr.","F");
    TCanvas*c1 = new TCanvas();
    c1->Divide(uStackSystematic.size());
    auto i = 0;
    for ( auto fStack : uStackSystematic )  {
        c1->cd(i);
        gStyle->SetOptStat(0);
        gPad->SetLogx();
        uStackSystematic.at(i)->SetMinimum(0);
        uStackSystematic.at(i)->SetMaximum(uStackSystematic.at(i)->GetMaximum()*1.4);
        uStackSystematic.at(i)->Draw("");
        lLegend->Draw();
        c1->SaveAs((fFolder+TString("/plots/Full_2D_Sys_")+TString(Form("%i.pdf",i))).Data());
        gROOT->SetBatch(kFALSE);
        i++;
    }
    //
    TFile      *fOutput =   new TFile   (Form("%s%s",fFolder.Data(),"/2D_Systematic.root"),"recreate");
    uBuildSystematicError(hStandard,hVariations,(fFolder+TString("/plots/BinByBinCheck/2D/")).Data(),"2D",fRelevantVariations,true)->Write();
    fOutput->Close();
    //
}
//
void
GeneralAnalysis
( ) {
    //
    // --------- RECOVER CALCULATED UNCERTAINTIES
    //
    SetStyle();
    fSetAllBins ();
    gROOT   ->  ProcessLine ( Form(".! mkdir -p %s",(TString(Form(kAnalysis_Systemt_Dir,"yield"))).Data()) );
    gROOT   ->  ProcessLine ( Form(".! mkdir -p %s",(TString(Form(kSystematicsPlot,"yield"))).Data()) );
    //
    // --------- BUILDING CONSTANT UNCERTAINTIES HISTOGRAMS
    //
    // ---------  --------- BRANCHING RATIO
    //
    TH1F   *h1DBranchingRatio     =   new TH1F    ("h1DBranchingRatio","h1DBranchingRatio",nBinPT1D,fArrPT1D);
    for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )    {
        h1DBranchingRatio->SetBinContent(iPT1D+1,kSysLow_BR);
    }
    h1DBranchingRatio           ->  SetLineStyle(1);
    h1DBranchingRatio           ->  SetLineColor(colors[1]);
    //
    TH2F   *h2DBranchingRatio   =   new TH2F    ("h2DBranchingRatio","h2DBranchingRatio",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )    {
        for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )    {
            h2DBranchingRatio->SetBinContent(iPT2D+1,jPT2D+1,kSysLow_BR*2);
        }
    }
    //
    // ---------  --------- GLOBAL TRACKING
    //
    TH1F   *hRTBranchingRatio     =   new TH1F    ("hRTBranchingRatio","hRTBranchingRatio",2,0,2);
    for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )    {
        hRTBranchingRatio->SetBinContent(iPT1D+1,kSysLow_BR);
    }
    //
    TH1F   *h1DITSTPCMatch_     =   new TH1F    ("h1DITSTPCMatch_","h1DITSTPCMatch_",nBinPT1D,fArrPT1D);
    for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )    {
        h1DITSTPCMatch_->SetBinContent(iPT1D+1,.08);
    }
    h1DITSTPCMatch_           ->  SetLineStyle(1);
    h1DITSTPCMatch_           ->  SetLineColor(colors[5]);
    //
    TH2F   *h2DITSTPCMatch_   =   new TH2F    ("h2DITSTPCMatch_","h2DITSTPCMatch_",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )    {
        for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )    {
            h2DITSTPCMatch_->SetBinContent(iPT2D+1,jPT2D+1,.16);
        }
    }
    //
    TH1F   *hRTITSTPCMatch_     =   new TH1F    ("hRTITSTPCMatch_","hRTITSTPCMatch_",2,0,2);
    hRTITSTPCMatch_->SetBinContent(1,0.);
    hRTITSTPCMatch_->SetBinContent(2,0.);
    //
    // --------- RECOVERING CALCULATED UNCERTAINTIES HISTOGRAMS
    //
    //
    // ---------  --------- SIGNAL EXTRACTION
    TFile      *f1DSignalExtraction     =   new TFile(Form("%s%s",Form(kAnalysis_SgExSys_Dir,"yield/Systematics/Standard/"),"/1D_Systematic.root"));
    TH1F       *h1DSignalExtraction     =   (TH1F*)(f1DSignalExtraction->Get("hRAW_1D"));
    h1DSignalExtraction                 ->  SetLineStyle(1);
    h1DSignalExtraction                 ->  SetLineColor(colors[2]);
    //
    TFile      *f2DSignalExtraction     =   new TFile(Form("%s%s",Form(kAnalysis_SgExSys_Dir,"yield/Systematics/Standard/"),"/2D_Systematic.root"));
    TH2F       *h2DSignalExtraction     =   (TH2F*)(f2DSignalExtraction->Get("hRAW_2D"));
    //
    TFile      *fRTSignalExtraction     =   new TFile(Form("%s%s",Form(kAnalysis_SgExSys_Dir,"yield/Systematics/Standard/"),"/RT_Systematic.root"));
    TH1F       *hRTSignalExtraction     =   (TH1F*)(fRTSignalExtraction->Get("hStandard"));
    //
    // ---------  --------- PID
    TFile      *f1DParticleIdentif_     =   new TFile(Form("%s%s",Form(kAnalysis_Systemt_Dir,"yield"),"/PID/1D_Systematic.root"));
    TH1F       *h1DParticleIdentif_     =   (TH1F*)(f1DParticleIdentif_->Get("hRES_1D_Stat"));
    h1DParticleIdentif_                 ->  SetLineStyle(1);
    h1DParticleIdentif_                 ->  SetLineColor(colors[3]);
    //
    TFile      *f2DParticleIdentif_     =   new TFile(Form("%s%s",Form(kAnalysis_Systemt_Dir,"yield"),"/PID/2D_Systematic.root"));
    TH2F       *h2DParticleIdentif_     =   (TH2F*)(f2DParticleIdentif_->Get("h2DStandard"));
    //
    TFile      *fRTParticleIdentif_     =   new TFile(Form("%s%s",Form(kAnalysis_Systemt_Dir,"yield"),"/PID/RT_Systematic.root"));
    TH1F       *hRTParticleIdentif_     =   (TH1F*)(fRTParticleIdentif_->Get("hStandard"));
    //
    // ---------  --------- ANALYSIS CUTS
    TFile      *f1DAnalysisCuts____     =   new TFile(Form("%s%s",Form(kAnalysis_Systemt_Dir,"yield"),"/TRK/1D_Systematic.root"));
    TH1F       *h1DAnalysisCuts____     =   (TH1F*)(f1DAnalysisCuts____->Get("hRES_1D_Stat"));
    h1DAnalysisCuts____                 ->  SetLineStyle(1);
    h1DAnalysisCuts____                 ->  SetLineColor(colors[4]);
    //
    TFile      *f2DAnalysisCuts____     =   new TFile(Form("%s%s",Form(kAnalysis_Systemt_Dir,"yield"),"/TRK/2D_Systematic.root"));
    TH2F       *h2DAnalysisCuts____     =   (TH2F*)(f2DAnalysisCuts____->Get("h2DStandard"));
    //
    TFile      *fRTAnalysisCuts____     =   new TFile(Form("%s%s",Form(kAnalysis_Systemt_Dir,"yield"),"/TRK/RT_Systematic.root"));
    TH1F       *hRTAnalysisCuts____     =   (TH1F*)(fRTAnalysisCuts____->Get("hStandard"));
    //
    // --------- BUILDING TOTAL UNCERTAINTY HISTOGRAM
    //
    TH1F   *h1DTotalSystematic     =   new TH1F    ("h1DTotalSystematic","h1DTotalSystematic",nBinPT1D,fArrPT1D);
    for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )    {
        h1DTotalSystematic->SetBinContent(iPT1D+1,SquareSum( {h1DBranchingRatio->GetBinContent(iPT1D+1), h1DSignalExtraction->GetBinContent(iPT1D+1), h1DParticleIdentif_->GetBinContent(iPT1D+1), h1DAnalysisCuts____->GetBinContent(iPT1D+1), h1DITSTPCMatch_->GetBinContent(iPT1D+1)} ));
    }
    h1DTotalSystematic->SetLineWidth(2);
    h1DTotalSystematic->SetLineColor(kBlack);
    h1DTotalSystematic->SetMinimum(.0);
    h1DTotalSystematic->SetMaximum(1.3*h1DTotalSystematic->GetMaximum());
    //
    TH2F   *h2DTotalSystematic     =   new TH2F    ("h2DTotalSystematic","h2DTotalSystematic",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )    {
        for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )    {
            h2DTotalSystematic->SetBinContent(iPT2D+1,jPT2D+1,SquareSum( {h2DBranchingRatio->GetBinContent(iPT2D+1,jPT2D+1), h2DSignalExtraction->GetBinContent(iPT2D+1,jPT2D+1), h2DParticleIdentif_->GetBinContent(iPT2D+1,jPT2D+1), h2DAnalysisCuts____->GetBinContent(iPT2D+1,jPT2D+1), h2DITSTPCMatch_->GetBinContent(iPT2D+1,jPT2D+1) } ));
        }
    }
    //
    TH1F   *hRTTotalSystematic     =   new TH1F    ("hRTTotalSystematic","hRTTotalSystematic",2,0,2);
    for ( Int_t iPTRT = 0; iPTRT < 2; iPTRT++ )    {
        hRTTotalSystematic->SetBinContent(iPTRT+1,SquareSum( { iPTRT == 0? kSysLow_BR : 0, hRTSignalExtraction->GetBinContent(iPTRT+1), hRTParticleIdentif_->GetBinContent(iPTRT+1), hRTAnalysisCuts____->GetBinContent(iPTRT+1)} ));
    }
    //
    // --------- FINAL CANVAS 1D
    //
    gROOT->SetBatch();
    //
    TLegend    *lLegend                 =   new TLegend(0.17,0.87,0.42,0.72);
    lLegend     ->  SetFillColorAlpha(kWhite,0.);
    lLegend     ->  SetLineColorAlpha(kWhite,0.);
    lLegend     ->  SetNColumns(2);
    lLegend     ->  AddEntry(h1DTotalSystematic,    "Total","L");
    lLegend     ->  AddEntry(h1DBranchingRatio,     "Branching Ratio","L");
    lLegend     ->  AddEntry(h1DSignalExtraction,   "Signal Extraction","L");
    lLegend     ->  AddEntry(h1DParticleIdentif_,   "PID","L");
    lLegend     ->  AddEntry(h1DITSTPCMatch_,       "Global Tracking","L");
    lLegend     ->  AddEntry(h1DAnalysisCuts____,   "Analysis cuts","L");
    //
    TCanvas    *cFullSyst               =   new TCanvas();
    gPad                                ->  SetLogx();
    h1DTotalSystematic                  ->  Draw();
    h1DBranchingRatio                   ->  Draw("same");
    h1DSignalExtraction                 ->  Draw("same");
    h1DParticleIdentif_                 ->  Draw("same");
    h1DAnalysisCuts____                 ->  Draw("same");
    h1DITSTPCMatch_                     ->  Draw("same");
    lLegend                             ->  Draw("same");
    cFullSyst                           ->  SaveAs(Form("%s/1DFullSyst.pdf",Form(kSystematicsPlot,"yield")));
    delete cFullSyst;
    //
    // --------- FINAL CANVAS 2D
    //
    gROOT->SetBatch();
    //
    for ( Int_t iPT = 0; iPT < nBinPT2D; iPT++ )    {
        cFullSyst                       =   new TCanvas();
        gPad                                ->  SetLogx();
        //
        auto hTotal = h2DTotalSystematic  ->  ProjectionY(Form("2DFL_%i",iPT),iPT+1,iPT+1);
        hTotal  ->SetLineWidth(2);
        hTotal  ->SetLineColor(kBlack);
        hTotal  ->SetMinimum(.0);
        hTotal  ->SetMaximum(1.3*hTotal->GetMaximum());
        hTotal  ->Draw();
        //
        auto hBranch=   h2DBranchingRatio   ->  ProjectionY(Form("2DBR_%i",iPT),iPT+1,iPT+1);
        hBranch ->  SetLineStyle(1);
        hBranch ->  SetLineColor(colors[1]);
        hBranch ->  Draw("SAME");
        //
        auto hSigEx =   h2DSignalExtraction ->  ProjectionY(Form("2DSE_%i",iPT),iPT+1,iPT+1);
        hSigEx  ->  SetLineStyle(1);
        hSigEx  ->  SetLineColor(colors[2]);
        hSigEx  ->  Draw("SAME");
        //
        auto hPID =   h2DParticleIdentif_ ->  ProjectionY(Form("2DPD_%i",iPT),iPT+1,iPT+1);
        hPID  ->  SetLineStyle(1);
        hPID  ->  SetLineColor(colors[3]);
        hPID  ->  Draw("SAME");
        //
        auto hAnCut =   h2DAnalysisCuts____ ->  ProjectionY(Form("2DAC_%i",iPT),iPT+1,iPT+1);
        hAnCut  ->  SetLineStyle(1);
        hAnCut  ->  SetLineColor(colors[4]);
        hAnCut  ->  Draw("SAME");
        //
        auto hITSPC =   h2DITSTPCMatch_ ->  ProjectionY(Form("2DIT_%i",iPT),iPT+1,iPT+1);
        hITSPC  ->  SetLineStyle(1);
        hITSPC  ->  SetLineColor(colors[5]);
        hITSPC  ->  Draw("SAME");
        //
        lLegend                                 ->  Draw("same");
        cFullSyst                       ->  SaveAs(Form("%s/2DFullSyst_%i.pdf",Form(kSystematicsPlot,"yield"),iPT));
        delete cFullSyst;
    }
    //
    TFile      *fIn_Data        =   new TFile   (Form(kASigExtp_FitCheckRst,"Yield"));
    TH1F       *h1DYieldHisto   =   (TH1F*)(fIn_Data->Get("hRES_1D_Stat"));
    //
    for ( Int_t iX = 0; iX < h1DYieldHisto->GetNbinsX(); iX++ )   {
        auto fIgnore = .08;
        h1DYieldHisto->SetBinError(iX+1,sqrt( h1DTotalSystematic->GetBinContent(iX+1)*h1DTotalSystematic->GetBinContent(iX+1) - fIgnore*fIgnore )*h1DYieldHisto->GetBinContent(iX+1));
    }
    //
    TH2F       *h2DYieldHisto   =   new TH2F("h2DYieldHisto","h2DYieldHisto",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    for ( Int_t iX = 0; iX < h2DYieldHisto->GetNbinsX(); iX++ )   {
        auto fSlice =   (TH1F*)(fIn_Data->Get(Form("hRES_2D_Cond1_Stat_%i",iX)));
        for ( Int_t iY = 0; iY < h2DYieldHisto->GetNbinsY(); iY++ )   {
            auto fIgnore = .16;
            h2DYieldHisto->SetBinContent    (iX+1,iY+1,fSlice->GetBinContent(iY+1));
            h2DYieldHisto->SetBinError      (iX+1,iY+1,sqrt( h2DTotalSystematic->GetBinContent(iX+1,iY+1)*h2DTotalSystematic->GetBinContent(iX+1,iY+1) - fIgnore*fIgnore )*fSlice->GetBinContent(iY+1));
        }
    }
    //
    auto    k1DStdIntegralErr  =   0.;
    auto    k2DStdIntegralErr  =   0.;
    auto    k1DStdIntegral     =   h1DYieldHisto->IntegralAndError(-1,10000,k1DStdIntegralErr,"width");
    auto    k2DStdIntegral     =   h2DYieldHisto->IntegralAndError(-1,10000,-1,10000,k2DStdIntegralErr,"width");
            k1DStdIntegralErr /=   k1DStdIntegral;
            k2DStdIntegralErr /=   k2DStdIntegral;
    //
    hRTTotalSystematic->GetYaxis()->SetTitle("Systematic uncertainty (%)");
    hRTTotalSystematic->GetXaxis()->SetNdivisions(2);
    hRTTotalSystematic->GetXaxis()->SetBinLabel(hRTTotalSystematic->GetXaxis()->FindBin(0.5),"#LT Y_{2#phi} #GT / #LT Y_{1#phi} #GT");
    hRTTotalSystematic->GetXaxis()->SetBinLabel(hRTTotalSystematic->GetXaxis()->FindBin(1.5),"#LT Y_{2#phi} #GT / #LT Y_{1#phi} #GT^{2}");
    hRTTotalSystematic->GetXaxis()->LabelsOption("h");
    hRTTotalSystematic->SetMarkerColor(colors[3]);
    hRTTotalSystematic->SetLineWidth(3);
    hRTTotalSystematic->SetMarkerStyle(markers[3]);
    hRTTotalSystematic->SetBinError  (1,0);
    hRTTotalSystematic->SetBinError  (2,0);
    hRTTotalSystematic->Scale(100);
    //
    TH1F   *hSimple     =   new TH1F("hLinear",     "", 2,  0,  2);
    hSimple->SetMarkerColor(colors[2]);
    hSimple->SetLineWidth(3);
    hSimple->SetMarkerStyle(markers[4]);
    hSimple->SetBinContent(1,k1DStdIntegralErr+k2DStdIntegralErr);
    hSimple->SetBinContent(2,2*k1DStdIntegralErr+k2DStdIntegralErr);
    hSimple->SetBinError  (1,0);
    hSimple->SetBinError  (2,0);
    hSimple->Scale(100);
    //
    TH1F   *hSquare     =   new TH1F("hSquare",     "", 2,  0,  2);
    hSquare->SetMarkerColor(colors[1]);
    hSquare->SetLineWidth(3);
    hSquare->SetMarkerStyle(markers[5]);
    hSquare->SetBinContent(1,SquareSum( {k1DStdIntegralErr,k2DStdIntegralErr} ));
    hSquare->SetBinContent(2,SquareSum( {2*k1DStdIntegralErr,k2DStdIntegralErr} ));
    hSquare->SetBinError  (1,0);
    hSquare->SetBinError  (2,0);
    hSquare->Scale(100);
    //
    auto fMaximum  = max ( hSimple->GetMaximum(), hRTTotalSystematic->GetMaximum() );
    hRTTotalSystematic->SetMaximum(fMaximum*1.3);
    //
    TCanvas    *c1 = new TCanvas("Ratio");
    //
    auto lLegen2     =   new TLegend(0.18,0.82,0.33,0.72);
    lLegen2     ->  AddEntry( hRTTotalSystematic,    "Ratio err.", "P" );
    lLegen2     ->  AddEntry( hSimple,      "Linear err.", "P" );
    lLegen2     ->  AddEntry( hSquare,      "Square err.", "P" );
    //
    hRTTotalSystematic->Draw("][ EP MIN0");
    hSimple->Draw("SAME EP ][");
    hSquare->Draw("SAME EP ][");
    lLegen2->Draw("same");
    //
    c1          ->  SaveAs(Form("%s/RTFullSyst.pdf",Form(kSystematicsPlot,"yield")));
    delete c1;
    
    
    
    /*
    
    
    auto    iTer    =   0;
    std::vector<Float_t>   kSimpleRatio;
    std::vector<Float_t>   kSquareRatio;
    auto    k1DStdIntegralErr  =   0.;
    auto    k2DStdIntegralErr  =   0.;
    auto    k1DStdIntegral     =   h1DStandard->IntegralAndError(-1,10000,k1DStdIntegralErr,"width");
    auto    k2DStdIntegral     =   h2DStandard->IntegralAndError(-1,10000,-1,10000,k2DStdIntegralErr,"width");
            k1DStdIntegralErr /=   k1DStdIntegral;
            k2DStdIntegralErr /=   k2DStdIntegral;
    //
    for ( auto hVariation : h1DVariations )   {
        //
        if ( h1DEfficiency ) h1DVariations.at(iTer)->Divide(h1DEfficiency);
        if ( h2DEfficiency ) h2DVariations.at(iTer)->Divide(h2DEfficiency);
        auto    k1Dintegral     =   h1DVariations.at(iTer)->Integral("width");
        auto    k2Dintegral     =   h2DVariations.at(iTer)->Integral("width");
        //
        kSimpleRatio.push_back((k1DStdIntegral*k2Dintegral)/(k1Dintegral*k2DStdIntegral)-1);
        kSquareRatio.push_back((k1DStdIntegral*k1DStdIntegral*k2Dintegral)/(k1Dintegral*k1Dintegral*k2DStdIntegral)-1);
        iTer++;
    }
    //
    TH1F       *hSimpleRatioError   =   uBuildTH1F(kSimpleRatio,2000,0,-0.5,0.5);
    auto        fSimpleRatioError   =   0.;
    TH1F       *hSquareRatioError   =   uBuildTH1F(kSquareRatio,2000,0,-0.5,0.5);
    auto        fSquareRatioError   =   0.;
    fSimpleRatioError   +=  fabs(hSimpleRatioError->GetMean());
    fSimpleRatioError   +=  hSimpleRatioError->GetRMS();
    fSquareRatioError   +=  fabs(hSquareRatioError->GetMean());
    fSquareRatioError   +=  hSquareRatioError->GetRMS();
    //
    TH1F   *hStandard   =   new TH1F("hStandard",   "", 2,  0,  2);
    hStandard->GetYaxis()->SetTitle("Systematic uncertainty (%)");
    hStandard->GetXaxis()->SetNdivisions(2);
    hStandard->GetXaxis()->SetBinLabel(hStandard->GetXaxis()->FindBin(0.5),"#LT Y_{2#phi} #GT / #LT Y_{1#phi} #GT");
    hStandard->GetXaxis()->SetBinLabel(hStandard->GetXaxis()->FindBin(1.5),"#LT Y_{2#phi} #GT / #LT Y_{1#phi} #GT^{2}");
    hStandard->GetXaxis()->LabelsOption("h");
    hStandard->SetMarkerColor(colors[3]);
    hStandard->SetLineWidth(3);
    hStandard->SetMarkerStyle(markers[3]);
    hStandard->SetBinContent(1,fSimpleRatioError);
    hStandard->SetBinContent(2,fSquareRatioError);
    hStandard->SetBinError  (1,0);
    hStandard->SetBinError  (2,0);
    hStandard->Scale(100);
    TH1F   *hSimple     =   new TH1F("hLinear",     "", 2,  0,  2);
    hSimple->SetMarkerColor(colors[2]);
    hSimple->SetLineWidth(3);
    hSimple->SetMarkerStyle(markers[4]);
    hSimple->SetBinContent(1,k1DStdIntegralErr+k2DStdIntegralErr);
    hSimple->SetBinContent(2,2*k1DStdIntegralErr+k2DStdIntegralErr);
    hSimple->SetBinError  (1,0);
    hSimple->SetBinError  (2,0);
    hSimple->Scale(100);
    TH1F   *hSquare     =   new TH1F("hSquare",     "", 2,  0,  2);
    hSquare->SetMarkerColor(colors[1]);
    hSquare->SetLineWidth(3);
    hSquare->SetMarkerStyle(markers[5]);
    hSquare->SetBinContent(1,SquareSum( {k1DStdIntegralErr,k2DStdIntegralErr} ));
    hSquare->SetBinContent(2,SquareSum( {2*k1DStdIntegralErr,k2DStdIntegralErr} ));
    hSquare->SetBinError  (1,0);
    hSquare->SetBinError  (2,0);
    hSquare->Scale(100);
    //
    auto fMaximum  = 100 * max ( 2*k1DStdIntegralErr+k2DStdIntegralErr, max (fSimpleRatioError, fSquareRatioError ) );
    hStandard->SetMaximum(fMaximum*1.3);
    //
    gROOT->SetBatch();
    TCanvas    *c1 = new TCanvas("Ratio");
    //
    TLegend    *lLegend = new TLegend(0.18,0.82,0.33,0.72);
    lLegend     ->  AddEntry( hStandard,    "Ratio err.", "P" );
    lLegend     ->  AddEntry( hSimple,      "Linear err.", "P" );
    lLegend     ->  AddEntry( hSquare,      "Square err.", "P" );
    //
    hStandard->Draw("][ EP MIN0");
    hSimple->Draw("SAME EP ][");
    hSquare->Draw("SAME EP ][");
    lLegend->Draw("same");
    //
    cFullSyst                       ->  SaveAs(Form("%s/RTFullSyst_%i.pdf",Form(kSystematicsPlot,"yield"),iPT));
    delete c1;
    */
    
    
    
    
    //
    gROOT->SetBatch(kFALSE);
    //
    TFile      *fOutput =   new TFile( Form("%s/Full_Systematics.root",(TString(Form(kAnalysis_Systemt_Dir,"yield"))).Data()), "RECREATE" );
    //
    h1DTotalSystematic->Write();
    h2DTotalSystematic->Write();
    hRTTotalSystematic->Write();
    //
    fOutput->Close();
}
//




















    /*
    auto    uWhichIsToBeconsidered  =   uBuildSystematicStack(h2D_Syst[0],uVariation2);
    //
    TH2F* fCheckMap = (TH2F*)h2D_Syst[0]->Clone();
    TCanvas*    cCanvas =   new TCanvas("cCanvas","cCanvas",1500,1200);
    cCanvas->Divide(5,4);
    cCanvas->cd(1);
    gPad->SetLogx();
    uBuildSystematicStack(h1D_Syst[0],uVariations)->Draw("HIST F");
    for ( Int_t iVar = 0; iVar < uWhichIsToBeconsidered.size(); iVar++ )    {
        cCanvas->cd(iVar+2);
        gPad->SetLogx();
        for ( Int_t jVar = 0; jVar < uWhichIsToBeconsidered.size(); jVar++ )    {
            auto fHistArray = uWhichIsToBeconsidered.at(iVar)->GetStack();
            auto fYValue = ((TH1F*)fHistArray->At(0))->GetBinContent(jVar+1);
            fYValue += ((TH1F*)fHistArray->At(1))->GetBinContent(jVar+1);
            fCheckMap->SetBinContent(iVar+1,jVar+1,fYValue);
        }
        auto f2Dhist = uBuildSystAndStatStack(h2D_Syst[0],uWhichIsToBeconsidered.at(iVar)(;
        f2Dhist->Draw("HIST F");
    }
    //
    TCanvas    *hSystErrMap =   new TCanvas();
    gPad->SetLogx();
    gPad->SetLogy();
    fCheckMap->Draw("COLZ");
    //
    uBuildSystAndStatStack(h1D_Syst[0],uBuildSystematicStack(h1D_Syst[0],uVariations));
    //
    TFile  *fOut    =   new TFile("Test.root","RECREATE");
    fCheckMap->Write();
    fOut->Close();
    //
     */
/*
    
    // Creating the histograms-------------------------------------------------------------------------------
    //
    hName               =   Form("h1D_Syst_Bin_StSy");
    hTitle              =   Form("Fractional variation of raw yield for bin of PT [%.2f#;%.2f]",fArrPT1D[0],fArrPT1D[nBinPT1D]);
    TH1F   *h1D_Syst_Bin_StSy   =   new TH1F (hName,hTitle,1000,-.5,.5);
    //
    hName               =   Form("h2D_Syst_Bin_StSy");
    hTitle              =   Form("Fractional variation of raw yield for bin of PT [%.2f#;%.2f]",fArrPT1D[0],fArrPT1D[nBinPT1D]);
    TH1F   *h2D_Syst_Bin_StSy   =   new TH1F (hName,hTitle,400,-2.,2.);
    //
    hName               =   Form("h1D_Stat_Bin");
    hTitle              =   Form("h1D_Stat_Bin");
    TH1F   *h1D_Stat_Bin    =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    //
    TH1F  **h1D_Syst_Bin    =   new TH1F   *[nBinPT1D+1];
    hName               =   Form("h1D_Syst_Bin_PT_%.2f_%.2f",fArrPT1D[0],fArrPT1D[nBinPT1D]);
    hTitle              =   Form("Fractional variation of raw yield for bin of PT [%.2f#;%.2f]",fArrPT1D[0],fArrPT1D[nBinPT1D]);
    h1D_Syst_Bin[0]  =   new TH1F (hName,hTitle,1000,-.5,.5);
    h1D_Syst_Bin[0]  ->SetTitle(Form("Fractional variation of raw yield for bin %i",-1));
    h1D_Syst_Bin[0]  ->GetXaxis()->SetTitle("Fractional variation");
    for ( Int_t iAll = 1; iAll <= nBinPT1D; iAll++ )
    {
        hName               =   Form("h1D_Syst_Bin_PT_%.2f_%.2f",fArrPT1D[iAll-1],fArrPT1D[iAll]);
        hTitle              =   Form("Fractional variation of raw yield for bin of PT [%.2f#;%.2f]",fArrPT1D[iAll-1],fArrPT1D[iAll]);
        h1D_Syst_Bin[iAll]  =   new TH1F (hName,hTitle,1000,-.5,.5);
        h1D_Syst_Bin[iAll]  ->SetTitle(Form("Fractional variation of raw yield for bin %i",iAll));
        h1D_Syst_Bin[iAll]  ->GetXaxis()->SetTitle("Fractional variation");
    }
    //
    TH2F   *hCheckFull1D        =   new TH2F("hCheckFull1D","hCheckFull1D",nBinPT1D,fArrPT1D,nBinSyst,fArrSyst);
    //
    TH1F   *h2D_Syst_Bin_All;
    hName                       =   Form("h2D_Syst_Bin_PT_%.2f_%.2f_%.2f_%.2f",fArrPT2D[0],fArrPT2D[nBinPT2D],fArrPT2D[0],fArrPT2D[nBinPT2D]);
    hTitle                      =   Form("Fractional variation of raw yield for bin of PT [%.2f#;%.2f] [%.2f#;%.2f]",fArrPT2D[0],fArrPT2D[nBinPT2D],fArrPT2D[0],fArrPT2D[nBinPT2D]);
    h2D_Syst_Bin_All            =   new TH1F    (hName,hTitle,400,-2.,2.);
    h2D_Syst_Bin_All            ->  SetTitle(Form("Fractional variation of raw yield for bin %i",-1));
    h2D_Syst_Bin_All            ->  GetXaxis()  ->  SetTitle("Fractional variation");
    //
    TH2F   *h2D_Stat_Bin;
    hName                       =   Form("h2D_Stat_Bin");
    hTitle                      =   Form("h2D_Stat_Bin");
    h2D_Stat_Bin                =   new TH2F    (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    //
    TH1F ***h2D_Syst_Bin        =   new TH1F  **[nBinPT2D];
    TH2F  **hCheckFull2D        =   new TH2F   *[nBinPT2D];
    for ( Int_t iAll = 0; iAll < nBinPT2D; iAll++ ) {
        h2D_Syst_Bin[iAll]      =   new TH1F   *[nBinPT2D];
        hName                   =   Form("hCheckFull2D_%i",iAll);
        hCheckFull2D[iAll]      =   new TH2F(hName,hName,nBinPT2D,fArrPT2D,nBinSyst,fArrSyst);
        for ( Int_t jAll = 0; jAll < nBinPT2D; jAll++ ) {
            hName                       =   Form("h2D_Syst_Bin_PT_%.2f_%.2f_%.2f_%.2f",fArrPT2D[iAll],fArrPT2D[iAll+1],fArrPT2D[jAll],fArrPT2D[jAll+1]);
            hTitle                      =   Form("Fractional variation of raw yield for bin of PT [%.2f#;%.2f] [%.2f#;%.2f]",fArrPT2D[iAll],fArrPT2D[iAll+1],fArrPT2D[jAll],fArrPT2D[jAll+1]);
            h2D_Syst_Bin[iAll][jAll]    =   new TH1F (hName,hTitle,1000,-.5,.5);
            h2D_Syst_Bin[iAll][jAll]    ->  SetTitle(Form("Fractional variation of raw yield for bin %i",iAll+1));
            h2D_Syst_Bin[iAll][jAll]    ->  GetXaxis()  ->  SetTitle("Fractional variation");
        }
    }
    //------------//
    //  ANALYSIS  //
    //------------//
    
    // Output File for Fit Check
    TFile*  outFileFit  =   new TFile("./result/yield/ExtractionSystematics/ExtractionSystematics_CheckRatioAndBins.root","recreate");
    
    //------ 1D Histograms ------//
    //
    TGraphAsymmErrors     **g1D_Stat            =   new TGraphAsymmErrors  *[nOptions+1];
    for ( Int_t iTer = 0; iTer <= nOptions; iTer++ )    {
        auto fCheck = new TH1F (*h1D_Syst[iTer]);
        fCheck->Divide(h1D_Syst[iTer],h1D_Syst[0]);
        fCheck->Write();
        g1D_Stat[iTer]     =   fTH1_to_TGAsymmErrors(h1D_Syst[iTer]);
    }
    //
    TGraphAsymmErrors     **g1D_Stat_VarErr     =   new TGraphAsymmErrors  *[nOptions+1];
    for ( Int_t iTer = 1; iTer <= nOptions; iTer++ )    {
        g1D_Stat_VarErr[iTer]     =   new TGraphAsymmErrors();
    }
    //
    TGraphAsymmErrors  *g1D_Stat_Err        =   new TGraphAsymmErrors   ();
    TGraphAsymmErrors  *g1D_Syst_Err        =   new TGraphAsymmErrors   ();
    for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )  {
        //
        auto    fStandard       =   g1D_Stat[0]         ->GetPointY     (iPT1D);
        auto    fStdError       =   g1D_Stat[0]         ->GetErrorYhigh (iPT1D);
        h1D_Stat_Bin            ->  SetBinContent(iPT1D+1,fStdError/fStandard);
        //
        auto    fXPoint         =   fArrPT1D[iPT1D] + .5*(fArrPT1D[iPT1D+1] - fArrPT1D[iPT1D]);
        auto    fXError         =   .5*(fArrPT1D[iPT1D+1] - fArrPT1D[iPT1D]);
        auto    fYPoint         =   0.;
        auto    fYErrorHig      =   g1D_Stat[0]         ->GetErrorYhigh (iPT1D) / (fStandard);
        auto    fYErrorLow      =   g1D_Stat[0]         ->GetErrorYlow  (iPT1D) / (fStandard);
        g1D_Stat_Err            ->  SetPoint        (iPT1D, fXPoint,    fYPoint);
        g1D_Stat_Err            ->  SetPointError   (iPT1D, fXError,    fXError,    fYErrorHig, fYErrorLow);
        g1D_Syst_Err            ->  SetPoint        (iPT1D, fXPoint,    fYPoint);
        if ( iPT1D == 0  || iPT1D == 19 ) g1D_Syst_Err            ->  SetPointError   (iPT1D, fXError,    fXError,    0.040, 0.040);
        if ( iPT1D >= 1  && iPT1D <= 8  ) g1D_Syst_Err            ->  SetPointError   (iPT1D, fXError,    fXError,    0.010, 0.010);
        if ( iPT1D >= 9  && iPT1D <= 14 ) g1D_Syst_Err            ->  SetPointError   (iPT1D, fXError,    fXError,    0.015, 0.015);
        if ( iPT1D >= 15 && iPT1D <= 18 ) g1D_Syst_Err            ->  SetPointError   (iPT1D, fXError,    fXError,    0.035, 0.035);
        //
        for ( Int_t iTer = 1; iTer <= nOptions; iTer++ )    {
            auto    fVariatin       =   g1D_Stat[iTer]          ->GetPointY     (iPT1D);
            //
            auto    fYErrVrHig      =   g1D_Stat[iTer]         ->GetErrorYhigh (iPT1D) / (fStandard);
            auto    fYErrVrLow      =   g1D_Stat[iTer]         ->GetErrorYlow  (iPT1D) / (fStandard);
            //
            auto    fYFrac          =   fVariatin / fStandard - 1;
            auto    fYFracEHig      =   sqrt( fabs(fYErrorHig*fYErrorHig - fYErrVrHig*fYErrVrHig) );
            auto    fYFracELow      =   sqrt( fabs(fYErrorLow*fYErrorLow - fYErrVrLow*fYErrVrLow) );
            g1D_Stat_VarErr[iTer]   ->  SetPoint        (iPT1D, fXPoint,    fYFrac);
            g1D_Stat_VarErr[iTer]   ->  SetPointError   (iPT1D, fXError,    fXError,    fYFracEHig, fYFracELow);
            //
            if ( fVariatin - fStandard == 0 ) continue;
            if ( fBarlowCheck(fStandard,max(fStandard*fYErrorHig,fStandard*fYErrorLow),fVariatin,max(fStandard*fYErrVrHig,fStandard*fYErrVrLow)) )   continue;
            h1D_Syst_Bin[iPT1D+1]   ->  Fill(fYFrac);
            h1D_Syst_Bin[0]         ->  Fill(fYFrac);
            hCheckFull1D            ->  Fill(fXPoint,fYFrac);
        }
    }
    //
    h1D_Syst_Bin_StSy->Write();
    for ( Int_t iPT1D = 0; iPT1D <= nBinPT1D; iPT1D++ )  {
        h1D_Syst_Bin[iPT1D] ->  Write();
    }
    //
    //------ 2D Histograms ------//
    //
    TGraphAsymmErrors    ***g2D_Stat            =   new TGraphAsymmErrors  **[nOption2+1];
    for ( Int_t iTer = 0; iTer <= nOption2; iTer++ )    {
        auto fCheck = new TH2F (*h2D_Syst[iTer]);
        fCheck->Divide(h2D_Syst[iTer],h2D_Syst[0]);
        fCheck->Write();
        g2D_Stat[iTer]     =   fTH2_to_TGAsymmErrors(h2D_Syst[iTer]);
    }
    //
    TGraphAsymmErrors    ***g2D_Stat_VarErr        =   new TGraphAsymmErrors **[nBinPT2D];
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )    {
        g2D_Stat_VarErr[iPT2D]      =   new TGraphAsymmErrors   *[nOption2+1];
        for ( Int_t iTer = 1; iTer <= nOption2; iTer++ )    {
            g2D_Stat_VarErr[iPT2D][iTer]   =   new TGraphAsymmErrors();
        }
    }
    //
    TH2F*   h2D_Syst_Ful3   =   new TH2F    ("h2D_Syst_Ful3","h2D_Syst_Full_Averaged",      nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    TGraphAsymmErrors     **g2D_Stat_Err        =   new TGraphAsymmErrors  *[nBinPT2D];
    TGraphAsymmErrors     **g2D_Syst_Err        =   new TGraphAsymmErrors  *[nBinPT2D];
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
        g2D_Stat_Err[iPT2D]     =   new TGraphAsymmErrors();
        g2D_Syst_Err[iPT2D]     =   new TGraphAsymmErrors();
        for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )  {
            //
            auto    fStandard       =   g2D_Stat[0][iPT2D]         ->GetPointY     (jPT2D);
            auto    fStdError       =   g2D_Stat[0][iPT2D]         ->GetErrorYhigh (jPT2D);
            h2D_Stat_Bin            ->  SetBinContent(iPT2D+1,jPT2D+1,fStdError/fStandard);
            //
            auto    fXPoint         =   fArrPT2D[jPT2D] + .5*(fArrPT2D[jPT2D+1] - fArrPT2D[jPT2D]);
            auto    fXError         =   .5*(fArrPT2D[jPT2D+1] - fArrPT2D[jPT2D]);
            auto    fYPoint         =   0.;
            auto    fYErrorHig      =   g2D_Stat[0][iPT2D]          ->GetErrorYhigh (jPT2D) / (fStandard);
            auto    fYErrorLow      =   g2D_Stat[0][iPT2D]          ->GetErrorYlow  (jPT2D) / (fStandard);
            g2D_Stat_Err[iPT2D]     ->  SetPoint        (jPT2D, fXPoint,    fYPoint);
            g2D_Stat_Err[iPT2D]     ->  SetPointError   (jPT2D, fXError,    fXError,    fYErrorHig, fYErrorLow);
            g2D_Syst_Err[iPT2D]     ->  SetPoint        (jPT2D, fXPoint,    fYPoint);
            g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.050, 0.050);
            if ( iPT2D == 0  || jPT2D == 0  )       g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.075, 0.075);
            else if ( iPT2D == 9  || jPT2D == 9  )  g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.150, 0.150);
            else if ( iPT2D == 1  || jPT2D == 1  )  g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.075, 0.075);
            else if ( iPT2D == 4  || jPT2D == 4  )  g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.075, 0.075);
            else if ( iPT2D == 6  || jPT2D == 6  )  g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.075, 0.075);
            else if ( iPT2D == 8  || jPT2D == 8  )  g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.075, 0.075);
            //
            // Specific Points
            //
            if ( iPT2D == 0  && jPT2D == 9  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.200, 0.200);
            if ( iPT2D == 9  && jPT2D == 0  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.200, 0.200);
            //
            if ( iPT2D == 0  && jPT2D == 6  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.150, 0.150);
            if ( iPT2D == 6  && jPT2D == 0  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.150, 0.150);
            //
            if ( iPT2D == 2  && jPT2D == 9  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.220, 0.220);
            if ( iPT2D == 9  && jPT2D == 2  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.220, 0.220);
            //
            if ( iPT2D == 3  && jPT2D == 9  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.070, 0.070);
            if ( iPT2D == 9  && jPT2D == 3  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.070, 0.070);
            //
            if ( iPT2D == 5  && jPT2D == 9  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.050, 0.050);
            if ( iPT2D == 9  && jPT2D == 5  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.050, 0.050);
            //
            if ( iPT2D == 1  && jPT2D == 2  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.040, 0.040);
            if ( iPT2D == 2  && jPT2D == 1  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.040, 0.040);
            //
            if ( iPT2D == 1  && jPT2D == 3  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.040, 0.040);
            if ( iPT2D == 3  && jPT2D == 1  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.040, 0.040);
            //
            if ( iPT2D == 1  && jPT2D == 4  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.040, 0.040);
            if ( iPT2D == 4  && jPT2D == 1  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.040, 0.040);
            //
            if ( iPT2D == 7  && jPT2D == 9  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.090, 0.090);
            if ( iPT2D == 9  && jPT2D == 7  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.090, 0.090);
            //
            if ( iPT2D == 2  && jPT2D == 5  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.030, 0.030);
            if ( iPT2D == 5  && jPT2D == 2  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.030, 0.030);
            //
            if ( iPT2D == 2  && jPT2D == 0  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.120, 0.120);
            if ( iPT2D == 0  && jPT2D == 2  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.120, 0.120);
            //
            // Diagonal
            if ( iPT2D == 0  && jPT2D == 0  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.400, 0.400);
            if ( iPT2D == 1  && jPT2D == 1  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.040, 0.040);
            if ( iPT2D == 2  && jPT2D == 2  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.070, 0.070);
            if ( iPT2D == 3  && jPT2D == 3  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.050, 0.050);
            if ( iPT2D == 4  && jPT2D == 4  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.100, 0.100);
            if ( iPT2D == 5  && jPT2D == 5  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.080, 0.080);
            if ( iPT2D == 6  && jPT2D == 6  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.140, 0.140);
            if ( iPT2D == 7  && jPT2D == 7  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.050, 0.050);
            if ( iPT2D == 8  && jPT2D == 8  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.100, 0.100);
            if ( iPT2D == 9  && jPT2D == 9  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.130, 0.130);
            //
            //
            auto fVal = g2D_Syst_Err[iPT2D] -> GetErrorYlow(jPT2D);
            h2D_Syst_Ful3->SetBinContent( iPT2D+1, jPT2D+1, fVal) ;//
            //
            for ( Int_t iTer = 1; iTer <= nOption2; iTer++ )    {
                auto    fVariatin       =   g2D_Stat[iTer][iPT2D]       ->GetPointY     (jPT2D);
                auto    fVarError       =   g2D_Stat[iTer][iPT2D]       ->GetErrorYhigh (jPT2D);
                //
                auto    fYErrVrHig      =   g2D_Stat[iTer][iPT2D]       ->GetErrorYhigh (jPT2D) / (fStandard);
                auto    fYErrVrLow      =   g2D_Stat[iTer][iPT2D]       ->GetErrorYlow  (jPT2D) / (fStandard);
                //
                auto    fYFrac          =   fVariatin / fStandard - 1;
                auto    fYFracEHig      =   sqrt( fabs(fYErrorHig*fYErrorHig - fYErrVrHig*fYErrVrHig) );
                auto    fYFracELow      =   sqrt( fabs(fYErrorLow*fYErrorLow - fYErrVrLow*fYErrVrLow) );
                g2D_Stat_VarErr[iPT2D][iTer]    ->  SetPoint        (jPT2D, fXPoint,    fYFrac);
                g2D_Stat_VarErr[iPT2D][iTer]    ->  SetPointError   (jPT2D, fXError,    fXError,    fYFracEHig, fYFracELow);
                //
                if ( fVariatin - fStandard == 0 ) continue;
                if ( fBarlowCheck(fStandard,max(fStandard*fYErrorHig,fStandard*fYErrorLow),fVariatin,max(fStandard*fYErrVrHig,fStandard*fYErrVrLow)) )   continue;
                h2D_Syst_Bin[iPT2D][jPT2D]  ->  Fill(fYFrac);
                h2D_Syst_Bin_All            ->  Fill(fYFrac);
                hCheckFull2D[iPT2D]         ->  Fill(fXPoint,fYFrac);
            }
        }
    }
    //
    h2D_Syst_Bin_StSy->Write();
    h2D_Syst_Bin_All->Write();
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
        for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )  {
            h2D_Syst_Bin[iPT2D][jPT2D] ->  Write();
        }
    }
    //
    //------ Ratio Histograms ------//
    //
    TH1F * hRatio1DVar  =   new TH1F("hRatio1DVar",     "<Y^{SE}_{#phi}>",          nOptions,0,nOptions);
    TH1F * hRatio2DVar  =   new TH1F("hRatio2DVar",     "<Y^{SE}_{#phi#phi}>",      nOptions,0,nOptions);
    TH1F * hRatio1D     =   new TH1F("hRatio1D",        "",                         50,-.05,.05);
    TH1F * hRatio2D     =   new TH1F("hRatio2D",        "",                         50,-.05,.05);
    hRatio1DVar         -> GetXaxis() -> SetTitle("Variation");
    hRatio2DVar         -> GetXaxis() -> SetTitle("Variation");
    hRatio1D            -> GetXaxis() -> SetTitle("Variation");
    hRatio2D            -> GetXaxis() -> SetTitle("Variation");
    hRatio1DVar         -> GetXaxis() -> SetTitleOffset(1.5);
    hRatio2DVar         -> GetXaxis() -> SetTitleOffset(1.5);
    hRatio1DVar         -> GetYaxis() -> SetTitle("Fractional Deviation");
    hRatio2DVar         -> GetYaxis() -> SetTitle("Fractional Deviation");
    auto iBin = 1;
    for ( auto iName : sOptions )   {
        hRatio1DVar        ->GetXaxis()    ->  SetBinLabel(iBin,iName);
        hRatio2DVar        ->GetXaxis()    ->  SetBinLabel(iBin,iName);
        iBin++;
    }
    
    TH1F * fCheckRatio  =   new TH1F("fCheckRatio", "<Y_{#phi#phi}> / <Y_{#phi}>",nOptions,0,nOptions);
    TH1F * fCheckRatio2 =   new TH1F("fCheckRatio2","<Y_{#phi#phi}> / <Y_{#phi}>^{2}",nOptions,0,nOptions);
    fCheckRatio     -> GetXaxis() -> SetTitle("Variation");
    fCheckRatio2    -> GetXaxis() -> SetTitle("Variation");
    fCheckRatio     -> GetXaxis() -> SetTitleOffset(1.5);
    fCheckRatio2    -> GetXaxis() -> SetTitleOffset(1.5);
    fCheckRatio     -> GetYaxis() -> SetTitle("Fractional Deviation");
    fCheckRatio2    -> GetYaxis() -> SetTitle("Fractional Deviation");
    fCheckRatio->GetXaxis()->LabelsOption("v");
    fCheckRatio2->GetXaxis()->LabelsOption("v");
    TH1F * fCheckRati_  =   new TH1F("fCheckRati_", "",50,-.05,.05);
    TH1F * fCheckRati_2 =   new TH1F("fCheckRati_2","",50,-.05,.05);
    fCheckRati_     -> GetXaxis() -> SetTitle("Fractional Deviation");
    fCheckRati_2    -> GetXaxis() -> SetTitle("Fractional Deviation");
    for ( Int_t iTer = 1; iTer <= nOptions; iTer++ )  {
        auto f1DRefErr =  0.;
        auto f1DRefHst = new TH1F(*h1D_Syst[0]);
        f1DRefHst->Divide(hEFF_1D);
        auto f1DRefInt = f1DRefHst->IntegralAndError(-1,1000,f1DRefErr,"width");
        auto f1DTstErr =  0.;
        auto f1DTstHst = new TH1F(*h1D_Syst[iTer]);
        f1DTstHst->Divide(hEFF_1D);
        auto f1DTstInt = f1DTstHst->IntegralAndError(-1,1000,f1DTstErr,"width");
        auto f2DRefErr =  0.;
        auto f2DRefHst = new TH2F(*h2D_Syst[0]);
        f2DRefHst->Divide(hEFF_2D);
        auto f2DRefInt = f2DRefHst->IntegralAndError(-1,1000,-1,1000,f2DRefErr,"width");
        auto f2DTstErr =  0.;
        auto f2DTstHst = new TH2F(*h2D_Syst[iTer]);
        f2DTstHst->Divide(hEFF_2D);
        auto f2DTstInt = f2DTstHst->IntegralAndError(-1,1000,-1,1000,f2DTstErr,"width");
        auto fRatio1Ref =   f2DRefInt/f1DRefInt;
        auto fRatio1Trg =   f2DTstInt/f1DTstInt;
        auto fRatio2Ref =   f2DRefInt/(f1DRefInt*f1DRefInt);
        auto fRatio2Trg =   f2DTstInt/(f1DTstInt*f1DTstInt);
        fCheckRatio     ->  SetBinContent       (iTer,  fRatio1Trg/fRatio1Ref -1.);
        fCheckRatio2    ->  SetBinContent       (iTer,  fRatio2Trg/fRatio2Ref -1.);
        fCheckRati_     ->  Fill                (fRatio1Trg/fRatio1Ref -1.);
        fCheckRati_2    ->  Fill                (fRatio2Trg/fRatio2Ref -1.);
    }
    TCanvas * c1 = new TCanvas("","",1600,1600);
    c1->Divide(2,2);
    c1->cd(1);
    gStyle->SetOptStat(0);
    fCheckRatio->Draw();
    c1->cd(2);
    gStyle->SetOptStat(0);
    fCheckRatio2->Draw();
    c1->cd(3);
    gStyle->SetOptStat(0);
    fCheckRati_->Draw();
    c1->cd(4);
    gStyle->SetOptStat(0);
    fCheckRati_2->Draw();
    c1->SaveAs("./result/yield/ExtractionSystematics/1D_2D.pdf");
    delete c1;
    TH1F * fChec2Ratio  =   new TH1F("fCheckRatio", "",nOption2-nOptions,nOptions,nOption2);
    TH1F * fChec2Ratio2 =   new TH1F("fCheckRatio2","",nOption2-nOptions,nOptions,nOption2);
    TH1F * fChec2Rati_  =   new TH1F("fCheckRati_", "",50,-.05,.05);
    TH1F * fChec2Rati_2 =   new TH1F("fCheckRati_2","",50,-.05,.05);
    for ( Int_t iTer = nOptions+1; iTer <= nOption2; iTer++ )  {
        auto f1DRefErr =  0.;
        auto f1DRefInt = h1D_Syst[0]->IntegralAndError(-1,1000,f1DRefErr,"width");
        auto f1DTstErr =  0.;
        auto f1DTstInt = h1D_Syst[0]->IntegralAndError(-1,1000,f1DTstErr,"width");
        auto f2DRefErr =  0.;
        auto f2DRefInt = h2D_Syst[0]->IntegralAndError(-1,1000,-1,1000,f2DRefErr,"width");
        auto f2DTstErr =  0.;
        auto f2DTstInt = h2D_Syst[iTer]->IntegralAndError(-1,1000,-1,1000,f2DTstErr,"width");
        auto fRatio1Ref =   f2DRefInt/f1DRefInt;
        auto fRatio1Trg =   f2DTstInt/f1DTstInt;
        auto fRatio2Ref =   f2DRefInt/(f1DRefInt*f1DRefInt);
        auto fRatio2Trg =   f2DTstInt/(f1DTstInt*f1DTstInt);
        fChec2Ratio     ->  SetBinContent       (iTer-nOptions,  fRatio1Trg/fRatio1Ref -1.);
        fChec2Ratio2    ->  SetBinContent       (iTer-nOptions,  fRatio2Trg/fRatio2Ref -1.);
        fChec2Rati_     ->  Fill                (fRatio1Trg/fRatio1Ref -1.);
        fChec2Rati_2    ->  Fill                (fRatio2Trg/fRatio2Ref -1.);
    }
    TCanvas * c2 = new TCanvas("","",1600,1600);
    c2->Divide(2,2);
    c2->cd(1);
    gStyle->SetOptStat(0);
    fChec2Ratio->Draw();
    c2->cd(2);
    gStyle->SetOptStat(0);
    fChec2Ratio2->Draw();
    c2->cd(3);
    gStyle->SetOptStat(0);
    fChec2Rati_->Draw();
    c2->cd(4);
    gStyle->SetOptStat(0);
    fChec2Rati_2->Draw();
    c2->SaveAs("./result/yield/ExtractionSystematics/2D.pdf");
    delete c2;
    //
    gROOT->SetBatch(true);
    // Output File for Fit Check
    TFile*  outFileFi2  =   new TFile("./result/yield/ExtractionSystematics/ExtractionSystematics_MeanAndRMS.root","recreate");
    //
    TH1F*   h1D_Syst_Mean   =   new TH1F    ("h1D_Syst_Mean","h1D_Syst_Mean",nBinPT1D,fArrPT1D);
    TH1F*   h1D_Syst_RMS_   =   new TH1F    ("h1D_Syst_RMS_","h1D_Syst_RMS_",nBinPT1D,fArrPT1D);
    TH1F*   h1D_Syst_Full   =   new TH1F    ("h1D_Syst_Full","h1D_Syst_Full",nBinPT1D,fArrPT1D);
    TH1F*   h1D_Syst_Ful2   =   new TH1F    ("h1D_Syst_Ful2","h1D_Syst_Full",nBinPT1D,fArrPT1D);
    for ( Int_t iPT1D = 1; iPT1D <= nBinPT1D; iPT1D++ )  {
        h1D_Syst_Mean   ->SetBinContent (iPT1D,  fabs( h1D_Syst_Bin[iPT1D] ->  GetMean() ) );
        h1D_Syst_RMS_   ->SetBinContent (iPT1D,  h1D_Syst_Bin[iPT1D] ->  GetRMS());
        h1D_Syst_Full   ->SetBinContent (iPT1D,  h1D_Syst_Bin[iPT1D] ->  GetRMS() + fabs(h1D_Syst_Bin[iPT1D] ->  GetMean() ));
        h1D_Syst_Ful2   ->SetBinContent (iPT1D,  h1D_Syst[0]->GetBinContent(iPT1D));
        h1D_Syst_Ful2   ->SetBinError   (iPT1D,  h1D_Syst[0]->GetBinContent(iPT1D)*g1D_Syst_Err->GetErrorYlow(iPT1D));
    }
    h1D_Stat_Bin    ->  Write();
    h1D_Syst_Mean   ->  Write();
    h1D_Syst_RMS_   ->  Write();
    h1D_Syst_Full   ->  Write();
    h1D_Syst_Ful2   ->  Write();
    //
    TH2F*   h2D_Syst_Mean   =   new TH2F    ("h2D_Syst_Mean","h2D_Syst_Mean",               nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    TH2F*   h2D_Syst_RMS_   =   new TH2F    ("h2D_Syst_RMS_","h2D_Syst_RMS_",               nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    TH2F*   h2D_Syst_Full   =   new TH2F    ("h2D_Syst_Full","h2D_Syst_Full",               nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    TH2F*   h2D_Syst_Ful2   =   new TH2F    ("h2D_Syst_Ful2","h2D_Syst_Full",               nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
        for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )  {
            h2D_Syst_Mean   ->SetBinContent (iPT2D+1,jPT2D+1, fabs( h2D_Syst_Bin[iPT2D][jPT2D] ->  GetMean()) );
            h2D_Syst_RMS_   ->SetBinContent (iPT2D+1,jPT2D+1, h2D_Syst_Bin[iPT2D][jPT2D] ->  GetRMS());
            h2D_Syst_Full   ->SetBinContent (iPT2D+1,jPT2D+1, h2D_Syst_Bin[iPT2D][jPT2D] ->  GetRMS() + fabs(h2D_Syst_Bin[iPT2D][jPT2D] ->  GetMean() ));
            h2D_Syst_Ful2   ->SetBinContent (iPT2D+1,jPT2D+1, h2D_Syst[0]->GetBinContent(iPT2D+1,jPT2D+1));
            h2D_Syst_Ful2   ->SetBinError   (iPT2D+1,jPT2D+1, h2D_Syst[0]->GetBinContent(iPT2D+1,jPT2D+1)*g2D_Syst_Err[iPT2D]->GetErrorYlow(jPT2D));
        }
    }
    SetAxis(h2D_Stat_Bin,"PT 2D");
    h2D_Stat_Bin    ->  Write();
    SetAxis(h2D_Syst_Mean,"PT 2D");
    h2D_Syst_Mean   ->  Write();
    SetAxis(h2D_Syst_RMS_,"PT 2D");
    h2D_Syst_RMS_   ->  Write();
    SetAxis(h2D_Syst_Full,"PT 2D");
    h2D_Syst_Full   ->  Write();
    SetAxis(h2D_Syst_Ful2,"PT 2D");
    h2D_Syst_Ful2   ->  Write();
    SetAxis(h2D_Syst_Ful3,"PT 2D");
    h2D_Syst_Ful3   ->  Write();
    
    TCanvas* cCompare = new TCanvas("","",1800,600);
    cCompare->Divide(3,1);
    cCompare->cd(1);
    gPad->SetLogx();
    gPad->SetLogy();
    h2D_Syst_Full->Scale(100);
    h2D_Syst_Full->DrawCopy("colz text28");
    cCompare->cd(2);
    gPad->SetLogx();
    gPad->SetLogy();
    h2D_Syst_Ful3->Scale(100);
    h2D_Syst_Ful3->Draw("colz text28");
    cCompare->cd(3);
    gPad->SetLogx();
    gPad->SetLogy();
    h2D_Syst_Full->Add(h2D_Syst_Ful3,-1.);
    h2D_Syst_Full->Draw("colz text28");
    cCompare->SaveAs("dddd.pdf");
    delete cCompare;
    //
    TH1F * fChec3Rati_  =   new TH1F("Ratio Errors", "<Y_{#phi#phi}> / <Y_{#phi}>",500,0.,5.);
    fChec3Rati_->GetXaxis()->SetTitle("Relative Uncertainty (%)");
    TH1F * fChec4Rati_  =   new TH1F("Square Combine Errors", "",500,0.,5.);
    TH1F * fChec5Rati_  =   new TH1F("Ratio Errors", "<Y_{#phi#phi}> / <Y_{#phi}>^{2}",500,0.,5.);
    fChec5Rati_->GetXaxis()->SetTitle("Relative Uncertainty (%)");
    TH1F * fChec6Rati_  =   new TH1F("Square Combine Errors", "",500,0.,5.);
    TH1F * fChec7Rati_  =   new TH1F("Linear Combine Errors", "",500,0.,5.);
    TH1F * fChec8Rati_  =   new TH1F("Linear Combine Errors", "",500,0.,5.);
    TLegend * L1 = new TLegend(0.9,0.9,0.5,0.8);
    L1->AddEntry(fChec3Rati_,"Ratio Errors");
    L1->AddEntry(fChec4Rati_,"Square Combined Errors");
    L1->AddEntry(fChec7Rati_,"Linear Combined Errors");
    TCanvas * c3 = new TCanvas("","",1000,500);
    c3->Divide(2,1);
    c3->cd(1);
    gStyle->SetOptStat(0);
    fChec3Rati_->Fill(100*(fCheckRati_->  GetRMS() + fabs(fCheckRati_ ->  GetMean() )));
    fChec3Rati_->SetLineColor(kRed);
    fChec3Rati_->Draw("");
    auto h1D_Syst_Ful2_Err  = 0.;
    auto h1D_Syst_Ful2_Hst = new TH1F(*h1D_Syst_Ful2);
    h1D_Syst_Ful2_Hst->Divide(hEFF_1D);
    auto h1D_Syst_Ful2_Int  = h1D_Syst_Ful2_Hst->IntegralAndError(-1,1000,h1D_Syst_Ful2_Err,"width");
    auto h2D_Syst_Ful2_Err  = 0.;
    auto h2D_Syst_Ful2_Hst = new TH2F(*h2D_Syst_Ful2);
    h2D_Syst_Ful2_Hst->Divide(hEFF_2D);
    auto h2D_Syst_Ful2_Int  = h2D_Syst_Ful2_Hst->IntegralAndError(-1,1000,-1,1000,h2D_Syst_Ful2_Err,"width");
    auto hTarget            = sqrt( (h1D_Syst_Ful2_Err*h1D_Syst_Ful2_Err)/(h1D_Syst_Ful2_Int*h1D_Syst_Ful2_Int) + (h2D_Syst_Ful2_Err*h2D_Syst_Ful2_Err)/(h2D_Syst_Ful2_Int*h2D_Syst_Ful2_Int) );
    fChec4Rati_->Fill(100*hTarget);
    fChec4Rati_->SetLineColor(kBlue);
    fChec4Rati_->Draw("SAME");
    hTarget  = (h1D_Syst_Ful2_Err)/(h1D_Syst_Ful2_Int) + (h2D_Syst_Ful2_Err)/(h2D_Syst_Ful2_Int);
    fChec7Rati_->SetLineColor(kGreen-2);
    fChec7Rati_->Fill(100*hTarget);
    fChec7Rati_->Draw("SAME");
    L1->Draw("same");
    c3->cd(2);
    gStyle->SetOptStat(0);
    fChec5Rati_->Fill(100*(fCheckRati_2->  GetRMS() + fabs(fCheckRati_2 ->  GetMean() )) );
    fChec5Rati_->SetLineColor(kRed);
    fChec5Rati_->Draw("");
    hTarget            = sqrt( 4*(h1D_Syst_Ful2_Err*h1D_Syst_Ful2_Err)/(h1D_Syst_Ful2_Int*h1D_Syst_Ful2_Int) + (h2D_Syst_Ful2_Err*h2D_Syst_Ful2_Err)/(h2D_Syst_Ful2_Int*h2D_Syst_Ful2_Int) );
    fChec6Rati_->Fill(100*hTarget);
    fChec6Rati_->SetLineColor(kBlue);
    fChec6Rati_->Draw("SAME");
    hTarget  = 2*(h1D_Syst_Ful2_Err)/(h1D_Syst_Ful2_Int) + (h2D_Syst_Ful2_Err)/(h2D_Syst_Ful2_Int);
    fChec8Rati_->SetLineColor(kGreen-2);
    fChec8Rati_->Fill(100*hTarget);
    fChec8Rati_->Draw("SAME");
    L1->Draw("same");
    c3->SaveAs("./result/yield/ExtractionSystematics/12D_Check.pdf");
    delete c3;
    //
    TLatex         *latext              =   new TLatex();
    TCanvas        *cDrawComparison     =   new TCanvas("cDrawComparison","");
    TLegend        *cComparisonLegend   =   new TLegend(0.15,0.75,0.25,0.85);
    //
    gStyle                              ->  SetOptStat(0);
    //
    h1D_Stat_Bin                        ->  SetLineWidth(3);
    h1D_Stat_Bin                        ->  SetLineColor(kBlue);
    h1D_Stat_Bin                        ->  SetMinimum(0);
    h1D_Stat_Bin                        ->  SetMaximum(0.07);
    h1D_Syst_RMS_                       ->  SetLineWidth(3);
    h1D_Syst_RMS_                       ->  SetLineColor(kRed);
    //
    cComparisonLegend                   ->  SetLineColorAlpha(1,0.);
    cComparisonLegend                   ->  AddEntry(h1D_Stat_Bin,"Stat.","L");
    cComparisonLegend                   ->  AddEntry(h1D_Syst_RMS_,"Syst.","L");
    //
    h1D_Stat_Bin                        ->  Draw();
    h1D_Syst_RMS_                       ->  Draw("SAME");
    cComparisonLegend                   ->  Draw("SAME");
    cDrawComparison                     ->  SaveAs("./result/yield/ExtractionSystematics/_Syst_Stat_overimp.pdf");
    cDrawComparison                     ->  SaveAs("./result/yield/ExtractionSystematics/_Syst_Stat_overimp.png");
    cDrawComparison                     ->  Write();
    //
    delete cDrawComparison;
    delete cComparisonLegend;
    //
                    cDrawComparison     =   new TCanvas("cDrawComparison","");
    gStyle                              ->  SetOptStat(0);
    //
    TF1            *fFlatDist           =   new TF1("fFlatDist","pol0",-100.,100.);
    //
    h1D_Syst_Mean                       ->  SetLineWidth(3);
    h1D_Syst_Mean                       ->  SetLineColor(kBlue);
    h1D_Syst_Mean                       ->  Fit(fFlatDist,"IMREQ0S");
    //
                    cComparisonLegend   =   new TLegend(0.15,0.75,0.35,0.85);
    cComparisonLegend                   ->  SetLineColorAlpha(1,0.);
    cComparisonLegend                   ->  AddEntry(h1D_Syst_Mean,"Mean of Bin Distr.","L");
    //
    h1D_Syst_Mean                       ->  Draw();
    cComparisonLegend                   ->  Draw("SAME");
    fFlatDist                           ->  Draw("SAME");
    latext                              ->  DrawLatexNDC(0.6, 0.83, Form("MEAN:  %3f",fFlatDist->GetParameter(0)));
    latext                              ->  DrawLatexNDC(0.6, 0.75, Form("ERROR: %3f",fFlatDist->GetParError(0)));
    cDrawComparison                     ->  SaveAs("./result/yield/ExtractionSystematics/_Mean_1D.pdf");
    cDrawComparison                     ->  SaveAs("./result/yield/ExtractionSystematics/_Mean_1D.png");
    cDrawComparison                     ->  Write();
    //
    
    
    //
    TCanvas                *cDrawCollection = new TCanvas("","",1600,1600);
    cDrawCollection->Divide(4,3);
    
    auto Check              =   fMultipleError(g1D_Stat_Err,g1D_Syst_Err,g1D_Stat_VarErr,1,nOptions,sOptions);
    gPad->SetLogx(true);
    auto fMultiGrap1        =   (TMultiGraph*)Check ->  GetPrimitive   ("cDrawAllGraphs");
    fMultiGrap1             ->  SetMaximum(+0.15);
    fMultiGrap1             ->  SetMinimum(-0.15);
    SetAxis(fMultiGrap1,"PT 1D");
    fMultiGrap1             ->  GetYaxis()->SetTitle("Fractional Variation");
    fMultiGrap1             ->  SetTitle(Form("PID Systematic in 1D"));
    Check                   ->  SaveAs("./result/yield/ExtractionSystematics/SE_ERROR_1D.pdf");
    Check                   ->  SaveAs("./result/yield/ExtractionSystematics/SE_ERROR_1D.png");
    Check                   ->  Write();
    cDrawCollection         ->  cd(1);
    Check                   ->  DrawClonePad();
    gPad->SetLogx(true);
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
        auto Check2D        =   fMultipleError  (g2D_Stat_Err[iPT2D],g2D_Syst_Err[iPT2D],g2D_Stat_VarErr[iPT2D],1,nOption2,sOption2);
        gPad->SetLogx(true);
        auto fMultiGrap2    =   (TMultiGraph*)Check2D             ->  GetPrimitive   ("cDrawAllGraphs");
        fMultiGrap2         ->  SetMaximum(+0.5);
        fMultiGrap2         ->  SetMinimum(-0.5);
        SetAxis(fMultiGrap2,"PT DD");
        fMultiGrap2         ->  GetYaxis()->SetTitle("Fractional Variation");
        fMultiGrap2             ->  SetTitle(Form("PID Systematic in 2D, PT %.1f-%.1f",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1]));
        Check2D             ->  SaveAs  (Form("./result/yield/ExtractionSystematics/SE_ERROR_2D_PT_Bin_%.1f_%.1f.pdf",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1]));
        Check2D             ->  SaveAs  (Form("./result/yield/ExtractionSystematics/SE_ERROR_2D_PT_Bin_%.1f_%.1f.png",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1]));
        Check2D             ->  Write   ();
        cDrawCollection     ->  cd(2+iPT2D);
        gPad->SetLogx(true);
        Check2D             ->  DrawClonePad();
    }
    cDrawCollection                   ->  SaveAs("./result/yield/ExtractionSystematics/SE_ERROR_FULL.pdf");
    cDrawCollection                   ->  SaveAs("./result/yield/ExtractionSystematics/SE_ERROR_FULL.png");
    //
    delete  cDrawCollection;
    delete  cComparisonLegend;
    //
    cDrawCollection = new TCanvas("","",1600,1600);
    gPad->SetLogx(true);
    cDrawCollection->Divide(4,3);
    //
    THStack  *hSystStack_1D =   new THStack("","");
    h1D_Syst_Mean           ->  SetFillColorAlpha(2,0.1);
    h1D_Syst_Mean           ->  SetLineColorAlpha(2,1.0);
    h1D_Syst_Mean           ->  SetLineWidth(1.2);
    hSystStack_1D           ->  Add(h1D_Syst_Mean);
    h1D_Syst_RMS_           ->  SetFillColorAlpha(4,0.1);
    h1D_Syst_RMS_           ->  SetLineColorAlpha(4,1.0);
    h1D_Syst_RMS_           ->  SetLineWidth(1.2);
    hSystStack_1D           ->  Add(h1D_Syst_RMS_);
    //
    cComparisonLegend   =   new TLegend(0.15,0.75,0.35,0.85);
    cComparisonLegend       ->  SetLineColorAlpha(1,0.);
    cComparisonLegend       ->  SetFillColorAlpha(1,0.);
    
    TMultiGraph    *cDrawAllGraphs      =   new TMultiGraph("cDrawAllGraphs","");
    //
    g1D_Stat_Err                              ->  SetFillColorAlpha(kGray,0.75);
    g1D_Syst_Err                              ->  SetFillColorAlpha(kGray+2,0.5);
    //
    cDrawAllGraphs                      ->  Add         (g1D_Stat_Err,      "AE2");
    cDrawAllGraphs                      ->  Add         (g1D_Syst_Err,      "AE2");
    //
    //
    cComparisonLegend       ->  AddEntry(g1D_Stat_Err,"Stat. Err.","F");
    cComparisonLegend       ->  AddEntry(g1D_Syst_Err,"Syst. Err.","F");
    cComparisonLegend       ->  AddEntry(h1D_Syst_Mean,"Mean Contr.","F");
    cComparisonLegend       ->  AddEntry(h1D_Syst_RMS_,"RMS Contr.","F");
    //
    TCanvas     *cStackShow =   new TCanvas();
    gPad->SetLogx(true);
    cDrawAllGraphs     ->  Draw("ALP");
    cDrawAllGraphs     ->  GetYaxis()->SetTitle("Fractional Variation");
    cDrawAllGraphs     ->  SetMinimum(0.0);
    cDrawAllGraphs     ->  SetMaximum(0.1);
    SetAxis(cDrawAllGraphs,"PT 1D");
    hSystStack_1D           ->  Draw("same");
    cComparisonLegend       ->  Draw("same");
    
    cDrawCollection         ->  cd(1);
    cStackShow              ->  DrawClonePad();
    //
    cStackShow              ->  SaveAs("./result/yield/ExtractionSystematics/SE_ERROR_SYST_1D.pdf");
    cStackShow              ->  SaveAs("./result/yield/ExtractionSystematics/SE_ERROR_SYST_1D.png");
    //
    delete  cStackShow;
    delete  cComparisonLegend;
    //
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
        //
        cComparisonLegend   =   new TLegend(0.15,0.75,0.35,0.85);
        cComparisonLegend       ->  SetLineColorAlpha(1,0.);
        cComparisonLegend       ->  SetFillColorAlpha(1,0.);
        //
        auto    hSlice_Mean     =   h2D_Syst_Mean->ProjectionY("Mean",iPT2D+1,iPT2D+1);
        auto    hSlice_RMS_     =   h2D_Syst_RMS_->ProjectionY("RMS_",iPT2D+1,iPT2D+1);
        //
        THStack  *hSystStack_2D =   new THStack("","");
        hSlice_Mean             ->  SetFillColorAlpha(2,0.1);
        hSlice_Mean             ->  SetLineColorAlpha(2,1.0);
        hSlice_Mean             ->  SetLineWidth(1.2);
        hSystStack_2D           ->  Add(hSlice_Mean);
        hSlice_RMS_             ->  SetFillColorAlpha(4,0.1);
        hSlice_RMS_             ->  SetLineColorAlpha(4,1.0);
        hSlice_RMS_             ->  SetLineWidth(1.2);
        hSystStack_2D           ->  Add(hSlice_RMS_);
        //
        TMultiGraph    *cDrawAllGraphs      =   new TMultiGraph("cDrawAllGraphs","");
        //
        g2D_Stat_Err[iPT2D]                              ->  SetFillColorAlpha(kGray,0.75);
        g2D_Syst_Err[iPT2D]                              ->  SetFillColorAlpha(kGray+2,0.5);
        //
        cDrawAllGraphs                      ->  Add         (g2D_Stat_Err[iPT2D],      "AE2");
        cDrawAllGraphs                      ->  Add         (g2D_Syst_Err[iPT2D],      "AE2");
        //
        //
        cComparisonLegend       ->  AddEntry(g2D_Stat_Err[iPT2D],"Stat. Err.","F");
        cComparisonLegend       ->  AddEntry(g2D_Syst_Err[iPT2D],"Syst. Err.","F");
        cComparisonLegend       ->  AddEntry(hSlice_Mean,"Mean Contr.","F");
        cComparisonLegend       ->  AddEntry(hSlice_RMS_,"RMS Contr.","F");
        //
                    cStackShow  =   new TCanvas("StackShow","StackShow");
        gPad->SetLogx(true);
        cDrawAllGraphs     ->  Draw("ALP");
        cDrawAllGraphs     ->  GetYaxis()->SetTitle("Fractional Variation");
        cDrawAllGraphs     ->  SetMinimum(0.0);
        cDrawAllGraphs     ->  SetMaximum(0.5);
        SetAxis(cDrawAllGraphs,"PT DD");
        hSystStack_2D           ->  Draw("same");
        cComparisonLegend       ->  Draw("same");
        //
        cDrawCollection         ->  cd(2+iPT2D);
        cStackShow              ->  DrawClonePad();
        //
        cStackShow             ->  SaveAs  (Form("./result/yield/ExtractionSystematics/SE_ERROR_SYST_2D_PT_Bin_%.1f_%.1f.pdf",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1]));
        cStackShow             ->  SaveAs  (Form("./result/yield/ExtractionSystematics/SE_ERROR_SYST_2D_PT_Bin_%.1f_%.1f.png",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1]));
        //
        delete cStackShow;
        delete cComparisonLegend;
    }
    cDrawCollection                   ->  SaveAs("./result/yield/ExtractionSystematics/SE_ERROR_SYST_FULL.pdf");
    cDrawCollection                   ->  SaveAs("./result/yield/ExtractionSystematics/SE_ERROR_SYST_FULL.png");
    //
    delete  cDrawCollection;
    //
    cDrawCollection = new TCanvas("","",1600,1600);
    //
    gPad->SetLogx(true);
    hCheckFull1D->Draw("colz");
    cDrawCollection                   ->  SaveAs("./result/yield/ExtractionSystematics/SE_ERROR_SYST_FULL_1D.pdf");
    cDrawCollection                   ->  SaveAs("./result/yield/ExtractionSystematics/SE_ERROR_SYST_FULL_1D.png");
    //
    delete  cDrawCollection;
    //
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
        cDrawCollection = new TCanvas("","",1600,1600);
        //
        gPad->SetLogx(true);
        hCheckFull2D[iPT2D]->Draw("colz");
        cDrawCollection                   ->  SaveAs(Form("./result/yield/ExtractionSystematics/SE_ERROR_SYST_FULL_2D_%i.pdf",iPT2D));
        cDrawCollection                   ->  SaveAs(Form("./result/yield/ExtractionSystematics/SE_ERROR_SYST_FULL_2D_%i.png",iPT2D));
        //
        delete  cDrawCollection;
    }
    //
    // >-> Close input File
    //
    outFileFit->Close();
    outFileFi2->Close();
    //
    for ( Int_t iTer = 0; iTer <= nOptions; iTer++ )    {
        insFileH1D[iTer]->Close();
    }
    //
    for ( Int_t iTer = 0; iTer < nOption2; iTer++ )    {
        insFileH2D[iTer]->Close();
    }
    //
    gROOT->SetBatch(false);
    return;
}
*/
