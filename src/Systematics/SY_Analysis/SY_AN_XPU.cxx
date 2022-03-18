// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../../inc/AliAnalysisPhiPair.h"
#include "../GeneralAnalysis.cxx"
#include "RooMsgService.h"

void
SY_AN_XPU
( TString fOption = "Yield", TString kFolder = "_p_p__7TeV" )    {
    //
    //-----------------------------//
    //  Setting general analysis   //
    //-----------------------------//
    //
    //  Option chosing
    if ( !fChooseOption(fOption) ) return;
    //
    //  Generating the binning array--------------------------------------------------------------------------
    fSetAllBins();
    //
    if ( kDoYield ) {
        //
        //  --- Recovering Standard Analysis
        TFile*  insFile_Data_YL     =   new TFile   (Form(kASigExtp_FitCheckRst,(TString("Yield")+kFolder).Data()));
        auto        hSyst1D_STD     =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "h1D_Nraw_stat",        "hSyst1D_STD" );
        auto        hSyst2D_STD     =   uLoadHistograms<0,TH2F> ( insFile_Data_YL,  "anSS2D_",              "hSyst2D_STD" );
        auto        hSysURT_STD     =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "hXD_Nyld_stat",        "hSystRT_STD" );
        auto        hSystRT_STD     =   new TH1F( "hSystRT_STD", "hSystRT_STD", 6, 0, 6 );
        auto        k1DYield_V      =   hSysURT_STD->GetBinContent  (1);
        auto        k1DYield_E      =   hSysURT_STD->GetBinError    (1);
        auto        k1DYield_R      =   k1DYield_E / k1DYield_V;
        auto        k2DYield_V      =   hSysURT_STD->GetBinContent  (2);
        auto        k2DYield_E      =   hSysURT_STD->GetBinError    (2);
        auto        k2DYield_R      =   k2DYield_E / k2DYield_V;
        hSystRT_STD ->  SetBinContent   ( 1, k1DYield_V );
        hSystRT_STD ->  SetBinError     ( 1, k1DYield_E );
        hSystRT_STD ->  SetBinContent   ( 2, k2DYield_V );
        hSystRT_STD ->  SetBinError     ( 2, k2DYield_E );
        hSystRT_STD ->  SetBinContent   ( 3, (k2DYield_V/(k1DYield_V)) );
        hSystRT_STD ->  SetBinError     ( 3, (k2DYield_V/(k1DYield_V))*SquareSum( {k2DYield_R,k1DYield_R} ) );
        hSystRT_STD ->  SetBinContent   ( 4, (k2DYield_V/(k1DYield_V*k1DYield_V)) );
        hSystRT_STD ->  SetBinError     ( 4, (k2DYield_V/(k1DYield_V*k1DYield_V))*SquareSum( {k2DYield_R,k1DYield_R,k1DYield_R} ) );
        hSystRT_STD ->  SetBinContent   ( 5, fSigmaPhiValue(k1DYield_V,k2DYield_V) );
        hSystRT_STD ->  SetBinError     ( 5, fSigmaPhiError(k1DYield_V,k2DYield_V,k1DYield_E,k2DYield_E) );
        hSystRT_STD ->  SetBinContent   ( 6, fGammaPhiValue(k1DYield_V,k2DYield_V) );
        hSystRT_STD ->  SetBinError     ( 6, fGammaPhiError(k1DYield_V,k2DYield_V,k1DYield_E,k2DYield_E) );
        hSyst1D_STD ->  SetMarkerStyle  ( uGetMarker(1) );
        hSyst2D_STD ->  SetMarkerStyle  ( uGetMarker(1) );
        hSystRT_STD ->  SetMarkerStyle  ( uGetMarker(1) );
        hSyst1D_STD ->  SetMarkerColor  ( uGetColor(1) );
        hSyst2D_STD ->  SetMarkerColor  ( uGetColor(1) );
        hSystRT_STD ->  SetMarkerColor  ( uGetColor(1) );
        hSyst1D_STD ->  SetLineColor    ( uGetColor(1) );
        hSyst2D_STD ->  SetLineColor    ( uGetColor(1) );
        hSystRT_STD ->  SetLineColor    ( uGetColor(1) );
        //
        //  --- Recovering High Rate Analysis
        TFile*  insFile_Data_HR     =   new TFile   (Form(kASigExtp_FitCheckRst,(TString("Yield")+kFolder+TString("/PileupTest/HighRate/")).Data()));
        auto        hSyst1D_PHR     =   uLoadHistograms<0,TH1F> ( insFile_Data_HR,  "h1D_Nraw_stat",        "hSyst1D_PHR" );
        auto        hSyst2D_PHR     =   uLoadHistograms<0,TH2F> ( insFile_Data_HR,  "anSS2D_",              "hSyst2D_PHR" );
        auto        hSysURT_PHR     =   uLoadHistograms<0,TH1F> ( insFile_Data_HR,  "hXD_Nyld_stat",        "hSystRT_PHR" );
        auto        hSystRT_PHR     =   new TH1F( "hSystRT_PHR", "hSystRT_PHR", 6, 0, 6 );
                    k1DYield_V      =   hSysURT_PHR->GetBinContent  (1);
                    k1DYield_E      =   hSysURT_PHR->GetBinError    (1);
                    k1DYield_R      =   k1DYield_E / k1DYield_V;
                    k2DYield_V      =   hSysURT_PHR->GetBinContent  (2);
                    k2DYield_E      =   hSysURT_PHR->GetBinError    (2);
                    k2DYield_R      =   k2DYield_E / k2DYield_V;
        hSystRT_PHR ->  SetBinContent   ( 1, k1DYield_V );
        hSystRT_PHR ->  SetBinError     ( 1, k1DYield_E );
        hSystRT_PHR ->  SetBinContent   ( 2, k2DYield_V );
        hSystRT_PHR ->  SetBinError     ( 2, k2DYield_E );
        hSystRT_PHR ->  SetBinContent   ( 3, (k2DYield_V/(k1DYield_V)) );
        hSystRT_PHR ->  SetBinError     ( 3, (k2DYield_V/(k1DYield_V))*SquareSum( {k2DYield_R,k1DYield_R} ) );
        hSystRT_PHR ->  SetBinContent   ( 4, (k2DYield_V/(k1DYield_V*k1DYield_V)) );
        hSystRT_PHR ->  SetBinError     ( 4, (k2DYield_V/(k1DYield_V*k1DYield_V))*SquareSum( {k2DYield_R,k1DYield_R,k1DYield_R} ) );
        hSystRT_PHR ->  SetBinContent   ( 5, fSigmaPhiValue(k1DYield_V,k2DYield_V) );
        hSystRT_PHR ->  SetBinError     ( 5, fSigmaPhiError(k1DYield_V,k2DYield_V,k1DYield_E,k2DYield_E) );
        hSystRT_PHR ->  SetBinContent   ( 6, fGammaPhiValue(k1DYield_V,k2DYield_V) );
        hSystRT_PHR ->  SetBinError     ( 6, fGammaPhiError(k1DYield_V,k2DYield_V,k1DYield_E,k2DYield_E) );
        hSyst1D_PHR ->  SetMarkerStyle  ( uGetMarker(1) );
        hSyst2D_PHR ->  SetMarkerStyle  ( uGetMarker(1) );
        hSystRT_PHR ->  SetMarkerStyle  ( uGetMarker(1) );
        hSyst1D_PHR ->  SetMarkerColor  ( uGetColor(1) );
        hSyst2D_PHR ->  SetMarkerColor  ( uGetColor(1) );
        hSystRT_PHR ->  SetMarkerColor  ( uGetColor(1) );
        hSyst1D_PHR ->  SetLineColor    ( uGetColor(1) );
        hSyst2D_PHR ->  SetLineColor    ( uGetColor(1) );
        hSystRT_PHR ->  SetLineColor    ( uGetColor(1) );
        hSyst1D_PHR ->  SetMarkerStyle  ( uGetMarker(2) );
        hSyst2D_PHR ->  SetMarkerStyle  ( uGetMarker(2) );
        hSystRT_PHR ->  SetMarkerStyle  ( uGetMarker(2) );
        hSyst1D_PHR ->  SetMarkerColor  ( uGetColor(2) );
        hSyst2D_PHR ->  SetMarkerColor  ( uGetColor(2) );
        hSystRT_PHR ->  SetMarkerColor  ( uGetColor(2) );
        hSyst1D_PHR ->  SetLineColor    ( uGetColor(2) );
        hSyst2D_PHR ->  SetLineColor    ( uGetColor(2) );
        hSystRT_PHR ->  SetLineColor    ( uGetColor(2) );
        //
        //  --- Recovering Low Rate Analysis
        TFile*  insFile_Data_LR     =   new TFile   (Form(kASigExtp_FitCheckRst,(TString("Yield")+kFolder+TString("/PileupTest/Low_Rate/")).Data()));
        auto        hSyst1D_PLR     =   uLoadHistograms<0,TH1F> ( insFile_Data_LR,  "h1D_Nraw_stat",        "hSyst1D_PLR" );
        auto        hSyst2D_PLR     =   uLoadHistograms<0,TH2F> ( insFile_Data_LR,  "anSS2D_",              "hSyst2D_PLR" );
        auto        hSysURT_PLR     =   uLoadHistograms<0,TH1F> ( insFile_Data_LR,  "hXD_Nyld_stat",        "hSystRT_PLR" );
        auto        hSystRT_PLR     =   new TH1F( "hSystRT_PLR", "hSystRT_PLR", 6, 0, 6 );
                    k1DYield_V      =   hSysURT_PLR->GetBinContent  (1);
                    k1DYield_E      =   hSysURT_PLR->GetBinError    (1);
                    k1DYield_R      =   k1DYield_E / k1DYield_V;
                    k2DYield_V      =   hSysURT_PLR->GetBinContent  (2);
                    k2DYield_E      =   hSysURT_PLR->GetBinError    (2);
                    k2DYield_R      =   k2DYield_E / k2DYield_V;
        hSystRT_PLR ->  SetBinContent   ( 1, k1DYield_V );
        hSystRT_PLR ->  SetBinError     ( 1, k1DYield_E );
        hSystRT_PLR ->  SetBinContent   ( 2, k2DYield_V );
        hSystRT_PLR ->  SetBinError     ( 2, k2DYield_E );
        hSystRT_PLR ->  SetBinContent   ( 3, (k2DYield_V/(k1DYield_V)) );
        hSystRT_PLR ->  SetBinError     ( 3, (k2DYield_V/(k1DYield_V))*SquareSum( {k2DYield_R,k1DYield_R} ) );
        hSystRT_PLR ->  SetBinContent   ( 4, (k2DYield_V/(k1DYield_V*k1DYield_V)) );
        hSystRT_PLR ->  SetBinError     ( 4, (k2DYield_V/(k1DYield_V*k1DYield_V))*SquareSum( {k2DYield_R,k1DYield_R,k1DYield_R} ) );
        hSystRT_PLR ->  SetBinContent   ( 5, fSigmaPhiValue(k1DYield_V,k2DYield_V) );
        hSystRT_PLR ->  SetBinError     ( 5, fSigmaPhiError(k1DYield_V,k2DYield_V,k1DYield_E,k2DYield_E) );
        hSystRT_PLR ->  SetBinContent   ( 6, fGammaPhiValue(k1DYield_V,k2DYield_V) );
        hSystRT_PLR ->  SetBinError     ( 6, fGammaPhiError(k1DYield_V,k2DYield_V,k1DYield_E,k2DYield_E) );
        hSyst1D_PLR ->  SetMarkerStyle  ( uGetMarker(1) );
        hSyst2D_PLR ->  SetMarkerStyle  ( uGetMarker(1) );
        hSystRT_PLR ->  SetMarkerStyle  ( uGetMarker(1) );
        hSyst1D_PLR ->  SetMarkerColor  ( uGetColor(1) );
        hSyst2D_PLR ->  SetMarkerColor  ( uGetColor(1) );
        hSystRT_PLR ->  SetMarkerColor  ( uGetColor(1) );
        hSyst1D_PLR ->  SetLineColor    ( uGetColor(1) );
        hSyst2D_PLR ->  SetLineColor    ( uGetColor(1) );
        hSystRT_PLR ->  SetLineColor    ( uGetColor(1) );
        hSyst1D_PLR ->  SetMarkerStyle  ( uGetMarker(3) );
        hSyst2D_PLR ->  SetMarkerStyle  ( uGetMarker(3) );
        hSystRT_PLR ->  SetMarkerStyle  ( uGetMarker(3) );
        hSyst1D_PLR ->  SetMarkerColor  ( uGetColor(3) );
        hSyst2D_PLR ->  SetMarkerColor  ( uGetColor(3) );
        hSystRT_PLR ->  SetMarkerColor  ( uGetColor(3) );
        hSyst1D_PLR ->  SetLineColor    ( uGetColor(3) );
        hSyst2D_PLR ->  SetLineColor    ( uGetColor(3) );
        hSystRT_PLR ->  SetLineColor    ( uGetColor(3) );
        //
        //  --- Make Output Directory
        gROOT->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_Systemt_Dir,(TString("Yield")+kFolder).Data()))+TString("XPU")).Data()));
        //
        TCanvas*    cDrawComparison_1D  =   new TCanvas("cDrawComparison_1D","cDrawComparison_1D",900,1300);
        //
        auto    kRatio_PHR  =   (TH1F*)(hSyst1D_PHR->Clone());
        kRatio_PHR  ->  Divide( hSyst1D_STD );
        auto    kRatio_PLR  =   (TH1F*)(hSyst1D_PLR->Clone());
        kRatio_PLR  ->  Divide( hSyst1D_STD );
        kRatio_PHR  ->  SetMaximum( 1.1 );
        kRatio_PHR  ->  SetMinimum( 0.9 );
        //
        TLegend*    cLegendComparison_1D    =   new TLegend( 0.18, 0.75, 0.48, 0.88 );
        cLegendComparison_1D    -> AddEntry( hSyst1D_STD, "Standard",   "EP" );
        cLegendComparison_1D    -> AddEntry( hSyst1D_PHR, "High Rate",  "EP" );
        cLegendComparison_1D    -> AddEntry( hSyst1D_PLR, "Low Rate",   "EP" );
        //
        TPad*   kUpperPlot  =   new TPad("kUpperPlot", "kUpperPlot", 0, 0.3, 1, 1.0);
        gStyle      ->  SetOptStat(0);
        kUpperPlot  ->  SetBottomMargin(0);
        kUpperPlot  -> SetLogy();
        kUpperPlot  ->  Draw();
        kUpperPlot  ->  cd();
        hSyst1D_STD -> Draw("SAME");
        hSyst1D_PLR -> Draw("SAME");
        hSyst1D_PHR -> Draw("SAME");
        cLegendComparison_1D->Draw("SAME");
        //
        cDrawComparison_1D-> cd();
        TPad*   kLowerPlot  =   new TPad("kLowerPlot", "kLowerPlot", 0, 0.0, 1, 0.3);
        kLowerPlot      ->  SetGridy();
        gStyle          ->  SetOptStat(0);
        gPad            ->  SetGridy();
        kLowerPlot->SetTopMargin(0);
        kLowerPlot->Draw();
        kLowerPlot->cd();
        kRatio_PHR->Draw("SAME");
        kRatio_PLR->Draw("SAME");
        //
        cDrawComparison_1D  ->  cd();
        //
        cDrawComparison_1D  ->  SaveAs(Form("%s/hPileUpTest_1D.pdf",(TString(Form(kAnalysis_Systemt_Dir,(TString("Yield")+kFolder).Data()))+TString("XPU")).Data()));
        //
        TCanvas*    cDrawComparison_RT  =   new TCanvas("cDrawComparison_1D","cDrawComparison_1D",900,1300);
        //
        kRatio_PHR  =   (TH1F*)(hSystRT_PHR->Clone());
        kRatio_PHR  ->  Divide( hSystRT_STD );
        kRatio_PLR  =   (TH1F*)(hSystRT_PLR->Clone());
        kRatio_PLR  ->  Divide( hSystRT_STD );
        kRatio_PHR  ->  SetMaximum( 1.5 );
        kRatio_PHR  ->  SetMinimum( 0.5 );
        kRatio_PHR  ->  SetTitle("");
        kRatio_PHR  ->  GetXaxis()  ->  SetBinLabel(1,"#frac{dN_{#phi}}{dy}");
        kRatio_PHR  ->  GetXaxis()  ->  SetBinLabel(2,"#frac{dN_{#phi#phi}}{dy}");
        kRatio_PHR  ->  GetXaxis()  ->  SetBinLabel(3,"#frac{#LT Y_{#phi#phi} #GT }{ #LT Y_{#phi} #GT }");
        kRatio_PHR  ->  GetXaxis()  ->  SetBinLabel(4,"#frac{#LT Y_{#phi#phi} #GT }{ #LT Y_{#phi} #GT^{2} }");
        kRatio_PHR  ->  GetXaxis()  ->  SetBinLabel(5,"#sigma_{#phi}");
        kRatio_PHR  ->  GetXaxis()  ->  SetBinLabel(6,"#gamma_{#phi}");
        kRatio_PHR  ->  GetXaxis()  ->  SetLabelSize(0.11);
        kRatio_PHR  ->  GetXaxis()  ->  SetLabelOffset(0.03);
        //
        kUpperPlot  =   new TPad("kUpperPlot", "kUpperPlot", 0, 0.35, 1, 1.0);
        gStyle      ->  SetOptStat(0);
        kUpperPlot  ->  SetBottomMargin(0);
        kUpperPlot  -> SetLogy();
        kUpperPlot  ->  Draw();
        kUpperPlot  ->  cd();
        hSystRT_STD -> Draw("SAME");
        hSystRT_PLR -> Draw("SAME");
        hSystRT_PHR -> Draw("SAME");
        cLegendComparison_1D->Draw("SAME");
        //
        cDrawComparison_RT-> cd();
        kLowerPlot  =   new TPad("kLowerPlot", "kLowerPlot", 0, 0.0, 1, 0.35);
        kLowerPlot  ->  SetGridy();
        kLowerPlot  ->  SetBottomMargin(0.29);
        gStyle      ->  SetOptStat(0);
        gPad        ->  SetGridy();
        kLowerPlot->SetTopMargin(0);
        kLowerPlot->Draw();
        kLowerPlot->cd();
        kRatio_PHR->Draw("SAME");
        kRatio_PLR->SetBinContent   (1, kRatio_PLR->GetBinContent(1));
        kRatio_PLR->SetBinContent   (5, kRatio_PLR->GetBinContent(5));
        kRatio_PHR->SetBinContent   (1, kRatio_PHR->GetBinContent(1));
        kRatio_PHR->SetBinContent   (5, kRatio_PHR->GetBinContent(5));
        kRatio_PHR->Draw("SAME");
        kRatio_PLR->Draw("SAME");
        //
        cDrawComparison_RT  ->  cd();
        //
        cDrawComparison_RT  ->  SaveAs(Form("%s/hPileUpTest_RT.pdf",(TString(Form(kAnalysis_Systemt_Dir,(TString("Yield")+kFolder).Data()))+TString("XPU")).Data()));
        //
        delete  cDrawComparison_RT;
        //
    }
}
