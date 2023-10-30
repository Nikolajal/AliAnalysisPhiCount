// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../../inc/AliAnalysisPhiPair.h"
#include "../GeneralAnalysis.cxx"
#include "RooMsgService.h"
//!
template <  typename TH_1D = TH1F,
            typename TH_2D = TH2F,
            typename TH_RT = TH1F,
            typename TH_PT = TH1F >
std::tuple< TString, TH_1D*, TH_2D*, TH_RT*, TH_PT* >
uLoadCustomSystematics
( TString kFolder, TString kType, TString kCurrentSystematic, TString kNewTag, Int_t iColor, Int_t iLineStyle ) {
    //! Book result variable
    std::tuple< TString, TH_1D*, TH_2D*, TH_RT*, TH_PT* > fResult;
    get<0>(fResult)     = kNewTag;
    //!
    TFile  *InputFile   = new TFile( Form( kAnalysis_Systemt_Dir + TString("/")+kCurrentSystematic+TString("/FullSystematics.root"), (kType+kFolder).Data()) );
    gROOT->cd();
    get<1>(fResult)     = uLoadHistograms<0,TH1F> ( InputFile, "k1DSystematics",    TString("h1D_Syst_") + kCurrentSystematic );
    get<2>(fResult)     = uLoadHistograms<0,TH2F> ( InputFile, "k2DSystematics",    TString("h2D_Syst_") + kCurrentSystematic );
    get<3>(fResult)     = uLoadHistograms<0,TH1F> ( InputFile, "kPTSystematics",    TString("hPT_Syst_") + kCurrentSystematic );
    get<4>(fResult)     = uLoadHistograms<0,TH1F> ( InputFile, "kRatioSystematics", TString("hRT_Syst_") + kCurrentSystematic );
    get<1>(fResult)     -> SetLineStyle(iLineStyle);
    get<2>(fResult)     -> SetLineStyle(iLineStyle);
    get<3>(fResult)     -> SetLineStyle(iLineStyle);
    get<4>(fResult)     -> SetLineStyle(iLineStyle);
    get<1>(fResult)     -> SetLineColor( uGetColor(iColor) );
    get<2>(fResult)     -> SetLineColor( uGetColor(iColor) );
    get<3>(fResult)     -> SetLineColor( uGetColor(iColor) );
    get<4>(fResult)     -> SetLineColor( uGetColor(iColor) );
    //!
    InputFile   -> Close();
    return fResult;
}
//!
template <  typename TH_1D = TH1F,
            typename TH_2D = TH2F,
            typename TH_RT = TH1F,
            typename TH_PT = TH1F >
std::tuple< TString, TH_1D*, TH_2D*, TH_RT*, TH_PT* >
uGenerateFlatSystematics
( TString kFolder, TString kType, TString kCurrentSystematic, TString kNewTag, Float_t k1DValue, Float_t k2DValue, Float_t kPTValue, Float_t k1D_Yield, Float_t k2D_Yield, Int_t iColor, Int_t iLineStyle ) {
    //! Book result variable
    std::tuple< TString, TH_1D*, TH_2D*, TH_RT*, TH_PT* > fResult;
    get<0>(fResult)     = kNewTag;
    //!
    TFile  *InputFile   = new TFile( Form( kAnalysis_Systemt_Dir + TString("/")+kCurrentSystematic+TString("/FullSystematics.root"), ( kType + kFolder).Data()) );
    gROOT->cd();
    get<1>(fResult)     = uLoadHistograms<0,TH1F> ( InputFile, "k1DSystematics",    TString("h1D_Syst_") + kNewTag );
    get<2>(fResult)     = uLoadHistograms<0,TH2F> ( InputFile, "k2DSystematics",    TString("h2D_Syst_") + kNewTag );
    get<3>(fResult)     = uLoadHistograms<0,TH1F> ( InputFile, "kPTSystematics",    TString("hPT_Syst_") + kNewTag );
    get<4>(fResult)     = uLoadHistograms<0,TH1F> ( InputFile, "kRatioSystematics", TString("hRT_Syst_") + kNewTag );
    get<1>(fResult)     = uFlatDistribution(get<1>(fResult),k1DValue);
    get<2>(fResult)     = uFlatDistribution(get<2>(fResult),k2DValue);
    get<3>(fResult)     = uFlatDistribution(get<3>(fResult),kPTValue);
    get<4>(fResult)     = uFlatDistribution(get<4>(fResult),0);
    get<4>(fResult)->SetBinContent(1,k1DValue);
    get<4>(fResult)->SetBinContent(2,k2DValue);
    get<4>(fResult)->SetBinContent(3,0);
    get<4>(fResult)->SetBinContent(4,k1DValue);
    auto P1_Sys = max( fabs( 1- fSigmaPhiValue( k1D_Yield*(1+k1DValue), k2D_Yield*(1+k2DValue) )/fSigmaPhiValue( k1D_Yield, k2D_Yield ) ), fabs( 1- fSigmaPhiValue( k1D_Yield*(1-k1DValue), k2D_Yield*(1-k2DValue) )/fSigmaPhiValue( k1D_Yield, k2D_Yield ) )  );
    get<4>(fResult)->SetBinContent(5,P1_Sys);
    auto P2_Sys = max( fabs( 1- fGammaPhiValue( k1D_Yield*(1+k1DValue), k2D_Yield*(1+k2DValue) )/fGammaPhiValue( k1D_Yield, k2D_Yield ) ), fabs( 1- fGammaPhiValue( k1D_Yield*(1-k1DValue), k2D_Yield*(1-k2DValue) )/fGammaPhiValue( k1D_Yield, k2D_Yield ) )  );
    get<4>(fResult)->SetBinContent(6,P2_Sys );
    
    get<1>(fResult)     -> SetLineStyle(iLineStyle);
    get<2>(fResult)     -> SetLineStyle(iLineStyle);
    get<3>(fResult)     -> SetLineStyle(iLineStyle);
    get<4>(fResult)     -> SetLineStyle(iLineStyle);
    get<1>(fResult)     -> SetLineColor( uGetColor(iColor) );
    get<2>(fResult)     -> SetLineColor( uGetColor(iColor) );
    get<3>(fResult)     -> SetLineColor( uGetColor(iColor) );
    get<4>(fResult)     -> SetLineColor( uGetColor(iColor) );
    //!
    InputFile   -> Close();
    return fResult;
}
//!
template <  typename TH_1D = TH1F,
            typename TH_2D = TH2F,
            typename TH_RT = TH1F,
            typename TH_PT = TH1F >
std::tuple< TString, TH_1D*, TH_2D*, TH_RT*, TH_PT* >
uLoadExternalSystematics
( TString kFile, std::tuple< TString, TString, TString, TString> kHistogramTags, TString kNewTag, TH1F* k1D, TH2F* k2D, Int_t iColor, Int_t iLineStyle ) {
    //! Book result variable
    std::tuple< TString, TH_1D*, TH_2D*, TH_RT*, TH_PT* > fResult;
    get<0>(fResult)     = kNewTag;
    //!
    TFile  *InputFile   = new TFile( kFile );
    TString kFolder = "_p_p__5TeV";
    TString kCurrentSystematic = "TRK";
    TString kType = "Yield";
    TFile  *InputFil2   = new TFile( Form( kAnalysis_Systemt_Dir + TString("/")+kCurrentSystematic+TString("/FullSystematics.root"), ( kType + kFolder).Data()) );
    gROOT->cd();
    get<1>(fResult)     = uLoadHistograms<0,TH1F> ( InputFile, get<0>(kHistogramTags),  TString("h1D_Syst_") + kNewTag );
    get<2>(fResult)     = uLoadHistograms<0,TH2F> ( InputFile, get<1>(kHistogramTags),  TString("h2D_Syst_") + kNewTag );
    get<3>(fResult)     = uLoadHistograms<0,TH1F> ( InputFil2, "kPTSystematics",    TString("hPT_Syst_") + kNewTag );
    get<4>(fResult)     = uLoadHistograms<0,TH1F> ( InputFil2, "kRatioSystematics", TString("hRT_Syst_") + kNewTag );
    get<1>(fResult)     -> SetLineStyle(iLineStyle);
    get<2>(fResult)     -> SetLineStyle(iLineStyle);
    get<3>(fResult)     -> SetLineStyle(iLineStyle);
    get<4>(fResult)     -> SetLineStyle(iLineStyle);
    get<1>(fResult)     -> SetLineColor( uGetColor(iColor) );
    get<2>(fResult)     -> SetLineColor( uGetColor(iColor) );
    get<3>(fResult)     -> SetLineColor( uGetColor(iColor) );
    get<4>(fResult)     -> SetLineColor( uGetColor(iColor) );
    //!
    InputFile   -> Close();
    return fResult;
}
//!
template <  typename TH_1D = TH1F,
            typename TH_2D = TH2F,
            typename TH_PT = TH1F,
            typename TH_RT = TH1F >
void
uCombineSystematics
( std::vector<std::tuple< TString, TH_1D*, TH_2D*, TH_PT*, TH_RT* >> &kArray_Sys ) {
    //!
    //! Redefine vector
    std::vector<TH_1D*> hSyst1D_Array;
    std::vector<TH_2D*> hSyst2D_Array;
    std::vector<TH_PT*> hSystPT_Array;
    std::vector<TH_RT*> hSystRT_Array;
    for ( auto kCurrentTarget : kArray_Sys ) {
        hSyst1D_Array.push_back(get<1>(kCurrentTarget));
        hSyst2D_Array.push_back(get<2>(kCurrentTarget));
        hSystPT_Array.push_back(get<3>(kCurrentTarget));
        hSystRT_Array.push_back(get<4>(kCurrentTarget));
    }
    //!  --- Combine All systematics
    auto        hSyst1D_ALL     =   uSumSystErrors<true,TH_1D>( hSyst1D_Array );
    auto        hSyst2D_ALL     =   uSumSystErrors<true,TH_2D>( hSyst2D_Array );
    auto        hSystPT_ALL     =   uSumSystErrors<true,TH_PT>( hSystPT_Array );
    auto        hSystRT_ALL     =   uSumSystErrors<true,TH_RT>( hSystRT_Array );
    hSyst1D_ALL ->  SetNameTitle("hFullSystematics1D","");
    hSyst2D_ALL ->  SetNameTitle("hFullSystematics2D","");
    hSystPT_ALL ->  SetNameTitle("hFullSystematicsPT","");
    hSystRT_ALL ->  SetNameTitle("hFullSystematicsRT","");
    //!
    hSyst1D_ALL ->  SetLineStyle(1);
    hSyst2D_ALL ->  SetLineStyle(1);
    hSystPT_ALL ->  SetLineStyle(1);
    hSystRT_ALL ->  SetLineStyle(1);
    hSyst1D_ALL ->  SetLineWidth(3);
    hSyst2D_ALL ->  SetLineWidth(3);
    hSystPT_ALL ->  SetLineWidth(3);
    hSystRT_ALL ->  SetLineWidth(3);
    hSyst1D_ALL ->  SetLineColor( uGetColor(0) );
    hSyst2D_ALL ->  SetLineColor( uGetColor(0) );
    hSystPT_ALL ->  SetLineColor( uGetColor(0) );
    hSystRT_ALL ->  SetLineColor( uGetColor(0) );
    hSyst1D_ALL ->  GetYaxis()  ->  SetTitle("Fractional Uncertainty (%)");
    hSyst2D_ALL ->  GetZaxis()  ->  SetTitle("Fractional Uncertainty (%)");
    hSystPT_ALL ->  GetYaxis()  ->  SetTitle("Fractional Uncertainty (%)");
    hSystRT_ALL ->  GetYaxis()  ->  SetTitle("Fractional Uncertainty (%)");
    //
    push_to_front(kArray_Sys,std::tuple< TString, TH_1D*, TH_2D*, TH_RT*, TH_PT* > ( "Full", hSyst1D_ALL, hSyst2D_ALL, hSystPT_ALL, hSystRT_ALL));
}
//!
template <  typename THTarget = TH1F >
void
uPlotSystCheck
( THTarget* kStatReference, std::vector<THTarget*> kArray_Sys, TLegend* kLegend, TString kSaveName, Bool_t kLogScale ) {
    TCanvas*    cDraw   = new TCanvas("cDraw","cDraw",1000,1000);
    gStyle      -> SetOptStat(0);
    gPad        -> SetTopMargin(0.04);
    gPad        -> SetBottomMargin(0.14);
    gPad        -> SetLeftMargin(0.11);
    gPad        -> SetRightMargin(0.04);
    if (kLogScale) gPad        -> SetLogx();
    kStatReference = uMakeRelativeUncertainty(kStatReference,true);
    kStatReference->SetFillColorAlpha(kGray+1,0.33);
    kStatReference->SetLineColorAlpha(kWhite,0.00);
    for ( auto kCurrentSystematics : kArray_Sys ) kCurrentSystematics->Scale(100.);
    for ( auto kCurrentSystematics : kArray_Sys ) kCurrentSystematics->GetYaxis()->SetTitleOffset(1.4);
    for ( auto kCurrentSystematics : kArray_Sys ) kCurrentSystematics->GetYaxis()->SetTitleSize(0.04);
    kArray_Sys.at(0)->SetMaximum(max(1.3*kArray_Sys.at(0)->GetMaximum(),1.3*kStatReference->GetMaximum()));
    kArray_Sys.at(0)->SetMinimum(0);
    kArray_Sys.at(0)->SetLineWidth(3);
    for ( auto kCurrentSystematics : kArray_Sys ) kCurrentSystematics->Draw("SAME HIST");
    kStatReference->Draw("SAME HIST");
    kLegend->Draw("SAME HIST");
    cDraw->SaveAs(kSaveName);
}
//!
template <  typename TH_1D = TH1F,
            typename TH_2D = TH2F,
            typename TH_PT = TH1F,
            typename TH_RT = TH1F >
void
uPlotFullSystematics
 ( std::tuple< TString, TH_1D*, TH_2D*, TH_PT*, TH_RT* > kStatReference, std::vector<std::tuple< TString, TH_1D*, TH_2D*, TH_PT*, TH_RT* >> kArray_Sys, TString kFolder, TString kType = "Yield", TString kSubFolder = "", Bool_t kLogScale = true ) {
    //! Building output directory
    gROOT->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_Systemt_Dir + kSubFolder + TString("/plots/"), (kType+kFolder).Data()) ));
    //! Build legend
    get<1>(kStatReference)->SetFillColorAlpha(kGray+1,0.33);
    get<1>(kStatReference)->SetLineColorAlpha(kWhite,0.00);
    TLegend*    cLegendUncertainties    =   new TLegend(0.18,0.95,0.88,0.80);
    cLegendUncertainties    ->  SetNColumns(4);
    cLegendUncertainties    ->  SetFillColorAlpha(kWhite,0.0);
    cLegendUncertainties    ->  SetLineColorAlpha(kWhite,0.0);
    cLegendUncertainties    ->  AddEntry( get<1>(kStatReference), get<0>(kStatReference), "F" );
    for ( auto kCurrentSystematics : kArray_Sys ) cLegendUncertainties    ->  AddEntry( get<1>(kCurrentSystematics), get<0>(kCurrentSystematics), "L" );
    //!
    std::vector<TH_1D*> hSyst1D_Array;
    std::vector<TH_2D*> hSyst2D_Array;
    std::vector<TH_PT*> hSystPT_Array;
    std::vector<TH_RT*> hSystRT_Array;
    for ( auto kCurrentTarget : kArray_Sys ) {
        hSyst1D_Array.push_back(get<1>(kCurrentTarget));
        hSyst2D_Array.push_back(get<2>(kCurrentTarget));
        hSystPT_Array.push_back(get<3>(kCurrentTarget));
        hSystRT_Array.push_back(get<4>(kCurrentTarget));
    }
    //! Build canvases
    //!
    //!     1D
    uPlotSystCheck( get<1>(kStatReference), hSyst1D_Array, cLegendUncertainties,(Form(kAnalysis_Systemt_Dir + kSubFolder + TString("/plots/hSyst1D.pdf"), (kType+kFolder).Data())), true );
    //!
    //!     2D
    auto iTer = -1;
    for ( auto iXproj = 0; iXproj < get<2>(kStatReference)->GetNbinsX(); iXproj++ ) {
        iTer++;
        auto kCurrentSlice_Stat = get<2>(kStatReference)->ProjectionX(TString(get<2>(kStatReference)->GetName())+TString(Form("%i",iXproj)),iTer+1,iTer+1);
        std::vector<TH1D*> kSystSlices;
        for ( auto kCurrent_Sys : kArray_Sys ) {
            kSystSlices.push_back( get<2>(kCurrent_Sys)->ProjectionX(TString(get<2>(kCurrent_Sys)->GetName())+TString(Form("%i",iXproj)),iTer+1,iTer+1) );
        }
        uPlotSystCheck( kCurrentSlice_Stat, kSystSlices, cLegendUncertainties,(Form(kAnalysis_Systemt_Dir + kSubFolder + TString(Form("/plots/hSyst2D_%i.pdf",iTer)), (kType+kFolder).Data())), true );
    }
    //!
    //!     PT
    //!
    //!     YIELDS
    uPlotSystCheck( get<4>(kStatReference), hSystRT_Array, cLegendUncertainties,(Form(kAnalysis_Systemt_Dir + kSubFolder + TString("/plots/hSystRT.pdf"), (kType+kFolder).Data())), false );
}
//!
void
SY_AN_All
( TString fOption = "Yield mult", TString kFolder = "_p_p__5TeV" )    {
    //!
    //! Set-Up
    //! --- Option chosing
    if ( !fChooseOption(fOption) ) return;
    //! --- Style
    SetStyle();
    //! --- Binning
    fSetAllBins();
    //!
    //! Yield analysis
    if ( kDoYield ) {
        //!
        //! Building output directory
        gROOT->  ProcessLine(Form(".! mkdir -p %s",Form(kSystematicsPlot,       (TString("Yield")+kFolder).Data())));
        //!
        //! Retrieve Statistics
        TFile*  insFile_Data_YL     =   new TFile   (Form(kASigExtp_FitCheckRst,(TString("Yield")+kFolder).Data()));
        gROOT->cd();
        auto        hSysU1D_STD     =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "h1D_Nres_stat",        "hSysU1D_STD" );
        auto        hSysU2D_STD     =   uLoadHistograms<0,TH2F> ( insFile_Data_YL,  "h2D_Nres_stat",        "hSysU2D_STD" );
        auto        hSysUPT_STD     =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "h2D_MeanPT_stat",      "hSysUPT_STD" );
        auto        hSysURT_STD     =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "hXD_Nfqs_stat",        "hSysURT_STD" );
        auto        hSysEXT_YLD     =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "hxD_YL_Extrapol_Syst", "hSysEXT_YLD" );
        auto        hSysEXT_MPT     =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "hxD_PT_Extrapol_Syst", "hSysEXT_MPT" );
        auto k1D_Yield = hSysURT_STD->GetBinContent(1);
        auto k2D_Yield = hSysURT_STD->GetBinContent(2);
        //!
        //! Retrieve All Systematics
        std::vector<std::tuple< TString, TH1F*, TH2F*, TH1F*, TH1F* >> kAllSysArray;
        //! 
        kAllSysArray.push_back(uLoadCustomSystematics   (kFolder, "Yield", "SEX",   "Sig. Extr.",   1,1));
        kAllSysArray.push_back(uLoadCustomSystematics   (kFolder, "Yield", "PID",   "PID",          2,1));
        kAllSysArray.push_back(uLoadCustomSystematics   (kFolder, "Yield", "TRK",   "Track Cuts",   3,1));
        kAllSysArray.push_back(uGenerateFlatSystematics (kFolder, "Yield", "TRK",   "Br. Ratio",    0.01,0.02,0.,k1D_Yield,k2D_Yield,4,1));
        kAllSysArray.push_back(uGenerateFlatSystematics (kFolder, "Yield", "TRK",   "Normalisation",0.025,0.025,0.,k1D_Yield,k2D_Yield,5,1));
        kAllSysArray.push_back(uLoadCustomSystematics   (kFolder, "Yield", "TKN",   "Tracking",     6,1));
        kAllSysArray.push_back(uLoadCustomSystematics   (kFolder, "Yield", "HIN",   "Had. Int.",    7,1));
        kAllSysArray.push_back(uLoadCustomSystematics   (kFolder, "Yield", "MBD",   "Mat. Bud.",    8,1));
        kAllSysArray.push_back(uLoadCustomSystematics   (kFolder, "Yield", "EXT",   "Extrap.",      9,1));
        uCombineSystematics( kAllSysArray );
        uPlotFullSystematics( {"Stat",hSysU1D_STD,hSysU2D_STD,hSysUPT_STD,hSysURT_STD}, kAllSysArray, kFolder);
        //!
        TFile*  OutputFile    =   new TFile(Form("%s/FullSystematics.root",Form(kAnalysis_Systemt_Dir,  (TString("Yield")+kFolder).Data())),"RECREATE");
        get<1>(kAllSysArray.at(0))->Scale(1./100);
        //get<2>(kAllSysArray.at(0))->Scale(1./100);
        get<3>(kAllSysArray.at(0))->Scale(1./100);
        get<4>(kAllSysArray.at(0))->Scale(1./100);
        get<1>(kAllSysArray.at(0))->Write();
        get<2>(kAllSysArray.at(0))->Write();
        get<3>(kAllSysArray.at(0))->Write();
        get<4>(kAllSysArray.at(0))->Write();
        //!
        OutputFile->Close();
    }
    if ( kDoMultiplicity ) {
        //!
        //! Building output directory
        gROOT->  ProcessLine(Form(".! mkdir -p %s",Form(kSystematicsPlot, (TString("Multiplicity")+kFolder).Data())));
        //!
        //! Retrieve Statistics
        TFile*  insFile_Data_MT     =   new TFile   (Form(kASigExtp_FitCheckRst,(TString("Multiplicity")+kFolder).Data()));
        gROOT->cd();
        auto    hSysU1D_STD         = uLoadHistograms<1,TH1F> ( insFile_Data_MT,  "h1D_Nres_stat_MT_%i", "hSysU1D_STD_MT_%i" );
        auto    hSysU2D_STD         = uLoadHistograms<1,TH2F> ( insFile_Data_MT,  "h2D_Nres_stat_MT_%i", "hSysU2D_STD_MT_%i" );
        auto    g1D_Nres_stat_Mult  = uLoadHistograms<0,TGraphErrors> ( insFile_Data_MT,    "g1D_Nres_Stat_Mult",          "g1D_Nres_stat_Mult" );
        auto    g2D_Nres_stat_Mult  = uLoadHistograms<0,TGraphErrors> ( insFile_Data_MT,    "g2D_Nres_Stat_Mult",          "g2D_Nres_stat_Mult" );
        //auto        hSysUPT_STD     =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "h2D_MeanPT_stat",      "hSysUPT_STD" );
        auto        hSysURT_STD     =   uLoadHistograms<1,TH1F> ( insFile_Data_MT,  "hXD_Nfqs_stat_MT_%i",        "hSysURT_STD" );
        //auto        hSysEXT_YLD     =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "hxD_YL_Extrapol_Syst", "hSysEXT_YLD" );
        //auto        hSysEXT_MPT     =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "hxD_PT_Extrapol_Syst", "hSysEXT_MPT" );
        //!
        for ( Int_t iMlt = 0; iMlt <= nBinMult; iMlt++ ) {
            auto k1D_Yield = hSysURT_STD.at(iMlt)->GetBinContent(1);
            auto k2D_Yield = hSysURT_STD.at(iMlt) ->GetBinContent(2);
            //!
            //! Retrieve All Systematics
            std::vector<std::tuple< TString, TH1F*, TH2F*, TH1F*, TH1F* >> kAllSysArray;
            kAllSysArray.push_back(uLoadCustomSystematics   (kFolder, "Multiplicity", Form("/MLT_%i/SEX",iMlt),    "Sig. Extr.",   1,1));
            kAllSysArray.push_back(uLoadCustomSystematics   (kFolder, "Multiplicity", Form("/MLT_%i/PID/",iMlt),                "PID",          2,1));
            kAllSysArray.push_back(uLoadCustomSystematics   (kFolder, "Multiplicity", Form("/MLT_%i/TRK/",iMlt),                "Track Cuts",   3,1));
            kAllSysArray.push_back(uGenerateFlatSystematics (kFolder, "Yield", "TRK", "Br. Ratio",    0.01,0.02,0.,k1D_Yield,k2D_Yield,4,1));
            kAllSysArray.push_back(uGenerateFlatSystematics (kFolder, "Yield", "TRK", "Normalisation",0.025,0.025,0.,k1D_Yield,k2D_Yield,5,1));
            kAllSysArray.push_back(uLoadCustomSystematics   (kFolder, "Yield", "TKN", "Tracking",    6,1));
            kAllSysArray.push_back(uLoadCustomSystematics   (kFolder, "Yield", "HIN", "Had. Int.",   7,1));
            kAllSysArray.push_back(uLoadCustomSystematics   (kFolder, "Yield", "MBD", "Mat. Bud.",   8,1));
            uCombineSystematics( kAllSysArray );
            uPlotFullSystematics( {"Stat",hSysU1D_STD.at(iMlt),hSysU2D_STD.at(iMlt),hSysU1D_STD.at(iMlt),hSysURT_STD.at(iMlt)}, kAllSysArray, kFolder, "Multiplicity", Form("/MLT_%i/",iMlt) );
            //!
            TFile*  OutputFile    =   new TFile(Form("%s/FullSystematics.root",Form(kAnalysis_Systemt_Dir+TString(Form("/MLT_%i/",iMlt)),(TString("Multiplicity")+kFolder).Data())),"RECREATE");
            get<1>(kAllSysArray.at(0))->Scale(1./100);
            //get<2>(kAllSysArray.at(0))->Scale(1./100);
            get<3>(kAllSysArray.at(0))->Scale(1./100);
            get<4>(kAllSysArray.at(0))->Scale(1./100);
            get<1>(kAllSysArray.at(0))->Write();
            get<2>(kAllSysArray.at(0))->Write();
            get<3>(kAllSysArray.at(0))->Write();
            get<4>(kAllSysArray.at(0))->Write();
            //!
            OutputFile->Close();
        }
    }
}


/*
 
 //
 //  --- --- Signal Extraction
 TFile      *fIn_Syst_SEX    =   new TFile( Form( kAnalysis_Systemt_Dir + TString("/SignalExtraction/FullSystematics.root"),  (TString("Yield") + kFolder).Data()) );
 auto        hSyst1D_SEX     =   uLoadHistograms<0,TH1F> ( fIn_Syst_SEX,  "k1DSystematics",      "hSyst1D_SEX" );
 auto        hSyst2D_SEX     =   uLoadHistograms<0,TH2F> ( fIn_Syst_SEX,  "k2DSystematics",      "hSyst2D_SEX" );
 auto        hSysUPT_SEX     =   uLoadHistograms<0,TH1F> ( fIn_Syst_SEX,  "kPTSystematics",      "hSysUPT_SEX" );
 auto        hSystRT_SEX     =   uLoadHistograms<0,TH1F> ( fIn_Syst_SEX,  "kRatioSystematics",   "hSystRT_SEX" );
 auto        hSystPT_SEX     =   hSyst2D_SEX->ProjectionX("hSystPT_SEX",1,1);
 //
 hSyst1D_SEX ->  Scale(100);
 hSyst2D_SEX ->  Scale(100);
 hSystRT_SEX ->  Scale(100);
 hSystPT_SEX ->  Scale(100);
 hSyst1D_SEX ->  SetLineStyle(1);
 hSyst2D_SEX ->  SetLineStyle(1);
 hSystRT_SEX ->  SetLineStyle(1);
 hSystPT_SEX ->  SetLineStyle(1);
 hSyst1D_SEX ->  SetLineColor( uGetColor(1) );
 hSyst2D_SEX ->  SetLineColor( uGetColor(1) );
 hSystRT_SEX ->  SetLineColor( uGetColor(1) );
 hSystPT_SEX ->  SetLineColor( uGetColor(1) );
 //
 for ( Int_t iTer = 1; iTer <= hSystPT_SEX->GetNbinsX(); iTer++ )  {
     hSystPT_SEX ->  SetBinContent   ( iTer, (hSysUPT_SEX->GetBinContent(iTer+2)) );
     hSystPT_SEX ->  SetBinError     ( iTer, 0. );
 }
 //
 //  --- --- PID
 TFile      *fIn_Syst_PID    =   new TFile( Form( kAnalysis_Systemt_Dir + TString("/PID/FullSystematics.root"),  (TString("Yield") + kFolder).Data()) );
 auto        hSyst1D_PID     =   uLoadHistograms<0,TH1F> ( fIn_Syst_PID,  "k1DSystematics",      "hSyst1D_PID" );
 auto        hSyst2D_PID     =   uLoadHistograms<0,TH2F> ( fIn_Syst_PID,  "k2DSystematics",      "hSyst2D_PID" );
 auto        hSysUPT_PID     =   uLoadHistograms<0,TH1F> ( fIn_Syst_PID,  "kPTSystematics",      "hSysUPT_PID" );
 auto        hSystRT_PID     =   uLoadHistograms<0,TH1F> ( fIn_Syst_PID,  "kRatioSystematics",   "hSystRT_PID" );
 auto        hSystPT_PID     =   hSyst2D_PID->ProjectionX("hSystPT_PID",1,1);
 //
 hSyst1D_PID ->  Scale(100);
 hSyst2D_PID ->  Scale(100);
 hSystPT_PID ->  Scale(100);
 hSystRT_PID ->  Scale(100);
 hSyst1D_PID ->  SetLineStyle(1);
 hSyst2D_PID ->  SetLineStyle(1);
 hSystPT_PID ->  SetLineStyle(1);
 hSystRT_PID ->  SetLineStyle(1);
 hSyst1D_PID ->  SetLineColor( uGetColor(2) );
 hSyst2D_PID ->  SetLineColor( uGetColor(2) );
 hSystPT_PID ->  SetLineColor( uGetColor(2) );
 hSystRT_PID ->  SetLineColor( uGetColor(2) );
 //
 for ( Int_t iTer = 1; iTer <= hSystPT_PID->GetNbinsX(); iTer++ )  {
     hSystPT_PID ->  SetBinContent   ( iTer, (hSysUPT_PID->GetBinContent(iTer+2)) );
     hSystPT_PID ->  SetBinError     ( iTer, 0. );
 }
 //
 //  --- --- Track Selection
 TFile      *fIn_Syst_TRK    =   new TFile( Form( kAnalysis_Systemt_Dir + TString("/TRK/FullSystematics.root"),  (TString("Yield") + kFolder).Data()) );
 auto        hSyst1D_TRK     =   uLoadHistograms<0,TH1F> ( fIn_Syst_TRK,  "k1DSystematics",      "hSyst1D_TRK" );
 auto        hSyst2D_TRK     =   uLoadHistograms<0,TH2F> ( fIn_Syst_TRK,  "k2DSystematics",      "hSyst2D_TRK" );
 auto        hSysUPT_TRK     =   uLoadHistograms<0,TH1F> ( fIn_Syst_TRK,  "kPTSystematics",      "hSysUPT_TRK" );
 auto        hSystRT_TRK     =   uLoadHistograms<0,TH1F> ( fIn_Syst_TRK,  "kRatioSystematics",   "hSystRT_TRK" );
 auto        hSystPT_TRK     =   hSyst2D_TRK->ProjectionX("hSystPT_TRK",1,1);
 //
 hSyst1D_TRK ->  Scale(100);
 hSyst2D_TRK ->  Scale(100);
 hSystPT_TRK ->  Scale(100);
 hSystRT_TRK ->  Scale(100);
 hSyst1D_TRK ->  SetLineStyle(1);
 hSyst2D_TRK ->  SetLineStyle(1);
 hSystPT_TRK ->  SetLineStyle(1);
 hSystRT_TRK ->  SetLineStyle(1);
 hSyst1D_TRK ->  SetLineColor( uGetColor(3) );
 hSyst2D_TRK ->  SetLineColor( uGetColor(3) );
 hSystPT_TRK ->  SetLineColor( uGetColor(3) );
 hSystRT_TRK ->  SetLineColor( uGetColor(3) );
 //
 for ( Int_t iTer = 1; iTer <= hSystPT_TRK->GetNbinsX(); iTer++ )  {
     hSystPT_TRK ->  SetBinContent   ( iTer, (hSysUPT_TRK->GetBinContent(iTer+2)) );
     hSystPT_TRK ->  SetBinError     ( iTer, 0. );
 }
 //
 //  --- Calculate All A Priori Systematics
 //
 //  --- --- Utility Yields
 auto    k1D_Yield   =   hSysURT_STD ->  GetBinContent   (1);
 auto    k2D_Yield   =   hSysURT_STD ->  GetBinContent   (2);
 auto    k1D_Error   =   hSysURT_STD ->  GetBinError     (1);
 auto    k2D_Error   =   hSysURT_STD ->  GetBinError     (2);
 //
 //  --- --- Branching Ratio
 auto        hSyst1D_BrR     =   uLoadHistograms<0,TH1F> ( fIn_Syst_SEX,  "k1DSystematics",      "hSyst1D_BrR" );
 auto        hSyst2D_BrR     =   uLoadHistograms<0,TH2F> ( fIn_Syst_SEX,  "k2DSystematics",      "hSyst2D_BrR" );
 auto        hSystRT_BrR     =   uLoadHistograms<0,TH1F> ( fIn_Syst_SEX,  "kRatioSystematics",   "hSystRT_BrR" );
 for ( Int_t iBin = 1; iBin <= hSyst1D_BrR->GetNbinsX(); iBin++ )  hSyst1D_BrR->SetBinContent( iBin, kSysHig_BR );
 for ( Int_t iBin = 1; iBin <= hSyst2D_BrR->GetNbinsX(); iBin++ ) for ( Int_t jBin = 1; jBin <= hSyst2D_BrR->GetNbinsY(); jBin++ ) hSyst2D_BrR->SetBinContent( iBin, jBin, kSysHig_BR*2  );
 hSystRT_BrR ->  SetBinContent( 1, kSysHig_BR );
 hSystRT_BrR ->  SetBinContent( 2, kSysHig_BR*kSysHig_BR );
 hSystRT_BrR ->  SetBinContent( 3, kSysHig_BR );
 hSystRT_BrR ->  SetBinContent( 4, 0. );
 hSystRT_BrR ->  SetBinContent( 5, fabs( 1- fSigmaPhiValue( k1D_Yield*(1+kSysHig_BR), k2D_Yield*(1+kSysHig_BR) )/fSigmaPhiValue( k1D_Yield, k2D_Yield ) ) );
 hSystRT_BrR ->  SetBinContent( 6, fabs( 1- fGammaPhiValue( k1D_Yield*(1+kSysHig_BR), k2D_Yield*(1+kSysHig_BR) )/fGammaPhiValue( k1D_Yield, k2D_Yield ) ) );
 //
 hSyst1D_BrR ->  Scale(100);
 hSyst2D_BrR ->  Scale(100);
 hSystRT_BrR ->  Scale(100);
 hSyst1D_BrR ->  SetLineStyle(1);
 hSyst2D_BrR ->  SetLineStyle(1);
 hSystRT_BrR ->  SetLineStyle(1);
 hSyst1D_BrR ->  SetLineColor( uGetColor(4) );
 hSyst2D_BrR ->  SetLineColor( uGetColor(4) );
 hSystRT_BrR ->  SetLineColor( uGetColor(4) );
 //
 //  --- --- Tracking Efficiency
 auto        hSyst1D_ITP     =   uLoadHistograms<0,TH1F> ( fIn_Syst_SEX,  "k1DSystematics",      "hSyst1D_ITP" );
 auto        hSyst2D_ITP     =   uLoadHistograms<0,TH2F> ( fIn_Syst_SEX,  "k2DSystematics",      "hSyst2D_ITP" );
 auto        hSystRT_ITP     =   uLoadHistograms<0,TH1F> ( fIn_Syst_SEX,  "kRatioSystematics",   "hSystRT_ITP" );
 for ( Int_t iBin = 1; iBin <= hSyst1D_ITP->GetNbinsX(); iBin++ )  hSyst1D_ITP->SetBinContent( iBin, kTRerr );
 for ( Int_t iBin = 1; iBin <= hSyst2D_ITP->GetNbinsX(); iBin++ ) for ( Int_t jBin = 1; jBin <= hSyst2D_ITP->GetNbinsY(); jBin++ ) hSyst2D_ITP->SetBinContent( iBin, jBin, kTRerr*2 );
 hSystRT_ITP ->  SetBinContent( 1, kTRerr );
 hSystRT_ITP ->  SetBinContent( 2, kTRerr*2 );
 hSystRT_ITP ->  SetBinContent( 3, kTRerr );
 hSystRT_ITP ->  SetBinContent( 4, 0. );
 hSystRT_ITP ->  SetBinContent( 5, fabs( 1- fSigmaPhiValue( k1D_Yield*(1+kTRerr), k2D_Yield*(1+kTRerr) )/fSigmaPhiValue( k1D_Yield, k2D_Yield ) ) );
 hSystRT_ITP ->  SetBinContent( 6, fabs( 1- fGammaPhiValue( k1D_Yield*(1+kTRerr), k2D_Yield*(1+kTRerr) )/fGammaPhiValue( k1D_Yield, k2D_Yield ) ) );
 //
 hSyst1D_ITP ->  Scale(100);
 hSyst2D_ITP ->  Scale(100);
 hSystRT_ITP ->  Scale(100);
 hSyst1D_ITP ->  SetLineStyle(1);
 hSyst2D_ITP ->  SetLineStyle(1);
 hSystRT_ITP ->  SetLineStyle(1);
 hSyst1D_ITP ->  SetLineColor( uGetColor(5) );
 hSyst2D_ITP ->  SetLineColor( uGetColor(5) );
 hSystRT_ITP ->  SetLineColor( uGetColor(5) );
 //
 //  --- --- Extrapolation Uncertainty
 auto        hSystRT_EXT     =   uLoadHistograms<0,TH1F> ( fIn_Syst_SEX,     "kRatioSystematics",    "hSystRT_EXT" );
 auto        hSystPT_EXT     =   hSyst2D_ITP->ProjectionX("hSystPT_EXT",1,1);
 //
 auto        kSystErrYL1D    =   1.;
 auto        kSystIntYL1D    =   hSysEXT_YLD ->  IntegralAndError( 1,  1, kSystErrYL1D, "width" );
 auto        kSystErrYL2D    =   1.;
 auto        kSystIntYL2D    =   hSysEXT_YLD ->  IntegralAndError( 3, -1, kSystErrYL2D, "width" );
 kSystErrYL1D     /=   kSystIntYL1D;
 kSystErrYL2D     /=   kSystIntYL2D;
 //
 hSystRT_EXT ->  SetBinContent( 1, kSystErrYL1D );
 hSystRT_EXT ->  SetBinContent( 2, kSystErrYL2D );
 hSystRT_EXT ->  SetBinContent( 3, SquareSum( {kSystErrYL2D,kSystErrYL1D} ) );
 hSystRT_EXT ->  SetBinContent( 4, SquareSum( {kSystErrYL2D,kSystErrYL1D,kSystErrYL1D} ) );
 hSystRT_EXT ->  SetBinContent( 5, fSigmaPhiError(k1D_Yield,k2D_Yield,k1D_Yield*kSystErrYL1D,k2D_Yield*kSystErrYL2D) / fSigmaPhiValue(k1D_Yield,k2D_Yield) );
 hSystRT_EXT ->  SetBinContent( 6, fGammaPhiError(k1D_Yield,k2D_Yield,k1D_Yield*kSystErrYL1D,k2D_Yield*kSystErrYL2D) / fGammaPhiValue(k1D_Yield,k2D_Yield) );
 //
 hSystRT_EXT ->  Scale(100);
 hSystRT_EXT ->  SetLineStyle(1);
 hSystRT_EXT ->  SetLineColor( uGetColor(6) );
 //
 for ( Int_t iTer = 1; iTer <= hSystPT_EXT->GetNbinsX(); iTer++ )  {
     hSystPT_EXT ->  SetBinContent   ( iTer, hSysEXT_MPT->GetBinError(iTer+2)/(hSysEXT_MPT->GetBinContent(iTer+2)) );
     hSystPT_EXT ->  SetBinError     ( iTer, 0. );
 }
 //
 hSystPT_EXT ->  Scale(100);
 hSystPT_EXT ->  SetLineStyle(1);
 hSystPT_EXT ->  SetLineColor( uGetColor(6) );
 //
 //  --- --- Normalisation Uncertainty
 auto        hSystRT_NMP     =   uLoadHistograms<0,TH1F> ( fIn_Syst_SEX,     "kRatioSystematics",    "hSystRT_NMP" );
 auto        hSystRT_NMM     =   uLoadHistograms<0,TH1F> ( fIn_Syst_SEX,     "kRatioSystematics",    "hSystRT_NMM" );
 //
 hSystRT_NMP ->  SetBinContent( 1, +kSysHig_TE );
 hSystRT_NMP ->  SetBinContent( 2, +kSysHig_TE );
 hSystRT_NMP ->  SetBinContent( 3, 0. );
 hSystRT_NMP ->  SetBinContent( 4, +kSysHig_TE );
 hSystRT_NMP ->  SetBinContent( 5, +fabs(1-fSigmaPhiValue(k1D_Yield*(1+kSysHig_TE),k2D_Yield*(1+kSysHig_TE)) / fSigmaPhiValue(k1D_Yield,k2D_Yield)) );
 hSystRT_NMP ->  SetBinContent( 6, +fabs(1-fGammaPhiValue(k1D_Yield*(1+kSysHig_TE),k2D_Yield*(1+kSysHig_TE)) / fGammaPhiValue(k1D_Yield,k2D_Yield)) );
 hSystRT_NMM ->  SetBinContent( 1, -kSysLow_TE  );
 hSystRT_NMM ->  SetBinContent( 2, -kSysLow_TE  );
 hSystRT_NMM ->  SetBinContent( 3, 0. );
 hSystRT_NMM ->  SetBinContent( 4, -kSysLow_TE  );
 hSystRT_NMM ->  SetBinContent( 5, -fabs(1-fSigmaPhiValue(k1D_Yield*(1-kSysLow_TE),k2D_Yield*(1-kSysLow_TE)) / fSigmaPhiValue(k1D_Yield,k2D_Yield)) );
 hSystRT_NMM ->  SetBinContent( 6, -fabs(1-fGammaPhiValue(k1D_Yield*(1-kSysLow_TE),k2D_Yield*(1-kSysLow_TE)) / fGammaPhiValue(k1D_Yield,k2D_Yield)) );
 //
 hSystRT_NMP ->  Scale(100);
 hSystRT_NMP ->  SetLineColorAlpha( 0., 0. );
 hSystRT_NMP ->  SetFillColorAlpha( uGetColor(7), 0.33 );
 hSystRT_NMM ->  Scale(100);
 hSystRT_NMM ->  SetLineColorAlpha( 0., 0. );
 hSystRT_NMM ->  SetFillColorAlpha( uGetColor(7), 0.33 );
 //
 //  --- --- Reference Statistical
 for ( Int_t iBin = 1; iBin <= hSysU1D_STD->GetNbinsX(); iBin++  )    {
     auto    kBinValue   =   hSysU1D_STD->GetBinContent  (iBin);
     auto    kBinError   =   hSysU1D_STD->GetBinError    (iBin);
     hSysU1D_STD->SetBinContent  ( iBin, 100*kBinError/kBinValue );
     hSysU1D_STD->SetBinError    ( iBin, 0. );
 }
 hSysU1D_STD->SetFillColorAlpha(kGray+1,0.33);
 hSysU1D_STD->SetLineColorAlpha(kWhite,0.00);
 for ( Int_t iBin = 1; iBin <= hSysUPT_STD->GetNbinsX(); iBin++  )    {
     auto    kBinValue   =   hSysUPT_STD->GetBinContent  (iBin);
     auto    kBinError   =   hSysUPT_STD->GetBinError    (iBin);
     hSysUPT_STD->SetBinContent  ( iBin, 100*kBinError/kBinValue );
     hSysUPT_STD->SetBinError    ( iBin, 0. );
 }
 hSysUPT_STD->SetFillColorAlpha(kGray+1,0.33);
 hSysUPT_STD->SetLineColorAlpha(kWhite,0.00);
 for ( Int_t iBin = 1; iBin <= hSysU2D_STD->GetNbinsX(); iBin++  )    {
     for ( Int_t jBin = 1; jBin <= hSysU2D_STD->GetNbinsY(); jBin++  )    {
         auto    kBinValue   =   hSysU2D_STD->GetBinContent  (iBin,jBin);
         auto    kBinError   =   hSysU2D_STD->GetBinError    (iBin,jBin);
         hSysU2D_STD->SetBinContent  ( iBin,jBin, 100*kBinError/kBinValue );
         hSysU2D_STD->SetBinError    ( iBin,jBin, 0. );
     }
 }
 hSysU2D_STD->SetFillColorAlpha(kGray+1,0.33);
 hSysU2D_STD->SetLineColorAlpha(kWhite,0.00);
 auto    hSystRT_STD =   (TH1F*)(hSystRT_ITP->Clone());
 hSystRT_STD->SetBinContent  ( 1, 100*k1D_Error/k1D_Yield );
 hSystRT_STD->SetBinError    ( 1, 0. );
 hSystRT_STD->SetBinContent  ( 2, 100*k2D_Error/k2D_Yield );
 hSystRT_STD->SetBinError    ( 2, 0. );
 hSystRT_STD->SetBinContent  ( 3, 100*SquareSum( {k2D_Error/k2D_Yield,k1D_Error/k1D_Yield} ) );
 hSystRT_STD->SetBinError    ( 3, 0. );
 hSystRT_STD->SetBinContent  ( 4, 100*SquareSum( {k2D_Error/k2D_Yield,k1D_Error/k1D_Yield,k1D_Error/k1D_Yield} ) );
 hSystRT_STD->SetBinError    ( 4, 0. );
 hSystRT_STD->SetBinContent  ( 5, 100*(fSigmaPhiError(k1D_Yield,k2D_Yield,k1D_Error,k2D_Error) / fSigmaPhiValue(k1D_Yield,k2D_Yield)) );
 hSystRT_STD->SetBinError    ( 5, 0. );
 hSystRT_STD->SetBinContent  ( 6, 100*(fGammaPhiError(k1D_Yield,k2D_Yield,k1D_Error,k2D_Error) / fGammaPhiValue(k1D_Yield,k2D_Yield)) );
 hSystRT_STD->SetBinError    ( 6, 0. );
 hSystRT_STD->SetFillColorAlpha(kGray+1,0.33);
 hSystRT_STD->SetLineColorAlpha(kWhite,0.00);
 //
 //  --- Combine All systematics
 auto        hSyst1D_ALL     =   uSumSystErrors<true,TH1F>( {hSyst1D_SEX,hSyst1D_PID,hSyst1D_TRK,hSyst1D_BrR,hSyst1D_ITP} );
 auto        hSyst2D_ALL     =   uSumSystErrors<true,TH2F>( {hSyst2D_SEX,hSyst2D_PID,hSyst2D_TRK,hSyst2D_BrR,hSyst2D_ITP} );
 auto        hSystPT_ALL     =   uSumSystErrors<true,TH1D>( {hSystPT_SEX,hSystPT_PID,hSystPT_TRK,hSystPT_EXT} );
 auto        hSystRT_ALL     =   uSumSystErrors<true,TH1F>( {hSystRT_SEX,hSystRT_PID,hSystRT_TRK,hSystRT_BrR,hSystRT_ITP,hSystRT_EXT} );
 hSyst1D_ALL ->  SetNameTitle("hFullSystematics1D","");
 hSyst2D_ALL ->  SetNameTitle("hFullSystematics2D","");
 hSystPT_ALL ->  SetNameTitle("hFullSystematicsPT","");
 hSystRT_ALL ->  SetNameTitle("hFullSystematicsRT","");
 //
 hSyst1D_ALL ->  SetMaximum( 1.3* max( hSyst1D_ALL->GetMaximum(), hSysU1D_STD->GetMaximum() ) );
 hSyst2D_ALL ->  SetMaximum( 1.3* max( hSyst2D_ALL->GetMaximum(), hSysU2D_STD->GetMaximum() ) );
 hSystPT_ALL ->  SetMaximum( 1.3* max( hSystPT_ALL->GetMaximum(), hSysUPT_STD->GetMaximum() ) );
 hSystRT_ALL ->  SetMaximum( 1.3* max( hSystRT_ALL->GetMaximum(), hSystRT_STD->GetMaximum() ) );
 hSyst1D_ALL ->  SetMinimum(0);
 hSyst2D_ALL ->  SetMinimum(0);
 hSystPT_ALL ->  SetMinimum(0);
 hSystRT_ALL ->  SetMinimum( -5. );
 hSyst1D_ALL ->  SetLineStyle(1);
 hSyst2D_ALL ->  SetLineStyle(1);
 hSystPT_ALL ->  SetLineStyle(1);
 hSystRT_ALL ->  SetLineStyle(1);
 hSyst1D_ALL ->  SetLineWidth(3);
 hSyst2D_ALL ->  SetLineWidth(3);
 hSystPT_ALL ->  SetLineWidth(3);
 hSystRT_ALL ->  SetLineWidth(3);
 hSyst1D_ALL ->  SetLineColor( uGetColor(0) );
 hSyst2D_ALL ->  SetLineColor( uGetColor(0) );
 hSystPT_ALL ->  SetLineColor( uGetColor(0) );
 hSystRT_ALL ->  SetLineColor( uGetColor(0) );
 hSyst1D_ALL ->  GetYaxis()  ->  SetTitle("Fractional Uncertainty (%)");
 hSyst2D_ALL ->  GetZaxis()  ->  SetTitle("Fractional Uncertainty (%)");
 hSystPT_ALL ->  GetYaxis()  ->  SetTitle("Fractional Uncertainty (%)");
 hSystRT_ALL ->  GetYaxis()  ->  SetTitle("Fractional Uncertainty (%)");
 //
 TLegend*    cLegendUncertainties    =   new TLegend(0.18,0.88,0.88,0.75);
 cLegendUncertainties    ->  SetNColumns(4);
 cLegendUncertainties    ->  SetFillColorAlpha(kWhite,0.0);
 cLegendUncertainties    ->  SetLineColorAlpha(kWhite,0.0);
 cLegendUncertainties    ->  AddEntry( hSysU1D_STD, "Stat", "F" );
 cLegendUncertainties    ->  AddEntry( hSyst1D_ALL, "Full (No Norm)", "L" );
 cLegendUncertainties    ->  AddEntry( hSyst1D_SEX, "Sig. Extr", "L" );
 cLegendUncertainties    ->  AddEntry( hSyst1D_PID, "PID", "L" );
 cLegendUncertainties    ->  AddEntry( hSyst1D_TRK, "Track Cuts", "L" );
 cLegendUncertainties    ->  AddEntry( hSyst1D_BrR, "Br. Ratio", "L" );
 cLegendUncertainties    ->  AddEntry( hSyst1D_ITP, "Tracking", "L" );
 //
 gROOT->SetBatch();
 //
 TCanvas*    cDrawFinalUncertainties_1D  =   new TCanvas("cDrawFinalUncertainties_1D","cDrawFinalUncertainties_1D",1000,1000);
 //
 gPad->SetLogx();
 hSyst1D_ALL->Draw("SAME");
 hSyst1D_TRK->Draw("SAME");
 hSyst1D_PID->Draw("SAME");
 hSyst1D_SEX->Draw("SAME");
 hSyst1D_BrR->Draw("SAME");
 hSyst1D_ITP->Draw("SAME");
 hSysU1D_STD->Draw("SAME");
 cLegendUncertainties->Draw("SAME");
 //
 delete cDrawFinalUncertainties_1D;
 //
 for ( Int_t iPT = 0; iPT < hSysU2D_STD->GetNbinsX(); iPT++ )    {
     TCanvas*    cDrawFinalUncertainties_2D  =   new TCanvas("cDrawFinalUncertainties_2D","cDrawFinalUncertainties_2D",1000,1000);
     gPad                                ->  SetLogx();
     //
     auto        hSlice_ALL              =   hSyst2D_ALL ->  ProjectionY(Form("2DFL_%i",iPT),iPT+1,iPT+1);
     hSlice_ALL  ->  SetMinimum( 0 );
     hSlice_ALL  ->  SetMaximum( hSyst2D_ALL->GetMaximum() );
     hSlice_ALL  ->  SetLineWidth(3);
     hSlice_ALL  ->  Draw();
     uLatex      ->  SetTextSize(0.035);
     uLatex      ->  DrawLatexNDC(0.18,0.05,Form("#it{p}_{T,#phi_{1}} (GeV/#it{c}) #in [%.1f;%.1f]",fArrPT2D[iPT],fArrPT2D[iPT+1]));
     cLegendUncertainties->Draw("SAME");
     //
     auto        hSlice_STD              =   hSysU2D_STD ->  ProjectionY(Form("2DFL_STD_%i",iPT),iPT+1,iPT+1);
     hSlice_STD  ->  Draw("SAME");
     //
     auto        hSlice_TRK              =   hSyst2D_TRK ->  ProjectionY(Form("2DFL_TRK_%i",iPT),iPT+1,iPT+1);
     hSlice_TRK  ->  Draw("SAME");
     //
     auto        hSlice_PID              =   hSyst2D_PID ->  ProjectionY(Form("2DFL_PID_%i",iPT),iPT+1,iPT+1);
     hSlice_PID  ->  Draw("SAME");
     //
     auto        hSlice_SEX              =   hSyst2D_SEX ->  ProjectionY(Form("2DFL_SEX_%i",iPT),iPT+1,iPT+1);
     hSlice_SEX  ->  Draw("SAME");
     //
     auto        hSlice_BrR              =   hSyst2D_BrR ->  ProjectionY(Form("2DFL_BrR_%i",iPT),iPT+1,iPT+1);
     hSlice_BrR  ->  Draw("SAME");
     //
     auto        hSlice_ITP              =   hSyst2D_ITP ->  ProjectionY(Form("2DFL_ITP_%i",iPT),iPT+1,iPT+1);
     hSlice_ITP  ->  Draw("SAME");
     //
     cDrawFinalUncertainties_2D  ->  SaveAs(Form(kSystematicsPlot + TString(Form("/hFullSyst_2D_%i.pdf",iPT)), (TString("Yield")+kFolder).Data()));
     //
     delete      cDrawFinalUncertainties_2D;
 }
 //
 cLegendUncertainties->AddEntry( hSystRT_EXT, "Extrapolation" );
 //
 TCanvas*    cDrawFinalUncertainties_PT  =   new TCanvas("cDrawFinalUncertainties_PT","cDrawFinalUncertainties_PT",1000,1000);
 //
 //  !TODO: FIX IN CREATING FILE
 //  --- --- --- --- --- --- --- ---
 uSetHisto(hSystPT_ALL,"SPT MPT 12D");
 hSystPT_ALL ->  GetYaxis()  ->  SetTitle("Fractional Uncertainty (%)");
 //  --- --- --- --- --- --- --- ---
 hSystPT_ALL->Draw("SAME HIST");
 hSystPT_EXT->Draw("SAME");
 hSystPT_SEX->Draw("SAME");
 hSystPT_TRK->Draw("SAME");
 hSystPT_PID->Draw("SAME");
 hSysUPT_STD->Draw("SAME");
 cLegendUncertainties->Draw("SAME");
 //
 cDrawFinalUncertainties_PT  ->  SaveAs(Form(kSystematicsPlot + TString("/hFullSyst_PT.pdf"), (TString("Yield")+kFolder).Data()));
 //
 delete cDrawFinalUncertainties_PT;
 //
 cLegendUncertainties    ->  AddEntry( hSystRT_NMP, "Normalisation" );
 cLegendUncertainties    ->  SetNColumns(3);
 //
 TCanvas*    cDrawFinalUncertainties_RT  =   new TCanvas("cDrawFinalUncertainties_RT","cDrawFinalUncertainties_RT",1000,1000);
 //
 hSystRT_ALL->Draw("SAME HIST");
 hSystRT_EXT->Draw("SAME HIST");
 hSystRT_TRK->Draw("SAME HIST");
 hSystRT_PID->Draw("SAME HIST");
 hSystRT_SEX->Draw("SAME HIST");
 hSystRT_BrR->Draw("SAME HIST");
 hSystRT_ITP->Draw("SAME HIST");
 hSystRT_STD->Draw("SAME HIST");
 hSystRT_NMM->Draw("SAME HIST");
 hSystRT_NMP->Draw("SAME HIST");
 cLegendUncertainties->Draw("SAME");
 //
 cDrawFinalUncertainties_RT  ->  SaveAs(Form(kSystematicsPlot + TString("/hFullSyst_RT.pdf"), (TString("Yield")+kFolder).Data()));
 //
 delete cDrawFinalUncertainties_RT;
 //
 gROOT->SetBatch(kFALSE);
 //
 TFile*  fOutFile    =   new TFile(Form("%s/FullSystematics.root",Form(kAnalysis_Systemt_Dir,  (TString("Yield")+kFolder).Data())),"RECREATE");
 hSyst1D_ALL->Scale(0.01);
 hSyst2D_ALL->Scale(0.01);
 hSystPT_ALL->Scale(0.01);
 hSystRT_ALL->Scale(0.01);
 hSyst1D_ALL->Write();
 hSyst2D_ALL->Write();
 hSystPT_ALL->Write();
 hSystRT_ALL->Write();
 //
 fOutFile->Close();
}
 */
