// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../../inc/AliAnalysisPhiPair.h"
#include "../GeneralAnalysis.cxx"
#include "RooMsgService.h"

void
SystematicsGeneralAnalysis
( TString fOption = "full", TString kFolder = "_p_p__5TeV", TString kSource = "PID" )    {
    //
    //  Setting general analysis
    //
    //  --- Option chosing
    if ( !fChooseOption(fOption) ) return;
    //
    //  --- Systematics options
    std::vector<TString> kOptions1D;
    std::vector<TString> kOptions2D;
    if( kSource.Contains("PID") ) {
        kOptions1D = kSyst_PID_XD_Options;
        kOptions2D = kSyst_PID_XD_Options;
    }
    if( kSource.Contains("TRK") ) {
        kOptions1D = kSyst_TRK_XD_Options;
        kOptions2D = kSyst_TRK_XD_Options;
    }
    //
    //  Generating the binning array--------------------------------------------------------------------------
    SetStyle();
    fSetAllBins();
    //
    if ( kDoYield ) {
        //  --- Recovering Standard Analysis
        TFile*  insFile_Data_YL     =   new TFile   (Form(kASigExtp_FitCheckRst,(TString("Yield")+kFolder).Data()));
        //
        auto        h1D_Nres_stat   =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "h1D_Nres_stat" );
        auto        h2D_Nres_stat   =   uLoadHistograms<0,TH2F> ( insFile_Data_YL,  "h2D_Nres_stat" );
        auto        h2D_Next_stat   =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "h2D_Nres_stat_stat_0" );
        //
        std::vector<TH1F*>  k1D_Variations;
        std::vector<TFile*> k1D_VarInFiles;
        for ( auto kCurrent_Syst : kOptions1D ) {
            push_to_front( k1D_VarInFiles, new TFile ( Form(kASigExtp_FitCheckRst,(TString("Yield")+kFolder+TString("/Systematics/")+kSource+TString("/")+(kCurrent_Syst)).Data(),"1D") ) );
            //
            auto    kCurrent_Variation  =   uLoadHistograms<0,TH1F> ( k1D_VarInFiles.at(0),  "h1D_Nres_stat", kCurrent_Syst );
            kCurrent_Variation  ->  SetName(kCurrent_Syst.Data());
            k1D_Variations.push_back( kCurrent_Variation );
        }
        //
        std::vector<TH2F*>  k2D_Variations;
        std::vector<TFile*> k2D_VarInFiles;
        for ( auto kCurrent_Syst : kOptions2D ) {
            push_to_front( k2D_VarInFiles, new TFile ( Form(kASigExtp_FitCheckRst,(TString("Yield")+kFolder+TString("/Systematics/")+kSource+TString("/")+(kCurrent_Syst)).Data(),"2D") ) );
            auto    kCurrent_Variation  =   uLoadHistograms<0,TH2F> ( k2D_VarInFiles.at(0),  "h2D_Nres_stat", kCurrent_Syst );
            kCurrent_Variation  ->  SetName(kCurrent_Syst.Data());
            k2D_Variations.push_back( kCurrent_Variation );
        }
        //
        fSetAllCustomFunctions();
        uSysEvaluate_Extrapolation_Custom1D  ( h1D_Nres_stat, k1D_Variations, h2D_Nres_stat, h2D_Next_stat, k2D_Variations, TString(Form(kAnalysis_Systemt_Dir,(TString("Yield")+kFolder).Data()))+kSource, "", false );
    }
    if ( kDoMultiplicity ) {
        //  --- Recovering Standard Analysis
        TFile      *insFile_Data_ML     = new TFile   (Form(kASigExtp_FitCheckRst,(TString("Multiplicity")+kFolder).Data()));
        auto        h1D_Nres_stat_MT    = uLoadHistograms<1,TH1F> ( insFile_Data_ML,  "h1D_Nres_stat_MT_%i" );
        auto        h2D_Nres_stat_MT    = uLoadHistograms<1,TH2F> ( insFile_Data_ML,  "h2D_Nres_stat_MT_%i" );
        auto        h2D_Next_stat_MT    = uLoadHistograms<1,TH1F> ( insFile_Data_ML,  "h2D_Nres_stat_MT_%i_stat_0" );
        //
        //  Loop on mult bins
        for ( int iMult = 0; iMult <= nBinMult; iMult++ ) {
            auto        h1D_Nres_stat   =   h1D_Nres_stat_MT.at(iMult);
            auto        h2D_Nres_stat   =   h2D_Nres_stat_MT.at(iMult);
            auto        h2D_Next_stat   =   h2D_Next_stat_MT.at(iMult);
            //
            std::vector<TH1F*>  k1D_Variations;
            std::vector<TFile*> k1D_VarInFiles;
            for ( auto kCurrent_Syst : kOptions1D ) {
                push_to_front( k1D_VarInFiles, new TFile ( Form(kASigExtp_FitCheckRst,(TString("Multiplicity")+kFolder+TString("/Systematics/")+kSource+TString("/")+(kCurrent_Syst)).Data(),"1D") ) );
                //
                auto    kCurrent_Variation  =   uLoadHistograms<0,TH1F> ( k1D_VarInFiles.at(0),  Form("h1D_Nres_stat_MT_%i",iMult), Form("1D_%s_MT_%i",kCurrent_Syst.Data(),iMult) );
                k1D_Variations.push_back( kCurrent_Variation );
            }
            //
            std::vector<TH2F*>  k2D_Variations;
            std::vector<TFile*> k2D_VarInFiles;
            for ( auto kCurrent_Syst : kOptions2D ) {
                push_to_front( k2D_VarInFiles, new TFile ( Form(kASigExtp_FitCheckRst,(TString("Multiplicity")+kFolder+TString("/Systematics/")+kSource+TString("/")+(kCurrent_Syst)).Data(),"2D") ) );
                auto    kCurrent_Variation  =   uLoadHistograms<0,TH2F> ( k2D_VarInFiles.at(0),  Form("h2D_Nraw_stat_MT_%i",iMult), Form("2D_%s_MT_%i",kCurrent_Syst.Data(),iMult) );
                k2D_Variations.push_back( kCurrent_Variation );
            }
            //
            fSetAllCustomFunctions();
            uSysEvaluate_Extrapolation_Custom1D  ( h1D_Nres_stat, k1D_Variations, h2D_Nres_stat, h2D_Next_stat, k2D_Variations, TString(Form(kAnalysis_Systemt_Dir+TString(Form("/MLT_%i/",iMult)),(TString("Multiplicity")+kFolder).Data()))+kSource, "", false );
            //
            //  --- Close all files
            for ( auto kCurrentFile : k1D_VarInFiles ) kCurrentFile->Close();
            for ( auto kCurrentFile : k2D_VarInFiles ) kCurrentFile->Close();
        }
    }
}
