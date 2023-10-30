// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../../inc/AliAnalysisPhiPair.h"
#include "../GeneralAnalysis.cxx"
#include "RooMsgService.h"

/*
void
SY_AN_PID
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
    SetStyle();
    fSetAllBins();
    //
    if ( kDoYield ) {
        //  --- Recovering Standard Analysis
        TFile*  insFile_Data_YL     =   new TFile   (Form(kASigExtp_FitCheckRst,(TString("Yield")+kFolder).Data()));
        //
        auto        h1D_Nraw_stat   =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "h1D_Nres_stat" );
        auto        h2D_Nraw_stat   =   uLoadHistograms<0,TH2F> ( insFile_Data_YL,  "h2D_Nres_stat" );
        //
        std::vector<TH1F*>  k1D_Variations;
        std::vector<TFile*> k1D_VarInFiles;
        for ( auto kCurrent_Syst : kSyst_PID_XD_Options ) {
            push_to_front( k1D_VarInFiles, new TFile ( Form(kASigExtp_FitCheckRst,(TString("Yield")+kFolder+TString("/Systematics/PID/")+(kCurrent_Syst)).Data(),"1D") ) );
            //
            auto    kCurrent_Variation  =   uLoadHistograms<0,TH1F> ( k1D_VarInFiles.at(0),  "h1D_Nres_stat", kCurrent_Syst );
            kCurrent_Variation  ->  SetName(kCurrent_Syst.Data());
            k1D_Variations.push_back( kCurrent_Variation );
        }
        //
        std::vector<TH2F*>  k2D_Variations;
        std::vector<TFile*> k2D_VarInFiles;
        for ( auto kCurrent_Syst : kSyst_PID_XD_Options ) {
            push_to_front( k2D_VarInFiles, new TFile ( Form(kASigExtp_FitCheckRst,(TString("Yield")+kFolder+TString("/Systematics/PID/")+(kCurrent_Syst)).Data(),"2D") ) );
            auto    kCurrent_Variation  =   uLoadHistograms<0,TH2F> ( k2D_VarInFiles.at(0),  "h2D_Nres_stat", kCurrent_Syst );
            kCurrent_Variation  ->  SetName(kCurrent_Syst.Data());
            k2D_Variations.push_back( kCurrent_Variation );
        }
        //
        fSetAllCustomFunctions();
        uSysEvaluate_Extrapolation_Custom1D  ( h1D_Nraw_stat, k1D_Variations, h2D_Nraw_stat, k2D_Variations, TString(Form(kAnalysis_Systemt_Dir,(TString("Yield")+kFolder).Data()))+TString("PID"), "", false );
    }
    if ( kDoMultiplicity ) {
        
    }
}
*/
