// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../../inc/AliAnalysisPhiPair.h"
#include "../GeneralAnalysis.cxx"
#include "RooMsgService.h"
//!
void
SystematicsGeneralAnalysisLoaded
 ( TString fOption = "full", TString kFolder = "_p_p__5TeV", bool kMatBud = true, bool bTK = true )    {
     //
     //  Setting general analysis
     //
     //  --- Option chosing
     if ( !fChooseOption(fOption) ) return;
     //
     //  Generating the binning array--------------------------------------------------------------------------
     SetStyle();
     fSetAllBins();
     //
     auto kFolderExternal = TString("/Users/nrubini/PythiaCustomVx/out_8308_0_12345_7000_2212_2212.root");
    auto kMatBud1D = "hError1D_MB";
    auto kHadInt1D = "hError1D_HI";
    auto kMatBud2D = "hError2D_MB";
    auto kHadInt2D = "hError2D_HI";
    auto kTracki1D = "hError1D_TK";
    auto kTracki2D = "hError2D_TK";
    TString kSource = bTK ? "TKN" : (kMatBud? "MBD" : "HIN");
     if ( kDoYield ) {
         //  --- Recovering Standard Analysis
         TFile*  insFile_Data_YL     =   new TFile   (Form(kASigExtp_FitCheckRst,(TString("Yield")+kFolder).Data()));
         //
         auto        h1D_Nres_stat   =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "h1D_Nres_stat" );
         auto        h2D_Nres_stat   =   uLoadHistograms<0,TH2F> ( insFile_Data_YL,  "h2D_Nres_stat" );
         auto        h2D_Next_stat   =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "h2D_Nres_stat_stat_0" );
         //
         std::vector<TH1F*>  k1D_Variations;
         for ( Int_t iTer = 0; iTer < 100; iTer++ ) {
             auto kCurrentVariation = uRandomisePoints( fSetSystErrors( h1D_Nres_stat, kFolderExternal, bTK ? kTracki1D : (kMatBud? kMatBud1D : kHadInt1D) ));
             kCurrentVariation->SetName(Form("%s%i",kSource.Data(),iTer));
             k1D_Variations.push_back( kCurrentVariation );
         }
         //
         std::vector<TH2F*>  k2D_Variations;
         for ( Int_t iTer = 0; iTer < 100; iTer++ ) {
             auto kCurrentVariation = uRandomisePoints( fSetSystErrors( h2D_Nres_stat, kFolderExternal, bTK ? kTracki2D : (kMatBud? kMatBud2D : kHadInt2D) ));
             kCurrentVariation->SetName(Form("%s%i",kSource.Data(),iTer));
             k2D_Variations.push_back( kCurrentVariation );
         }
         //
         fSetAllCustomFunctions();
         uSysEvaluate_Extrapolation_Custom1D  ( h1D_Nres_stat, k1D_Variations, h2D_Nres_stat, h2D_Next_stat, k2D_Variations, TString(Form(kAnalysis_Systemt_Dir,(TString("Yield")+kFolder).Data()))+kSource, "", "", true );
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
             for ( Int_t iTer = 0; iTer < 100; iTer++ ) {
                 auto kCurrentVariation = uRandomisePoints( fSetSystErrors( h1D_Nres_stat, kFolderExternal, bTK ? kTracki1D : (kMatBud? kMatBud1D : kHadInt1D) ));
                 kCurrentVariation->SetName(Form("%s%i",kSource.Data(),iTer));
                 k1D_Variations.push_back( kCurrentVariation );
             }
             //
             std::vector<TH2F*>  k2D_Variations;
             for ( Int_t iTer = 0; iTer < 100; iTer++ ) {
                 auto kCurrentVariation = uRandomisePoints( fSetSystErrors( h2D_Nres_stat, kFolderExternal, bTK ? kTracki2D : (kMatBud? kMatBud2D : kHadInt2D) ));
                 kCurrentVariation->SetName(Form("%s%i",kSource.Data(),iTer));
                 k2D_Variations.push_back( kCurrentVariation );
             }
             //
             fSetAllCustomFunctions();
             uSysEvaluate_Extrapolation_Custom1D  ( h1D_Nres_stat, k1D_Variations, h2D_Nres_stat, h2D_Next_stat, k2D_Variations, TString(Form(kAnalysis_Systemt_Dir+TString(Form("/MLT_%i/",iMult)),(TString("Multiplicity")+kFolder).Data()))+kSource, "", "", true );
             //
         }
     }
 }
