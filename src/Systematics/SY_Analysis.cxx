#include "../../inc/AliAnalysisPhiPair.h"
#include "./GeneralAnalysis.cxx"
// !TODO: Nothing

void
SY_Analysis
 ( TString fOption = "Yield", TString fType = "SEX" )    {
     fSetAllBins();
     if ( fType.Contains("SEX") ) {
         //  --- Load Files
         TFile*  insFile_Data_YL     =   new TFile   (Form(kASigExtr_FitCheckRst,"Yield"));
         TFile*  insFile_Effc_YL     =   new TFile   (Form(kAnalysis_MCTruthHist,"Yield"));
         //
         //  --- Load Histograms
         auto        fHEventCount    =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "fQC_Event_Enum_FLL" );
         auto        h1D_Nraw        =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "h1D_Nraw_" );
         auto        h1D_Nrec        =   uLoadHistograms<0,TH1F> ( insFile_Effc_YL,  "h1D_Nrec" );
         auto        h1D_Ngen        =   uLoadHistograms<0,TH1F> ( insFile_Effc_YL,  "h1D_Ngen" );
         auto        h2D_Nraw        =   uLoadHistograms<0,TH2F> ( insFile_Data_YL,  "anSS2D_" );
         auto        h1D_Nrec_2Db    =   uLoadHistograms<0,TH1F> ( insFile_Effc_YL,  "h1D_Nrec_2Db" );
         auto        h1D_Ngen_2Db    =   uLoadHistograms<0,TH1F> ( insFile_Effc_YL,  "h1D_Ngen_2Db" );
         //
         //  --- Minimum Bias Normalisation
         auto        kN_Trg          =   (fHEventCount->GetBinContent(kEventCount::kTrigger));
         auto        kN_Vtx          =   (fHEventCount->GetBinContent(kEventCount::kVertex));
         auto        kN_MB           =   (fHEventCount->GetBinContent(kEventCount::kVertex10));
         Double_t    f1DCorrection   =   (1./kBR)        *(1./kN_MB) *(kTriggerEff/1.)   *(kN_Vtx/kN_Trg);
         Double_t    f2DCorrection   =   (1./(kBR*kBR))  *(1./kN_MB) *(kTriggerEff/1.)   *(kN_Vtx/kN_Trg);
         if ( kEnergy == 5 ) {   f1DCorrection   *=  kTriggerEff15n17pq; f2DCorrection   *=  kTriggerEff15n17pq; }
         if ( kEnergy == 7 ) {   f1DCorrection   *=  kTriggerEff10bcdef; f2DCorrection   *=  kTriggerEff10bcdef; }
         //
                     h1D_Nraw        ->  Scale( 1., "width" );
         TH1F*       h1D_Nraw_stat   =   uEfficiencyCorrection1D ( h1D_Nraw, h1D_Nrec, h1D_Ngen, f1DCorrection );
         SetAxis(h1D_Nraw_stat,"PT1D");
         //
         std::vector<TH1F*>  k1D_Variations;
         std::vector<TFile*> k1D_InputFiles;
         for ( auto kCurrent_Option : kSyst_SEX_1D_Options ) {
             auto   kCurrent_File       =   new TFile   ( Form(kASigExtr_FitChkRstSY,"Yield/Systematics/","1D",kCurrent_Option.Data()) );
             auto   h1D_Current_Nraw    =   uLoadHistograms<0,TH1F> ( kCurrent_File, Form("h1D_Nraw_%s",kCurrent_Option.Data()) );
             k1D_InputFiles.push_back( kCurrent_File );
             h1D_Current_Nraw           ->  Scale( 1., "width" );
             TH1F*       h1D_Current_Nraw_stat  =   uEfficiencyCorrection1D ( h1D_Current_Nraw, h1D_Nrec, h1D_Ngen, f1DCorrection );
             h1D_Current_Nraw_stat      ->  SetName( kCurrent_Option.Data() );
             SetAxis(h1D_Current_Nraw_stat,"PT1D");
             k1D_Variations.push_back( h1D_Current_Nraw_stat );
         }
         GeneralAnalysis( h1D_Nraw_stat, k1D_Variations, Form(kAnalysis_SigExtr_Dir,"Yield/Systematics") );
         //
                    h2D_Nraw        ->  Scale( 1., "width" );
         TH2F*      h2D_Nraw_stat   =   uEfficiencyCorrection2D ( h2D_Nraw, h1D_Nrec_2Db, h1D_Ngen_2Db, f2DCorrection );
         SetAxis(h2D_Nraw_stat,"PT2D");
         //
         std::vector<TH2F*>  k2D_Variations;
         std::vector<TFile*> k2D_InputFiles;
         for ( auto kCurrent_Option : kSyst_SEX_2D_Options ) {
             auto   kCurrent_File       =   new TFile   ( Form(kASigExtr_FitChkRstSY,"Yield/Systematics/","2D",kCurrent_Option.Data()) );
             auto   h2D_Current_Nraw    =   uLoadHistograms<0,TH2F> ( kCurrent_File, Form("anSS2D_%s",kCurrent_Option.Data()) );
             k2D_InputFiles.push_back( kCurrent_File );
             h2D_Current_Nraw           ->  Scale( 1., "width" );
             TH2F*       h2D_Current_Nraw_stat  =   uEfficiencyCorrection2D ( h2D_Current_Nraw, h1D_Nrec_2Db, h1D_Ngen_2Db, f2DCorrection );
             h2D_Current_Nraw_stat      ->  SetName( kCurrent_Option.Data() );
             SetAxis(h2D_Current_Nraw_stat,"PT2D");
             k2D_Variations.push_back( h2D_Current_Nraw_stat );
         }
         GeneralAnalysis( h2D_Nraw_stat, k2D_Variations, Form(kAnalysis_SigExtr_Dir,"Yield/Systematics") );
         //
         GeneralAnalysis( h1D_Nraw_stat, k1D_Variations, h2D_Nraw_stat, k2D_Variations, Form(kAnalysis_SigExtr_Dir,"Yield/Systematics") );
         //
         for ( auto kFile : k1D_InputFiles ) kFile->Close();
         for ( auto kFile : k2D_InputFiles ) kFile->Close();
         //
         insFile_Data_YL->Close();
         insFile_Effc_YL->Close();
     }
}
