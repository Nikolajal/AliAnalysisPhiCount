// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../../inc/AliAnalysisPhiPair.h"
#include "../GeneralAnalysis.cxx"
#include "RooMsgService.h"

void
SY_AN_SigExtraction
( TString fOption = "yield mult", TString kFolder="_p_p__5TeV")    {
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
        TFile*  insFile_Data_YL     =   new TFile   (Form(kASigExtr_FitCheckRst,(TString("Yield")+kFolder).Data()));
        TFile*  insFile_Dat2_YL     =   new TFile   (Form(kASigExtp_FitCheckRst,(TString("Yield")+kFolder).Data()));
        TFile*  insFile_Effc_YL     =   new TFile   (Form(kAnalysis_MCTruthHist,(TString("Yield")+kFolder).Data()));
        //
        auto        fHEventCount    =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "fQC_Event_Enum_FLL" );
        auto        h1D_Nraw        =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "h1D_Nraw_" );
        auto        h1D_Nrec        =   uLoadHistograms<0,TH1F> ( insFile_Effc_YL,  "h1D_Nrec" );
        auto        h1D_Ngen        =   uLoadHistograms<0,TH1F> ( insFile_Effc_YL,  "h1D_Ngen" );
        auto        h2D_Nraw        =   uLoadHistograms<0,TH2F> ( insFile_Data_YL,  "anSS2D_" );
        auto        h1D_Nrec_2Db    =   uLoadHistograms<0,TH1F> ( insFile_Effc_YL,  "h1D_Nrec_2Db" );
        auto        h1D_Ngen_2Db    =   uLoadHistograms<0,TH1F> ( insFile_Effc_YL,  "h1D_Ngen_2Db" );
        auto        h2D_Next_stat   =   uLoadHistograms<0,TH1F> ( insFile_Dat2_YL,  "h2D_Nres_stat_stat_0" );
        //
        //      TODO: Make a separate function to calculate normalisation
        //  --- Minimum Bias Normalisation
        auto        kN_PU           =   (fHEventCount->GetBinContent(kEventCount::kPU_MB));
        auto        kN_Trg          =   -kN_PU +(fHEventCount->GetBinContent(kEventCount::kTrigger));
        auto        kN_Vtx          =   -kN_PU +(fHEventCount->GetBinContent(kEventCount::kVertex));
        auto        kN_MB           =   -kN_PU +(fHEventCount->GetBinContent(kEventCount::kVertex10));
        Double_t    f1DCorrection   =   (1./kBR)        *(1./kN_MB) *(kTriggerEff/1.)   *(kN_Vtx/kN_Trg);
        Double_t    f2DCorrection   =   (1./(kBR*kBR))  *(1./kN_MB) *(kTriggerEff/1.)   *(kN_Vtx/kN_Trg);
        if ( is_pp_anl && ( kEnergy == 5 ) )    {   f1DCorrection   *=  kTriggerEff15n17pq; f2DCorrection   *=  kTriggerEff15n17pq; }
        if ( is_pp_anl && ( kEnergy == 7 ) )    {   f1DCorrection   *=  kTriggerEff10bcdef; f2DCorrection   *=  kTriggerEff10bcdef; }
        //
                    h1D_Nraw        ->  Scale(1.,"width");
        //
        TH1F*       h1D_Nraw_stat   =   uEfficiencyCorrection1D ( h1D_Nraw, h1D_Nrec, h1D_Ngen, f1DCorrection );
        SetAxis(h1D_Nraw_stat,"PT1D");
        h1D_Nraw_stat   ->  SetName("h1D_Nraw_stat");
        //
        std::vector<TH1F*>  k1D_Variations;
        std::vector<TFile*> k1D_VarInFiles;
        for ( auto kCurrent_Syst : kSyst_SEX_1D_Options ) {
            hName           =   Form("h1D_Nraw_%s",kCurrent_Syst.Data());
            push_to_front( k1D_VarInFiles, new TFile ( Form(kASigExtr_FitChkRstSY,(TString("Yield")+kFolder+TString("/Systematics/")).Data(),(kCurrent_Syst).Data(),"1D") ) );
            auto    kCurrent_Variation  =   uEfficiencyCorrection1D ( (TH1F*)((k1D_VarInFiles.at(0))->Get(hName)), h1D_Nrec, h1D_Ngen, f1DCorrection );
            kCurrent_Variation  ->  Scale(1.,"width");
            kCurrent_Variation  ->  SetName(kCurrent_Syst.Data());
            k1D_Variations.push_back( kCurrent_Variation );
        }
        //
                    h2D_Nraw        ->  Scale(1.,"width");
        //
        TH2F*       h2D_Nres_stat   =   uEfficiencyCorrection2D ( h2D_Nraw, h1D_Nrec_2Db, h1D_Ngen_2Db, f2DCorrection );
        SetAxis(h2D_Nres_stat,"PT2D");
        h2D_Nres_stat   ->  SetName("h2D_Nraw_stat");
        //
        std::vector<TH2F*>  k2D_Variations;
        std::vector<TFile*> k2D_VarInFiles;
        for ( auto kCurrent_Syst : kSyst_SEX_2D_Options ) {
            hName           =   Form("anSS2D_%s",kCurrent_Syst.Data());
            push_to_front( k2D_VarInFiles, new TFile ( Form(kASigExtr_FitChkRstSY,(TString("Yield")+kFolder+TString("/Systematics/")).Data(),(kCurrent_Syst).Data(),"2D") ) );
            auto    kCurrent_Variation  =   uEfficiencyCorrection2D ( (TH2F*)((k2D_VarInFiles.at(0))->Get(hName)), h1D_Nrec_2Db, h1D_Ngen_2Db, f2DCorrection );
            kCurrent_Variation  ->  Scale(1.,"width");
            kCurrent_Variation  ->  SetName(kCurrent_Syst.Data());
            k2D_Variations.push_back( kCurrent_Variation );
        }
        //
        fSetAllCustomFunctions();
        uSysEvaluate_Extrapolation_Custom1D  ( h1D_Nraw_stat, k1D_Variations, h2D_Nres_stat, h2D_Next_stat, k2D_Variations, TString(Form(kAnalysis_Systemt_Dir,(TString("Yield")+kFolder).Data()))+TString("SEX"), "", false );
    }
    if ( kDoMultiplicity ) {
        //  --- Recovering Standard Analysis
        TFile*  insFile_Data_ML     =   new TFile   (Form(kASigExtr_FitCheckRst,(TString("Multiplicity")+kFolder).Data()));
        TFile*  insFile_Dat2_YL     =   new TFile   (Form(kASigExtp_FitCheckRst,(TString("Multiplicity")+kFolder).Data()));
        TFile*  insFile_Effc_ML     =   new TFile   (Form(kAnalysis_MCTruthHist,(TString("Multiplicity")+kFolder).Data()));
        //
        auto        fHEventCount    =   uLoadHistograms<0,TH1F> ( insFile_Data_ML,  "fQC_Event_Enum_FLL"    );
        auto        fHEventCntMlt   =   uLoadHistograms<0,TH1F> ( insFile_Data_ML,  "fQC_Event_Enum_V0M"    );
        auto        h1D_Nraw_MT     =   uLoadHistograms<1,TH1F> ( insFile_Data_ML,  "h1D_Nraw_MT_%i",       "h1D_Nraw_MT_%i" );
        auto        h1D_Nrec        =   uLoadHistograms<0,TH1F> ( insFile_Effc_ML,  "h1D_Nrec" );
        auto        h1D_Ngen        =   uLoadHistograms<0,TH1F> ( insFile_Effc_ML,  "h1D_Ngen" );
        auto        h2D_Nraw_MT     =   uLoadHistograms<1,TH2F> ( insFile_Data_ML,  "anSS2D_MT_%i",         "h2D_Nraw_MT_%i" );
        auto        h1D_Nrec_2Db    =   uLoadHistograms<0,TH1F> ( insFile_Effc_ML,  "h1D_Nrec_2Db" );
        auto        h1D_Ngen_2Db    =   uLoadHistograms<0,TH1F> ( insFile_Effc_ML,  "h1D_Ngen_2Db" );
        auto        h2D_Next_stat_MT=   uLoadHistograms<1,TH1F> ( insFile_Dat2_YL,  "h2D_Nres_stat_MT_%i_stat_0" );
        //
        Int_t   iMult = -1;
        for ( auto kTarget : h1D_Nraw_MT )   {
            iMult++;
            //
            auto    f1DCorrection    =   (1.)/(fEvaluateINELgt0(iMult-1,fHEventCntMlt) * kBR );
            auto    f2DCorrection    =   (1.)/(fEvaluateINELgt0(iMult-1,fHEventCntMlt) * kBR * kBR );
            if ( is_pb_anl ) {
                f1DCorrection   =   (1.)/(fEvaluateINELgt0(iMult-1,fHEventCntMlt,false) * kBR ) * 2;
                f2DCorrection   =   (1.)/(fEvaluateINELgt0(iMult-1,fHEventCntMlt,false) * kBR * kBR ) * 2;
            }
            //
                        h1D_Nraw_MT.at(iMult)        ->  Scale(1.,"width");
            //
            TH1F*       h1D_Nraw_stat   =   uEfficiencyCorrection1D ( h1D_Nraw_MT.at(iMult), h1D_Nrec, h1D_Ngen, f1DCorrection );
            SetAxis(h1D_Nraw_stat,"PT1D");
            h1D_Nraw_stat   ->  SetName(Form("h1D_Nraw_stat_MT_%i",iMult));
            //
            std::vector<TH1F*>  k1D_Variations;
            std::vector<TFile*> k1D_VarInFiles;
            for ( auto kCurrent_Syst : kSyst_SEX_1D_Options ) {
                hName           =   Form("h1D_Nraw_%s",kCurrent_Syst.Data());
                push_to_front( k1D_VarInFiles, new TFile ( Form(kASigExtr_FitChkRstSY,(TString("Multiplicity")+kFolder+TString("/Systematics/")).Data(),(kCurrent_Syst).Data(),Form("1D_MLT_%i",iMult)) ) );
                auto    kCurrent_Variation  =   uEfficiencyCorrection1D ( (TH1F*)((k1D_VarInFiles.at(0))->Get(hName)), h1D_Nrec, h1D_Ngen, f1DCorrection );
                kCurrent_Variation  ->  Scale(1.,"width");
                kCurrent_Variation  ->  SetName(kCurrent_Syst.Data());
                k1D_Variations.push_back( kCurrent_Variation );
            }
            //
                        h2D_Nraw_MT.at(iMult)       ->  Scale(1.,"width");
            //
            TH2F*       h2D_Nraw_stat   =   uEfficiencyCorrection2D ( h2D_Nraw_MT.at(iMult), h1D_Nrec_2Db, h1D_Ngen_2Db, f2DCorrection );
            SetAxis(h2D_Nraw_stat,"PT2D");
            h2D_Nraw_stat   ->  SetName(Form("h2D_Nraw_stat_MT_%i",iMult));
            //
            std::vector<TH2F*>  k2D_Variations;
            std::vector<TFile*> k2D_VarInFiles;
            for ( auto kCurrent_Syst : kSyst_SEX_2D_Options ) {
                hName           =   Form("anSS2D_%s",kCurrent_Syst.Data());
                push_to_front( k2D_VarInFiles, new TFile ( Form(kASigExtr_FitChkRstSY,(TString("Multiplicity")+kFolder+TString("/Systematics/")).Data(),(kCurrent_Syst).Data(),Form("2D_MLT_%i",iMult)) ) );
                auto    kCurrent_Variation  =   uEfficiencyCorrection2D ( (TH2F*)((k2D_VarInFiles.at(0))->Get(hName)), h1D_Nrec_2Db, h1D_Ngen_2Db, f2DCorrection );
                kCurrent_Variation  ->  Scale(1.,"width");
                kCurrent_Variation  ->  SetName(kCurrent_Syst.Data());
                k2D_Variations.push_back( kCurrent_Variation );
            }
            //
            fSetAllCustomFunctions();
            uSysEvaluate_Extrapolation_Custom1D  ( h1D_Nraw_stat, k1D_Variations, h2D_Nraw_stat, h2D_Next_stat_MT.at(iMult), k2D_Variations, TString(Form(kAnalysis_Systemt_Dir,(TString("Multiplicity")+kFolder).Data()))+TString(Form("/MLT_%i/SEX/",iMult)), "", false );
        }
    }
}
