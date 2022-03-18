// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../../inc/AliAnalysisPhiPair.h"
#include "../GeneralAnalysis.cxx"
#include "RooMsgService.h"

void
SY_AN_SigExtraction
( TString fOption, TString kFolder )    {
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
        //  --- Recovering Standard Analysis
        TFile*  insFile_Data_YL     =   new TFile   (Form(kASigExtr_FitCheckRst,(TString("Yield")+kFolder).Data()));
        TFile*  insFile_Effc_YL     =   new TFile   (Form(kAnalysis_MCTruthHist,(TString("Yield")+kFolder).Data()));
        //
        auto        fHEventCount    =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "fQC_Event_Enum_FLL" );
        auto        h1D_Nraw        =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "h1D_Nraw_" );
        auto        h1D_Nrec        =   uLoadHistograms<0,TH1F> ( insFile_Effc_YL,  "h1D_Nrec" );
        auto        h1D_Ngen        =   uLoadHistograms<0,TH1F> ( insFile_Effc_YL,  "h1D_Ngen" );
        auto        h2D_Nraw        =   uLoadHistograms<0,TH2F> ( insFile_Data_YL,  "anSS2D_" );
        auto        h1D_Nrec_2Db    =   uLoadHistograms<0,TH1F> ( insFile_Effc_YL,  "h1D_Nrec_2Db" );
        auto        h1D_Ngen_2Db    =   uLoadHistograms<0,TH1F> ( insFile_Effc_YL,  "h1D_Ngen_2Db" );
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
        TH2F*       h2D_Nraw_stat   =   uEfficiencyCorrection2D ( h2D_Nraw, h1D_Nrec_2Db, h1D_Ngen_2Db, f2DCorrection );
        SetAxis(h2D_Nraw_stat,"PT2D");
        h2D_Nraw_stat   ->  SetName("h2D_Nraw_stat");
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
        uSysEvaluate_Extrapolation_Custom1D  ( h1D_Nraw_stat, k1D_Variations, h2D_Nraw_stat, k2D_Variations, TString(Form(kAnalysis_Systemt_Dir,(TString("Yield")+kFolder).Data()))+TString("SignalExtraction"), "", false );
    }
    if ( kDoMultiplicity ) {
        //  --- Recovering Standard Analysis
        TFile*  insFile_Data_ML     =   new TFile   (Form(kASigExtr_FitCheckRst,(TString("Multiplicity")+kFolder).Data()));
        TFile*  insFile_Effc_ML     =   new TFile   (Form(kAnalysis_MCTruthHist,(TString("Multiplicity")+kFolder).Data()));
        //
        auto        fHEventCount    =   uLoadHistograms<0,TH1F> ( insFile_Data_ML,  "fQC_Event_Enum_FLL"     );
        auto        fHEventCntMlt   =   uLoadHistograms<0,TH1F> ( insFile_Data_ML,  "fQC_Event_Enum_V0M"     );
        auto        h1D_Nraw_MT     =   uLoadHistograms<1,TH1F> ( insFile_Data_ML,  "h1D_Nraw_MT_%i",       "h1D_Nraw_MT_%i" );
        auto        h1D_Nrec        =   uLoadHistograms<0,TH1F> ( insFile_Effc_ML,  "h1D_Nrec" );
        auto        h1D_Ngen        =   uLoadHistograms<0,TH1F> ( insFile_Effc_ML,  "h1D_Ngen" );
        auto        h2D_Nraw_MT     =   uLoadHistograms<1,TH2F> ( insFile_Data_ML,  "anSS2D_MT_%i",         "h2D_Nraw_MT_%i" );
        auto        h1D_Nrec_2Db    =   uLoadHistograms<0,TH1F> ( insFile_Effc_ML,  "h1D_Nrec_2Db" );
        auto        h1D_Ngen_2Db    =   uLoadHistograms<0,TH1F> ( insFile_Effc_ML,  "h1D_Ngen_2Db" );
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
            uSysEvaluate_Extrapolation_Custom1D  ( h1D_Nraw_stat, k1D_Variations, h2D_Nraw_stat, k2D_Variations, TString(Form(kAnalysis_Systemt_Dir,(TString("Multiplicity")+kFolder).Data()))+TString(Form("SignalExtraction/MLT_%i",iMult)), "", false );
        }
        //
        /*
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
            push_to_front( k1D_VarInFiles, new TFile ( Form(kASigExtr_FitChkRstSY,(TString("Multiplicity")+kFolder+TString("/Systematics/")).Data(),(kCurrent_Syst).Data(),"1D") ) );
            auto    kCurrent_Variation  =   uEfficiencyCorrection1D ( (TH1F*)((k1D_VarInFiles.at(0))->Get(hName)), h1D_Nrec, h1D_Ngen, f1DCorrection );
            kCurrent_Variation  ->  Scale(1.,"width");
            kCurrent_Variation  ->  SetName(kCurrent_Syst.Data());
            k1D_Variations.push_back( kCurrent_Variation );
        }
        //
                    h2D_Nraw        ->  Scale(1.,"width");
        //
        TH2F*       h2D_Nraw_stat   =   uEfficiencyCorrection2D ( h2D_Nraw, h1D_Nrec_2Db, h1D_Ngen_2Db, f2DCorrection );
        SetAxis(h2D_Nraw_stat,"PT2D");
        h2D_Nraw_stat   ->  SetName("h2D_Nraw_stat");
        //
        std::vector<TH2F*>  k2D_Variations;
        std::vector<TFile*> k2D_VarInFiles;
        for ( auto kCurrent_Syst : kSyst_SEX_2D_Options ) {
            hName           =   Form("anSS2D_%s",kCurrent_Syst.Data());
            push_to_front( k2D_VarInFiles, new TFile ( Form(kASigExtr_FitChkRstSY,(TString("Multiplicity")+kFolder+TString("/Systematics/")).Data(),(kCurrent_Syst).Data(),"2D") ) );
            auto    kCurrent_Variation  =   uEfficiencyCorrection2D ( (TH2F*)((k2D_VarInFiles.at(0))->Get(hName)), h1D_Nrec_2Db, h1D_Ngen_2Db, f2DCorrection );
            kCurrent_Variation  ->  Scale(1.,"width");
            kCurrent_Variation  ->  SetName(kCurrent_Syst.Data());
            k2D_Variations.push_back( kCurrent_Variation );
        }
        //
        fSetAllCustomFunctions();
        uSysEvaluate_Extrapolation_Custom1D  ( h1D_Nraw_stat, k1D_Variations, h2D_Nraw_stat, k2D_Variations, TString(Form(kAnalysis_Systemt_Dir,(TString("Multiplicity")+kFolder).Data()))+TString("SignalExtraction"), "", false );
        
        
        //
        //  --- Output File for Fit Check
        TFile*  outFile_Chck_ML     =   new TFile(Form(kASigExtp_FitCheckPlt,(TString("Multiplicity")+kFolder).Data()),"recreate");
        //
        //  --- Build output structure
        Int_t   iMult = -1;
        std::vector<TString>    kPlotFolder1D;
        std::vector<TString>    kPlotFolder2D;
        std::vector<TH1F*>      h1D_Nres_stat_MT;
        std::vector<TH1F*>      h1D_Nres_syst_MT;
        std::vector<TH2F*>      h2D_Nraw_stat_MT;
        std::vector<TH2F*>      h2D_Nraw_syst_MT;
        std::vector<std::map<TString,std::tuple<Float_t,Float_t,Float_t>>>                  kMult1D_Results;
        std::vector<std::vector<std::map<TString,std::tuple<Float_t,Float_t,Float_t>>>>     kMult2D_Results;
        TString kPlotFolder1D_Full  =   TString( Form(kASigExtp_Plot_Direct,(TString("Multiplicity")+kFolder).Data()) ) +TString(Form("/1D/"));
        TString kPlotFolder2D_Full  =   TString( Form(kASigExtp_Plot_Direct,(TString("Multiplicity")+kFolder).Data()) ) +TString(Form("/1D/"));
        gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,(TString("Multiplicity")+kFolder).Data()))+TString(Form("/1D/")));
        gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,(TString("Multiplicity")+kFolder).Data()))+TString(Form("/2D/")));
        for ( auto kTarget : h1D_Nraw_MT )   {
            iMult++;
            //
            fSetAllCustomFunctions();
            //
            auto    f1DCorrection    =   (1.)/(fEvaluateINELgt0(iMult-1,fHEventCntMlt) * kBR );
            auto    f2DCorrection    =   (1.)/(fEvaluateINELgt0(iMult-1,fHEventCntMlt) * kBR * kBR );
            if ( is_pb_anl ) {
                f1DCorrection   =   (1.)/(fEvaluateINELgt0(iMult-1,fHEventCntMlt,false) * kBR ) * 2;
                f2DCorrection   =   (1.)/(fEvaluateINELgt0(iMult-1,fHEventCntMlt,false) * kBR * kBR ) * 2;
            }
            //
            gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,(TString("Multiplicity")+kFolder).Data()))+TString(Form("/MLT_%i/1D/",iMult)));
            gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,(TString("Multiplicity")+kFolder).Data()))+TString(Form("/MLT_%i/2D/",iMult)));
            gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,(TString("Multiplicity")+kFolder).Data()))+TString(Form("/MLT_%i/Full/",iMult)));
            kPlotFolder1D.push_back( TString( Form(kASigExtp_Plot_Direct,(TString("Multiplicity")+kFolder).Data()) ) +TString(Form("/MLT_%i/1D/",iMult)) );
            kPlotFolder2D.push_back( TString( Form(kASigExtp_Plot_Direct,(TString("Multiplicity")+kFolder).Data()) ) +TString(Form("/MLT_%i/2D/",iMult)) );
            //
            h1D_Nraw_MT.at(iMult)        ->  Scale(1.,"width");
            //
            h1D_Nres_stat_MT.push_back( uEfficiencyCorrection1D ( h1D_Nraw_MT.at(iMult), h1D_Nrec, h1D_Ngen, f1DCorrection ) );
            h1D_Nres_syst_MT.push_back( fSetSystErrors ( h1D_Nres_stat_MT.at(iMult), "" ) );
            SetAxis(h1D_Nres_stat_MT.at(iMult),"PT1D");
            SetAxis(h1D_Nres_syst_MT.at(iMult),"PT1D");
            h1D_Nres_stat_MT.at(iMult)   ->  SetName(Form("h1D_Nres_stat_MT_%i",iMult));
            h1D_Nres_syst_MT.at(iMult)   ->  SetName(Form("h1D_Nres_syst_MT_%i",iMult));
            //
            kMult1D_Results.push_back( uMeasureFullYield( h1D_Nres_stat_MT.at(iMult),h1D_Nres_syst_MT.at(iMult),kStandardSystematicFitFunctions,kStatEvalCycles,kPlotFolder1D.at(iMult),"1D") );
            //
            h2D_Nraw_MT.at(iMult)        ->  Scale(1.,"width");
            h2D_Nraw_stat_MT.push_back( uEfficiencyCorrection2D ( h2D_Nraw_MT.at(iMult), h1D_Nrec_2Db, h1D_Ngen_2Db, f2DCorrection ) );
            h2D_Nraw_syst_MT.push_back( fSetSystErrors ( h2D_Nraw_stat_MT.at(iMult), "" ) );
            h2D_Nraw_stat_MT.at(iMult)   ->  SetName(Form("h2D_Nraw_stat_MT_%i",iMult));
            h2D_Nraw_syst_MT.at(iMult)   ->  SetName(Form("h2D_Nraw_syst_MT_%i",iMult));
            //
            kMult2D_Results.push_back( uMeasureFullYield2D( h2D_Nraw_stat_MT.at(iMult),h2D_Nraw_syst_MT.at(iMult),kStandardSystematicFitFunctions,kStatEvalCycles,kPlotFolder2D.at(iMult), "2D_%i" ) );
            //
        }
        //
        //  --- Build output plots
        //
        TGraphErrors*   g1D_Nres_Mult  =   new TGraphErrors ( );
        g1D_Nres_Mult->SetNameTitle( "g1D_Nres_Mult", "g1D_Nres_Mult" );
        TGraphErrors*   g2D_Nres_Mult  =   new TGraphErrors ( );
        g2D_Nres_Mult->SetNameTitle( "g2D_Nres_Mult", "g2D_Nres_Mult" );
        //
        auto iTer = -1;
        for ( auto kCurrent_1D_Result : kMult1D_Results )    {
            iTer++;
            if ( iTer <= 0 ) continue;
            g1D_Nres_Mult->SetPoint        ( iTer-1, fArrRMlt[iTer], get<0>( kCurrent_1D_Result["YL_FLL"] ) );
            g1D_Nres_Mult->SetPointError   ( iTer-1, 0,              get<1>( kCurrent_1D_Result["YL_FLL"] ) );
        }
        iTer = -1;
        for ( auto kCurrent_2D_Result : kMult2D_Results )    {
            iTer++;
            if ( iTer <= 0 ) continue;
            g2D_Nres_Mult->SetPoint        ( iTer-1, fArrRMlt[iTer], get<0>( kCurrent_2D_Result.at(0)["YL_FLL"] ) );
            g2D_Nres_Mult->SetPointError   ( iTer-1, 0,              get<1>( kCurrent_2D_Result.at(0)["YL_FLL"] ) );
        }
        //
        // --- Printing to Plots
        SetStyle();
        gROOT       ->  SetBatch(kTRUE);
        //
        auto cDraw1DMultNorm    =   new TCanvas( "cDraw1DMultNorm", "cDraw1DMultNorm", 1500, 1500 );
        //
        g1D_Nres_Mult      ->  Draw("EP");
        //
        cDraw1DMultNorm         ->  SaveAs( kPlotFolder1D_Full + TString("/Test.pdf") );
        cDraw1DMultNorm         ->  SaveAs( kPlotFolder1D_Full + TString("/Test.eps") );
        //
        delete  cDraw1DMultNorm;
        //
        gROOT           ->  SetBatch(kFALSE);
        //
        //  --- Output
        TFile *outFile_Rslt_ML      =   new TFile   (Form(kASigExtp_FitCheckRst,(TString("Multiplicity")+kFolder).Data()),"recreate");
        //
        g1D_Nres_Mult->Write();
        g2D_Nres_Mult->Write();
        for ( auto kSave : h1D_Nres_stat_MT )   kSave ->  Write();
        for ( auto kSave : h1D_Nres_syst_MT )   kSave ->  Write();
        for ( auto kSave : h2D_Nraw_stat_MT )   kSave ->  Write();
        for ( auto kSave : h2D_Nraw_syst_MT )   kSave ->  Write();
        //
        outFile_Rslt_ML->Close();
        outFile_Chck_ML->Close();
        insFile_Effc_ML->Close();
        insFile_Data_ML->Close();
         */
    }
}
