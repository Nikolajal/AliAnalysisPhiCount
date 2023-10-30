#include "../../inc/AliAnalysisPhiPair.h"

void AN_SigCorrections   ( TString fOption = "all", TString kFolder = "_p_p__7TeV", Bool_t fSilent = true )    {
    // --- --- --- --- --- --- --- SET-UP --- --- --- --- --- --- --- --- --- --- ---
    //
    //  --- INFO on Set-up variables
    fChooseOption(fOption);
    //
    //  --- Silencing warnings for smoother running
    if ( fSilent )  {
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
        RooMsgService::instance().setSilentMode(fSilent);
        gErrorIgnoreLevel   =   kWarning;
    }
    //
    //  --- Setting the input datastructure
    fSetAllBins();
    //
    //  --- Utility variables
    Int_t fTotalCount, fProgrCount;
    //
    // --- YIELD ANALYSIS
    if ( kDoYield ) {
        //  --- Load Files
        TFile*  insFile_Data_YL     =   new TFile   (Form(kASigExtr_FitCheckRst,(TString("Yield")+kFolder).Data()));
        TFile*  insFile_Effc_YL     =   new TFile   (Form(kAnalysis_MCTruthHist,(TString("Yield")+kFolder).Data()));
        //
        //  --- Load Histograms
        auto        fHEventCount    =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "fQC_Event_Enum_FLL"    );
        auto        h1D_Nraw        =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "h1D_Nraw_"             );
        auto        h1D_Nrec        =   uLoadHistograms<0,TH1F> ( insFile_Effc_YL,  "h1D_Nrec"              );
        auto        h1D_Ngen        =   uLoadHistograms<0,TH1F> ( insFile_Effc_YL,  "h1D_Ngen"              );
        auto        h2D_Nraw        =   uLoadHistograms<0,TH2F> ( insFile_Data_YL,  "anSS2D_",              "h2D_Nraw_" );
        auto        h1D_Nrec_2Db    =   uLoadHistograms<0,TH1F> ( insFile_Effc_YL,  "h1D_Nrec_2Db"          );
        auto        h1D_Ngen_2Db    =   uLoadHistograms<0,TH1F> ( insFile_Effc_YL,  "h1D_Ngen_2Db"          );
        TH1D*       hxD_YL_Extrapol_Syst =   new TH1D ( "hxD_YL_Extrapol_Syst", "hxD_YL_Extrapol_Syst", nBinPT2D+3, 0, nBinPT2D+3 );
        TH1D*       hxD_PT_Extrapol_Syst =   new TH1D ( "hxD_PT_Extrapol_Syst", "hxD_PT_Extrapol_Syst", nBinPT2D+3, 0, nBinPT2D+3 );
        //
        //  --- Build output structure
        gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_SigExtp_Dir,(TString("Yield")+kFolder).Data())));
        gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,(TString("Yield")+kFolder).Data())));
        gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,(TString("Yield")+kFolder).Data()))+TString("/1D"));
        gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,(TString("Yield")+kFolder).Data()))+TString("/2D"));
        gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,(TString("Yield")+kFolder).Data()))+TString("/Full"));
        auto        kPlotFolder1D   =   TString( Form(kASigExtp_Plot_Direct,(TString("Yield")+kFolder).Data()) ) + TString( "/1D/" );
        auto        kPlotFolder2D   =   TString( Form(kASigExtp_Plot_Direct,(TString("Yield")+kFolder).Data()) ) + TString( "/2D/" );
        //  --- --- Build Output Histogram
        TH1F*   hXD_Nyld_stat       =   new TH1F ( "hXD_Nyld_stat", "hXD_Nyld_stat", 5, 0.5, 5.5 ); // Yield 1D, Yield 2D, PT 1D, EXT 1D, EXT 2D
        TH1F*   hXD_Nyld_syst       =   new TH1F ( "hXD_Nyld_syst", "hXD_Nyld_syst", 5, 0.5, 5.5 ); // --
        TH1F*   hXD_Nfqs_stat       =   new TH1F ( "hXD_Nfqs_stat", "hXD_Nfqs_stat", 6, 0.5, 6.5 ); // Yield 1D, Yield 2D, R1, R2, P1, P2
        TH1F*   hXD_Nfqs_syst       =   new TH1F ( "hXD_Nfqs_syst", "hXD_Nfqs_syst", 6, 0.5, 6.5 ); // --
        //
        //  --- Output File for Fit Check
        TFile*  outFile_Chck_YL     =   new TFile(Form(kASigExtp_FitCheckPlt,(TString("Yield")+kFolder).Data()),"recreate");
        //
        //  --- Minimum Bias Normalisation
        Double_t    f1DCorrection   =   -1;
        Double_t    f2DCorrection   =   -1;
        if ( is_pp_anl ) {
            auto        kN_PU           =   (fHEventCount->GetBinContent(kEventCount::kPU_MB));
            auto        kN_Trg          =   (fHEventCount->GetBinContent(kEventCount::kTrigger));
            auto        kN_Vtx          =   (fHEventCount->GetBinContent(kEventCount::kVertex));
            auto        kN_MB           =   -kN_PU +(fHEventCount->GetBinContent(kEventCount::kVertex10));
            f1DCorrection               =   (1./kBR)        *(1./kN_MB) *(kN_Vtx/kN_Trg);
            f2DCorrection               =   (1./(kBR*kBR))  *(1./kN_MB) *(kN_Vtx/kN_Trg);
            if ( fabs( kEnergy - 5  ) < 0.1 ) { f1DCorrection   *=  kTriggerEff15n17pq; f2DCorrection   *=  kTriggerEff15n17pq; }
            if ( fabs( kEnergy - 7  ) < 0.1 ) { f1DCorrection   *=  kTriggerEff10bcdef; f2DCorrection   *=  kTriggerEff10bcdef; }
            if ( fabs( kEnergy - 13 ) < 0.1 ) { f1DCorrection   *=  kTriggerEff13TeV;   f2DCorrection   *=  kTriggerEff13TeV;   }
        } else if ( is_pb_anl ) {
            auto        kN_PU           =   (fHEventCount->GetBinContent(kEventCount::kPU_MB));
            auto        kN_Trg          =   (fHEventCount->GetBinContent(kEventCount::kTrigger));
            auto        kN_Vtx          =   (fHEventCount->GetBinContent(kEventCount::kVertex));
            auto        kN_MB           =   -kN_PU +(fHEventCount->GetBinContent(kEventCount::kVertex10));
            f1DCorrection               =   (1./kBR)        *(1./kN_MB) *(kN_Vtx/kN_Trg) *(1./0.5); // y normalisation 0.5
            f2DCorrection               =   (1./(kBR*kBR))  *(1./kN_MB) *(kN_Vtx/kN_Trg) *(1./0.5); // y normalisation 0.5
            if ( fabs( kEnergy - 5  ) < 0.1 ) {   f1DCorrection   *=  0.978;     f2DCorrection       *=  0.978; }
        }
        //
                    h1D_Nraw        ->  Scale(1.,"width");
        //
        TH1F*       h1D_Nraw_stat   =   uEfficiencyCorrection1D ( h1D_Nraw, h1D_Nrec, h1D_Ngen, f1DCorrection );
        TH1F*       h1D_Nraw_syst   =   fSetSystErrors ( h1D_Nraw_stat, Form("%s/FullSystematics.root",Form(kAnalysis_Systemt_Dir,  (TString("Yield")+kFolder).Data())), "hFullSystematics1D" );
        SetAxis(h1D_Nraw_stat,"PT1D");
        SetAxis(h1D_Nraw_syst,"PT1D");
        h1D_Nraw_stat   ->  SetName("h1D_Nres_stat");
        h1D_Nraw_syst   ->  SetName("h1D_Nres_syst");
        //
        fSetAllCustomFunctions();
        //
        auto        k1D_Results     =   uMeasureFullYield(h1D_Nraw_stat,h1D_Nraw_syst,kStandardSystematicFitFunctions,kStatEvalCycles,kPlotFolder1D,"1D");
                    hXD_Nyld_stat   ->  SetBinContent   ( 1, get<0>( k1D_Results["YL_FLL"] ) );
                    hXD_Nyld_stat   ->  SetBinError     ( 1, get<1>( k1D_Results["YL_FLL"] ) );
                    hXD_Nyld_syst   ->  SetBinContent   ( 1, get<0>( k1D_Results["YL_FLL"] ) );
                    hXD_Nyld_syst   ->  SetBinError     ( 1, get<2>( k1D_Results["YL_FLL"] ) );
        //
                    hXD_Nyld_stat   ->  SetBinContent   ( 3, get<0>( k1D_Results["PT_FLL"] ) );
                    hXD_Nyld_stat   ->  SetBinError     ( 3, get<1>( k1D_Results["PT_FLL"] ) );
                    hXD_Nyld_syst   ->  SetBinContent   ( 3, get<0>( k1D_Results["PT_FLL"] ) );
                    hXD_Nyld_syst   ->  SetBinError     ( 3, get<2>( k1D_Results["PT_FLL"] ) );
        //
                    hXD_Nyld_stat   ->  SetBinContent   ( 4, get<0>( k1D_Results["YL_EXT"] ) );
                    hXD_Nyld_stat   ->  SetBinError     ( 4, get<1>( k1D_Results["YL_EXT"] ) );
                    hXD_Nyld_syst   ->  SetBinContent   ( 4, get<0>( k1D_Results["YL_EXT"] ) );
                    hXD_Nyld_syst   ->  SetBinError     ( 4, get<2>( k1D_Results["YL_EXT"] ) );
        //
                    hxD_YL_Extrapol_Syst    ->  SetBinContent   ( 1, get<0>( k1D_Results["YL_FIT"] ) );
                    hxD_YL_Extrapol_Syst    ->  SetBinError     ( 1, get<1>( k1D_Results["YL_FIT"] ) );
                    hxD_PT_Extrapol_Syst    ->  SetBinContent   ( 1, get<0>( k1D_Results["PT_FIT"] ) );
                    hxD_PT_Extrapol_Syst    ->  SetBinError     ( 1, get<1>( k1D_Results["PT_FIT"] ) );
        //
                    h2D_Nraw        ->  Scale(1.,"width");
        //
        auto        h2D_Nraw_stat   =   uEfficiencyCorrection2D ( h2D_Nraw, h1D_Nrec_2Db, h1D_Ngen_2Db, f2DCorrection );
        TH2F*       h2D_Nraw_syst   =   fSetSystErrors ( h2D_Nraw_stat, Form("%s/FullSystematics.root",Form(kAnalysis_Systemt_Dir,  (TString("Yield")+kFolder).Data())), "hFullSystematics2D" );
                    h2D_Nraw_stat   ->  SetName("h2D_Nres_stat");
                    h2D_Nraw_syst   ->  SetName("h2D_Nres_syst");
        //
        auto        k2D_Results     =   uMeasureFullYield2D( h2D_Nraw_stat, h2D_Nraw_syst, kStandardSystematicFitFunctions,kStatEvalCycles,kPlotFolder2D,"2D_%i");
                    hXD_Nyld_stat   ->  SetBinContent   ( 2, get<0>( k2D_Results.at(0)["YL_FLL"] ) );
                    hXD_Nyld_stat   ->  SetBinError     ( 2, get<1>( k2D_Results.at(0)["YL_FLL"] ) );
                    hXD_Nyld_syst   ->  SetBinContent   ( 2, get<0>( k2D_Results.at(0)["YL_FLL"] ) );
                    hXD_Nyld_syst   ->  SetBinError     ( 2, get<2>( k2D_Results.at(0)["YL_FLL"] ) );
        //
                    hXD_Nyld_stat   ->  SetBinContent   ( 5, get<0>( k2D_Results.at(1)["YL_EXT"] ) );
                    hXD_Nyld_stat   ->  SetBinError     ( 5, get<1>( k2D_Results.at(1)["YL_EXT"] ) );
                    hXD_Nyld_syst   ->  SetBinContent   ( 5, get<0>( k2D_Results.at(1)["YL_EXT"] ) );
                    hXD_Nyld_syst   ->  SetBinError     ( 5, get<2>( k2D_Results.at(1)["YL_EXT"] ) );
        //
        //  --- Conditional Spectra
        std::vector<TH1D*>  kStat_Array;
        std::vector<TH1D*>  kSyst_Array;
        for ( Int_t iBin = 1; iBin <= h2D_Nraw_stat->GetNbinsX(); iBin++ )   kStat_Array.push_back( (TH1D*)(h2D_Nraw_stat->ProjectionX(Form("%s_stat_%i",h2D_Nraw_stat->GetName(),iBin),iBin,iBin))->Clone() );
        for ( Int_t iBin = 1; iBin <= h2D_Nraw_syst->GetNbinsX(); iBin++ )   kSyst_Array.push_back( (TH1D*)(h2D_Nraw_syst->ProjectionX(Form("%s_syst_%i",h2D_Nraw_syst->GetName(),iBin),iBin,iBin))->Clone() );
        //
        auto    iBin    = -1;
        TH1D*   h2D_LowPT_Extrap_stat   =   (TH1D*)(h2D_Nraw_stat->ProjectionX(Form("%s_stat_%i",h2D_Nraw_stat->GetName(),0),1,1))->Clone();
        TH1D*   h2D_LowPT_Extrap_syst   =   (TH1D*)(h2D_Nraw_syst->ProjectionX(Form("%s_syst_%i",h2D_Nraw_syst->GetName(),0),1,1))->Clone();
        TH1D*   h2D_MeanPT_stat         =   new TH1D ( "h2D_MeanPT_stat", "h2D_MeanPT_stat", nBinPT2D+1, fArrPT2D_Comp);
        TH1D*   h2D_MeanPT_syst         =   new TH1D ( "h2D_MeanPT_syst", "h2D_MeanPT_syst", nBinPT2D+1, fArrPT2D_Comp);
        //
        for ( auto kResult : k2D_Results )    {
            iBin++;
            //
            if ( iBin <= 0 ) continue;
            //
            auto kMPTValue  =   get<0>(kResult["PT_FLL"]);
            auto kMPT_stat  =   get<1>(kResult["PT_FLL"]);
            auto kMPT_syst  =   get<2>(kResult["PT_FLL"]);
            h2D_MeanPT_stat         ->  SetBinContent   ( iBin, kMPTValue );
            h2D_MeanPT_stat         ->  SetBinError     ( iBin, kMPT_stat );
            h2D_MeanPT_syst         ->  SetBinContent   ( iBin, kMPTValue );
            h2D_MeanPT_syst         ->  SetBinError     ( iBin, kMPT_syst );
            //
            auto kBinVYLFT  =   get<0>(kResult["YL_FIT"]);
            auto kBin_YLFT  =   get<1>(kResult["YL_FIT"]);
            auto kBinVPTFT  =   get<0>(kResult["PT_FIT"]);
            auto kBin_PTFT  =   get<1>(kResult["PT_FIT"]);
            hxD_YL_Extrapol_Syst    ->  SetBinContent   ( iBin+2, kBinVYLFT );
            hxD_YL_Extrapol_Syst    ->  SetBinError     ( iBin+2, kBin_YLFT );
            hxD_PT_Extrapol_Syst    ->  SetBinContent   ( iBin+2, kBinVPTFT );
            hxD_PT_Extrapol_Syst    ->  SetBinError     ( iBin+2, kBin_PTFT );
            //
            if ( iBin <= 1 ) continue;
            //
            auto kBinValue  =   get<0>(kResult["YL_EXT"]);
            auto kBin_stat  =   get<1>(kResult["YL_EXT"]);
            auto kBin_syst  =   get<2>(kResult["YL_EXT"]);
            h2D_LowPT_Extrap_stat   ->  SetBinContent   ( iBin-1, kBinValue );
            h2D_LowPT_Extrap_stat   ->  SetBinError     ( iBin-1, kBin_stat );
            h2D_LowPT_Extrap_syst   ->  SetBinContent   ( iBin-1, kBinValue );
            h2D_LowPT_Extrap_syst   ->  SetBinError     ( iBin-1, kBin_syst );
            //
        }
        h2D_LowPT_Extrap_stat       ->  Scale( 1., "width" );
        h2D_LowPT_Extrap_syst       ->  Scale( 1., "width" );
        push_to_front( kStat_Array, h2D_LowPT_Extrap_stat );
        push_to_front( kSyst_Array, h2D_LowPT_Extrap_syst );
        //
        h2D_MeanPT_syst = fSetSystErrors ( h2D_MeanPT_syst, Form("%s/FullSystematics.root",Form(kAnalysis_Systemt_Dir,  (TString("Yield")+kFolder).Data())), "hFullSystematicsPT" );
        //
        //  --- Final Quantities
        hXD_Nfqs_stat   =   uPlotDerivedQuantitiesRaw( {hXD_Nyld_stat->GetBinContent(1),hXD_Nyld_stat->GetBinError(1),hXD_Nyld_stat->GetBinError(1)}, {hXD_Nyld_stat->GetBinContent(2),hXD_Nyld_stat->GetBinError(2),hXD_Nyld_stat->GetBinError(2)} );
        hXD_Nfqs_syst   =   hXD_Nfqs_stat;
        hXD_Nfqs_syst   =   fSetSystErrors ( hXD_Nfqs_syst, Form("%s/FullSystematics.root",Form(kAnalysis_Systemt_Dir,  (TString("Yield")+kFolder).Data())), "hFullSystematicsRT" );
        hXD_Nfqs_stat   ->  SetName("hXD_Nfqs_stat");
        hXD_Nfqs_syst   ->  SetName("hXD_Nfqs_syst");
        //
        //  --- Output
        TFile *outFile_Rslt_YL      =   new TFile   (Form(kASigExtp_FitCheckRst,(TString("Yield")+kFolder).Data()),"recreate");
        hXD_Nyld_stat               ->  Write();
        hXD_Nyld_syst               ->  Write();
        h1D_Nraw_stat               ->  Write();
        h1D_Nraw_syst               ->  Write();
        h2D_Nraw_stat               ->  Write();
        h2D_Nraw_syst               ->  Write();
        hxD_YL_Extrapol_Syst        ->  Write();
        hxD_PT_Extrapol_Syst        ->  Write();
        hXD_Nfqs_syst               ->  Write();
        hXD_Nfqs_stat               ->  Write();
        h2D_MeanPT_stat             ->  Write();
        h2D_MeanPT_syst             ->  Write();
        for ( auto kSave : kStat_Array ) kSave   ->  Write();
        for ( auto kSave : kSyst_Array ) kSave   ->  Write();
        //
        // --- Printing to Plots
        SetStyle();
        gROOT       ->  SetBatch(kTRUE);
        //
        auto cDrawResult = uPlotSpectrum(h1D_Nraw_stat,h1D_Nraw_syst,"SPT 1D T");
        cDrawResult ->  SaveAs(kPlotFolder1D+TString("/BeauSpectrum_YLD_1D.pdf"));
        cDrawResult ->  SaveAs(kPlotFolder1D+TString("/BeauSpectrum_YLD_1D.eps"));
        delete      cDrawResult;
        //
        cDrawResult = uPlotSpectrum(h2D_MeanPT_stat,h2D_MeanPT_syst,"1D R MPT ");
        cDrawResult ->  SaveAs(kPlotFolder2D+TString("/BeauSpectrum_MPT_2D.pdf"));
        cDrawResult ->  SaveAs(kPlotFolder2D+TString("/BeauSpectrum_MPT_2D.eps"));
        delete      cDrawResult;
        //
        auto iPT = -1;
        for ( auto kSave : kStat_Array )    {
            iPT++;
            cDrawResult = uPlotSpectrum(kStat_Array.at(iPT),kSyst_Array.at(iPT),"SPT 12D T ");
            uLatex      ->  DrawLatexNDC( 0.28, 0.83, Form("#it{p}_{T,#phi_{2}} #in [%.2f-%.2f] GeV/#it{c}",fArrPT2D_Comp[iPT],fArrPT2D_Comp[iPT+1]));
            cDrawResult ->  SaveAs(kPlotFolder2D+TString(Form("/BeauSpectrum_YLD_2D_PT%i.pdf",iPT)));
            cDrawResult ->  SaveAs(kPlotFolder2D+TString(Form("/BeauSpectrum_YLD_2D_PT%i.eps",iPT)));
            delete      cDrawResult;
        }
        //
        cDrawResult = uPlotSpectrum(h2D_Nraw_stat,h2D_Nraw_syst,"SPT 2D T");
        cDrawResult ->  SaveAs(kPlotFolder2D+TString("/BeauSpectrum_YLD_2D.pdf"));
        cDrawResult ->  SaveAs(kPlotFolder2D+TString("/BeauSpectrum_YLD_2D.eps"));
        delete      cDrawResult;
        //
        gROOT           ->  SetBatch(kFALSE);
        //
        //  --- Close Files
        //
        outFile_Rslt_YL->Close();
        outFile_Chck_YL->Close();
        insFile_Effc_YL->Close();
        insFile_Data_YL->Close();
    }
    // --- MULTIPLICITY ANALYSIS
    if ( kDoMultiplicity )  {
        //  --- Load Files
        TFile*  insFile_Data_ML     =   new TFile   (Form(kASigExtr_FitCheckRst,(TString("Multiplicity")+kFolder).Data()));
        TFile*  insFile_Effc_ML     =   new TFile   (Form(kAnalysis_MCTruthHist,(TString("Multiplicity")+kFolder).Data()));
        //
        //  --- Load Histograms
        auto        fHEventCount    =   uLoadHistograms<0,TH1F> ( insFile_Data_ML,  "fQC_Event_Enum_FLL"     );
        auto        fHEventCntMlt   =   uLoadHistograms<0,TH1F> ( insFile_Data_ML,  "fQC_Event_Enum_V0M"     );
        auto        h1D_Nraw_MT     =   uLoadHistograms<1,TH1F> ( insFile_Data_ML,  "h1D_Nraw_MT_%i",       "h1D_Nraw_MT_%i" );
        auto        h1D_Nrec        =   uLoadHistograms<0,TH1F> ( insFile_Effc_ML,  "h1D_Nrec" );
        auto        h1D_Ngen        =   uLoadHistograms<0,TH1F> ( insFile_Effc_ML,  "h1D_Ngen" );
        auto        h2D_Nraw_MT     =   uLoadHistograms<1,TH2F> ( insFile_Data_ML,  "anSS2D_MT_%i",         "h2D_Nraw_MT_%i" );
        auto        h1D_Nrec_2Db    =   uLoadHistograms<0,TH1F> ( insFile_Effc_ML,  "h1D_Nrec_2Db" );
        auto        h1D_Ngen_2Db    =   uLoadHistograms<0,TH1F> ( insFile_Effc_ML,  "h1D_Ngen_2Db" );
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
        std::vector<TH1D*>      h2D_MnPt_stat_MT;
        std::vector<TH1D*>      h2D_MnPt_syst_MT;
        std::vector<TH1F*>      hXD_Nfqs_stat;
        std::vector<TH1F*>      hXD_Nfqs_syst;
        std::vector<std::vector<TH1D*>> h2D_Slcs_stat_MT;
        std::vector<std::vector<TH1D*>> h2D_Slcs_syst_MT;
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
            //  --- Multiplicity Normalisation
            Double_t    f1DCorrection   =   -1;
            Double_t    f2DCorrection   =   -1;
            if ( is_pp_anl ) {
                f1DCorrection   =   (1.)/(fEvaluateINELgt0(iMult-1,fHEventCntMlt) * kBR );
                f2DCorrection   =   (1.)/(fEvaluateINELgt0(iMult-1,fHEventCntMlt) * kBR * kBR );
            } else if ( is_pb_anl ) {
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
            h1D_Nres_syst_MT.push_back( fSetSystErrors ( h1D_Nres_stat_MT.at(iMult), Form("%s/FullSystematics.root",Form(kAnalysis_Systemt_Dir+TString(Form("/MLT_%i/",iMult)),(TString("Multiplicity")+kFolder).Data())), "hFullSystematics1D" ));
            SetAxis(h1D_Nres_stat_MT.at(iMult),"PT1D");
            SetAxis(h1D_Nres_syst_MT.at(iMult),"PT1D");
            h1D_Nres_stat_MT.at(iMult)   ->  SetName(Form("h1D_Nres_stat_MT_%i",iMult));
            h1D_Nres_syst_MT.at(iMult)   ->  SetName(Form("h1D_Nres_syst_MT_%i",iMult));
            //
            kMult1D_Results.push_back( uMeasureFullYield( h1D_Nres_stat_MT.at(iMult),h1D_Nres_syst_MT.at(iMult),kStandardSystematicFitFunctions,kStatEvalCycles,kPlotFolder1D.at(iMult),"1D") );
            //
            h2D_Nraw_MT.at(iMult)        ->  Scale(1.,"width");
            h2D_Nraw_stat_MT.push_back( uEfficiencyCorrection2D ( h2D_Nraw_MT.at(iMult), h1D_Nrec_2Db, h1D_Ngen_2Db, f2DCorrection ) );
                                       h2D_Nraw_syst_MT.push_back( fSetSystErrors ( h2D_Nraw_stat_MT.at(iMult), Form("%s/FullSystematics.root",Form(kAnalysis_Systemt_Dir+TString(Form("/MLT_%i/",iMult)),(TString("Multiplicity")+kFolder).Data())), "hFullSystematics2D"));
            h2D_Nraw_stat_MT.at(iMult)   ->  SetName(Form("h2D_Nres_stat_MT_%i",iMult));
            h2D_Nraw_syst_MT.at(iMult)   ->  SetName(Form("h2D_Nres_syst_MT_%i",iMult));
            //
            kMult2D_Results.push_back( uMeasureFullYield2D( h2D_Nraw_stat_MT.at(iMult),h2D_Nraw_syst_MT.at(iMult),kStandardSystematicFitFunctions,kStatEvalCycles,kPlotFolder2D.at(iMult), "2D_%i" ) );
            //
            //  --- Conditional Spectra
            std::vector<TH1D*>  k2DProjectionArrayStat;
            std::vector<TH1D*>  k2DProjectionArraySyst;
            for ( Int_t iBin = 1; iBin <= h2D_Nraw_stat_MT.at(iMult)->GetNbinsX(); iBin++ ) k2DProjectionArrayStat.push_back( (TH1D*)(h2D_Nraw_stat_MT.at(iMult)->ProjectionX(Form("%s_stat_%i",h2D_Nraw_stat_MT.at(iMult)->GetName(),iBin),iBin,iBin))->Clone() );
            for ( Int_t iBin = 1; iBin <= h2D_Nraw_syst_MT.at(iMult)->GetNbinsX(); iBin++ ) k2DProjectionArraySyst.push_back( (TH1D*)(h2D_Nraw_syst_MT.at(iMult)->ProjectionX(Form("%s_syst_%i",h2D_Nraw_syst_MT.at(iMult)->GetName(),iBin),iBin,iBin))->Clone() );
            //
            auto    iBin    = -1;
            TH1D*   h2D_LowPT_Extrap_stat   =   (TH1D*)(h2D_Nraw_stat_MT.at(iMult)->ProjectionX(Form("%s_stat_%i",h2D_Nraw_stat_MT.at(iMult)->GetName(),0),1,1))->Clone();
            TH1D*   h2D_LowPT_Extrap_syst   =   (TH1D*)(h2D_Nraw_syst_MT.at(iMult)->ProjectionX(Form("%s_syst_%i",h2D_Nraw_syst_MT.at(iMult)->GetName(),0),1,1))->Clone();
            TH1D*   h2D_MeanPT_stat         =   new TH1D ( Form("h2D_MeanPT_stat_MT_%i",iMult), "h2D_MeanPT_stat", nBinPT2D+1, fArrPT2D_Comp);
            TH1D*   h2D_MeanPT_syst         =   new TH1D ( Form("h2D_MeanPT_syst_MT_%i",iMult), "h2D_MeanPT_syst", nBinPT2D+1, fArrPT2D_Comp);
            //
            for ( auto kResult : kMult2D_Results.at(iMult) )    {
                iBin++;
                //
                if ( iBin <= 0 ) continue;
                //
                auto kMPTValue  =   get<0>(kResult["PT_FLL"]);
                auto kMPT_stat  =   get<1>(kResult["PT_FLL"]);
                auto kMPT_syst  =   get<2>(kResult["PT_FLL"]);
                h2D_MeanPT_stat         ->  SetBinContent   ( iBin, kMPTValue );
                h2D_MeanPT_stat         ->  SetBinError     ( iBin, kMPT_stat );
                h2D_MeanPT_syst         ->  SetBinContent   ( iBin, kMPTValue );
                h2D_MeanPT_syst         ->  SetBinError     ( iBin, kMPT_syst );
                //
                if ( iBin <= 1 ) continue;
                //
                auto kBinValue  =   get<0>(kResult["YL_EXT"]);
                auto kBin_stat  =   get<1>(kResult["YL_EXT"]);
                auto kBin_syst  =   get<2>(kResult["YL_EXT"]);
                h2D_LowPT_Extrap_stat   ->  SetBinContent   ( iBin-1, kBinValue );
                h2D_LowPT_Extrap_stat   ->  SetBinError     ( iBin-1, kBin_stat );
                h2D_LowPT_Extrap_syst   ->  SetBinContent   ( iBin-1, kBinValue );
                h2D_LowPT_Extrap_syst   ->  SetBinError     ( iBin-1, kBin_syst );
                //
            }
            h2D_LowPT_Extrap_stat       ->  Scale( 1., "width" );
            h2D_LowPT_Extrap_syst       ->  Scale( 1., "width" );
            push_to_front( k2DProjectionArrayStat, h2D_LowPT_Extrap_stat );
            push_to_front( k2DProjectionArraySyst, h2D_LowPT_Extrap_syst );
            //
            fSetSystErrors ( h2D_MeanPT_syst, Form("%s/FullSystematics.root",Form(kAnalysis_Systemt_Dir+TString(Form("/MLT_%i/",iMult)),(TString("Multiplicity")+kFolder).Data())), "hFullSystematicsPT");
            //
            h2D_MnPt_stat_MT.push_back(h2D_MeanPT_stat);
            h2D_MnPt_syst_MT.push_back(h2D_MeanPT_syst);
            h2D_Slcs_stat_MT.push_back(k2DProjectionArrayStat);
            h2D_Slcs_syst_MT.push_back(k2DProjectionArraySyst);
        }
        //
        //  --- Build output plots
        TGraphErrors*   g1D_Nres_Stat_Mult   =   new TGraphErrors ( );
        TGraphErrors*   g2D_Nres_Stat_Mult   =   new TGraphErrors ( );
        TGraphErrors*   gR1_Nres_Stat_Mult   =   new TGraphErrors ( );
        TGraphErrors*   gR2_Nres_Stat_Mult   =   new TGraphErrors ( );
        TGraphErrors*   gP1_Nres_Stat_Mult   =   new TGraphErrors ( );
        TGraphErrors*   gP2_Nres_Stat_Mult   =   new TGraphErrors ( );
        TGraphErrors*   g1D_Nres_Syst_Mult   =   new TGraphErrors ( );
        TGraphErrors*   g2D_Nres_Syst_Mult   =   new TGraphErrors ( );
        TGraphErrors*   gR1_Nres_Syst_Mult   =   new TGraphErrors ( );
        TGraphErrors*   gR2_Nres_Syst_Mult   =   new TGraphErrors ( );
        TGraphErrors*   gP1_Nres_Syst_Mult   =   new TGraphErrors ( );
        TGraphErrors*   gP2_Nres_Syst_Mult   =   new TGraphErrors ( );
        g1D_Nres_Stat_Mult   ->  SetNameTitle( "g1D_Nres_Stat_Mult", "g1D_Nres_Stat_Mult" );
        g2D_Nres_Stat_Mult   ->  SetNameTitle( "g2D_Nres_Stat_Mult", "g2D_Nres_Stat_Mult" );
        gR1_Nres_Stat_Mult   ->  SetNameTitle( "gR1_Nres_Stat_Mult", "gR1_Nres_Stat_Mult" );
        gR2_Nres_Stat_Mult   ->  SetNameTitle( "gR2_Nres_Stat_Mult", "gR2_Nres_Stat_Mult" );
        gP1_Nres_Stat_Mult   ->  SetNameTitle( "gP1_Nres_Stat_Mult", "gP1_Nres_Stat_Mult" );
        gP2_Nres_Stat_Mult   ->  SetNameTitle( "gP2_Nres_Stat_Mult", "gP2_Nres_Stat_Mult" );
        g1D_Nres_Stat_Mult   ->  SetNameTitle( "g1D_Nres_Stat_Mult", "g1D_Nres_Stat_Mult" );
        g2D_Nres_Stat_Mult   ->  SetNameTitle( "g2D_Nres_Stat_Mult", "g2D_Nres_Stat_Mult" );
        g1D_Nres_Syst_Mult   ->  SetNameTitle( "g1D_Nres_Syst_Mult", "g1D_Nres_Syst_Mult" );
        g2D_Nres_Syst_Mult   ->  SetNameTitle( "g2D_Nres_Syst_Mult", "g2D_Nres_Syst_Mult" );
        gR1_Nres_Syst_Mult   ->  SetNameTitle( "gR1_Nres_Syst_Mult", "gR1_Nres_Syst_Mult" );
        gR2_Nres_Syst_Mult   ->  SetNameTitle( "gR2_Nres_Syst_Mult", "gR2_Nres_Syst_Mult" );
        gP1_Nres_Syst_Mult   ->  SetNameTitle( "gP1_Nres_Syst_Mult", "gP1_Nres_Syst_Mult" );
        gP2_Nres_Syst_Mult   ->  SetNameTitle( "gP2_Nres_Syst_Mult", "gP2_Nres_Syst_Mult" );
        g1D_Nres_Syst_Mult   ->  SetNameTitle( "g1D_Nres_Syst_Mult", "g1D_Nres_Syst_Mult" );
        g2D_Nres_Syst_Mult   ->  SetNameTitle( "g2D_Nres_Syst_Mult", "g2D_Nres_Syst_Mult" );
        //
        auto iTer = -1;
        for ( auto kCurrent_1D_Result : kMult1D_Results )    {
            iTer++;
            auto    kYield_1D   =   get<0>( kCurrent_1D_Result["YL_FLL"] );
            auto    kError_1D   =   get<1>( kCurrent_1D_Result["YL_FLL"] );
            auto    kYield_2D   =   get<0>( ((kMult2D_Results.at(iTer)).at(0))["YL_FLL"] );
            auto    kError_2D   =   get<1>( ((kMult2D_Results.at(iTer)).at(0))["YL_FLL"] );
            //
            hXD_Nfqs_stat.push_back( new TH1F ( Form("hXD_Nfqs_stat_MT_%i",iTer), "hXD_Nfqs_stat", 6, 0.5, 6.5 ) );
            hXD_Nfqs_syst.push_back( new TH1F ( Form("hXD_Nfqs_syst_MT_%i",iTer), "hXD_Nfqs_syst", 6, 0.5, 6.5 ) ); // Yield 1D, Yield 2D, R1, R2, P1, P2
            //
            hXD_Nfqs_stat.at(iTer) = uPlotDerivedQuantitiesRaw( {kYield_1D, kError_1D, 0.}, {kYield_2D, kError_2D, 0.} );
            hXD_Nfqs_syst.at(iTer) = hXD_Nfqs_stat.at(iTer);
            hXD_Nfqs_syst.at(iTer) = fSetSystErrors ( hXD_Nfqs_syst.at(iTer), Form("%s/FullSystematics.root",Form(kAnalysis_Systemt_Dir+TString(Form("/MLT_%i/",iTer)),(TString("Multiplicity")+kFolder).Data())), "hFullSystematicsRT" );
            hXD_Nfqs_stat.at(iTer) -> SetName(Form("hXD_Nfqs_stat_MT_%i",iTer));
            hXD_Nfqs_syst.at(iTer) -> SetName(Form("hXD_Nfqs_syst_MT_%i",iTer));
            auto    kErrSy_1D   =   hXD_Nfqs_syst.at(iTer)->GetBinError(1);
            auto    kErrSy_2D   =   hXD_Nfqs_syst.at(iTer)->GetBinError(2);
            auto    kErRel_1D   =   kError_1D/kYield_1D;
            auto    kErRel_2D   =   kError_2D/kYield_2D;
            auto    kErRSys1D   =   kErrSy_1D/kYield_1D;
            auto    kErRSys2D   =   kErrSy_2D/kYield_2D;
            //
            if ( iTer <= 0 ) continue;
            //
            g1D_Nres_Stat_Mult->SetPoint        ( iTer-1, fArrRMlt[iTer], kYield_1D );
            g1D_Nres_Stat_Mult->SetPointError   ( iTer-1, 0,              kError_1D );
            //
            g2D_Nres_Stat_Mult->SetPoint        ( iTer-1, fArrRMlt[iTer], kYield_2D );
            g2D_Nres_Stat_Mult->SetPointError   ( iTer-1, 0,              kError_2D );
            //
            gR1_Nres_Stat_Mult->SetPoint        ( iTer-1, fArrRMlt[iTer], kYield_2D/kYield_1D );
            gR1_Nres_Stat_Mult->SetPointError   ( iTer-1, 0,              (kYield_2D/kYield_1D)*SquareSum( { kErRel_2D, kErRel_1D } ) );
            //
            gR2_Nres_Stat_Mult->SetPoint        ( iTer-1, fArrRMlt[iTer], kYield_2D/(kYield_1D*kYield_1D) );
            gR2_Nres_Stat_Mult->SetPointError   ( iTer-1, 0,              (kYield_2D/(kYield_1D*kYield_1D))*SquareSum( { 2*kErRel_2D, kErRel_1D } ) );
            //
            gP1_Nres_Stat_Mult->SetPoint        ( iTer-1, fArrRMlt[iTer], fSigmaPhiValue ( kYield_1D, kYield_2D ) );
            gP1_Nres_Stat_Mult->SetPointError   ( iTer-1, 0,              fSigmaPhiError ( kYield_1D, kYield_2D, kError_1D, kError_2D ) );
            //
            gP2_Nres_Stat_Mult->SetPoint        ( iTer-1, fArrRMlt[iTer], fGammaPhiValue ( kYield_1D, kYield_2D ) );
            gP2_Nres_Stat_Mult->SetPointError   ( iTer-1, 0,              fGammaPhiError ( kYield_1D, kYield_2D, kError_1D, kError_2D ) );
            //
            g1D_Nres_Syst_Mult->SetPoint        ( iTer-1, fArrRMlt[iTer], kYield_1D );
            g1D_Nres_Syst_Mult->SetPointError   ( iTer-1, 0,              kErrSy_1D );
            //
            g2D_Nres_Syst_Mult->SetPoint        ( iTer-1, fArrRMlt[iTer], kYield_2D );
            g2D_Nres_Syst_Mult->SetPointError   ( iTer-1, 0,              kErrSy_2D );
            //
            gR1_Nres_Syst_Mult->SetPoint        ( iTer-1, fArrRMlt[iTer], kYield_2D/kYield_1D );
            gR1_Nres_Syst_Mult->SetPointError   ( iTer-1, 0,              (kYield_2D/kYield_1D)*SquareSum( { kErRSys2D, kErRSys1D } ) );
            //
            gR2_Nres_Syst_Mult->SetPoint        ( iTer-1, fArrRMlt[iTer], kYield_2D/(kYield_1D*kYield_1D) );
            gR2_Nres_Syst_Mult->SetPointError   ( iTer-1, 0,              (kYield_2D/(kYield_1D*kYield_1D))*SquareSum( { 2*kErRSys2D, kErRSys1D } ) );
            //
            gP1_Nres_Syst_Mult->SetPoint        ( iTer-1, fArrRMlt[iTer], fSigmaPhiValue ( kYield_1D, kYield_2D ) );
            gP1_Nres_Syst_Mult->SetPointError   ( iTer-1, 0,              fSigmaPhiError ( kYield_1D, kYield_2D, kErrSy_1D, kErrSy_2D ) );
            //
            gP2_Nres_Syst_Mult->SetPoint        ( iTer-1, fArrRMlt[iTer], fGammaPhiValue ( kYield_1D, kYield_2D ) );
            gP2_Nres_Syst_Mult->SetPointError   ( iTer-1, 0,              fGammaPhiError ( kYield_1D, kYield_2D, kErrSy_1D, kErrSy_2D ) );
        }
        //  --- Output
        TFile *outFile_Rslt_ML      =   new TFile   (Form(kASigExtp_FitCheckRst,(TString("Multiplicity")+kFolder).Data()),"recreate");
        //
        g1D_Nres_Stat_Mult->Write();
        g2D_Nres_Stat_Mult->Write();
        gR1_Nres_Stat_Mult->Write();
        gR2_Nres_Stat_Mult->Write();
        gP1_Nres_Stat_Mult->Write();
        gP2_Nres_Stat_Mult->Write();
        g1D_Nres_Syst_Mult->Write();
        g2D_Nres_Syst_Mult->Write();
        gR1_Nres_Syst_Mult->Write();
        gR2_Nres_Syst_Mult->Write();
        gP1_Nres_Syst_Mult->Write();
        gP2_Nres_Syst_Mult->Write();
        for ( auto kSave : h1D_Nres_stat_MT )   kSave ->  Write();
        for ( auto kSave : h1D_Nres_syst_MT )   kSave ->  Write();
        for ( auto kSave : h2D_Nraw_stat_MT )   kSave ->  Write();
        for ( auto kSave : h2D_Nraw_syst_MT )   kSave ->  Write();
        for ( auto kSave : h2D_MnPt_stat_MT )   kSave ->  Write();
        for ( auto kSave : h2D_MnPt_syst_MT )   kSave ->  Write();
        for ( auto kSaveArray : h2D_Slcs_stat_MT ) for ( auto kSave : kSaveArray ) kSave ->  Write();
        for ( auto kSaveArray : h2D_Slcs_syst_MT ) for ( auto kSave : kSaveArray ) kSave ->  Write();
        for ( auto kSave : hXD_Nfqs_stat )   kSave ->  Write();
        for ( auto kSave : hXD_Nfqs_syst )   kSave ->  Write();
        //
        outFile_Rslt_ML->Close();
        outFile_Chck_ML->Close();
        insFile_Effc_ML->Close();
        insFile_Data_ML->Close();
    }
    // --- CORRELATION ANALYSIS
    /* if ( kDoCorrelation && false) {
        //  --- Load Files
        TFile*  insFile_Data_CR     =   new TFile   (Form(kASigExtr_FitCheckRst,(TString("Correlation")+kFolder).Data()));
        TFile*  insFile_Effc_CR     =   new TFile   (Form(kAnalysis_MCTruthHist,(TString("Correlation")+kFolder).Data()));
        //
        //  --- Load Histograms
        auto        fHEventCount    =   uLoadHistograms<0,TH1F> ( insFile_Data_CR, "fQC_Event_Enum_FLL" );
        auto        h1D_Nrec_2Db    =   uLoadHistograms<0,TH1F> ( insFile_Effc_CR, "h1D_Nrec_2Db" );
        auto        h1D_Ngen_2Db    =   uLoadHistograms<0,TH1F> ( insFile_Effc_CR, "h1D_Ngen_2Db" );
        auto        h2D_Nraw_CR     =   uLoadHistograms<1,TH2F> ( insFile_Data_CR, "anSS2D_CR_%i" );
        //
        //  --- Output File for Fit Check
        TFile*  outFile_Chck_CR     =   new TFile(Form(kASigExtp_FitCheckPlt,(TString("Correlation")+kFolder).Data()),"recreate");
        //
        //  --- Minimum Bias Normalisation
        Double_t    f1DCorrection   =   -1;
        Double_t    f2DCorrection   =   -1;
        if ( is_pp_anl ) {
            auto        kN_PU           =   (fHEventCount->GetBinContent(kEventCount::kPU_MB));
            auto        kN_Trg          =   (fHEventCount->GetBinContent(kEventCount::kTrigger));
            auto        kN_Vtx          =   (fHEventCount->GetBinContent(kEventCount::kVertex));
            auto        kN_MB           =   -kN_PU +(fHEventCount->GetBinContent(kEventCount::kVertex10));
            f1DCorrection               =   (1./kBR)        *(1./kN_MB) *(kN_Vtx/kN_Trg);
            f2DCorrection               =   (1./(kBR*kBR))  *(1./kN_MB) *(kN_Vtx/kN_Trg);
            if ( fabs( kEnergy - 5  ) < 0.1 ) { f1DCorrection   *=  kTriggerEff15n17pq; f2DCorrection   *=  kTriggerEff15n17pq; }
            if ( fabs( kEnergy - 7  ) < 0.1 ) { f1DCorrection   *=  kTriggerEff10bcdef; f2DCorrection   *=  kTriggerEff10bcdef; }
            if ( fabs( kEnergy - 13 ) < 0.1 ) { f1DCorrection   *=  kTriggerEff13TeV;   f2DCorrection   *=  kTriggerEff13TeV;   }
        } else if ( is_pb_anl ) {
            auto        kN_PU           =   (fHEventCount->GetBinContent(kEventCount::kPU_MB));
            auto        kN_Trg          =   (fHEventCount->GetBinContent(kEventCount::kTrigger));
            auto        kN_Vtx          =   (fHEventCount->GetBinContent(kEventCount::kVertex));
            auto        kN_MB           =   -kN_PU +(fHEventCount->GetBinContent(kEventCount::kVertex10));
            f1DCorrection               =   (1./kBR)        *(1./kN_MB) *(kN_Vtx/kN_Trg) *(1./0.5); // y normalisation 0.5
            f2DCorrection               =   (1./(kBR*kBR))  *(1./kN_MB) *(kN_Vtx/kN_Trg) *(1./0.5); // y normalisation 0.5
            if ( fabs( kEnergy - 5  ) < 0.1 ) {   f1DCorrection   *=  0.978;     f2DCorrection       *=  0.978; }
        }
        //  --- Build output structure
        Int_t   iCrPh = -1;
        std::vector<TString>    kPlotFolder1D;
        std::vector<TString>    kPlotFolder2D;
        std::vector<TH1F*>      h1D_Nres_stat_CR;
        std::vector<TH1F*>      h1D_Nres_syst_CR;
        std::vector<TH2F*>      h2D_Nraw_stat_CR;
        std::vector<TH2F*>      h2D_Nraw_syst_CR;
        std::vector<std::map<TString,std::tuple<Float_t,Float_t,Float_t>>>                  kMult1D_Results;
        std::vector<std::vector<std::map<TString,std::tuple<Float_t,Float_t,Float_t>>>>     kMult2D_Results;
        TString kPlotFolder1D_Full  =   TString( Form(kASigExtp_Plot_Direct,(TString("Correlation")+kFolder).Data()) ) +TString(Form("/1D/"));
        TString kPlotFolder2D_Full  =   TString( Form(kASigExtp_Plot_Direct,(TString("Correlation")+kFolder).Data()) ) +TString(Form("/2D/"));
        TString kPlotFolderFL_Full  =   TString( Form(kASigExtp_Plot_Direct,(TString("Correlation")+kFolder).Data()) ) +TString(Form("/Full/"));
        gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,(TString("Correlation")+kFolder).Data()))+TString(Form("/1D/")));
        gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,(TString("Correlation")+kFolder).Data()))+TString(Form("/2D/")));
        for ( auto kTarget : h2D_Nraw_CR )   {
            iCrPh++;
            //
            fSetAllCustomFunctions();
            //
            gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,(TString("Correlation")+kFolder).Data()))+TString(Form("/CRL_%i/1D/",iCrPh)));
            gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,(TString("Correlation")+kFolder).Data()))+TString(Form("/CRL_%i/2D/",iCrPh)));
            gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,(TString("Correlation")+kFolder).Data()))+TString(Form("/CRL_%i/Full/",iCrPh)));
            kPlotFolder1D.push_back( TString( Form(kASigExtp_Plot_Direct,(TString("Correlation")+kFolder).Data()) ) +TString(Form("/CRL_%i/1D/",iCrPh)) );
            kPlotFolder2D.push_back( TString( Form(kASigExtp_Plot_Direct,(TString("Correlation")+kFolder).Data()) ) +TString(Form("/CRL_%i/2D/",iCrPh)) );
            //
            h2D_Nraw_CR.at(iCrPh)        ->  Scale(1.,"width");
            h2D_Nraw_stat_CR.push_back( uEfficiencyCorrection2D ( h2D_Nraw_CR.at(iCrPh), h1D_Nrec_2Db, h1D_Ngen_2Db, f2DCorrection ) );
            h2D_Nraw_syst_CR.push_back( fSetSystErrors ( h2D_Nraw_stat_CR.at(iCrPh), "" ) );
            h2D_Nraw_stat_CR.at(iCrPh)   ->  SetName(Form("h2D_Nraw_stat_CR_%i",iCrPh));
            h2D_Nraw_syst_CR.at(iCrPh)   ->  SetName(Form("h2D_Nraw_syst_CR_%i",iCrPh));
            //
            kMult2D_Results.push_back( uMeasureFullYield2D( h2D_Nraw_stat_CR.at(iCrPh),h2D_Nraw_syst_CR.at(iCrPh),kStandardSystematicFitFunctions,kStatEvalCycles,kPlotFolder2D.at(iCrPh), "2D_%i" ) );
            //
        }
        //
        //  --- Final Result Histograms
        //
        hName   =   Form( "h1D_PhiCorr_Stat" );
        hTitle  =   Form( "h1D_PhiCorr_Stat" );
        TH1F*   h1D_PhiCorr_Stat    =   new TH1F ( hName, hTitle, nBinCrPh, fArrCrPh );
        SetAxis( h1D_PhiCorr_Stat, "CR" );
        //
        hName   =   Form( "h1D_PhiCorr_Syst" );
        hTitle  =   Form( "h1D_PhiCorr_Syst" );
        TH1F*   h1D_PhiCorr_Syst    =   new TH1F ( hName, hTitle, nBinCrPh, fArrCrPh );
        SetAxis( h1D_PhiCorr_Syst, "CR" );
        //
        auto iBin = -1;
        for ( auto kCurrent_2D_Result : kMult2D_Results ) {
            iBin++;
            auto    kYield_2D   =   get<0>( kCurrent_2D_Result.at(0)["YL_FLL"] );
            auto    kError_2D   =   get<1>( kCurrent_2D_Result.at(0)["YL_FLL"] );
            auto    kErSys_2D   =   get<2>( kCurrent_2D_Result.at(0)["YL_FLL"] );
            auto    kErRel_2D   =   kError_2D/kYield_2D;
            //
            cout << "iBin: " << iBin << endl;
            cout << "kYield_2D: " << kYield_2D << endl;
            cout << "kError_2D: " << kError_2D << endl;
            cout << "kErSys_2D: " << kErSys_2D << endl;
            cout << "kErRel_2D: " << kErRel_2D << endl;
            //
            h1D_PhiCorr_Stat->SetBinContent ( iBin+1, kYield_2D );
            h1D_PhiCorr_Stat->SetBinError   ( iBin+1, kError_2D );
            h1D_PhiCorr_Syst->SetBinContent ( iBin+1, kYield_2D );
            h1D_PhiCorr_Syst->SetBinError   ( iBin+1, kErSys_2D );
        }
        //
        fStopTimer("Correlation Analysis Signal Correction");
        //
        //  --- Output
        TFile *outFile_Rslt_CR      =   new TFile   (Form(kASigExtp_FitCheckRst,(TString("Correlation")+kFolder).Data()),"recreate");
        //
        h1D_PhiCorr_Stat->Scale(1.,"width");
        h1D_PhiCorr_Syst->Scale(1.,"width");
        //
        h1D_PhiCorr_Stat->Write();
        h1D_PhiCorr_Syst->Write();
        //
        uPlotCorrelation ( h1D_PhiCorr_Stat, h1D_PhiCorr_Syst, kPlotFolderFL_Full )->Write();
        //
        outFile_Rslt_CR->Close();
        outFile_Chck_CR->Close();
        insFile_Effc_CR->Close();
        insFile_Data_CR->Close();
    }*/
}
