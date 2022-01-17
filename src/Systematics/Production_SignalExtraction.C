// File for 1-Dimensional Analysis:
// !TODO: 1. Move the results folder in the signal extraction folder 
#include "../../inc/AliAnalysisPhiPair.h"
#include "RooMsgService.h"

void
Production_SignalExtraction
 ( TString fOption = "", Bool_t fSilent = true )  {
    //
    //-----------------------------//
    //  Setting general analysis   //
    //-----------------------------//
    //
    //  Verbose option
    if ( fSilent )  {
        gErrorIgnoreLevel = kWarning;
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
        RooMsgService::instance().setSilentMode(fSilent);
    }
    //
    //  Option chosing
    if ( !fChooseOption(fOption) ) return;
    //
    //  Generating the binning array--------------------------------------------------------------------------
    fSetBinPT1D();
    fSetBinIM1D();
    fSetBinPT2D();
    fSetBinIM2D();
    fSetBinRap_();
    fSetBinMult();
    fSetBinNTup();
    //
    //  Utility variables-------------------------------------------------------------------------------------
    Int_t fTotalCount, fProgrCount;
    //
    TFile*  insFile_DT_Slp     =   new TFile   ( Form( kMassResolution_Anal,  "yield" ) );
    TH1F       *hSlop_Reference;
    TH1F       *hSlop_Referen2D;
    hName                       =   "hSigmaCnt_1D";
    hSlop_Reference             =   (TH1F*)(insFile_DT_Slp->Get(hName));
    hName                       =   "hSigmaCnt_1D_in_2D_bin";
    hSlop_Referen2D             =   (TH1F*)(insFile_DT_Slp->Get(hName));
    //
    if ( kDoYield ) {
        //
        //  Retrieving PreProcessed data histograms & reference resolution------------------------------------
        //
        TFile*  insFile_RawData     =   new TFile   ( Form( kAnalysis_InvMassHist,  "yield" ) );
        TFile*  insFile_RefSlop     =   new TFile   ( Form( kMassResolution_Anal,  "yield" ) );
        //
        //  RawData Histograms
        //  >   Definition
        TH1F**  hREC_1D_in_PT           =   new TH1F*   [nBinPT1D];
        TH1F**  hREC_1D_in_2D_bin_in_PT =   new TH1F*   [nBinPT2D];
        TH2F*** hREC_2D_in_PT           =   new TH2F**  [nBinPT2D];
        //  >   Recovering
        for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )  {
            hName                       =   Form("hREC_1D_in_PT_%i",iPT1D);
            hREC_1D_in_PT[iPT1D]        =   (TH1F*)(insFile_RawData->Get(hName));
        }
        for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
            hName                           =   Form("hREC_1D_in_PT_2D_bin_%i",iPT2D);
            hREC_1D_in_2D_bin_in_PT[iPT2D]  =   (TH1F*)(insFile_RawData->Get(hName));
            hREC_2D_in_PT[iPT2D]            =   new TH2F       *[nBinPT2D];
            
            for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )  {
                hName                           =   Form("hREC_2D_in_PT_%i_%i",iPT2D,jPT2D);
                hREC_2D_in_PT[iPT2D][jPT2D]     =  (TH2F*)(insFile_RawData->Get(hName));
            }
        }
        //
        //  Fitting histograms with systematic variations-----------------------------------------------------
        //  >   Full analysis
        fStartTimer("Full Systematics Yield Extraction Evaluation");
        //
        //  >   1D Analysis
        fStartTimer("Systematics Yield Extraction Evaluation 1D");
        //  >   >   Loop over Systematical variations
        fProgrCount =   0;
        fTotalCount =   kSyst_SEX_1D_Options.size();
        for ( auto kCurrent_Syst : kSyst_SEX_1D_Options ) {
            //  >   Progress print
            fPrintLoopTimer("Systematics Yield Extraction Evaluation 1D",fProgrCount,fTotalCount,1);
            //  >   Make fit check directory and file
            gROOT                       ->  ProcessLine (Form(".! mkdir -p %s/Systematics/ExtractionCheck/%s/1D/",(Form(kAnalysis_SigExtr_Dir,"yield")),kCurrent_Syst.Data()));
            TFile      *fCheckFit       =   new TFile   (Form("%s/Systematics/ExtractionCheck/%s/1D/CheckFitResults_%s.root",(Form(kAnalysis_SigExtr_Dir,"yield")),kCurrent_Syst.Data(),kCurrent_Syst.Data()),"recreate");
            //  >   Fit the Model
            auto    fFitResults_1DYield =   FitModel    (hREC_1D_in_PT,hSlop_Reference,Form("%s/Systematics/ExtractionCheck/%s/1D/",(Form(kAnalysis_SigExtr_Dir,"yield")),kCurrent_Syst.Data()),kCurrent_Syst.Data(),kCurrent_Syst.Data());
            //  >   Progressive Count
            fProgrCount++;
            //  >   Save output
            TFile      *fResultFit       =  new TFile   (Form("%s/Systematics/ExtractionCheck/%s/1D/FitResults_%s.root",(Form(kAnalysis_SigExtr_Dir,"yield")),kCurrent_Syst.Data(),kCurrent_Syst.Data()),"recreate");
            // !TODO: Re-arrange this
            for ( auto hSave : fFitResults_1DYield )    {
                if ( strncmp(hSave->GetName(),"hRAW_1D",7) == 0 )   hSave->Scale(1.,"width");
                if ( strncmp(hSave->GetName(),"h1D_anBB",7) == 0 )  hSave->Scale(1.,"width");
                hSave   ->  Write();
            }
            //  >   Wrap up, closing files
            fResultFit->Close();
            fCheckFit->Close();
        }
        fStopTimer("Systematics Yield Extraction Evaluation 1D");
        //
        //  >   2D Analysis
        fStartTimer("Systematics Yield Extraction Evaluation 2D");
        //  >   >   Loop over Systematical variations
        fProgrCount =   0;
        fTotalCount =   kSyst_SEX_2D_Options.size();
        for ( auto kCurrent_Syst : kSyst_SEX_2D_Options ) {
            //  >   Progress print
            fPrintLoopTimer("Systematics Yield Extraction Evaluation 2D",fProgrCount,fTotalCount,1);
            //  >   Make fit check directory and file
            gROOT                       ->  ProcessLine (Form(".! mkdir -p %s/Systematics/ExtractionCheck/%s/2D/",(Form(kAnalysis_SigExtr_Dir,"yield")),kCurrent_Syst.Data()));
            TFile      *fCheckFit       =   new TFile   (Form("%s/Systematics/ExtractionCheck/%s/2D/CheckFitResults_%s.root",(Form(kAnalysis_SigExtr_Dir,"yield")),kCurrent_Syst.Data(),kCurrent_Syst.Data()),"recreate");
            //  >   Fit the Model
            std::vector<TH1F*>  f1DCheck;
            auto    fFitResults_2DYield = FitModel(hREC_1D_in_2D_bin_in_PT,hSlop_Referen2D,hREC_2D_in_PT,f1DCheck,Form("%s/Systematics/ExtractionCheck/%s/2D/",(Form(kAnalysis_SigExtr_Dir,"yield")),kCurrent_Syst.Data()),kCurrent_Syst.Data(),kCurrent_Syst.Data());
            //  >   Progressive Count
            fProgrCount++;
            //  >   Save output
            TFile      *fResultFit       =   new TFile   (Form("%s/Systematics/ExtractionCheck/%s/2D/FitResults_%s.root",(Form(kAnalysis_SigExtr_Dir,"yield")),kCurrent_Syst.Data(),kCurrent_Syst.Data()),"recreate");
            // !TODO: Re-arrange this
            for ( auto hSave : f1DCheck )    {
                if ( strncmp(hSave->GetName(),"hRAW_1D_in_2D_bin",17) == 0 )    hSave->Scale(1.,"width");
                if ( strncmp(hSave->GetName(),"2Dbin_anBB",17) == 0 )           hSave->Scale(1.,"width");
                hSave   ->  Write();
            }
            for ( auto hSave : fFitResults_2DYield )    {
                if ( strncmp(hSave->GetName(),"hRAW_2D",7) == 0 )   hSave->Scale(1.,"width");
                if ( strncmp(hSave->GetName(),"anSB2D",6) == 0 )    hSave->Scale(1.,"width");
                if ( strncmp(hSave->GetName(),"anBS2D",6) == 0 )    hSave->Scale(1.,"width");
                if ( strncmp(hSave->GetName(),"anBB2D",6) == 0 )    hSave->Scale(1.,"width");
                hSave   ->  Write();
            }
            //  >   Wrap up, closing files
            fResultFit->Close();
            fCheckFit->Close();
        }
        fStopTimer("Systematics Yield Extraction Evaluation 2D");
        //
        //  Wrapping up, closing files------------------------------------------------------------------------
        insFile_RawData->Close();
        insFile_RefSlop->Close();
        fStopTimer("Full Systematics Yield Extraction Evaluation");
        //
    }
    if ( kDoMultiplicity ) {
        //
        // Retrieving PreProcessed data histograms & reference resolution------------------------------------
        //
        TFile*  insFile_RawData     =   new TFile   ( Form( kAnalysis_InvMassHist,  "multiplicity" ) );
        TFile*  insFile_RefSlop     =   new TFile   ( Form( kMassResolution_Anal,   "multiplicity" ) );
        //
        //  RawData Histograms
        //  >   Definition
        TH1F*** hREC_1D_in_PT_in_MT             =   new TH1F**  [nBinMult];
        TH1F*** hREC_1D_in_2D_bin_in_PT_in_MT   =   new TH1F**  [nBinMult];
        TH2F****hREC_2D_in_PT_in_MT             =   new TH2F*** [nBinMult];
        //  >   Recovering
        for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )  {
            hREC_1D_in_PT_in_MT[iMult]              =   new TH1F*   [nBinPT1D];
            hREC_1D_in_2D_bin_in_PT_in_MT[iMult]    =   new TH1F*   [nBinPT2D];
            hREC_2D_in_PT_in_MT[iMult]              =   new TH2F**  [nBinPT2D];
            for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )  {
                hName                                   =   Form("hREC_1D_in_MT_PT_%i_%i",iMult,iPT1D);
                hREC_1D_in_PT_in_MT[iMult][iPT1D]       =   (TH1F*)(insFile_RawData->Get(hName));
            }
            for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
                hName                                       =   Form("hREC_1D_in_PT_2D_bin_%i_in_MT_%i",iPT2D,iMult);
                hREC_1D_in_2D_bin_in_PT_in_MT[iMult][iPT2D] =   (TH1F*)(insFile_RawData->Get(hName));
                hREC_2D_in_PT_in_MT[iMult][iPT2D]           =   new TH2F       *[nBinPT2D];
                
                for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )  {
                    hName                                       =   Form("hREC_2D_in_PT_%i_%i_in_MT_%i",iPT2D,jPT2D,iMult);
                    hREC_2D_in_PT_in_MT[iMult][iPT2D][jPT2D]    =  (TH2F*)(insFile_RawData->Get(hName));
                }
            }
        }
        //
        //  Fitting histograms with systematic variations-----------------------------------------------------
        //  >   Full analysis
        fStartTimer("Full Systematics Multiplicity Yield Extraction Evaluation");
        //
        //  >   1D Analysis
        fStartTimer("Systematics Multiplicity Yield Extraction Evaluation 1D");
        //  >   >   Loop over Systematical variations
        fProgrCount =   0;
        fTotalCount =   (nBinMult+1)*kSyst_SEX_1D_Options.size();
        for ( auto kCurrent_Syst : kSyst_SEX_1D_Options ) {
            for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )  {
                //  >   Progress print
                fPrintLoopTimer("Systematics Multiplicity Yield Extraction Evaluation 1D",fProgrCount,fTotalCount,1);
                //  >   Make fit check directory and file
                gROOT                       ->  ProcessLine (Form(".! mkdir -p %s/Systematics/ExtractionCheck/%s/1D/MLT_%i/",(Form(kAnalysis_SigExtr_Dir,"multiplicity")),kCurrent_Syst.Data(),iMult));
                TFile      *fCheckFit       =   new TFile   (Form("%s/Systematics/ExtractionCheck/%s/1D/MLT_%i/CheckFitResults_%s.root",(Form(kAnalysis_SigExtr_Dir,"multiplicity")),kCurrent_Syst.Data(),iMult,kCurrent_Syst.Data()),"recreate");
                //  >   Fit the Model
                auto    fFitResults_1DYield =   FitModel    (hREC_1D_in_PT_in_MT[iMult],hSlop_Reference,Form("%s/Systematics/ExtractionCheck/%s/1D/MLT_%i/",(Form(kAnalysis_SigExtr_Dir,"multiplicity")),kCurrent_Syst.Data(),iMult),kCurrent_Syst.Data(),kCurrent_Syst.Data());
                //  >   Progressive Count
                fProgrCount++;
                //  >   Save output
                TFile      *fResultFit       =  new TFile   (Form("%s/Systematics/ExtractionCheck/%s/1D/MLT_%i/FitResults_%s.root",(Form(kAnalysis_SigExtr_Dir,"multiplicity")),kCurrent_Syst.Data(),iMult,kCurrent_Syst.Data()),"recreate");
                // !TODO: Re-arrange this
                for ( auto hSave : fFitResults_1DYield )    {
                    auto    kSaveName   =   hSave-> GetName();
                    hSave->SetName(Form("%s_MLT_%i",kSaveName,iMult));
                    if ( strncmp(kSaveName, "hRAW_1D",  7 ) == 0 )  hSave->Scale(1.,"width");
                    if ( strncmp(kSaveName, "h1D_anBB", 7 ) == 0 )  hSave->Scale(1.,"width");
                    hSave   ->  Write();
                }
                //  >   Wrap up, closing files
                fResultFit->Close();
                fCheckFit->Close();
            }
        }
        fStopTimer("Systematics Multiplicity Yield Extraction Evaluation 1D");
        //
        //  >   2D Analysis
        fStartTimer("Systematics Multiplicity Yield Extraction Evaluation 2D");
        //  >   >   Loop over Systematical variations
        fProgrCount =   0;
        fTotalCount =   (nBinMult+1)*kSyst_SEX_2D_Options.size();
        for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )  {
            for ( auto kCurrent_Syst : kSyst_SEX_2D_Options ) {
                //  >   Progress print
                fPrintLoopTimer("Systematics Multiplicity Yield Extraction Evaluation 2D",fProgrCount,fTotalCount,1);
                //  >   Make fit check directory and file
                gROOT                       ->  ProcessLine (Form(".! mkdir -p %s/Systematics/ExtractionCheck/%s/2D/MLT_%i/",(Form(kAnalysis_SigExtr_Dir,"multiplicity")),kCurrent_Syst.Data(),iMult));
                TFile      *fCheckFit       =   new TFile   (Form("%s/Systematics/ExtractionCheck/%s/2D/MLT_%i/CheckFitResults_%s.root",(Form(kAnalysis_SigExtr_Dir,"multiplicity")),kCurrent_Syst.Data(),iMult,kCurrent_Syst.Data()),"recreate");
                //  >   Fit the Model
                std::vector<TH1F*>  f1DCheck;
                auto    fFitResults_2DYield = FitModel(hREC_1D_in_2D_bin_in_PT_in_MT[iMult],hSlop_Referen2D,hREC_2D_in_PT_in_MT[iMult],f1DCheck,Form("%s/Systematics/ExtractionCheck/%s/2D/MLT_%i/",(Form(kAnalysis_SigExtr_Dir,"multiplicity")),kCurrent_Syst.Data(),iMult),kCurrent_Syst.Data(),kCurrent_Syst.Data());
                //  >   Progressive Count
                fProgrCount++;
                //  >   Save output
                TFile      *fResultFit       =   new TFile   (Form("%s/Systematics/ExtractionCheck/%s/2D/MLT_%i/FitResults_%s.root",(Form(kAnalysis_SigExtr_Dir,"multiplicity")),kCurrent_Syst.Data(),iMult,kCurrent_Syst.Data()),"recreate");
                // !TODO: Re-arrange this
                for ( auto hSave : f1DCheck )    {
                    auto    kSaveName   =   hSave-> GetName();
                    hSave->SetName(Form("%s_MLT_%i",kSaveName,iMult));
                    if ( strncmp(hSave->GetName(),"hRAW_1D_in_2D_bin",17) == 0 )    hSave->Scale(1.,"width");
                    if ( strncmp(hSave->GetName(),"2Dbin_anBB",17) == 0 )           hSave->Scale(1.,"width");
                    hSave   ->  Write();
                }
                for ( auto hSave : fFitResults_2DYield )    {
                    auto    kSaveName   =   hSave-> GetName();
                    hSave->SetName(Form("%s_MLT_%i",kSaveName,iMult));
                    if ( strncmp(hSave->GetName(),"hRAW_2D",7) == 0 )   hSave->Scale(1.,"width");
                    if ( strncmp(hSave->GetName(),"anSB2D",6) == 0 )    hSave->Scale(1.,"width");
                    if ( strncmp(hSave->GetName(),"anBS2D",6) == 0 )    hSave->Scale(1.,"width");
                    if ( strncmp(hSave->GetName(),"anBB2D",6) == 0 )    hSave->Scale(1.,"width");
                    hSave   ->  Write();
                }
                //  >   Wrap up, closing files
                fResultFit->Close();
                fCheckFit->Close();
            }
        }
        fStopTimer("Systematics Multiplicity Yield Extraction Evaluation 2D");
        //
        //  Wrapping up, closing files------------------------------------------------------------------------
        insFile_RawData->Close();
        insFile_RefSlop->Close();
        fStopTimer("Full Systematics Multiplicity Yield Extraction Evaluation");
        //
    }
    /*
    
        // Recovering the histograms-------------------------------------------------------------------------------

        // >-> YIELD ANALYSIS //

        // >->-->-> 1-Dimension analysis //
        //
        //  Declaring all histograms
        //
        TH1F       *hREC_1D;
        //
        //  Defining cumulative histogram over measurable pT
        //
        hName       =   "hREC_1D";
        hREC_1D     =   (TH1F*)(insFile_DT_Yield->Get(hName));
        hName                       =   "hSigmaCnt_1D";
        hSlop_Reference             =   (TH1F*)(insFile_DT_Slp->Get(hName));
        hName                       =   "hSigmaCnt_1D_in_2D_bin";
        hSlop_Referen2D             =   (TH1F*)(insFile_DT_Slp->Get(hName));
        //
        //  Defining pT-Differential histograms over measurable pT
        //
        for ( Int_t iHisto = 0; iHisto < nBinPT1D; iHisto++ )
        {
            hName                   =   Form("hREC_1D_in_PT_%i",iHisto);
            hREC_1D_in_PT[iHisto]   =   (TH1F*)(insFile_DT_Yield->Get(hName));
        }

        // >->-->-> 2-Dimension analysis //
        //
        //  Declaring all histograms
        //
        TH2F       *hREC_2D;
        TH1F      **hREC_1D_in_2D_bin_in_PT        = new TH1F     *[nBinPT2D];
        TH2F     ***hREC_2D_in_PT               = new TH2F    **[nBinPT2D];
        //
        //  Defining cumulative histogram over measurable pT
        //
        hName       =   "hREC_2D";
        hREC_2D     =   (TH2F*)(insFile_DT_Yield->Get(hName));
        //
        //  Defining pT-Differential histograms over measurable pT
        //
        for ( Int_t iHisto = 0; iHisto < nBinPT2D; iHisto++ )
        {
            hName = Form("hREC_1D_in_2D_bin_in_PT_%i",iHisto);
            hREC_1D_in_2D_bin_in_PT[iHisto]    =   (TH1F*)(insFile_DT_Yield->Get(hName));
            hREC_2D_in_PT[iHisto]           =   new TH2F       *[nBinPT2D];
            
            for ( Int_t jHisto = 0; jHisto < nBinPT2D; jHisto++ )
            {
                hName = Form("hREC_2D_in_PT_%i_%i",iHisto,jHisto);
                hREC_2D_in_PT[iHisto][jHisto]    = (TH2F*)(insFile_DT_Yield->Get(hName));
            }
        }
        //
        //-------------------------//
        //  Systematics spectra    //
        //-------------------------//
        //
        //  Making utility variables
        Int_t fTotalCount, fProgrCount;
        //
        // Total Fit number abnd progressive
        fTotalCount = nOptions;
        fProgrCount = 0.;
        //
        // Starting cycle
        //
        fStartTimer("Systematics Yield Extraction Evaluation 1D");
        //
        for ( Int_t iSys = 0; iSys < nOptions; iSys++ ) {
            //
            gROOT                       ->  ProcessLine (Form(".! mkdir -p %s/Systematics/ExtractionCheck/%s/1D/",(Form(kAnalysis_SigExtr_Dir,"yield")),sOptions.at(iSys).Data()));
            TFile      *fCheckFit       =   new TFile   (Form("%s/Systematics/ExtractionCheck/%s/1D/CheckFitResults_%s.root",(Form(kAnalysis_SigExtr_Dir,"yield")),sOptions.at(iSys).Data(),sOptions.at(iSys).Data()),"recreate");
            //
            //>>    Fit the Model
            auto    fFitResults_1DYield = FitModel        (hREC_1D_in_PT,hSlop_Reference,Form("%s/Systematics/ExtractionCheck/%s/1D/",(Form(kAnalysis_SigExtr_Dir,"yield")),sOptions.at(iSys).Data()),sOptions.at(iSys).Data(),sOptions.at(iSys).Data());
            //
            //>>    Progressive Count
            fProgrCount++;
            //
            TFile      *fResultFit       =   new TFile   (Form("%s/Systematics/ExtractionCheck/%s/1D/FitResults_%s.root",(Form(kAnalysis_SigExtr_Dir,"yield")),sOptions.at(iSys).Data(),sOptions.at(iSys).Data()),"recreate");
            for ( auto hSave : fFitResults_1DYield )    {
                if ( strncmp(hSave->GetName(),"hRAW_1D",7) == 0 )   hSave->Scale(1.,"width");
                if ( strncmp(hSave->GetName(),"h1D_anBB",7) == 0 )  hSave->Scale(1.,"width");
                hSave   ->  Write();
            }
            //
            fResultFit->Close();
            fCheckFit->Close();
            //
            //>>    Print Progress
            fPrintLoopTimer("Systematics Yield Extraction Evaluation 1D",fProgrCount,fTotalCount,1);
        }
        //
        fStopTimer("Systematics Yield Extraction Evaluation 1D");
        //
        // Total Fit number abnd progressive
        fTotalCount = nOption2;
        fProgrCount = 0.;
        //
        fStartTimer("Systematics Yield Extraction Evaluation 2D");
        //
        for ( Int_t iSys = 0; iSys < nOption2; iSys++ ) {
            //
            gROOT                       ->  ProcessLine (Form(".! mkdir -p %s/Systematics/ExtractionCheck/%s/2D/",(Form(kAnalysis_SigExtr_Dir,"yield")),kCurrent_Syst.Data()));
            TFile      *fCheckFit       =   new TFile   (Form("%s/Systematics/ExtractionCheck/%s/2D/CheckFitResults_%s.root",(Form(kAnalysis_SigExtr_Dir,"yield")),sOption2.at(iSys).Data(),sOption2.at(iSys).Data()),"recreate");
            //
            //>>    Fit
            std::vector<TH1F*>  f1DCheck;
            auto    fFitResults_2DYield = FitModel(hREC_1D_in_2D_bin_in_PT,hSlop_Referen2D,hREC_2D_in_PT,f1DCheck,Form("%s/Systematics/ExtractionCheck/%s/2D/",(Form(kAnalysis_SigExtr_Dir,"yield")),sOption2.at(iSys).Data()),sOption2.at(iSys).Data(),sOption2.at(iSys).Data());
            //
            //Progressive Count
            fProgrCount++;
            //
            //
            TFile      *fResultFit       =   new TFile   (Form("%s/Systematics/ExtractionCheck/%s/2D/FitResults_%s.root",(Form(kAnalysis_SigExtr_Dir,"yield")),sOption2.at(iSys).Data(),sOption2.at(iSys).Data()),"recreate");
            for ( auto hSave : f1DCheck )    {
                if ( strncmp(hSave->GetName(),"hRAW_1D_in_2D_bin",17) == 0 )    hSave->Scale(1.,"width");
                if ( strncmp(hSave->GetName(),"2Dbin_anBB",17) == 0 )           hSave->Scale(1.,"width");
                hSave   ->  Write();
            }
            for ( auto hSave : fFitResults_2DYield )    {
                if ( strncmp(hSave->GetName(),"hRAW_2D",7) == 0 )   hSave->Scale(1.,"width");
                if ( strncmp(hSave->GetName(),"anSB2D",6) == 0 )    hSave->Scale(1.,"width");
                if ( strncmp(hSave->GetName(),"anBS2D",6) == 0 )    hSave->Scale(1.,"width");
                if ( strncmp(hSave->GetName(),"anBB2D",6) == 0 )    hSave->Scale(1.,"width");
                hSave   ->  Write();
            }
            //
            fResultFit->Close();
            fCheckFit->Close();
            //
            //>>    Print Progress
            fPrintLoopTimer("Systematics Yield Extraction Evaluation 2D",fProgrCount,fTotalCount,1);
        }
        fStopTimer("Systematics Yield Extraction Evaluation 2D");
        insFile_DT_Yield->Close();
    }
     */
}
