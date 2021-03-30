// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../inc/AliAnalysisPhiPair.h"
#include "RooMsgService.h"

void SE_Production ( bool fSilent = true, TString fOption = "" )
{
    //---------------------//
    //  Setting up input   //
    //---------------------//
    
    //>-> OPTIONS
    
    // Silencing warnings for smoother running
    if ( fSilent )
    {
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
        RooMsgService::instance().setSilentMode(fSilent);
    }
    fChooseOption(fOption);
    
    // Retrieving PreProcessed data histograms
    TFile*  insFile_DT_Yield            =   new TFile   (fYldPreProc);
    
    // Recovering the histograms-------------------------------------------------------------------------------

    // >-> YIELD ANALYSIS //

    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F       *hREC_1D;
    TH1F      **hREC_1D_in_PT               = new TH1F     *[nBinPT1D];
    //
    //  Defining cumulative histogram over measurable pT
    //
    hName       =   "hREC_1D";
    hREC_1D     =   (TH1F*)(insFile_DT_Yield->Get(hName));
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
    TH1F      **hREC_1D_in_PT_2D_bin        = new TH1F     *[nBinPT2D];
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
        hName = Form("hREC_1D_in_PT_2D_bin_%i",iHisto);
        hREC_1D_in_PT_2D_bin[iHisto]    =   (TH1F*)(insFile_DT_Yield->Get(hName));
        hREC_2D_in_PT[iHisto]           =   new TH2F       *[nBinPT2D];
        
        for ( Int_t jHisto = 0; jHisto < nBinPT2D; jHisto++ )
        {
            hName = Form("hREC_2D_in_PT_%i_%i",iHisto,jHisto);
            hREC_2D_in_PT[iHisto][jHisto]    = (TH2F*)(insFile_DT_Yield->Get(hName));
        }
    }
    //
    //---------------------//
    //  Setting up output  //
    //---------------------//
    
    // Creating the histograms-------------------------------------------------------------------------------

    // Generating the binning array--------------------------------------------------------------------------
    fSetBinPT1D();
    fSetBinIM1D();
    fSetBinPT2D();
    fSetBinIM2D();
    fSetBinRap_();
    fSetBinMult();
    fSetBinNTup();
    
    // Creating the histograms-------------------------------------------------------------------------------

    // >-> YIELD ANALYSIS //

    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F       *hRAW_1D;
    TH1F       *hRAW_1D_in_2D_bin;
    //
    //  Defining yield histogram over measurable pT
    //
    hName       =   Form("hRAW_1D");
    hTitle      =   Form("hRAW_1D");
    hRAW_1D     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hRAW_1D,"PT 1D");
    //
    hName       =   Form("hRAW_1D_in_2D_bin");
    hTitle      =   Form("hRAW_1D_in_2D_bin");
    hRAW_1D_in_2D_bin     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hRAW_1D_in_2D_bin,"PT 1D");
    //
    // >->-->-> 2-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH2F       *hRAW_2D;
    //
    //  Defining yield histogram over measurable pT
    //
    hName       =   Form("hRAW_2D");
    hTitle      =   Form("hRAW_2D");
    hRAW_2D     =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    SetAxis(hRAW_2D,"PT 2D");
    //
    //-------------------------//
    //  Filling output objects //
    //-------------------------//
    //
    //  Making utility variables
    Int_t fTotalCount, fProgrCount;
    //
    fStartTimer("Systematics Yield Extraction Evaluation");
         
    // Total Fit number abnd progressive
    fTotalCount = nOptions*nBinPT1D+nOption2*(nBinPT2D*(nBinPT2D+1));
    fProgrCount = 0;
    //
    // Starting cycle
    //
    for ( Int_t iSys = 0; iSys < nOptions; iSys++ )
    {
        TFile      *outCheckSE_Syst =   new TFile (Form("result/Syst_SE/FitCheck_1D_%s.root",sOptions[iSys].c_str()),"recreate");
        cout << Form("[INFO] Starting yield analysis in 1D, mode: %s",sOptions[iSys].c_str()) << endl;
        for ( Int_t iFit = 0; iFit < nBinPT1D; iFit++ )
        {
            //Progressive Count
            fProgrCount++;
            //
            // Print Progress
            fPrintLoopTimer("Systematics Yield Extraction Evaluation",fProgrCount,fTotalCount,1);
            //
            // Fit
            auto fResults       =   FitModel(hREC_1D_in_PT[iFit],"1D",true,iFit,1,sOptions[iSys].c_str());
            
            if ( !fResults ) continue;
            
            // Building N_Raw histogram
            auto N_Raw      = static_cast<RooRealVar*>(fResults->floatParsFinal().at(Signal));
            hRAW_1D->SetBinContent      (iFit+1,N_Raw->getVal());
            hRAW_1D->SetBinError        (iFit+1,N_Raw->getError());
        }
        //
        outCheckSE_Syst->Close();
        
        TFile      *outFinalResults =   new TFile (Form("result/Syst_SE/Result_1D_%s.root",sOptions[iSys].c_str()),"recreate");
        
        hRAW_1D->Write();
        
        outFinalResults->Close();
    }
    
    for ( Int_t iSys = 0; iSys < nOption2; iSys++ )
    {
        // Storing Results for shape
        //
        RooFitResult  **fShapeStore =   new RooFitResult   *[nBinPT2D];
        //
        TFile      *outCheckSE_Syst =   new TFile (Form("result/Syst_SE/FitCheck_2D_%s.root",sOption2[iSys].c_str()),"recreate");
        cout << Form("[INFO] Starting yield analysis in 1D in 2D bins, mode: %s",sOption2[iSys].c_str()) << endl;
        for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ )
        {
            //Progressive Count
            fProgrCount++;
            
            // Print Progress
            fPrintLoopTimer("Systematics Yield Extraction Evaluation",fProgrCount,fTotalCount,1);
            
            // Fit
            fShapeStore[iFit]   =   FitModel(hREC_1D_in_PT_2D_bin[iFit],"2D_bin",true,iFit,2,sOption2[iSys].c_str());
            
            if ( !fShapeStore[iFit]  ) continue;
            
            // Building N_Raw histogram
            auto N_Raw      = static_cast<RooRealVar*>(fShapeStore[iFit] ->floatParsFinal().at(Signal));
            hRAW_1D_in_2D_bin->SetBinContent      (iFit+1,N_Raw->getVal());
            hRAW_1D_in_2D_bin->SetBinError        (iFit+1,N_Raw->getError());
        }
        //
        cout << Form("[INFO] Starting yield analysis in 2D, mode: %s",sOption2[iSys].c_str()) << endl;
        for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ ) {
            for ( Int_t jFit = 0; jFit < nBinPT2D; jFit++ ) {
                //Progressive Count
                fProgrCount++;
                    
                // Print Progress
                fPrintLoopTimer("Systematics Yield Extraction Evaluation",fProgrCount,fTotalCount,1);
                    
                if ( !fShapeStore[iFit] ) continue;
                if ( !fShapeStore[jFit] ) continue;
                
                // Fit
                auto fResults       =   FitModel(hREC_2D_in_PT[iFit][jFit],fShapeStore[iFit],fShapeStore[jFit],"2D",true,iFit,jFit,sOption2[iSys].c_str());
                
                if ( !fResults ) continue;
                
                // Building N_Raw histogram
                auto N_Raw      = static_cast<RooRealVar*>(fResults->floatParsFinal().at(SignlSignl));
                hRAW_2D->SetBinContent      (iFit+1,jFit+1,N_Raw->getVal());
                hRAW_2D->SetBinError        (iFit+1,jFit+1,N_Raw->getError());
            }
        }
        //
        delete[] fShapeStore;
        fStopTimer("Systematics Yield Extraction Evaluation");
        //
        outCheckSE_Syst->Close();
        
        TFile      *outFinalResults =   new TFile (Form("result/Syst_SE/Result_2D_%s.root",sOption2[iSys].c_str()),"recreate");
        hRAW_1D_in_2D_bin->Write();
        hRAW_2D->Write();
        
        outFinalResults->Close();
    }
    
    insFile_DT_Yield->Close();
}
