// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../inc/AliAnalysisPhiPair.h"
#include "GeneralAnalysis.cxx"
#include "RooMsgService.h"

void
Analysis_SignalExtraction
( bool fSilent = false)    {
    
}



    /*
     //
     //  YIELD Efficiency
     TFile      *insFile_EF_Yield    =   new TFile   ( Form(kAnalysis_MCTruthHist,"yield") );
     //
     TH1F       *hStandard;
     hName                           =   Form("hRAW_1D");
     hStandard                       =   (TH1F*)((new TFile(Form(kASigExtr_FitCheckRst,"Yield")))->Get(hName));
     //
     std::vector<TH1F*>  fVariations;
     hName           =   Form("hRAW_1D");
     for ( Int_t iTer = 0; iTer < nOptions; iTer++ )    {
         auto    fInput  =   new TFile   (Form("%s/ExtractionCheck/%s/1D/FitResults_%s.root",Form(kAnalysis_SgExSys_Dir,"yield"),sOptions.at(iTer).Data(),sOptions.at(iTer).Data()));
         auto    fTarget =   new TH1F( *(TH1F*)(fInput->Get(hName)) );
         fTarget->SetName(sOptions.at(iTer).Data());
         fVariations.push_back(fTarget);
         delete fInput;
         delete fTarget;
     }
     //
     GeneralAnalysis(hStandard,fVariations,TString(Form(kAnalysis_SgExSys_Dir,"yield")));
     //
    
    TFile **insFileH1D  =   new TFile*  [nOptions+1];
    TFile **insFileH2D  =   new TFile*  [nOption2];
    
    //Recovering histograms
    TH1F**  h1D_Syst    =   new TH1F*   [nOptions+1];
    TH2F**  h2D_Syst    =   new TH2F*   [nOption2+1];
    
    insFileH1D[0]   =   new TFile   (Form(kASigExtr_FitCheckRst,"yield"));
    hName           =   Form("hRAW_1D");
    h1D_Syst[0]     =   (TH1F*)(insFileH1D[0]->Get(hName));
    hName           =   Form("hRAW_2D");
    h2D_Syst[0]     =   (TH2F*)(insFileH1D[0]->Get(hName));
    for ( Int_t iTer = 1; iTer <= nOptions; iTer++ )
    {
        insFileH1D[iTer]=   new TFile   (Form("result/yield/ExtractionSystematics/ExtractionCheck/%s/1D/FitResults_%s.root",sOptions.at(iTer-1).Data(),sOptions.at(iTer-1).Data()));
        hName           =   Form("hRAW_1D");
        h1D_Syst[iTer]  =   (TH1F*)(insFileH1D[iTer]->Get(hName));
        //h1D_Syst[iTer]      ->Scale(1.,"width");
    }
    for ( Int_t iTer = 0; iTer < nOption2; iTer++ )
    {
        insFileH2D[iTer]=   new TFile   (Form("result/yield/ExtractionSystematics/ExtractionCheck/%s/2D/FitResults_%s.root",sOption2.at(iTer).Data(),sOption2.at(iTer).Data()));
        hName           =   Form("hRAW_2D");
        h2D_Syst[iTer+1]=   (TH2F*)(insFileH2D[iTer]->Get(hName));
        //h2D_Syst[iTer+1]    ->Scale(1.,"width");
    }
     */
/*
void Analysis_SignalExtraction ( bool fSilent = false)
{
    
    TFile*  insFile_EF_Yield            =   new TFile   (Form(kAnalysis_MCTruthHist,"yield"));
    
    TFile **insFileH1D  =   new TFile*  [nOptions+1];
    TFile **insFileH2D  =   new TFile*  [nOption2];
    
    //Recovering histograms
    TH1F**  h1D_Syst    =   new TH1F*   [nOptions+1];
    TH2F**  h2D_Syst    =   new TH2F*   [nOption2+1];
    
    insFileH1D[0]   =   new TFile   (Form(kASigExtr_FitCheckRst,"yield"));
    hName           =   Form("hRAW_1D");
    h1D_Syst[0]     =   (TH1F*)(insFileH1D[0]->Get(hName));
    hName           =   Form("hRAW_2D");
    h2D_Syst[0]     =   (TH2F*)(insFileH1D[0]->Get(hName));
    for ( Int_t iTer = 1; iTer <= nOptions; iTer++ )
    {
        insFileH1D[iTer]=   new TFile   (Form("result/yield/ExtractionSystematics/ExtractionCheck/%s/1D/FitResults_%s.root",sOptions.at(iTer-1).Data(),sOptions.at(iTer-1).Data()));
        hName           =   Form("hRAW_1D");
        h1D_Syst[iTer]  =   (TH1F*)(insFileH1D[iTer]->Get(hName));
        //h1D_Syst[iTer]      ->Scale(1.,"width");
    }
    for ( Int_t iTer = 0; iTer < nOption2; iTer++ )
    {
        insFileH2D[iTer]=   new TFile   (Form("result/yield/ExtractionSystematics/ExtractionCheck/%s/2D/FitResults_%s.root",sOption2.at(iTer).Data(),sOption2.at(iTer).Data()));
        hName           =   Form("hRAW_2D");
        h2D_Syst[iTer+1]=   (TH2F*)(insFileH2D[iTer]->Get(hName));
        //h2D_Syst[iTer+1]    ->Scale(1.,"width");
    }
    //
    TH1F        *hEFF_1D;
    TH2F        *hEFF_2D;
    //
    hName       =   "hEFF_1D";
    hEFF_1D     =   (TH1F*)(insFile_EF_Yield->Get(hName));
    //
    hName       =   "hEFF_2D_fr_1D";
    hEFF_2D     =   (TH2F*)(insFile_EF_Yield->Get(hName));
    
    
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
        for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )  {
            hEFF_2D->SetBinError( iPT2D, jPT2D, 0);
        }
    }
    for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )  {
        hEFF_1D->SetBinError( iPT1D, 0);
    }
    
    //---------------------//
    //  Setting up output  //
    //---------------------//
    
    //--- Generating the binning array ---//
    fSetBinPT1D();
    fSetBinIM1D();
    fSetBinPT2D();
    fSetBinIM2D();
    fSetBinSyst();
    
    // Creating the histograms-------------------------------------------------------------------------------
    //
    hName               =   Form("h1D_Syst_Bin_StSy");
    hTitle              =   Form("Fractional variation of raw yield for bin of PT [%.2f#;%.2f]",fArrPT1D[0],fArrPT1D[nBinPT1D]);
    TH1F   *h1D_Syst_Bin_StSy   =   new TH1F (hName,hTitle,1000,-.5,.5);
    //
    hName               =   Form("h2D_Syst_Bin_StSy");
    hTitle              =   Form("Fractional variation of raw yield for bin of PT [%.2f#;%.2f]",fArrPT1D[0],fArrPT1D[nBinPT1D]);
    TH1F   *h2D_Syst_Bin_StSy   =   new TH1F (hName,hTitle,400,-2.,2.);
    //
    hName               =   Form("h1D_Stat_Bin");
    hTitle              =   Form("h1D_Stat_Bin");
    TH1F   *h1D_Stat_Bin    =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    //
    TH1F  **h1D_Syst_Bin    =   new TH1F   *[nBinPT1D+1];
    hName               =   Form("h1D_Syst_Bin_PT_%.2f_%.2f",fArrPT1D[0],fArrPT1D[nBinPT1D]);
    hTitle              =   Form("Fractional variation of raw yield for bin of PT [%.2f#;%.2f]",fArrPT1D[0],fArrPT1D[nBinPT1D]);
    h1D_Syst_Bin[0]  =   new TH1F (hName,hTitle,1000,-.5,.5);
    h1D_Syst_Bin[0]  ->SetTitle(Form("Fractional variation of raw yield for bin %i",-1));
    h1D_Syst_Bin[0]  ->GetXaxis()->SetTitle("Fractional variation");
    for ( Int_t iAll = 1; iAll <= nBinPT1D; iAll++ )
    {
        hName               =   Form("h1D_Syst_Bin_PT_%.2f_%.2f",fArrPT1D[iAll-1],fArrPT1D[iAll]);
        hTitle              =   Form("Fractional variation of raw yield for bin of PT [%.2f#;%.2f]",fArrPT1D[iAll-1],fArrPT1D[iAll]);
        h1D_Syst_Bin[iAll]  =   new TH1F (hName,hTitle,1000,-.5,.5);
        h1D_Syst_Bin[iAll]  ->SetTitle(Form("Fractional variation of raw yield for bin %i",iAll));
        h1D_Syst_Bin[iAll]  ->GetXaxis()->SetTitle("Fractional variation");
    }
    //
    TH2F   *hCheckFull1D        =   new TH2F("hCheckFull1D","hCheckFull1D",nBinPT1D,fArrPT1D,nBinSyst,fArrSyst);
    //
    TH1F   *h2D_Syst_Bin_All;
    hName                       =   Form("h2D_Syst_Bin_PT_%.2f_%.2f_%.2f_%.2f",fArrPT2D[0],fArrPT2D[nBinPT2D],fArrPT2D[0],fArrPT2D[nBinPT2D]);
    hTitle                      =   Form("Fractional variation of raw yield for bin of PT [%.2f#;%.2f] [%.2f#;%.2f]",fArrPT2D[0],fArrPT2D[nBinPT2D],fArrPT2D[0],fArrPT2D[nBinPT2D]);
    h2D_Syst_Bin_All            =   new TH1F    (hName,hTitle,400,-2.,2.);
    h2D_Syst_Bin_All            ->  SetTitle(Form("Fractional variation of raw yield for bin %i",-1));
    h2D_Syst_Bin_All            ->  GetXaxis()  ->  SetTitle("Fractional variation");
    //
    TH2F   *h2D_Stat_Bin;
    hName                       =   Form("h2D_Stat_Bin");
    hTitle                      =   Form("h2D_Stat_Bin");
    h2D_Stat_Bin                =   new TH2F    (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    //
    TH1F ***h2D_Syst_Bin        =   new TH1F  **[nBinPT2D];
    TH2F  **hCheckFull2D        =   new TH2F   *[nBinPT2D];
    for ( Int_t iAll = 0; iAll < nBinPT2D; iAll++ ) {
        h2D_Syst_Bin[iAll]      =   new TH1F   *[nBinPT2D];
        hName                   =   Form("hCheckFull2D_%i",iAll);
        hCheckFull2D[iAll]      =   new TH2F(hName,hName,nBinPT2D,fArrPT2D,nBinSyst,fArrSyst);
        for ( Int_t jAll = 0; jAll < nBinPT2D; jAll++ ) {
            hName                       =   Form("h2D_Syst_Bin_PT_%.2f_%.2f_%.2f_%.2f",fArrPT2D[iAll],fArrPT2D[iAll+1],fArrPT2D[jAll],fArrPT2D[jAll+1]);
            hTitle                      =   Form("Fractional variation of raw yield for bin of PT [%.2f#;%.2f] [%.2f#;%.2f]",fArrPT2D[iAll],fArrPT2D[iAll+1],fArrPT2D[jAll],fArrPT2D[jAll+1]);
            h2D_Syst_Bin[iAll][jAll]    =   new TH1F (hName,hTitle,1000,-.5,.5);
            h2D_Syst_Bin[iAll][jAll]    ->  SetTitle(Form("Fractional variation of raw yield for bin %i",iAll+1));
            h2D_Syst_Bin[iAll][jAll]    ->  GetXaxis()  ->  SetTitle("Fractional variation");
        }
    }
    //------------//
    //  ANALYSIS  //
    //------------//
    
    // Output File for Fit Check
    TFile*  outFileFit  =   new TFile("./result/yield/ExtractionSystematics/ExtractionSystematics_CheckRatioAndBins.root","recreate");
    
    //------ 1D Histograms ------//
    //
    TGraphAsymmErrors     **g1D_Stat            =   new TGraphAsymmErrors  *[nOptions+1];
    for ( Int_t iTer = 0; iTer <= nOptions; iTer++ )    {
        auto fCheck = new TH1F (*h1D_Syst[iTer]);
        fCheck->Divide(h1D_Syst[iTer],h1D_Syst[0]);
        fCheck->Write();
        g1D_Stat[iTer]     =   fTH1_to_TGAsymmErrors(h1D_Syst[iTer]);
    }
    //
    TGraphAsymmErrors     **g1D_Stat_VarErr     =   new TGraphAsymmErrors  *[nOptions+1];
    for ( Int_t iTer = 1; iTer <= nOptions; iTer++ )    {
        g1D_Stat_VarErr[iTer]     =   new TGraphAsymmErrors();
    }
    //
    TGraphAsymmErrors  *g1D_Stat_Err        =   new TGraphAsymmErrors   ();
    TGraphAsymmErrors  *g1D_Syst_Err        =   new TGraphAsymmErrors   ();
    for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )  {
        //
        auto    fStandard       =   g1D_Stat[0]         ->GetPointY     (iPT1D);
        auto    fStdError       =   g1D_Stat[0]         ->GetErrorYhigh (iPT1D);
        h1D_Stat_Bin            ->  SetBinContent(iPT1D+1,fStdError/fStandard);
        //
        auto    fXPoint         =   fArrPT1D[iPT1D] + .5*(fArrPT1D[iPT1D+1] - fArrPT1D[iPT1D]);
        auto    fXError         =   .5*(fArrPT1D[iPT1D+1] - fArrPT1D[iPT1D]);
        auto    fYPoint         =   0.;
        auto    fYErrorHig      =   g1D_Stat[0]         ->GetErrorYhigh (iPT1D) / (fStandard);
        auto    fYErrorLow      =   g1D_Stat[0]         ->GetErrorYlow  (iPT1D) / (fStandard);
        g1D_Stat_Err            ->  SetPoint        (iPT1D, fXPoint,    fYPoint);
        g1D_Stat_Err            ->  SetPointError   (iPT1D, fXError,    fXError,    fYErrorHig, fYErrorLow);
        g1D_Syst_Err            ->  SetPoint        (iPT1D, fXPoint,    fYPoint);
        if ( iPT1D == 0  || iPT1D == 19 ) g1D_Syst_Err            ->  SetPointError   (iPT1D, fXError,    fXError,    0.040, 0.040);
        if ( iPT1D >= 1  && iPT1D <= 8  ) g1D_Syst_Err            ->  SetPointError   (iPT1D, fXError,    fXError,    0.010, 0.010);
        if ( iPT1D >= 9  && iPT1D <= 14 ) g1D_Syst_Err            ->  SetPointError   (iPT1D, fXError,    fXError,    0.015, 0.015);
        if ( iPT1D >= 15 && iPT1D <= 18 ) g1D_Syst_Err            ->  SetPointError   (iPT1D, fXError,    fXError,    0.035, 0.035);
        //
        for ( Int_t iTer = 1; iTer <= nOptions; iTer++ )    {
            auto    fVariatin       =   g1D_Stat[iTer]          ->GetPointY     (iPT1D);
            //
            auto    fYErrVrHig      =   g1D_Stat[iTer]         ->GetErrorYhigh (iPT1D) / (fStandard);
            auto    fYErrVrLow      =   g1D_Stat[iTer]         ->GetErrorYlow  (iPT1D) / (fStandard);
            //
            auto    fYFrac          =   fVariatin / fStandard - 1;
            auto    fYFracEHig      =   sqrt( fabs(fYErrorHig*fYErrorHig - fYErrVrHig*fYErrVrHig) );
            auto    fYFracELow      =   sqrt( fabs(fYErrorLow*fYErrorLow - fYErrVrLow*fYErrVrLow) );
            g1D_Stat_VarErr[iTer]   ->  SetPoint        (iPT1D, fXPoint,    fYFrac);
            g1D_Stat_VarErr[iTer]   ->  SetPointError   (iPT1D, fXError,    fXError,    fYFracEHig, fYFracELow);
            //
            if ( fVariatin - fStandard == 0 ) continue;
            if ( fBarlowCheck(fStandard,max(fStandard*fYErrorHig,fStandard*fYErrorLow),fVariatin,max(fStandard*fYErrVrHig,fStandard*fYErrVrLow)) )   continue;
            h1D_Syst_Bin[iPT1D+1]   ->  Fill(fYFrac);
            h1D_Syst_Bin[0]         ->  Fill(fYFrac);
            hCheckFull1D            ->  Fill(fXPoint,fYFrac);
        }
    }
    //
    h1D_Syst_Bin_StSy->Write();
    for ( Int_t iPT1D = 0; iPT1D <= nBinPT1D; iPT1D++ )  {
        h1D_Syst_Bin[iPT1D] ->  Write();
    }
    //
    //------ 2D Histograms ------//
    //
    TGraphAsymmErrors    ***g2D_Stat            =   new TGraphAsymmErrors  **[nOption2+1];
    for ( Int_t iTer = 0; iTer <= nOption2; iTer++ )    {
        auto fCheck = new TH2F (*h2D_Syst[iTer]);
        fCheck->Divide(h2D_Syst[iTer],h2D_Syst[0]);
        fCheck->Write();
        g2D_Stat[iTer]     =   fTH2_to_TGAsymmErrors(h2D_Syst[iTer]);
    }
    //
    TGraphAsymmErrors    ***g2D_Stat_VarErr        =   new TGraphAsymmErrors **[nBinPT2D];
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )    {
        g2D_Stat_VarErr[iPT2D]      =   new TGraphAsymmErrors   *[nOption2+1];
        for ( Int_t iTer = 1; iTer <= nOption2; iTer++ )    {
            g2D_Stat_VarErr[iPT2D][iTer]   =   new TGraphAsymmErrors();
        }
    }
    //
    TH2F*   h2D_Syst_Ful3   =   new TH2F    ("h2D_Syst_Ful3","h2D_Syst_Full_Averaged",      nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    TGraphAsymmErrors     **g2D_Stat_Err        =   new TGraphAsymmErrors  *[nBinPT2D];
    TGraphAsymmErrors     **g2D_Syst_Err        =   new TGraphAsymmErrors  *[nBinPT2D];
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
        g2D_Stat_Err[iPT2D]     =   new TGraphAsymmErrors();
        g2D_Syst_Err[iPT2D]     =   new TGraphAsymmErrors();
        for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )  {
            //
            auto    fStandard       =   g2D_Stat[0][iPT2D]         ->GetPointY     (jPT2D);
            auto    fStdError       =   g2D_Stat[0][iPT2D]         ->GetErrorYhigh (jPT2D);
            h2D_Stat_Bin            ->  SetBinContent(iPT2D+1,jPT2D+1,fStdError/fStandard);
            //
            auto    fXPoint         =   fArrPT2D[jPT2D] + .5*(fArrPT2D[jPT2D+1] - fArrPT2D[jPT2D]);
            auto    fXError         =   .5*(fArrPT2D[jPT2D+1] - fArrPT2D[jPT2D]);
            auto    fYPoint         =   0.;
            auto    fYErrorHig      =   g2D_Stat[0][iPT2D]          ->GetErrorYhigh (jPT2D) / (fStandard);
            auto    fYErrorLow      =   g2D_Stat[0][iPT2D]          ->GetErrorYlow  (jPT2D) / (fStandard);
            g2D_Stat_Err[iPT2D]     ->  SetPoint        (jPT2D, fXPoint,    fYPoint);
            g2D_Stat_Err[iPT2D]     ->  SetPointError   (jPT2D, fXError,    fXError,    fYErrorHig, fYErrorLow);
            g2D_Syst_Err[iPT2D]     ->  SetPoint        (jPT2D, fXPoint,    fYPoint);
            g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.050, 0.050);
            if ( iPT2D == 0  || jPT2D == 0  )       g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.075, 0.075);
            else if ( iPT2D == 9  || jPT2D == 9  )  g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.150, 0.150);
            else if ( iPT2D == 1  || jPT2D == 1  )  g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.075, 0.075);
            else if ( iPT2D == 4  || jPT2D == 4  )  g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.075, 0.075);
            else if ( iPT2D == 6  || jPT2D == 6  )  g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.075, 0.075);
            else if ( iPT2D == 8  || jPT2D == 8  )  g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.075, 0.075);
            //
            // Specific Points
            //
            if ( iPT2D == 0  && jPT2D == 9  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.200, 0.200);
            if ( iPT2D == 9  && jPT2D == 0  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.200, 0.200);
            //
            if ( iPT2D == 0  && jPT2D == 6  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.150, 0.150);
            if ( iPT2D == 6  && jPT2D == 0  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.150, 0.150);
            //
            if ( iPT2D == 2  && jPT2D == 9  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.220, 0.220);
            if ( iPT2D == 9  && jPT2D == 2  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.220, 0.220);
            //
            if ( iPT2D == 3  && jPT2D == 9  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.070, 0.070);
            if ( iPT2D == 9  && jPT2D == 3  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.070, 0.070);
            //
            if ( iPT2D == 5  && jPT2D == 9  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.050, 0.050);
            if ( iPT2D == 9  && jPT2D == 5  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.050, 0.050);
            //
            if ( iPT2D == 1  && jPT2D == 2  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.040, 0.040);
            if ( iPT2D == 2  && jPT2D == 1  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.040, 0.040);
            //
            if ( iPT2D == 1  && jPT2D == 3  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.040, 0.040);
            if ( iPT2D == 3  && jPT2D == 1  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.040, 0.040);
            //
            if ( iPT2D == 1  && jPT2D == 4  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.040, 0.040);
            if ( iPT2D == 4  && jPT2D == 1  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.040, 0.040);
            //
            if ( iPT2D == 7  && jPT2D == 9  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.090, 0.090);
            if ( iPT2D == 9  && jPT2D == 7  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.090, 0.090);
            //
            if ( iPT2D == 2  && jPT2D == 5  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.030, 0.030);
            if ( iPT2D == 5  && jPT2D == 2  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.030, 0.030);
            //
            if ( iPT2D == 2  && jPT2D == 0  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.120, 0.120);
            if ( iPT2D == 0  && jPT2D == 2  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.120, 0.120);
            //
            // Diagonal
            if ( iPT2D == 0  && jPT2D == 0  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.400, 0.400);
            if ( iPT2D == 1  && jPT2D == 1  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.040, 0.040);
            if ( iPT2D == 2  && jPT2D == 2  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.070, 0.070);
            if ( iPT2D == 3  && jPT2D == 3  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.050, 0.050);
            if ( iPT2D == 4  && jPT2D == 4  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.100, 0.100);
            if ( iPT2D == 5  && jPT2D == 5  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.080, 0.080);
            if ( iPT2D == 6  && jPT2D == 6  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.140, 0.140);
            if ( iPT2D == 7  && jPT2D == 7  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.050, 0.050);
            if ( iPT2D == 8  && jPT2D == 8  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.100, 0.100);
            if ( iPT2D == 9  && jPT2D == 9  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.130, 0.130);
            //
            //
            auto fVal = g2D_Syst_Err[iPT2D] -> GetErrorYlow(jPT2D);
            h2D_Syst_Ful3->SetBinContent( iPT2D+1, jPT2D+1, fVal) ;//
            //
            for ( Int_t iTer = 1; iTer <= nOption2; iTer++ )    {
                auto    fVariatin       =   g2D_Stat[iTer][iPT2D]       ->GetPointY     (jPT2D);
                auto    fVarError       =   g2D_Stat[iTer][iPT2D]       ->GetErrorYhigh (jPT2D);
                //
                auto    fYErrVrHig      =   g2D_Stat[iTer][iPT2D]       ->GetErrorYhigh (jPT2D) / (fStandard);
                auto    fYErrVrLow      =   g2D_Stat[iTer][iPT2D]       ->GetErrorYlow  (jPT2D) / (fStandard);
                //
                auto    fYFrac          =   fVariatin / fStandard - 1;
                auto    fYFracEHig      =   sqrt( fabs(fYErrorHig*fYErrorHig - fYErrVrHig*fYErrVrHig) );
                auto    fYFracELow      =   sqrt( fabs(fYErrorLow*fYErrorLow - fYErrVrLow*fYErrVrLow) );
                g2D_Stat_VarErr[iPT2D][iTer]    ->  SetPoint        (jPT2D, fXPoint,    fYFrac);
                g2D_Stat_VarErr[iPT2D][iTer]    ->  SetPointError   (jPT2D, fXError,    fXError,    fYFracEHig, fYFracELow);
                //
                if ( fVariatin - fStandard == 0 ) continue;
                if ( fBarlowCheck(fStandard,max(fStandard*fYErrorHig,fStandard*fYErrorLow),fVariatin,max(fStandard*fYErrVrHig,fStandard*fYErrVrLow)) )   continue;
                h2D_Syst_Bin[iPT2D][jPT2D]  ->  Fill(fYFrac);
                h2D_Syst_Bin_All            ->  Fill(fYFrac);
                hCheckFull2D[iPT2D]         ->  Fill(fXPoint,fYFrac);
            }
        }
    }
    //
    h2D_Syst_Bin_StSy->Write();
    h2D_Syst_Bin_All->Write();
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
        for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )  {
            h2D_Syst_Bin[iPT2D][jPT2D] ->  Write();
        }
    }
    //
    //------ Ratio Histograms ------//
    //
    TH1F * hRatio1DVar  =   new TH1F("hRatio1DVar",     "<Y^{SE}_{#phi}>",          nOptions,0,nOptions);
    TH1F * hRatio2DVar  =   new TH1F("hRatio2DVar",     "<Y^{SE}_{#phi#phi}>",      nOptions,0,nOptions);
    TH1F * hRatio1D     =   new TH1F("hRatio1D",        "",                         50,-.05,.05);
    TH1F * hRatio2D     =   new TH1F("hRatio2D",        "",                         50,-.05,.05);
    hRatio1DVar         -> GetXaxis() -> SetTitle("Variation");
    hRatio2DVar         -> GetXaxis() -> SetTitle("Variation");
    hRatio1D            -> GetXaxis() -> SetTitle("Variation");
    hRatio2D            -> GetXaxis() -> SetTitle("Variation");
    hRatio1DVar         -> GetXaxis() -> SetTitleOffset(1.5);
    hRatio2DVar         -> GetXaxis() -> SetTitleOffset(1.5);
    hRatio1DVar         -> GetYaxis() -> SetTitle("Fractional Deviation");
    hRatio2DVar         -> GetYaxis() -> SetTitle("Fractional Deviation");
    auto iBin = 1;
    for ( auto iName : sOptions )   {
        hRatio1DVar        ->GetXaxis()    ->  SetBinLabel(iBin,iName);
        hRatio2DVar        ->GetXaxis()    ->  SetBinLabel(iBin,iName);
        iBin++;
    }
    
    TH1F * fCheckRatio  =   new TH1F("fCheckRatio", "<Y_{#phi#phi}> / <Y_{#phi}>",nOptions,0,nOptions);
    TH1F * fCheckRatio2 =   new TH1F("fCheckRatio2","<Y_{#phi#phi}> / <Y_{#phi}>^{2}",nOptions,0,nOptions);
    fCheckRatio     -> GetXaxis() -> SetTitle("Variation");
    fCheckRatio2    -> GetXaxis() -> SetTitle("Variation");
    fCheckRatio     -> GetXaxis() -> SetTitleOffset(1.5);
    fCheckRatio2    -> GetXaxis() -> SetTitleOffset(1.5);
    fCheckRatio     -> GetYaxis() -> SetTitle("Fractional Deviation");
    fCheckRatio2    -> GetYaxis() -> SetTitle("Fractional Deviation");
    fCheckRatio->GetXaxis()->LabelsOption("v");
    fCheckRatio2->GetXaxis()->LabelsOption("v");
    TH1F * fCheckRati_  =   new TH1F("fCheckRati_", "",50,-.05,.05);
    TH1F * fCheckRati_2 =   new TH1F("fCheckRati_2","",50,-.05,.05);
    fCheckRati_     -> GetXaxis() -> SetTitle("Fractional Deviation");
    fCheckRati_2    -> GetXaxis() -> SetTitle("Fractional Deviation");
    for ( Int_t iTer = 1; iTer <= nOptions; iTer++ )  {
        auto f1DRefErr =  0.;
        auto f1DRefHst = new TH1F(*h1D_Syst[0]);
        f1DRefHst->Divide(hEFF_1D);
        auto f1DRefInt = f1DRefHst->IntegralAndError(-1,1000,f1DRefErr,"width");
        auto f1DTstErr =  0.;
        auto f1DTstHst = new TH1F(*h1D_Syst[iTer]);
        f1DTstHst->Divide(hEFF_1D);
        auto f1DTstInt = f1DTstHst->IntegralAndError(-1,1000,f1DTstErr,"width");
        auto f2DRefErr =  0.;
        auto f2DRefHst = new TH2F(*h2D_Syst[0]);
        f2DRefHst->Divide(hEFF_2D);
        auto f2DRefInt = f2DRefHst->IntegralAndError(-1,1000,-1,1000,f2DRefErr,"width");
        auto f2DTstErr =  0.;
        auto f2DTstHst = new TH2F(*h2D_Syst[iTer]);
        f2DTstHst->Divide(hEFF_2D);
        auto f2DTstInt = f2DTstHst->IntegralAndError(-1,1000,-1,1000,f2DTstErr,"width");
        auto fRatio1Ref =   f2DRefInt/f1DRefInt;
        auto fRatio1Trg =   f2DTstInt/f1DTstInt;
        auto fRatio2Ref =   f2DRefInt/(f1DRefInt*f1DRefInt);
        auto fRatio2Trg =   f2DTstInt/(f1DTstInt*f1DTstInt);
        fCheckRatio     ->  SetBinContent       (iTer,  fRatio1Trg/fRatio1Ref -1.);
        fCheckRatio2    ->  SetBinContent       (iTer,  fRatio2Trg/fRatio2Ref -1.);
        fCheckRati_     ->  Fill                (fRatio1Trg/fRatio1Ref -1.);
        fCheckRati_2    ->  Fill                (fRatio2Trg/fRatio2Ref -1.);
    }
    TCanvas * c1 = new TCanvas("","",1600,1600);
    c1->Divide(2,2);
    c1->cd(1);
    gStyle->SetOptStat(0);
    fCheckRatio->Draw();
    c1->cd(2);
    gStyle->SetOptStat(0);
    fCheckRatio2->Draw();
    c1->cd(3);
    gStyle->SetOptStat(0);
    fCheckRati_->Draw();
    c1->cd(4);
    gStyle->SetOptStat(0);
    fCheckRati_2->Draw();
    c1->SaveAs("./result/yield/ExtractionSystematics/1D_2D.pdf");
    delete c1;
    TH1F * fChec2Ratio  =   new TH1F("fCheckRatio", "",nOption2-nOptions,nOptions,nOption2);
    TH1F * fChec2Ratio2 =   new TH1F("fCheckRatio2","",nOption2-nOptions,nOptions,nOption2);
    TH1F * fChec2Rati_  =   new TH1F("fCheckRati_", "",50,-.05,.05);
    TH1F * fChec2Rati_2 =   new TH1F("fCheckRati_2","",50,-.05,.05);
    for ( Int_t iTer = nOptions+1; iTer <= nOption2; iTer++ )  {
        auto f1DRefErr =  0.;
        auto f1DRefInt = h1D_Syst[0]->IntegralAndError(-1,1000,f1DRefErr,"width");
        auto f1DTstErr =  0.;
        auto f1DTstInt = h1D_Syst[0]->IntegralAndError(-1,1000,f1DTstErr,"width");
        auto f2DRefErr =  0.;
        auto f2DRefInt = h2D_Syst[0]->IntegralAndError(-1,1000,-1,1000,f2DRefErr,"width");
        auto f2DTstErr =  0.;
        auto f2DTstInt = h2D_Syst[iTer]->IntegralAndError(-1,1000,-1,1000,f2DTstErr,"width");
        auto fRatio1Ref =   f2DRefInt/f1DRefInt;
        auto fRatio1Trg =   f2DTstInt/f1DTstInt;
        auto fRatio2Ref =   f2DRefInt/(f1DRefInt*f1DRefInt);
        auto fRatio2Trg =   f2DTstInt/(f1DTstInt*f1DTstInt);
        fChec2Ratio     ->  SetBinContent       (iTer-nOptions,  fRatio1Trg/fRatio1Ref -1.);
        fChec2Ratio2    ->  SetBinContent       (iTer-nOptions,  fRatio2Trg/fRatio2Ref -1.);
        fChec2Rati_     ->  Fill                (fRatio1Trg/fRatio1Ref -1.);
        fChec2Rati_2    ->  Fill                (fRatio2Trg/fRatio2Ref -1.);
    }
    TCanvas * c2 = new TCanvas("","",1600,1600);
    c2->Divide(2,2);
    c2->cd(1);
    gStyle->SetOptStat(0);
    fChec2Ratio->Draw();
    c2->cd(2);
    gStyle->SetOptStat(0);
    fChec2Ratio2->Draw();
    c2->cd(3);
    gStyle->SetOptStat(0);
    fChec2Rati_->Draw();
    c2->cd(4);
    gStyle->SetOptStat(0);
    fChec2Rati_2->Draw();
    c2->SaveAs("./result/yield/ExtractionSystematics/2D.pdf");
    delete c2;
    //
    gROOT->SetBatch(true);
    // Output File for Fit Check
    TFile*  outFileFi2  =   new TFile("./result/yield/ExtractionSystematics/ExtractionSystematics_MeanAndRMS.root","recreate");
    //
    TH1F*   h1D_Syst_Mean   =   new TH1F    ("h1D_Syst_Mean","h1D_Syst_Mean",nBinPT1D,fArrPT1D);
    TH1F*   h1D_Syst_RMS_   =   new TH1F    ("h1D_Syst_RMS_","h1D_Syst_RMS_",nBinPT1D,fArrPT1D);
    TH1F*   h1D_Syst_Full   =   new TH1F    ("h1D_Syst_Full","h1D_Syst_Full",nBinPT1D,fArrPT1D);
    TH1F*   h1D_Syst_Ful2   =   new TH1F    ("h1D_Syst_Ful2","h1D_Syst_Full",nBinPT1D,fArrPT1D);
    for ( Int_t iPT1D = 1; iPT1D <= nBinPT1D; iPT1D++ )  {
        h1D_Syst_Mean   ->SetBinContent (iPT1D,  fabs( h1D_Syst_Bin[iPT1D] ->  GetMean() ) );
        h1D_Syst_RMS_   ->SetBinContent (iPT1D,  h1D_Syst_Bin[iPT1D] ->  GetRMS());
        h1D_Syst_Full   ->SetBinContent (iPT1D,  h1D_Syst_Bin[iPT1D] ->  GetRMS() + fabs(h1D_Syst_Bin[iPT1D] ->  GetMean() ));
        h1D_Syst_Ful2   ->SetBinContent (iPT1D,  h1D_Syst[0]->GetBinContent(iPT1D));
        h1D_Syst_Ful2   ->SetBinError   (iPT1D,  h1D_Syst[0]->GetBinContent(iPT1D)*g1D_Syst_Err->GetErrorYlow(iPT1D));
    }
    h1D_Stat_Bin    ->  Write();
    h1D_Syst_Mean   ->  Write();
    h1D_Syst_RMS_   ->  Write();
    h1D_Syst_Full   ->  Write();
    h1D_Syst_Ful2   ->  Write();
    //
    TH2F*   h2D_Syst_Mean   =   new TH2F    ("h2D_Syst_Mean","h2D_Syst_Mean",               nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    TH2F*   h2D_Syst_RMS_   =   new TH2F    ("h2D_Syst_RMS_","h2D_Syst_RMS_",               nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    TH2F*   h2D_Syst_Full   =   new TH2F    ("h2D_Syst_Full","h2D_Syst_Full",               nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    TH2F*   h2D_Syst_Ful2   =   new TH2F    ("h2D_Syst_Ful2","h2D_Syst_Full",               nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
        for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )  {
            h2D_Syst_Mean   ->SetBinContent (iPT2D+1,jPT2D+1, fabs( h2D_Syst_Bin[iPT2D][jPT2D] ->  GetMean()) );
            h2D_Syst_RMS_   ->SetBinContent (iPT2D+1,jPT2D+1, h2D_Syst_Bin[iPT2D][jPT2D] ->  GetRMS());
            h2D_Syst_Full   ->SetBinContent (iPT2D+1,jPT2D+1, h2D_Syst_Bin[iPT2D][jPT2D] ->  GetRMS() + fabs(h2D_Syst_Bin[iPT2D][jPT2D] ->  GetMean() ));
            h2D_Syst_Ful2   ->SetBinContent (iPT2D+1,jPT2D+1, h2D_Syst[0]->GetBinContent(iPT2D+1,jPT2D+1));
            h2D_Syst_Ful2   ->SetBinError   (iPT2D+1,jPT2D+1, h2D_Syst[0]->GetBinContent(iPT2D+1,jPT2D+1)*g2D_Syst_Err[iPT2D]->GetErrorYlow(jPT2D));
        }
    }
    SetAxis(h2D_Stat_Bin,"PT 2D");
    h2D_Stat_Bin    ->  Write();
    SetAxis(h2D_Syst_Mean,"PT 2D");
    h2D_Syst_Mean   ->  Write();
    SetAxis(h2D_Syst_RMS_,"PT 2D");
    h2D_Syst_RMS_   ->  Write();
    SetAxis(h2D_Syst_Full,"PT 2D");
    h2D_Syst_Full   ->  Write();
    SetAxis(h2D_Syst_Ful2,"PT 2D");
    h2D_Syst_Ful2   ->  Write();
    SetAxis(h2D_Syst_Ful3,"PT 2D");
    h2D_Syst_Ful3   ->  Write();
    
    TCanvas* cCompare = new TCanvas("","",1800,600);
    cCompare->Divide(3,1);
    cCompare->cd(1);
    gPad->SetLogx();
    gPad->SetLogy();
    h2D_Syst_Full->Scale(100);
    h2D_Syst_Full->DrawCopy("colz text28");
    cCompare->cd(2);
    gPad->SetLogx();
    gPad->SetLogy();
    h2D_Syst_Ful3->Scale(100);
    h2D_Syst_Ful3->Draw("colz text28");
    cCompare->cd(3);
    gPad->SetLogx();
    gPad->SetLogy();
    h2D_Syst_Full->Add(h2D_Syst_Ful3,-1.);
    h2D_Syst_Full->Draw("colz text28");
    cCompare->SaveAs("dddd.pdf");
    delete cCompare;
    //
    TH1F * fChec3Rati_  =   new TH1F("Ratio Errors", "<Y_{#phi#phi}> / <Y_{#phi}>",500,0.,5.);
    fChec3Rati_->GetXaxis()->SetTitle("Relative Uncertainty (%)");
    TH1F * fChec4Rati_  =   new TH1F("Square Combine Errors", "",500,0.,5.);
    TH1F * fChec5Rati_  =   new TH1F("Ratio Errors", "<Y_{#phi#phi}> / <Y_{#phi}>^{2}",500,0.,5.);
    fChec5Rati_->GetXaxis()->SetTitle("Relative Uncertainty (%)");
    TH1F * fChec6Rati_  =   new TH1F("Square Combine Errors", "",500,0.,5.);
    TH1F * fChec7Rati_  =   new TH1F("Linear Combine Errors", "",500,0.,5.);
    TH1F * fChec8Rati_  =   new TH1F("Linear Combine Errors", "",500,0.,5.);
    TLegend * L1 = new TLegend(0.9,0.9,0.5,0.8);
    L1->AddEntry(fChec3Rati_,"Ratio Errors");
    L1->AddEntry(fChec4Rati_,"Square Combined Errors");
    L1->AddEntry(fChec7Rati_,"Linear Combined Errors");
    TCanvas * c3 = new TCanvas("","",1000,500);
    c3->Divide(2,1);
    c3->cd(1);
    gStyle->SetOptStat(0);
    fChec3Rati_->Fill(100*(fCheckRati_->  GetRMS() + fabs(fCheckRati_ ->  GetMean() )));
    fChec3Rati_->SetLineColor(kRed);
    fChec3Rati_->Draw("");
    auto h1D_Syst_Ful2_Err  = 0.;
    auto h1D_Syst_Ful2_Hst = new TH1F(*h1D_Syst_Ful2);
    h1D_Syst_Ful2_Hst->Divide(hEFF_1D);
    auto h1D_Syst_Ful2_Int  = h1D_Syst_Ful2_Hst->IntegralAndError(-1,1000,h1D_Syst_Ful2_Err,"width");
    auto h2D_Syst_Ful2_Err  = 0.;
    auto h2D_Syst_Ful2_Hst = new TH2F(*h2D_Syst_Ful2);
    h2D_Syst_Ful2_Hst->Divide(hEFF_2D);
    auto h2D_Syst_Ful2_Int  = h2D_Syst_Ful2_Hst->IntegralAndError(-1,1000,-1,1000,h2D_Syst_Ful2_Err,"width");
    auto hTarget            = sqrt( (h1D_Syst_Ful2_Err*h1D_Syst_Ful2_Err)/(h1D_Syst_Ful2_Int*h1D_Syst_Ful2_Int) + (h2D_Syst_Ful2_Err*h2D_Syst_Ful2_Err)/(h2D_Syst_Ful2_Int*h2D_Syst_Ful2_Int) );
    fChec4Rati_->Fill(100*hTarget);
    fChec4Rati_->SetLineColor(kBlue);
    fChec4Rati_->Draw("SAME");
    hTarget  = (h1D_Syst_Ful2_Err)/(h1D_Syst_Ful2_Int) + (h2D_Syst_Ful2_Err)/(h2D_Syst_Ful2_Int);
    fChec7Rati_->SetLineColor(kGreen-2);
    fChec7Rati_->Fill(100*hTarget);
    fChec7Rati_->Draw("SAME");
    L1->Draw("same");
    c3->cd(2);
    gStyle->SetOptStat(0);
    fChec5Rati_->Fill(100*(fCheckRati_2->  GetRMS() + fabs(fCheckRati_2 ->  GetMean() )) );
    fChec5Rati_->SetLineColor(kRed);
    fChec5Rati_->Draw("");
    hTarget            = sqrt( 4*(h1D_Syst_Ful2_Err*h1D_Syst_Ful2_Err)/(h1D_Syst_Ful2_Int*h1D_Syst_Ful2_Int) + (h2D_Syst_Ful2_Err*h2D_Syst_Ful2_Err)/(h2D_Syst_Ful2_Int*h2D_Syst_Ful2_Int) );
    fChec6Rati_->Fill(100*hTarget);
    fChec6Rati_->SetLineColor(kBlue);
    fChec6Rati_->Draw("SAME");
    hTarget  = 2*(h1D_Syst_Ful2_Err)/(h1D_Syst_Ful2_Int) + (h2D_Syst_Ful2_Err)/(h2D_Syst_Ful2_Int);
    fChec8Rati_->SetLineColor(kGreen-2);
    fChec8Rati_->Fill(100*hTarget);
    fChec8Rati_->Draw("SAME");
    L1->Draw("same");
    c3->SaveAs("./result/yield/ExtractionSystematics/12D_Check.pdf");
    delete c3;
    //
    TLatex         *latext              =   new TLatex();
    TCanvas        *cDrawComparison     =   new TCanvas("cDrawComparison","");
    TLegend        *cComparisonLegend   =   new TLegend(0.15,0.75,0.25,0.85);
    //
    gStyle                              ->  SetOptStat(0);
    //
    h1D_Stat_Bin                        ->  SetLineWidth(3);
    h1D_Stat_Bin                        ->  SetLineColor(kBlue);
    h1D_Stat_Bin                        ->  SetMinimum(0);
    h1D_Stat_Bin                        ->  SetMaximum(0.07);
    h1D_Syst_RMS_                       ->  SetLineWidth(3);
    h1D_Syst_RMS_                       ->  SetLineColor(kRed);
    //
    cComparisonLegend                   ->  SetLineColorAlpha(1,0.);
    cComparisonLegend                   ->  AddEntry(h1D_Stat_Bin,"Stat.","L");
    cComparisonLegend                   ->  AddEntry(h1D_Syst_RMS_,"Syst.","L");
    //
    h1D_Stat_Bin                        ->  Draw();
    h1D_Syst_RMS_                       ->  Draw("SAME");
    cComparisonLegend                   ->  Draw("SAME");
    cDrawComparison                     ->  SaveAs("./result/yield/ExtractionSystematics/_Syst_Stat_overimp.pdf");
    cDrawComparison                     ->  SaveAs("./result/yield/ExtractionSystematics/_Syst_Stat_overimp.png");
    cDrawComparison                     ->  Write();
    //
    delete cDrawComparison;
    delete cComparisonLegend;
    //
                    cDrawComparison     =   new TCanvas("cDrawComparison","");
    gStyle                              ->  SetOptStat(0);
    //
    TF1            *fFlatDist           =   new TF1("fFlatDist","pol0",-100.,100.);
    //
    h1D_Syst_Mean                       ->  SetLineWidth(3);
    h1D_Syst_Mean                       ->  SetLineColor(kBlue);
    h1D_Syst_Mean                       ->  Fit(fFlatDist,"IMREQ0S");
    //
                    cComparisonLegend   =   new TLegend(0.15,0.75,0.35,0.85);
    cComparisonLegend                   ->  SetLineColorAlpha(1,0.);
    cComparisonLegend                   ->  AddEntry(h1D_Syst_Mean,"Mean of Bin Distr.","L");
    //
    h1D_Syst_Mean                       ->  Draw();
    cComparisonLegend                   ->  Draw("SAME");
    fFlatDist                           ->  Draw("SAME");
    latext                              ->  DrawLatexNDC(0.6, 0.83, Form("MEAN:  %3f",fFlatDist->GetParameter(0)));
    latext                              ->  DrawLatexNDC(0.6, 0.75, Form("ERROR: %3f",fFlatDist->GetParError(0)));
    cDrawComparison                     ->  SaveAs("./result/yield/ExtractionSystematics/_Mean_1D.pdf");
    cDrawComparison                     ->  SaveAs("./result/yield/ExtractionSystematics/_Mean_1D.png");
    cDrawComparison                     ->  Write();
    //
    
    
    //
    TCanvas                *cDrawCollection = new TCanvas("","",1600,1600);
    cDrawCollection->Divide(4,3);
    
    auto Check              =   fMultipleError(g1D_Stat_Err,g1D_Syst_Err,g1D_Stat_VarErr,1,nOptions,sOptions);
    gPad->SetLogx(true);
    auto fMultiGrap1        =   (TMultiGraph*)Check ->  GetPrimitive   ("cDrawAllGraphs");
    fMultiGrap1             ->  SetMaximum(+0.15);
    fMultiGrap1             ->  SetMinimum(-0.15);
    SetAxis(fMultiGrap1,"PT 1D");
    fMultiGrap1             ->  GetYaxis()->SetTitle("Fractional Variation");
    fMultiGrap1             ->  SetTitle(Form("PID Systematic in 1D"));
    Check                   ->  SaveAs("./result/yield/ExtractionSystematics/SE_ERROR_1D.pdf");
    Check                   ->  SaveAs("./result/yield/ExtractionSystematics/SE_ERROR_1D.png");
    Check                   ->  Write();
    cDrawCollection         ->  cd(1);
    Check                   ->  DrawClonePad();
    gPad->SetLogx(true);
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
        auto Check2D        =   fMultipleError  (g2D_Stat_Err[iPT2D],g2D_Syst_Err[iPT2D],g2D_Stat_VarErr[iPT2D],1,nOption2,sOption2);
        gPad->SetLogx(true);
        auto fMultiGrap2    =   (TMultiGraph*)Check2D             ->  GetPrimitive   ("cDrawAllGraphs");
        fMultiGrap2         ->  SetMaximum(+0.5);
        fMultiGrap2         ->  SetMinimum(-0.5);
        SetAxis(fMultiGrap2,"PT DD");
        fMultiGrap2         ->  GetYaxis()->SetTitle("Fractional Variation");
        fMultiGrap2             ->  SetTitle(Form("PID Systematic in 2D, PT %.1f-%.1f",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1]));
        Check2D             ->  SaveAs  (Form("./result/yield/ExtractionSystematics/SE_ERROR_2D_PT_Bin_%.1f_%.1f.pdf",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1]));
        Check2D             ->  SaveAs  (Form("./result/yield/ExtractionSystematics/SE_ERROR_2D_PT_Bin_%.1f_%.1f.png",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1]));
        Check2D             ->  Write   ();
        cDrawCollection     ->  cd(2+iPT2D);
        gPad->SetLogx(true);
        Check2D             ->  DrawClonePad();
    }
    cDrawCollection                   ->  SaveAs("./result/yield/ExtractionSystematics/SE_ERROR_FULL.pdf");
    cDrawCollection                   ->  SaveAs("./result/yield/ExtractionSystematics/SE_ERROR_FULL.png");
    //
    delete  cDrawCollection;
    delete  cComparisonLegend;
    //
    cDrawCollection = new TCanvas("","",1600,1600);
    gPad->SetLogx(true);
    cDrawCollection->Divide(4,3);
    //
    THStack  *hSystStack_1D =   new THStack("","");
    h1D_Syst_Mean           ->  SetFillColorAlpha(2,0.1);
    h1D_Syst_Mean           ->  SetLineColorAlpha(2,1.0);
    h1D_Syst_Mean           ->  SetLineWidth(1.2);
    hSystStack_1D           ->  Add(h1D_Syst_Mean);
    h1D_Syst_RMS_           ->  SetFillColorAlpha(4,0.1);
    h1D_Syst_RMS_           ->  SetLineColorAlpha(4,1.0);
    h1D_Syst_RMS_           ->  SetLineWidth(1.2);
    hSystStack_1D           ->  Add(h1D_Syst_RMS_);
    //
    cComparisonLegend   =   new TLegend(0.15,0.75,0.35,0.85);
    cComparisonLegend       ->  SetLineColorAlpha(1,0.);
    cComparisonLegend       ->  SetFillColorAlpha(1,0.);
    
    TMultiGraph    *cDrawAllGraphs      =   new TMultiGraph("cDrawAllGraphs","");
    //
    g1D_Stat_Err                              ->  SetFillColorAlpha(kGray,0.75);
    g1D_Syst_Err                              ->  SetFillColorAlpha(kGray+2,0.5);
    //
    cDrawAllGraphs                      ->  Add         (g1D_Stat_Err,      "AE2");
    cDrawAllGraphs                      ->  Add         (g1D_Syst_Err,      "AE2");
    //
    //
    cComparisonLegend       ->  AddEntry(g1D_Stat_Err,"Stat. Err.","F");
    cComparisonLegend       ->  AddEntry(g1D_Syst_Err,"Syst. Err.","F");
    cComparisonLegend       ->  AddEntry(h1D_Syst_Mean,"Mean Contr.","F");
    cComparisonLegend       ->  AddEntry(h1D_Syst_RMS_,"RMS Contr.","F");
    //
    TCanvas     *cStackShow =   new TCanvas();
    gPad->SetLogx(true);
    cDrawAllGraphs     ->  Draw("ALP");
    cDrawAllGraphs     ->  GetYaxis()->SetTitle("Fractional Variation");
    cDrawAllGraphs     ->  SetMinimum(0.0);
    cDrawAllGraphs     ->  SetMaximum(0.1);
    SetAxis(cDrawAllGraphs,"PT 1D");
    hSystStack_1D           ->  Draw("same");
    cComparisonLegend       ->  Draw("same");
    
    cDrawCollection         ->  cd(1);
    cStackShow              ->  DrawClonePad();
    //
    cStackShow              ->  SaveAs("./result/yield/ExtractionSystematics/SE_ERROR_SYST_1D.pdf");
    cStackShow              ->  SaveAs("./result/yield/ExtractionSystematics/SE_ERROR_SYST_1D.png");
    //
    delete  cStackShow;
    delete  cComparisonLegend;
    //
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
        //
        cComparisonLegend   =   new TLegend(0.15,0.75,0.35,0.85);
        cComparisonLegend       ->  SetLineColorAlpha(1,0.);
        cComparisonLegend       ->  SetFillColorAlpha(1,0.);
        //
        auto    hSlice_Mean     =   h2D_Syst_Mean->ProjectionY("Mean",iPT2D+1,iPT2D+1);
        auto    hSlice_RMS_     =   h2D_Syst_RMS_->ProjectionY("RMS_",iPT2D+1,iPT2D+1);
        //
        THStack  *hSystStack_2D =   new THStack("","");
        hSlice_Mean             ->  SetFillColorAlpha(2,0.1);
        hSlice_Mean             ->  SetLineColorAlpha(2,1.0);
        hSlice_Mean             ->  SetLineWidth(1.2);
        hSystStack_2D           ->  Add(hSlice_Mean);
        hSlice_RMS_             ->  SetFillColorAlpha(4,0.1);
        hSlice_RMS_             ->  SetLineColorAlpha(4,1.0);
        hSlice_RMS_             ->  SetLineWidth(1.2);
        hSystStack_2D           ->  Add(hSlice_RMS_);
        //
        TMultiGraph    *cDrawAllGraphs      =   new TMultiGraph("cDrawAllGraphs","");
        //
        g2D_Stat_Err[iPT2D]                              ->  SetFillColorAlpha(kGray,0.75);
        g2D_Syst_Err[iPT2D]                              ->  SetFillColorAlpha(kGray+2,0.5);
        //
        cDrawAllGraphs                      ->  Add         (g2D_Stat_Err[iPT2D],      "AE2");
        cDrawAllGraphs                      ->  Add         (g2D_Syst_Err[iPT2D],      "AE2");
        //
        //
        cComparisonLegend       ->  AddEntry(g2D_Stat_Err[iPT2D],"Stat. Err.","F");
        cComparisonLegend       ->  AddEntry(g2D_Syst_Err[iPT2D],"Syst. Err.","F");
        cComparisonLegend       ->  AddEntry(hSlice_Mean,"Mean Contr.","F");
        cComparisonLegend       ->  AddEntry(hSlice_RMS_,"RMS Contr.","F");
        //
                    cStackShow  =   new TCanvas("StackShow","StackShow");
        gPad->SetLogx(true);
        cDrawAllGraphs     ->  Draw("ALP");
        cDrawAllGraphs     ->  GetYaxis()->SetTitle("Fractional Variation");
        cDrawAllGraphs     ->  SetMinimum(0.0);
        cDrawAllGraphs     ->  SetMaximum(0.5);
        SetAxis(cDrawAllGraphs,"PT DD");
        hSystStack_2D           ->  Draw("same");
        cComparisonLegend       ->  Draw("same");
        //
        cDrawCollection         ->  cd(2+iPT2D);
        cStackShow              ->  DrawClonePad();
        //
        cStackShow             ->  SaveAs  (Form("./result/yield/ExtractionSystematics/SE_ERROR_SYST_2D_PT_Bin_%.1f_%.1f.pdf",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1]));
        cStackShow             ->  SaveAs  (Form("./result/yield/ExtractionSystematics/SE_ERROR_SYST_2D_PT_Bin_%.1f_%.1f.png",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1]));
        //
        delete cStackShow;
        delete cComparisonLegend;
    }
    cDrawCollection                   ->  SaveAs("./result/yield/ExtractionSystematics/SE_ERROR_SYST_FULL.pdf");
    cDrawCollection                   ->  SaveAs("./result/yield/ExtractionSystematics/SE_ERROR_SYST_FULL.png");
    //
    delete  cDrawCollection;
    //
    cDrawCollection = new TCanvas("","",1600,1600);
    //
    gPad->SetLogx(true);
    hCheckFull1D->Draw("colz");
    cDrawCollection                   ->  SaveAs("./result/yield/ExtractionSystematics/SE_ERROR_SYST_FULL_1D.pdf");
    cDrawCollection                   ->  SaveAs("./result/yield/ExtractionSystematics/SE_ERROR_SYST_FULL_1D.png");
    //
    delete  cDrawCollection;
    //
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
        cDrawCollection = new TCanvas("","",1600,1600);
        //
        gPad->SetLogx(true);
        hCheckFull2D[iPT2D]->Draw("colz");
        cDrawCollection                   ->  SaveAs(Form("./result/yield/ExtractionSystematics/SE_ERROR_SYST_FULL_2D_%i.pdf",iPT2D));
        cDrawCollection                   ->  SaveAs(Form("./result/yield/ExtractionSystematics/SE_ERROR_SYST_FULL_2D_%i.png",iPT2D));
        //
        delete  cDrawCollection;
    }
    //
    // >-> Close input File
    //
    outFileFit->Close();
    outFileFi2->Close();
    //
    for ( Int_t iTer = 0; iTer <= nOptions; iTer++ )    {
        insFileH1D[iTer]->Close();
    }
    //
    for ( Int_t iTer = 0; iTer < nOption2; iTer++ )    {
        insFileH2D[iTer]->Close();
    }
    //
    gROOT->SetBatch(false);
    return;
}

*/
