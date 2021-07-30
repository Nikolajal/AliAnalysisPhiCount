// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../inc/AliAnalysisPhiPair.h"
#include "RooMsgService.h"

void Analysis_PID ( bool fSilent = false)
{
    TFile **insFileH1D  =   new TFile*  [nPIDFiles+1];
    
    //Recovering histograms
    TH1F**  h1D_Syst    =   new TH1F*   [nPIDFiles+1];
    
    insFileH1D[0]   =   new TFile   (Form("result/yield/PIDSystematics/Standard/SCHistograms.root"));
    hName           =   Form("hRES_1D_Stat");
    h1D_Syst[0]     =   (TH1F*)(insFileH1D[0]->Get(hName));
    
    for ( Int_t iTer = 1; iTer <= nPIDFiles; iTer++ )
    {
        insFileH1D[iTer]=   new TFile   (Form("result/yield/PIDSystematics/PID_%i/SCHistograms.root",iTer));
        hName           =   Form("hRES_1D_Stat");
        h1D_Syst[iTer]  =   (TH1F*)(insFileH1D[iTer]->Get(hName));
        //h1D_Syst[iTer]      ->Scale(1.,"width");
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
    TFile*  outFileFit  =   new TFile("./result/yield/PIDSystematics/ExtractionSystematics_CheckRatioAndBins.root","recreate");
    
    //------ 1D Histograms ------//
    //
    TGraphAsymmErrors     **g1D_Stat            =   new TGraphAsymmErrors  *[nPIDFiles+1];
    for ( Int_t iTer = 0; iTer <= nPIDFiles; iTer++ )    {
        auto fCheck = new TH1F (*h1D_Syst[iTer]);
        fCheck->Divide(h1D_Syst[iTer],h1D_Syst[0]);
        fCheck->Write();
        g1D_Stat[iTer]     =   fTH1_to_TGAsymmErrors(h1D_Syst[iTer]);
    }
    //
    TGraphAsymmErrors     **g1D_Stat_VarErr     =   new TGraphAsymmErrors  *[nPIDFiles+1];
    for ( Int_t iTer = 1; iTer <= nPIDFiles; iTer++ )    {
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
        if ( iPT1D == 0  || iPT1D == 19 ) g1D_Syst_Err            ->  SetPointError   (iPT1D, fXError,    fXError,    0.020, 0.020);
        if ( iPT1D >= 1  && iPT1D <= 1  ) g1D_Syst_Err            ->  SetPointError   (iPT1D, fXError,    fXError,    0.015, 0.015);
        if ( iPT1D >= 2  && iPT1D <= 2  ) g1D_Syst_Err            ->  SetPointError   (iPT1D, fXError,    fXError,    0.010, 0.010);
        if ( iPT1D >= 3  && iPT1D <= 18 ) g1D_Syst_Err            ->  SetPointError   (iPT1D, fXError,    fXError,    0.005, 0.005);
        //
        for ( Int_t iTer = 1; iTer <= nPIDFiles; iTer++ )    {
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
    
    TGraphAsymmErrors    ***g2D_Stat            =   new TGraphAsymmErrors  **[nPIDFiles+1];
    for ( Int_t iTer = 0; iTer <= nPIDFiles; iTer++ )    {
        g2D_Stat[iTer] = new TGraphAsymmErrors* [nBinPT2D];
        for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )    {
            hName           =   Form("hRES_2D_Cond1_Stat_%i",iPT2D);
            g2D_Stat[iTer][iPT2D]  =   fTH1_to_TGAsymmErrors((TH1F*)(insFileH1D[iTer]->Get(hName)));
            
            /*
            auto fCheck = new TH2F (*h2D_Syst[iTer]);
            fCheck->Divide(h2D_Syst[iTer],h2D_Syst[0]);
            fCheck->Write();
            g2D_Stat[iTer]     =   fTH2_to_TGAsymmErrors(h2D_Syst[iTer]);
            fTH1_to_TGAsymmErrors
            */
            
        }
    }
    //
    TGraphAsymmErrors    ***g2D_Stat_VarErr        =   new TGraphAsymmErrors **[nBinPT2D];
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )    {
        g2D_Stat_VarErr[iPT2D]      =   new TGraphAsymmErrors   *[nPIDFiles+1];
        for ( Int_t iTer = 1; iTer <= nPIDFiles; iTer++ )    {
            g2D_Stat_VarErr[iPT2D][iTer]   =   new TGraphAsymmErrors();
        }
    }
    //
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
            if ( iPT2D == 0  || jPT2D == 0  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.050, 0.050);
            if ( iPT2D == 9  || jPT2D == 9  ) g2D_Syst_Err[iPT2D]            ->  SetPointError   (jPT2D, fXError,    fXError,    0.200, 0.200);
            //
            for ( Int_t iTer = 1; iTer <= nPIDFiles; iTer++ )    {
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
    TH1F * fCheckRatio  =   new TH1F("fCheckRatio", "",nPIDFiles,0,nPIDFiles);
    TH1F * fCheckRatio2 =   new TH1F("fCheckRatio2","",nPIDFiles,0,nPIDFiles);
    auto iBin = 1;
    for ( auto iName : sOptions )   {
        fCheckRatio     ->GetXaxis()    ->  SetBinLabel(iBin,iName);
        fCheckRatio2    ->GetXaxis()    ->  SetBinLabel(iBin,iName);
        iBin++;
    }
    fCheckRatio->GetXaxis()->LabelsOption("v");
    fCheckRatio2->GetXaxis()->LabelsOption("v");
    TH1F * fCheckRati_  =   new TH1F("fCheckRati_", "",50,-.05,.05);
    TH1F * fCheckRati_2 =   new TH1F("fCheckRati_2","",50,-.05,.05);
    for ( Int_t iTer = 1; iTer <= nPIDFiles; iTer++ )  {
        auto f1DRefErr =  0.;
        auto f1DRefInt = h1D_Syst[0]->IntegralAndError(-1,1000,f1DRefErr,"width");
        auto f1DTstErr =  0.;
        auto f1DTstInt = h1D_Syst[iTer]->IntegralAndError(-1,1000,f1DTstErr,"width");
        auto f2DRefErr =  0.;
        auto f2DRefInt = 0.;
        auto f2DTstErr =  0.;
        auto f2DTstInt = 0.;
        for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )  {
            auto fErr1 = 0.;
            auto fErr2 = 0.;
            hName           =   Form("hRES_2D_Cond1_Stat_%i",jPT2D);
            f2DRefInt += ((TH1F*)(insFileH1D[0]   ->Get(hName)))->IntegralAndError(-1,1000,fErr1,"width");
            f2DRefErr += fErr1;
            f2DTstInt += ((TH1F*)(insFileH1D[iTer]->Get(hName)))->IntegralAndError(-1,1000,fErr2,"width");
            f2DTstErr += fErr2;
        }
        auto fRatio1Ref =   f2DRefInt/f1DRefInt;
        auto fRatio1Trg =   f2DTstInt/f1DTstInt;
        auto fRatio2Ref =   f2DRefInt/(f1DRefInt*f1DRefInt);
        auto fRatio2Trg =   f2DTstInt/(f1DTstInt*f1DTstInt);
        fCheckRatio     ->  SetBinContent       (iTer,  fRatio1Trg/fRatio1Ref -1.);
        fCheckRatio2    ->  SetBinContent       (iTer,  fRatio2Trg/fRatio2Ref -1.);
        fCheckRati_     ->  Fill                (fRatio1Trg/fRatio1Ref -1.);
        fCheckRati_2    ->  Fill                (fRatio2Trg/fRatio2Ref -1.);
    }
    TCanvas * c1 = new TCanvas("","",1000,1000);
    c1->Divide(2,2);
    c1->cd(1);
    gStyle->SetOptStat(111111);
    fCheckRatio->Draw();
    c1->cd(2);
    gStyle->SetOptStat(111111);
    fCheckRatio2->Draw();
    c1->cd(3);
    gStyle->SetOptStat(111111);
    fCheckRati_->Draw();
    c1->cd(4);
    gStyle->SetOptStat(111111);
    fCheckRati_2->Draw();
    c1->SaveAs("./result/yield/PIDSystematics/1D_2D.pdf");
    delete c1;
    /*
    TH1F * fChec2Ratio  =   new TH1F("fCheckRatio", "",nPIDFiles-nPIDFiles,nPIDFiles,nPIDFiles);
    TH1F * fChec2Ratio2 =   new TH1F("fCheckRatio2","",nPIDFiles-nPIDFiles,nPIDFiles,nPIDFiles);
    TH1F * fChec2Rati_  =   new TH1F("fCheckRati_", "",50,-.05,.05);
    TH1F * fChec2Rati_2 =   new TH1F("fCheckRati_2","",50,-.05,.05);
    for ( Int_t iTer = nPIDFiles+1; iTer <= nPIDFiles; iTer++ )  {
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
        fChec2Ratio     ->  SetBinContent       (iTer-nPIDFiles,  fRatio1Trg/fRatio1Ref -1.);
        fChec2Ratio2    ->  SetBinContent       (iTer-nPIDFiles,  fRatio2Trg/fRatio2Ref -1.);
        fChec2Rati_     ->  Fill                (fRatio1Trg/fRatio1Ref -1.);
        fChec2Rati_2    ->  Fill                (fRatio2Trg/fRatio2Ref -1.);
    }
    TCanvas * c2 = new TCanvas("","",1000,1000);
    c2->Divide(2,2);
    c2->cd(1);
    gStyle->SetOptStat(111111);
    fChec2Ratio->Draw();
    c2->cd(2);
    gStyle->SetOptStat(111111);
    fChec2Ratio2->Draw();
    c2->cd(3);
    gStyle->SetOptStat(111111);
    fChec2Rati_->Draw();
    c2->cd(4);
    gStyle->SetOptStat(111111);
    fChec2Rati_2->Draw();
    c2->SaveAs("./result/yield/PIDSystematics/2D.pdf");
    delete c2;
     */
    //
    gROOT->SetBatch(true);
    // Output File for Fit Check
    TFile*  outFileFi2  =   new TFile("./result/yield/PIDSystematics/ExtractionSystematics_MeanAndRMS.root","recreate");
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
    
    TH2F*   h2D_Syst_Mean   =   new TH2F    ("h2D_Syst_Mean","h2D_Syst_Mean",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    TH2F*   h2D_Syst_RMS_   =   new TH2F    ("h2D_Syst_RMS_","h2D_Syst_RMS_",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    TH2F*   h2D_Syst_Full   =   new TH2F    ("h2D_Syst_Full","h2D_Syst_Full",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    TH2F*   h2D_Syst_Ful2   =   new TH2F    ("h2D_Syst_Ful2","h2D_Syst_Full",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
        for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )  {
            h2D_Syst_Mean   ->SetBinContent (iPT2D+1,jPT2D+1, fabs( h2D_Syst_Bin[iPT2D][jPT2D] ->  GetMean()) );
            h2D_Syst_RMS_   ->SetBinContent (iPT2D+1,jPT2D+1, h2D_Syst_Bin[iPT2D][jPT2D] ->  GetRMS());
            h2D_Syst_Full   ->SetBinContent (iPT2D+1,jPT2D+1, h2D_Syst_Bin[iPT2D][jPT2D] ->  GetRMS() + fabs(h2D_Syst_Bin[iPT2D][jPT2D] ->  GetMean() ));
            
            hName           =   Form("hRES_2D_Cond1_Stat_%i",iPT2D);
            h2D_Syst_Ful2   ->SetBinContent (iPT2D+1,jPT2D+1, ((TH1F*)(insFileH1D[0]   ->Get(hName)))->GetBinContent(jPT2D+1));
            h2D_Syst_Ful2   ->SetBinError   (iPT2D+1,jPT2D+1, ((TH1F*)(insFileH1D[0]   ->Get(hName)))->GetBinContent(jPT2D+1)*g2D_Syst_Err[iPT2D]->GetErrorYlow(jPT2D));
        }
    }
    h2D_Stat_Bin    ->  Write();
    h2D_Syst_Mean   ->  Write();
    h2D_Syst_RMS_   ->  Write();
    h2D_Syst_Full   ->  Write();
    h2D_Syst_Ful2   ->  Write();
    //
    TH1F * fChec3Rati_  =   new TH1F("Ratio Errors", "<Y_{#phi#phi}> / <Y_{#phi}>",500,0.,.05);
    TH1F * fChec4Rati_  =   new TH1F("Square Combine Errors", "",500,0.,.05);
    TH1F * fChec5Rati_  =   new TH1F("Ratio Errors", "<Y_{#phi#phi}> / <Y_{#phi}>^{2}",500,0.,.05);
    TH1F * fChec6Rati_  =   new TH1F("Square Combine Errors", "",500,0.,.05);
    TH1F * fChec7Rati_  =   new TH1F("Linear Combine Errors", "",500,0.,.05);
    TH1F * fChec8Rati_  =   new TH1F("Linear Combine Errors", "",500,0.,.05);
    TLegend * L1 = new TLegend(0.9,0.9,0.5,0.8);
    L1->AddEntry(fChec3Rati_,"Ratio Errors");
    L1->AddEntry(fChec4Rati_,"Square Combined Errors");
    L1->AddEntry(fChec7Rati_,"Linear Combined Errors");
    TCanvas * c3 = new TCanvas("","",1000,500);
    c3->Divide(2,1);
    c3->cd(1);
    gStyle->SetOptStat(0);
    fChec3Rati_->Fill(fCheckRati_->  GetRMS() + fabs(fCheckRati_ ->  GetMean() ));
    fChec3Rati_->SetLineColor(kRed);
    fChec3Rati_->Draw("");
    auto h1D_Syst_Ful2_Err  = 0.;
    auto h1D_Syst_Ful2_Int  = h1D_Syst_Ful2->IntegralAndError(-1,1000,h1D_Syst_Ful2_Err,"width");
    auto h2D_Syst_Ful2_Err  = 0.;
    auto h2D_Syst_Ful2_Int  = h2D_Syst_Ful2->IntegralAndError(-1,1000,-1,1000,h2D_Syst_Ful2_Err,"width");
    auto hTarget            = sqrt( (h1D_Syst_Ful2_Err*h1D_Syst_Ful2_Err)/(h1D_Syst_Ful2_Int*h1D_Syst_Ful2_Int) + (h2D_Syst_Ful2_Err*h2D_Syst_Ful2_Err)/(h2D_Syst_Ful2_Int*h2D_Syst_Ful2_Int) );
    fChec4Rati_->Fill(hTarget);
    fChec4Rati_->SetLineColor(kBlue);
    fChec4Rati_->Draw("SAME");
    hTarget  = (h1D_Syst_Ful2_Err)/(h1D_Syst_Ful2_Int) + (h2D_Syst_Ful2_Err)/(h2D_Syst_Ful2_Int);
    fChec7Rati_->SetLineColor(kGreen-2);
    fChec7Rati_->Fill(hTarget);
    fChec7Rati_->Draw("SAME");
    L1->Draw("same");
    c3->cd(2);
    gStyle->SetOptStat(0);
    fChec5Rati_->Fill(fCheckRati_2->  GetRMS() + fabs(fCheckRati_2 ->  GetMean() ));
    fChec5Rati_->SetLineColor(kRed);
    fChec5Rati_->Draw("");
    hTarget            = sqrt( 4*(h1D_Syst_Ful2_Err*h1D_Syst_Ful2_Err)/(h1D_Syst_Ful2_Int*h1D_Syst_Ful2_Int) + (h2D_Syst_Ful2_Err*h2D_Syst_Ful2_Err)/(h2D_Syst_Ful2_Int*h2D_Syst_Ful2_Int) );
    fChec6Rati_->Fill(hTarget);
    fChec6Rati_->SetLineColor(kBlue);
    fChec6Rati_->Draw("SAME");
    hTarget  = 2*(h1D_Syst_Ful2_Err)/(h1D_Syst_Ful2_Int) + (h2D_Syst_Ful2_Err)/(h2D_Syst_Ful2_Int);
    fChec8Rati_->SetLineColor(kGreen-2);
    fChec8Rati_->Fill(hTarget);
    fChec8Rati_->Draw("SAME");
    L1->Draw("same");
    c3->SaveAs("./result/yield/PIDSystematics/12D_Check.pdf");
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
    cDrawComparison                     ->  SaveAs("./result/yield/PIDSystematics/_Syst_Stat_overimp.pdf");
    cDrawComparison                     ->  SaveAs("./result/yield/PIDSystematics/_Syst_Stat_overimp.png");
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
    cDrawComparison                     ->  SaveAs("./result/yield/PIDSystematics/_Mean_1D.pdf");
    cDrawComparison                     ->  SaveAs("./result/yield/PIDSystematics/_Mean_1D.png");
    cDrawComparison                     ->  Write();
    //
    
    
    //
    TCanvas                *cDrawCollection = new TCanvas("","",1000,1000);
    cDrawCollection->Divide(4,3);
    
    auto Check              =   fMultipleError(g1D_Stat_Err,g1D_Syst_Err,g1D_Stat_VarErr,1,nPIDFiles,sOptiPID);
    gPad->SetLogx(true);
    auto fMultiGrap1        =   (TMultiGraph*)Check ->  GetPrimitive   ("cDrawAllGraphs");
    fMultiGrap1             ->  SetMaximum(+0.025);
    fMultiGrap1             ->  SetMinimum(-0.025);
    SetAxis(fMultiGrap1,"PT 1D");
    fMultiGrap1             ->  GetYaxis()->SetTitle("Fractional Variation");
    fMultiGrap1             ->  SetTitle(Form("PID Systematic in 1D"));
    Check                   ->  SaveAs("./result/yield/PIDSystematics/SE_ERROR_1D.pdf");
    Check                   ->  SaveAs("./result/yield/PIDSystematics/SE_ERROR_1D.png");
    Check                   ->  Write();
    cDrawCollection         ->  cd(1);
    Check                   ->  DrawClonePad();
    gPad->SetLogx(true);
    
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
        auto Check2D        =   fMultipleError  (g2D_Stat_Err[iPT2D],g2D_Syst_Err[iPT2D],g2D_Stat_VarErr[iPT2D],1,nPIDFiles,sOptiPID);
        gPad->SetLogx(true);
        auto fMultiGrap2    =   (TMultiGraph*)Check2D             ->  GetPrimitive   ("cDrawAllGraphs");
        fMultiGrap2         ->  SetMaximum(+0.5);
        fMultiGrap2         ->  SetMinimum(-0.5);
        SetAxis(fMultiGrap2,"PT DD");
        fMultiGrap2         ->  GetYaxis()->SetTitle("Fractional Variation");
        fMultiGrap2             ->  SetTitle(Form("PID Systematic in 2D, PT %.1f-%.1f",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1]));
        Check2D             ->  SaveAs  (Form("./result/yield/PIDSystematics/SE_ERROR_2D_PT_Bin_%.1f_%.1f.pdf",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1]));
        Check2D             ->  SaveAs  (Form("./result/yield/PIDSystematics/SE_ERROR_2D_PT_Bin_%.1f_%.1f.png",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1]));
        Check2D             ->  Write   ();
        cDrawCollection     ->  cd(2+iPT2D);
        gPad->SetLogx(true);
        Check2D             ->  DrawClonePad();
    }
     
    cDrawCollection                   ->  SaveAs("./result/yield/PIDSystematics/SE_ERROR_FULL.pdf");
    cDrawCollection                   ->  SaveAs("./result/yield/PIDSystematics/SE_ERROR_FULL.png");
    //
    delete  cDrawCollection;
    delete  cComparisonLegend;
    //
    cDrawCollection = new TCanvas("","",1000,1000);
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
    SetAxis(cDrawAllGraphs,"PT DD");
    hSystStack_1D           ->  Draw("same");
    cComparisonLegend       ->  Draw("same");
    
    cDrawCollection         ->  cd(1);
    cStackShow              ->  DrawClonePad();
    //
    cStackShow              ->  SaveAs("./result/yield/PIDSystematics/SE_ERROR_SYST_1D.pdf");
    cStackShow              ->  SaveAs("./result/yield/PIDSystematics/SE_ERROR_SYST_1D.png");
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
        cStackShow             ->  SaveAs  (Form("./result/yield/PIDSystematics/SE_ERROR_SYST_2D_PT_Bin_%.1f_%.1f.pdf",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1]));
        cStackShow             ->  SaveAs  (Form("./result/yield/PIDSystematics/SE_ERROR_SYST_2D_PT_Bin_%.1f_%.1f.png",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1]));
        //
        delete cStackShow;
        delete cComparisonLegend;
    }
    cDrawCollection                   ->  SaveAs("./result/yield/PIDSystematics/SE_ERROR_SYST_FULL.pdf");
    cDrawCollection                   ->  SaveAs("./result/yield/PIDSystematics/SE_ERROR_SYST_FULL.png");
    //
    delete  cDrawCollection;
    //
    cDrawCollection = new TCanvas("","",1000,1000);
    //
    gPad->SetLogx(true);
    hCheckFull1D->Draw("colz");
    cDrawCollection                   ->  SaveAs("./result/yield/PIDSystematics/SE_ERROR_SYST_FULL_1D.pdf");
    cDrawCollection                   ->  SaveAs("./result/yield/PIDSystematics/SE_ERROR_SYST_FULL_1D.png");
    //
    delete  cDrawCollection;
    //
    
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
        cDrawCollection = new TCanvas("","",1000,1000);
        //
        gPad->SetLogx(true);
        hCheckFull2D[iPT2D]->Draw("colz");
        cDrawCollection                   ->  SaveAs(Form("./result/yield/PIDSystematics/SE_ERROR_SYST_FULL_2D_%i.pdf",iPT2D));
        cDrawCollection                   ->  SaveAs(Form("./result/yield/PIDSystematics/SE_ERROR_SYST_FULL_2D_%i.png",iPT2D));
        //
        delete  cDrawCollection;
    }
    //
    // >-> Close input File
    //
    outFileFit->Close();
    outFileFi2->Close();
    //
    for ( Int_t iTer = 0; iTer <= nPIDFiles; iTer++ )    {
        insFileH1D[iTer]->Close();
    }
    //
    gROOT->SetBatch(false);
    return;
}
