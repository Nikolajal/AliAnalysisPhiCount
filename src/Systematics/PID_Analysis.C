// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../inc/AliAnalysisPhiPair.h"
#include "RooMsgService.h"

void PID_Analysis ( bool fSilent = false)
{
    TFile             **insFileHXD      =   new TFile              *[nPIDFiles+1];
    
    //Recovering histograms
    TGraphAsymmErrors **g1D_Stat        =   new TGraphAsymmErrors  *[nPIDFiles+1];
    TGraphAsymmErrors***g2D_Stat        =   new TGraphAsymmErrors **[nPIDFiles+1];
    
    for ( Int_t iFls = 0; iFls <= nPIDFiles; iFls++ )    {
        insFileHXD[iFls]    =   new TFile(Form("./result/Syst_PID_%i/SCHistograms.root",iFls));
        hName               =   Form("gRES_1D_Stat");
        g1D_Stat[iFls]      =   (TGraphAsymmErrors*)(insFileHXD[iFls]->Get(hName));
        
        g2D_Stat[iFls]      =   new TGraphAsymmErrors   *[nBinPT2D];
        for ( Int_t iPT2D = 0; iPT2D <= nBinPT2D; iPT2D++ )    {
            hName                   =   Form("gRES_2D_Stat_%i",iPT2D);
            g2D_Stat[iFls][iPT2D]   =   (TGraphAsymmErrors*)(insFileHXD[iFls]->Get(hName));
        }
    }
    
    //---------------------//
    //  Setting up output  //
    //---------------------//
    
    //--- Generating the binning array ---//
    fSetBinPT1D();
    fSetBinIM1D();
    fSetBinPT2D();
    fSetBinIM2D();
    
    // Creating the histograms-------------------------------------------------------------------------------
    //
    hName               =   Form("h1D_Syst_Bin_StSy");
    hTitle              =   Form("Percentage variation of raw yield for bin of PT [%.2f#;%.2f]",fArrPT1D[0],fArrPT1D[nBinPT1D]);
    TH1F   *h1D_Syst_Bin_StSy   =   new TH1F (hName,hTitle,100,-.5,.5);
    //
    hName               =   Form("h2D_Syst_Bin_StSy");
    hTitle              =   Form("Percentage variation of raw yield for bin of PT [%.2f#;%.2f]",fArrPT1D[0],fArrPT1D[nBinPT1D]);
    TH1F   *h2D_Syst_Bin_StSy   =   new TH1F (hName,hTitle,400,-2.,2.);
    //
    hName               =   Form("h1D_Stat_Bin");
    hTitle              =   Form("h1D_Stat_Bin");
    TH1F   *h1D_Stat_Bin    =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    //
    TH1F  **h1D_Syst_Bin    =   new TH1F   *[nBinPT1D+1];
    hName               =   Form("h1D_Syst_Bin_PT_%.2f_%.2f",fArrPT1D[0],fArrPT1D[nBinPT1D]);
    hTitle              =   Form("Percentage variation of raw yield for bin of PT [%.2f#;%.2f]",fArrPT1D[0],fArrPT1D[nBinPT1D]);
    h1D_Syst_Bin[0]  =   new TH1F    (hName,hTitle,400,-2.,2.);
    h1D_Syst_Bin[0]  ->SetTitle(Form("Percentage variation of raw yield for bin %i",-1));
    h1D_Syst_Bin[0]  ->GetXaxis()->SetTitle("Percentage variation");
    for ( Int_t iAll = 1; iAll <= nBinPT1D; iAll++ )
    {
        hName               =   Form("h1D_Syst_Bin_PT_%.2f_%.2f",fArrPT1D[iAll-1],fArrPT1D[iAll]);
        hTitle              =   Form("Percentage variation of raw yield for bin of PT [%.2f#;%.2f]",fArrPT1D[iAll-1],fArrPT1D[iAll]);
        h1D_Syst_Bin[iAll]  =   new TH1F    (hName,hTitle,400,-2.,2.);
        h1D_Syst_Bin[iAll]  ->SetTitle(Form("Percentage variation of raw yield for bin %i",iAll));
        h1D_Syst_Bin[iAll]  ->GetXaxis()->SetTitle("Percentage variation");
    }
    //
    TH1F   *h2D_Syst_Bin_All;
    hName                       =   Form("h2D_Syst_Bin_PT_%.2f_%.2f_%.2f_%.2f",fArrPT2D[0],fArrPT2D[nBinPT2D],fArrPT2D[0],fArrPT2D[nBinPT2D]);
    hTitle                      =   Form("Percentage variation of raw yield for bin of PT [%.2f#;%.2f] [%.2f#;%.2f]",fArrPT2D[0],fArrPT2D[nBinPT2D],fArrPT2D[0],fArrPT2D[nBinPT2D]);
    h2D_Syst_Bin_All            =   new TH1F    (hName,hTitle,100,0.,2.);
    h2D_Syst_Bin_All            ->  SetTitle(Form("Percentage variation of raw yield for bin %i",-1));
    h2D_Syst_Bin_All            ->  GetXaxis()  ->  SetTitle("Percentage variation");
    //
    TH2F   *h2D_Stat_Bin;
    hName                       =   Form("h2D_Stat_Bin");
    hTitle                      =   Form("h2D_Stat_Bin");
    h2D_Stat_Bin                =   new TH2F    (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    //
    TH1F ***h2D_Syst_Bin    =   new TH1F  **[nBinPT2D];
    for ( Int_t iAll = 0; iAll < nBinPT2D; iAll++ ) {
        h2D_Syst_Bin[iAll]  =   new TH1F   *[nBinPT2D];
        for ( Int_t jAll = 0; jAll < nBinPT2D; jAll++ ) {
            hName                       =   Form("h2D_Syst_Bin_PT_%.2f_%.2f_%.2f_%.2f",fArrPT2D[iAll],fArrPT2D[iAll+1],fArrPT2D[jAll],fArrPT2D[jAll+1]);
            hTitle                      =   Form("Percentage variation of raw yield for bin of PT [%.2f#;%.2f] [%.2f#;%.2f]",fArrPT2D[iAll],fArrPT2D[iAll+1],fArrPT2D[jAll],fArrPT2D[jAll+1]);
            h2D_Syst_Bin[iAll][jAll]    =   new TH1F    (hName,hTitle,400,-2.,2.);
            h2D_Syst_Bin[iAll][jAll]    ->  SetTitle(Form("Percentage variation of raw yield for bin %i",iAll+1));
            h2D_Syst_Bin[iAll][jAll]    ->  GetXaxis()  ->  SetTitle("Percentage variation");
        }
    }
    //------------//
    //  ANALYSIS  //
    //------------//
    
    // Output File for Fit Check
    TFile*  outFileFit  =   new TFile("./result/Syst_PID/Syst_PID_CheckRatioAndBins.root","recreate");
    
    //------ 1D Histograms ------//
    //
    h1D_Syst_Bin[0]         ->  Fill(0.,1.*nPIDFiles);
    for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )  {
        h1D_Syst_Bin[iPT1D+1]   ->  Fill(0);
    }
    for ( Int_t iTer = 1; iTer <= nPIDFiles; iTer++ )    {
        for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )  {
            auto    fStandard       =   g1D_Stat[0]         ->GetPointY     (iPT1D);
            auto    fStdError       =   g1D_Stat[0]         ->GetErrorYhigh  (iPT1D);
            auto    fVariatin       =   g1D_Stat[iTer]      ->GetPointY     (iPT1D);
            auto    fVarError       =   g1D_Stat[iTer]      ->GetErrorYhigh  (iPT1D);
            auto    fStdVarRt       =   fVariatin/fStandard;
            auto    fSVRError       =   (fVariatin/fStandard)*sqrt(fStdError*fStdError/(fStandard*fStandard)+fVariatin*fVariatin/(fVarError*fVarError));
            if ( fBarlowCheck(fStandard,fStdError,fVariatin,fVarError) )   continue;
            h1D_Syst_Bin[iPT1D+1]   ->  Fill(1-fStdVarRt);
            h1D_Syst_Bin[0]         ->  Fill(1-fStdVarRt);
            //h1D_Syst_Bin_StSy       ->  Fill(min(0.,fabs(h1D_Systematics_Util->GetBinContent(iPT1D)-1)-(h1D_Syst[0]->GetBinError(iPT1D))/(h1D_Syst[0]->GetBinContent(iPT1D))));
        }
    }
    //
    for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )  {
        auto    fStandard       =   g1D_Stat[0]         ->GetPointY     (iPT1D);
        auto    fStdError       =   g1D_Stat[0]         ->GetErrorYhigh (iPT1D);
        h1D_Stat_Bin            ->  SetBinContent(iPT1D+1,fStdError/fStandard);
    }
    //
    h1D_Syst_Bin_StSy->Write();
    for ( Int_t iPT1D = 0; iPT1D <= nBinPT1D; iPT1D++ )  {
        h1D_Syst_Bin[iPT1D] ->  Write();
    }
    //
    
    //------ 2D Histograms ------//
    //
    /*
    TH2F   *h2D_Systematics_Util    =   new TH2F("2D","2D",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    for ( Int_t iTer = 1; iTer < nOption2; iTer++ )    {
        h2D_Systematics_Util        ->  Divide(h2D_Syst[iTer],h2D_Syst[0]);
        h2D_Systematics_Util        ->  SetName(Form("hDivide_2D_%s",sOption2[iTer-1].c_str()));
        for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
            for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )  {
                auto    fStandard       =   h2D_Syst[0]         ->GetBinContent (iPT2D+1,jPT2D+1);
                auto    fStdError       =   h2D_Syst[0]         ->GetBinError   (iPT2D+1,jPT2D+1);
                auto    fVariatin       =   h2D_Syst[iTer]      ->GetBinContent (iPT2D+1,jPT2D+1);
                auto    fVarError       =   h2D_Syst[iTer]      ->GetBinError   (iPT2D+1,jPT2D+1);
                auto    fStdVarRt       =   h2D_Systematics_Util->GetBinContent (iPT2D+1,jPT2D+1);
                auto    fSVRError       =   h2D_Systematics_Util->GetBinError   (iPT2D+1,jPT2D+1);
                h2D_Syst_Bin[iPT2D][jPT2D]  ->  Fill(fStdVarRt);
                h2D_Syst_Bin_All            ->  Fill(fStdVarRt);
                h2D_Syst_Bin_StSy           ->  Fill(min(0.,fabs(h2D_Systematics_Util->GetBinContent(iPT2D+1,jPT2D+1)-1)-(h2D_Syst[0]->GetBinError(iPT2D+1,jPT2D+1))/(h2D_Syst[0]->GetBinContent(iPT2D+1,jPT2D+1))));
            }
        }
        h2D_Systematics_Util->Write();
    }
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
        for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )  {
            auto    fStandard       =   h2D_Syst[0]         ->GetBinContent (iPT2D+1,jPT2D+1);
            auto    fStdError       =   h2D_Syst[0]         ->GetBinError   (iPT2D+1,jPT2D+1);
            h2D_Stat_Bin                ->  SetBinContent(iPT2D+1,jPT2D+1,fStdError/fStandard);
        }
    }
    */
    h2D_Syst_Bin_StSy->Write();
    h2D_Syst_Bin_All->Write();
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
        for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )  {
            h2D_Syst_Bin[iPT2D][jPT2D] ->  Write();
        }
    }
    //
    // Output File for Fit Check
    TFile*  outFileFi2  =   new TFile("./result/Syst_PID/Syst_PIS_MeanAndRMS.root","recreate");
    //
    //--
    //  Individual Histograms
    //--
    //
    TH1F*   h1D_Syst_Mean   =   new TH1F    ("h1D_Syst_Mean","h1D_Syst_Mean",nBinPT1D,fArrPT1D);
    TH1F*   h1D_Syst_RMS_   =   new TH1F    ("h1D_Syst_RMS_","h1D_Syst_RMS_",nBinPT1D,fArrPT1D);
    TH1F*   h1D_Syst_Full   =   new TH1F    ("h1D_Syst_Full","h1D_Syst_Full",nBinPT1D,fArrPT1D);
    for ( Int_t iPT1D = 1; iPT1D <= nBinPT1D; iPT1D++ )  {
        h1D_Syst_Mean   ->SetBinContent(iPT1D,fabs(h1D_Syst_Bin[iPT1D] ->  GetMean() ));
        h1D_Syst_RMS_   ->SetBinContent(iPT1D,h1D_Syst_Bin[iPT1D] ->  GetRMS());
        h1D_Syst_Full   ->SetBinContent(iPT1D,h1D_Syst_Bin[iPT1D] ->  GetRMS() + fabs(h1D_Syst_Bin[iPT1D] ->  GetMean() ));
    }
    //
    TCanvas        *cDrawComparison     =   new TCanvas("cDrawComparison","");
    gStyle                              ->  SetOptStat(0);
    //
    h1D_Stat_Bin                        ->  SetLineWidth(3);
    h1D_Stat_Bin                        ->  SetLineColor(kBlue);
    h1D_Syst_Full                       ->  SetLineWidth(3);
    h1D_Syst_Full                       ->  SetLineColor(kRed);
    //
    TLegend        *cComparisonLegend   =   new TLegend(0.15,0.75,0.25,0.85);
    cComparisonLegend                   ->  SetLineColorAlpha(1,0.);
    cComparisonLegend                   ->  AddEntry(h1D_Stat_Bin,"Stat.","L");
    cComparisonLegend                   ->  AddEntry(h1D_Syst_Full,"Syst.","L");
    //
    h1D_Stat_Bin                        ->  Draw();
    h1D_Syst_Full                       ->  Draw("SAME");
    cComparisonLegend                   ->  Draw("SAME");
    cDrawComparison                     ->  Write();
    //
    TCanvas        *cDrawCompariso2     =   new TCanvas("cDrawCompariso2","");
    cDrawCompariso2                     ->  SetLogy();
    //
    TLegend        *cComparisonLegen2   =   new TLegend();
    TMultiGraph    *cComparisonGraphs   =   new TMultiGraph();
    for ( Int_t iFls = 0; iFls <= nPIDFiles; iFls++ )    {
        cComparisonGraphs               ->  Add(g1D_Stat[iFls],"LEP");
        cComparisonLegen2               ->  AddEntry(g1D_Stat[iFls]);
    }
    //
    cComparisonGraphs                   ->  Draw("ALP");
    cComparisonLegen2                   ->  Draw("SAME");
    cDrawCompariso2                     ->  Write();
    //
    h1D_Stat_Bin    ->  Write();
    h1D_Syst_Mean   ->  Write();
    h1D_Syst_RMS_   ->  Write();
    h1D_Syst_Full   ->  Write();
    //
    TH2F*   h2D_Syst_Mean   =   new TH2F    ("h2D_Syst_Mean","h2D_Syst_Mean",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    TH2F*   h2D_Syst_RMS_   =   new TH2F    ("h2D_Syst_RMS_","h2D_Syst_RMS_",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    TH2F*   h2D_Syst_Full   =   new TH2F    ("h2D_Syst_Full","h2D_Syst_Full",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
        for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )  {
            h2D_Syst_Mean   ->SetBinContent(iPT2D+1,jPT2D+1,fabs(h2D_Syst_Bin[iPT2D][jPT2D] ->  GetMean() ));
            h2D_Syst_RMS_   ->SetBinContent(iPT2D+1,jPT2D+1,h2D_Syst_Bin[iPT2D][jPT2D] ->  GetRMS());
            h2D_Syst_Full   ->SetBinContent(iPT2D+1,jPT2D+1,h2D_Syst_Bin[iPT2D][jPT2D] ->  GetRMS() + fabs(h2D_Syst_Bin[iPT2D][jPT2D] ->  GetMean() ));
        }
    }
    h2D_Stat_Bin    ->  Write();
    h2D_Syst_Mean   ->  Write();
    h2D_Syst_RMS_   ->  Write();
    h2D_Syst_Full   ->  Write();
    //
    // >-> Close input File
    //
    for ( Int_t iFls = 0; iFls <= nPIDFiles; iFls++ )    {
        insFileHXD[iFls]    ->  Close();
    }
    //
    outFileFit->Close();
    outFileFi2->Close();
}
