// File for final combinations of results:
// !TODO: N/A
#include "../../inc/AliAnalysisPhiPair.h"

void
FinalPlots
() {
    //
    //  Recovering calculated Systematic Errors
    TFile*  insFile_SY_Yield    =   new TFile   (Form("%s%s",Form(kAnalysis_Systemt_Dir,"Yield_7TeV"),"/Full_Systematics.root"));
    //
    TH1F*       hFullSyst;
    hName       =   "hCalculateFull_Sng";
    hFullSyst   =   (TH1F*)(insFile_SY_Yield->Get(hName));
    //
    //  Recovering calculated Systematic Results
    TFile*  insFile_DT_Yield    =   new TFile   (Form(kASigExtp_FitCheckRst,"Yield_7TeV"));
    TFile*  insFile_DT_Yiel2    =   new TFile   (Form(kASigExtp_FitCheckRst,"Yield"));
    TFile*  insFile_DT_Mult_    =   new TFile   (Form(kASigExtp_FitCheckRst,"Multiplicity"));
    //
    TH1F   *hMultipleResults_Stat,*hMultipleResults_Syst,*hMultipleResult2_Stat,*hMultipleResult2_Syst;
    hName       =   "hMultipleResults_Stat";
    hMultipleResults_Stat   =   (TH1F*)(insFile_DT_Yield->Get(hName));
    //
    hName       =   "hMultipleResults_Syst";
    hMultipleResults_Syst   =   (TH1F*)(insFile_DT_Yield->Get(hName));
    //
    hName       =   "hMultipleResults_Stat";
    hMultipleResult2_Stat   =   (TH1F*)(insFile_DT_Yiel2->Get(hName));
    //
    hName       =   "hMultipleResults_Syst";
    hMultipleResult2_Syst   =   (TH1F*)(insFile_DT_Yiel2->Get(hName));
    //
    std::vector<TString> kShowWhat = {"Y_1","Y_2","R_1","R_2","P_1","P_2"};
    for ( Int_t iTer = 0; iTer < 6; iTer++ )    {
        hMultipleResults_Syst->SetBinError  ( iTer+1, hMultipleResults_Stat->GetBinContent(iTer+1)*hFullSyst->GetBinContent(iTer+1) );
        
        auto k1Yield    =   hMultipleResults_Stat->GetBinContent(1);
        auto k2Yield    =   hMultipleResults_Stat->GetBinContent(2);
        auto kHigNorm   =   hMultipleResults_Stat->GetBinContent(iTer+1)*kSysHig_TE;
        auto kLowNorm   =   hMultipleResults_Stat->GetBinContent(iTer+1)*kSysLow_TE;
        if ( iTer +1 == 3 ) {
            kHigNorm    =   0;
            kLowNorm    =   0;
        }
        if ( iTer +1 == 4 ) {
            kHigNorm    =   k2Yield/(k1Yield*(1+kSysLow_TE));
            kLowNorm    =   k2Yield/(k1Yield*(1+kSysHig_TE));
        }
        if ( iTer +1 == 5 ) {
            kHigNorm = fabs( fSigmaPhiValue(k1Yield*(1+kSysHig_TE),k2Yield*(1+kSysHig_TE)) - hMultipleResults_Stat->GetBinContent(5) );
            kLowNorm = fabs( fSigmaPhiValue(k1Yield*(1-kSysLow_TE),k2Yield*(1-kSysLow_TE)) - hMultipleResults_Stat->GetBinContent(5) );
        }
        if ( iTer +1 == 6 ) {
            kHigNorm = fabs( fGammaPhiValue(k1Yield*(1-kSysLow_TE),k2Yield*(1-kSysLow_TE)) - hMultipleResults_Stat->GetBinContent(6) );
            kLowNorm = fabs( fGammaPhiValue(k1Yield*(1+kSysHig_TE),k2Yield*(1+kSysHig_TE)) - hMultipleResults_Stat->GetBinContent(6) );
        }
        cout << "-------------------------------" << endl;
        cout << "-       " << "       +" << hMultipleResults_Stat->GetBinError(iTer+1) << " +" << hMultipleResults_Stat->GetBinContent(iTer+1)*hFullSyst->GetBinContent(iTer+1) << " +" << kHigNorm << endl;
        cout << "- " << kShowWhat.at(iTer).Data() << " : " << hMultipleResults_Stat->GetBinContent(iTer+1) << endl;
        cout << "-       " << "       -" << hMultipleResults_Stat->GetBinError(iTer+1) << " -" << hMultipleResults_Stat->GetBinContent(iTer+1)*hFullSyst->GetBinContent(iTer+1) << " -" << kLowNorm << endl;
        cout << "-------------------------------" << endl;
    }
    //
    TCanvas*    cDrawResults    =   new TCanvas();
    //
    //  7TeV
    //  STAT
    hMultipleResults_Stat->SetMarkerStyle(fGetMarker(3));
    hMultipleResults_Stat->SetMarkerColor(fGetColor(4));
    hMultipleResults_Stat->SetFillColor(fGetColor(4));
    hMultipleResults_Stat->SetLineColorAlpha(fGetColor(3),1.);
    //  SYST
    hMultipleResults_Syst->SetMarkerStyle(fGetMarker(3));
    hMultipleResults_Syst->SetMarkerColor(fGetColor(5));
    hMultipleResults_Syst->SetFillColorAlpha(fGetColor(5),0.3);
    hMultipleResults_Syst->SetLineColorAlpha(0.,0.);
    hMultipleResults_Stat->Draw("E1P MIN0");
    hMultipleResults_Syst->Draw("E2 SAME MIN0");
    //
    //  5TeV
    //  STAT
    hMultipleResult2_Stat->SetMarkerStyle(fGetMarker(3));
    hMultipleResult2_Stat->SetMarkerColor(fGetColor(3));
    hMultipleResult2_Stat->SetFillColor(fGetColor(3));
    hMultipleResult2_Stat->SetLineColorAlpha(fGetColor(3),1.);
    hMultipleResult2_Stat->Draw("E1P SAME MIN0");
    //
    //  TLegend();
    TLegend*    lLegend =   new TLegend(0.7,0.6,0.88,0.85);
    lLegend->AddEntry(hMultipleResults_Stat,"7TeV Stat","EP");
    lLegend->AddEntry(hMultipleResults_Syst,"7TeV Syst","F");
    lLegend->AddEntry(hMultipleResult2_Stat,"5TeV Stat","EP");
    lLegend->Draw("SAME");
    //
    cDrawResults->SaveAs("TEST.pdf");
    for ( Int_t iTer = 0; iTer < 6; iTer++ ) {
        hMultipleResults_Syst->GetXaxis()->SetRange(iTer+1,iTer+1);
        hMultipleResults_Syst->SetMaximum(hMultipleResults_Syst->GetBinContent(iTer+1) + 3*hMultipleResults_Syst->GetBinError(iTer+1));
        hMultipleResults_Syst->SetMinimum(hMultipleResults_Syst->GetBinContent(iTer+1) - 3*hMultipleResults_Syst->GetBinError(iTer+1));
        hMultipleResults_Syst->Draw("E2 P");
        hMultipleResults_Stat->Draw("E1 SAME");
        hMultipleResult2_Stat->Draw("E1P SAME");
        if ( iTer == 0 ) uLatex->DrawLatexNDC(0.18,0.80,"#frac{dN_{#phi}}{dy}");
        if ( iTer == 1 ) uLatex->DrawLatexNDC(0.18,0.80,"#frac{dN_{#phi#phi}}{dy}");
        if ( iTer == 2 ) uLatex->DrawLatexNDC(0.18,0.80,"#frac{#LT Y_{#phi#phi} #GT}{#LT Y_{#phi} #GT}");
        if ( iTer == 3 ) uLatex->DrawLatexNDC(0.18,0.80,"#frac{#LT Y_{#phi#phi} #GT}{#LT Y_{#phi} #GT^{2}}");
        if ( iTer == 4 ) uLatex->DrawLatexNDC(0.18,0.80,"#sigma^{2}_{#phi}");
        if ( iTer == 5 ) uLatex->DrawLatexNDC(0.18,0.80,"#gamma_{#phi}");
        lLegend->Draw("SAME");
        cDrawResults->SaveAs(Form("TEST_%i.pdf",iTer));
    }
    for ( Int_t iTer = 0; iTer < 6; iTer++ ) {
        TH1F * hCurrent_Plot;
        if ( iTer == 0 ) hCurrent_Plot = (TH1F*)(insFile_DT_Mult_->Get("hShow1D"));
        if ( iTer == 1 ) hCurrent_Plot = (TH1F*)(insFile_DT_Mult_->Get("hShow2D"));
        if ( iTer == 2 ) hCurrent_Plot = (TH1F*)(insFile_DT_Mult_->Get("hShowR1"));
        if ( iTer == 3 ) hCurrent_Plot = (TH1F*)(insFile_DT_Mult_->Get("hShowR2"));
        if ( iTer == 4 ) hCurrent_Plot = (TH1F*)(insFile_DT_Mult_->Get("hShowP1"));
        if ( iTer == 5 ) hCurrent_Plot = (TH1F*)(insFile_DT_Mult_->Get("hShowP2"));
        
        hCurrent_Plot->SetMinimum(0);
        hCurrent_Plot->Draw();
        
        if ( iTer == 0 ) uLatex->DrawLatexNDC(0.18,0.80,"#frac{dN_{#phi}}{dy}");
        if ( iTer == 1 ) uLatex->DrawLatexNDC(0.18,0.80,"#frac{dN_{#phi#phi}}{dy}");
        if ( iTer == 2 ) uLatex->DrawLatexNDC(0.18,0.80,"#frac{#LT Y_{#phi#phi} #GT}{#LT Y_{#phi} #GT}");
        if ( iTer == 3 ) uLatex->DrawLatexNDC(0.18,0.80,"#frac{#LT Y_{#phi#phi} #GT}{#LT Y_{#phi} #GT^{2}}");
        if ( iTer == 4 ) uLatex->DrawLatexNDC(0.18,0.80,"#sigma^{2}_{#phi}");
        if ( iTer == 5 ) uLatex->DrawLatexNDC(0.18,0.80,"#gamma_{#phi}");
        
        cDrawResults->SaveAs(Form("TEST_show_%i.pdf",iTer));
    }
    delete  cDrawResults;
}
