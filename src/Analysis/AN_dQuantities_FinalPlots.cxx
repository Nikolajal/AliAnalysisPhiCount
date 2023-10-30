#include "../../inc/AliAnalysisPhiPair.h"

template <typename TH1X_Target = TH1F>
TCanvas *
uPlotPtSpectra2D(std::vector<TH1X_Target *> hTargetList, std::vector<TH1X_Target *> hTargetListSystematics, Float_t kMultiplicationFactor = 15, Float_t kEnergy = 5.02, Bool_t kLogScale = true, TString tSavePlot = "")
{
    //! Silencing TCanvas pop-up
    gROOT->SetBatch(kTRUE);
    TCanvas *cPlotPtSpectra = new TCanvas("", "", 1200, 1000);
    gStyle->SetOptStat(0);
    gPad->SetTopMargin(0.04);
    gPad->SetBottomMargin(0.14);
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.04);
    if (kLogScale)
        gPad->SetLogy();
    //!
    TLegend *lPlotLegend = new TLegend(0.19, 0.82, 0.58, 0.90);
    lPlotLegend->SetNColumns(max(1, (int)(hTargetList.size() / 3. + 1.)));
    lPlotLegend->SetFillColorAlpha(kWhite, 0.);
    lPlotLegend->SetLineColorAlpha(kWhite, 0.);
    lPlotLegend->SetTextSize(.02);
    //!
    //! --- Range of plot
    auto iTer = -1;
    Double_t kMaximum = 0;
    Double_t kMinimum = 1e9;
    for (auto kCurrentTarget : hTargetList)
    {
        iTer++;
        auto kCurrentMaximum = kLogScale ? kCurrentTarget->GetMaximum() * pow(kMultiplicationFactor, iTer) : kCurrentTarget->GetMaximum() + kMultiplicationFactor * iTer;
        kMaximum = kMaximum > kCurrentMaximum ? kMaximum : kCurrentMaximum;
        kMinimum = kMinimum < kCurrentTarget->GetMinimum() ? kMinimum : kCurrentTarget->GetMinimum();
    }
    if (kLogScale)
        kMaximum = 1.e+8;
    else
        kMaximum = 2.5;
    if (kLogScale)
        kMinimum = 1.e-6;
    else
        kMinimum = 0.5;
    //!
    iTer = kLogScale ? -2 : -1;
    for (auto kCurrentTarget : hTargetList)
    {
        iTer++;
        if (iTer == -1)
        {
            continue;
        }
        auto kCurrentHistogram = (TH1X_Target *)kCurrentTarget->Clone(Form("tmp_%i", iTer));
        TH1X_Target *kCurrentSystemErr;
        if (kLogScale)
            kCurrentSystemErr = (TH1X_Target *)hTargetListSystematics.at(iTer + 1)->Clone(Form("tmpsys_%i", iTer));
        else
            kCurrentSystemErr = (TH1X_Target *)hTargetListSystematics.at(iTer)->Clone(Form("tmpsys_%i", iTer));
        //! Scale
        if (kLogScale)
        {
            kCurrentHistogram->Scale(pow(kMultiplicationFactor, iTer));
        }
        else
        {
            kCurrentHistogram = uShift(kCurrentHistogram, kMultiplicationFactor * iTer);
        }
        kCurrentHistogram->SetMaximum(kMaximum);
        kCurrentHistogram->SetMinimum(kMinimum);
        if (kLogScale)
        {
            kCurrentSystemErr->Scale(pow(kMultiplicationFactor, iTer));
        }
        else
        {
            kCurrentSystemErr = uShift(kCurrentSystemErr, kMultiplicationFactor * iTer);
        }
        kCurrentSystemErr->SetMaximum(kMaximum);
        kCurrentSystemErr->SetMinimum(kMinimum);
        //! Style
        kCurrentSystemErr->SetTitle("");
        kCurrentHistogram->SetLineColor(kRainbowColor[iTer]);
        kCurrentHistogram->SetMarkerColor(kRainbowColor[iTer]);
        kCurrentHistogram->SetMarkerStyle(kRainbowMarker[iTer]);
        kCurrentHistogram->SetMarkerSize(2);
        kCurrentSystemErr->SetLineColor(kRainbowColor[iTer]);
        kCurrentSystemErr->SetMarkerColor(kRainbowColor[iTer]);
        kCurrentSystemErr->SetMarkerStyle(kRainbowMarker[iTer]);
        kCurrentSystemErr->SetFillColor(kRainbowColor[iTer] - 8);
        kCurrentSystemErr->SetMarkerSize(2);
        //! Axes
        kCurrentHistogram->GetYaxis()->SetNdivisions(8);
        //! Plot
        lPlotLegend->AddEntry(kCurrentHistogram, Form("%s", kCurrentTarget->GetTitle()));
        kCurrentSystemErr->Draw("SAME PE2");
        kCurrentHistogram->Draw("SAME PE X0");
    }
    uLatex->SetTextFont(60);
    uLatex->SetTextSize(0.05);
    uLatex->DrawLatexNDC(0.65, 0.88, "ALICE");
    uLatex->SetTextFont(42);
    uLatex->SetTextSize(0.04);
    uLatex->DrawLatexNDC(0.65, 0.83, Form("pp #sqrt{#it{s}}= %.2f TeV", kEnergy));
    uLatex->DrawLatexNDC(0.65, 0.78, "#phi #rightarrow K^{+}K^{-}, |#it{y} | < 0.5");
    uLatex->SetTextSize(0.035);
    uLatex->DrawLatexNDC(0.19, 0.91, Form("Conditional p_{T} class"));
    uLatex->DrawLatexNDC(0.19, 0.22, Form("Uncertainties:"));
    uLatex->DrawLatexNDC(0.19, 0.17, Form("Stat. (bars), Syst. (boxes)"));
    lPlotLegend->Draw("same");
    //!
    if (!tSavePlot.IsNull())
        cPlotPtSpectra->SaveAs(tSavePlot + TString(".pdf"));
    if (!tSavePlot.IsNull())
        cPlotPtSpectra->SaveAs(tSavePlot + TString(".eps"));
    //!
    gROOT->SetBatch(kFALSE);
    return cPlotPtSpectra;
}

template <typename TH1X_Target = TH1F>
TCanvas *
uPlotPtSpectra2DRatio(std::vector<TH1X_Target *> hTargetList, std::vector<TH1X_Target *> hTargetListSystematics, Int_t kMultiplicationFactor = 15, Float_t kEnergy = 5.02, Bool_t kLogScale = true, TString tSavePlot = "")
{
    //! Silencing TCanvas pop-up
    gROOT->SetBatch(kFALSE);
    TCanvas *cPlotPtSpectra = new TCanvas("", "", 1200, 1000);
    gStyle->SetOptStat(0);
    gPad->SetTopMargin(0.04);
    gPad->SetBottomMargin(0.14);
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.04);
    if (kLogScale)
        gPad->SetLogy();
    //!
    TLegend *lPlotLegend = new TLegend(0.19, 0.17, 0.64, 0.25);
    lPlotLegend->SetNColumns(max(1, (int)(hTargetList.size() / 3)));
    lPlotLegend->SetFillColorAlpha(kWhite, 0.);
    lPlotLegend->SetLineColorAlpha(kWhite, 0.);
    lPlotLegend->SetTextSize(.02);
    //!
    //! --- Range of plot
    Int_t iTer = -2;
    Double_t kMaximum = 0;
    Double_t kMinimum = 1e9;
    for (auto kCurrentTarget : hTargetList)
    {
        iTer++;
        if (iTer == -1)
        {
            continue;
        }
        auto kCurrentMaximum = kLogScale ? kCurrentTarget->GetMaximum() * pow(kMultiplicationFactor, iTer) : kCurrentTarget->GetMaximum() + kMultiplicationFactor * iTer;
        kMaximum = kMaximum > kCurrentMaximum ? kMaximum : kCurrentMaximum;
        kMinimum = kMinimum < kCurrentTarget->GetMinimum() ? kMinimum : kCurrentTarget->GetMinimum();
    }
    if (kLogScale)
        kMaximum *= pow(kMultiplicationFactor, 3);
    else
        kMaximum += kMultiplicationFactor * 2;
    if (kLogScale)
        kMinimum /= pow(kMultiplicationFactor, 4);
    else
        kMinimum -= kMultiplicationFactor * 2;
    //!
    iTer = -2;
    for (auto kCurrentTarget : hTargetList)
    {
        iTer++;
        if (iTer == -1)
        {
            continue;
        }
        auto kCurrentHistogram = (TH1X_Target *)kCurrentTarget->Clone(Form("tmp_%i", iTer));
        auto kCurrentSystemErr = (TH1X_Target *)hTargetListSystematics.at(iTer + 1)->Clone(Form("tmpsys_%i", iTer));
        //! Scale
        if (kLogScale)
        {
            kCurrentHistogram->Scale(pow(kMultiplicationFactor, iTer));
        }
        else
        {
            kCurrentHistogram = uShift(kCurrentHistogram, kMultiplicationFactor * iTer);
        }
        kCurrentHistogram->SetMaximum(kMaximum);
        kCurrentHistogram->SetMinimum(kMinimum);
        if (kLogScale)
        {
            kCurrentSystemErr->Scale(pow(kMultiplicationFactor, iTer));
        }
        else
        {
            kCurrentSystemErr = uShift(kCurrentSystemErr, kMultiplicationFactor * iTer);
        }
        kCurrentSystemErr->SetMaximum(kMaximum);
        kCurrentSystemErr->SetMinimum(kMinimum);
        //! Style
        kCurrentSystemErr->SetTitle("");
        kCurrentHistogram->SetLineColor(kRainbowColor[iTer]);
        kCurrentHistogram->SetMarkerColor(kRainbowColor[iTer]);
        kCurrentHistogram->SetMarkerStyle(kRainbowMarker[iTer]);
        kCurrentHistogram->SetMarkerSize(2);
        kCurrentSystemErr->SetLineColor(kRainbowColor[iTer]);
        kCurrentSystemErr->SetMarkerColor(kRainbowColor[iTer]);
        kCurrentSystemErr->SetMarkerStyle(kRainbowMarker[iTer]);
        kCurrentSystemErr->SetFillColor(kRainbowColor[iTer] - 8);
        kCurrentSystemErr->SetMarkerSize(2);
        //! Axes
        kCurrentHistogram->GetYaxis()->SetNdivisions(8);
        //! Plot
        lPlotLegend->AddEntry(kCurrentHistogram, Form("%s", kCurrentTarget->GetTitle()));
        kCurrentSystemErr->Draw("SAME PE2");
        kCurrentHistogram->Draw("SAME PE X0");
    }
    uLatex->SetTextFont(60);
    uLatex->SetTextSize(0.05);
    uLatex->DrawLatexNDC(0.58, 0.90, "ALICE");
    uLatex->SetTextFont(42);
    uLatex->SetTextSize(0.04);
    uLatex->DrawLatexNDC(0.58, 0.85, Form("pp #sqrt{#it{s}}= %.2f TeV", kEnergy));
    uLatex->DrawLatexNDC(0.58, 0.80, "#phi #rightarrow K^{+}K^{-}, |#it{y} | < 0.5");
    uLatex->SetTextSize(0.035);
    uLatex->DrawLatexNDC(0.19, 0.27, Form("Multiplicity class"));
    uLatex->SetTextSize(0.035);
    uLatex->DrawLatexNDC(0.62, 0.23, Form("Uncertainties:"));
    uLatex->DrawLatexNDC(0.62, 0.18, Form("Stat. (bars), Syst. (boxes)"));
    lPlotLegend->Draw("same");
    //!
    if (!tSavePlot.IsNull())
        cPlotPtSpectra->SaveAs(tSavePlot + TString(".pdf"));
    if (!tSavePlot.IsNull())
        cPlotPtSpectra->SaveAs(tSavePlot + TString(".eps"));
    //!
    gROOT->SetBatch(kFALSE);
    return cPlotPtSpectra;
}

std::vector<TString> kClassesLabel = {"I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "XVII", "XVIII", "IX", "XX"};
std::vector<TString> kMultClassesLabel = {"I+II", "III+IV", "V+VI", "VII+VIII", "IX+X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "XVII", "XVIII", "IX", "XX"};

void AN_dQuantities_FinalPlots()
{
    //!
    //! SET-UP
    fSetAllBins();
    SetStyle();
    //!
    //! Load all needed histograms
    TString kFolder = "";
    kFolder = "_p_p__7TeV";
    TFile *insFile_Data_YL_p_p__07 = new TFile(Form(kASigExtp_FitCheckRst, (TString("Yield") + kFolder).Data()));
    TFile *insFile_Syst_YL_p_p__07 = new TFile(Form("%s/FullSystematics.root", Form(kAnalysis_Systemt_Dir, (TString("Yield") + kFolder).Data())));
    kFolder = "_p_p__5TeV";
    TFile *insFile_Data_YL_p_p__05 = new TFile(Form(kASigExtp_FitCheckRst, (TString("Yield") + kFolder).Data()));
    TFile *insFile_Data_MT_p_p__05 = new TFile(Form(kASigExtp_FitCheckRst, (TString("Multiplicity") + kFolder).Data()));
    TFile *insFile_Syst_YL_p_p__05 = new TFile(Form("%s/FullSystematics.root", Form(kAnalysis_Systemt_Dir, (TString("Yield") + kFolder).Data())));
    TFile *insFile_Syst_MT_p_p__05 = new TFile(Form("%s/FullSystematics.root", Form(kAnalysis_Systemt_Dir, (TString("Multiplicity") + kFolder).Data())));
    //!
    auto hXD_Nyld_stat_pp_05 = uLoadHistograms<0, TH1F>(insFile_Data_YL_p_p__05, "hXD_Nfqs_stat", "hXD_Nfqs_stat_p_p__05");
    auto hXD_Nyld_stat_pp_07 = uLoadHistograms<0, TH1F>(insFile_Data_YL_p_p__07, "hXD_Nfqs_stat", "hXD_Nfqs_stat_p_p__07");
    auto h2D_MeanPT_stat_pp_05 = uLoadHistograms<0, TH1F>(insFile_Data_YL_p_p__05, "h2D_MeanPT_stat", "h2D_MeanPT_stat_pp_05");
    auto h2D_MeanPT_syst_pp_05 = uLoadHistograms<0, TH1F>(insFile_Data_YL_p_p__05, "h2D_MeanPT_syst", "h2D_MeanPT_syst_pp_05");
    auto h2D_MeanPT_stat_pp_07 = uLoadHistograms<0, TH1F>(insFile_Data_YL_p_p__07, "h2D_MeanPT_stat", "h2D_MeanPT_stat_pp_07");
    auto h2D_MeanPT_syst_pp_07 = uLoadHistograms<0, TH1F>(insFile_Data_YL_p_p__07, "h2D_MeanPT_syst", "h2D_MeanPT_syst_pp_07");
    auto h2D_Nres_stat_pp_05 = uLoadHistograms<1, TH1F>(insFile_Data_YL_p_p__05, "h2D_Nres_stat_stat_%i", "h2D_Nres_stat_stat_%i_p_p__05");
    auto h2D_Nres_stat_pp_07 = uLoadHistograms<1, TH1F>(insFile_Data_YL_p_p__07, "h2D_Nres_stat_stat_%i", "h2D_Nres_stat_stat_%i_p_p__07");
    auto h2D_Nres_syst_pp_05 = uLoadHistograms<1, TH1F>(insFile_Data_YL_p_p__05, "h2D_Nres_syst_syst_%i", "h2D_Nres_syst_syst_%i_p_p__05");
    auto h2D_Nres_syst_pp_07 = uLoadHistograms<1, TH1F>(insFile_Data_YL_p_p__07, "h2D_Nres_syst_syst_%i", "h2D_Nres_syst_syst_%i_p_p__07");
    auto hXD_Efrc_syst_pp_05 = uLoadHistograms<0, TH1F>(insFile_Syst_YL_p_p__05, "hFullSystematicsRT", "hFullSystematicsRT_p_p__07");
    auto hXD_Efrc_syst_pp_07 = uLoadHistograms<0, TH1F>(insFile_Syst_YL_p_p__07, "hFullSystematicsRT", "hFullSystematicsRT_p_p__05");
    auto g1D_Nres_stat_Mult_pp_05 = uLoadHistograms<0, TGraphErrors>(insFile_Data_MT_p_p__05, "g1D_Nres_Stat_Mult", "g1D_Nres_stat_Mult");
    auto g2D_Nres_stat_Mult_pp_05 = uLoadHistograms<0, TGraphErrors>(insFile_Data_MT_p_p__05, "g2D_Nres_Stat_Mult", "g2D_Nres_stat_Mult");
    auto gR2_Nres_stat_Mult_pp_05 = uLoadHistograms<0, TGraphErrors>(insFile_Data_MT_p_p__05, "gR2_Nres_Stat_Mult", "gR2_Nres_stat_Mult");
    auto gP1_Nres_stat_Mult_pp_05 = uLoadHistograms<0, TGraphErrors>(insFile_Data_MT_p_p__05, "gP1_Nres_Stat_Mult", "gP1_Nres_stat_Mult");
    auto gP2_Nres_stat_Mult_pp_05 = uLoadHistograms<0, TGraphErrors>(insFile_Data_MT_p_p__05, "gP2_Nres_Stat_Mult", "gP2_Nres_stat_Mult");
    auto g1D_Nres_syst_Mult_pp_05 = uLoadHistograms<0, TGraphErrors>(insFile_Data_MT_p_p__05, "g1D_Nres_Syst_Mult", "g1D_Nres_syst_Mult");
    auto g2D_Nres_syst_Mult_pp_05 = uLoadHistograms<0, TGraphErrors>(insFile_Data_MT_p_p__05, "g2D_Nres_Syst_Mult", "g2D_Nres_syst_Mult");
    auto gR2_Nres_syst_Mult_pp_05 = uLoadHistograms<0, TGraphErrors>(insFile_Data_MT_p_p__05, "gR2_Nres_Syst_Mult", "gR2_Nres_syst_Mult");
    auto gP2_Nres_syst_Mult_pp_05 = uLoadHistograms<0, TGraphErrors>(insFile_Data_MT_p_p__05, "gP2_Nres_Syst_Mult", "gP2_Nres_syst_Mult");
    auto h2D_MeanPT_stat_MT_pp_05 = uLoadHistograms<1, TH1F>(insFile_Data_MT_p_p__05, "h2D_MeanPT_stat_MT_%i", "h2D_MeanPT_stat_pp_05_MT_%i");
    auto h2D_MeanPT_syst_MT_pp_05 = uLoadHistograms<1, TH1F>(insFile_Data_MT_p_p__05, "h2D_MeanPT_syst_MT_%i", "h2D_MeanPT_syst_pp_05_MT_%i");
    //!
    //! Load all MC histograms
    TFile *fComparison01 = new TFile("Pythia8Monash.root");
    auto gP800_N1 = uLoadHistograms<0, TGraphErrors>(fComparison01, "g1DYield", "gP800_N1");
    auto gP800_N2 = uLoadHistograms<0, TGraphErrors>(fComparison01, "g2DYield", "gP800_N2");
    auto gP800_R1 = uLoadHistograms<0, TGraphErrors>(fComparison01, "gR1Yield", "gP800_R1");
    auto gP800_R2 = uLoadHistograms<0, TGraphErrors>(fComparison01, "gR2Yield", "gP800_R2");
    auto gP800_P1 = uLoadHistograms<0, TGraphErrors>(fComparison01, "gP1Yield", "gP800_P1");
    auto gP800_P2 = uLoadHistograms<0, TGraphErrors>(fComparison01, "gP2Yield", "gP800_P2");
    gP800_N1->SetLineColor(kBlue);
    gP800_N1->SetLineWidth(3);
    gP800_N2->SetLineColor(kBlue);
    gP800_N2->SetLineWidth(3);
    gP800_R1->SetLineColor(kBlue);
    gP800_R1->SetLineWidth(3);
    gP800_R2->SetLineColor(kBlue);
    gP800_R2->SetLineWidth(3);
    gP800_P1->SetLineColor(kBlue);
    gP800_P1->SetLineWidth(3);
    gP800_P2->SetLineColor(kBlue);
    gP800_P2->SetLineWidth(3);
    TFile *fComparison02 = new TFile("Pythia8MonashRopes.root");
    auto gP808_N1 = uLoadHistograms<0, TGraphErrors>(fComparison02, "g1DYield", "gP800_N1");
    auto gP808_N2 = uLoadHistograms<0, TGraphErrors>(fComparison02, "g2DYield", "gP800_N2");
    auto gP808_R1 = uLoadHistograms<0, TGraphErrors>(fComparison02, "gR1Yield", "gP800_R1");
    auto gP808_R2 = uLoadHistograms<0, TGraphErrors>(fComparison02, "gR2Yield", "gP800_R2");
    auto gP808_P1 = uLoadHistograms<0, TGraphErrors>(fComparison02, "gP1Yield", "gP800_P1");
    auto gP808_P2 = uLoadHistograms<0, TGraphErrors>(fComparison02, "gP2Yield", "gP800_P2");
    gP808_N1->SetLineColor(kBlue);
    gP808_N1->SetLineWidth(3);
    gP808_N1->SetLineStyle(3);
    gP808_N2->SetLineColor(kBlue);
    gP808_N2->SetLineWidth(3);
    gP808_N2->SetLineStyle(3);
    gP808_R1->SetLineColor(kBlue);
    gP808_R1->SetLineWidth(3);
    gP808_R1->SetLineStyle(3);
    gP808_R2->SetLineColor(kBlue);
    gP808_R2->SetLineWidth(3);
    gP808_R2->SetLineStyle(3);
    gP808_P1->SetLineColor(kBlue);
    gP808_P1->SetLineWidth(3);
    gP808_P1->SetLineStyle(3);
    gP808_P2->SetLineColor(kBlue);
    gP808_P2->SetLineWidth(3);
    gP808_P2->SetLineStyle(3);
    //!
    std::vector<TH1F *> h2D_MeanPT_stat;
    h2D_MeanPT_stat.push_back(h2D_MeanPT_stat_pp_05);
    h2D_MeanPT_stat.push_back(h2D_MeanPT_stat_pp_07);
    std::vector<TH1F *> h2D_MeanPT_syst;
    h2D_MeanPT_syst.push_back(h2D_MeanPT_syst_pp_05);
    h2D_MeanPT_syst.push_back(h2D_MeanPT_syst_pp_07);
    //!
    auto iTer = -1;
    for (auto kCurrentTarget : h2D_Nres_stat_pp_05)
    {
        iTer++;
        uSetHisto(kCurrentTarget, "SPT 12D");
        kCurrentTarget->SetTitle(Form("%s (x%i^{%i})", kClassesLabel.at(iTer).Data(), 30, iTer));
    }
    iTer = -1;
    for (auto kCurrentTarget : h2D_Nres_syst_pp_05)
    {
        iTer++;
        uSetHisto(kCurrentTarget, "SPT 12D");
    }
    iTer = -1;
    for (auto kCurrentTarget : h2D_Nres_stat_pp_07)
    {
        iTer++;
        uSetHisto(kCurrentTarget, "SPT 12D");
        kCurrentTarget->SetTitle(Form("%s (x%i^{%i})", kClassesLabel.at(iTer).Data(), 30, iTer));
    }
    iTer = -1;
    for (auto kCurrentTarget : h2D_Nres_syst_pp_07)
    {
        iTer++;
        uSetHisto(kCurrentTarget, "SPT 12D");
    }
    iTer = -1;
    for (auto kCurrentTarget : h2D_MeanPT_stat)
    {
        iTer++;
        uSetHisto(kCurrentTarget, "SPT 12D");
    }
    h2D_MeanPT_stat.at(0)->SetTitle("5 TeV");
    h2D_MeanPT_stat.at(1)->SetTitle("7 TeV");
    iTer = -1;
    for (auto kCurrentTarget : h2D_MeanPT_syst)
    {
        iTer++;
        uSetHisto(kCurrentTarget, "MPT 12D");
    }
    iTer = -2;
    for (auto kCurrentTarget : h2D_MeanPT_stat_MT_pp_05)
    {
        iTer++;
        if (iTer == -1)
        {
            continue;
        }
        uSetHisto(kCurrentTarget, "MPT 12D");
        kCurrentTarget->SetTitle(Form("%s (+%i*%i)", kMultClassesLabel.at(iTer).Data(), 2, iTer));
    }
    iTer = -1;
    for (auto kCurrentTarget : h2D_MeanPT_syst_MT_pp_05)
    {
        iTer++;
        uSetHisto(kCurrentTarget, "MPT 12D");
    }
    //!
    TGraphErrors *gInclusive_stat_pp_07_Y1 = new TGraphErrors();
    gInclusive_stat_pp_07_Y1->SetPoint(0, 4.60, hXD_Nyld_stat_pp_07->GetBinContent(1));
    gInclusive_stat_pp_07_Y1->SetPointError(0, 0.1, hXD_Nyld_stat_pp_07->GetBinError(1));
    TGraphErrors *gInclusive_syst_pp_07_Y1 = new TGraphErrors();
    gInclusive_syst_pp_07_Y1->SetPoint(0, 4.60, hXD_Nyld_stat_pp_07->GetBinContent(1));
    gInclusive_syst_pp_07_Y1->SetPointError(0, 0.1, hXD_Nyld_stat_pp_07->GetBinContent(1) * hXD_Efrc_syst_pp_07->GetBinContent(1));
    TGraphErrors *gInclusive_stat_pp_07_Y2 = new TGraphErrors();
    gInclusive_stat_pp_07_Y2->SetPoint(0, 4.60, 10 * hXD_Nyld_stat_pp_07->GetBinContent(2));
    gInclusive_stat_pp_07_Y2->SetPointError(0, 0.1, 10 * hXD_Nyld_stat_pp_07->GetBinError(2));
    TGraphErrors *gInclusive_syst_pp_07_Y2 = new TGraphErrors();
    gInclusive_syst_pp_07_Y2->SetPoint(0, 4.60, 10 * hXD_Nyld_stat_pp_07->GetBinContent(2));
    gInclusive_syst_pp_07_Y2->SetPointError(0, 0.1, 10 * hXD_Nyld_stat_pp_07->GetBinContent(2) * hXD_Efrc_syst_pp_07->GetBinContent(2));
    //!
    TGraphErrors *gInclusive_stat_pp_07_R2 = new TGraphErrors();
    gInclusive_stat_pp_07_R2->SetPoint(0, 4.60, hXD_Nyld_stat_pp_07->GetBinContent(4));
    gInclusive_stat_pp_07_R2->SetPointError(0, 0.1, hXD_Nyld_stat_pp_07->GetBinError(4));
    TGraphErrors *gInclusive_syst_pp_07_R2 = new TGraphErrors();
    gInclusive_syst_pp_07_R2->SetPoint(0, 4.60, hXD_Nyld_stat_pp_07->GetBinContent(4));
    gInclusive_syst_pp_07_R2->SetPointError(0, 0.1, hXD_Nyld_stat_pp_07->GetBinContent(4) * hXD_Efrc_syst_pp_07->GetBinContent(4));
    TGraphErrors *gInclusive_stat_pp_07_P2 = new TGraphErrors();
    gInclusive_stat_pp_07_P2->SetPoint(0, 4.60, 10 * hXD_Nyld_stat_pp_07->GetBinContent(6));
    gInclusive_stat_pp_07_P2->SetPointError(0, 0.1, 10 * hXD_Nyld_stat_pp_07->GetBinError(6));
    TGraphErrors *gInclusive_syst_pp_07_P2 = new TGraphErrors();
    gInclusive_syst_pp_07_P2->SetPoint(0, 4.60, 10 * hXD_Nyld_stat_pp_07->GetBinContent(6));
    gInclusive_syst_pp_07_P2->SetPointError(0, 0.1, 10 * hXD_Nyld_stat_pp_07->GetBinContent(6) * hXD_Efrc_syst_pp_07->GetBinContent(6));
    TCanvas *c1 = new TCanvas();
    hXD_Nyld_stat_pp_07->Draw();
    //!
    //! Make final plots
    //! --- pT spectra
    uPlotPtSpectra2D(h2D_Nres_stat_pp_05, h2D_Nres_syst_pp_05, 30, 5.02, true, "/Users/nrubini/Analysis/ALICE/PWG-LF/PAG-RSN/_1020_Phi_Pair/result/FinalResults/h2D_Nres_pp_05");
    uPlotPtSpectra2D(h2D_Nres_stat_pp_07, h2D_Nres_syst_pp_07, 30, 7.00, true, "/Users/nrubini/Analysis/ALICE/PWG-LF/PAG-RSN/_1020_Phi_Pair/result/FinalResults/h2D_Nres_pp_07");
    //!
    //! --- Mean pT spectra
    uPlotPtSpectra2D(h2D_MeanPT_stat, h2D_MeanPT_syst, 0.5, 5.02, false, "/Users/nrubini/Analysis/ALICE/PWG-LF/PAG-RSN/_1020_Phi_Pair/result/FinalResults/h2D_MeanPT");
    uPlotPtSpectra2DRatio(h2D_MeanPT_stat_MT_pp_05, h2D_MeanPT_syst_MT_pp_05, 1, 5.02, false, "/Users/nrubini/Analysis/ALICE/PWG-LF/PAG-RSN/_1020_Phi_Pair/result/FinalResults/h2D_MeanPT_MT");
    //!
    //! --- Yields
    TCanvas *cPlotFinalQ1 = new TCanvas("", "", 1200, 1000);
    gStyle->SetOptStat(0);
    gPad->SetTopMargin(0.04);
    gPad->SetBottomMargin(0.14);
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.04);
    //!
    TLegend *lPlotLegend = new TLegend(0.18, 0.73, 0.58, 0.83);
    lPlotLegend->SetNColumns(2);
    lPlotLegend->SetFillColorAlpha(kWhite, 0.);
    lPlotLegend->SetLineColorAlpha(kWhite, 0.);
    lPlotLegend->SetTextSize(.03);
    lPlotLegend->AddEntry(g1D_Nres_stat_Mult_pp_05, "pp, #sqrt{s}=5.02 TeV", "P");
    lPlotLegend->AddEntry(gInclusive_stat_pp_07_Y1, "#LT Y_{#phi} #GT", "L");
    lPlotLegend->AddEntry(gInclusive_stat_pp_07_Y1, "pp, #sqrt{s}=7.00 TeV", "P");
    lPlotLegend->AddEntry(g2D_Nres_stat_Mult_pp_05, "#LT Y_{#phi#phi} #GT (x10)", "L");
    TLegend *lPlotLegMC0 = new TLegend(0.18, 0.73, 0.58, 0.63);
    lPlotLegMC0->SetFillColorAlpha(kWhite, 0.);
    lPlotLegMC0->SetLineColorAlpha(kWhite, 0.);
    lPlotLegMC0->SetTextSize(.03);
    lPlotLegMC0->AddEntry(gP800_R2, "Pythia 8 Monash", "L");
    lPlotLegMC0->AddEntry(gP808_R2, "Pythia 8 Monash Ropes", "L");
    //!
    TGraph *hDummyDrawRange = new TGraph();
    hDummyDrawRange->SetPoint(0, 0.5, 0);
    hDummyDrawRange->SetPoint(1, 19.5, 0);
    //!
    TMultiGraph *mg2 = new TMultiGraph();
    uScale(g2D_Nres_stat_Mult_pp_05, 10);
    uScale(g2D_Nres_syst_Mult_pp_05, 10);
    uScale(gP800_N2, 10);
    uScale(gP808_N2, 10);
    uSetXerror(g1D_Nres_stat_Mult_pp_05, 0.1);
    uSetXerror(g1D_Nres_syst_Mult_pp_05, 0.1);
    uSetXerror(g2D_Nres_stat_Mult_pp_05, 0.1);
    uSetXerror(g2D_Nres_syst_Mult_pp_05, 0.1);
    //!
    gInclusive_stat_pp_07_Y1->SetMarkerStyle(kRainbowMarker[1]);
    gInclusive_stat_pp_07_Y1->SetMarkerColor(kRainbowColor[4]);
    gInclusive_stat_pp_07_Y1->SetLineColor(kRainbowColor[4]);
    gInclusive_stat_pp_07_Y1->SetMarkerSize(2);
    gInclusive_syst_pp_07_Y1->SetMarkerStyle(kRainbowMarker[1]);
    gInclusive_syst_pp_07_Y1->SetMarkerColor(kRainbowColor[4]);
    gInclusive_syst_pp_07_Y1->SetLineColor(kRainbowColor[4] - 8);
    gInclusive_syst_pp_07_Y1->SetMarkerSize(2);
    //!
    gInclusive_stat_pp_07_Y2->SetMarkerStyle(kRainbowMarker[1]);
    gInclusive_stat_pp_07_Y2->SetMarkerColor(kRainbowColor[0]);
    gInclusive_stat_pp_07_Y2->SetLineColor(kRainbowColor[0]);
    gInclusive_stat_pp_07_Y2->SetMarkerSize(2);
    gInclusive_syst_pp_07_Y2->SetMarkerStyle(kRainbowMarker[1]);
    gInclusive_syst_pp_07_Y2->SetMarkerColor(kRainbowColor[0]);
    gInclusive_syst_pp_07_Y2->SetLineColor(kRainbowColor[0] - 8);
    gInclusive_syst_pp_07_Y2->SetMarkerSize(2);
    //!
    g1D_Nres_stat_Mult_pp_05->SetMarkerStyle(kRainbowMarker[0]);
    g1D_Nres_stat_Mult_pp_05->SetMarkerColor(kRainbowColor[4]);
    g1D_Nres_stat_Mult_pp_05->SetLineColor(kRainbowColor[4]);
    g1D_Nres_stat_Mult_pp_05->SetMarkerSize(2);
    g1D_Nres_syst_Mult_pp_05->SetMarkerStyle(kRainbowMarker[0]);
    g1D_Nres_syst_Mult_pp_05->SetMarkerColor(kRainbowColor[4]);
    g1D_Nres_syst_Mult_pp_05->SetLineColor(kRainbowColor[4] - 8);
    g1D_Nres_syst_Mult_pp_05->SetMarkerSize(2);
    //!
    g2D_Nres_stat_Mult_pp_05->SetMarkerStyle(kRainbowMarker[0]);
    g2D_Nres_stat_Mult_pp_05->SetMarkerColor(kRainbowColor[0]);
    g2D_Nres_stat_Mult_pp_05->SetLineColor(kRainbowColor[0]);
    g2D_Nres_stat_Mult_pp_05->SetMarkerSize(2);
    g2D_Nres_syst_Mult_pp_05->SetMarkerStyle(kRainbowMarker[0]);
    g2D_Nres_syst_Mult_pp_05->SetMarkerColor(kRainbowColor[0]);
    g2D_Nres_syst_Mult_pp_05->SetLineColor(kRainbowColor[0] - 8);
    g2D_Nres_syst_Mult_pp_05->SetMarkerSize(2);
    //!
    mg2->Add(g1D_Nres_syst_Mult_pp_05, "PE5");
    mg2->Add(g2D_Nres_syst_Mult_pp_05, "PE5");
    mg2->Add(gInclusive_syst_pp_07_Y1, "PE5");
    mg2->Add(gInclusive_syst_pp_07_Y2, "PE5");
    mg2->Add(g1D_Nres_stat_Mult_pp_05, "PE");
    mg2->Add(g2D_Nres_stat_Mult_pp_05, "PE");
    // mg2->Add(gP1_Nres_stat_Mult_pp_05, "PE");
    mg2->Add(gInclusive_stat_pp_07_Y1, "PE");
    mg2->Add(gInclusive_stat_pp_07_Y2, "PE");
    mg2->Add(hDummyDrawRange, "");
    mg2->Add(gP800_N1, "L");
    mg2->Add(gP800_N2, "L");
    mg2->Add(gP808_N1, "L");
    mg2->Add(gP808_N2, "L");
    mg2->Draw("SAME");
    mg2->GetXaxis()->SetRangeUser(0.,20.);
    mg2->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    mg2->GetYaxis()->SetTitle("d#it{N}/(d#it{y})");
    lPlotLegend->Draw("SAME");
    lPlotLegMC0->Draw("SAME");
    uLatex->SetTextFont(60);
    uLatex->SetTextSize(0.05);
    uLatex->DrawLatexNDC(0.19, 0.88, "ALICE");
    uLatex->SetTextFont(42);
    uLatex->SetTextSize(0.04);
    uLatex->DrawLatexNDC(0.19, 0.83, "#phi #rightarrow K^{+}K^{-}, |#it{y} | < 0.5");
    uLatex->SetTextSize(0.035);
    uLatex->DrawLatexNDC(0.62, 0.23, Form("Uncertainties:"));
    uLatex->DrawLatexNDC(0.62, 0.18, Form("Stat. (bars), Syst. (boxes)"));
    // draw an axis on the right side
    /*
    TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(), gPad->GetUymax(),0,0.015,510,"+L");
    axis->SetTextFont(gStyle->GetTextFont());
    axis->SetLabelFont(gStyle->GetLabelFont("y"));
    axis->Draw();
     */
    //!
    cPlotFinalQ1->SaveAs("/Users/nrubini/Analysis/ALICE/PWG-LF/PAG-RSN/_1020_Phi_Pair/result/FinalResults/hFullQ1.pdf");
    //!
    TCanvas *cPlotFinalQ2 = new TCanvas("", "", 1200, 1000);
    gStyle->SetOptStat(0);
    gPad->SetTopMargin(0.04);
    gPad->SetBottomMargin(0.14);
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.04);
    //!
    TLegend *lPlotLegen2 = new TLegend(0.50, 0.79, 0.92, 0.93);
    lPlotLegen2->SetNColumns(2);
    lPlotLegen2->SetFillColorAlpha(kWhite, 0.);
    lPlotLegen2->SetLineColorAlpha(kWhite, 0.);
    lPlotLegen2->SetTextSize(.03);
    TLegend *lPlotLegMC1 = new TLegend(0.50, 0.79, 0.92, 0.63);
    lPlotLegMC1->SetFillColorAlpha(kWhite, 0.);
    lPlotLegMC1->SetLineColorAlpha(kWhite, 0.);
    lPlotLegMC1->SetTextSize(.03);
    //!
    TMultiGraph *mg1 = new TMultiGraph();
    uScale(gP2_Nres_stat_Mult_pp_05, 10);
    uScale(gP2_Nres_syst_Mult_pp_05, 10);
    uSetXerror(gR2_Nres_stat_Mult_pp_05, 0.1);
    uSetXerror(gR2_Nres_syst_Mult_pp_05, 0.1);
    uSetXerror(gP2_Nres_stat_Mult_pp_05, 0.1);
    uSetXerror(gP2_Nres_syst_Mult_pp_05, 0.1);
    uScale(gP800_P2, 10);
    uScale(gP808_P2, 10);
    //!
    gInclusive_stat_pp_07_R2->SetMarkerStyle(kRainbowMarker[1]);
    gInclusive_stat_pp_07_R2->SetMarkerColor(kRainbowColor[4]);
    gInclusive_stat_pp_07_R2->SetLineColor(kRainbowColor[4]);
    gInclusive_stat_pp_07_R2->SetMarkerSize(2);
    gInclusive_syst_pp_07_R2->SetMarkerStyle(kRainbowMarker[1]);
    gInclusive_syst_pp_07_R2->SetMarkerColor(kRainbowColor[4]);
    gInclusive_syst_pp_07_R2->SetLineColor(kRainbowColor[4] - 8);
    gInclusive_syst_pp_07_R2->SetMarkerSize(2);
    //!
    gInclusive_stat_pp_07_P2->SetMarkerStyle(kRainbowMarker[1]);
    gInclusive_stat_pp_07_P2->SetMarkerColor(kRainbowColor[0]);
    gInclusive_stat_pp_07_P2->SetLineColor(kRainbowColor[0]);
    gInclusive_stat_pp_07_P2->SetMarkerSize(2);
    gInclusive_syst_pp_07_P2->SetMarkerStyle(kRainbowMarker[1]);
    gInclusive_syst_pp_07_P2->SetMarkerColor(kRainbowColor[0]);
    gInclusive_syst_pp_07_P2->SetLineColor(kRainbowColor[0] - 8);
    gInclusive_syst_pp_07_P2->SetMarkerSize(2);
    //!
    gR2_Nres_stat_Mult_pp_05->SetMarkerStyle(kRainbowMarker[0]);
    gR2_Nres_stat_Mult_pp_05->SetMarkerColor(kRainbowColor[4]);
    gR2_Nres_stat_Mult_pp_05->SetLineColor(kRainbowColor[4]);
    gR2_Nres_stat_Mult_pp_05->SetMarkerSize(2);
    gR2_Nres_syst_Mult_pp_05->SetMarkerStyle(kRainbowMarker[0]);
    gR2_Nres_syst_Mult_pp_05->SetMarkerColor(kRainbowColor[4]);
    gR2_Nres_syst_Mult_pp_05->SetLineColor(kRainbowColor[4] - 8);
    gR2_Nres_syst_Mult_pp_05->SetMarkerSize(2);
    //!
    gP2_Nres_stat_Mult_pp_05->SetMarkerStyle(kRainbowMarker[0]);
    gP2_Nres_stat_Mult_pp_05->SetMarkerColor(kRainbowColor[0]);
    gP2_Nres_stat_Mult_pp_05->SetLineColor(kRainbowColor[0]);
    gP2_Nres_stat_Mult_pp_05->SetMarkerSize(2);
    gP2_Nres_syst_Mult_pp_05->SetMarkerStyle(kRainbowMarker[0]);
    gP2_Nres_syst_Mult_pp_05->SetMarkerColor(kRainbowColor[0]);
    gP2_Nres_syst_Mult_pp_05->SetLineColor(kRainbowColor[0] - 8);
    gP2_Nres_syst_Mult_pp_05->SetMarkerSize(2);

    // MArker 4 per 7TeV
    //!
    mg1->Add(gR2_Nres_syst_Mult_pp_05, "PE5");
    mg1->Add(gP2_Nres_syst_Mult_pp_05, "PE5");
    mg1->Add(gInclusive_syst_pp_07_R2, "PE5");
    mg1->Add(gInclusive_syst_pp_07_P2, "PE5");
    mg1->Add(gR2_Nres_stat_Mult_pp_05, "PE");
    mg1->Add(gP2_Nres_stat_Mult_pp_05, "PE");
    mg1->Add(gInclusive_stat_pp_07_R2, "PE");
    mg1->Add(gInclusive_stat_pp_07_P2, "PE");
    mg1->Add(gP800_R2, "L");
    mg1->Add(gP800_P2, "L");
    mg1->Add(gP808_R2, "L");
    mg1->Add(gP808_P2, "L");
    mg1->Draw("A");
    mg1->SetMinimum(-0.5);
    mg1->SetMaximum(+2.2);
    mg1->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    mg1->GetYaxis()->SetTitle("Value");
    lPlotLegen2->Draw("SAME");
    lPlotLegMC1->Draw("SAME");
    uLatex->SetTextFont(60);
    uLatex->SetTextSize(0.05);
    uLatex->DrawLatexNDC(0.24, 0.88, "ALICE");
    uLatex->SetTextFont(42);
    uLatex->SetTextSize(0.04);
    uLatex->DrawLatexNDC(0.24, 0.83, "#phi #rightarrow K^{+}K^{-}, |#it{y} | < 0.5");
    uLatex->SetTextSize(0.035);
    uLatex->DrawLatexNDC(0.19, 0.23, Form("Uncertainties:"));
    uLatex->DrawLatexNDC(0.19, 0.18, Form("Stat. (bars), Syst. (boxes)"));
    //!
    TLine *kLine = new TLine();
    kLine->SetLineWidth(3);
    kLine->SetLineStyle(2);
    kLine->SetLineColor(kRainbowColor[4]);
    kLine->DrawLineNDC(0.16, 0.44370370, 0.96, 0.44370370);
    lPlotLegen2->AddEntry(g1D_Nres_stat_Mult_pp_05, "R_{#phi}, #sqrt{s}=5 TeV", "P");
    lPlotLegen2->AddEntry(g2D_Nres_stat_Mult_pp_05, "#gamma_{#phi} (x10), #sqrt{s}=5 TeV", "P");
    lPlotLegen2->AddEntry(gInclusive_stat_pp_07_Y1, "R_{#phi}, #sqrt{s}=7 TeV", "P");
    lPlotLegen2->AddEntry(gInclusive_stat_pp_07_Y2, "#gamma_{#phi} (x10), #sqrt{s}=7 TeV", "P");
    lPlotLegen2->AddEntry((TLine *)(kLine->Clone()), "Pois. Lim.", "L");
    kLine->SetLineColor(kRainbowColor[0]);
    kLine->DrawLineNDC(0.16, 0.29185185, 0.96, 0.29185185);
    lPlotLegen2->AddEntry((TLine *)(kLine->Clone()), "Pois. Lim.", "L");
    lPlotLegMC1->AddEntry(gP800_R2, "Pythia 8 Monash", "L");
    lPlotLegMC1->AddEntry(gP808_R2, "Pythia 8 Monash Ropes", "L");
    //!
    cPlotFinalQ2->SaveAs("/Users/nrubini/Analysis/ALICE/PWG-LF/PAG-RSN/_1020_Phi_Pair/result/FinalResults/hFullQ2.pdf");

    /*
    TCanvas* cPlotFinalQ2 = new TCanvas("","",1200,1000);
    gStyle      -> SetOptStat(0);
    gPad        -> SetTopMargin(0.04);
    gPad        -> SetBottomMargin(0.14);
    gPad        -> SetLeftMargin(0.16);
    gPad        -> SetRightMargin(0.04);
    //!
    TLegend* lPlotLegend2    = new TLegend( 0.19, 0.17, 0.64, 0.25);
    lPlotLegend->SetNColumns(1);
    lPlotLegend->SetFillColorAlpha(kWhite,0.);
    lPlotLegend->SetLineColorAlpha(kWhite,0.);
    lPlotLegend->SetTextSize(.02);
    //!
    TMultiGraph *mg2 = new TMultiGraph();
    uScale(gP2_Nres_stat_Mult_pp_05,10);
    uScale(gP2_Nres_syst_Mult_pp_05,10);
    uSetXerror(gR2_Nres_stat_Mult_pp_05,0.1);
    uSetXerror(gR2_Nres_syst_Mult_pp_05,0.1);
    uSetXerror(gP2_Nres_stat_Mult_pp_05,0.1);
    uSetXerror(gP2_Nres_syst_Mult_pp_05,0.1);
    //!
    gR2_Nres_stat_Mult_pp_05->SetMarkerStyle(kRainbowMarker[0]);
    gR2_Nres_stat_Mult_pp_05->SetMarkerColor(kRainbowColor[0]);
    gR2_Nres_stat_Mult_pp_05->SetLineColor  (kRainbowColor[0]);
    gR2_Nres_stat_Mult_pp_05->SetMarkerSize(2);
    gR2_Nres_syst_Mult_pp_05->SetMarkerStyle(kRainbowMarker[0]);
    gR2_Nres_syst_Mult_pp_05->SetMarkerColor(kRainbowColor[0]);
    gR2_Nres_syst_Mult_pp_05->SetLineColor  (kRainbowColor[0]-8);
    gR2_Nres_syst_Mult_pp_05->SetMarkerSize(2);
    //!
    gP2_Nres_stat_Mult_pp_05->SetMarkerStyle(kRainbowMarker[1]);
    gP2_Nres_stat_Mult_pp_05->SetMarkerColor(kRainbowColor[5]);
    gP2_Nres_stat_Mult_pp_05->SetLineColor  (kRainbowColor[5]);
    gP2_Nres_stat_Mult_pp_05->SetMarkerSize(2);
    gP2_Nres_syst_Mult_pp_05->SetMarkerStyle(kRainbowMarker[1]);
    gP2_Nres_syst_Mult_pp_05->SetMarkerColor(kRainbowColor[5]);
    gP2_Nres_syst_Mult_pp_05->SetLineColor  (kRainbowColor[5]-8);
    gP2_Nres_syst_Mult_pp_05->SetMarkerSize(2);
    //!
    mg2->Add(gR2_Nres_syst_Mult_pp_05, "PE5");
    mg2->Add(gP2_Nres_syst_Mult_pp_05, "PE5");
    mg2->Add(gR2_Nres_stat_Mult_pp_05, "PE");
    mg2->Add(gP2_Nres_stat_Mult_pp_05, "PE");
    mg2->Draw("ALP");
    mg2->SetMinimum(-0.5);
    mg2->SetMaximum(+2.0);
    mg2->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    uLatex->DrawLatexNDC(0.62, 0.23,Form("Uncertainties:"));
    uLatex->DrawLatexNDC(0.62, 0.18,Form("Stat. (bars), Syst. (boxes)"));
    //!
    cPlotFinalQ2->SaveAs("/Users/nrubini/Analysis/ALICE/PWG-LF/PAG-RSN/_1020_Phi_Pair/result/FinalResults/hFullQ2.pdf");
    */
    //
    //  --- ----    ----    A WORKING MESS BUT A MESS STILL
    //
    // --- YIELD ANALYSIS
    /*
    if ( kDoYield ) {
        //  --- Load Files
        kFolder = "_p_p__7TeV";
        TFile*  insFile_Data_YL_pp_07   =   new TFile   ( Form(kASigExtp_FitCheckRst,(TString("Yield")+kFolder).Data()) );
        TFile*  insFile_Syst_YL_pp_07   =   new TFile   ( Form("%s/FullSystematics.root",Form(kAnalysis_Systemt_Dir,  (TString("Yield")+kFolder).Data())) );
        kFolder = "_p_p__5TeV";
        TFile*  insFile_Data_YL_pp_05   =   new TFile   ( Form(kASigExtp_FitCheckRst,(TString("Yield")+kFolder).Data()) );
        kFolder = "_p_p__13TeV";
        TFile*  insFile_Data_YL_pp_13   =   new TFile   ( Form(kASigExtp_FitCheckRst,(TString("Yield")+kFolder).Data()) );
        //TFile*  insFile_Syst_YL_pp_05   =   new TFile   ( Form("%s/FullSystematics.root",Form(kAnalysis_Systemt_Dir,  (TString("Yield")+kFolder).Data())) );
        kFolder = "_p_Pb_5TeV";
        TFile*  insFile_Data_YL_Pb_05   =   new TFile   ( Form(kASigExtp_FitCheckRst,(TString("Yield")+kFolder).Data()) );
        //TFile*  insFile_Syst_YL_Pb_05   =   new TFile   ( Form("%s/FullSystematics.root",Form(kAnalysis_Systemt_Dir,  (TString("Yield")+kFolder).Data())) );
        //
        //  --- Load Histograms
        auto    hXD_Nyld_stat_pp_05     =   uLoadHistograms<0,TH1F> ( insFile_Data_YL_pp_05,    "hXD_Nyld_stat",        "hXD_Nyld_stat_pp_05" );
        auto    hXD_Nyld_stat_pp_07     =   uLoadHistograms<0,TH1F> ( insFile_Data_YL_pp_07,    "hXD_Nyld_stat",        "hXD_Nyld_stat_pp_07" );
        auto    hXD_Nyld_stat_pp_13     =   uLoadHistograms<0,TH1F> ( insFile_Data_YL_pp_13,    "hXD_Nyld_stat",        "hXD_Nyld_stat_pp_13" );
        auto    hXD_Nyld_stat_Pb_05     =   uLoadHistograms<0,TH1F> ( insFile_Data_YL_Pb_05,    "hXD_Nyld_stat",        "hXD_Nyld_stat_Pb_05" );
        auto    hXD_Efrc_syst_pp_07     =   uLoadHistograms<0,TH1F> ( insFile_Syst_YL_pp_07,    "hFullSystematicsRT",   "hFullSystematicsRT_pp_07" );
        //auto    hXD_Efrc_syst_pp_05     =   uLoadHistograms<0,TH1F> ( insFile_Data_YL_pp_05,    "hXD_Nyld_stat",    "hXD_Nyld_stat_pp_05" );
        //auto    hXD_Efrc_syst_Pb_05     =   uLoadHistograms<0,TH1F> ( insFile_Data_YL_Pb_05,    "hXD_Nyld_stat",    "hXD_Nyld_stat_Pb_05" );
        //
        //  --- Output directory
        TString kPlotDirectory          =   Form(kDIR_FinalResultsPlots,(TString("Yield")).Data());
        gROOT   ->  ProcessLine(Form(".! mkdir -p %s",kPlotDirectory.Data()));
        //
        TGraphErrors*   g1DResults_pp_07, *g2DResults_pp_07, *gR1Results_pp_07, *gR2Results_pp_07, *gP1Results_pp_07, *gP2Results_pp_07, *g1DResults_pp_05, *g2DResults_pp_05, *gR1Results_pp_05, *gR2Results_pp_05, *gP1Results_pp_05, *gP2Results_pp_05, *g1DResults_pp_13, *g2DResults_pp_13, *gR1Results_pp_13, *gR2Results_pp_13, *gP1Results_pp_13, *gP2Results_pp_13, *g1DResults_Pb_05, *g2DResults_Pb_05, *gR1Results_Pb_05, *gR2Results_Pb_05, *gP1Results_Pb_05, *gP2Results_Pb_05;
        g1DResults_pp_07    =   new TGraphErrors();
        g2DResults_pp_07    =   new TGraphErrors();
        gR1Results_pp_07    =   new TGraphErrors();
        gR2Results_pp_07    =   new TGraphErrors();
        gP1Results_pp_07    =   new TGraphErrors();
        gP2Results_pp_07    =   new TGraphErrors();
        g1DResults_pp_05    =   new TGraphErrors();
        g2DResults_pp_05    =   new TGraphErrors();
        gR1Results_pp_05    =   new TGraphErrors();
        gR2Results_pp_05    =   new TGraphErrors();
        gP1Results_pp_05    =   new TGraphErrors();
        gP2Results_pp_05    =   new TGraphErrors();
        g1DResults_Pb_05    =   new TGraphErrors();
        g2DResults_Pb_05    =   new TGraphErrors();
        gR1Results_Pb_05    =   new TGraphErrors();
        gR2Results_Pb_05    =   new TGraphErrors();
        gP1Results_Pb_05    =   new TGraphErrors();
        gP2Results_Pb_05    =   new TGraphErrors();
        g1DResults_pp_13    =   new TGraphErrors();
        g2DResults_pp_13    =   new TGraphErrors();
        gR1Results_pp_13    =   new TGraphErrors();
        gR2Results_pp_13    =   new TGraphErrors();
        gP1Results_pp_13    =   new TGraphErrors();
        gP2Results_pp_13    =   new TGraphErrors();
        g1DResults_pp_07    ->  SetMarkerStyle( uGetMarker(1) );
        g2DResults_pp_07    ->  SetMarkerStyle( uGetMarker(1) );
        gR1Results_pp_07    ->  SetMarkerStyle( uGetMarker(1) );
        gR2Results_pp_07    ->  SetMarkerStyle( uGetMarker(1) );
        gP1Results_pp_07    ->  SetMarkerStyle( uGetMarker(1) );
        gP2Results_pp_07    ->  SetMarkerStyle( uGetMarker(1) );
        g1DResults_pp_05    ->  SetMarkerStyle( uGetMarker(2) );
        g2DResults_pp_05    ->  SetMarkerStyle( uGetMarker(2) );
        gR1Results_pp_05    ->  SetMarkerStyle( uGetMarker(2) );
        gR2Results_pp_05    ->  SetMarkerStyle( uGetMarker(2) );
        gP1Results_pp_05    ->  SetMarkerStyle( uGetMarker(2) );
        gP2Results_pp_05    ->  SetMarkerStyle( uGetMarker(2) );
        g1DResults_Pb_05    ->  SetMarkerStyle( uGetMarker(3) );
        g2DResults_Pb_05    ->  SetMarkerStyle( uGetMarker(3) );
        gR1Results_Pb_05    ->  SetMarkerStyle( uGetMarker(3) );
        gR2Results_Pb_05    ->  SetMarkerStyle( uGetMarker(3) );
        gP1Results_Pb_05    ->  SetMarkerStyle( uGetMarker(3) );
        gP2Results_Pb_05    ->  SetMarkerStyle( uGetMarker(3) );
        g1DResults_pp_13    ->  SetMarkerStyle( uGetMarker(4) );
        g2DResults_pp_13    ->  SetMarkerStyle( uGetMarker(4) );
        gR1Results_pp_13    ->  SetMarkerStyle( uGetMarker(4) );
        gR2Results_pp_13    ->  SetMarkerStyle( uGetMarker(4) );
        gP1Results_pp_13    ->  SetMarkerStyle( uGetMarker(4) );
        gP2Results_pp_13    ->  SetMarkerStyle( uGetMarker(4) );
        g1DResults_pp_07    ->  SetMarkerColor( uGetColor(1) );
        g2DResults_pp_07    ->  SetMarkerColor( uGetColor(1) );
        gR1Results_pp_07    ->  SetMarkerColor( uGetColor(1) );
        gR2Results_pp_07    ->  SetMarkerColor( uGetColor(1) );
        gP1Results_pp_07    ->  SetMarkerColor( uGetColor(1) );
        gP2Results_pp_07    ->  SetMarkerColor( uGetColor(1) );
        g1DResults_pp_05    ->  SetMarkerColor( uGetColor(2) );
        g2DResults_pp_05    ->  SetMarkerColor( uGetColor(2) );
        gR1Results_pp_05    ->  SetMarkerColor( uGetColor(2) );
        gR2Results_pp_05    ->  SetMarkerColor( uGetColor(2) );
        gP1Results_pp_05    ->  SetMarkerColor( uGetColor(2) );
        gP2Results_pp_05    ->  SetMarkerColor( uGetColor(2) );
        g1DResults_Pb_05    ->  SetMarkerColor( uGetColor(3) );
        g2DResults_Pb_05    ->  SetMarkerColor( uGetColor(3) );
        gR1Results_Pb_05    ->  SetMarkerColor( uGetColor(3) );
        gR2Results_Pb_05    ->  SetMarkerColor( uGetColor(3) );
        gP1Results_Pb_05    ->  SetMarkerColor( uGetColor(3) );
        gP2Results_Pb_05    ->  SetMarkerColor( uGetColor(3) );
        g1DResults_pp_13    ->  SetMarkerColor( uGetColor(4) );
        g2DResults_pp_13    ->  SetMarkerColor( uGetColor(4) );
        gR1Results_pp_13    ->  SetMarkerColor( uGetColor(4) );
        gR2Results_pp_13    ->  SetMarkerColor( uGetColor(4) );
        gP1Results_pp_13    ->  SetMarkerColor( uGetColor(4) );
        gP2Results_pp_13    ->  SetMarkerColor( uGetColor(4) );
        g1DResults_pp_07    ->  SetLineColor( uGetColor(1) );
        g2DResults_pp_07    ->  SetLineColor( uGetColor(1) );
        gR1Results_pp_07    ->  SetLineColor( uGetColor(1) );
        gR2Results_pp_07    ->  SetLineColor( uGetColor(1) );
        gP1Results_pp_07    ->  SetLineColor( uGetColor(1) );
        gP2Results_pp_07    ->  SetLineColor( uGetColor(1) );
        g1DResults_pp_05    ->  SetLineColor( uGetColor(2) );
        g2DResults_pp_05    ->  SetLineColor( uGetColor(2) );
        gR1Results_pp_05    ->  SetLineColor( uGetColor(2) );
        gR2Results_pp_05    ->  SetLineColor( uGetColor(2) );
        gP1Results_pp_05    ->  SetLineColor( uGetColor(2) );
        gP2Results_pp_05    ->  SetLineColor( uGetColor(2) );
        g1DResults_Pb_05    ->  SetLineColor( uGetColor(3) );
        g2DResults_Pb_05    ->  SetLineColor( uGetColor(3) );
        gR1Results_Pb_05    ->  SetLineColor( uGetColor(3) );
        gR2Results_Pb_05    ->  SetLineColor( uGetColor(3) );
        gP1Results_Pb_05    ->  SetLineColor( uGetColor(3) );
        gP2Results_Pb_05    ->  SetLineColor( uGetColor(3) );
        g1DResults_pp_13    ->  SetLineColor( uGetColor(4) );
        g2DResults_pp_13    ->  SetLineColor( uGetColor(4) );
        gR1Results_pp_13    ->  SetLineColor( uGetColor(4) );
        gR2Results_pp_13    ->  SetLineColor( uGetColor(4) );
        gP1Results_pp_13    ->  SetLineColor( uGetColor(4) );
        gP2Results_pp_13    ->  SetLineColor( uGetColor(4) );
        g1DResults_pp_07    ->  SetPoint( 0, 7, hXD_Nyld_stat_pp_07->GetBinContent(1) );
        g2DResults_pp_07    ->  SetPoint( 0, 7, hXD_Nyld_stat_pp_07->GetBinContent(2) );
        gR1Results_pp_07    ->  SetPoint( 0, 7, hXD_Nyld_stat_pp_07->GetBinContent(2)/(hXD_Nyld_stat_pp_07->GetBinContent(1)) );
        gR2Results_pp_07    ->  SetPoint( 0, 7, hXD_Nyld_stat_pp_07->GetBinContent(2)/(hXD_Nyld_stat_pp_07->GetBinContent(1)*hXD_Nyld_stat_pp_07->GetBinContent(1)) );
        gP1Results_pp_07    ->  SetPoint( 0, 7, fSigmaPhiValue(hXD_Nyld_stat_pp_07->GetBinContent(1),hXD_Nyld_stat_pp_07->GetBinContent(2)) );
        gP2Results_pp_07    ->  SetPoint( 0, 7, fGammaPhiValue(hXD_Nyld_stat_pp_07->GetBinContent(1),hXD_Nyld_stat_pp_07->GetBinContent(2)) );
        g1DResults_pp_07    ->  SetPointError( 0, 0, hXD_Nyld_stat_pp_07->GetBinError(1) );
        g2DResults_pp_07    ->  SetPointError( 0, 0, hXD_Nyld_stat_pp_07->GetBinError(2) );
        gR1Results_pp_07    ->  SetPointError( 0, 0, ( SquareSum( { hXD_Nyld_stat_pp_07->GetBinError(1)/(hXD_Nyld_stat_pp_07->GetBinContent(1)), hXD_Nyld_stat_pp_07->GetBinError(2)/(hXD_Nyld_stat_pp_07->GetBinContent(2)) }) )*hXD_Nyld_stat_pp_07->GetBinContent(2)/(hXD_Nyld_stat_pp_07->GetBinContent(1)) );
        gR2Results_pp_07    ->  SetPointError( 0, 0, ( SquareSum( { hXD_Nyld_stat_pp_07->GetBinError(1)/(hXD_Nyld_stat_pp_07->GetBinContent(1)), hXD_Nyld_stat_pp_07->GetBinError(1)/(hXD_Nyld_stat_pp_07->GetBinContent(1)), hXD_Nyld_stat_pp_07->GetBinError(2)/(hXD_Nyld_stat_pp_07->GetBinContent(2)) }) )*hXD_Nyld_stat_pp_07->GetBinContent(2)/(hXD_Nyld_stat_pp_07->GetBinContent(1)*hXD_Nyld_stat_pp_07->GetBinContent(1)) );
        gP1Results_pp_07    ->  SetPointError( 0, 0, fSigmaPhiError(hXD_Nyld_stat_pp_07->GetBinContent(1),hXD_Nyld_stat_pp_07->GetBinContent(2),hXD_Nyld_stat_pp_07->GetBinError(1),hXD_Nyld_stat_pp_07->GetBinError(2)) );
        gP2Results_pp_07    ->  SetPointError( 0, 0, fGammaPhiError(hXD_Nyld_stat_pp_07->GetBinContent(1),hXD_Nyld_stat_pp_07->GetBinContent(2),hXD_Nyld_stat_pp_07->GetBinError(1),hXD_Nyld_stat_pp_07->GetBinError(2))  );
        g1DResults_pp_05    ->  SetPoint( 0, 5, hXD_Nyld_stat_pp_05->GetBinContent(1) );
        g2DResults_pp_05    ->  SetPoint( 0, 5, hXD_Nyld_stat_pp_05->GetBinContent(2) );
        gR1Results_pp_05    ->  SetPoint( 0, 5, hXD_Nyld_stat_pp_05->GetBinContent(2)/(hXD_Nyld_stat_pp_05->GetBinContent(1)) );
        gR2Results_pp_05    ->  SetPoint( 0, 5, hXD_Nyld_stat_pp_05->GetBinContent(2)/(hXD_Nyld_stat_pp_05->GetBinContent(1)*hXD_Nyld_stat_pp_05->GetBinContent(1)) );
        gP1Results_pp_05    ->  SetPoint( 0, 5, fSigmaPhiValue(hXD_Nyld_stat_pp_05->GetBinContent(1),hXD_Nyld_stat_pp_05->GetBinContent(2)) );
        gP2Results_pp_05    ->  SetPoint( 0, 5, fGammaPhiValue(hXD_Nyld_stat_pp_05->GetBinContent(1),hXD_Nyld_stat_pp_05->GetBinContent(2)) );
        g1DResults_pp_05    ->  SetPointError( 0, 0, hXD_Nyld_stat_pp_05->GetBinError(1) );
        g2DResults_pp_05    ->  SetPointError( 0, 0, hXD_Nyld_stat_pp_05->GetBinError(2) );
        gR1Results_pp_05    ->  SetPointError( 0, 0, ( SquareSum( { hXD_Nyld_stat_pp_05->GetBinError(1)/(hXD_Nyld_stat_pp_05->GetBinContent(1)), hXD_Nyld_stat_pp_05->GetBinError(2)/(hXD_Nyld_stat_pp_05->GetBinContent(2)) }) )*hXD_Nyld_stat_pp_05->GetBinContent(2)/(hXD_Nyld_stat_pp_05->GetBinContent(1)) );
        gR2Results_pp_05    ->  SetPointError( 0, 0, ( SquareSum( { hXD_Nyld_stat_pp_05->GetBinError(1)/(hXD_Nyld_stat_pp_05->GetBinContent(1)), hXD_Nyld_stat_pp_05->GetBinError(1)/(hXD_Nyld_stat_pp_05->GetBinContent(1)), hXD_Nyld_stat_pp_05->GetBinError(2)/(hXD_Nyld_stat_pp_05->GetBinContent(2)) }) )*hXD_Nyld_stat_pp_05->GetBinContent(2)/(hXD_Nyld_stat_pp_05->GetBinContent(1)*hXD_Nyld_stat_pp_05->GetBinContent(1)) );
        gP1Results_pp_05    ->  SetPointError( 0, 0, fSigmaPhiError(hXD_Nyld_stat_pp_05->GetBinContent(1),hXD_Nyld_stat_pp_05->GetBinContent(2),hXD_Nyld_stat_pp_05->GetBinError(1),hXD_Nyld_stat_pp_05->GetBinError(2)) );
        gP2Results_pp_05    ->  SetPointError( 0, 0, fGammaPhiError(hXD_Nyld_stat_pp_05->GetBinContent(1),hXD_Nyld_stat_pp_05->GetBinContent(2),hXD_Nyld_stat_pp_05->GetBinError(1),hXD_Nyld_stat_pp_05->GetBinError(2))  );
        g1DResults_pp_13    ->  SetPoint( 0, 13, hXD_Nyld_stat_pp_13->GetBinContent(1) );
        g2DResults_pp_13    ->  SetPoint( 0, 13, hXD_Nyld_stat_pp_13->GetBinContent(2) );
        gR1Results_pp_13    ->  SetPoint( 0, 13, hXD_Nyld_stat_pp_13->GetBinContent(2)/(hXD_Nyld_stat_pp_13->GetBinContent(1)) );
        gR2Results_pp_13    ->  SetPoint( 0, 13, hXD_Nyld_stat_pp_13->GetBinContent(2)/(hXD_Nyld_stat_pp_13->GetBinContent(1)*hXD_Nyld_stat_pp_13->GetBinContent(1)) );
        gP1Results_pp_13    ->  SetPoint( 0, 13, fSigmaPhiValue(hXD_Nyld_stat_pp_13->GetBinContent(1),hXD_Nyld_stat_pp_13->GetBinContent(2)) );
        gP2Results_pp_13    ->  SetPoint( 0, 13, fGammaPhiValue(hXD_Nyld_stat_pp_13->GetBinContent(1),hXD_Nyld_stat_pp_13->GetBinContent(2)) );
        g1DResults_pp_13    ->  SetPointError( 0, 0, hXD_Nyld_stat_pp_13->GetBinError(1) );
        g2DResults_pp_13    ->  SetPointError( 0, 0, hXD_Nyld_stat_pp_13->GetBinError(2) );
        gR1Results_pp_13    ->  SetPointError( 0, 0, ( SquareSum( { hXD_Nyld_stat_pp_13->GetBinError(1)/(hXD_Nyld_stat_pp_13->GetBinContent(1)), hXD_Nyld_stat_pp_13->GetBinError(2)/(hXD_Nyld_stat_pp_13->GetBinContent(2)) }) )*hXD_Nyld_stat_pp_13->GetBinContent(2)/(hXD_Nyld_stat_pp_13->GetBinContent(1)) );
        gR2Results_pp_13    ->  SetPointError( 0, 0, ( SquareSum( { hXD_Nyld_stat_pp_13->GetBinError(1)/(hXD_Nyld_stat_pp_13->GetBinContent(1)), hXD_Nyld_stat_pp_13->GetBinError(1)/(hXD_Nyld_stat_pp_13->GetBinContent(1)), hXD_Nyld_stat_pp_13->GetBinError(2)/(hXD_Nyld_stat_pp_13->GetBinContent(2)) }) )*hXD_Nyld_stat_pp_13->GetBinContent(2)/(hXD_Nyld_stat_pp_13->GetBinContent(1)*hXD_Nyld_stat_pp_13->GetBinContent(1)) );
        gP1Results_pp_13    ->  SetPointError( 0, 0, fSigmaPhiError(hXD_Nyld_stat_pp_13->GetBinContent(1),hXD_Nyld_stat_pp_13->GetBinContent(2),hXD_Nyld_stat_pp_13->GetBinError(1),hXD_Nyld_stat_pp_13->GetBinError(2)) );
        gP2Results_pp_13    ->  SetPointError( 0, 0, fGammaPhiError(hXD_Nyld_stat_pp_13->GetBinContent(1),hXD_Nyld_stat_pp_13->GetBinContent(2),hXD_Nyld_stat_pp_13->GetBinError(1),hXD_Nyld_stat_pp_13->GetBinError(2))  );
        g1DResults_Pb_05    ->  SetPoint( 0, 5, hXD_Nyld_stat_Pb_05->GetBinContent(1) );
        g2DResults_Pb_05    ->  SetPoint( 0, 5, hXD_Nyld_stat_Pb_05->GetBinContent(2) );
        gR1Results_Pb_05    ->  SetPoint( 0, 5, hXD_Nyld_stat_Pb_05->GetBinContent(2)/(hXD_Nyld_stat_Pb_05->GetBinContent(1)) );
        gR2Results_Pb_05    ->  SetPoint( 0, 5, hXD_Nyld_stat_Pb_05->GetBinContent(2)/(hXD_Nyld_stat_Pb_05->GetBinContent(1)*hXD_Nyld_stat_Pb_05->GetBinContent(1)) );
        gP1Results_Pb_05    ->  SetPoint( 0, 5, fSigmaPhiValue(hXD_Nyld_stat_Pb_05->GetBinContent(1),hXD_Nyld_stat_Pb_05->GetBinContent(2)) );
        gP2Results_Pb_05    ->  SetPoint( 0, 5, fGammaPhiValue(hXD_Nyld_stat_Pb_05->GetBinContent(1),hXD_Nyld_stat_Pb_05->GetBinContent(2)) );
        g1DResults_Pb_05    ->  SetPointError( 0, 0, hXD_Nyld_stat_Pb_05->GetBinError(1) );
        g2DResults_Pb_05    ->  SetPointError( 0, 0, hXD_Nyld_stat_Pb_05->GetBinError(2) );
        gR1Results_Pb_05    ->  SetPointError( 0, 0, ( SquareSum( { hXD_Nyld_stat_Pb_05->GetBinError(1)/(hXD_Nyld_stat_Pb_05->GetBinContent(1)), hXD_Nyld_stat_Pb_05->GetBinError(2)/(hXD_Nyld_stat_Pb_05->GetBinContent(2)) }) )*hXD_Nyld_stat_Pb_05->GetBinContent(2)/(hXD_Nyld_stat_Pb_05->GetBinContent(1)) );
        gR2Results_Pb_05    ->  SetPointError( 0, 0, ( SquareSum( { hXD_Nyld_stat_Pb_05->GetBinError(1)/(hXD_Nyld_stat_Pb_05->GetBinContent(1)), hXD_Nyld_stat_Pb_05->GetBinError(1)/(hXD_Nyld_stat_Pb_05->GetBinContent(1)), hXD_Nyld_stat_Pb_05->GetBinError(2)/(hXD_Nyld_stat_Pb_05->GetBinContent(2)) }) )*hXD_Nyld_stat_Pb_05->GetBinContent(2)/(hXD_Nyld_stat_Pb_05->GetBinContent(1)*hXD_Nyld_stat_Pb_05->GetBinContent(1)) );
        gP1Results_Pb_05    ->  SetPointError( 0, 0, fSigmaPhiError(hXD_Nyld_stat_Pb_05->GetBinContent(1),hXD_Nyld_stat_Pb_05->GetBinContent(2),hXD_Nyld_stat_Pb_05->GetBinError(1),hXD_Nyld_stat_Pb_05->GetBinError(2)) );
        gP2Results_Pb_05    ->  SetPointError( 0, 0, fGammaPhiError(hXD_Nyld_stat_Pb_05->GetBinContent(1),hXD_Nyld_stat_Pb_05->GetBinContent(2),hXD_Nyld_stat_Pb_05->GetBinError(1),hXD_Nyld_stat_Pb_05->GetBinError(2))  );
        //
        TLegend*    cDrawLegend =   new TLegend( 0.15,0.55,0.45,0.9);
        cDrawLegend     ->  AddEntry( g1DResults_pp_05, "pp @5.02TeV", "EP" );
        cDrawLegend     ->  AddEntry( g1DResults_pp_07, "pp @7.00TeV", "EP" );
        cDrawLegend     ->  AddEntry( g1DResults_pp_13, "pp @13.0TeV", "EP" );
        cDrawLegend     ->  AddEntry( g1DResults_Pb_05, "pPb @5.02TeV", "EP" );
        cDrawLegend     ->  SetFillColorAlpha(0.,0.);
        cDrawLegend     ->  SetLineColorAlpha(0.,0.);
        //
        TH1F*   hDrawRange  = new TH1F( "hDrawRange", "hDrawRange", 1, 0, 15 );
        hDrawRange->GetXaxis()->SetTitle("#sqrt{s} (TeV)");
        //
        TCanvas*    cDrawResults    =   new TCanvas( "cDrawResults", "cDrawResults", 1500, 1500 );
        cDrawResults    ->  Divide(3,2);
        gStyle->SetOptStat(0000);
        //
        cDrawResults    ->  cd(1);
        hDrawRange->SetTitle( "#LT Y_{#phi} #GT" );
        hDrawRange->SetMinimum(0);
        hDrawRange->SetMaximum(0.2);
        hDrawRange->DrawCopy("");
        g1DResults_pp_07->Draw("EP SAME");
        g1DResults_pp_05->Draw("EP SAME");
        g1DResults_pp_13->Draw("EP SAME");
        g1DResults_Pb_05->Draw("EP SAME");
        cDrawLegend->Draw("SAME");
        //
        cDrawResults    ->  cd(2);
        hDrawRange->SetTitle( "#LT Y_{#phi#phi} #GT / #LT Y_{#phi} #GT" );
        hDrawRange->SetMinimum(0);
        hDrawRange->SetMaximum(0.1);
        hDrawRange->DrawCopy("");
        gR1Results_pp_07->Draw("EP SAME");
        gR1Results_pp_05->Draw("EP SAME");
        gR1Results_pp_13->Draw("EP SAME");
        gR1Results_Pb_05->Draw("EP SAME");
        //
        cDrawResults    ->  cd(3);
        hDrawRange->SetTitle( "#sigma_{#phi}" );
        hDrawRange->SetMinimum(0);
        hDrawRange->SetMaximum(0.2);
        hDrawRange->DrawCopy("");
        gP1Results_pp_07->Draw("EP SAME");
        gR1Results_pp_05->Draw("EP SAME");
        gR1Results_pp_13->Draw("EP SAME");
        gP1Results_Pb_05->Draw("EP SAME");
        //
        cDrawResults    ->  cd(4);
        hDrawRange->SetTitle( "#LT Y_{#phi#phi} #GT" );
        hDrawRange->SetMinimum(0);
        hDrawRange->SetMaximum(0.02);
        hDrawRange->DrawCopy("");
        g2DResults_pp_07->Draw("EP SAME");
        g2DResults_pp_05->Draw("EP SAME");
        g2DResults_pp_13->Draw("EP SAME");
        g2DResults_Pb_05->Draw("EP SAME");
        //
        cDrawResults    ->  cd(5);
        hDrawRange->SetTitle( "#LT Y_{#phi#phi} #GT / #LT Y_{#phi} #GT^{2}" );
        hDrawRange->SetBinContent(1,0.5);
        hDrawRange->SetMinimum(0);
        hDrawRange->SetMaximum(2.5);
        hDrawRange->DrawCopy("");
        gR2Results_pp_07->Draw("EP SAME");
        gR2Results_pp_05->Draw("EP SAME");
        gR2Results_pp_13->Draw("EP SAME");
        gR2Results_Pb_05->Draw("EP SAME");
        //
        cDrawResults    ->  cd(6);
        hDrawRange->SetTitle( "#gamma_{#phi}" );
        hDrawRange->SetBinContent(1,0.);
        hDrawRange->SetMinimum(-0.05);
        hDrawRange->SetMaximum(+0.15);
        hDrawRange->DrawCopy("");
        gP2Results_pp_07->Draw("EP SAME");
        gP2Results_pp_05->Draw("EP SAME");
        gP2Results_pp_13->Draw("EP SAME");
        gP2Results_Pb_05->Draw("EP SAME");
        //
        cDrawResults    ->  SaveAs(kPlotDirectory+TString("EnergyDependence.pdf"));
        //
        delete cDrawResults;
        //
        //  --- Final Results
        TFile*  outFile   =   new TFile   ( Form(kDIR_FinalResults,(TString("Yield")+kFolder).Data()) );
        //
        std::vector<TString> kShowWhat = {"Y_1","Y_2","R_1","R_2","P_1","P_2"};
        for ( Int_t iTer = 0; iTer < 6; iTer++ )    {
            auto k1Yield    =   hXD_Nyld_stat_pp_07->GetBinContent(1);
            auto k2Yield    =   hXD_Nyld_stat_pp_07->GetBinContent(2);
            auto k1YieldSt  =   hXD_Nyld_stat_pp_07->GetBinError(1);
            auto k2YieldSt  =   hXD_Nyld_stat_pp_07->GetBinError(2);
            auto kHigNorm   =   hXD_Nyld_stat_pp_07->GetBinContent(iTer+1)*kSysHig_TE;
            auto kLowNorm   =   hXD_Nyld_stat_pp_07->GetBinContent(iTer+1)*kSysLow_TE;
            auto kCurrent   =   0.;
            auto kCurrStat  =   0.;
            switch ( iTer ) {
                case 0:
                    kCurrent    =   k1Yield;
                    kHigNorm    =   kCurrent*kSysHig_TE;
                    kLowNorm    =   kCurrent*kSysLow_TE;
                    kCurrStat   =   k1YieldSt;
                    break;
                case 1:
                    kCurrent    =   k2Yield;
                    kHigNorm    =   kCurrent*kSysHig_TE;
                    kLowNorm    =   kCurrent*kSysLow_TE;
                    kCurrStat   =   k2YieldSt;
                    break;
                case 2:
                    kCurrent    =   k2Yield/k1Yield;
                    kHigNorm    =   0;
                    kLowNorm    =   0;
                    kCurrStat   =   (k2Yield/k1Yield)*SquareSum( {k1YieldSt/k1Yield, k2YieldSt/k2Yield} );
                    break;
                case 3:
                    kCurrent    =   k2Yield/( k1Yield * k1Yield );
                    kHigNorm    =   k2Yield/(k1Yield*(1+kSysLow_TE));
                    kLowNorm    =   k2Yield/(k1Yield*(1+kSysHig_TE));
                    kCurrStat   =   (k2Yield/(k1Yield*k1Yield))*SquareSum( {4*k1YieldSt/k1Yield, k2YieldSt/k2Yield} );
                    break;
                case 4:
                    kCurrent    =   fSigmaPhiValue(k1Yield,k2Yield);
                    kHigNorm    =   fabs( fSigmaPhiValue(k1Yield*(1+kSysHig_TE),k2Yield*(1+kSysHig_TE)) - kCurrent );
                    kLowNorm    =   fabs( fSigmaPhiValue(k1Yield*(1-kSysLow_TE),k2Yield*(1-kSysLow_TE)) - kCurrent );
                    kCurrStat   =   fSigmaPhiError(k1Yield,k2Yield,k1YieldSt,k2YieldSt);
                    break;
                case 5:
                    kCurrent    =   fGammaPhiValue(k1Yield,k2Yield) ;
                    kHigNorm    =   fabs( fGammaPhiValue(k1Yield*(1-kSysLow_TE),k2Yield*(1-kSysLow_TE)) - kCurrent );
                    kLowNorm    =   fabs( fGammaPhiValue(k1Yield*(1+kSysHig_TE),k2Yield*(1+kSysHig_TE)) - kCurrent );
                    kCurrStat   =   fGammaPhiError(k1Yield,k2Yield,k1YieldSt,k2YieldSt);
                    break;
            }
            cout << "-------------------------------" << endl;
            cout << "-       " << "       +" << kCurrStat << " +" << kCurrent*hXD_Efrc_syst_pp_07->GetBinContent(iTer+1) << " +" << kHigNorm << endl;
            cout << "- " << kShowWhat.at(iTer).Data() << " : " << kCurrent << endl;
            cout << "-       " << "       -" << kCurrStat << " -" << kCurrent*hXD_Efrc_syst_pp_07->GetBinContent(iTer+1) << " -" << kLowNorm << endl;
            cout << "-------------------------------" << endl;
        }
        /*
        cout << "-------------------------------" << endl;
        cout << "-       " << "       +" << hXD_Nyld_stat->GetBinError(3) << " +" << hXD_Nyld_stat->GetBinContent(3)*hReferenc2->GetBinContent(1) << endl;
        cout << "- " << "PT" << " : " << hXD_Nyld_stat->GetBinContent(3) << endl;
        cout << "-       " << "       -" << hXD_Nyld_stat->GetBinError(3) << " -" << hXD_Nyld_stat->GetBinContent(3)*hReferenc2->GetBinContent(1) << endl;
        cout << "-------------------------------" << endl;

        //
    }
*/
    //
    /*
    if ( kDoMultiplicity ) {
        //  --- Load Files
        kFolder = "_p_p__5TeV";
        TFile*  insFile_Data_ML_pp  =   new TFile   ( Form(kASigExtp_FitCheckRst,(TString("Multiplicity")+kFolder).Data()) );
        TFile*  insFile_Syst_ML_pp  =   new TFile   ( Form("%s/FullSystematics.root",Form(kAnalysis_Systemt_Dir,  (TString("Multiplicity")+kFolder).Data())) );
        kFolder = "_p_p__13TeV";
        TFile*  insFile_Data_ML_p3  =   new TFile   ( Form(kASigExtp_FitCheckRst,(TString("Multiplicity")+kFolder).Data()) );
        TFile*  insFile_Syst_ML_p3  =   new TFile   ( Form("%s/FullSystematics.root",Form(kAnalysis_Systemt_Dir,  (TString("Multiplicity")+kFolder).Data())) );
        kFolder = "_p_Pb_5TeV";
        TFile*  insFile_Data_ML_Pb  =   new TFile   ( Form(kASigExtp_FitCheckRst,(TString("Multiplicity")+kFolder).Data()) );
        TFile*  insFile_Syst_ML_Pb  =   new TFile   ( Form("%s/FullSystematics.root",Form(kAnalysis_Systemt_Dir,  (TString("Multiplicity")+kFolder).Data())) );
        //
        //  --- Load Histograms
        auto        g1DResults_pp   =   uLoadHistograms<0,TGraphErrors> ( insFile_Data_ML_pp,  "g1D_Nres_Mult",    "g1DResults_pp" );
        auto        g2DResults_pp   =   uLoadHistograms<0,TGraphErrors> ( insFile_Data_ML_pp,  "g2D_Nres_Mult",    "g2DResults_pp" );
        auto        gR1Results_pp   =   uLoadHistograms<0,TGraphErrors> ( insFile_Data_ML_pp,  "g2D_Nres_Mult",    "gR1Results_pp" );
        auto        gR2Results_pp   =   uLoadHistograms<0,TGraphErrors> ( insFile_Data_ML_pp,  "g2D_Nres_Mult",    "gR1Results_pp" );
        auto        gP1Results_pp   =   uLoadHistograms<0,TGraphErrors> ( insFile_Data_ML_pp,  "g2D_Nres_Mult",    "gP1Results_pp" );
        auto        gP2Results_pp   =   uLoadHistograms<0,TGraphErrors> ( insFile_Data_ML_pp,  "g2D_Nres_Mult",    "gP2Results_pp" );
        auto        g1DResults_p3   =   uLoadHistograms<0,TGraphErrors> ( insFile_Data_ML_p3,  "g1D_Nres_Mult",    "g1DResults_p3" );
        auto        g2DResults_p3   =   uLoadHistograms<0,TGraphErrors> ( insFile_Data_ML_p3,  "g2D_Nres_Mult",    "g2DResults_p3" );
        auto        gR1Results_p3   =   uLoadHistograms<0,TGraphErrors> ( insFile_Data_ML_p3,  "g2D_Nres_Mult",    "gR1Results_p3" );
        auto        gR2Results_p3   =   uLoadHistograms<0,TGraphErrors> ( insFile_Data_ML_p3,  "g2D_Nres_Mult",    "gR1Results_p3" );
        auto        gP1Results_p3   =   uLoadHistograms<0,TGraphErrors> ( insFile_Data_ML_p3,  "g2D_Nres_Mult",    "gP1Results_p3" );
        auto        gP2Results_p3   =   uLoadHistograms<0,TGraphErrors> ( insFile_Data_ML_p3,  "g2D_Nres_Mult",    "gP2Results_p3" );
        auto        g1DResults_Pb   =   uLoadHistograms<0,TGraphErrors> ( insFile_Data_ML_Pb,  "g1D_Nres_Mult",    "g1DResults_Pb" );
        auto        g2DResults_Pb   =   uLoadHistograms<0,TGraphErrors> ( insFile_Data_ML_Pb,  "g2D_Nres_Mult",    "g2DResults_Pb" );
        auto        gR1Results_Pb   =   uLoadHistograms<0,TGraphErrors> ( insFile_Data_ML_Pb,  "g1D_Nres_Mult",    "gR1Results_Pb" );
        auto        gR2Results_Pb   =   uLoadHistograms<0,TGraphErrors> ( insFile_Data_ML_Pb,  "g2D_Nres_Mult",    "gR2Results_Pb" );
        auto        gP1Results_Pb   =   uLoadHistograms<0,TGraphErrors> ( insFile_Data_ML_Pb,  "g1D_Nres_Mult",    "gP1Results_Pb" );
        auto        gP2Results_Pb   =   uLoadHistograms<0,TGraphErrors> ( insFile_Data_ML_Pb,  "g2D_Nres_Mult",    "gP2Results_Pb" );
        //
        TFile*  insFile_MCPr_YL_pp_05_EPOS      =   new TFile   ( "/Users/nikolajal/alice/AliAnalysisPhiCount/result/Yield_p_p__7TeV/MC_Production/ComparePlots_EPOS_5TeV_test.root" );
        TFile*  insFile_MCPr_YL_pp_05_PYTHIA6   =   new TFile   ( "/Users/nikolajal/alice/AliAnalysisPhiCount/result/Yield_p_p__7TeV/MC_Production/ComparePlots_PYTHIA6_5TeV.root" );
        TFile*  insFile_MCPr_YL_pp_05_PYTHIA8   =   new TFile   ( "/Users/nikolajal/alice/AliAnalysisPhiCount/result/Yield_p_p__7TeV/MC_Production/ComparePlots_PYTHIA8_5TeV.root" );
        auto    h1D_Ntru_Mult_EPOS      =   uLoadHistograms<0,TH1F> ( insFile_MCPr_YL_pp_05_EPOS,    "h1D_Ntru_Mult",   "h1D_Ntru_Mult_EPOS" );
        auto    h2D_Ntru_Mult_EPOS      =   uLoadHistograms<0,TH1F> ( insFile_MCPr_YL_pp_05_EPOS,    "h2D_Ntru_Mult",   "h2D_Ntru_Mult_EPOS" );
        auto    hR1_Ntru_Mult_EPOS      =   uLoadHistograms<0,TH1F> ( insFile_MCPr_YL_pp_05_EPOS,    "hR1_Ntru_Mult",   "hR1_Ntru_Mult_EPOS" );
        auto    hR2_Ntru_Mult_EPOS      =   uLoadHistograms<0,TH1F> ( insFile_MCPr_YL_pp_05_EPOS,    "hR2_Ntru_Mult",   "hR2_Ntru_Mult_EPOS" );
        auto    hP1_Ntru_Mult_EPOS      =   uLoadHistograms<0,TH1F> ( insFile_MCPr_YL_pp_05_EPOS,    "hP1_Ntru_Mult",   "hP1_Ntru_Mult_EPOS" );
        auto    hP2_Ntru_Mult_EPOS      =   uLoadHistograms<0,TH1F> ( insFile_MCPr_YL_pp_05_EPOS,    "hP2_Ntru_Mult",   "hP2_Ntru_Mult_EPOS" );
        auto    h1D_Ntru_Mult_P6      =   uLoadHistograms<0,TH1F> ( insFile_MCPr_YL_pp_05_PYTHIA6,    "h1D_Ntru_Mult",   "h1D_Ntru_Mult_P6" );
        auto    h2D_Ntru_Mult_P6      =   uLoadHistograms<0,TH1F> ( insFile_MCPr_YL_pp_05_PYTHIA6,    "h2D_Ntru_Mult",   "h2D_Ntru_Mult_P6" );
        auto    hR1_Ntru_Mult_P6      =   uLoadHistograms<0,TH1F> ( insFile_MCPr_YL_pp_05_PYTHIA6,    "hR1_Ntru_Mult",   "hR1_Ntru_Mult_P6" );
        auto    hR2_Ntru_Mult_P6      =   uLoadHistograms<0,TH1F> ( insFile_MCPr_YL_pp_05_PYTHIA6,    "hR2_Ntru_Mult",   "hR2_Ntru_Mult_P6" );
        auto    hP1_Ntru_Mult_P6      =   uLoadHistograms<0,TH1F> ( insFile_MCPr_YL_pp_05_PYTHIA6,    "hP1_Ntru_Mult",   "hP1_Ntru_Mult_P6" );
        auto    hP2_Ntru_Mult_P6      =   uLoadHistograms<0,TH1F> ( insFile_MCPr_YL_pp_05_PYTHIA6,    "hP2_Ntru_Mult",   "hP2_Ntru_Mult_P6" );
        auto    h1D_Ntru_Mult_P8      =   uLoadHistograms<0,TH1F> ( insFile_MCPr_YL_pp_05_PYTHIA8,    "h1D_Ntru_Mult",   "h1D_Ntru_Mult_P8" );
        auto    h2D_Ntru_Mult_P8     =   uLoadHistograms<0,TH1F> ( insFile_MCPr_YL_pp_05_PYTHIA8,    "h2D_Ntru_Mult",   "h2D_Ntru_Mult_P8" );
        auto    hR1_Ntru_Mult_P8      =   uLoadHistograms<0,TH1F> ( insFile_MCPr_YL_pp_05_PYTHIA8,    "hR1_Ntru_Mult",   "hR1_Ntru_Mult_P8" );
        auto    hR2_Ntru_Mult_P8      =   uLoadHistograms<0,TH1F> ( insFile_MCPr_YL_pp_05_PYTHIA8,    "hR2_Ntru_Mult",   "hR2_Ntru_Mult_P8" );
        auto    hP1_Ntru_Mult_P8      =   uLoadHistograms<0,TH1F> ( insFile_MCPr_YL_pp_05_PYTHIA8,    "hP1_Ntru_Mult",   "hP1_Ntru_Mult_P8" );
        auto    hP2_Ntru_Mult_P8      =   uLoadHistograms<0,TH1F> ( insFile_MCPr_YL_pp_05_PYTHIA8,    "hP2_Ntru_Mult",   "hP2_Ntru_Mult_P8" );
        h1D_Ntru_Mult_EPOS    ->SetLineColor( uGetColor(4) );
        h2D_Ntru_Mult_EPOS    ->SetLineColor( uGetColor(4) );
        hR1_Ntru_Mult_EPOS    ->SetLineColor( uGetColor(4) );
        hR2_Ntru_Mult_EPOS    ->SetLineColor( uGetColor(4) );
        hP1_Ntru_Mult_EPOS    ->SetLineColor( uGetColor(4) );
        hP2_Ntru_Mult_EPOS    ->SetLineColor( uGetColor(4) );
        h1D_Ntru_Mult_P6    ->SetLineColor( uGetColor(5) );
        h2D_Ntru_Mult_P6    ->SetLineColor( uGetColor(5) );
        hR1_Ntru_Mult_P6    ->SetLineColor( uGetColor(5) );
        hR2_Ntru_Mult_P6    ->SetLineColor( uGetColor(5) );
        hP1_Ntru_Mult_P6    ->SetLineColor( uGetColor(5) );
        hP2_Ntru_Mult_P6    ->SetLineColor( uGetColor(5) );
        h1D_Ntru_Mult_P8       ->SetLineColor( uGetColor(6) );
        h2D_Ntru_Mult_P8       ->SetLineColor( uGetColor(6) );
        hR1_Ntru_Mult_P8       ->SetLineColor( uGetColor(6) );
        hR2_Ntru_Mult_P8       ->SetLineColor( uGetColor(6) );
        hP1_Ntru_Mult_P8       ->SetLineColor( uGetColor(6) );
        hP2_Ntru_Mult_P8       ->SetLineColor( uGetColor(6) );
        //
        //  --- Output directory
        TString kPlotDirectory          =   Form(kDIR_FinalResultsPlots,(TString("Multiplicity")).Data());
        gROOT   ->  ProcessLine(Form(".! mkdir -p %s",kPlotDirectory.Data()));
        //
        g1DResults_pp   ->  SetMarkerStyle( uGetMarker(2) );
        g2DResults_pp   ->  SetMarkerStyle( uGetMarker(2) );
        gR1Results_pp   ->  SetMarkerStyle( uGetMarker(2) );
        gR2Results_pp   ->  SetMarkerStyle( uGetMarker(2) );
        gP1Results_pp   ->  SetMarkerStyle( uGetMarker(2) );
        gP2Results_pp   ->  SetMarkerStyle( uGetMarker(2) );
        g1DResults_Pb   ->  SetMarkerStyle( uGetMarker(3) );
        g2DResults_Pb   ->  SetMarkerStyle( uGetMarker(3) );
        gR1Results_Pb   ->  SetMarkerStyle( uGetMarker(3) );
        gR2Results_Pb   ->  SetMarkerStyle( uGetMarker(3) );
        gP1Results_Pb   ->  SetMarkerStyle( uGetMarker(3) );
        gP2Results_Pb   ->  SetMarkerStyle( uGetMarker(3) );
        g1DResults_p3   ->  SetMarkerStyle( uGetMarker(4) );
        g2DResults_p3   ->  SetMarkerStyle( uGetMarker(4) );
        gR1Results_p3   ->  SetMarkerStyle( uGetMarker(4) );
        gR2Results_p3   ->  SetMarkerStyle( uGetMarker(4) );
        gP1Results_p3   ->  SetMarkerStyle( uGetMarker(4) );
        gP2Results_p3   ->  SetMarkerStyle( uGetMarker(4) );
        g1DResults_pp   ->  SetMarkerColor( uGetColor(2) );
        g2DResults_pp   ->  SetMarkerColor( uGetColor(2) );
        gR1Results_pp   ->  SetMarkerColor( uGetColor(2) );
        gR2Results_pp   ->  SetMarkerColor( uGetColor(2) );
        gP1Results_pp   ->  SetMarkerColor( uGetColor(2) );
        gP2Results_pp   ->  SetMarkerColor( uGetColor(2) );
        g1DResults_Pb   ->  SetMarkerColor( uGetColor(3) );
        g2DResults_Pb   ->  SetMarkerColor( uGetColor(3) );
        gR1Results_Pb   ->  SetMarkerColor( uGetColor(3) );
        gR2Results_Pb   ->  SetMarkerColor( uGetColor(3) );
        gP1Results_Pb   ->  SetMarkerColor( uGetColor(3) );
        gP2Results_Pb   ->  SetMarkerColor( uGetColor(3) );
        g1DResults_p3   ->  SetMarkerColor( uGetColor(4) );
        g2DResults_p3   ->  SetMarkerColor( uGetColor(4) );
        gR1Results_p3   ->  SetMarkerColor( uGetColor(4) );
        gR2Results_p3   ->  SetMarkerColor( uGetColor(4) );
        gP1Results_p3   ->  SetMarkerColor( uGetColor(4) );
        gP2Results_p3   ->  SetMarkerColor( uGetColor(4) );
        g1DResults_pp   ->  SetLineColor( uGetColor(2) );
        g2DResults_pp   ->  SetLineColor( uGetColor(2) );
        gR1Results_pp   ->  SetLineColor( uGetColor(2) );
        gR2Results_pp   ->  SetLineColor( uGetColor(2) );
        gP1Results_pp   ->  SetLineColor( uGetColor(2) );
        gP2Results_pp   ->  SetLineColor( uGetColor(2) );
        g1DResults_Pb   ->  SetLineColor( uGetColor(3) );
        g2DResults_Pb   ->  SetLineColor( uGetColor(3) );
        gR1Results_Pb   ->  SetLineColor( uGetColor(3) );
        gR2Results_Pb   ->  SetLineColor( uGetColor(3) );
        gP1Results_Pb   ->  SetLineColor( uGetColor(3) );
        gP2Results_Pb   ->  SetLineColor( uGetColor(3) );
        g1DResults_p3   ->  SetLineColor( uGetColor(4) );
        g2DResults_p3   ->  SetLineColor( uGetColor(4) );
        gR1Results_p3   ->  SetLineColor( uGetColor(4) );
        gR2Results_p3   ->  SetLineColor( uGetColor(4) );
        gP1Results_p3   ->  SetLineColor( uGetColor(4) );
        gP2Results_p3   ->  SetLineColor( uGetColor(4) );
        //
        TLegend*    cDrawLegend =   new TLegend( 0.15,0.55,0.45,0.9);
        cDrawLegend     ->  SetFillColorAlpha(0.,0.);
        cDrawLegend     ->  SetLineColorAlpha(0.,0.);
        cDrawLegend     ->  AddEntry( g1DResults_pp, "pp @5.02TeV", "EP" );
        cDrawLegend     ->  AddEntry( g1DResults_p3, "pp @13.0TeV", "EP" );
        cDrawLegend     ->  AddEntry( g1DResults_Pb, "pPb @5.02TeV", "EP" );
       // cDrawLegend     ->  AddEntry( h1D_Ntru_Mult_EPOS, "EPOS", "L" );
       // cDrawLegend     ->  AddEntry( h1D_Ntru_Mult_P6, "PYTHIA6", "L" );
       // cDrawLegend     ->  AddEntry( h1D_Ntru_Mult_P8, "PYTHIA8", "L" );
        //
        for ( Int_t iTer = 0; iTer < g1DResults_pp->GetN(); iTer++ )    {
            auto    kXMean__    =   g1DResults_pp->GetPointX(iTer);
            auto    k1DYield    =   g1DResults_pp->GetPointY(iTer);
            auto    k2DYield    =   g2DResults_pp->GetPointY(iTer);
            auto    k1DError    =   g1DResults_pp->GetEY()[iTer];
            auto    k2DError    =   g2DResults_pp->GetEY()[iTer];
            auto    k1DRlErr    =   k1DError/k1DYield;
            auto    k2DRlErr    =   k2DError/k2DYield;
            gR1Results_pp->SetPoint         ( iTer, kXMean__,    k2DYield/k1DYield );
            gR1Results_pp->SetPointError    ( iTer, 0,           (k2DYield/k1DYield)*SquareSum( {k2DRlErr,k1DRlErr} ) );
            gR2Results_pp->SetPoint         ( iTer, kXMean__,    k2DYield/(k1DYield*k1DYield) );
            gR2Results_pp->SetPointError    ( iTer, 0,           (k2DYield/(k1DYield*k1DYield))*SquareSum( {k2DRlErr,k1DRlErr,k1DRlErr} ) );
            gP1Results_pp->SetPoint         ( iTer, kXMean__,    fSigmaPhiValue(k1DYield,k2DYield) );
            gP1Results_pp->SetPointError    ( iTer, 0,           fSigmaPhiError(k1DYield,k2DYield,k1DError,k2DError) );
            gP2Results_pp->SetPoint         ( iTer, kXMean__,    fGammaPhiValue(k1DYield,k2DYield) );
            gP2Results_pp->SetPointError    ( iTer, 0,           fGammaPhiError(k1DYield,k2DYield,k1DError,k2DError) );
        }
        for ( Int_t iTer = 0; iTer < g1DResults_p3->GetN(); iTer++ )    {
            auto    kXMean__    =   g1DResults_p3->GetPointX(iTer);
            auto    k1DYield    =   g1DResults_p3->GetPointY(iTer);
            auto    k2DYield    =   g2DResults_p3->GetPointY(iTer);
            auto    k1DError    =   g1DResults_p3->GetEY()[iTer];
            auto    k2DError    =   g2DResults_p3->GetEY()[iTer];
            auto    k1DRlErr    =   k1DError/k1DYield;
            auto    k2DRlErr    =   k2DError/k2DYield;
            gR1Results_p3->SetPoint         ( iTer, kXMean__,    k2DYield/k1DYield );
            gR1Results_p3->SetPointError    ( iTer, 0,           (k2DYield/k1DYield)*SquareSum( {k2DRlErr,k1DRlErr} ) );
            gR2Results_p3->SetPoint         ( iTer, kXMean__,    k2DYield/(k1DYield*k1DYield) );
            gR2Results_p3->SetPointError    ( iTer, 0,           (k2DYield/(k1DYield*k1DYield))*SquareSum( {k2DRlErr,k1DRlErr,k1DRlErr} ) );
            gP1Results_p3->SetPoint         ( iTer, kXMean__,    fSigmaPhiValue(k1DYield,k2DYield) );
            gP1Results_p3->SetPointError    ( iTer, 0,           fSigmaPhiError(k1DYield,k2DYield,k1DError,k2DError) );
            gP2Results_p3->SetPoint         ( iTer, kXMean__,    fGammaPhiValue(k1DYield,k2DYield) );
            gP2Results_p3->SetPointError    ( iTer, 0,           fGammaPhiError(k1DYield,k2DYield,k1DError,k2DError) );
        }
        for ( Int_t iTer = 0; iTer < g1DResults_Pb->GetN(); iTer++ )    {
            auto    kXMean__    =   g1DResults_Pb->GetPointX(iTer);
            auto    k1DYield    =   g1DResults_Pb->GetPointY(iTer);
            auto    k2DYield    =   g2DResults_Pb->GetPointY(iTer);
            auto    k1DError    =   g1DResults_Pb->GetEY()[iTer];
            auto    k2DError    =   g2DResults_Pb->GetEY()[iTer];
            auto    k1DRlErr    =   k1DError/k1DYield;
            auto    k2DRlErr    =   k2DError/k2DYield;
            gR1Results_Pb->SetPoint         ( iTer, kXMean__,    k2DYield/k1DYield );
            gR1Results_Pb->SetPointError    ( iTer, 0,           (k2DYield/k1DYield)*SquareSum( {k2DRlErr,k1DRlErr} ) );
            gR2Results_Pb->SetPoint         ( iTer, kXMean__,    k2DYield/(k1DYield*k1DYield) );
            gR2Results_Pb->SetPointError    ( iTer, 0,           (k2DYield/(k1DYield*k1DYield))*SquareSum( {k2DRlErr,k1DRlErr,k1DRlErr} ) );
            gP1Results_Pb->SetPoint         ( iTer, kXMean__,    fSigmaPhiValue(k1DYield,k2DYield) );
            gP1Results_Pb->SetPointError    ( iTer, 0,           fSigmaPhiError(k1DYield,k2DYield,k1DError,k2DError) );
            gP2Results_Pb->SetPoint         ( iTer, kXMean__,    fGammaPhiValue(k1DYield,k2DYield) );
            gP2Results_Pb->SetPointError    ( iTer, 0,           fGammaPhiError(k1DYield,k2DYield,k1DError,k2DError) );
        }
        //
        //g1DResults_pp->Fit("pol1");
        g1DResults_p3->Fit("pol1");
        //auto    pol1_pp = g1DResults_pp->GetFunction("pol1");
        auto    pol1_p3 = g1DResults_p3->GetFunction("pol1");
        //pol1_pp ->  SetLineStyle(3);
        pol1_p3 ->  SetLineStyle(3);
        //pol1_pp ->  SetLineColor(uGetColor(2));
        pol1_p3 ->  SetLineColor(uGetColor(4));
        //
        //g2DResults_pp->Fit("pol2");
        g2DResults_p3->Fit("pol2");
        //auto    pol2_pp = g2DResults_pp->GetFunction("pol2");
        auto    pol2_p3 = g2DResults_p3->GetFunction("pol2");
        //pol2_pp ->  SetLineStyle(3);
        pol2_p3 ->  SetLineStyle(3);
        //pol2_pp ->  SetLineColor(uGetColor(2));
        pol2_p3 ->  SetLineColor(uGetColor(4));
        //
        TF1*    fR1 =   new TF1("fR1","([0]+[1]*x+[2]*x*x)/(([3]+[4]*x))",0,100);
        TF1*    fR2 =   new TF1("fR2","([0]+[1]*x+[2]*x*x)/(([3]+[4]*x)*([3]+[4]*x))",0,100);
        TF1*    fP1 =   new TF1("fP1","2*([0]+[1]*x+[2]*x*x)+([3]+[4]*x)-(([3]+[4]*x)*([3]+[4]*x))",0,100);
        TF1*    fP2 =   new TF1("fP2","2*([0]+[1]*x+[2]*x*x)/(([3]+[4]*x))-([3]+[4]*x)",0,100);
        //
        fR1 ->  SetParameter(0,pol2_p3->GetParameter(0));
        fR1 ->  SetParameter(1,pol2_p3->GetParameter(1));
        fR1 ->  SetParameter(2,pol2_p3->GetParameter(2));
        fR1 ->  SetParameter(3,pol1_p3->GetParameter(0));
        fR1 ->  SetParameter(4,pol1_p3->GetParameter(1));
        fR1 ->  SetLineStyle(3);
        fR1 ->  SetLineColor(uGetColor(4));
        //
        fR2 ->  SetParameter(0,pol2_p3->GetParameter(0));
        fR2 ->  SetParameter(1,pol2_p3->GetParameter(1));
        fR2 ->  SetParameter(2,pol2_p3->GetParameter(2));
        fR2 ->  SetParameter(3,pol1_p3->GetParameter(0));
        fR2 ->  SetParameter(4,pol1_p3->GetParameter(1));
        fR2 ->  SetLineStyle(3);
        fR2 ->  SetLineColor(uGetColor(4));
        //
        fP1 ->  SetParameter(0,pol2_p3->GetParameter(0));
        fP1 ->  SetParameter(1,pol2_p3->GetParameter(1));
        fP1 ->  SetParameter(2,pol2_p3->GetParameter(2));
        fP1 ->  SetParameter(3,pol1_p3->GetParameter(0));
        fP1 ->  SetParameter(4,pol1_p3->GetParameter(1));
        fP1 ->  SetLineStyle(3);
        fP1 ->  SetLineColor(uGetColor(4));
        //
        fP2 ->  SetParameter(0,pol2_p3->GetParameter(0));
        fP2 ->  SetParameter(1,pol2_p3->GetParameter(1));
        fP2 ->  SetParameter(2,pol2_p3->GetParameter(2));
        fP2 ->  SetParameter(3,pol1_p3->GetParameter(0));
        fP2 ->  SetParameter(4,pol1_p3->GetParameter(1));
        fP2 ->  SetLineStyle(3);
        fP2 ->  SetLineColor(uGetColor(4));
        //
        TH1F*   hDrawRange  = new TH1F( "hDrawRange", "hDrawRange", 1, 0, 45 );
        hDrawRange->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
        //
        TCanvas*    cDrawResults    =   new TCanvas( "cDrawResults", "cDrawResults", 1500, 1500 );
        cDrawResults    ->  Divide(3,2);
        gStyle->SetOptStat(0000);
        //
        cDrawResults    ->  cd(1);
        hDrawRange->SetTitle( "#LT Y_{#phi} #GT" );
        hDrawRange->SetMinimum(0);
        hDrawRange->SetMaximum(0.4);
        hDrawRange->DrawCopy("");
        g1DResults_pp->Draw("EP SAME");
        g1DResults_p3->Draw("EP SAME");
        g1DResults_Pb->Draw("EP SAME");
        cDrawLegend->Draw("SAME");
        //h1D_Ntru_Mult_EPOS->Draw("SAME");
        //h1D_Ntru_Mult_P6->Draw("SAME");
        //h1D_Ntru_Mult_P8->Draw("SAME");
        //
        cDrawResults    ->  cd(2);
        hDrawRange->SetTitle( "#LT Y_{#phi#phi} #GT / #LT Y_{#phi} #GT" );
        hDrawRange->SetMinimum(0);
        hDrawRange->SetMaximum(0.2);
        hDrawRange->DrawCopy("");
        gR1Results_pp->Draw("EP SAME");
        gR1Results_p3->Draw("EP SAME");
        gR1Results_Pb->Draw("EP SAME");
        fR1->  Draw("same");
        //hR1_Ntru_Mult_EPOS->Draw("SAME");
        //hR1_Ntru_Mult_P6->Draw("SAME");
        //hR1_Ntru_Mult_P8->Draw("SAME");
        //
        cDrawResults    ->  cd(3);
        hDrawRange->SetTitle( "#sigma_{#phi}" );
        hDrawRange->SetMinimum(0);
        hDrawRange->SetMaximum(0.3);
        hDrawRange->DrawCopy("");
        gP1Results_pp->Draw("EP SAME");
        gP1Results_p3->Draw("EP SAME");
        gP1Results_Pb->Draw("EP SAME");
        fP1->Draw("SAME");
        //hP1_Ntru_Mult_EPOS->Draw("SAME");
        //hP1_Ntru_Mult_P6->Draw("SAME");
        //hP1_Ntru_Mult_P8->Draw("SAME");
        //
        cDrawResults    ->  cd(4);
        hDrawRange->SetTitle( "#LT Y_{#phi#phi} #GT" );
        hDrawRange->SetMinimum(0);
        hDrawRange->SetMaximum(0.05);
        hDrawRange->DrawCopy("");
        g2DResults_pp->Draw("EP SAME");
        g2DResults_p3->Draw("EP SAME");
        g2DResults_Pb->Draw("EP SAME");
        //h2D_Ntru_Mult_EPOS->Draw("SAME");
        //h2D_Ntru_Mult_P6->Draw("SAME");
        //h2D_Ntru_Mult_P8->Draw("SAME");
        //
        cDrawResults    ->  cd(5);
        hDrawRange->SetTitle( "#LT Y_{#phi#phi} #GT / #LT Y_{#phi} #GT^{2}" );
        hDrawRange->SetBinContent(1,0.5);
        hDrawRange->SetMinimum(0);
        hDrawRange->SetMaximum(2.5);
        hDrawRange->DrawCopy("");
        gR2Results_pp->Draw("EP SAME");
        gR2Results_p3->Draw("EP SAME");
        gR2Results_Pb->Draw("EP SAME");
        fR2->Draw("SAME");
       // hR2_Ntru_Mult_EPOS->Draw("SAME");
        //hR2_Ntru_Mult_P6->Draw("SAME");
        //hR2_Ntru_Mult_P8->Draw("SAME");
        //
        cDrawResults    ->  cd(6);
        hDrawRange->SetTitle( "#gamma_{#phi}" );
        hDrawRange->SetBinContent(1,0.);
        hDrawRange->SetMinimum(-0.25);
        hDrawRange->SetMaximum(+0.25);
        hDrawRange->DrawCopy("");
        gP2Results_pp->Draw("EP SAME");
        gP2Results_p3->Draw("EP SAME");
        gP2Results_Pb->Draw("EP SAME");
        fP2->Draw("same");
        //hP2_Ntru_Mult_EPOS->Draw("SAME");
        //hP2_Ntru_Mult_P6->Draw("SAME");
        //hP2_Ntru_Mult_P8->Draw("SAME");
        //
        cDrawResults    ->  SaveAs( kPlotDirectory+ TString("MultiplicityDependence.pdf"));
        //
        delete cDrawResults;
    }*/
}
