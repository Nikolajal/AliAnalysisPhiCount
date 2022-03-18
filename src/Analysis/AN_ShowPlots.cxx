#include "../../inc/AliAnalysisPhiPair.h"

void AN_ShowPlots ( TString fOption = "all", TString kFolder = "_p_p__7TeV", std::vector<TString> kComparisonMCTagList = { "PYTHIA6_7TeV", "PYTHIA8_7TeV_CR", "PYTHIA8_7TeV_ROPES" }, std::vector<TString> kComparisonLegendList = { "Pythia6 Perugia 2011", "Pythia8 CR", "Pythia8 Ropes" } )    {
    // --- --- --- --- --- --- --- SET-UP --- --- --- --- --- --- --- --- --- --- ---
    //
    //  --- INFO on Set-up variables
    fChooseOption(fOption);
    //
    //  --- Setting the input datastructure
    fSetAllBins();
    //
    //  --- Setting the style
    SetStyle();
    //
    // --- YIELD ANALYSIS
    if ( kDoYield ) {
        //
        //  --- Load Data Files
        TFile*  insFile_Data_YL = new TFile ( Form(kASigExtp_FitCheckRst,(TString("Yield")+kFolder).Data()) );
        TFile*  insFile_Syst_YL = new TFile ( Form("%s/FullSystematics.root",Form(kAnalysis_Systemt_Dir,  (TString("Yield")+kFolder).Data())) );
        //
        //  --- Load Data Histograms
        auto    hXD_Nyld_stat   = uLoadHistograms<0,TH1F> ( insFile_Data_YL, "hXD_Nyld_stat" );
        auto    hXD_Nyld_syst   = uLoadHistograms<0,TH1F> ( insFile_Data_YL, "hXD_Nyld_syst" );
        auto    hXD_Nfqs_stat   = uLoadHistograms<0,TH1F> ( insFile_Data_YL, "hXD_Nfqs_stat" );
        auto    hXD_Nfqs_syst   = uLoadHistograms<0,TH1F> ( insFile_Data_YL, "hXD_Nfqs_syst" );
        auto    h1D_Nres_stat   = uLoadHistograms<0,TH1F> ( insFile_Data_YL, "h1D_Nres_stat" );
        auto    h1D_Nres_syst   = uLoadHistograms<0,TH1F> ( insFile_Data_YL, "h1D_Nres_syst" );
        auto    h2D_Nres_stat   = uLoadHistograms<1,TH1F> ( insFile_Data_YL, "h2D_Nres_stat_stat_%i" );
        auto    h2D_Nres_syst   = uLoadHistograms<1,TH1F> ( insFile_Data_YL, "h2D_Nres_syst_syst_%i" );
        auto    h2D_MeanPT_stat = uLoadHistograms<0,TH1F> ( insFile_Data_YL, "h2D_MeanPT_stat" );
        auto    h2D_MeanPT_syst = uLoadHistograms<0,TH1F> ( insFile_Data_YL, "h2D_MeanPT_syst" );
        //
        // --- MC Elements
        std::vector<TFile*> insFile_MntC_YL;
        std::vector<TH1F*> hMPT_NTru;
        std::vector<TH1F*> hMPT_NTru_Fin;
        std::vector<TH1F*> h1D_Ntru_Fin;
        std::vector<TH2F*> h2D_Ntru_Fin;
        std::vector<TH2F*> h2D_Ntru_Fin_Nrm;
        std::vector<TH1F*> h1D_Ntru;
        std::vector<TH2F*> h2D_Ntru;
        std::vector<TH1F*> hFullQuantities;
        TLegend*           lMCLegend_1D;
        TLegend*           lMCLegend_2D;
        TLegend*           lMCLegend_PT;
        //
        if ( kComparisonMCTagList.size() ) {
            // --- Load MC Files
            for ( auto kCurrent_MC_File : kComparisonMCTagList ) insFile_MntC_YL.push_back( new TFile   ( Form(kProduction_MC_Ofl,(TString("Yield")+kFolder).Data(),kCurrent_MC_File.Data()) ) );
            //
            // --- Build Legend
            lMCLegend_1D    = new TLegend( 0.55, 0.70, 0.85, 0.55 );
            lMCLegend_2D    = new TLegend( 0.55, 0.65, 0.85, 0.50 );
            lMCLegend_PT    = new TLegend( 0.65, 0.03, 0.88, 0.18 );
            //
            //  --- Load MC Histograms
            auto iClr = 0;
            for ( auto kFile : insFile_MntC_YL )   {
                iClr++;
                //  --- Production Histograms
                //  --- --- 1D
                h1D_Ntru                    .push_back( (TH1F*)((( kFile )->Get("h1D_Ntru"))        -> Clone(Form("h1D_Ntru_%s",kFile->GetName())) ));
                h1D_Ntru.at( iClr-1 )       -> SetLineWidth( 2 );
                h1D_Ntru.at( iClr-1 )       -> SetLineColor( uGetColor( iClr ) );
                h1D_Ntru.at( iClr-1 )       -> SetFillColorAlpha( uGetColor( iClr ), 0.33 );
                h1D_Ntru_Fin                .push_back( (TH1F*)((( kFile )->Get("h1D_Ntru_Fin"))    -> Clone(Form("h1D_Ntru_Fin_%s",kFile->GetName())) ));
                h1D_Ntru_Fin.at( iClr-1 )   -> SetLineWidth( 2 );
                h1D_Ntru_Fin.at( iClr-1 )   -> SetLineColor( uGetColor( iClr ) );
                h1D_Ntru_Fin.at( iClr-1 )   -> SetFillColorAlpha( uGetColor( iClr ), 0.33 );
                //
                //  --- --- 2D
                h2D_Ntru                    .push_back( (TH2F*)((( kFile )->Get("h2D_Ntru"))        -> Clone(Form("h2D_Ntru_%s",kFile->GetName())) ));
                h2D_Ntru.at( iClr-1 )       -> SetLineWidth( 2 );
                h2D_Ntru.at( iClr-1 )       -> SetLineColor( uGetColor( iClr ) );
                h2D_Ntru.at( iClr-1 )       -> SetFillColorAlpha( uGetColor( iClr ), 0.33 );
                h2D_Ntru_Fin                .push_back( (TH2F*)((( kFile )->Get("h2D_Ntru_Fin"))    -> Clone(Form("h2D_Ntru_Fin_%s",kFile->GetName())) ));
                h2D_Ntru_Fin.at( iClr-1 )   -> SetLineWidth( 2 );
                h2D_Ntru_Fin.at( iClr-1 )   -> SetLineColor( uGetColor( iClr ) );
                h2D_Ntru_Fin.at( iClr-1 )   -> SetFillColorAlpha( uGetColor( iClr ), 0.33 );
                h2D_Ntru_Fin_Nrm            .push_back( (TH2F*)((( kFile )->Get("h2D_Ntru_Fin_Nrm"))-> Clone(Form("h2D_Ntru_Fin_Nrm_%s",kFile->GetName())) ));
                h2D_Ntru_Fin_Nrm.at( iClr-1 )-> SetLineWidth( 2 );
                h2D_Ntru_Fin_Nrm.at( iClr-1 )-> SetLineColor( uGetColor( iClr ) );
                h2D_Ntru_Fin_Nrm.at( iClr-1 )-> SetFillColorAlpha( uGetColor( iClr ), 0.33 );
                //
                //  --- --- MPT
                hMPT_NTru                   .push_back( (TH1F*)((( kFile )->Get("hMPT_NTru"))       -> Clone(Form("hMPT_NTru_%s",kFile->GetName())) ));
                hMPT_NTru.at( iClr-1 )      -> SetLineWidth( 2 );
                hMPT_NTru.at( iClr-1 )      -> SetLineColor( uGetColor( iClr ) );
                hMPT_NTru.at( iClr-1 )      -> SetFillColorAlpha( uGetColor( iClr ), 0.33 );
                hMPT_NTru_Fin               .push_back( (TH1F*)((( kFile )->Get("hMPT_NTru_Fin"))   -> Clone(Form("hMPT_NTru_Fin_%s",kFile->GetName())) ));
                hMPT_NTru_Fin.at( iClr-1 )  -> SetLineWidth( 2 );
                hMPT_NTru_Fin.at( iClr-1 )  -> SetLineColor( uGetColor( iClr ) );
                hMPT_NTru_Fin.at( iClr-1 )  -> SetFillColorAlpha( uGetColor( iClr ), 0.33 );
                //
                //  --- --- Other
                hFullQuantities             .push_back( (TH1F*)((( kFile )->Get("hFullQuantities")) -> Clone(Form("hFullQuantities_%s",kFile->GetName())) ));
                hFullQuantities.at( iClr-1 )-> SetLineWidth( 2 );
                hFullQuantities.at( iClr-1 )-> SetLineColor( uGetColor( iClr ) );
                hFullQuantities.at( iClr-1 )-> SetFillColorAlpha( uGetColor( iClr ), 0.33 );
                //
                lMCLegend_1D                -> AddEntry( hMPT_NTru.at( iClr-1 ), kComparisonLegendList.at( iClr-1 ), "L" );
                lMCLegend_2D                -> AddEntry( hMPT_NTru.at( iClr-1 ), kComparisonLegendList.at( iClr-1 ), "L" );
                lMCLegend_PT                -> AddEntry( hMPT_NTru.at( iClr-1 ), kComparisonLegendList.at( iClr-1 ), "L" );
            }
        }
        //
        //  --- Output directory
        TString kPlotDirectory          =   Form(kDIR_ShowPlots,(TString("Yield")+kFolder).Data());
        gROOT   ->  ProcessLine(Form(".! mkdir -p %s",kPlotDirectory.Data()));
        //
        //  --- Yield Plots
        //  --- --- Batch Mode
        gROOT   ->  SetBatch( kTRUE );
        //
        //  --- --- Final quantities
        uSetHisto( hXD_Nfqs_stat,       "SPT STAT 1D" );
        uSetHisto( hXD_Nfqs_syst,       "SPT SYST 1D" );
        //
        TLegend*    lLegend         = new TLegend( 0.43, 0.03, 0.63, 0.20 );
        lLegend     ->  SetLineColorAlpha(0.,0.);
        lLegend     ->  SetFillColorAlpha(0.,0.);
        lLegend     ->  AddEntry    (hXD_Nfqs_stat,"Measured","P");
        lLegend     ->  AddEntry    (hXD_Nfqs_stat,"Stat. Uncert.","EL");
        lLegend     ->  AddEntry    (hXD_Nfqs_syst,"Syst. Uncert.","F");
        //
        TCanvas*    cDrawResults    = new TCanvas( "cDrawResults", "cDrawResults", 1200, 1000 );
        gPad                    -> SetLogy();
        //
        //  --- Upper Plot
        TPad*   kUpperPlot  =   new TPad("kUpperPlot_RT", "kUpperPlot", 0, 0.3, 1, 1.0);
        gStyle      -> SetOptStat(0);
        kUpperPlot  -> SetBottomMargin(0);
        kUpperPlot  -> SetFillColorAlpha( 0., 0. );
        kUpperPlot  -> SetLogy();
        kUpperPlot  -> cd();
        //
        hXD_Nfqs_stat           -> Draw("PE2 SAME");
        hXD_Nfqs_syst           -> Draw("E1 SAME");
        if ( kComparisonMCTagList.size() ) {
            for ( auto kCurrent_FQ : hFullQuantities ) kCurrent_FQ-> Draw("SAME EP");
        }
        lLegend                 -> Draw("SAME");
        lMCLegend_PT            -> Draw("SAME");
        //
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.18, 0.83,"ALICE Preliminary");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        uLatex->DrawLatexNDC(0.18, 0.78,"pp #sqrt{#it{s}}= 7 TeV");
        uLatex->DrawLatexNDC(0.18, 0.73,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
        //
        //  --- Lower Plot
        TPad*   kLowerPlot  =   new TPad("kLowerPlot_RT", "kLowerPlot", 0, 0.0, 1, 0.3);
        kLowerPlot  -> SetGridy();
        kLowerPlot  -> SetTopMargin(0);
        kLowerPlot  -> SetFillColorAlpha( 0., 0. );
        kLowerPlot  -> cd();
        //
        if ( kComparisonMCTagList.size() ) {
            auto kFullUncertainties = uSumErrors( hXD_Nfqs_stat, hXD_Nfqs_syst );
            auto kUtilityPlot       = uScale(kFullUncertainties,0.,-2.);
            kUtilityPlot    -> GetYaxis() -> SetTitle("Ratio Model / Data ");
            kUtilityPlot    -> GetYaxis() -> SetTitleSize(0.10);
            kUtilityPlot    -> GetYaxis() -> SetTitleOffset(0.3);
            kUtilityPlot    -> SetMaximum( 1.9 );
            kUtilityPlot    -> SetMinimum( 0.1 );
            kUtilityPlot    -> Draw();
            for ( auto kCurrent_FQ : hFullQuantities ) {
                auto kCurrentHist   =   (TH1F*)kCurrent_FQ -> Clone(Form("tmp_%s",kCurrent_FQ->GetName()));
                kCurrentHist    -> Divide( kFullUncertainties );
                kCurrentHist    -> DrawCopy("SAME E2");
                kCurrentHist    -> SetFillColorAlpha( 0.,0. );
                kCurrentHist    -> DrawCopy("SAME HIST");
            }
        }
        //
        cDrawResults    -> cd();
        kUpperPlot      -> Draw();
        kLowerPlot      -> Draw();
        cDrawResults    -> SaveAs( kPlotDirectory + TString("FullQuantitites.pdf") );
        cDrawResults    -> SaveAs( kPlotDirectory + TString("FullQuantitites.eps") );
        delete lLegend;
        delete kUpperPlot;
        delete kLowerPlot;
        delete cDrawResults;
        //
        //  --- --- Full 1D Spectrum //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        TH1F*   h1D_Nres_stat_full      = new TH1F( "h1D_Nres_stat_full", "h1D_Nres_stat_full", nBinPT1D+1, fArrPT1D_Comp );
        TH1F*   h1D_Nres_syst_full      = new TH1F( "h1D_Nres_syst_full", "h1D_Nres_syst_full", nBinPT1D+1, fArrPT1D_Comp );
        //
        h1D_Nres_stat_full  -> SetBinContent( 1, hXD_Nyld_stat->GetBinContent(4) / (fArrPT1D_Comp[1]-fArrPT1D_Comp[0]) );
        h1D_Nres_syst_full  -> SetBinContent( 1, hXD_Nyld_stat->GetBinContent(4) / (fArrPT1D_Comp[1]-fArrPT1D_Comp[0]) );
        h1D_Nres_stat_full  -> SetBinError  ( 1, hXD_Nyld_stat->GetBinError  (4) / (fArrPT1D_Comp[1]-fArrPT1D_Comp[0]) );
        h1D_Nres_syst_full  -> SetBinError  ( 1, hXD_Nyld_syst->GetBinError  (4) / (fArrPT1D_Comp[1]-fArrPT1D_Comp[0]) );
        h1D_Nres_syst_full  -> SetMaximum( h1D_Nres_syst->GetMaximum()*2.5 );
        h1D_Nres_syst_full  -> SetMinimum( h1D_Nres_syst->GetMinimum()*0.5 );
        uSetHisto( h1D_Nres_stat,       "SPT STAT 1D" );
        uSetHisto( h1D_Nres_syst,       "SPT SYST 1D" );
        uSetHisto( h1D_Nres_stat_full,  "SPT STAT 1D" );
        uSetHisto( h1D_Nres_syst_full,  "SPT SYST 1D" );
        h1D_Nres_stat_full  -> SetMarkerStyle( uGetMarker(8) );
        h1D_Nres_syst_full  -> SetMarkerStyle( uGetMarker(8) );
        //
        lLegend     = new TLegend( 0.18, 0.03, 0.38, 0.20 );
        lLegend     ->  SetLineColorAlpha(0.,0.);
        lLegend     ->  SetFillColorAlpha(0.,0.);
        lLegend     ->  AddEntry    (h1D_Nres_stat,"Measured","P");
        lLegend     ->  AddEntry    (h1D_Nres_stat_full,"Extrapolated","P");
        lLegend     ->  AddEntry    (h1D_Nres_stat,"Stat. Uncert.","EL");
        lLegend     ->  AddEntry    (h1D_Nres_syst,"Syst. Uncert.","F");
        //
        cDrawResults    = new TCanvas( "cDrawResults", "cDrawResults", 1200, 1000 );
        //
        //  --- Upper Plot
        kUpperPlot  =   new TPad("kUpperPlot_1D", "kUpperPlot", 0, 0.3, 1, 1.0);
        gStyle      -> SetOptStat(0);
        kUpperPlot  -> SetBottomMargin(0);
        kUpperPlot  -> SetFillColorAlpha( 0., 0. );
        kUpperPlot  -> SetLogy();
        kUpperPlot  -> cd();
        //
        h1D_Nres_syst_full  -> Draw("PE2 SAME");
        h1D_Nres_stat_full  -> Draw("E1 SAME");
        h1D_Nres_syst       -> Draw("PE2 SAME");
        h1D_Nres_stat       -> Draw("E1 SAME");
        if ( kComparisonMCTagList.size() ) {
            for ( auto kCurrent_1D : h1D_Ntru_Fin ) {
                kCurrent_1D-> DrawCopy("SAME E3");
                kCurrent_1D-> SetFillColorAlpha( 0.,0. );
                kCurrent_1D-> DrawCopy("SAME HIST L");
            }
        }
        lLegend             -> Draw("SAME");
        lMCLegend_1D        -> Draw("SAME");
        //
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.55, 0.83,"ALICE Preliminary");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        uLatex->DrawLatexNDC(0.55, 0.78,"pp #sqrt{#it{s}}= 7 TeV");
        uLatex->DrawLatexNDC(0.55, 0.73,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
        //
        //  --- Lower Plot
        kLowerPlot  =   new TPad("kLowerPlot_1D", "kLowerPlot", 0, 0.0, 1, 0.3);
        kLowerPlot  -> SetGridy();
        kLowerPlot  -> SetTopMargin(0);
        kLowerPlot  -> SetFillColorAlpha( 0., 0. );
        kLowerPlot  -> cd();
        //
        if ( kComparisonMCTagList.size() ) {
            auto    h1D_Nres_stat_full_utl  = (TH1F*)h1D_Nres_stat_full->Clone("h1D_Nres_stat_full_utl");
            auto    h1D_Nres_syst_full_utl  = (TH1F*)h1D_Nres_stat_full->Clone("h1D_Nres_syst_full_utl");
            for ( Int_t iBin = 1; iBin <= h1D_Nres_stat->GetNbinsX(); iBin++ ) {
                h1D_Nres_stat_full_utl->SetBinContent(iBin+1,h1D_Nres_stat->GetBinContent(iBin));
                h1D_Nres_syst_full_utl->SetBinContent(iBin+1,h1D_Nres_stat->GetBinContent(iBin));
                h1D_Nres_stat_full_utl->SetBinError(iBin+1,h1D_Nres_stat->GetBinError(iBin));
                h1D_Nres_syst_full_utl->SetBinError(iBin+1,h1D_Nres_syst->GetBinError(iBin));
            }
            auto kFullUncertainties = uSumErrors( h1D_Nres_stat_full_utl, h1D_Nres_syst_full_utl );
            auto kUtilityPlot       = uScale(kFullUncertainties,0.,-2.);
            kUtilityPlot    -> GetYaxis() -> SetTitle("Ratio Model / Data ");
            kUtilityPlot    -> GetYaxis() -> SetTitleSize(0.10);
            kUtilityPlot    -> GetYaxis() -> SetTitleOffset(0.3);
            kUtilityPlot    -> SetMaximum( 1.9 );
            kUtilityPlot    -> SetMinimum( 0.1 );
            kUtilityPlot    -> Draw();
            for ( auto kCurrent_1D : h1D_Ntru ) {
                auto kCurrentHist   =   (TH1F*)kCurrent_1D -> Clone(Form("tmp_%s",kCurrent_1D->GetName()));
                kCurrentHist    -> Divide( kFullUncertainties );
                kCurrentHist    -> DrawCopy("SAME E2");
                kCurrentHist    -> SetFillColorAlpha( 0.,0. );
                kCurrentHist    -> DrawCopy("SAME HIST");
            }
        }
        //
        cDrawResults    -> cd();
        kUpperPlot      -> Draw();
        kLowerPlot      -> Draw();
        cDrawResults    -> SaveAs( kPlotDirectory + TString("Yield_1D.pdf") );
        cDrawResults    -> SaveAs( kPlotDirectory + TString("Yield_1D.eps") );
        delete cDrawResults;    //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        //
        for ( Int_t iPT2D = 1; iPT2D < h2D_Nres_stat.size(); iPT2D++ ) {
            //
            //  --- --- Conditional 2D Spectrum //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            TH1F*   h2D_Nres_stat_full      = new TH1F( Form("h2D_Nres_stat_full_%i",iPT2D), Form("h2D_Nres_stat_full_%i",iPT2D), nBinPT2D+1, fArrPT2D_Comp );
            TH1F*   h2D_Nres_syst_full      = new TH1F( Form("h2D_Nres_syst_full_%i",iPT2D), Form("h2D_Nres_syst_full_%i",iPT2D), nBinPT2D+1, fArrPT2D_Comp );
            //
            h2D_Nres_stat_full  -> SetBinContent( 1, h2D_Nres_stat.at(0)->GetBinContent(iPT2D) );
            h2D_Nres_syst_full  -> SetBinContent( 1, h2D_Nres_stat.at(0)->GetBinContent(iPT2D) );
            h2D_Nres_stat_full  -> SetBinError  ( 1, h2D_Nres_stat.at(0)->GetBinError  (iPT2D) );
            h2D_Nres_syst_full  -> SetBinError  ( 1, h2D_Nres_syst.at(0)->GetBinError  (iPT2D) );
            h2D_Nres_syst_full  -> SetMaximum( h2D_Nres_syst.at(iPT2D)->GetMaximum()*2.5 );
            h2D_Nres_syst_full  -> SetMinimum( h2D_Nres_syst.at(iPT2D)->GetMinimum()*0.5 );
            uSetHisto( h2D_Nres_stat_full,      " SPT STAT 12D " );
            uSetHisto( h2D_Nres_syst_full,      " SPT SYST 12D " );
            uSetHisto( h2D_Nres_stat.at(iPT2D), " SPT STAT 12D " );
            uSetHisto( h2D_Nres_syst.at(iPT2D), " SPT SYST 12D " );
            h2D_Nres_stat_full  -> SetMarkerStyle( uGetMarker(8) );
            h2D_Nres_syst_full  -> SetMarkerStyle( uGetMarker(8) );
            //
            TCanvas*    cDrawResults    = new TCanvas( "cDrawResults", "cDrawResults", 1200, 1000 );
            //
            //  --- Upper Plot
            TPad*   kUpperPlot  =   new TPad("kUpperPlot_2D", "kUpperPlot", 0, 0.3, 1, 1.0);
            gStyle      -> SetOptStat(0);
            kUpperPlot  -> SetBottomMargin(0);
            kUpperPlot  -> SetFillColorAlpha( 0., 0. );
            kUpperPlot  -> SetLogy();
            kUpperPlot  -> cd();
            //
            h2D_Nres_syst_full      -> Draw("SAME PE2");
            h2D_Nres_stat_full      -> Draw("SAME  E1");
            h2D_Nres_syst.at(iPT2D) -> Draw("SAME PE2");
            h2D_Nres_stat.at(iPT2D) -> Draw("SAME  E1");
            if ( kComparisonMCTagList.size() ) {
                for ( auto kCurrent_2D : h2D_Ntru_Fin_Nrm ) {
                    auto kCurrentHist  =   kCurrent_2D -> ProjectionY("tmp", iPT2D+1 , iPT2D+1 );
                    kCurrentHist   ->  SetLineWidth( 2 );
                    kCurrentHist   ->  DrawCopy("SAME E3");
                    kCurrentHist   ->  SetFillColorAlpha(0.,0.);
                    kCurrentHist   ->  DrawCopy("SAME HIST L");
                }
            }
            lLegend                 -> Draw("SAME");
            lMCLegend_2D            -> Draw("SAME");
            //
            uLatex->SetTextFont(60);
            uLatex->SetTextSize(0.05);
            uLatex->DrawLatexNDC(0.55, 0.83,"ALICE Preliminary");
            uLatex->SetTextFont(42);
            uLatex->SetTextSize(0.04);
            uLatex->DrawLatexNDC(0.55, 0.78,"pp #sqrt{#it{s}}= 7 TeV");
            uLatex->DrawLatexNDC(0.55, 0.73,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
            uLatex->DrawLatexNDC(0.55, 0.68,Form("#it{p}_{T,#phi_{2}} [%.2f;%.2f] (GeV/#it{c})",fArrPT2D_Comp[iPT2D],fArrPT2D_Comp[iPT2D+1]));
            //
            //  --- Lower Plot
            TPad*   kLowerPlot  =   new TPad("kLowerPlot_2D", "kLowerPlot", 0, 0.0, 1, 0.3);
            kLowerPlot  -> SetGridy();
            kLowerPlot  -> SetTopMargin(0);
            kLowerPlot  -> SetFillColorAlpha( 0., 0. );
            kLowerPlot  -> cd();
            //
            if ( kComparisonMCTagList.size() ) {
                auto    h2D_Nres_stat_full_utl  = (TH1F*)h2D_Nres_stat_full->Clone("h2D_Nres_stat_full_utl");
                auto    h2D_Nres_syst_full_utl  = (TH1F*)h2D_Nres_stat_full->Clone("h2D_Nres_syst_full_utl");
                for ( Int_t iBin = 1; iBin <= h2D_Nres_stat.at(iPT2D)->GetNbinsX(); iBin++ ) {
                    h2D_Nres_stat_full_utl->SetBinContent(iBin+1,h2D_Nres_stat.at(iPT2D)->GetBinContent(iBin));
                    h2D_Nres_syst_full_utl->SetBinContent(iBin+1,h2D_Nres_stat.at(iPT2D)->GetBinContent(iBin));
                    h2D_Nres_stat_full_utl->SetBinError(iBin+1,h2D_Nres_stat.at(iPT2D)->GetBinError(iBin));
                    h2D_Nres_syst_full_utl->SetBinError(iBin+1,h2D_Nres_syst.at(iPT2D)->GetBinError(iBin));
                }
                auto kFullUncertainties = uSumErrors( h2D_Nres_stat_full_utl, h2D_Nres_syst_full_utl );
                auto kUtilityPlot       = uScale(kFullUncertainties,0.,-2.);
                kUtilityPlot    -> GetYaxis() -> SetTitle("Ratio Model / Data ");
                kUtilityPlot    -> GetYaxis() -> SetTitleSize(0.10);
                kUtilityPlot    -> GetYaxis() -> SetTitleOffset(0.3);
                kUtilityPlot    -> SetMaximum( 1.9 );
                kUtilityPlot    -> SetMinimum( 0.1 );
                kUtilityPlot    -> Draw();
                for ( auto kCurrent_2D : h2D_Ntru ) {
                    auto kCurrentHist   =   kCurrent_2D -> ProjectionY(Form("tmp_%s",kCurrent_2D->GetName()), iPT2D+1 , iPT2D+1 );
                    kCurrentHist    -> Divide( kFullUncertainties );
                    kCurrentHist    -> DrawCopy("SAME E2");
                    kCurrentHist    -> SetFillColorAlpha( 0.,0. );
                    kCurrentHist    -> DrawCopy("SAME HIST");
                }
            }
            //
            cDrawResults    -> cd();
            kUpperPlot      -> Draw();
            kLowerPlot      -> Draw();
            cDrawResults    -> SaveAs( kPlotDirectory + TString(Form("Yield_2D_%i.pdf",iPT2D)) );
            cDrawResults    -> SaveAs( kPlotDirectory + TString(Form("Yield_2D_%i.eps",iPT2D)) );
            delete kUpperPlot;
            delete kLowerPlot;
            delete cDrawResults;    //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        }
        //
        //  --- --- Mean pT Spectrum //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        TH1F*   h2D_MeanPT_stat_full    = new TH1F( "h2D_MeanPT_stat_full", "h2D_MeanPT_stat_full", nBinPT2D+1, fArrPT2D_Comp );
        TH1F*   h2D_MeanPT_syst_full    = new TH1F( "h2D_MeanPT_syst_full", "h2D_MeanPT_syst_full", nBinPT2D+1, fArrPT2D_Comp );
        //
        h2D_MeanPT_stat_full  -> SetBinContent( 1, h2D_MeanPT_stat->GetBinContent(1) );
        h2D_MeanPT_syst_full  -> SetBinContent( 1, h2D_MeanPT_stat->GetBinContent(1) );
        h2D_MeanPT_stat_full  -> SetBinError  ( 1, h2D_MeanPT_stat->GetBinError  (1) );
        h2D_MeanPT_syst_full  -> SetBinError  ( 1, h2D_MeanPT_syst->GetBinError  (1) );
        h2D_MeanPT_syst       -> SetMaximum( h2D_MeanPT_syst->GetMaximum()*1.25 );
        h2D_MeanPT_syst       -> SetMinimum( h2D_MeanPT_syst->GetMinimum()*0.40 );
        uSetHisto( h2D_MeanPT_stat_full,  " MPT STAT 12D " );
        uSetHisto( h2D_MeanPT_syst_full,  " MPT SYST 12D " );
        uSetHisto( h2D_MeanPT_stat, " MPT STAT 12D " );
        uSetHisto( h2D_MeanPT_syst, " MPT SYST 12D " );
        h2D_MeanPT_stat_full    -> SetMarkerStyle( uGetMarker(8) );
        h2D_MeanPT_syst_full    -> SetMarkerStyle( uGetMarker(8) );
        //
                cDrawResults    = new TCanvas( "cDrawResults", "cDrawResults", 1200, 1000 );
        //
        //  --- Upper Plot
        kUpperPlot  =   new TPad("kUpperPlot_MPT", "kUpperPlot", 0, 0.3, 1, 1.0);
        gStyle      -> SetOptStat(0);
        kUpperPlot  -> SetBottomMargin(0);
        kUpperPlot  -> SetFillColorAlpha( 0., 0. );
        kUpperPlot  -> cd();
        //
        h2D_MeanPT_syst         -> Draw("SAME PE2");
        h2D_MeanPT_stat         -> Draw("SAME  E1");
        h2D_MeanPT_stat_full    -> Draw("SAME  E1");
        if ( kComparisonMCTagList.size() ) {
            for ( auto kCurrent_MPT : hMPT_NTru_Fin )   {
                kCurrent_MPT-> DrawCopy("SAME E3");
                kCurrent_MPT-> SetFillColorAlpha(0.,0.);
                kCurrent_MPT-> DrawCopy("SAME HIST L");
            }
        }
        lLegend                 -> Draw("SAME");
        lMCLegend_PT            -> Draw("SAME");
        //
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.38, 0.15,"ALICE Preliminary");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        uLatex->DrawLatexNDC(0.38, 0.10,"pp #sqrt{#it{s}}= 7 TeV");
        uLatex->DrawLatexNDC(0.38, 0.05,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
        //
        //  --- Lower Plot
        kLowerPlot  =   new TPad("kLowerPlot_MPT", "kLowerPlot", 0, 0.0, 1, 0.3);
        kLowerPlot  -> SetGridy();
        kLowerPlot  -> SetTopMargin(0);
        kLowerPlot  -> SetFillColorAlpha( 0., 0. );
        kLowerPlot  -> cd();
        //
        if ( kComparisonMCTagList.size() ) {
            auto    h2D_MeanPT_stat_full_utl  = (TH1F*)h2D_MeanPT_stat_full->Clone("h2D_MeanPT_stat_full_utl");
            auto    h2D_MeanPT_syst_full_utl  = (TH1F*)h2D_MeanPT_stat_full->Clone("h2D_MeanPT_syst_full_utl");
            for ( Int_t iBin = 1; iBin <= h2D_MeanPT_stat->GetNbinsX(); iBin++ ) {
                h2D_MeanPT_stat_full_utl->SetBinContent(iBin+1,h2D_MeanPT_stat->GetBinContent(iBin));
                h2D_MeanPT_syst_full_utl->SetBinContent(iBin+1,h2D_MeanPT_stat->GetBinContent(iBin));
                h2D_MeanPT_stat_full_utl->SetBinError(iBin+1,h2D_MeanPT_stat->GetBinError(iBin));
                h2D_MeanPT_syst_full_utl->SetBinError(iBin+1,h2D_MeanPT_syst->GetBinError(iBin));
            }
            auto kFullUncertainties = uSumErrors( h2D_MeanPT_stat_full_utl, h2D_MeanPT_syst_full_utl );
            auto kUtilityPlot       = uScale(kFullUncertainties,0.,-2.);
            kUtilityPlot    -> GetYaxis() -> SetTitle("Ratio Model / Data ");
            kUtilityPlot    -> GetYaxis() -> SetTitleSize(0.10);
            kUtilityPlot    -> GetYaxis() -> SetTitleOffset(0.3);
            kUtilityPlot    -> SetMaximum( 1.9 );
            kUtilityPlot    -> SetMinimum( 0.1 );
            kUtilityPlot    -> Draw();
            for ( auto kCurrent_MPT : hMPT_NTru ) {
                auto kCurrentHist   =   (TH1F*)kCurrent_MPT -> Clone(Form("tmp_%s",kCurrent_MPT->GetName()));
                kCurrentHist    -> Divide( kFullUncertainties );
                kCurrentHist   -> DrawCopy("SAME E2");
                kCurrentHist   -> SetFillColorAlpha(0.,0.);
                kCurrentHist   -> DrawCopy("SAME HIST");
            }
        }
        //
        cDrawResults    -> cd();
        kUpperPlot      -> Draw();
        kLowerPlot      -> Draw();
        cDrawResults    -> SaveAs( kPlotDirectory + TString(Form("MeanPT_2D_%i.pdf",0)) );
        cDrawResults    -> SaveAs( kPlotDirectory + TString(Form("MeanPT_2D_%i.eps",0)) );
        delete kUpperPlot;
        delete kLowerPlot;
        delete cDrawResults; //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        //
        //  --- --- Down to 0 pT conditional 2D Spectrum //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        TH1F*   h2D_Nres_stat_full      = new TH1F( "h2D_Nres_stat_full_0", "h2D_Nres_stat_full_0", nBinPT2D+1, fArrPT2D_Comp );
        TH1F*   h2D_Nres_syst_full      = new TH1F( "h2D_Nres_syst_full_0", "h2D_Nres_syst_full_0", nBinPT2D+1, fArrPT2D_Comp );
        //
        h2D_Nres_stat_full  -> SetBinContent( 1, hXD_Nyld_stat->GetBinContent(5) / (fArrPT2D_Comp[1]-fArrPT2D_Comp[0]) );
        h2D_Nres_syst_full  -> SetBinContent( 1, hXD_Nyld_stat->GetBinContent(5) / (fArrPT2D_Comp[1]-fArrPT2D_Comp[0]) );
        h2D_Nres_stat_full  -> SetBinError  ( 1, hXD_Nyld_stat->GetBinError  (5) / (fArrPT2D_Comp[1]-fArrPT2D_Comp[0]) );
        h2D_Nres_syst_full  -> SetBinError  ( 1, hXD_Nyld_syst->GetBinError  (5) / (fArrPT2D_Comp[1]-fArrPT2D_Comp[0]) );
        h2D_Nres_syst_full  -> SetMaximum( h2D_Nres_syst.at(0)->GetMaximum()*2.5 );
        h2D_Nres_syst_full  -> SetMinimum( h2D_Nres_syst.at(0)->GetMinimum()*0.5 );
        uSetHisto( h2D_Nres_stat_full,  " SPT STAT 12D " );
        uSetHisto( h2D_Nres_syst_full,  " SPT SYST 12D " );
        uSetHisto( h2D_Nres_stat.at(0), " SPT STAT 12D " );
        uSetHisto( h2D_Nres_syst.at(0), " SPT SYST 12D " );
        h2D_Nres_stat.at(0)  -> SetMarkerStyle( uGetMarker(8) );
        h2D_Nres_syst.at(0)  -> SetMarkerStyle( uGetMarker(8) );
        h2D_Nres_stat_full  -> SetMarkerStyle( uGetMarker(5) );
        h2D_Nres_syst_full  -> SetMarkerStyle( uGetMarker(5) );
        //
                    lLegend         = new TLegend( 0.18, 0.03, 0.38, 0.20 );
        lLegend     ->  SetLineColorAlpha(0.,0.);
        lLegend     ->  SetFillColorAlpha(0.,0.);
        lLegend     ->  AddEntry    (h2D_Nres_stat.at(0),"Extrapolated","P");
        lLegend     ->  AddEntry    (h2D_Nres_stat_full,"Double Extrapolated","P");
        lLegend     ->  AddEntry    (h2D_Nres_stat.at(0),"Stat. Uncert.","EL");
        lLegend     ->  AddEntry    (h2D_Nres_syst.at(0),"Syst. Uncert.","F");
        //
                cDrawResults    = new TCanvas( "cDrawResults", "cDrawResults", 1200, 1000 );
        //
        //  --- Upper Plot
        kUpperPlot  =   new TPad("kUpperPlot", "kUpperPlot", 0, 0.3, 1, 1.0);
        gStyle      -> SetOptStat(0);
        kUpperPlot  -> SetBottomMargin(0);
        kUpperPlot  -> SetFillColorAlpha( 0., 0. );
        kUpperPlot  -> SetLogy();
        kUpperPlot  -> cd();
        //
        h2D_Nres_syst_full      -> Draw("SAME PE2");
        h2D_Nres_stat_full      -> Draw("SAME  E1");
        h2D_Nres_syst.at(0)     -> Draw("SAME PE2");
        h2D_Nres_stat.at(0)     -> Draw("SAME  E1");
        if ( kComparisonMCTagList.size() ) {
            for ( auto kCurrent_2D : h2D_Ntru_Fin_Nrm ) {
                auto kCurrentHist  =   kCurrent_2D -> ProjectionY("tmp", 1 , 1 );
                kCurrentHist-> DrawCopy("SAME E3");
                kCurrentHist-> SetFillColorAlpha(0.,0.);
                kCurrentHist-> DrawCopy("SAME HIST L");
            }
        }
        lLegend                 -> Draw("SAME");
        lMCLegend_2D            -> Draw("SAME");
        //
        uLatex  -> SetTextFont  (60);
        uLatex  -> SetTextSize  (0.05);
        uLatex  -> DrawLatexNDC ( 0.55, 0.83, "ALICE Preliminary" );
        uLatex  -> SetTextFont  (42);
        uLatex  -> SetTextSize  (0.04);
        uLatex  -> DrawLatexNDC (0.55, 0.78, "pp #sqrt{#it{s}}= 7 TeV" );
        uLatex  -> DrawLatexNDC (0.55, 0.73, "#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5" );
        uLatex  -> DrawLatexNDC (0.55, 0.68, Form("#it{p}_{T,#phi_{2}} [%.2f;%.2f] (GeV/#it{c})",fArrPT2D_Comp[0],fArrPT2D_Comp[1]) );
        //
        //  --- Lower Plot
        kLowerPlot  =   new TPad("kLowerPlot", "kLowerPlot", 0, 0.0, 1, 0.3);
        kLowerPlot  -> SetGridy();
        kLowerPlot  -> SetTopMargin(0);
        kLowerPlot  -> SetFillColorAlpha( 0., 0. );
        kLowerPlot  -> cd();
        //
        if ( kComparisonMCTagList.size() ) {
            auto    h2D_Nres_stat_full_utl  = (TH1F*)h2D_Nres_stat_full->Clone("h2D_Nres_stat_full_utl");
            auto    h2D_Nres_syst_full_utl  = (TH1F*)h2D_Nres_stat_full->Clone("h2D_Nres_syst_full_utl");
            for ( Int_t iBin = 1; iBin <= h2D_Nres_stat.at(0)->GetNbinsX(); iBin++ ) {
                h2D_Nres_stat_full_utl->SetBinContent(iBin+1,h2D_Nres_stat.at(0)->GetBinContent(iBin));
                h2D_Nres_syst_full_utl->SetBinContent(iBin+1,h2D_Nres_stat.at(0)->GetBinContent(iBin));
                h2D_Nres_stat_full_utl->SetBinError(iBin+1,h2D_Nres_stat.at(0)->GetBinError(iBin));
                h2D_Nres_syst_full_utl->SetBinError(iBin+1,h2D_Nres_syst.at(0)->GetBinError(iBin));
            }
            auto kFullUncertainties = uSumErrors( h2D_Nres_stat_full_utl, h2D_Nres_syst_full_utl );
            auto kUtilityPlot       = uScale(kFullUncertainties,0.,-2.);
            kUtilityPlot    -> GetYaxis() -> SetTitle("Ratio Model / Data ");
            kUtilityPlot    -> GetYaxis() -> SetTitleSize(0.10);
            kUtilityPlot    -> GetYaxis() -> SetTitleOffset(0.3);
            kUtilityPlot    -> SetMaximum( 1.9 );
            kUtilityPlot    -> SetMinimum( 0.1 );
            kUtilityPlot    -> Draw();
            for ( auto kCurrent_2D : h2D_Ntru ) {
                auto kCurrentHist   =   kCurrent_2D -> ProjectionY(Form("tmp_%s",kCurrent_2D->GetName()), 1, 1 );
                kCurrentHist    -> Divide( kFullUncertainties );
                kCurrentHist    -> DrawCopy("SAME E2");
                kCurrentHist    -> SetFillColorAlpha( 0.,0. );
                kCurrentHist    -> DrawCopy("SAME HIST");
            }
        }
        //
        cDrawResults    -> cd();
        kUpperPlot      -> Draw();
        kLowerPlot      -> Draw();
        cDrawResults    -> SaveAs( kPlotDirectory + TString(Form("Yield_2D_%i.pdf",0)) );
        cDrawResults    -> SaveAs( kPlotDirectory + TString(Form("Yield_2D_%i.eps",0)) );
        delete kUpperPlot;
        delete kLowerPlot;
        delete cDrawResults; //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        //
        gROOT   ->  SetBatch( kFALSE );
    }
}
