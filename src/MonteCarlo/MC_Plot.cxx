#include "../../inc/AliAnalysisPhiPair.h"

void
MC_Plot
( std::vector<TString> fTags = {"Pythia6"}, std::vector<TString> fLegendTags = {"Pythia6"} , TString fOption = "yield", TString kFolder = "_p_p__7TeV"  ) {
    // --- --- --- --- --- --- --- SET-UP --- --- --- --- --- --- --- --- --- --- ---
    //
    // --- Retrieving Data
    TFile*      insFile_Data_YL    =   new TFile   ( Form(kASigExtp_FitCheckRst,(TString("Yield")+kFolder).Data()) );
    cout << Form(kASigExtp_FitCheckRst,(TString("Yield")+kFolder).Data()) << endl;
    //
    auto        h1D_Nres_stat   =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "h1D_Nres_stat",    "h1D_Nres_stat" );
    auto        h1D_Nres_syst   =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "h1D_Nres_syst",    "h1D_Nres_syst" );
    auto        h2D_Nres_stat   =   uLoadHistograms<0,TH2F> ( insFile_Data_YL,  "h2D_Nres_stat",    "h2D_Nres_stat" );
    auto        h2D_Nres_syst   =   uLoadHistograms<0,TH2F> ( insFile_Data_YL,  "h2D_Nres_syst",    "h2D_Nres_syst" );
    //
    // --- Retrieving MC Comparisons
    std::vector<TFile*> insFile_MntC_YL;
    for ( auto kCurrent_MC_File : fTags ) insFile_MntC_YL.push_back( new TFile   ( Form(kProduction_MC_Ofl,(TString("Yield")+kFolder).Data(),kCurrent_MC_File.Data()) ) );
    //
    std::vector<TH1F*> h1D_Ntru_Fin;
    std::vector<TH2F*> h2D_Ntru_Fin;
    std::vector<TH1F*> h1D_Ntru;
    std::vector<TH2F*> h2D_Ntru;
    std::vector<TH1F*> hFullQuantities;
    for ( auto kFile : insFile_MntC_YL )   {
        // --- Single Histograms
        h1D_Ntru_Fin        .push_back( (TH1F*)((( kFile )->Get("h1D_Ntru_Fin"))    ->Clone(Form("h1D_Ntru_Fin_%s",kFile->GetName())) ));
        h2D_Ntru_Fin        .push_back( (TH2F*)((( kFile )->Get("h2D_Ntru_Fin"))    ->Clone(Form("h2D_Ntru_Fin_%s",kFile->GetName())) ));
        h1D_Ntru            .push_back( (TH1F*)((( kFile )->Get("h1D_Ntru"))        ->Clone(Form("h1D_Ntru_%s",kFile->GetName())) ));
        h2D_Ntru            .push_back( (TH2F*)((( kFile )->Get("h2D_Ntru"))        ->Clone(Form("h2D_Ntru_%s",kFile->GetName())) ));
        hFullQuantities     .push_back( (TH1F*)((( kFile )->Get("hFullQuantities")) ->Clone(Form("hFullQuantities_%s",kFile->GetName())) ));
        //
    }
    //
    // --- --- --- --- --- --- --- OUTPUT --- --- --- --- --- --- --- --- --- --- ---
    //
    fSetAllBins();
    auto    cMC_Compare_1D  =uPlotYieldWithMC( h1D_Nres_stat, h1D_Nres_syst, h1D_Ntru, h1D_Ntru_Fin, fLegendTags, "SPT 1D" );
    uLatex->SetTextFont(60);
    uLatex->SetTextSize(0.05);
    uLatex->DrawLatexNDC(0.63, 0.83,"ALICE");
    uLatex->SetTextFont(42);
    uLatex->SetTextSize(0.04);
    uLatex->DrawLatexNDC(0.63, 0.78,"pp #sqrt{#it{s}}= 7 TeV");
    uLatex->DrawLatexNDC(0.63, 0.73,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
    cMC_Compare_1D  ->  SaveAs("/Users/nikolajal/alice/AliAnalysisPhiCount/result/Yield_p_p__7TeV/MC_Production/Plots/1D/Test1D.pdf");
    //
    /*
    std::vector<std::vector<TH1D*>> k2DTru;
    std::vector<std::vector<TH1D*>> k2DTru_Fin;
    auto    kUtilitySlice   =   h2D_Ntru_Fin.at(0)->ProjectionX("_tmp_utl_",1,1);
    for ( Int_t iBin = 1; iBin <= h2D_Nres_stat.at(0)->GetNbinsX()+1; iBin++ )  {
        auto    kFirstBin       =   kUtilitySlice->FindBin(h2D_Nres_stat.at(0)->GetBinLowEdge(iBin-1));
        auto    kLastBin        =   kUtilitySlice->FindBin(h2D_Nres_stat.at(0)->GetBinLowEdge(iBin));
        auto iTer = -1;
        std::vector<TH1D*> hTruUtility;
        std::vector<TH1D*> hTruUtility_Fin;
        for ( auto kCurrentTru : h2D_Ntru )  {
            iTer++;
            auto    k2DSlice_Fin    =   h2D_Ntru_Fin.at(0)->ProjectionX(Form("tmp_prj_iBin_%i",iTer),kFirstBin,kLastBin);
            auto    k2DSlice        =   h2D_Ntru    .at(0)->ProjectionX(Form("tmp_prj_iBin_%i_tru",iTer),iBin,iBin);
            hTruUtility             .push_back(k2DSlice);
            hTruUtility_Fin         .push_back(k2DSlice_Fin);
        }
        k2DTru.push_back(hTruUtility);
        k2DTru_Fin.push_back(hTruUtility_Fin);
    }
    //
    for ( Int_t iBin = 1; iBin <= h2D_Nres_stat.at(0)->GetNbinsX()+1; iBin++ )  {
        gROOT->SetBatch(true);
        auto    cMC_Compare_2D  =   uPlotYieldWithMC<TH1F,TH1D>( h2D_Nres_stat.at(iBin-1), h2D_Nres_syst.at(iBin-1), k2DTru.at(iBin-1), k2DTru.at(iBin-1), fLegendTags, "SPT 12D" );
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.63, 0.83,"ALICE");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        uLatex->DrawLatexNDC(0.63, 0.78,"pp #sqrt{#it{s}}= 7 TeV");
        uLatex->DrawLatexNDC(0.63, 0.73,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
        uLatex->DrawLatexNDC(0.28, 0.83,Form("#it{p}_{T,#phi_{2}} #in [%.2f-%.2f] GeV/#it{c}",fArrPT2D_Comp[iBin-1],fArrPT2D_Comp[iBin]));
        cMC_Compare_2D  ->  SaveAs(Form("/Users/nikolajal/alice/AliAnalysisPhiCount/result/Yield_p_p__7TeV/MC_Production/Plots/2D/Test2D_%i.pdf",iBin));
        delete cMC_Compare_2D;
        gROOT->SetBatch(false);
    }
     */
    //
    SetStyle();
    TH1F*   hPlotResults_Full_stat = new TH1F( "hPlotResults_Full_stat", "", 6, 0.5, 6.5 );
    TH1F*   hPlotResults_Full_syst = new TH1F( "hPlotResults_Full_syst", "", 6, 0.5, 6.5 );
    uSetHisto( hPlotResults_Full_stat, "SPT 12D STAT");
    uSetHisto( hPlotResults_Full_syst, "SPT 12D SYST");
    auto hPlotResults_Ratio_stat = (TH1F*)(hPlotResults_Full_stat->Clone("hPlotResults_Ratio_stat"));
    auto hPlotResults_Ratio_syst = (TH1F*)(hPlotResults_Full_syst->Clone("hPlotResults_Ratio_syst"));
    
    hPlotResults_Full_stat->GetXaxis()->SetBinLabel(1,"#frac{dN_{#phi}}{dy}");
    hPlotResults_Full_stat->GetXaxis()->SetBinLabel(2,"#frac{dN_{#phi#phi}}{dy}");
    hPlotResults_Full_stat->GetXaxis()->SetBinLabel(3,"#frac{#LT Y_{#phi#phi} #GT }{ #LT Y_{#phi} #GT }");
    hPlotResults_Full_stat->GetXaxis()->SetBinLabel(4,"#frac{#LT Y_{#phi#phi} #GT }{ #LT Y_{#phi} #GT^{2} }");
    hPlotResults_Full_stat->GetXaxis()->SetBinLabel(5,"#sigma_{#phi}");
    hPlotResults_Full_stat->GetXaxis()->SetBinLabel(6,"#gamma_{#phi}");
    
    hPlotResults_Ratio_syst->GetXaxis()->SetBinLabel(1,"#frac{dN_{#phi}}{dy}");
    hPlotResults_Ratio_syst->GetXaxis()->SetBinLabel(2,"#frac{dN_{#phi#phi}}{dy}");
    hPlotResults_Ratio_syst->GetXaxis()->SetBinLabel(3,"#frac{#LT Y_{#phi#phi} #GT }{ #LT Y_{#phi} #GT }");
    hPlotResults_Ratio_syst->GetXaxis()->SetBinLabel(4,"#frac{#LT Y_{#phi#phi} #GT }{ #LT Y_{#phi} #GT^{2} }");
    hPlotResults_Ratio_syst->GetXaxis()->SetBinLabel(5,"#sigma_{#phi}");
    hPlotResults_Ratio_syst->GetXaxis()->SetBinLabel(6,"#gamma_{#phi}");
    
    hPlotResults_Full_stat->SetBinContent(1,0.0318336);
    hPlotResults_Full_syst->SetBinContent(1,0.0318336);
    hPlotResults_Ratio_stat->SetBinContent(1,1);
    hPlotResults_Ratio_syst->SetBinContent(1,1);
    hPlotResults_Full_stat->SetBinError(1,.000174066);
    hPlotResults_Full_syst->SetBinError(1,.0026183);
    hPlotResults_Ratio_stat->SetBinError(1,.000174066/.0318336);
    hPlotResults_Ratio_syst->SetBinError(1,.002618300/.0318336);
    
    hPlotResults_Full_stat->SetBinContent(2,0.00143711);
    hPlotResults_Full_syst->SetBinContent(2,0.00143711);
    hPlotResults_Ratio_stat->SetBinContent(2,1);
    hPlotResults_Ratio_syst->SetBinContent(2,1);
    hPlotResults_Full_stat->SetBinError(2,.000269857);
    hPlotResults_Full_syst->SetBinError(2,.000226415);
    hPlotResults_Ratio_stat->SetBinError(2,.000269857/0.00143711);
    hPlotResults_Ratio_syst->SetBinError(2,.000226415/0.00143711);
    
    hPlotResults_Full_stat->SetBinContent(3,0.0451443);
    hPlotResults_Full_syst->SetBinContent(3,0.0451443);
    hPlotResults_Ratio_stat->SetBinContent(3,1);
    hPlotResults_Ratio_syst->SetBinContent(3,1);
    hPlotResults_Full_stat->SetBinError(3,.00848068);
    hPlotResults_Full_syst->SetBinError(3,.00475983);
    hPlotResults_Ratio_stat->SetBinError(3,.00848068/0.0451443);
    hPlotResults_Ratio_syst->SetBinError(3,.00475983/0.0451443);
    
    hPlotResults_Full_stat->SetBinContent(4,1.40491);
    hPlotResults_Full_syst->SetBinContent(4,1.40491);
    hPlotResults_Ratio_stat->SetBinContent(4,1);
    hPlotResults_Ratio_syst->SetBinContent(4,1);
    hPlotResults_Full_stat->SetBinError(4,0.291734);
    hPlotResults_Full_syst->SetBinError(4,0.117907);
    hPlotResults_Ratio_stat->SetBinError(4,0.291734/1.40491);
    hPlotResults_Ratio_syst->SetBinError(4,0.117907/1.40491);
    
    hPlotResults_Full_stat->SetBinContent(5,0.0339882);
    hPlotResults_Full_syst->SetBinContent(5,0.0339882);
    hPlotResults_Ratio_stat->SetBinContent(5,1);
    hPlotResults_Ratio_syst->SetBinContent(5,1);
    hPlotResults_Full_stat->SetBinError(5,0.000620131);
    hPlotResults_Full_syst->SetBinError(5,0.00304672);
    hPlotResults_Ratio_stat->SetBinError(5,0.000620131/0.0339882);
    hPlotResults_Ratio_syst->SetBinError(5,0.003046720/0.0339882);
    
    hPlotResults_Full_stat->SetBinContent(6,0.0581334);
    hPlotResults_Full_syst->SetBinContent(6,0.0581334);
    hPlotResults_Ratio_stat->SetBinContent(6,1);
    hPlotResults_Ratio_syst->SetBinContent(6,1);
    hPlotResults_Full_stat->SetBinError(6,0.0192981);
    hPlotResults_Full_syst->SetBinError(6,0.00898665);
    hPlotResults_Ratio_stat->SetBinError(6,0.0192981/0.0581334);
    hPlotResults_Ratio_syst->SetBinError(6,0.00898665/0.0581334);
    
    hPlotResults_Full_stat->GetXaxis()->SetTitle("");
    hPlotResults_Full_stat->GetYaxis()->SetTitle("");
    hPlotResults_Full_syst->GetXaxis()->SetTitle("");
    hPlotResults_Full_syst->GetYaxis()->SetTitle("");
    hPlotResults_Ratio_stat->GetXaxis()->SetTitle("");
    hPlotResults_Ratio_stat->GetYaxis()->SetTitle("");
    hPlotResults_Ratio_syst->GetXaxis()->SetTitle("");
    hPlotResults_Ratio_syst->GetYaxis()->SetTitle("");
    
    hPlotResults_Full_syst->SetMinimum(1.e-4);
    hPlotResults_Full_syst->SetMaximum(15.);
    
    TLegend*    lMCProduction2 =   new TLegend(0.18,0.55,0.5,0.85);
    lMCProduction2->SetNColumns(2);
    lMCProduction2->AddEntry( hPlotResults_Full_stat, "Stat",  "EP" );
    lMCProduction2->AddEntry( hPlotResults_Full_syst, "Syst",  "F" );
    
    
    //
    //  --- Final Canvas
    //
    //  --- --- Upper Plot
    TCanvas*    cPlotResults = new TCanvas("cPlotResults","cPlotResults",1000,1000);
    TPad*   kUpperPlot  =   new TPad("kUpperPlot", "kUpperPlot", 0, 0.35, 1, 1.0);
    kUpperPlot      ->  SetLogy();
    gStyle          ->  SetOptStat(0);
    kUpperPlot->SetBottomMargin(0);
    kUpperPlot->Draw();
    kUpperPlot->cd();
    hPlotResults_Full_syst->Draw("SAME PE2");
    hPlotResults_Full_stat->Draw("SAME E1");
    auto iTer = 0;
    for ( auto kCurrent_MC : hFullQuantities )  {
        kCurrent_MC     ->  SetLineColor( uGetColor( iTer+2 ) );
        kCurrent_MC  ->  SetLineWidth( 2 );
        kCurrent_MC  ->  Draw("SAME HIST EP ][");
        if ( iTer+1 > fLegendTags.size() )          lMCProduction2->AddEntry( kCurrent_MC, kCurrent_MC->GetName(),  "EP" );
        else if ( !fLegendTags.at(iTer).IsNull() )  lMCProduction2->AddEntry( kCurrent_MC, fLegendTags.at(iTer),    "EP" );
        else                                        lMCProduction2->AddEntry( kCurrent_MC, kCurrent_MC->GetName(),  "EP" );
        iTer++;
    }
    lMCProduction2->Draw("same");
    //
    //  --- --- Lower Plot
    cPlotResults-> cd();
    TPad*   kLowerPlot  =   new TPad("kLowerPlot", "kLowerPlot", 0, 0.0, 1, 0.35);
    gStyle          ->  SetOptStat(0);
    kLowerPlot->SetTopMargin(0);
    kLowerPlot->SetBottomMargin(0.29);
    kLowerPlot->Draw();
    kLowerPlot->cd();
    hPlotResults_Ratio_syst->GetXaxis()->SetLabelSize(0.11);
    hPlotResults_Ratio_syst->GetXaxis()->SetLabelOffset(0.03);
    hPlotResults_Ratio_syst->SetMaximum(2.);
    hPlotResults_Ratio_syst->SetMinimum(0.);
    hPlotResults_Ratio_syst->Draw("SAME PE2");
    hPlotResults_Ratio_stat->Draw("SAME E1");

    
    /*
    hPlotResult4->Divide(hPlotResults_Full_stat,hPlotResults_Full_stat);
    hPlotResult4  ->  SetBinError(1,SquareSum( { hPlotResults_Full_stat->GetBinError(1)/hPlotResults_Full_stat->GetBinContent(1),  hPlotResults_Full_syst->GetBinError(1)/hPlotResults_Full_stat->GetBinContent(1) } ));
    hPlotResult4  ->  SetBinError(2,SquareSum( { hPlotResults_Full_stat->GetBinError(2)/hPlotResults_Full_stat->GetBinContent(2),  hPlotResults_Full_syst->GetBinError(2)/hPlotResults_Full_stat->GetBinContent(2) } ));
    hPlotResult4  ->  SetBinError(3,SquareSum( { hPlotResults_Full_stat->GetBinError(3)/hPlotResults_Full_stat->GetBinContent(3),  hPlotResults_Full_syst->GetBinError(3)/hPlotResults_Full_stat->GetBinContent(3) } ));
    hPlotResult4  ->  SetBinError(4,SquareSum( { hPlotResults_Full_stat->GetBinError(4)/hPlotResults_Full_stat->GetBinContent(4),  hPlotResults_Full_syst->GetBinError(4)/hPlotResults_Full_stat->GetBinContent(4) } ));
    hPlotResult4  ->  SetBinError(5,SquareSum( { hPlotResults_Full_stat->GetBinError(5)/hPlotResults_Full_stat->GetBinContent(5),  hPlotResults_Full_syst->GetBinError(5)/hPlotResults_Full_stat->GetBinContent(5) } ));
    hPlotResult4  ->  SetBinError(6,SquareSum( { hPlotResults_Full_stat->GetBinError(6)/hPlotResults_Full_stat->GetBinContent(6),  hPlotResults_Full_syst->GetBinError(6)/hPlotResults_Full_stat->GetBinContent(6) } ));
    hPlotResult4->SetMinimum(0.0);
    hPlotResult4->SetMaximum(2.0);
    hPlotResult4->Draw("SAME");
    */
    //
    iTer = 0;
    for ( auto kCurrent_MC : hFullQuantities )  {
        auto    kCurrent_Ratio  =   ( TH1F* )( kCurrent_MC->Clone() );
        auto    kDivider        =   uScale( hPlotResults_Full_stat, 1, -2 );
        kCurrent_Ratio  ->  GetYaxis()  ->  SetTitle( "" );
        kCurrent_Ratio  ->  SetTitle( "" );
        kCurrent_Ratio  ->  Divide( kCurrent_MC, kDivider );
        kCurrent_Ratio  ->  SetMaximum( 2.0 );
        kCurrent_Ratio  ->  SetMinimum( 0.0 );
        kCurrent_Ratio  ->  SetLineColor( uGetColor( iTer+2 ) );
        kCurrent_Ratio  ->  SetLineWidth ( 3 );
        kCurrent_Ratio  ->  Draw("SAME");
        iTer++;
    }
    //
    cPlotResults->SaveAs(Form("/Users/nikolajal/alice/AliAnalysisPhiCount/result/Yield_p_p__7TeV/MC_Production/Plots/TestXD.pdf"));
}
