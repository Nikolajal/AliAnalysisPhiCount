#include "../inc/AliAnalysisPhiPair.h"
#include "/Users/nikolajal/alice/AliAnalysisPhiCount/AliAnalysisUtility/ReweightEfficiency.C"
#include "/Users/nikolajal/alice/AliAnalysisPhiCount/src/MonteCarlo/MC_Plot.cxx"
#include "/Users/nikolajal/alice/AliAnalysisPhiCount/src/MonteCarlo/MC_PlotGen.cxx"
#include "/Users/nikolajal/alice/AliAnalysisPhiCount/src/MonteCarlo/mergetrees.C"
// !TODO: Can there really be a TODO in a Test Macro?
//
TH1F* maketh1f ( TGraphErrors* gTarget, TString kName = "" ) {
    Float_t*    kBinning    =   new Float_t[int(1e6)];
    for ( Int_t iBin = 0; iBin < gTarget->GetN(); iBin++ )    kBinning[iBin] = gTarget->GetPointX(iBin);
    kBinning[gTarget->GetN()] = gTarget->GetPointY(gTarget->GetN()-1);
    TH1F*   hTEST = new TH1F( kName, kName, gTarget->GetN(), kBinning );
    for ( Int_t iBin = 0; iBin < gTarget->GetN(); iBin++ ) {
        hTEST->SetBinContent( iBin+1, gTarget->GetErrorX(iBin) );
        hTEST->SetBinError  ( iBin+1, gTarget->GetErrorY(iBin) );
    }
    return hTEST;
}
//
TH1F* maketh1f( TString kInputFile, TString kName = "" ) {
    return maketh1f( new TGraphErrors(kInputFile), kName );
}

TH2F* maketh2f ( TGraphAsymmErrors* gTarget, TString kName = "" ) {
    Float_t*    kXBinning   =   new Float_t[int(1e6)];
    Float_t*    kYBinning   =   new Float_t[int(1e6)];
    Int_t       nXbin       =   0;
    Int_t       nYbin       =   0;
    kXBinning[ nXbin ]  =   gTarget->GetPointX      (0);
    kYBinning[ nYbin ]  =   gTarget->GetErrorXlow   (0);
    nXbin++;
    nYbin++;
    for ( Int_t iBin = 0; iBin < gTarget->GetN(); iBin++ )  {
        auto    kCurrent_XBin   =   gTarget->GetPointX      (iBin);
        auto    kCurrent_YBin   =   gTarget->GetErrorXlow   (iBin);
        auto    kCurrentXbin    =   nXbin;
        auto    kCurrentYbin    =   nYbin;
        auto    FillX           =   true;
        auto    FillY           =   true;
        for ( Int_t xBin = 0; xBin < kCurrentXbin; xBin++ )    {
            if ( fabs( kXBinning[ xBin ]  - kCurrent_XBin ) < 1.e-4 ) { FillX = false; }
        }
        if ( FillX ) {
            kXBinning[ nXbin ]   =   kCurrent_XBin;
            nXbin++;
        }
        for ( Int_t yBin = 0; yBin < kCurrentYbin; yBin++ )    {
            if ( fabs( kYBinning[ yBin ] - kCurrent_YBin ) < 1.e-4 )  FillY = false;
        }
        if ( FillY ) {
            kYBinning[ nYbin ]   =   kCurrent_YBin;
            nYbin++;
        }
    }
    kXBinning[ nXbin ]   = gTarget->GetPointY    ( gTarget->GetN() -1 );
    kYBinning[ nYbin ]   = gTarget->GetErrorXhigh( gTarget->GetN() -1 );
    //
    TH2F*   hReturn = new TH2F( kName, kName, nXbin, kXBinning, nYbin, kYBinning );
    for ( Int_t iBin = 0; iBin < gTarget->GetN(); iBin++ )  {
        auto    kCurrent_XBin   =   gTarget->GetPointX      (iBin);
        auto    kCurrent_YBin   =   gTarget->GetErrorXlow   (iBin);
        auto xBin = hReturn ->  GetXaxis()  ->  FindBin( kCurrent_XBin );
        auto yBin = hReturn ->  GetYaxis()  ->  FindBin( kCurrent_YBin );
        hReturn ->  SetBinContent   ( xBin, yBin, gTarget->GetErrorYlow     (iBin) );
        hReturn ->  SetBinError     ( xBin, yBin, gTarget->GetErrorYhigh    (iBin) );
    }
    return hReturn;
}
//
TH2F* maketh2f( TString kInputFile, TString kName = "" ) {
    return maketh2f( new TGraphAsymmErrors(kInputFile), kName );
}
//
Float_t
uCalculateCorrectionFactor
 ( std::vector<Float_t> kYields, std::vector<Float_t> kEffics, std::vector<Float_t> kNevnts ) {
    //
    //  Calculate overall correction factor
    auto kCorr  =   0.;
    auto kNum_  =   0.;
    auto kDen_  =   0.;
    auto kFullEvs   =   0.;
    for ( auto kCurrent_Events : kNevnts )  kFullEvs += kCurrent_Events;
    for ( Int_t iTer = 0; iTer < kYields.size(); iTer++ ) {
        kNum_  +=   (kNevnts.at(iTer)/kFullEvs)*(kEffics.at(iTer));
        kDen_  +=   (kYields.at(iTer));
    }
    auto kScaleFactor   =   kNum_/kDen_;
    for ( Int_t iTer = 0; iTer < kYields.size(); iTer++ ) {
        kCorr   +=   kScaleFactor*((kYields.at(iTer))/(kEffics.at(iTer)));
    }
    return kCorr;
}
//
void
TestMacro
()  {
    /*
    auto _1D_Yield = 0.0321405;
    auto _2D_Yield = 0.00145206;
    auto _1D_Error = 0.000177788;
    auto _2D_Error = 0.000301294;
    
    TH1F* hSigma = new TH1F("hSigma","hSigma",10000,-0.1,0.3);
    TH1F* hGamma = new TH1F("hGamma","hGamma",10000,-0.1,0.3);
    TH1F* hCheck1D = new TH1F("hCheck1D","hCheck1D",10000,-0.1,0.3);
    TH1F* hCheck2D = new TH1F("hCheck2D","hCheck2D",10000,-0.1,0.3);
    
    for ( Int_t i = 0; i < (int)(1.e7); i++ ) {
        auto k1DY   =   uRandomGen->Gaus( _1D_Yield, _1D_Error );
        auto k2DY   =   uRandomGen->Gaus( _2D_Yield, _2D_Error );
        hSigma->Fill( fSigmaPhiValue( k1DY, k2DY ) );
        hGamma->Fill( fGammaPhiValue( k1DY, k2DY ) );
        hCheck1D->Fill( k1DY );
        hCheck2D->Fill( k2DY );
    }
    
    
    TCanvas* c1 = new TCanvas();
    c1->Divide(2,2);
    c1->cd(1);
    hSigma->Draw();
    uLatex->DrawLatexNDC(0.15,0.8,Form("Calc: %.5f",fSigmaPhiError( _1D_Yield, _2D_Yield, _1D_Error, _2D_Error )));
    uLatex->DrawLatexNDC(0.15,0.75,Form("E_{rel}: %.5f",fSigmaPhiError( _1D_Yield, _2D_Yield, _1D_Error, _2D_Error )/(fSigmaPhiValue( _1D_Yield, _2D_Yield ))));
    c1->cd(2);
    hGamma->Draw();
    uLatex->DrawLatexNDC(0.15,0.8,Form("Calc: %.5f",fGammaPhiError( _1D_Yield, _2D_Yield, _1D_Error, _2D_Error )));
    uLatex->DrawLatexNDC(0.15,0.75,Form("E_{rel}: %.5f",fGammaPhiError( _1D_Yield, _2D_Yield, _1D_Error, _2D_Error )/(fGammaPhiValue( _1D_Yield, _2D_Yield ))));
    c1->cd(3);
    hCheck1D->Draw();
    c1->cd(4);
    hCheck2D->Draw();
    return;
    
    
    */
    
    /*
     cout << "0. The MPI-based original Pythia 8 scheme." << endl;
     cout << "1. The new more QCD based scheme." << endl;
     cout << "2. The new gluon-move model." << endl;
     cout << "3. The SK I e^+ e^- CR model." << endl;
     cout << "4. The SK II e^+ e^- CR model." << endl;
     cout << "5. Mode 0 from https://arxiv.org/pdf/1505.01681.pdf." << endl;
     cout << "6. Mode 2 from https://arxiv.org/pdf/1505.01681.pdf." << endl;
     cout << "7. Mode 3 from https://arxiv.org/pdf/1505.01681.pdf." << endl;
     cout << "8. Custom mode 1: QCD based CR w/ ropes" << endl;
     cout << "9. Custom mode 2: QCD based CR w/o ropes" << endl;
     */
    
    //for ( Int_t iTer = 0; iTer <10; iTer++ ) mergetrees(Form("/Users/nikolajal/alice/GridDownload/outpythia_%i.root",iTer),Form("/Users/nikolajal/alice/GridDownload/Pythia_%i.root",iTer),iTer);
    //for ( Int_t iTer = 0; iTer <10; iTer++ ) MC_PlotGen(Form("/Users/nikolajal/alice/GridDownload/Pythia_%i.root",iTer),Form("Pythia8_%i",iTer),"yield",-1,"_p_p__7TeV",iTer);
    
    fSetAllBins();
    
    TFile*  kChekcpp5TeV    =   new TFile("/Users/nikolajal/alice/AliAnalysisPhiCount/result/Multiplicity_p_p__5TeV_test/SignalExtraction/FitResults.root");
    TFile*  kChekcpp5Te3    =   new TFile("/Users/nikolajal/alice/AliAnalysisPhiCount/result/Multiplicity_p_p__5TeV_test/PreProcessing/IM_MonteCarloTruth.root");
    TFile*  kChekcpp5Te2    =   new TFile("/Users/nikolajal/Downloads/spectrum_mul_evmix.root");
    
    auto        fHEventCount    =   uLoadHistograms<0,TH1F> ( kChekcpp5TeV,  "fQC_Event_Enum_FLL"   );
    auto        h1D_Eff         =   uLoadHistograms<0,TH1F> ( kChekcpp5Te3,  "h1D_Eff"   );
    auto        fHEventCntMlt   =   uLoadHistograms<0,TH1F> ( kChekcpp5TeV,  "fQC_Event_Enum_V0M"   );
    
    auto        h1D_Nraw_MT     =   uLoadHistograms<1,TH1F> ( kChekcpp5TeV,  "h1D_Nraw_MT_%i" );
    auto        h1D_Nraw_MT_1   =   uLoadHistograms<0,TH1F> ( kChekcpp5Te2,  "hrawBin_0_5" );
    /*
    auto        h1D_Nraw_MT_2   =   uLoadHistograms<0,TH1F> ( kChekcpp5Te2,  "h5_10" );
    auto        h1D_Nraw_MT_3   =   uLoadHistograms<0,TH1F> ( kChekcpp5Te2,  "h10_15" );
    auto        h1D_Nraw_MT_4   =   uLoadHistograms<0,TH1F> ( kChekcpp5Te2,  "h15_20" );
    auto        h1D_Nraw_MT_5   =   uLoadHistograms<0,TH1F> ( kChekcpp5Te2,  "h20_30" );
    auto        h1D_Nraw_MT_6   =   uLoadHistograms<0,TH1F> ( kChekcpp5Te2,  "h30_40" );
    auto        h1D_Nraw_MT_7   =   uLoadHistograms<0,TH1F> ( kChekcpp5Te2,  "h40_50" );
    auto        h1D_Nraw_MT_8   =   uLoadHistograms<0,TH1F> ( kChekcpp5Te2,  "h50_70" );
    auto        h1D_Nraw_MT_9   =   uLoadHistograms<0,TH1F> ( kChekcpp5Te2,  "h70_100" );
     */
    auto        h1D_Nraw_MT_0   =   uLoadHistograms<0,TH1F> ( kChekcpp5Te2,  "hrawBin_0_100" );
    auto        h1D_Eff_0       =   uLoadHistograms<0,TH1F> ( kChekcpp5Te2,  "Efficiency_0_100" );
    
    auto iMult = 0;
    auto hPrevMeasurement = h1D_Eff_0;
    auto hThisMeasurement = h1D_Eff ;//h1D_Nraw_MT.at(iMult);
    
    
    auto    f1DCorrection    =   (1.)/(fEvaluateINELgt0(iMult-1,fHEventCntMlt) );
    auto    f2DCorrection    =   (1.)/(fEvaluateINELgt0(iMult-1,fHEventCntMlt) * kBR * kBR );
    //hThisMeasurement->Scale(1.,"width");
    //hThisMeasurement->Scale(f1DCorrection);
    hThisMeasurement->Scale(0.01);
    
    auto dump = new TCanvas();
    //
    hThisMeasurement    ->  SetTitle("");
    hThisMeasurement    ->  GetYaxis()  ->  SetTitleOffset(0.92);
    hThisMeasurement    ->  SetMarkerStyle(uGetMarker(3));
    hPrevMeasurement    ->  SetMarkerStyle(uGetMarker(4));
    //
    auto    hThisMeasurementRatio =   (TH1F*)(hThisMeasurement->Clone("test"));
    auto    hPrevMeasurementRatio =   (TH1F*)(hPrevMeasurement->Clone("tes2"));
    //
    uSetHisto( hThisMeasurementRatio, "SPT 1D" );
    uSetHisto( hPrevMeasurementRatio, "SPT 1D" );
    //
    fSetAllCustomFunctions();
    //hThisMeasurement->Fit( fLevyTsallis );
    //
    delete dump;
    
    TCanvas*        cDrawEfficiencies   =   new TCanvas("cDrawEfficiencies","cDrawEfficiencies",1300,800);
    //
    TLegend*        lEfficiencies   =   new TLegend(0.625,0.88,0.88,0.7);
    lEfficiencies   ->  SetNColumns(2);
    lEfficiencies   ->  SetFillColorAlpha(0.,0.);
    lEfficiencies   ->  SetLineColorAlpha(0.,0.);
    lEfficiencies->AddEntry(hThisMeasurement,"This Work","EP");
    lEfficiencies->AddEntry(hPrevMeasurement,"Prev Work","EP");
    //
    TPad*   kUpperPlot  =   new TPad("kUpperPlot", "kUpperPlot", 0, 0.33, 1, 1.0);
    gStyle  ->  SetOptStat(0);
    //kUpperPlot  ->  SetLogy();
    kUpperPlot  ->SetBottomMargin(0);
    kUpperPlot  ->Draw();
    kUpperPlot  ->cd();
    hThisMeasurement->GetYaxis()->SetTitleSize(0.05);
    hThisMeasurement->Draw("SAME []");
    hPrevMeasurement->Draw("same[]");
    lEfficiencies->Draw("SAME");
    //
    cDrawEfficiencies-> cd();
    TPad*   kLowerPlot  =   new TPad("kLowerPlot", "kLowerPlot", 0, 0.0, 1, 0.33);
    gStyle          ->  SetOptStat(0);
    kLowerPlot->SetTopMargin(0);
    kLowerPlot->SetBottomMargin(0.2);
    kLowerPlot->Draw();
    kLowerPlot->cd();
    kLowerPlot->SetGridy();
    //
    hThisMeasurementRatio   ->  Divide(hThisMeasurement,hPrevMeasurement);
    hPrevMeasurementRatio   ->  Divide(hThisMeasurement,hPrevMeasurement);
    hThisMeasurementRatio   ->  SetMarkerStyle(uGetMarker(3));
    hPrevMeasurementRatio   ->  SetMarkerStyle(uGetMarker(4));
    hThisMeasurementRatio   ->  SetMaximum(1.5);
    hThisMeasurementRatio   ->  SetMinimum(0.5);
    hThisMeasurementRatio   ->  GetXaxis()->SetTitleSize(0.1);
    hThisMeasurementRatio   ->  GetXaxis()->SetTitleOffset(0.72);
    hThisMeasurementRatio   ->  GetYaxis()->SetTitle("Ratio Nicola / Sushanta");
    hThisMeasurementRatio   ->  GetYaxis()->SetTitleSize(0.08);
    hThisMeasurementRatio   ->  GetYaxis()->SetTitleOffset(0.25);
    hThisMeasurementRatio   ->  Draw("SAME[]");
    hPrevMeasurementRatio   ->  Draw("SAME[]");
    //
    kUpperPlot->cd();
    
    return;
    
    
    MC_Plot({/*"HERWIG",*/"EPOS","Pythia6","Pythia8_0"/*,"Pythia8_1","Pythia8_2","Pythia8_3","Pythia8_4","Pythia8_5","Pythia8_6","Pythia8_7","Pythia8_8","Pythia8_9"*/},{/*"HERWIG",*/"EPOS LHC","Pythia6 Perugia 2011","Pythia8 Monasch 2013"/*,"More QCD","Gluon-move model","SK I e^+ e^- CR model","SK II e^+ e^- CR model"*/,"Pythia8 CR Mode 0","CR Mode 2","Pythia8 CR Mode 3","Pythia8 Ropes"/*,"QCD based CR w/o ropes"*/},"yield","_p_p__7TeV");
    
    return;
    
    /*
    auto dump = new TCanvas();
    TFile* kInput = new TFile("/Users/nikolajal/alice/AliAnalysisPhiCount/result/Yield_p_p__7TeV/PreviousResults/HEPData-ins1762364-v1-root.root");
    //
    TH1F*   fHEventCount    =   (TH1F*) kInput->Get( "Table 4/Hist1D_y1" );
    TH1F*   fHEventCoun2    =   (TH1F*) kInput->Get( "Table 4/Hist1D_y1_e1" );
    TH1F*   fHEventCoun3    =   (TH1F*) kInput->Get( "Table 4/Hist1D_y1_e2" );
    for ( Int_t iBin = 1; iBin <= fHEventCount->GetNbinsX(); iBin++ ) fHEventCount->SetBinError(iBin, SquareSum({fHEventCoun2->GetBinContent(iBin),fHEventCoun3->GetBinContent(iBin)}) );
    //
    TFile* kINFILE_3  =   new TFile   ("/Users/nikolajal/alice/AliAnalysisPhiCount/result/Yield_p_p__7TeV/SignalExtrapolation/FitResults.root");
    //
    auto    hResStat        =   uLoadHistograms<0,TH1F>( kINFILE_3, "h1D_Nres_stat", "h1D_Nraw_stat" );
    auto    hResSyst        =   uLoadHistograms<0,TH1F>( kINFILE_3, "h1D_Nres_syst", "h1D_Nraw_syst" );
    auto    hSumErr =   (TH1F*)(uSumErrors(hResStat,hResSyst)->Clone());
    //
    hSumErr->SetTitle("");
    hSumErr->GetYaxis()->SetTitleOffset(0.92);
    hSumErr         ->SetMarkerStyle(uGetMarker(3));
    fHEventCount    ->SetMarkerStyle(uGetMarker(4));
    //
    auto    hNewRat =   (TH1F*)(hSumErr->Clone());
    auto    hOldRat =   (TH1F*)(fHEventCount->Clone());
    //
    uSetHisto( hNewRat, "SPT 1D" );
    uSetHisto( fHEventCount, "SPT 1D" );
    fSetAllCustomFunctions();
    hSumErr->Fit(fLevyTsallis);
    delete dump;
    //
    TCanvas*        cDrawEfficiencies   =   new TCanvas("cDrawEfficiencies","cDrawEfficiencies",1300,800);
    //
    TLegend*        lEfficiencies   =   new TLegend(0.625,0.88,0.88,0.7);
    lEfficiencies   ->  SetNColumns(2);
    lEfficiencies   ->  SetFillColorAlpha(0.,0.);
    lEfficiencies   ->  SetLineColorAlpha(0.,0.);
    lEfficiencies->AddEntry(hSumErr,"This Work","EP");
    lEfficiencies->AddEntry(hOldRat,"Prev Work","EP");
    //
    TPad*   kUpperPlot  =   new TPad("kUpperPlot", "kUpperPlot", 0, 0.33, 1, 1.0);
    gStyle  ->  SetOptStat(0);
    kUpperPlot  ->  SetLogy();
    kUpperPlot  ->SetBottomMargin(0);
    kUpperPlot  ->Draw();
    kUpperPlot  ->cd();
    hSumErr->GetYaxis()->SetTitleSize(0.05);
    hSumErr->Draw("SAME []");
    fHEventCount->Draw("same[]");
    lEfficiencies->Draw("SAME");
    //
    cDrawEfficiencies-> cd();
    TPad*   kLowerPlot  =   new TPad("kLowerPlot", "kLowerPlot", 0, 0.0, 1, 0.33);
    gStyle          ->  SetOptStat(0);
    kLowerPlot->SetTopMargin(0);
    kLowerPlot->SetBottomMargin(0.2);
    kLowerPlot->Draw();
    kLowerPlot->cd();
    //
    hNewRat->Divide(fLevyTsallis);
    hOldRat->Divide(fLevyTsallis);
    hNewRat         ->SetMarkerStyle(uGetMarker(3));
    hOldRat         ->SetMarkerStyle(uGetMarker(4));
    hNewRat->SetMaximum(1.35);
    hNewRat->SetMinimum(0.65);
    hNewRat->GetXaxis()->SetTitleSize(0.1);
    hNewRat->GetXaxis()->SetTitleOffset(0.72);
    hNewRat->GetYaxis()->SetTitle("Ratio to Levy Fit");
    hNewRat->GetYaxis()->SetTitleSize(0.08);
    hNewRat->GetYaxis()->SetTitleOffset(0.25);
    hNewRat->Draw("SAME[]");
    hOldRat->Draw("SAME[]");
    //
    kUpperPlot->cd();
    //
    //delete cDrawEfficiencies;

    //
    
    /*
    SetStyle();
    TH1F*   hPlotResults = new TH1F( "hPlotResults", "", 6, 0.5, 6.5 );
    TH1F*   hPlotResult2 = new TH1F( "hPlotResult2", "", 6, 0.5, 6.5 );
    TH1F*   hPlotResult3 = new TH1F( "hPlotResult3", "", 6, 0.5, 6.5 );
    
    hPlotResults->GetXaxis()->SetBinLabel(1,"#frac{dN_{#phi}}{dy}");
    hPlotResults->GetXaxis()->SetBinLabel(2,"#frac{dN_{#phi#phi}}{dy}");
    hPlotResults->GetXaxis()->SetBinLabel(3,"#frac{#LT Y_{#phi#phi} #GT }{ #LT Y_{#phi} #GT }");
    hPlotResults->GetXaxis()->SetBinLabel(4,"#frac{#LT Y_{#phi#phi} #GT }{ #LT Y_{#phi} #GT^{2} }");
    hPlotResults->GetXaxis()->SetBinLabel(5,"#sigma_{#phi}");
    hPlotResults->GetXaxis()->SetBinLabel(6,"#gamma_{#phi}");
    
    hPlotResult2->GetXaxis()->SetBinLabel(1,"#frac{dN_{#phi}}{dy}");
    hPlotResult2->GetXaxis()->SetBinLabel(2,"#frac{dN_{#phi#phi}}{dy}");
    hPlotResult2->GetXaxis()->SetBinLabel(3,"#frac{#LT Y_{#phi#phi} #GT }{ #LT Y_{#phi} #GT }");
    hPlotResult2->GetXaxis()->SetBinLabel(4,"#frac{#LT Y_{#phi#phi} #GT }{ #LT Y_{#phi} #GT^{2} }");
    hPlotResult2->GetXaxis()->SetBinLabel(5,"#sigma_{#phi}");
    hPlotResult2->GetXaxis()->SetBinLabel(6,"#gamma_{#phi}");
    
    hPlotResults->SetBinContent(1,0.0318336);
    hPlotResult2->SetBinContent(1,0.0318336);
    hPlotResults->SetBinError(1,.000174066);
    hPlotResult2->SetBinError(1,.0026183);
    
    hPlotResults->SetBinContent(2,0.00143711);
    hPlotResult2->SetBinContent(2,0.00143711);
    hPlotResults->SetBinError(2,.000269857);
    hPlotResult2->SetBinError(2,.000226415);
    
    hPlotResults->SetBinContent(3,0.0451443);
    hPlotResult2->SetBinContent(3,0.0451443);
    hPlotResults->SetBinError(3,.00848068);
    hPlotResult2->SetBinError(3,.00475983);
    
    hPlotResults->SetBinContent(4,1.41813);
    hPlotResult2->SetBinContent(4,1.41813);
    hPlotResults->SetBinError(4,.268094);
    hPlotResult2->SetBinError(4,.0799365);
    
    hPlotResults->SetBinContent(5,0.0336945);
    hPlotResult2->SetBinContent(5,0.0336945);
    hPlotResults->SetBinError(5,.000565412);
    hPlotResult2->SetBinError(5,.00313393);
    
    hPlotResults->SetBinContent(6,0.0584549);
    hPlotResult2->SetBinContent(6,0.0584549);
    hPlotResults->SetBinError(6,.0176219);
    hPlotResult2->SetBinError(6,.00670411);
    
    uSetHisto( hPlotResults, "SPT 12D STAT");
    uSetHisto( hPlotResult2, "SPT 12D SYST");
    hPlotResults->GetXaxis()->SetTitle("");
    hPlotResults->GetYaxis()->SetTitle("");
    hPlotResult2->GetXaxis()->SetTitle("");
    hPlotResult2->GetYaxis()->SetTitle("");
    
    hPlotResult2->SetMinimum(1.e-4);
    hPlotResult2->SetMaximum(5.);
    
    TCanvas*    cPlotResults = new TCanvas("cPlotResults","cPlotResults",1000,1000);
    TPad*   kUpperPlot  =   new TPad("kUpperPlot", "kUpperPlot", 0, 0.35, 1, 1.0);
    kUpperPlot      ->  SetLogy();
    //kUpperPlot      ->  SetGridx();
    //kUpperPlot      ->  SetGridy();
    gStyle          ->  SetOptStat(0);
    kUpperPlot->SetBottomMargin(0);
    kUpperPlot->Draw();
    kUpperPlot->cd();
    
    hPlotResult2->Draw("SAME");
    hPlotResults->Draw("SAME E1");
    
    cPlotResults-> cd();
    TPad*   kLowerPlot  =   new TPad("kLowerPlot", "kLowerPlot", 0, 0.0, 1, 0.35);
    //kLowerPlot      ->  SetGridx();
    //kLowerPlot      ->  SetGridy();
    gStyle          ->  SetOptStat(0);
    kLowerPlot->SetTopMargin(0);
    kLowerPlot->SetBottomMargin(0.29);
    kLowerPlot->Draw();
    kLowerPlot->cd();
    hPlotResults->GetXaxis()->SetLabelSize(0.11);
    hPlotResults->GetXaxis()->SetLabelOffset(0.03);
    
    auto hPlotResult4 = (TH1F*)(hPlotResults->Clone("hPlotResult4"));
    
    hPlotResult4->Divide(hPlotResults,hPlotResults);
    
    hPlotResult4->SetMinimum(0.5);
    hPlotResult4->SetMaximum(1.5);
    
    hPlotResult4->Draw("SAME");

    
    //cPlotResults;
    */
    /*
    fSetAllBins();
    //
    //  --- EFF
    TFile* kINFILE_1  =   new TFile   ("/Users/nikolajal/alice/AliAnalysisPhiCount/result/Yield_p_p__7TeV/PreProcessing/IM_MonteCarloTruth.root");
    //
    auto    hEffNew         =   uLoadHistograms<0,TH1F>( kINFILE_1, "h1D_Eff", "h1D_Eff_New" );
    auto    hRecNew         =   uLoadHistograms<0,TH1F>( kINFILE_1, "h1D_Nrec", "h1D_rec_New" );
    auto    hGenNew         =   uLoadHistograms<0,TH1F>( kINFILE_1, "h1D_Ngen", "h1D_gen_New" );
    auto    hGenNew_2Db     =   uLoadHistograms<0,TH1F>( kINFILE_1, "h1D_Ngen_2Db", "h1D_gen_New_2Db" );
    auto    hTruNew         =   uLoadHistograms<0,TH1F>( kINFILE_1, "h1D_Ntru", "h1D_tru_New" );
    //
    //  --- EFF
    TFile* kINFILE_3  =   new TFile   ("/Users/nikolajal/alice/AliAnalysisPhiCount/result/Yield_p_p__7TeV/SignalExtrapolation/FitResults.root");
    //
    auto    hResStat        =   uLoadHistograms<0,TH1F>( kINFILE_3, "h1D_Nraw_stat", "h1D_Nraw_stat" );
    auto    hResSyst        =   uLoadHistograms<0,TH1F>( kINFILE_3, "h1D_Nraw_syst", "h1D_Nraw_syst" );
    //
    TFile* kOUTILE_1  =   new TFile   ("/Users/nikolajal/alice/AliAnalysisPhiCount/test_Anders.root","recreate");
    //
    fSetAllFunctions();
    ReweightEfficiency( uSumErrors( hResStat, hResSyst ), fLevyTsallis );
    //
    TFile* kOUTILE_2  =   new TFile   ("/Users/nikolajal/alice/AliAnalysisPhiCount/test_NewMet.root","recreate");
    //
    hEffNew ->Write();
    hRawNew ->Write();
    hRecNew ->Write();
    hGenNew ->Write();
    hGenNew_2Db ->Write();
    hTruNew ->Write();
    //
    kOUTILE_2   ->  Close();
    kOUTILE_1   ->  Close();
    kINFILE_3   ->  Close();
    kINFILE_1   ->  Close();
    //
    //TFile*          kOutput = new TFile( "./out.root", "RECREATE" );
    //
    //maketh1f(kMPT,"kMPT")->Write();
    //maketh1f("/Users/nikolajal/alice/AliAnalysisPhiCount/AliAnalysisRivetPhiCount/rivet-plots/ALICE_2022_Test/ptspec.0333.1D.dat","kS1D")->Write();
    //maketh1f(kS2D,"kS2D")->Write();
    //maketh1f("/Users/nikolajal/alice/AliAnalysisPhiCount/AliAnalysisRivetPhiCount/rivet-plots/ALICE_2022_Test/prdprb.0333.XD.dat","kSPR")->Write();
    //maketh1f("/Users/nikolajal/alice/AliAnalysisPhiCount/AliAnalysisRivetPhiCount/rivet-plots/ALICE_2022_Test/phidif.0333.1D.dat","kPHI")->Write();
    //maketh2f("/Users/nikolajal/alice/AliAnalysisPhiCount/AliAnalysisRivetPhiCount/rivet-plots/ALICE_2022_Test/ptspec.0333.2D_Rivet.dat","kS2D")->Write();
    //
    kOutput->Close();
    return;
     */
    //
    /*
    for ( Int_t i = 0; i < nBinPT2D; i++ )  {
        for ( Int_t j = i; j < nBinPT2D; j+=2 )  {
            cout << "\\begin{figure}[!h]" << endl;
            cout << "\\centering" << endl;
            cout << Form("\\includegraphics[width=0.9\\linewidth]{../result/Yield/SignalExtraction/Plots/2D/PT_%2.1f_%2.1f__%2.1f_%2.1f_.pdf}",fArrPT2D[i],fArrPT2D[i+1],fArrPT2D[j],fArrPT2D[j+1]) << endl;
            cout << Form("\\includegraphics[width=0.9\\linewidth]{../result/Yield/SignalExtraction/Plots/2D/PT_%2.1f_%2.1f__%2.1f_%2.1f_.pdf}",fArrPT2D[i],fArrPT2D[i+1],fArrPT2D[j+1],fArrPT2D[j+2]) << endl;
            cout << "\\end{figure}[!h]" << endl;
        }
    }
    
    return;
    TFile*  kINFILE_HI   =   new TFile   ("/Users/nikolajal/alice/AliAnalysisPhiCount/result/Reference_HI_MB.root");
    TProfile * hTest = (TProfile*)(kINFILE_HI->Get("hHISysErr"));
    TH1F* hTest2 = new TH1F("","",nBinPT1D,fArrPT1D);
    
    uRebin(hTest2,hTest);
    
    hTest2->Draw();
    hTest->Draw("SAME");
    return;
    
    std::vector<TString> fInputFiles;
    for ( auto kOption : kSyst_SEX_1D_Options ) {
        fInputFiles.push_back( TString( Form( "/Users/nikolajal/alice/AliAnalysisPhiCount/result/Yield/SignalExtraction/Systematics/ExtractionCheck/%s/1D/FitResults_%s.root" , kOption.Data(), kOption.Data() ) ) );
    }
    uCompareResultsTH1F( fInputFiles, "hRAW_1D", kSyst_SEX_1D_Legend );
    
    return;
    
    TFile* kINFILE_1  =   new TFile   ("/Users/nikolajal/alice/AliAnalysisPhiCount/result/Yield/PreProcessing/IM_MonteCarloTruth.root");
    TFile* kINFILE_2  =   new TFile   ("/Users/nikolajal/alice/AliAnalysisPhiCount/result/Yield/Systematics/MaterialBudget/PreProcessing/IM_MonteCarloTruth_902.root");
    TFile* kINFILE_3  =   new TFile   ("/Users/nikolajal/alice/AliAnalysisPhiCount/result/Yield/Systematics/MaterialBudget/PreProcessing/IM_MonteCarloTruth_900.root");
    TFile* kINFILE_4  =   new TFile   ("/Users/nikolajal/alice/AliAnalysisPhiCount/result/Yield/Systematics/MaterialBudget/PreProcessing/IM_MonteCarloTruth_900.root");
    
    TH1F*   hEFF_FLL    =   (TH1F*)(kINFILE_1->Get("hEFF_1D"));
    TH1F*   hEFF_200    =   (TH1F*)(kINFILE_2->Get("hEFF_1D"));
    TH1F*   hEFF_201    =   (TH1F*)(kINFILE_3->Get("hEFF_1D"));
    
    TH1F*   hREC_FLL    =   (TH1F*)(kINFILE_1->Get("hREC_Rw_1D"));
    TH1F*   hREC_200    =   (TH1F*)(kINFILE_2->Get("hREC_Rw_1D"));
    TH1F*   hREC_201    =   (TH1F*)(kINFILE_3->Get("hREC_Rw_1D"));
    TH1F*   hREC_202    =   (TH1F*)(kINFILE_4->Get("hREC_Rw_1D"));
    TH1F*   hGEN_FLL    =   (TH1F*)(kINFILE_1->Get("hGEN_Rw_1D"));
    TH1F*   hGEN_200    =   (TH1F*)(kINFILE_2->Get("hGEN_Rw_1D"));
    TH1F*   hGEN_201    =   (TH1F*)(kINFILE_3->Get("hGEN_Rw_1D"));
    TH1F*   hGEN_202    =   (TH1F*)(kINFILE_4->Get("hGEN_Rw_1D"));
    
    auto hRatio1 = (TH1F*)hEFF_FLL->Clone();
    auto hRatio2 = (TH1F*)hEFF_FLL->Clone();
    auto hRatio3 = (TH1F*)hEFF_FLL->Clone();
    
    auto hRatio4 = (TH1F*)hREC_FLL->Clone();
    auto hRatio5 = (TH1F*)hREC_FLL->Clone();
    auto hRatio6 = (TH1F*)hREC_FLL->Clone();
    auto hRatio7 = (TH1F*)hREC_FLL->Clone();
    auto hRatio8 = (TH1F*)hREC_FLL->Clone();
    auto hRatio9 = (TH1F*)hREC_FLL->Clone();
    
    hRatio1->Divide(hEFF_FLL,hEFF_200,1.,1.,"b");
    hRatio2->Divide(hEFF_FLL,hEFF_201,1.,1.,"b");
    hRatio3->Divide(hEFF_200,hEFF_201,1.,1.,"b");
    
    hRatio4->Divide(hREC_FLL,hGEN_FLL,1.,1.,"b");
    hRatio5->Divide(hREC_200,hGEN_200,1.,1.,"b");
    hRatio6->Divide(hREC_201,hGEN_201,1.,1.,"b");
    
    hRatio7->Divide(hRatio4,hRatio5,1.,1.,"b");
    hRatio8->Divide(hRatio4,hRatio6,1.,1.,"b");
    hRatio9->Divide(hRatio5,hRatio6,1.,1.,"b");

    uOffset(hRatio1,-1,true);
    uOffset(hRatio2,-1,true);
    //uOffset(hRatio3,-1,true);
    
    TH1F* hReference1D = new TH1F("hReference1D","",nBinPT1D,fArrPT1D);
    TH1F* hReference2D = new TH1F("hReference2D","",nBinPT2D,fArrPT2D);
    
    
    TF1* fCorr = new TF1("fCorr","1 + [0] * TMath::Exp( [1] + x * [2] )",0,100);
    fCorr->SetParameter(0,+1.00);
    fCorr->SetParameter(1,-0.01);
    fCorr->SetParameter(2,+0.03);
    
    TCanvas* c1 = new TCanvas("","",1600,400);
    c1->Divide(3,2);
    c1->cd(1);
    gPad->SetLogx();
    hRatio7->Fit(fCorr,"IMEQS","R",0.,5.);
    hRatio7->DrawCopy("HIST");
    fCorr->DrawCopy("SAME");
    for ( Int_t i = 1; i <= nBinPT1D; i++ ) {
        auto kRefErr = fCorr->Integral( fArrPT1D[i-1], fArrPT1D[i] ) / ( fArrPT1D[i] - fArrPT1D[i-1] );
        hReference1D->SetBinContent( i, kRefErr-1 );
    }
    c1->cd(4);
    gPad->SetLogx();
    hReference1D->DrawCopy("HIST");
    fCorr->SetParameter(0,+1.00);
    fCorr->SetParameter(1,-0.01);
    fCorr->SetParameter(2,+0.03);
    c1->cd(2);
    gPad->SetLogx();
    hRatio8->Fit(fCorr,"IMEQS","R",0.,5.);
    hRatio8->DrawCopy("HIST");
    fCorr->DrawCopy("SAME");
    for ( Int_t i = 1; i <= nBinPT1D; i++ ) {
        auto kRefErr = fCorr->Integral( fArrPT1D[i-1], fArrPT1D[i] ) / ( fArrPT1D[i] - fArrPT1D[i-1] );
        hReference1D->SetBinContent( i, kRefErr-1 );
    }
    c1->cd(5);
    gPad->SetLogx();
    hReference1D->DrawCopy("HIST");
    fCorr->SetParameter(0,+1.00);
    fCorr->SetParameter(1,-0.01);
    fCorr->SetParameter(2,+0.03);
    c1->cd(3);
    gPad->SetLogx();
    hRatio9->Fit(fCorr,"IMEQS","R",0.,5.);
    hRatio9->DrawCopy("HIST");
    fCorr->DrawCopy("SAME");
    for ( Int_t i = 1; i <= nBinPT1D; i++ ) {
        auto kRefErr = fCorr->Integral( fArrPT1D[i-1], fArrPT1D[i] ) / ( fArrPT1D[i] - fArrPT1D[i-1] );
        hReference1D->SetBinContent( i, kRefErr-1 );
    }
    c1->cd(6);
    gPad->SetLogx();
    hReference1D->DrawCopy("HIST");
    
    return;
    
    gROOT->ProcessLine(".! $ROOT_SYS hadd -f ./result/MCout.root ./result/MC_Production/outGeneratorMC_00*.root");
    auto kREBIN = 1;
    TFile* kINFILE2  =   new TFile   ("/Users/nikolajal/alice/AliAnalysisPhiCount/result/MCout.root");
    auto hPhi1 = (TProfile*)(kINFILE2->Get("hYPhiPro1"));
    auto hPhi2 = (TProfile*)(kINFILE2->Get("hYPhiPro2"));
    hPhi1->Rebin(kREBIN);
    hPhi2->Rebin(kREBIN);
    cout << hPhi1->Integral("width") << endl;
    auto hPhiGamma  = new TH1F("hPhiGamma", "hPhiGamma",    (int)(100./kREBIN),0.,100);
    auto hPhiSigma  = new TH1F("hPhiSigma", "hPhiSigma",    (int)(100./kREBIN),0.,100);
    auto hPhiR2     = new TH1F("hPhiR2",    "hPhiR2",       (int)(100./kREBIN),0.,100);
    auto hPhiR1     = new TH1F("hPhiR1",    "hPhiR1",       (int)(100./kREBIN),0.,100);
    
    for ( Int_t iTer = 1; iTer <= (int)(100./kREBIN); iTer++ )  {
        auto nPhi1 = hPhi1->GetBinContent(iTer);
        auto nPhi2 = hPhi2->GetBinContent(iTer);
        auto nPhE1 = hPhi1->GetBinError(iTer);
        auto nPhE2 = hPhi2->GetBinError(iTer);
        if ( nPhi1 == 0 ) continue;
        if ( nPhi2 == 0 ) continue;
        hPhiGamma->SetBinContent( iTer, fGammaPhiValue( nPhi1, nPhi2 ) );
        hPhiGamma->SetBinError  ( iTer, fGammaPhiError( nPhi1, nPhi2, nPhE1, nPhE2 ) );
        hPhiSigma->SetBinContent( iTer, fSigmaPhiValue( nPhi1, nPhi2 ) );
        hPhiSigma->SetBinError  ( iTer, fSigmaPhiError( nPhi1, nPhi2, nPhE1, nPhE2 ) );
        hPhiR2->SetBinContent   ( iTer, nPhi2/(nPhi1*nPhi1) );
        hPhiR2->SetBinError     ( iTer, (nPhi2/(nPhi1*nPhi1))*sqrt( nPhE2*nPhE2/(nPhi2*nPhi2) + 4*nPhE1*nPhE1/(nPhi1*nPhi1) ) );
        hPhiR1->SetBinContent   ( iTer, nPhi2/(nPhi1) );
        hPhiR1->SetBinError     ( iTer, (nPhi2/(nPhi1))*sqrt( nPhE2*nPhE2/(nPhi2*nPhi2) + nPhE1*nPhE1/(nPhi1*nPhi1) ) );
    }
    
    hPhi1->SetLineColor(2);
    hPhi1->SetLineStyle(5);
    hPhi1->SetLineWidth(3);
    
    hPhi2->SetLineColor(2);
    hPhi2->SetLineStyle(5);
    hPhi2->SetLineWidth(3);
    
    hPhiR1->SetLineColor(2);
    hPhiR1->SetLineStyle(5);
    hPhiR1->SetLineWidth(3);
    
    hPhiR2->SetLineColor(2);
    hPhiR2->SetLineStyle(5);
    hPhiR2->SetLineWidth(3);
    
    hPhiGamma->SetLineColor(2);
    hPhiGamma->SetLineStyle(5);
    hPhiGamma->SetLineWidth(3);
    
    hPhiSigma->SetLineColor(2);
    hPhiSigma->SetLineStyle(5);
    hPhiSigma->SetLineWidth(3);
    
    TFile* kINFILE  =   new TFile   (Form(kASigExtp_FitCheckRst,"Multiplicity"));
    float    kNCH[]  =   {5.398,9.44,13.13,17.58,24};
    char*    kName[] =   {"1D","2D","R1","R2","P1","P2"};
    char*    kLabel[]=    {"#frac{dN_{#phi}}{dy}","#frac{dN_{#phi#phi}}{dy}","#frac{#LT Y_{#phi#phi} #GT}{#LT Y_{#phi} #GT}","#frac{#LT Y_{#phi#phi} #GT}{#LT Y_{#phi} #GT^{2}}","#sigma^{2}_{#phi}","#gamma_{#phi}"};
    
    
    std::vector<TGraphErrors*> fOutput;
    for ( int i = 0; i < 6; i++ )   {
        auto fCurrentTH1F = (TH1D*)(kINFILE->Get(Form("hShow%s",kName[i])));
        fOutput.push_back(new TGraphErrors());
        for ( int j = 0; j < 5; j++ )   {
            fOutput.at(i)->SetPoint(j,kNCH[j],fCurrentTH1F->GetBinContent(j+1));
            fOutput.at(i)->SetPointError(j,0.,fCurrentTH1F->GetBinError(j+1));
        }
        fOutput.at(i)->SetPoint(5,-1,fCurrentTH1F->GetBinContent(6));
        fOutput.at(i)->SetPointError(5,0.,fCurrentTH1F->GetBinError(6));
        fOutput.at(i)->SetMinimum(0);
        fOutput.at(i)->SetMarkerStyle(4);
        fOutput.at(i)->SetMarkerColor(4);
        fOutput.at(i)->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
        if ( i == 5 ) fOutput.at(i)->SetMaximum(0.06);
        if ( i == 2 ) fOutput.at(i)->SetMaximum(0.08);
        if ( i == 3 ) fOutput.at(i)->SetMinimum(0.40);
        if ( i == 3 ) fOutput.at(i)->SetMaximum(1.50);
    }
    /*
    fOutput.push_back(new TGraphErrors());
    for ( int j = 0; j < 5; j++ )   {
        double fY1D, fY2D, fDump;
        fOutput.at(0)->GetPoint(j,fDump,fY1D);
        fOutput.at(4)->GetPoint(j,fDump,fY2D);
        
        fOutput.at(6)->SetPoint(j,kNCH[j],fY1D/fY2D);
        fOutput.at(6)->SetPointError(j,0.,0);
        fOutput.at(6)->SetMarkerStyle(markers[2+(6%4)]);
        fOutput.at(6)->SetMarkerColor(colors[1+(6%3)]);
    }
     */
    /*
    TCanvas* cTest = new TCanvas("","",1000,700);
    cTest->Divide(3,2);
    auto iTer = 1;
    for ( auto kGraph : fOutput )    {
        cTest->cd(iTer);
        kGraph->Draw("APE");
        uLatex->DrawLatexNDC(0.18,0.80,kLabel[iTer-1]);
        if ( iTer == 1 )    {
            hPhi1->DrawCopy("same HIST L ");
            hPhi1->SetFillColorAlpha(3,0.33);
            hPhi1->DrawCopy("same E3L ");
        }
        if ( iTer == 2 )    {
            hPhi2->DrawCopy("same HIST L ");
            hPhi2->SetFillColorAlpha(3,0.33);
            hPhi2->DrawCopy("same E3L ");
        }
        if ( iTer == 3 )    {
            hPhiR1->DrawCopy("same HIST L ");
            hPhiR1->SetFillColorAlpha(3,0.33);
            hPhiR1->DrawCopy("same E3L ");
        }
        if ( iTer == 4 )    {
            hPhiR2->DrawCopy("same HIST L ");
            hPhiR2->SetFillColorAlpha(3,0.33);
            hPhiR2->DrawCopy("same E3L ");
        }
        if ( iTer == 5 )    {
            hPhiSigma->DrawCopy("same HIST L ");
            hPhiSigma->SetFillColorAlpha(3,0.33);
            hPhiSigma->DrawCopy("same E3L ");
        }
        if ( iTer == 6 )    {
            hPhiGamma->DrawCopy("same HIST L ");
            hPhiGamma->SetFillColorAlpha(3,0.33);
            hPhiGamma->DrawCopy("same E3L ");
        }
        iTer++;
    }
    /*
    auto fResult = fOutput.at(0)->Fit("pol1","IMEQ0S");
    auto k1D_0  =   fResult->GetParams()[0];
    auto k1D_1  =   fResult->GetParams()[1];
    fResult = fOutput.at(1)->Fit("pol2","IMEQ0S");
    auto k2D_0  =   fResult->GetParams()[0];
    auto k2D_1  =   fResult->GetParams()[1];
    auto k2D_2  =   fResult->GetParams()[2];
    TH1F* cPrediction = new TH1F("","",30,0,30);
    for ( int i = 1; i <= 30; i++ )  {
        cPrediction->SetBinContent(i,fGammaPhiValue( k1D_0 + k1D_1*i, k2D_0 + k2D_1*i + k2D_2*i*i ) );
    }
    cPrediction->SetLineColor(colors[1]);
    cTest->cd(6);
    cPrediction->SetMinimum(0);
    cPrediction->Draw("HIST MIN0");
    fOutput.at(5)->Draw("EP SAME");
    cTest->cd(1);
    TMultiGraph* cmulti = new TMultiGraph();
    cmulti->Add(fOutput.at(0),"EP");
    cmulti->Add(fOutput.at(4),"EP");
    cmulti->Draw("APE");
     */
}
