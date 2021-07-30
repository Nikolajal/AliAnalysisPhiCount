// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../inc/AliAnalysisPhiPair.h"

void MakeSyst ()    {
    
    fSetAllBins();
    
    TH1F    *hTotalSystematics  =   new TH1F("Syst_",       "Total",                nBinPT1D,fArrPT1D);
    TH1F    *hSignalExtraction  =   new TH1F("SE_1D",       "Signal Extraction",    nBinPT1D,fArrPT1D);
    TH1F    *hPIDSelection      =   new TH1F("PID1D",       "PID",                  nBinPT1D,fArrPT1D);
    TH1F    *hTrackCuts         =   new TH1F("hTrackCuts",  "Track Cuts",           nBinPT1D,fArrPT1D);
    TH1F    *hTracking          =   new TH1F("hTracking",   "Tracking",             nBinPT1D,fArrPT1D);
    
    for ( Int_t iPT1D = 1; iPT1D <= nBinPT1D; iPT1D++ )  {
        // Signal Extraction
        if ( iPT1D == 1  || iPT1D == 20 )   hSignalExtraction       ->  SetBinContent   (iPT1D, 0.040);
        if ( iPT1D >= 2  && iPT1D <= 9  )   hSignalExtraction       ->  SetBinContent   (iPT1D, 0.010);
        if ( iPT1D >= 10 && iPT1D <= 15 )   hSignalExtraction       ->  SetBinContent   (iPT1D, 0.015);
        if ( iPT1D >= 16 && iPT1D <= 19 )   hSignalExtraction       ->  SetBinContent   (iPT1D, 0.035);
        // PID
        if ( iPT1D == 1  || iPT1D == 20 )   hPIDSelection           ->  SetBinContent   (iPT1D, 0.020);
        if ( iPT1D >= 2  && iPT1D <= 2  )   hPIDSelection           ->  SetBinContent   (iPT1D, 0.015);
        if ( iPT1D >= 3  && iPT1D <= 3  )   hPIDSelection           ->  SetBinContent   (iPT1D, 0.010);
        if ( iPT1D >= 4  && iPT1D <= 19 )   hPIDSelection           ->  SetBinContent   (iPT1D, 0.005);
        //  Track cuts
        if ( iPT1D != 0  || iPT1D != 1  )   hTrackCuts              ->  SetBinContent   (iPT1D, 0.030);
        //  Tracking
        if ( iPT1D != 0  || iPT1D != 1  )   hTracking               ->  SetBinContent   (iPT1D, 0.060);
        // Total
        hTotalSystematics   ->  SetBinContent   (iPT1D, SquareSum({hSignalExtraction->GetBinContent(iPT1D),hPIDSelection->GetBinContent(iPT1D),hTrackCuts->GetBinContent(iPT1D),hTracking->GetBinContent(iPT1D)}));
    }
    
    TCanvas*c1 = new TCanvas();
    gStyle->SetOptStat(0);
    gPad->SetLogx();
    
    // 2D Analysis
    
    //  Total
    hTotalSystematics   ->  SetMinimum(0.0);
    hTotalSystematics   ->  SetMaximum(0.12);
    hTotalSystematics   ->  SetLineColor(kRed);
    hTotalSystematics   ->  SetLineWidth(3);
    hTotalSystematics   ->  SetLineStyle(1);
    hTotalSystematics   ->  Draw("");
    //  Signal Extraction
    hSignalExtraction   ->  SetLineColor(kRainbowColor[1]);
    hSignalExtraction   ->  SetLineWidth(3);
    hSignalExtraction   ->  SetLineStyle(4);
    hSignalExtraction   ->  Draw("SAME");
    //  PID Selection
    hPIDSelection       ->  SetLineColor(kRainbowColor[2]);
    hPIDSelection       ->  SetLineWidth(3);
    hPIDSelection       ->  SetLineStyle(5);
    hPIDSelection       ->  Draw("SAME");
    //  Track Cuts
    hTrackCuts          ->  SetLineColor(kRainbowColor[3]);
    hTrackCuts          ->  SetLineWidth(3);
    hTrackCuts          ->  SetLineStyle(6);
    hTrackCuts          ->  Draw("SAME");
    //  Tracking
    hTracking           ->  SetLineColor(kRainbowColor[4]);
    hTracking           ->  SetLineWidth(3);
    hTracking           ->  SetLineStyle(7);
    hTracking           ->  Draw("SAME");
    
    TH2F    *hTotalSystematics  =   new TH2F("Syst_",       "Total",                nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    TH2F    *hSignalExtraction  =   new TH2F("SE_1D",       "Signal Extraction",    nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    TH2F    *hPIDSelection      =   new TH2F("PID1D",       "PID",                  nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    TH2F    *hTrackCuts         =   new TH2F("hTrackCuts",  "Track Cuts",           nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    TH2F    *hTracking          =   new TH2F("hTracking",   "Tracking",             nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    
    for ( Int_t iPT2D = 1; iPT2D <= nBinPT2D; iPT2D++ )  {
        for ( Int_t jPT2D = 1; jPT2D <= nBinPT2D; jPT2D++ )  {
            // Signal Extraction
            if ( iPT2D == 1  || iPT2D == 20 )   hSignalExtraction       ->  SetBinContent   (iPT2D,jPT2D, 0.040);
            if ( iPT2D >= 2  && iPT2D <= 9  )   hSignalExtraction       ->  SetBinContent   (iPT2D,jPT2D, 0.010);
            if ( iPT2D >= 10 && iPT2D <= 15 )   hSignalExtraction       ->  SetBinContent   (iPT2D,jPT2D, 0.015);
            if ( iPT2D >= 16 && iPT2D <= 19 )   hSignalExtraction       ->  SetBinContent   (iPT2D,jPT2D, 0.035);
            // PID
            if ( iPT2D == 1  || iPT2D == 20 )   hPIDSelection           ->  SetBinContent   (iPT2D,jPT2D, 0.020);
            if ( iPT2D >= 2  && iPT2D <= 2  )   hPIDSelection           ->  SetBinContent   (iPT2D,jPT2D, 0.015);
            if ( iPT2D >= 3  && iPT2D <= 3  )   hPIDSelection           ->  SetBinContent   (iPT2D,jPT2D, 0.010);
            if ( iPT2D >= 4  && iPT2D <= 19 )   hPIDSelection           ->  SetBinContent   (iPT2D,jPT2D, 0.005);
            //  Track cuts
            if ( iPT2D != 0  || iPT2D != 1  )   hTrackCuts              ->  SetBinContent   (iPT2D,jPT2D, 0.030);
            //  Tracking
            if ( iPT2D != 0  || iPT2D != 1  )   hTracking               ->  SetBinContent   (iPT2D,jPT2D, 0.060);
            // Total
            hTotalSystematics   ->  SetBinContent   (iPT2D, jPT2D, SquareSum({hSignalExtraction->GetBinContent(iPT2D),hPIDSelection->GetBinContent(iPT2D),hTrackCuts->GetBinContent(iPT2D),hTracking->GetBinContent(iPT2D)}));
        }
    }
    
    TCanvas*c2 = new TCanvas();
    c2->Divide(2,5);
    gStyle->SetOptStat(0);
    gPad->SetLogx();
    
    //  Total
    hTotalSystematics   ->  SetMinimum(0.0);
    hTotalSystematics   ->  SetMaximum(0.12);
    hTotalSystematics   ->  SetLineColor(kRed);
    hTotalSystematics   ->  SetLineWidth(3);
    hTotalSystematics   ->  SetLineStyle(1);
    hTotalSystematics   ->  Draw("");
    //  Signal Extraction
    hSignalExtraction   ->  SetLineColor(kRainbowColor[1]);
    hSignalExtraction   ->  SetLineWidth(3);
    hSignalExtraction   ->  SetLineStyle(4);
    hSignalExtraction   ->  Draw("SAME");
    //  PID Selection
    hPIDSelection       ->  SetLineColor(kRainbowColor[2]);
    hPIDSelection       ->  SetLineWidth(3);
    hPIDSelection       ->  SetLineStyle(5);
    hPIDSelection       ->  Draw("SAME");
    //  Track Cuts
    hTrackCuts          ->  SetLineColor(kRainbowColor[3]);
    hTrackCuts          ->  SetLineWidth(3);
    hTrackCuts          ->  SetLineStyle(6);
    hTrackCuts          ->  Draw("SAME");
    //  Tracking
    hTracking           ->  SetLineColor(kRainbowColor[4]);
    hTracking           ->  SetLineWidth(3);
    hTracking           ->  SetLineStyle(7);
    hTracking           ->  Draw("SAME");
    
    gPad->BuildLegend();
}
