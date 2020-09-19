// Global Functions and data structures file
// !TODO: General housekeeping and new manage rules
// !TODO: Headers clean-up
// !TODO: 2D plot options
// !TODO: Transfer to RooFit the Fit PT Count
// !TODO: Legend as in ...Count.C the slice plots

#ifndef SETFUNCTIONS_H
#define SETFUNCTIONS_H
#include "SetValues.h"
#include "SpectraUtils.C"

// C++
#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <chrono>

// ROOT
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TFile.h"

// RooFit
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooAbsData.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooGlobalFunc.h"
#include "RooMsgService.h"

// RooFitFunction
#include "RooChebychev.h"
#include "RooArgusBG.h"
#include "RooBreitWigner.h"
#include "RooExponential.h"
#include "RooVoigtian.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooUniform.h"

using namespace std;
using namespace RooFit;

//--------------------------------------------//
//  Utilities for histogram creations         //
//--------------------------------------------//

int             fLegendSelect       ( string fOption )
{
    if ( !fOption.compare("InvMass1D") )   return 1;
    if ( !fOption.compare("xInvMass2D") )  return 2;
    if ( !fOption.compare("yInvMass2D") )  return 2;
    else return -1;
}

void            fLegendMaker        ( RooPlot * fRooPlot, const char * fSelect, TLegend * fLegend )
{
    switch (fLegendSelect(fSelect))
    {
        case 1:
            fLegend                     ->SetFillColor(kWhite);
            fLegend                     ->SetLineColor(kWhite);
            fLegend                     ->AddEntry(fRooPlot->findObject("RooData"), "Data",                 "EP");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooSS"),   "Fit (Sig)",            "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooBB"),   "Fit (Bkg)",            "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooMod"),  "Fit (Model)",          "L");
            break;
        case 2:
            fLegend                     ->SetFillColor(kWhite);
            fLegend                     ->SetLineColor(kWhite);
            fLegend                     ->AddEntry(fRooPlot->findObject("RooData"), "Data",                 "EP");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooSS"),   "Fit (Sig #times Sig)", "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooBS"),   "Fit (Bkg #times Sig)", "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooSB"),   "Fit (Sig #times Bkg)", "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooBB"),   "Fit (Bkg #times Bkg)", "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooMod"),  "Fit (Model)",          "L");
            break;
        default:
            cout << "Improper option, no changes made" << endl;
            break;
    }
}

int             fAxisSelect         ( string fOption )
{
    if ( !fOption.compare("InvMass1D") )   return 1;
    if ( !fOption.compare("xInvMass2D") )  return 2;
    if ( !fOption.compare("yInvMass2D") )  return 3;
    else return -1;
}

void            fAxisMaker          ( RooPlot * fRooPlot, const char * fSelect )
{
    switch (fAxisSelect(fSelect))
    {
        case 1:
            fRooPlot                    ->GetXaxis()->SetTitle("m_{K_{+}K_{-}} (GeV/c^{2})");
            break;
        case 2:
            fRooPlot                    ->GetXaxis()->SetTitle("m^{x}_{K_{+}K_{-}} (GeV/c^{2})");
            break;
        case 3:
            fRooPlot                    ->GetXaxis()->SetTitle("m^{y}_{K_{+}K_{-}} (GeV/c^{2})");
            break;
        default:
            cout << "Improper option, no changes made" << endl;
            break;
    }
}

int             fPlotterSelect      ( string fOption )
{
    if ( !fOption.compare("InvMass1D") )   return 1;
    if ( !fOption.compare("xInvMass2D") )  return 2;
    if ( !fOption.compare("yInvMass2D") )  return 2;
    else return -1;
}

void            fRooPlotPlotter     ( RooPlot * fRooPlot, const char * fSelect, RooAddPdf fModel , RooDataHist * fData )
{
    switch (fPlotterSelect(fSelect))
    {
        case 1:
            fData                           ->plotOn(fRooPlot,      MarkerColor(38),                MarkerStyle(26),    Name("RooData"));
            fModel                          .plotOn (fRooPlot,      LineColor(4),                   LineStyle(kDashed), Name("RooMod"));
            fModel                          .plotOn (fRooPlot,      Components("fBkg"),             LineStyle(kDashed), LineColor(38),      Name("RooBB"));
            fModel                          .plotOn (fRooPlot,      Components("fSig"),             LineColor(2),       Name("RooSS"));
            break;
        case 2:
            fData                           ->plotOn(fRooPlot,      CutRange("fDrawRange"),         MarkerColor(38),    MarkerStyle(26) ,   Name("RooData"));
            fModel                          .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  LineColor(4),       LineStyle(kDashed), Name("RooMod"));
            fModel                          .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  Components("fBkg"), LineStyle(kDashed), LineColor(38),      Name("RooBB"));
            fModel                      .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  Components("fSigSig"),LineColor(2),       Name("RooSS"));
            fModel                      .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  Components("fSigBkg"),LineStyle(kDashed), LineColor(33),    Name("RooSB"));
            fModel                      .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  Components("fBkgSig"),LineStyle(kDashed), LineColor(36),    Name("RooBS"));
            break;
        default:
            cout << "Improper option, no changes made" << endl;
            break;
    }
}

void            fRooPlotMaker       ( RooPlot * fRooPlot, TLegend * fLegend, RooAddPdf fModel , RooDataHist * fData, const char * fSelect )
{
    fRooPlotPlotter(fRooPlot,fSelect,fModel,fData);
    fLegendMaker(fRooPlot,fSelect,fLegend);
    fAxisMaker(fRooPlot,fSelect);
}

void            setMarker           ( TH1F * hTarget, Int_t iType = 0 )
{
    hTarget->SetMarkerColor(kColor[iType]);
    hTarget->SetMarkerStyle(kStyle[iType]);
    hTarget->SetMarkerSize(kWidth[iType]);
}

void            setLine             ( TH1F * hTarget, Int_t iType = 0 )
{
    hTarget->SetLineColor(kColor[iType]);
    hTarget->SetLineStyle(kStyle[iType]);
    hTarget->SetLineWidth(kWidth[iType]);
}

void            setMarker           ( TH1D * hTarget, Int_t iType = 0 )
{
    hTarget->SetMarkerColor(kColor[iType]);
    hTarget->SetMarkerStyle(kStyle[iType]);
    hTarget->SetMarkerSize(kWidth[iType]);
}

void            setMarker           ( TGraphAsymmErrors * hTarget, Int_t iType = 0 )
{
    hTarget->SetMarkerColor(kColor[iType]);
    hTarget->SetMarkerStyle(kStyle[iType]);
    hTarget->SetMarkerSize(kWidth[iType]);
}

void            setLine             ( TH1D * hTarget, Int_t iType = 0 )
{
    hTarget->SetLineColor(kColor[iType]);
    hTarget->SetLineStyle(kStyle[iType]);
    hTarget->SetLineWidth(kWidth[iType]);
}

void            TH1_SetDescription ( TH1F * hTarget )
{
    hTarget->GetXaxis()->SetTitle("p_{T} #phi (GeV/c)");
    hTarget->GetYaxis()->SetTitle("#frac{d^{2}N_{#phi}}{dydp_{T}}(GeV/c)^{-1}");
    hTarget->GetXaxis()->SetTitleOffset(1.15);
    hTarget->GetYaxis()->SetTitleOffset(1.15);
}

void            TH2_SetDescription ( TH2F * hTarget )
{
    hTarget->GetXaxis()->SetTitle("p_{T} #phi_{1} (GeV/c)");
    hTarget->GetYaxis()->SetTitle("p_{T} #phi_{2} (GeV/c)");
    hTarget->GetXaxis()->SetTitleOffset(1.5);
    hTarget->GetYaxis()->SetTitleOffset(1.5);
}


// initialised values for fits in 2D
Float_t fLevyPar1 [] = {5.45,   5.05,   5.55,   5.85,   6.15,   6.10,   6.20,   5.2,    5.2,    4.1};
Float_t fLevyPar2 [] = {.238,   .234,   .262,   .283,   .302,   .315,   .347,   .344,   .385,   0.41};
Float_t fLevyPar3 [] = {0.69e-3,0.69e-3,0.69e-3,0.58e-3,0.51e-3,0.37e-3,.193e-3,83.2e-6,16.9e-6,1.01e-6};

//--------------------------------------------//
//  Signal Extraction from Invariant Mass hst //
//--------------------------------------------//

enum            fitresults1D
{
    Background, Signal
};

enum            fitresults2D
{

    BackgBackg, BackgSignl, SignlBackg, SignlSignl
};

void            SetBoundaries   ( string fOption, Double_t &aValMin, Double_t &aValMax )
{
    aValMin = 0.99;
    aValMax = 1.05;
    if ( fOption.find("RA") != -1 )
    {
        aValMin =   0.990;
        aValMax =   1.040;
    }
    if ( fOption.find("RB") != -1 )
    {
        aValMin =   0.990;
        aValMax =   1.060;
    }
    if ( fOption.find("RC") != -1 )
    {
        aValMin =   0.990;
        aValMax =   1.070;
    }
    if ( fOption.find("RD") != -1 )
    {
        aValMin =   0.995;
        aValMax =   1.040;
    }
    if ( fOption.find("RE") != -1 )
    {
        aValMin =   0.995;
        aValMax =   1.050;
    }
    if ( fOption.find("RF") != -1 )
    {
        aValMin =   0.995;
        aValMax =   1.060;
    }
    if ( fOption.find("RG") != -1 )
    {
        aValMin =   0.995;
        aValMax =   1.070;
    }
    if ( fOption.find("RH") != -1 )
    {
        aValMin =   1.000;
        aValMax =   1.040;
    }
    if ( fOption.find("RI") != -1 )
    {
        aValMin =   1.000;
        aValMax =   1.050;
    }
    if ( fOption.find("RJ") != -1 )
    {
        aValMin =   1.000;
        aValMax =   1.060;
    }
    if ( fOption.find("RK") != -1 )
    {
        aValMin =   1.000;
        aValMax =   1.070;
    }
}

RooFitResult*   FitModel        ( TH1F * THdata, const char* fName = "", Bool_t fSaveToFile = false, Int_t PTindex = -1, Int_t PTDimension = 1, string fOption = "" )
{
    // Silencing TCanvas Pop-Up
    gROOT->SetBatch();
    
    Bool_t  fCheb3  =   false;
    if  ( fOption.find("CH3") != -1 )   fCheb3  =   true;
    
    Bool_t  fCheb5  =   false;
    if  ( fOption.find("CH5") != -1 )   fCheb5  =   true;
    
    Bool_t  fWidth  =   false;
    if  ( fOption.find("W") != -1 )     fWidth  =   true;
    
    Bool_t  fMass_  =   false;
    if  ( fOption.find("M") != -1 )     fMass_  =   true;
    
    Double_t fInvMassValMax, fInvMassValMin;
    SetBoundaries(fOption,fInvMassValMin,fInvMassValMax);
    
    // Global Variables
    Int_t nEntries      = THdata->GetEntries();
    RooRealVar InvMass  = RooRealVar        ("InvMass","InvMass",fInvMassValMin,fInvMassValMax);
    RooDataHist* data   = new RooDataHist   (fName,fName,InvMass,Import(*THdata));
    Int_t kNCycle       = 5;
    
    // Background PDF Coefficients
    RooRealVar ch0      = RooRealVar        ("ch0","ch0"      ,0.5,-1,1);//,0.5,-1,1);
    RooRealVar ch1      = RooRealVar        ("ch1","ch1"      ,-0.1,-1,1);//,-0.1,-1,1);
    RooRealVar ch2      = RooRealVar        ("ch2","ch2"      ,0.01,-1,1);//,0.01,-1,1);
    RooRealVar ch3      = RooRealVar        ("ch3","ch3"      ,-0.05,-1,1);//,-0.05,-1,1);
    
    RooRealVar ch4, ch5;
    if ( fCheb3 && !fCheb5 )    ch4     = RooRealVar        ("ch4","ch4"        ,0.);
    else                        ch4     = RooRealVar        ("ch4","ch4"        ,0.2,-1,1);
    
    if ( fCheb5 )               ch5     = RooRealVar        ("ch5","ch5"        ,0.,-1,1);
    else                        ch5     = RooRealVar        ("ch5","ch5"        ,0.);
        
    //Signal
    RooRealVar sMass, sWidt, sSlop;
    if ( fWidth )               sWidt   = RooRealVar        ("sWidt","sWidt"    ,kPWid);
    else                        sWidt   = RooRealVar        ("sWidt","sWidt"    ,kPWid,kPWid*0.99,kPWid*1.01);
    
    if ( fMass_ )               sMass   = RooRealVar        ("sMass","sMass"    ,kPMas);
    else                        sMass   = RooRealVar        ("sMass","sMass"    ,kPMas,kPMas*0.99,kPMas*1.01);
    
    if ( bPythiaTest )          sSlop   = RooRealVar        ("sSlop","sSlop"    ,0.);
    else                        sSlop   = RooRealVar        ("sSlop","sSlop"    ,0.5,0.,1.);
    
    // Coefficients
    RooRealVar nSS      = RooRealVar        ("1nSS","1nSS"      ,0.5*nEntries,0.,nEntries);
    RooRealVar nBB      = RooRealVar        ("1nBB","1nBB"      ,0.5*nEntries,0.,nEntries);
    
    // PDFs
    RooVoigtian     fSig= RooVoigtian      ("fSig","fSig"      ,InvMass,sMass,sWidt,sSlop);
    RooChebychev    fBkg= RooChebychev     ("fBkg","fBkg"      ,InvMass,RooArgSet(ch0,ch1,ch2,ch3,ch4,ch5));
    RooAddPdf       fMod= RooAddPdf        ("fMod","fMod"      ,RooArgList(fBkg,fSig),RooArgList(nBB,nSS));
    
    
    RooFitResult* result;
    for ( Int_t iCycle = 0; iCycle < kNCycle; iCycle++ )
    {
        result = fMod.fitTo(*data,Extended(kTRUE),SumW2Error(kTRUE),Save(),Range(""),NormRange(""));
    }
    
    if ( fSaveToFile )
    {
        hName                       = "InvMass";
        hTitle                      = "Invariant Mass of Kaons in pT 0-6 GeV";
        if ( PTindex != -1 && !bPythiaTest )
        {
            if ( PTDimension == 1 ) hTitle = Form("Invariant Mass of Kaons in pT %.1f-%.1f GeV",fArrPT1D[PTindex],fArrPT1D[PTindex+1]);
            if ( PTDimension == 2 ) hTitle = Form("Invariant Mass of Kaons in pT %.1f-%.1f GeV",fArrPT2D[PTindex],fArrPT2D[PTindex+1]);
        }
        if ( PTindex != -1 &&  bPythiaTest )
        {
            if ( PTDimension == 1 ) hTitle = Form("Invariant Mass of Kaons in pT %.1f-%.1f GeV for MC",fArrPT1D[PTindex],fArrPT1D[PTindex+1]);
            if ( PTDimension == 2 ) hTitle = Form("Invariant Mass of Kaons in pT %.1f-%.1f GeV for MC",fArrPT2D[PTindex],fArrPT2D[PTindex+1]);
        }
        
        TCanvas * fSaveToCanvas;
        if ( PTDimension == 1 )fSaveToCanvas   =   new TCanvas(
                                                Form("PT_%.1f_%.1f_1D_%s",fArrPT1D[PTindex],fArrPT1D[PTindex+1],fName),
                                                Form("PT_%.1f_%.1f_1D_%s",fArrPT1D[PTindex],fArrPT1D[PTindex+1],fName)
                                                );
        
        if ( PTDimension == 2 )fSaveToCanvas   =   new TCanvas(
                                                Form("PT_%.1f_%.1f_2D_%s",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fName),
                                                Form("PT_%.1f_%.1f_2D_%s",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fName)
                                                );
        
        RooPlot * fSaveToFrame      = InvMass.frame(Name(hName),Title(hTitle));
        TLegend * fLegend           = new TLegend   (0.12,0.60,0.45,0.85);
        
        fRooPlotMaker(fSaveToFrame,fLegend,fMod,data,"InvMass1D");
        
        fSaveToFrame                ->Draw("same");
        fLegend                     ->Draw("same");
        fSaveToCanvas               ->Write ();
        delete fSaveToCanvas;
    }
    
    // Un-Silencing TCanvas Pop-Up
    gROOT->SetBatch(false);
    
    return result;
}

RooFitResult*   FitModel        ( TH2F * THdata, RooFitResult * fFitShapeX, RooFitResult * fFitShapeY, string fHistName = "", Bool_t fSaveToFile = false, Int_t PTindex = -1, Int_t PTjndex = -1, string fOption = "" )
{
    // Silencing TCanvas Pop-Up
    gROOT->SetBatch();
    
    Bool_t  fBackg  =   false;
    if  ( fOption.find("BK") != -1 )     fBackg  =   true;
    
    Double_t fInvMassValMax, fInvMassValMin;
    SetBoundaries(fOption,fInvMassValMin,fInvMassValMax);
    
    // Global Variables
    Int_t nEntries      = THdata->GetEntries();
    RooRealVar varx     = RooRealVar        ("xInvMass2D","xInvMass2D",fInvMassValMin,fInvMassValMax);
    RooRealVar vary     = RooRealVar        ("yInvMass2D","yInvMass2D",fInvMassValMin,fInvMassValMax);
    RooDataHist* data   = new RooDataHist   (fHistName.c_str(),fHistName.c_str(),RooArgList(varx,vary),Import(*THdata));
    Int_t kNCycle       = 5;
    
    RooArgSet  *utilityx    =   new RooArgSet(fFitShapeX->floatParsFinal(),fFitShapeX->constPars());
    RooArgSet  *utilityy    =   new RooArgSet(fFitShapeY->floatParsFinal(),fFitShapeY->constPars());
    
    
    // Background
    RooRealVar ch0x, ch1x, ch2x, ch3x, ch4x, ch5x, ch0y, ch1y, ch2y, ch3y, ch4y, ch5y;
    
    if ( fBackg )
    {
        ch0x     = RooRealVar ("ch0x","ch0x"     ,utilityx->getRealValue("ch0",0),utilityx->getRealValue("ch0",0)*0.75,utilityx->getRealValue("ch0",0)*1.25);
        ch1x     = RooRealVar ("ch1x","ch1x"     ,utilityx->getRealValue("ch1",0),utilityx->getRealValue("ch1",0)*0.75,utilityx->getRealValue("ch1",0)*1.25);
        ch2x     = RooRealVar ("ch2x","ch2x"     ,utilityx->getRealValue("ch2",0),utilityx->getRealValue("ch2",0)*0.75,utilityx->getRealValue("ch2",0)*1.25);
        ch3x     = RooRealVar ("ch3x","ch3x"     ,utilityx->getRealValue("ch3",0),utilityx->getRealValue("ch3",0)*0.75,utilityx->getRealValue("ch3",0)*1.25);
        ch4x     = RooRealVar ("ch4x","ch4x"     ,utilityx->getRealValue("ch4",0),utilityx->getRealValue("ch4",0)*0.75,utilityx->getRealValue("ch4",0)*1.25);
        ch5x     = RooRealVar ("ch5x","ch5x"     ,utilityx->getRealValue("ch5",0),utilityx->getRealValue("ch5",0)*0.75,utilityx->getRealValue("ch5",0)*1.25);
        ch0y     = RooRealVar ("ch0y","ch0y"     ,utilityy->getRealValue("ch0",0),utilityy->getRealValue("ch0",0)*0.75,utilityy->getRealValue("ch0",0)*1.25);
        ch1y     = RooRealVar ("ch1y","ch1y"     ,utilityy->getRealValue("ch1",0),utilityy->getRealValue("ch1",0)*0.75,utilityy->getRealValue("ch1",0)*1.25);
        ch2y     = RooRealVar ("ch2y","ch2y"     ,utilityy->getRealValue("ch2",0),utilityy->getRealValue("ch2",0)*0.75,utilityy->getRealValue("ch2",0)*1.25);
        ch3y     = RooRealVar ("ch3y","ch3y"     ,utilityy->getRealValue("ch3",0),utilityy->getRealValue("ch3",0)*0.75,utilityy->getRealValue("ch3",0)*1.25);
        ch4y     = RooRealVar ("ch4y","ch4y"     ,utilityy->getRealValue("ch4",0),utilityy->getRealValue("ch4",0)*0.75,utilityy->getRealValue("ch4",0)*1.25);
        ch5y     = RooRealVar ("ch5y","ch5y"     ,utilityy->getRealValue("ch5",0),utilityy->getRealValue("ch5",0)*0.75,utilityy->getRealValue("ch5",0)*1.25);
    }
    else
    {
        ch0x     = RooRealVar ("ch0x","ch0x"     ,utilityx->getRealValue("ch0",0));
        ch1x     = RooRealVar ("ch1x","ch1x"     ,utilityx->getRealValue("ch1",0));
        ch2x     = RooRealVar ("ch2x","ch2x"     ,utilityx->getRealValue("ch2",0));
        ch3x     = RooRealVar ("ch3x","ch3x"     ,utilityx->getRealValue("ch3",0));
        ch4x     = RooRealVar ("ch4x","ch4x"     ,utilityx->getRealValue("ch4",0));
        ch5x     = RooRealVar ("ch5x","ch5x"     ,utilityx->getRealValue("ch5",0));
        ch0y     = RooRealVar ("ch0y","ch0y"     ,utilityy->getRealValue("ch0",0));
        ch1y     = RooRealVar ("ch1y","ch1y"     ,utilityy->getRealValue("ch1",0));
        ch2y     = RooRealVar ("ch2y","ch2y"     ,utilityy->getRealValue("ch2",0));
        ch3y     = RooRealVar ("ch3y","ch3y"     ,utilityy->getRealValue("ch3",0));
        ch4y     = RooRealVar ("ch4y","ch4y"     ,utilityy->getRealValue("ch4",0));
        ch5y     = RooRealVar ("ch5y","ch5y"     ,utilityy->getRealValue("ch5",0));
    }
    
    //Signal
    RooRealVar pMassx   = RooRealVar ("pMassx","pMassx" ,utilityx->getRealValue("sMass",0));
    RooRealVar pWidthx  = RooRealVar ("pWidtx","pWidtx" ,utilityx->getRealValue("sWidt",0));
    RooRealVar pSlopex  = RooRealVar ("pSlopx","pSlopx" ,utilityx->getRealValue("sSlop",0));
    RooRealVar pMassy   = RooRealVar ("pMassy","pMassy" ,utilityy->getRealValue("sMass",0));
    RooRealVar pWidthy  = RooRealVar ("pWidty","pWidty" ,utilityy->getRealValue("sWidt",0));
    RooRealVar pSlopey  = RooRealVar ("pSlopy","pSlopy" ,utilityy->getRealValue("sSlop",0));
    
    // Coefficients
    auto fS__x  =   utilityx->getRealValue("1nSS",0);
    auto fB__x  =   utilityx->getRealValue("1nBB",0);
    auto fS__y  =   utilityy->getRealValue("1nSS",0);
    auto fB__y  =   utilityy->getRealValue("1nBB",0);
    auto fTot_  =   fS__x*fS__y+fB__x*fB__y+fS__x*fB__y+fB__x*fS__y;
    RooRealVar n0       = RooRealVar ("anSS2D","anSS2D" ,nEntries*(fS__x*fS__y)/fTot_,0.,nEntries);
    RooRealVar n1       = RooRealVar ("anBB2D","anBB2D" ,nEntries*(fB__x*fB__y)/fTot_,0.,nEntries);
    RooRealVar n2       = RooRealVar ("anBS2D","anBS2D" ,nEntries*(fB__x*fS__y)/fTot_,0.,nEntries);
    RooRealVar n3       = RooRealVar ("anSB2D","anSB2D" ,nEntries*(fS__x*fB__y)/fTot_,0.,nEntries);
    
    // PDFs
    RooChebychev        fBkgx ("fBkgx","fBkgx"          ,varx,RooArgSet(ch0x,ch1x,ch2x,ch3x,ch4x,ch5x));
    RooVoigtian         fSigx ("fSigx","fSigx"          ,varx,pMassx,pWidthx,pSlopex);
    RooChebychev        fBkgy ("fBkgy","fBkgy"          ,vary,RooArgSet(ch0y,ch1y,ch2y,ch3y,ch4y,ch5x));
    RooVoigtian         fSigy ("fSigy","fSigy"          ,vary,pMassy,pWidthy,pSlopey);
    RooProdPdf          fBB   ("fBkg","fBkg"            ,fBkgx,fBkgy);
    RooProdPdf          fSB   ("fSigBkg","fSBWBkg"      ,fSigx,fBkgy);
    RooProdPdf          fBS   ("fBkgSig","fBkgSig"      ,fBkgx,fSigy);
    RooProdPdf          fSS   ("fSigSig","fSigSig"      ,fSigx,fSigy);
    RooAddPdf           fMod  ("fMod2D","fMod2D"        ,RooArgList(fBB,fSS,fSB,fBS),RooArgList(n1,n0,n3,n2));
    

    RooFitResult* FitResults;
    for ( Int_t iCycle = 0; iCycle < kNCycle; iCycle++ )
    {
       FitResults = fMod.fitTo(*data,Extended(kTRUE),SumW2Error(kTRUE),Save());
    }
    
    // Save to file
    if ( fSaveToFile )
    {
        int         nBinsPrint      =   3;
        double      dIncrement      =   (fMaxIM2D-fMinIM2D)/nBinsPrint;
        TLatex*     latext          =   new TLatex();
        TCanvas*    cTotal          =   new TCanvas("","",0,45,1440,855);
                    cTotal          ->  SetTitle(Form("Slices of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]));
                    cTotal          ->  SetName(Form("PT_%.1f_%.1f__%.1f_%.1f_%s",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1],fHistName.c_str()));
                    cTotal          ->  Divide(nBinsPrint,2);
        
                            varx.setRange("fDrawRange",fMinIM2D,fMaxIM2D);
                            vary.setRange("fDrawRange",fMinIM2D,fMaxIM2D);
        for (int i = 0; i < nBinsPrint; i++)
        {
            hName                       = "Slice of 2D Invariant Mass of Kaons";
            hTitle                      = "Slice of 2D Invariant Mass of Kaons";
            if ( PTindex != -1 && !bPythiaTest ) hTitle = Form("Slice of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]);
            if ( PTindex != -1 &&  bPythiaTest ) hTitle = Form("Slice of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV for MC",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]);
            
            TCanvas * fSaveToCanvas =   new TCanvas(
                                                    Form("xInvMass_%.3f_%.3f_PTx_%.3f_%.3f_PTy_%.3f_%.3f_%s",fMinIM2D+dIncrement*i,fMinIM2D+dIncrement*(i+1),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1],fHistName.c_str()),
                                                    Form("xInvMass_%.3f_%.3f_PTx_%.3f_%.3f_PTy_%.3f_%.3f",fMinIM2D+dIncrement*i,fMinIM2D+dIncrement*(i+1),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1])
                                                    );
            
            RooPlot * fSaveToFrame  =   vary.frame(Name(hName),Title(hTitle));
            TLegend * fLegend           = new TLegend   (0.12,0.60,0.45,0.85);

                            varx.setRange("fDrawRange",fMinIM2D+i*dIncrement,fMinIM2D+(i+1)*dIncrement);
                            vary.setRange("fDrawRange",fMinIM2D,fMaxIM2D);

            fRooPlotMaker(fSaveToFrame,fLegend,fMod,data,"yInvMass2D");
            
            cTotal->cd( i+1 );
            fSaveToFrame                ->Draw("same");
            fLegend                     ->Draw("same");
            latext                      ->DrawLatexNDC(0.5, 0.8, Form("%.3f < m^{x}_{K_{+}K_{-}} < %.3f",fMinIM2D+dIncrement*i,fMinIM2D+dIncrement*(i+1)));
            fSaveToCanvas->cd();
            fSaveToFrame                ->Draw("same");
            fLegend                     ->Draw("same");
            latext                      ->DrawLatexNDC(0.5, 0.8, Form("%.3f < m^{x}_{K_{+}K_{-}} < %.3f",fMinIM2D+dIncrement*i,fMinIM2D+dIncrement*(i+1)));
            fSaveToCanvas               ->Write();
            delete fSaveToCanvas;
        }
                                        varx.setRange("fDrawRange",fMinIM2D,fMaxIM2D);
                                        vary.setRange("fDrawRange",fMinIM2D,fMaxIM2D);
        for (int i = 0; i < nBinsPrint; i++)
        {
            hName                       = "Slice of 2D Invariant Mass of Kaons";
            hTitle                      = "Slice of 2D Invariant Mass of Kaons";
            if ( PTindex != -1 && !bPythiaTest ) hTitle = Form("Slice of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]);
            if ( PTindex != -1 &&  bPythiaTest ) hTitle = Form("Slice of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV for MC",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]);
            
            TCanvas * fSaveToCanvas =   new TCanvas(
                                                    Form("yInvMass_%.3f_%.3f_PTx_%.3f_%.3f_PTy_%.3f_%.3f_%s",fMinIM2D+dIncrement*i,fMinIM2D+dIncrement*(i+1),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1],fHistName.c_str()),
                                                    Form("yInvMass_%.3f_%.3f_PTx_%.3f_%.3f_PTy_%.3f_%.3f",fMinIM2D+dIncrement*i,fMinIM2D+dIncrement*(i+1),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1])
                                                    );
            
            RooPlot * fSaveToFrame      =   varx.frame(Name(hName),Title(hTitle));
            TLegend * fLegend           = new TLegend   (0.12,0.60,0.45,0.85);
            
                                        varx.setRange("fDrawRange",fMinIM2D,fMaxIM2D);
                                        vary.setRange("fDrawRange",fMinIM2D+i*dIncrement,fMinIM2D+(i+1)*dIncrement);
                                                                            
            fRooPlotMaker(fSaveToFrame,fLegend,fMod,data,"xInvMass2D");
            
            cTotal->cd( i+1 +3 );
            fSaveToFrame                ->Draw("same");
            fLegend                     ->Draw("same");
            latext                      ->DrawLatexNDC(0.5, 0.8, Form("%.2f < m^{y}_{K_{+}K_{-}} < %.2f",fMinIM2D+dIncrement*i,fMinIM2D+dIncrement*(i+1)));
            
            fSaveToCanvas->cd();
            fSaveToFrame                ->Draw("same");
            fLegend                     ->Draw("same");
            latext                      ->DrawLatexNDC(0.5, 0.8, Form("%.2f < m^{y}_{K_{+}K_{-}} < %.2f",fMinIM2D+dIncrement*i,fMinIM2D+dIncrement*(i+1)));
            fSaveToCanvas               ->Write();
            delete fSaveToCanvas;
        }
                                        varx.setRange("fDrawRange",fMinIM2D,fMaxIM2D);
                                        vary.setRange("fDrawRange",fMinIM2D,fMaxIM2D);
        cTotal ->Write();
        delete cTotal;
    }
    
    // Un-Silencing TCanvas Pop-Up
    gROOT->SetBatch(false);
    
    // Fit
    return FitResults;
}

//--------------------------------------------//
//  Final Yield Extraction and extrapolation  //
//--------------------------------------------//

// - // Errors and management

TGraphAsymmErrors * SetTGrAsymE1D   ( TH1  * hTarget )
{
    TGraphAsymmErrors * gUtility  =   new TGraphAsymmErrors(hTarget);
    gUtility->RemovePoint(0);
    gUtility->RemovePoint(0);
    gUtility->RemovePoint(0);
    gUtility->RemovePoint(0);
    TGraphAsymmErrors * gResult___  =   new TGraphAsymmErrors();
    return gResult___;
}

TGraphAsymmErrors * SetTGrAsymE2D   ( TH1  * hTarget )
{
    TGraphAsymmErrors * gResult___  =   new TGraphAsymmErrors(hTarget);
    gResult___->RemovePoint(0);
    gResult___->RemovePoint(0);
    return gResult___;
}

TGraphAsymmErrors * AddTGrAsymESQ   ( TGraphAsymmErrors * gTarget, TGraphAsymmErrors * gTargetAdd )
{
    TGraphAsymmErrors * gResult___  =   new TGraphAsymmErrors();
    
    // Running on all points
    for ( Int_t iPoint = 0; iPoint < gTarget->GetMaxSize(); iPoint++ )
    {
        // Generating TGraph Points
        Double_t        XPoint__    =   (gTarget->  GetX())     [iPoint];
        Double_t        XPointEl    =   (gTarget->  GetEXlow()) [iPoint];
        Double_t        XPointEh    =   (gTarget->  GetEXhigh())[iPoint];
        Double_t        YPoint__    =   (gTarget->  GetY())     [iPoint];
        Double_t        YPointEl    =   pow((gTarget->  GetEYlow()) [iPoint],2) + pow((gTargetAdd->  GetEYlow()) [iPoint],2);
        Double_t        YPointEh    =   pow((gTarget->  GetEYhigh())[iPoint],2) + pow((gTargetAdd->  GetEYhigh())[iPoint],2);
        
        // Setting Final values
        gResult___  ->  SetPointX       (iPoint-3,XPoint__);
        gResult___  ->  SetPointEXhigh  (iPoint-3,XPointEl);
        gResult___  ->  SetPointEXlow   (iPoint-3,XPointEh);
        gResult___  ->  SetPointY       (iPoint-3,YPoint__);
        gResult___  ->  SetPointEYhigh  (iPoint-3,sqrt(YPointEl));
        gResult___  ->  SetPointEYlow   (iPoint-3,sqrt(YPointEh));
    }
    return gResult___;
}

TGraphAsymmErrors * AddTGrAsymESM   ( TGraphAsymmErrors * gTarget, TGraphAsymmErrors * gTargetAdd )
{
    TGraphAsymmErrors * gResult___  =   new TGraphAsymmErrors();
    
    // Running on all points
    for ( Int_t iPoint = 0; iPoint < gTarget->GetMaxSize(); iPoint++ )
    {
        // Generating TGraph Points
        Double_t        XPoint__    =   (gTarget->  GetX())     [iPoint];
        Double_t        XPointEl    =   (gTarget->  GetEXlow()) [iPoint];
        Double_t        XPointEh    =   (gTarget->  GetEXhigh())[iPoint];
        Double_t        YPoint__    =   (gTarget->  GetY())     [iPoint];
        Double_t        YPointEl    =   (gTarget->  GetEYlow()) [iPoint] + (gTargetAdd->  GetEYlow()) [iPoint];
        Double_t        YPointEh    =   (gTarget->  GetEYhigh())[iPoint] + (gTargetAdd->  GetEYhigh())[iPoint];
        
        // Setting Final values
        gResult___  ->  SetPointX       (iPoint-3,XPoint__);
        gResult___  ->  SetPointEXhigh  (iPoint-3,XPointEl);
        gResult___  ->  SetPointEXlow   (iPoint-3,XPointEh);
        gResult___  ->  SetPointY       (iPoint-3,YPoint__);
        gResult___  ->  SetPointEYhigh  (iPoint-3,sqrt(YPointEl));
        gResult___  ->  SetPointEYlow   (iPoint-3,sqrt(YPointEh));
    }
    return gResult___;
}

void                SetUnCorrltdE   ( TGraphAsymmErrors * gTarget )
{
    TGraphAsymmErrors * gResult___  =   new TGraphAsymmErrors();
    
    Double_t*       FYPntErrH_  =   gTarget->GetEYhigh();
    Double_t*       FYPntErrL_  =   gTarget->GetEYlow();
    
    // Running on all points
    for ( Int_t iPoint = 0; iPoint < gTarget->GetMaxSize(); iPoint++ )
    {
        Double_t        YPointEl    =   pow(FYPntErrH_[iPoint],2);
        Double_t        YPointEh    =   pow(FYPntErrL_[iPoint],2);
        
        // Event Trigger Efficiency
        YPointEh                    +=  pow(kEventEfficienERP,2);
        YPointEl                    +=  pow(kEventEfficienERM,2);
        
        // Setting Final values
        gResult___  ->  SetPointEYhigh  (iPoint,YPointEl);
        gResult___  ->  SetPointEYlow   (iPoint,YPointEh);
    }
}

void                SetCorrelatdE   ( TGraphAsymmErrors * gTarget )
{
    TGraphAsymmErrors * gResult___  =   new TGraphAsymmErrors();
    
    Double_t*       FYPntErrH_  =   gTarget->GetEYhigh();
    Double_t*       FYPntErrL_  =   gTarget->GetEYlow();
    
    // Running on all points
    for ( Int_t iPoint = 0; iPoint < gTarget->GetMaxSize(); iPoint++ )
    {
        Double_t        YPointEl    =   pow(FYPntErrH_[iPoint],2);
        Double_t        YPointEh    =   pow(FYPntErrL_[iPoint],2);
        
        YPointEh                    +=  0;
        YPointEl                    +=  0;
        
        // Setting Final values
        gResult___  ->  SetPointEYhigh  (iPoint,YPointEl);
        gResult___  ->  SetPointEYlow   (iPoint,YPointEh);
    }
}

TGraphAsymmErrors * SetSystErrors   ( TH1  * hTarget, string fOption )
{
    TGraphAsymmErrors * gResult___  =   new TGraphAsymmErrors();
    
    Bool_t              fCorrelatd  =   false;
    Bool_t              fUnCorrltd  =   false;
    
    for ( Char_t iOption:fOption )
    {
        if ( fOption.empty() )      break;
        if ( iOption == 'C' )       fCorrelatd = kTRUE;
        if ( iOption == 'U' )       fUnCorrltd = kTRUE;
        if ( iOption == '1' )       gResult___  =   SetTGrAsymE1D( hTarget );
        if ( iOption == '2' )       gResult___  =   SetTGrAsymE2D( hTarget );
    }
    
    if  ( fCorrelatd )              SetCorrelatdE(gResult___);
    if  ( fUnCorrltd )              SetUnCorrltdE(gResult___);
    
    return  gResult___;
}

TGraphAsymmErrors * SetSystErrors   ( TH1F * hTarget )
{
    TGraphAsymmErrors * _Return = new TGraphAsymmErrors(hTarget);
    
    Int_t iPoint = 0;
    for ( Int_t iHisto = 0; iHisto < nBinPT1D; iHisto++ )
    {
        auto fBinContent    =   hTarget->GetBinContent  (iHisto+1);
        if ( fBinContent == 0 )
        {
            _Return ->  RemovePoint(0);
            continue;
        }
        auto fBinError__    =   hTarget->GetBinError    (iHisto+1);
        auto ErrorYhigh = fBinContent*sqrt(pow(fBinError__/fBinContent,2)+ pow(kRapidityInterERP/kRapidityInterval,2)+pow(kEventEfficienERP/kEventEfficiency_,2)+pow(kBranchingRatiERP/kBranchingRatio__,2)+pow(kVertexEfficieERP/kVertexEfficiency,2)+pow(kSignalExtractERP/kSignalExtraction,2)+pow(kTrackingEfficERP/kTrackingEfficien,2)+pow(kParticleIdentERP/kParticleIdentifi,2));
        auto ErrorYlow_ = fBinContent*sqrt(pow(fBinError__/fBinContent,2)+ pow(kRapidityInterERM/kRapidityInterval,2)+pow(kEventEfficienERM/kEventEfficiency_,2)+pow(kBranchingRatiERM/kBranchingRatio__,2)+pow(kVertexEfficieERM/kVertexEfficiency,2)+pow(kSignalExtractERM/kSignalExtraction,2)+pow(kTrackingEfficERM/kTrackingEfficien,2)+pow(kParticleIdentERM/kParticleIdentifi,2));
        _Return ->  SetPointEYhigh  (iPoint,ErrorYhigh);
        _Return ->  SetPointEYlow   (iPoint,ErrorYlow_);
        iPoint++;
    }
    
    return _Return;
}

TGraphAsymmErrors * SetSystErrors   ( TH1D * hTarget )
{
    TGraphAsymmErrors * _Return = new TGraphAsymmErrors(hTarget);
    
    Int_t iPoint = 0;
    for ( Int_t iHisto = 0; iHisto < nBinPT2D; iHisto++ )
    {
        auto fBinContent    =   hTarget->GetBinContent  (iHisto+1);
        if ( fBinContent == 0 )
        {
            _Return ->  RemovePoint(0);
            continue;
        }
        auto fBinError__    =   hTarget->GetBinError    (iHisto+1);
        auto ErrorYhigh = fBinContent*sqrt(pow(fBinError__/fBinContent,2)+ pow(kRapidityInterERP/kRapidityInterval,2)+pow(kEventEfficienERP/kEventEfficiency_,2)+pow(2*kBranchingRatiERP/kBranchingRatio__,2)+pow(kVertexEfficieERP/kVertexEfficiency,2)+pow(kSignalExtractERP/kSignalExtraction,2)+pow(kTrackingEfficERP/kTrackingEfficien,2)+pow(kParticleIdentERP/kParticleIdentifi,2));
        auto ErrorYlow_ = fBinContent*sqrt(pow(fBinError__/fBinContent,2)+ pow(kRapidityInterERM/kRapidityInterval,2)+pow(kEventEfficienERM/kEventEfficiency_,2)+pow(2*kBranchingRatiERM/kBranchingRatio__,2)+pow(kVertexEfficieERM/kVertexEfficiency,2)+pow(kSignalExtractERM/kSignalExtraction,2)+pow(kTrackingEfficERM/kTrackingEfficien,2)+pow(kParticleIdentERM/kParticleIdentifi,2));
        _Return ->  SetPointEYhigh  (iPoint,ErrorYhigh);
        _Return ->  SetPointEYlow   (iPoint,ErrorYlow_);
        iPoint++;
    }
    return _Return;
}

void                SetSystErr1D_   ( TGraphAsymmErrors * gTarget )
{
    Double_t*           FXPoints__  =   gTarget->GetX();
    Double_t*           FXPntErrH_  =   gTarget->GetEXhigh();
    Double_t*           FXPntErrL_  =   gTarget->GetEXlow();
    Double_t*           FYPoints__  =   gTarget->GetY();
    Double_t*           FYPntErrH_  =   gTarget->GetEYhigh();
    Double_t*           FYPntErrL_  =   gTarget->GetEYlow();
    
    for ( Int_t iHisto = 0; iHisto < gTarget->GetMaxSize(); iHisto++ )
    {
        Double_t    EYlow   =   (FYPntErrH_[iHisto]/FYPoints__[iHisto])+kRapidityInterERM+kEventEfficienERM+kBranchingRatiERM+kVertexEfficieERM+kSignalExtractERM+kTrackingEfficERM+kParticleIdentERM;
        Double_t    EYhig   =   (FYPntErrL_[iHisto]/FYPoints__[iHisto])+kRapidityInterERP+kEventEfficienERP+kBranchingRatiERP+kVertexEfficieERP+kSignalExtractERP+kTrackingEfficERP+kParticleIdentERP;
        
        gTarget->SetPointEYlow  (iHisto,EYlow*FYPoints__[iHisto]);
        gTarget->SetPointEYhigh (iHisto,EYhig*FYPoints__[iHisto]);
    }
}


TGraphAsymmErrors * SetRandPoints   ( TGraphAsymmErrors * gTarget )
{
    TGraphAsymmErrors * gResult___  =   new TGraphAsymmErrors();
                    
    Double_t*           FXPoints__  =   gTarget->GetX();
    Double_t*           FXPntErrH_  =   gTarget->GetEXhigh();
    Double_t*           FXPntErrL_  =   gTarget->GetEXlow();
    Double_t*           FYPoints__  =   gTarget->GetY();
    Double_t*           FYPntErrH_  =   gTarget->GetEYhigh();
    Double_t*           FYPntErrL_  =   gTarget->GetEYlow();
    
    // Running on all points
    for ( Int_t iPoint = 0; iPoint < gTarget->GetMaxSize(); iPoint++ )
    {
        gResult___  ->  SetPointX       (iPoint,FXPoints__[iPoint]);
        gResult___  ->  SetPointEXhigh  (iPoint,FXPntErrH_[iPoint]);
        gResult___  ->  SetPointEXlow   (iPoint,FXPntErrL_[iPoint]);
    
        // Generating Shuffling
        Double_t        LowBound    =   (FYPoints__[iPoint] - FYPntErrL_[iPoint]);
        Double_t        HigBound    =   (FYPoints__[iPoint] + FYPntErrH_[iPoint]);
        Double_t        MidPoint;
        Bool_t kCheck = true;
        while ( kCheck )
        {
            MidPoint    =   fRandomGen  ->  Gaus( FYPoints__[iPoint], HigBound );
            if ( MidPoint >= LowBound && MidPoint <= HigBound ) kCheck = false;
        }
        
        // Setting Final values
        gResult___  ->  SetPointY       (iPoint,MidPoint);
        gResult___  ->  SetPointEYhigh  (iPoint,HigBound - MidPoint);
        gResult___  ->  SetPointEYlow   (iPoint,MidPoint - LowBound);
    }
    return gResult___;
}

// - // Integrals

Float_t             TGrAEIntegral    ( TGraphAsymmErrors * gTarget, Float_t & _ErrorHigh, Float_t & _ErrorLow_, Int_t kMinBin = 0, Int_t kMaxBin = 0 )
{
    Float_t         _Return     = 0;
    Float_t         fErrorHigh  = 0;
    Float_t         fErrorLow_  = 0;
    
    Double_t*       FXPoints__  =   gTarget->GetX();
    Double_t*       FXPntErrH_  =   gTarget->GetEXhigh();
    Double_t*       FXPntErrL_  =   gTarget->GetEXlow();
    Double_t*       FYPoints__  =   gTarget->GetY();
    Double_t*       FYPntErrH_  =   gTarget->GetEYhigh();
    Double_t*       FYPntErrL_  =   gTarget->GetEYlow();
    
    // Running on requested points
    for ( Int_t iPoint = kMinBin; iPoint < kMaxBin; iPoint++ )
    {
        _Return     += FYPoints__[iPoint]*(FXPntErrH_[iPoint]+FXPntErrL_[iPoint]);
        fErrorHigh  += FYPntErrH_[iPoint]*FYPntErrH_[iPoint];
        fErrorLow_  += FYPntErrL_[iPoint]*FYPntErrL_[iPoint];
    }
    
    _ErrorHigh = sqrt(fErrorHigh);
    _ErrorLow_ = sqrt(fErrorLow_);
    return _Return;
}

// - // Fits

Double_t            fLevyFunc1D     ( Double_t * fVar, Double_t * fParams )
{
    Double_t    fPT     = fVar[0];
    Double_t    fMass   = fParams[0];
    Double_t    fEnne   = fParams[1];
    Double_t    fSlop   = fParams[2];
    Double_t    fdNdY   = fParams[3];
    
    Double_t    fNum1   = (fEnne-1)*(fEnne-2);
    Double_t    fDen1   = fEnne*fSlop*(fEnne*fSlop+fMass*(fEnne-2));
    Double_t    fFac1   = fNum1/fDen1;
    
    Double_t    fMasT   = sqrt(fMass*fMass+fPT*fPT);
    Double_t    fNum2   = fMasT - fMass;
    Double_t    fDen2   = fEnne*fSlop;
    Double_t    fFac2   = TMath::Power((1 + fNum2/fDen2),(-fEnne));
    
    return      fPT*fdNdY*fFac1*fFac2;
}

TF1 *               fLevyFit1D      = new TF1 ("fLevyFunc1D",fLevyFunc1D,fMinPT1D,fMaxPT1D,4);

void                SetLevyTsalPT1D ( )
{
    // - // Setting up Fit parameters
    
    // Mass
    fLevyFit1D  ->  SetParLimits(0,kPMas,kPMas);
    fLevyFit1D  ->  SetParameter(0,kPMas);
    
    // n-Parameter
    fLevyFit1D  ->  SetParLimits(1,6.0,7.4);
    fLevyFit1D  ->  SetParameter(1,6.7); // 6.7
    
    // T-Parameter
    fLevyFit1D  ->  SetParLimits(2,.25,.30);
    fLevyFit1D  ->  SetParameter(2,.272); // .272
    
    // dN/dy
    fLevyFit1D  ->  SetParLimits(3,0.028,0.036);
    fLevyFit1D  ->  SetParameter(3,0.032);
}
 
void                SetLevyTsalPT2D ( )
{
    // - // Setting up Fit parameters
    
    // Mass
    fLevyFit1D  ->  SetParLimits(0,kPMas,kPMas);
    fLevyFit1D  ->  SetParameter(0,kPMas);
    
    // n-Parameter
    fLevyFit1D  ->  SetParLimits(1,3.5,7.5);
    fLevyFit1D  ->  SetParameter(1,4.5); // 6.7
    
    // T-Parameter
    fLevyFit1D  ->  SetParLimits(2,.18,.45);
    fLevyFit1D  ->  SetParameter(2,.272); // .272
    
    // dN/dy
    fLevyFit1D  ->  SetParLimits(3,0.5e-6,1.e-3);
    fLevyFit1D  ->  SetParameter(3,1.e-6);
}


void                SetLevyTsalPT2D ( Int_t iSet )
{
    // - // Setting up Fit parameters
    
    // Mass
    fLevyFit1D  ->  SetParLimits(0,kPMas,kPMas);
    fLevyFit1D  ->  SetParameter(0,kPMas);
    
    // n-Parameter
    fLevyFit1D  ->  SetParLimits(1,3.5,7.5);
    fLevyFit1D  ->  SetParameter(1,fLevyPar1[iSet]); // 6.7
    
    // T-Parameter
    fLevyFit1D  ->  SetParLimits(2,.18,.45);
    fLevyFit1D  ->  SetParameter(2,fLevyPar2[iSet]); // .272
    
    // dN/dy
    fLevyFit1D  ->  SetParLimits(3,1.e-8,1.e-3);
    fLevyFit1D  ->  SetParameter(3,fLevyPar3[iSet]);
}

Float_t             LevyTsalPT1D    ( TGraphAsymmErrors * gTarget, Float_t & _Error, TFitResultPtr & _FitRst, Float_t kMinVal = 0, Float_t kMaxVal = 0, Bool_t fSaveToFile = false )
{
    Float_t _Return;
    Float_t Cycle1, Cycle2, Cycle3;
    _FitRst     =   gTarget     ->  Fit(fLevyFit1D,"MREQ0S","EX0",0.4,10.);
    
    if ( fSaveToFile )
    {
        TCanvas * cSave = new TCanvas();
        gStyle->SetOptStat(0);
        gPad->SetLogy();
        gTarget->Draw();
        fLevyFit1D->Draw("same");
        cSave->Write();
        delete cSave;
    }
    
    _Return     =   fLevyFit1D ->Integral(kMinVal,kMaxVal);
    _Error      =   fLevyFit1D ->IntegralError(kMinVal,kMaxVal);
    return _Return;
}

// Final Extraction

Float_t             EvlIntegral     ( TGraphAsymmErrors * gTarget, Float_t & _ErrorHigh, Float_t & _ErrorLow_, Int_t kMinBin = 0, Int_t kMaxBin = nBinPT1D-4, string fName = "", Int_t kFitPrep = -1 )
{
    // Speeding multiple fits
    gROOT->SetBatch();
    
    Float_t _Return, _Integral, _IntErrHig, _IntErrLow, _IntFit_, _IntErrFit;
    TFitResultPtr FUtility;
    
    // Integrating in the known region
    _Integral                   =   TGrAEIntegral(gTarget,_IntErrHig,_IntErrLow,kMinBin,kMaxBin);
    
    //Extrapolating from the Fit
    _IntFit_                    =   LevyTsalPT1D(gTarget,_IntErrFit,FUtility,0,0.4);
    
    Int_t kPoints   = 3e3;
    
    // - // Initialising the random point mover
    TGraphAsymmErrors*gUtility  =   new TGraphAsymmErrors();
    Double_t * xUtility = new Double_t [kPoints];
    Double_t * xUtilit1 = new Double_t [kPoints];
    Double_t * xUtilit2 = new Double_t [kPoints];
    Double_t * xUtilit3 = new Double_t [kPoints];
    Double_t * wUtility = new Double_t [kPoints];
    _IntFit_                    =   LevyTsalPT1D(gTarget,_IntErrFit,FUtility,0,0.4,true);
    
    for ( Int_t iTer = 0; iTer < kPoints; iTer++ )
    {
        // Inizializzare il fit sempre uguale
        if( kFitPrep < 0  )     SetLevyTsalPT1D();
        if( kFitPrep >= 0 )     SetLevyTsalPT2D(kFitPrep);
        gUtility                =   SetRandPoints(gTarget);
        SetSystErr1D_(gUtility);
        _IntFit_                =   LevyTsalPT1D(gUtility,_IntErrFit,FUtility,0,0.4);
        xUtilit1[iTer]          =   fLevyFit1D->GetParameter(1);
        xUtilit2[iTer]          =   fLevyFit1D->GetParameter(2);
        xUtilit3[iTer]          =   fLevyFit1D->GetParameter(3);
        xUtility[iTer]          =   _IntFit_;
        wUtility[iTer]          =   1;
        if ( iTer%(kPoints/10) == 0 )cout << "iTer:" << iTer << endl;
    }
    
    TVectorD fUtility (kPoints,xUtility);
    TVectorD fUtilit1 (kPoints,xUtilit1);
    TVectorD fUtilit2 (kPoints,xUtilit2);
    TVectorD fUtilit3 (kPoints,xUtilit3);
    
    TH1F * hUtility             =   new TH1F (Form("gs_%s",fName.c_str()),"",75,fUtility.Min(),fUtility.Max());
    hUtility->FillN(kPoints,xUtility,wUtility);
    
    TH1F * hUtilit1             =   new TH1F (Form("L1_%s",fName.c_str()),"",75,fUtilit1.Min(),fUtilit1.Max());
    hUtilit1->FillN(kPoints,xUtilit1,wUtility);
    
    TH1F * hUtilit2             =   new TH1F (Form("L2_%s",fName.c_str()),"",75,fUtilit2.Min(),fUtilit2.Max());
    hUtilit2->FillN(kPoints,xUtilit2,wUtility);
    
    TH1F * hUtilit3             =   new TH1F (Form("L3_%s",fName.c_str()),"",75,fUtilit3.Min(),fUtilit3.Max());
    hUtilit3->FillN(kPoints,xUtilit3,wUtility);
    
    hUtility->Fit("gaus","Q");
    _IntFit_                    =   hUtility->GetFunction("gaus")->GetParameter(1);
    _IntErrFit                  =   hUtility->GetFunction("gaus")->GetParameter(2);
    hUtility                    ->Write();
    hUtilit1                    ->Write();
    hUtilit2                    ->Write();
    hUtilit3                    ->Write();
    
    _Return                     =   _Integral + _IntFit_;
    _ErrorHigh                  =   sqrt( pow(_IntErrFit,2) + pow(_IntErrHig,2) );
    _ErrorLow_                  =   sqrt( pow(_IntErrFit,2) + pow(_IntErrLow,2) );
    
    // Speeding multiple fits
    gROOT->SetBatch(false);
    return _Return;
}

Float_t **          EvlIntegral     ( TH2F * hTarget, Float_t **& _ErrorHig, Float_t **& _ErrorLow )
{
    Float_t ** _Return__    = new Float_t *[2];
    _Return__[0]            = new Float_t [nBinPT2D+1];
    _Return__[1]            = new Float_t [nBinPT2D+1];
    _ErrorHig               = new Float_t *[2];
    _ErrorHig[0]            = new Float_t [nBinPT2D+1];
    _ErrorHig[1]            = new Float_t [nBinPT2D+1];
    _ErrorLow               = new Float_t *[2];
    _ErrorLow[0]            = new Float_t [nBinPT2D+1];
    _ErrorLow[1]            = new Float_t [nBinPT2D+1];
    TH1D *      hSliceFX_   = new TH1D("","",nBinPT2D,fArrPT2D);
    TH1D *      hSliceFY_   = new TH1D("","",nBinPT2D,fArrPT2D);
    TF1  **     fLevyMem_   = new TF1 * [2];
    fLevyMem_[0]            = new TF1 [nBinPT2D+1];
    fLevyMem_[1]            = new TF1 [nBinPT2D+1];
    
    TCanvas * cDrawAllX  = new TCanvas("cDrawAllX","cDrawAllX");
    cDrawAllX->Divide(5,2);
    TCanvas * cDrawAllY  = new TCanvas("cDrawAllY","cDrawAllY");
    cDrawAllY->Divide(5,2);
    TCanvas * cDrawCHXY  = new TCanvas("cDrawCHXY","cDrawCHXY");
    cDrawCHXY->Divide(5,2);
    
    Int_t iHisto = 0;
    for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ )
    {
        if ( fArrPT2D[iFit+1] <= 0.41 ) continue;
        
        // Prepping the Function
        SetLevyTsalPT2D(iHisto);
        
        // X-Projection
        hName = Form("XProjection_PT_%.1f_%.1f",fArrPT2D[iFit],fArrPT2D[iFit+1]);
        TGraphAsymmErrors * g1D_ResX    =   new TGraphAsymmErrors(static_cast<TH1D*>(hTarget->ProjectionX(hName,iFit+1,iFit+1)));
        g1D_ResX                        ->  RemovePoint(0);
        g1D_ResX                        ->  RemovePoint(0);
        g1D_ResX                        ->  SetTitle(Form("Slice in #phi_{1} P_{T} from %.1f to %.1f GeV",fArrPT2D[iFit],fArrPT2D[iFit+1]));
        g1D_ResX                        ->  GetXaxis()  ->  SetTitle("#phi_{2} P_{T} GeV");
        g1D_ResX                        ->  GetYaxis()  ->  SetTitle("#frac{d^{3}N #phi_{1} }{dydp_{T}d#phi_{2}}(GeV/c)^{-1}");
        
        _Return__[0][iHisto+1]  =   EvlIntegral     (g1D_ResX,_ErrorHig[0][iHisto+1],_ErrorLow[0][iHisto+1],0,nBinPT2D-2,Form("X_%i",iFit),iHisto);
        hSliceFY_               ->  SetBinContent   (iFit+1,static_cast<Double_t>(_Return__[0][iHisto+1]));
        hSliceFY_               ->  SetBinError     (iFit+1,static_cast<Double_t>(_ErrorHig[0][iHisto+1]));
        fLevyMem_[0][iHisto+1]  =   *fLevyFit1D;
        
        // Prepping the Function
        SetLevyTsalPT2D(iHisto);
        
        // Y Projection
        hName = Form("YProjection_PT_%.1f_%.1f",fArrPT2D[iFit],fArrPT2D[iFit+1]);
        TGraphAsymmErrors * g1D_ResY    =   new TGraphAsymmErrors(static_cast<TH1D*>(hTarget->ProjectionY(hName,iFit+1,iFit+1)));
        g1D_ResY                        ->  RemovePoint(0);
        g1D_ResY                        ->  RemovePoint(0);
        g1D_ResY                        ->  SetTitle(Form("Slice in #phi_{2} P_{T} from %.1f to %.1f GeV",fArrPT2D[iFit],fArrPT2D[iFit+1]));
        g1D_ResY                        ->  GetXaxis()  ->  SetTitle("#phi_{1} P_{T} GeV");
        g1D_ResY                        ->  GetYaxis()  ->  SetTitle("#frac{d^{3}N #phi_{2} }{dydp_{T}d#phi_{1}}(GeV/c)^{-1}");
        _Return__[1][iHisto+1]  =   EvlIntegral     (g1D_ResY,_ErrorHig[1][iHisto+1],_ErrorLow[1][iHisto+1],0,nBinPT2D-2,Form("Y_%i",iFit),iHisto);
        hSliceFX_               ->  SetBinContent   (iFit+1,static_cast<Double_t>(_Return__[1][iHisto+1]));
        hSliceFX_               ->  SetBinError     (iFit+1,static_cast<Double_t>(_ErrorHig[1][iHisto+1]));
        fLevyMem_[1][iHisto+1]  =   *fLevyFit1D;
        
        cDrawAllX               ->  cd(iHisto+1);
        gPad                    ->  SetLogy();
        g1D_ResX                ->  Draw("AP");
        fLevyMem_[0][iHisto+1]  .   Draw("same");
    
        cDrawAllY               ->  cd(iHisto+1);
        gPad                    ->  SetLogy();
        g1D_ResY                ->  Draw("AP");
        fLevyMem_[1][iHisto+1]  .   Draw("same");
        
        
        cDrawCHXY               ->  cd(iHisto+1);
        cDrawCHXY               ->  SetLogy();
        hName       = Form("XProjection_PT_%.1f_%.1f",fArrPT2D[iFit],fArrPT2D[iFit+1]);
        TH1D * X__  = new TH1D(*static_cast<TH1D*>(hTarget->ProjectionX(hName,iFit+1,iFit+1)));
        hName = Form("YProjection_PT_%.1f_%.1f",fArrPT2D[iFit],fArrPT2D[iFit+1]);
        TH1D * Y__  = new TH1D(*static_cast<TH1D*>(hTarget->ProjectionY(hName,iFit+1,iFit+1)));
        X__ ->Draw("");
        Y__ ->Draw("same");
        /*
        hFitY[iHisto]->SetNameTitle(hName,hName);
        hFitY[iHisto]->GetYaxis()->SetTitle("#frac{d^{2}N_{#phi}}{dydp_{T}}(GeV/c)^{-1}");
        hFitY[iHisto]->GetXaxis()->SetTitleOffset(1.);
        hFitY[iHisto]->GetYaxis()->SetTitleOffset(1.4);
         */
        
        cout << "iHisto:" << iHisto << endl;
        
        iHisto++;
    }
    
    SetLevyTsalPT2D();
    TGraphAsymmErrors*      gSliceFX_       =   SetSystErrors(hSliceFX_);
    _Return__[0][0]     =   EvlIntegral     (gSliceFX_,_ErrorHig[0][0],_ErrorLow[0][0],3,nBinPT2D-2,Form("AllX"));
    fLevyMem_[0][0]     =   *fLevyFit1D;
    TGraphAsymmErrors*      gSliceFY_       =   SetSystErrors(hSliceFY_);
    _Return__[1][0]     =   EvlIntegral     (gSliceFX_,_ErrorHig[1][0],_ErrorLow[1][0],3,nBinPT2D-2,Form("AllY"));
    fLevyMem_[1][0]     =   *fLevyFit1D;
    
    
    TCanvas * cDrawComX  = new TCanvas("cDrawComX","cDrawComX");
    gPad                    ->  SetLogy();
    gSliceFX_               ->  Draw("");
    fLevyMem_[0][0]         .   Draw("same");
    
    TCanvas * cDrawComY  = new TCanvas("cDrawComY","cDrawComY");
    gPad                    ->  SetLogy();
    gSliceFY_               ->  Draw("");
    fLevyMem_[1][0]         .   Draw("same");
    
    cDrawAllX->Write();
    cDrawAllY->Write();
    cDrawCHXY->Write();
    return _Return__;
     
}

// Mean PT calculation

Float_t             EvlMeanPT__     ( TGraphAsymmErrors * gTarget, Float_t & _ErrorHigh, Float_t & _ErrorLow_, Int_t kMinBin = 0, Int_t kMaxBin = 1.e3 )
{
    TFitResultPtr   FUtilFit__;
    Float_t         FUtilError;
    Float_t         FTotWeight  =   0;
    Float_t         FTotError_  =   0;
    Float_t         FTotYPoint  =   0;
    Float_t         FFinalMean  =   0;
    Int_t           FkPointFit  =   2;
    Int_t           FkPointCmp  =   0;
    Double_t*       FFitErrors  =   new Double_t [FkPointFit];
    Double_t*       FFitXPoint  =   new Double_t [FkPointFit];
    Double_t*       FFitYPoint  =   new Double_t [FkPointFit];
    Double_t*       FXComplPnt  =   new Double_t [FkPointFit+gTarget->GetMaxSize()];
    Double_t*       FXComplPEh  =   new Double_t [FkPointFit+gTarget->GetMaxSize()];
    Double_t*       FXComplPEl  =   new Double_t [FkPointFit+gTarget->GetMaxSize()];
    Double_t*       FYComplPnt  =   new Double_t [FkPointFit+gTarget->GetMaxSize()];
    Double_t*       FYComplPEh  =   new Double_t [FkPointFit+gTarget->GetMaxSize()];
    Double_t*       FYComplPEl  =   new Double_t [FkPointFit+gTarget->GetMaxSize()];
    Double_t*       F_ComplWei  =   new Double_t [FkPointFit+gTarget->GetMaxSize()];
    Double_t*       F_ComplErr  =   new Double_t [FkPointFit+gTarget->GetMaxSize()];
    Double_t*       FXPoints__  =   gTarget->GetX();
    Double_t*       FXPntErrH_  =   gTarget->GetEXhigh();
    Double_t*       FXPntErrL_  =   gTarget->GetEXlow();
    Double_t*       FYPoints__  =   gTarget->GetY();
    Double_t*       FYPntErrH_  =   gTarget->GetEYhigh();
    Double_t*       FYPntErrL_  =   gTarget->GetEYlow();
    
    // Recovering Errors from Fit
    LevyTsalPT1D    (gTarget,FUtilError,FUtilFit__);
    
    // Initialising Fit Points
    for ( Int_t iPoint = 0; iPoint < FkPointFit; iPoint++ )
    {
        FFitXPoint[iPoint]      =   iPoint*(0.4)/(FkPointFit) + (0.4)/(2*FkPointFit);
        FFitYPoint[iPoint]      =   fLevyFit1D->Eval(FFitXPoint[iPoint]);
    }
    
    FUtilFit__      ->  GetConfidenceIntervals(FkPointFit, 1, 1, FFitXPoint, FFitErrors, 0.683, false);
    
    // Evaluating all weights
    for ( Int_t iPoint = 0; iPoint <= FkPointFit + gTarget->GetMaxSize() +10; iPoint++ )
    {
        
        if ( FXPoints__[iPoint-FkPointFit] <= FXPoints__[iPoint-1-FkPointFit] && iPoint > FkPointFit) break;
        
        if ( iPoint < FkPointFit )
        {
            FXComplPnt[iPoint]  =   FFitXPoint[iPoint];
            FXComplPEh[iPoint]  =   (0.4)/(2*FkPointFit);
            FXComplPEl[iPoint]  =   (0.4)/(2*FkPointFit);
            FYComplPnt[iPoint]  =   FFitYPoint[iPoint];
            FYComplPEh[iPoint]  =   FFitErrors[iPoint];
            FYComplPEl[iPoint]  =   FFitErrors[iPoint];
        }
        else
        {
            FXComplPnt[iPoint]  =   FXPoints__[iPoint-FkPointFit];
            FXComplPEh[iPoint]  =   FXPntErrH_[iPoint-FkPointFit];
            FXComplPEl[iPoint]  =   FXPntErrL_[iPoint-FkPointFit];
            FYComplPnt[iPoint]  =   FYPoints__[iPoint-FkPointFit];
            FYComplPEh[iPoint]  =   FYPntErrH_[iPoint-FkPointFit];
            FYComplPEl[iPoint]  =   FYPntErrL_[iPoint-FkPointFit];
        }
        FTotYPoint  +=  FYComplPnt[iPoint];
        FkPointCmp++;
    }
    
    for ( Int_t iPoint = 0; iPoint < FkPointCmp; iPoint++ )
    {
        auto    FPart1____  =   FXComplPnt[iPoint]*FXComplPnt[iPoint]/(FTotYPoint*FTotYPoint);
        auto    FPart2____  =   FYComplPnt[iPoint]*FYComplPnt[iPoint]*(FXComplPnt[iPoint]/FXComplPEh[iPoint])*(FXComplPnt[iPoint]/FXComplPEh[iPoint]);
        auto    FPart3____  =   (FTotYPoint-FYComplPnt[iPoint])*(FTotYPoint-FYComplPnt[iPoint]);
        auto    FPart4____  =   FYComplPEh[iPoint]*FYComplPEh[iPoint]/(FTotYPoint*FTotYPoint);
        F_ComplErr[iPoint]  =   FPart1____*(FPart3____*FPart4____);
        FTotError_          +=  F_ComplErr[iPoint];
        F_ComplWei[iPoint]  =   FYComplPnt[iPoint]*(FXComplPEl[iPoint]+FXComplPEh[iPoint]);
        FTotWeight          +=  F_ComplWei[iPoint];
    }
    
    for ( Int_t iPoint = 0; iPoint < FkPointCmp; iPoint++ )
    {
        FFinalMean  +=  F_ComplWei[iPoint]*FXComplPnt[iPoint]/FTotWeight;
    }
    
    _ErrorHigh  =   sqrt(FTotError_);
    _ErrorLow_  =   sqrt(FTotError_);
    return FFinalMean;
}

Float_t **          EvlMeanPT__     ( TH2F * hTarget, Float_t **& _ErrorHig, Float_t **& _ErrorLow )
{
    Float_t ** _Return__    = new Float_t *[2];
    _Return__[0]            = new Float_t [nBinPT2D+1];
    _Return__[1]            = new Float_t [nBinPT2D+1];
    _ErrorHig               = new Float_t *[2];
    _ErrorHig[0]            = new Float_t [nBinPT2D+1];
    _ErrorHig[1]            = new Float_t [nBinPT2D+1];
    _ErrorLow               = new Float_t *[2];
    _ErrorLow[0]            = new Float_t [nBinPT2D+1];
    _ErrorLow[1]            = new Float_t [nBinPT2D+1];
    TH1D *      hSliceFX_   = new TH1D("","",nBinPT2D,fArrPT2D);
    TH1D *      hSliceFY_   = new TH1D("","",nBinPT2D,fArrPT2D);
    TF1  **     fLevyMem_   = new TF1 * [2];
    fLevyMem_[0]            = new TF1 [nBinPT2D+1];
    fLevyMem_[1]            = new TF1 [nBinPT2D+1];
    
    TCanvas * cDrawAllX  = new TCanvas("cDrawAllX","cDrawAllX");
    cDrawAllX->Divide(3,3);
    TCanvas * cDrawAllY  = new TCanvas("cDrawAllY","cDrawAllY");
    cDrawAllY->Divide(3,3);
    TCanvas * cDrawCHXY  = new TCanvas("cDrawCHXY","cDrawCHXY");
    cDrawCHXY->Divide(3,3);
    
    Int_t iHisto = 0;
    for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ )
    {
        if ( fArrPT2D[iFit+1] <= 0.41 ) continue;
        
        cout << "Case: " << iHisto << endl;
        
        // Prepping the Function
        SetLevyTsalPT2D();
        
        // X-Projection
        hName = Form("XProjection_PT_%.1f_%.1f",fArrPT2D[iFit],fArrPT2D[iFit+1]);
        TGraphAsymmErrors * g1D_ResX = SetSystErrors(static_cast<TH1D*>(hTarget->ProjectionX(hName,iFit+1,iFit+1)));
        _Return__[0][iHisto]    =   EvlIntegral     (g1D_ResX,_ErrorHig[0][iHisto],_ErrorLow[0][iHisto]);
        hSliceFY_               ->  SetBinContent   (iFit+1,static_cast<Double_t>(_Return__[0][iHisto]));
        hSliceFY_               ->  SetBinContent   (iFit+1,static_cast<Double_t>(_ErrorHig[0][iHisto]));
        fLevyMem_[0][iHisto]    =   *fLevyFit1D;
        
        
        cout << "X- p3: " << fLevyMem_[0][iHisto].GetParameter(3) << " INT:  " << _Return__[0][iHisto] << endl;
        cout << "X- p3E:" << fLevyMem_[0][iHisto].GetParError(3) << " INTE: " << _ErrorHig[0][iHisto] << endl;
        
        // Prepping the Function
        SetLevyTsalPT2D();
        
        // Y Projection
        hName = Form("YProjection_PT_%.1f_%.1f",fArrPT2D[iFit],fArrPT2D[iFit+1]);
        TGraphAsymmErrors * g1D_ResY = SetSystErrors(static_cast<TH1D*>(hTarget->ProjectionY(hName,iFit+1,iFit+1)));
        _Return__[1][iHisto]    =   EvlIntegral     (g1D_ResY,_ErrorHig[1][iHisto],_ErrorLow[1][iHisto]);
        hSliceFX_               ->  SetBinContent   (iFit+1,static_cast<Double_t>(_Return__[1][iHisto]));
        hSliceFX_               ->  SetBinContent   (iFit+1,static_cast<Double_t>(_ErrorHig[1][iHisto]));
        fLevyMem_[1][iHisto]    =   *fLevyFit1D;
        
        cout << "Y- p3: " << fLevyMem_[1][iHisto].GetParameter(3) << " INT:  " << _Return__[1][iHisto] << endl;
        cout << "Y- p3E:" << fLevyMem_[1][iHisto].GetParError(3) << " INTE: " << _ErrorHig[1][iHisto] << endl;
        
        iHisto++;
        
        if ( iHisto+1 > 9 )     continue;
        cDrawAllX               ->  cd(iHisto+1);
        gPad                    ->  SetLogy();
        g1D_ResX                ->  Draw("");
        fLevyMem_[0][iHisto]    .   Draw("same");
    
        cDrawAllY               ->  cd(iHisto+1);
        gPad                    ->  SetLogy();
        g1D_ResY                ->  Draw("");
        fLevyMem_[1][iHisto]    .   Draw("same");
        
        cDrawCHXY               ->  cd(iHisto+1);
        cDrawCHXY               ->  SetLogy();
        static_cast<TH1D*>(hTarget->ProjectionX(hName,iFit+1,iFit+1))   ->Draw("same");
        static_cast<TH1D*>(hTarget->ProjectionY(hName,iFit+1,iFit+1))   ->Draw("same");
        
    }
    
    // Prepping the Function
    SetLevyTsalPT2D();
    _Return__[0][iHisto]    =   EvlIntegral     (SetSystErrors(hSliceFX_),_ErrorHig[0][iHisto],_ErrorLow[0][iHisto]);
    
    // Prepping the Function
    SetLevyTsalPT2D();
    _Return__[0][iHisto]    =   EvlIntegral     (SetSystErrors(hSliceFY_),_ErrorHig[0][iHisto],_ErrorLow[0][iHisto]);
    
    cDrawAllX->Write();
    cDrawAllY->Write();
    cDrawCHXY->Write();
    return _Return__;
     
}

void                TGrCompare1D    ( TGraphAsymmErrors * gDataStat, TH1F * h1D_Tru_P6, TH1F * h1D_Tru_P8, TGraphAsymmErrors * gCompare = nullptr )
{
    // Setting Canvas
    TCanvas *       cCompare1D  =   new TCanvas ("cCompare1D","cCompare1D");
    gPad->SetLogy();
    
    // Setting Standard Colors
    setMarker               (gCompare,3);
    setLine                 (h1D_Tru_P6,1);
    setLine                 (h1D_Tru_P8,2);
    
    // Setting Statistical error
    gDataStat->SetMarkerStyle(22);
    gDataStat->SetMarkerColor(38);
    gDataStat->SetFillColorAlpha(33,0.33);
    gDataStat->SetLineColorAlpha(33,0.66);
    gDataStat->GetYaxis()->SetRangeUser(4e-6,4e-2);
    gDataStat->GetXaxis()->SetRangeUser(0,1e1);
    
    // Setting Legend
    TLegend * lLegend1          =   new TLegend(0.65,0.65,0.85,0.85);
    lLegend1                    ->SetFillColor(kWhite);
    lLegend1                    ->SetLineColor(kWhite);
    lLegend1                    ->AddEntry(gDataStat,    "Data",             "P");
    lLegend1                    ->AddEntry(gDataStat,    "Statistical",      "EF");
    lLegend1                    ->AddEntry(gCompare,    "Previous Paper",   "EP");
    lLegend1                    ->AddEntry(h1D_Tru_P6,  "Pythia 6",         "EP");
    lLegend1                    ->AddEntry(h1D_Tru_P8,  "Pythia 8",         "EP");
    
    
    // Setting Multigraph
    TMultiGraph *   mCompare1D  =   new TMultiGraph();
    mCompare1D      ->Add   (gDataStat,"AP35");
    mCompare1D      ->Add   (gCompare,"P");
    
    // Drawing
    mCompare1D      ->Draw  ("AP");
    h1D_Tru_P6      ->Draw  ("same HIST L");
    h1D_Tru_P8      ->Draw  ("same HIST L");
    lLegend1        ->Draw  ("same");
    
    cCompare1D->Write();
}


//------------------------------//
//                              //
//      Trsnsv. Mom. FITs       //
//                              //
//------------------------------//

Double_t        fLevyFunc2D ( Double_t * fVar, Double_t * fParams )
{
    Double_t * fParx = new Double_t [4];
    Double_t * fPary = new Double_t [4];
    Double_t * fVarx = new Double_t [1];
    Double_t * fVary = new Double_t [1];
    
    fVarx = &fVar[0];
    fVary = &fVar[1];
    for ( int i = 0; i < 4; i++ )
    {
        fParx[i] = fParams[i];
        fPary[i] = fParams[i+4];
    }
    return  fPary[0]*(LevyTsallis_Func(fVarx,fParx)*LevyTsallis_Func(fVary,fParx));
}

Double_t        fFlatFuncXD ( Double_t * fVar, Double_t * fParams )
{
    return fParams[0];
}

TF2 * fLevyFit2D = new TF2 ("fLevyFunc2D",fLevyFunc2D,fMinPT2D,fMaxPT2D,fMinPT2D,fMaxPT2D,8);

TF1 * fFlatFit1D = new TF1 ("fLevyFunc1D",fFlatFuncXD,fMinPT1D,fMaxPT1D,1);

TF2 * fFlatFit2D = new TF2 ("fLevyFunc2D",fFlatFuncXD,fMinPT2D,fMaxPT2D,fMinPT2D,fMaxPT2D,1);

void fSliceCheck ( TH1F * hCheck1, TH1F * hCheck2, const char* hName1 = "1", const char* hName2 = "2", const char* fName = "-1", const char* fTitle = "" )
{
    TLatex * latext             =   new TLatex();
    
    TH1F hUtility1 = *hCheck1;
    TH1F hUtility2 = *hCheck2;
    hUtility1.SetMarkerStyle(29);
    hUtility1.SetMarkerColor(2);
    hUtility1.SetName("1");
    hUtility2.SetMarkerStyle(33);
    hUtility2.SetMarkerColor(4);
    hUtility2.SetName("2");
    
    TCanvas * fSaveToCanvas;
    fSaveToCanvas               =   new TCanvas(fName,fTitle);
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    hUtility1                   .Draw("same");
    hUtility2                   .Draw("same");
    TLegend * fLegend           = new TLegend   ();
    fLegend                     ->SetFillColor(kWhite);
    fLegend                     ->SetLineColor(kWhite);
    fLegend                     ->AddEntry("1",    hName1,      "EP");
    fLegend                     ->AddEntry("2",    hName2,      "EP");
    fLegend                     ->Draw("same");
    fSaveToCanvas->Write();
    
    delete fSaveToCanvas;
}

void fSliceCheck ( TH1D * hCheck1, TH1D * hCheck2, const char* hName1 = "1", const char* hName2 = "2", const char* fName = "-1", const char* fTitle = "")
{
    TLatex * latext             =   new TLatex();
    
    TH1D hUtility1 = *hCheck1;
    TH1D hUtility2 = *hCheck2;
    TH1D hUtility3 = *hCheck1;
    hUtility1.SetName("1");
    hUtility2.SetName("2");
    hUtility1.SetMarkerStyle(29);
    hUtility1.SetMarkerColor(2);
    hUtility3.SetMarkerStyle(29);
    hUtility3.SetMarkerColor(2);
    hUtility2.SetMarkerStyle(33);
    hUtility2.SetMarkerColor(4);
    
    TCanvas * fSaveToCanvas;
    fSaveToCanvas               =   new TCanvas(fName,fTitle);
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    hUtility1                   .Draw("same");
    hUtility2                   .Draw("same");
    TLegend * fLegend           = new TLegend   ();
    fLegend                     ->SetFillColor(kWhite);
    fLegend                     ->SetLineColor(kWhite);
    fLegend                     ->AddEntry("1",    hName1,      "EP");
    fLegend                     ->AddEntry("2",    hName2,      "EP");
    fLegend                     ->Draw("same");
    fSaveToCanvas->Write();
    
    delete fSaveToCanvas;
}

void fSliceCheckPT2D ( TH2F * hCheck1, TH2F * hCheck2, const char* hName1 = "1", const char* hName2 = "2" )
{
    TH1D * hUtility1X;
    TH1D * hUtility2X;
    TH1D * hUtility1Y;
    TH1D * hUtility2Y;
    for ( int iChk = 0; iChk < nBinPT2D; iChk++ )
    {
        // X-Projection
        hUtility1X  = (hCheck1->ProjectionX(Form("%s_X_Profile1_PT_%.2f_%.2f",hName1,fArrPT2D[iChk],fArrPT2D[iChk+1]),iChk+1,iChk+1));
        hUtility2X  = (hCheck2->ProjectionX(Form("%s_X_Profile2_PT_%.2f_%.2f",hName2,fArrPT2D[iChk],fArrPT2D[iChk+1]),iChk+1,iChk+1));
        hUtility1X->SetTitle(Form("X_%s",hName1));
        hUtility2X->SetTitle(Form("X_%s",hName2));
        fSliceCheck(hUtility1X,hUtility2X,hName1,hName2,Form("%s_%s_X_Profile_PT_%.2f_%.2f",hName1,hName2,fArrPT2D[iChk],fArrPT2D[iChk+1]),Form("X Profile in pT %.2f to %.2f GeV",fArrPT2D[iChk],fArrPT2D[iChk+1]));
                                         
        // Y-Projection
        hUtility1Y  = (hCheck1->ProjectionY(Form("%s_Y_Profile1_PT_%.2f_%.2f",hName1,fArrPT2D[iChk],fArrPT2D[iChk+1]),iChk+1,iChk+1));
        hUtility2Y  = (hCheck2->ProjectionY(Form("%s_Y_Profile2_PT_%.2f_%.2f",hName2,fArrPT2D[iChk],fArrPT2D[iChk+1]),iChk+1,iChk+1));
        hUtility1Y->SetTitle(Form("Y_%s",hName1));
        hUtility2Y->SetTitle(Form("Y_%s",hName2));
        fSliceCheck(hUtility1Y,hUtility2Y,hName1,hName2,Form("%s_%s_Y_Profile_PT_%.2f_%.2f",hName1,hName2,fArrPT2D[iChk],fArrPT2D[iChk+1]),Form("Y Profile in pT %.2f to %.2f GeV",fArrPT2D[iChk],fArrPT2D[iChk+1]));
    }
}

void fSlicePT2D ( TH2F * hCheck1, const char* hName1 = "1")
{
    TH1D * hUtility1X;
    TH1D * hUtility1Y;
    for ( int iChk = 0; iChk < nBinPT2D; iChk++ )
    {
        // X-Projection
        hUtility1X  = (hCheck1->ProjectionX(Form("%s_X_Profile1_PT_%.2f_%.2f",hName1,fArrPT2D[iChk],fArrPT2D[iChk+1]),iChk+1,iChk+1));
        hUtility1X->SetTitle(Form("X_%s",hName1));
        hUtility1X->SetMaximum(1.5);
        hUtility1X->SetMinimum(0.5);
        hUtility1X->Write();
                                         
        // Y-Projection
        hUtility1Y  = (hCheck1->ProjectionY(Form("%s_Y_Profile1_PT_%.2f_%.2f",hName1,fArrPT2D[iChk],fArrPT2D[iChk+1]),iChk+1,iChk+1));
        hUtility1Y->SetTitle(Form("Y_%s",hName1));
        hUtility1Y->SetMaximum(1.5);
        hUtility1Y->SetMinimum(0.5);
        hUtility1Y->Write();
    }
}

void fMosaicCanvas ( TH1F ** hData, const char*  fDTOpt = "", const char*  fName = "", Int_t xDiv = 1, Int_t yDiv = 1, Bool_t fSaveToFile = false, Bool_t fLogy = false, TF1 ** hFunc = nullptr, const char*  fFTOpt = "" )
{
    TCanvas * cResult = new TCanvas (fName,fName);
    cResult->Divide(xDiv,yDiv);
    for ( Int_t iHisto = 1; iHisto <= xDiv*yDiv; iHisto++ )
    {
        cResult->cd(iHisto);
        if ( fLogy ) gPad->SetLogy();
        
        // Data Draw
        if ( !hData ) continue;
        if ( !hData[iHisto-1] ) continue;
        hData[iHisto-1]->SetMarkerColor(1);
        hData[iHisto-1]->SetMarkerStyle(33);
        hData[iHisto-1]->Draw(fDTOpt);
        
        // Function Draw
        if ( !hFunc ) continue;
        if ( !hFunc[iHisto-1] ) continue;
        hFunc[iHisto-1]->SetLineColor(2);
        hFunc[iHisto-1]->Draw(Form("same %s",fFTOpt));
    }
    if ( fSaveToFile )
    {
        cResult->Write();
        cResult->SaveAs(Form("%s.pdf",fName));
        cResult->SaveAs(Form("%s.png",fName));
        delete cResult;
    }
}

void fMosaicCanvas ( TH1D ** hData, const char*  fDTOpt = "", const char*  fName = "", Int_t xDiv = 1, Int_t yDiv = 1, Bool_t fSaveToFile = false, Bool_t fLogy = false, TF1 ** hFunc = nullptr, const char*  fFTOpt = "" )
{
    TCanvas * cResult = new TCanvas (fName,fName);
    cResult->Divide(xDiv,yDiv);
    for ( Int_t iHisto = 1; iHisto <= xDiv*yDiv; iHisto++ )
    {
        cResult->cd(iHisto);
        if ( fLogy ) gPad->SetLogy();
        
        // Data Draw
        if ( !hData ) continue;
        if ( !hData[iHisto-1] ) continue;
        hData[iHisto-1]->SetMarkerColor(1);
        hData[iHisto-1]->SetMarkerStyle(33);
        hData[iHisto-1]->Draw(fDTOpt);
        
        // Function Draw
        if ( !hFunc ) continue;
        if ( !hFunc[iHisto-1] ) continue;
        hFunc[iHisto-1]->SetLineColor(2);
        hFunc[iHisto-1]->Draw(Form("same %s",fFTOpt));
    }
    if ( fSaveToFile )
    {
        cResult->Write();
        cResult->SaveAs(Form("%s.pdf",fName));
        cResult->SaveAs(Form("%s.png",fName));
        delete cResult;
    }
}

void fMosaicCanvas ( TH1F ** hData, const char*  fDTOpt = "", const char*  fName = "", Int_t xDiv = 1, Int_t yDiv = 1, Bool_t fSaveToFile = false, Bool_t fLogy = false, TH1F ** hDat2 = nullptr, const char*  fFTOpt = "" )
{
    TCanvas * cResult = new TCanvas (fName,fName);
    cResult->Divide(xDiv,yDiv);
    for ( Int_t iHisto = 1; iHisto <= xDiv*yDiv; iHisto++ )
    {
        cResult->cd(iHisto);
        if ( fLogy ) gPad->SetLogy();
        
        // Data Draw
        if ( !hData ) continue;
        if ( !hData[iHisto-1] ) continue;
        hData[iHisto-1]->SetMarkerColor(1);
        hData[iHisto-1]->SetMarkerStyle(33);
        hData[iHisto-1]->Draw(fDTOpt);
        
        // Data 2 Draw
        if ( !hDat2 ) continue;
        if ( !hDat2[iHisto-1] ) continue;
        hDat2[iHisto-1]->SetMarkerColor(2);
        hDat2[iHisto-1]->SetMarkerStyle(33);
        hDat2[iHisto-1]->Draw(Form("same %s",fFTOpt));
    }
    if ( fSaveToFile )
    {
        cResult->Write();
        cResult->SaveAs(Form("%s.pdf",fName));
        cResult->SaveAs(Form("%s.png",fName));
        delete cResult;
    }
}

void fMosaicCanvas ( TH1D ** hData, const char*  fDTOpt = "", const char*  fName = "", Int_t xDiv = 1, Int_t yDiv = 1, Bool_t fSaveToFile = false, Bool_t fLogy = false, TH1D ** hDat2 = nullptr, const char*  fFTOpt = "" )
{
    TCanvas * cResult = new TCanvas (fName,fName);
    cResult->Divide(xDiv,yDiv);
    for ( Int_t iHisto = 1; iHisto <= xDiv*yDiv; iHisto++ )
    {
        cResult->cd(iHisto);
        if ( fLogy ) gPad->SetLogy();
               
        // Data Draw
        if ( !hData ) continue;
        if ( !hData[iHisto-1] ) continue;
        hData[iHisto-1]->SetMarkerColor(1);
        hData[iHisto-1]->SetMarkerStyle(33);
        hData[iHisto-1]->Draw(fDTOpt);
               
        // Data 2 Draw
        if ( !hDat2 ) continue;
        if ( !hDat2[iHisto-1] ) continue;
        hDat2[iHisto-1]->SetMarkerColor(2);
        hDat2[iHisto-1]->SetMarkerStyle(33);
        hDat2[iHisto-1]->Draw(Form("same %s",fFTOpt));
    }
    if ( fSaveToFile )
    {
        cResult->Write();
        cResult->SaveAs(Form("%s.pdf",fName));
        cResult->SaveAs(Form("%s.png",fName));
        delete cResult;
    }
}

void fMosaicCanvas ( TH2F * hData, const char*  fDTOpt = "", const char*  fName = "", Int_t xDiv = 1, Int_t yDiv = 1, Bool_t fSaveToFile = false, Bool_t fLogy = false, Bool_t fPT2D = false, Double_t fMax = -1., Double_t fMin = -1. )
{
    // Building utility histograms
    TH1D ** hSlcX;
    if ( fPT2D )    hSlcX = new TH1D * [nBinPT2D];
    else            hSlcX = new TH1D * [xDiv*yDiv];
    TH1D ** hSlcY;
    if ( fPT2D )    hSlcY = new TH1D * [nBinPT2D];
    else            hSlcY = new TH1D * [xDiv*yDiv];
    TH1D ** fAmbg   = nullptr;
    
    auto iHist2 = 0;
    auto nTotal = 0;
    if ( fPT2D )    nTotal = nBinPT2D;
    else            nTotal = xDiv*yDiv;
    for ( Int_t iHisto = 0; iHisto < nTotal; iHisto++ )
    {
        if ( fPT2D && fArrPT2D[iHisto+1] <= 0.41 ) continue;
        
        // X-Projection
        if ( fPT2D )    hName = Form("XProjection_PT_%.1f_%.1f_%s",fArrPT2D[iHisto],fArrPT2D[iHisto+1],fName);
        else            hName = Form("XProjection_Bin_%i_%i_%s",iHisto,iHisto+1,fName);
        if ( fPT2D )    hSlcX[iHist2] = static_cast<TH1D*>(hData->ProjectionX("",iHisto+1,iHisto+1));
        else            hSlcX[iHisto] = static_cast<TH1D*>(hData->ProjectionX("",iHisto+1,iHisto+1));
        if ( fPT2D )    hSlcX[iHist2]->SetNameTitle(hName,hName);
        else            hSlcX[iHisto]->SetNameTitle(hName,hName);
        if ( fPT2D && fMax != -1. )  hSlcX[iHist2]->SetMaximum(fMax);
        if ( fPT2D && fMin != -1. )  hSlcX[iHist2]->SetMinimum(fMin);
        if ( !fPT2D&& fMax != -1. )  hSlcX[iHisto]->SetMaximum(fMax);
        if ( !fPT2D&& fMin != -1. )  hSlcX[iHisto]->SetMinimum(fMin);
        
        // Y Projection
        if ( fPT2D )    hName = Form("YProjection_PT_%.1f_%.1f_%s",fArrPT2D[iHisto],fArrPT2D[iHisto+1],fName);
        else            hName = Form("YProjection_Bin_%i_%i_%s",iHisto,iHisto+1,fName);
        if ( fPT2D )    hSlcY[iHist2] = static_cast<TH1D*>(hData->ProjectionY("",iHisto+1,iHisto+1));
        else            hSlcY[iHisto] = static_cast<TH1D*>(hData->ProjectionY("",iHisto+1,iHisto+1));
        if ( fPT2D )    hSlcY[iHist2]->SetNameTitle(hName,hName);
        else            hSlcY[iHisto]->SetNameTitle(hName,hName);
        if ( fPT2D && fMax != -1. )  hSlcY[iHist2]->SetMaximum(fMax);
        if ( fPT2D && fMin != -1. )  hSlcY[iHist2]->SetMinimum(fMin);
        if ( !fPT2D&& fMax != -1. )  hSlcY[iHisto]->SetMaximum(fMax);
        if ( !fPT2D&& fMin != -1. )  hSlcY[iHisto]->SetMinimum(fMin);
        
        iHist2++;
    }
    
    fMosaicCanvas(hSlcX,"",Form("XProjection_%s",fName),xDiv,yDiv,fSaveToFile,fLogy,fAmbg);
    fMosaicCanvas(hSlcY,"",Form("YProjection_%s",fName),xDiv,yDiv,fSaveToFile,fLogy,fAmbg);
    return;
}
/*
TF1 FitModelPT1D (TH1D * data, Int_t nEntries_DT, TF1 var, Bool_t fSaveToFile = false,  const char * fName = "" )
{
    // Scaling the data
    data    ->Scale(1./nEntries_DT);
    
    // Setting up the Fit
    fLevyFit1D->SetParLimits(0,kPMas,kPMas);
    fLevyFit1D->SetParameter(0,kPMas);
    fLevyFit1D->SetParLimits(1,2.+1.e-9,1.e1);
    fLevyFit1D->SetParameter(1,4.1);
    fLevyFit1D->SetParLimits(2,0.,1.);
    fLevyFit1D->SetParameter(2,0.3);
    fLevyFit1D->SetParLimits(3,1e-9,1.);
    fLevyFit1D->SetParameter(3,0.004);
    
    // Fitting
    data    ->Fit(fLevyFit1D,"IMREQ0","",0.4,10.);
    
    // Final Values
    double fInt, fErr;
    
    // Recovering the known yield
    fInt =  data                    ->IntegralAndError(3,nBinPT2D,fErr,"width");
    fErr *= fErr;
    fInt += fLevyFit1D              ->Integral(0.,0.4,1e-12);
    fErr += TMath::Power(fLevyFit1D  ->IntegralError(0.,0.4),2);
    
    // Save the results to visualise Fit goodness
    if ( fSaveToFile )
    {
        data->SetName(Form("THF_%s",fName));
        data->SetLineColor(kBlue);
        data->SetMarkerColor(kBlack);
        data->SetMarkerStyle(33);
        data->Write();
        fLevyFit1D->SetName(Form("FIT_%s",fName));
        fLevyFit1D->Write();
        TCanvas * c1 = new TCanvas(Form("CNV_%s",fName),fName);
        gPad->SetLogy();
        fLevyFit1D->SetRange(0.,10.);
        data->Draw("same");
        fLevyFit1D->Draw("same");
        c1->Write();
        delete c1;
    }
    
    var = *fLevyFit1D;
    
    // Results
    fLevyFit1D->SetParameter(3,fInt);
    fLevyFit1D->SetParError(3,sqrt(fErr));
    
    return var;
}

TH2F * HistoModel2D (RooFitResult * utilityx, RooFitResult * utilityy, RooFitResult * input, RooRealVar varx, RooRealVar vary, char *  hName)
{
    // Background
    RooRealVar ch0x     = RooRealVar ("ch0x","ch0x"     ,static_cast<RooRealVar*>(utilityx->floatParsFinal().at(ChebyPar0_))->getVal());
    RooRealVar ch1x     = RooRealVar ("ch1x","ch1x"     ,static_cast<RooRealVar*>(utilityx->floatParsFinal().at(ChebyPar1_))->getVal());
    RooRealVar ch2x     = RooRealVar ("ch2x","ch2x"     ,static_cast<RooRealVar*>(utilityx->floatParsFinal().at(ChebyPar2_))->getVal());
    RooRealVar ch3x     = RooRealVar ("ch3x","ch3x"     ,static_cast<RooRealVar*>(utilityx->floatParsFinal().at(ChebyPar3_))->getVal());
    RooRealVar ch4x     = RooRealVar ("ch4x","ch4x"     ,static_cast<RooRealVar*>(utilityx->floatParsFinal().at(ChebyPar4_))->getVal());
    RooRealVar ch0y     = RooRealVar ("ch0y","ch0y"     ,static_cast<RooRealVar*>(utilityy->floatParsFinal().at(ChebyPar0_))->getVal());
    RooRealVar ch1y     = RooRealVar ("ch1y","ch1y"     ,static_cast<RooRealVar*>(utilityy->floatParsFinal().at(ChebyPar1_))->getVal());
    RooRealVar ch2y     = RooRealVar ("ch2y","ch2y"     ,static_cast<RooRealVar*>(utilityy->floatParsFinal().at(ChebyPar2_))->getVal());
    RooRealVar ch3y     = RooRealVar ("ch3y","ch3y"     ,static_cast<RooRealVar*>(utilityy->floatParsFinal().at(ChebyPar3_))->getVal());
    RooRealVar ch4y     = RooRealVar ("ch4y","ch4y"     ,static_cast<RooRealVar*>(utilityy->floatParsFinal().at(ChebyPar4_))->getVal());
    
    //Signal
    RooRealVar pMassx   = RooRealVar ("pMassx","pMassx" ,static_cast<RooRealVar*>(utilityx->floatParsFinal().at(PhiMass_1D))->getVal());
    RooRealVar pWidthx  = RooRealVar ("pWidtx","pWidtx" ,static_cast<RooRealVar*>(utilityx->floatParsFinal().at(PhiWidth1D))->getVal());
    RooRealVar pSlopex;
    if ( !bPythiaTest )  pSlopex = RooRealVar ("pSlopx","pSlopx" ,static_cast<RooRealVar*>(utilityy->floatParsFinal().at(PhiSlope1D))->getVal());
    RooRealVar pMassy   = RooRealVar ("pMassy","pMassy" ,static_cast<RooRealVar*>(utilityy->floatParsFinal().at(PhiMass_1D))->getVal());
    RooRealVar pWidthy  = RooRealVar ("pWidty","pWidty" ,static_cast<RooRealVar*>(utilityy->floatParsFinal().at(PhiWidth1D))->getVal());
    RooRealVar pSlopey;
    if ( !bPythiaTest )  pSlopey = RooRealVar ("pSlopx","pSlopx" ,static_cast<RooRealVar*>(utilityy->floatParsFinal().at(PhiSlope1D))->getVal());
    
    // Coefficients
    RooRealVar n0       = RooRealVar ("nSS2D","nSS2D"   ,static_cast<RooRealVar*>(input->floatParsFinal().at(SignlSignl))->getVal());
    RooRealVar n1       = RooRealVar ("nBB2D","nBB2D"   ,static_cast<RooRealVar*>(input->floatParsFinal().at(BackgBackg))->getVal());
    RooRealVar n2       = RooRealVar ("nBS2D","nBS2D"   ,static_cast<RooRealVar*>(input->floatParsFinal().at(BackgSignl))->getVal());
    RooRealVar n3       = RooRealVar ("nSB2D","nSB2D"   ,static_cast<RooRealVar*>(input->floatParsFinal().at(SignlBackg))->getVal());
    
    // PDFs
    RooChebychev        fBkgx ("fBkgx","fBkgx"          ,varx,RooArgSet(ch0x,ch1x,ch2x,ch3x,ch4x));
    RooBreitWigner      fSigx ("fSigx","fSigx"          ,varx,pMassx,pWidthx);
    RooChebychev        fBkgy ("fBkgy","fBkgy"          ,vary,RooArgSet(ch0y,ch1y,ch2y,ch3y,ch4y));
    RooBreitWigner      fSigy ("fSigy","fSigy"          ,vary,pMassy,pWidthy);
    RooProdPdf          fBB   ("fBB2D","fBB2D"          ,fBkgx,fBkgy);
    RooProdPdf          fSB   ("fSB2D","fSB2D"          ,fSigx,fBkgy);
    RooProdPdf          fBS   ("fBS2D","fBS2D"          ,fBkgx,fSigy);
    RooProdPdf          fSS   ("fSS2D","fSS2D"          ,fSigx,fSigy);
    RooAddPdf           fMod  ("fMod2D","fMod2D"        ,RooArgList(fBB,fSS,fSB,fBS),RooArgList(n1,n0,n3,n2));
    
    auto return_histo = (TH2F *) fMod.createHistogram(hName,varx,Binning(nBinIM2D,fMinIM2D,fMaxIM2D),YVar(vary,Binning(nBinIM2D,fMinIM2D,fMaxIM2D)));
    return_histo->SetName(hName);
    return return_histo;
}
*/

 // Utility

void                SetLevyTsalis   ( )
{
    // - // Setting up Fit parameters
    
    // Mass
    fLevyFit1D  ->  SetParLimits(0,kPMas,kPMas);
    fLevyFit1D  ->  SetParameter(0,kPMas);
    
    // n-Parameter
    fLevyFit1D  ->  SetParLimits(1,2.1,7.5);
    fLevyFit1D  ->  SetParameter(1,6.); // 6.7
    
    // T-Parameter
    fLevyFit1D  ->  SetParLimits(2,.21,.750);
    fLevyFit1D  ->  SetParameter(2,.272); // .272
    
    // dN/dy
    fLevyFit1D  ->  SetParLimits(3,1.e-7,1.e-1);
    fLevyFit1D  ->  SetParameter(3,0.032);
    
    if ( bPythiaTest )
    {
        // Mass
        fLevyFit1D  ->  SetParLimits(0,kPMas,kPMas);
        fLevyFit1D  ->  SetParameter(0,kPMas);
           
        // n-Parameter
        fLevyFit1D  ->  SetParLimits(1,2.1,7.5);
        fLevyFit1D  ->  SetParameter(1,4.); // 6.7
           
        // T-Parameter
        fLevyFit1D  ->  SetParLimits(2,.01,.750);
        fLevyFit1D  ->  SetParameter(2,.21); // .272
           
        // dN/dy
        fLevyFit1D  ->  SetParLimits(3,1.e-9,1.);
        fLevyFit1D  ->  SetParameter(3,.04);
    }
}

Double_t            FuncIntegrals   ( TF1 * aFunction, Double_t aLowBound, Double_t aHigBound, string aOption = "", Double_t aEpsilon = 1.e-3 )
{
    Double_t    fResult = 0;
    Int_t       fiTer,  fnPowr;
    Bool_t      fWidth, fPower;
    
    if ( aOption.find("W") != -1 )  fWidth  =   true;
    else                                fWidth  =   false;
    
    if ( aOption.find("x^") != -1 )   { fPower  =   true;   fnPowr  =   aOption.at(aOption.find("x^")+2)-'0'; }
    else                                fPower  =   false;
    
    // Starting iterations at 0
    fiTer = 0;
    while ( true )
    {
        // Determining the pT point to calculate
        auto fPoint = aLowBound + aEpsilon*fiTer;
        
        // Exiting the loop if out of bound
        if ( fPoint >= aHigBound ) break;
        
        Double_t        fCycleAdd   =   (aFunction->Eval( fPoint ));
        
        if ( fWidth )   fCycleAdd  *=   ( aEpsilon );
        if ( fPower )   fCycleAdd  *=   pow( ( fPoint ) , fnPowr );
        
        fResult += fCycleAdd;
        fiTer++;
    }
    
    return fResult;
}

Double_t            HistIntegrals   ( TH1F* aHistogrm, string aOption = "", Double_t aLowBound = 0, Double_t aHigBound = 0 )
{
    Double_t    fResult = 0;
    Int_t       fnPowr;
    Bool_t      fWidth, fPower, fError;
    
    if ( aOption.find("W") != -1 )      fWidth  =   true;
    else                                fWidth  =   false;
    
    if ( aOption.find("x^") != -1 )   { fPower  =   true;   fnPowr  =   aOption.at(aOption.find("x^")+2)-'0'; }
    else                                fPower  =   false;
    
    if ( aOption.find("E") != -1 )      fError  =   true;
    else                                fError  =   false;
    
    // Starting iterations at 0
    for ( Int_t fiTer = 1; fiTer <= aHistogrm->GetNbinsX(); fiTer++ )
    {
        // Determining the pT point to calculate
        auto fPoint =   aHistogrm->GetBinCenter(fiTer);
        
        // Exiting the loop if out of bound
        if ( aHistogrm->GetBinLowEdge(fiTer)   <= aLowBound && aLowBound != aHigBound ) continue;
        if ( aHistogrm->GetBinLowEdge(fiTer+1) >= aHigBound && aLowBound != aHigBound ) break;
        
        Double_t        fCycleAdd   =   ( aHistogrm->GetBinContent(fiTer) );
        if ( fError )   fCycleAdd   =   ( aHistogrm->GetBinError(fiTer) );
        if ( fWidth )   fCycleAdd  *=   ( aHistogrm->GetBinWidth(fiTer) );
        if ( fPower )   fCycleAdd  *=   pow( ( aHistogrm->GetBinCenter(fiTer) ) , fnPowr );
        if ( fError )   fCycleAdd  *=   fCycleAdd;
        
        fResult += fCycleAdd;
    }
    
    if ( fError )   return sqrt(fResult);
    return fResult;
}

Double_t            HistIntegrals   ( TH1D* aHistogrm, string aOption = "", Double_t aLowBound = 0, Double_t aHigBound = 0 )
{
    
    Double_t    fResult = 0;
    Int_t       fnPowr;
    Bool_t      fWidth, fPower, fError;
    
    if ( aOption.find("W") != -1 )      fWidth  =   true;
    else                                fWidth  =   false;
    
    if ( aOption.find("x^") != -1 )   { fPower  =   true;   fnPowr  =   aOption.at(aOption.find("x^")+2)-'0'; }
    else                                fPower  =   false;
    
    if ( aOption.find("E") != -1 )      fError  =   true;
    else                                fError  =   false;
    
    // Starting iterations at 0
    for ( Int_t fiTer = 1; fiTer <= aHistogrm->GetNbinsX(); fiTer++ )
    {
        // Determining the pT point to calculate
        auto fPoint =   aHistogrm->GetBinCenter(fiTer);
        
        // Exiting the loop if out of bound
        if ( aHistogrm->GetBinLowEdge(fiTer)   <= aLowBound && aLowBound != aHigBound ) continue;
        if ( aHistogrm->GetBinLowEdge(fiTer+1) >= aHigBound && aLowBound != aHigBound ) break;
        
        Double_t        fCycleAdd   =   ( aHistogrm->GetBinContent(fiTer) );
        if ( fError )   fCycleAdd   =   ( aHistogrm->GetBinError(fiTer) );
        if ( fWidth )   fCycleAdd  *=   ( aHistogrm->GetBinWidth(fiTer) );
        if ( fPower )   fCycleAdd  *=   pow( ( aHistogrm->GetBinCenter(fiTer) ) , fnPowr );
        if ( fError )   fCycleAdd  *=   fCycleAdd;
        
        fResult += fCycleAdd;
    }
    
    if ( fError )   return sqrt(fResult);
    return fResult;
}

TH1F *              SetRandPoints   ( TH1F * aTarget )
{
    TH1F *  fReturn =   new TH1F (*aTarget);
    Double_t    fError,    fValue;
    for ( Int_t iBin = 1; iBin <= aTarget->GetNbinsX(); iBin++ )
    {
        fError  =   aTarget ->GetBinError   (iBin);
        fValue  =   aTarget ->GetBinContent (iBin);
        fReturn ->  SetBinContent           (iBin,   fValue + fRandomGen->Gaus(0.,fError) );
        fReturn ->  SetBinError             (iBin,   fError );
    }
    return fReturn;
}

TH1D *              SetRandPoints   ( TH1D * aTarget )
{
    TH1D *  fReturn =   new TH1D (*aTarget);
    Double_t    fError,    fValue;
    for ( Int_t iBin = 1; iBin <= aTarget->GetNbinsX(); iBin++ )
    {
        fError  =   aTarget ->GetBinError   (iBin);
        fValue  =   aTarget ->GetBinContent (iBin);
        fReturn ->  SetBinContent           (iBin,   fValue + fRandomGen->Gaus(0.,fError) );
        fReturn ->  SetBinError             (iBin,   fError );
    }
    return fReturn;
}

TH1F *              SetSystErrorsh  ( TH1F * aTarget ) // Togliere l'h
{
    TH1F *  fReturn =   new TH1F (*aTarget);
    Double_t    fError,    fValue;
    for ( Int_t iBin = 1; iBin <= aTarget->GetNbinsX(); iBin++ )
    {
        fError  =   aTarget ->GetBinError   (iBin);
        fValue  =   aTarget ->GetBinContent (iBin);
        fReturn ->  SetBinContent           (iBin,   fValue );
        fReturn ->  SetBinError             (iBin,   fError );//+ kSystematicalErrP*fValue );
    }
    return fReturn;
}

TH1D *              SetSystErrorsh  ( TH1D * aTarget ) // Togliere l'h
{
    TH1D *  fReturn =   new TH1D (*aTarget);
    Double_t    fError,    fValue;
    for ( Int_t iBin = 1; iBin <= aTarget->GetNbinsX(); iBin++ )
    {
        fError  =   aTarget ->GetBinError   (iBin);
        fValue  =   aTarget ->GetBinContent (iBin);
        fReturn ->  SetBinContent           (iBin,   fValue );
        fReturn ->  SetBinError             (iBin,   fError );//+ kSystematicalErrP*fValue );
    }
    return fReturn;
}

 // Measurements

Double_t *          EvalStErInteg   ( TH1F * aTarget,   string aName__ = "",  Bool_t aSaveFit = false )
{
    Int_t       fPoints =   1.e3;
    Double_t *  fReturn = new Double_t [2];
    Double_t *  fUtilPt = new Double_t [fPoints];
    Double_t *  fUtilMn = new Double_t [fPoints];
    Double_t *  fUtilP1 = new Double_t [fPoints];
    Double_t *  fUtilP2 = new Double_t [fPoints];
    Double_t *  fUtilP3 = new Double_t [fPoints];
    Double_t *  fUtilW0 = new Double_t [fPoints];
    
    // Speeding multiple fits
    gROOT->SetBatch();
    
    TH1F *  fRandPt,    *fTotal_;
    for ( Int_t iTer = 0; iTer < fPoints; iTer++ )
    {
        fRandPt =   SetRandPoints(aTarget);
        fTotal_ =   SetSystErrorsh(fRandPt);
        SetLevyTsalis();
        fRandPt->Fit(fLevyFit1D,"IMREQ0S","",0.4,10.);
        if ( bPythiaTest ) fRandPt->Fit(fLevyFit1D,"IMREQ0S","",0.4,1.6);
        fUtilPt[iTer]   =   fLevyFit1D->Moment(1,0.0,0.4);
        fUtilMn[iTer]   =   fLevyFit1D->Integral(0.0,0.4);
        fUtilP1[iTer]   =   fLevyFit1D->GetParameter(1);
        fUtilP2[iTer]   =   fLevyFit1D->GetParameter(2);
        fUtilP3[iTer]   =   fLevyFit1D->GetParameter(3);
        fUtilW0[iTer]   =   1;
        
        if ( aSaveFit )
        {
            TCanvas * c1    =   new TCanvas();
            gPad->SetLogy();
            fTotal_     ->Draw();
            fLevyFit1D  ->Draw("same");
            c1          ->Write();
            delete  c1;
        }
    }

    // Speeding multiple fits
    gROOT->SetBatch(false);

    TVectorD fUtilVP    (fPoints,fUtilPt);
    TVectorD fUtilVM    (fPoints,fUtilMn);
    TVectorD fUtilV1    (fPoints,fUtilP1);
    TVectorD fUtilV2    (fPoints,fUtilP2);
    TVectorD fUtilV3    (fPoints,fUtilP3);
    
    TH1F *  fUtilHP     =   new TH1F    (Form("gp_%s",aName__.c_str()),"",25,fUtilVP.Min(),fUtilVP.Max());
    fUtilHP             ->  FillN       (fPoints,fUtilPt,fUtilW0);
    fUtilHP             ->  Write();
    
    TH1F *  fUtilHM     =   new TH1F    (Form("gm_%s",aName__.c_str()),"",25,fUtilVM.Min(),fUtilVM.Max());
    fUtilHM             ->  FillN       (fPoints,fUtilMn,fUtilW0);
    fUtilHM             ->  Write();
    
    TH1F *  fUtilH1     =   new TH1F    (Form("P1_%s",aName__.c_str()),"",25,fUtilV1.Min(),fUtilV1.Max());
    fUtilH1             ->  FillN       (fPoints,fUtilP1,fUtilW0);
    fUtilH1             ->  Write();
    
    TH1F *  fUtilH2     =   new TH1F    (Form("P2_%s",aName__.c_str()),"",25,fUtilV2.Min(),fUtilV2.Max());
    fUtilH2             ->  FillN       (fPoints,fUtilP2,fUtilW0);
    fUtilH2             ->  Write();
    
    TH1F *  fUtilH3     =   new TH1F    (Form("P3_%s",aName__.c_str()),"",25,fUtilV3.Min(),fUtilV3.Max());
    fUtilH3             ->  FillN       (fPoints,fUtilP3,fUtilW0);
    fUtilH3             ->  Write();
    
    fUtilHM             ->  Fit ("gaus","IMREQ0S");
    fUtilHP             ->  Fit ("gaus","IMREQ0S");
    fReturn[0]  =   fUtilHM      ->  GetFunction ("gaus")->GetParameter(2);
    fReturn[1]  =   fUtilHP      ->  GetFunction ("gaus")->GetParameter(2);
    return fReturn;
}

Double_t *          EvalStErInteg   ( TH1D * aTarget,   string aName__ = "",  Bool_t aSaveFit = false )
{
    Int_t       fPoints =   1.e2;
    Double_t *  fReturn = new Double_t [2];
    Double_t *  fUtilPt = new Double_t [fPoints];
    Double_t *  fUtilMn = new Double_t [fPoints];
    Double_t *  fUtilP1 = new Double_t [fPoints];
    Double_t *  fUtilP2 = new Double_t [fPoints];
    Double_t *  fUtilP3 = new Double_t [fPoints];
    Double_t *  fUtilW0 = new Double_t [fPoints];
    
    // Speeding multiple fits
    gROOT->SetBatch();
    
    TH1D *  fRandPt,    *fTotal_;
    for ( Int_t iTer = 0; iTer < fPoints; iTer++ )
    {
        fRandPt =   SetRandPoints(aTarget);
        fTotal_ =   SetSystErrorsh(fRandPt);
        SetLevyTsalis();
        fRandPt->Fit(fLevyFit1D,"IMREQ0S","",0.4,10.);
        if ( bPythiaTest ) fRandPt->Fit(fLevyFit1D,"IMREQ0S","",0.4,1.6);
        fUtilPt[iTer]   =   fLevyFit1D->Moment(1,0.0,0.4);
        fUtilMn[iTer]   =   fLevyFit1D->Integral(0.0,0.4);
        fUtilP1[iTer]   =   fLevyFit1D->GetParameter(1);
        fUtilP2[iTer]   =   fLevyFit1D->GetParameter(2);
        fUtilP3[iTer]   =   fLevyFit1D->GetParameter(3);
        fUtilW0[iTer]   =   1;
        
        if ( aSaveFit )
        {
            TCanvas * c1    =   new TCanvas();
            gPad->SetLogy();
            fTotal_     ->Draw();
            fLevyFit1D  ->Draw("same");
            c1          ->Write();
            delete  c1;
        }
    }

    // Speeding multiple fits
    gROOT->SetBatch(false);

    TVectorD fUtilVP    (fPoints,fUtilPt);
    TVectorD fUtilVM    (fPoints,fUtilMn);
    TVectorD fUtilV1    (fPoints,fUtilP1);
    TVectorD fUtilV2    (fPoints,fUtilP2);
    TVectorD fUtilV3    (fPoints,fUtilP3);
    TVectorD fUtilV0    (fPoints,fUtilW0);
    
    TH1F *  fUtilHP     =   new TH1F    (Form("gp_%s",aName__.c_str()),"",25,fUtilVP.Min(),fUtilVP.Max());
    fUtilHP             ->  FillN       (fPoints,fUtilPt,fUtilW0);
    fUtilHP             ->  Write();
    
    TH1F *  fUtilHM     =   new TH1F    (Form("gm_%s",aName__.c_str()),"",25,fUtilVM.Min(),fUtilVM.Max());
    fUtilHM             ->  FillN       (fPoints,fUtilMn,fUtilW0);
    fUtilHM             ->  Write();
    
    TH1F *  fUtilH1     =   new TH1F    (Form("P1_%s",aName__.c_str()),"",25,fUtilV1.Min(),fUtilV1.Max());
    fUtilH1             ->  FillN       (fPoints,fUtilP1,fUtilW0);
    fUtilH1             ->  Write();
    
    TH1F *  fUtilH2     =   new TH1F    (Form("P2_%s",aName__.c_str()),"",25,fUtilV2.Min(),fUtilV2.Max());
    fUtilH2             ->  FillN       (fPoints,fUtilP2,fUtilW0);
    fUtilH2             ->  Write();
    
    TH1F *  fUtilH3     =   new TH1F    (Form("P3_%s",aName__.c_str()),"",25,fUtilV3.Min(),fUtilV3.Max());
    fUtilH3             ->  FillN       (fPoints,fUtilP3,fUtilW0);
    fUtilH3             ->  Write();
    
    fUtilHM             ->  Fit ("gaus","IMREQ0S");
    fUtilHP             ->  Fit ("gaus","IMREQ0S");
    fReturn[0]  =   fUtilHM      ->  GetFunction ("gaus")->GetParameter(2);
    fReturn[1]  =   fUtilHP      ->  GetFunction ("gaus")->GetParameter(2);
    return fReturn;
}

Double_t *          ExtrapolateVl   ( TH1F * aTarget,   string aName__ = "", Bool_t aSaveFit = false )
{
    Double_t *  fReturn =   new Double_t    [6];    // Result of the Process
                                                    // 1. Mean Value    2. Stat Err     3. Syst Err     4. Mean PT
    // Prepping the FIT
    SetLevyTsalis();
    aTarget                     ->  Fit         (fLevyFit1D,"IMREQ0S","",0.4,10.);
    if ( bPythiaTest ) aTarget  ->  Fit         (fLevyFit1D,"IMREQ0S","",0.4,1.6);
    if (aSaveFit)
    {
        aTarget     ->Write();
        fLevyFit1D  ->Write();
        TCanvas * fSaveToCanvas = new TCanvas(aName__.c_str());
        aTarget     ->Draw();
        if ( bPythiaTest )  fLevyFit1D->SetRange(0.,2.);
        fLevyFit1D  ->Draw("SAME");
        fSaveToCanvas->Write();
        delete fSaveToCanvas;
    }
    Double_t    fHIntWidth      =   HistIntegrals( aTarget,     "W" );
    Double_t    fHIntWidtE      =   HistIntegrals( aTarget,     "WE" );
    Double_t    fHIntWidx1      =   HistIntegrals( aTarget,     "W x^1" );
    Double_t    fHIntWidxE      =   HistIntegrals( aTarget,     "WE x^1" );
    Double_t    fFIntWidth      =   FuncIntegrals( fLevyFit1D,  0.,0.4,"W");
    Double_t    fFIntWidx1      =   FuncIntegrals( fLevyFit1D,  0.,0.4,"W x^1");
    Double_t *  fStatErr        =   EvalStErInteg(aTarget,aName__);
    
    // Mean Value
    fReturn[0]  =   fHIntWidth + fFIntWidth;
    
    //Stat Error
    fReturn[1]  =   fStatErr[0] + fHIntWidtE;
    cout << fStatErr[0] << " " << fHIntWidtE << endl;
    
    //Syst Error
    fReturn[2]  =   fReturn[1] + kSystematicalErrP*fReturn[0];
    
    // Mean PT
    fReturn[3]  =   (fHIntWidx1+fFIntWidx1)/(fHIntWidth+fFIntWidth);
    
    //Stat Error
    fReturn[4]  =   fStatErr[1] + fHIntWidxE;
    
    //Syst Error
    fReturn[5]  =   0;

    return fReturn;
}

Double_t *          ExtrapolateVl   ( TH1D * aTarget,   string aName__ = "", Bool_t aSaveFit = false )
{
    Double_t *  fReturn =   new Double_t    [6];    // Result of the Process
                                                    // 1. Mean Value    2. Stat Err     3. Syst Err     4. Mean PT
    // Prepping the FIT
    SetLevyTsalis();
    aTarget                     ->  Fit         (fLevyFit1D,"IMREQ0S","",0.4,10.);
    if ( bPythiaTest ) aTarget  ->  Fit         (fLevyFit1D,"IMREQ0S","",0.4,1.6);
    if (aSaveFit)
    {
        aTarget     ->Write();
        fLevyFit1D  ->Write();
        TCanvas * fSaveToCanvas = new TCanvas(aName__.c_str());
        aTarget     ->Draw();
        if ( bPythiaTest )  fLevyFit1D->SetRange(0.,2.);
        fLevyFit1D  ->Draw("SAME");
        fSaveToCanvas->Write();
        delete fSaveToCanvas;
    }
    Double_t    fHIntWidth      =   HistIntegrals( aTarget,     "W" );
    Double_t    fHIntWidtE      =   HistIntegrals( aTarget,     "WE" );
    Double_t    fHIntWidx1      =   HistIntegrals( aTarget,     "W x^1" );
    Double_t    fHIntWidxE      =   HistIntegrals( aTarget,     "WE x^1" );
    Double_t    fFIntWidth      =   FuncIntegrals( fLevyFit1D,  0.,0.4,"W");
    Double_t    fFIntWidx1      =   FuncIntegrals( fLevyFit1D,  0.,0.4,"W x^1");
    Double_t *  fStatErr        =   EvalStErInteg(aTarget,aName__);
    
    // Mean Value
    fReturn[0]  =   fHIntWidth + fFIntWidth;
    
    //Stat Error
    fReturn[1]  =   fStatErr[0] + fHIntWidtE;
    cout << fStatErr[0] << " " << fHIntWidtE << endl;
    
    //Syst Error
    fReturn[2]  =   fReturn[1] + kSystematicalErrP*fReturn[0];
    
    // Mean PT
    fReturn[3]  =   (fHIntWidx1+fFIntWidx1)/(fHIntWidth+fFIntWidth);
    
    //Stat Error
    fReturn[4]  =   fStatErr[1] + fHIntWidxE;
    
    //Syst Error
    fReturn[5]  =   0;

    return fReturn;
}

Double_t ***        ExtrapolateVl   ( TH2F * aTarget )
{
    Double_t***     fReturn =   new Double_t**  [nBinPT2D+1];
    TF1  **         fLevyMm =   new TF1 *       [nBinPT2D+1];
    for ( Int_t iTer = 0; iTer < nBinPT2D+1; iTer++ )
    {
        fReturn[iTer]       =   new Double_t*   [2];
        fReturn[iTer][0]    =   new Double_t    [6];
        fReturn[iTer][1]    =   new Double_t    [6];
        fLevyMm[iTer]       =   new TF1         [2];
    }
    
    TH1D *      hSliceFX    =   new TH1D("hSliceFX_","hSliceFX_",nBinPT2D,fArrPT2D);
    TH1D *      hSliceFY    =   new TH1D("hSliceFY_","hSliceFY_",nBinPT2D,fArrPT2D);
    
    TCanvas *   fDrawAllX   =   new TCanvas("cDrawAllX","cDrawAllX");
    fDrawAllX               ->  Divide(5,2);
    TCanvas *   fDrawAllY   =   new TCanvas("cDrawAllY","cDrawAllY");
    fDrawAllY               ->  Divide(5,2);
    TCanvas *   fDrawFull   =   new TCanvas("fDrawFull","fDrawFull");
    fDrawFull               ->  Divide(2,1);
    
    Int_t iHisto = 0;
    for ( Int_t iFit = 1; iFit <= nBinPT2D-2; iFit++ )
    {
        // X-Projection
        fDrawAllX           ->  cd(iFit);
        gStyle              ->  SetOptStat(0);
        gPad                ->  SetLogy();
        hName = Form("XProjection_PT_%.1f_%.1f",fArrPT2D[iFit+1],fArrPT2D[iFit+2]);
        TH1D *  h1D_ResX    =   new TH1D    (*(aTarget->ProjectionX(hName,iFit+2,iFit+2)));
        h1D_ResX            ->  SetTitle(Form("Slice in #phi_{1} P_{T} from %.1f to %.1f GeV",fArrPT2D[iFit+1],fArrPT2D[iFit+2]));
        h1D_ResX            ->  GetXaxis()  ->  SetTitle("#phi_{2} P_{T} GeV");
        h1D_ResX            ->  GetYaxis()  ->  SetTitle("#frac{d^{3}N #phi_{1} }{dydp_{T}d#phi_{2}}(GeV/c)^{-1}");
        h1D_ResX            ->  Fit(fLevyFit1D,"IMREQ0S","",0.4,10.);
        fLevyMm[iFit][0]    =   *fLevyFit1D;
        h1D_ResX            ->  Draw();
        fLevyMm[iFit][0]    .   Draw("same");
        
        fReturn[iFit][0]    =   ExtrapolateVl   (h1D_ResX,Form("X_%d",iFit),true);
        hSliceFY            ->  SetBinContent   (iFit+2,fReturn[iFit][0][0]);
        hSliceFY            ->  SetBinError     (iFit+2,fReturn[iFit][0][1]);
    
        // Y-Projection
        fDrawAllY           ->  cd(iFit);
        gStyle              ->  SetOptStat(0);
        gPad                ->  SetLogy();
        hName = Form("XProjection_PT_%.1f_%.1f",fArrPT2D[iFit+1],fArrPT2D[iFit+2]);
        TH1D *  h1D_ResY    =   new TH1D    (*(aTarget->ProjectionY(hName,iFit+2,iFit+2)));
        h1D_ResY            ->  SetTitle(Form("Slice in #phi_{1} P_{T} from %.1f to %.1f GeV",fArrPT2D[iFit+1],fArrPT2D[iFit+2]));
        h1D_ResY            ->  GetXaxis()  ->  SetTitle("#phi_{2} P_{T} GeV");
        h1D_ResY            ->  GetYaxis()  ->  SetTitle("#frac{d^{3}N #phi_{1} }{dydp_{T}d#phi_{2}}(GeV/c)^{-1}");
        h1D_ResY            ->  Fit(fLevyFit1D,"IMREQ0S","",0.4,10.);
        fLevyMm[iFit][1]    =   *fLevyFit1D;
        h1D_ResY            ->  Draw();
        fLevyMm[iFit][1]    .   Draw("same");
        
        fReturn[iFit][1]    =   ExtrapolateVl   (h1D_ResY,Form("Y_%d",iFit),true);
        hSliceFX            ->  SetBinContent   (iFit+2,fReturn[iFit][1][0]);
        hSliceFX            ->  SetBinError     (iFit+2,fReturn[iFit][1][1]);
    }
    
    fDrawFull           ->  cd(0);
    gStyle              ->  SetOptStat(0);
    gPad                ->  SetLogy();
    hSliceFY            ->  SetTitle("Slice in #phi_{1} P_{T} Y");
    hSliceFY            ->  GetXaxis()  ->  SetTitle("#phi_{2} P_{T} GeV");
    hSliceFY            ->  GetYaxis()  ->  SetTitle("#frac{d^{3}N #phi_{1} }{dydp_{T}d#phi_{2}}(GeV/c)^{-1}");
    hSliceFY            ->  Fit(fLevyFit1D,"IMREQ0S","",0.4,10.);
    fLevyMm[0][0]       =   *fLevyFit1D;
    hSliceFY            ->  Draw();
    fLevyMm[0][0]       .   Draw("same");
    
    fReturn[0][0]       =   ExtrapolateVl   (hSliceFY,"Y_Full",true);
    
    fDrawFull           ->  cd(1);
    gStyle              ->  SetOptStat(0);
    gPad                ->  SetLogy();
    hSliceFX            ->  SetTitle("Slice in #phi_{1} P_{T} X");
    hSliceFX            ->  GetXaxis()  ->  SetTitle("#phi_{2} P_{T} GeV");
    hSliceFX            ->  GetYaxis()  ->  SetTitle("#frac{d^{3}N #phi_{1} }{dydp_{T}d#phi_{2}}(GeV/c)^{-1}");
    hSliceFX            ->  Fit(fLevyFit1D,"IMREQ0S","",0.4,10.);
    fLevyMm[0][1]       =   *fLevyFit1D;
    hSliceFX            ->  Draw();
    fLevyMm[0][1]       .   Draw("same");
    
    fReturn[0][1]       =   ExtrapolateVl   (hSliceFX,"X_Full",true);
    
    fDrawFull->Write();
    fDrawAllX->Write();
    fDrawAllY->Write();
    delete fDrawAllX;
    delete fDrawAllY;
    return fReturn;
}
 
TGraphErrors **     BuildTGraphEr   ( Double_t * aExtrapolateVl1D, Double_t *** aExtrapolateVl2D )
{
    TGraphErrors ** fResult =   new TGraphErrors* [16];
    for ( Int_t iPnt = 0; iPnt < 16; iPnt++ )
    {
        fResult[iPnt]   =   new TGraphErrors();
    }
    
    // Final Results yields 1D
    fResult[0]      ->  SetNameTitle    ("1D_dN_Stat","1D_dN_Stat");
    fResult[0]      ->  SetPoint        (0,1,aExtrapolateVl1D[0]);
    fResult[0]      ->  SetPointError   (0,0,aExtrapolateVl1D[1]);
    fResult[1]      ->  SetNameTitle    ("1D_dN_Syst","1D_dN_Syst");
    fResult[1]      ->  SetPoint        (0,1,aExtrapolateVl1D[0]);
    fResult[1]      ->  SetPointError   (0,0,aExtrapolateVl1D[2]);
    
    // Final Results yields 2D
    fResult[2]      ->  SetNameTitle    ("2D_dN_Stat","2D_dN_Stat");
    fResult[2]      ->  SetPoint        (0,2,aExtrapolateVl2D[0][0][0]);
    fResult[2]      ->  SetPointError   (0,0,aExtrapolateVl2D[0][0][1]);
    fResult[3]      ->  SetNameTitle    ("2D_dN_Syst","2D_dN_Syst");
    fResult[3]      ->  SetPoint        (0,2,aExtrapolateVl2D[0][0][0]);
    fResult[3]      ->  SetPointError   (0,0,aExtrapolateVl2D[0][0][2]);
    /*
    fResult[1]      ->  SetNameTitle    ("2D_dN_Stat","2D_dN_Stat");
    fResult[2]      ->  SetPoint        (0,2,aExtrapolateVl2D[0][0][0]);
    fResult[2]      ->  SetPointError   (0,0,aExtrapolateVl2D[0][1][1]);
    fResult[1]      ->  SetNameTitle    ("2D_dN_Syst","2D_dN_Syst");
    fResult[3]      ->  SetPoint        (0,2,aExtrapolateVl2D[0][0][0]);
    fResult[3]      ->  SetPointError   (0,0,aExtrapolateVl2D[0][1][2]);
    */
     
    // Final Results yields 2D - XY Proj
    for ( Int_t iPnt = 1; iPnt <= 10; iPnt++ )
    {
        fResult[4]   ->  SetNameTitle   ("2D_dN_Stat_X","2D_dN_Stat_X");
        fResult[4]   ->  SetPoint       (iPnt-1,.5*fArrPT2D[iPnt+1]+.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][0][0]);
        fResult[5]   ->  SetPoint       (iPnt-1,.5*fArrPT2D[iPnt+1]+.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][1][0]);
        fResult[5]   ->  SetNameTitle   ("2D_dN_Stat_Y","2D_dN_Stat_Y");
        fResult[4]   ->  SetPointError  (iPnt-1,.5*fArrPT2D[iPnt+1]-.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][0][1]);
        fResult[5]   ->  SetPointError  (iPnt-1,.5*fArrPT2D[iPnt+1]-.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][1][1]);
        fResult[6]   ->  SetNameTitle   ("2D_dN_Syst_X","2D_dN_Syst_X");
        fResult[6]   ->  SetPoint       (iPnt-1,.5*fArrPT2D[iPnt+1]+.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][0][0]);
        fResult[7]   ->  SetPoint       (iPnt-1,.5*fArrPT2D[iPnt+1]+.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][1][0]);
        fResult[7]   ->  SetNameTitle   ("2D_dN_Syst_Y","2D_dN_Syst_Y");
        fResult[6]   ->  SetPointError  (iPnt-1,.5*fArrPT2D[iPnt+1]-.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][0][2]);
        fResult[7]   ->  SetPointError  (iPnt-1,.5*fArrPT2D[iPnt+1]-.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][1][2]);
    }
    
    // Final Results PT 1D
    fResult[8]      ->  SetNameTitle    ("1D_PT_Stat","1D_PT_Stat");
    fResult[8]      ->  SetPoint        (0,1,aExtrapolateVl1D[3]);
    fResult[8]      ->  SetPointError   (0,0,aExtrapolateVl1D[4]);
    fResult[9]      ->  SetNameTitle    ("1D_PT_Syst","1D_PT_Syst");
    fResult[9]      ->  SetPoint        (0,1,aExtrapolateVl1D[3]);
    fResult[9]      ->  SetPointError   (0,0,aExtrapolateVl1D[5]);
    
    // Final Results PT 2D
    fResult[10]     ->  SetNameTitle    ("2D_PT_Stat","2D_PT_Stat");
    fResult[10]     ->  SetPoint        (0,2,aExtrapolateVl2D[0][0][3]);
    fResult[10]     ->  SetPointError   (0,0,aExtrapolateVl2D[0][0][4]);
    fResult[11]     ->  SetNameTitle    ("2D_PT_Syst","2D_PT_Syst");
    fResult[11]     ->  SetPoint        (0,2,aExtrapolateVl2D[0][0][3]);
    fResult[11]     ->  SetPointError   (0,0,aExtrapolateVl2D[0][0][5]);
    
    // Final Results PT 2D - XY Proj
    for ( Int_t iPnt = 1; iPnt <= 10; iPnt++ )
    {
        fResult[12]   ->  SetNameTitle   ("2D_PT_Stat_X","2D_PT_Stat_X");
        fResult[12]   ->  SetPoint       (iPnt-1,  .5*fArrPT2D[iPnt+1]+.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][0][3]);
        fResult[13]   ->  SetPoint       (iPnt-1,  .5*fArrPT2D[iPnt+1]+.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][1][3]);
        fResult[13]   ->  SetNameTitle   ("2D_PT_Stat_Y","2D_PT_Stat_Y");
        fResult[12]   ->  SetPointError  (iPnt-1,  .5*fArrPT2D[iPnt+1]-.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][0][4]);
        fResult[13]   ->  SetPointError  (iPnt-1,  .5*fArrPT2D[iPnt+1]-.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][1][4]);
        fResult[14]   ->  SetNameTitle   ("2D_PT_Syst_X","2D_PT_Syst_X");
        fResult[14]   ->  SetPoint       (iPnt-1,.5*fArrPT2D[iPnt+1]+.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][0][3]);
        fResult[15]   ->  SetPoint       (iPnt-1,.5*fArrPT2D[iPnt+1]+.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][1][3]);
        fResult[15]   ->  SetNameTitle   ("2D_PT_Syst_Y","2D_PT_Syst_Y");
        fResult[14]   ->  SetPointError  (iPnt-1,.5*fArrPT2D[iPnt+1]-.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][0][5]);
        fResult[15]   ->  SetPointError  (iPnt-1,.5*fArrPT2D[iPnt+1]-.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][1][5]);
    }
    
    return fResult;
}

#endif
