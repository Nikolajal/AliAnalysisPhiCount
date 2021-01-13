// Global Values and constants file
// !TODO: All set!

#ifndef ALIANALYSISUTILITY_H
#define ALIANALYSISUTILITY_H

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
#include "TRandom.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TBenchmark.h"

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
#include "RooClassFactory.h"

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

//------------------------------//
//      GLOBAL VARIABLES        //
//------------------------------//
//
// Random Generator
//
TRandom    *fRandomGen      =   new TRandom();
//
// Benchmark
//
TBenchmark *fBenchmark      =   new TBenchmark();
//
// Title and Name for histograms
//
auto        hName           =   "Name";
auto        hTitle          =   "Title";
//
//---------------------------------------//
//- Physics is by defualt in GeV, we    -//
//- define some constants to evaluate   -//
//- results in other units              -//
//---------------------------------------//
//
auto const  KeV             =   1e6;
auto const  MeV             =   1e3;
auto const  GeV             =   1;
auto const  TeV             =   1e-3;
//
//--------------------------------//
//      BENCHMARK UTILITIES       //
//--------------------------------//
//
TString     fMSG_PrintTimer =   "[INFO] Event # %4.f %s | %02.0f %% | %2.2f %s events/s | Time: %02.0f:%02.0f | ETA: %02.0f:%02.0f \n";
//
//_____________________________________________________________________________
//
void    fStartTimer             ( TString fTimerName )    {
    fBenchmark->Start(fTimerName.Data());
    printf("[INFO] Starting %s \n", fTimerName.Data());
    fflush(stdout);
}
//
//_____________________________________________________________________________
//
void    fStopTimer              ( TString fTimerName )     {
    fBenchmark->Stop(fTimerName.Data());
    printf("[INFO] Stopping %s \n", fTimerName.Data());
    Float_t fElapsedS   = (int)(fBenchmark->GetRealTime(fTimerName.Data()));
    Float_t fElapsedM   = (int)(fElapsedS/60.);
    printf("[INFO] It took %02.0f:%02.0f \n",   fElapsedM,  fElapsedS);
    fflush(stdout);
}
//
//_____________________________________________________________________________
//
void    fPrintLoopTimer         ( TString fTimerName, Int_t iEvent, Int_t nEntries, Int_t iPrintInterval )   {
    if ( iEvent%iPrintInterval != 0 || iEvent == 0 ) return;
    
    // Suffix for events
    TString     fSuffix =   "";
    Int_t       fSfxCor =   iPrintInterval;
    if ( iPrintInterval == 1000 )       {
        fSuffix =   "k";
        fSfxCor =   iPrintInterval%1000;
    }
    if ( iPrintInterval == 1000000 )    {
        fSuffix =   "mln";
        fSfxCor =   iPrintInterval%1000000;
    }
    if ( iPrintInterval == 1000000000 ) {
        fSuffix =   "mld";
        fSfxCor =   iPrintInterval%1000000000;
    }
    
    // Stopping timer
    fBenchmark->Stop(fTimerName.Data());
    
    // Evaluating informations
    Float_t fFraction   =   (float)iEvent/((float)nEntries);
    Float_t fElapsedS   =   (float)(fBenchmark->GetRealTime(fTimerName.Data()));
    Float_t fElapsedM   =   (float)(fElapsedS/60.);
    Float_t fSpeedvsS   =   (float)iEvent/((float)fElapsedS);
    Float_t fSpeedvsM   =   (float)fSpeedvsS*60.;
    Float_t fEta____S   =   (float)nEntries/((float)fSpeedvsS) - (float)fElapsedS;
    Float_t fEta____M   =   (float)(fEta____S/60.);
    Float_t fPrintEvt   =   (float)iEvent*(float)fSfxCor/((float)iPrintInterval);
    
    // Printing
    "[INFO] Event # %4.f %s | %02.0f %% | %2.2f %s events/s | Time: %02.0f:%02.0f | ETA: %02.0f:%02.0f \n";
    printf(fMSG_PrintTimer.Data(),  fPrintEvt,  fSuffix.Data(), 100.*fFraction, fSpeedvsS,  fSuffix.Data(), fElapsedM,  (int)fElapsedS%60,  fEta____M,  (int)fEta____S%60);
    fflush(stdout);
    
    // Resuming timer
    fBenchmark->Start(fTimerName.Data());
}
//
//------------------------------//
//    HISTOGRAM UTILITIES       //
//------------------------------//
//
template < class Tclass >
bool    fIsWorthFitting         ( Tclass * aTarget )    {
    if ( aTarget->GetEntries() <= 2.*aTarget->GetNbinsX()*aTarget->GetNbinsY()*aTarget->GetNbinsZ() )
    {
        cout << "[WARNING] Skipping empty or scarsely populated histogram!" << endl;
        return false;
    }
    return true;
}
//
//_____________________________________________________________________________

#endif
