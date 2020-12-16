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

// Random Generator
TRandom    *fRandomGen      =   new TRandom();

// Benchmark
TBenchmark *fBenchmark      =   new TBenchmark();

// Title and Name for histograms
auto        hName           =   "Name";
auto        hTitle          =   "Title";

//---------------------------------------//
//- Physics is by defualt in GeV, we    -//
//- define some constants to evaluate   -//
//- results in other units              -//
//---------------------------------------//

auto const  KeV             =   1e6;
auto const  MeV             =   1e3;
auto const  GeV             =   1;
auto const  TeV             =   1e-3;

//--------------------------------//
//      BENCHMARK UTILITIES       //
//--------------------------------//

void    fStartTimer( string fTimerName )    {
    fBenchmark->Start(fTimerName.c_str());
    cout << "[INFO] Starting " << fTimerName.c_str() << endl;
}

void    fStopTimer( string fTimerName )     {
    fBenchmark->Stop(fTimerName.c_str());
    cout << "[INFO] Stopping " << fTimerName.c_str() << endl;
    Float_t elapsed = fBenchmark->GetRealTime(fTimerName.c_str());
    printf("[INFO] It took %02.0f:%02.0f \n",elapsed / 60., ((int)(elapsed) % 60)*1. );
}

void    fPrintLoopTimer( string fTimerName, Int_t iEvent, Int_t nEntries, Int_t iPrintInterval )   {
    if ( iEvent%iPrintInterval != 0 || iEvent == 0 ) return;
    fBenchmark->Stop(fTimerName.c_str());
    Float_t frac = (Float_t)iEvent / (Float_t)nEntries;
    Float_t elapsed = fBenchmark->GetRealTime(fTimerName.c_str());
    Float_t speed = (Float_t)iEvent / elapsed;
    Float_t eta = (Float_t)nEntries / speed - elapsed;
    printf("[INFO] Event # %4.d mln | %02.0f %% | %7.f events/s | Time: %02.0f:%02.0f | ETA: %02.0f:%02.0f \n", iEvent/iPrintInterval, 100. * frac, speed, elapsed / 60., ((int)(elapsed) % 60)*1. , eta / 60., ((int)(eta) % 60)*1. );
    fflush(stdout);
    fBenchmark->Start(fTimerName.c_str());
}

#endif
