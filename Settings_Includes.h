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

// RooFitFunction
#include "RooArgusBG.h"
#include "RooBreitWigner.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooUniform.h"

using namespace std;
using namespace RooFit;
using namespace std::chrono;

auto const outMC    = "outGeneratorMC.root";
auto const outPP    = "outPreProcessing.root";
auto const outDF    = "outDataFormat.root";
auto const PTreeName= "PythiaTree";
auto const PEvents  = 1e5;
auto const nBins    = 6e2;
auto const minBound = 0.98;
auto const maxBound = 1.18;


typedef struct
{
    Int_t nKaon, particleID[1024], mother1[1024], mother2[1024], motherID[1024];
    Float_t px[1024], py[1024], pz[1024], pT[1024], e[1024];
} EVKAON;

typedef struct
{
    Int_t   nKaonCouple, iKaon[1024], jKaon[1024];
    Bool_t  bPhi[1024];
    Float_t InvMass[1024], px[1024], py[1024], pz[1024], pT[1024], e[1024];
} EVKAONCOUPLE;

bool checkBool (int a, int b, int c, int d)
{
    if ( a == b ) return false;
    if ( a == c ) return false;
    if ( a == d ) return false;
    if ( b == c ) return false;
    if ( b == d ) return false;
    if ( c == d ) return false;
    return true;
}
