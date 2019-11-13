#ifndef SETFUNCTIONS_H
#define SETFUNCTIONS_H
#include "SetValues.h"

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
#include "RooChebychev.h"
#include "RooArgusBG.h"
#include "RooBreitWigner.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooUniform.h"

using namespace std;
using namespace RooFit;
using namespace std::chrono;

enum fitresults1D
{
    ich0,
    ich1,
    ich2,
    ich3,
    ich4,
    inBkg,
    inSig,
    ipMass,
    ipWidth
};
enum fitresults2D
{
    inBB,
    inBS,
    inSB,
    inSS,
};

RooFitResult *  FitModel1D (RooDataHist * data, RooRealVar var)
{
    // Variables
    Int_t nEntries      = data->sumEntries();
    
    // Background
    RooRealVar ch0      = RooRealVar ("ch0","ch0"       ,0.1,-1.,1.);
    RooRealVar ch1      = RooRealVar ("ch1","ch1"       ,0.1,-1.,1.);
    RooRealVar ch2      = RooRealVar ("ch2","ch2"       ,0.1,-1.,1.);
    RooRealVar ch3      = RooRealVar ("ch3","ch3"       ,0.1,-1.,1.);
    RooRealVar ch4      = RooRealVar ("ch4","ch4"       ,0.1,-1.,1.);
    
    //Signal
    RooRealVar pMass    = RooRealVar ("pMass","pMass"   ,1.020,1.010,1.030);
    RooRealVar pWidth   = RooRealVar ("pWidth","pWidth" ,0.005,0.0005,0.01);
    
    // Coefficients
    RooRealVar n0       = RooRealVar ("nSS","nSS"       ,0.5*nEntries,0.,nEntries);
    RooRealVar n1       = RooRealVar ("nBB","nBB"       ,0.5*nEntries,0.,nEntries);
    
    // PDFs
    RooBreitWigner      fSig ("fSig","fSig"             ,var,pMass,pWidth);
    RooChebychev        fBkg ("fBkg","fBkg"             ,var,RooArgSet(ch0,ch1,ch2,ch3,ch4));
    RooAddPdf           fMod ("fMod","fMod"             ,RooArgList(fBkg,fSig),RooArgList(n1,n0));
    RooGenericPdf       fBkg2("fBkg2","max(0,@0*@1-@2*@1-@3)",RooArgSet(var,ch0,ch1,ch2,ch3,ch4));
    RooAddPdf           fMod2("fMod2","fMod2"           ,RooArgList(fBkg2,fSig),RooArgList(n1,n0));
    
    // Fit
    if (BKG2) return    fMod2.fitTo(*data,Extended(kTRUE),SumW2Error(kTRUE),Save());
    return              fMod.fitTo(*data,Extended(kTRUE),SumW2Error(kTRUE),Save());
}

TH1F * HistoModel   (RooFitResult * input, RooRealVar var, char *  hName)
{
    // Background
    RooRealVar ch0      = RooRealVar ("ch0","ch0"       ,static_cast<RooRealVar*>(input->floatParsFinal().at(ich0))->getVal());
    RooRealVar ch1      = RooRealVar ("ch1","ch1"       ,static_cast<RooRealVar*>(input->floatParsFinal().at(ich1))->getVal());
    RooRealVar ch2      = RooRealVar ("ch2","ch2"       ,static_cast<RooRealVar*>(input->floatParsFinal().at(ich2))->getVal());
    RooRealVar ch3      = RooRealVar ("ch3","ch3"       ,static_cast<RooRealVar*>(input->floatParsFinal().at(ich3))->getVal());
    RooRealVar ch4      = RooRealVar ("ch4","ch4"       ,static_cast<RooRealVar*>(input->floatParsFinal().at(ich4))->getVal());
    
    //Signal
    RooRealVar pMass    = RooRealVar ("pMass","pMass"   ,static_cast<RooRealVar*>(input->floatParsFinal().at(ipMass))->getVal());
    RooRealVar pWidth   = RooRealVar ("pWidth","pWidth" ,static_cast<RooRealVar*>(input->floatParsFinal().at(ipWidth))->getVal());
    
    // Coefficients
    RooRealVar n0       = RooRealVar ("nSS","nSS"       ,static_cast<RooRealVar*>(input->floatParsFinal().at(inSig))->getVal());
    RooRealVar n1       = RooRealVar ("nBB","nBB"       ,static_cast<RooRealVar*>(input->floatParsFinal().at(inBkg))->getVal());
    
    // PDFs
    RooBreitWigner      fSig ("fSig","fSig"             ,var,pMass,pWidth);
    RooChebychev        fBkg ("fBkg","fBkg"             ,var,RooArgSet(ch0,ch1,ch2,ch3,ch4));
    RooAddPdf           fMod ("fMod","fMod"             ,RooArgList(fBkg,fSig),RooArgList(n1,n0));
    RooGenericPdf       fBkg2("fBkg2","max(0,@0*@1-@2*@1-@3)",RooArgSet(var,ch0,ch1,ch2,ch3,ch4));
    RooAddPdf           fMod2("fMod2","fMod2"           ,RooArgList(fBkg2,fSig),RooArgList(n1,n0));
    
    if (BKG2) return (TH1F *)   fMod2.createHistogram(hName,var,Binning(nBinIM1D,fMinIM1D,fMaxIM1D));
    return (TH1F *)             fMod.createHistogram(hName,var,Binning(nBinIM1D,fMinIM1D,fMaxIM1D));
}

RooFitResult *  FitModel2D (RooFitResult * utilityx, RooFitResult * utilityy, RooDataHist * data, RooRealVar varx, RooRealVar vary)
{
    // Variables
    Int_t nEntries      = data->sumEntries();
    
    // Background
    RooRealVar ch0x     = RooRealVar ("ch0x","ch0x"     ,static_cast<RooRealVar*>(utilityx->floatParsFinal().at(ich0))->getVal());
    RooRealVar ch1x     = RooRealVar ("ch1x","ch1x"     ,static_cast<RooRealVar*>(utilityx->floatParsFinal().at(ich1))->getVal());
    RooRealVar ch2x     = RooRealVar ("ch2x","ch2x"     ,static_cast<RooRealVar*>(utilityx->floatParsFinal().at(ich2))->getVal());
    RooRealVar ch3x     = RooRealVar ("ch3x","ch3x"     ,static_cast<RooRealVar*>(utilityx->floatParsFinal().at(ich3))->getVal());
    RooRealVar ch4x     = RooRealVar ("ch4x","ch4x"     ,static_cast<RooRealVar*>(utilityx->floatParsFinal().at(ich4))->getVal());
    RooRealVar ch0y     = RooRealVar ("ch0y","ch0y"     ,static_cast<RooRealVar*>(utilityy->floatParsFinal().at(ich0))->getVal());
    RooRealVar ch1y     = RooRealVar ("ch1y","ch1y"     ,static_cast<RooRealVar*>(utilityy->floatParsFinal().at(ich1))->getVal());
    RooRealVar ch2y     = RooRealVar ("ch2y","ch2y"     ,static_cast<RooRealVar*>(utilityy->floatParsFinal().at(ich2))->getVal());
    RooRealVar ch3y     = RooRealVar ("ch3y","ch3y"     ,static_cast<RooRealVar*>(utilityy->floatParsFinal().at(ich3))->getVal());
    RooRealVar ch4y     = RooRealVar ("ch4y","ch4y"     ,static_cast<RooRealVar*>(utilityy->floatParsFinal().at(ich4))->getVal());
    
    //Signal
    RooRealVar pMassx   = RooRealVar ("pMassx","pMassx"  ,static_cast<RooRealVar*>(utilityx->floatParsFinal().at(ipMass))->getVal());
    RooRealVar pWidthx  = RooRealVar ("pWidthx","pWidthx",static_cast<RooRealVar*>(utilityx->floatParsFinal().at(ipWidth))->getVal());
    RooRealVar pMassy   = RooRealVar ("pMassy","pMassy"  ,static_cast<RooRealVar*>(utilityy->floatParsFinal().at(ipMass))->getVal());
    RooRealVar pWidthy  = RooRealVar ("pWidthy","pWidthy",static_cast<RooRealVar*>(utilityy->floatParsFinal().at(ipWidth))->getVal());
    
    // Coefficients
    RooRealVar n0       = RooRealVar ("nSS2D","nSS2D"   ,0.25*nEntries,0.,nEntries);
    RooRealVar n1       = RooRealVar ("nBB2D","nBB2D"   ,0.25*nEntries,0.,nEntries);
    RooRealVar n2       = RooRealVar ("nBS2D","nBS2D"   ,0.25*nEntries,0.,nEntries);
    RooRealVar n3       = RooRealVar ("nSB2D","nSB2D"   ,0.25*nEntries,0.,nEntries);
    
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
    
    // Fit
    return              fMod.fitTo(*data,Extended(kTRUE),SumW2Error(kTRUE),Save());
}

TH2F * HistoModel   (RooFitResult * utilityx, RooFitResult * utilityy, RooFitResult * input, RooRealVar varx, RooRealVar vary, char *  hName)
{
    // Background
    RooRealVar ch0x     = RooRealVar ("ch0x","ch0x"     ,static_cast<RooRealVar*>(utilityx->floatParsFinal().at(ich0))->getVal());
    RooRealVar ch1x     = RooRealVar ("ch1x","ch1x"     ,static_cast<RooRealVar*>(utilityx->floatParsFinal().at(ich1))->getVal());
    RooRealVar ch2x     = RooRealVar ("ch2x","ch2x"     ,static_cast<RooRealVar*>(utilityx->floatParsFinal().at(ich2))->getVal());
    RooRealVar ch3x     = RooRealVar ("ch3x","ch3x"     ,static_cast<RooRealVar*>(utilityx->floatParsFinal().at(ich3))->getVal());
    RooRealVar ch4x     = RooRealVar ("ch4x","ch4x"     ,static_cast<RooRealVar*>(utilityx->floatParsFinal().at(ich4))->getVal());
    RooRealVar ch0y     = RooRealVar ("ch0y","ch0y"     ,static_cast<RooRealVar*>(utilityy->floatParsFinal().at(ich0))->getVal());
    RooRealVar ch1y     = RooRealVar ("ch1y","ch1y"     ,static_cast<RooRealVar*>(utilityy->floatParsFinal().at(ich1))->getVal());
    RooRealVar ch2y     = RooRealVar ("ch2y","ch2y"     ,static_cast<RooRealVar*>(utilityy->floatParsFinal().at(ich2))->getVal());
    RooRealVar ch3y     = RooRealVar ("ch3y","ch3y"     ,static_cast<RooRealVar*>(utilityy->floatParsFinal().at(ich3))->getVal());
    RooRealVar ch4y     = RooRealVar ("ch4y","ch4y"     ,static_cast<RooRealVar*>(utilityy->floatParsFinal().at(ich4))->getVal());
    
    //Signal
    RooRealVar pMassx   = RooRealVar ("pMassx","pMassx"  ,static_cast<RooRealVar*>(utilityx->floatParsFinal().at(ipMass))->getVal());
    RooRealVar pWidthx  = RooRealVar ("pWidthx","pWidthx",static_cast<RooRealVar*>(utilityx->floatParsFinal().at(ipWidth))->getVal());
    RooRealVar pMassy   = RooRealVar ("pMassy","pMassy"  ,static_cast<RooRealVar*>(utilityy->floatParsFinal().at(ipMass))->getVal());
    RooRealVar pWidthy  = RooRealVar ("pWidthy","pWidthy",static_cast<RooRealVar*>(utilityy->floatParsFinal().at(ipWidth))->getVal());
    
    // Coefficients
    RooRealVar n0       = RooRealVar ("nSS2D","nSS2D"   ,static_cast<RooRealVar*>(input->floatParsFinal().at(inSS))->getVal());
    RooRealVar n1       = RooRealVar ("nBB2D","nBB2D"   ,static_cast<RooRealVar*>(input->floatParsFinal().at(inBB))->getVal());
    RooRealVar n2       = RooRealVar ("nBS2D","nBS2D"   ,static_cast<RooRealVar*>(input->floatParsFinal().at(inBS))->getVal());
    RooRealVar n3       = RooRealVar ("nSB2D","nSB2D"   ,static_cast<RooRealVar*>(input->floatParsFinal().at(inSB))->getVal());
    
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
    
    return (TH2F *)     fMod.createHistogram(hName,varx,Binning(nBinIM2D,fMinIM2D,fMaxIM2D),YVar(vary,Binning(nBinIM2D,fMinIM2D,fMaxIM2D)));
}

#endif
