// Global file
//
// !TODO: 1.
//
#ifndef ALIANALYSISPHIPAIR_H
#define ALIANALYSISPHIPAIR_H
//
//>>    Analysis Utility
#include "../AliAnalysisUtility/AliAnalysisUtility.h"
//
//>>    Analysis constants and binning utilities
#include "./AliAnalysisPhiPair_Constants.h"
//
//------------------------------//
//      GLOBAL VARIABLES        //
//------------------------------//
//
//>>    Performance Regulation Values
auto const  kCPU_use                =   1;
auto const  kNCycle_                =   3;
auto        kStatEvalCycles         =   1000;
auto const  kPrintIntervalPP        =   1000000;
auto const  k2DErrorLimit           =   10.;
auto const  kCPUStrategy            =   1;
auto const  kFitOffset              =   true;
auto const  kFitInitHesse           =   true;
auto const  kFitMinos               =   true;
auto const  kFitMinuitStrategy      =   2;
auto const  kSaveToFile             =   true;
//
//_____________________________________________________________________________
//
//>>    Analysis Values
auto const  bPythiaTest             =   kFALSE;
auto const  kOnlyTrue               =   false;
//
//_____________________________________________________________________________
//
//*****                 !TODO: CLEAN
std::vector<std::pair<TF1*,std::vector<float>>> fSystFitFunctions  =   {  {fMTExponential,{1.2,1.4,1.6,2.0}}, {fBoseEinstein,{1.2,1.4,1.6,2.0}}, {fBoltzmann,{1.2,1.4,1.6,2.0}}, {fPowerLaw,{2.0,2.8,4.0}}, {fLevyTsallis,{1.2,1.4,1.6,2.0,2.8,4.0}} };//, {fBGBlastWave,{1.2,1.4,1.6,2.0,2.8,4.0}} };
//
//*****

Float_t
fMeasureMeanPT
 ( TF1 * fLowFit, TGraphAsymmErrors * gTotal, Bool_t fReFit = false );//! TODO: TO BE CLEAN
template < class Tclass >
Float_t                 fEvaluateMeanPT
 ( Tclass* fTotalUncertaintySpectrum, TF1 * fLowPTModel = nullptr, Bool_t kRequestReFit = false )   {
    //
    //  Result
    Float_t fResult = 0;
    Float_t fDenominator = 0;
    //
    //  Re-Fit if requested
    if ( kRequestReFit && fLowPTModel )    fTotalUncertaintySpectrum->Fit(fLowPTModel, "IMEQS");
    //
    //  Storing necessary information
    std::vector<float>  fBinIntegrals;
    std::vector<float>  fBinWidth;
    std::vector<float>  fBinMeanPT;
    //
    //  Extrapolated Contribution
    if ( fLowPTModel )     {
        fBinIntegrals.  push_back( fLowPTModel->Integral( 0., fTotalUncertaintySpectrum->GetBinLowEdge(1) ) );
        fBinWidth.      push_back( 1. );
        fBinMeanPT.     push_back( fLowPTModel->Moment( 1, 0, fTotalUncertaintySpectrum->GetBinLowEdge(1) ) );
    }
    //
    //  Integrated Contribution
    for ( Int_t iBin = 1; iBin <= fTotalUncertaintySpectrum->GetNbinsX(); iBin++ )   {
        fBinIntegrals.  push_back( fTotalUncertaintySpectrum->GetBinContent(iBin) );
        fBinWidth.      push_back( fTotalUncertaintySpectrum->GetBinWidth(iBin) );
        fBinMeanPT.     push_back( fTotalUncertaintySpectrum->GetBinCenter(iBin) );
    }
    //
    //  Calculate Mean PT
    for ( Int_t iBin = 0; iBin < fBinWidth.size(); iBin++ )   {
        fResult +=  fBinIntegrals.at(iBin)*fBinWidth.at(iBin)*fBinMeanPT.at(iBin);
        fDenominator    += fBinIntegrals.at(iBin)*fBinWidth.at(iBin);
    }
    return fResult/fDenominator;
}
template < class Tclass >
void                    fEvaluateMeanPT
 ( Tclass* fTotalUncertaintySpectrum, Double_t &kError )   {
    //
    //  Result
    Float_t fResult = 0;
    Float_t fDenominator = 0;
    //
    //  Storing necessary information
    std::vector<float>  fBinIntegrals;
    std::vector<float>  fBinWidth;
    std::vector<float>  fBinMeanPT;
    //
    //  Integrated Contribution
    for ( Int_t iBin = 1; iBin <= fTotalUncertaintySpectrum->GetNbinsX(); iBin++ )   {
        fBinIntegrals.  push_back( fTotalUncertaintySpectrum->GetBinContent(iBin) + fTotalUncertaintySpectrum->GetBinError(iBin) );
        fBinWidth.      push_back( fTotalUncertaintySpectrum->GetBinWidth(iBin) );
        fBinMeanPT.     push_back( fTotalUncertaintySpectrum->GetBinCenter(iBin) );
    }
    //
    //  Calculate Mean PT
    for ( Int_t iBin = 0; iBin < fBinWidth.size(); iBin++ )   {
        fResult +=  fBinIntegrals.at(iBin)*fBinWidth.at(iBin)*fBinMeanPT.at(iBin);
        fDenominator    += fBinIntegrals.at(iBin)*fBinWidth.at(iBin);
    }
    //
    auto    kUpperLimit =   fResult/fDenominator;
    //
    fResult = 0;
    fDenominator = 0;
    //
    //  Storing necessary information
    fBinIntegrals.clear();
    fBinWidth.clear();
    fBinMeanPT.clear();
    //
    //  Integrated Contribution
    for ( Int_t iBin = 1; iBin <= fTotalUncertaintySpectrum->GetNbinsX(); iBin++ )   {
        fBinIntegrals.  push_back( fTotalUncertaintySpectrum->GetBinContent(iBin) - fTotalUncertaintySpectrum->GetBinError(iBin) );
        fBinWidth.      push_back( fTotalUncertaintySpectrum->GetBinWidth(iBin) );
        fBinMeanPT.     push_back( fTotalUncertaintySpectrum->GetBinCenter(iBin) );
    }
    //
    //  Calculate Mean PT
    for ( Int_t iBin = 0; iBin < fBinWidth.size(); iBin++ )   {
        fResult +=  fBinIntegrals.at(iBin)*fBinWidth.at(iBin)*fBinMeanPT.at(iBin);
        fDenominator    += fBinIntegrals.at(iBin)*fBinWidth.at(iBin);
    }
    //
    auto    kLowerLimit =   fResult/fDenominator;
    //
    kError  =   fabs ( kUpperLimit - kLowerLimit ) / 2. ;
}
//
//-------------------------------------//
//      Analysis Fit Functions         //
//-------------------------------------//
//
//  --  --  General Analysis Functions  --  --  //
//
void
SetBoundaries
 ( TString fOption, Double_t &aValMin, Double_t &aValMax )   {
    aValMin = 0.998;
    aValMax = 1.065;
    //----
    if ( fOption.Contains("RA") )
    {
        aValMin =   0.996;
        aValMax =   1.059;
    }
    if ( fOption.Contains("RB") )
    {
        aValMin =   0.996;
        aValMax =   1.062;
    }
    if ( fOption.Contains("RC") )
    {
        aValMin =   0.996;
        aValMax =   1.065;
    }
    if ( fOption.Contains("RD") )
    {
        aValMin =   0.996;
        aValMax =   1.068;
    }
    if ( fOption.Contains("RE") )
    {
        aValMin =   0.996;
        aValMax =   1.071;
    }
    //----
    if ( fOption.Contains("RF") )
    {
        aValMin =   0.998;
        aValMax =   1.059;
    }
    if ( fOption.Contains("RG") )
    {
        aValMin =   0.998;
        aValMax =   1.062;
    }
    if ( fOption.Contains("RH") )
    {
        aValMin =   0.998;
        aValMax =   1.068;
    }
    if ( fOption.Contains("RI") )
    {
        aValMin =   0.998;
        aValMax =   1.071;
    }
    //----
    if ( fOption.Contains("RJ") )
    {
        aValMin =   1.000;
        aValMax =   1.059;
    }
    if ( fOption.Contains("RK") )
    {
        aValMin =   1.000;
        aValMax =   1.062;
    }
    if ( fOption.Contains("RL") )
    {
        aValMin =   1.000;
        aValMax =   1.065;
    }
    if ( fOption.Contains("RM") )
    {
        aValMin =   1.000;
        aValMax =   1.068;
    }
    if ( fOption.Contains("RN") )
    {
        aValMin =   1.000;
        aValMax =   1.071;
    }
}
//
//  --  --  General Customisation Functions  --  --  //
//
int
fLegendSelect
 ( string fOption )  {
    if ( !fOption.compare("InvMass1D") )    return 1;
    if ( !fOption.compare("RapTru") )       return 1;
    if ( !fOption.compare("Rap") )          return 3;
    if ( !fOption.compare("xInvMass2D") )   return 2;
    if ( !fOption.compare("yInvMass2D") )   return 2;
    else return -1;
}
void
fLegendMaker
 ( RooPlot * fRooPlot, const char * fSelect, TLegend * fLegend )    {
    fLegend                     ->SetFillColorAlpha(kWhite,0.);
    fLegend                     ->SetLineColorAlpha(kWhite,0.);
    fLegend                     ->SetTextSize(.03);
    switch (fLegendSelect(fSelect))
    {
        case 1:
            fLegend                     ->AddEntry(fRooPlot->findObject("RooData"), "stat errors",          "EP");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooSS"),   "Voigtian peak",        "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooBB"),   "#splitline{Chebychev}{background}", "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooMod"),  "Full model",           "L");
            break;
        case 2:
            fLegend                     ->AddEntry(fRooPlot->findObject("RooData"), "stat errors",          "EP");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooSS"),   "Voigt #times Voigt", "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooBS"),   "Cheby #times Voigt", "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooSB"),   "Voigt #times Cheby", "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooBB"),   "Cheby #times Cheby", "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooMod"),  "Full model",          "L");
            break;
        case 3:
            fLegend                     ->AddEntry(fRooPlot->findObject("RooData"), "Data",                 "EP");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooSS"),   "Fit (Sig)",            "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooBB"),   "Fit (Bkg Lin)",        "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooB2"),   "Fit (Bkg Gaus)",       "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooMod"),  "Fit (Model)",          "L");
            break;
        default:
            cout << "Improper option, no changes made" << endl;
            break;
    }
}
int
fAxisSelect
 ( string fOption ){
    if ( !fOption.compare("InvMass1D") )    return 1;
    if ( !fOption.compare("xInvMass2D") )   return 2;
    if ( !fOption.compare("yInvMass2D") )   return 3;
    if ( !fOption.compare("Rap") )          return 4;
    if ( !fOption.compare("RapTru") )       return 4;
    else return -1;
}
void
fAxisMaker
 ( RooPlot * fRooPlot, const char * fSelect ){
    switch (fAxisSelect(fSelect))
    {
        case 1:
            fRooPlot                    ->GetXaxis()->SetTitle("M_{KK} (GeV/#it{c}^{2})");
            break;
        case 2:
            fRooPlot                    ->GetXaxis()->SetTitle("M_{KK,#phi_{1}} (GeV/#it{c}^{2})");
            break;
        case 3:
            fRooPlot                    ->GetXaxis()->SetTitle("M_{KK,#phi_{2}} (GeV/#it{c}^{2})");
            break;
        case 4:
            fRooPlot                    ->GetXaxis()->SetTitle("|#Delta y| #phi_{1,2}");
            break;
        default:
            cout << "Improper option, no changes made" << endl;
            break;
    }
}
int
fPlotterSelect
 ( string fOption ){
    if ( !fOption.compare("InvMass1D") )    return 1;
    if ( !fOption.compare("RapTru") )       return 1;
    if ( !fOption.compare("Rap") )          return 3;
    if ( !fOption.compare("xInvMass2D") )   return 2;
    if ( !fOption.compare("yInvMass2D") )   return 2;
    else return -1;
}
void
fRooPlotPlotter
 ( RooPlot * fRooPlot, const char * fSelect, RooAddPdf fModel , RooDataHist * fData ){
    switch (fPlotterSelect(fSelect))
    {
        case 1:
            fData                           ->plotOn(fRooPlot,      MarkerColor(colors[0]),                MarkerStyle(markers[2]),    Name("RooData"));
            fModel                          .plotOn (fRooPlot,      LineColor(colors[2]),                   LineStyle(kSolid), Name("RooMod"));
            fModel                          .plotOn (fRooPlot,      Components("fBkg"),             LineStyle(kDashed), LineColor(colors[2]),      Name("RooBB"));
            fModel                          .plotOn (fRooPlot,      Components("fSig"),             LineColor(colors[1]),       Name("RooSS"));
            fData                           ->plotOn(fRooPlot,      MarkerColor(colors[0]),                MarkerStyle(markers[2]),    Name("RooData"));
            break;
        case 2:
            fData                           ->plotOn(fRooPlot,      CutRange("fDrawRange"),         MarkerColor(colors[0]),    MarkerStyle(markers[2]) ,   Name("RooData"));
            fModel                          .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  LineColor(colors[2]),       LineStyle(kSolid), Name("RooMod"));
            fModel                          .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  Components("fBkg"), LineStyle(kDashed), LineColor(colors[4]),      Name("RooBB"));
            fModel                      .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  Components("fSigSig"),LineColor(colors[1]),       Name("RooSS"));
            fModel                      .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  Components("fSigBkg"),LineStyle(kDashed), LineColor(colors[5]),    Name("RooSB"));
            fModel                      .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  Components("fBkgSig"),LineStyle(kDashed), LineColor(colors[6]),    Name("RooBS"));
            fData                           ->plotOn(fRooPlot,      CutRange("fDrawRange"),         MarkerColor(colors[0]),    MarkerStyle(markers[2]) ,   Name("RooData"));
            break;
        case 3:
            fData                           ->plotOn(fRooPlot,      MarkerColor(colors[0]),                MarkerStyle(markers[2]),    Name("RooData"));
            fModel                          .plotOn (fRooPlot,      LineColor(colors[1]),                   LineStyle(kSolid), Name("RooMod"));
            fModel                          .plotOn (fRooPlot,      Components("fBkg"),             LineStyle(kDashed), LineColor(colors[2]),      Name("RooBB"));
            fModel                          .plotOn (fRooPlot,      Components("fBk2"),             LineStyle(kDashed), LineColor(33),      Name("RooB2"));
            fModel                          .plotOn (fRooPlot,      Components("fSig"),             LineColor(colors[2]),       Name("RooSS"));
            fData                           ->plotOn(fRooPlot,      MarkerColor(colors[0]),                MarkerStyle(markers[2]),    Name("RooData"));
            break;
        default:
            cout << "Improper option, no changes made" << endl;
            break;
    }
}
void
fRooPlotMaker
 ( RooPlot * fRooPlot, TLegend * fLegend, RooAddPdf fModel , RooDataHist * fData, const char * fSelect ){
    fRooPlot->SetTitle("");
    fRooPlotPlotter(fRooPlot,fSelect,fModel,fData);
    fLegendMaker(fRooPlot,fSelect,fLegend);
    fAxisMaker(fRooPlot,fSelect);
}
//
//  --  --  1D Analysis Functions  --  --  //
//
RooFitResult*
fFitCoreModel
 ( TH1 * THdata, TH1F* hSlopReference, TString fName = "", TString fOption = "", Int_t PTindex = -1, Int_t PTDimension = 1, TString fPathToSave = "./result/SEFitCheck" )    {
    //
    //>>    Silencing TCanvas Pop-Up and speeding Fit procedures
    gROOT->SetBatch(kTRUE);
    //
    //>>    Selecting all Options
    Bool_t  fDegree2,   fDegree4,   fLosWidt,   fUsePoly,   fAlwaysFit, fSaveToFile,    fUseHighRes,    fUseLowRes,    fUseFreeRes;
    //
    fDegree2    =   fOption.Contains("DG2", TString::kIgnoreCase);
    fDegree4    =   fOption.Contains("DG4", TString::kIgnoreCase);
    fLosWidt    =   fOption.Contains("WDT", TString::kIgnoreCase);
    fUsePoly    =   fOption.Contains("POL", TString::kIgnoreCase);
    fUseHighRes =   fOption.Contains("RSH", TString::kIgnoreCase);
    fUseLowRes  =   fOption.Contains("RSL", TString::kIgnoreCase);
    fUseFreeRes =   fOption.Contains("RSX", TString::kIgnoreCase);
    fAlwaysFit  =   fOption.Contains("FIT", TString::kIgnoreCase);
    fSaveToFile =   fOption.Contains("SAVE",TString::kIgnoreCase);
    //
    //>>    Check there is a reasonable amount of entries form the general method, or if it has been overrided
    if ( !fIsWorthFitting( THdata ) || fAlwaysFit ) return nullptr;
    //
    //>>    Global Variables
    Double_t                        fInvMassValMax, fInvMassValMin;
    SetBoundaries(  fOption,        fInvMassValMin, fInvMassValMax  );
    Int_t           nEntries    =   THdata->GetEntries();
    RooRealVar      InvMass     =   RooRealVar        ("InvMass",   "InvMass",      fInvMassValMin, fInvMassValMax                  );
    RooDataHist*    data        =   new RooDataHist   ("Data",      "Data",         InvMass,        Import(*THdata)                 );
    RooDataHist*    dataLoose   =   new RooDataHist   ("DataLoose", "DataLoose",    InvMass,        Import(*fLooseErrors(THdata))   );
    Int_t           kNCycle     =   kNCycle_;
    Float_t         fRescaleRes =   1.;
    if ( fUseHighRes )  fRescaleRes = 1.15;
    if ( fUseLowRes  )  fRescaleRes = 0.95;
    //
    //>>    Background PDF Coefficients
    RooRealVar ch1, ch2, ch3, ch4;
                                ch1     =   RooRealVar      ("ch1",     "ch1",      0.5,    -1,         1   );
                                ch2     =   RooRealVar      ("ch2",     "ch2",      -.06,   -1,         1   );
    if ( fDegree2 && !fDegree4 )ch3     =   RooRealVar      ("ch3",     "ch3",      0.  );
    else                        ch3     =   RooRealVar      ("ch3",     "ch3",      0.01,   -1,         1   );
    if ( fDegree4 )             ch4     =   RooRealVar      ("ch4",     "ch4",      0.,     -1,         1   );
    else                        ch4     =   RooRealVar      ("ch4",     "ch4",      0.  );
    //
    //>>    Signal PDF Coefficients
    RooRealVar sMass, sWidt, sSlop;
                                sMass   =   RooRealVar      ("bMass",   "bMass",    kPhiMesonMass_,  kPhiMesonMass_*0.5,  kPhiMesonMass_*1.5);
    if ( fLosWidt )             sWidt   =   RooRealVar      ("bWidt",   "bWidt",    kPhiMesonWidth,  kPhiMesonWidth*0.5,  kPhiMesonWidth*1.5);
    else                        sWidt   =   RooRealVar      ("bWidt",   "bWidt",    kPhiMesonWidth);
    if ( bPythiaTest )          sSlop   =   RooRealVar      ("bSlop",   "bSlop",    0.);
    else if ( !fUseFreeRes )    sSlop   =   RooRealVar      ("bSlop",   "bSlop",    fRescaleRes*hSlopReference->GetBinContent(PTindex+1)/1000.);
    else                        sSlop   =   RooRealVar      ("bSlop",   "bSlop",    fRescaleRes*hSlopReference->GetBinContent(PTindex+1)/1000.);
    //
    //>>    Normalisation Coefficients
    RooRealVar  nSS,    nBB;
                                nSS     =   RooRealVar      ("anSS",    "anSS",     .5*nEntries,    0., nEntries);
                                nBB     =   RooRealVar      ("anBB",    "anBB",     .5*nEntries,    0., nEntries);
    //
    //>>    Building the PDFs
    RooVoigtian                 fSig    =   RooVoigtian     ("fSig",    "fSig",     InvMass,    sMass,  sWidt,  sSlop);
    RooChebychev                fBkg    =   RooChebychev    ("fBkg",    "fBkg",     InvMass,    RooArgSet(ch1,ch2,ch3,ch4));
    RooPolynomial               fBkgPol =   RooPolynomial   ("fBkg",    "fBkg",     InvMass,    RooArgSet(ch1,ch2,ch3,ch4));
    RooAddPdf                  *fMod;
    if ( fUsePoly )   {
        fMod                            =   new RooAddPdf   ("fMod",    "fMod",     RooArgList(fBkgPol,fSig),   RooArgList(nBB,nSS));
        ch1.removeRange();
        ch2.removeRange();
        ch3.removeRange();
        ch4.removeRange();
    }   else    {
        fMod                            =   new RooAddPdf   ("fMod",    "fMod",     RooArgList(fBkg,fSig),      RooArgList(nBB,nSS));
    }
    //
    RooFitResult* fFitResults;
    fFitResults      =   fMod->fitTo(*dataLoose,Save(),NumCPU(kCPU_use,kCPUStrategy));
    for ( Int_t iCycle = 0; iCycle < kNCycle; iCycle++ )    {
        if ( !kOnlyTrue )   {
            fFitResults      =   fMod->fitTo(*data,Extended(kTRUE),Save(),NumCPU(kCPU_use,kCPUStrategy),Offset(kFitOffset),Strategy(kFitMinuitStrategy),InitialHesse(kFitInitHesse),Minos(kFitMinos));
            auto N_Raw  =   static_cast<RooRealVar*>(fFitResults ->floatParsFinal().find("anSS"));
            if ( fIsResultAcceptable(N_Raw->getVal(),N_Raw->getError()) ) break;
        }
        if ( kOnlyTrue )    {
            fFitResults      =   fSig.fitTo(*data,Extended(kTRUE),Save(),NumCPU(kCPU_use,kCPUStrategy),Offset(kFitOffset),Strategy(kFitMinuitStrategy),InitialHesse(kFitInitHesse),Minos(kFitMinos));
        }
    }
    //
    // Modify Raw with missing signal
    
    //
    if ( fSaveToFile || kSaveToFile )
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
        
        SetStyle();
        
        TCanvas * fSaveToCanvas;
        if ( PTDimension == 1 )fSaveToCanvas   =   new TCanvas(
                                                Form("PT_%.1f_%.1f_1D_%s",fArrPT1D[PTindex],fArrPT1D[PTindex+1],fName.Data()),
                                                Form("PT_%.1f_%.1f_1D_%s",fArrPT1D[PTindex],fArrPT1D[PTindex+1],fName.Data())
                                                );
        
        if ( PTDimension == 2 )fSaveToCanvas   =   new TCanvas(
                                                Form("PT_%.1f_%.1f_2D_%s",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fName.Data()),
                                                Form("PT_%.1f_%.1f_2D_%s",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fName.Data())
                                                );
        
        RooPlot * fSaveToFrame      = InvMass.frame(Name(hName),Title(hTitle));
        TLegend * fLegend           = new TLegend   (0.18,0.85,0.33,0.60);
        
        fRooPlotMaker(fSaveToFrame,fLegend,*fMod,data,"InvMass1D");
        
        if ( PTDimension == 1 )fSaveToFrame                ->GetYaxis()->SetTitle(Form("Counts/( %.1f MeV/#it{c}^{2} )",1000*kBinningPrecision1D));
        if ( PTDimension == 2 )fSaveToFrame                ->GetYaxis()->SetTitle(Form("Counts/( %.1f MeV/#it{c}^{2} )",1000*kBinningPrecision2D));
        
        auto fMaximum   =   fSaveToFrame->GetMaximum();
        fSaveToFrame    ->  SetMaximum(fMaximum*1.30);
        
        fSaveToFrame                ->Draw("same");
        fLegend                     ->Draw("same");
        
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.59, 0.83,"ALICE Performance");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        uLatex->DrawLatexNDC(0.59, 0.77,"pp #sqrt{#it{s}}= 7 TeV");
        if ( PTDimension == 1 ) uLatex->DrawLatexNDC(0.59, 0.72,Form("%.2f < #it{p}_{T} < %.2f GeV/#it{c}",fArrPT1D[PTindex],fArrPT1D[PTindex+1]));
        if ( PTDimension == 2 ) uLatex->DrawLatexNDC(0.59, 0.72,Form("%.2f < #it{p}_{T} < %.2f GeV/#it{c}",fArrPT2D[PTindex],fArrPT2D[PTindex+1]));
        uLatex->DrawLatexNDC(0.59, 0.67,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
        
        fSaveToCanvas               ->Write ();
        if ( PTDimension == 1 )fSaveToCanvas               ->SaveAs(Form("%s/PT_%.1f_%.1f_1D_%s.pdf",fPathToSave.Data(),fArrPT1D[PTindex],fArrPT1D[PTindex+1],fName.Data()));
        if ( PTDimension == 2 )fSaveToCanvas               ->SaveAs(Form("%s/PT_%.1f_%.1f_1D_%s.pdf",fPathToSave.Data(),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fName.Data()));
        if ( PTDimension == 1 )fSaveToCanvas               ->SaveAs(Form("%s/PT_%.1f_%.1f_1D_%s.eps",fPathToSave.Data(),fArrPT1D[PTindex],fArrPT1D[PTindex+1],fName.Data()));
        if ( PTDimension == 2 )fSaveToCanvas               ->SaveAs(Form("%s/PT_%.1f_%.1f_1D_%s.eps",fPathToSave.Data(),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fName.Data()));
        delete fSaveToCanvas;
    }
    
    // Un-Silencing TCanvas Pop-Up
    gROOT->SetBatch(kFALSE);
    
    return fFitResults;
}
//
std::vector<TH1F*>
FitModel
 ( TH1F **hTarget, TH1F* hSlopReference, RooFitResult** &fFitresultsStore, Int_t fDimension, TString fOption = "", TString fTargetPath = "./result/SEFitCheck", TString fNameFile = "" )  {
    std::vector<TH1F*>  fResults;
    if ( !fFitresultsStore )  {  fFitresultsStore = new RooFitResult*[((fDimension == 1)? nBinPT1D : nBinPT2D)]; }
    for ( Int_t iFit = 0; iFit < ((fDimension == 1)? nBinPT1D : nBinPT2D); iFit++ )
    {
        //>>    Fit to model
        fFitresultsStore[iFit]      =   fFitCoreModel(hTarget[iFit],hSlopReference,fNameFile,fOption,iFit,fDimension,fTargetPath);
        //
        //>>    SegFault Protection
        if ( !fFitresultsStore[iFit] ) continue;
        //
        //>>    Building Raw Count histograms
        if ( iFit == 0 )   {
            for ( auto fCoeff : fFitresultsStore[iFit]->floatParsFinal() )   {
                auto N_Raw      = static_cast<RooRealVar*>(fCoeff);
                if  ( fDimension == 1)  fResults.push_back( new TH1F(Form("%s%s_%s","h1D_",    N_Raw->GetName(), fNameFile.Data()),  Form("%s%s","h1D_",     N_Raw->GetName()),nBinPT1D,fArrPT1D) );
                else                    fResults.push_back( new TH1F(Form("%s%s_%s","h2Dbin_", N_Raw->GetName(), fNameFile.Data()),  Form("%s%s","h2Dbin_",  N_Raw->GetName()),nBinPT2D,fArrPT2D) );
            }
            for ( auto fCoeff : fFitresultsStore[iFit]->constPars() )   {
                auto N_Raw      = static_cast<RooRealVar*>(fCoeff);
                if  ( fDimension == 1)  fResults.push_back( new TH1F(Form("%s%s_%s","h1D_",    N_Raw->GetName(), fNameFile.Data()),  Form("%s%s","h1D_",     N_Raw->GetName()),nBinPT1D,fArrPT1D) );
                else                    fResults.push_back( new TH1F(Form("%s%s_%s","h2Dbin_", N_Raw->GetName(), fNameFile.Data()),  Form("%s%s","h2Dbin_",  N_Raw->GetName()),nBinPT2D,fArrPT2D) );
            }
            if  ( fDimension == 1)  {
                fResults.push_back( new TH1F("hRAW_1D","hRAW_1D",nBinPT1D,fArrPT1D) );
                fResults.push_back( new TH1F("hCOR_1D","hCOR_1D",nBinPT1D,fArrPT1D) );
            }   else    {
                fResults.push_back( new TH1F("hRAW_1D_in_2D_bin","hRAW_1D_in_2D_bin",nBinPT1D,fArrPT1D) );
                fResults.push_back( new TH1F("hCOR_1D_in_2D_bin","hCOR_1D_in_2D_bin",nBinPT1D,fArrPT1D) );
            }
        }
        //
        //>>    Filling Raw Count Histograms
        Int_t   iTer = 0;
        Float_t fCount, fCounE, fMean, fWidth;
        Double_t    fMin,   fMax;
        for ( auto fCoeff : fFitresultsStore[iFit]->floatParsFinal() )   {
            auto N_Raw      = static_cast<RooRealVar*>(fCoeff);
            fResults.at(iTer)->SetBinContent          (iFit+1,N_Raw->getVal());
            fResults.at(iTer)->SetBinError            (iFit+1,N_Raw->getError());
            if ( strncmp(N_Raw->GetName(),"anSS",4) == 0 )  {
                fCount = N_Raw->getVal();
                fCounE = N_Raw->getError();
            }
            if ( strncmp(N_Raw->GetName(),"bMass",5) == 0 ) fMean   =   N_Raw->getVal();
            if ( strncmp(N_Raw->GetName(),"bWidt",5) == 0 ) fWidth  =   N_Raw->getVal();
            iTer++;
        }
        for ( auto fCoeff : fFitresultsStore[iFit]->constPars() )   {
            auto N_Raw      = static_cast<RooRealVar*>(fCoeff);
            fResults.at(iTer)->SetBinContent          (iFit+1,N_Raw->getVal());
            fResults.at(iTer)->SetBinError            (iFit+1,N_Raw->getError());
            if ( strncmp(N_Raw->GetName(),"bMass",5) == 0 ) fMean   =   N_Raw->getVal();
            if ( strncmp(N_Raw->GetName(),"bWidt",5) == 0 ) fWidth  =   N_Raw->getVal();
            iTer++;
        }
        //
        RooRealVar      vIntU   =   RooRealVar      ("vIntU",   "vIntU",    2*kKaonMass,  1000 );
        RooRealVar      vMean   =   RooRealVar      ("vMean",   "vMean",    fMean   );
        RooRealVar      vStdv   =   RooRealVar      ("vStdv",   "vStdv",    fWidth  );
        //
        SetBoundaries(fOption,fMin,fMax);
        //
        vIntU.setRange("Full",2*kKaonMass,1000);
        vIntU.setRange("Meas",fMin,fMax);
        //
        RooBreitWigner  hUtil   =   RooBreitWigner  ("hUtil",   "hUtil",    vIntU,  vMean,  vStdv);
        //
        Float_t         kCorr   =   hUtil.analyticalIntegral(1,"Meas")/hUtil.analyticalIntegral(1,"Full");
        //
        fResults.at(iTer)->SetBinContent          (iFit+1,fCount/kCorr);
        fResults.at(iTer)->SetBinError            (iFit+1,fCounE/kCorr);
        iTer++;
        //
        fResults.at(iTer)->SetBinContent          (iFit+1,kCorr);
        fResults.at(iTer)->SetBinError            (iFit+1,0);
        iTer++;
    }
    //
    return  fResults;
}
std::vector<TH1F*>
FitModel
 ( TH1F **hTarget, TH1F* hSlopReference, TString fTargetPath = "./result/SEFitCheck", TString fNameFile = "", TString fOption = "" )  {
    return FitModel(hTarget,hSlopReference,NULL_ROOFITPTR2,1,fOption,fTargetPath,fNameFile);
}
std::vector<TH1F*>
FitModel
 ( TH1F **hTarget, TH1F* hSlopReference, RooFitResult** &fFitresultsStore, TString fTargetPath = "./result/SEFitCheck", TString fNameFile = "", TString fOption = "" )  {
    return FitModel(hTarget,hSlopReference,fFitresultsStore,2,fOption,fTargetPath,fNameFile);
}
//
//  --  --  2D Analysis Functions  --  --  //
//
RooFitResult*
fFitCoreModel
 ( TH2F * THdata, TH1F* hSlopReference, RooFitResult * fFitShapeX, RooFitResult * fFitShapeY, TString fName = "", TString fOption = "", Int_t PTindex = -1, Int_t PTjndex = -1, TString fPathToSave = "./result/SEFitCheck" )    {
    //
    //>>    Silencing TCanvas Pop-Up and speeding Fit procedures
    gROOT->SetBatch(kTRUE);
    //
    //>>    Selecting all Options
    Bool_t  fFreeBkg,   fUsePoly,   fAlwaysFit, fSaveToFile;
    //
    fFreeBkg    =   fOption.Contains("BKG", TString::kIgnoreCase);
    fUsePoly    =   fOption.Contains("POL", TString::kIgnoreCase);
    fAlwaysFit  =   fOption.Contains("FIT", TString::kIgnoreCase);
    fSaveToFile =   fOption.Contains("SAVE",TString::kIgnoreCase);
    //
    //>>    Check there is a reasonable amount of entries form the general method, or if it has been overrided
    if ( !fIsWorthFitting( THdata ) || fAlwaysFit ) return nullptr;
    //
    //>>    Global Variables
    Double_t                        fInvMassValMax, fInvMassValMin;
    SetBoundaries(  fOption,        fInvMassValMin, fInvMassValMax  );
    Int_t           nEntries    =   THdata->GetEntries();
    RooRealVar      xInvMass    =   RooRealVar        ("xInvMass",  "xInvMass",     fInvMassValMin, fInvMassValMax                                  );
    RooRealVar      yInvMass    =   RooRealVar        ("yInvMass",  "yInvMass",     fInvMassValMin, fInvMassValMax                                  );
    RooDataHist*    data        =   new RooDataHist   ("Data",      "Data",         RooArgList(xInvMass,yInvMass),  Import(*THdata)                 );
    RooDataHist*    dataLoose   =   new RooDataHist   ("DataLoose", "DataLoose",    RooArgList(xInvMass,yInvMass),  Import(*fLooseErrors(THdata))   );
    Int_t kNCycle               =   kNCycle_;
    //
    //>>    Recovering the previous fits
    RooArgSet      *fXShapes    =   new RooArgSet(fFitShapeX->floatParsFinal(),fFitShapeX->constPars());
    RooArgSet      *fYShapes    =   new RooArgSet(fFitShapeY->floatParsFinal(),fFitShapeY->constPars());
    //
    //>>    Background
    RooRealVar  ch1x, ch2x, ch3x, ch4x, ch1y, ch2y, ch3y, ch4y;
    //
    if ( fFreeBkg )
    {
        ch1x     = RooRealVar ("ch1x","ch1x"     ,fXShapes->getRealValue("ch1",0),fXShapes->getRealValue("ch1",0)-0.1,fXShapes->getRealValue("ch1",0)+0.1);
        ch2x     = RooRealVar ("ch2x","ch2x"     ,fXShapes->getRealValue("ch2",0),fXShapes->getRealValue("ch2",0)-0.1,fXShapes->getRealValue("ch2",0)+0.1);
        ch3x     = RooRealVar ("ch3x","ch3x"     ,fXShapes->getRealValue("ch3",0),fXShapes->getRealValue("ch3",0)-0.1,fXShapes->getRealValue("ch3",0)+0.1);
        ch4x     = RooRealVar ("ch4x","ch4x"     ,fXShapes->getRealValue("ch4",0),fXShapes->getRealValue("ch4",0)-0.1,fXShapes->getRealValue("ch4",0)+0.1);
        ch1y     = RooRealVar ("ch1y","ch1y"     ,fYShapes->getRealValue("ch1",0),fYShapes->getRealValue("ch1",0)-0.1,fYShapes->getRealValue("ch1",0)+0.1);
        ch2y     = RooRealVar ("ch2y","ch2y"     ,fYShapes->getRealValue("ch2",0),fYShapes->getRealValue("ch2",0)-0.1,fYShapes->getRealValue("ch2",0)+0.1);
        ch3y     = RooRealVar ("ch3y","ch3y"     ,fYShapes->getRealValue("ch3",0),fYShapes->getRealValue("ch3",0)-0.1,fYShapes->getRealValue("ch3",0)+0.1);
        ch4y     = RooRealVar ("ch4y","ch4y"     ,fYShapes->getRealValue("ch4",0),fYShapes->getRealValue("ch4",0)-0.1,fYShapes->getRealValue("ch4",0)+0.1);
    }
    else
    {
        ch1x     = RooRealVar ("ch1x","ch1x"     ,fXShapes->getRealValue("ch1",0));
        ch2x     = RooRealVar ("ch2x","ch2x"     ,fXShapes->getRealValue("ch2",0));
        ch3x     = RooRealVar ("ch3x","ch3x"     ,fXShapes->getRealValue("ch3",0));
        ch4x     = RooRealVar ("ch4x","ch4x"     ,fXShapes->getRealValue("ch4",0));
        ch1y     = RooRealVar ("ch1y","ch1y"     ,fYShapes->getRealValue("ch1",0));
        ch2y     = RooRealVar ("ch2y","ch2y"     ,fYShapes->getRealValue("ch2",0));
        ch3y     = RooRealVar ("ch3y","ch3y"     ,fYShapes->getRealValue("ch3",0));
        ch4y     = RooRealVar ("ch4y","ch4y"     ,fYShapes->getRealValue("ch4",0));
    }
    
    //Signal
    RooRealVar pMassx   = RooRealVar ("pMassx","pMassx" ,fXShapes->getRealValue("bMass",0));
    RooRealVar pWidthx  = RooRealVar ("pWidtx","pWidtx" ,fXShapes->getRealValue("bWidt",0));
    RooRealVar pSlopex  = RooRealVar ("pSlopx","pSlopx" ,fXShapes->getRealValue("bSlop",0));
    RooRealVar pMassy   = RooRealVar ("pMassy","pMassy" ,fYShapes->getRealValue("bMass",0));
    RooRealVar pWidthy  = RooRealVar ("pWidty","pWidty" ,fYShapes->getRealValue("bWidt",0));
    RooRealVar pSlopey  = RooRealVar ("pSlopy","pSlopy" ,fYShapes->getRealValue("bSlop",0));
    
    // Coefficients
    auto fS__x  =   fXShapes->getRealValue("anSS",0);
    auto fB__x  =   fXShapes->getRealValue("anBB",0);
    auto fS__y  =   fYShapes->getRealValue("anSS",0);
    auto fB__y  =   fYShapes->getRealValue("anBB",0);
    auto fTot_  =   fS__x*fS__y+fB__x*fB__y+fS__x*fB__y+fB__x*fS__y;
    
    RooRealVar n0       = RooRealVar ("anSS2D","anSS2D" ,nEntries*(fS__x*fS__y)/fTot_,0.,nEntries);
    RooRealVar n1       = RooRealVar ("anBB2D","anBB2D" ,nEntries*(fB__x*fB__y)/fTot_,0.,nEntries);
    RooRealVar n2       = RooRealVar ("anBS2D","anBS2D" ,nEntries*(fB__x*fS__y)/fTot_,0.,nEntries);
    RooRealVar n3       = RooRealVar ("anSB2D","anSB2D" ,nEntries*(fS__x*fB__y)/fTot_,0.,nEntries);
    
    // PDFs
    RooChebychev        fBkgx   ("fBkgx","fBkgx"        ,xInvMass,RooArgSet(ch1x,ch2x,ch3x,ch4x));
    RooPolynomial       fBk2x   ("fBkgx","fBkgx"        ,xInvMass,RooArgSet(ch1x,ch2x,ch3x,ch4x));
    RooVoigtian         fSigx   ("fSigx","fSigx"        ,xInvMass,pMassx,pWidthx,pSlopex);
    RooChebychev        fBkgy   ("fBkgy","fBkgy"        ,yInvMass,RooArgSet(ch1y,ch2y,ch3y,ch4y));
    RooPolynomial       fBk2y   ("fBkgy","fBkgy"        ,yInvMass,RooArgSet(ch1y,ch2y,ch3y,ch4y));
    RooVoigtian         fSigy   ("fSigy","fSigy"        ,yInvMass,pMassy,pWidthy,pSlopey);
    RooProdPdf         *fBB;
    if ( fUsePoly )   {
        fBB= new RooProdPdf         ("fBkg","fBkg"          ,fBk2x,fBk2y);
    }   else    {
        fBB= new RooProdPdf         ("fBkg","fBkg"          ,fBkgx,fBkgy);
    }
    RooProdPdf         *fSB;
    if ( fUsePoly )   {
        fSB= new RooProdPdf         ("fSigBkg","fSigBkg"    ,fSigx,fBk2y);
    }   else    {
        fSB= new RooProdPdf         ("fSigBkg","fSigBkg"    ,fSigx,fBkgy);
    }
    RooProdPdf         *fBS;
    if ( fUsePoly )   {
        fBS= new RooProdPdf         ("fBkgSig","fBkgSig"          ,fBk2x,fSigy);
    }   else    {
        fBS= new RooProdPdf         ("fBkgSig","fBkgSig"          ,fBkgx,fSigy);
    }
    RooProdPdf          fSS     ("fSigSig","fSigSig"    ,fSigx,fSigy);
    RooAddPdf           fMod    ("fMod2D","fMod2D"      ,RooArgList(*fBB,fSS,*fSB,*fBS),RooArgList(n1,n0,n3,n2));
    
    RooFitResult* fFitResults;
    fFitResults      =   fMod.fitTo(*dataLoose,Save(),NumCPU(kCPU_use,kCPUStrategy));
    for ( Int_t iCycle = 0; iCycle < kNCycle; iCycle++ )
    {
        fFitResults      =   fMod.fitTo(*data,Extended(kTRUE),Save(),NumCPU(kCPU_use,kCPUStrategy),Offset(kFitOffset),Strategy(kFitMinuitStrategy),InitialHesse(kFitInitHesse),Minos(kFitMinos));
        auto N_Raw      =   static_cast<RooRealVar*>(fFitResults ->floatParsFinal().find("anSS2D"));
        if ( fIsResultAcceptable(N_Raw->getVal(),N_Raw->getError(),k2DErrorLimit) ) break;
    }
    // Save to file
    if ( fSaveToFile || kSaveToFile )
    {
        SetStyle();
        
        int         nBinsPrint      =   4;
        double      dIncrement      =   0.014;//(fInvMassValMax-fInvMassValMin)/nBinsPrint;
        TCanvas*    cTotal          =   new TCanvas("","",0,45,1800,1400);
                    cTotal          ->  SetTitle(Form("Slices of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]));
                    cTotal          ->  SetName(Form("PT_%.1f_%.1f__%.1f_%.1f_%s",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1],fName.Data()));
                    cTotal          ->  Divide(2,nBinsPrint);
        
                            xInvMass.setRange("fDrawRange",fInvMassValMin,fInvMassValMax);
                            yInvMass.setRange("fDrawRange",fInvMassValMin,fInvMassValMax);
        
        for (int i = 0; i < nBinsPrint; i++)
        {
            hName                       = "Slice of 2D Invariant Mass of Kaons";
            hTitle                      = "Slice of 2D Invariant Mass of Kaons";
            if ( PTindex != -1 && !bPythiaTest ) hTitle = Form("Slice of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]);
            if ( PTindex != -1 &&  bPythiaTest ) hTitle = Form("Slice of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV for MC",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]);
            
            TCanvas * fSaveToCanvas =   new TCanvas(
                                                    Form("xInvMass_%.3f_%.3f_PTx_%.3f_%.3f_PTy_%.3f_%.3f_%s",fInvMassValMin+dIncrement*i,fInvMassValMin+dIncrement*(i+1),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1],fName.Data()),
                                                    Form("xInvMass_%.3f_%.3f_PTx_%.3f_%.3f_PTy_%.3f_%.3f",fInvMassValMin+dIncrement*i,fInvMassValMin+dIncrement*(i+1),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1])
                                                    );
            
            RooPlot * fSaveToFrame  =   yInvMass.frame(Name(hName),Title(hTitle));
            TLegend * fLegend           = new TLegend   (0.18,0.85,0.33,0.60);

                            xInvMass.setRange("fDrawRange",fInvMassValMin+i*dIncrement,fInvMassValMin+(i+1)*dIncrement);
                            yInvMass.setRange("fDrawRange",fInvMassValMin,fInvMassValMax);

            fRooPlotMaker(fSaveToFrame,fLegend,fMod,data,"yInvMass2D");
            
            fSaveToCanvas->cd();
            
            fSaveToFrame                ->GetYaxis()->SetTitle(Form("Counts/( %.1f MeV/#it{c}^{2} )",1000*kBinningPrecision2D));
            
            auto fMaximum   =   fSaveToFrame->GetMaximum();
            fSaveToFrame    ->  SetMaximum(fMaximum*1.30);
            
            fSaveToFrame                ->Draw("same");
            fLegend                     ->Draw("same");
            
            uLatex->SetTextFont(60);
            uLatex->SetTextSize(0.05);
            uLatex->DrawLatexNDC(0.50, 0.83,"ALICE Performance");
            uLatex->SetTextFont(42);
            uLatex->SetTextSize(0.035);
            uLatex->DrawLatexNDC(0.50, 0.77,"pp #sqrt{#it{s}}= 7 TeV, #phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
            uLatex->DrawLatexNDC(0.50, 0.72,Form("%.2f < #it{p}_{T,#phi_{1}} < %.2f GeV/#it{c}",fArrPT2D[PTindex],fArrPT2D[PTindex+1]));
            uLatex->DrawLatexNDC(0.51, 0.67,Form("%.2f < #it{p}_{T,#phi_{2}} < %.2f GeV/#it{c}",fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]));
            uLatex->DrawLatexNDC(0.50, 0.62,Form("%.3f < M_{KK,#phi_{1}} < %.3f GeV/#it{c}^{2}",fInvMassValMin+dIncrement*i,fInvMassValMin+dIncrement*(i+1)));
            
            fSaveToCanvas               ->Write();
            fSaveToCanvas               ->SaveAs(Form("%s/PT_%.1f_%.1f__%.1f_%.1f_%s_%i.pdf",fPathToSave.Data(),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1],fName.Data(),i));
            fSaveToCanvas               ->SaveAs(Form("%s/PT_%.1f_%.1f__%.1f_%.1f_%s_%i.eps",fPathToSave.Data(),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1],fName.Data(),i));
            
            cTotal->cd( 2*i+1 );
            fSaveToFrame                ->Draw("same");
            if ( i == 0 )   {
                fLegend                     ->Draw("same");
                uLatex->SetTextFont(60);
                uLatex->SetTextSize(0.05);
                uLatex->DrawLatexNDC(0.55, 0.83,"ALICE Performance");
                uLatex->SetTextFont(42);
                uLatex->SetTextSize(0.04);
                uLatex->DrawLatexNDC(0.55, 0.77,"pp #sqrt{#it{s}}= 7 TeV, #phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
                uLatex->DrawLatexNDC(0.55, 0.72,Form("%.2f < #it{p}_{T,#phi_{1}} < %.2f GeV/#it{c}",fArrPT2D[PTindex],fArrPT2D[PTindex+1]));
                uLatex->DrawLatexNDC(0.56, 0.67,Form("%.2f < #it{p}_{T,#phi_{2}} < %.2f GeV/#it{c}",fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]));
                uLatex->DrawLatexNDC(0.55, 0.62,Form("%.3f < M_{KK,#phi_{1}} < %.3f GeV/#it{c}^{2}",fInvMassValMin+dIncrement*i,fInvMassValMin+dIncrement*(i+1)));
            }   else    {
                uLatex->DrawLatexNDC(0.55, 0.83,Form("%.3f < M_{KK,#phi_{1}} < %.3f GeV/#it{c}^{2}",fInvMassValMin+dIncrement*i,fInvMassValMin+dIncrement*(i+1)));
            }
            
            delete fSaveToCanvas;
        }
                                        xInvMass.setRange("fDrawRange",fInvMassValMin,fInvMassValMax);
                                        yInvMass.setRange("fDrawRange",fInvMassValMin,fInvMassValMax);
        for (int i = 0; i < nBinsPrint; i++)
        {
            hName                       = "Slice of 2D Invariant Mass of Kaons";
            hTitle                      = "Slice of 2D Invariant Mass of Kaons";
            if ( PTindex != -1 && !bPythiaTest ) hTitle = Form("Slice of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]);
            if ( PTindex != -1 &&  bPythiaTest ) hTitle = Form("Slice of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV for MC",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]);
            
            TCanvas * fSaveToCanvas =   new TCanvas(
                                                    Form("yInvMass_%.3f_%.3f_PTx_%.3f_%.3f_PTy_%.3f_%.3f_%s",fInvMassValMin+dIncrement*i,fInvMassValMin+dIncrement*(i+1),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1],fName.Data()),
                                                    Form("yInvMass_%.3f_%.3f_PTx_%.3f_%.3f_PTy_%.3f_%.3f",fInvMassValMin+dIncrement*i,fInvMassValMin+dIncrement*(i+1),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1])
                                                    );
            
            RooPlot * fSaveToFrame      =   xInvMass.frame(Name(hName),Title(hTitle));
            TLegend * fLegend           = new TLegend   (0.18,0.85,0.33,0.60);
            
                                        xInvMass.setRange("fDrawRange",fInvMassValMin,fInvMassValMax);
                                        yInvMass.setRange("fDrawRange",fInvMassValMin+i*dIncrement,fInvMassValMin+(i+1)*dIncrement);
                                                                            
            fRooPlotMaker(fSaveToFrame,fLegend,fMod,data,"xInvMass2D");
            
            fSaveToCanvas->cd();
            
            fSaveToFrame                ->GetYaxis()->SetTitle(Form("Counts/( %.1f MeV/#it{c}^{2} )",1000*kBinningPrecision2D));
            
            auto fMaximum   =   fSaveToFrame->GetMaximum();
            fSaveToFrame    ->  SetMaximum(fMaximum*1.30);
            
            fSaveToFrame                ->Draw("same");
            fLegend                     ->Draw("same");
            
            uLatex->SetTextFont(60);
            uLatex->SetTextSize(0.05);
            uLatex->DrawLatexNDC(0.50, 0.83,"ALICE Performance");
            uLatex->SetTextFont(42);
            uLatex->SetTextSize(0.035);
            uLatex->DrawLatexNDC(0.50, 0.77,"pp #sqrt{#it{s}}= 7 TeV, #phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
            uLatex->DrawLatexNDC(0.50, 0.72,Form("%.2f < #it{p}_{T,#phi_{1}} < %.2f GeV/#it{c}",fArrPT2D[PTindex],fArrPT2D[PTindex+1]));
            uLatex->DrawLatexNDC(0.51, 0.67,Form("%.2f < #it{p}_{T,#phi_{2}} < %.2f GeV/#it{c}",fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]));
            uLatex->DrawLatexNDC(0.50, 0.62,Form("%.3f < M_{KK,#phi_{2}} < %.3f GeV/#it{c}^{2}",fInvMassValMin+dIncrement*i,fInvMassValMin+dIncrement*(i+1)));
            
            fSaveToCanvas               ->Write();
            fSaveToCanvas               ->SaveAs(Form("%s/PT_%.1f_%.1f__%.1f_%.1f_%s_%i.pdf",fPathToSave.Data(),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1],fName.Data(),i+nBinsPrint));
            fSaveToCanvas               ->SaveAs(Form("%s/PT_%.1f_%.1f__%.1f_%.1f_%s_%i.eps",fPathToSave.Data(),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1],fName.Data(),i+nBinsPrint));
            
            cTotal->cd( 2*i+2 );
            fSaveToFrame                ->Draw("same");
            uLatex->DrawLatexNDC(0.50, 0.83,Form("%.3f < M_{KK,#phi_{2}} < %.3f GeV/#it{c}^{2}",fInvMassValMin+dIncrement*i,fInvMassValMin+dIncrement*(i+1)));
            
            delete fSaveToCanvas;
        }
                                        xInvMass.setRange("fDrawRange",fInvMassValMin,fInvMassValMax);
                                        yInvMass.setRange("fDrawRange",fInvMassValMin,fInvMassValMax);
        cTotal ->Write();
        cTotal               ->SaveAs(Form("%s/PT_%.1f_%.1f__%.1f_%.1f_%s.pdf",fPathToSave.Data(),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1],fName.Data()));
        cTotal               ->SaveAs(Form("%s/PT_%.1f_%.1f__%.1f_%.1f_%s.eps",fPathToSave.Data(),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1],fName.Data()));
        delete cTotal;
    }
    
    // Un-Silencing TCanvas Pop-Up
    gROOT->SetBatch(false);
    
    // Fit
    return fFitResults;
}
//
std::vector<TH2F*>
FitModel
 ( TH1F  **hShapeFit, TH1F* hSlopReference, TH2F***hTarget, RooFitResult*** &fFitresultsStore = NULL_ROOFITPTR3, TString fOption = "", TString fTargetPath = "./result/SEFitCheck", TString fNameFile = "", std::vector<TH1F*> &f1Din2DbinCheck = NULL_VECTOR )  {
    std::vector<TH2F*>  fResults;
    RooFitResult      **fShapeStore = nullptr;
    if ( !fFitresultsStore )    {
        fFitresultsStore = new RooFitResult**[nBinPT2D];
        for ( Int_t iTer = 0; iTer < nBinPT2D; iTer++ ) {
            fFitresultsStore[iTer]  =    new RooFitResult    *[nBinPT2D];
        }
    }
    //
    //>>    Recover 1D Shapes
    f1Din2DbinCheck     =   FitModel(hShapeFit,hSlopReference,fShapeStore,2,fOption,fTargetPath,fNameFile);
    //
    //>>    Fit 2D Histograms
    for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ ) {
        for ( Int_t jFit = iFit; jFit < nBinPT2D; jFit++ ) {
            //
            //>>    Protection Against SegFault
            if ( !fShapeStore[iFit] ) continue;
            if ( !fShapeStore[jFit] ) continue;
            //
            //>>    Fit
            fFitresultsStore[iFit][jFit]       =  fFitCoreModel(hTarget[iFit][jFit],hSlopReference,fShapeStore[iFit],fShapeStore[jFit],fNameFile,fOption,iFit,jFit,fTargetPath);
            //
            //>>    Protection Against SegFault
            if ( !fFitresultsStore[iFit][jFit] ) continue;
            //
            //>>    Building Raw Count histograms
            if ( iFit == 0 && jFit == 0 )   {
                Int_t   iTer = 0;
                for ( auto fCoeff : fFitresultsStore[iFit][jFit]->floatParsFinal() )   {
                    auto N_Raw      = static_cast<RooRealVar*>(fCoeff);
                    fResults.push_back( new TH2F(Form("%s_%s",N_Raw->GetName(),fNameFile.Data()),N_Raw->GetName(),nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D) );
                    if ( strncmp(N_Raw->GetName(),"anSS2D",6) == 0 ) fResults.at(iTer)->SetName("hRAW_2D");
                    iTer++;
                }
                for ( auto fCoeff : fFitresultsStore[iFit][jFit]->constPars() )   {
                    auto N_Raw      = static_cast<RooRealVar*>(fCoeff);
                    fResults.push_back( new TH2F(Form("%s_%s",N_Raw->GetName(),fNameFile.Data()),N_Raw->GetName(),nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D) );
                }
            }
            //
            //>>    Filling Raw Count Histograms
            Int_t   iTer = 0;
            for ( auto fCoeff : fFitresultsStore[iFit][jFit]->floatParsFinal() )   {
                auto N_Raw      = static_cast<RooRealVar*>(fCoeff);
                fResults.at(iTer)->SetBinContent          (iFit+1,jFit+1,N_Raw->getVal());
                fResults.at(iTer)->SetBinError            (iFit+1,jFit+1,N_Raw->getError());
                fResults.at(iTer)->SetBinContent          (jFit+1,iFit+1,N_Raw->getVal());
                fResults.at(iTer)->SetBinError            (jFit+1,iFit+1,N_Raw->getError());
                if ( iFit == jFit ) {
                    fResults.at(iTer)->SetBinContent          (iFit+1,jFit+1,2.*N_Raw->getVal());
                    fResults.at(iTer)->SetBinError            (iFit+1,jFit+1,2.*N_Raw->getError());
                }
                iTer++;
            }
        }
    }
    return fResults;
}
std::vector<TH2F*>
FitModel
 ( TH2F***hTarget, TH1F* hSlopReference,  RooFitResult**  fShapeStore, RooFitResult*** &fFitresultsStore = NULL_ROOFITPTR3, TString fOption = "", TString fTargetPath = "./result/SEFitCheck", TString fNameFile = "", std::vector<TH1F*> &f1Din2DbinCheck = NULL_VECTOR )  {
    std::vector<TH2F*>  fResults;
    if ( !fFitresultsStore )    {
        fFitresultsStore = new RooFitResult**[nBinPT2D];
        for ( Int_t iTer = 0; iTer < nBinPT2D; iTer++ ) {
            fFitresultsStore[iTer]  =    new RooFitResult    *[nBinPT2D];
        }
    }
    //
    //>>    Fit 2D Histograms
    for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ ) {
        for ( Int_t jFit = 0; jFit < nBinPT2D; jFit++ ) {
            //
            //>>    Protection Against SegFault
            if ( !fShapeStore[iFit] ) continue;
            if ( !fShapeStore[jFit] ) continue;
            //
            //>>    Fit
            fFitresultsStore[iFit][jFit]       =  fFitCoreModel(hTarget[iFit][jFit],hSlopReference,fShapeStore[iFit],fShapeStore[jFit],fNameFile,fOption,iFit,jFit,fTargetPath);
            //
            //>>    Protection Against SegFault
            if ( !fFitresultsStore[iFit][jFit] ) continue;
            //
            //>>    Building Raw Count histograms
            if ( iFit == 0 && jFit == 0 )   {
                Int_t   iTer = 0;
                for ( auto fCoeff : fFitresultsStore[iFit][jFit]->floatParsFinal() )   {
                    auto N_Raw      = static_cast<RooRealVar*>(fCoeff);
                    fResults.push_back( new TH2F(N_Raw->GetName(),N_Raw->GetName(),nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D) );
                    if ( strncmp(N_Raw->GetName(),"anSS2D",6) == 0 ) fResults.at(iTer)->SetName("hRAW_2D");
                    iTer++;
                }
            }
            //
            //>>    Filling Raw Count Histograms
            Int_t   iTer = 0;
            for ( auto fCoeff : fFitresultsStore[iFit][jFit]->floatParsFinal() )   {
                auto N_Raw      = static_cast<RooRealVar*>(fCoeff);
                fResults.at(iTer)->SetBinContent          (iFit+1,jFit+1,N_Raw->getVal());
                fResults.at(iTer)->SetBinError            (iFit+1,jFit+1,N_Raw->getError());
                iTer++;
            }
        }
    }
    return fResults;
}
std::vector<TH2F*>
FitModel                        ( TH1F  **hShapeFit, TH1F* hSlopReference, TH2F***hTarget, TString fTargetPath = "./result/SEFitCheck", TString fNameFile = "", TString fOption = "",  std::vector<TH1F*> &f1Din2DbinCheck = NULL_VECTOR )  {
    return  FitModel(hShapeFit,hSlopReference,hTarget,NULL_ROOFITPTR3,fOption,fTargetPath,fNameFile);
}
std::vector<TH2F*>
FitModel
 ( TH1F  **hShapeFit, TH1F* hSlopReference, TH2F***hTarget, std::vector<TH1F*> &f1Din2DbinCheck, TString fTargetPath = "./result/SEFitCheck", TString fNameFile = "", TString fOption = "" )  {
    return  FitModel(hShapeFit,hSlopReference,hTarget,NULL_ROOFITPTR3,fOption,fTargetPath,fNameFile,f1Din2DbinCheck);
}
std::vector<TH2F*>
FitModel
 ( TH2F***hTarget, TH1F* hSlopReference, RooFitResult**  fShapeStore, std::vector<TH1F*> &f1Din2DbinCheck, TString fTargetPath = "./result/SEFitCheck", TString fNameFile = "" )  {
    return  FitModel(hTarget,hSlopReference,fShapeStore,NULL_ROOFITPTR3,"",fTargetPath,fNameFile,f1Din2DbinCheck);
}
std::vector<TH2F*>
FitModel
 ( TH2F***hTarget, TH1F* hSlopReference, RooFitResult**  fShapeStore, std::vector<TH1F*> &f1Din2DbinCheck, TString fTargetPath = "./result/SEFitCheck", TString fNameFile = "", TString fOption = "" )  {
    return  FitModel(hTarget,hSlopReference,fShapeStore,NULL_ROOFITPTR3,fOption,fTargetPath,fNameFile,f1Din2DbinCheck);
}
//
//_____________________________________________________________________________

//------------------------------//
//    !TODO: Clean the mess below      //
//------------------------------//
//
template < class Tclass >
void    SetAxis             ( Tclass * aTarget, string aOption = "" )   {
    if ( aOption.find("IM") != -1 )
    {
        if ( aOption.find("2D") != -1 )
        {
            // X-Axis formatting
            aTarget->GetXaxis()->SetTitle("M_{KK,#phi_{1}} (GeV/#it{c}^{2})");
            aTarget->GetXaxis()->SetTitleOffset(1.15);
            
            // Y-Axis formatting
            aTarget->GetYaxis()->SetTitle("M_{KK,#phi_{2}} (GeV/#it{c}^{2})");
            aTarget->GetYaxis()->SetTitleOffset(1.15);
        }
        else if ( aOption.find("1D") != -1 )
        {
            // X-Axis formatting
            aTarget->GetXaxis()->SetTitle("M_{KK} (GeV/#it{c}^{2})");
            aTarget->GetXaxis()->SetTitleOffset(1.15);
        }
    }
    else if ( aOption.find("PT") != -1 )
    {
        if ( aOption.find("DD") != -1 )
        {
            // X-Axis formatting
            aTarget->GetXaxis()->SetTitle("#it{p}_{T,#phi_{2}} (GeV/#it{c})");
            aTarget->GetXaxis()->SetTitleOffset(1.15);
            
            // Y-Axis formatting
            aTarget->GetYaxis()->SetTitle("#frac{d^{2}N_{#phi#phi}}{dydp_{T}#phi_{2}}(GeV/c)^{-1}");
            aTarget->GetYaxis()->SetTitleOffset(1.15);
        }
        else if ( aOption.find("2D") != -1 )
        {
            // X-Axis formatting
            aTarget->GetXaxis()->SetTitle("#it{p}_{T,#phi_{1}} (GeV/#it{c})");
            aTarget->GetXaxis()->SetTitleOffset(1.15);
                
            // Y-Axis formatting
            aTarget->GetYaxis()->SetTitle("#it{p}_{T,#phi_{2}} (GeV/#it{c})");
            aTarget->GetYaxis()->SetTitleOffset(1.15);
                
            // Z-Axis formatting
            //aTarget->GetZaxis()->SetTitle("#frac{d^{3}N_{#phi#phi}}{dydp_{T}#phi_{1}dp_{T}#phi_{2}}(GeV/c)^{-1}");
            //aTarget->GetZaxis()->SetTitleOffset(1.15);
        }
        else if ( aOption.find("1D") != -1 )
        {
            // X-Axis formatting
            aTarget->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            aTarget->GetXaxis()->SetTitleOffset(1.15);
            
            // Y-Axis formatting
            aTarget->GetYaxis()->SetTitle("#frac{d^{2}N_{#phi}}{dydp_{T}}(GeV/c)^{-1}");
            aTarget->GetYaxis()->SetTitleOffset(1.15);
        }
    }
}
//
//--
//
template < class Tclass >
Double_t        HistIntegrals               ( Tclass* aTarget, string aOption = "", Double_t aLowBound = 0, Double_t aHigBound = 0 )    {
    
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
    for ( Int_t fiTer = 1; fiTer <= aTarget->GetNbinsX(); fiTer++ )
    {
        // Determining the pT point to calculate
        auto fPoint =   aTarget->GetBinCenter(fiTer);
        
        // Exiting the loop if out of bound
        if ( aTarget->GetBinLowEdge(fiTer)   <= aLowBound && aLowBound != aHigBound ) continue;
        if ( aTarget->GetBinLowEdge(fiTer+1) >= aHigBound && aLowBound != aHigBound ) break;
        
        Double_t        fCycleAdd   =   ( aTarget->GetBinContent(fiTer) );
        if ( fError )   fCycleAdd   =   ( aTarget->GetBinError(fiTer) );
        if ( fWidth )   fCycleAdd  *=   ( aTarget->GetBinWidth(fiTer) );
        if ( fPower )   fCycleAdd  *=   pow( ( aTarget->GetBinCenter(fiTer) ) , fnPowr );
        if ( fError )   fCycleAdd  *=   fCycleAdd;
        fResult += fCycleAdd;
    }
    
    if ( fError )   return sqrt(fResult);
    return fResult;
}
//
//------------------------------//
//    ANALYSISI SEPCIFIC Fncs   //
//------------------------------//
//
//_____________________________________________________________________________
//
Int_t kColor[] = {38,kBlue,kBlue+3,46,38};
Int_t kStyle[] = {26,9,10,25,22};
Int_t kWidth[] = {1,3,3,1,1};
//
//_____________________________________________________________________________
//
void                fChooseOption                   ( TString fOption ) {
    if ( fOption.IsNull() ) cout << "[INFO] No option chosen, standard all inclusive analysis enabled" <<endl;
    else    {
        kDoMultiplicity         =   false;
        kDoYield                =   false;
        kDoTrigger              =   false;
        kDoRapidity             =   false;
        if ( fOption.Contains("Multiplicity",TString::kIgnoreCase) )    { kDoMultiplicity = true;   cout << "[INFO] Multiplicity option chosen" <<endl;}
        if ( fOption.Contains("Yield",TString::kIgnoreCase) )           { kDoYield = true;          cout << "[INFO] Yield option chosen" <<endl;}
        if ( fOption.Contains("Trigger",TString::kIgnoreCase) )         { kDoTrigger = true;        cout << "[INFO] Trigger option chosen" <<endl;}
        if ( fOption.Contains("Rapidity",TString::kIgnoreCase) )        { kDoRapidity = true;       cout << "[INFO] Rapidity option chosen" <<endl;}
    }
}
//
//_____________________________________________________________________________
//
Double_t
fGammaPhiValue
 ( Double_t fYieldPhi, Double_t fYieldPhiPhi )  {
    return  2*fYieldPhiPhi/fYieldPhi -fYieldPhi;
}
Double_t
fGammaPhiError
 ( Double_t fYieldPhi, Double_t fYieldPhiPhi, Double_t fErrorPhi, Double_t fErrorPhiPhi)  {
    auto    fPar1   =   2*fErrorPhiPhi/fYieldPhi;
    auto    fPar2   =   (2*fYieldPhiPhi/(fYieldPhi*fYieldPhi)+1)*fErrorPhi;
    return  fPar1 + fPar2;
}
Double_t
fSigmaPhiValue
 ( Double_t fYieldPhi, Double_t fYieldPhiPhi )  {
    return  2*fYieldPhiPhi + fYieldPhi - fYieldPhi*fYieldPhi;
}
Double_t
fSigmaPhiError
 ( Double_t fYieldPhi, Double_t fYieldPhiPhi, Double_t fErrorPhi, Double_t fErrorPhiPhi)  {
    return SquareSum( { 2*fErrorPhiPhi, (-1+fYieldPhi)*fErrorPhi } );
}
//
//_____________________________________________________________________________
//
void                fSetFunction                    ( TF1* fFitFunction = fLevyTsallis, Double_t fIntegral = 0.032 ) {
    fSetAllFunctions();
    //
    //____________________________________________Mass
    auto    nMass   =   fFitFunction    ->GetParNumber("Mass");
    if ( nMass != -1 )   {
        fFitFunction    ->  FixParameter (nMass,kPhiMesonMass_);
    }
    //
    //____________________________________________dN/dy
    auto    ndNdy   =   fFitFunction    ->GetParNumber("dN_dy");
    if ( ndNdy != -1 )   {
        fFitFunction    ->  SetParLimits(ndNdy,1.e-12,1.);
        fFitFunction    ->  SetParameter(ndNdy,fIntegral);
    }
    //
    //____________________________________________dN/dy
    auto    nEnne   =   fFitFunction    ->GetParNumber("n");
    if ( nEnne != -1 )   {
        fFitFunction    ->  SetParLimits(nEnne,0.,100.);
        fFitFunction    ->  SetParameter(nEnne,7.18);
        if ( strncmp(fFitFunction->GetName(),"LevyTsallis",11) == 0 )   {
            fFitFunction  ->  SetParLimits(nEnne,2.1,15.);
            fFitFunction  ->  SetParameter(nEnne,7.18);
        }
        if ( strncmp(fFitFunction->GetName(),"PowerLaw",13) == 0 )          fFitFunction  ->  SetParameter(nEnne,13);
    }
    //
    //____________________________________________dN/dy
    auto    nSlop   =   fFitFunction    ->GetParNumber("T");
    if ( nSlop != -1 )   {
        fFitFunction    ->  SetParLimits(nSlop,0.05,30);
        fFitFunction    ->  SetParameter(nSlop,.330);
        if ( strncmp(fFitFunction->GetName(),"LevyTsallis",11) == 0 )       fFitFunction  ->  SetParameter(nSlop,.300);
        if ( strncmp(fFitFunction->GetName(),"MTExponential",13) == 0 )     fFitFunction  ->  SetParameter(nSlop,.560);
        if ( strncmp(fFitFunction->GetName(),"PTExponential",13) == 0 )     fFitFunction  ->  SetParameter(nSlop,.560);
        if ( strncmp(fFitFunction->GetName(),"FermiDirac",13) == 0 )        fFitFunction  ->  SetParameter(nSlop,.330);
        if ( strncmp(fFitFunction->GetName(),"Boltzmann",13) == 0 )         fFitFunction  ->  SetParameter(nSlop,.350);
        if ( strncmp(fFitFunction->GetName(),"BoseEinstein",13) == 0 )      fFitFunction  ->  SetParameter(nSlop,.350);
        if ( strncmp(fFitFunction->GetName(),"PowerLaw",13) == 0 )          fFitFunction  ->  SetParameter(nSlop,.550);
    }
}
void                fSetFunction                    ( TF2* fFitFunction , Double_t fIntegral = 3.e3 ) {
    fSetAllFunctions();
    //
    //____________________________________________Mass
    auto    xMass   =   fFitFunction    ->GetParNumber("XMass");
    auto    yMass   =   fFitFunction    ->GetParNumber("YMass");
    fFitFunction    ->  FixParameter (xMass,kPhiMesonMass_);
    fFitFunction    ->  FixParameter (yMass,kPhiMesonMass_);
    //
    //____________________________________________dN/dy
    auto    XdNdy   =   fFitFunction    ->GetParNumber("XdN_dy");
    auto    YdNdy   =   fFitFunction    ->GetParNumber("YdN_dy");
    auto    _dNdy   =   fFitFunction    ->GetParNumber("dN_dy");
    fFitFunction    ->  SetParLimits(XdNdy,0.,1.e6);
    fFitFunction    ->  SetParameter(XdNdy,fIntegral);
    fFitFunction    ->  SetParLimits(YdNdy,0.,1.e6);
    fFitFunction    ->  SetParameter(YdNdy,fIntegral);
    fFitFunction    ->  SetParLimits(_dNdy,0.,1.e6);
    fFitFunction    ->  SetParameter(_dNdy,fIntegral);
    //
    //____________________________________________dN/dy
    auto    XEnne   =   fFitFunction    ->GetParNumber("Xn");
    auto    YEnne   =   fFitFunction    ->GetParNumber("Yn");
    fFitFunction    ->  SetParLimits(XEnne,0.,100.);
    fFitFunction    ->  SetParameter(XEnne,7.18);
    fFitFunction    ->  SetParLimits(YEnne,0.,100.);
    fFitFunction    ->  SetParameter(YEnne,7.18);
    //
    //____________________________________________dN/dy
    auto    xSlop   =   fFitFunction    ->GetParNumber("XT");
    auto    ySlop   =   fFitFunction    ->GetParNumber("YT");
    fFitFunction    ->  SetParLimits(xSlop,0.05,30);
    fFitFunction    ->  SetParameter(xSlop,.330);
    fFitFunction    ->  SetParLimits(ySlop,0.05,30);
    fFitFunction    ->  SetParameter(ySlop,.330);
}

//
//_____________________________________________________________________________
//
void                fFitLevyTsalis                  ( TGraphAsymmErrors* gToBeFitted ) {
    return;
}
//
//_____________________________________________________________________________
//
template < class Tclass >
Tclass*
fSetSystErrors
 ( Tclass*                          hStatistics ) {
    auto        fResult         =   (Tclass*)(hStatistics->Clone());
    //
    TFile      *fReference      =   new TFile( Form("%s/Full_Systematics.root",(TString(Form(kAnalysis_Systemt_Dir,"yield"))).Data()) );
    Tclass     *hReference      =   (Tclass*)(fReference->Get("h1DTotalSystematic"));
    //
    for ( Int_t iBin = 0; iBin < fResult->GetNbinsX(); iBin++ ) {
        auto    fBinContent     =   fResult->GetBinContent( iBin );
        fResult->SetBinError( iBin, fBinContent*0.1);//(hReference->GetBinContent(iBin)));
    }
    fReference->Close();
    //
    return  fResult;
}
template < class Tclass >
std::vector<Tclass*>
fSetSystErrors
 ( std::vector<Tclass*>             hStatistics ) {
    std::vector<Tclass*> fResult;
    //
    TFile      *fReference      =   new TFile( Form("%s/Full_Systematics.root",(TString(Form(kAnalysis_Systemt_Dir,"yield"))).Data()) );
    TH2F       *hReference      =   (TH2F*)(fReference->Get("h2DTotalSystematic"));
    //
    auto jBin = 0;
    for ( auto& iHisto : hStatistics )  {
        auto    fCurrentHisto   =   new Tclass(*iHisto);
        for ( Int_t iBin = 0; iBin < iHisto->GetNbinsX(); iBin++ ) {
            auto    fBinContent     =   fCurrentHisto->GetBinContent( iBin+1 );
            fCurrentHisto->SetBinError( iBin+1, fBinContent*0.1);//(hReference->GetBinContent(jBin+1,iBin+1)));
        }
        jBin++;
        fResult.push_back(fCurrentHisto);
    }
    return fResult;
}
//
template < class Tclass >
Tclass*
fSetSystErrors
( Tclass* hTarget,                  Int_t jBin ) {
    Tclass* hResult         =   (Tclass*)(hTarget->Clone());
    TFile*  fReference      =   new TFile( Form("%s/Full_Systematics.root",(TString(Form(kAnalysis_Systemt_Dir,"yield"))).Data()) );
    TH2F*   hReference      =   (TH2F*)(fReference->Get("h2DTotalSystematic"));
    for ( Int_t iBin = 0; iBin < hTarget->GetNbinsX(); iBin++ ) {
        hResult->SetBinError    ( iBin+1, (hTarget->GetBinContent( iBin+1 ))*0.1);//(hReference->GetBinContent(jBin+1,iBin+1)) );
    }
    return hResult;
}
//
//_____________________________________________________________________________
//
//!    TODO: Generalise and move to alianalysis utility, add function from outside
template < class Tclass >
std::vector<float>
fEvaluateExtrapolationUncertainty
 ( bool fIsConditional, Tclass* fStaticUncertaintySpectrum, Tclass* fVariableUncertaintySpectrum, Float_t fIntegral, Double_t fMaximumFitRange = fMaxPT1D, Double_t fMinimumFitRange = fMinPT1D, TString fName = "", TString fFolder = "", Int_t kEvaluationCycles = kStatEvalCycles )  {
    //
    //  Result container
    //  Format: [0] Stat Extrap Unc //! TODO: update
    std::vector<float> fResult;
    //
    //  Build the Total Uncertainty Spectrum
    Tclass*     fTotalUncertaintySpectrum   =   fSumErrors( fStaticUncertaintySpectrum, fVariableUncertaintySpectrum );
                fTotalUncertaintySpectrum   ->  SetName("TotalUncertaintySpectrum");
    //
    //  Measure the baselines
    fTotalUncertaintySpectrum       ->  Fit( fLevyTsallis, "Q",         "R", fMinimumFitRange,  fMaximumFitRange );
    fTotalUncertaintySpectrum       ->  Fit( fLevyTsallis, "QMES",      "R", fMinimumFitRange,  fMaximumFitRange );
    fTotalUncertaintySpectrum       ->  Fit( fLevyTsallis, "QMESI",     "R", fMinimumFitRange,  fMaximumFitRange );
    //
    auto    fStandardIntegral       =   fTotalUncertaintySpectrum->Integral("width");
    auto    fStandardExtrapolation  =   fLevyTsallis->Integral(0.,fMinPT1D);
    auto    fStandardYield          =   fStandardIntegral+fStandardExtrapolation;
    auto    fStandardMeanPT         =   fEvaluateMeanPT( fTotalUncertaintySpectrum, fLevyTsallis );
    //
    //  Canvas to draw the Results checks
    SetStyle();
    auto    cDrawResults            =   new TCanvas("cDrawResult","cDrawResult",1800,600);
    gStyle                          ->  SetOptStat(0);
    cDrawResults                    ->  Divide(3,1);
    cDrawResults                    ->  cd(1);
    gPad                            ->  SetLogy();
    //
    //  Draw Spectra
    fTotalUncertaintySpectrum       ->  SetMarkerStyle(markers[2]);
    fTotalUncertaintySpectrum       ->  SetMarkerSize(2);
    fTotalUncertaintySpectrum       ->  SetMarkerColor(colors[0]);
    fTotalUncertaintySpectrum       ->  SetLineColor(kBlack);
    fTotalUncertaintySpectrum       ->  DrawClone("EP");
    //
    //  Generate utility histogram to evaluate statistical error
    std::vector<float> fYieldVariation;
    std::vector<float> fPTVariation;
    //
    //  Loop requested times
    //      Protection against unreasonable values
    if ( kEvaluationCycles < 2 )    {
        kEvaluationCycles = kStatEvalCycles;
        cout << "[WARNING] Evaluation cycles must be at least 2, set to default value";
        if ( kEvaluationCycles < 2 ) {
            kEvaluationCycles = 10;
            cout << "[WARNING] Default evaluation cycles is invalid, set to 10";
        }
    }
    for ( Int_t iFit = 0; iFit < kEvaluationCycles; iFit++ )  {
        //
        //  Dump Canvas
        auto    cDump   =   new TCanvas();
        //  Randomise points
        auto    fCurrentUtilitySpectrum =   uRandomizePoints( fStaticUncertaintySpectrum, fVariableUncertaintySpectrum );
        //
        //  //! TODO: Looks like a cheap fix
        //
        //  Set Standard Fit
        fSetFunction        (fLevyTsallis,  fCurrentUtilitySpectrum->Integral());
        fCurrentUtilitySpectrum    ->  Fit (fLevyTsallis,  "Q");
        fCurrentUtilitySpectrum    ->  Fit (fLevyTsallis,  "QEMS",  "", fMinimumFitRange,fMaximumFitRange);
        fCurrentUtilitySpectrum    ->  Fit (fLevyTsallis,  "QEMSI", "", fMinimumFitRange,fMaximumFitRange);
        //
        delete  cDump;
        //
        fYieldVariation.push_back( fStandardIntegral + fLevyTsallis->Integral(0.,fMinPT1D) );
        //
        fPTVariation.push_back( fEvaluateMeanPT( fTotalUncertaintySpectrum, fLevyTsallis ) );
        //
        //  Draw Function
        cDrawResults->cd(1);
        fCurrentUtilitySpectrum->GetFunction("LevyTsallis")->DrawCopy("SAME");
    }
    //
    cDrawResults                ->  cd(2);
    //
    //  Build the Variation histogram, fixed to 50 bins
    auto    hYieldVariations    =   uBuildTH1F( fYieldVariation, 50 );
    hYieldVariations            ->  SetNameTitle( Form("hYieldVariations_%s",fName.Data()), "dN/dy variations" );
    //
    //  Gauss Fit
    hYieldVariations            ->  Fit( "gaus", "IEMQS" );
    auto fErrorMean             =   hYieldVariations->GetFunction("gaus")->GetParameter(1);
    auto fErrorSTDV             =   hYieldVariations->GetFunction("gaus")->GetParameter(2);
    auto fErrorFULL             =   fabs( fErrorMean - fStandardYield ) + fErrorSTDV;
    //
    //  Save in Results
    fResult.push_back(fErrorFULL);
    fResult.push_back(fErrorFULL);
    //
    //  Draw the results
    hYieldVariations            ->  SetMaximum( 1.4*hYieldVariations->GetMaximum() );
    hYieldVariations            ->  Draw();
    uLatex                      ->  DrawLatexNDC(0.20,0.83,Form("#sigma_{Mean}: %.2f %s",   100*fabs(fErrorMean-fStandardYield)/(fStandardYield),"%"));
    uLatex                      ->  DrawLatexNDC(0.60,0.83,Form("#sigma_{Tot}: %.2f %s",    100*(fabs(fErrorMean-fStandardYield)+fErrorSTDV)/(fStandardYield),"%"));
    uLatex                      ->  DrawLatexNDC(0.23,0.77,Form("#sigma_{Dev}: %.2f %s",    100*(fErrorSTDV)/(fStandardYield),"%"));
    //
    cDrawResults                ->  cd(3);
    //
    //  Build the Variation histogram, fixed to 50 bins
    auto    hPTVariation        =   uBuildTH1F( fPTVariation, 50 );
    hPTVariation                ->  SetNameTitle( Form("hYieldVariations_%s",fName.Data()), "dN/dy variations" );
    //
    //  Gauss Fit
    hPTVariation                ->  Fit( "gaus", "IEMQS" );
    fErrorMean                  =   hPTVariation->GetFunction("gaus")->GetParameter(1);
    fErrorSTDV                  =   hPTVariation->GetFunction("gaus")->GetParameter(2);
    fErrorFULL                  =   fabs( fErrorMean - fStandardMeanPT ) + fErrorSTDV;
    //
    //  Save in Results
    fResult.push_back(fErrorFULL);
    fResult.push_back(fErrorFULL);
    //
    //  Draw the results
    hPTVariation                ->  SetMaximum( 1.4*hPTVariation->GetMaximum() );
    hPTVariation                ->  Draw();
    uLatex                      ->  DrawLatexNDC(0.20,0.83,Form("#sigma_{Mean}: %.2f %s",   100*fabs(fErrorMean-fStandardMeanPT)/(fStandardMeanPT),"%"));
    uLatex                      ->  DrawLatexNDC(0.60,0.83,Form("#sigma_{Tot}: %.2f %s",    100*(fabs(fErrorMean-fStandardMeanPT)+fErrorSTDV)/(fStandardMeanPT),"%"));
    uLatex                      ->  DrawLatexNDC(0.23,0.77,Form("#sigma_{Dev}: %.2f %s",    100*(fErrorSTDV)/(fStandardMeanPT),"%"));
    //
    cDrawResults             ->  cd(1);
    fTotalUncertaintySpectrum->Draw("same EP");
    //
    cDrawResults->Write();
    if ( fIsConditional )   cDrawResults->SaveAs(Form("%s/2D/ErrorFits_%s.pdf",fFolder.Data(),fName.Data()));
    else                    cDrawResults->SaveAs(Form("%s/1D/ErrorFits_%s.pdf",fFolder.Data(),fName.Data()));
    delete cDrawResults;
    //
    return fResult;
}
//
template < class Tclass >
std::vector<float>
fEvaluateExtrapolationUncertainty
 ( bool fIsConditional, Tclass* fStaticUncertaintySpectrum, Tclass* fVariableUncertaintySpectrum, std::vector<std::pair<TF1*,std::vector<float>>> fSystFunc, Double_t fMaximumFitRange = fMaxPT1D, Double_t fMinimumFitRange = fMinPT1D, TString fName = "", TString fFolder = "" )  {
    //
    //  Result container
    //  Format: [0] Stat Extrap Unc //! TODO: update
    std::vector<float> fResult;
    //
    //  Build the Total Uncertainty Spectrum
    Tclass*     fTotalUncertaintySpectrum   =   fSumErrors( fStaticUncertaintySpectrum, fVariableUncertaintySpectrum );
                fTotalUncertaintySpectrum   ->  SetName("TotalUncertaintySpectrum");
    //
    //  Measure the baselines
    fTotalUncertaintySpectrum       ->  Fit( fLevyTsallis, "Q",         "R", fMinimumFitRange,  fMaximumFitRange );
    fTotalUncertaintySpectrum       ->  Fit( fLevyTsallis, "QMES",      "R", fMinimumFitRange,  fMaximumFitRange );
    fTotalUncertaintySpectrum       ->  Fit( fLevyTsallis, "QMESI",     "R", fMinimumFitRange,  fMaximumFitRange );
    //
    auto    fStandardIntegral       =   fTotalUncertaintySpectrum->Integral("width");
    auto    fStandardExtrapolation  =   fLevyTsallis->Integral(0.,fMinPT1D);
    auto    fStandardMeanPT         =   fEvaluateMeanPT( fTotalUncertaintySpectrum );
    //
    //  Canvas to draw the Results checks
    SetStyle();
    auto    cDrawResults            =   new TCanvas("cDrawResult","cDrawResult",1800,600);
    gStyle                          ->  SetOptStat(0);
    cDrawResults                    ->  Divide(3,1);
    cDrawResults                    ->  cd(1);
    gPad                            ->  SetLogy();
    //
    //  Draw Spectra
    fTotalUncertaintySpectrum       ->  SetMarkerStyle(markers[2]);
    fTotalUncertaintySpectrum       ->  SetMarkerSize(2);
    fTotalUncertaintySpectrum       ->  SetMarkerColor(colors[0]);
    fTotalUncertaintySpectrum       ->  SetLineColor(kBlack);
    fTotalUncertaintySpectrum       ->  GetXaxis()  ->  SetRangeUser(0.,4.0);
    fTotalUncertaintySpectrum       ->  DrawClone("EP");
    //
    //  Generate utility histogram to evaluate statistical error
    std::vector<float> fYieldVariation;
    std::vector<float> fPTVariation;
    TLegend    *lFitFunctions   =   new TLegend(0.35,0.15,0.6,0.3);
    lFitFunctions   ->  SetLineColorAlpha(kWhite,0.);
    lFitFunctions   ->  SetFillColorAlpha(kWhite,0.);
    auto iTer = 0;
    auto jTer = 0;
    for ( auto iFuncRange : fSystFunc )  {
        //
        //  Prepping Fit Function
        fSetAllFunctions();
        auto    iFunc   =   iFuncRange.first;
        iFunc->SetLineColor( fGetRainbowColor(iTer,true) );
        //
        //  Prepping Associated Ranges
        for ( auto iRange : iFuncRange.second )  {
            fSetFunction(iFunc);
            //
            fTotalUncertaintySpectrum      ->  Fit(iFunc,          "EMRQS","", fMinPT1D,fMaxPT1D);
            fTotalUncertaintySpectrum      ->  Fit(iFunc,          "EMRQS","", fMinPT1D,fMaxPT1D);
            fTotalUncertaintySpectrum      ->  Fit(iFunc,          "EMRQS","", fMinimumFitRange,iRange);
            fTotalUncertaintySpectrum      ->  Fit(iFunc,          "EMRQSI","",fMinimumFitRange,iRange);
            //
            fYieldVariation.push_back   ( fStandardIntegral + iFunc  ->Integral(0.,fMinPT1D) );
            fPTVariation.push_back      ( fEvaluateMeanPT( fTotalUncertaintySpectrum, iFunc ) );
            //
            //  Draw Function
            iFunc->SetRange(0.,iRange);
            cDrawResults            ->  cd(1);
            iFunc->DrawCopy("same");
            jTer++;
        }
        //
        //  TLegend
        lFitFunctions    ->  AddEntry(iFunc,iFunc->GetName(),"L");
        iTer++;
    }
    //
    cDrawResults                ->  cd(2);
    //
    //  Build the Variation histogram, fixed to 50 bins
    auto    hYieldVariations    =   uBuildTH1F( fYieldVariation, 10 );
    hYieldVariations            ->  SetNameTitle( Form("hYieldVariations_%s",fName.Data()), "dN/dy variations" );
    //
    //  Gauss Fit
    auto fErrorMean             =   hYieldVariations->GetMean();
    auto fErrorSTDV             =   hYieldVariations->GetRMS();
    auto fErrorFULL             =   fabs( fErrorMean - ( fStandardExtrapolation + fStandardIntegral ) ) + fErrorSTDV;
    //
    //  Save in Results
    fResult.push_back(fErrorFULL);
    fResult.push_back(fErrorFULL);
    //
    //  Draw the results
    hYieldVariations            ->  SetMaximum( 1.4*hYieldVariations->GetMaximum() );
    hYieldVariations            ->  Draw();
    uLatex                      ->  DrawLatexNDC(0.20,0.83,Form("#sigma_{Mean}: %.2f %s",   100*fabs(fErrorMean-fStandardExtrapolation-fStandardIntegral)/(fStandardExtrapolation+fStandardIntegral),"%"));
    uLatex                      ->  DrawLatexNDC(0.60,0.83,Form("#sigma_{Tot}: %.2f %s",    100*(fabs(fErrorMean-fStandardExtrapolation-fStandardIntegral)+fErrorSTDV)/(fStandardExtrapolation+fStandardIntegral),"%"));
    uLatex                      ->  DrawLatexNDC(0.23,0.77,Form("#sigma_{Dev}: %.2f %s",    100*(fErrorSTDV)/(fStandardExtrapolation+fStandardIntegral),"%"));
    //
    cDrawResults                ->  cd(3);
    //
    //  Build the Variation histogram, fixed to 50 bins
    auto    hPTVariation        =   uBuildTH1F( fPTVariation, 10 );
    hPTVariation                ->  SetNameTitle( Form("hYieldVariations_%s",fName.Data()), "dN/dy variations" );
    //
    //  Gauss Fit
    fErrorMean                  =   hPTVariation->GetMean();
    fErrorSTDV                  =   hPTVariation->GetRMS();
    fErrorFULL                  =   fabs( fErrorMean - fStandardMeanPT ) + fErrorSTDV;
    //
    //  Save in Results
    fResult.push_back(fErrorFULL);
    fResult.push_back(fErrorFULL);
    //
    //  Draw the results
    hPTVariation                ->  SetMaximum( 1.4*hPTVariation->GetMaximum() );
    hPTVariation                ->  Draw();
    uLatex                      ->  DrawLatexNDC(0.20,0.83,Form("#sigma_{Mean}: %.2f %s",   100*fabs(fErrorMean-fStandardMeanPT)/(fStandardMeanPT),"%"));
    uLatex                      ->  DrawLatexNDC(0.60,0.83,Form("#sigma_{Tot}: %.2f %s",    100*(fabs(fErrorMean-fStandardMeanPT)+fErrorSTDV)/(fStandardMeanPT),"%"));
    uLatex                      ->  DrawLatexNDC(0.23,0.77,Form("#sigma_{Dev}: %.2f %s",    100*(fErrorSTDV)/(fStandardMeanPT),"%"));
    //
    cDrawResults                ->  cd(1);
    lFitFunctions               ->  Draw("same");
    fTotalUncertaintySpectrum   ->  Draw("same EP");
    //
    cDrawResults->Write();
    if ( fIsConditional )   cDrawResults->SaveAs(Form("%s/2D/ErrorFits_%s.pdf",fFolder.Data(),fName.Data()));
    else                    cDrawResults->SaveAs(Form("%s/1D/ErrorFits_%s.pdf",fFolder.Data(),fName.Data()));
    delete cDrawResults;
    //
    return fResult;
}
//
//_____________________________________________________________________________
//
template < class Tclass >
Double_t*
fExtrapolateModel
 ( bool fIsConditional, Tclass* gStatistics, Tclass* gSystematics, Double_t fIntegral = 0.032, TString fName = "ExtrapolateSignal", Double_t fMaximumFitRange = fMaxPT1D, Double_t fMinimumFitRange = fMinPT1D, TString fFolder = ""  )    {
    //  Optimisation mode
    gROOT->SetBatch(true);
    //
    //  Result format: Integral, Stat err low, Stat err high, Syst err low, syst err high, Mean pT, Stat err low, Stat err high, Syst err low, syst err high
    Double_t   *fResult     =   new Double_t    [13];
    //
    fLevyTsallis->Draw();
    //  Initialising the Fit Function
    fSetFunction(fLevyTsallis,fIntegral);
    fLevyTsallis->SetLineColor(kRed);
    fLevyTsallis->SetRange(0.,fMaxPT1D);
    //
    //  Setting the Fit Range
    Double_t    fMaxFitInt  =   TMath::Min( fMaximumFitRange, (double)fMaxPT1D );
    Double_t    fMinFitInt  =   TMath::Max( fMinimumFitRange, (double)fMinPT1D );
    //
    //  Generating a Total Error Spectra to fit and extrapolating at low pT
    Tclass      *gTotal      =   new Tclass(*(fSumErrors(gStatistics,gSystematics)));
    //
    //  Fitting a first time to evaluate integral in non-measured region
    gTotal                              ->  Fit(fLevyTsallis,"EMRQS","",fMinFitInt,fMaxFitInt);
    gTotal                              ->  Fit(fLevyTsallis,"EMRQSI","",fMinFitInt,fMaxFitInt);
    fResult[0]                          =   fLevyTsallis  ->Integral(0.,fMinPT1D);
    //
    fResult[5]                          =   fEvaluateMeanPT(gTotal,fLevyTsallis);
    //
    //  TCanvas w/ options
    TCanvas                *cDrawResult =   new TCanvas(Form("%s_%s",gStatistics->GetName(),fName.Data()),"",1600,800);
    cDrawResult->Divide(2,1);
    gStyle  ->SetOptStat(0);
    gTotal->SetMarkerStyle(markers[2]);
    gTotal->SetMarkerColor(colors[0]);
    gTotal->SetLineColor(kBlack);
    //
    //  Draw Spectra w/ function
    cDrawResult->cd(1);
    gPad    ->SetLogy();
    gTotal->Draw();
    fLevyTsallis->Draw("same");
    //
    //  Write info on the Extrapolation:
    TLatex  *fText   =   new TLatex();
    fText   ->  DrawLatexNDC(0.5,0.825,Form("1/N_{ev}dN/dy in [%.1f;%.1f]: %.7f",0.,fMinPT1D,fResult[0]));
    fText   ->  DrawLatexNDC(0.5,0.750,Form("Function: %s",fLevyTsallis->GetName()));
    fText   ->  DrawLatexNDC(0.5,0.700,Form("Fit range: [%.1f;%.1f]",fMinFitInt,fMaxFitInt));
    //
    cDrawResult->cd(2);
    gPad    ->SetLogy();
    auto gTotalCloseUp = (Tclass*)gTotal->Clone();
    gTotalCloseUp->GetXaxis()->SetRangeUser(0.01,2.5);
    gTotalCloseUp->Draw();
    //fLevyTsallis->Draw("same");
    //
    //  Save To File
    cDrawResult->Write();
    //
    //  Graphical Check
    if ( fIsConditional )   cDrawResult->   SaveAs(Form("%s/2D/EvaluationFit_%s.pdf",fFolder.Data(),fName.Data()));
    else                    cDrawResult->   SaveAs(Form("%s/1D/EvaluationFit_%s.pdf",fFolder.Data(),fName.Data()));
    delete cDrawResult;
    //
    uRandomGen->SetSeed(1);
    auto fStatResults                   =   fEvaluateExtrapolationUncertainty(fIsConditional,gSystematics,gStatistics,fIntegral,fMinFitInt,fMaxFitInt,TString("Stat_")+fName,fFolder);
    uRandomGen->SetSeed(1);
    auto fSystResults                   =   fEvaluateExtrapolationUncertainty(fIsConditional,gStatistics,gSystematics,fIntegral,fMinFitInt,fMaxFitInt,TString("Syst_")+fName,fFolder);
    //
    auto fSystResuFit                   =   fEvaluateExtrapolationUncertainty(fIsConditional,gStatistics,gSystematics,fSystFitFunctions,fMinFitInt,fMaxFitInt,TString("SFit_")+fName,fFolder);
    //
    fResult[1] = fResult[2] = fStatResults.at(0);
    fResult[3] = fResult[4] = SquareSum( { fSystResults.at(0), fSystResuFit.at(0) } );
    fResult[6] = fResult[7] = fStatResults.at(2);
    fResult[8] = fResult[9] = SquareSum( { fSystResults.at(2), fSystResuFit.at(2) } );
    //_____________________________________
    //
    //  End Optimisation mode
    gROOT->SetBatch(false);
    //
    return fResult;
}
Double_t*           fExtrapolateModel               ( std::vector<TH1F*>  gStatistics, std::vector<TH1F*>  gSystematics, Double_t fIntegral = 0.032, TString fName = "ExtrapolateSignal", Double_t fMaximumFitRange = fMaxPT1D, Double_t fMinimumFitRange = fMinPT1D, TString fFolder = ""  )    {
    //  Optimisation mode
    gROOT->SetBatch(true);
    //
    //  Result format: Integral, Stat err low, Stat err high, Syst err low, syst err high, Mean pT, Stat err low, Stat err high, Syst err low, syst err high
    Double_t   *fResult     =   new Double_t    [10];
    //
    for ( int i(0); i < 10; i++ ) fResult[i] = 0;
    //  End Optimisation mode
    gROOT->SetBatch(false);
    //
    return fResult;
}
//
//_____________________________________________________________________________
//
Double_t*           fIntegrateModel                 ( std::vector<TH1F*>  gStatistics, std::vector<TH1F*>  gSystematics, TString fName = "IntegrateSignal", TString fFolder = ""  )     {
    //  Optimisation mode
    gROOT->SetBatch(true);
    //
    //  Result format: Integral, Stat err low, Stat err high, Syst err low, syst err high, Mean pT, Stat err low, Stat err high, Syst err low, syst err high
    Double_t   *fResult     =   new Double_t    [10];
    //
    for ( int i(0); i < 10; i++ ) fResult[i] = 0;
    //
    //  End Optimisation mode
    gROOT->SetBatch(false);
    //
    return fResult;
}
//
//_____________________________________________________________________________
//
Double_t*           fMeasureFullYield               ( TH1F* gStatistics, TH1F* gSystematics, TString fName = "MeasureFullYield", TString fFolder = ""  )     {
    // Optimisation mode
    gROOT->SetBatch(true);
    //
    // Result format:  Integral, Stat err low, Stat err high, Syst err low, syst err high, Mean pT, Stat err low, Stat err high, Syst err low, syst err high
    Double_t   *fResult             =   new Double_t        [16];
    //
    bool fIsConditional = fName.Contains("2D");
    //
    auto        fIntegralStat       =   0.;
    auto        fIntegralSyst       =   0.;
    auto        fIntegral           =   gStatistics->IntegralAndError(-1,1000,fIntegralStat,"width");
                                        gSystematics->IntegralAndError(-1,1000,fIntegralSyst,"width");
    auto        fExtrapolResults    =   fExtrapolateModel   (fIsConditional,gStatistics,gSystematics,fIntegral,fName,fMinPT1D,fMaxPT1D,fFolder);
    Double_t    fIntegralMeanPTStat, fIntegralMeanPTSyst;
    fEvaluateMeanPT( gStatistics, fIntegralMeanPTStat );
    fEvaluateMeanPT( gSystematics, fIntegralMeanPTSyst );
    //
    //  Mean Value of Result
    fResult[0]  =   fIntegral +   fExtrapolResults[0];
    //
    //  Statistical Error of Result
    fResult[1]  =   SquareSum( { fIntegralStat, fExtrapolResults[1] } );
    fResult[2]  =   SquareSum( { fIntegralStat, fExtrapolResults[2] } );
    //
    //  Systematical Error of Result
    fResult[3]  =   SquareSum( { fIntegralSyst, fExtrapolResults[3] } );
    fResult[4]  =   SquareSum( { fIntegralSyst, fExtrapolResults[4] } );
    //
    //  Mean Value of pT
    fResult[5]  =   fExtrapolResults[5];
    //
    //  Statistical Error of Result
    fResult[6]  =   SquareSum( { fIntegralStat, fExtrapolResults[1] } );
    fResult[7]  =   SquareSum( { fIntegralStat, fExtrapolResults[2] } );
    //
    //  Systematical Error of Result
    fResult[8]  =   SquareSum( { fIntegralSyst, fExtrapolResults[3] } );
    fResult[9]  =   SquareSum( { fIntegralSyst, fExtrapolResults[4] } );
    //
    //  Extrapolation only
    fResult[10] =   fExtrapolResults[0];
    fResult[11] =   fExtrapolResults[1];
    fResult[12] =   fExtrapolResults[3];
    //
    //  Extrapolation only
    fResult[13] =   fIntegral;
    fResult[14] =   fIntegralStat;
    fResult[15] =   fIntegralSyst;
    //
    // End Optimisation mode
    gROOT->SetBatch(false);
    //
    return fResult;
}
Double_t*           fMeasureFullYield               ( std::vector<TH1F*>  gStatistics, std::vector<TH1F*>  gSystematics, TString fName = "MeasureFullYield", TString fFolder = "" )     {
    // Optimisation mode
    gROOT->SetBatch(true);
    //
    // Result format:  Integral, Stat err low, Stat err high, Syst err low, syst err high, Mean pT, Stat err low, Stat err high, Syst err low, syst err high
    Double_t   *fResult             =   new Double_t        [10];
    //
    auto        fIntegralResults    =   fIntegrateModel     (gStatistics,gSystematics,fName,fFolder);
    auto        fExtrapolResults    =   fExtrapolateModel   (gStatistics,gSystematics,fIntegralResults[0],fName,fMinPT1D,fMaxPT1D,fFolder);
    //
    //  Mean Value of Result
    fResult[0]  =   fIntegralResults[0] +   fExtrapolResults[0];
    //
    //  Statistical Error of Result
    fResult[1]  =   fIntegralResults[1] +   fExtrapolResults[1];
    fResult[2]  =   fIntegralResults[2] +   fExtrapolResults[2];
    //
    //  Systematical Error of Result
    fResult[3]  =   fIntegralResults[3] +   fExtrapolResults[4];
    fResult[4]  =   fIntegralResults[3] +   fExtrapolResults[4];
    //
    //  Mean Value of pT
    fResult[5]  =   fExtrapolResults[5]/fResult[0];
    //
    //  Statistical Error of Result
    fResult[6]  =   fResult[5]*( fExtrapolResults[6]/fExtrapolResults[5] + fResult[1]/fResult[0] );
    fResult[7]  =   fResult[5]*( fExtrapolResults[7]/fExtrapolResults[5] + fResult[2]/fResult[0] );
    //
    //  Systematical Error of Result
    fResult[8]  =   fResult[5]*( fExtrapolResults[8]/fExtrapolResults[5] + fResult[3]/fResult[0] );
    fResult[9]  =   fResult[5]*( fExtrapolResults[9]/fExtrapolResults[5] + fResult[4]/fResult[0] );
    //
    // End Optimisation mode
    gROOT->SetBatch(false);
    //
    return fResult;
}
//
//_____________________________________________________________________________
//
Double_t            fEvaluateINELgt0                ( Int_t iMultBin, TH1  *hMultCounter)  {
    Double_t    fResult =   0;
    auto    kUtilCount  =   (TH1F*)(hMultCounter->Clone());
    for ( Int_t iTer = 1; iTer <= kUtilCount->GetNbinsX(); iTer++ ) {
        auto    kBinCenter  =   hMultCounter->GetBinCenter(iTer);
        if ( kBinCenter < kMltTrgECls[0] || kBinCenter > kMltTrgECls[nMltTrgECls-1] ) continue;
        auto    kBinContent =   hMultCounter->GetBinContent(iTer);
        for ( Int_t iTe2 = 1; iTe2 < nMltTrgECls; iTe2++ )  {
            if ( kBinCenter < kMltTrgECls[iTe2] )   {
                kUtilCount->SetBinContent(iTer,kBinContent/kMultTrgEff[iTe2-1]);
                break;
            }
        }
    }
    if ( iMultBin < 0 || iMultBin > nBinMult )  {
        for ( Int_t iTer = 0; iTer < nBinMult; iTer++ ) fResult += fEvaluateINELgt0(iTer,kUtilCount);
    } else {
        fResult =   (kUtilCount->Integral(2+fArrMult[iMultBin],1+fArrMult[iMultBin+1]));
    }
    return      fResult;
}
//
//__________________________________________________________________________
//
void        fSetPhiCandidate                    ( TTree* TPhiCandidate, Struct_PhiEfficiency &evPhiEfficiency )    {
    TPhiCandidate-> SetBranchAddress    ("EventMask",       &evPhiEfficiency.EventMask);
    TPhiCandidate-> SetBranchAddress    ("TrueEventMask",   &evPhiEfficiency.TrueEventMask);
    TPhiCandidate-> SetBranchAddress    ("Multiplicity",    &evPhiEfficiency.Multiplicity);
    TPhiCandidate-> SetBranchAddress    ("nPhi",            &evPhiEfficiency.nPhi);
    TPhiCandidate-> SetBranchAddress    ("Px",              &evPhiEfficiency.Px);
    TPhiCandidate-> SetBranchAddress    ("Py",              &evPhiEfficiency.Py);
    TPhiCandidate-> SetBranchAddress    ("Pz",              &evPhiEfficiency.Pz);
    TPhiCandidate-> SetBranchAddress    ("Selection",       &evPhiEfficiency.Selection);
}
//
//_____________________________________________________________________________
//
void        fSetKaonCandidate                    ( TTree* TKaonCandidate, Struct_KaonEfficiency &evKaonEfficiency )    {
    TKaonCandidate->SetBranchAddress    ("EventMask",       &evKaonEfficiency.EventMask);
    TKaonCandidate->SetBranchAddress    ("TrueEventMask",   &evKaonEfficiency.TrueEventMask);
    TKaonCandidate->SetBranchAddress    ("Multiplicity",    &evKaonEfficiency.Multiplicity);
    TKaonCandidate->SetBranchAddress    ("nKaon",           &evKaonEfficiency.nKaon);
    TKaonCandidate->SetBranchAddress    ("Px",              &evKaonEfficiency.Px);
    TKaonCandidate->SetBranchAddress    ("Py",              &evKaonEfficiency.Py);
    TKaonCandidate->SetBranchAddress    ("Pz",              &evKaonEfficiency.Pz);
    TKaonCandidate->SetBranchAddress    ("Selection",       &evKaonEfficiency.Selection);
}
//
//_____________________________________________________________________________
//
Bool_t      fSetCandidates                      ( TTree* TPhiCnd, Struct_PhiEfficiency &fTargetPhi, TTree* TKaonCn, Struct_KaonEfficiency &fTargetKaon )    {
    if ( !TPhiCnd && !TKaonCn ) return false;
    if ( TPhiCnd )  fSetPhiCandidate(TPhiCnd,fTargetPhi);
    if ( TKaonCn )  fSetKaonCandidate(TKaonCn,fTargetKaon);
    return true;
}
//
//_____________________________________________________________________________
//
void        fSetPhiCandidate                    ( TTree* TPhiCandidate, Struct_PhiCandidate &evPhiCandidate )    {
    TPhiCandidate-> SetBranchAddress    ("EventMask",       &evPhiCandidate.EventMask);
    TPhiCandidate-> SetBranchAddress    ("Multiplicity",    &evPhiCandidate.Multiplicity);
    TPhiCandidate-> SetBranchAddress    ("nPhi",            &evPhiCandidate.nPhi);
    TPhiCandidate-> SetBranchAddress    ("Px",              &evPhiCandidate.Px);
    TPhiCandidate-> SetBranchAddress    ("Py",              &evPhiCandidate.Py);
    TPhiCandidate-> SetBranchAddress    ("Pz",              &evPhiCandidate.Pz);
    TPhiCandidate-> SetBranchAddress    ("InvMass",         &evPhiCandidate.InvMass);
    TPhiCandidate-> SetBranchAddress    ("iKaon",           &evPhiCandidate.iKaon);
    TPhiCandidate-> SetBranchAddress    ("jKaon",           &evPhiCandidate.jKaon);
    TPhiCandidate-> SetBranchAddress    ("TrueInvMass",     &evPhiCandidate.TrueInvMass);
}
//
//_____________________________________________________________________________
//
void        fSetKaonCandidate                    ( TTree* TKaonCandidate, Struct_KaonCandidate &evKaonCandidate )    {
    TKaonCandidate-> SetBranchAddress   ("EventMask",       &evKaonCandidate.EventMask);
    TKaonCandidate-> SetBranchAddress   ("Multiplicity",    &evKaonCandidate.Multiplicity);
    TKaonCandidate-> SetBranchAddress   ("nKaon",           &evKaonCandidate.nKaon);
    TKaonCandidate-> SetBranchAddress   ("Px",              &evKaonCandidate.Px);
    TKaonCandidate-> SetBranchAddress   ("Py",              &evKaonCandidate.Py);
    TKaonCandidate-> SetBranchAddress   ("Pz",              &evKaonCandidate.Pz);
    TKaonCandidate-> SetBranchAddress   ("Charge",          &evKaonCandidate.Charge);
    TKaonCandidate-> SetBranchAddress   ("TOFSigma",        &evKaonCandidate.SigmaTOF);
    TKaonCandidate-> SetBranchAddress   ("TPCSigma",        &evKaonCandidate.SigmaTPC);
}
//
//_____________________________________________________________________________
//
Bool_t      fSetCandidates                      ( TTree* TPhiCnd, Struct_PhiCandidate &fTargetPhi, TTree* TKaonCn, Struct_KaonCandidate &fTargetKaon )    {
    if ( !TPhiCnd && !TKaonCn ) return false;
    if ( TPhiCnd )  fSetPhiCandidate(TPhiCnd,fTargetPhi);
    if ( TKaonCn )  fSetKaonCandidate(TKaonCn,fTargetKaon);
    return true;
}
//
//_____________________________________________________________________________
//




Float_t             fMeasureMeanPT                  ( TF1 * fLowFit, TH1D * gTotal, Bool_t fReFit = false )   {
    Float_t fResult = 0.;
    if ( fReFit )   {
        gTotal->Fit(fLowFit);
    }
    std::vector<float>  fIntegral;
    std::vector<float>  fMeanPT;
    std::vector<float>  fWidth;
    auto    nEntries        =   gTotal->GetNbinsX();
    fIntegral.push_back(fLowFit->Integral(0.,gTotal->GetBinLowEdge(1)));
    fMeanPT.push_back(  fLowFit->Moment(1,0.,gTotal->GetBinLowEdge(1)));
    fWidth.push_back(   1.);
    for ( Int_t iBin = 1; iBin <= nEntries; iBin++ )   {
        fWidth.push_back(   gTotal->GetBinWidth(iBin) );
        fIntegral.push_back(gTotal->GetBinContent(iBin));
        fMeanPT.push_back(  fLowFit->Moment(1,gTotal->GetBinLowEdge(iBin),gTotal->GetBinLowEdge(iBin+1)));
    }
    fIntegral.push_back(fLowFit->Integral(gTotal->GetBinLowEdge(nEntries+1),1.e2));
    fMeanPT.push_back(  fLowFit->Moment(1,gTotal->GetBinLowEdge(nEntries+1),1.e2));
    fWidth.push_back(   1.);
    for ( Int_t iTer = 0; iTer < fIntegral.size(); iTer++ ) {
        fResult +=  fIntegral.at(iTer)*fMeanPT.at(iTer)*fWidth.at(iTer);
    }
    return fResult;
    /*
    
    //--------------------
    //  !TODO: Scorporare la funzione che calcola mean pT
    //  Evaluating the Mean pT
    if ( !fIsConditional )  {
        fResult[5]                          =   (fMinPT1D)*fFitFunc->Moment(1,0.,fMinPT1D)*(fResult[0]);
        for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )  {
            fResult[5]                     +=   (fArrPT1D[iPT1D+1] - fArrPT1D[iPT1D])*(fFitFunc->Moment(1,fArrPT1D[iPT1D],fArrPT1D[iPT1D+1]))*(gTotal->GetPointY(iPT1D));
        }
    }   else    {
        fResult[5]                          =   (fMinPT2D)*fFitFunc->Moment(1,0.,fMinPT2D)*(fResult[0]);
        for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
            fResult[5]                     +=   (fArrPT2D[iPT2D+1] - fArrPT2D[iPT2D])*(fFitFunc->Moment(1,fArrPT2D[iPT2D],fArrPT2D[iPT2D+1]))*(gTotal->GetPointY(iPT2D));
        }
    }
    //--------------------
     */
}
Float_t             fMeasureMeanPT                  ( TF1 * fLowFit, TH1F * gTotal, Bool_t fReFit = false )   {
    Float_t fResult = 0.;
    if ( fReFit )   {
        gTotal->Fit(fLowFit);
    }
    std::vector<float>  fIntegral;
    std::vector<float>  fMeanPT;
    std::vector<float>  fWidth;
    auto    nEntries        =   gTotal->GetNbinsX();
    fIntegral.push_back(fLowFit->Integral(0.,gTotal->GetBinLowEdge(1)));
    fMeanPT.push_back(  fLowFit->Moment(1,0.,gTotal->GetBinLowEdge(1)));
    fWidth.push_back(   1.);
    for ( Int_t iBin = 1; iBin <= nEntries; iBin++ )   {
        fWidth.push_back(   gTotal->GetBinWidth(iBin) );
        fIntegral.push_back(gTotal->GetBinContent(iBin));
        fMeanPT.push_back(  fLowFit->Moment(1,gTotal->GetBinLowEdge(iBin),gTotal->GetBinLowEdge(iBin+1)));
    }
    fIntegral.push_back(fLowFit->Integral(gTotal->GetBinLowEdge(nEntries+1),1.e2));
    fMeanPT.push_back(  fLowFit->Moment(1,gTotal->GetBinLowEdge(nEntries+1),1.e2));
    fWidth.push_back(   1.);
    for ( Int_t iTer = 0; iTer < fIntegral.size(); iTer++ ) {
        fResult +=  fIntegral.at(iTer)*fMeanPT.at(iTer)*fWidth.at(iTer);
    }
    return fResult;
    /*
    
    //--------------------
    //  !TODO: Scorporare la funzione che calcola mean pT
    //  Evaluating the Mean pT
    if ( !fIsConditional )  {
        fResult[5]                          =   (fMinPT1D)*fFitFunc->Moment(1,0.,fMinPT1D)*(fResult[0]);
        for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )  {
            fResult[5]                     +=   (fArrPT1D[iPT1D+1] - fArrPT1D[iPT1D])*(fFitFunc->Moment(1,fArrPT1D[iPT1D],fArrPT1D[iPT1D+1]))*(gTotal->GetPointY(iPT1D));
        }
    }   else    {
        fResult[5]                          =   (fMinPT2D)*fFitFunc->Moment(1,0.,fMinPT2D)*(fResult[0]);
        for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
            fResult[5]                     +=   (fArrPT2D[iPT2D+1] - fArrPT2D[iPT2D])*(fFitFunc->Moment(1,fArrPT2D[iPT2D],fArrPT2D[iPT2D+1]))*(gTotal->GetPointY(iPT2D));
        }
    }
    //--------------------
     */
}

Double_t            fIntegralOverHalf               ( std::vector<TH1F*>  hTarget, Double_t &fError = fDumpVar_Double_t  ) {
    Double_t    fResult =   0.;
    for ( Int_t iTer = 0; iTer < hTarget.size(); iTer++ )   {
        for ( Int_t jTer = iTer+1; jTer < hTarget.size(); jTer++ )   {
            fResult             +=  hTarget.at(iTer)->GetBinContent(jTer);
            fError              +=  hTarget.at(iTer)->GetBinError(jTer)*hTarget.at(iTer)->GetBinError(jTer);
        }
    }
    sqrt(fError);
    return fResult;
}
Double_t            fIntegralDiagonal               ( std::vector<TH1F*>  hTarget, Double_t &fError = fDumpVar_Double_t  ) {
    Double_t    fResult =   0.;
    for ( Int_t iTer = 0; iTer < hTarget.size(); iTer++ )   {
        fResult             +=  hTarget.at(iTer)->GetBinContent(iTer);
        fError              +=  hTarget.at(iTer)->GetBinError(iTer)*hTarget.at(iTer)->GetBinError(iTer);
    }
    sqrt(fError);
    return fResult;
}
Double_t            fIntegralDiagHalf               ( std::vector<TH1F*>  hTarget, Double_t &fError = fDumpVar_Double_t  ) {
    Double_t    fResult =   0.;
    fResult +=  fIntegralOverHalf( hTarget, fDumpVar_Double_t );
    fError  +=  fDumpVar_Double_t*fDumpVar_Double_t;
    fResult +=  fIntegralDiagonal( hTarget, fDumpVar_Double_t );
    fError  +=  fDumpVar_Double_t*fDumpVar_Double_t;
    sqrt(fError);
    return  fResult;
}
Double_t            fIntegralCorrHalf               ( std::vector<TH1F*>  hTarget, Double_t &fError = fDumpVar_Double_t  ) {
    Double_t    fResult =   0.;
    fResult +=  2*fIntegralOverHalf( hTarget, fDumpVar_Double_t );
    fError  +=  4*fDumpVar_Double_t*fDumpVar_Double_t;
    fResult +=  fIntegralDiagonal( hTarget, fDumpVar_Double_t );
    fError  +=  fDumpVar_Double_t*fDumpVar_Double_t;
    sqrt(fError);
    return  fResult;
}

std::pair<std::vector<float>,TH1F*>     uIntegralError                  ( std::vector<TH1F*>  hTargetStat, std::vector<TH1F*>  hTargetSyst )  {
    //
    
    // auto fError = ; 
    
    //
    //  Result Container
    std::pair<std::vector<float>,TH1F*> fResult;
    //
    //  Generate utility histogram to evaluate error
    std::vector<float>      uIntegralVar;
    std::vector<float>*     uIntegralVarCond = new std::vector<float> [hTargetStat.size()];
    std::vector<float>      uExtrapolVar;
    std::vector<float>      uExtrapolVa2;
    std::vector<float>*     uExtrapolVarCond = new std::vector<float> [hTargetStat.size()];
    //
    //  Evaluate Central Value
    TH1F   *hConditionalYieldStat   =   new TH1F("hConditionalYieldStat","hConditionalYieldStat",nBinPT2D,fArrPT2D);
    TH1F   *hConditionalYieldSyst   =   new TH1F("hConditionalYieldSyst","hConditionalYieldSyst",nBinPT2D,fArrPT2D);
    TH1F   *hConditionalYieldINT_   =   new TH1F("hConditionalYieldINT_","hConditionalYieldINT_",nBinPT2D,fArrPT2D);
    TH1F   *hConditionalYieldEXT_   =   new TH1F("hConditionalYieldEXT_","hConditionalYieldEXT_",nBinPT2D,fArrPT2D);
    TH1F   *hConditionalYieldTest   =   new TH1F("hConditionalYieldTest","hConditionalYieldTest",nBinPT2D,fArrPT2D);
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
        fSetAllFunctions();
        fSetFunction(fLevyTsallis);
        auto    fNewTarget  =   fSumErrors(hTargetStat.at(iPT2D),hTargetSyst.at(iPT2D));
        fNewTarget          ->  Fit(fLevyTsallis,   "EMRQS","",     fMinPT1D,fMaxPT1D);
        fNewTarget          ->  Fit(fLevyTsallis,   "EMRQSI","",    fMinPT1D,fMaxPT1D);
        //
        TCanvas*cTest = new TCanvas();
        fNewTarget->DrawCopy();
        fLevyTsallis->DrawCopy("same");
        cTest->SaveAs(Form("./result/Yield/ExtrapolateCheck/ttt_%i.pdf",iPT2D));
        delete cTest;
        //
        //  Integral
        auto    fIntegral   =   fNewTarget->Integral("width");
        //
        //  Extrapolation
        auto    fExtraplt   =   fLevyTsallis->Integral(0.,fMinPT1D);
        auto    fExtraplE   =   fLevyTsallis->IntegralError(0.,fMinPT1D);
        //
        //  Set Values
        hConditionalYieldStat->SetBinContent(iPT2D+1,fIntegral+fExtraplt);
        hConditionalYieldSyst->SetBinContent(iPT2D+1,fIntegral+fExtraplt);
        hConditionalYieldINT_->SetBinContent(iPT2D+1,fIntegral);
        hConditionalYieldEXT_->SetBinContent(iPT2D+1,fExtraplt);
        hConditionalYieldEXT_->SetBinError  (iPT2D+1,fExtraplE);
        //
        //  Set test
        hConditionalYieldTest->SetBinContent(iPT2D+1,fNewTarget->GetBinContent  (iPT2D+1));
        hConditionalYieldTest->SetBinError  (iPT2D+1,fNewTarget->GetBinError    (iPT2D+1));
    }
    //
    //  Save to File
    hConditionalYieldStat->Write();
    hConditionalYieldSyst->Write();
    hConditionalYieldINT_->Write();
    hConditionalYieldEXT_->Write();
    hConditionalYieldTest->Write();
    //
    //  Evaluate Variations
    auto cCheckMultiFit = new TCanvas();
    gStyle->SetOptStat(0);
    cCheckMultiFit->Divide(4,3);
    TH1F   *hConditionalUtil   =   new TH1F("hConditionalUtil","hConditionalUtil",nBinPT2D,fArrPT2D);
    for ( Int_t iTer = 0; iTer < kStatEvalCycles; iTer++ )  {
        //
        //  Randomise points Symmetrically
        auto fNewTarget = uRandomizePointsSymm(hTargetStat,hTargetSyst);
        //
        //  Loop over Randomised Conditional Spectra
        auto        jTer        =   0;
        Double_t    fIntegral   =   0.;
        Double_t    fExtraplt   =   0.;
        for ( auto hHistogram : fNewTarget )    {
            //
            //  Recover Integral of current Spectra & save it
            auto fCurrentIntegral = hHistogram->Integral("width");
            uIntegralVarCond[jTer].push_back(fCurrentIntegral);
            //
            //  Cumulative Integral
            fIntegral   +=  fCurrentIntegral*(fArrPT2D[jTer+1]-fArrPT2D[jTer]);
            //
            //  Prepping the Fit
            fSetAllFunctions();
            fSetFunction(fLevyTsallis);
            hHistogram->Fit(fLevyTsallis,"IQMRE0S");
            auto    fCurrentExtraplt    =   fLevyTsallis->Integral(0.,fMinPT1D);
            uExtrapolVarCond[jTer].push_back(fCurrentExtraplt);
            //
            //  Cumulative Extrapolation
            fExtraplt   +=  fCurrentExtraplt;
            //
            //  Make the Conditional Yield Histogram
            hConditionalUtil->SetBinContent (jTer+1, fCurrentIntegral+fCurrentExtraplt);
            //hConditionalUtil->SetBinError   (jTer+1, fError[jTer]);
            hConditionalYieldTest->SetBinContent( jTer+1, hHistogram->GetBinContent(jTer));
            //
            //  Graphics
            cCheckMultiFit->cd( jTer + 1 );
            gPad->SetLogy();
            if ( iTer == 0 )                    fSumErrors(hTargetStat.at(jTer),hTargetSyst.at(jTer))->Draw("SAME");
            fLevyTsallis->DrawCopy("same");
            if ( iTer == kStatEvalCycles-1 )    fSumErrors(hTargetStat.at(jTer),hTargetSyst.at(jTer))->Draw("SAME");
            //
            jTer++;
        }
        //
        uIntegralVar.push_back( fIntegral );
        //
        fSetAllFunctions();
        fSetFunction(fLevyTsallis);
        hConditionalUtil->Fit(fLevyTsallis,"IQMRE0S");
        cCheckMultiFit->cd( 11 );
        gPad->SetLogy();
        if ( iTer == 0 )                    hConditionalUtil->Draw("SAME");
        fLevyTsallis->DrawCopy("same");
        if ( iTer == kStatEvalCycles-1 )    hConditionalUtil->Draw("SAME");
        //
        uExtrapolVar.push_back( fLevyTsallis->Integral(0.,fMinPT1D) );
        //
        fSetAllFunctions();
        fSetFunction(fLevyTsallis);
        hConditionalYieldTest->Fit(fLevyTsallis,"IQMRE0S","",fMinPT1D,2.);
        cCheckMultiFit->cd( 12 );
        gPad->SetLogy();
        if ( iTer == 0 )                    hConditionalYieldTest->Draw("SAME");
        fLevyTsallis->DrawCopy("same");
        if ( iTer == kStatEvalCycles-1 )    hConditionalYieldTest->Draw("SAME");
        //
        uExtrapolVa2.push_back( fLevyTsallis->Integral(0.,fMinPT1D) );
    }
    cCheckMultiFit->SaveAs("./result/Yield/ExtrapolateCheck/cCheckMultiFit.pdf");
    delete cCheckMultiFit;
    //
    //  Analysis of Integral Variations
    auto cCheckMultiIntegral    = new TCanvas();
    gStyle->SetOptStat(0);
    cCheckMultiIntegral->Divide(4,3);
    //  First-Differential Spectrum
    cCheckMultiIntegral->cd(1);
    //  Generating the Utility Histogram
    auto    hIntegral       =   uBuildTH1F( uIntegralVar, 50, -hConditionalYieldINT_->Integral("width") );
    hIntegral               ->  SetNameTitle( "hIntegral", "dN/dy variations" );
    fSetAllFunctions();
    fSetFunction(fGauss);
    fGauss->SetParameter(0,hIntegral->Integral());
    fGauss->FixParameter(1,0);
    fGauss->SetParameter(2,hIntegral->GetRMS());
    hIntegral               ->  Fit(fGauss,"IMRE0S");
    
    hIntegral               ->  Draw();
    hIntegral               ->  Write();
    fGauss                  ->  DrawCopy("SAME");
    fGauss                  ->  Write();
    
    for ( Int_t iTer = 0; iTer < hTargetStat.size(); iTer++ )    {
        cCheckMultiIntegral->cd(iTer+2);
        
        auto    hIntegralCond       =   uBuildTH1F( uIntegralVarCond[iTer], 50, -hConditionalYieldINT_->GetBinContent(iTer+1) );
        hIntegralCond               ->  SetNameTitle( Form("hIntegral_%i",iTer), "dN/dy variations" );
        hIntegralCond               ->  Fit(fGauss,"IMRE0S");
        
        hIntegralCond               ->  Draw();
        hIntegralCond               ->  Write();
        fGauss                      ->  DrawCopy("SAME");
        fGauss                      ->  Write();
    }
    //
    cCheckMultiIntegral->SaveAs(Form("result/yield/Integral.pdf"));
    delete cCheckMultiIntegral;
    //
    //
    fSetAllFunctions();
    fSetFunction(fLevyTsallis);
    hConditionalYieldStat->Fit(fLevyTsallis,"IQMRE0S");
    //
    //  Analysis of Extrapolation Variations
    auto cCheckMultiExtrap    = new TCanvas();
    gStyle->SetOptStat(0);
    cCheckMultiExtrap->Divide(4,3);
    //  Second-Differential Spectrum
    cCheckMultiExtrap->cd(1);
    //  Generating the Utility Histogram
    auto    hExtrap         =   uBuildTH1F( uExtrapolVar, 50, -fLevyTsallis->Integral(0.,fMinPT1D) );
    hExtrap                 ->  SetNameTitle( "hExtrap", "dN/dy variations" );
    fSetAllFunctions();
    fSetFunction(fGauss);
    fGauss->SetParameter(0,hExtrap->Integral());
    fGauss->FixParameter(1,0);
    fGauss->SetParameter(2,hExtrap->GetRMS());
    hExtrap                 ->  Fit(fGauss,"IMRE0S");
    
    hExtrap                 ->  Draw();
    hExtrap                 ->  Write();
    fGauss                  ->  DrawCopy("SAME");
    fGauss                  ->  Write();
    
    cCheckMultiExtrap->cd(12);
    //  Generating the Utility Histogram
    auto    hExtra2         =   uBuildTH1F( uExtrapolVa2, 50, -fLevyTsallis->Integral(0.,fMinPT1D) );
    hExtra2                 ->  SetNameTitle( "hExtrap", "dN/dy variations" );
    fSetAllFunctions();
    fSetFunction(fGauss);
    fGauss->SetParameter(0,hExtra2->Integral());
    fGauss->FixParameter(1,0);
    fGauss->SetParameter(2,hExtra2->GetRMS());
    hExtrap                 ->  Fit(fGauss,"IMRE0S");
    
    hExtra2                 ->  Draw();
    hExtra2                 ->  Write();
    fGauss                  ->  DrawCopy("SAME");
    fGauss                  ->  Write();
    
    for ( Int_t iTer = 0; iTer < hTargetStat.size(); iTer++ )    {
        cCheckMultiExtrap->cd(iTer+2);
        
        auto    hExtrapCond         =   uBuildTH1F( uExtrapolVarCond[iTer], 50, -hConditionalYieldEXT_->GetBinContent(iTer+1) );
        hExtrapCond                 ->  SetNameTitle( Form("hExtrap_%i",iTer), "dN/dy variations" );
        hExtrapCond                 ->  Fit(fGauss,"IMRE0S");
        
        hExtrapCond                 ->  Draw();
        hExtrapCond                 ->  Write();
        fGauss                      ->  DrawCopy("SAME");
        fGauss                      ->  Write();
    }
    //
    cCheckMultiExtrap->SaveAs(Form("result/yield/Extrap.pdf"));
    delete cCheckMultiExtrap;
    
    
    
    /*
    // RESULT CHECK
    
    auto cDrawResult = new TCanvas();
    cDrawResult->Divide(4,3);

    cDrawResult->cd(1);
    auto    hStatIntegral   =   uBuildTH1F( uIntegralVar, 50, -hConditionalYieldINT_->Integral("width") );
    hStatIntegral           ->  SetNameTitle( "  ssdd", "dN/dy variations" );
    auto fErr = 0.;
    auto    jTer    =   0;
    for ( auto hIntegrand : hTargetStat )    {
        auto gfs = 0.;
        hIntegrand->IntegralAndError(-1,10000,gfs,"width");
        fErr += gfs*gfs;
    }
    hStatIntegral->Fit("gaus","IMREQ0S");
    auto fResultsss = hStatIntegral->GetFunction("gaus");
    auto siggg = fResultsss->GetParameter(2);
    hStatIntegral->Draw();
    fResultsss->Draw("same");
    //uLatex->DrawLatexNDC(0.60,0.80,Form("OLD : %.4f",1.e6*fIntegralStat));
    //uLatex->DrawLatexNDC(0.60,0.80,Form("OLD : %.4f",1.e6*sqrt(fErr)));
    uLatex->DrawLatexNDC(0.60,0.75,Form("GSS : %.4f",1.e6*siggg));
    uLatex->DrawLatexNDC(0.60,0.70,Form("RMS : %.4f",1.e6*hStatIntegral->GetRMS()));
    
    hStatIntegral->SetName("NEW_-1");
    hStatIntegral->Write();
    for ( Int_t iTer = 0; iTer < hTargetStat.size(); iTer++ )    {
        cDrawResult->cd(iTer+2);
        uCleanOutsiders(  uIntegralVarCond[iTer] );
        auto    hConditional   =   uBuildTH1F( uIntegralVarCond[iTer], 50, -hConditionalYieldINT_->GetBinContent(iTer+1) );
        hConditional           ->  SetNameTitle( Form("NEW_%i",iTer), "dN/dy variations" );
        hConditional->Fit("gaus","IMREQ0S");
        hConditional->Draw();
        fResultsss = hConditional->GetFunction("gaus");
        siggg = fResultsss->GetParameter(2);
        fResultsss->Draw("same");
        //uLatex->DrawLatexNDC(0.60,0.80,Form("OLD : %.4f",1.e6*hRES_2D_Cond2_Stat_INT->GetBinError(iTer+1)));
        uLatex->DrawLatexNDC(0.60,0.75,Form("GSS : %.4f",1.e6*siggg));
        uLatex->DrawLatexNDC(0.60,0.70,Form("RMS : %.4f",1.e6*hConditional->GetRMS()));
        hConditional->Write();
    }
    
    cDrawResult->SaveAs(Form("result/yield/Integral.pdf"));
    delete cDrawResult;
    
    cDrawResult = new TCanvas();
    cDrawResult->Divide(4,3);
    cDrawResult->cd(1);
    
    fSetAllFunctions();
    fSetFunction(fLevyTsallis);
    hConditionalYieldStat->Fit(fLevyTsallis,"IQMRE0S");
    
    auto hStatExtrapol   =   uBuildTH1F( uExtrapolVar, 50, -fLevyTsallis->Integral(0.,fMinPT1D) );
    hStatExtrapol           ->  SetNameTitle( "  hStatExtrapol", "dN/dy variations" );
    hStatExtrapol->Fit("gaus","IMRQE0S");
    hStatExtrapol->Draw();
    auto fResulFits = hStatExtrapol->GetFunction("gaus");
    siggg = fResulFits->GetParameter(2);
    fResulFits->DrawCopy("same");
    //uLatex->DrawLatexNDC(0.60,0.80,Form("OLD : %.4f",1.e6*fExtrapolResults[1]));
    uLatex->DrawLatexNDC(0.60,0.75,Form("GSS : %.4f",1.e6*siggg));
    uLatex->DrawLatexNDC(0.60,0.70,Form("RMS : %.4f",1.e6*hStatIntegral->GetRMS()));
    
    hStatIntegral->SetName("ENEW_-1");
    hStatIntegral->Write();
    for ( Int_t iTer = 0; iTer < hTargetStat.size(); iTer++ )    {
        cDrawResult->cd(iTer+2);
        uCleanOutsiders(  uExtrapolVarCond[iTer] );
        auto    hConditional   =   uBuildTH1F( uExtrapolVarCond[iTer], 50, -hConditionalYieldEXT_->GetBinContent(iTer+1) );
        hConditional           ->  SetNameTitle( Form("ENEW_%i",iTer), "dN/dy variations" );
        hConditional->Fit("gaus","IMRE0QS");
        hConditional->Draw();
        fResultsss = hConditional->GetFunction("gaus");
        siggg = fResultsss->GetParameter(2);
        fResultsss->Draw("same");
        //uLatex->DrawLatexNDC(0.60,0.80,Form("OLD : %.4f",1.e6*hRES_2D_Cond2_Stat_EXT->GetBinError(iTer+1)));
        uLatex->DrawLatexNDC(0.60,0.75,Form("GSS : %.4f",1.e6*siggg));
        uLatex->DrawLatexNDC(0.60,0.70,Form("RMS : %.4f",1.e6*hConditional->GetRMS()));
        
        hConditional->Write();
    }
    
    cDrawResult->SaveAs(Form("result/yield/Extrapol.pdf"));
    //delete cDrawResult;
    */
    
    return fResult;
}

void
fCalculateSystematics
(){
    
}

#endif
