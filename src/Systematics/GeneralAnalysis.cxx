// File for 1-Dimensional Analysis:
// !TODO: Set Global Counter to avoid memory loss due to overlap in names of th**
// !TODO: Make copies for TH1F/TH1D etc.
#include "../../inc/AliAnalysisPhiPair.h"
#include "RooMsgService.h"
//
void
GeneralAnalysis
( TH1F* hStandard, std::vector<TH1F*> hVariations, TString fFolder, Bool_t kNoBarlowCheck = false  ) {
    //
    gErrorIgnoreLevel   =   kWarning;
    //
    // --------- FIND THE SYSTEMATICAL RELEVANT VARIATIONS
    //
    SetStyle();
    gROOT   ->  ProcessLine ( Form(".! mkdir -p %s",(fFolder+TString("/plots/BarlowCheck/1D/")).Data()) );
    gROOT   ->  ProcessLine ( Form(".! mkdir -p %s",(fFolder+TString("/plots/BinByBinCheck/1D/")).Data()) );
    auto fRelevantVariations    =  uIsRelevantVariation(hStandard,hVariations,(fFolder+TString("/plots/BarlowCheck/1D/")).Data(),"1D");
    //
    if ( kNoBarlowCheck )   {
        for ( auto fRelvenant : fRelevantVariations ) fRelvenant = true;
    }
    // --------- USE THE SYSTEMATICAL RELEVANT VARIATIONS TO DETERMINE THE UNCERTAINTY
    //
    //for ( [...] ) Loop to exclude non Barlow variations
    //
    auto    uStackSystematic    =   uBuildSystematicStack(hStandard,hVariations,fRelevantVariations);
    //
    gROOT->SetBatch();
    TLegend    *lLegend =   new TLegend(0.2,0.85,0.4,0.7);
    TCanvas*c1 = new TCanvas();
    gStyle->SetOptStat(0);
    gPad->SetLogx();
    uStackSystematic->SetMinimum(0);
    uStackSystematic->SetMaximum(uStackSystematic->GetMaximum()*1.3);
    uStackSystematic->Draw("");
    uStackSystematic->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    uStackSystematic->GetYaxis()->SetTitleOffset(1.5);
    uStackSystematic->GetYaxis()->SetTitle("Systematic Uncertainty (%)");
    lLegend->AddEntry((TH1F*)uStackSystematic->GetHists()->At(0),"#mu contr.","F");
    lLegend->AddEntry((TH1F*)uStackSystematic->GetHists()->At(1),"#sigma contr.","F");
    lLegend->Draw();
    c1->SaveAs((fFolder+TString("/plots/Full_1D_Sys.pdf")).Data());
    gROOT->SetBatch(kFALSE);
    //
    TFile      *fOutput =   new TFile   (Form("%s%s",fFolder.Data(),"/1D_Systematic.root"),"recreate");
    uBuildSystematicError(hStandard,hVariations,(fFolder+TString("/plots/BinByBinCheck/1D/")).Data(),"1D",fRelevantVariations)->Write();
    fOutput->Close();
    //
    gErrorIgnoreLevel   =   kInfo;
    //
}
//
void
GeneralAnalysis
( TH2F* hStandard, std::vector<TH2F*> hVariations, TString fFolder, Bool_t kNoBarlowCheck = false ) {
    //
    gErrorIgnoreLevel   =   kWarning;
    //
    // --------- FIND THE SYSTEMATICAL RELEVANT VARIATIONS
    //
    SetStyle();
    gROOT   ->  ProcessLine ( Form(".! mkdir -p %s",(fFolder+TString("/plots/BarlowCheck/2D/")).Data()) );
    gROOT   ->  ProcessLine ( Form(".! mkdir -p %s",(fFolder+TString("/plots/BinByBinCheck/2D/")).Data()) );
    auto fRelevantVariations    =  uIsRelevantVariation(hStandard,hVariations,(fFolder+TString("/plots/BarlowCheck/2D/")).Data(),"2D",true);
    //
    if ( kNoBarlowCheck )   {
        for ( auto fRelvenant : fRelevantVariations ) fRelvenant = true;
    }
    //
    // --------- USE THE SYSTEMATICAL RELEVANT VARIATIONS TO DETERMINE THE UNCERTAINTY
    //
    //for ( [...] ) Loop to exclude non Barlow variations
    //
    auto    uStackSystematic    =   uBuildSystematicStack(hStandard,hVariations,fRelevantVariations);
    //
    gROOT->SetBatch();
    TLegend    *lLegend =   new TLegend(0.2,0.85,0.4,0.7);
    lLegend->AddEntry((TH1F*)uStackSystematic.at(0)->GetHists()->At(0),"#mu contr.","F");
    lLegend->AddEntry((TH1F*)uStackSystematic.at(0)->GetHists()->At(1),"#sigma contr.","F");
    TCanvas*c1 = new TCanvas();
    c1->Divide(uStackSystematic.size());
    auto i = 0;
    for ( auto fStack : uStackSystematic )  {
        c1->cd(i);
        gStyle->SetOptStat(0);
        gPad->SetLogx();
        uStackSystematic.at(i)->SetMinimum(0);
        uStackSystematic.at(i)->SetMaximum(uStackSystematic.at(i)->GetMaximum()*1.4);
        uStackSystematic.at(i)->Draw("");
        uStackSystematic.at(i)->GetXaxis()->SetTitle("#it{p}_{T,#phi_{1}} (GeV/#it{c})");
        uStackSystematic.at(i)->GetYaxis()->SetTitleOffset(1.5);
        uStackSystematic.at(i)->GetYaxis()->SetTitle("Systematic Uncertainty (%)");
        uLatex->DrawLatexNDC(0.5,0.82,Form("#it{p}_{T,#phi_{2}} (GeV/#it{c}) #in [%.1f;%.1f]",fArrPT2D[i],fArrPT2D[i+1]));
        lLegend->Draw();
        c1->SaveAs((fFolder+TString("/plots/Full_2D_Sys_")+TString(Form("%i.pdf",i))).Data());
        gROOT->SetBatch(kFALSE);
        i++;
    }
    //
    TFile      *fOutput =   new TFile   (Form("%s%s",fFolder.Data(),"/2D_Systematic.root"),"recreate");
    uBuildSystematicError(hStandard,hVariations,(fFolder+TString("/plots/BinByBinCheck/2D/")).Data(),"2D",fRelevantVariations,true)->Write();
    fOutput->Close();
    //
    gErrorIgnoreLevel   =   kInfo;
    //
}
//
void
GeneralAnalysis
( TH1F* hStandard, std::vector<TH1F*> hVariations, TH2F* h2Standard, std::vector<TH2F*> h2Variations, TString fFolder = "" ) {
    uEvaluateRatioError(hStandard,hVariations,h2Standard,h2Variations,fFolder);
    //
    //! TODO: (1) This is a temporary fix, please address it AYEC.
    //! **** (1)
    TFile*  fInputData  =   new TFile   ( Form(kASigExtp_FitCheckRst,"Yield/Systematics/Standard/") );
    TH1F*   fUtilityC1  =   (TH1F*)(fInputData->Get("hRES_2D_Cond2_Stat"));
    TH1F*   fUtilityC2  =   (TH1F*)(fInputData->Get("hRES_2D_Cond2_Syst"));
    //! **** (1)
    //
    gROOT->SetBatch();
    //
    gErrorIgnoreLevel   =   kFatal;
    //
    fSetAllBins();
    //
    gROOT    ->  ProcessLine(Form(".! mkdir -p %s/plots/1D/",fFolder.Data()));
    gROOT    ->  ProcessLine(Form(".! mkdir -p %s/plots/2D/",fFolder.Data()));
    //
    fStartTimer("Systematic uncertainties determination: 1D");
    auto    k1DStandardExtrap   =   fExtrapolateModel(false,hStandard,fSetSystErrors(hStandard),0.032,"SystEval",fMinPT1D,fMaxPT1D,fFolder+TString("/plots/"));
    auto    k1DStandardIntegr   =   hStandard->Integral("width");
    auto    k1DFullYieldStand   =   k1DStandardIntegr + k1DStandardExtrap[0];
    //
    std::vector<Float_t>   k1DVariationExtrap;
    std::vector<Float_t>   k1DVariationIntegr;
    std::vector<Float_t>   k1DVariationFull__;
    std::vector<Float_t>   k1DVariatinExtVanl;
    std::vector<Float_t>   k1DVariatinIntVanl;
    int iTer = 0;
    for ( auto  h1DVarHisto : hVariations ) {
        fPrintLoopTimer("Systematic uncertainties determination: 1D",iTer+1,hVariations.size()+1,1);
        auto    k1DVariatinExtrap   =   fExtrapolateModel(false,h1DVarHisto,fSetSystErrors(h1DVarHisto),0.032,"SystEval",fMinPT1D,fMaxPT1D,fFolder+TString("/plots/"));
        auto    k1DVariatinIntegr   =   h1DVarHisto->Integral("width");
        k1DVariatinIntVanl.push_back( k1DVariatinIntegr );
        k1DVariatinExtVanl.push_back( k1DVariatinExtrap[0] );
        k1DVariationIntegr.push_back( 1 - (k1DVariatinIntegr/k1DStandardIntegr) );
        k1DVariationExtrap.push_back( 1 - (k1DVariatinExtrap[0]/k1DStandardExtrap[0]) );
        k1DVariationFull__.push_back( 1 - ( (k1DVariatinIntegr+k1DVariatinExtrap[0]) / (k1DStandardIntegr+k1DStandardExtrap[0]) ) );
        iTer++;
    }
    //
    //      1D Integral Error
    TH1F       *h1DIntegralError    =   uBuildTH1F(k1DVariationIntegr,2000,0,-0.5,0.5);
    auto        f1DIntegralError    =   0.;
    f1DIntegralError               +=   fabs(h1DIntegralError->GetMean());
    f1DIntegralError               +=   h1DIntegralError->GetRMS();
    //
    //      1D Extrapolation Error
    TH1F       *h1DExtrapolError    =   uBuildTH1F(k1DVariationExtrap,2000,0,-0.5,0.5);
    auto        f1DExtrapolError    =   0.;
    f1DExtrapolError               +=   fabs(h1DExtrapolError->GetMean());
    f1DExtrapolError               +=   h1DExtrapolError->GetRMS();
    //
    //      1D Calculated ( INT + EXT ) Error
    TH1F       *h1DCalculatError    =   uBuildTH1F(k1DVariationFull__,2000,0,-0.5,0.5);
    auto        f1DCalculatError    =   0.;
    f1DCalculatError               +=   fabs(h1DCalculatError->GetMean());
    f1DCalculatError               +=   h1DCalculatError->GetRMS();
    //
    fStopTimer("Systematic uncertainties determination: 1D");
    fStartTimer("Systematic uncertainties determination: 2D");
    //
    hName   =   Form("h2DStandard_Extrap");
    hTitle  =   Form("h2DStandard_Extrap");
    TH1F       *h2DStandard_Extrap_Stat =   new TH1F(hName,hTitle,nBinPT2D,fArrPT2D);
    TH1F       *h2DStandard_Extrap_Syst =   new TH1F(hName,hTitle,nBinPT2D,fArrPT2D);
    //
    for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ ) {
        //
        auto hTarget    =   h2Standard->ProjectionX(Form("Proj_STD_%i",iFit+1),iFit+1,iFit+1);
        //
        auto fResults   =   fExtrapolateModel(true,hTarget,fSetSystErrors(hTarget,iFit),0.0005,Form("SystEval2D_%i",iFit),fMinPT2D,fMaxPT2D,fFolder+TString("/plots/"));
        //
        h2DStandard_Extrap_Stat  ->  SetBinContent   ( iFit+1, fResults[0] );
        h2DStandard_Extrap_Syst  ->  SetBinContent   ( iFit+1, fResults[0] );
        //
        //! **** (1)
        h2DStandard_Extrap_Stat  ->  SetBinError     ( iFit+1, fUtilityC1->GetBinError(iFit+1) );
        h2DStandard_Extrap_Syst  ->  SetBinError     ( iFit+1, fUtilityC2->GetBinError(iFit+1) );
        //! **** (1)
    }
    //
    auto    k2DStandardExtrap_Ut    =   fExtrapolateModel(true,h2DStandard_Extrap_Stat,h2DStandard_Extrap_Syst,0.0005,Form("SystEval2D_2D"),fMinPT2D,fMaxPT2D,fFolder+TString("/plots/"));
    auto    k2DStandardExtrap       =   k2DStandardExtrap_Ut[0] + 2*h2DStandard_Extrap_Stat->Integral("width");
    auto    k2DStandardIntegr       =   h2Standard->Integral("width");
    auto    k2DFullYieldStand       =   k2DStandardIntegr + k2DStandardExtrap;
    //
    std::vector<Float_t>   k2DVariationExtrap;
    std::vector<Float_t>   k2DVariationIntegr;
    std::vector<Float_t>   k2DVariationFull__;
    std::vector<Float_t>   k2DVariationFullR1;
    std::vector<Float_t>   k2DVariationFullR2;
    std::vector<Float_t>   k2DVariationFullP1;
    std::vector<Float_t>   k2DVariationFullP2;
    iTer = 0;
    for ( auto  h2DVarHisto : h2Variations ) {
        fPrintLoopTimer("Systematic uncertainties determination: 2D",iTer+1,h2Variations.size()+1,1);
        if ( iTer+1 > hVariations.size() ) continue;
        for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ ) {
            //
            auto hTarget    =   h2DVarHisto->ProjectionX(Form("Proj_STD_%i",iFit+1),iFit+1,iFit+1);
            //
            auto fResults   =   fExtrapolateModel(true,hTarget,fSetSystErrors(hTarget,iFit),0.0005,Form("SystEval2D_%i",iFit),fMinPT2D,fMaxPT2D,fFolder+TString("/plots/"));
            //
            h2DStandard_Extrap_Stat  ->  SetBinContent   ( iFit+1, fResults[0] );
            h2DStandard_Extrap_Syst  ->  SetBinContent   ( iFit+1, fResults[0] );
            //! **** (1)
            h2DStandard_Extrap_Stat  ->  SetBinError     ( iFit+1, fUtilityC1->GetBinError(iFit+1) );
            h2DStandard_Extrap_Syst  ->  SetBinError     ( iFit+1, fUtilityC2->GetBinError(iFit+1) );
            //! **** (1)
        }
        auto    k2DVariatinExtrap_Ut    =   fExtrapolateModel(true,h2DStandard_Extrap_Stat,h2DStandard_Extrap_Syst,0.0005,Form("SystEval2D_2D"),fMinPT2D,fMaxPT2D,fFolder+TString("/plots/"));
        auto    k2DVariatinExtrap       =   k2DVariatinExtrap_Ut[0] + 2*h2DStandard_Extrap_Stat->Integral("width");
        auto    k2DVariatinIntegr       =   h2DVarHisto->Integral("width");
        k2DVariationIntegr.push_back( 1 - (k2DVariatinIntegr/k2DStandardIntegr) );
        k2DVariationExtrap.push_back( 1 - (k2DVariatinExtrap/k2DStandardExtrap) );
        k2DVariationFull__.push_back( 1 - ( (k2DVariatinIntegr+k2DVariatinExtrap) / (k2DStandardIntegr+k2DStandardExtrap) ) );
        auto    k2DFullYieldVarit       =   k2DVariatinIntegr + k2DVariatinExtrap;
        auto    k1DVariatinExtrap       =   k1DVariatinExtVanl.at(iTer);
        auto    k1DVariatinIntegr       =   k1DVariatinIntVanl.at(iTer);
        auto    k1DFullYieldVarit       =   k1DVariatinIntegr + k1DVariatinExtrap;
        //
        //  Ratios and composite quantitites;
        auto    kR1 =   (k2DFullYieldVarit*k1DFullYieldStand)/(k1DFullYieldVarit*k2DFullYieldStand);
        auto    kR2 =   (k2DFullYieldVarit*k1DFullYieldStand*k1DFullYieldStand)/(k2DFullYieldStand*k1DFullYieldVarit*k1DFullYieldVarit);
        auto    kP1 =   (fSigmaPhiValue(k1DFullYieldVarit,k2DFullYieldVarit))/(fSigmaPhiValue(k1DFullYieldStand,k2DFullYieldStand));
        auto    kP2 =   (fGammaPhiValue(k1DFullYieldVarit,k2DFullYieldVarit))/(fGammaPhiValue(k1DFullYieldStand,k2DFullYieldStand));
        k2DVariationFullR1.push_back( 1 - kR1 );
        k2DVariationFullR2.push_back( 1 - kR2 );
        k2DVariationFullP1.push_back( 1 - kP1 );
        k2DVariationFullP2.push_back( 1 - kP2 );
        iTer++;
    }
    //
    //      2D Integral Error
    TH1F       *h2DIntegralError    =   uBuildTH1F(k2DVariationIntegr,2000,0,-0.5,0.5);
    auto        f2DIntegralError    =   0.;
    f2DIntegralError               +=   fabs(h2DIntegralError->GetMean());
    f2DIntegralError               +=   h2DIntegralError->GetRMS();
    //
    //      2D Extrapolation Error
    TH1F       *h2DExtrapolError    =   uBuildTH1F(k2DVariationExtrap,2000,0,-0.5,0.5);
    auto        f2DExtrapolError    =   0.;
    f2DExtrapolError               +=   fabs(h2DExtrapolError->GetMean());
    f2DExtrapolError               +=   h2DExtrapolError->GetRMS();
    //
    //      2D Calculated Error
    TH1F       *h2DCalculatError    =   uBuildTH1F(k2DVariationFull__,2000,0,-0.5,0.5);
    auto        f2DCalculatError    =   0.;
    f2DCalculatError               +=   fabs(h2DCalculatError->GetMean());
    f2DCalculatError               +=   h2DCalculatError->GetRMS();
    //
    //      Ratio 1 Calculated Error
    TH1F       *hR1CalculatError    =   uBuildTH1F(k2DVariationFullR1,2000,0,-0.5,0.5);
    auto        fR1CalculatError    =   0.;
    fR1CalculatError               +=   fabs(hR1CalculatError->GetMean());
    fR1CalculatError               +=   hR1CalculatError->GetRMS();
    //
    //      Ratio 2 Calculated Error
    TH1F       *hR2CalculatError    =   uBuildTH1F(k2DVariationFullR2,2000,0,-0.5,0.5);
    auto        fR2CalculatError    =   0.;
    fR2CalculatError               +=   fabs(hR2CalculatError->GetMean());
    fR2CalculatError               +=   hR2CalculatError->GetRMS();
    //
    //      Parameter 1 Calculated Error
    TH1F       *hP1CalculatError    =   uBuildTH1F(k2DVariationFullP1,2000,0,-0.5,0.5);
    auto        fP1CalculatError    =   0.;
    fP1CalculatError               +=   fabs(hP1CalculatError->GetMean());
    fP1CalculatError               +=   hP1CalculatError->GetRMS();
    //
    //      Parameter 2 Calculated Error
    TH1F       *hP2CalculatError    =   uBuildTH1F(k2DVariationFullP2,2000,0,-0.5,0.5);
    auto        fP2CalculatError    =   0.;
    fP2CalculatError               +=   fabs(hP2CalculatError->GetMean());
    fP2CalculatError               +=   hP2CalculatError->GetRMS();
    //
    fStopTimer("Systematic uncertainties determination: 2D");
    //
    TH1F   *hCalculateFull_Sng   =   new TH1F("hCalculateFull_Sng",   "", 6,  0,  6);
    hCalculateFull_Sng->GetYaxis()->SetTitle("Systematic uncertainty (%)");
    hCalculateFull_Sng->GetXaxis()->SetNdivisions(6);
    hCalculateFull_Sng->GetXaxis()->SetBinLabel(hCalculateFull_Sng->GetXaxis()->FindBin(0.5),"#LT Y_{#phi} #GT");
    hCalculateFull_Sng->GetXaxis()->SetBinLabel(hCalculateFull_Sng->GetXaxis()->FindBin(1.5),"#LT Y_{#phi#phi} #GT");
    hCalculateFull_Sng->GetXaxis()->SetBinLabel(hCalculateFull_Sng->GetXaxis()->FindBin(2.5),"#frac{ #LT Y_{#phi#phi} #GT }{ #LT Y_{#phi} #GT  }");
    hCalculateFull_Sng->GetXaxis()->SetBinLabel(hCalculateFull_Sng->GetXaxis()->FindBin(3.5),"#frac{ #LT Y_{#phi#phi} #GT }{ #LT Y_{#phi} #GT^{2}  }");
    hCalculateFull_Sng->GetXaxis()->SetBinLabel(hCalculateFull_Sng->GetXaxis()->FindBin(4.5),"#sigma^{2}_{#phi}");
    hCalculateFull_Sng->GetXaxis()->SetBinLabel(hCalculateFull_Sng->GetXaxis()->FindBin(5.5),"#gamma_{#phi}");
    hCalculateFull_Sng->GetXaxis()->LabelsOption("h");
    hCalculateFull_Sng->GetXaxis()->SetLabelSize(0.04);
    hCalculateFull_Sng->SetMarkerColor(colors[4]);
    hCalculateFull_Sng->SetLineWidth(3);
    hCalculateFull_Sng->SetMarkerStyle(markers[5]);
    hCalculateFull_Sng->SetBinContent(1,f1DCalculatError);
    hCalculateFull_Sng->SetBinContent(2,f2DCalculatError);
    hCalculateFull_Sng->SetBinContent(3,fR1CalculatError);
    hCalculateFull_Sng->SetBinContent(4,fR2CalculatError);
    hCalculateFull_Sng->SetBinContent(5,fP1CalculatError);
    hCalculateFull_Sng->SetBinContent(6,fP2CalculatError);
    hCalculateFull_Sng->SetBinError  (1,0);
    hCalculateFull_Sng->SetBinError  (2,0);
    hCalculateFull_Sng->SetBinError  (3,0);
    hCalculateFull_Sng->SetBinError  (4,0);
    hCalculateFull_Sng->SetBinError  (5,0);
    hCalculateFull_Sng->SetBinError  (6,0);
    hCalculateFull_Sng->Scale(100);
    //
    TH1F   *hIntegral_Sng   =   new TH1F("hIntegral_Sng",   "", 6,  0,  6);
    hIntegral_Sng->SetMarkerColor(colors[1]);
    hIntegral_Sng->SetLineWidth(3);
    hIntegral_Sng->SetMarkerStyle(markers[2]);
    hIntegral_Sng->SetBinContent(1,f1DIntegralError);
    hIntegral_Sng->SetBinContent(2,f2DIntegralError);
    hIntegral_Sng->SetBinError  (1,0);
    hIntegral_Sng->SetBinError  (2,0);
    hIntegral_Sng->SetBinError  (3,0);
    hIntegral_Sng->SetBinError  (4,0);
    hIntegral_Sng->SetBinError  (5,0);
    hIntegral_Sng->SetBinError  (6,0);
    hIntegral_Sng->SetMinimum   (0.);
    hIntegral_Sng->Scale(100);
    //
    TH1F   *hCombined_Sng   =   new TH1F("hCombined_Sng",   "", 6,  0,  6);
    hCombined_Sng->SetMarkerColor(colors[2]);
    hCombined_Sng->SetLineWidth(3);
    hCombined_Sng->SetMarkerStyle(markers[3]);
    hCombined_Sng->SetBinContent(1,sqrt(f1DIntegralError*f1DIntegralError*k1DStandardIntegr*k1DStandardIntegr + f1DExtrapolError*f1DExtrapolError*k1DStandardExtrap[0]*k1DStandardExtrap[0]) / (k1DStandardIntegr+k1DStandardExtrap[0]));
    hCombined_Sng->SetBinContent(2,sqrt(f2DIntegralError*f2DIntegralError*k2DStandardIntegr*k2DStandardIntegr + f2DExtrapolError*f2DExtrapolError*k2DStandardExtrap*k2DStandardExtrap) / (k2DStandardIntegr+k2DStandardExtrap));
    hCombined_Sng->SetBinContent(3,sqrt(f1DCalculatError*f1DCalculatError+f2DCalculatError*f2DCalculatError));
    hCombined_Sng->SetBinContent(4,sqrt(4*f1DCalculatError*f1DCalculatError+f2DCalculatError*f2DCalculatError));
    hCombined_Sng->SetBinContent(5,fSigmaPhiError(k1DFullYieldStand,0.5*k2DFullYieldStand,k1DFullYieldStand*f1DCalculatError,0.5*k2DFullYieldStand*f2DCalculatError)/fSigmaPhiValue(k1DFullYieldStand,0.5*k2DFullYieldStand));
    hCombined_Sng->SetBinContent(6,fGammaPhiError(k1DFullYieldStand,0.5*k2DFullYieldStand,k1DFullYieldStand*f1DCalculatError,0.5*k2DFullYieldStand*f2DCalculatError)/fGammaPhiValue(k1DFullYieldStand,0.5*k2DFullYieldStand));
    hCombined_Sng->SetBinError  (1,0);
    hCombined_Sng->SetBinError  (2,0);
    hCombined_Sng->SetBinError  (3,0);
    hCombined_Sng->SetBinError  (4,0);
    hCombined_Sng->SetBinError  (5,0);
    hCombined_Sng->SetBinError  (6,0);
    hCombined_Sng->Scale(100);
    //
    TH1F   *hExtrapolate_Sng   =   new TH1F("hExtrapolate_Sng",   "", 6,  0,  6);
    hExtrapolate_Sng->SetMarkerColor(colors[3]);
    hExtrapolate_Sng->SetLineWidth(3);
    hExtrapolate_Sng->SetMarkerStyle(markers[4]);
    hExtrapolate_Sng->SetBinContent(1,f1DExtrapolError);
    hExtrapolate_Sng->SetBinContent(2,f2DExtrapolError);
    hExtrapolate_Sng->SetBinError  (1,0);
    hExtrapolate_Sng->SetBinError  (2,0);
    hExtrapolate_Sng->SetBinError  (3,0);
    hExtrapolate_Sng->SetBinError  (4,0);
    hExtrapolate_Sng->SetBinError  (5,0);
    hExtrapolate_Sng->SetBinError  (6,0);
    hExtrapolate_Sng->Scale(100);
    //
    hCalculateFull_Sng->SetMaximum( 1.25*max( hCalculateFull_Sng->GetMaximum(), max ( hCombined_Sng->GetMaximum(), max ( hExtrapolate_Sng->GetMaximum() , hIntegral_Sng->GetMaximum() ) ) ) );
    //
    TLegend*    lLegend =   new TLegend(0.2,0.77,0.7,0.92);
    lLegend     ->  SetFillColorAlpha(0.,0.);
    lLegend     ->  SetLineColorAlpha(0.,0.);
    lLegend     ->  SetNColumns(4);
    lLegend     ->  AddEntry(hIntegral_Sng,"Int Err","P");
    lLegend     ->  AddEntry(hExtrapolate_Sng,"Ext Err","P");
    lLegend     ->  AddEntry(hCombined_Sng,"Sq(I+E) Err","P");
    lLegend     ->  AddEntry(hCalculateFull_Sng,"Calc Err","P");
    //
    TCanvas*    c1  = new TCanvas("2","22",1200,1200);
    hCalculateFull_Sng->Draw("EP SAME MIN0");
    hIntegral_Sng->Draw("EP SAME");
    hExtrapolate_Sng->Draw("EP SAME");
    hCombined_Sng->Draw("EP SAME");
    lLegend->Draw("same");
    c1->SaveAs(Form("%s/plots/Ratio_Analysis.pdf",fFolder.Data()));
    delete c1;
    //
    TFile      *fOutput =   new TFile   (Form("%s%s",fFolder.Data(),"/RT_Systematic.root"),"recreate");
    hCalculateFull_Sng->Scale(0.01);
    hCalculateFull_Sng->Write();
    fOutput->Close();
    //
    gROOT->SetBatch(kFALSE);
    //
    gErrorIgnoreLevel   =   kPrint;
    //
    return;
}
//
void
GeneralAnalysis
( ) {
    //
    // --------- RECOVER CALCULATED UNCERTAINTIES
    //
    SetStyle();
    fSetAllBins ();
    gROOT   ->  ProcessLine ( Form(".! mkdir -p %s",(TString(Form(kAnalysis_Systemt_Dir,"yield"))).Data()) );
    gROOT   ->  ProcessLine ( Form(".! mkdir -p %s",(TString(Form(kSystematicsPlot,"yield"))).Data()) );
    //
    // --------- BUILDING CONSTANT UNCERTAINTIES HISTOGRAMS
    //
    // ---------  --------- BRANCHING RATIO
    //
    TH1F   *h1DBranchingRatio     =   new TH1F    ("h1DBranchingRatio","h1DBranchingRatio",nBinPT1D,fArrPT1D);
    for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )    {
        h1DBranchingRatio->SetBinContent(iPT1D+1,kSysLow_BR);
    }
    h1DBranchingRatio           ->  SetLineStyle(1);
    h1DBranchingRatio           ->  SetLineColor(colors[1]);
    //
    TH2F   *h2DBranchingRatio   =   new TH2F    ("h2DBranchingRatio","h2DBranchingRatio",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )    {
        for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )    {
            h2DBranchingRatio->SetBinContent(iPT2D+1,jPT2D+1,kSysLow_BR*2);
        }
    }
    //
    TH1F   *hRTBranchingRatio     =   new TH1F    ("hRTBranchingRatio","hRTBranchingRatio",6,0,6);
    //
    TH1F   *hRTNormalisationH     =   new TH1F    ("hRTNormalisation_","hRTNormalisation_",6,0,6);
    TH1F   *hRTNormalisationL     =   new TH1F    ("hRTNormalisation_","hRTNormalisation_",6,0,6);
    hRTNormalisationH->SetBinContent(1,+1.*kSysHig_TE);
    hRTNormalisationH->SetBinContent(2,+1.*kSysHig_TE);
    hRTNormalisationH->SetBinContent(3,+1.*0.);
    hRTNormalisationH->SetBinContent(4,+0.0451296/1.41657);
    hRTNormalisationH->SetBinContent(5,+0.00246015/0.0349741);
    hRTNormalisationH->SetBinContent(6,+0.00116127/0.0604572);
    hRTNormalisationL->SetBinContent(1,-1.*kSysLow_TE);
    hRTNormalisationL->SetBinContent(2,-1.*kSysLow_TE);
    hRTNormalisationL->SetBinContent(3,-1.*0.);
    hRTNormalisationL->SetBinContent(4,-0.0435496/1.41657);
    hRTNormalisationL->SetBinContent(5,-0.00119453/0.0349741);
    hRTNormalisationL->SetBinContent(6,-0.00239997/0.0604572);
    hRTNormalisationH                 ->  SetLineWidth(1);
    hRTNormalisationH                 ->  SetLineStyle(1);
    hRTNormalisationH                 ->  SetLineColor(colors[6]);
    hRTNormalisationH                 ->  SetFillColorAlpha(colors[6],0.5);
    hRTNormalisationL                 ->  SetLineWidth(1);
    hRTNormalisationL                 ->  SetLineStyle(1);
    hRTNormalisationL                 ->  SetLineColor(colors[6]);
    hRTNormalisationL                 ->  SetFillColorAlpha(colors[6],0.5);
    hRTNormalisationH                 ->  Scale(100);
    hRTNormalisationL                 ->  Scale(100);
    //
    // --------- RECOVERING CALCULATED UNCERTAINTIES HISTOGRAMS
    //
    //
    // ---------  --------- SIGNAL EXTRACTION
    TFile      *f1DSignalExtraction     =   new TFile(Form("%s%s",Form(kAnalysis_SgExSys_Dir,"yield/Systematics/Standard/"),"/1D_Systematic.root"));
    TH1F       *h1DSignalExtraction     =   (TH1F*)(f1DSignalExtraction->Get("hRAW_1D"));
    h1DSignalExtraction                 ->  SetLineStyle(1);
    h1DSignalExtraction                 ->  SetLineColor(colors[2]);
    //
    TFile      *f2DSignalExtraction     =   new TFile(Form("%s%s",Form(kAnalysis_SgExSys_Dir,"yield/Systematics/Standard/"),"/2D_Systematic.root"));
    TH2F       *h2DSignalExtraction     =   (TH2F*)(f2DSignalExtraction->Get("hRAW_2D"));
    //
    TFile      *fRTSignalExtraction     =   new TFile(Form("%s%s",Form(kAnalysis_SgExSys_Dir,"yield/Systematics/Standard/"),"/RT_Systematic.root"));
    TH1F       *hRTSignalExtraction     =   (TH1F*)(fRTSignalExtraction->Get("hCalculateFull_Sng"));
    hRTSignalExtraction                 ->  SetLineWidth(1);
    hRTSignalExtraction                 ->  SetLineStyle(1);
    hRTSignalExtraction                 ->  SetLineColor(colors[2]);
    //
    // ---------  --------- PID
    TFile      *f1DParticleIdentif_     =   new TFile(Form("%s%s",Form(kAnalysis_Systemt_Dir,"yield"),"/PID/1D_Systematic.root"));
    TH1F       *h1DParticleIdentif_     =   (TH1F*)(f1DParticleIdentif_->Get("hRES_1D_Stat"));
    h1DParticleIdentif_                 ->  SetLineStyle(1);
    h1DParticleIdentif_                 ->  SetLineColor(colors[3]);
    //
    TFile      *f2DParticleIdentif_     =   new TFile(Form("%s%s",Form(kAnalysis_Systemt_Dir,"yield"),"/PID/2D_Systematic.root"));
    TH2F       *h2DParticleIdentif_     =   (TH2F*)(f2DParticleIdentif_->Get("h2DStandard"));
    //
    TFile      *fRTParticleIdentif_     =   new TFile(Form("%s%s",Form(kAnalysis_Systemt_Dir,"yield"),"/PID/RT_Systematic.root"));
    TH1F       *hRTParticleIdentif_     =   (TH1F*)(fRTParticleIdentif_->Get("hCalculateFull_Sng"));
    hRTParticleIdentif_                 ->  SetLineWidth(1);
    hRTParticleIdentif_                 ->  SetLineStyle(1);
    hRTParticleIdentif_                 ->  SetLineColor(colors[3]);
    //
    // ---------  --------- ANALYSIS CUTS
    TFile      *f1DAnalysisCuts____     =   new TFile(Form("%s%s",Form(kAnalysis_Systemt_Dir,"yield"),"/TRK/1D_Systematic.root"));
    TH1F       *h1DAnalysisCuts____     =   (TH1F*)(f1DAnalysisCuts____->Get("hRES_1D_Stat"));
    h1DAnalysisCuts____                 ->  SetLineStyle(1);
    h1DAnalysisCuts____                 ->  SetLineColor(colors[4]);
    //
    TFile      *f2DAnalysisCuts____     =   new TFile(Form("%s%s",Form(kAnalysis_Systemt_Dir,"yield"),"/TRK/2D_Systematic.root"));
    TH2F       *h2DAnalysisCuts____     =   (TH2F*)(f2DAnalysisCuts____->Get("h2DStandard"));
    //
    TFile      *fRTAnalysisCuts____     =   new TFile(Form("%s%s",Form(kAnalysis_Systemt_Dir,"yield"),"/TRK/RT_Systematic.root"));
    TH1F       *hRTAnalysisCuts____     =   (TH1F*)(fRTAnalysisCuts____->Get("hCalculateFull_Sng"));
    hRTAnalysisCuts____                 ->  SetLineWidth(1);
    hRTAnalysisCuts____                 ->  SetLineStyle(1);
    hRTAnalysisCuts____                 ->  SetLineColor(colors[4]);
    //
    // ---------  --------- GLOBAL TRACKING
    TFile      *f1DGlobalTracking__     =   new TFile(Form("%s%s",Form(kAnalysis_Systemt_Dir,"yield"),"/GTK/1D_Systematic.root"));
    TH1F       *h1DGlobalTracking__     =   (TH1F*)(f1DGlobalTracking__->Get("hRES_1D_Stat"));
    h1DGlobalTracking__                 ->  SetLineStyle(1);
    h1DGlobalTracking__                 ->  SetLineColor(colors[5]);
    //
    TFile      *f2DGlobalTracking__     =   new TFile(Form("%s%s",Form(kAnalysis_Systemt_Dir,"yield"),"/GTK/2D_Systematic.root"));
    TH2F       *h2DGlobalTracking__     =   (TH2F*)(f2DGlobalTracking__->Get("h2DStandard"));
    //
    TFile      *fRTGlobalTracking__     =   new TFile(Form("%s%s",Form(kAnalysis_Systemt_Dir,"yield"),"/GTK/RT_Systematic.root"));
    TH1F       *hRTGlobalTracking__     =   (TH1F*)(fRTGlobalTracking__->Get("hCalculateFull_Sng"));
    hRTGlobalTracking__                 ->  SetLineWidth(1);
    hRTGlobalTracking__                 ->  SetLineStyle(1);
    hRTGlobalTracking__                 ->  SetLineColor(colors[5]);
    //
    // --------- BUILDING TOTAL UNCERTAINTY HISTOGRAM
    //
    TH1F   *h1DTotalSystematic     =   new TH1F    ("h1DTotalSystematic","h1DTotalSystematic",nBinPT1D,fArrPT1D);
    for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )    {
        h1DTotalSystematic->SetBinContent(iPT1D+1,SquareSum( {h1DBranchingRatio->GetBinContent(iPT1D+1), h1DSignalExtraction->GetBinContent(iPT1D+1), h1DParticleIdentif_->GetBinContent(iPT1D+1), h1DAnalysisCuts____->GetBinContent(iPT1D+1), h1DGlobalTracking__->GetBinContent(iPT1D+1)} ));
    }
    h1DTotalSystematic->SetLineWidth(2);
    h1DTotalSystematic->SetLineColor(kBlack);
    h1DTotalSystematic->SetMinimum(.0);
    h1DTotalSystematic->SetMaximum(1.3*h1DTotalSystematic->GetMaximum());
    h1DTotalSystematic->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h1DTotalSystematic->GetYaxis()->SetTitleOffset(1.5);
    h1DTotalSystematic->GetYaxis()->SetTitle("Systematic Uncertainty (%)");
    //
    TH2F   *h2DTotalSystematic     =   new TH2F    ("h2DTotalSystematic","h2DTotalSystematic",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )    {
        for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )    {
            h2DTotalSystematic->SetBinContent(iPT2D+1,jPT2D+1,SquareSum( {h2DBranchingRatio->GetBinContent(iPT2D+1,jPT2D+1), h2DSignalExtraction->GetBinContent(iPT2D+1,jPT2D+1), h2DParticleIdentif_->GetBinContent(iPT2D+1,jPT2D+1), h2DAnalysisCuts____->GetBinContent(iPT2D+1,jPT2D+1), h2DGlobalTracking__->GetBinContent(iPT2D+1,jPT2D+1) } ));
        }
    }
    //
    TH1F   *hRTTotalSystematic     =   (TH1F*)(hRTSignalExtraction->Clone());
    hRTTotalSystematic->GetXaxis()->SetBinLabel(hRTTotalSystematic->GetXaxis()->FindBin(4.5),"#sigma^{2}_{#phi}");
    for ( Int_t iPTRT = 0; iPTRT < 6; iPTRT++ )    {
        hRTBranchingRatio->SetBinContent(iPTRT+1, ((0.005)/(0.08))*hRTGlobalTracking__->GetBinContent(iPTRT+1));
        hRTBranchingRatio->SetBinError(iPTRT+1, 0);
        hRTTotalSystematic->SetBinContent(iPTRT+1,SquareSum( { hRTSignalExtraction->GetBinContent(iPTRT+1), hRTParticleIdentif_->GetBinContent(iPTRT+1), hRTAnalysisCuts____->GetBinContent(iPTRT+1), hRTGlobalTracking__->GetBinContent(iPTRT+1), hRTBranchingRatio->GetBinContent(iPTRT+1)} ));
    }
    hRTTotalSystematic->SetLineWidth(2);
    hRTTotalSystematic->SetLineColor(kBlack);
    hRTTotalSystematic->Scale(100);
    hRTTotalSystematic->SetMaximum(1.4*hRTTotalSystematic->GetMaximum());
    hRTTotalSystematic->SetMinimum(1.4*hRTNormalisationL->GetMinimum());
    //
    // --------- FINAL CANVAS 1D
    //
    gROOT->SetBatch();
    //
    TLegend    *lLegend                 =   new TLegend(0.17,0.87,0.82,0.76);
    lLegend     ->  SetFillColorAlpha(kWhite,0.);
    lLegend     ->  SetLineColorAlpha(kWhite,0.);
    lLegend     ->  SetNColumns(4);
    lLegend     ->  AddEntry(h1DTotalSystematic,    "Total (No Norm.)","L");
    lLegend     ->  AddEntry(h1DBranchingRatio,     "Branching Ratio","L");
    lLegend     ->  AddEntry(h1DSignalExtraction,   "Signal Extraction","L");
    lLegend     ->  AddEntry(h1DParticleIdentif_,   "PID","L");
    lLegend     ->  AddEntry(h1DGlobalTracking__,   "Global Tracking","L");
    lLegend     ->  AddEntry(h1DAnalysisCuts____,   "Analysis cuts","L");
    //
    TCanvas    *cFullSyst               =   new TCanvas();
    gPad                                ->  SetLogx();
    h1DTotalSystematic                  ->  Draw();
    h1DBranchingRatio                   ->  Draw("same");
    h1DSignalExtraction                 ->  Draw("same");
    h1DParticleIdentif_                 ->  Draw("same");
    h1DAnalysisCuts____                 ->  Draw("same");
    h1DGlobalTracking__                 ->  Draw("same");
    lLegend                             ->  Draw("same");
    cFullSyst                           ->  SaveAs(Form("%s/1DFullSyst.pdf",Form(kSystematicsPlot,"yield")));
    delete cFullSyst;
    //
    // --------- FINAL CANVAS 2D
    //
    gROOT->SetBatch();
    //
    for ( Int_t iPT = 0; iPT < nBinPT2D; iPT++ )    {
        cFullSyst                           =   new TCanvas();
        gPad                                ->  SetLogx();
        //
        auto hTotal = h2DTotalSystematic    ->  ProjectionY(Form("2DFL_%i",iPT),iPT+1,iPT+1);
        hTotal  ->SetLineWidth(2);
        hTotal  ->SetLineColor(kBlack);
        hTotal  ->SetMinimum(.0);
        hTotal  ->SetMaximum(1.3*hTotal->GetMaximum());
        hTotal  ->GetXaxis()->SetTitle("#it{p}_{T,#phi_{1}} (GeV/#it{c})");
        hTotal  ->GetYaxis()->SetTitleOffset(1.5);
        hTotal  ->GetYaxis()->SetTitle("Systematic Uncertainty (%)");
        hTotal  ->Draw();
        uLatex  ->DrawLatexNDC(0.15,0.05,Form("#it{p}_{T,#phi_{2}} (GeV/#it{c}) #in [%.1f;%.1f]",fArrPT2D[iPT],fArrPT2D[iPT+1]));
        //
        auto hBranch=   h2DBranchingRatio   ->  ProjectionY(Form("2DBR_%i",iPT),iPT+1,iPT+1);
        hBranch ->  SetLineStyle(1);
        hBranch ->  SetLineColor(colors[1]);
        hBranch ->  Draw("SAME");
        //
        auto hSigEx =   h2DSignalExtraction ->  ProjectionY(Form("2DSE_%i",iPT),iPT+1,iPT+1);
        hSigEx  ->  SetLineStyle(1);
        hSigEx  ->  SetLineColor(colors[2]);
        hSigEx  ->  Draw("SAME");
        //
        auto hPID =   h2DParticleIdentif_   ->  ProjectionY(Form("2DPD_%i",iPT),iPT+1,iPT+1);
        hPID  ->  SetLineStyle(1);
        hPID  ->  SetLineColor(colors[3]);
        hPID  ->  Draw("SAME");
        //
        auto hAnCut =   h2DAnalysisCuts____ ->  ProjectionY(Form("2DAC_%i",iPT),iPT+1,iPT+1);
        hAnCut  ->  SetLineStyle(1);
        hAnCut  ->  SetLineColor(colors[4]);
        hAnCut  ->  Draw("SAME");
        //
        auto hITSPC =   h2DGlobalTracking__ ->  ProjectionY(Form("2DIT_%i",iPT),iPT+1,iPT+1);
        hITSPC  ->  SetLineStyle(1);
        hITSPC  ->  SetLineColor(colors[5]);
        hITSPC  ->  Draw("SAME");
        //
        lLegend                             ->  Draw("same");
        cFullSyst                           ->  SaveAs(Form("%s/2DFullSyst_%i.pdf",Form(kSystematicsPlot,"yield"),iPT));
        delete cFullSyst;
    }
    //
    TH1F   *hRTBaseLineInfo             =   new TH1F    ("hRTBranchingRatio","hRTBranchingRatio",6,0,6);
    hRTBaseLineInfo                     ->  SetLineWidth(1);
    hRTBaseLineInfo                     ->  SetLineStyle(3);
    hRTBaseLineInfo                     ->  SetLineColor(colors[0]);
    //
    cFullSyst               =   new TCanvas();
    hRTSignalExtraction                 ->  Scale(100.);
    hRTParticleIdentif_                 ->  Scale(100.);
    hRTAnalysisCuts____                 ->  Scale(100.);
    hRTGlobalTracking__                 ->  Scale(100.);
    hRTBranchingRatio                   ->  Scale(100.);
    hRTTotalSystematic                  ->  Draw("SAME MIN0");
    hRTSignalExtraction                 ->  Draw("SAME");
    hRTParticleIdentif_                 ->  Draw("SAME");
    hRTAnalysisCuts____                 ->  Draw("SAME");
    hRTGlobalTracking__                 ->  Draw("SAME");
    hRTBranchingRatio                   ->  Draw("SAME");
    hRTNormalisationH                   ->  Draw("SAME HIST B");
    hRTNormalisationL                   ->  Draw("SAME HIST B");
    hRTBaseLineInfo                     ->  Draw("SAME");
    lLegend                             ->  AddEntry( hRTNormalisationL, "Normalisation", "F" );
    lLegend                             ->  Draw("SAME");
    cFullSyst                           ->  SaveAs(Form("%s/RTFullSyst.pdf",Form(kSystematicsPlot,"yield")));
    delete cFullSyst;
    //
    gROOT->SetBatch(kFALSE);
    //
    TFile      *fOutput =   new TFile( Form("%s/Full_Systematics.root",(TString(Form(kAnalysis_Systemt_Dir,"yield"))).Data()), "RECREATE" );
    h1DTotalSystematic->Write();
    h2DTotalSystematic->Write();
    hRTTotalSystematic->Scale(0.01);
    hRTTotalSystematic->Write();
    fOutput->Close();
}
//
