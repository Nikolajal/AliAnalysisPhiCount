#include "../../inc/AliAnalysisPhiPair.h"

//  RopeFragPars::
// Initial step size for a calculation.
const double DELTAA = 0.1;
// Convergence criterion for a calculation.
const double ACONV = 0.001;
// Low z cut-off in fragmentation function.
const double ZCUT = 1.0e-4;
Float_t
uCalculateAlpha
( Float_t kRho, Float_t kX, Float_t kY ) {
    return (1+2*kRho*kX+9*kY+6*kX*kY*kRho+3*kX*kX*kY*kRho*kRho)/(2+kRho);
}
Double_t
fragf(double z, double a, double b, double mT2) {
  if (z < ZCUT) return 0.0;
  return pow(1 - z, a) * exp(-b * mT2 / z) / z;
}
Double_t
trapIntegrate
 ( double a, double b, double mT2, double sOld, int n ) {

  // Compute the nth correction to the integral of fragfunc between 0 and 1
  // using extended trapezoidal rule.
  if (n == 1) return 0.5 * (fragf(0.0, a, b, mT2) + fragf(1.0, a, b, mT2));
  // We want 2^(n-2) interior points (intp). Use bitwise shift to speed up.
  int intp = 1;
  intp <<= n - 2;
  double deltaz = 1.0 / double(intp);
  double z = 0.5 * deltaz;
  double sum = 0.0;
  // Do the integral.
  for (int i = 0; i < intp; ++i, z += deltaz) sum += fragf( z, a, b, mT2);
  return 0.5 * (sOld + sum / double(intp));

}
Double_t
integrateFragFun(double a, double b, double mT2) {

  // Using Simpson's rule to integrate the Lund fragmentation function.
  double nextIter, nextComb;
  double thisComb = 0.0, thisIter = 0.0;
  // The target error on the integral should never be changed.
  double error = 1.0e-2;

  // 20 is the max number of iterations, 3 is min. Should not be changed.
  for (int i = 1; i < 20; ++i) {
    nextIter = trapIntegrate( a, b, mT2, thisIter, i);
    nextComb = (4.0 * nextIter - thisIter) / 3.0;
    if (i > 3 && abs(nextComb - thisComb) < error * abs(nextComb))
      return nextComb;
    thisIter = nextIter;
    thisComb = nextComb;
  }
  return 0.0;

}
Double_t
aEffective
 (double aOrig, double thisb, double mT2, Float_t kDefaultB) {

  // Calculate initial normalization constants.
  double N    = integrateFragFun(aOrig, kDefaultB, mT2);
  double NEff = integrateFragFun(aOrig, thisb, mT2);
  int    s    = (N < NEff) ? -1 : 1;
  double da   = DELTAA;
  double aNew = aOrig - s * da;

  // Iterate until we meet preset convergence criterion.
  do {
    // Calculate normalization with current a.
    NEff = integrateFragFun(aNew, thisb, mT2);
    if ( ((N < NEff) ? -1 : 1) != s ) {
      s = (N < NEff) ? -1 : 1;
      // If we have crossed over the solution, decrease
      // the step size and turn around.
      da /= 10.0;
    }
    aNew -= s * da;
    if (aNew < 0.0) {aNew = 0.1; break;}
    if (aNew > 2.0) {aNew = 2.0; break;}
  } while (da > ACONV);
  return aNew;

}
Float_t
uCalculateA
( Float_t kRe_CalcB, Float_t mT2, Bool_t isDiquark, Float_t kDefaultA, Float_t kDefaultAdiq, Float_t kDefaultB ) {
    return aEffective( kDefaultA + ( isDiquark? kDefaultAdiq : 0. ), kRe_CalcB, mT2, kDefaultB );
}
std::map<TString,Float_t>
uCalculateStringEnhancement
( Float_t kRelativeTension = 1., Float_t kUsedBeta = 0.2, Float_t kUsedKappa = 0.2 ) {
    //! Output
    std::map<TString,Float_t> kResultMap;
    //!
    //cout << "[INFO] Requested a custom string tension enhacnement of: " << kRelativeTension << endl;
    if ( kRelativeTension < 0. ) {
        cout << "[ERROR] Invalid query, enhancement factor must be positive. Nothing done." << endl;
        return kResultMap;
    }
    /*
    if ( kRelativeTension == 1. ) {
        cout << "[WARNING] Invalid query, enhancement factor must be different than 1. Nothing done." << endl;
        return kResultMap;
    }
     */
    //! Defining h and 1/h for utility
    const Float_t kEnhancement    = kRelativeTension;
    const Float_t kInvEnhancement = 1./kRelativeTension;
    //!
    //! --- --- Standard Parameter
    //!
    //! Kappa: A base value of the string tension can be added, and modified along with other parameters, to allow for studies of exotic quark production in the Rope model.
    //! StringFlav:kappa
    const Float_t kDefaultKappa = 0.2;
    //! Rho: Suppression of s quark production relative to u or d type production
    //! StringFlav:probStoUD
    const Float_t kDefaultRho   = 0.217;
    //! Xi: Suppression of diquark production relative to quark production
    //! StringFlav:probQQtoQ
    const Float_t kDefaultXi    = 0.081;
    //! X: Suppression of strange diquark production relative to light diquark production
    //! StringFlav:probSQtoQQ
    const Float_t kDefaultX     = 0.915;
    //! Y: Suppression of spin 1 diquarks relative to spin 0 diquarks
    //! StringFlav:probQQ1toQQ0
    const Float_t kDefaultY     = 0.0275;
    //! Sigma: Width of the transverse momentum distribution in string break-ups
    //! StringPT:sigma
    const Float_t kDefaultSigma = 0.335;
    //! Alpha: Derived parameter for Xi calculation
    //! -
    const Float_t kDefaultAlpha = uCalculateAlpha( kDefaultRho, kDefaultX, kDefaultY );
    //! Beta: In the current implementation of the rope model, the theoretical ignorance about baryon production has been parameterized, assuming that the parameter StringFlav:probQQtoQ will factorize into two parts, one which will scale with effective string tension, one which will not. This parameter controls how large a fraction of the parameter will scale with string tension.
    //! Ropewalk:beta
    const Float_t kDefaultBeta  = 0.2;
    //! A: The a parameter of the Lund symmetric fragmentation function.
    //! StringZ:aLund
    const Float_t kDefaultA     = 0.680;
    //! Adiq: The a parameter of the Lund symmetric fragmentation function.
    //! StringZ:aExtraDiquark
    const Float_t kDefaultAdiq  = 0.970;
    //! B: The a parameter of the Lund symmetric fragmentation function.
    //! StringZ:bLund
    const Float_t kDefaultB     = 0.980;
    //!
    //! --- --- Re-calculated parameters
    //!
    const Float_t kRe_CalcKappa = kUsedKappa*kEnhancement;
    const Float_t kRe_CalcBeta  = kUsedBeta;
    const Float_t kRe_CalcRho   = pow( kDefaultRho,     kInvEnhancement );
    const Float_t kRe_CalcX     = pow( kDefaultX,       kInvEnhancement );
    const Float_t kRe_CalcY     = pow( kDefaultY,       kInvEnhancement );
    const Float_t kRe_CalcSigma = pow( kDefaultSigma,   kInvEnhancement );
    const Float_t kRe_CalcAlpha = uCalculateAlpha( kRe_CalcRho, kRe_CalcX, kRe_CalcY );
    const Float_t kRe_CalcXi_U  = kRe_CalcAlpha*kRe_CalcBeta*pow( kDefaultXi/(kDefaultAlpha*kRe_CalcBeta), kInvEnhancement );
    const Float_t kRe_CalcXi    = min( 1., 1.*max( kRe_CalcXi_U, kDefaultXi ) );
    const Float_t kRe_CalcB_U   = kDefaultB*(2+kRe_CalcRho)/(2+kDefaultRho);
    const Float_t kRe_CalcB     = min( 2., 1.*max( kRe_CalcB_U, kDefaultXi ) );
    const Float_t kRe_CalcA     = uCalculateA( kRe_CalcB, 1.0, false, kDefaultA, kDefaultAdiq, kDefaultB );
    const Float_t kRe_CalcAdiq  = uCalculateA( kRe_CalcB, 1.0, true, kDefaultA, kDefaultAdiq, kDefaultB ) - kRe_CalcA;
    //!
    //! --- --- Verbose info
    //!
    //cout << "[WARNING] !! Preliminary implementation of custom string tension !! Please use caution and cross-check results" << endl;
    //cout << "[INFO] Def.    K:"     << kUsedKappa       << "        Re-Calc. K:"    << kRe_CalcKappa    << endl;
    //cout << "[INFO] Def.    Beta:"  << kDefaultBeta     << "     Re-Calc. Beta:"    << kRe_CalcBeta     << endl;
    //cout << "[INFO] Def.    Rho:"   << kDefaultRho      << "      Re-Calc. Rho:"    << kRe_CalcRho      << endl;
    //cout << "[INFO] Def.    X:"     << kDefaultX        << "        Re-Calc. X:"    << kRe_CalcX        << endl;
    //cout << "[INFO] Def.    Y:"     << kDefaultY        << "        Re-Calc. Y:"    << kRe_CalcY        << endl;
    //cout << "[INFO] Def.    Sigma:" << kDefaultSigma    << "    Re-Calc. Sigma:"    << kRe_CalcSigma    << endl;
    //cout << "[INFO] Def*.   Alpha:" << kDefaultAlpha    << "    Re-Calc*.Alpha:"    << kRe_CalcAlpha    << endl;
    //cout << "[INFO] Def.    Xi:"    << kDefaultXi       << "       Re-Calc. Xi:"    << kRe_CalcXi       << "    No-Cut Xi:" << kRe_CalcXi_U << endl;
    //cout << "[INFO] Def.    A:"     << kDefaultA        << "        Re-Calc. A:"    << kRe_CalcA        << endl;
    //cout << "[INFO] Def.    Adiq:"  << kDefaultAdiq     << "     Re-Calc. Adiq:"    << kRe_CalcAdiq     << endl;
    //cout << "[INFO] Def.    B:"     << kDefaultB        << "        Re-Calc. B:"    << kRe_CalcB        << "     No-Cut B:" << kRe_CalcB_U << endl;
    //!
    //! --- --- Setting parameters
    //!
    kResultMap["kappa"] = kRe_CalcKappa;
    kResultMap["Rho"]   = kRe_CalcRho;
    kResultMap["Xi_U"]  = kRe_CalcXi_U;
    kResultMap["Xi"]    = kRe_CalcXi;
    kResultMap["X"]     = kRe_CalcX;
    kResultMap["Y"]     = kRe_CalcY;
    kResultMap["Sigma"] = kRe_CalcSigma;
    kResultMap["Beta"]  = kRe_CalcBeta;
    kResultMap["A"]     = kRe_CalcA;
    kResultMap["Adiq"]  = kRe_CalcAdiq;
    kResultMap["B_U"]   = kRe_CalcB_U;
    kResultMap["B"]     = kRe_CalcB;
    //!
    return kResultMap;
}

void AN_ShowPlots ( TString fOption = "all", TString kFolder = "_p_p__5TeV" , std::vector<TString> kComparisonMCTagList = { /*"Pythia8X_0_0", "Pythia8X_0_1", "Pythia8X0_"*/ }, std::vector<TString> kComparisonLegendList = { /*"P8M13R1", "P8M13R2", "P8M13RHome"*/ } )    {
    // --- --- --- --- --- --- --- SET-UP --- --- --- --- --- --- --- --- --- --- ---
    //
    //  --- INFO on Set-up variables
    fChooseOption(fOption);
    //
    //  --- Setting the input datastructure
    fSetAllBins();
    //
    //  --- Setting the style
    SetStyle();
    //
    // --- YIELD ANALYSIS
    if ( kDoYield ) {
        //
        //  --- Load Data Files
        TFile*  insFile_Data_YL = new TFile ( Form(kASigExtp_FitCheckRst,(TString("Yield")+kFolder).Data()) );
        TFile*  insFile_Syst_YL = new TFile ( Form("%s/FullSystematics.root",Form(kAnalysis_Systemt_Dir,  (TString("Yield")+kFolder).Data())) );
        TFile*  insFile_Data_CP = new TFile ( Form(kASigExtp_FitCheckRst,(TString("Correlation")+kFolder).Data()) );
        TFile*  insFile_Data_C2 = new TFile ( Form(kASigExtp_FitCheckRst,(TString("Correlation")+kFolder+TString("_InvMassTest")).Data()) );
        TFile*  outFile_Dump    = new TFile ( "./temp.root", "RECREATE" );
        //
        //  --- Load Data Histograms
        auto    hXD_Nyld_stat   = uLoadHistograms<0,TH1F> ( insFile_Data_YL, "hXD_Nyld_stat" );
        auto    hXD_Nyld_syst   = uLoadHistograms<0,TH1F> ( insFile_Data_YL, "hXD_Nyld_syst" );
        auto    hXD_Nfqs_stat   = uLoadHistograms<0,TH1F> ( insFile_Data_YL, "hXD_Nfqs_stat" );
        auto    hXD_Nfqs_syst   = uLoadHistograms<0,TH1F> ( insFile_Data_YL, "hXD_Nfqs_syst" );
        auto    h1D_Nres_stat   = uLoadHistograms<0,TH1F> ( insFile_Data_YL, "h1D_Nres_stat" );
        auto    h1D_Nres_syst   = uLoadHistograms<0,TH1F> ( insFile_Data_YL, "h1D_Nres_syst" );
        auto    h2D_Nres_stat   = uLoadHistograms<1,TH1F> ( insFile_Data_YL, "h2D_Nres_stat_stat_%i" );
        auto    h2D_Nres_syst   = uLoadHistograms<1,TH1F> ( insFile_Data_YL, "h2D_Nres_syst_syst_%i" );
        auto    h2D_MeanPT_stat = uLoadHistograms<0,TH1F> ( insFile_Data_YL, "h2D_MeanPT_stat" );
        auto    h2D_MeanPT_syst = uLoadHistograms<0,TH1F> ( insFile_Data_YL, "h2D_MeanPT_syst" );
        auto    h2D_PhiCorrelat = uLoadHistograms<0,TH1F> ( insFile_Data_CP, "h1D_PhiCorr_Stat" );
        auto    h2D_PhiCorrela2 = uLoadHistograms<0,TH1F> ( insFile_Data_C2, "h1D_PhiCorr_Stat" );
        //
        // --- MC Elements
        std::vector<TFile*> insFile_MntC_YL;
        std::vector<TH1F*> hMPT_NTru;
        std::vector<TH1F*> hMPT_NTru_Fin;
        std::vector<TH1F*> h1D_Ntru_Fin;
        std::vector<TH2F*> h2D_Ntru_Fin;
        std::vector<TH2F*> h2D_Ntru_Fin_Nrm;
        std::vector<TH1F*> h1D_Ntru;
        std::vector<TH2F*> h2D_Ntru;
        std::vector<TH1F*> hFullQuantities;
        std::vector<TH1F*> hProductionProb;
        std::vector<TH1F*> h2D_PhiCorrelation;
        TLegend*           lMCLegend_1D;
        TLegend*           lMCLegend_2D;
        TLegend*           lMCLegend_PT;
        TLegend*           lMCLegend_PP;
        TLegend*           lMCLegend_CR;
        //
        if ( kComparisonMCTagList.size() ) {
            // --- Load MC Files
            for ( auto kCurrent_MC_File : kComparisonMCTagList ) insFile_MntC_YL.push_back( new TFile   ( Form(kProduction_MC_Ofl,(TString("Yield")+kFolder).Data(),kCurrent_MC_File.Data()) ) );
            //
            // --- Build Legend
            lMCLegend_1D    = new TLegend( 0.500, 0.700, 0.880, 0.550 );
            lMCLegend_2D    = new TLegend( 0.550, 0.650, 0.880, 0.450 );
            lMCLegend_PT    = new TLegend( 0.600, 0.030, 0.880, 0.220 );
            lMCLegend_PP    = new TLegend( 0.175, 0.175, 0.525, 0.375 );
            lMCLegend_CR    = new TLegend( 0.675, 0.675, 0.875, 0.875 );
            lMCLegend_CR    -> SetNColumns(2);
            //
            //  --- Load MC Histograms
            auto iClr = 0;
            for ( auto kFile : insFile_MntC_YL )   {
                iClr++;
                //  --- Production Histograms
                //  --- --- 1D
                h1D_Ntru                    .push_back( (TH1F*)((( kFile )->Get("h1D_Ntru"))        -> Clone(Form("h1D_Ntru_%s",kFile->GetName())) ));
                h1D_Ntru.at( iClr-1 )       -> SetLineWidth( 2 );
                h1D_Ntru.at( iClr-1 )       -> SetLineColor( uGetColor( iClr+1 ) );
                h1D_Ntru.at( iClr-1 )       -> SetFillColor( uGetFillColor( iClr+1 ) );
                h1D_Ntru_Fin                .push_back( (TH1F*)((( kFile )->Get("h1D_Ntru_Fin"))    -> Clone(Form("h1D_Ntru_Fin_%s",kFile->GetName())) ));
                h1D_Ntru_Fin.at( iClr-1 )   -> SetLineWidth( 2 );
                h1D_Ntru_Fin.at( iClr-1 )   -> SetLineColor( uGetColor( iClr+1 ) );
                h1D_Ntru_Fin.at( iClr-1 )   -> SetFillColor( uGetFillColor( iClr+1 ) );
                //
                //  --- --- 2D
                h2D_Ntru                    .push_back( (TH2F*)((( kFile )->Get("h2D_Ntru"))        -> Clone(Form("h2D_Ntru_%s",kFile->GetName())) ));
                h2D_Ntru.at( iClr-1 )       -> SetLineWidth( 2 );
                h2D_Ntru.at( iClr-1 )       -> SetLineColor( uGetColor( iClr+1 ) );
                h2D_Ntru.at( iClr-1 )       -> SetFillColorAlpha( uGetColor( iClr+1 ), 0.33 );
                h2D_Ntru_Fin                .push_back( (TH2F*)((( kFile )->Get("h2D_Ntru_Fin"))    -> Clone(Form("h2D_Ntru_Fin_%s",kFile->GetName())) ));
                h2D_Ntru_Fin.at( iClr-1 )   -> SetLineWidth( 2 );
                h2D_Ntru_Fin.at( iClr-1 )   -> SetLineColor( uGetColor( iClr+1 ) );
                h2D_Ntru_Fin.at( iClr-1 )   -> SetFillColor( uGetFillColor( iClr+1 ) );
                h2D_Ntru_Fin_Nrm            .push_back( (TH2F*)((( kFile )->Get("h2D_Ntru_Fin_Nrm"))-> Clone(Form("h2D_Ntru_Fin_Nrm_%s",kFile->GetName())) ));
                h2D_Ntru_Fin_Nrm.at( iClr-1 )-> SetLineWidth( 2 );
                h2D_Ntru_Fin_Nrm.at( iClr-1 )-> SetLineColor( uGetColor( iClr+1 ) );
                h2D_Ntru_Fin_Nrm.at( iClr-1 )-> SetFillColor( uGetFillColor( iClr+1 ) );
                //
                //  --- --- MPT
                hMPT_NTru                   .push_back( (TH1F*)((( kFile )->Get("hMPT_NTru"))       -> Clone(Form("hMPT_NTru_%s",kFile->GetName())) ));
                hMPT_NTru.at( iClr-1 )      -> SetLineWidth( 2 );
                hMPT_NTru.at( iClr-1 )      -> SetLineColor( uGetColor( iClr+1 ) );
                hMPT_NTru.at( iClr-1 )      -> SetFillColor( uGetFillColor( iClr+1 ) );
                hMPT_NTru_Fin               .push_back( (TH1F*)((( kFile )->Get("hMPT_NTru_Fin"))   -> Clone(Form("hMPT_NTru_Fin_%s",kFile->GetName())) ));
                hMPT_NTru_Fin.at( iClr-1 )  -> SetLineWidth( 2 );
                hMPT_NTru_Fin.at( iClr-1 )  -> SetLineColor( uGetColor( iClr+1 ) );
                hMPT_NTru_Fin.at( iClr-1 )  -> SetFillColor( uGetFillColor( iClr+1 ) );
                //
                //  --- --- Other
                hFullQuantities             .push_back( (TH1F*)((( kFile )->Get("hFullQuantities")) -> Clone(Form("hFullQuantities_%s",kFile->GetName())) ));
                hFullQuantities.at( iClr-1 )-> SetLineWidth( 2 );
                hFullQuantities.at( iClr-1 )-> SetLineColor( uGetColor( iClr+1 ) );
                hFullQuantities.at( iClr-1 )-> SetFillColor( uGetFillColor( iClr+1 ) );
                //
                hProductionProb             .push_back( (TH1F*)((( kFile )->Get("hProductionProb")) -> Clone(Form("hProductionProb_%s",kFile->GetName())) ));
                hProductionProb.at( iClr-1 )-> SetLineWidth( 2 );
                hProductionProb.at( iClr-1 )-> SetLineColor( uGetColor( iClr+1 ) );
                hProductionProb.at( iClr-1 )-> SetFillColor( uGetFillColor( iClr+1 ) );
                //
                h2D_PhiCorrelation             .push_back( (TH1F*)((( kFile )->Get("hPhiCorrelation")) -> Clone(Form("hPhiCorrelation_%s",kFile->GetName())) ));
                h2D_PhiCorrelation.at( iClr-1 )-> SetLineWidth( 2 );
                h2D_PhiCorrelation.at( iClr-1 )-> SetLineColor( uGetColor( iClr+1 ) );
                h2D_PhiCorrelation.at( iClr-1 )-> SetFillColor( uGetFillColor( iClr+1 ) );
                //
                lMCLegend_1D                -> AddEntry( hMPT_NTru.at( iClr-1 ), kComparisonLegendList.at( iClr-1 ), "L" );
                lMCLegend_2D                -> AddEntry( hMPT_NTru.at( iClr-1 ), kComparisonLegendList.at( iClr-1 ), "L" );
                lMCLegend_PT                -> AddEntry( hMPT_NTru.at( iClr-1 ), kComparisonLegendList.at( iClr-1 ), "L" );
                lMCLegend_PP                -> AddEntry( hMPT_NTru.at( iClr-1 ), kComparisonLegendList.at( iClr-1 ), "L" );
                lMCLegend_CR                -> AddEntry( hMPT_NTru.at( iClr-1 ), kComparisonLegendList.at( iClr-1 ), "L" );
            }
        }
        //
        //  --- Output directory
        TString kPlotDirectory          =   Form(kDIR_ShowPlots,(TString("Yield")+kFolder).Data());
        gROOT   ->  ProcessLine(Form(".! mkdir -p %s",kPlotDirectory.Data()));
        //
        //  --- Yield Plots
        //  --- --- Batch Mode
        gROOT   ->  SetBatch( kTRUE );
        //
        //  --- --- Production Probability
        TCanvas*    cDrawResults    = new TCanvas( "cDrawResults", "cDrawResults", 1200, 1000 );
        gPad                    -> SetLogy();
        //
        TH1F* hUtility = new TH1F("","", 5, 0., +5 );
        hUtility->SetMaximum(1.);
        hUtility->SetMinimum(1.e-4);
        hUtility->GetXaxis()->SetNdivisions(5);
        hUtility->GetXaxis()->SetTitle("Number of #phi mesons produced");
        hUtility->GetYaxis()->SetTitle("Event frequency");
        hUtility->Draw();
        if ( kComparisonMCTagList.size() ) {
            for ( auto kCurrent_PP : hProductionProb ) {
                kCurrent_PP-> SetMarkerStyle(0);
                kCurrent_PP-> ResetAttFill();
                kCurrent_PP-> DrawCopy("SAME HIST L");
            }
        }
        if ( kComparisonMCTagList.size() ) lMCLegend_PP    ->  Draw("SAME");
        //
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.55, 0.83,"ALICE Simulation");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        uLatex->DrawLatexNDC(0.68, 0.78,"pp #sqrt{#it{s}}= 7 TeV");
        uLatex->DrawLatexNDC(0.74, 0.73,"#phi, |#it{y}|<0.5");
        //
        cDrawResults    -> SaveAs( kPlotDirectory + TString("hProductionProb.pdf") );
        cDrawResults    -> SaveAs( kPlotDirectory + TString("hProductionProb.eps") );
        delete cDrawResults;
        //
        //  --- --- Final quantities
        uSetHisto( hXD_Nfqs_stat,       "FNL STAT 1D" );
        uSetHisto( hXD_Nfqs_syst,       "FNL SYST 1D" );
        //
        TLegend*    lLegend         = new TLegend( 0.43, 0.03, 0.62, 0.25 );
        lLegend     ->  SetLineColorAlpha(0.,0.);
        lLegend     ->  SetFillColorAlpha(0.,0.);
        lLegend     ->  AddEntry    (hXD_Nfqs_stat,"Measured","P");
        lLegend     ->  AddEntry    (hXD_Nfqs_stat,"Stat. Uncert.","EL");
        lLegend     ->  AddEntry    (hXD_Nfqs_syst,"Syst. Uncert.","F");
        //
                cDrawResults    = new TCanvas( "cDrawResults", "cDrawResults", 1200, 1000 );
        gPad                    -> SetLogy();
        //
        //  --- Upper Plot
        TPad*   kUpperPlot  =   new TPad("kUpperPlot_RT", "kUpperPlot", 0, 0.40, 1, 1.0);
        gStyle      -> SetOptStat(0);
        kUpperPlot  -> SetTopMargin(0.05);
        kUpperPlot  -> SetBottomMargin(0);
        kUpperPlot  -> SetFillColorAlpha( 0., 0. );
        kUpperPlot  -> SetLogy();
        kUpperPlot  -> cd();
        //
        hXD_Nfqs_syst   -> GetYaxis() -> SetTitleSize(0.06);
        hXD_Nfqs_syst   -> GetYaxis() -> SetTitleOffset(0.98);
        hXD_Nfqs_syst   -> GetYaxis() -> SetLabelSize(0.055);
        hXD_Nfqs_syst   -> Draw("SAME PE2");
        hXD_Nfqs_stat   -> Draw("SAME PE X0");
        if ( kComparisonMCTagList.size() ) {
            for ( auto kCurrent_FQ : hFullQuantities )  {
                kCurrent_FQ-> SetMarkerStyle(0);
                kCurrent_FQ-> ResetAttFill();
                kCurrent_FQ-> Draw("SAME EP");
            }
        }
        lLegend                 -> Draw("SAME");
        if ( kComparisonMCTagList.size() ) lMCLegend_PT            -> Draw("SAME");
        //
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.07);
        uLatex->DrawLatexNDC(0.18, 0.86,"ALICE Preliminary");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.18, 0.80,"pp #sqrt{#it{s}}= 7 TeV");
        uLatex->DrawLatexNDC(0.18, 0.74,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
        //
        //  --- Lower Plot
        TPad*   kLowerPlot  =   new TPad("kLowerPlot_RT", "kLowerPlot", 0, 0.0, 1, 0.40);
        kLowerPlot  -> SetTopMargin(0);
        kLowerPlot  -> SetBottomMargin(0.3);
        kLowerPlot  -> SetFillColorAlpha( 0., 0. );
        kLowerPlot  -> cd();
        //
        auto kUtility_syst  = uScale( hXD_Nfqs_syst, -1 );
        kUtility_syst -> SetName("tmp_syst_1D");
        uSetHisto( kUtility_syst, " FNL SYST 12D " );
        kUtility_syst -> SetFillColor   ( uGetFillColor(1) );
        kUtility_syst -> SetMarkerStyle     ( 0 );
        kUtility_syst -> GetYaxis() -> SetTitle("Model / Data");
        kUtility_syst -> GetYaxis() -> SetTitleSize(0.12);
        kUtility_syst -> GetYaxis() -> SetTitleOffset(0.55);
        kUtility_syst -> GetYaxis() -> SetLabelSize(0.09);
        kUtility_syst -> GetYaxis() -> SetNdivisions(7);
        kUtility_syst -> GetXaxis() -> SetLabelSize(0.14);
        kUtility_syst -> GetXaxis() -> SetLabelOffset(0.0115);
        kUtility_syst -> GetXaxis() -> SetTitleSize(0.01);
        kUtility_syst -> GetXaxis() -> SetTitleOffset(1.0);
        kUtility_syst -> SetMaximum( 1.7 );
        kUtility_syst -> SetMinimum( 0.3 );
        kUtility_syst -> Draw("SAME PE2");
        //
        auto kUtility_stat = uScale( hXD_Nfqs_stat, -1 );
        kUtility_syst -> SetName("tmp_stat_1D");
        kUtility_stat -> Draw("SAME PE X0");
        //
        if ( kComparisonMCTagList.size() ) {
            auto kUtilityPlot       = uScale( hXD_Nfqs_stat, 1., -2.);
            for ( auto kCurrent_FQ : hFullQuantities ) {
                auto kCurrentHist   =   (TH1F*)kCurrent_FQ -> Clone(Form("tmp_%s",kCurrent_FQ->GetName()));
                kCurrentHist    -> Divide( kUtilityPlot );
                kCurrentHist    -> DrawCopy("SAME EP");
            }
        }
        //
        cDrawResults    -> cd();
        kUpperPlot      -> Draw();
        kLowerPlot      -> Draw();
        cDrawResults    -> SaveAs( kPlotDirectory + TString("FullQuantitites.pdf") );
        cDrawResults    -> SaveAs( kPlotDirectory + TString("FullQuantitites.eps") );
        delete lLegend;
        delete kUpperPlot;
        delete kLowerPlot;
        delete cDrawResults;
        //
        //  --- --- Full 1D Spectrum //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        TH1F*   h1D_Nres_stat_full      = new TH1F( "h1D_Nres_stat_full", "h1D_Nres_stat_full", nBinPT1D+1, fArrPT1D_Comp );
        TH1F*   h1D_Nres_syst_full      = new TH1F( "h1D_Nres_syst_full", "h1D_Nres_syst_full", nBinPT1D+1, fArrPT1D_Comp );
        //
        h1D_Nres_stat_full  -> SetBinContent( 1, hXD_Nyld_stat->GetBinContent(4) / (fArrPT1D_Comp[1]-fArrPT1D_Comp[0]) );
        h1D_Nres_syst_full  -> SetBinContent( 1, hXD_Nyld_stat->GetBinContent(4) / (fArrPT1D_Comp[1]-fArrPT1D_Comp[0]) );
        h1D_Nres_stat_full  -> SetBinError  ( 1, hXD_Nyld_stat->GetBinError  (4) / (fArrPT1D_Comp[1]-fArrPT1D_Comp[0]) );
        h1D_Nres_syst_full  -> SetBinError  ( 1, hXD_Nyld_syst->GetBinError  (4) / (fArrPT1D_Comp[1]-fArrPT1D_Comp[0]) );
        h1D_Nres_syst_full  -> SetMaximum   ( h1D_Nres_syst->GetMaximum()*2.5 );
        h1D_Nres_syst_full  -> SetMinimum   ( h1D_Nres_syst->GetMinimum()*0.5 );
        uSetHisto( h1D_Nres_stat,       "SPT STAT 1D" );
        uSetHisto( h1D_Nres_syst,       "SPT SYST 1D" );
        uSetHisto( h1D_Nres_stat_full,  "SPT STAT 1D VAR1 " );
        uSetHisto( h1D_Nres_syst_full,  "SPT SYST 1D VAR1 " );
        //
        lLegend     = new TLegend( 0.20, 0.03, 0.38, 0.25 );
        lLegend     ->  SetLineColorAlpha(0.,0.);
        lLegend     ->  SetFillColorAlpha(0.,0.);
        lLegend     ->  AddEntry    (h1D_Nres_stat,"Measured","P");
        lLegend     ->  AddEntry    (h1D_Nres_stat_full,"Extrapolated","P");
        lLegend     ->  AddEntry    (h1D_Nres_stat,"Stat. Uncert.","EL");
        lLegend     ->  AddEntry    (h1D_Nres_syst,"Syst. Uncert.","F");
        //
        cDrawResults    = new TCanvas( "cDrawResults", "cDrawResults", 1200, 1000 );
        //
        //  --- Upper Plot
        kUpperPlot  =   new TPad("kUpperPlot_1D", "kUpperPlot", 0, 0.35, 1, 1.0);
        gStyle      -> SetOptStat(0);
        kUpperPlot  -> SetTopMargin(0.05);
        kUpperPlot  -> SetBottomMargin(0);
        kUpperPlot  -> SetFillColorAlpha( 0., 0. );
        kUpperPlot  -> SetLogy();
        kUpperPlot  -> cd();
        //
        h1D_Nres_syst_full  -> GetYaxis() -> SetTitleSize(0.06);
        h1D_Nres_syst_full  -> GetYaxis() -> SetTitleOffset(0.98);
        h1D_Nres_syst_full  -> GetYaxis() -> SetLabelSize(0.055);
        h1D_Nres_syst_full  -> Draw("SAME PE2");
        h1D_Nres_stat_full  -> Draw("SAME PE X0");
        h1D_Nres_syst       -> Draw("SAME PE2");
        h1D_Nres_stat       -> Draw("SAME PE X0");
        lLegend             -> Draw("SAME");
        if ( kComparisonMCTagList.size() ) lMCLegend_1D        -> Draw("SAME");
        if ( kComparisonMCTagList.size() ) {
            for ( auto kCurrent_1D : h1D_Ntru_Fin ) {
                kCurrent_1D-> DrawCopy("SAME E3");
                kCurrent_1D-> SetFillColor( 0. );
                kCurrent_1D-> ResetAttFill();
                kCurrent_1D-> DrawCopy("SAME HIST L");
            }
        }
        //
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.55, 0.83,"ALICE Preliminary");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        uLatex->DrawLatexNDC(0.55, 0.78,"pp #sqrt{#it{s}}= 7 TeV");
        uLatex->DrawLatexNDC(0.55, 0.73,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
        //
        //  --- Lower Plot
        kLowerPlot  =   new TPad("kLowerPlot_1D", "kLowerPlot", 0, 0.0, 1, 0.35);
        kLowerPlot  -> SetTopMargin(0);
        kLowerPlot  -> SetBottomMargin(0.3);
        kLowerPlot  -> SetFillColorAlpha( 0., 0. );
        kLowerPlot  -> cd();
        //
        auto    h1D_Nres_divd_nfll_utl  = (TH1F*)h1D_Nres_stat_full->Clone("h1D_Nres_stat_nfll_utl");
        auto    h1D_Nres_stat_nfll_utl  = (TH1F*)h1D_Nres_stat->Clone("h1D_Nres_stat_nfll_utl");
        auto    h1D_Nres_syst_nfll_utl  = (TH1F*)h1D_Nres_syst->Clone("h1D_Nres_syst_nfll_utl");
        auto    h1D_Nres_stat_full_utl  = (TH1F*)h1D_Nres_stat_full->Clone("h1D_Nres_stat_full_utl");
        auto    h1D_Nres_syst_full_utl  = (TH1F*)h1D_Nres_syst_full->Clone("h1D_Nres_syst_full_utl");
        for ( Int_t iBin = 1; iBin <= nBinPT1D; iBin++ ) {
            h1D_Nres_divd_nfll_utl -> SetBinContent (iBin+1,h1D_Nres_stat_nfll_utl->GetBinContent   (iBin));
            h1D_Nres_divd_nfll_utl -> SetBinError   (iBin+1,h1D_Nres_stat_nfll_utl->GetBinError     (iBin));
        }
        h1D_Nres_divd_nfll_utl          = uScale( h1D_Nres_divd_nfll_utl, 1., -2 );
        h1D_Nres_stat_nfll_utl          = uScale( h1D_Nres_stat_nfll_utl, -1. );
        h1D_Nres_syst_nfll_utl          = uScale( h1D_Nres_syst_nfll_utl, -1. );
        h1D_Nres_stat_full_utl          = uScale( h1D_Nres_stat_full_utl, -1. );
        h1D_Nres_syst_full_utl          = uScale( h1D_Nres_syst_full_utl, -1. );
        h1D_Nres_syst_full_utl -> SetMarkerStyle     ( 0 );
        h1D_Nres_syst_full_utl -> GetYaxis() -> SetTitle("Model / Data");
        h1D_Nres_syst_full_utl -> GetYaxis() -> SetTitleSize(0.14);
        h1D_Nres_syst_full_utl -> GetYaxis() -> SetTitleOffset(0.45);
        h1D_Nres_syst_full_utl -> GetYaxis() -> SetLabelSize(0.09);
        h1D_Nres_syst_full_utl -> GetYaxis() -> SetNdivisions(8);
        h1D_Nres_syst_full_utl -> GetXaxis() -> SetLabelSize(0.09);
        h1D_Nres_syst_full_utl -> GetXaxis() -> SetLabelSize(0.09);
        h1D_Nres_syst_full_utl -> GetXaxis() -> SetTitleSize(0.12);
        h1D_Nres_syst_full_utl -> GetXaxis() -> SetTitleOffset(1.);
        h1D_Nres_syst_full_utl -> SetMaximum( 1.9 );
        h1D_Nres_syst_full_utl -> SetMinimum( 0.0 );
        h1D_Nres_syst_full_utl -> Draw("SAME PE2");
        h1D_Nres_stat_full_utl -> Draw("SAME PE X0");
        h1D_Nres_syst_nfll_utl -> Draw("SAME PE2");
        h1D_Nres_stat_nfll_utl -> Draw("SAME PE X0");
        //
        if ( kComparisonMCTagList.size() ) {
            for ( auto kCurrent_1D : h1D_Ntru ) {
                auto kCurrentHist   =   (TH1F*)kCurrent_1D -> Clone(Form("tmp_%s",kCurrent_1D->GetName()));
                kCurrentHist    -> Divide( h1D_Nres_divd_nfll_utl );
                kCurrentHist    -> DrawCopy("SAME E2");
                kCurrentHist    -> ResetAttFill();
                kCurrentHist    -> SetMarkerStyle(0);
                kCurrentHist    -> DrawCopy("SAME HIST");
            }
        }
        //
        cDrawResults    -> cd();
        kUpperPlot      -> Draw();
        kLowerPlot      -> Draw();
        cDrawResults    -> SaveAs( kPlotDirectory + TString("Yield_1D.pdf") );
        cDrawResults    -> SaveAs( kPlotDirectory + TString("Yield_1D.eps") );
        delete cDrawResults;    //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        //
        for ( Int_t iPT2D = 1; iPT2D < h2D_Nres_stat.size(); iPT2D++ ) {
            //
            //  --- --- Conditional 2D Spectrum //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            TH1F*   h2D_Nres_stat_full      = new TH1F( Form("h2D_Nres_stat_full_%i",iPT2D), Form("h2D_Nres_stat_full_%i",iPT2D), nBinPT2D+1, fArrPT2D_Comp );
            TH1F*   h2D_Nres_syst_full      = new TH1F( Form("h2D_Nres_syst_full_%i",iPT2D), Form("h2D_Nres_syst_full_%i",iPT2D), nBinPT2D+1, fArrPT2D_Comp );
            //
            h2D_Nres_stat_full  -> SetBinContent( 1, h2D_Nres_stat.at(0)->GetBinContent(iPT2D) );
            h2D_Nres_syst_full  -> SetBinContent( 1, h2D_Nres_stat.at(0)->GetBinContent(iPT2D) );
            h2D_Nres_stat_full  -> SetBinError  ( 1, h2D_Nres_stat.at(0)->GetBinError  (iPT2D) );
            h2D_Nres_syst_full  -> SetBinError  ( 1, h2D_Nres_syst.at(0)->GetBinError  (iPT2D) );
            h2D_Nres_syst_full  -> SetMaximum( h2D_Nres_syst.at(iPT2D)->GetMaximum()*2.5 );
            h2D_Nres_syst_full  -> SetMinimum( h2D_Nres_syst.at(iPT2D)->GetMinimum()*0.5 );
            uSetHisto( h2D_Nres_stat.at(iPT2D), " SPT STAT 12D " );
            uSetHisto( h2D_Nres_syst.at(iPT2D), " SPT SYST 12D " );
            uSetHisto( h2D_Nres_stat_full,      " SPT STAT 12D VAR1 " );
            uSetHisto( h2D_Nres_syst_full,      " SPT SYST 12D VAR1 " );
            //
            TCanvas*    cDrawResults    = new TCanvas( "cDrawResults", "cDrawResults", 1200, 1000 );
            //
            //  --- Upper Plot
            TPad*   kUpperPlot  =   new TPad("kUpperPlot_2D", "kUpperPlot", 0, 0.35, 1, 1.0);
            gStyle      -> SetOptStat(0);
            kUpperPlot  -> SetTopMargin(0.05);
            kUpperPlot  -> SetBottomMargin(0);
            kUpperPlot  -> SetFillColorAlpha( 0., 0. );
            kUpperPlot  -> SetLogy();
            kUpperPlot  -> cd();
            //
            h2D_Nres_syst_full      -> GetYaxis() -> SetTitleSize(0.06);
            h2D_Nres_syst_full      -> GetYaxis() -> SetTitleOffset(0.98);
            h2D_Nres_syst_full      -> GetYaxis() -> SetLabelSize(0.055);
            h2D_Nres_syst_full      -> Draw("SAME PE2");
            h2D_Nres_stat_full      -> Draw("SAME PE X0");
            h2D_Nres_syst.at(iPT2D) -> Draw("SAME PE2");
            h2D_Nres_stat.at(iPT2D) -> Draw("SAME PE X0");
            lLegend                 -> Draw("SAME");
            if ( kComparisonMCTagList.size() ) lMCLegend_2D            -> Draw("SAME");
            if ( kComparisonMCTagList.size() ) {
                for ( auto kCurrent_2D : h2D_Ntru_Fin_Nrm ) {
                    auto kCurrentHist  =   kCurrent_2D -> ProjectionY("tmp", iPT2D+1 , iPT2D+1 );
                    kCurrentHist    -> SetLineWidth( 2 );
                    kCurrentHist    -> DrawCopy("SAME E3");
                    kCurrentHist    -> SetMarkerStyle(0);
                    kCurrentHist    -> ResetAttFill();
                    kCurrentHist    -> DrawCopy("SAME HIST L");
                }
            }
            //
            uLatex->SetTextFont(60);
            uLatex->SetTextSize(0.05);
            uLatex->DrawLatexNDC(0.55, 0.88,"ALICE Preliminary");
            uLatex->SetTextFont(42);
            uLatex->SetTextSize(0.04);
            uLatex->DrawLatexNDC(0.55, 0.83,"pp #sqrt{#it{s}}= 7 TeV");
            uLatex->DrawLatexNDC(0.55, 0.78,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
            uLatex->DrawLatexNDC(0.55, 0.73,Form("#it{p}_{T,#phi_{2}} [%.2f;%.2f] (GeV/#it{c})",fArrPT2D_Comp[iPT2D],fArrPT2D_Comp[iPT2D+1]));
            //
            //  --- Lower Plot
            TPad*   kLowerPlot  =   new TPad("kLowerPlot_2D", "kLowerPlot", 0, 0.0, 1, 0.35);
            kLowerPlot  -> SetTopMargin(0);
            kLowerPlot  -> SetBottomMargin(0.3);
            kLowerPlot  -> SetFillColorAlpha( 0., 0. );
            kLowerPlot  -> cd();
            //
            auto    h2D_Nres_stat_nfll_utl  = (TH1F*)h2D_Nres_stat.at(iPT2D)->Clone("h2D_Nres_stat_nfll_utl");
            auto    h2D_Nres_syst_nfll_utl  = (TH1F*)h2D_Nres_syst.at(iPT2D)->Clone("h2D_Nres_syst_nfll_utl");
            auto    h2D_Nres_stat_full_utl  = (TH1F*)h2D_Nres_stat_full->Clone("h2D_Nres_stat_full_utl");
            auto    h2D_Nres_syst_full_utl  = (TH1F*)h2D_Nres_syst_full->Clone("h2D_Nres_syst_full_utl");
            //auto    kUtilityPlot            = uScale( h2D_Nres_stat_full_utl, 1., -2 );
            h2D_Nres_stat_nfll_utl          = uScale( h2D_Nres_stat_nfll_utl, -1. );
            h2D_Nres_syst_nfll_utl          = uScale( h2D_Nres_syst_nfll_utl, -1. );
            h2D_Nres_stat_full_utl          = uScale( h2D_Nres_stat_full_utl, -1. );
            h2D_Nres_syst_full_utl          = uScale( h2D_Nres_syst_full_utl, -1. );
            h2D_Nres_syst_full_utl -> SetMarkerStyle     ( 0 );
            h2D_Nres_syst_full_utl -> GetYaxis() -> SetTitle("Model / Data");
            h2D_Nres_syst_full_utl -> GetYaxis() -> SetTitleSize(0.14);
            h2D_Nres_syst_full_utl -> GetYaxis() -> SetTitleOffset(0.45);
            h2D_Nres_syst_full_utl -> GetYaxis() -> SetLabelSize(0.09);
            h2D_Nres_syst_full_utl -> GetYaxis() -> SetNdivisions(8);
            h2D_Nres_syst_full_utl -> GetXaxis() -> SetLabelSize(0.09);
            h2D_Nres_syst_full_utl -> GetXaxis() -> SetLabelSize(0.09);
            h2D_Nres_syst_full_utl -> GetXaxis() -> SetTitleSize(0.12);
            h2D_Nres_syst_full_utl -> GetXaxis() -> SetTitleOffset(1.);
            h2D_Nres_syst_full_utl -> SetMaximum( 1.9 );
            h2D_Nres_syst_full_utl -> SetMinimum( 0.0 );
            h2D_Nres_syst_full_utl -> Draw("SAME PE2");
            h2D_Nres_stat_full_utl -> Draw("SAME PE X0");
            h2D_Nres_syst_nfll_utl -> Draw("SAME PE2");
            h2D_Nres_stat_nfll_utl -> Draw("SAME PE X0");
            //
            if ( kComparisonMCTagList.size() ) {
                for ( auto kCurrent_2D : h2D_Ntru ) {
                    auto kCurrentHist   =   kCurrent_2D -> ProjectionY(Form("tmp_%s",kCurrent_2D->GetName()), iPT2D+1 , iPT2D+1 );
                    //kCurrentHist    -> Divide( kUtilityPlot );
                    kCurrentHist    -> DrawCopy("SAME E2");
                    kCurrentHist    -> ResetAttFill();
                    kCurrentHist    -> DrawCopy("SAME HIST");
                }
                kLowerPlot->SetLogy();
            }
            //
            cDrawResults    -> cd();
            kUpperPlot      -> Draw();
            kLowerPlot      -> Draw();
            cDrawResults    -> SaveAs( kPlotDirectory + TString(Form("Yield_2D_%i.pdf",iPT2D)) );
            cDrawResults    -> SaveAs( kPlotDirectory + TString(Form("Yield_2D_%i.eps",iPT2D)) );
            delete kUpperPlot;
            delete kLowerPlot;
            delete cDrawResults;    //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        }
        delete lLegend;
        //
        //  --- --- Mean pT Spectrum //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        TH1F*   h2D_MeanPT_stat_full    = new TH1F( "h2D_MeanPT_stat_full", "h2D_MeanPT_stat_full", nBinPT2D+1, fArrPT2D_Comp );
        TH1F*   h2D_MeanPT_syst_full    = new TH1F( "h2D_MeanPT_syst_full", "h2D_MeanPT_syst_full", nBinPT2D+1, fArrPT2D_Comp );
        //
        h2D_MeanPT_stat_full  -> SetBinContent( 1, h2D_MeanPT_stat->GetBinContent(1) );
        h2D_MeanPT_syst_full  -> SetBinContent( 1, h2D_MeanPT_stat->GetBinContent(1) );
        h2D_MeanPT_stat_full  -> SetBinError  ( 1, h2D_MeanPT_stat->GetBinError  (1) );
        h2D_MeanPT_syst_full  -> SetBinError  ( 1, h2D_MeanPT_syst->GetBinError  (1) );
        h2D_MeanPT_stat_full -> SetMaximum( h2D_MeanPT_syst->GetMaximum()*1.25 );
        h2D_MeanPT_stat_full -> SetMinimum( h2D_MeanPT_syst->GetMinimum()*0.40 );
        h2D_MeanPT_syst -> GetXaxis() ->SetRangeUser (0.5,100);
        h2D_MeanPT_stat -> GetXaxis() ->SetRangeUser (0.5,100);
        uSetHisto( h2D_MeanPT_stat, " MPT STAT 12D " );
        uSetHisto( h2D_MeanPT_syst, " MPT SYST 12D " );
        uSetHisto( h2D_MeanPT_stat_full,  " MPT STAT 12D VAR1 " );
        uSetHisto( h2D_MeanPT_syst_full,  " MPT SYST 12D VAR1 " );
        //
                cDrawResults    = new TCanvas( "cDrawResults", "cDrawResults", 1200, 1000 );
        //
        //  --- Upper Plot
        kUpperPlot  =   new TPad("kUpperPlot_MPT", "kUpperPlot", 0, 0.35, 1, 1.0);
        gStyle      -> SetOptStat(0);
        kUpperPlot  -> SetBottomMargin(0);
        kUpperPlot  -> SetFillColorAlpha( 0., 0. );
        kUpperPlot  -> cd();
        //
                    lLegend         = new TLegend( 0.600, 0.230, 0.880, 0.430 );
        lLegend     ->  SetLineColorAlpha(0.,0.);
        lLegend     ->  SetFillColorAlpha(0.,0.);
        lLegend     ->  AddEntry    (h2D_MeanPT_stat,       "Measured","P");
        lLegend     ->  AddEntry    (h2D_MeanPT_stat_full,  "Extrapolated","P");
        lLegend     ->  AddEntry    (h2D_MeanPT_stat,       "Stat. Uncert.","EL");
        lLegend     ->  AddEntry    (h2D_MeanPT_syst,       "Syst. Uncert.","F");
        //
        h2D_MeanPT_stat_full    -> GetYaxis() -> SetTitleSize(0.06);
        h2D_MeanPT_stat_full    -> GetYaxis() -> SetTitleOffset(0.98);
        h2D_MeanPT_stat_full    -> GetYaxis() -> SetLabelSize(0.055);
        h2D_MeanPT_stat_full    -> Draw("SAME PE X0");
        h2D_MeanPT_syst_full    -> Draw("SAME PE2");
        h2D_MeanPT_syst         -> Draw("SAME PE2");
        h2D_MeanPT_stat         -> Draw("SAME PE X0");
        lLegend                 -> Draw("SAME");
        if ( kComparisonMCTagList.size() ) lMCLegend_PT            -> Draw("SAME");
        if ( kComparisonMCTagList.size() ) {
            for ( auto kCurrent_MPT : hMPT_NTru )   {
                kCurrent_MPT-> DrawCopy("SAME E3");
                kCurrent_MPT    -> ResetAttFill();
                kCurrent_MPT-> DrawCopy("SAME HIST L");
            }
        }
        //
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.38, 0.15,"ALICE Preliminary");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        uLatex->DrawLatexNDC(0.38, 0.10,"pp #sqrt{#it{s}}= 7 TeV");
        uLatex->DrawLatexNDC(0.38, 0.05,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
        //
        //  --- Lower Plot
        kLowerPlot  =   new TPad("kLowerPlot_MPT", "kLowerPlot", 0, 0.0, 1, 0.35);
        kLowerPlot  -> SetTopMargin(0);
        kLowerPlot  -> SetBottomMargin(0.3);
        kLowerPlot  -> SetFillColorAlpha( 0., 0. );
        kLowerPlot  -> cd();
        //
        auto    h2D_MeanPT_stat_nfll_utl    = (TH1F*)h2D_MeanPT_stat->Clone("h2D_MeanPT_stat_nfll_utl");
        auto    h2D_MeanPT_syst_nfll_utl    = (TH1F*)h2D_MeanPT_syst->Clone("h2D_MeanPT_syst_nfll_utl");
        auto    h2D_MeanPT_stat_full_utl    = (TH1F*)h2D_MeanPT_stat_full->Clone("h2D_MeanPT_stat_full_utl");
        auto    h2D_MeanPT_syst_full_utl    = (TH1F*)h2D_MeanPT_syst_full->Clone("h2D_MeanPT_syst_full_utl");
        h2D_MeanPT_stat_nfll_utl            = uScale( h2D_MeanPT_stat_nfll_utl, -1. );
        h2D_MeanPT_syst_nfll_utl            = uScale( h2D_MeanPT_syst_nfll_utl, -1. );
        h2D_MeanPT_stat_full_utl            = uScale( h2D_MeanPT_stat_full_utl, -1. );
        h2D_MeanPT_syst_full_utl            = uScale( h2D_MeanPT_syst_full_utl, -1. );
        h2D_MeanPT_syst_full_utl -> SetMarkerStyle     ( 0 );
        h2D_MeanPT_syst_full_utl -> GetYaxis() -> SetTitle("Model / Data");
        h2D_MeanPT_syst_full_utl -> GetYaxis() -> SetTitleSize(0.14);
        h2D_MeanPT_syst_full_utl -> GetYaxis() -> SetTitleOffset(0.45);
        h2D_MeanPT_syst_full_utl -> GetYaxis() -> SetLabelSize(0.09);
        h2D_MeanPT_syst_full_utl -> GetYaxis() -> SetNdivisions(8);
        h2D_MeanPT_syst_full_utl -> GetXaxis() -> SetLabelSize(0.09);
        h2D_MeanPT_syst_full_utl -> GetXaxis() -> SetLabelSize(0.09);
        h2D_MeanPT_syst_full_utl -> GetXaxis() -> SetTitleSize(0.12);
        h2D_MeanPT_syst_full_utl -> GetXaxis() -> SetTitleOffset(1.);
        h2D_MeanPT_syst_full_utl -> SetMaximum( 1.1 );
        h2D_MeanPT_syst_full_utl -> SetMinimum( 0.9 );
        h2D_MeanPT_syst_full_utl -> Draw("SAME PE2");
        h2D_MeanPT_stat_full_utl -> Draw("SAME PE X0");
        h2D_MeanPT_syst_nfll_utl -> Draw("SAME PE2");
        h2D_MeanPT_stat_nfll_utl -> Draw("SAME PE X0");
        //
        if ( kComparisonMCTagList.size() ) {
            for ( auto kCurrent_MPT : hMPT_NTru ) {
                auto kCurrentHist   =   (TH1F*)kCurrent_MPT -> Clone(Form("tmp_%s",kCurrent_MPT->GetName()));
                //kCurrentHist    -> Divide( kUtilityPlot );
                kCurrentHist    -> DrawCopy("SAME E2");
                kCurrentHist    -> ResetAttFill();
                kCurrentHist    -> SetMarkerStyle(0);
                kCurrentHist    -> DrawCopy("SAME HIST");
            }
        }
        //
        cDrawResults    -> cd();
        kUpperPlot      -> Draw();
        kLowerPlot      -> Draw();
        cDrawResults    -> SaveAs( kPlotDirectory + TString(Form("MeanPT_2D_%i.pdf",0)) );
        cDrawResults    -> SaveAs( kPlotDirectory + TString(Form("MeanPT_2D_%i.eps",0)) );
        delete kUpperPlot;
        delete kLowerPlot;
        delete cDrawResults; //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        //
        //  --- --- Correlation Results //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        TF1* guasflat = new TF1( "flatgaus", "[0]*exp(-0.5*((x-[1])/[2])**2)+[3]", -500, 500 );
        
        h2D_PhiCorrelat -> Scale(1.e6);
        h2D_PhiCorrelat -> SetMaximum(8.0);
        h2D_PhiCorrelat -> SetMinimum(0.1);
        uSetHisto( h2D_PhiCorrelat, " CRPH " );
        h2D_PhiCorrela2 -> Scale(1.e6);
        h2D_PhiCorrela2 -> SetMaximum(8.0);
        h2D_PhiCorrela2 -> SetMinimum(0.1);
        uSetHisto( h2D_PhiCorrela2, " CRPH " );
        h2D_PhiCorrela2 -> SetLineColor(kBlue);
        h2D_PhiCorrela2 -> SetMarkerColor(kBlue);
        
        guasflat->FixParameter(1,0);
        guasflat->SetParameter(2,40);
        guasflat->SetParameter(3,3);
        h2D_PhiCorrela2 ->  Fit(guasflat, "IMS");
        h2D_PhiCorrelat ->  Fit(guasflat, "IMS");
                cDrawResults    = new TCanvas( "cDrawResults", "cDrawResults", 1200, 1000 );
        //
        //  --- Upper Plot
        kUpperPlot  =   new TPad("kUpperPlot_MPT", "kUpperPlot", 0, 0.35, 1, 1.0);
        gStyle      -> SetOptStat(0);
        kUpperPlot  -> SetBottomMargin(0);
        kUpperPlot  -> SetFillColorAlpha( 0., 0. );
        kUpperPlot  -> cd();
        //
                    lLegend         = new TLegend( 0.500, 0.675, 0.675, 0.575 );
        lLegend     ->  SetLineColorAlpha(0.,0.);
        lLegend     ->  SetFillColorAlpha(0.,0.);
        lLegend     ->  AddEntry    (h2D_PhiCorrelat,       "Measured","P");
        lLegend     ->  AddEntry    (h2D_PhiCorrelat,       "Stat. Uncert.","EL");
        lLegend     ->  AddEntry    (h2D_PhiCorrela2,       "Low InvMass","EL");
        //
        h2D_PhiCorrelat     ->  GetYaxis() -> SetTitleSize(0.06);
        h2D_PhiCorrelat     ->  GetYaxis() -> SetTitleOffset(0.98);
        h2D_PhiCorrelat     ->  GetYaxis() -> SetLabelSize(0.055);
        h2D_PhiCorrelat     ->  Draw("SAME PE");
        h2D_PhiCorrela2     ->  Draw("SAME PE");
        lLegend             ->  Draw("SAME");
        if ( kComparisonMCTagList.size() ) lMCLegend_CR -> Draw("SAME");
        if ( kComparisonMCTagList.size() ) {
            for ( auto kCurrent_CrPh : h2D_PhiCorrelation )   {
                kCurrent_CrPh -> Scale(1.e6);
                kCurrent_CrPh -> DrawCopy("SAME E3");
                kCurrent_CrPh -> ResetAttFill();
                kCurrent_CrPh -> DrawCopy("SAME HIST L");
            }
        }
        //
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.500, 0.825,"ALICE");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        uLatex->DrawLatexNDC(0.500, 0.775,"pp #sqrt{#it{s}}= 5 TeV");
        uLatex->DrawLatexNDC(0.500, 0.725,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
        //
        //  --- Lower Plot
        kLowerPlot  =   new TPad("kLowerPlot_MPT", "kLowerPlot", 0, 0.0, 1, 0.35);
        kLowerPlot  -> SetTopMargin(0);
        kLowerPlot  -> SetBottomMargin(0.3);
        kLowerPlot  -> SetFillColorAlpha( 0., 0. );
        kLowerPlot  -> cd();
        //
        auto    h2D_PhiCorrelat_utl =   (TH1F*)h2D_PhiCorrelat->Clone("h2D_PhiCorrelat_utl");
        auto    kUtilityPlot        =   uScale( h2D_PhiCorrelat_utl, 1., -2 );
        h2D_PhiCorrelat_utl         =   uScale( h2D_PhiCorrelat_utl, -1. );
        h2D_PhiCorrelat_utl         -> SetMaximum( 1.58 );
        h2D_PhiCorrelat_utl         -> SetMinimum( 0.18 );
        uSetHisto( h2D_PhiCorrelat_utl, " CRPH BTM " );
        h2D_PhiCorrelat_utl         ->  Draw("SAME PEX0");
        //
        if ( kComparisonMCTagList.size() ) {
            for ( auto kCurrentPhCr : h2D_PhiCorrelation ) {
                auto kCurrentHist   =   (TH1F*)kCurrentPhCr -> Clone(Form("tmp_%s",kCurrentPhCr->GetName()));
                kCurrentHist    -> Divide( kUtilityPlot );
                kCurrentHist    -> DrawCopy("SAME E2");
                kCurrentHist    -> ResetAttFill();
                kCurrentHist    -> SetMarkerStyle(0);
                kCurrentHist    -> DrawCopy("SAME HIST");
            }
        }
        //
        cDrawResults    -> cd();
        kUpperPlot      -> Draw();
        kLowerPlot      -> Draw();
        cDrawResults    -> SaveAs( kPlotDirectory + TString(Form("PhiCorrelation.pdf")) );
        cDrawResults    -> SaveAs( kPlotDirectory + TString(Form("PhiCorrelation.eps")) );
        delete kUpperPlot;
        delete kLowerPlot;
        delete cDrawResults; //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        //
        std::map<TString,TH1F*> kHistograms;
        auto kValueMap = uCalculateStringEnhancement(0.);
        for ( auto const& [key, val] : kValueMap ) {
            kHistograms[key] = new TH1F( key.Data(), key.Data(), 100, 0., 10. );
            kHistograms[key]->GetXaxis()->SetTitle("Enhancement factor of string tension (h)");
            kHistograms[key]->GetYaxis()->SetTitle("Parameter value");
        }
        for ( Int_t iTer = 0; iTer < 100; iTer++ ) {
            kValueMap = uCalculateStringEnhancement(iTer/10.);
            for ( auto const& [key, val] : kValueMap ) {
                kHistograms[key]->SetBinContent(iTer+1,val);
            }
        }
        TFile* kTEST = new TFile( kPlotDirectory + TString("/test.root"), "RECREATE" );
        for ( auto const& [key, val] : kHistograms ) {
            val->Write();
        }
        kTEST->Close();

        gROOT   ->  SetBatch( kFALSE );
    }
    //
    /*
    // --- YIELD ANALYSIS
    if ( kDoMultiplicity ) {
        //
        //  --- Load Data Files
        TFile*  insFile_Data_YL = new TFile ( Form(kASigExtp_FitCheckRst,(TString("Multiplicity")+kFolder).Data()) );
        TFile*  insFile_Syst_YL = new TFile ( Form("%s/FullSystematics.root",Form(kAnalysis_Systemt_Dir,  (TString("Multiplicity")+kFolder).Data())) );
        TFile*  insFile_Data_CP = new TFile ( Form(kASigExtp_FitCheckRst,(TString("Correlation")+kFolder).Data()) );
        //
        //  --- Load Data Histograms
        auto    g1D_Nres_Mult   = uLoadHistograms<0,TGraphErrors> ( insFile_Data_YL, "g1D_Nres_Mult" );
        auto    g2D_Nres_Mult   = uLoadHistograms<0,TGraphErrors> ( insFile_Data_YL, "g2D_Nres_Mult" );
        auto    gR1_Nres_Mult   = uLoadHistograms<0,TGraphErrors> ( insFile_Data_YL, "gR1_Nres_Mult" );
        auto    gR2_Nres_Mult   = uLoadHistograms<0,TGraphErrors> ( insFile_Data_YL, "gR2_Nres_Mult" );
        auto    gP1_Nres_Mult   = uLoadHistograms<0,TGraphErrors> ( insFile_Data_YL, "gP1_Nres_Mult" );
        auto    gP2_Nres_Mult   = uLoadHistograms<0,TGraphErrors> ( insFile_Data_YL, "gP2_Nres_Mult" );
        //
        // --- MC Elements
        std::vector<TFile*> insFile_MntC_YL;
        std::vector<TH1F*> hMPT_NTru;
        std::vector<TH1F*> hMPT_NTru_Fin;
        std::vector<TH1F*> h1D_Ntru_Fin;
        std::vector<TH2F*> h2D_Ntru_Fin;
        std::vector<TH2F*> h2D_Ntru_Fin_Nrm;
        std::vector<TH1F*> h1D_Ntru;
        std::vector<TH2F*> h2D_Ntru;
        std::vector<TH1F*> hFullQuantities;
        std::vector<TH1F*> hProductionProb;
        std::vector<TH1F*> h2D_PhiCorrelation;
        TLegend*           lMCLegend_1D;
        TLegend*           lMCLegend_2D;
        TLegend*           lMCLegend_PT;
        TLegend*           lMCLegend_PP;
        TLegend*           lMCLegend_CR;
        //
        if ( kComparisonMCTagList.size() ) {
            // --- Load MC Files
            for ( auto kCurrent_MC_File : kComparisonMCTagList ) insFile_MntC_YL.push_back( new TFile   ( Form(kProduction_MC_Ofl,(TString("Multiplicity")+kFolder).Data(),kCurrent_MC_File.Data()) ) );
            //
            // --- Build Legend
            lMCLegend_1D    = new TLegend( 0.500, 0.700, 0.880, 0.550 );
            lMCLegend_2D    = new TLegend( 0.550, 0.650, 0.880, 0.450 );
            lMCLegend_PT    = new TLegend( 0.600, 0.030, 0.880, 0.220 );
            lMCLegend_PP    = new TLegend( 0.175, 0.175, 0.525, 0.375 );
            lMCLegend_CR    = new TLegend( 0.675, 0.675, 0.875, 0.875 );
            lMCLegend_CR    -> SetNColumns(2);
            //
            //  --- Load MC Histograms
            auto iClr = 0;
            for ( auto kFile : insFile_MntC_YL )   {
                iClr++;
                //  --- Production Histograms
                //  --- --- 1D
                h1D_Ntru                    .push_back( (TH1F*)((( kFile )->Get("h1D_Ntru"))        -> Clone(Form("h1D_Ntru_%s",kFile->GetName())) ));
                h1D_Ntru.at( iClr-1 )       -> SetLineWidth( 2 );
                h1D_Ntru.at( iClr-1 )       -> SetLineColor( uGetColor( iClr+1 ) );
                h1D_Ntru.at( iClr-1 )       -> SetFillColor( uGetFillColor( iClr+1 ) );
                h1D_Ntru_Fin                .push_back( (TH1F*)((( kFile )->Get("h1D_Ntru_Fin"))    -> Clone(Form("h1D_Ntru_Fin_%s",kFile->GetName())) ));
                h1D_Ntru_Fin.at( iClr-1 )   -> SetLineWidth( 2 );
                h1D_Ntru_Fin.at( iClr-1 )   -> SetLineColor( uGetColor( iClr+1 ) );
                h1D_Ntru_Fin.at( iClr-1 )   -> SetFillColor( uGetFillColor( iClr+1 ) );
                //
                //  --- --- 2D
                h2D_Ntru                    .push_back( (TH2F*)((( kFile )->Get("h2D_Ntru"))        -> Clone(Form("h2D_Ntru_%s",kFile->GetName())) ));
                h2D_Ntru.at( iClr-1 )       -> SetLineWidth( 2 );
                h2D_Ntru.at( iClr-1 )       -> SetLineColor( uGetColor( iClr+1 ) );
                h2D_Ntru.at( iClr-1 )       -> SetFillColorAlpha( uGetColor( iClr+1 ), 0.33 );
                h2D_Ntru_Fin                .push_back( (TH2F*)((( kFile )->Get("h2D_Ntru_Fin"))    -> Clone(Form("h2D_Ntru_Fin_%s",kFile->GetName())) ));
                h2D_Ntru_Fin.at( iClr-1 )   -> SetLineWidth( 2 );
                h2D_Ntru_Fin.at( iClr-1 )   -> SetLineColor( uGetColor( iClr+1 ) );
                h2D_Ntru_Fin.at( iClr-1 )   -> SetFillColor( uGetFillColor( iClr+1 ) );
                h2D_Ntru_Fin_Nrm            .push_back( (TH2F*)((( kFile )->Get("h2D_Ntru_Fin_Nrm"))-> Clone(Form("h2D_Ntru_Fin_Nrm_%s",kFile->GetName())) ));
                h2D_Ntru_Fin_Nrm.at( iClr-1 )-> SetLineWidth( 2 );
                h2D_Ntru_Fin_Nrm.at( iClr-1 )-> SetLineColor( uGetColor( iClr+1 ) );
                h2D_Ntru_Fin_Nrm.at( iClr-1 )-> SetFillColor( uGetFillColor( iClr+1 ) );
                //
                //  --- --- MPT
                hMPT_NTru                   .push_back( (TH1F*)((( kFile )->Get("hMPT_NTru"))       -> Clone(Form("hMPT_NTru_%s",kFile->GetName())) ));
                hMPT_NTru.at( iClr-1 )      -> SetLineWidth( 2 );
                hMPT_NTru.at( iClr-1 )      -> SetLineColor( uGetColor( iClr+1 ) );
                hMPT_NTru.at( iClr-1 )      -> SetFillColor( uGetFillColor( iClr+1 ) );
                hMPT_NTru_Fin               .push_back( (TH1F*)((( kFile )->Get("hMPT_NTru_Fin"))   -> Clone(Form("hMPT_NTru_Fin_%s",kFile->GetName())) ));
                hMPT_NTru_Fin.at( iClr-1 )  -> SetLineWidth( 2 );
                hMPT_NTru_Fin.at( iClr-1 )  -> SetLineColor( uGetColor( iClr+1 ) );
                hMPT_NTru_Fin.at( iClr-1 )  -> SetFillColor( uGetFillColor( iClr+1 ) );
                //
                //  --- --- Other
                hFullQuantities             .push_back( (TH1F*)((( kFile )->Get("hFullQuantities")) -> Clone(Form("hFullQuantities_%s",kFile->GetName())) ));
                hFullQuantities.at( iClr-1 )-> SetLineWidth( 2 );
                hFullQuantities.at( iClr-1 )-> SetLineColor( uGetColor( iClr+1 ) );
                hFullQuantities.at( iClr-1 )-> SetFillColor( uGetFillColor( iClr+1 ) );
                //
                hProductionProb             .push_back( (TH1F*)((( kFile )->Get("hProductionProb")) -> Clone(Form("hProductionProb_%s",kFile->GetName())) ));
                hProductionProb.at( iClr-1 )-> SetLineWidth( 2 );
                hProductionProb.at( iClr-1 )-> SetLineColor( uGetColor( iClr+1 ) );
                hProductionProb.at( iClr-1 )-> SetFillColor( uGetFillColor( iClr+1 ) );
                //
                h2D_PhiCorrelation             .push_back( (TH1F*)((( kFile )->Get("hPhiCorrelation")) -> Clone(Form("hPhiCorrelation_%s",kFile->GetName())) ));
                h2D_PhiCorrelation.at( iClr-1 )-> SetLineWidth( 2 );
                h2D_PhiCorrelation.at( iClr-1 )-> SetLineColor( uGetColor( iClr+1 ) );
                h2D_PhiCorrelation.at( iClr-1 )-> SetFillColor( uGetFillColor( iClr+1 ) );
                //
                lMCLegend_1D                -> AddEntry( hMPT_NTru.at( iClr-1 ), kComparisonLegendList.at( iClr-1 ), "L" );
                lMCLegend_2D                -> AddEntry( hMPT_NTru.at( iClr-1 ), kComparisonLegendList.at( iClr-1 ), "L" );
                lMCLegend_PT                -> AddEntry( hMPT_NTru.at( iClr-1 ), kComparisonLegendList.at( iClr-1 ), "L" );
                lMCLegend_PP                -> AddEntry( hMPT_NTru.at( iClr-1 ), kComparisonLegendList.at( iClr-1 ), "L" );
                lMCLegend_CR                -> AddEntry( hMPT_NTru.at( iClr-1 ), kComparisonLegendList.at( iClr-1 ), "L" );
            }
        }
        //
        //  --- Output directory
        TString kPlotDirectory          =   Form(kDIR_ShowPlots,(TString("Multiplicity")+kFolder).Data());
        gROOT   ->  ProcessLine(Form(".! mkdir -p %s",kPlotDirectory.Data()));
        //
        //  --- Yield Plots
        //  --- --- Batch Mode
        gROOT   ->  SetBatch( kTRUE );
        //
        //  --- --- Production Probability
        TCanvas*    cDrawResults    = new TCanvas( "cDrawResults", "cDrawResults", 1200, 1000 );
        gPad                    -> SetLogy();
        //
        TH1F* hUtility = new TH1F("","", 5, 0., +5 );
        hUtility->SetMaximum(1.);
        hUtility->SetMinimum(1.e-4);
        hUtility->GetXaxis()->SetNdivisions(5);
        hUtility->GetXaxis()->SetTitle("Number of #phi mesons produced");
        hUtility->GetYaxis()->SetTitle("Event frequency");
        hUtility->Draw();
        if ( kComparisonMCTagList.size() ) {
            for ( auto kCurrent_PP : hProductionProb ) {
                kCurrent_PP-> SetMarkerStyle(0);
                kCurrent_PP-> ResetAttFill();
                kCurrent_PP-> DrawCopy("SAME HIST L");
            }
        }
        if ( kComparisonMCTagList.size() ) lMCLegend_PP    ->  Draw("SAME");
        //
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.55, 0.83,"ALICE Simulation");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        uLatex->DrawLatexNDC(0.68, 0.78,"pp #sqrt{#it{s}}= 7 TeV");
        uLatex->DrawLatexNDC(0.74, 0.73,"#phi, |#it{y}|<0.5");
        //
        cDrawResults    -> SaveAs( kPlotDirectory + TString("hProductionProb.pdf") );
        cDrawResults    -> SaveAs( kPlotDirectory + TString("hProductionProb.eps") );
        delete cDrawResults;
        //
        //  --- --- Final quantities
        uSetHisto( hXD_Nfqs_stat,       "FNL STAT 1D" );
        uSetHisto( hXD_Nfqs_syst,       "FNL SYST 1D" );
        //
        TLegend*    lLegend         = new TLegend( 0.43, 0.03, 0.62, 0.25 );
        lLegend     ->  SetLineColorAlpha(0.,0.);
        lLegend     ->  SetFillColorAlpha(0.,0.);
        lLegend     ->  AddEntry    (hXD_Nfqs_stat,"Measured","P");
        lLegend     ->  AddEntry    (hXD_Nfqs_stat,"Stat. Uncert.","EL");
        lLegend     ->  AddEntry    (hXD_Nfqs_syst,"Syst. Uncert.","F");
        //
                cDrawResults    = new TCanvas( "cDrawResults", "cDrawResults", 1200, 1000 );
        gPad                    -> SetLogy();
        //
        //  --- Upper Plot
        TPad*   kUpperPlot  =   new TPad("kUpperPlot_RT", "kUpperPlot", 0, 0.40, 1, 1.0);
        gStyle      -> SetOptStat(0);
        kUpperPlot  -> SetTopMargin(0.05);
        kUpperPlot  -> SetBottomMargin(0);
        kUpperPlot  -> SetFillColorAlpha( 0., 0. );
        kUpperPlot  -> SetLogy();
        kUpperPlot  -> cd();
        //
        hXD_Nfqs_syst   -> GetYaxis() -> SetTitleSize(0.06);
        hXD_Nfqs_syst   -> GetYaxis() -> SetTitleOffset(0.98);
        hXD_Nfqs_syst   -> GetYaxis() -> SetLabelSize(0.055);
        hXD_Nfqs_syst   -> Draw("SAME PE2");
        hXD_Nfqs_stat   -> Draw("SAME PE X0");
        if ( kComparisonMCTagList.size() ) {
            for ( auto kCurrent_FQ : hFullQuantities )  {
                kCurrent_FQ-> SetMarkerStyle(0);
                kCurrent_FQ-> ResetAttFill();
                kCurrent_FQ-> Draw("SAME EP");
            }
        }
        lLegend                 -> Draw("SAME");
        if ( kComparisonMCTagList.size() ) lMCLegend_PT            -> Draw("SAME");
        //
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.07);
        uLatex->DrawLatexNDC(0.18, 0.86,"ALICE Preliminary");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.18, 0.80,"pp #sqrt{#it{s}}= 7 TeV");
        uLatex->DrawLatexNDC(0.18, 0.74,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
        //
        //  --- Lower Plot
        TPad*   kLowerPlot  =   new TPad("kLowerPlot_RT", "kLowerPlot", 0, 0.0, 1, 0.40);
        kLowerPlot  -> SetTopMargin(0);
        kLowerPlot  -> SetBottomMargin(0.3);
        kLowerPlot  -> SetFillColorAlpha( 0., 0. );
        kLowerPlot  -> cd();
        //
        auto kUtility_syst  = uScale( hXD_Nfqs_syst, -1 );
        kUtility_syst -> SetName("tmp_syst_1D");
        uSetHisto( kUtility_syst, " FNL SYST 12D " );
        kUtility_syst -> SetFillColor   ( uGetFillColor(1) );
        kUtility_syst -> SetMarkerStyle     ( 0 );
        kUtility_syst -> GetYaxis() -> SetTitle("Model / Data");
        kUtility_syst -> GetYaxis() -> SetTitleSize(0.12);
        kUtility_syst -> GetYaxis() -> SetTitleOffset(0.55);
        kUtility_syst -> GetYaxis() -> SetLabelSize(0.09);
        kUtility_syst -> GetYaxis() -> SetNdivisions(7);
        kUtility_syst -> GetXaxis() -> SetLabelSize(0.14);
        kUtility_syst -> GetXaxis() -> SetLabelOffset(0.0115);
        kUtility_syst -> GetXaxis() -> SetTitleSize(0.01);
        kUtility_syst -> GetXaxis() -> SetTitleOffset(1.0);
        kUtility_syst -> SetMaximum( 1.7 );
        kUtility_syst -> SetMinimum( 0.3 );
        kUtility_syst -> Draw("SAME PE2");
        //
        auto kUtility_stat = uScale( hXD_Nfqs_stat, -1 );
        kUtility_syst -> SetName("tmp_stat_1D");
        kUtility_stat -> Draw("SAME PE X0");
        //
        if ( kComparisonMCTagList.size() ) {
            auto kUtilityPlot       = uScale( hXD_Nfqs_stat, 1., -2.);
            for ( auto kCurrent_FQ : hFullQuantities ) {
                auto kCurrentHist   =   (TH1F*)kCurrent_FQ -> Clone(Form("tmp_%s",kCurrent_FQ->GetName()));
                kCurrentHist    -> Divide( kUtilityPlot );
                kCurrentHist    -> DrawCopy("SAME EP");
            }
        }
        //
        cDrawResults    -> cd();
        kUpperPlot      -> Draw();
        kLowerPlot      -> Draw();
        cDrawResults    -> SaveAs( kPlotDirectory + TString("FullQuantitites.pdf") );
        cDrawResults    -> SaveAs( kPlotDirectory + TString("FullQuantitites.eps") );
        delete lLegend;
        delete kUpperPlot;
        delete kLowerPlot;
        delete cDrawResults;
        //
        //  --- --- Full 1D Spectrum //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        TH1F*   h1D_Nres_stat_full      = new TH1F( "h1D_Nres_stat_full", "h1D_Nres_stat_full", nBinPT1D+1, fArrPT1D_Comp );
        TH1F*   h1D_Nres_syst_full      = new TH1F( "h1D_Nres_syst_full", "h1D_Nres_syst_full", nBinPT1D+1, fArrPT1D_Comp );
        //
        h1D_Nres_stat_full  -> SetBinContent( 1, hXD_Nyld_stat->GetBinContent(4) / (fArrPT1D_Comp[1]-fArrPT1D_Comp[0]) );
        h1D_Nres_syst_full  -> SetBinContent( 1, hXD_Nyld_stat->GetBinContent(4) / (fArrPT1D_Comp[1]-fArrPT1D_Comp[0]) );
        h1D_Nres_stat_full  -> SetBinError  ( 1, hXD_Nyld_stat->GetBinError  (4) / (fArrPT1D_Comp[1]-fArrPT1D_Comp[0]) );
        h1D_Nres_syst_full  -> SetBinError  ( 1, hXD_Nyld_syst->GetBinError  (4) / (fArrPT1D_Comp[1]-fArrPT1D_Comp[0]) );
        h1D_Nres_syst_full  -> SetMaximum   ( h1D_Nres_syst->GetMaximum()*2.5 );
        h1D_Nres_syst_full  -> SetMinimum   ( h1D_Nres_syst->GetMinimum()*0.5 );
        uSetHisto( h1D_Nres_stat,       "SPT STAT 1D" );
        uSetHisto( h1D_Nres_syst,       "SPT SYST 1D" );
        uSetHisto( h1D_Nres_stat_full,  "SPT STAT 1D VAR1 " );
        uSetHisto( h1D_Nres_syst_full,  "SPT SYST 1D VAR1 " );
        //
        lLegend     = new TLegend( 0.20, 0.03, 0.38, 0.25 );
        lLegend     ->  SetLineColorAlpha(0.,0.);
        lLegend     ->  SetFillColorAlpha(0.,0.);
        lLegend     ->  AddEntry    (h1D_Nres_stat,"Measured","P");
        lLegend     ->  AddEntry    (h1D_Nres_stat_full,"Extrapolated","P");
        lLegend     ->  AddEntry    (h1D_Nres_stat,"Stat. Uncert.","EL");
        lLegend     ->  AddEntry    (h1D_Nres_syst,"Syst. Uncert.","F");
        //
        cDrawResults    = new TCanvas( "cDrawResults", "cDrawResults", 1200, 1000 );
        //
        //  --- Upper Plot
        kUpperPlot  =   new TPad("kUpperPlot_1D", "kUpperPlot", 0, 0.35, 1, 1.0);
        gStyle      -> SetOptStat(0);
        kUpperPlot  -> SetTopMargin(0.05);
        kUpperPlot  -> SetBottomMargin(0);
        kUpperPlot  -> SetFillColorAlpha( 0., 0. );
        kUpperPlot  -> SetLogy();
        kUpperPlot  -> cd();
        //
        h1D_Nres_syst_full  -> GetYaxis() -> SetTitleSize(0.06);
        h1D_Nres_syst_full  -> GetYaxis() -> SetTitleOffset(0.98);
        h1D_Nres_syst_full  -> GetYaxis() -> SetLabelSize(0.055);
        h1D_Nres_syst_full  -> Draw("SAME PE2");
        h1D_Nres_stat_full  -> Draw("SAME PE X0");
        h1D_Nres_syst       -> Draw("SAME PE2");
        h1D_Nres_stat       -> Draw("SAME PE X0");
        lLegend             -> Draw("SAME");
        if ( kComparisonMCTagList.size() ) lMCLegend_1D        -> Draw("SAME");
        if ( kComparisonMCTagList.size() ) {
            for ( auto kCurrent_1D : h1D_Ntru_Fin ) {
                kCurrent_1D-> DrawCopy("SAME E3");
                kCurrent_1D-> SetFillColor( 0. );
                kCurrent_1D-> ResetAttFill();
                kCurrent_1D-> DrawCopy("SAME HIST L");
            }
        }
        //
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.55, 0.83,"ALICE Preliminary");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        uLatex->DrawLatexNDC(0.55, 0.78,"pp #sqrt{#it{s}}= 7 TeV");
        uLatex->DrawLatexNDC(0.55, 0.73,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
        //
        //  --- Lower Plot
        kLowerPlot  =   new TPad("kLowerPlot_1D", "kLowerPlot", 0, 0.0, 1, 0.35);
        kLowerPlot  -> SetTopMargin(0);
        kLowerPlot  -> SetBottomMargin(0.3);
        kLowerPlot  -> SetFillColorAlpha( 0., 0. );
        kLowerPlot  -> cd();
        //
        auto    h1D_Nres_divd_nfll_utl  = (TH1F*)h1D_Nres_stat_full->Clone("h1D_Nres_stat_nfll_utl");
        auto    h1D_Nres_stat_nfll_utl  = (TH1F*)h1D_Nres_stat->Clone("h1D_Nres_stat_nfll_utl");
        auto    h1D_Nres_syst_nfll_utl  = (TH1F*)h1D_Nres_syst->Clone("h1D_Nres_syst_nfll_utl");
        auto    h1D_Nres_stat_full_utl  = (TH1F*)h1D_Nres_stat_full->Clone("h1D_Nres_stat_full_utl");
        auto    h1D_Nres_syst_full_utl  = (TH1F*)h1D_Nres_syst_full->Clone("h1D_Nres_syst_full_utl");
        for ( Int_t iBin = 1; iBin <= nBinPT1D; iBin++ ) {
            h1D_Nres_divd_nfll_utl -> SetBinContent (iBin+1,h1D_Nres_stat_nfll_utl->GetBinContent   (iBin));
            h1D_Nres_divd_nfll_utl -> SetBinError   (iBin+1,h1D_Nres_stat_nfll_utl->GetBinError     (iBin));
        }
        h1D_Nres_divd_nfll_utl          = uScale( h1D_Nres_divd_nfll_utl, 1., -2 );
        h1D_Nres_stat_nfll_utl          = uScale( h1D_Nres_stat_nfll_utl, -1. );
        h1D_Nres_syst_nfll_utl          = uScale( h1D_Nres_syst_nfll_utl, -1. );
        h1D_Nres_stat_full_utl          = uScale( h1D_Nres_stat_full_utl, -1. );
        h1D_Nres_syst_full_utl          = uScale( h1D_Nres_syst_full_utl, -1. );
        h1D_Nres_syst_full_utl -> SetMarkerStyle     ( 0 );
        h1D_Nres_syst_full_utl -> GetYaxis() -> SetTitle("Model / Data");
        h1D_Nres_syst_full_utl -> GetYaxis() -> SetTitleSize(0.14);
        h1D_Nres_syst_full_utl -> GetYaxis() -> SetTitleOffset(0.45);
        h1D_Nres_syst_full_utl -> GetYaxis() -> SetLabelSize(0.09);
        h1D_Nres_syst_full_utl -> GetYaxis() -> SetNdivisions(8);
        h1D_Nres_syst_full_utl -> GetXaxis() -> SetLabelSize(0.09);
        h1D_Nres_syst_full_utl -> GetXaxis() -> SetLabelSize(0.09);
        h1D_Nres_syst_full_utl -> GetXaxis() -> SetTitleSize(0.12);
        h1D_Nres_syst_full_utl -> GetXaxis() -> SetTitleOffset(1.);
        h1D_Nres_syst_full_utl -> SetMaximum( 1.9 );
        h1D_Nres_syst_full_utl -> SetMinimum( 0.0 );
        h1D_Nres_syst_full_utl -> Draw("SAME PE2");
        h1D_Nres_stat_full_utl -> Draw("SAME PE X0");
        h1D_Nres_syst_nfll_utl -> Draw("SAME PE2");
        h1D_Nres_stat_nfll_utl -> Draw("SAME PE X0");
        //
        if ( kComparisonMCTagList.size() ) {
            for ( auto kCurrent_1D : h1D_Ntru ) {
                auto kCurrentHist   =   (TH1F*)kCurrent_1D -> Clone(Form("tmp_%s",kCurrent_1D->GetName()));
                kCurrentHist    -> Divide( h1D_Nres_divd_nfll_utl );
                kCurrentHist    -> DrawCopy("SAME E2");
                kCurrentHist    -> ResetAttFill();
                kCurrentHist    -> SetMarkerStyle(0);
                kCurrentHist    -> DrawCopy("SAME HIST");
            }
        }
        //
        cDrawResults    -> cd();
        kUpperPlot      -> Draw();
        kLowerPlot      -> Draw();
        cDrawResults    -> SaveAs( kPlotDirectory + TString("Yield_1D.pdf") );
        cDrawResults    -> SaveAs( kPlotDirectory + TString("Yield_1D.eps") );
        delete cDrawResults;    //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        //
        for ( Int_t iPT2D = 1; iPT2D < h2D_Nres_stat.size(); iPT2D++ ) {
            //
            //  --- --- Conditional 2D Spectrum //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            TH1F*   h2D_Nres_stat_full      = new TH1F( Form("h2D_Nres_stat_full_%i",iPT2D), Form("h2D_Nres_stat_full_%i",iPT2D), nBinPT2D+1, fArrPT2D_Comp );
            TH1F*   h2D_Nres_syst_full      = new TH1F( Form("h2D_Nres_syst_full_%i",iPT2D), Form("h2D_Nres_syst_full_%i",iPT2D), nBinPT2D+1, fArrPT2D_Comp );
            //
            h2D_Nres_stat_full  -> SetBinContent( 1, h2D_Nres_stat.at(0)->GetBinContent(iPT2D) );
            h2D_Nres_syst_full  -> SetBinContent( 1, h2D_Nres_stat.at(0)->GetBinContent(iPT2D) );
            h2D_Nres_stat_full  -> SetBinError  ( 1, h2D_Nres_stat.at(0)->GetBinError  (iPT2D) );
            h2D_Nres_syst_full  -> SetBinError  ( 1, h2D_Nres_syst.at(0)->GetBinError  (iPT2D) );
            h2D_Nres_syst_full  -> SetMaximum( h2D_Nres_syst.at(iPT2D)->GetMaximum()*2.5 );
            h2D_Nres_syst_full  -> SetMinimum( h2D_Nres_syst.at(iPT2D)->GetMinimum()*0.5 );
            uSetHisto( h2D_Nres_stat.at(iPT2D), " SPT STAT 12D " );
            uSetHisto( h2D_Nres_syst.at(iPT2D), " SPT SYST 12D " );
            uSetHisto( h2D_Nres_stat_full,      " SPT STAT 12D VAR1 " );
            uSetHisto( h2D_Nres_syst_full,      " SPT SYST 12D VAR1 " );
            //
            TCanvas*    cDrawResults    = new TCanvas( "cDrawResults", "cDrawResults", 1200, 1000 );
            //
            //  --- Upper Plot
            TPad*   kUpperPlot  =   new TPad("kUpperPlot_2D", "kUpperPlot", 0, 0.35, 1, 1.0);
            gStyle      -> SetOptStat(0);
            kUpperPlot  -> SetTopMargin(0.05);
            kUpperPlot  -> SetBottomMargin(0);
            kUpperPlot  -> SetFillColorAlpha( 0., 0. );
            kUpperPlot  -> SetLogy();
            kUpperPlot  -> cd();
            //
            h2D_Nres_syst_full      -> GetYaxis() -> SetTitleSize(0.06);
            h2D_Nres_syst_full      -> GetYaxis() -> SetTitleOffset(0.98);
            h2D_Nres_syst_full      -> GetYaxis() -> SetLabelSize(0.055);
            h2D_Nres_syst_full      -> Draw("SAME PE2");
            h2D_Nres_stat_full      -> Draw("SAME PE X0");
            h2D_Nres_syst.at(iPT2D) -> Draw("SAME PE2");
            h2D_Nres_stat.at(iPT2D) -> Draw("SAME PE X0");
            lLegend                 -> Draw("SAME");
            if ( kComparisonMCTagList.size() ) lMCLegend_2D            -> Draw("SAME");
            if ( kComparisonMCTagList.size() ) {
                for ( auto kCurrent_2D : h2D_Ntru_Fin_Nrm ) {
                    auto kCurrentHist  =   kCurrent_2D -> ProjectionY("tmp", iPT2D+1 , iPT2D+1 );
                    kCurrentHist    -> SetLineWidth( 2 );
                    kCurrentHist    -> DrawCopy("SAME E3");
                    kCurrentHist    -> SetMarkerStyle(0);
                    kCurrentHist    -> ResetAttFill();
                    kCurrentHist    -> DrawCopy("SAME HIST L");
                }
            }
            //
            uLatex->SetTextFont(60);
            uLatex->SetTextSize(0.05);
            uLatex->DrawLatexNDC(0.55, 0.88,"ALICE Preliminary");
            uLatex->SetTextFont(42);
            uLatex->SetTextSize(0.04);
            uLatex->DrawLatexNDC(0.55, 0.83,"pp #sqrt{#it{s}}= 7 TeV");
            uLatex->DrawLatexNDC(0.55, 0.78,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
            uLatex->DrawLatexNDC(0.55, 0.73,Form("#it{p}_{T,#phi_{2}} [%.2f;%.2f] (GeV/#it{c})",fArrPT2D_Comp[iPT2D],fArrPT2D_Comp[iPT2D+1]));
            //
            //  --- Lower Plot
            TPad*   kLowerPlot  =   new TPad("kLowerPlot_2D", "kLowerPlot", 0, 0.0, 1, 0.35);
            kLowerPlot  -> SetTopMargin(0);
            kLowerPlot  -> SetBottomMargin(0.3);
            kLowerPlot  -> SetFillColorAlpha( 0., 0. );
            kLowerPlot  -> cd();
            //
            auto    h2D_Nres_stat_nfll_utl  = (TH1F*)h2D_Nres_stat.at(iPT2D)->Clone("h2D_Nres_stat_nfll_utl");
            auto    h2D_Nres_syst_nfll_utl  = (TH1F*)h2D_Nres_syst.at(iPT2D)->Clone("h2D_Nres_syst_nfll_utl");
            auto    h2D_Nres_stat_full_utl  = (TH1F*)h2D_Nres_stat_full->Clone("h2D_Nres_stat_full_utl");
            auto    h2D_Nres_syst_full_utl  = (TH1F*)h2D_Nres_syst_full->Clone("h2D_Nres_syst_full_utl");
            //auto    kUtilityPlot            = uScale( h2D_Nres_stat_full_utl, 1., -2 );
            h2D_Nres_stat_nfll_utl          = uScale( h2D_Nres_stat_nfll_utl, -1. );
            h2D_Nres_syst_nfll_utl          = uScale( h2D_Nres_syst_nfll_utl, -1. );
            h2D_Nres_stat_full_utl          = uScale( h2D_Nres_stat_full_utl, -1. );
            h2D_Nres_syst_full_utl          = uScale( h2D_Nres_syst_full_utl, -1. );
            h2D_Nres_syst_full_utl -> SetMarkerStyle     ( 0 );
            h2D_Nres_syst_full_utl -> GetYaxis() -> SetTitle("Model / Data");
            h2D_Nres_syst_full_utl -> GetYaxis() -> SetTitleSize(0.14);
            h2D_Nres_syst_full_utl -> GetYaxis() -> SetTitleOffset(0.45);
            h2D_Nres_syst_full_utl -> GetYaxis() -> SetLabelSize(0.09);
            h2D_Nres_syst_full_utl -> GetYaxis() -> SetNdivisions(8);
            h2D_Nres_syst_full_utl -> GetXaxis() -> SetLabelSize(0.09);
            h2D_Nres_syst_full_utl -> GetXaxis() -> SetLabelSize(0.09);
            h2D_Nres_syst_full_utl -> GetXaxis() -> SetTitleSize(0.12);
            h2D_Nres_syst_full_utl -> GetXaxis() -> SetTitleOffset(1.);
            h2D_Nres_syst_full_utl -> SetMaximum( 1.9 );
            h2D_Nres_syst_full_utl -> SetMinimum( 0.0 );
            h2D_Nres_syst_full_utl -> Draw("SAME PE2");
            h2D_Nres_stat_full_utl -> Draw("SAME PE X0");
            h2D_Nres_syst_nfll_utl -> Draw("SAME PE2");
            h2D_Nres_stat_nfll_utl -> Draw("SAME PE X0");
            //
            if ( kComparisonMCTagList.size() ) {
                for ( auto kCurrent_2D : h2D_Ntru ) {
                    auto kCurrentHist   =   kCurrent_2D -> ProjectionY(Form("tmp_%s",kCurrent_2D->GetName()), iPT2D+1 , iPT2D+1 );
                    //kCurrentHist    -> Divide( kUtilityPlot );
                    kCurrentHist    -> DrawCopy("SAME E2");
                    kCurrentHist    -> ResetAttFill();
                    kCurrentHist    -> DrawCopy("SAME HIST");
                }
                kLowerPlot->SetLogy();
            }
            //
            cDrawResults    -> cd();
            kUpperPlot      -> Draw();
            kLowerPlot      -> Draw();
            cDrawResults    -> SaveAs( kPlotDirectory + TString(Form("Yield_2D_%i.pdf",iPT2D)) );
            cDrawResults    -> SaveAs( kPlotDirectory + TString(Form("Yield_2D_%i.eps",iPT2D)) );
            delete kUpperPlot;
            delete kLowerPlot;
            delete cDrawResults;    //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        }
        delete lLegend;
        //
        //  --- --- Mean pT Spectrum //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        TH1F*   h2D_MeanPT_stat_full    = new TH1F( "h2D_MeanPT_stat_full", "h2D_MeanPT_stat_full", nBinPT2D+1, fArrPT2D_Comp );
        TH1F*   h2D_MeanPT_syst_full    = new TH1F( "h2D_MeanPT_syst_full", "h2D_MeanPT_syst_full", nBinPT2D+1, fArrPT2D_Comp );
        //
        h2D_MeanPT_stat_full  -> SetBinContent( 1, h2D_MeanPT_stat->GetBinContent(1) );
        h2D_MeanPT_syst_full  -> SetBinContent( 1, h2D_MeanPT_stat->GetBinContent(1) );
        h2D_MeanPT_stat_full  -> SetBinError  ( 1, h2D_MeanPT_stat->GetBinError  (1) );
        h2D_MeanPT_syst_full  -> SetBinError  ( 1, h2D_MeanPT_syst->GetBinError  (1) );
        h2D_MeanPT_stat_full -> SetMaximum( h2D_MeanPT_syst->GetMaximum()*1.25 );
        h2D_MeanPT_stat_full -> SetMinimum( h2D_MeanPT_syst->GetMinimum()*0.40 );
        h2D_MeanPT_syst -> GetXaxis() ->SetRangeUser (0.5,100);
        h2D_MeanPT_stat -> GetXaxis() ->SetRangeUser (0.5,100);
        uSetHisto( h2D_MeanPT_stat, " MPT STAT 12D " );
        uSetHisto( h2D_MeanPT_syst, " MPT SYST 12D " );
        uSetHisto( h2D_MeanPT_stat_full,  " MPT STAT 12D VAR1 " );
        uSetHisto( h2D_MeanPT_syst_full,  " MPT SYST 12D VAR1 " );
        //
                cDrawResults    = new TCanvas( "cDrawResults", "cDrawResults", 1200, 1000 );
        //
        //  --- Upper Plot
        kUpperPlot  =   new TPad("kUpperPlot_MPT", "kUpperPlot", 0, 0.35, 1, 1.0);
        gStyle      -> SetOptStat(0);
        kUpperPlot  -> SetBottomMargin(0);
        kUpperPlot  -> SetFillColorAlpha( 0., 0. );
        kUpperPlot  -> cd();
        //
                    lLegend         = new TLegend( 0.600, 0.230, 0.880, 0.430 );
        lLegend     ->  SetLineColorAlpha(0.,0.);
        lLegend     ->  SetFillColorAlpha(0.,0.);
        lLegend     ->  AddEntry    (h2D_MeanPT_stat,       "Measured","P");
        lLegend     ->  AddEntry    (h2D_MeanPT_stat_full,  "Extrapolated","P");
        lLegend     ->  AddEntry    (h2D_MeanPT_stat,       "Stat. Uncert.","EL");
        lLegend     ->  AddEntry    (h2D_MeanPT_syst,       "Syst. Uncert.","F");
        //
        h2D_MeanPT_stat_full    -> GetYaxis() -> SetTitleSize(0.06);
        h2D_MeanPT_stat_full    -> GetYaxis() -> SetTitleOffset(0.98);
        h2D_MeanPT_stat_full    -> GetYaxis() -> SetLabelSize(0.055);
        h2D_MeanPT_stat_full    -> Draw("SAME PE X0");
        h2D_MeanPT_syst_full    -> Draw("SAME PE2");
        h2D_MeanPT_syst         -> Draw("SAME PE2");
        h2D_MeanPT_stat         -> Draw("SAME PE X0");
        lLegend                 -> Draw("SAME");
        if ( kComparisonMCTagList.size() ) lMCLegend_PT            -> Draw("SAME");
        if ( kComparisonMCTagList.size() ) {
            for ( auto kCurrent_MPT : hMPT_NTru )   {
                kCurrent_MPT-> DrawCopy("SAME E3");
                kCurrent_MPT    -> ResetAttFill();
                kCurrent_MPT-> DrawCopy("SAME HIST L");
            }
        }
        //
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.38, 0.15,"ALICE Preliminary");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        uLatex->DrawLatexNDC(0.38, 0.10,"pp #sqrt{#it{s}}= 7 TeV");
        uLatex->DrawLatexNDC(0.38, 0.05,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
        //
        //  --- Lower Plot
        kLowerPlot  =   new TPad("kLowerPlot_MPT", "kLowerPlot", 0, 0.0, 1, 0.35);
        kLowerPlot  -> SetTopMargin(0);
        kLowerPlot  -> SetBottomMargin(0.3);
        kLowerPlot  -> SetFillColorAlpha( 0., 0. );
        kLowerPlot  -> cd();
        //
        auto    h2D_MeanPT_stat_nfll_utl    = (TH1F*)h2D_MeanPT_stat->Clone("h2D_MeanPT_stat_nfll_utl");
        auto    h2D_MeanPT_syst_nfll_utl    = (TH1F*)h2D_MeanPT_syst->Clone("h2D_MeanPT_syst_nfll_utl");
        auto    h2D_MeanPT_stat_full_utl    = (TH1F*)h2D_MeanPT_stat_full->Clone("h2D_MeanPT_stat_full_utl");
        auto    h2D_MeanPT_syst_full_utl    = (TH1F*)h2D_MeanPT_syst_full->Clone("h2D_MeanPT_syst_full_utl");
        h2D_MeanPT_stat_nfll_utl            = uScale( h2D_MeanPT_stat_nfll_utl, -1. );
        h2D_MeanPT_syst_nfll_utl            = uScale( h2D_MeanPT_syst_nfll_utl, -1. );
        h2D_MeanPT_stat_full_utl            = uScale( h2D_MeanPT_stat_full_utl, -1. );
        h2D_MeanPT_syst_full_utl            = uScale( h2D_MeanPT_syst_full_utl, -1. );
        h2D_MeanPT_syst_full_utl -> SetMarkerStyle     ( 0 );
        h2D_MeanPT_syst_full_utl -> GetYaxis() -> SetTitle("Model / Data");
        h2D_MeanPT_syst_full_utl -> GetYaxis() -> SetTitleSize(0.14);
        h2D_MeanPT_syst_full_utl -> GetYaxis() -> SetTitleOffset(0.45);
        h2D_MeanPT_syst_full_utl -> GetYaxis() -> SetLabelSize(0.09);
        h2D_MeanPT_syst_full_utl -> GetYaxis() -> SetNdivisions(8);
        h2D_MeanPT_syst_full_utl -> GetXaxis() -> SetLabelSize(0.09);
        h2D_MeanPT_syst_full_utl -> GetXaxis() -> SetLabelSize(0.09);
        h2D_MeanPT_syst_full_utl -> GetXaxis() -> SetTitleSize(0.12);
        h2D_MeanPT_syst_full_utl -> GetXaxis() -> SetTitleOffset(1.);
        h2D_MeanPT_syst_full_utl -> SetMaximum( 1.1 );
        h2D_MeanPT_syst_full_utl -> SetMinimum( 0.9 );
        h2D_MeanPT_syst_full_utl -> Draw("SAME PE2");
        h2D_MeanPT_stat_full_utl -> Draw("SAME PE X0");
        h2D_MeanPT_syst_nfll_utl -> Draw("SAME PE2");
        h2D_MeanPT_stat_nfll_utl -> Draw("SAME PE X0");
        //
        if ( kComparisonMCTagList.size() ) {
            for ( auto kCurrent_MPT : hMPT_NTru ) {
                auto kCurrentHist   =   (TH1F*)kCurrent_MPT -> Clone(Form("tmp_%s",kCurrent_MPT->GetName()));
                //kCurrentHist    -> Divide( kUtilityPlot );
                kCurrentHist    -> DrawCopy("SAME E2");
                kCurrentHist    -> ResetAttFill();
                kCurrentHist    -> SetMarkerStyle(0);
                kCurrentHist    -> DrawCopy("SAME HIST");
            }
        }
        //
        cDrawResults    -> cd();
        kUpperPlot      -> Draw();
        kLowerPlot      -> Draw();
        cDrawResults    -> SaveAs( kPlotDirectory + TString(Form("MeanPT_2D_%i.pdf",0)) );
        cDrawResults    -> SaveAs( kPlotDirectory + TString(Form("MeanPT_2D_%i.eps",0)) );
        delete kUpperPlot;
        delete kLowerPlot;
        delete cDrawResults; //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        //
        //  --- --- Correlation Results //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        h2D_PhiCorrelat -> Scale(1.e6);
        h2D_PhiCorrelat -> SetMaximum(8.0);
        h2D_PhiCorrelat -> SetMinimum(0.1);
        uSetHisto( h2D_PhiCorrelat, " CRPH " );
                cDrawResults    = new TCanvas( "cDrawResults", "cDrawResults", 1200, 1000 );
        //
        //  --- Upper Plot
        kUpperPlot  =   new TPad("kUpperPlot_MPT", "kUpperPlot", 0, 0.35, 1, 1.0);
        gStyle      -> SetOptStat(0);
        kUpperPlot  -> SetBottomMargin(0);
        kUpperPlot  -> SetFillColorAlpha( 0., 0. );
        kUpperPlot  -> cd();
        //
                    lLegend         = new TLegend( 0.500, 0.675, 0.675, 0.575 );
        lLegend     ->  SetLineColorAlpha(0.,0.);
        lLegend     ->  SetFillColorAlpha(0.,0.);
        lLegend     ->  AddEntry    (h2D_PhiCorrelat,       "Measured","P");
        lLegend     ->  AddEntry    (h2D_PhiCorrelat,       "Stat. Uncert.","EL");
        //
        h2D_PhiCorrelat     ->  GetYaxis() -> SetTitleSize(0.06);
        h2D_PhiCorrelat     ->  GetYaxis() -> SetTitleOffset(0.98);
        h2D_PhiCorrelat     ->  GetYaxis() -> SetLabelSize(0.055);
        h2D_PhiCorrelat     ->  Draw("SAME PE");
        lLegend             ->  Draw("SAME");
        if ( kComparisonMCTagList.size() ) lMCLegend_CR -> Draw("SAME");
        if ( kComparisonMCTagList.size() ) {
            for ( auto kCurrent_CrPh : h2D_PhiCorrelation )   {
                kCurrent_CrPh -> Scale(1.e6);
                kCurrent_CrPh -> DrawCopy("SAME E3");
                kCurrent_CrPh -> ResetAttFill();
                kCurrent_CrPh -> DrawCopy("SAME HIST L");
            }
        }
        //
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.500, 0.825,"ALICE");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        uLatex->DrawLatexNDC(0.500, 0.775,"pp #sqrt{#it{s}}= 7 TeV");
        uLatex->DrawLatexNDC(0.500, 0.725,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
        //
        //  --- Lower Plot
        kLowerPlot  =   new TPad("kLowerPlot_MPT", "kLowerPlot", 0, 0.0, 1, 0.35);
        kLowerPlot  -> SetTopMargin(0);
        kLowerPlot  -> SetBottomMargin(0.3);
        kLowerPlot  -> SetFillColorAlpha( 0., 0. );
        kLowerPlot  -> cd();
        //
        auto    h2D_PhiCorrelat_utl =   (TH1F*)h2D_PhiCorrelat->Clone("h2D_PhiCorrelat_utl");
        auto    kUtilityPlot        =   uScale( h2D_PhiCorrelat_utl, 1., -2 );
        h2D_PhiCorrelat_utl         =   uScale( h2D_PhiCorrelat_utl, -1. );
        h2D_PhiCorrelat_utl         -> SetMaximum( 1.58 );
        h2D_PhiCorrelat_utl         -> SetMinimum( 0.18 );
        uSetHisto( h2D_PhiCorrelat_utl, " CRPH BTM " );
        h2D_PhiCorrelat_utl         ->  Draw("SAME PEX0");
        //
        if ( kComparisonMCTagList.size() ) {
            for ( auto kCurrentPhCr : h2D_PhiCorrelation ) {
                auto kCurrentHist   =   (TH1F*)kCurrentPhCr -> Clone(Form("tmp_%s",kCurrentPhCr->GetName()));
                kCurrentHist    -> Divide( kUtilityPlot );
                kCurrentHist    -> DrawCopy("SAME E2");
                kCurrentHist    -> ResetAttFill();
                kCurrentHist    -> SetMarkerStyle(0);
                kCurrentHist    -> DrawCopy("SAME HIST");
            }
        }
        //
        cDrawResults    -> cd();
        kUpperPlot      -> Draw();
        kLowerPlot      -> Draw();
        cDrawResults    -> SaveAs( kPlotDirectory + TString(Form("PhiCorrelation.pdf")) );
        cDrawResults    -> SaveAs( kPlotDirectory + TString(Form("PhiCorrelation.eps")) );
        delete kUpperPlot;
        delete kLowerPlot;
        delete cDrawResults; //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        //
        std::map<TString,TH1F*> kHistograms;
        auto kValueMap = uCalculateStringEnhancement(0.);
        for ( auto const& [key, val] : kValueMap ) {
            kHistograms[key] = new TH1F( key.Data(), key.Data(), 100, 0., 10. );
            kHistograms[key]->GetXaxis()->SetTitle("Enhancement factor of string tension (h)");
            kHistograms[key]->GetYaxis()->SetTitle("Parameter value");
        }
        for ( Int_t iTer = 0; iTer < 100; iTer++ ) {
            kValueMap = uCalculateStringEnhancement(iTer/10.);
            for ( auto const& [key, val] : kValueMap ) {
                kHistograms[key]->SetBinContent(iTer+1,val);
            }
        }
        TFile* kTEST = new TFile( kPlotDirectory + TString("/test.root"), "RECREATE" );
        for ( auto const& [key, val] : kHistograms ) {
            val->Write();
        }
        kTEST->Close();
        
        gROOT   ->  SetBatch( kFALSE );
    }*/
}

