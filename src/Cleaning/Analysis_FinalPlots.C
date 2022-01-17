#include "../../inc/AliAnalysisPhiPair.h"
// !TODO: [INFO] About trees in input

void Analysis_FinalPlots ()
{
    //---------------------//
    //  Setting up input   //
    //---------------------//
    
    // Retrieving PreProcessed data histograms
    TFile  *insFile_DT_Yield            =   new TFile   (fYldSigCorr);
    TFile  *insFile_DT_Mult             =   new TFile   (fMltSigCorr);
    //TFile **fEfficiency                 =   new TFile   [nPeriodFiles];
    
    // Previous Publications checks
    TFile*  insFile_PrevPap             =   new TFile   ("./result/HEPData-ins1762364-v1-Table_4.root");
    
    // Recovering the histograms-------------------------------------------------------------------------------

    // >-> YIELD ANALYSIS //
    //
    TGraphAsymmErrors                  *g1D_Res_Stat,      *g1D_Res_Syst;
    TGraphAsymmErrors                 **g2D_Res_Stat,     **g2D_Res_Syst;
    TGraphAsymmErrors                  *fYield_Stat,       *fYield_Syst;
    //
    hName                               =   Form("gRES_1D_Stat");
    g1D_Res_Stat                        =   (TGraphAsymmErrors*)(insFile_DT_Yield->Get(hName));
    //
    hName                               =   Form("gRES_1D_Syst");
    g1D_Res_Syst                        =   (TGraphAsymmErrors*)(insFile_DT_Yield->Get(hName));
    //
    g2D_Res_Stat                        =   new TGraphAsymmErrors*  [nBinPT2D];
    g2D_Res_Syst                        =   new TGraphAsymmErrors*  [nBinPT2D];
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ ) {
        hName                           =   Form("gRES_2D_Stat_%i",iPT2D);
        g2D_Res_Stat[iPT2D]             =   (TGraphAsymmErrors*)(insFile_DT_Yield->Get(hName));
        hName                           =   Form("gRES_2D_Syst_%i",iPT2D);
        g2D_Res_Syst[iPT2D]             =   (TGraphAsymmErrors*)(insFile_DT_Yield->Get(hName));
    }
    //
    hName                               =   Form("fYield_Stat");
    fYield_Stat                         =   (TGraphAsymmErrors*)(insFile_DT_Yield->Get(hName));
    //
    hName                               =   Form("fYield_Syst");
    fYield_Syst                         =   (TGraphAsymmErrors*)(insFile_DT_Yield->Get(hName));
    
    
    
    /*
    // >-> MULTIPLICITY ANALYSIS //
    //
    TGraphAsymmErrors **fYield_Stat_in_MT  =   new TGraphAsymmErrors  *[nBinMult];
    TGraphAsymmErrors **fYield_Syst_in_MT  =   new TGraphAsymmErrors  *[nBinMult];
    //
    for ( Int_t iMult = 0; iMult < nBinMult; iMult++ )   {
        hName                               =   Form("fYield_Stat_%i",iMult);
        fYield_Stat_in_MT[iMult]            =   (TGraphAsymmErrors*)(insFile_DT_Mult->Get(hName));
        //
        hName                               =   Form("fYield_Syst_%i",iMult);
        fYield_Syst_in_MT[iMult]            =   (TGraphAsymmErrors*)(insFile_DT_Mult->Get(hName));
    }
     */
    //
    //---------------------//
    //  Setting up output  //
    //---------------------//
    
    // Generating the binning array--------------------------------------------------------------------------
    fSetBinPT1D();
    fSetBinIM1D();
    fSetBinPT2D();
    fSetBinIM2D();
    fSetBinRap_();
    fSetBinMult();
    fSetBinNTup();
    Int_t       U_AccCand[1024];
    Int_t       U_nAccept;
    
    // Creating the histograms-------------------------------------------------------------------------------
    //
    //
    //-------------------------//
    //  Filling output objects //
    //-------------------------//
    //
    //  >>->> YIELD ANALYSIS //
    //
    TFile   *fCheck =   new TFile("eeeee.root","recreate");
    TCanvas    *cInclusiveSpectrum      =   fPublishSpectrum(g1D_Res_Stat,g1D_Res_Syst);
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ ) {
        TCanvas     *cConditionalSpectrum   =    fPublishSpectrum(g2D_Res_Stat[iPT2D],g2D_Res_Syst[iPT2D]);
    }
    fCheck->Close();
}
    /*
    //  >>->>->> Recovering Basic values
    Double_t    fY_phi          =   fYield_Stat->GetPointY(0);
    Double_t    fY_phi_stat     =   fYield_Stat->GetErrorY(0);
    Double_t    fY_phi_syst     =   fYield_Syst->GetErrorY(0);
    Double_t    fY_phiphi       =   fYield_Stat->GetPointY(1);
    Double_t    fY_phiphi_stat  =   fYield_Stat->GetErrorY(1);
    Double_t    fY_phiphi_syst  =   fYield_Syst->GetErrorY(1);
    //
    //  >>->>->> Calculating GammaPhi
    //
    fGammaPhi_Stat->SetPoint        (0, 1,  fGammaPhiValue(fY_phi,fY_phiphi));
    fGammaPhi_Stat->SetPointError   (0, 0,  0,  fGammaPhiError(fY_phi,fY_phiphi,fY_phi_stat,fY_phiphi_stat),    fGammaPhiError(fY_phi,fY_phiphi,fY_phi_stat,fY_phiphi_stat));
    fGammaPhi_Syst->SetPoint        (0, 1,  fGammaPhiValue(fY_phi,fY_phiphi));
    fGammaPhi_Syst->SetPointError   (0, 0,  0,  fGammaPhiError(fY_phi,fY_phiphi,fY_phi_stat,fY_phiphi_stat),    fGammaPhiError(fY_phi,fY_phiphi,fY_phi_stat,fY_phiphi_stat));
    //
     */
    
    
    //  >>->> MULTIPLICITY ANALYSIS //
    //
    /*
    for ( Int_t iMult = 0; iMult < nBinMult; iMult++ )   {
        //  >>->>->> Recovering Basic values
        Double_t    fY_phi_MT          =   fYield_Stat_in_MT[iMult]->GetPointY(0);
        Double_t    fY_phi_stat_MT     =   fYield_Stat_in_MT[iMult]->GetErrorY(0);
        Double_t    fY_phi_syst_MT     =   fYield_Stat_in_MT[iMult]->GetErrorY(0);
        Double_t    fY_phiphi_MT       =   fYield_Stat_in_MT[iMult]->GetPointY(1);
        Double_t    fY_phiphi_stat_MT  =   fYield_Stat_in_MT[iMult]->GetErrorY(1);
        Double_t    fY_phiphi_syst_MT  =   fYield_Stat_in_MT[iMult]->GetErrorY(1);
        //
        //  >>->>->> Calculating GammaPhi
        //
        fGammaPhi_Stat_in_MT->SetPoint(iMult,iMult,fGammaPhiValue(fY_phi_MT,fY_phiphi_MT));
        fGammaPhi_Stat_in_MT->SetPointError(iMult,0,0,fGammaPhiError(fY_phi_MT,fY_phiphi_MT,fY_phi_stat_MT,fY_phiphi_stat_MT),fGammaPhiError(fY_phi_MT,fY_phiphi_MT,fY_phi_stat_MT,fY_phiphi_stat_MT));
        fGammaPhi_Syst_in_MT->SetPoint(iMult,iMult,fGammaPhiValue(fY_phi_MT,fY_phiphi_MT));
        fGammaPhi_Syst_in_MT->SetPointError(iMult,0,0,fGammaPhiError(fY_phi_MT,fY_phiphi_MT,fY_phi_stat_MT,fY_phiphi_stat_MT),fGammaPhiError(fY_phi_MT,fY_phiphi_MT,fY_phi_stat_MT,fY_phiphi_stat_MT));
        //
    }
    //
    std::vector<TH1D*>  fVector_Full_2D_Conditional_Spectra;
    //fVector_Full_2D_Conditional_Spectra.push_back(Final);
    for ( Int_t iPT2D = 1; iPT2D <= nBinPT2D; iPT2D++ ) {
        TH1D    *   fConditionalYield   =   hRES_2D->ProjectionX(Form("hRES_2D_PorjX_%i",iPT2D),iPT2D+1,iPT2D+1);
        if ( fConditionalYield->GetEntries() <= 0 ) continue;
        fVector_Full_2D_Conditional_Spectra.push_back(fConditionalYield);
    }*/
    //
    //-------------------------//
    //  Analysis               //
    //-------------------------//
    //
    /*
    TCanvas        *fGammaPhi       =   new TCanvas();
    TLegend        *fLegendGPhi     =   new TLegend();
    TMultiGraph    *fGammaPhiFull   =   new TMultiGraph();
    //
    fGammaPhi_Syst  ->  GetHistogram()  ->  SetMaximum(0.1);
    fGammaPhi_Syst  ->  GetHistogram()  ->  SetMinimum(-0.05);
    fGammaPhi_Stat  ->  SetMarkerStyle(33);
    fGammaPhi_Stat  ->  SetMarkerColor(2);
    fGammaPhi_Stat  ->  SetMarkerSize(2);
    fGammaPhi_Stat  ->  SetDrawOption("EP");
    
    fGammaPhi_Syst  ->  SetMarkerStyle(33);
    fGammaPhi_Syst  ->  SetMarkerColor(2);
    fGammaPhi_Syst  ->  SetMarkerSize(2);
    fGammaPhi_Syst  ->  SetDrawOption("EP");
    //
    fGammaPhiFull   ->Add           (fGammaPhi_Syst);
    fGammaPhiFull   ->Add           (fGammaPhi_Stat);
    //
    fLegendGPhi     ->AddEntry      (fGammaPhi_Syst,    "Data",         "P");
    fLegendGPhi     ->AddEntry      (fGammaPhi_Syst,    "Syst. Err.",   "F");
    fLegendGPhi     ->AddEntry      (fGammaPhi_Stat,    "Stat. Err.",   "E");
    //
    fGammaPhiFull   ->Draw      ("ALP");
    fLegendGPhi     ->Draw      ("SAME");
    fGammaPhi       ->SaveAs    ("fGammaPhi.pdf");
    //
    //--
    //
    TCanvas        *fGammaPhi_in_MT =   new TCanvas();
    TLegend        *fLegendGPhiMT   =   new TLegend();
    TMultiGraph    *fGammaPhiFull_in_MT   =   new TMultiGraph();
    //
    fGammaPhi_Syst_in_MT    ->  GetHistogram()  ->  SetMaximum(0.1);
    fGammaPhi_Syst_in_MT    ->  GetHistogram()  ->  SetMinimum(-0.05);
    fGammaPhi_Stat_in_MT    ->  SetMarkerStyle(33);
    fGammaPhi_Stat_in_MT    ->  SetMarkerColor(2);
    fGammaPhi_Stat_in_MT    ->  SetMarkerSize(2);
    fGammaPhi_Stat_in_MT    ->  SetDrawOption("EP");
    
    fGammaPhi_Syst_in_MT    ->  SetMarkerStyle(33);
    fGammaPhi_Syst_in_MT    ->  SetMarkerColor(2);
    fGammaPhi_Syst_in_MT    ->  SetMarkerSize(2);
    fGammaPhi_Syst_in_MT    ->  SetDrawOption("EP");
    //
    fGammaPhiFull_in_MT     ->Add           (fGammaPhi_Syst_in_MT);
    fGammaPhiFull_in_MT     ->Add           (fGammaPhi_Stat_in_MT);
    //
    fLegendGPhiMT           ->AddEntry      (fGammaPhi_Syst,    "Data",         "P");
    fLegendGPhiMT           ->AddEntry      (fGammaPhi_Syst,    "Syst. Err.",   "F");
    fLegendGPhiMT           ->AddEntry      (fGammaPhi_Stat,    "Stat. Err.",   "E");
    //
    fGammaPhiFull_in_MT     ->Draw      ("ALP");
    fLegendGPhiMT           ->Draw      ("SAME");
    fGammaPhi_in_MT         ->SaveAs    ("fGammaPhi_in_MT.pdf");
    //
    //--
    //
    TCanvas        *fGammaPhi_All   =   new TCanvas();
    TLegend        *fLegendGPhiAll  =   new TLegend();
    
    fLegendGPhiAll          ->Draw      ("SAME");
    fGammaPhi_All           ->SaveAs    ("fGammaPhi_All.pdf");
    //
    TCanvas        *fCanvas_Conditional_Spectra    =   new TCanvas();
    TLegend        *fLegend_Conditional_Spectra    =   new TLegend(0.7,0.4,0.9,0.9);
    gPad                    ->SetLogy();
    gStyle                  ->SetOptStat(0);
    
    Int_t   fIterator   =   3;
    fSetRainbowStyle(fVector_Full_2D_Conditional_Spectra);
    for ( TH1D*&    fConditionalYield   :   fVector_Full_2D_Conditional_Spectra )  {
        fLegend_Conditional_Spectra ->AddEntry(fConditionalYield,   Form("#phi p_{T} #in [%2.1f-%2.1f]",fArrPT2D[fIterator-1],fArrPT2D[fIterator]), "EP");
        fConditionalYield->Draw("SAME EP");
        fConditionalYield->Draw("SAME HIST L");
        fIterator++;
    }
    
    fLegend_Conditional_Spectra ->Draw      ("SAME");
    fCanvas_Conditional_Spectra ->SaveAs    ("fCanvas_Conditional_Spectra.pdf");
    //
    //fCheckPublishedResults(hRES_1D,hTRU,gTRU)->Draw();*/

/*

//-// Recovering Data histograms

//------ 1D Histograms Recovery ------//

hName                           =   "h1D_Raw";                              // Name of the histogram in the preprocessed file
TH1F *  h1D_Raw                 =   (TH1F*)(insFile_DT->Get(hName));
hName                           =   "hNP_1D_Eff_PT_S";                      // Name of the histogram in the preprocessed file
TH1F *  h1D_Eff                 =   (TH1F*)(insFile_EF->Get(hName));
hName                           =   "Entry_DT";                             // Name of the histogram in the preprocessed file
TH1F *  h1D_EDT                 =   (TH1F*)(insFile_DT->Get(hName));
hName                           =   "hNP_1D_Tru_PT_S";                      // Name of the histogram in the preprocessed file
TH1F *  h1D_Tru                 =   (TH1F*)(insFile_EF->Get(hName));
hName                           =   "hNP_1D_Rec_PT_S";                      // Name of the histogram in the preprocessed file
TH1F *  h1D_Rec                 =   (TH1F*)(insFile_EF->Get(hName));
hName                           =   "hNP_1D_Tru_PT_S";                      // Name of the histogram in the preprocessed file
TH1F *  h1D_Tru_P6              =   (TH1F*)(insFile_MC_P6->Get(hName));
hName                           =   "hNP_1D_Tru_PT_S";                      // Name of the histogram in the preprocessed file
TH1F *  h1D_Tru_P8              =   (TH1F*)(insFile_MC_P8->Get(hName));
hName                           =   "hNP_1D_Rec_PT_S";                      // Name of the histogram in the preprocessed file
TH1F *  h1D_Rec_P6              =   (TH1F*)(insFile_MC_P6->Get(hName));
hName                           =   "hNP_1D_Rec_PT_S";                      // Name of the histogram in the preprocessed file
TH1F *  h1D_Rec_P8              =   (TH1F*)(insFile_MC_P8->Get(hName));

//------ 2D Histograms Recovery ------//

hName                           =   "h2D_Raw";                              // Name of the histogram in the preprocessed file
TH2F *  h2D_Raw                 =   (TH2F*)(insFile_DT->Get(hName));
hName                           =   "hNP_2D_Eff_X2_S";                      // Name of the histogram in the preprocessed file
TH2F *  h2D_Eff                 =   (TH2F*)(insFile_EF->Get(hName));
hName                           =   "hNP_2D_Tru_PT_S";                      // Name of the histogram in the preprocessed file
TH2F *  h2D_Tru                 =   (TH2F*)(insFile_EF->Get(hName));
hName                           =   "hNP_2D_Rec_PT_S";                      // Name of the histogram in the preprocessed file
TH2F *  h2D_Rec                 =   (TH2F*)(insFile_EF->Get(hName));
hName                           =   "hNP_2D_Tru_PT_S";                      // Name of the histogram in the preprocessed file
TH2F *  h2D_Tru_P6              =   (TH2F*)(insFile_MC_P6->Get(hName));
hName                           =   "hNP_2D_Tru_PT_S";                      // Name of the histogram in the preprocessed file
TH2F *  h2D_Tru_P8              =   (TH2F*)(insFile_MC_P8->Get(hName));
hName                           =   "hNP_2D_Rec_PT_S";                      // Name of the histogram in the preprocessed file
TH2F *  h2D_Rec_P6              =   (TH2F*)(insFile_MC_P6->Get(hName));
hName                           =   "hNP_2D_Rec_PT_S";                      // Name of the histogram in the preprocessed file
TH2F *  h2D_Rec_P8              =   (TH2F*)(insFile_MC_P8->Get(hName));

//------ TGraphAsymmErrors Rec ------//
hName                           =   "Table 2/Graph1D_y1";                   // Name of the histogram in the preprocessed file
TGraphAsymmErrors *  gCheck_    =   (TGraphAsymmErrors*)(insCheck_->Get(hName));

//---------------------//
//  Setting up output  //
//---------------------//

// Creating the histograms-------------------------------------------------------------------------------

//--- Generating the binning array ---//
fSetBinPT1D();
fSetBinIM1D();
fSetBinPT2D();
fSetBinIM2D();

//---------- 1D Histograms ----------//

hName   = "h1D_Res";
hTitle  = "Number of #phi in |y|<5";
TH1F*   h1D_Res     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
TH1_SetDescription (h1D_Res);

//---------- 2D Histograms ----------//

hName   = "h2D_Res";
hTitle  = "Number of #phi in |y|<5";
TH2F*   h2D_Res     =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
TH2_SetDescription (h2D_Res);
*/


/*------------*/
/*  ANALYSIS  */
/*------------*/
/*
// Output File for Fit Check
TFile*  outFile_FT  =   new TFile(fAnlResHist,"recreate");

//------- Histograms Scaling  -------//

//         N_raw    1     1    Eff    1     1     1     1     1
// N_res = ----- X --- X --- X --- X --- X --- X --- X --- X ---
//          eff    DpT   Dy    Evn   BR    Vrt   Sig   Trk   PID

// Error propagation
h1D_Res->Sumw2();
h2D_Res->Sumw2();

// Producing of N_res by applying corrections to N_raw
h1D_Res->Divide(h1D_Raw,h1D_Eff,1.,1.,"");
h2D_Res->Divide(h2D_Raw,h2D_Eff,1.,1.,"");

// Scaling in pT
h1D_Raw->Scale(1.,"width");
h2D_Raw->Scale(1.,"width");
h1D_Res->Scale(1.,"width");
h2D_Res->Scale(1.,"width");

// Scaling in Rapidity Interval
h1D_Res->Scale(1./kRapidityInterval);
h2D_Res->Scale(1./kRapidityInterval);

// Scaling for events corrected for efficiency
h1D_Raw->Scale(kEventEfficiency_/(h1D_EDT->GetBinContent(1)));
h2D_Raw->Scale(kEventEfficiency_/(h1D_EDT->GetBinContent(1)));
h1D_Res->Scale(kEventEfficiency_/(h1D_EDT->GetBinContent(1)));
h2D_Res->Scale(kEventEfficiency_/(h1D_EDT->GetBinContent(1)));

// Scaling for Branching Ratio
h1D_Res->Scale(1./kBranchingRatio__);
h2D_Res->Scale(1./(kBranchingRatio__*kBranchingRatio__));

if ( bPythiaTest )
{
    // Scaling for MC biases
    h1D_Raw->Scale(1./(kEventEfficiency_*kPythia1DEfficien));
    h2D_Raw->Scale(1./(kEventEfficiency_*kPythia2DEfficien));
    h1D_Res->Scale(1./(kEventEfficiency_*kPythia1DEfficien));
    h2D_Res->Scale(1./(kEventEfficiency_*kPythia2DEfficien));
}
else
{
    // Scaling for Vertex Selection Efficiency
    h1D_Res->Scale(1./kVertexEfficiency);
    h2D_Res->Scale(1./kVertexEfficiency);
}



//------- Histograms Fitting  -------//

SetLevyTsalis();
gCheck_->Fit(fLevyFit1D,"IMRE0S","",0.4,6.0);
cout << "M:" << fLevyFit1D->Mean(0.,10.) << endl;
TCanvas * c133 = new TCanvas();
//gPad->SetLogy();
//gCheck_->Draw();
//fLevyFit1D->Draw("same");
h1D_Res->Divide(fLevyFit1D);
h1D_Res->Draw();
c133->SaveAs("./graphs/eeee.pdf");

// 1D PT Eval
Double_t *  fResults1D = new Double_t [6];
gSystem->RedirectOutput("/dev/null");
fResults1D = ExtrapolateVl( h1D_Res );
gSystem->RedirectOutput(0,0);

// 2D PT Eval
Double_t***     fResults2D =   new Double_t**  [nBinPT2D+1];
for ( Int_t iTer = 0; iTer < nBinPT2D+1; iTer++ )
{
    fResults2D[iTer]        =   new Double_t*   [2];
    fResults2D[iTer][0]     =   new Double_t    [6];
    fResults2D[iTer][1]     =   new Double_t    [6];
}
gSystem->RedirectOutput("/dev/null");
fResults2D = ExtrapolateVl( h2D_Res );
gSystem->RedirectOutput(0,0);

Double_t errREC, errRAW;
Double_t intREC = h1D_Rec->IntegralAndError(5,100,errREC,"width");
Double_t intRAW = h1D_Raw->IntegralAndError(5,100,errRAW,"width");
cout << "1D:" << endl;
cout << "INT_Rec: " << intREC               << "+-" << errREC                                       << endl;
cout << "INT_Raw: " << intRAW               << "+-" << errRAW                                       << endl;
cout << "Raw/Rec: " << intRAW/intREC        << "+-" << (1/intREC)*(errRAW+errREC*(intRAW/intREC))   << endl;
cout << endl;
intREC = h1D_Tru->IntegralAndError(5,100,errREC,"width");
intRAW = h1D_Res->IntegralAndError(5,100,errRAW,"width");
cout << "1D:" << endl;
cout << "INT_Tru: " << intREC               << "+-" << errREC                                       << endl;
cout << "INT_Res: " << intRAW               << "+-" << errRAW                                       << endl;
cout << "Res/Tru: " << intRAW/intREC        << "+-" << (1/intREC)*(errRAW+errREC*(intRAW/intREC))   << endl;
cout << endl;
intREC = h1D_Tru->IntegralAndError(-1,100,errREC,"width");
cout << "1D:" << endl;
cout << "INT_Tru: " << intREC               << "+-" << errREC                                                   << endl;
cout << "INT_Res: " << fResults1D[0]        << "+-" << fResults1D[1] << " +-" << fResults1D[2]                  << endl;
cout << "Res/Tru: " << fResults1D[0]/intREC << "+-" << (1/intREC)*(fResults1D[1]+errREC*(fResults1D[0]/intREC)) << endl;
cout << endl;
cout << endl;
intREC = h2D_Rec->IntegralAndError(5,100,5,100,errREC,"width");
intRAW = h2D_Raw->IntegralAndError(5,100,5,100,errRAW,"width");
cout << "2D:" << endl;
cout << "INT_Rec: " << intREC               << "+-" << errREC                                       << endl;
cout << "INT_Raw: " << intRAW               << "+-" << errRAW                                       << endl;
cout << "Raw/Rec: " << intRAW/intREC        << "+-" << (1/intREC)*(errRAW+errREC*(intRAW/intREC))   << endl;
cout << endl;
intREC = h2D_Tru->IntegralAndError(5,100,5,100,errREC,"width");
intRAW = h2D_Res->IntegralAndError(5,100,5,100,errRAW,"width");
cout << "2D:" << endl;
cout << "INT_Tru: " << intREC               << "+-" << errREC                                       << endl;
cout << "INT_Res: " << intRAW               << "+-" << errRAW                                       << endl;
cout << "Res/Tru: " << intRAW/intREC        << "+-" << (1/intREC)*(errRAW+errREC*(intRAW/intREC))   << endl;
cout << endl;
intREC = h2D_Tru->IntegralAndError(-1,100,-1,100,errREC,"width");
cout << "2DX:" << endl;
cout << "INT_Tru: " << intREC                       << "+-" << errREC                                                               << endl;
cout << "INT_Res: " << fResults2D[0][0][0]          << "+-" << fResults2D[0][0][1]  << " +-" << fResults2D[0][0][2]                 << endl;
cout << "Res/Tru: " << fResults2D[0][0][0]/intREC   << "+-" << (1/intREC)*(fResults2D[0][0][1]+errREC*(fResults2D[0][0][0]/intREC))<< endl;
cout << "2DY:" << endl;
cout << "INT_Tru: " << intREC                       << "+-" << errREC                                                               << endl;
cout << "INT_Res: " << fResults2D[0][1][0]          << "+-" << fResults2D[0][1][1] << " +-" << fResults2D[0][1][2]                  << endl;
cout << "Res/Tru: " << fResults2D[0][1][0]/intREC   << "+-" << (1/intREC)*(fResults2D[0][1][1]+errREC*(fResults2D[0][1][0]/intREC))<< endl;
cout << endl;
auto err1p  = sqrt(fResults1D[2]*fResults1D[2]+fResults1D[0]*kEventEfficienERP*fResults1D[0]*kEventEfficienERP);
auto err1m  = sqrt(fResults1D[2]*fResults1D[2]+fResults1D[0]*kEventEfficienERM*fResults1D[0]*kEventEfficienERM);
cout << "1D YIELD:" << endl;
cout << "YIELD: " << fResults1D[0]    << "+-" <<  fResults1D[1]   << "+" <<   err1p  << "-" <<  err1m  << endl;
cout << "MeanP: " << fResults1D[3]    << "+-" <<  fResults1D[4]   << endl;
cout << endl;
auto errxp  =   sqrt(fResults2D[0][0][2]*fResults2D[0][0][2]+fResults2D[0][0][0]*kEventEfficienERP*fResults2D[0][0][0]*kEventEfficienERP);
auto errxm  =   sqrt(fResults2D[0][0][2]*fResults2D[0][0][2]+fResults2D[0][0][0]*kEventEfficienERM*fResults2D[0][0][0]*kEventEfficienERM);
auto erryp  =   sqrt(fResults2D[0][1][2]*fResults2D[0][1][2]+fResults2D[0][1][0]*kEventEfficienERP*fResults2D[0][1][0]*kEventEfficienERP);
auto errym  =   sqrt(fResults2D[0][1][2]*fResults2D[0][1][2]+fResults2D[0][1][0]*kEventEfficienERM*fResults2D[0][1][0]*kEventEfficienERM);
cout << "2DX YIELD:" << endl;
cout << "YIELD: " << fResults2D[0][0][0]    << "+-" <<  fResults2D[0][0][1]   << "+" <<  errxp    << "-" <<   errxm  << endl;
cout << "MeanP: " << fResults2D[0][0][3]    << "+-" <<  fResults2D[0][0][4]   << endl;
cout << endl;
cout << "2DY YIELD:" << endl;
cout << "YIELD: " << fResults2D[0][1][0]    << "+-" <<  fResults2D[0][1][1]   << "+" <<   erryp   << "-" <<   errym << endl;
cout << "MeanP: " << fResults2D[0][1][3]    << "+-" <<  fResults2D[0][1][4]   << endl;
cout << endl;
cout << endl;
cout << endl;
cout << "P61D: " << h1D_Tru_P6->Integral("width") << endl;
cout << "P81D: " << h1D_Tru_P8->Integral("width") << endl;
cout << "P62D: " << h2D_Tru_P6->Integral("width") << endl;
cout << "P82D: " << h2D_Tru_P8->Integral("width") << endl;

auto valy   =   fResults2D[0][0][0]/(fResults1D[0]*fResults1D[0]);
auto valx   =   fResults2D[0][1][0]/(fResults1D[0]*fResults1D[0]);
auto errx   =   sqrt(fResults2D[0][0][1]*fResults2D[0][1][1]/(fResults2D[0][0][0]*fResults2D[0][0][0])+4*(fResults1D[1]*fResults1D[1])/(fResults1D[0]*fResults1D[0]));
auto erry   =   sqrt(fResults2D[0][1][1]*fResults2D[0][1][1]/(fResults2D[0][1][0]*fResults2D[0][1][0])+4*(fResults1D[1]*fResults1D[1])/(fResults1D[0]*fResults1D[0]));
auto erxP   =   sqrt((pow(kSyst_SigExtr2D,2)+pow(kEventEfficienERP,2)+(pow((2*fResults1D[2]),2)-pow(2*kSyst_SigExtr1D*fResults1D[0],2))/pow((fResults1D[0]),2)*pow(1+(pow(fResults2D[0][0][0]/fResults1D[0],2)),2)));
auto eryP   =   sqrt((pow(kSyst_SigExtr2D,2)+pow(kEventEfficienERP,2)+(pow((2*fResults1D[2]),2)-pow(2*kSyst_SigExtr1D*fResults1D[0],2))/pow((fResults1D[0]),2)*pow(1+(pow(fResults2D[0][1][0]/fResults1D[0],2)),2)));
auto erxM   =   sqrt((pow(kSyst_SigExtr2D,2)+pow(kEventEfficienERM,2)+(pow((2*fResults1D[2]),2)-pow(2*kSyst_SigExtr1D*fResults1D[0],2))/pow((fResults1D[0]),2)*pow(1+(pow(fResults2D[0][0][0]/fResults1D[0],2)),2)));
auto eryM   =   sqrt((pow(kSyst_SigExtr2D,2)+pow(kEventEfficienERM,2)+(pow((2*fResults1D[2]),2)-pow(2*kSyst_SigExtr1D*fResults1D[0],2))/pow((fResults1D[0]),2)*pow(1+(pow(fResults2D[0][1][0]/fResults1D[0],2)),2)));
cout << endl;
cout << endl;
cout << endl;
cout << "2DX/1D^2:" << valx << "+-" << errx << "+" << valx*erxP << "-" << valx*erxM << endl;
cout << endl;
cout << "2DY/1D^2:" << valy << "+-" << erry << "+" << valx*erxP << "-" << valx*erxM << endl;
cout << endl;
cout << "2DY/1D^2 P6:" << (h2D_Tru_P6->Integral("width"))/(pow(h1D_Tru_P6->Integral("width"),2)) << endl;
cout << "2DY/1D^2 P8:" << (h2D_Tru_P8->Integral("width"))/(pow(h1D_Tru_P8->Integral("width"),2)) << endl;
cout << endl;
auto errsigstat =   sqrt(2*errx*errx + pow(1-2*fResults1D[0],2)*pow(fResults1D[1],2));
auto errsigstatP=   sqrt(2*valx*erxP*valx*erxP + pow(1-2*fResults1D[0],2)*(pow(fResults1D[2],2)+kEventEfficienERP*kEventEfficienERP));
auto errsigstatM=   sqrt(2*valx*erxM*valx*erxM + pow(1-2*fResults1D[0],2)*(pow(fResults1D[2],2)+kEventEfficienERM*kEventEfficienERM));
cout << "Sigma:" << 2*fResults2D[0][0][0] + fResults1D[0] -fResults1D[0]*fResults1D[0] << "+-" << errsigstat  <<"+" << errsigstatP << "-" << errsigstatM << endl;
cout << "Sigma6:" << 2*(h2D_Tru_P6->Integral("width")) + h1D_Tru_P6->Integral("width") - (pow(h1D_Tru_P6->Integral("width"),2)) << endl;
cout << "Sigma8:" << 2*(h2D_Tru_P8->Integral("width")) + h1D_Tru_P8->Integral("width") - (pow(h1D_Tru_P8->Integral("width"),2)) << endl;
cout << endl;
cout << endl;
cout << endl;
auto poisson = 2*fResults2D[0][0][0]/fResults1D[0]-fResults1D[0];
auto poissonstat =   sqrt(4*pow(poisson*fResults1D[1]/fResults1D[0],2) + 4*pow(poisson*fResults2D[0][0][1]/fResults2D[0][0][0],2) + pow(fResults1D[1],2));
auto poissonsyst=   sqrt(poisson*poisson*(kBranchingRatio__*kBranchingRatio__+kSyst_TrackEff1D*kSyst_TrackEff1D+kSyst_PID1D*kSyst_PID1D+kSyst_SigExtr1D*kSyst_SigExtr1D+kSyst_SigExtr2D*kSyst_SigExtr2D)+pow(fResults1D[2],2));
auto errsigst2tM=   sqrt(poisson*poisson*(kBranchingRatio__*kBranchingRatio__+kSyst_TrackEff1D*kSyst_TrackEff1D+kSyst_PID1D*kSyst_PID1D+kSyst_SigExtr1D*kSyst_SigExtr1D+kSyst_SigExtr2D*kSyst_SigExtr2D)+pow(fResults1D[2],2));
cout << "Sigma:" <<  poisson << "+-" << poissonstat  <<"+" << poissonsyst << "-" << poissonsyst << endl;
cout << "Sigma6:" << 2*(h2D_Tru_P6->Integral("width")/h1D_Tru_P6->Integral("width")) - h1D_Tru_P6->Integral("width") << endl;
cout << "Sigma8:" << 2*(h2D_Tru_P8->Integral("width")/h1D_Tru_P8->Integral("width")) - h1D_Tru_P8->Integral("width") << endl;

//---- Graphics of Presentation -----//

TGraphAsymmErrors   *g1D_Res_Stat    =   new TGraphAsymmErrors(h1D_Res);
TGraphAsymmErrors   *g1D_Res_Syst    =   new TGraphAsymmErrors(SetSystErrorsh(h1D_Res));
for ( Int_t iTer = 0; iTer < 4; iTer++ )
{
    g1D_Res_Stat->RemovePoint(0);
    g1D_Res_Syst->RemovePoint(0);
}

TGrCompare1D(g1D_Res_Stat,g1D_Res_Syst,h1D_Tru_P6,h1D_Tru_P8,gCheck_);
TGrCompare2D(h2D_Res,h2D_Tru_P6,h2D_Tru_P8);
TGraphAEGeneratorPT(fResults2D,h2D_Tru_P6,h2D_Tru_P8);

TCanvas *c____1 = new TCanvas("c1", "Intensity of LED 1",0,  0, 800, 600);
c____1->SetMargin(0.2,0.1,0.9,0.9);

TLegend        *lLegend1D           =   new TLegend(0.6,0.45,0.9,0.6);
TMultiGraph    *gFinal1D            =   new TMultiGraph();
TGraphErrors   *gFinal1DDataStat    =   new TGraphErrors();
gFinal1DDataStat->SetPoint(0,1,fResults1D[0]);
gFinal1DDataStat->SetPointError(0,0,fResults1D[1]);
TGraphAsymmErrors   *gFinal1DDataSyst    =   new TGraphAsymmErrors();
gFinal1DDataSyst->SetPoint(0,1,fResults1D[0]);
gFinal1DDataSyst->SetPointError(0,0.025,0.025,sqrt(fResults1D[2]*fResults1D[2]+fResults1D[0]*kEventEfficienERM*fResults1D[0]*kEventEfficienERM),sqrt(fResults1D[2]*fResults1D[2]+fResults1D[0]*kEventEfficienERP*fResults1D[0]*kEventEfficienERP));
TGraphErrors   *gFinal1DPythia6_    =   new TGraphErrors();
gFinal1DPythia6_->SetPoint(0,1.05,h1D_Tru_P6->Integral("width"));
TGraphErrors   *gFinal1DPythia8_    =   new TGraphErrors();
gFinal1DPythia8_->SetPoint(0,1.1,h1D_Tru_P8->Integral("width"));

gFinal1DDataStat->SetMarkerStyle(20);
gFinal1DDataStat->SetMarkerColor(46);
gFinal1DDataStat->SetLineColor(kRed);
gFinal1DDataSyst->SetMarkerStyle(20);
gFinal1DDataSyst->SetMarkerColor(46);
gFinal1DDataSyst->SetLineColor(kRed);
gFinal1DPythia6_->SetMarkerStyle(26);
gFinal1DPythia6_->SetMarkerColor(kBlue);
gFinal1DPythia8_->SetMarkerStyle(32);
gFinal1DPythia8_->SetMarkerColor(kBlue+3);

gFinal1D->Add(gFinal1DDataSyst,"EP5");
gFinal1D->Add(gFinal1DDataStat,"EP");
gFinal1D->Add(gFinal1DPythia6_,"P");
gFinal1D->Add(gFinal1DPythia8_,"P");
gFinal1D->GetYaxis()->SetTitle("#frac{dN_{#phi}}{dy}");
gFinal1D->GetXaxis()->SetLimits(0.85,1.5);
gFinal1D->SetMaximum(0.041);
gFinal1D->SetMinimum(0.025);
TAxis *xAxis1  =   new TAxis(*gFinal1D->GetYaxis());

lLegend1D->AddEntry(gFinal1DDataStat,"Stat. Err.","EP");
lLegend1D->AddEntry(gFinal1DDataSyst,"Syst. Err.","F");
lLegend1D->AddEntry(gFinal1DPythia6_,"Pythia 6","P");
lLegend1D->AddEntry(gFinal1DPythia8_,"Pythia 8","P");

gFinal1D->Draw("PA");
lLegend1D->Draw("same");
xAxis1->Draw();
c____1->Write();
c____1->SaveAs("./graphs/c____1.pdf");
c____1->SaveAs("./graphs/c____1.png");

TCanvas *c____2 = new TCanvas("c22", "Intensity of LED 1",0,  0, 800, 600);
c____2->SetMargin(0.2,0.1,0.9,0.9);

TMultiGraph    *gFinal2D            =   new TMultiGraph();
TGraphErrors   *gFinal2DDataStat    =   new TGraphErrors();
gFinal2DDataStat->SetPoint(0,1,fResults2D[0][0][0]);
gFinal2DDataStat->SetPointError(0,0,fResults2D[0][0][1]);
TGraphAsymmErrors   *gFinal2DDataSyst    =   new TGraphAsymmErrors();
gFinal2DDataSyst->SetPoint(0,1,fResults2D[0][0][0]);
gFinal2DDataSyst->SetPointError(0,0.025,0.025,sqrt(fResults2D[0][0][2]*fResults2D[0][0][2]+fResults2D[0][0][0]*kEventEfficienERM*fResults2D[0][0][0]*kEventEfficienERM),sqrt(fResults2D[0][0][2]*fResults2D[0][0][2]+fResults2D[0][0][0]*kEventEfficienERP*fResults2D[0][0][0]*kEventEfficienERP));
TGraphErrors   *gFinal2DPythia6_    =   new TGraphErrors();
gFinal2DPythia6_->SetPoint(0,1.05,h2D_Tru_P6->Integral("width"));
TGraphErrors   *gFinal2DPythia8_    =   new TGraphErrors();
gFinal2DPythia8_->SetPoint(0,1.1,h2D_Tru_P8->Integral("width"));

gFinal2DDataStat->SetMarkerStyle(20);
gFinal2DDataStat->SetMarkerColor(46);
gFinal2DDataStat->SetLineColor(kRed);
gFinal2DDataSyst->SetMarkerStyle(20);
gFinal2DDataSyst->SetMarkerColor(46);
gFinal2DDataSyst->SetLineColor(kRed);
gFinal2DPythia6_->SetMarkerStyle(26);
gFinal2DPythia6_->SetMarkerColor(kBlue);
gFinal2DPythia8_->SetMarkerStyle(32);
gFinal2DPythia8_->SetMarkerColor(kBlue+3);

gFinal2D->Add(gFinal2DDataSyst,"EP5");
gFinal2D->Add(gFinal2DDataStat,"EP");
gFinal2D->Add(gFinal2DPythia6_,"P");
gFinal2D->Add(gFinal2DPythia8_,"P");
gFinal2D->GetYaxis()->SetTitle("#frac{dN_{#phi#phi}}{dy}");
gFinal2D->GetXaxis()->SetLimits(0.85,1.5);
gFinal2D->SetMaximum(0.0021);
gFinal2D->SetMinimum(0.0008);
TAxis *xAxis2  =   new TAxis(*gFinal2D->GetYaxis());

gFinal2D->Draw("PA");
lLegend1D->Draw("same");
xAxis2->Draw();
c____2->Write();
c____2->SaveAs("./graphs/c____2.pdf");
c____2->SaveAs("./graphs/c____2.png");

TCanvas *c____3 = new TCanvas("c32", "Intensity of LED 1",0,  0, 800, 600);
c____3->SetMargin(0.2,0.1,0.9,0.9);

TMultiGraph    *gFinal3D            =   new TMultiGraph();
TGraphErrors   *gFinal3DDataStat    =   new TGraphErrors();
gFinal3DDataStat->SetPoint(0,1,poisson);
gFinal3DDataStat->SetPointError(0,0,poissonstat);
TGraphAsymmErrors   *gFinal3DDataSyst    =   new TGraphAsymmErrors();
gFinal3DDataSyst->SetPoint(0,1,poisson);
gFinal3DDataSyst->SetPointError(0,0.025,0.025,poissonsyst,poissonsyst);
TGraphErrors   *gFinal3DPythia6_    =   new TGraphErrors();
gFinal3DPythia6_->SetPoint(0,1.05,2*(h2D_Tru_P6->Integral("width")/h1D_Tru_P6->Integral("width")) - h1D_Tru_P6->Integral("width"));
TGraphErrors   *gFinal3DPythia8_    =   new TGraphErrors();
gFinal3DPythia8_->SetPoint(0,1.1,2*(h2D_Tru_P8->Integral("width")/h1D_Tru_P8->Integral("width")) - h1D_Tru_P8->Integral("width"));
TGraph * hPoissonLine = new TGraph();
hPoissonLine->SetPoint(0,0,0);
hPoissonLine->SetPoint(1,3,0);

lLegend1D->AddEntry(hPoissonLine,"Pois. Distr.","L");

gFinal3DDataStat->SetMarkerStyle(20);
gFinal3DDataStat->SetMarkerColor(46);
gFinal3DDataStat->SetLineColor(kRed);
gFinal3DDataSyst->SetMarkerStyle(20);
gFinal3DDataSyst->SetMarkerColor(46);
gFinal3DDataSyst->SetLineColor(kRed);
gFinal3DPythia6_->SetMarkerStyle(26);
gFinal3DPythia6_->SetMarkerColor(kBlue);
gFinal3DPythia8_->SetMarkerStyle(32);
gFinal3DPythia8_->SetMarkerColor(kBlue+3);

gFinal3D->Add(gFinal3DDataSyst,"EP5");
gFinal3D->Add(gFinal3DDataStat,"EP");
gFinal3D->Add(gFinal3DPythia6_,"P");
gFinal3D->Add(gFinal3DPythia8_,"P");
gFinal3D->Add(hPoissonLine,"L");
gFinal3D->GetYaxis()->SetTitle("#sigma^{2}_{#phi}/#mu_{#phi} -1");
gFinal3D->GetXaxis()->SetLimits(0.85,1.5);
gFinal3D->SetMaximum(0.1);
gFinal3D->SetMinimum(-0.05);
TAxis *xAxis3  =   new TAxis(*gFinal3D->GetYaxis());

gFinal3D->Draw("PA");
lLegend1D->Draw("same");
xAxis3->Draw();
c____3->Write();
c____3->SaveAs("./graphs/c____3.pdf");
c____3->SaveAs("./graphs/c____3.png");

TCanvas *c____4 = new TCanvas("c3244", "Intensity of LED 1",0,  0, 800, 600);
c____4->SetMargin(0.2,0.1,0.9,0.9);

TMultiGraph    *gFinal4D            =   new TMultiGraph();
TGraphErrors   *gFinal4DDataStat    =   new TGraphErrors();
gFinal4DDataStat->SetPoint(0,1,valx);
gFinal4DDataStat->SetPointError(0,0,errx);
TGraphAsymmErrors   *gFinal4DDataSyst    =   new TGraphAsymmErrors();
gFinal4DDataSyst->SetPoint(0,1,valx);
gFinal4DDataSyst->SetPointError(0,0.025,0.025,erxM,erxP);
TGraphErrors   *gFinal4DPythia6_    =   new TGraphErrors();
gFinal4DPythia6_->SetPoint(0,1.05,(h2D_Tru_P6->Integral("width")/pow(h1D_Tru_P6->Integral("width"),2)));
TGraphErrors   *gFinal4DPythia8_    =   new TGraphErrors();
gFinal4DPythia8_->SetPoint(0,1.1,(h2D_Tru_P8->Integral("width")/pow(h1D_Tru_P8->Integral("width"),2)));
TGraph * hPoissonLin2 = new TGraph();
hPoissonLin2->SetPoint(0,0,0.5);
hPoissonLin2->SetPoint(1,3,0.5);

gFinal4DDataStat->SetMarkerStyle(20);
gFinal4DDataStat->SetMarkerColor(46);
gFinal4DDataStat->SetLineColor(kRed);
gFinal4DDataSyst->SetMarkerStyle(20);
gFinal4DDataSyst->SetMarkerColor(46);
gFinal4DDataSyst->SetLineColor(kRed);
gFinal4DPythia6_->SetMarkerStyle(26);
gFinal4DPythia6_->SetMarkerColor(kBlue);
gFinal4DPythia8_->SetMarkerStyle(32);
gFinal4DPythia8_->SetMarkerColor(kBlue+3);

gFinal4D->Add(gFinal4DDataSyst,"EP5");
gFinal4D->Add(gFinal4DDataStat,"EP");
gFinal4D->Add(gFinal4DPythia6_,"P");
gFinal4D->Add(gFinal4DPythia8_,"P");
gFinal4D->Add(hPoissonLin2,"L");
gFinal4D->GetYaxis()->SetTitle("#frac{dN_{#phi#phi}}{dy}/(#frac{dN_{#phi}}{dy})^{2}");
gFinal4D->GetXaxis()->SetLimits(0.85,1.5);
gFinal4D->SetMaximum(1.6);
gFinal4D->SetMinimum(0.4);
TAxis *xAxis4  =   new TAxis(*gFinal4D->GetYaxis());

gFinal4D->Draw("PA");
lLegend1D->Draw("same");
xAxis4->Draw();
c____4->Write();
c____4->SaveAs("./graphs/c____4.pdf");
c____4->SaveAs("./graphs/c____4.png");

ratioplot(h1D_Res,h1D_Tru_P6,h1D_Tru_P8,"1D");
TH1D * hSlice6 = new TH1D ("hSlice6","hSlice6",nBinPT2D,fArrPT2D);
TH1D * hSlice8 = new TH1D ("hSlice8","hSlice8",nBinPT2D,fArrPT2D);
for ( Int_t iTer = 0; iTer < nBinPT2D; iTer++ )
{
    auto    hProjXPythia6   =   h2D_Tru_P6->ProjectionX(Form("P6X_%i",iTer),iTer+1,iTer+1);
    auto    hProjXPythia8   =   h2D_Tru_P8->ProjectionX(Form("P8X_%i",iTer),iTer+1,iTer+1);
    auto    hProjXDataset   =   h2D_Res->ProjectionX(Form("DTX_%i",iTer),iTer+1,iTer+1);
    auto    hProjYPythia6   =   h2D_Tru_P6->ProjectionY(Form("P6Y_%i",iTer),iTer+1,iTer+1);
    auto    hProjYPythia8   =   h2D_Tru_P8->ProjectionY(Form("P8Y_%i",iTer),iTer+1,iTer+1);
    auto    hProjYDataset   =   h2D_Res->ProjectionY(Form("DTY_%i",iTer),iTer+1,iTer+1);
    hSlice6->SetBinContent(iTer+1,hProjXPythia6->Integral("width"));
    hSlice8->SetBinContent(iTer+1,hProjXPythia8->Integral("width"));
    if ( iTer > 1 )
    {
        ratioplot(hProjXDataset,hProjXPythia6,hProjXPythia8,Form("2DX_%i",iTer+1));
        ratioplot(hProjYDataset,hProjYPythia6,hProjYPythia8,Form("2DY_%i",iTer+1));
    }
}

ratioplot(h1D_Raw,h1D_Rec_P6,h1D_Rec_P8,"1D_Raw");
for ( Int_t iTer = 0; iTer < nBinPT2D; iTer++ )
{
    auto    hProjXPythia6_   =   h2D_Rec_P6->ProjectionX(Form("P6X_%i_Raw",iTer),iTer+1,iTer+1);
    auto    hProjXPythia8_   =   h2D_Rec_P8->ProjectionX(Form("P8X_%i_Raw",iTer),iTer+1,iTer+1);
    auto    hProjXDataset_   =   h2D_Raw->ProjectionX(Form("DTX_%i_Raw",iTer),iTer+1,iTer+1);
    auto    hProjYPythia6_   =   h2D_Rec_P6->ProjectionY(Form("P6Y_%i_Raw",iTer),iTer+1,iTer+1);
    auto    hProjYPythia8_   =   h2D_Rec_P8->ProjectionY(Form("P8Y_%i_Raw",iTer),iTer+1,iTer+1);
    auto    hProjYDataset_   =   h2D_Raw->ProjectionY(Form("DTY_%i_Raw",iTer),iTer+1,iTer+1);
    if ( iTer > 1 )
    {
        ratioplot(hProjXDataset_,hProjXPythia6_,hProjXPythia8_,Form("2DX_%i_Raw",iTer+1));
        ratioplot(hProjYDataset_,hProjYPythia6_,hProjYPythia8_,Form("2DY_%i_Raw",iTer+1));
    }
}


TCanvas *cCompareFinal  =   new TCanvas();
TH1D * hData__ = TH1DGeneratorYield(fResults2D,false);
hData__->GetXaxis()->SetTitle("P_{T}#phi_{2} (GeV/c)");
hData__->GetYaxis()->SetTitle("#frac{d^{2}N#phi#phi}{dydp_{T}#phi_{2}}(GeV/c)^{-1}");
TLegend * lLegend1          =   new TLegend(0.65,0.65,0.85,0.85);
lLegend1                    ->SetFillColor(kWhite);
lLegend1                    ->SetLineColor(kWhite);
lLegend1                    ->AddEntry(hSlice8,"Res","L");
lLegend1                    ->AddEntry(hData__,"MC","EP");
gPad->SetLogy();
hSlice8->SetLineColor(kBlue);
hSlice8->Draw("HIST L");
hData__->SetMarkerStyle(25);
hData__->SetMarkerColor(kRed);
hData__->Draw("SAME EP");
lLegend1->Draw("SAME");
cCompareFinal->SaveAs("./graphs/cCompareFinal.pdf");
ratioplot(fResults2D,hSlice6,hSlice8);
 
//---------------------//
// Output and wrap up  //-------------------------------------------------------------------------------
//---------------------//

// Output File for Results
TFile*  outFile_RS  =   new TFile(fAnlResults,"recreate");

h1D_Tru->Write();
h1D_Res->Write();
h2D_Tru->Write();
h2D_Res->Write();

// Closing opened files
outFile_FT->Close();
outFile_RS->Close();
insFile_DT->Close();
insFile_EF->Close();

}
//

//_____________________________________________________________________________
//                                                  //  !TODO: Make the kMarkerstyle in a kArray
//                                                  //  Relates to global variable kRainbowColor
template < typename Tclass >
void                    fSetRainbowStyle            ( std::vector<Tclass**> fInputObjectLits, Int_t fStartIteratorAt = 0 )  {
    TCanvas*    fResult     =   new TCanvas();
    Int_t       fIterator   =   fStartIteratorAt;
    for ( Tclass*&fObjectToPlot   :   fInputObjectLits )  {
        if ( fIterator == 12 ) fIterator = 0;
        fObjectToPlot->SetMarkerStyle(20+fIterator >= 31 ? 20+fIterator+1 : 20+fIterator );
        fObjectToPlot->SetMarkerColor(kRainbowColor[fIterator]);
        fObjectToPlot->SetLineWidth(2.);
        fObjectToPlot->SetLineStyle(9);
        fObjectToPlot->SetLineColor(kRainbowColor[fIterator]);
        fIterator++;
    }
    return fResult;
}
//
//_____________________________________________________________________________
//
template < typename Tclass >
void                    fSetUniformBinning          ( Tclass *fArrBin, Tclass fMinBin, Tclass fMaxBin, Int_t fNPoints )  {
    for (int i = 0; i <= fNPoints; i++ )
    {
        fArrBin[i] = fMinBin+(i)*(fMaxBin - fMinBin)/(static_cast<Float_t>(fNPoints));
    }
}
////_____________________________________________________________________________
//
void
uSetHisto
( TH1* hTarget, TString fOption ){
    hTarget->SetTitle("");
    if ( fOption.Contains("1D") || fOption.Contains("12D") )   {
        if ( fOption.Contains("EFF") ) {
            // Y-Axis limits
            hTarget->Scale(100.);
            hTarget->SetMaximum(100.);
            hTarget->SetMinimum(0.);
            // Y-Axis title
            hTarget->GetYaxis()->SetTitle("Efficiency #times Acceptance (%)");
            // Marker style
            hTarget->SetMarkerStyle(kMarkers[5]);
            // Colour scheme
            hTarget->SetMarkerColor(kColors[2]);
            hTarget->SetLineColor(kColors[2]);
            hTarget->SetFillColorAlpha(kFillColors[2],0.33);
            //
            hTarget->SetOption("EP L");
            if ( fOption.Contains("EFF2") ) {
                hTarget->SetMarkerStyle(kMarkers[9]);
                hTarget->SetMarkerColor(kColors[1]);
                hTarget->SetLineColor(kColors[1]);
                hTarget->SetFillColorAlpha(kFillColors[1],0.33);
                hTarget->SetOption("EP L SAME");
            }
            if ( fOption.Contains("SL") ) {
                // Y-Axis title
                hTarget->GetYaxis()->SetTitle("Signal Loss (%)");
                uOffset(hTarget,-100);
                hTarget->SetMaximum(10.);
                hTarget->SetMinimum(-1.);
            }
        }
        if ( fOption.Contains("SPT") ) {
            
             // Preferred kColors and kMarkers
             //const Int_t kFillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
             //const Int_t kColors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
             //const Int_t kMarkers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};
            
            hTarget->SetMarkerStyle(kMarkers[3]);
            hTarget->SetMarkerColor(kColors[2]);
            hTarget->SetMarkerSize(1);
            hTarget->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            hTarget->GetYaxis()->SetTitleOffset(1.5);
            hTarget->GetYaxis()->SetTitle("1/N_{ev}dN^{2}/(dydp_{T})");
            if ( fOption.Contains("12D") )  {
                hTarget->GetYaxis()->SetTitle("1/N_{ev}dN^{3}/(dydp_{T,#phi_{1}}dp_{T,#phi_{2}})");
                hTarget->GetXaxis()->SetTitle("#it{p}_{T,#phi_{1}} (GeV/#it{c})");
            }
        }
        if ( fOption.Contains("MPT") ) {
            hTarget->SetMarkerStyle(kMarkers[3]);
            hTarget->SetMarkerColor(kColors[2]);
            hTarget->SetMarkerSize(1);
            hTarget->GetXaxis()->SetTitleOffset(1.1);
            hTarget->GetXaxis()->SetTitle("#it{p}_{T,#phi_{1}} (GeV/#it{c})");
            hTarget->GetYaxis()->SetTitleOffset(1.5);
            hTarget->GetYaxis()->SetTitle("#LT #it{p}_{T,#phi_{2}} #GT (GeV/#it{c})");
        }
        if ( fOption.Contains("STAT") ) {
            hTarget->SetLineColor(kColors[2]);
            hTarget->SetFillColorAlpha(kColors[2],0.33);
            hTarget->SetOption("PE2");
        }
        if ( fOption.Contains("SYST") ) {
            hTarget->SetLineColor(kColors[3]);
            hTarget->SetFillColorAlpha(kColors[3],0.);
            hTarget->SetOption("PE2");
        }
    } else if ( fOption.Contains("2D") )   {
        
    } else if ( fOption.Contains("3D") )   {
        cout << " Buu 3D " << endl;
    } else  {
        cout << " CANT GUESS DIMENSION " << endl;
    }
}
//
TCanvas*
uPlotSpectrum
( TH1* hTarget, TH1* hTrSyst, TString fOption = "" ){
    //
    SetStyle();
    //
    TCanvas    *cDrawResult =   new TCanvas();
    gStyle->SetOptStat(0);
    if ( fOption.Contains("SPT") ) gPad->SetLogy();
    uSetHisto(hTarget,fOption + TString(" STAT"));
    uSetHisto(hTrSyst,fOption + TString(" SYST"));
    hTarget->SetMaximum(1.25*max(hTarget->GetMaximum(),hTrSyst->GetMaximum()));
    hTarget->SetMinimum(0.75*min(hTarget->GetMinimum(),hTrSyst->GetMinimum()));
    if ( fOption.Contains("SPT") )  hTarget->SetMaximum(2.0*max(hTarget->GetMaximum(),hTrSyst->GetMaximum()));
    if ( fOption.Contains("SPT") )  hTarget->SetMinimum(0.5*min(hTarget->GetMinimum(),hTrSyst->GetMinimum()));
    //
    TLegend    *lLegend;
    if ( fOption.Contains("R") )    lLegend =   new TLegend(0.65,0.35,0.85,0.5);
        else                        lLegend =   new TLegend(0.2,0.35,0.4,0.5);
    lLegend     ->  SetFillColorAlpha(0.,0.);
    lLegend     ->  AddEntry    (hTarget,"Data","P");
    lLegend     ->  AddEntry    (hTarget,"Stat","F");
    lLegend     ->  AddEntry    (hTrSyst,"Syst","F");
    //
    hTarget->Draw();
    hTrSyst->Draw("SAME E1");
    lLegend->Draw("SAME");
    //
    if ( fOption.Contains("R") )    {
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.65, 0.3,"ALICE");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        uLatex->DrawLatexNDC(0.65, 0.25,"pp #sqrt{#it{s}}= 7 TeV");
        uLatex->DrawLatexNDC(0.65, 0.2,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
    } else  {
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.20, 0.3,"ALICE");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        uLatex->DrawLatexNDC(0.20, 0.25,"pp #sqrt{#it{s}}= 7 TeV");
        uLatex->DrawLatexNDC(0.20, 0.2,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
    }
    //
    return cDrawResult;
}
//
void
uSetHisto
( TGraphMultiErrors* hTarget, TString fOption = "" ){
    hTarget->SetTitle("");
    hTarget->GetYaxis()->SetTitle("1/N_{ev}dN/dy");
    //hTarget->GetXaxis()->SetNdivisions(2);
    //hTarget->GetXaxis()->SetBinLabel(hTarget->GetXaxis()->FindBin(1),"#LT Y_{1#phi} #GT");
    //hTarget->GetXaxis()->SetBinLabel(hTarget->GetXaxis()->FindBin(2),"#LT Y_{2#phi} #GT");
    //hTarget->GetXaxis()->LabelsOption("h");
    hTarget->GetXaxis()->SetLabelSize(0.08);
    hTarget->SetLineColorAlpha(0.,0.);
    hTarget->SetMarkerStyle(21);
    hTarget->SetMarkerColor(kRed);
    hTarget->GetAttLine(0)->SetLineColor(38);
    hTarget->GetAttLine(1)->SetLineColor(46);
    hTarget->GetAttFill(0)->SetFillColorAlpha(38,0.33);
    hTarget->GetAttFill(1)->SetFillColorAlpha(46,0.33);
}
//
//  --  --  Data Handling Utilities  --  --  //
//
//_____________________________________________________________________________
//
TGraphAsymmErrors      *fSumErrors                  ( TGraphAsymmErrors* gBasic, TGraphAsymmErrors* gAddition )    {
    //  Checking the consistency of TGraphs
    Int_t   fNPoints =   gBasic ->  GetN();
    if  ( fNPoints  != gAddition ->  GetN() )   {
        cout << "[ERROR] Systematics and Statistics do not have the same number of points! Skipping this one..." << endl;
        return nullptr;
    }
    //
    TGraphAsymmErrors  *fResult =   new TGraphAsymmErrors(*gBasic);
    for ( Int_t iFit = 0; iFit < fNPoints; iFit++ ) {
        auto    fXErrBsicLow    =   ( gBasic ->  GetErrorXlow(iFit) );
        auto    fXErrBsicHigh   =   ( gBasic ->  GetErrorXhigh(iFit) );
        auto    fYErrBsicLow    =   ( gBasic ->  GetErrorYlow(iFit) );
        auto    fYErrBsicHigh   =   ( gBasic ->  GetErrorYhigh(iFit) );
        auto    fYErrAddtLow    =   ( gAddition ->  GetErrorYlow(iFit) );
        auto    fYErrAddtHigh   =   ( gAddition ->  GetErrorYhigh(iFit) );
        fResult ->  SetPointError(iFit,fXErrBsicLow,fXErrBsicHigh,sqrt(fYErrBsicLow*fYErrBsicLow+fYErrAddtLow*fYErrAddtLow),sqrt(fYErrBsicHigh*fYErrBsicHigh+fYErrAddtHigh*fYErrAddtHigh));
    }
    //
    return  fResult;
}
//
template < class Tclass >
Tclass                   *fSumErrors                  ( Tclass* gBasic, Tclass* gAddition )    {
    Tclass  *fResult =   new Tclass(*gBasic);
    for ( Int_t iBin = 0; iBin < gBasic->GetNbinsX(); iBin++ ) {
        fResult ->  SetBinError( iBin, SquareSum( { gBasic->GetBinError(iBin), gAddition->GetBinError(iBin) } ) );
    }
    //
    return  fResult;
}
//
//_____________________________________________________________________________
//

TGraphAsymmErrors      *fScaleWithError             ( TGraphAsymmErrors* gBasic, Double_t fScale, Double_t fScaleErrHigh = 0., Double_t fScaleErrLow = 0. )    {
    TGraphAsymmErrors  *fResult =   new TGraphAsymmErrors(*gBasic);
    for ( Int_t iFit = 0; iFit < fResult->GetN(); iFit++ ) {
        auto    fYValue         =   ( gBasic ->  GetPointY(iFit) );
        auto    fYErrBsicLow    =   ( gBasic ->  GetErrorYlow(iFit) );
        auto    fYErrBsicHigh   =   ( gBasic ->  GetErrorYhigh(iFit) );
        fResult ->  SetPointY       ( iFit, fYValue*fScale);
        fResult ->  SetPointEYhigh  ( iFit, (fYValue*fScale)*sqrt(fYErrBsicHigh*fYErrBsicHigh/(fYValue*fYValue) + fScaleErrHigh*fScaleErrHigh/(fScale*fScale)));
        fResult ->  SetPointEYlow   ( iFit, (fYValue*fScale)*sqrt(fYErrBsicLow*fYErrBsicLow/(fYValue*fYValue) + fScaleErrLow*fScaleErrLow/(fScale*fScale)));
    }
    //
    return  fResult;
}
TH1F                   *fScaleWithError             ( TH1F* gBasic, Double_t fScale, Double_t fScaleError = 0. )    {
    TH1F  *fResult =   new TH1F(*gBasic);
    for ( Int_t iBin = 1; iBin <= fResult->GetNbinsX(); iBin++ ) {
        fResult->SetBinContent( iBin, gBasic->GetBinContent(iBin)/fScale );
        fResult->SetBinError( iBin, (gBasic->GetBinContent(iBin)/fScale)*SquareSum( { gBasic->GetBinError( iBin )/gBasic->GetBinContent(iBin), fScaleError/fScale } ) );
    }
    return  fResult;
}
//
//_____________________________________________________________________________
//
TGraphAsymmErrors      *uRandomizePoints            ( TGraphAsymmErrors* gStatic, TGraphAsymmErrors* gMoveable )    {
    //  Checking the consistency of TGraphs
    Int_t   fNPoints =   gStatic ->  GetN();
    if  ( fNPoints  != gMoveable ->  GetN() )
    {
        cout << "[ERROR] Systematics and Statistics do not have the same number of points! Skipping this one..." << endl;
        return nullptr;
    }
    //
    TGraphAsymmErrors  *fResult =   new TGraphAsymmErrors(*gStatic);
    for ( Int_t iFit = 0; iFit < fNPoints; iFit++ ) {
        auto    fXValue         =   ( gStatic ->  GetPointX(iFit) );
        auto    fYValue         =   ( gStatic ->  GetPointY(iFit) );
        auto    fYErrMvblLow    =   ( gMoveable ->  GetErrorYlow(iFit) );
        auto    fYErrMvblHigh   =   ( gMoveable ->  GetErrorYhigh(iFit) );
        auto    fIsFluctLow     =   true;
        auto    fYNewValue      =   fYValue;
        ( uRandomGen  ->  Uniform (0.,1.) ) > 0.5? fIsFluctLow = true : fIsFluctLow = false;
        if ( fIsFluctLow )  {
            fYNewValue  -= fabs(uRandomGen  ->  Gaus(fYValue,fYErrMvblLow) - fYValue);
        }   else    {
            fYNewValue  += fabs(uRandomGen  ->  Gaus(fYValue,fYErrMvblHigh) - fYValue);
        }
        fResult->SetPoint(iFit,fXValue,fYNewValue);
    }
    return  fSumErrors(fResult,gMoveable);
}
std::vector<TH1F*>      uRandomizePoints            ( std::vector<TH1F*>  gStatic, std::vector<TH1F*>  gMoveable )    {
    std::vector<TH1F*>  fResult;
    auto iTer = 0;
    for ( auto hStatic : gStatic ) {
        fResult.push_back(uRandomizePoints((TH1F*)hStatic->Clone(),(TH1F*)gMoveable.at(iTer)->Clone()));
        iTer++;
    }
    return  fResult;
}
std::vector<TH1F*>      uRandomizePointsSymm        ( std::vector<TH1F*>  gStatic, std::vector<TH1F*>  gMoveable )    {
    std::vector<TH1F*>  fResult;
    for ( auto hTarget : gStatic )  {
        fResult.push_back(new TH1F(*hTarget));
    }
    auto nHisto =   gStatic.size();
    auto nBins  =   gStatic.at(0)->GetNbinsX();
    for ( Int_t iTer = 1; iTer <= nHisto; iTer++ ) {
        for ( Int_t jTer = iTer; jTer <= nBins; jTer++ ) {
            auto    fBinContent =   gStatic.at(iTer-1)    ->GetBinContent (jTer);
            auto    fStatError  =   gStatic.at(iTer-1)    ->GetBinError   (jTer);
            auto    fMoveError  =   gMoveable.at(iTer-1)  ->GetBinError   (jTer);
            auto    fYNewValue  =   max( 0., uRandomGen -> Gaus(fBinContent,fMoveError) );
            auto    fYNewError  =   sqrt( fStatError*fStatError + fMoveError*fMoveError );
            fResult.at(iTer-1)->SetBinContent ( jTer, fYNewValue );
            fResult.at(iTer-1)->SetBinError   ( jTer, fYNewError );
            fResult.at(jTer-1)->SetBinContent ( iTer, fYNewValue );
            fResult.at(jTer-1)->SetBinError   ( iTer, fYNewError );
        }
    }
    return  fResult;
}
std::vector<TH1D*>      uRandomizePointsSymm        ( std::vector<TH1D*>  gStatic, std::vector<TH1D*>  gMoveable )    {
    std::vector<TH1D*>  fResult;
    for ( auto hTarget : gStatic )  {
        fResult.push_back(new TH1D(*hTarget));
    }
    auto nHisto =   gStatic.size();
    auto nBins  =   gStatic.at(0)->GetNbinsX();
    for ( Int_t iTer = 1; iTer <= nHisto; iTer++ ) {
        for ( Int_t jTer = iTer; jTer <= nBins; jTer++ ) {
            auto    fBinContent =   gStatic.at(iTer-1)    ->GetBinContent (jTer);
            auto    fStatError  =   gStatic.at(iTer-1)    ->GetBinError   (jTer);
            auto    fMoveError  =   gMoveable.at(iTer-1)  ->GetBinError   (jTer);
            auto    fYNewValue  =   max( 0., uRandomGen -> Gaus(fBinContent,fMoveError) );
            auto    fYNewError  =   sqrt( fStatError*fStatError + fMoveError*fMoveError );
            fResult.at(iTer-1)->SetBinContent ( jTer, fYNewValue );
            fResult.at(iTer-1)->SetBinError   ( jTer, fYNewError );
            fResult.at(jTer-1)->SetBinContent ( iTer, fYNewValue );
            fResult.at(jTer-1)->SetBinError   ( iTer, fYNewError );
        }
    }
    return  fResult;
}
//
//
//_____________________________________________________________________________
//
TH1F                    *fEfficiencycorrection       ( TH1   *fToBeCorrected, TH1    *fAccepted,  TH1   *fTotal,    Double_t fScale = 1. )  {
    TH1F   *fEfficiency =   (TH1F*)fAccepted->Clone();
    TH1F   *fResult     =   (TH1F*)fToBeCorrected->Clone();
    fEfficiency         ->  Divide(fAccepted,fTotal,1.,1.,"b");
    fResult             ->  Divide(fToBeCorrected,fEfficiency,fScale);
    return  fResult;
}
//
//_____________________________________________________________________________
//
TGraphAsymmErrors      *fEfficiencycorrection       ( TGraphAsymmErrors *fToBeCorrected, TH1    *fAccepted,  TH1    *fTotal,    Double_t fScale = 1. )  {
    //  Copying accordingly the TH1*
    TGraphAsymmErrors  *fResult =   new TGraphAsymmErrors(*fToBeCorrected);
    //
    fResult->Divide(fAccepted,fTotal,"cl=0.683 b(1,1) mode");
    //
    Int_t   fNPoints =   fResult ->  GetN();
    for ( Int_t iFit = 0; iFit < fNPoints; iFit++ ) {
        auto    fYValueResult   =   ( fToBeCorrected ->  GetPointY(iFit) );
        auto    fYErrorRstHigh  =   ( fToBeCorrected ->  GetErrorYhigh(iFit) );
        auto    fYErrorRstLow   =   ( fToBeCorrected ->  GetErrorYlow(iFit) );
        auto    fYValueEffic    =   ( fResult ->  GetPointY(iFit) );
        auto    fYErrorEffHigh  =   ( fResult ->  GetErrorYhigh(iFit) );
        auto    fYErrorEffLow   =   ( fResult ->  GetErrorYlow(iFit) );
        fResult ->  SetPointY       (iFit,fScale*fYValueResult/fYValueEffic);
        fResult ->  SetPointEYlow   (iFit,(fScale*fYValueResult/fYValueEffic)*sqrt(fYErrorEffHigh*fYErrorEffHigh/(fYValueEffic*fYValueEffic) + fYErrorRstHigh*fYErrorRstHigh/(fYValueResult*fYValueResult) ));
        fResult ->  SetPointEYhigh  (iFit,(fScale*fYValueResult/fYValueEffic)*sqrt(fYErrorEffLow*fYErrorEffLow/(fYValueEffic*fYValueEffic) + fYErrorRstLow*fYErrorRstLow/(fYValueResult*fYValueResult) ));
    }
    //
    return  fResult;
}
//
//_____________________________________________________________________________
//                                                          // !TODO: To Be Implemented (...)
std::vector<TGraphAsymmErrors*> fEfficiencycorrection       ( TH2   *fToBeCorrected, TH2    *fAccepted,  TH2    *fTotal,    Double_t fScale = 1. )  {
    return std::vector<TGraphAsymmErrors*>();
}
//
std::vector<TH1F*> fEfficiencycorrection       ( TH2   *fToBeCorrected, TH1    *fAccepted,  TH1    *fTotal,    Double_t fScale = 1. )  {
    std::vector<TH1F*> fResult;
    if ( !fToBeCorrected )  { cout << "No fToBeCorrected" << endl; return fResult; }
    if ( !fAccepted )  { cout << "No fAccepted" << endl; return fResult; }
    if ( !fTotal )  { cout << "No fTotal" << endl; return fResult; }
    TH1F   *fEfficiency =   (TH1F*)fAccepted->Clone();
    fEfficiency         ->  Divide(fAccepted,fTotal,1.,1.,"b");
    for ( Int_t iHisto = 1; iHisto <= fToBeCorrected->GetNbinsY(); iHisto++ )    {
        auto    fConditional    =   fToBeCorrected->ProjectionY(Form("dd_%i",iHisto),iHisto,iHisto);
        TH1F*    fCorrCondit    =   fEfficiencycorrection( fConditional, fAccepted, fTotal, fScale );
        fResult.push_back( fScaleWithError( fCorrCondit, fEfficiency->GetBinContent(iHisto), fEfficiency->GetBinError(iHisto) ) );
    }
    return fResult;
}
//
//_____________________________________________________________________________
//
TGraphAsymmErrors*              fEfficiencycorrection       ( TGraphAsymmErrors    *fToBeCorrected,    TGraphAsymmErrors   *fEfficiency,     Double_t fScale = 1. )  {
    //  Checking the consistency of TGraphs
    Int_t   fNPoints =   fToBeCorrected ->  GetN();
    if  ( fNPoints  != fEfficiency ->  GetN() )   {
        cout << "[ERROR] Systematics and Statistics do not have the same number of points! Skipping this one..." << endl;
        return nullptr;
    }
    TGraphAsymmErrors  *fResult =   new TGraphAsymmErrors(*fToBeCorrected);
    for ( Int_t iHisto = 1; iHisto <= fNPoints; iHisto++ )    {
        auto    fTargetY        =   fToBeCorrected->GetPointY(iHisto-1);
        auto    fEfficnY        =   fEfficiency->GetPointY(iHisto-1);
        auto    fTargetYlow     =   fToBeCorrected->GetErrorYlow(iHisto-1);
        auto    fEfficnYlow     =   fEfficiency->GetErrorYlow(iHisto-1);
        auto    fTargetYhigh    =   fToBeCorrected->GetErrorYlow(iHisto-1);
        auto    fEfficnYhigh    =   fEfficiency->GetErrorYlow(iHisto-1);
        fResult->SetPointY      (iHisto-1,  fTargetY*fEfficnY*fScale);
        fResult->SetPointEYlow  (iHisto-1,  sqrt(fTargetYlow*fTargetYlow+fEfficnYlow*fEfficnYlow));
        fResult->SetPointEYhigh (iHisto-1,  sqrt(fTargetYhigh*fTargetYhigh+fEfficnYhigh*fEfficnYhigh));
    }
    return fResult;
}
//
//_____________________________________________________________________________
//
//  --  --  Specific purpose histogram generation  --  --  //
//
//_____________________________________________________________________________
//>>    LEGACY, TO BE CHECKED AGAIN
//
TGraphAsymmErrors*
fTH1_to_TGAsymmErrors
 ( TH1*  hTarget )    {
    TGraphAsymmErrors  *fResult                     =   new TGraphAsymmErrors();
    Int_t       fNBins  =   hTarget->GetNbinsX();
    for ( Int_t iBin = 1; iBin <= fNBins; iBin++ )  {
        auto    fXValue =   hTarget         ->GetBinCenter(iBin);
        auto    fYValue =   hTarget         ->GetBinContent(iBin);
        auto    fXError =   0.5*( hTarget->GetBinLowEdge(iBin+1) - hTarget->GetBinLowEdge(iBin) );
        auto    fYError =   hTarget         ->GetBinError(iBin);
        fResult         ->  SetPoint        (iBin-1,fXValue,fYValue);
        fResult         ->  SetPointError   (iBin-1,fXError,fXError,fYError,fYError);
    }
    fResult             ->SetMaximum(hTarget->GetMaximum());
    fResult             ->SetMinimum(hTarget->GetMinimum());
    return fResult;
}
//
TGraphAsymmErrors**
fTH2_to_TGAsymmErrors
 ( TH2*  hTarget )    {
    Int_t       fNBinsY =   hTarget->GetNbinsY();
    TGraphAsymmErrors **fResult                     =   new TGraphAsymmErrors*[fNBinsY];
    for ( Int_t iBin = 1; iBin <= fNBinsY; iBin++ )  {
        fResult[iBin-1] =   fTH1_to_TGAsymmErrors(hTarget ->  ProjectionX(Form("%i",iBin),iBin,iBin));
    }
    return fResult;
}
//
TCanvas                *uPlotReferenceValue         ( TGraphAsymmErrors*  hMeasured,    Float_t fReference,   Float_t fRefError, TString fLabel = "PDG Value" )    {
    TCanvas    *fResult         =   new TCanvas();
    //
    auto        fNPoints        =   hMeasured->GetN();
    auto        fXLow           =   hMeasured->GetPointX(0);
    auto        fXLowErr        =   hMeasured->GetErrorXlow(0);
    auto        fXHig           =   hMeasured->GetPointX(fNPoints-1);
    auto        fXHigErr        =   hMeasured->GetErrorXhigh(fNPoints-1);
                hMeasured       ->  SetMaximum(max((float)(hMeasured->GetMaximum()+0.25*(hMeasured->GetMaximum()-hMeasured->GetMinimum())), fReference+fRefError*5));
                hMeasured       ->  SetMinimum(min((float)(hMeasured->GetMinimum()-0.25*(hMeasured->GetMaximum()-hMeasured->GetMinimum())), fReference-fRefError*5));
    //
                hMeasured       ->  Draw("APE");
    //
    TLegend    *fLegend         =   new TLegend(0.6,0.9,0.9,0.75);
                fLegend         ->  SetLineColorAlpha(1,0.);
                fLegend         ->  SetFillColorAlpha(1,0.);
                fLegend         ->  SetNColumns(2);
    //
    TH1F      **fPlotReference  =   new TH1F   *[7];
    for ( Int_t iTer = 0; iTer < 7; iTer ++ )   {
        fPlotReference[iTer]    =   new TH1F    (Form("%i",iTer),Form("%i",iTer),1,(fXLow-fXLowErr*1.5),(fXHig+fXHigErr*1.5));
        fPlotReference[iTer]    ->  SetBinContent(1, fReference + ( iTer -3 )*fRefError );
        fPlotReference[iTer]    ->  SetLineStyle(10);
        fPlotReference[iTer]    ->  SetLineWidth(3);
        fPlotReference[iTer]    ->  SetLineColorAlpha(920+(4-fabs(iTer-3)),0.75);
        fPlotReference[iTer]    ->  Draw("same ][");
    }
    fPlotReference[3]           ->  SetLineColorAlpha(2,0.75);
    fPlotReference[3]           ->  SetLineWidth(5);
    //
    fLegend                     ->  AddEntry( fPlotReference[3],    fLabel.Data(),  "L");
    fLegend                     ->  AddEntry( fPlotReference[2],    "#pm 1 #sigma", "L");
    fLegend                     ->  AddEntry( fPlotReference[1],    "#pm 2 #sigma", "L");
    fLegend                     ->  AddEntry( fPlotReference[0],    "#pm 3 #sigma", "L");
    fLegend                     ->  Draw("same");
    //
    return fResult;
}
TCanvas                *uPlotReferenceValue         ( TH1*  hMeasured,    Float_t fReference,   Float_t fRefError, TString fLabel = "PDG Value"  )    {
    return uPlotReferenceValue( fTH1_to_TGAsymmErrors(hMeasured), fReference, fRefError, fLabel );
}
//
//
TCanvas*
uPlotEfficiencies
 ( std::vector<TH1F*> hTarget, std::vector<TString> fLegend = {} )  {
    TCanvas        *cDrawEff    =   new TCanvas("","",1200,1200);
    gStyle          ->  SetOptStat(0);
    gPad            ->  SetLogx();
    gPad            ->  SetGridy();
    //
    TLegend*        lEfficiencies   =   new TLegend(0.625,0.88,0.88,0.7);
    lEfficiencies   ->  SetNColumns(2);
    lEfficiencies   ->  SetFillColorAlpha(0.,0.);
    lEfficiencies   ->  SetLineColorAlpha(0.,0.);
    //
    auto iTer = 0;
    for ( auto k1D_Eff : hTarget )  {
        uSetHisto( k1D_Eff, "EFF 1D" );
        if ( iTer != 0 )    k1D_Eff ->  SetMarkerStyle ( uGetMarker(4) );
        k1D_Eff ->  SetMarkerColor ( uGetColor(iTer) );
        k1D_Eff ->  SetLineColor ( uGetColor(iTer) );
        k1D_Eff ->  Draw( "SAME" );
        if ( iTer+1 > fLegend.size() )          lEfficiencies->AddEntry( k1D_Eff, k1D_Eff->GetName(),   "EP" );
        else if ( !fLegend.at(iTer).IsNull() )  lEfficiencies->AddEntry( k1D_Eff, fLegend.at(iTer),     "EP" );
        else                                    lEfficiencies->AddEntry( k1D_Eff, k1D_Eff->GetName(),   "EP" );
        iTer++;
    }
    lEfficiencies->Draw("SAME");
    //
    return cDrawEff;
}
//
*/
