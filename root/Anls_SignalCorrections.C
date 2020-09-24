// File for 1-Dimensional Analysis:
// !TODO: All Set!


#include "../inc/SetValues.h"
#include "../inc/SetFunctions.h"
#include "RooMsgService.h"

void Anls_SignalCorrections ( bool fSilent = false )
{
    //---------------------//
    //  Setting up input   //-------------------------------------------------------------------------------
    //---------------------//
    
    //-// OPTIONS
    
    // Silencing warnings for smoother
    if ( fSilent )
    {
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
        RooMsgService::instance().setSilentMode(fSilent);
    }
    
    // Opening Data and Efficiencies File
    TFile*  insFile_DT              =   new TFile   (fFitResults);
    TFile*  insFile_EF              =   new TFile   (fEfficiHist);
    TFile*  insCheck_               =   new TFile   ("./result/HEPData-ins1182213-v1-root.root");
    TFile*  insFile_MC_P6           =   new TFile   ("./result/MCTruth_Pythia6.root");
    TFile*  insFile_MC_P8           =   new TFile   ("./result/MCTruth_Pythia8.root");
    
    // Fit Variables for Roofit
    RooRealVar  xTransverseMom1D("xTransverseMom1D","xTransverseMom1D", fMinPT1D,fMaxPT1D);
    RooRealVar  xTransverseMom2D("xTransverseMom2D","xTransverseMom2D", fMinPT2D,fMaxPT2D);
    RooRealVar  yTransverseMom2D("yTransverseMom2D","yTransverseMom2D", fMinPT2D,fMaxPT2D);
    
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
    
    
    /*------------*/
    /*  ANALYSIS  */
    /*------------*/
    
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
    cout << "1D YIELD:" << endl;
    cout << "YIELD: " << fResults1D[0]    << "+-" <<  fResults1D[1]   << "+" <<   fResults1D[2]+kEventEfficienERP*fResults1D[0]  << "-" <<   fResults1D[2]+kEventEfficienERM*fResults1D[0]  << endl;
    cout << "MeanP: " << fResults1D[3]    << "+-" <<  fResults1D[4]   << endl;
    cout << endl;
    cout << "2DX YIELD:" << endl;
    cout << "YIELD: " << fResults2D[0][0][0]    << "+-" <<  fResults2D[0][0][1]   << "+" <<   fResults2D[0][0][2] + fResults2D[0][0][0]*kEventEfficienERP   << "-" <<   fResults2D[0][0][2] + fResults2D[0][0][0]*kEventEfficienERM   << endl;
    cout << "MeanP: " << fResults2D[0][0][3]    << "+-" <<  fResults2D[0][0][4]   << endl;
    cout << endl;
    cout << "2DY YIELD:" << endl;
    cout << "YIELD: " << fResults2D[0][1][0]    << "+-" <<  fResults2D[0][1][1]   << "+" <<   fResults2D[0][1][2] + fResults2D[0][1][0]*kEventEfficienERP   << "-" <<   fResults2D[0][1][2] + fResults2D[0][1][0]*kEventEfficienERM   << endl;
    cout << "MeanP: " << fResults2D[0][1][3]    << "+-" <<  fResults2D[0][1][4]   << endl;
    cout << endl;
    cout << endl;
    cout << endl;
    auto valy   =   fResults2D[0][0][0]/(fResults1D[0]*fResults1D[0]);
    auto valx   =   fResults2D[0][1][0]/(fResults1D[0]*fResults1D[0]);
    auto errx   =   sqrt(fResults2D[0][0][1]*fResults2D[0][1][1]/(fResults2D[0][0][0]*fResults2D[0][0][0])+2*(fResults1D[1]*fResults1D[1])/(fResults1D[0]*fResults1D[0]));
    auto erry   =   sqrt(fResults2D[0][1][1]*fResults2D[0][1][1]/(fResults2D[0][1][0]*fResults2D[0][1][0])+2*(fResults1D[1]*fResults1D[1])/(fResults1D[0]*fResults1D[0]));
    auto ersx   =   sqrt((pow(2*kSyst_SigExtr1D+kSyst_SigExtr2D,2)+pow(kEventEfficienERP,2)));
    auto ersy   =   sqrt((pow(2*kSyst_SigExtr1D+kSyst_SigExtr2D,2)+pow(kEventEfficienERP,2)));
    auto er2x   =   sqrt((pow(2*kSyst_SigExtr1D+kSyst_SigExtr2D,2)+pow(kEventEfficienERM,2)));
    auto er2y   =   sqrt((pow(2*kSyst_SigExtr1D+kSyst_SigExtr2D,2)+pow(kEventEfficienERM,2)));
    cout << "2DX/1D^2:" << valx << "+-" << valx*errx << "+" << valx*ersx << "-" << valx*er2x << endl;
    cout << endl;
    cout << "2DY/1D^2:" << valy << "+-" << valy*erry << "+" << valy*ersy << "-" << valy*er2y << endl;
    cout << endl;
    cout << endl;
    cout << endl;
    
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
    TGraphAEGeneratorPT(fResults2D);
    
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
