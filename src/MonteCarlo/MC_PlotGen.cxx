#include "../../inc/AliAnalysisPhiPair.h"

void
MC_PlotGen_
( TString fFileName = "/Users/nrubini/Analysis/ALICE/Simulation/Pythia8PhiPair/FinalResult_29_7000.root", TString kProductionTag = "PythiaX29", TString fOption = "yield", Int_t nEventsCut = -1, TString kFolder = "_p_p__7TeV", Int_t iMode = 0, Int_t kEnergy = 5000 );

void
MC_PlotGen2
() {
    /*
    MC_PlotGen_("/Users/nrubini/Analysis/ALICE/Simulation/Pythia8PhiPair/5000/0/FinalResults_0_5000.root","Pythia8X0","yield",-1.,"_p_p__5TeV",0);
    MC_PlotGen_("/Users/nrubini/Analysis/ALICE/Simulation/Pythia8PhiPair/5000/8/FinalResults_8_5000.root","Pythia8X8","yield",-1.,"_p_p__5TeV",8);
    MC_PlotGen_("/Users/nrubini/Analysis/ALICE/Simulation/Pythia8PhiPair/5000/23/FinalResults_23_5000.root","Pythia8X23","yield",-1.,"_p_p__5TeV",23);
    MC_PlotGen_("/Users/nrubini/Analysis/ALICE/Simulation/Pythia8PhiPair/5000/24/FinalResults_24_5000.root","Pythia8X24","yield",-1.,"_p_p__5TeV",24);
    MC_PlotGen_("/Users/nrubini/Analysis/ALICE/Simulation/Pythia8PhiPair/5000/25/FinalResults_25_5000.root","Pythia8X25","yield",-1.,"_p_p__5TeV",25);
    MC_PlotGen_("/Users/nrubini/Analysis/ALICE/Simulation/Pythia8PhiPair/5000/26/FinalResults_26_5000.root","Pythia8X26","yield",-1.,"_p_p__5TeV",26);
    MC_PlotGen_("/Users/nrubini/Analysis/ALICE/Simulation/Pythia8PhiPair/5000/27/FinalResults_27_5000.root","Pythia8X27","yield",-1.,"_p_p__5TeV",27);
    MC_PlotGen_("/Users/nrubini/Analysis/ALICE/Simulation/Pythia8PhiPair/5000/28/FinalResults_28_5000.root","Pythia8X28","yield",-1.,"_p_p__5TeV",28);
    MC_PlotGen_("/Users/nrubini/Analysis/ALICE/Simulation/Pythia8PhiPair/5000/29/FinalResults_29_5000.root","Pythia8X29","yield",-1.,"_p_p__5TeV",29);
    MC_PlotGen_("/Users/nrubini/Analysis/ALICE/Simulation/Pythia8PhiPair/5000/29/FinalResults_29_5000.root","Pythia8X29","yield",-1.,"_p_p__5TeV",29);
    */
    MC_PlotGen_("/Users/nrubini/Simulation/Macros/Pythia8PhiPair/FinalResult_0_7000.root","Pythia8X0_","yield",-1.,"_p_p__7TeV",0,7000);
    MC_PlotGen_("/Users/nrubini/Simulation/Macros/Pythia8PhiPair/FinalResult_8_7000.root","Pythia8X8_","yield",-1.,"_p_p__7TeV",8,7000);
    MC_PlotGen_("/Users/nrubini/Simulation/Macros/Pythia8PhiPair/5000/0/FinalResults_0_5000.root","Pythia8X0_","yield",-1.,"_p_p__5TeV",0,5000);
    MC_PlotGen_("/Users/nrubini/Simulation/Macros/Pythia8PhiPair/5000/8/FinalResults_8_5000.root","Pythia8X8_","yield",-1.,"_p_p__5TeV",8,5000);
}

void
MC_PlotGen
() {
    TString kFolder = "/Users/nrubini/Simulation/Pythia8/DataSet/7000/";
    //MC_PlotGenTrue_( kFolder + TString("Pythia8_07TeV_Monash_2013.root"),       "Pythia8X_0_0","yield",-1,"_p_p__7TeV");
    //MC_PlotGenTrue_( kFolder + TString("Pythia8_Monash13_07TeV.root"),          "Pythia8X_0_1","yield",-1,"_p_p__7TeV");
    //MC_PlotGenTrue_( kFolder + TString("Pythia8_07TeV_Monash_2013_Ropes.root"), "Pythia8X_8_0","yield",-1,"_p_p__7TeV");
    MC_PlotGen_( kFolder + TString("AnalysisResults-5.root"),          "Pythia8X_8_1","yield",-1,"_p_p__7TeV",-1.);
}

void
MC_PlotGen_
( TString fFileName = "/Users/nrubini/Analysis/ALICE/Simulation/Pythia8PhiPair/FinalResult_29_7000.root", TString kProductionTag = "PythiaX29", TString fOption = "yield", Int_t nEventsCut = -1, TString kFolder = "_p_p__7TeV", Int_t iMode = 29, Int_t kEnergy = 5000 ) {
    // --- --- --- --- --- --- --- SET-UP --- --- --- --- --- --- --- --- --- --- ---
    //
    //  Check Correct Syntax
    if ( fFileName == "" )  {
        cout << "[ERROR] Must Specify an input root file" << endl;
        return;
    }
    //
    // --- INFO on Set-up variables
    if ( nEventsCut != -1 ) cout << "[INFO] Choosing to limit the datasample to " << nEventsCut << " events" <<endl;
    fChooseOption(fOption);
    //
    // --- Retrieving Event data
    TFile*      insFile_Data    =   new TFile   ( fFileName );
    //
    //
    auto kDebugLimit = 1.;
    //
    TTree*      TPhiEfficiency;
    TList*      fQCOutputList;
    TH1D*       fHEventCount;
    TH1D*       fQC_Event_Enum_E05;
    TH1D*       fQC_Event_Enum_V0M;
    Struct_MCParticle    evPhiEfficiency;
    if ( iMode >= 0 ) {
        // --- Retrieving TTree
        TPhiEfficiency  =   (TTree*)insFile_Data->Get( Form("Prt_E%i_M%i",kEnergy,iMode) );
        //
        // --- Retrieving Utility Histograms
        fHEventCount        =   (TH1D*) insFile_Data->Get( "fQC_Event_Enum_FLL" );
        fQC_Event_Enum_E05  =   (TH1D*) insFile_Data->Get( "fQC_Event_Enum_FLL" );
        fQC_Event_Enum_V0M  =   (TH1D*) insFile_Data->Get( "fQC_Event_Enum_FLL" );
        //
        // --- Setting the input datastructure
        fSetPhiCandidateMC(TPhiEfficiency,evPhiEfficiency);
    } else {
        // --- Retrieving TTree
        TPhiEfficiency  =   (TTree*)insFile_Data->Get( "PhiCandidate_" );
        //
        // --- Retrieving Utility Histograms
        fQCOutputList       =   (TList*)insFile_Data->Get( "fQCOutputList_PhiCount_STD" );
        fHEventCount        =   (TH1D*) fQCOutputList->FindObject( "fQC_Event_Enum_FLL" );
        fQC_Event_Enum_E05  =   (TH1D*) fQCOutputList->FindObject( "fQC_Event_Enum_E10" );
        fQC_Event_Enum_V0M  =   (TH1D*) fQCOutputList->FindObject( "fQC_Event_Enum_V0M" );
        //
        // --- Setting the input datastructure
        fSetPhiCandidateMC(TPhiEfficiency,evPhiEfficiency);
    }
    //
    // --- Setting the output datastructure
    fSetAllBins();
    auto    kUniformBinning_200MeV =   uUniformBinning( 200*MeV, 0., 50. );
    auto    kUniformBinning_400MeV =   uUniformBinning( 400*MeV, 0., 50. );
    //
    //  --- Result Histograms
    //  --- --- 1D Histograms
    TH1D*       h1D_Ntru_Fin;
    TH1D*       h1D_Ntru;
    //  --- --- 2D Histograms
    TH2D*       h2D_Ntru_Fin;
    TH2D*       h2D_Ntru_Fin_Nrm;
    TH2D*       h2D_Ntru;
    TProfile*   hMPT_NTru;
    TProfile*   hMPT_NTru_Fin;
    //  --- --- Combination Histograms
    TH1D*       h1D_Ntru_Mult;
    TH1D*       h2D_Ntru_Mult;
    TH1D*       hR1_Ntru_Mult;
    TH1D*       hR2_Ntru_Mult;
    TH1D*       hP1_Ntru_Mult;
    TH1D*       hP2_Ntru_Mult;
    TH1D*       hFullQuantities;
    TH1D*       hPhiCorrelation;
    TH1D*       hPhiCorrelation_Full;
    TH1D*       hPhiCorrelation_FBkg;
    TH1D*       hPhiCorrlRapSig;
    TH1D*       hPhiCorrlRapBkg;
    TH1D*       hProductionProb;
    TH1D*       hMultCorr;
    TH1D*       hMultUCor;
    TH2D*       hDPhiDyCh;
    TH2D*       hDPhiDyCh_Bkg;
    TH2D*       hDPhiDyCh_Tst;
    TH2D*       hDPhiDyCh_Eta;
    TH2D*       hDPhiDyCh_Phi;
    TH2D*       hDPhiDyCh_Eta_Bkg;
    TH2D*       hDPhiDyCh_Phi_Bkg;
    TH2D*       hDPhiDyCh_Rap_Eta_Bkg;
    TH2D*       hDPhiDyCh_Rap_Phi_Bkg;
    TH2D*       hDPhiDyCh_Rap_Eta;
    TH2D*       hDPhiDyCh_Rap_Phi;
    //
    hName           =   Form( "h1D_Ntru_Fin" );
    hTitle          =   Form( "True Kaons from phi" );
    h1D_Ntru_Fin    =   new TH1D ( hName, hTitle, kUniformBinning_200MeV.first , kUniformBinning_200MeV.second );
    SetAxis( h1D_Ntru_Fin, "PT 1D" );
    //
    hName           =   Form( "h2D_Ntru_Fin" );
    hTitle          =   Form( "True Kaons from phi" );
    h2D_Ntru_Fin    =   new TH2D ( hName, hTitle, kUniformBinning_200MeV.first , kUniformBinning_200MeV.second, kUniformBinning_200MeV.first , kUniformBinning_200MeV.second );
    SetAxis( h2D_Ntru_Fin, "PT 2D" );
    //
    hName           =   Form( "h1D_Ntru" );
    hTitle          =   Form( "True Kaons from phi" );
    h1D_Ntru        =   new TH1D ( hName, hTitle, nBinPT1D+1, fArrPT1D_Comp );
    SetAxis( h1D_Ntru, "PT 1D" );
    //
    hName           =   Form( "h2D_Ntru" );
    hTitle          =   Form( "True Kaons from phi" );
    h2D_Ntru        =   new TH2D ( hName, hTitle, nBinPT2D+1, fArrPT2D_Comp, nBinPT2D+1, fArrPT2D_Comp );
    SetAxis( h2D_Ntru, "PT 2D" );
    //
    hName           =   Form( "h2D_Ntru_Fin_Nrm" );
    hTitle          =   Form( "True Kaons from phi" );
    h2D_Ntru_Fin_Nrm=   new TH2D ( hName, hTitle, nBinPT2D+1, fArrPT2D_Comp, kUniformBinning_200MeV.first , kUniformBinning_200MeV.second );
    SetAxis( h2D_Ntru_Fin_Nrm, "PT 2D" );
    //
    hName           =   Form( "h1D_Ntru_Mult" );
    hTitle          =   Form( "True Kaons from phi in mult E10" );
    h1D_Ntru_Mult   =   new TH1D ( hName, hTitle, 5000,   0., 5000.);
    SetAxis( h1D_Ntru_Mult, "PT 1D" );
    //
    hName           =   Form( "h2D_Ntru_Mult" );
    hTitle          =   Form( "True Kaons from phi in mult E10" );
    h2D_Ntru_Mult   =   new TH1D ( hName, hTitle, 5000,   0., 5000.);
    SetAxis( h2D_Ntru_Mult, "PT 2D" );
    //
    hName           =   Form( "hR1_Ntru_Mult" );
    hTitle          =   Form( "True Kaons from phi in mult E10" );
    hR1_Ntru_Mult   =   new TH1D ( hName, hTitle, 5000,   0., 5000.);
    SetAxis( hR1_Ntru_Mult, "PT 1D" );
    //
    hName           =   Form( "hR2_Ntru_Mult" );
    hTitle          =   Form( "True Kaons from phi in mult E10" );
    hR2_Ntru_Mult   =   new TH1D ( hName, hTitle, 5000,   0., 5000.);
    SetAxis( hR2_Ntru_Mult, "PT 1D" );
    //
    hName           =   Form( "hP1_Ntru_Mult" );
    hTitle          =   Form( "True Kaons from phi in mult E10" );
    hP1_Ntru_Mult   =   new TH1D ( hName, hTitle, 5000,   0., 5000.);
    SetAxis( hP1_Ntru_Mult, "PT 1D" );
    //
    hName           =   Form( "hP2_Ntru_Mult" );
    hTitle          =   Form( "True Kaons from phi in mult E10" );
    hP2_Ntru_Mult   =   new TH1D ( hName, hTitle, 5000,   0., 5000.);
    SetAxis( hP2_Ntru_Mult, "PT 1D" );
    //
    hName           =   Form( "hMPT_NTru" );
    hTitle          =   Form( "True Kaons from phi in hMPT_NTru" );
    hMPT_NTru   =   new TProfile ( hName, hTitle, nBinPT2D+1, fArrPT2D_Comp );
    SetAxis( hMPT_NTru, "PT 1D" );
    //
    hName           =   Form( "hMPT_NTru_Fin" );
    hTitle          =   Form( "True Kaons from phi in hMPT_NTru" );
    hMPT_NTru_Fin   =   new TProfile ( hName, hTitle, kUniformBinning_400MeV.first, kUniformBinning_400MeV.second );
    SetAxis( hMPT_NTru_Fin, "PT 1D" );
    //
    hName           =   Form( "hFullQuantities" );
    hTitle          =   Form( "hFullQuantities" );
    hFullQuantities =   new TH1D ( hName, hTitle, 6, +0.5, +6.5 );
    //
    hName           =   Form( "hProductionProb" );
    hTitle          =   Form( "hProductionProb" );
    hProductionProb =   new TH1D ( hName, hTitle, 7, -0.5, +6.5 );
    //
    hName           =   Form( "hPhiCorrelation" );
    hTitle          =   Form( "hPhiCorrelation" );
    hPhiCorrelation =   new TH1D ( hName, hTitle, nBinCrPh,  fArrCrPh);
    SetAxis( hPhiCorrelation, "PT 1D" );
    //
    hName           =   Form( "hPhiCorrelation_Full" );
    hTitle          =   Form( "hPhiCorrelation_Full" );
    hPhiCorrelation_Full =   new TH1D ( hName, hTitle, 100, -4., 4. );
    SetAxis( hPhiCorrelation_Full, "PT 1D" );
    //
    hName           =   Form( "hPhiCorrelation_FBkg" );
    hTitle          =   Form( "hPhiCorrelation_FBkg" );
    hPhiCorrelation_FBkg =   new TH1D ( hName, hTitle, 100, -4., 4. );
    SetAxis( hPhiCorrelation_FBkg, "PT 1D" );
    //
    hName           =   Form( "hPhiCorrlRapSig" );
    hTitle          =   Form( "hPhiCorrlRapSig" );
    hPhiCorrlRapSig =   new TH1D ( hName, hTitle, 100, -1., 1. );
    SetAxis( hPhiCorrlRapSig, "PT 1D" );
    //
    hName           =   Form( "hPhiCorrlRapBkg" );
    hTitle          =   Form( "hPhiCorrlRapBkg" );
    hPhiCorrlRapBkg =   new TH1D ( hName, hTitle, 100, -1., 1. );
    SetAxis( hPhiCorrlRapBkg, "PT 1D" );
    //
    hName           =   Form( "hMultCorr" );
    hTitle          =   Form( "True Kaons from phi" );
    hMultCorr    =   new TH1D ( hName, hTitle, kUniformBinning_200MeV.first , kUniformBinning_200MeV.second );
    SetAxis( hMultCorr, "PT 1D" );
    //
    hName           =   Form( "hMultUCor" );
    hTitle          =   Form( "True Kaons from phi" );
    hMultUCor    =   new TH1D ( hName, hTitle, kUniformBinning_200MeV.first , kUniformBinning_200MeV.second );
    SetAxis( hMultUCor, "PT 1D" );
    //
    hName           =   Form( "hDPhiDyCh" );
    hTitle          =   Form( "True Kaons from phi" );
    hDPhiDyCh    =   new TH2D ( hName, hTitle, 100, -105, 255, 100, -1.6, 1.6 );
    hDPhiDyCh   ->  GetXaxis()  ->  SetTitle("#Delta_{#varphi} (deg.)");
    hDPhiDyCh   ->  GetYaxis()  ->  SetTitle("#Delta_{#eta}");
    //
    hName           =   Form( "hDPhiDyCh_Bkg" );
    hTitle          =   Form( "True Kaons from phi" );
    hDPhiDyCh_Bkg    =   new TH2D ( hName, hTitle, 100, -105, 255, 100, -1.6, 1.6 );
    hDPhiDyCh_Bkg   ->  GetXaxis()  ->  SetTitle("#Delta_{#varphi} (deg.)");
    hDPhiDyCh_Bkg   ->  GetYaxis()  ->  SetTitle("#Delta_{#eta}");
    //
    hName           =   Form( "hDPhiDyCh_Rap_Phi" );
    hTitle          =   Form( "True Kaons from phi" );
    hDPhiDyCh_Rap_Phi    =   new TH2D ( hName, hTitle, 100, -105, 255, 100, -1.0, 1.0 );
    hDPhiDyCh_Rap_Phi   ->  GetXaxis()  ->  SetTitle("#Delta_{#varphi} (deg.)");
    hDPhiDyCh_Rap_Phi   ->  GetYaxis()  ->  SetTitle("#Delta_{y}");
    //
    hName           =   Form( "hDPhiDyCh_Rap_Phi_Bkg" );
    hTitle          =   Form( "True Kaons from phi" );
    hDPhiDyCh_Rap_Phi_Bkg    =   new TH2D ( hName, hTitle, 100, -105, 255, 100, -1.0, 1.0 );
    hDPhiDyCh_Rap_Phi_Bkg   ->  GetXaxis()  ->  SetTitle("#Delta_{#varphi} (deg.)");
    hDPhiDyCh_Rap_Phi_Bkg   ->  GetYaxis()  ->  SetTitle("#Delta_{y}");
    //
    hName           =   Form( "hDPhiDyCh_Rap_Eta" );
    hTitle          =   Form( "True Kaons from phi" );
    hDPhiDyCh_Rap_Eta    =   new TH2D ( hName, hTitle, 100, -1.6, 1.6, 100, -1.0, 1.0 );
    hDPhiDyCh_Rap_Eta   ->  GetXaxis()  ->  SetTitle("#Delta_{#eta}");
    hDPhiDyCh_Rap_Eta   ->  GetYaxis()  ->  SetTitle("#Delta_{y}");
    //
    hName           =   Form( "hDPhiDyCh_Rap_Eta_Bkg" );
    hTitle          =   Form( "True Kaons from phi" );
    hDPhiDyCh_Rap_Eta_Bkg    =   new TH2D ( hName, hTitle, 100, -1.6, 1.6, 100, -1.0, 1.0 );
    hDPhiDyCh_Rap_Eta_Bkg   ->  GetXaxis()  ->  SetTitle("#Delta_{#eta}");
    hDPhiDyCh_Rap_Eta_Bkg   ->  GetYaxis()  ->  SetTitle("#Delta_{y}");
    //
    //hDPhiDyCh_Eta
    //
    hName           =   Form( "hDPhiDyCh_Eta" );
    hTitle          =   Form( "hDPhiDyCh_Eta" );
    hDPhiDyCh_Eta    =   new TH2D ( hName, hTitle, 100, -0.8, 0.8, 100, -0.8, 0.8 );
    hDPhiDyCh_Eta   ->  GetXaxis()  ->  SetTitle("#eta_{1}");
    hDPhiDyCh_Eta   ->  GetYaxis()  ->  SetTitle("#eta_{2}");
    //
    hName           =   Form( "hDPhiDyCh_Phi" );
    hTitle          =   Form( "hDPhiDyCh_Phi" );
    hDPhiDyCh_Phi    =   new TH2D ( hName, hTitle, 100, -180, 180, 100, -180, 180 );
    hDPhiDyCh_Phi   ->  GetXaxis()  ->  SetTitle("#varphi_{1}");
    hDPhiDyCh_Phi   ->  GetYaxis()  ->  SetTitle("#varphi_{2}");
    //
    hName           =   Form( "hDPhiDyCh_Eta_Bkg" );
    hTitle          =   Form( "hDPhiDyCh_Eta_Bkg" );
    hDPhiDyCh_Eta_Bkg    =   new TH2D ( hName, hTitle, 100, -0.8, 0.8, 100, -0.8, 0.8 );
    hDPhiDyCh_Eta_Bkg   ->  GetXaxis()  ->  SetTitle("#eta_{1}");
    hDPhiDyCh_Eta_Bkg   ->  GetYaxis()  ->  SetTitle("#eta_{2}");
    //
    hName           =   Form( "hDPhiDyCh_Phi_Bkg" );
    hTitle          =   Form( "hDPhiDyCh_Phi_Bkg" );
    hDPhiDyCh_Phi_Bkg    =   new TH2D ( hName, hTitle, 100, -180, 180, 100, -180, 180 );
    hDPhiDyCh_Phi_Bkg   ->  GetXaxis()  ->  SetTitle("#varphi_{1}");
    hDPhiDyCh_Phi_Bkg   ->  GetYaxis()  ->  SetTitle("#varphi_{2}");
    //
    hName           =   Form( "hPhiCorrelation_IM_Sig" );
    hTitle          =   Form( "hPhiCorrelation_IM_Sig" );
    TH1D *hPhiCorrelation_IM_Sig =   new TH1D ( hName, hTitle, 1800,2.,20. );
    SetAxis( hPhiCorrelation_IM_Sig, "IM 1D" );
    hPhiCorrelation_IM_Sig->GetXaxis()->SetTitleOffset(0.9);
    hPhiCorrelation_IM_Sig->GetXaxis()->SetTitle("M_{#phi#phi} (GeV/#it{c}^{2})");
    //
    hName           =   Form( "hPhiCorrelation_IM_Bkg" );
    hTitle          =   Form( "hPhiCorrelation_IM_Bkg" );
    TH1D *hPhiCorrelation_IM_Bkg =   new TH1D ( hName, hTitle, 1800,2.,20. );
    SetAxis( hPhiCorrelation_IM_Bkg, "IM 1D" );
    hPhiCorrelation_IM_Bkg->GetXaxis()->SetTitleOffset(0.9);
    hPhiCorrelation_IM_Bkg->GetXaxis()->SetTitle("M_{#phi#phi} (GeV/#it{c}^{2})");
    //
    hName           =   Form( "hPhiCorrelation_IM_Rat" );
    hTitle          =   Form( "hPhiCorrelation_IM_Rat" );
    TH1D *hPhiCorrelation_IM_Rat =   new TH1D ( hName, hTitle, 1800,2.,20. );
    SetAxis( hPhiCorrelation_IM_Rat, "IM 1D" );
    hPhiCorrelation_IM_Rat->GetXaxis()->SetTitleOffset(0.9);
    hPhiCorrelation_IM_Rat->GetXaxis()->SetTitle("M_{#phi#phi} (GeV/#it{c}^{2})");
    //
    //  --- Applying limitation sample cuts
    auto nEvents = (!TPhiEfficiency) ? 0 : ( nEventsCut == -1.? TPhiEfficiency->GetEntries() : nEventsCut);
    nEvents *= kDebugLimit;
    if ( nEvents > 0 )  fStartTimer("Comparison Plots Production");
    //
    Struct_MCParticle fCurrent_MCTrue;
    Struct_MCParticle fPrevious_MCTrue;
    //
    fCurrent_MCTrue.nPart   =   0;
    fPrevious_MCTrue.nPart  =   0;
    for ( Int_t iEvent = 0; iEvent < nEvents; iEvent++ )    {
        TPhiEfficiency->GetEntry(iEvent);
        fPrintLoopTimer("Comparison Plots Production",iEvent,nEvents,kPrintIntervalPP);
        //
        //  --- Loop over candidates
        fCurrent_MCTrue.nPart    =   0;
        for ( Int_t iPhi = 0; iPhi < evPhiEfficiency.nPart; iPhi++ )  {
            TLorentzVector  kTLVUtility;   // TODO: Update to new lorentzvector class
            kTLVUtility.SetXYZM(evPhiEfficiency.Px[iPhi],evPhiEfficiency.Py[iPhi],evPhiEfficiency.Pz[iPhi],kPhiMesonMass_);
            //  --- Load Useful Variables for histograms filling
            fCurrent_MCTrue.pT          [ fCurrent_MCTrue.nPart ]    =   kTLVUtility.Pt();
            fCurrent_MCTrue.Px          [ fCurrent_MCTrue.nPart ]    =   kTLVUtility.Px();
            fCurrent_MCTrue.Py          [ fCurrent_MCTrue.nPart ]    =   kTLVUtility.Py();
            fCurrent_MCTrue.Pz          [ fCurrent_MCTrue.nPart ]    =   kTLVUtility.Pz();
            fCurrent_MCTrue.Eta         [ fCurrent_MCTrue.nPart ]    =   kTLVUtility.Eta();
            fCurrent_MCTrue.Rap         [ fCurrent_MCTrue.nPart ]    =   kTLVUtility.Rapidity();
            fCurrent_MCTrue.Phi         [ fCurrent_MCTrue.nPart ]    =   kTLVUtility.Phi() * 360 / ( TMath::Pi() * 2 );
            fCurrent_MCTrue.iPT1D       [ fCurrent_MCTrue.nPart ]    =   fGetBinPT1D( kTLVUtility.Pt() );
            fCurrent_MCTrue.iPT2D       [ fCurrent_MCTrue.nPart ]    =   fGetBinPT2D( kTLVUtility.Pt() );
            if ( is_pp_anl )    fCurrent_MCTrue.kHasRap     [ fCurrent_MCTrue.nPart ]    =   fabs( kTLVUtility.Rapidity() ) < 0.5;
            if ( is_pb_anl )    fCurrent_MCTrue.kHasRap     [ fCurrent_MCTrue.nPart ]    =   ( kTLVUtility.Rapidity() ) < 0.035 && ( kTLVUtility.Rapidity() ) > -0.465;
            fCurrent_MCTrue.nPart++;
        }
        //
        hProductionProb ->Fill( fCurrent_MCTrue.nPart );
        //
        //  --- Discarding if no valid candidates are found
        if ( fCurrent_MCTrue.nPart < 1 ) continue;
        //
        //  --- Event variables
        fCurrent_MCTrue.IsMB            =   true;//( evPhiEfficiency.TrueEventMask == 0 );
        fCurrent_MCTrue.Eta_05          =   evPhiEfficiency.Eta_05;
        //fCurrent_MCTrue.iMult           =   -1;//fGetBinMult(evPhiEfficiency.Eta_10);
        //fCurrent_MCTrue.kHasMult        =   -1;//fCurrent_MCTrue.iMult != -1;
        //
        for ( Int_t iPhi = 0; iPhi < fCurrent_MCTrue.nPart; iPhi++ )  {
            //
            if ( fCurrent_MCTrue.kHasRap[iPhi] && fCurrent_MCTrue.IsMB )  { //  --- Mid-Rapidity Analyses
                if ( kDoYield ) {                                           //  --- YIELD ANALYSIS
                    h1D_Ntru        ->  Fill( fCurrent_MCTrue.pT[iPhi] );
                    h1D_Ntru_Fin    ->  Fill( fCurrent_MCTrue.pT[iPhi] );
                    hFullQuantities ->  Fill( 1 );
                    h1D_Ntru_Mult   ->  Fill( fCurrent_MCTrue.Eta_05 );
                }
                for ( Int_t jPhi = 0; jPhi < fPrevious_MCTrue.nPart; jPhi++ )  {
                    //
                    if ( fPrevious_MCTrue.kHasRap[jPhi] && fPrevious_MCTrue.IsMB )  { //  --- Mid-Rapidity Analyses
                        if ( kDoYield ) {                                           //  --- YIELD ANALYSIS
                            hDPhiDyCh_Bkg   -> Fill( fGetDltCrPh( fCurrent_MCTrue.Phi[iPhi] - fPrevious_MCTrue.Phi[jPhi] ), fCurrent_MCTrue.Eta[iPhi] - fPrevious_MCTrue.Eta[jPhi] );
                            hDPhiDyCh_Eta_Bkg   -> Fill( fCurrent_MCTrue.Eta[iPhi], fPrevious_MCTrue.Eta[jPhi] );
                            hDPhiDyCh_Phi_Bkg   -> Fill( fCurrent_MCTrue.Phi[iPhi], fPrevious_MCTrue.Phi[jPhi] );
                            hDPhiDyCh_Rap_Eta_Bkg   -> Fill( fCurrent_MCTrue.Eta[iPhi] - fPrevious_MCTrue.Eta[jPhi], fCurrent_MCTrue.Rap[iPhi] - fPrevious_MCTrue.Rap[jPhi] );
                            hDPhiDyCh_Rap_Phi_Bkg   -> Fill( fGetDltCrPh( fCurrent_MCTrue.Phi[iPhi] - fPrevious_MCTrue.Phi[jPhi] ), fCurrent_MCTrue.Rap[iPhi] - fPrevious_MCTrue.Rap[jPhi] );
                            hPhiCorrlRapBkg ->  Fill( fCurrent_MCTrue.Rap[iPhi] - fPrevious_MCTrue.Rap[jPhi] );
                            TLorentzVector  kTLVUtility_1;
                            TLorentzVector  kTLVUtility_2;
                            TLorentzVector  kTLVUtility_3;
                            kTLVUtility_1.SetXYZM(fCurrent_MCTrue.Px[iPhi],fCurrent_MCTrue.Py[iPhi],fCurrent_MCTrue.Pz[iPhi],kPhiMesonMass_);
                            kTLVUtility_2.SetXYZM(fCurrent_MCTrue.Px[jPhi],fCurrent_MCTrue.Py[jPhi],fCurrent_MCTrue.Pz[jPhi],kPhiMesonMass_);
                            kTLVUtility_3   =   kTLVUtility_1 + kTLVUtility_2;
                            hPhiCorrelation_IM_Bkg->Fill( kTLVUtility_3.Mag() );
                        }
                    }
                    hPhiCorrelation_FBkg->Fill( fCurrent_MCTrue.Rap[iPhi] - fPrevious_MCTrue.Rap[jPhi] );
                }
            }
            //
            //  --- Speed-up protection
            if ( fCurrent_MCTrue.nPart < 2 ) continue;
            for ( Int_t jPhi = iPhi+1; jPhi < fCurrent_MCTrue.nPart; jPhi++ )  {
                //
                if ( fCurrent_MCTrue.kHasRap[iPhi] && fCurrent_MCTrue.kHasRap[jPhi] && fCurrent_MCTrue.IsMB )  {    //  --- Mid-Rapidity Analyses
                    //  --- TEST
                    auto    kDeltaPhi   =   ( fCurrent_MCTrue.Phi[iPhi] - fCurrent_MCTrue.Phi[jPhi]  );
                    kDeltaPhi = kDeltaPhi < -180 ? kDeltaPhi + 360 : kDeltaPhi > 180 ? kDeltaPhi -360 : kDeltaPhi;
                    //
                    if ( kDoYield ) {   //  --- YIELD ANALYSIS
                        hPhiCorrelation ->  Fill( fGetDltCrPh( fCurrent_MCTrue.Phi[iPhi] - fCurrent_MCTrue.Phi[jPhi] ) );
                        hPhiCorrlRapSig ->  Fill( fCurrent_MCTrue.Rap[iPhi] - fCurrent_MCTrue.Rap[jPhi] );
                        h2D_Ntru        ->  Fill( fCurrent_MCTrue.pT[iPhi], fCurrent_MCTrue.pT[jPhi] );
                        h2D_Ntru_Fin    ->  Fill( fCurrent_MCTrue.pT[iPhi], fCurrent_MCTrue.pT[jPhi] );
                        h2D_Ntru_Fin_Nrm->  Fill( fCurrent_MCTrue.pT[iPhi], fCurrent_MCTrue.pT[jPhi] );
                        hMPT_NTru       ->  Fill( fCurrent_MCTrue.pT[iPhi], fCurrent_MCTrue.pT[jPhi] );
                        hMPT_NTru_Fin   ->  Fill( fCurrent_MCTrue.pT[iPhi], fCurrent_MCTrue.pT[jPhi] );
                        h2D_Ntru_Mult   ->  Fill( fCurrent_MCTrue.Eta_10 );
                        hFullQuantities ->  Fill( 2 );
                        hDPhiDyCh       ->  Fill( fGetDltCrPh( fCurrent_MCTrue.Phi[iPhi] - fCurrent_MCTrue.Phi[jPhi] ), fCurrent_MCTrue.Eta[iPhi] - fCurrent_MCTrue.Eta[jPhi] );
                        hDPhiDyCh_Rap_Phi ->  Fill( fGetDltCrPh( fCurrent_MCTrue.Phi[iPhi] - fCurrent_MCTrue.Phi[jPhi] ), fCurrent_MCTrue.Rap[iPhi] - fCurrent_MCTrue.Rap[jPhi] );
                        hDPhiDyCh_Rap_Eta ->  Fill( fCurrent_MCTrue.Eta[iPhi] - fCurrent_MCTrue.Eta[jPhi], fCurrent_MCTrue.Rap[iPhi] - fCurrent_MCTrue.Rap[jPhi] );
                        hDPhiDyCh_Eta   ->  Fill( fCurrent_MCTrue.Eta[iPhi], fCurrent_MCTrue.Eta[jPhi] );
                        hDPhiDyCh_Phi   ->  Fill( fCurrent_MCTrue.Phi[iPhi], fCurrent_MCTrue.Phi[jPhi] );
                        //if ( fCurrent_MCTrue.nPart ) hPhiCorrelatio2 ->  Fill( fGetDltCrPh( fCurrent_MCTrue.Phi[iPhi] - fCurrent_MCTrue.Phi[jPhi] ) );
                        //hMultCorr       ->  Fill(  );
                        TLorentzVector  kTLVUtility_1;
                        TLorentzVector  kTLVUtility_2;
                        TLorentzVector  kTLVUtility_3;
                        kTLVUtility_1.SetXYZM(fCurrent_MCTrue.Px[iPhi],fCurrent_MCTrue.Py[iPhi],fCurrent_MCTrue.Pz[iPhi],kPhiMesonMass_);
                        kTLVUtility_2.SetXYZM(fCurrent_MCTrue.Px[jPhi],fCurrent_MCTrue.Py[jPhi],fCurrent_MCTrue.Pz[jPhi],kPhiMesonMass_);
                        kTLVUtility_3   =   kTLVUtility_1 + kTLVUtility_2;
                        hPhiCorrelation_IM_Sig->Fill( kTLVUtility_3.Mag() );
                   }
                }
                hPhiCorrelation_Full->Fill( fCurrent_MCTrue.Rap[iPhi] - fCurrent_MCTrue.Rap[jPhi] );
            }
        }
        fPrevious_MCTrue = fCurrent_MCTrue;
        //
    }
    if ( nEvents > 0 )  fStopTimer("Comparison Plots Production");
    //
    // --- --- --- --- --- --- --- OUTPUT --- --- --- --- --- --- --- --- --- --- ---
    //
    auto    kEventNormalisation =   fHEventCount->GetBinContent(2);
    h1D_Ntru_Fin    ->  Scale( 1., "width" );
    h1D_Ntru_Fin    ->  Scale( 1./kEventNormalisation );
    h2D_Ntru_Fin    ->  Scale( 1., "width" );
    h2D_Ntru_Fin    ->  Scale( 1./kEventNormalisation );
    h2D_Ntru_Fin_Nrm->  Scale( 1., "width" );
    h2D_Ntru_Fin_Nrm->  Scale( 1./kEventNormalisation );
    h1D_Ntru        ->  Scale( 1., "width" );
    h1D_Ntru        ->  Scale( 1./kEventNormalisation );
    h2D_Ntru        ->  Scale( 1., "width" );
    h2D_Ntru        ->  Scale( 1./kEventNormalisation );
    hPhiCorrelation ->  Scale( 1., "width" );
    hPhiCorrelation ->  Scale( 1./kEventNormalisation );
    hProductionProb ->  SetBinContent( 1, kEventNormalisation-hProductionProb->Integral() );
    hProductionProb ->  Scale( 1./kEventNormalisation );
    hFullQuantities ->  Scale( 1./kEventNormalisation );
    h1D_Ntru_Mult   ->  Divide(fQC_Event_Enum_E05);
    h2D_Ntru_Mult   ->  Divide(fQC_Event_Enum_E05);
    hR1_Ntru_Mult   ->  Divide(h2D_Ntru_Mult,h1D_Ntru_Mult);
    hR2_Ntru_Mult   ->  Divide(hR1_Ntru_Mult,h1D_Ntru_Mult);
    for ( Int_t iBin = 1; iBin <= h1D_Ntru_Mult->GetNbinsX(); iBin++ )    {
        if ( h1D_Ntru_Mult->GetBinContent(iBin) == 0 ) continue;
        if ( h2D_Ntru_Mult->GetBinContent(iBin) == 0 ) continue;
        hP1_Ntru_Mult->SetBinContent(iBin, fSigmaPhiValue( h1D_Ntru_Mult->GetBinContent(iBin), h2D_Ntru_Mult->GetBinContent(iBin) ) );
        hP2_Ntru_Mult->SetBinContent(iBin, fGammaPhiValue( h1D_Ntru_Mult->GetBinContent(iBin), h2D_Ntru_Mult->GetBinContent(iBin) ) );
    }
    auto    k1DYield    =   hFullQuantities->GetBinContent(1);
    auto    k2DYield    =   hFullQuantities->GetBinContent(2);
    auto    k1DError    =   hFullQuantities->GetBinError(1);
    auto    k2DError    =   hFullQuantities->GetBinError(2);
    auto    k1DErrRl    =   k1DError / k1DYield;
    auto    k2DErrRl    =   k2DError / k2DYield;
    hFullQuantities ->  SetBinContent   ( 3, k2DYield/k1DYield );
    hFullQuantities ->  SetBinError     ( 3, (k2DYield/k1DYield)*SquareSum( { k2DErrRl, k1DErrRl } ) );
    hFullQuantities ->  SetBinContent   ( 4, (k2DYield/(k1DYield*k1DYield)) );
    hFullQuantities ->  SetBinError     ( 4, (k2DYield/(k1DYield*k1DYield))*SquareSum( { k2DErrRl, k1DErrRl, k1DErrRl } ) );
    hFullQuantities ->  SetBinContent   ( 5, fSigmaPhiValue(k1DYield,k2DYield) );
    hFullQuantities ->  SetBinError     ( 5, fSigmaPhiError(k1DYield,k2DYield,k1DError,k2DError) );
    hFullQuantities ->  SetBinContent   ( 6, fGammaPhiValue(k1DYield,k2DYield) );
    hFullQuantities ->  SetBinError     ( 6, fGammaPhiError(k1DYield,k2DYield,k1DError,k2DError) );
    //
    hDPhiDyCh       ->  Scale( 1., "width" );
    hDPhiDyCh       ->  Scale( 1./kEventNormalisation );
    hDPhiDyCh_Bkg   ->  Scale( 1., "width" );
    hDPhiDyCh_Bkg   ->  Scale( 1./kEventNormalisation );
    auto hDPhiDyCh_Rat = (TH2F*)(hDPhiDyCh->Clone("hDPhiDyCh_Rat"));
    hDPhiDyCh_Rat       -> Divide( hDPhiDyCh, hDPhiDyCh_Bkg );
    //
    double kSigInt = hDPhiDyCh       ->  Integral( 50, 100, 90, 100 );
    double kBkgInt = hDPhiDyCh_Bkg   ->  Integral( 50, 100, 90, 100 );
    //
    auto hDPhiDyCh_Sub = (TH2F*)(hDPhiDyCh->Clone("hDPhiDyCh_Sub"));
    hDPhiDyCh_Sub       -> Add( hDPhiDyCh_Bkg, -kSigInt/kBkgInt );
    //
    auto hDPhiDyCh_Eta_Rat  = (TH2F*)(hDPhiDyCh_Eta->Clone("hDPhiDyCh_Eta_Rat"));
    hDPhiDyCh_Eta_Rat       -> Divide( hDPhiDyCh_Eta, hDPhiDyCh_Eta_Bkg );
    auto hDPhiDyCh_Phi_Rat  = (TH2F*)(hDPhiDyCh_Phi->Clone("hDPhiDyCh_Phi_Rat"));
    hDPhiDyCh_Phi_Rat       -> Divide( hDPhiDyCh_Phi, hDPhiDyCh_Phi_Bkg );
    //
    auto hDPhiDyCh_Rap_Eta_Rat  = (TH2F*)(hDPhiDyCh_Rap_Eta->Clone("hDPhiDyCh_Rap_Eta_Rat"));
    hDPhiDyCh_Rap_Eta_Rat       -> Divide( hDPhiDyCh_Rap_Eta, hDPhiDyCh_Rap_Eta_Bkg );
    auto hDPhiDyCh_Rap_Phi_Rat  = (TH2F*)(hDPhiDyCh_Rap_Phi->Clone("hDPhiDyCh_Rap_Phi_Rat"));
    hDPhiDyCh_Rap_Phi_Rat       -> Divide( hDPhiDyCh_Rap_Phi, hDPhiDyCh_Rap_Phi_Bkg );
    //
    auto   hPhiCorrlRapRat = (TH2F*)(hPhiCorrlRapSig->Clone("hPhiCorrlRapRat"));
    hPhiCorrlRapRat       -> Divide( hPhiCorrlRapSig, hPhiCorrlRapBkg );
    //
    auto hPhiCorrelation_FRat  = (TH2F*)(hPhiCorrelation_Full->Clone("hPhiCorrelation_FRat"));
    hPhiCorrelation_FRat       -> Divide( hPhiCorrelation_Full, hPhiCorrelation_FBkg );
    //
    hPhiCorrelation_IM_Rat  = (TH1D*)(hPhiCorrelation_IM_Sig->Clone("hPhiCorrelation_IM_Rat"));
    hPhiCorrelation_IM_Rat       -> Divide( hPhiCorrelation_IM_Sig, hPhiCorrelation_IM_Bkg );
    //
    /*
    TF1*    kCustom = new TF1( "kCustom", "[3]*([0]*exp(-0.5*((x-[1])/[2])**2)+(1-[0]))" );
    TF1*    kCusto2 = new TF1( "kCusto2", "[0]*exp(-0.5*((x-[1])/[2])**2)+[3]*exp(-0.5*((x-[4])/[5])**2)+[6]" );
    kCustom->SetParLimits(0,0.,1.);
    kCustom->SetParameter(0,0.1);
    kCustom->FixParameter(1,0.);
    kCusto2->FixParameter(1,0.);
    kCusto2->FixParameter(4,0.);
    //hPhiCorrlRapRat->Fit(kCustom);
    //hPhiCorrelation->Fit(kCustom);
    //hPhiCorrelation_FRat->Fit(kCusto2);
     */
    //
    //  --- --- YIELD ANALYSIS
    if ( kDoYield ) {
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kProduction_MC_Dir,(TString("Yield")+kFolder).Data())));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kProduction_MC_Plot,(TString("Yield")+kFolder).Data()))+TString("/1D/")).Data()));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kProduction_MC_Plot,(TString("Yield")+kFolder).Data()))+TString("/2D/")).Data()));
        TFile *outFile_Yield    =   new TFile   (Form(kProduction_MC_Ofl,(TString("Yield")+kFolder).Data(),kProductionTag.Data()),"recreate");
        //
        // --- Saving to File
        fHEventCount            ->  Write();
        h1D_Ntru_Fin            ->  Write();
        h2D_Ntru_Fin_Nrm        ->  Write();
        h2D_Ntru_Fin            ->  Write();
        h1D_Ntru                ->  Write();
        h2D_Ntru                ->  Write();
        h1D_Ntru_Mult           ->  Write();
        h2D_Ntru_Mult           ->  Write();
        hR1_Ntru_Mult           ->  Write();
        hR2_Ntru_Mult           ->  Write();
        hP1_Ntru_Mult           ->  Write();
        hP2_Ntru_Mult           ->  Write();
        hMPT_NTru               ->  Write();
        hMPT_NTru_Fin           ->  Write();
        hFullQuantities         ->  Write();
        hPhiCorrelation         ->  Write();
        hProductionProb         ->  Write();
        hMultCorr               ->  Write();
        hMultUCor               ->  Write();
        hDPhiDyCh               ->  Write();
        hDPhiDyCh_Bkg           ->  Write();
        hDPhiDyCh_Rat           ->  Write();
        hDPhiDyCh_Sub           ->  Write();
        hDPhiDyCh_Rap_Eta_Bkg   ->  Write();
        hDPhiDyCh_Rap_Phi_Bkg   ->  Write();
        hDPhiDyCh_Rap_Eta       ->  Write();
        hDPhiDyCh_Rap_Phi       ->  Write();
        hDPhiDyCh_Rap_Eta_Rat   ->  Write();
        hDPhiDyCh_Rap_Phi_Rat   ->  Write();
        hDPhiDyCh_Eta           ->  Write();
        hDPhiDyCh_Phi           ->  Write();
        hDPhiDyCh_Eta_Bkg       ->  Write();
        hDPhiDyCh_Phi_Bkg       ->  Write();
        hDPhiDyCh_Eta_Rat       ->  Write();
        hDPhiDyCh_Phi_Rat       ->  Write();
        hPhiCorrlRapSig         ->  Write();
        hPhiCorrlRapBkg         ->  Write();
        hPhiCorrlRapRat         ->  Write();
        hPhiCorrelation_Full    ->  Write();
        hPhiCorrelation_FBkg    ->  Write();
        hPhiCorrelation_FRat    ->  Write();
        hPhiCorrelation_IM_Sig  ->  Write();
        hPhiCorrelation_IM_Bkg  ->  Write();
        hPhiCorrelation_IM_Rat  ->  Write();
        //
        outFile_Yield           ->  Close();
    }
}
