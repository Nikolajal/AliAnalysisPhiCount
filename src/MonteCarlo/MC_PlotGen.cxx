#include "../../inc/AliAnalysisPhiPair.h"

void
MC_PlotGen
( TString fFileName = "/Volumes/NRUBINI_DATASTASH/Dataset/_Sim/p_p__7TeV/LHC14j4_STD.root", TString kProductionTag = "Pythia6", TString fOption = "yield", Int_t nEventsCut = -1, TString kFolder = "_p_p__7TeV", Int_t iMode = 0 ) {
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
    
    auto kDebugLimit = 1.0;
    
    // --- Retrieving TTree
    //TTree*      TPhiEfficiency  =   (TTree*)insFile_Data->Get( Form("Phi_SX_E7000_M%i",iMode) );
    //TTree*      TPhiEfficiency  =   (TTree*)insFile_Data->Get( fPhiCandidateEff_Tree );
    TTree*      TPhiEfficiency  =   (TTree*)insFile_Data->Get( "PhiCandidate_" );
    //
    // --- Retrieving Utility Histograms
    //TH1D*       fHEventCount    =   (TH1D*) insFile_Data->Get( "kEventCount" );
    TList*      fQCOutputList       =   (TList*)insFile_Data->Get( "fQCOutputList_PhiCount_STD" );
    TH1D*       fHEventCount        =   (TH1D*) fQCOutputList->FindObject( "fQC_Event_Enum_FLL" );
    TH1D*       fQC_Event_Enum_E05  =   (TH1D*) fQCOutputList->FindObject( "fQC_Event_Enum_E10" );
    TH1D*       fQC_Event_Enum_V0M  =   (TH1D*) fQCOutputList->FindObject( "fQC_Event_Enum_V0M" );
    //
    // --- Setting the input datastructure
    Struct_MCParticle    evPhiEfficiency;
    fSetPhiCandidateMC(TPhiEfficiency,evPhiEfficiency);
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
    hFullQuantities =   new TH1D ( hName, hTitle, 6, 0.5, 6.5 );
    //
    hName           =   Form( "hPhiCorrelation" );
    hTitle          =   Form( "hPhiCorrelation" );
    hPhiCorrelation =   new TH1D ( hName, hTitle, 100, fMinCrPh, fMaxCrPh );
    SetAxis( hPhiCorrelation, "PT 1D" );
    //
    /*
    hName           =   Form( "hPhiCorrelatio2" );
    hTitle          =   Form( "hPhiCorrelatio2" );
    hPhiCorrelatio2 =   new TH1D ( hName, hTitle, 100, fMinCrPh, fMaxCrPh );
    SetAxis( hPhiCorrelatio2, "PT 1D" );
     */
    //
    //  --- Applying limitation sample cuts
    auto nEvents = (!TPhiEfficiency) ? 0 : ( nEventsCut == -1.? TPhiEfficiency->GetEntries() : nEventsCut);
    nEvents *= kDebugLimit;
    if ( nEvents > 0 )  fStartTimer("Comparison Plots Production");
    //
    for ( Int_t iEvent = 0; iEvent < nEvents; iEvent++ )    {
        TPhiEfficiency->GetEntry(iEvent);
        fPrintLoopTimer("Comparison Plots Production",iEvent,nEvents,kPrintIntervalPP);
        //
        //  --- Loop over candidates
        Struct_MCParticle fCurrent_MCTrue;
        fCurrent_MCTrue.nPart    =   0;
        for ( Int_t iPhi = 0; iPhi < evPhiEfficiency.nPart; iPhi++ )  {
            TLorentzVector  kTLVUtility;   // TODO: Update to new lorentzvector class
            kTLVUtility.SetXYZM(evPhiEfficiency.Px[iPhi],evPhiEfficiency.Py[iPhi],evPhiEfficiency.Pz[iPhi],kPhiMesonMass_);
            //  --- Load Useful Variables for histograms filling
            fCurrent_MCTrue.pT          [ fCurrent_MCTrue.nPart ]    =   kTLVUtility.Pt();
            fCurrent_MCTrue.Rap         [ fCurrent_MCTrue.nPart ]    =   kTLVUtility.Rapidity();
            fCurrent_MCTrue.Phi         [ fCurrent_MCTrue.nPart ]    =   kTLVUtility.Phi() * 360 / ( TMath::Pi() * 2 );
            fCurrent_MCTrue.iPT1D       [ fCurrent_MCTrue.nPart ]    =   fGetBinPT1D( kTLVUtility.Pt() );
            fCurrent_MCTrue.iPT2D       [ fCurrent_MCTrue.nPart ]    =   fGetBinPT2D( kTLVUtility.Pt() );
            if ( is_pp_anl )    fCurrent_MCTrue.kHasRap     [ fCurrent_MCTrue.nPart ]    =   fabs( kTLVUtility.Rapidity() ) < 0.5;
            if ( is_pb_anl )    fCurrent_MCTrue.kHasRap     [ fCurrent_MCTrue.nPart ]    =   ( kTLVUtility.Rapidity() ) < 0.035 && ( kTLVUtility.Rapidity() ) > -0.465;
            fCurrent_MCTrue.nPart++;
        }
        //
        //  --- Discarding if no valid candidates are found
        if ( fCurrent_MCTrue.nPart < 1 ) continue;
        //
        //  --- Event variables
        fCurrent_MCTrue.IsMB            =   true;//( evPhiEfficiency.TrueEventMask == 0 );
        fCurrent_MCTrue.Eta_10          =   evPhiEfficiency.Eta_10;
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
                    h1D_Ntru_Mult   ->  Fill( fCurrent_MCTrue.Eta_10 );
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
                        h2D_Ntru        ->  Fill( fCurrent_MCTrue.pT[iPhi], fCurrent_MCTrue.pT[jPhi] );
                        h2D_Ntru_Fin    ->  Fill( fCurrent_MCTrue.pT[iPhi], fCurrent_MCTrue.pT[jPhi] );
                        h2D_Ntru_Fin_Nrm->  Fill( fCurrent_MCTrue.pT[iPhi], fCurrent_MCTrue.pT[jPhi] );
                        hMPT_NTru       ->  Fill( fCurrent_MCTrue.pT[iPhi], fCurrent_MCTrue.pT[jPhi] );
                        hMPT_NTru_Fin   ->  Fill( fCurrent_MCTrue.pT[iPhi], fCurrent_MCTrue.pT[jPhi] );
                        h2D_Ntru_Mult   ->  Fill( fCurrent_MCTrue.Eta_10 );
                        hFullQuantities ->  Fill( 2 );
                        //if ( fCurrent_MCTrue.nPart ) hPhiCorrelatio2 ->  Fill( fGetDltCrPh( fCurrent_MCTrue.Phi[iPhi] - fCurrent_MCTrue.Phi[jPhi] ) );
                   }
                }
            }
        }
    }
    if ( nEvents > 0 )  fStopTimer("Comparison Plots Production");
    //
    // --- --- --- --- --- --- --- OUTPUT --- --- --- --- --- --- --- --- --- --- ---
    //
    //auto    kEventNormalisation =   fHEventCount->GetEntries();
    auto    kEventNormalisation =   kDebugLimit*fHEventCount->GetBinContent(1);
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
        //
        outFile_Yield           ->  Close();
    }
}
