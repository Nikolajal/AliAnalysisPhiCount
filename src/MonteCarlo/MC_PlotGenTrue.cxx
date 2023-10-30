#include "../../inc/AliAnalysisPhiPair.h"

void
MC_PlotGenTrue_
 ( TString fFileName = "/Users/nrubini/Analysis/ALICE/PWG-LF/PAG-RSN/_1020_Phi_Pair/result/Yield_p_p__7TeV/MC_Production/7000/DataSet/Pythia6_MB_7TeV_Perugia_2011.root", TString kProductionTag = "Pythia6X0", TString fOption = "yield", Int_t nEventsCut = 30000000, TString kFolder = "_p_p__7TeV" );


void
MC_PlotGenTrue_
( TString fFileName = "/Users/nrubini/Analysis/ALICE/PWG-LF/PAG-RSN/_1020_Phi_Pair/result/Yield_p_p__7TeV/MC_Production/7000/DataSet/Pythia6_MB_7TeV_Perugia_2011.root", TString kProductionTag = "Pythia6X0", TString fOption = "yield", Int_t nEventsCut = -1, TString kFolder = "_p_p__7TeV" ) {
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
    // --- Retrieving TTree
    TTree*      TPhiEfficiency  =   (TTree*)insFile_Data->Get( "PhiCandidate_" );
    //
    // --- Retrieving Utility Histograms
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
    TH1D*       hProductionProb;
    TH1D*       hMultCorr;
    TH1D*       hMultUCor;
    TH2D*       hDPhiDyCh;
    TH2D*       hDPhiDyCh_Bkg;
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
    h1D_Ntru_Mult   =   new TH1D ( hName, hTitle, 500,   0.,  500.);
    SetAxis( h1D_Ntru_Mult, "PT 1D" );
    //
    hName           =   Form( "h2D_Ntru_Mult" );
    hTitle          =   Form( "True Kaons from phi in mult E10" );
    h2D_Ntru_Mult   =   new TH1D ( hName, hTitle, 500,   0.,  500.);
    SetAxis( h2D_Ntru_Mult, "PT 2D" );
    //
    hName           =   Form( "hR1_Ntru_Mult" );
    hTitle          =   Form( "True Kaons from phi in mult E10" );
    hR1_Ntru_Mult   =   new TH1D ( hName, hTitle, 500,   0.,  500.);
    SetAxis( hR1_Ntru_Mult, "PT 1D" );
    //
    hName           =   Form( "hR2_Ntru_Mult" );
    hTitle          =   Form( "True Kaons from phi in mult E10" );
    hR2_Ntru_Mult   =   new TH1D ( hName, hTitle, 500,   0.,  500.);
    SetAxis( hR2_Ntru_Mult, "PT 1D" );
    //
    hName           =   Form( "hP1_Ntru_Mult" );
    hTitle          =   Form( "True Kaons from phi in mult E10" );
    hP1_Ntru_Mult   =   new TH1D ( hName, hTitle, 500,   0.,  500.);
    SetAxis( hP1_Ntru_Mult, "PT 1D" );
    //
    hName           =   Form( "hP2_Ntru_Mult" );
    hTitle          =   Form( "True Kaons from phi in mult E10" );
    hP2_Ntru_Mult   =   new TH1D ( hName, hTitle, 500,   0.,  500.);
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
    hPhiCorrelation =   new TH1D ( hName, hTitle, 100, fMinCrPh, fMaxCrPh );
    SetAxis( hPhiCorrelation, "PT 1D" );
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
    hDPhiDyCh    =   new TH2D ( hName, hTitle, 100, -105, 255, 100, -1.0, 1.0 );
    SetAxis( hDPhiDyCh, "PT 2D" );
    //
    hName           =   Form( "hDPhiDyCh_Bkg" );
    hTitle          =   Form( "True Kaons from phi" );
    hDPhiDyCh_Bkg    =   new TH2D ( hName, hTitle, 100, -105, 255, 100, -1.0, 1.0 );
    SetAxis( hDPhiDyCh_Bkg, "PT 2D" );
    //
    //  --- Applying limitation sample cuts
    auto nEvents = (!TPhiEfficiency) ? 0 : ( nEventsCut == -1.? TPhiEfficiency->GetEntries() : nEventsCut);
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
                            hDPhiDyCh_Bkg   -> Fill( fGetDltCrPh( fCurrent_MCTrue.Phi[iPhi] - fPrevious_MCTrue.Phi[jPhi] ), fCurrent_MCTrue.Rap[iPhi] - fPrevious_MCTrue.Rap[jPhi] );
                        }
                    }
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
                        h2D_Ntru_Mult   ->  Fill( fCurrent_MCTrue.Eta_05 );
                        hFullQuantities ->  Fill( 2 );
                        hDPhiDyCh       ->  Fill( fGetDltCrPh( fCurrent_MCTrue.Phi[iPhi] - fCurrent_MCTrue.Phi[jPhi] ), fCurrent_MCTrue.Rap[iPhi] - fCurrent_MCTrue.Rap[jPhi] );
                        //if ( fCurrent_MCTrue.nPart ) hPhiCorrelatio2 ->  Fill( fGetDltCrPh( fCurrent_MCTrue.Phi[iPhi] - fCurrent_MCTrue.Phi[jPhi] ) );
                        //hMultCorr       ->  Fill(  );
                   }
                }
            }
        }
        fPrevious_MCTrue = fCurrent_MCTrue;
        //
    }
    if ( nEvents > 0 )  fStopTimer("Comparison Plots Production");
    //
    // --- --- --- --- --- --- --- OUTPUT --- --- --- --- --- --- --- --- --- --- ---
    //
    auto    kTotalNormalisation = fHEventCount->GetBinContent(1);
    auto    kEventNormalisation = nEventsCut < 0 ? kTotalNormalisation : kTotalNormalisation*(nEventsCut)/(TPhiEfficiency->GetEntries());
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
    uSetHisto( hFullQuantities,       "FNL STAT 1D" );
    //
    hDPhiDyCh       ->  Scale( 1., "width" );
    hDPhiDyCh       ->  Scale( 1./kEventNormalisation );
    hDPhiDyCh_Bkg   ->  Scale( 1., "width" );
    hDPhiDyCh_Bkg   ->  Scale( 1./kEventNormalisation );
    //
    double kSigInt = hDPhiDyCh       ->  Integral( 50, 100, 90, 100 );
    double kBkgInt = hDPhiDyCh_Bkg   ->  Integral( 50, 100, 90, 100 );
    //
    auto hDPhiDyCh_2 = (TH2F*)(hDPhiDyCh->Clone("hDPhiDyCh_2"));
    hDPhiDyCh_2       -> Add( hDPhiDyCh_Bkg, -kSigInt/kBkgInt );
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
        hDPhiDyCh_2             ->  Write();
        hDPhiDyCh_Bkg           ->  Write();
        //
        outFile_Yield           ->  Close();
    }
}
