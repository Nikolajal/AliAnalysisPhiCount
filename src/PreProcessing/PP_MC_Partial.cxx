#include "../../inc/AliAnalysisPhiPair.h"
// !TODO: AXIS and TITLES conforming to ALICE guidelines
// !TODO: Multiplicity comparison of efficiencies

void PP_MC_Partial ( TString fFileName, TString fOption, Int_t nEventsCut = -1, TString kFolder = "" )  {
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
    TTree*      TPhiEfficiency  =   (TTree*)insFile_Data->Get( fPhiCandidateEff_Tree +TString("") );
    TTree*      TPhiCandidate   =   (TTree*)insFile_Data->Get( fPhiCandidate_Tree +TString("") );
    //
    // --- Retrieving Utility Histograms
    TList*      fQCOutputList   =   (TList*)insFile_Data->Get( "fQCOutputList" );
    TH1D*       fHEventCount    =   (TH1D*) fQCOutputList->FindObject( "fQC_Event_Enum_FLL" );
    TH1D*       fHEvCountMlt    =   (TH1D*) fQCOutputList->FindObject( "fQC_Event_Enum_V0M" );
    //
    // --- Setting the input datastructure
    Struct_PhiCandidate     evPhiCandidate;
    Struct_PhiEfficiency    evPhiEfficiency;
    Struct_KaonCandidate    kfix;           //TODO: Clean the fix
    Struct_KaonEfficiency   kfi2;           //TODO: Clean the fix
    if ( !fSetCandidates(TPhiCandidate,evPhiCandidate,nullptr,kfix) )   cout << "WW" << endl;    // TODO: Fix the TrueInvMass in case of Data load
    if ( !fSetCandidates(TPhiEfficiency,evPhiEfficiency,nullptr,kfi2) ) cout << "rW" << endl;
    //
    // --- Setting the output datastructure
    fSetAllBins();
    auto    kUniformBinning =   uUniformBinning( 10*MeV, fMinPT1D, fMaxPT1D );
    //
    TH1F*       h1D_Nrec;
    TH1F*       h1D_Ngen;
    TH1F*       h1D_Ntru;
    TH1F*       h1D_Nrec_2Db;
    TH1F*       h1D_Ngen_2Db;
    TH1F*       h1D_Ntru_2Db;
    TH1F*       h1D_Nrec_Fin;
    TH1F*       h1D_Ngen_Fin;
    TH1F*       h1D_Ntru_Fin;
    TH1F*       h1D_Nrec_CrPh;
    TH1F*       h1D_Ngen_CrPh;
    TH1F*       h1D_Ntru_CrPh;
    TH2F*       h2D_Nrec;
    TH2F*       h2D_Ngen;
    TH2F*       h2D_Ntru;
    std::vector<TH1F*>  h1D_TruInvMass;
    std::vector<TH1F*>  h1D_RecInvMass;
    std::vector<TH1F*>  h1D_InvMassRes;
    std::vector<TH1F*>  h1D_TruInvMass_2Db;
    std::vector<TH1F*>  h1D_RecInvMass_2Db;
    std::vector<TH1F*>  h1D_InvMassRes_2Db;
    //
    hName           =   Form( "h1D_Nrec" );
    hTitle          =   Form( "Recordable Kaons from phi" );
    h1D_Nrec        =   new TH1F ( hName, hTitle, nBinPT1D, fArrPT1D );
    SetAxis( h1D_Nrec, "PT 1D" );
    //
    hName           =   Form( "h1D_Ngen" );
    hTitle          =   Form( "Generated Kaons from phi" );
    h1D_Ngen        =   new TH1F ( hName, hTitle, nBinPT1D, fArrPT1D );
    SetAxis( h1D_Ngen, "PT 1D" );
    //
    hName           =   Form( "h1D_Ntru" );
    hTitle          =   Form( "True Kaons from phi" );
    h1D_Ntru        =   new TH1F ( hName, hTitle, nBinPT1D, fArrPT1D );
    SetAxis( h1D_Ntru, "PT 1D" );
    //
    hName           =   Form( "h1D_Nrec_2Db" );
    hTitle          =   Form( "Recordable Kaons from phi" );
    h1D_Nrec_2Db    =   new TH1F ( hName, hTitle, nBinPT2D, fArrPT2D );
    SetAxis( h1D_Nrec_2Db, "PT 1D" );
    //
    hName           =   Form( "h1D_Ngen_2Db" );
    hTitle          =   Form( "Generated Kaons from phi" );
    h1D_Ngen_2Db    =   new TH1F ( hName, hTitle, nBinPT2D, fArrPT2D );
    SetAxis( h1D_Ngen_2Db, "PT 1D" );
    //
    hName           =   Form( "h1D_Ntru_2Db" );
    hTitle          =   Form( "True Kaons from phi" );
    h1D_Ntru_2Db    =   new TH1F ( hName, hTitle, nBinPT2D, fArrPT2D );
    SetAxis( h1D_Ntru_2Db, "PT 1D" );
    //
    hName           =   Form( "h1D_Nrec_Fin" );
    hTitle          =   Form( "Recordable Kaons from phi" );
    h1D_Nrec_Fin    =   new TH1F ( hName, hTitle, kUniformBinning.first , kUniformBinning.second );
    SetAxis( h1D_Nrec_Fin, "PT 1D" );
    //
    hName           =   Form( "h1D_Ngen_Fin" );
    hTitle          =   Form( "Generated Kaons from phi" );
    h1D_Ngen_Fin    =   new TH1F ( hName, hTitle, kUniformBinning.first , kUniformBinning.second );
    SetAxis( h1D_Ngen_Fin, "PT 1D" );
    //
    hName           =   Form( "h1D_Ntru_Fin" );
    hTitle          =   Form( "True Kaons from phi" );
    h1D_Ntru_Fin    =   new TH1F ( hName, hTitle, kUniformBinning.first , kUniformBinning.second );
    SetAxis( h1D_Ntru_Fin, "PT 1D" );
    //
    hName           =   Form( "h1D_Nrec_CrPh" );
    hTitle          =   Form( "True Kaons from phi" );
    h1D_Nrec_CrPh   =   new TH1F ( hName, hTitle, nBinCrPh, fArrCrPh );
    SetAxis( h1D_Nrec_CrPh, "PT 1D" );
    //
    hName           =   Form( "h1D_Ngen_CrPh" );
    hTitle          =   Form( "True Kaons from phi" );
    h1D_Ngen_CrPh   =   new TH1F ( hName, hTitle, nBinCrPh, fArrCrPh );
    SetAxis( h1D_Ngen_CrPh, "PT 1D" );
    //
    hName           =   Form( "h1D_Ntru_CrPh" );
    hTitle          =   Form( "True Kaons from phi" );
    h1D_Ntru_CrPh   =   new TH1F ( hName, hTitle, nBinCrPh, fArrCrPh );
    SetAxis( h1D_Ntru_CrPh, "PT 1D" );
    //
    hName           =   Form( "h2D_Nrec" );
    hTitle          =   Form( "Recordable Kaons from phi" );
    h2D_Nrec        =   new TH2F ( hName, hTitle, nBinPT2D, fArrPT2D, nBinPT2D, fArrPT2D );
    SetAxis( h1D_Nrec, "PT 2D" );
    //
    hName           =   Form( "h2D_Ngen" );
    hTitle          =   Form( "Generated Kaons from phi" );
    h2D_Ngen        =   new TH2F ( hName, hTitle, nBinPT2D, fArrPT2D, nBinPT2D, fArrPT2D );
    SetAxis( h1D_Ngen, "PT 2D" );
    //
    hName           =   Form( "h2D_Ntru" );
    hTitle          =   Form( "True Kaons from phi" );
    h2D_Ntru        =   new TH2F ( hName, hTitle, nBinPT2D, fArrPT2D, nBinPT2D, fArrPT2D );
    SetAxis( h1D_Ntru, "PT 2D" );
    //
    for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ ) {
        TH1F*                       h1D_TruInvMass_PT_Util;
        hName                       =   Form( "h1D_TruInvMass_PT_%i", iPT1D );
        hTitle                      =   Form( "m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c", fArrPT1D[iPT1D], fArrPT1D[iPT1D+1] );
        h1D_TruInvMass_PT_Util      =   new TH1F ( hName, hTitle, nBinIMTR, fArrIMTR );
        SetAxis( h1D_TruInvMass_PT_Util, "IM 1D" );
        h1D_TruInvMass.push_back( h1D_TruInvMass_PT_Util );
        //
        TH1F*                       h1D_RecInvMass_PT_Util;
        hName                       =   Form( "h1D_RecInvMass_PT_%i", iPT1D );
        hTitle                      =   Form( "m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c", fArrPT1D[iPT1D], fArrPT1D[iPT1D+1] );
        h1D_RecInvMass_PT_Util      =   new TH1F ( hName, hTitle, nBinIMRC, fArrIMRC );
        SetAxis( h1D_RecInvMass_PT_Util, "IM 1D" );
        h1D_RecInvMass.push_back( h1D_RecInvMass_PT_Util );
        //
        TH1F*                       h1D_InvMassRes_PT_Util;
        hName                       =   Form( "h1D_InvMassRes_PT_%i", iPT1D );
        hTitle                      =   Form( "m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c", fArrPT1D[iPT1D], fArrPT1D[iPT1D+1] );
        h1D_InvMassRes_PT_Util      =   new TH1F ( hName, hTitle, nBinIMRs, fArrIMRs );
        SetAxis( h1D_InvMassRes_PT_Util, "IM 1D" );
        h1D_InvMassRes.push_back( h1D_InvMassRes_PT_Util );
    }
    //
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ ) {
        TH1F*                       h1D_TruInvMass_PT_Util;
        hName                       =   Form( "h1D_TruInvMass_2Db_PT_%i", iPT2D );
        hTitle                      =   Form( "m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c", fArrPT2D[iPT2D], fArrPT2D[iPT2D+1] );
        h1D_TruInvMass_PT_Util      =   new TH1F ( hName, hTitle, nBinIMTR, fArrIMTR );
        SetAxis( h1D_TruInvMass_PT_Util, "IM 1D" );
        h1D_TruInvMass_2Db.push_back( h1D_TruInvMass_PT_Util );
        //
        TH1F*                       h1D_RecInvMass_PT_Util;
        hName                       =   Form( "h1D_RecInvMass_2Db_PT_%i", iPT2D );
        hTitle                      =   Form( "m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c", fArrPT2D[iPT2D], fArrPT2D[iPT2D+1] );
        h1D_RecInvMass_PT_Util      =   new TH1F ( hName, hTitle, nBinIMRC, fArrIMRC );
        SetAxis( h1D_RecInvMass_PT_Util, "IM 1D" );
        h1D_RecInvMass_2Db.push_back( h1D_RecInvMass_PT_Util );
        //
        TH1F*                       h1D_InvMassRes_PT_Util;
        hName                       =   Form( "h1D_InvMassRes_2Db_PT_%i", iPT2D );
        hTitle                      =   Form( "m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c", fArrPT2D[iPT2D], fArrPT2D[iPT2D+1] );
        h1D_InvMassRes_PT_Util      =   new TH1F ( hName, hTitle, nBinIMRs, fArrIMRs );
        SetAxis( h1D_InvMassRes_PT_Util, "IM 1D" );
        h1D_InvMassRes_2Db.push_back( h1D_InvMassRes_PT_Util );
    }
    //
    // --- --- --- --- --- --- --- ANALYSIS --- --- --- --- --- --- --- --- --- --- -
    //
    //  --- Applying limitation sample cuts
    Int_t nEvents = (!TPhiCandidate) ? 0 : ( nEventsCut == -1.? TPhiCandidate->GetEntries() : nEventsCut);
    if ( nEvents > 0 )  fStartTimer("Resolution Analysis");
    //
    for ( Int_t iEvent = 0; iEvent < nEvents; iEvent++ )    {
        TPhiCandidate->GetEntry(iEvent);
        fPrintLoopTimer("Resolution Analysis",iEvent,nEvents,kPrintIntervalPP);
        //
        //  --- Evaluate the Event is of interest for the analysis
        //
        //  --- Discarding Pile-up events
        if ( kCheckPileUp   &&  kDoYield           &&  fCheckMask(evPhiCandidate.EventMask,1) ) continue;
        if ( kCheckPileUp   &&  kDoMultiplicity    &&  ( fCheckMask(evPhiCandidate.EventMask,1) || fCheckMask(evPhiCandidate.EventMask,2) )) continue;
        //
        //  --- Loop over candidates
        Struct_PhiCandidate fCurrent_Candidates;
        fCurrent_Candidates.nPhi    =   0;
        for ( Int_t iPhi = 0; iPhi < evPhiCandidate.nPhi; iPhi++ )  {
            TLorentzVector  kTLVUtility;   // TODO: Update to new lorentzvector class
            kTLVUtility.SetXYZM(evPhiCandidate.Px[iPhi],evPhiCandidate.Py[iPhi],evPhiCandidate.Pz[iPhi],evPhiCandidate.InvMass[iPhi]);
            if ( !fAcceptCandidate(evPhiCandidate.InvMass[iPhi],kTLVUtility.Pt()) ) continue;
            //  --- Load Useful Variables for histograms filling
            fCurrent_Candidates.InvMass     [ fCurrent_Candidates.nPhi ]    =   evPhiCandidate.InvMass[iPhi];
            fCurrent_Candidates.TrueInvMass [ fCurrent_Candidates.nPhi ]    =   evPhiCandidate.TrueInvMass[iPhi];
            fCurrent_Candidates.pT          [ fCurrent_Candidates.nPhi ]    =   kTLVUtility.Pt();
            fCurrent_Candidates.Rap         [ fCurrent_Candidates.nPhi ]    =   kTLVUtility.Rapidity();
            fCurrent_Candidates.iPT1D       [ fCurrent_Candidates.nPhi ]    =   fGetBinPT1D( kTLVUtility.Pt() );
            fCurrent_Candidates.iPT2D       [ fCurrent_Candidates.nPhi ]    =   fGetBinPT2D( kTLVUtility.Pt() );
            if ( is_pp_anl )    fCurrent_Candidates.kHasRap     [ fCurrent_Candidates.nPhi ]    =   fabs( kTLVUtility.Rapidity() ) < 0.5;
            if ( is_pb_anl )    fCurrent_Candidates.kHasRap     [ fCurrent_Candidates.nPhi ]    =   min( 0., kTLVUtility.Rapidity() +0.465 ) < 0.5;
            fCurrent_Candidates.nPhi++;
        }
        //
        //  --- Discarding if no valid candidates are found
        if ( fCurrent_Candidates.nPhi == 0 ) continue;
        //
        //  --- Ordering in pT the candidates
        if ( fCurrent_Candidates.nPhi >  1 )    uOrderPTCandidates(fCurrent_Candidates);
        //
        //  --- Event variables
        fCurrent_Candidates.Multiplicity    =   evPhiCandidate.Multiplicity;
        fCurrent_Candidates.iMult           =   fGetBinMult(evPhiCandidate.Multiplicity);
        fCurrent_Candidates.kHasMult        =   fCurrent_Candidates.iMult != -1;
        //
        for ( Int_t iPhi = 0; iPhi < fCurrent_Candidates.nPhi; iPhi++ )  {
            //
            //  --- Selecting valid candidates
            if ( !fAcceptCandidate( fCurrent_Candidates, iPhi ) ) continue;
            //
            if ( fCurrent_Candidates.kHasRap[iPhi] )  {       //  --- Mid-Rapidity Analyses
                if ( kDoYield ) {       //  --- YIELD ANALYSIS
                    if ( fCurrent_Candidates.TrueInvMass[iPhi] != 0 )   {
                        h1D_TruInvMass      .at( fCurrent_Candidates.iPT1D[iPhi] )  ->  Fill    ( fCurrent_Candidates.TrueInvMass[iPhi] );
                        h1D_RecInvMass      .at( fCurrent_Candidates.iPT1D[iPhi] )  ->  Fill    ( fCurrent_Candidates.InvMass[iPhi] );
                        h1D_InvMassRes      .at( fCurrent_Candidates.iPT1D[iPhi] )  ->  Fill    ( fCurrent_Candidates.TrueInvMass[iPhi] - fCurrent_Candidates.InvMass[iPhi] );
                        h1D_TruInvMass_2Db  .at( fCurrent_Candidates.iPT2D[iPhi] )  ->  Fill    ( fCurrent_Candidates.TrueInvMass[iPhi] );
                        h1D_RecInvMass_2Db  .at( fCurrent_Candidates.iPT2D[iPhi] )  ->  Fill    ( fCurrent_Candidates.InvMass[iPhi] );
                        h1D_InvMassRes_2Db  .at( fCurrent_Candidates.iPT2D[iPhi] )  ->  Fill    ( fCurrent_Candidates.TrueInvMass[iPhi] - fCurrent_Candidates.InvMass[iPhi] );
                    }
                } if ( kDoMultiplicity )    {   //  --- MULTIPLICITY ANALYSIS
                    
                }
            }
            //
            //  --- Speed-up protection
            if ( fCurrent_Candidates.nPhi < 2 ) continue;
            for ( Int_t jPhi = iPhi+1; jPhi < fCurrent_Candidates.nPhi; jPhi++ )  {
                //
                //  --- Selecting valid candidates
                if ( !fAcceptCandidate( fCurrent_Candidates, iPhi ) ) continue;
                //
                if ( fCurrent_Candidates.kHasRap[iPhi] && fCurrent_Candidates.kHasRap[jPhi] )  {    //  --- Mid-Rapidity Analyses
                    if ( kDoYield ) {   //  --- YIELD ANALYSIS
                        
                    } if ( kDoMultiplicity )    {   //  --- MULTIPLICITY ANALYSIS
                        
                    }
                }
            }
        }
    }
    if ( nEvents > 0 )  fStopTimer("Resolution Analysis");
    //
    //  --- Applying limitation sample cuts
         nEvents = (!TPhiEfficiency) ? 0 : ( nEventsCut == -1.? TPhiEfficiency->GetEntries() : nEventsCut);
    if ( nEvents > 0 )  fStartTimer("Efficiency Analysis");
    //
    for ( Int_t iEvent = 0; iEvent < nEvents; iEvent++ )    {
        TPhiEfficiency->GetEntry(iEvent);
        fPrintLoopTimer("Efficiency Analysis",iEvent,nEvents,kPrintIntervalPP);
        //
        //  --- Evaluate the Event is of interest for the analysis
        //
        //  --- Loop over candidates
        Struct_PhiEfficiency fCurrent_MCTrue;
        fCurrent_MCTrue.nPhi    =   0;
        for ( Int_t iPhi = 0; iPhi < evPhiEfficiency.nPhi; iPhi++ )  {
            TLorentzVector  kTLVUtility;   // TODO: Update to new lorentzvector class
            kTLVUtility.SetXYZM(evPhiEfficiency.Px[iPhi],evPhiEfficiency.Py[iPhi],evPhiEfficiency.Pz[iPhi],kPhiMesonMass_);
            //  --- Load Useful Variables for histograms filling
            fCurrent_MCTrue.pT          [ fCurrent_MCTrue.nPhi ]    =   kTLVUtility.Pt();
            fCurrent_MCTrue.Rap         [ fCurrent_MCTrue.nPhi ]    =   kTLVUtility.Rapidity();
            fCurrent_MCTrue.Phi         [ fCurrent_MCTrue.nPhi ]    =   kTLVUtility.Phi() * 360 / ( TMath::Pi() * 2 );
            fCurrent_MCTrue.iPT1D       [ fCurrent_MCTrue.nPhi ]    =   fGetBinPT1D( kTLVUtility.Pt() );
            fCurrent_MCTrue.iPT2D       [ fCurrent_MCTrue.nPhi ]    =   fGetBinPT2D( kTLVUtility.Pt() );
            fCurrent_MCTrue.kIsGen      [ fCurrent_MCTrue.nPhi ]    =   evPhiEfficiency.Selection[iPhi] > 0;
            fCurrent_MCTrue.kIsRec      [ fCurrent_MCTrue.nPhi ]    =   evPhiEfficiency.Selection[iPhi] > 1;
            if ( is_pp_anl )    fCurrent_MCTrue.kHasRap     [ fCurrent_MCTrue.nPhi ]    =   fabs( kTLVUtility.Rapidity() ) < 0.5;
            if ( is_pb_anl )    fCurrent_MCTrue.kHasRap     [ fCurrent_MCTrue.nPhi ]    =   min( 0., kTLVUtility.Rapidity() +0.465 ) < 0.5;
            fCurrent_MCTrue.nPhi++;
        }
        //
        //  --- Discarding if no valid candidates are found
        if ( fCurrent_MCTrue.nPhi == 0 ) continue;
        //
        //  --- Event variables
        fCurrent_MCTrue.IsMB            =   ( evPhiEfficiency.TrueEventMask == 0 );
        fCurrent_MCTrue.Multiplicity    =   evPhiEfficiency.Multiplicity;
        fCurrent_MCTrue.iMult           =   fGetBinMult(evPhiEfficiency.Multiplicity);
        fCurrent_MCTrue.kHasMult        =   fCurrent_MCTrue.iMult != -1;
        //
        for ( Int_t iPhi = 0; iPhi < fCurrent_MCTrue.nPhi; iPhi++ )  {
            //
            if ( fCurrent_MCTrue.kHasRap[iPhi] && fCurrent_MCTrue.IsMB )  {       //  --- Mid-Rapidity Analyses
                if ( kDoYield ) {       //  --- YIELD ANALYSIS
                    h1D_Ntru        ->  Fill( fCurrent_MCTrue.pT[iPhi] );
                    h1D_Ntru_2Db    ->  Fill( fCurrent_MCTrue.pT[iPhi] );
                    h1D_Ntru_Fin    ->  Fill( fCurrent_MCTrue.pT[iPhi] );
                    if ( fCurrent_MCTrue.kIsGen[iPhi] ) {
                        h1D_Ngen        ->  Fill( fCurrent_MCTrue.pT[iPhi] );
                        h1D_Ngen_2Db    ->  Fill( fCurrent_MCTrue.pT[iPhi] );
                        h1D_Ngen_Fin    ->  Fill( fCurrent_MCTrue.pT[iPhi] );
                    } if ( fCurrent_MCTrue.kIsRec[iPhi] ) {
                        h1D_Nrec        ->  Fill( fCurrent_MCTrue.pT[iPhi] );
                        h1D_Nrec_2Db    ->  Fill( fCurrent_MCTrue.pT[iPhi] );
                        h1D_Nrec_Fin    ->  Fill( fCurrent_MCTrue.pT[iPhi] );
                    }
                } if ( kDoMultiplicity )    {   //  --- MULTIPLICITY ANALYSIS
                    
                }
            }
            //
            //  --- Speed-up protection
            if ( fCurrent_MCTrue.nPhi < 2 ) continue;
            for ( Int_t jPhi = iPhi+1; jPhi < fCurrent_MCTrue.nPhi; jPhi++ )  {
                //
                if ( fCurrent_MCTrue.kHasRap[iPhi] && fCurrent_MCTrue.kHasRap[jPhi] && fCurrent_MCTrue.IsMB )  {    //  --- Mid-Rapidity Analyses
                    //  --- TEST
                    auto    kDeltaPhi   =   ( fCurrent_MCTrue.Phi[iPhi] - fCurrent_MCTrue.Phi[jPhi]  );
                    kDeltaPhi = kDeltaPhi < -180 ? kDeltaPhi + 360 : kDeltaPhi > 180 ? kDeltaPhi -360 : kDeltaPhi;
                    //
                    if ( kDoYield ) {   //  --- YIELD ANALYSIS
                        h2D_Ntru        ->  Fill( fCurrent_MCTrue.pT[iPhi], fCurrent_MCTrue.pT[jPhi] );
                        if ( fCurrent_MCTrue.kIsGen[iPhi] ) {
                            h2D_Ngen        ->  Fill( fCurrent_MCTrue.pT[iPhi], fCurrent_MCTrue.pT[jPhi] );
                        } if ( fCurrent_MCTrue.kIsRec[iPhi] ) {
                            h2D_Nrec        ->  Fill( fCurrent_MCTrue.pT[iPhi], fCurrent_MCTrue.pT[jPhi] );
                        }
                    } if ( kDoMultiplicity )    {   //  --- MULTIPLICITY ANALYSIS
                        
                    } if ( kDoCorrelation )    {    //  --- CORRELATION ANALYSIS
                        h1D_Ntru_CrPh   ->  Fill ( kDeltaPhi );
                        if ( fCurrent_MCTrue.kIsGen[iPhi] ) {
                            h1D_Ngen_CrPh   ->  Fill ( kDeltaPhi );
                        } if ( fCurrent_MCTrue.kIsRec[iPhi] ) {
                            h1D_Nrec_CrPh   ->  Fill ( kDeltaPhi );
                        }
                   }
                }
            }
        }
    }
    if ( nEvents > 0 )  fStopTimer("Efficiency Analysis");
    //
    // --- --- --- --- --- --- --- OUTPUT --- --- --- --- --- --- --- --- --- --- ---
    //
    auto    kEventNormalisation =   fHEventCount->GetBinContent(1);
    //
    //  --- --- YIELD ANALYSIS
    if ( kDoYield ) {
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_PreProc_Dir,(TString("Yield")+kFolder).Data())));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Yield")+kFolder).Data()))+TString("/Plots/1D/")).Data()));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Yield")+kFolder).Data()))+TString("/Plots/2D/")).Data()));
        TFile *outFile_Yield    =   new TFile   (Form(kAnalysis_MCTruthHist,(TString("Yield")+kFolder).Data()),"recreate");
        //
        // --- Saving to File
        fHEventCount    ->  Write();
        fHEvCountMlt    ->  Write();
        h1D_Nrec        ->  Write();
        h1D_Ngen        ->  Write();
        h1D_Ntru        ->  Write();
        h1D_Nrec_2Db    ->  Write();
        h1D_Ngen_2Db    ->  Write();
        h1D_Ntru_2Db    ->  Write();
        h1D_Nrec_Fin    ->  Write();
        h1D_Ngen_Fin    ->  Write();
        h1D_Ntru_Fin    ->  Write();
        h2D_Nrec        ->  Write();
        h2D_Ngen        ->  Write();
        h2D_Ntru        ->  Write();
        for ( auto kSave : h1D_TruInvMass )     kSave   ->  Write();
        for ( auto kSave : h1D_RecInvMass )     kSave   ->  Write();
        for ( auto kSave : h1D_InvMassRes )     kSave   ->  Write();
        for ( auto kSave : h1D_TruInvMass_2Db ) kSave   ->  Write();
        for ( auto kSave : h1D_RecInvMass_2Db ) kSave   ->  Write();
        for ( auto kSave : h1D_InvMassRes_2Db ) kSave   ->  Write();
        //
        // --- Printing to Plots
        //
        outFile_Yield           ->  Close();
    }
    //  --- --- MULTIPLICITY ANALYSIS
    if ( kDoMultiplicity ) {
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_PreProc_Dir,(TString("Multiplicity/")+kFolder).Data())));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Multiplicity/")+kFolder).Data()))+TString("/Plots/1D/")).Data()));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Multiplicity/")+kFolder).Data()))+TString("/Plots/2D/")).Data()));
        TFile *outFile_Yield    =   new TFile   (Form(kAnalysis_MCTruthHist,(TString("Multiplicity/")+kFolder).Data()),"recreate");
        //
        // --- Saving to File
        fHEventCount    ->  Write();
        fHEvCountMlt    ->  Write();
        h1D_Nrec        ->  Write();
        h1D_Ngen        ->  Write();
        h1D_Ntru        ->  Write();
        h1D_Nrec_2Db    ->  Write();
        h1D_Ngen_2Db    ->  Write();
        h1D_Ntru_2Db    ->  Write();
        h1D_Nrec_Fin    ->  Write();
        h1D_Ngen_Fin    ->  Write();
        h1D_Ntru_Fin    ->  Write();
        h2D_Nrec        ->  Write();
        h2D_Ngen        ->  Write();
        h2D_Ntru        ->  Write();
        for ( auto kSave : h1D_TruInvMass )     kSave   ->  Write();
        for ( auto kSave : h1D_RecInvMass )     kSave   ->  Write();
        for ( auto kSave : h1D_InvMassRes )     kSave   ->  Write();
        for ( auto kSave : h1D_TruInvMass_2Db ) kSave   ->  Write();
        for ( auto kSave : h1D_RecInvMass_2Db ) kSave   ->  Write();
        for ( auto kSave : h1D_InvMassRes_2Db ) kSave   ->  Write();
        //
        // --- Printing to Plots
        //
        outFile_Yield           ->  Close();
    }
    //  --- --- CORRELATION ANALYSIS
    if ( kDoCorrelation ) {
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_PreProc_Dir,(TString("Correlation/")+kFolder).Data())));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Correlation/")+kFolder).Data()))+TString("/Plots/1D/")).Data()));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Correlation/")+kFolder).Data()))+TString("/Plots/2D/")).Data()));
        TFile *outFile_Yield    =   new TFile   (Form(kAnalysis_MCTruthHist,(TString("Correlation/")+kFolder).Data()),"recreate");
        //
        //  --- Normalisation
        h1D_Nrec_CrPh   ->  Scale(1.,"width");
        h1D_Ngen_CrPh   ->  Scale(1.,"width");
        h1D_Ntru_CrPh   ->  Scale(1.,"width");
        h1D_Nrec_CrPh   ->  Scale(1./kEventNormalisation);
        h1D_Ngen_CrPh   ->  Scale(1./kEventNormalisation);
        h1D_Ntru_CrPh   ->  Scale(1./kEventNormalisation);
        
        //
        //  --- Saving to File
        fHEventCount    ->  Write();
        fHEvCountMlt    ->  Write();
        h1D_Nrec        ->  Write();
        h1D_Ngen        ->  Write();
        h1D_Ntru        ->  Write();
        h1D_Nrec_2Db    ->  Write();
        h1D_Ngen_2Db    ->  Write();
        h1D_Ntru_2Db    ->  Write();
        h1D_Nrec_Fin    ->  Write();
        h1D_Ngen_Fin    ->  Write();
        h1D_Ntru_Fin    ->  Write();
        h1D_Nrec_CrPh   ->  Write();
        h1D_Ngen_CrPh   ->  Write();
        h1D_Ntru_CrPh   ->  Write();
        h2D_Nrec        ->  Write();
        h2D_Ngen        ->  Write();
        h2D_Ntru        ->  Write();
        for ( auto kSave : h1D_TruInvMass )     kSave   ->  Write();
        for ( auto kSave : h1D_RecInvMass )     kSave   ->  Write();
        for ( auto kSave : h1D_InvMassRes )     kSave   ->  Write();
        for ( auto kSave : h1D_TruInvMass_2Db ) kSave   ->  Write();
        for ( auto kSave : h1D_RecInvMass_2Db ) kSave   ->  Write();
        for ( auto kSave : h1D_InvMassRes_2Db ) kSave   ->  Write();
        //
        // --- Printing to Plots
        //
        outFile_Yield           ->  Close();
    }
    //
    insFile_Data->Close();
}
