#include "../../inc/AliAnalysisPhiPair.h"
// !TODO: AXIS and TITLES conforming to ALICE guidelines
// !TODO: Make the InvMass distinction for all options

void PP_Data ( TString fFileName = "/Volumes/[HD][Nikolajal]_Toshiba 2/Dataset/tmp/LHC15TeV_DT/LHC15TeV_DT_STD.root", TString fOption = "mult", Int_t nEventsCut = -1., TString kFolder = "" )    {
    // --- --- --- --- --- --- --- SET-UP --- --- --- --- --- --- --- --- --- --- ---
    //
    //  Check Correct Syntax
    if ( fFileName == "" )  {
        cout << "[ERROR] Must Specify an input root file" << endl;
        cout << "[INFO] Usage ( \"Filename\", \"Option\" , \"Eventcut\" , \"Analysis\" ) " << endl;
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
    TTree*      TPhiCandidate   =   (TTree*)insFile_Data->Get( fPhiCandidate_Tree + TString("") );
    //
    // --- Retrieving Utility Histograms
    TList*      fQCOutputList   =   (TList*)insFile_Data->Get( "fQCOutputList" );
    TH1D*       fHEventCount    =   (TH1D*) fQCOutputList->FindObject( "fQC_Event_Enum_FLL" );
    TH1D*       fHEvCountMlt    =   (TH1D*) fQCOutputList->FindObject( "fQC_Event_Enum_V0M" );
    //
    // --- Setting the input datastructure
    Struct_PhiCandidate     evPhiCandidate;
    Struct_KaonCandidate    kfix;           //TODO: Clean the fix
    if ( !fSetCandidates(TPhiCandidate,evPhiCandidate,nullptr,kfix) )    return;    // TODO: Fix the TrueInvMass in case of Data load
    //
    // --- Setting the output datastructure
    fSetAllBins();
    //
    // ---  --- YIELD ANALYSIS
    //
    TH2F*   hTest = new TH2F("test","test",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    TH1F*               h1D_Nrec;
    std::vector<TH1F*>  h1D_Nrec_PT;
    //
    hName           =   Form( "h1D_Nrec" );
    hTitle          =   Form( "m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c", fMinPT1D, fMaxPT1D );
    h1D_Nrec        =   new TH1F ( hName, hTitle, nBinIM1D, fArrIM1D );
    SetAxis( h1D_Nrec, "IM 1D" );
    //
    for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ ) {
        TH1F*               h1D_Nrec_PT_Util;
        hName               =   Form( "h1D_Nrec_PT_%i", iPT1D );
        hTitle              =   Form( "m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c", fArrPT1D[iPT1D], fArrPT1D[iPT1D+1] );
        h1D_Nrec_PT_Util    =   new TH1F ( hName, hTitle, nBinIM1D, fArrIM1D );
        SetAxis( h1D_Nrec_PT_Util, "IM 1D" );
        h1D_Nrec_PT.push_back( h1D_Nrec_PT_Util );
    }
    //
    TH2F*                               h2D_Nrec;
    std::vector<TH1F*>                  h1D_Nrec_2Db_PT;
    std::vector<std::vector<TH2F*>>     h2D_Nrec_PT;
    //
    hName           =   Form( "h2D_Nrec" );
    hTitle          =   Form( "m_{K^{+}K^{-}} in p_{T} range (X)[%.2f#;%.2f] (Y)[%.2f#;%.2f] GeV/c", fMinPT2D, fMaxPT2D, fMinPT2D, fMaxPT2D );
    h2D_Nrec        =   new TH2F ( hName, hTitle, nBinIM2D, fArrIM2D, nBinIM2D, fArrIM2D );
    SetAxis( h2D_Nrec, "IM 2D" );
    //
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ ) {
        TH1F*                   h1D_Nrec_2Db_PT_Util;
        hName                   =   Form( "h1D_Nrec_2Db_PT_%i", iPT2D );
        hTitle                  =   Form( "m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c", fArrPT2D[iPT2D], fArrPT2D[iPT2D+1] );
        h1D_Nrec_2Db_PT_Util    =   new TH1F ( hName, hTitle, nBinIM2D, fArrIM2D );
        SetAxis( h1D_Nrec_2Db_PT_Util, "IM 1D" );
        h1D_Nrec_2Db_PT.push_back( h1D_Nrec_2Db_PT_Util );
        std::vector<TH2F*>  h2D_Nrec_PT_VecUtil;
        for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ ) {
            TH2F*                   h2D_Nrec_PT_Util;
            hName                   =   Form( "h2D_Nrec_PT_%i_%i", iPT2D, jPT2D );
            hTitle                  =   Form( "m_{K^{+}K^{-}} in p_{T} range (X)[%.2f#;%.2f] (Y)[%.2f#;%.2f] GeV/c", fArrPT2D[iPT2D], fArrPT2D[iPT2D+1], fArrPT2D[jPT2D], fArrPT2D[jPT2D+1] );
            h2D_Nrec_PT_Util        =   new TH2F ( hName, hTitle, nBinIM2D, fArrIM2D, nBinIM2D, fArrIM2D );
            SetAxis( h2D_Nrec_PT_Util, "IM 2D" );
            h2D_Nrec_PT_VecUtil.push_back( h2D_Nrec_PT_Util );
        }
        h2D_Nrec_PT.push_back(h2D_Nrec_PT_VecUtil);
    }
    //
    // ---  --- MULTIPLICITY ANALYSIS
    //
    //  TODO: Make a Function for building the histograms for multiple differentiations
    std::vector<std::vector<TH1F*>>     h1D_Nrec_MT_PT;
    //
    for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ ) {
        std::vector<TH1F*>  h1D_Nrec_MT_PT_VecUtil;
        for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ ) {
            TH1F*               h1D_Nrec_MT_PT_Util;
            hName               =   Form( "h1D_Nrec_MT_%i_PT_%i", iMult, iPT1D );
            hTitle              =   Form( "m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c, Multiplicity [%.2f#;%.2f]", fArrPT1D[iPT1D], fArrPT1D[iPT1D+1], fArrMult[iMult], fArrMult[iMult+1] );
            h1D_Nrec_MT_PT_Util =   new TH1F ( hName, hTitle, nBinIM1D, fArrIM1D );
            SetAxis( h1D_Nrec_MT_PT_Util, "IM 1D" );
            h1D_Nrec_MT_PT_VecUtil.push_back( h1D_Nrec_MT_PT_Util );
        }
        h1D_Nrec_MT_PT.push_back( h1D_Nrec_MT_PT_VecUtil );
    }
    //
    std::vector<std::vector<TH1F*>>                 h1D_Nrec_2Db_MT_PT;
    std::vector<std::vector<std::vector<TH2F*>>>    h2D_Nrec_MT_PT;
    //
    for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ ) {
        std::vector<TH1F*>              h1D_2Db_Nrec_PT_MT_VecUtil;
        std::vector<std::vector<TH2F*>> h2D_Nrec_PT_MT_VecUtil;
        for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ ) {
            TH1F*                   h1D_Nrec_2Db_MT_PT_Util;
            hName                   =   Form( "h1D_Nrec_2Db_MT_%i_PT_%i", iMult, iPT2D );
            hTitle                  =   Form( "m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c, Multiplicity [%.2f#;%.2f]", fArrPT2D[iPT2D], fArrPT2D[iPT2D+1], fArrMult[iMult], fArrMult[iMult+1] );
            h1D_Nrec_2Db_MT_PT_Util =   new TH1F ( hName, hTitle, nBinIM2D, fArrIM2D );
            SetAxis( h1D_Nrec_2Db_MT_PT_Util, "IM 1D" );
            h1D_2Db_Nrec_PT_MT_VecUtil.push_back( h1D_Nrec_2Db_MT_PT_Util );
            std::vector<TH2F*>  h2D_Nrec_PT_VecUtil;
            for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ ) {
                TH2F*                   h2D_Nrec_MT_PT_Util;
                hName                   =   Form( "h2D_Nrec_MT_%i_PT_%i_%i", iMult, iPT2D, jPT2D );
                hTitle                  =   Form( "m_{K^{+}K^{-}} in p_{T} range (X)[%.2f#;%.2f] (Y)[%.2f#;%.2f] GeV/c, Multiplicity [%.2f#;%.2f]", fArrPT2D[iPT2D], fArrPT2D[iPT2D+1], fArrPT2D[jPT2D], fArrPT2D[jPT2D+1], fArrMult[iMult], fArrMult[iMult+1] );
                h2D_Nrec_MT_PT_Util     =   new TH2F ( hName, hTitle, nBinIM2D, fArrIM2D, nBinIM2D, fArrIM2D );
                SetAxis( h2D_Nrec_MT_PT_Util, "IM 2D" );
                h2D_Nrec_PT_VecUtil.push_back( h2D_Nrec_MT_PT_Util );
            }
            h2D_Nrec_PT_MT_VecUtil.push_back( h2D_Nrec_PT_VecUtil );
        }
        h1D_Nrec_2Db_MT_PT.push_back(h1D_2Db_Nrec_PT_MT_VecUtil);
        h2D_Nrec_MT_PT.push_back(h2D_Nrec_PT_MT_VecUtil);
    }
    //
    // ---  --- CORRELATION ANALYSIS
    //
    //  TODO: Make a Function for building the histograms for multiple differentiations
    std::vector<std::vector<TH1F*>>     h1D_Nrec_CR_PT;
    //
    for ( Int_t iCrPh = 0; iCrPh < nBinCrPh; iCrPh++ ) {
        std::vector<TH1F*>  h1D_Nrec_CR_PT_VecUtil;
        for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ ) {
            TH1F*               h1D_Nrec_CR_PT_Util;
            hName               =   Form( "h1D_Nrec_CR_%i_PT_%i", iCrPh, iPT1D );
            hTitle              =   Form( "m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c, Multiplicity [%.2f#;%.2f]", fArrPT1D[iPT1D], fArrPT1D[iPT1D+1], fArrCrPh[iCrPh], fArrCrPh[iCrPh+1] );
            h1D_Nrec_CR_PT_Util =   new TH1F ( hName, hTitle, nBinIM1D, fArrIM1D );
            SetAxis( h1D_Nrec_CR_PT_Util, "IM 1D" );
            h1D_Nrec_CR_PT_VecUtil.push_back( h1D_Nrec_CR_PT_Util );
        }
        h1D_Nrec_CR_PT.push_back( h1D_Nrec_CR_PT_VecUtil );
    }
    //
    std::vector<std::vector<TH1F*>>                 h1D_Nrec_2Db_CR_PT;
    std::vector<std::vector<std::vector<TH2F*>>>    h2D_Nrec_CR_PT;
    //
    for ( Int_t iCrPh = 0; iCrPh < nBinCrPh; iCrPh++ ) {
        std::vector<TH1F*>              h1D_2Db_Nrec_PT_CR_VecUtil;
        std::vector<std::vector<TH2F*>> h2D_Nrec_PT_CR_VecUtil;
        for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ ) {
            TH1F*                   h1D_Nrec_2Db_CR_PT_Util;
            hName                   =   Form( "h1D_Nrec_2Db_CR_%i_PT_%i", iCrPh, iPT2D );
            hTitle                  =   Form( "m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c, Multiplicity [%.2f#;%.2f]", fArrPT2D[iPT2D], fArrPT2D[iPT2D+1], fArrCrPh[iCrPh], fArrCrPh[iCrPh+1] );
            h1D_Nrec_2Db_CR_PT_Util =   new TH1F ( hName, hTitle, nBinIM2D, fArrIM2D );
            SetAxis( h1D_Nrec_2Db_CR_PT_Util, "IM 1D" );
            h1D_2Db_Nrec_PT_CR_VecUtil.push_back( h1D_Nrec_2Db_CR_PT_Util );
            std::vector<TH2F*>  h2D_Nrec_PT_VecUtil;
            for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ ) {
                TH2F*                   h2D_Nrec_CR_PT_Util;
                hName                   =   Form( "h2D_Nrec_CR_%i_PT_%i_%i", iCrPh, iPT2D, jPT2D );
                hTitle                  =   Form( "m_{K^{+}K^{-}} in p_{T} range (X)[%.2f#;%.2f] (Y)[%.2f#;%.2f] GeV/c, Multiplicity [%.2f#;%.2f]", fArrPT2D[iPT2D], fArrPT2D[iPT2D+1], fArrPT2D[jPT2D], fArrPT2D[jPT2D+1], fArrCrPh[iCrPh], fArrCrPh[iCrPh+1] );
                h2D_Nrec_CR_PT_Util     =   new TH2F ( hName, hTitle, nBinIM2D, fArrIM2D, nBinIM2D, fArrIM2D );
                SetAxis( h2D_Nrec_CR_PT_Util, "IM 2D" );
                h2D_Nrec_PT_VecUtil.push_back( h2D_Nrec_CR_PT_Util );
            }
            h2D_Nrec_PT_CR_VecUtil.push_back( h2D_Nrec_PT_VecUtil );
        }
        h1D_Nrec_2Db_CR_PT.push_back(h1D_2Db_Nrec_PT_CR_VecUtil);
        h2D_Nrec_CR_PT.push_back(h2D_Nrec_PT_CR_VecUtil);
    }
    //
    //
    // --- --- --- --- --- --- --- ANALYSIS --- --- --- --- --- --- --- --- --- --- -
    //
    //  --- Applying limitation sample cuts
    Int_t nEvents = (!TPhiCandidate) ? 0 : ( nEventsCut == -1.? TPhiCandidate->GetEntries() : nEventsCut);
    if ( nEvents > 0 )  fStartTimer("Phi Yield Analysis");
    //
    //  --- Start Analysis cycle
    for ( Int_t iEvent = 0; iEvent < nEvents; ++iEvent )    {
        TPhiCandidate->GetEntry(iEvent);
        fPrintLoopTimer("Phi Yield Analysis",iEvent,nEvents,kPrintIntervalPP);
        //
        //  --- Evaluate the Event is of interest for the analysis
        //
        //  --- Discarding Pile-up events
        if ( kCheckPileUp   &&  kDoYield           &&    fCheckMask(evPhiCandidate.EventMask,1)                                             ) continue;
        if ( kCheckPileUp   &&  kDoMultiplicity    &&  ( fCheckMask(evPhiCandidate.EventMask,1) || fCheckMask(evPhiCandidate.EventMask,2) ) ) continue;
        //
        //  --- Loop over candidates
        Struct_PhiCandidate fCurrent_Candidates;
        fCurrent_Candidates.nPhi    =   0;
        for ( Int_t iPhi = 0; iPhi < evPhiCandidate.nPhi; ++iPhi )  {
            TLorentzVector  kTLVUtility;   // TODO: Update to new lorentzvector class
            kTLVUtility.SetXYZM(evPhiCandidate.Px[iPhi],evPhiCandidate.Py[iPhi],evPhiCandidate.Pz[iPhi],evPhiCandidate.InvMass[iPhi]);
            if ( !fAcceptCandidate(evPhiCandidate.InvMass[iPhi],kTLVUtility.Pt()) ) continue;
            //  --- Load Useful Variables for histograms filling
            fCurrent_Candidates.InvMass     [ fCurrent_Candidates.nPhi ]    =   evPhiCandidate.InvMass[iPhi];
            fCurrent_Candidates.iKaon       [ fCurrent_Candidates.nPhi ]    =   evPhiCandidate.iKaon[iPhi];
            fCurrent_Candidates.jKaon       [ fCurrent_Candidates.nPhi ]    =   evPhiCandidate.jKaon[iPhi];
            fCurrent_Candidates.pT          [ fCurrent_Candidates.nPhi ]    =   kTLVUtility.Pt();
            fCurrent_Candidates.Phi         [ fCurrent_Candidates.nPhi ]    =   kTLVUtility.Phi()  * 360 / ( TMath::Pi() * 2 );
            fCurrent_Candidates.Rap         [ fCurrent_Candidates.nPhi ]    =   kTLVUtility.Rapidity();
            fCurrent_Candidates.iPT1D       [ fCurrent_Candidates.nPhi ]    =   fGetBinPT1D( kTLVUtility.Pt() );
            fCurrent_Candidates.iPT2D       [ fCurrent_Candidates.nPhi ]    =   fGetBinPT2D( kTLVUtility.Pt() );
            if ( is_pp_anl )    fCurrent_Candidates.kHasRap     [ fCurrent_Candidates.nPhi ]    =   fabs( kTLVUtility.Rapidity() ) < 0.5;
            if ( is_pb_anl )    fCurrent_Candidates.kHasRap     [ fCurrent_Candidates.nPhi ]    =   ( kTLVUtility.Rapidity() +0.465 ) < 0.5 && ( kTLVUtility.Rapidity() +0.465 ) > 0;
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
        fCurrent_Candidates.iMult           =   1+fGetBinMult(evPhiCandidate.Multiplicity);
        fCurrent_Candidates.kHasMult        =   fCurrent_Candidates.iMult != -1;
        //
        for ( Int_t iPhi = 0; iPhi < fCurrent_Candidates.nPhi; iPhi++ )  {
            //
            //  --- Selecting valid candidates
            if ( !fAcceptCandidate( fCurrent_Candidates, iPhi ) ) continue;
            //
            if ( fCurrent_Candidates.kHasRap[iPhi] )  {       //  --- Mid-Rapidity Analyses
                if ( kDoYield || kDoCorrelation ) {       //  --- YIELD ANALYSIS and CORRELATION UTILITY
                    h1D_Nrec                                                                                    ->  Fill( fCurrent_Candidates.InvMass[iPhi] );
                    h1D_Nrec_PT                                         .at( fCurrent_Candidates.iPT1D[iPhi] )  ->  Fill( fCurrent_Candidates.InvMass[iPhi] );
                    h1D_Nrec_2Db_PT                                     .at( fCurrent_Candidates.iPT2D[iPhi] )  ->  Fill( fCurrent_Candidates.InvMass[iPhi] );
                } if ( kDoMultiplicity && fCurrent_Candidates.kHasMult )    {   //  --- MULTIPLICITY ANALYSIS
                    h1D_Nrec_MT_PT      .at(0)                          .at( fCurrent_Candidates.iPT1D[iPhi] )  ->  Fill( fCurrent_Candidates.InvMass[iPhi] );
                    h1D_Nrec_2Db_MT_PT  .at(0)                          .at( fCurrent_Candidates.iPT2D[iPhi] )  ->  Fill( fCurrent_Candidates.InvMass[iPhi] );
                    h1D_Nrec_MT_PT      .at(fCurrent_Candidates.iMult)  .at( fCurrent_Candidates.iPT1D[iPhi] )  ->  Fill( fCurrent_Candidates.InvMass[iPhi] );
                    h1D_Nrec_2Db_MT_PT  .at(fCurrent_Candidates.iMult)  .at( fCurrent_Candidates.iPT2D[iPhi] )  ->  Fill( fCurrent_Candidates.InvMass[iPhi] );
                }
            }
            //
            //  --- Speed-up protection
            if ( fCurrent_Candidates.nPhi < 2 ) continue;
            for ( Int_t jPhi = iPhi+1; jPhi < fCurrent_Candidates.nPhi; jPhi++ )  {
                //
                //  --- Selecting valid candidates
                if ( !fAcceptCandidate( fCurrent_Candidates, iPhi, jPhi ) ) continue;
                //
                if ( fCurrent_Candidates.kHasRap[iPhi] && fCurrent_Candidates.kHasRap[jPhi] )  {    //  --- Mid-Rapidity Analyses
                    if ( kDoYield ) {   //  --- YIELD ANALYSIS
                        h2D_Nrec                                                                                                                        ->  Fill( fCurrent_Candidates.InvMass[iPhi], fCurrent_Candidates.InvMass[jPhi] );
                        h2D_Nrec_PT                                     .at( fCurrent_Candidates.iPT2D[iPhi] )  .at( fCurrent_Candidates.iPT2D[jPhi] )  ->  Fill( fCurrent_Candidates.InvMass[iPhi], fCurrent_Candidates.InvMass[jPhi] );
                        hTest   -> Fill( fCurrent_Candidates.pT[iPhi], fCurrent_Candidates.pT[jPhi] );
                    } if ( kDoMultiplicity && fCurrent_Candidates.kHasMult )    {   //  --- MULTIPLICITY ANALYSIS
                        h2D_Nrec_MT_PT  .at(0)                          .at( fCurrent_Candidates.iPT2D[iPhi] )  .at( fCurrent_Candidates.iPT2D[jPhi] )  ->  Fill( fCurrent_Candidates.InvMass[iPhi], fCurrent_Candidates.InvMass[jPhi] );
                        h2D_Nrec_MT_PT  .at(fCurrent_Candidates.iMult)  .at( fCurrent_Candidates.iPT2D[iPhi] )  .at( fCurrent_Candidates.iPT2D[jPhi] )  ->  Fill( fCurrent_Candidates.InvMass[iPhi], fCurrent_Candidates.InvMass[jPhi] );
                        
                    } if ( kDoCorrelation ) {
                        auto    kDeltaPhi   =   ( fCurrent_Candidates.Phi[iPhi] - fCurrent_Candidates.Phi[jPhi]  );
                        auto    iCrPh       =   fGetBinCrPh( kDeltaPhi < -180 ? kDeltaPhi + 360 : kDeltaPhi > 180 ? kDeltaPhi -360 : kDeltaPhi );
                        h2D_Nrec_CR_PT      .at(iCrPh)                  .at( fCurrent_Candidates.iPT2D[iPhi] )  .at( fCurrent_Candidates.iPT2D[jPhi] )  ->  Fill( fCurrent_Candidates.InvMass[iPhi], fCurrent_Candidates.InvMass[jPhi] );
                    }
                }
            }
        }
    }
    if ( nEvents > 0 )  fStopTimer("Phi Yield Analysis");
    //
    // --- --- --- --- --- --- --- OUTPUT --- --- --- --- --- --- --- --- --- --- ---
    //
    //  --- --- YIELD ANALYSIS
    if ( kDoYield ) {
        auto    kFolder1D       =   TString(Form(kAnalysis_PreProc_Dir,(TString("Yield")+kFolder).Data()))+TString("/Plots/1D/");
        auto    kFolder2D       =   TString(Form(kAnalysis_PreProc_Dir,(TString("Yield")+kFolder).Data()))+TString("/Plots/2D/");
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_PreProc_Dir,(TString("Yield")+kFolder).Data())));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",kFolder1D.Data()));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",kFolder2D.Data()));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s/InvMass/",kFolder1D.Data()));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s/InvMass/",kFolder2D.Data()));
        TFile *outFile_Yield    =   new TFile   (Form(kAnalysis_InvMassHist,(TString("Yield")+kFolder).Data()),"recreate");
        //
        // --- Saving to File
        hTest                   ->  Write();
        fHEventCount            ->  Write();
        fHEvCountMlt            ->  Write();
        h1D_Nrec                ->  Write();
        h2D_Nrec                ->  Write();
        for ( auto kSave    : h1D_Nrec_PT )                                     kSave   ->  Write();
        for ( auto kSave    : h1D_Nrec_2Db_PT )                                 kSave   ->  Write();
        for ( auto kVecSave : h2D_Nrec_PT )     for ( auto kSave : kVecSave )   kSave   ->  Write();
        outFile_Yield->Close();
        //
        // --- Printing to Plots
        // TODO: Make the same for 1D, generate for Mult histos
        for ( auto kPlot : h1D_Nrec_PT )                                    if ( kPlot->GetEntries() != 0 ) uPlotInvMass(kPlot,kFolder1D+TString("/InvMass/"));
        //  TODO: make the 2D equivalent for 1D histos
        //for ( auto kPlot : h1D_Nrec_2Db_PT )                                if ( kPlot->GetEntries() != 0 ) uPlotInvMass(kPlot,kFolder2D+TString("/InvMass/"));
        for ( auto kVecPlot : h2D_Nrec_PT ) for ( auto kPlot : kVecPlot )   if ( kPlot->GetEntries() != 0 ) uPlotInvMass(kPlot,kFolder2D+TString("/InvMass/"));
    }
    //  --- --- MULTIPLICITY ANALYSIS
    if ( kDoMultiplicity ) {
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_PreProc_Dir,(TString("Multiplicity")+kFolder).Data())));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Multiplicity")+kFolder).Data()))+TString("/Plots/1D/")).Data()));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Multiplicity")+kFolder).Data()))+TString("/Plots/2D/")).Data()));
        TFile *outFile_Yield    =   new TFile   (Form(kAnalysis_InvMassHist,(TString("Multiplicity")+kFolder).Data()),"recreate");
        //
        // --- Saving to File
        fHEventCount            ->  Write();
        fHEvCountMlt            ->  Write();
        for ( auto kVecSave : h1D_Nrec_MT_PT )      for ( auto kSave    : kVecSave )    kSave   ->  Write();
        for ( auto kVecSave : h1D_Nrec_2Db_MT_PT )  for ( auto kSave    : kVecSave )    kSave   ->  Write();
        for ( auto kVe2Save : h2D_Nrec_MT_PT )      for ( auto kVecSave : kVe2Save )    for ( auto kSave : kVecSave )   kSave   ->  Write();
        //
        outFile_Yield->Close();
        //
        // --- Printing to Plots
    }
    //  --- --- CORRELATION ANALYSIS
    if ( kDoCorrelation ) {
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_PreProc_Dir,(TString("Correlation")+kFolder).Data())));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Correlation")+kFolder).Data()))+TString("/Plots/1D/")).Data()));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Correlation")+kFolder).Data()))+TString("/Plots/2D/")).Data()));
        TFile *outFile_Yield    =   new TFile   (Form(kAnalysis_InvMassHist,(TString("Correlation")+kFolder).Data()),"recreate");
        //
        // --- Saving to File
        fHEventCount            ->  Write();
        fHEvCountMlt            ->  Write();
        for ( auto kSave    : h1D_Nrec_PT )                                                                             kSave   ->  Write();
        for ( auto kSave    : h1D_Nrec_2Db_PT )                                                                         kSave   ->  Write();
        for ( auto kVe2Save : h2D_Nrec_CR_PT )      for ( auto kVecSave : kVe2Save )    for ( auto kSave : kVecSave )   kSave   ->  Write();
        //
        outFile_Yield->Close();
        //
        // --- Printing to Plots
    }
    //
    insFile_Data->Close();
}
