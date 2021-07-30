#include "../../inc/AliAnalysisPhiPair.h"
// !TODO: All Set!

void MaterialBudget ( string fFileNameHig = "", string fFileNameLow = "", Int_t nEventsCut = -1., TString fOption = "" )   {
    //---------------------//
    //  Setting up input   //
    //---------------------//
    //
    // >-> Initialisation warnings
    //
    if ( fFileNameHig == "" || fFileNameHig == "" )  {
        cout << "[WARNING] Must Specify two input root file" << endl;
        cout << "[INFO] Usage MaterialBudget.cxx(\"RootFile_High.root\",\"RootFile_Low.root\")" << endl;
        return;
    }
    if ( nEventsCut != -1 ) cout << "[WARNING] Choosing to limit the datasample to " << nEventsCut << " events" <<endl;
    fChooseOption(fOption);
    
    //Retrieving Event data
    TFile *insFileMC_Hig    =   new TFile   (fFileNameHig.c_str());
    TFile *insFileMC_Low    =   new TFile   (fFileNameLow.c_str());
    
    //Retrieving Event data TTree
    TTree   *TPhiCandidateH =   (TTree*)insFileMC_Hig->Get(fPhiCandidateEff_Tree+TString("_name"));
    TTree   *TKaonCandidatH =   nullptr;//(TTree*)insFileMC_Hig->Get(fKaonCandidateEff_Tree);
    TTree   *TPhiCandidateL =   (TTree*)insFileMC_Low->Get(fPhiCandidateEff_Tree+TString(/*"_name"*/"_PhiPhiAnalysis_0_3_5_3_PhiCandidate"));
    TTree   *TKaonCandidatL =   nullptr;//(TTree*)insFileMC_Low->Get(fKaonCandidateEff_Tree);
    
    // Define tree data structures
    Struct_PhiEfficiency    evPhiEfficiencyH;
    Struct_KaonEfficiency   evKaonEfficiencH;
    Struct_PhiEfficiency    evPhiEfficiencyL;
    Struct_KaonEfficiency   evKaonEfficiencL;

    // Setting the input Candidates in the Trees
    if ( !fSetCandidates(TPhiCandidateH,evPhiEfficiencyH,TKaonCandidatH,evKaonEfficiencH) ); //return;
    if ( !fSetCandidates(TPhiCandidateL,evPhiEfficiencyL,TKaonCandidatL,evKaonEfficiencL) ); //return;
    
    //---------------------//
    //  Setting up output  //
    //---------------------//
    
    // Generating the binning array--------------------------------------------------------------------------
    fSetAllBins ();
    auto    fNBinning = (int)((fMaxPT1D-fMinPT1D)*10.);
    Float_t *fUniformBinning100MeV = new Float_t[fNBinning+1];
    fSetUniformBinning(fUniformBinning100MeV,fMinPT1D,fMaxPT1D,fNBinning);
    Int_t       U_AccCand[1024];
    Int_t       U_nAccept,  U_nAccep2;
    
    // Creating the histograms-------------------------------------------------------------------------------

    // >> YIELD ANALYSIS //

    // >>-->> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F       *h1D_REC_Hig;
    TH1F       *h1D_REC_Low;
    TH1F       *h1D_Gen_Hig;
    TH1F       *h1D_Gen_Low;
    TH1F       *h1D_Eff_Hig;
    TH1F       *h1D_Eff_Low;
    //
    //  Defining Efficiency and check utilities
    //
    hName       =   Form("h1D_REC_Hig");
    hTitle      =   Form("h1D_REC_Hig");
    h1D_REC_Hig =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(h1D_REC_Hig,"PT 1D");
    //
    hName       =   Form("h1D_REC_Low");
    hTitle      =   Form("h1D_REC_Low");
    h1D_REC_Low =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(h1D_REC_Low,"PT 1D");
    //
    hName       =   Form("h1D_Gen_Hig");
    hTitle      =   Form("h1D_Gen_Hig");
    h1D_Gen_Hig =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(h1D_Gen_Hig,"PT 1D");
    //
    hName       =   Form("h1D_Gen_Low");
    hTitle      =   Form("h1D_Gen_Low");
    h1D_Gen_Low =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(h1D_Gen_Low,"PT 1D");
    //
    hName       =   Form("h1D_Eff_Hig");
    hTitle      =   Form("h1D_Eff_Hig");
    h1D_Eff_Hig =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(h1D_Eff_Hig,"PT 1D");
    //
    hName       =   Form("h1D_Eff_Low");
    hTitle      =   Form("h1D_Eff_Low");
    h1D_Eff_Low =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(h1D_Eff_Low,"PT 1D");

    //-------------------------//
    //  Filling output objects //
    //-------------------------//
    
    TString fNameStringHig  =   TString("Reading High Material File");
    fStartTimer(fNameStringHig);
    
    // Evaluating entries
    Int_t nEvents = (!TPhiCandidateH) ? 0 : ( nEventsCut == -1.? TPhiCandidateH->GetEntries() : ( nEventsCut > TPhiCandidateH->GetEntries() ? TPhiCandidateH->GetEntries() : nEventsCut ) );
    
    // Starting cycle
    for ( Int_t iEvent = 0; iEvent < nEvents; iEvent++ )    {
        // Recovering events
        TPhiCandidateH->GetEntry(iEvent);
        
        fPrintLoopTimer(fNameStringHig,iEvent,nEvents,kPrintIntervalPP);

        // Utilities
        TLorentzVector  LPhi_candidate1,    LPhi_candidate2;
        U_nAccept = 0;
        
        for ( Int_t iPhi = 0; iPhi < evPhiEfficiencyH.nPhi; iPhi++ ) {
            LPhi_candidate1.SetXYZM(evPhiEfficiencyH.Px[iPhi],evPhiEfficiencyH.Py[iPhi],evPhiEfficiencyH.Pz[iPhi],kPhiMesonMass_);
            if ( !fAcceptCandidate(kPhiMesonMass_,LPhi_candidate1.Pt()) ) continue;
            U_AccCand[U_nAccept] = iPhi;
            U_nAccept++;
        }
        for ( Int_t iPhi = 0; iPhi < U_nAccept; iPhi++ )    {
            // Must have at least 1 candidate
            if ( U_nAccept < 1 ) break;

            // Building First Candidate
            LPhi_candidate1.SetXYZM(evPhiEfficiencyH.Px[U_AccCand[iPhi]],evPhiEfficiencyH.Py[U_AccCand[iPhi]],evPhiEfficiencyH.Pz[U_AccCand[iPhi]],kPhiMesonMass_);

            // >> 1-Dimensional Analysis Fill
            //
            // >>-->>-->> True Phis
            //
            Int_t   iSelection          =   (int)evPhiEfficiencyH.Selection[U_AccCand[iPhi]];
            Bool_t  iIsGen              =   (iSelection >= 1);
            Bool_t  iIsRec              =   (iSelection >= 2);
            Float_t iTransMom           =   LPhi_candidate1.Pt();
            Bool_t  fHasRapidity        =   fabs(LPhi_candidate1.Rapidity()) <0.5;
            //
            if ( iIsGen && fHasRapidity )   {
                h1D_Gen_Hig                             ->  Fill(iTransMom);
            } if ( iIsRec&& fHasRapidity )  {
                h1D_REC_Hig                             ->  Fill(iTransMom);
            }
        }
    }
    
    fStopTimer(fNameStringHig);
    
    TString fNameStringLow  =   TString("Reading Low Material File");
    fStartTimer(fNameStringLow);
    
    // Evaluating entries
    nEvents = (!TPhiCandidateL) ? 0 : ( nEventsCut == -1.? TPhiCandidateL->GetEntries() : nEventsCut);
    
    // Starting cycle
    for ( Int_t iEvent = 0; iEvent < nEvents; iEvent++ )    {
        // Recovering events
        TPhiCandidateL->GetEntry(iEvent);
        
        fPrintLoopTimer(fNameStringLow,iEvent,nEvents,kPrintIntervalPP);

        // Utilities
        TLorentzVector  LPhi_candidate1,    LPhi_candidate2;
        U_nAccept = 0;
        
        for ( Int_t iPhi = 0; iPhi < evPhiEfficiencyL.nPhi; iPhi++ ) {
            LPhi_candidate1.SetXYZM(evPhiEfficiencyL.Px[iPhi],evPhiEfficiencyL.Py[iPhi],evPhiEfficiencyL.Pz[iPhi],kPhiMesonMass_);
            if ( !fAcceptCandidate(kPhiMesonMass_,LPhi_candidate1.Pt()) ) continue;
            U_AccCand[U_nAccept] = iPhi;
            U_nAccept++;
        }
        for ( Int_t iPhi = 0; iPhi < U_nAccept; iPhi++ )    {
            // Must have at least 1 candidate
            if ( U_nAccept < 1 ) break;

            // Building First Candidate
            LPhi_candidate1.SetXYZM(evPhiEfficiencyL.Px[U_AccCand[iPhi]],evPhiEfficiencyL.Py[U_AccCand[iPhi]],evPhiEfficiencyL.Pz[U_AccCand[iPhi]],kPhiMesonMass_);

            // >> 1-Dimensional Analysis Fill
            //
            // >>-->>-->> True Phis
            //
            Int_t   iSelection          =   (int)evPhiEfficiencyL.Selection[U_AccCand[iPhi]];
            Bool_t  iIsGen              =   (iSelection >= 1);
            Bool_t  iIsRec              =   (iSelection >= 2);
            Float_t iTransMom           =   LPhi_candidate1.Pt();
            Bool_t  fHasRapidity        =   fabs(LPhi_candidate1.Rapidity()) <0.5;
            //
            if ( iIsGen && fHasRapidity )   {
                h1D_Gen_Low                             ->  Fill(iTransMom);
            } if ( iIsRec && fHasRapidity )  {
                h1D_REC_Low                             ->  Fill(iTransMom);
            }
        }
    }
    
    fStopTimer(fNameStringLow);
    
    h1D_Eff_Hig                             ->Divide(h1D_REC_Hig,           h1D_Gen_Hig,            1.,1.,"b");
    h1D_Eff_Low                             ->Divide(h1D_REC_Low,           h1D_Gen_Low,            1.,1.,"b");
    
    TCanvas*    c1  =   new TCanvas("","",1000,500);
    c1->Divide(2,1);
    c1->cd(1);
    h1D_Eff_Low->SetLineColor(kRed);
    h1D_Eff_Low->SetMarkerColor(kRed);
    h1D_Eff_Low->DrawCopy();
    h1D_Eff_Hig->SetLineColor(kBlue);
    h1D_Eff_Hig->SetMarkerColor(kBlue);
    h1D_Eff_Hig->DrawCopy("same");
    c1->cd(2);
    TH1F * hHI = new TH1F(*h1D_Eff_Hig);
    h1D_Eff_Hig->Divide(h1D_Eff_Low,h1D_Eff_Hig);
    for ( Int_t iter = 1; iter <= h1D_Eff_Hig->GetNbinsX(); iter++ ) {
        hHI->SetBinContent(iter,100*fabs(h1D_Eff_Hig->GetBinContent(iter)-1));
    }
    hHI->Draw();
    gPad->BuildLegend();
}
    /*
    //--------------------------//
    // PostProcessin output obj //
    //--------------------------//
    //
    // >> YIELD ANALYSIS //
    //
    auto fNormEvent = fHEventCount->GetBinContent(1);
    hEFF_1D                             ->Divide(hREC_1D,           hGEN_1D,            1.,1.,"b");
    hEFF_1D_in_2Dbin                    ->Divide(hREC_1D_in_2Dbin,  hGEN_1D_in_2Dbin,   1.,1.,"b");
    hEFF_2D                             ->Divide(hREC_2D,           hGEN_2D,            1.,1.,"b");
    hREC_1D->Scale(1.,"width");
    hGEN_1D->Scale(1.,"width");
    hREC_Rw_1D->Scale(1.,"width");
    hGEN_Rw_1D->Scale(1.,"width");
    hTRU_1D->Scale(1.,"width");
    hTRU_RECVTX_1D->Scale(1.,"width");
    hTRU_ALLVTX_1D->Scale(1.,"width");
    hREC_1D_in_2Dbin->Scale(1.,"width");
    hGEN_1D_in_2Dbin->Scale(1.,"width");
    hTRU_1D_in_2Dbin->Scale(1.,"width");
    hREC_2D->Scale(1.,"width");
    hGEN_2D->Scale(1.,"width");
    hTRU_2D->Scale(1.,"width");
    hTRU_RECVTX_2D->Scale(1.,"width");
    hTRU_ALLVTX_2D->Scale(1.,"width");
    hREC_1D->Scale(1./fNormEvent);
    hGEN_1D->Scale(1./fNormEvent);
    hTRU_1D->Scale(1./fNormEvent);
    hTRU_RECVTX_1D->Scale(1./fNormEvent);
    hTRU_ALLVTX_1D->Scale(1./fNormEvent);
    hREC_1D_in_2Dbin->Scale(1./fNormEvent);
    hGEN_1D_in_2Dbin->Scale(1./fNormEvent);
    hTRU_1D_in_2Dbin->Scale(1./fNormEvent);
    hREC_2D->Scale(1./fNormEvent);
    hGEN_2D->Scale(1./fNormEvent);
    hTRU_2D->Scale(1./fNormEvent);
    hTRU_RECVTX_2D->Scale(1./fNormEvent);
    hTRU_ALLVTX_2D->Scale(1./fNormEvent);
    //
    
    // >> MULTIPLICITY ANALYSIS //
    //
    for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )
    {
        hEFF_1D_in_MT[iMult]               ->Divide(hREC_1D_in_MT[iMult],         hGEN_1D_in_MT[iMult],          1.,1.,"b");
        hEFF_1D_in_MT_in_2Dbin[iMult]      ->Divide(hREC_1D_in_MT_in_2Dbin[iMult],hGEN_1D_in_MT_in_2Dbin[iMult], 1.,1.,"b");
        hEFF_2D_in_MT[iMult]               ->Divide(hREC_2D_in_MT[iMult],         hGEN_2D_in_MT[iMult],          1.,1.,"b");
        hREC_1D_in_MT[iMult]->Scale(1.,"width");
        hGEN_1D_in_MT[iMult]->Scale(1.,"width");
        hTRU_1D_in_MT[iMult]->Scale(1.,"width");
        hREC_1D_in_MT_in_2Dbin[iMult]->Scale(1.,"width");
        hGEN_1D_in_MT_in_2Dbin[iMult]->Scale(1.,"width");
        hTRU_1D_in_MT_in_2Dbin[iMult]->Scale(1.,"width");
        hREC_2D_in_MT[iMult]->Scale(1.,"width");
        hGEN_2D_in_MT[iMult]->Scale(1.,"width");
        hTRU_2D_in_MT[iMult]->Scale(1.,"width");
        hREC_1D_in_MT[iMult]->Scale(1./fNormEvent);
        hGEN_1D_in_MT[iMult]->Scale(1./fNormEvent);
        hTRU_1D_in_MT[iMult]->Scale(1./fNormEvent);
        hREC_1D_in_MT_in_2Dbin[iMult]->Scale(1./fNormEvent);
        hGEN_1D_in_MT_in_2Dbin[iMult]->Scale(1./fNormEvent);
        hTRU_1D_in_MT_in_2Dbin[iMult]->Scale(1./fNormEvent);
        hREC_2D_in_MT[iMult]->Scale(1./fNormEvent);
        hGEN_2D_in_MT[iMult]->Scale(1./fNormEvent);
        hTRU_2D_in_MT[iMult]->Scale(1./fNormEvent);
    }
    //
    
    // >> RAPIDITY ANALYSIS //
    //
    for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ )
    {
        hEFF_1D_in_RP[iRap]               ->Divide(hREC_1D_in_RP[iRap],         hGEN_1D_in_RP[iRap],          1.,1.,"b");
        hEFF_1D_in_RP_in_2Dbin[iRap]      ->Divide(hREC_1D_in_RP_in_2Dbin[iRap],hGEN_1D_in_RP_in_2Dbin[iRap], 1.,1.,"b");
        hEFF_2D_in_RP[iRap]               ->Divide(hREC_2D_in_RP[iRap],         hGEN_2D_in_RP[iRap],          1.,1.,"b");
        hREC_1D_in_RP[iRap]->Scale(1.,"width");
        hGEN_1D_in_RP[iRap]->Scale(1.,"width");
        hTRU_1D_in_RP[iRap]->Scale(1.,"width");
        hREC_1D_in_RP_in_2Dbin[iRap]->Scale(1.,"width");
        hGEN_1D_in_RP_in_2Dbin[iRap]->Scale(1.,"width");
        hTRU_1D_in_RP_in_2Dbin[iRap]->Scale(1.,"width");
        hREC_2D_in_RP[iRap]->Scale(1.,"width");
        hGEN_2D_in_RP[iRap]->Scale(1.,"width");
        hTRU_2D_in_RP[iRap]->Scale(1.,"width");
        hREC_1D_in_RP[iRap]->Scale(1./fNormEvent);
        hGEN_1D_in_RP[iRap]->Scale(1./fNormEvent);
        hTRU_1D_in_RP[iRap]->Scale(1./fNormEvent);
        hREC_1D_in_RP_in_2Dbin[iRap]->Scale(1./fNormEvent);
        hGEN_1D_in_RP_in_2Dbin[iRap]->Scale(1./fNormEvent);
        hTRU_1D_in_RP_in_2Dbin[iRap]->Scale(1./fNormEvent);
        hREC_2D_in_RP[iRap]->Scale(1./fNormEvent);
        hGEN_2D_in_RP[iRap]->Scale(1./fNormEvent);
        hTRU_2D_in_RP[iRap]->Scale(1./fNormEvent);
    }
    //
    
    // >> TRIGGER ANALYSIS //
    
    //--------------------------//
    //  Printing output objects //
    //--------------------------//
    //
    // >> Trigger Analysis
    //
    if ( kDoTrigger )   {
        TFile *outFil1  =   new TFile   (fTrgPrePrMC,"recreate");
        //
        outFil1->Close();
    }
    //
    // >> Yield Analysis
    //
    if ( kDoYield || kDoRapidity ) {
        TFile *outFil2  =   new TFile   (fYldPrePrMC,"recreate");
        //
        hREC_1D->Write();
        hGEN_1D->Write();
        hREC_Rw_1D->Write();
        hGEN_Rw_1D->Write();
        hTRU_1D->Write();
        hTRU_RECVTX_1D->Write();
        hTRU_ALLVTX_1D->Write();
        hEFF_1D->Write();
        hREC_1D_in_2Dbin->Write();
        hGEN_1D_in_2Dbin->Write();
        hTRU_1D_in_2Dbin->Write();
        hEFF_1D_in_2Dbin->Write();
        hREC_2D->Write();
        hGEN_2D->Write();
        hTRU_2D->Write();
        hTRU_RECVTX_2D->Write();
        hTRU_ALLVTX_2D->Write();
        hEFF_2D->Write();
        //
        outFil2->Close();
    }
    //
    // >> Multiplicity Analysis
    //
    if ( kDoMultiplicity )  {
        TFile *outFil3  =   new TFile   (fMltPrePrMC,"recreate");
        //
        for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )
        {
            hREC_1D_in_MT[iMult]->Write();
            hGEN_1D_in_MT[iMult]->Write();
            hTRU_1D_in_MT[iMult]->Write();
            hEFF_1D_in_MT[iMult]->Write();
            hREC_1D_in_MT_in_2Dbin[iMult]->Write();
            hGEN_1D_in_MT_in_2Dbin[iMult]->Write();
            hTRU_1D_in_MT_in_2Dbin[iMult]->Write();
            hEFF_1D_in_MT_in_2Dbin[iMult]->Write();
            hREC_2D_in_MT[iMult]->Write();
            hGEN_2D_in_MT[iMult]->Write();
            hTRU_2D_in_MT[iMult]->Write();
            hEFF_2D_in_MT[iMult]->Write();
        }
        //
        outFil3->Close();
    }
    //
    if ( kDoRapidity )  {
        TFile *outFil4  =   new TFile   (fRapPrePrMC,"recreate");
        //
        for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ )
        {
            hREC_1D_in_RP[iRap]->Write();
            hGEN_1D_in_RP[iRap]->Write();
            hTRU_1D_in_RP[iRap]->Write();
            hEFF_1D_in_RP[iRap]->Write();
            hREC_1D_in_RP_in_2Dbin[iRap]->Write();
            hGEN_1D_in_RP_in_2Dbin[iRap]->Write();
            hTRU_1D_in_RP_in_2Dbin[iRap]->Write();
            hEFF_1D_in_RP_in_2Dbin[iRap]->Write();
            hREC_2D_in_RP[iRap]->Write();
            hGEN_2D_in_RP[iRap]->Write();
            hTRU_2D_in_RP[iRap]->Write();
            hEFF_2D_in_RP[iRap]->Write();
        }
        //
        outFil4->Close();
    }
    //
    // >-> Close input File
    //
    insFileMC->Close();
    //
}

*/
