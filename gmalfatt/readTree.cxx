void
readTree
() {
    TFile*  kInputFile  = new TFile("/Users/nikolajal/alice/AliAnalysisPhiCount/MALFATTORE_GIOVANNI_TEST_XXX_VACCINO_5G/AnalysisResults.MC.root");
    //
    auto    kTreeList   =   (TDirectoryFile*)(kInputFile->Get("deuthelium"));
    auto    kTreeRecn   =   (TTree*) kTreeList->FindObjectAny("fDeuteronsTree");
    auto    kTreeTrue   =   (TTree*) kTreeList->FindObjectAny("fMCDeuteronsTree");
    //
    Int_t       fnD_tru,        fPDGd[1024];
    Float_t     fDPx_tru[1024], fDPy_tru[1024], fDPz_tru[1024], fPT_rec[1024];
    Int_t       fnD_rec;
    Float_t     fDPx_rec[1024], fDPy_rec[1024], fDPz_rec[1024];
    //
    kTreeTrue->SetBranchAddress("fnDTrue",      &fnD_tru );
    kTreeTrue->SetBranchAddress("fDPx",         &fDPx_tru );
    kTreeTrue->SetBranchAddress("fDPy",         &fDPy_tru );
    kTreeTrue->SetBranchAddress("fDPz",         &fDPz_tru );
    //
    kTreeRecn->SetBranchAddress("fnD",          &fnD_rec );
    kTreeRecn->SetBranchAddress("fPDGd",        &fPDGd );
    kTreeRecn->SetBranchAddress("fDPx",         &fDPx_rec );
    kTreeRecn->SetBranchAddress("fDPy",         &fDPy_rec );
    kTreeRecn->SetBranchAddress("fDPz",         &fDPz_rec );
    kTreeRecn->SetBranchAddress("fDCAxyD",         &fPT_rec );
    //
    TH1F*   hTru    =   new TH1F("hTru","hTru",100,0,10);
    TH1F*   hRec    =   new TH1F("hRec","hRec",100,0,10);
    TH1F*   hEff    =   new TH1F("hEff","hEff",100,0,10);
    //
    for ( Int_t iEv = 0; iEv < kTreeTrue->GetEntries(); iEv++ ) {
        kTreeTrue->GetEntry(iEv);
        //
        for ( Int_t iDt = 0; iDt < fnD_tru; iDt++ ) {
            TLorentzVector  kCurrentDeuteron;
            kCurrentDeuteron.SetXYZM( fDPx_tru[iDt], fDPy_tru[iDt], fDPz_tru[iDt], 1.875 );
            //if ( fabs(kCurrentDeuteron.Rapidity()) < 0.5 ) continue;
            hTru->Fill( kCurrentDeuteron.Pt() );
        }
        //
    }
    //
    for ( Int_t iEv = 0; iEv < kTreeRecn->GetEntries(); iEv++ ) {
        kTreeRecn->GetEntry(iEv);
        //
        for ( Int_t iDt = 0; iDt < fnD_rec; iDt++ ) {
            TLorentzVector  kCurrentDeuteron;
            kCurrentDeuteron.SetXYZM( fDPx_rec[iDt], fDPy_rec[iDt], fDPz_rec[iDt], 1.875 );
            //if ( fabs(kCurrentDeuteron.Rapidity()) < 0.5 ) continue;
            if ( fPDGd[iDt] != 1000010020 ) continue;
            hRec->Fill( fPT_rec[iDt] );
        }
        //
    }
    //
    hEff->Divide( hRec, hTru, 1., 1., "b" );
    //
    TCanvas* c1 = new TCanvas();
    //
    gPad->SetGridy();
    hEff->Draw("same");
    c1->SaveAs("hEff.pdf");
    //
    delete c1;
    //
    c1 = new TCanvas();
    //
    hRec->Draw("same");
    c1->SaveAs("hRec.pdf");
    //
    delete c1;
    //
    c1 = new TCanvas();
    //
    hTru->Draw("same");
    c1->SaveAs("hTru.pdf");
    //
    delete c1;
}
