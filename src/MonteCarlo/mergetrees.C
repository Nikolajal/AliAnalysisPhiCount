void
mergetrees
 ( TString  fInputFile, TString fOutputFile, Int_t iMode = 0, Int_t kEnergy = 0 ) {
    TFile*  kInputFile  =   new TFile(fInputFile);
    if ( fOutputFile.IsNull() ) fOutputFile = fInputFile + TString(".merged.root");
    TFile*  kOutputFile  =   new TFile(fOutputFile,"RECREATE");
    TList*  Trees = new TList();
    auto iTer = 0;
    while ( true ) {
        iTer++;
        if ( !(kInputFile->Get(Form("Prt_S%i_E%i_M%i",iTer,kEnergy,iMode))) ) break;
        auto knewtarget = ((TTree*)(kInputFile->Get(Form("Prt_S%i_E%i_M%i",iTer,kEnergy,iMode))))->CloneTree();
        Trees->Add(knewtarget);
    }
    TTree*  kOutTree = TTree::MergeTrees(Trees);
    kOutTree->SetName(Form("Prt_SX_E%i_M%i",kEnergy,iMode));
    cout << kOutTree->GetEntries() << endl;
    (kInputFile->Get("kEventCount"))->Write();
    kOutTree->Write();
    kInputFile->Close();
    kOutputFile->Close();
}
