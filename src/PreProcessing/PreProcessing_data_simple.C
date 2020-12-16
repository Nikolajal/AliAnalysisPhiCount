#include "../../inc/AliAnalysisPhiPair.h"
// !TODO: [INFO] About trees in input

void PreProcessing_data ( string fFileName = "" )
{
    if ( fFileName == "" )
    {
        cout << "[WARNING] Must Specify an input root file" << endl;
        cout << "[INFO] Usage PreProcessing_data(\"Root_file_name.root\")" << endl;
        return;
    }
    
    //Retrieving Event data
    TFile *insFileDT        =   new TFile   (fFileName.c_str());
    
    //Retrieving Event data TTree
    TTree   *TPhiCandidate  =   (TTree*)insFileDT->Get(fPhiCandidate_Tree);
    TTree   *TKaonCandidate =   (TTree*)insFileDT->Get(fKaonCandidate_Tree);
    
    // Define tree data structures
    Struct_PhiCandidate     evPhiCandidate;
        TPhiCandidate-> SetBranchAddress    ("fMultiplicity",    &evPhiCandidate.Multiplicity);
        TPhiCandidate-> SetBranchAddress    ("nPhi",            &evPhiCandidate.nPhi);
        TPhiCandidate-> SetBranchAddress    ("Px",              &evPhiCandidate.Px);
        TPhiCandidate-> SetBranchAddress    ("Py",              &evPhiCandidate.Py);
        TPhiCandidate-> SetBranchAddress    ("Pz",              &evPhiCandidate.Pz);
        TPhiCandidate-> SetBranchAddress    ("InvMass",         &evPhiCandidate.InvMass);
        TPhiCandidate-> SetBranchAddress    ("iKaon",           &evPhiCandidate.iKaon);
        TPhiCandidate-> SetBranchAddress    ("jKaon",           &evPhiCandidate.jKaon);

    
    Int_t nEvents = TPhiCandidate->GetEntries();

    // Starting cycle
    for ( Int_t iEvent = 0; iEvent < (int)(nEvents*0.01); iEvent++ )
    {
        // Recovering events
        TPhiCandidate->GetEntry(iEvent);
        
        if ( iEvent%1000000 == 0 && iEvent != 0) fPrintLoopTimer("Analysis",iEvent,nEvents);
        
        if ( (int)(evPhiCandidate.nPhi) != 0 )cout << "[INFO] " << (int)(evPhiCandidate.nPhi) << endl;
        if ( (int)(evPhiCandidate.nPhi) >= 153 )
        {
            cout << "[INFO] Skipping overflow event" << endl;
            cout << "[INFO] " << (int)(evPhiCandidate.nPhi) << endl;
            continue;
        }
    }
    insFileDT->Close();
}
