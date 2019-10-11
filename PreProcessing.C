#include "SetValues.h"

int main ()
{
    //Retrieving Event/MC data
    TFile *inFile = new TFile(outMC);
    //TFile *inFile = new TFile(outDF);
    
    TTree *Ptree = (TTree*)inFile->Get(PTreeName);
    
    // Define some simple data structures
    EVKAON evKaon;
    EVKAONCOUPLE evKaonCouple;
    
    //Setting Branch Addresses
    // Filling 1st branch
    /*
    Ptree  ->SetBranchAdress    ("evKaon.nKaon",&evKaon.nKaon,"nKaon/I");
    Ptree  ->SetBranchAdress    ("evKaon.particleID",&evKaon.particleID,"particleID[nKaon]/I");
    Ptree  ->SetBranchAdress    ("evKaon.mother1",&evKaon.mother1,"mother1[nKaon]/I");
    Ptree  ->SetBranchAdress    ("evKaon.mother2",&evKaon.mother2,"mother2[nKaon]/I");
    Ptree  ->SetBranchAdress    ("evKaon.motherID",&evKaon.motherID,"motherID[nKaon]/I");
    Ptree  ->SetBranchAdress    ("evKaon.px",&evKaon.px,"px[nKaon]/F");
    Ptree  ->SetBranchAdress    ("evKaon.py",&evKaon.py,"py[nKaon]/F");
    Ptree  ->SetBranchAdress    ("evKaon.pz",&evKaon.pz,"pz[nKaon]/F");
    Ptree  ->SetBranchAdress    ("evKaon.pT",&evKaon.pT,"pT[nKaon]/F");
    Ptree  ->SetBranchAdress    ("evKaon.e",&evKaon.e,"e[nKaon]/F");
    */
    // Filling 2nd branch
    Ptree  ->SetBranchAddress    ("evKaonCouple.nKaonCouple",&evKaonCouple.nKaonCouple);
    Ptree  ->SetBranchAddress    ("evKaonCouple.iKaon",&evKaonCouple.iKaon);
    Ptree  ->SetBranchAddress    ("evKaonCouple.jKaon",&evKaonCouple.jKaon);
    Ptree  ->SetBranchAddress    ("evKaonCouple.bPhi",&evKaonCouple.bPhi);
    Ptree  ->SetBranchAddress    ("evKaonCouple.InvMass",&evKaonCouple.InvMass);
    Ptree  ->SetBranchAddress    ("evKaonCouple.pT",&evKaonCouple.pz);
    /*
    Ptree  ->SetBranchAdress    ("evKaonCouple.px",&evKaonCouple.px,"px[nKaonCouple]/F");
    Ptree  ->SetBranchAdress    ("evKaonCouple.py",&evKaonCouple.py,"py[nKaonCouple]/F");
    Ptree  ->SetBranchAdress    ("evKaonCouple.pz",&evKaonCouple.pz,"pz[nKaonCouple]/F");
    Ptree  ->SetBranchAdress    ("evKaonCouple.e",&evKaonCouple.e,"e[nKaonCouple]/F");
    */
    
    //1D Map of Phi invariant mass
    TH1F * Kaon_InvM1D = new TH1F("Kaon_InvMass1D","Kaon_InvMass1D",nBins,minBound,maxBound);
    Kaon_InvM1D->SetOption("colz");
    Kaon_InvM1D->SetStats(0);
    
    //2D Map of Phi invariant mass
    TH2F * Kaon_InvM2D  = new TH2F("Kaon_InvMass2D","Kaon_InvMass2D",nBins,minBound,maxBound,nBins,minBound,maxBound);
    Kaon_InvM2D->SetOption("colz");
    Kaon_InvM2D->SetStats(0);
    
    //2D Map of Phi invariant mass
    TH3F * Kaon_InvM3D  = new TH3F("Kaon_InvMass3D","Kaon_InvMass3D",nBins,minBound,maxBound,nBins,minBound,maxBound,nBins,minBound,maxBound);
    Kaon_InvM3D->SetOption("colz");
    Kaon_InvM3D->SetStats(0);

    
    //Starting loop on Events
    Int_t nEvents   = Ptree->GetEntries();
    Int_t nKaons    = 0;
    
    int Phi2Dcount = 0;
    int Phi2DcountDouble = 0;
    int Phi2DBackGround = 0;
    
    for (Int_t iEvent = 0; iEvent < nEvents; iEvent++)
    {
        Ptree->GetEntry(iEvent);
        nKaons    = evKaonCouple.nKaonCouple;
        for (Int_t iKaonLoop = 0; iKaonLoop < nKaons; iKaonLoop++)
        {
            Kaon_InvM1D->Fill(evKaonCouple.InvMass[iKaonLoop]);
            for (Int_t jKaonLoop = 0; jKaonLoop < evKaonCouple.nKaonCouple; jKaonLoop++)
            {
                if( checkBool(evKaonCouple.iKaon[iKaonLoop],evKaonCouple.jKaon[iKaonLoop],evKaonCouple.iKaon[jKaonLoop],evKaonCouple.jKaon[jKaonLoop]))
                {
                    Phi2DBackGround++;
                    Kaon_InvM2D->Fill(evKaonCouple.InvMass[iKaonLoop],evKaonCouple.InvMass[jKaonLoop]);
                    if ( evKaonCouple.bPhi[iKaonLoop] == 1 || evKaonCouple.bPhi[jKaonLoop] == 1 )
                    {
                        Phi2Dcount++;
                        Phi2DBackGround--;
                        if ( evKaonCouple.bPhi[iKaonLoop] == evKaonCouple.bPhi[jKaonLoop])
                        {
                            Phi2DcountDouble++;
                        }
                    }
                }
            }
        }
    }
    
    cout << "xy Phi: " << Phi2Dcount << endl;
    cout << "xy Phidouble: " << Phi2DcountDouble << endl;
    cout << "xy Phibkg: " << Phi2DBackGround << endl;
    
    TFile *outFile = new TFile(outPP,"recreate");
    Kaon_InvM1D-> Write();
    Kaon_InvM2D-> Write();
    //Kaon_InvM3D-> Write();
    outFile-> Close();
    return 0;
}
