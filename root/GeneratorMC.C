#include "../inc/SetValues.h"
// Pythia
#include "Pythia8/Pythia.h"

int main (int argc, char *argv[])
{
    // Define some simple data structures
    EVKAONCOUPLE evKaonCouple;
    EVPHI evPhi;

    //Initialisation of TTree
    TFile * outFile = new TFile(argv[1],"recreate");

    TTree * PtreeK2     = new   TTree   (PTreeNameK2,"A ROOT tree for pythia MC - Kaon");
    TTree * PtreePhi    = new   TTree   (PTreeNamePhi,"A ROOT tree for pythia MC - Phi");

    // Filling Kaon Couple TTree
    PtreeK2->Branch    ("evKaonCouple.nKaonCouple"  ,&evKaonCouple.nKaonCouple, "nKaonCouple/I");
    PtreeK2->Branch    ("evKaonCouple.iKaon"        ,&evKaonCouple.iKaon,       "iKaon[nKaonCouple]/I");
    PtreeK2->Branch    ("evKaonCouple.jKaon"        ,&evKaonCouple.jKaon,       "jKaon[nKaonCouple]/I");
    PtreeK2->Branch    ("evKaonCouple.bPhi"         ,&evKaonCouple.bPhi,        "bPhi[nKaonCouple]/O");
    PtreeK2->Branch    ("evKaonCouple.bRec"         ,&evKaonCouple.bRec,        "bRec[nKaonCouple]/O");
    PtreeK2->Branch    ("evKaonCouple.InvMass"      ,&evKaonCouple.InvMass,     "InvMass[nKaonCouple]/F");
    PtreeK2->Branch    ("evKaonCouple.pT"           ,&evKaonCouple.pT,          "pT[nKaonCouple]/F");
    
    // Filling Phi TTree
    PtreePhi->Branch    ("evPhi.nPhi"               ,&evPhi.nPhi,               "nPhi/I");
    PtreePhi->Branch    ("evPhi.pT"                 ,&evPhi.pT,                 "pT[nPhi]/F");
    
    /* PYTHIA INITIALISATION ( Add string pointer for ... ) */
    Pythia8::Pythia pythia;
    
    //Settings
    pythia.readString("SoftQCD:nonDiffractive = on");
    pythia.readString("ParticleDecays:limitTau0 = on");
    pythia.readString("333:mMin = 0.99");
    pythia.readString("333:mMax = 1.09");
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = 0");
    pythia.init();
    
    // save the ID of kaons here
    int nKaon, nPhi, kaonID[1024], phiID[1024];
    bool kaonRec[1024];
    
    for (int iEvent = 0; iEvent < PEvents; iEvent++)
    {
        pythia.next();
        nPhi    = 0;
        nKaon   = 0;
        for (int iParticle = 0; iParticle < pythia.event.size() ; iParticle++)
        {
            const auto particle = pythia.event[iParticle];
            kaonRec[nKaon] = 1;
            // Counting Phi
            if (particle.id() == 333)
            {
                phiID[nPhi] = iParticle;
                nPhi++;
            }
            
            //Skipping non-final particles
            if(!particle.isFinal())                 continue;
            
            //Skipping particles in Eta non-acceptance region
            if (abs(particle.eta()) > 0.8)          kaonRec[nKaon] = 0;
            
            //Skipping particles in pT non-acceptance region
            if (particle.pT() < 0.15)               kaonRec[nKaon] = 0;
            
            //Getting Kaons only
            if (particle.id() == 321 || particle.id() == -321)
            {
                kaonID[nKaon] = iParticle;
                nKaon++;
            }
        }
        // now I create Kaon+- pairs TTree entry
        evKaonCouple.nKaonCouple = 0;
        for (int iKaon = 0; iKaon < nKaon; iKaon++)
        {
            const auto Kaon1 = pythia.event[kaonID[iKaon]];
            
            for (int jKaon = (iKaon+1); jKaon < nKaon; jKaon++)
            {
                const auto Kaon2 = pythia.event[kaonID[jKaon]];
                
                if ( Kaon1.id() != Kaon2.id() )
                {
                    const auto pPhi = Kaon1.p() + Kaon2.p();
                    evKaonCouple.iKaon[evKaonCouple.nKaonCouple]     = iKaon;
                    evKaonCouple.jKaon[evKaonCouple.nKaonCouple]     = jKaon;
                    evKaonCouple.bPhi[evKaonCouple.nKaonCouple]      = (Kaon1.mother2() == 0 && Kaon1.mother1() == Kaon2.mother1() && (pythia.event[Kaon1.mother1()]).id() == 333);
                    evKaonCouple.bRec[evKaonCouple.nKaonCouple]      = (kaonRec[iKaon] && kaonRec[jKaon]);
                    evKaonCouple.pT[evKaonCouple.nKaonCouple]        = pPhi.pT();
                    evKaonCouple.InvMass[evKaonCouple.nKaonCouple]   = pPhi.mCalc();
                    
                    //Cut on Invariant Mass
                    if (evKaonCouple.InvMass[evKaonCouple.nKaonCouple] > maxBound) continue;
                    //Cut on Rapidity
                    if (abs(pPhi.rap()) >= 0.5) continue;
                    evKaonCouple.nKaonCouple++;
                }
            }
        }
        evPhi.nPhi = 0;
        for (int iPhi = 0; iPhi < nPhi; iPhi++)
        {
            const auto pPhi                 = pythia.event[phiID[iPhi]];
            evPhi.pT[evPhi.nPhi]            = pPhi.pT();
            
            //Cut on Rapidity
            if (abs((pPhi.p()).rap()) >= 0.5) continue;
            evPhi.nPhi++;
        }
        //Filling TTrees
        if(nPhi  > 0) PtreeK2->Fill();
        if(nKaon > 0) PtreePhi->Fill();
    }

    PtreeK2     ->Write();
    PtreePhi    ->Write();
    outFile     ->Close();
    return 0;
}
