#include "../inc/SetValues.h"
// Pythia
#include "Pythia8/Pythia.h"

int main ()
{
    // Define some simple data structures
    EVKAONCOUPLE evKaonS;
    EVKAONCOUPLE evKaonD;

    //Initialisation of TTree
    
    TFile * outFile     = new   TFile   (oFilePreBKG,   "recreate");
    TTree * PtreeKS     = new   TTree   (PTreeNameKS,   "A ROOT tree for pythia MC - Kaon Couples");
    TTree * PtreeKD     = new   TTree   (PTreeNameKD,   "A ROOT tree for pythia MC - Kaon Couples");
    
    // Filling Kaon Couple TTree
    PtreeKS->Branch    ("evKaonCouple.nKaonCouple"  ,&evKaonS.nKaonCouple, "nKaonCouple/I");
    PtreeKS->Branch    ("evKaonCouple.iKaon"        ,&evKaonS.iKaon,       "iKaon[nKaonCouple]/I");
    PtreeKS->Branch    ("evKaonCouple.jKaon"        ,&evKaonS.jKaon,       "jKaon[nKaonCouple]/I");
    PtreeKS->Branch    ("evKaonCouple.bPhi"         ,&evKaonS.bPhi,        "bPhi[nKaonCouple]/O");
    PtreeKS->Branch    ("evKaonCouple.bRec"         ,&evKaonS.bRec,        "bRec[nKaonCouple]/O");
    PtreeKS->Branch    ("evKaonCouple.InvMass"      ,&evKaonS.InvMass,     "InvMass[nKaonCouple]/F");
    PtreeKS->Branch    ("evKaonCouple.pT"           ,&evKaonS.pT,          "pT[nKaonCouple]/F");
    
    PtreeKD->Branch    ("evKaonCouple.nKaonCouple"  ,&evKaonD.nKaonCouple, "nKaonCouple/I");
    PtreeKD->Branch    ("evKaonCouple.iKaon"        ,&evKaonD.iKaon,       "iKaon[nKaonCouple]/I");
    PtreeKD->Branch    ("evKaonCouple.jKaon"        ,&evKaonD.jKaon,       "jKaon[nKaonCouple]/I");
    PtreeKD->Branch    ("evKaonCouple.bPhi"         ,&evKaonD.bPhi,        "bPhi[nKaonCouple]/O");
    PtreeKD->Branch    ("evKaonCouple.bRec"         ,&evKaonD.bRec,        "bRec[nKaonCouple]/O");
    PtreeKD->Branch    ("evKaonCouple.InvMass"      ,&evKaonD.InvMass,     "InvMass[nKaonCouple]/F");
    PtreeKD->Branch    ("evKaonCouple.pT"           ,&evKaonD.pT,          "pT[nKaonCouple]/F");
    
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
        nKaon   = 0;
        for (int iParticle = 0; iParticle < pythia.event.size() ; iParticle++)
        {
            const auto particle = pythia.event[iParticle];
            kaonRec[nKaon] = 1;
            
            //Skipping non-final particles
            if(!particle.isFinal())                 continue;
            
            //Skipping particles in Eta non-acceptance region
            if (abs(particle.eta()) > 0.8)          kaonRec[nKaon] = false;
            
            //Skipping particles in pT non-acceptance region
            if (particle.pT() < 0.15)               kaonRec[nKaon] = false;
            
            //Getting Kaons only
            if (particle.id() == 321 || particle.id() == -321)
            {
                kaonID[nKaon] = iParticle;
                nKaon++;
            }
        }
        
        // now I create Kaon++-- pairs TTree entry
        evKaonS.nKaonCouple = 0;
        evKaonD.nKaonCouple = 0;
        for (int iKaon = 0; iKaon < nKaon; iKaon++)
        {
            const auto Kaon1 = pythia.event[kaonID[iKaon]];
            for (int jKaon = (iKaon+1); jKaon < nKaon; jKaon++)
            {
                const auto Kaon2 = pythia.event[kaonID[jKaon]];
                const auto pPhi = Kaon1.p() + Kaon2.p();
                
                //Cut on Rapidity
                if (abs(pPhi.rap()) >= 0.5) continue;
                
                //Cut on Invariant Mass
                if (pPhi.mCalc() < fMinIM2D) continue;
                if (pPhi.mCalc() > fMaxIM2D) continue;
                if ( Kaon1.id() == Kaon2.id() )
                {
                    evKaonS.InvMass[evKaonS.nKaonCouple]   = pPhi.mCalc();
                    evKaonS.pT[evKaonS.nKaonCouple]        = pPhi.pT();
                    evKaonS.bRec[evKaonS.nKaonCouple]      = (kaonRec[iKaon] && kaonRec[jKaon]);
                    evKaonS.bPhi[evKaonS.nKaonCouple]      = 0;
                    evKaonS.iKaon[evKaonS.nKaonCouple]     = iKaon;
                    evKaonS.jKaon[evKaonS.nKaonCouple]     = jKaon;
                    evKaonS.nKaonCouple++;
                }
                if ( Kaon1.id() != Kaon2.id() )
                {
                    evKaonD.InvMass[evKaonD.nKaonCouple]   = pPhi.mCalc();
                    evKaonD.pT[evKaonD.nKaonCouple]        = pPhi.pT();
                    evKaonD.bRec[evKaonD.nKaonCouple]      = (kaonRec[iKaon] && kaonRec[jKaon]);
                    evKaonD.bPhi[evKaonD.nKaonCouple]      = (Kaon1.mother2() == 0 && Kaon1.mother1() == Kaon2.mother1() && (pythia.event[Kaon1.mother1()]).id() == 333);
                    evKaonD.iKaon[evKaonD.nKaonCouple]     = iKaon;
                    evKaonD.jKaon[evKaonD.nKaonCouple]     = jKaon;
                    evKaonD.nKaonCouple++;
                }
            }
        }
        PtreeKS->Fill();
        PtreeKD->Fill();
    }
    PtreeKS     ->Write();
    PtreeKD     ->Write();
    outFile     ->Close();
    return 0;
}
