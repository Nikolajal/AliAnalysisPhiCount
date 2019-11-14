#include "../inc/SetValues.h"
// Pythia
#include "Pythia8/Pythia.h"

//int main ()
int main (int argc, char *argv[])
{
    // Define some simple data structures
    EVKAONCOUPLE    evKaonS;
    EVKAONCOUPLE    evKaonD;
    EVPHI           evPhi;

    //Initialisation of TTree
    
    //TFile * outFile     = new   TFile   (oFilePreBKG,   "recreate");
    TFile * outFile     = new   TFile   (argv[1],       "recreate");
    TTree * PTreeKSig   = new   TTree   (PTreeKSigName, "A ROOT tree for pythia MC - Kaon++ Couples");
    TTree * PTreeKBkg   = new   TTree   (PTreeKBkgName, "A ROOT tree for pythia MC - Kaon+- Couples");
    TTree * PTreePTru   = new   TTree   (PTreePTruName, "A ROOT tree for pythia MC - Phi");
    
    // Filling Kaon Couple TTree
    PTreeKSig->Branch    ("evKaonCouple.nKaonCouple"  ,&evKaonS.nKaonCouple, "nKaonCouple/I");
    PTreeKSig->Branch    ("evKaonCouple.iKaon"        ,&evKaonS.iKaon,       "iKaon[nKaonCouple]/I");
    PTreeKSig->Branch    ("evKaonCouple.jKaon"        ,&evKaonS.jKaon,       "jKaon[nKaonCouple]/I");
    PTreeKSig->Branch    ("evKaonCouple.bPhi"         ,&evKaonS.bPhi,        "bPhi[nKaonCouple]/O");
    PTreeKSig->Branch    ("evKaonCouple.bRec"         ,&evKaonS.bRec,        "bRec[nKaonCouple]/O");
    PTreeKSig->Branch    ("evKaonCouple.bEta"         ,&evKaonS.bEta,        "bEta[nKaonCouple]/O");
    PTreeKSig->Branch    ("evKaonCouple.InvMass"      ,&evKaonS.InvMass,     "InvMass[nKaonCouple]/F");
    PTreeKSig->Branch    ("evKaonCouple.pT"           ,&evKaonS.pT,          "pT[nKaonCouple]/F");
    
    PTreeKBkg->Branch    ("evKaonCouple.nKaonCouple"  ,&evKaonD.nKaonCouple, "nKaonCouple/I");
    PTreeKBkg->Branch    ("evKaonCouple.iKaon"        ,&evKaonD.iKaon,       "iKaon[nKaonCouple]/I");
    PTreeKBkg->Branch    ("evKaonCouple.jKaon"        ,&evKaonD.jKaon,       "jKaon[nKaonCouple]/I");
    PTreeKBkg->Branch    ("evKaonCouple.bPhi"         ,&evKaonD.bPhi,        "bPhi[nKaonCouple]/O");
    PTreeKBkg->Branch    ("evKaonCouple.bRec"         ,&evKaonD.bRec,        "bRec[nKaonCouple]/O");
    PTreeKBkg->Branch    ("evKaonCouple.bEta"         ,&evKaonD.bEta,        "bEta[nKaonCouple]/O");
    PTreeKBkg->Branch    ("evKaonCouple.InvMass"      ,&evKaonD.InvMass,     "InvMass[nKaonCouple]/F");
    PTreeKBkg->Branch    ("evKaonCouple.pT"           ,&evKaonD.pT,          "pT[nKaonCouple]/F");
    
    // Filling Phi TTree
    PTreePTru->Branch    ("evPhi.nPhi"               ,&evPhi.nPhi,               "nPhi/I");
    PTreePTru->Branch    ("evPhi.bEta"               ,&evPhi.bEta,               "bEta[nPhi]/O");
    PTreePTru->Branch    ("evPhi.bRec"               ,&evPhi.bRec,               "bRec[nPhi]/O");
    PTreePTru->Branch    ("evPhi.pT"                 ,&evPhi.pT,                 "pT[nPhi]/F");
    
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
    int nKaon, nPhi, kaonID[1024], phiID[1024], phiRec[1024];
    bool kaonRec[1024];
    
    for (int iEvent = 0; iEvent < PEvents; iEvent++)
    {
        pythia.next();
        nKaon   = 0;
        nPhi    = 0;
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
        evPhi.nPhi = 1;
        for (int iKaon = 0; iKaon < nKaon; iKaon++)
        {
            const auto Kaon1 = pythia.event[kaonID[iKaon]];
            for (int jKaon = (iKaon+1); jKaon < nKaon; jKaon++)
            {
                const auto Kaon2 = pythia.event[kaonID[jKaon]];
                const auto pPhi = Kaon1.p() + Kaon2.p();
                
                //Cut on Invariant Mass
                if (pPhi.mCalc() < fMinIM2D) continue;
                if (pPhi.mCalc() > fMaxIM2D) continue;
                
                if ( Kaon1.id() == Kaon2.id() )
                {
                    evKaonS.InvMass[evKaonS.nKaonCouple]   = pPhi.mCalc();
                    evKaonS.pT[evKaonS.nKaonCouple]        = pPhi.pT();
                    evKaonS.bRec[evKaonS.nKaonCouple]      = (kaonRec[iKaon] && kaonRec[jKaon]);
                    evKaonS.bPhi[evKaonS.nKaonCouple]      = 0;
                    evKaonS.bEta[evKaonS.nKaonCouple]      = (abs(pPhi.rap()) >= 0.5);
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
                    evKaonD.bEta[evKaonS.nKaonCouple]      = (abs(pPhi.rap()) >= 0.5);
                    evKaonD.iKaon[evKaonD.nKaonCouple]     = iKaon;
                    evKaonD.jKaon[evKaonD.nKaonCouple]     = jKaon;
                    if (evKaonD.bRec[evKaonD.nKaonCouple] && evKaonD.bPhi[evKaonD.nKaonCouple])
                    {
                        phiRec[evPhi.nPhi] = Kaon1.mother1();
                        evPhi.nPhi++;
                    }
                    evKaonD.nKaonCouple++;
                }
            }
        }
        phiRec[0] = evPhi.nPhi;
        evPhi.nPhi = 0;
        for (int iPhi = 0; iPhi < nPhi; iPhi++)
        {
            const auto pPhi                 = pythia.event[phiID[iPhi]];
            evPhi.pT[evPhi.nPhi]            = pPhi.pT();
            evPhi.bEta[evPhi.nPhi]          = (abs((pPhi.p()).rap()) >= 0.5);
            evPhi.bRec[evPhi.nPhi]          = false;
            for ( Int_t iRec = 0; iRec < phiRec[0]; iRec++ )
            {
                if ( phiID[iRec] == phiRec[iRec] )
                {
                    evPhi.bRec[evPhi.nPhi] = true;
                    break;
                }
            }
            evPhi.nPhi++;
        }
        PTreeKSig     ->Fill();
        PTreeKBkg     ->Fill();
        PTreePTru     ->Fill();
    }
    PTreeKSig   ->Write();
    PTreeKBkg   ->Write();
    PTreePTru   ->Write();
    outFile     ->Close();
    return 0;
}
