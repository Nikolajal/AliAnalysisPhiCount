#include "SetValues.h"
// Pythia
#include "Pythia8/Pythia.h"

int main ()
{
    Int_t nPhi=0;
    
    // Define some simple data structures
    EVKAON evKaon;
    EVKAONCOUPLE evKaonCouple;
    
    //Initialisation of TTree
    TTree * Ptree       = new   TTree   ("PythiaTree","A ROOT tree for pythia MC");

    // Filling 1st branch
    /*
    Ptree  ->Branch    ("evKaon.nKaon",&evKaon.nKaon,"nKaon/I");
    Ptree  ->Branch    ("evKaon.particleID",&evKaon.particleID,"particleID[nKaon]/I");
    Ptree  ->Branch    ("evKaon.mother1",&evKaon.mother1,"mother1[nKaon]/I");
    Ptree  ->Branch    ("evKaon.mother2",&evKaon.mother2,"mother2[nKaon]/I");
    Ptree  ->Branch    ("evKaon.motherID",&evKaon.motherID,"motherID[nKaon]/I");
    Ptree  ->Branch    ("evKaon.px",&evKaon.px,"px[nKaon]/F");
    Ptree  ->Branch    ("evKaon.py",&evKaon.py,"py[nKaon]/F");
    Ptree  ->Branch    ("evKaon.pz",&evKaon.pz,"pz[nKaon]/F");
    Ptree  ->Branch    ("evKaon.pT",&evKaon.pT,"pT[nKaon]/F");
    Ptree  ->Branch    ("evKaon.e",&evKaon.e,"e[nKaon]/F");
    */
    // Filling 2nd branch
    Ptree  ->Branch    ("evKaonCouple.nKaonCouple",&evKaonCouple.nKaonCouple,"nKaonCouple/I");
    Ptree  ->Branch    ("evKaonCouple.iKaon",&evKaonCouple.iKaon,"iKaon[nKaonCouple]/I");
    Ptree  ->Branch    ("evKaonCouple.jKaon",&evKaonCouple.jKaon,"jKaon[nKaonCouple]/I");
    Ptree  ->Branch    ("evKaonCouple.bPhi",&evKaonCouple.bPhi,"bPhi[nKaonCouple]/O");
    Ptree  ->Branch    ("evKaonCouple.InvMass",&evKaonCouple.InvMass,"InvMass[nKaonCouple]/F");
    Ptree  ->Branch    ("evKaonCouple.pT",&evKaonCouple.pz,"pT[nKaonCouple]/F");
    /*
    Ptree  ->Branch    ("evKaonCouple.px",&evKaonCouple.px,"px[nKaonCouple]/F");
    Ptree  ->Branch    ("evKaonCouple.py",&evKaonCouple.py,"py[nKaonCouple]/F");
    Ptree  ->Branch    ("evKaonCouple.pz",&evKaonCouple.pz,"pz[nKaonCouple]/F");
    Ptree  ->Branch    ("evKaonCouple.e",&evKaonCouple.e,"e[nKaonCouple]/F");
    */
    
    /* PYTHIA INITIALISATION ( Add string pointer for ... ) */
    Pythia8::Pythia pythia;
    
    //Settings
    pythia.readString("SoftQCD:nonDiffractive = on");
    pythia.readString("ParticleDecays:limitTau0 = on");
    //pythia.readString("333:oneChannel = 1 1. 0 -321 321");
    pythia.init();
    
    // save the ID of kaons here
    int nKaon;
    int kaonID[1024];
    
    for (int iEvent = 0; iEvent < PEvents; iEvent++)
    {
        pythia.next();
        
        nKaon = 0;
        for (int iParticle = 0; iParticle < pythia.event.size() ; iParticle++)
        {
            const auto particle = pythia.event[iParticle];
            
            // Counting Phi
            if (particle.id() == 333) nPhi++;
            
            //Skipping non-final particles
            if(!particle.isFinal()) continue;
            
            //Getting Kaons only
            if (particle.id() == 321 || particle.id() == -321)
            {
                
                kaonID[nKaon] = iParticle;
                nKaon++;
            }
        }
        
        /*
        // now I save all the kaons that I have found in the tree
        evKaon.nKaon = nKaon;
        for (int iKaon = 0; iKaon < nKaon; iKaon++)
        {
            const auto particle = pythia.event[kaonID[iKaon]];
            evKaon.particleID[iKaon]   = particle.id();
            evKaon.mother1[iKaon]      = particle.mother1();
            evKaon.mother2[iKaon]      = particle.mother2();
            evKaon.motherID[iKaon]     = (pythia.event[particle.mother1()]).id();
            evKaon.px[iKaon]           = particle.px();
            evKaon.py[iKaon]           = particle.py();
            evKaon.pz[iKaon]           = particle.pz();
            evKaon.pT[iKaon]           = particle.px() + particle.py();
            evKaon.e[iKaon]            = particle.e();
        }
        */
        // now I create Kaon+- pairs branch
        evKaonCouple.nKaonCouple = 0;
        for (int iKaon = 0; iKaon < nKaon; iKaon++)
        {
            const auto Kaon1 = pythia.event[kaonID[iKaon]];
            
            for (int jKaon = (iKaon+1); jKaon < nKaon; jKaon++)
            {
                const auto Kaon2 = pythia.event[kaonID[jKaon]];
                
                if ( Kaon1.id() != Kaon2.id() )
                {
                    evKaonCouple.iKaon[evKaonCouple.nKaonCouple]     = iKaon;
                    evKaonCouple.jKaon[evKaonCouple.nKaonCouple]     = jKaon;
                    evKaonCouple.bPhi[evKaonCouple.nKaonCouple]      = (Kaon1.mother2() == 0 && Kaon1.mother1() == Kaon2.mother1() && (pythia.event[Kaon1.mother1()]).id() == 333);
                    evKaonCouple.pT[evKaonCouple.nKaonCouple]        = Kaon1.pT() + Kaon2.pT();
                    evKaonCouple.InvMass[evKaonCouple.nKaonCouple]   = Pythia8::m(Kaon1.p(),Kaon2.p());
                    /*
                    evKaonCouple.px[evKaonCouple.nKaonCouple]        = Kaon1.px() + Kaon2.px();
                    evKaonCouple.py[evKaonCouple.nKaonCouple]        = Kaon1.py() + Kaon2.py();
                    evKaonCouple.pz[evKaonCouple.nKaonCouple]        = Kaon1.pz() + Kaon2.pz();
                    evKaonCouple.e[evKaonCouple.nKaonCouple]         = Kaon1.e()  + Kaon2.e();
                    */
                    //Cut on Variant Mass
                    if (evKaonCouple.InvMass[evKaonCouple.nKaonCouple] > maxBound) continue;
                    
                    evKaonCouple.nKaonCouple++;
                }
            }
        }
        
        // Filling Tree
        Ptree->Fill();
    }

    TFile * outFile = new TFile(outMC,"recreate");
    Ptree->Write();
    outFile->Close();
    cout << nPhi << endl;
    return 0;
}
