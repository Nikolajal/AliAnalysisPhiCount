
// Pythia
#include "Pythia8/Pythia.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TBenchmark.h"

using namespace std;
using namespace ROOT;
using namespace Pythia8;
TBenchmark *fBenchmark      =   new TBenchmark();

typedef struct
{
    unsigned char nPhi,           iKaon[1024],    jKaon[1024];
    float Multiplicity,   Px[1024],       Py[1024],       Pz[1024],   InvMass[1024];
} Struct_PhiCandidate;

typedef struct
{
    unsigned char nKaon,          Charge[1024];
    char  SigmaTOF[1024], SigmaTPC[1024];
    float Multiplicity,   Px[1024],       Py[1024],       Pz[1024],   InvMass[1024];
} Struct_KaonCandidate;

typedef struct
{
    unsigned char nPhi,           Selection[1024];
    float Multiplicity,   Px[1024],       Py[1024],       Pz[1024],   InvMass[1024];
    bool  fTru,           fGen,           fRec;
} Struct_PhiEfficiency;

typedef struct
{
    unsigned char nKaon,          Charge[1024],   Selection[1024];
    float Multiplicity,   Px[1024],       Py[1024],       Pz[1024],   InvMass[1024];
    bool  ftru;
} Struct_KaonEfficiency;

//-// InvMass range Pythia MC
const   float   fMinIMMC  =   0.75;
const   float   fMaxIMMC  =   1.25;

//-// Tree Names
auto const  fPhiCandidate_Tree      =   "PhiCandidate";
auto const  fPhiCandidateEff_Tree   =   "PhiEfficiency";
auto const  fKaonCandidate_Tree     =   "KaonCandidate";
auto const  fKaonCandidateEff_Tree  =   "KaonEfficiency";

TString     fMSG_PrintTimer =   "[INFO] Event # %4.f %s | %02.0f %% | %2.2f %s events/s | Time: %02.0f:%02.0f | ETA: %02.0f:%02.0f \n";
//
//_____________________________________________________________________________
//
void    fStartTimer             ( TString fTimerName )    {
    fBenchmark->Start(fTimerName.Data());
    printf("[INFO] Starting %s \n", fTimerName.Data());
    fflush(stdout);
}
//
//_____________________________________________________________________________
//
void    fStopTimer              ( TString fTimerName )     {
    fBenchmark->Stop(fTimerName.Data());
    printf("[INFO] Stopping %s \n", fTimerName.Data());
    Float_t fElapsedS   = (float)(fBenchmark->GetRealTime(fTimerName.Data()));
    Float_t fElapsedM   = (Int_t)(fElapsedS/60.);
    printf("[INFO] It took %02.0f:%02.0f \n",   fElapsedM,  fElapsedS - 60.*fElapsedM);
    fflush(stdout);
}
//
//_____________________________________________________________________________
//
void    fPrintLoopTimer         ( TString fTimerName, Int_t iEvent, Int_t nEntries, Int_t iPrintInterval )   {
    if ( iEvent%iPrintInterval != 0 || iEvent == 0 ) return;
    
    // Suffix for events
    TString     fSuffix =   "";
    Int_t       fSfxCor =   iPrintInterval;
    if ( iPrintInterval/1e3 < 1e3         && iPrintInterval/1e3 >= 1 )       {
        fSuffix =   "k";
        fSfxCor =   (int)(iPrintInterval/1e3) + iPrintInterval%(int)1e3;
    }
    if ( iPrintInterval/1e6 < 1e6         && iPrintInterval/1e6 >= 1 )       {
        fSuffix =   "mln";
        fSfxCor =   (int)(iPrintInterval/1e6) + iPrintInterval%(int)1e6;
        cout << fSfxCor << endl;
    }
    if ( iPrintInterval/1e9 < 1e9         && iPrintInterval/1e9 >= 1 )       {
        fSuffix =   "mld";
        fSfxCor =   (int)(iPrintInterval/1e9) + iPrintInterval%(int)1e9;
        cout << fSfxCor << endl;
    }
    
    // Stopping timer
    fBenchmark->Stop(fTimerName.Data());
    
    // Evaluating informations
    Float_t fFraction   =   (float)iEvent/((float)nEntries);
    Float_t fElapsedS   =   (float)(fBenchmark->GetRealTime(fTimerName.Data()));
    Float_t fElapsedM   =   (Int_t)(fElapsedS/60.);
    Float_t fPrintEvt   =   (float)iEvent*(float)fSfxCor/((float)iPrintInterval);
    Float_t fSpeedvsS   =   fPrintEvt/fElapsedS;
    Float_t fEta____S   =   (Int_t)(fElapsedS*((float)nEntries/((float)iEvent) -1));
    Float_t fEta____M   =   (Int_t)(fEta____S/60.);
    
    // Printing
    "[INFO] Event # %4.f %s | %02.0f %% | %1.2f %s events/s | Time: %02.0f:%02.0f | ETA: %02.0f:%02.0f \n";
    printf(fMSG_PrintTimer.Data(),  fPrintEvt,  fSuffix.Data(), 100.*fFraction, fSpeedvsS,  fSuffix.Data(), fElapsedM,  fElapsedS -60.*fElapsedM,  fEta____M,  fEta____S -60.*fEta____M);
    fflush(stdout);
    
    // Resuming timer
    fBenchmark->Start(fTimerName.Data());
}


int main (int argc, char *argv[])
{
    // Check everything is good
    if (argc < 4)
    {
        cout << "ERROR: Insufficient parameters given!" << endl;
        cout << "Please use as: ./GeneratorMC [filename] [nevents] [seed]" << endl;
        return -1;
    }
    
    // Definition of number of events
    int   nEvents = atoi(argv[2]);
    
    // Output File
    TFile * outFile     = new   TFile   (Form("%s.root",argv[1]),   "recreate", "", 101);
    
    //Initialisation of TTree
    TTree * fPhiCandidate,  * fPhiEfficiency,   * fKaonCandidate,   * fKaonEfficiency;
    
    // Define tree data structures
    Struct_PhiCandidate     evPhiCandidate;
    Struct_PhiEfficiency    evPhiEfficiency;
    Struct_KaonCandidate    evKaonCandidate;
    Struct_KaonEfficiency   evKaonEfficiency;
    
    // PhiCandidate Tree Set-Up
    fPhiCandidate = new TTree   (fPhiCandidate_Tree,    "Data Tree for Phi Candidates");
    fPhiCandidate->Branch       ("Multiplicity",    &evPhiCandidate.Multiplicity,   "Multiplicity/F");
    fPhiCandidate->Branch       ("nPhi",            &evPhiCandidate.nPhi,           "nPhi/b");
    fPhiCandidate->Branch       ("Px",              &evPhiCandidate.Px,             "Px[nPhi]/F");
    fPhiCandidate->Branch       ("Py",              &evPhiCandidate.Py,             "Py[nPhi]/F");
    fPhiCandidate->Branch       ("Pz",              &evPhiCandidate.Pz,             "Pz[nPhi]/F");
    fPhiCandidate->Branch       ("InvMass",         &evPhiCandidate.InvMass,        "InvMass[nPhi]/F");
    fPhiCandidate->Branch       ("iKaon",           &evPhiCandidate.iKaon,          "iKaon[nPhi]/b");
    fPhiCandidate->Branch       ("jKaon",           &evPhiCandidate.jKaon,          "jKaon[nPhi]/b");
    
    // KaonCandidate Tree Set-Up
    fKaonCandidate = new TTree  (fKaonCandidate_Tree,   "Data Tree for Kaon Candidates");
    fKaonCandidate->Branch      ("Multiplicity",    &evKaonCandidate.Multiplicity,  "Multiplicity/F");
    fKaonCandidate->Branch      ("nKaon",           &evKaonCandidate.nKaon,         "nKaon/b");
    fKaonCandidate->Branch      ("Px",              &evKaonCandidate.Px,            "Px[nKaon]/F");
    fKaonCandidate->Branch      ("Py",              &evKaonCandidate.Py,            "Py[nKaon]/F");
    fKaonCandidate->Branch      ("Pz",              &evKaonCandidate.Pz,            "Pz[nKaon]/F");
    fKaonCandidate->Branch      ("Charge",          &evKaonCandidate.Charge,        "Charge[nKaon]/B");
    fKaonCandidate->Branch      ("TOFSigma",        &evKaonCandidate.SigmaTOF,      "TOFSigma[nKaon]/B");
    fKaonCandidate->Branch      ("TPCSigma",        &evKaonCandidate.SigmaTPC,      "TPCSigma[nKaon]/B");
    
    // PhiEfficiency Tree Set-Up
    fPhiEfficiency = new TTree  (fPhiCandidateEff_Tree, "MC Tree for Phi Efficiency");
    fPhiEfficiency->Branch      ("Multiplicity",    &evPhiEfficiency.Multiplicity,  "Multiplicity/F");
    fPhiEfficiency->Branch      ("nPhi",            &evPhiEfficiency.nPhi,          "nPhi/b");
    fPhiEfficiency->Branch      ("Px",              &evPhiEfficiency.Px,            "Px[nPhi]/F");
    fPhiEfficiency->Branch      ("Py",              &evPhiEfficiency.Py,            "Py[nPhi]/F");
    fPhiEfficiency->Branch      ("Pz",              &evPhiEfficiency.Pz,            "Pz[nPhi]/F");
    fPhiEfficiency->Branch      ("Selection",       &evPhiEfficiency.Selection,     "Selection[nPhi]/b");
    
    // KaonEfficiency Tree Set-Up
    fKaonEfficiency = new TTree (fKaonCandidateEff_Tree,"MC Tree for Kaon Efficiency");
    fKaonEfficiency->Branch     ("Multiplicity",    &evKaonEfficiency.Multiplicity, "Multiplicity/F");
    fKaonEfficiency->Branch     ("nKaon",           &evKaonEfficiency.nKaon,        "nKaon/b");
    fKaonEfficiency->Branch     ("Px",              &evKaonEfficiency.Px,           "Px[nKaon]/F");
    fKaonEfficiency->Branch     ("Py",              &evKaonEfficiency.Py,           "Py[nKaon]/F");
    fKaonEfficiency->Branch     ("Pz",              &evKaonEfficiency.Pz,           "Pz[nKaon]/F");
    fKaonEfficiency->Branch     ("Charge",          &evKaonEfficiency.Charge,       "Charge[nKaon]/B");
    fKaonEfficiency->Branch     ("Selection",       &evKaonEfficiency.Selection,    "Selection[nKaon]/b");
    
    TList  *fQCOutputList   = new TList();
    fQCOutputList   ->SetOwner(kTRUE);
    TH1D   *hEventCount =   new TH1D    ("fQC_Event_Enumerate","fQC_Event_Enumerate",17, 0.5, 17.5);
    fQCOutputList   ->Add(hEventCount);
    
    // PYTHIA INITIALISATION
    Pythia8::Pythia pythia;
    
    // Settings
    pythia.readString("SoftQCD:nonDiffractive = on");
    pythia.readString("ParticleDecays:limitTau0 = on");
    pythia.readString(Form("333:mMin = %f",fMinIMMC));
    pythia.readString(Form("333:mMax = %f",fMaxIMMC));
    pythia.readString("Random:setSeed = on");
    pythia.readString("Print:quiet = on");
    pythia.readString(Form("Random:seed = %i",atoi(argv[3])));
    pythia.settings.listAll();
    pythia.init();
    
    // Save the ID of kaons here
    //int nKaon, nPhi, nRecMistake, kaonID[1024], phiID[1024], phiRec[1024], RecMistake[1024];
    //bool kaonRec[1024];
    
    // Utility variables
    TLorentzVector  iKaon_p,    jKaon_p,    Phi_p;
    int           fMultiplicityCouter;
    
    fStartTimer("Production");
    
    // Cycling through events
    for ( int iEvent = 0; iEvent < nEvents; iEvent++ )
    {
        // Next event
        pythia.next();
        hEventCount->Fill(0);
        hEventCount->Fill(1);
        hEventCount->Fill(2);
        hEventCount->Fill(3);
        hEventCount->Fill(4);
        hEventCount->Fill(5);
        hEventCount->Fill(6);
        hEventCount->Fill(7);
        hEventCount->Fill(8);
        hEventCount->Fill(9);
        hEventCount->Fill(12);
        hEventCount->Fill(16);
        fPrintLoopTimer("Production",iEvent,nEvents,10000);
        
        // Set counters to zero
        evPhiCandidate.nPhi     =   0;
        evPhiEfficiency.nPhi    =   0;
        evKaonCandidate.nKaon   =   0;
        evKaonEfficiency.nKaon  =   0;
        fMultiplicityCouter     =   0;
        
        // Starting cycling through event particles
        for ( int iParticle = 0; iParticle < pythia.event.size() ; iParticle++ )
        {
            // Saving particles
            const auto Current_Particle = pythia.event[iParticle];
            
            // Storing True Phis
            if ( Current_Particle.id() == 333 )
            {
                evPhiEfficiency.Px[evPhiEfficiency.nPhi]        =   Current_Particle.px();
                evPhiEfficiency.Py[evPhiEfficiency.nPhi]        =   Current_Particle.py();
                evPhiEfficiency.Pz[evPhiEfficiency.nPhi]        =   Current_Particle.pz();
                evPhiEfficiency.Selection[evPhiEfficiency.nPhi] =   0;
                
                auto const Dau1                                 =   ( pythia.event[Current_Particle.daughter1()] );
                auto const Dau2                                 =   ( pythia.event[Current_Particle.daughter2()] );
                
                if  ( ( Current_Particle.daughterList().size() == 2 ) &&
                     ( Dau1.id() == -Dau2.id() ) &&
                     ( abs(Dau1.id()) == 321 ) )
                    evPhiEfficiency.Selection[evPhiEfficiency.nPhi]++;
                
                if  ( evPhiEfficiency.Selection[evPhiEfficiency.nPhi] == 1 &&
                     ( fabs(Dau1.eta()) < 0.8 ) &&
                     ( fabs(Dau2.eta()) < 0.8 ) &&
                     ( Dau1.pT() > 0.15 ) &&
                     ( Dau2.pT() > 0.15 ) )
                    evPhiEfficiency.Selection[evPhiEfficiency.nPhi]++;
                
                /*
                evPhi.ID        [evPhi.nPhi]    =   iParticle;
                evPhi.bEta      [evPhi.nPhi]    =   (fabs(particle.p().rap()) <= 0.5);
                evPhi.pT        [evPhi.nPhi]    =   particle.pT();
                evPhi.Dau1      [evPhi.nPhi]    =   particle.daughter1();
                evPhi.Dau2      [evPhi.nPhi]    =   particle.daughter2();
                auto const Dau1                 =   ( pythia.event[particle.daughter1()] );
                auto const Dau2                 =   ( pythia.event[particle.daughter2()] );
                
                evPhi.bKdc      [evPhi.nPhi]    =   ( ( particle.daughterList().size() == 2 ) &&
                                                ( Dau1.id() == -Dau2.id() ) &&
                                                ( abs(Dau1.id()) == 321 ) );
                
                evPhi.bRec      [evPhi.nPhi]    =   ( evPhi.bKdc[evPhi.nPhi] &&
                                                ( fabs(Dau1.eta()) < 0.8 ) &&
                                                ( fabs(Dau2.eta()) < 0.8 ) &&
                                                ( Dau1.pT() > 0.15 ) &&
                                                ( Dau2.pT() > 0.15 ) );
                */
                
                evPhiEfficiency.nPhi++;
            }
            
            //Skipping non-final particles
            if ( !Current_Particle.isFinal() )      continue;
            
            //Skipping non-charged particles
            if ( !Current_Particle.isCharged() )    continue;
            
            // Multiplicity evaluation
            fMultiplicityCouter++;
            
            // Storing Kaons
            if ( fabs(Current_Particle.id()) == 321 )
            {
                if ( !((fabs(Current_Particle.eta()) < 0.8) && (Current_Particle.pT() > 0.15)) ) continue;
                evKaonEfficiency.Px[evKaonEfficiency.nKaon]         =   Current_Particle.px();
                evKaonEfficiency.Py[evKaonEfficiency.nKaon]         =   Current_Particle.py();
                evKaonEfficiency.Pz[evKaonEfficiency.nKaon]         =   Current_Particle.pz();
                evKaonEfficiency.Charge[evKaonEfficiency.nKaon]     =   Current_Particle.charge();
                evKaonEfficiency.Selection[evKaonEfficiency.nKaon]  =   0;
    
                if ( !((fabs(Current_Particle.eta()) < 0.8) && (Current_Particle.pT() > 0.15)) )
                    evKaonEfficiency.Selection[evKaonEfficiency.nKaon]++;
                
                // VV DA IMPLEMENTARE VV
                /*
                if ( ( evKaon.Mom1[iKaon] == evKaon.Mom1[jKaon] ) &&
                    ( evKaon.Mom2[iKaon] == evKaon.Mom2[jKaon] ) &&
                    ( pythia.event[evKaon.Mom1[iKaon]]).id() == 333 ) &&
                    ( evKaon.Mom2[iKaon] == 0 ) )
                    evKaonEfficiency.Selection[evKaonEfficiency.nKaon]++;
                */
                /*
                evKaon.ID       [evKaon.nKaon]  =   iParticle;
                evKaon.bRec     [evKaon.nKaon]  =   (fabs(particle.eta()) < 0.8) && (particle.pT() > 0.15);
                evKaon.Mom1     [evKaon.nKaon]  =   particle.mother1();
                evKaon.Mom2     [evKaon.nKaon]  =   particle.mother2();
                evKaon.nKaon++;
                 */
                
                evKaonEfficiency.nKaon++;
            }
        }
        
        // Cycling through Kaons found
        for ( int iKaon = 0; iKaon < evKaonEfficiency.nKaon; iKaon++ )
        {
            // Storing first Kaon kinematics and sign
            iKaon_p.SetXYZM(evKaonEfficiency.Px[iKaon],evKaonEfficiency.Py[iKaon],evKaonEfficiency.Pz[iKaon],.493677);
            
            evKaonCandidate.Px[evKaonCandidate.nKaon]         =   evKaonEfficiency.Px[iKaon];
            evKaonCandidate.Py[evKaonCandidate.nKaon]         =   evKaonEfficiency.Py[iKaon];
            evKaonCandidate.Pz[evKaonCandidate.nKaon]         =   evKaonEfficiency.Pz[iKaon];
            evKaonCandidate.Charge[evKaonCandidate.nKaon]     =   evKaonEfficiency.Charge[iKaon];
        
            evKaonCandidate.nKaon++;
            
            for ( int jKaon = (iKaon+1); jKaon < evKaonEfficiency.nKaon; jKaon++ )
            {
                // Storing first Kaon kinematics and sign
                jKaon_p.SetXYZM(evKaonEfficiency.Px[jKaon],evKaonEfficiency.Py[jKaon],evKaonEfficiency.Pz[jKaon],.493677);
                
                // Exclude same sig Kaons
                if ( evKaonEfficiency.Charge[iKaon] == evKaonEfficiency.Charge[jKaon] ) continue;
                
                // Building the candidate Phi
                Phi_p   =   ( iKaon_p + jKaon_p );
                
                // Setting the Tree entry
                evPhiCandidate.Px[evPhiCandidate.nPhi]      =   Phi_p.Px();
                evPhiCandidate.Py[evPhiCandidate.nPhi]      =   Phi_p.Py();
                evPhiCandidate.Pz[evPhiCandidate.nPhi]      =   Phi_p.Pz();
                evPhiCandidate.InvMass[evPhiCandidate.nPhi] =   (Phi_p).Mag();
                evPhiCandidate.iKaon[evPhiCandidate.nPhi]   =   iKaon;
                evPhiCandidate.jKaon[evPhiCandidate.nPhi]   =   jKaon;
                
                evPhiCandidate.nPhi++;
            }
        }
        
        // Assigning Multiplicity
        evPhiCandidate.Multiplicity     =   (1./4.)*fMultiplicityCouter;
        evKaonCandidate.Multiplicity    =   (1./4.)*fMultiplicityCouter;
        evPhiEfficiency.Multiplicity    =   (1./4.)*fMultiplicityCouter;
        evKaonEfficiency.Multiplicity   =   (1./4.)*fMultiplicityCouter;
        
        if ( evPhiCandidate.nPhi <= 0 )  {
            hEventCount->Fill(10);
        }
        else    {
            fPhiCandidate   ->Fill();
            fPhiEfficiency  ->Fill();
            fKaonCandidate  ->Fill();
            fKaonEfficiency ->Fill();
        }
    }
    
    fStopTimer("Production");
    
    fQCOutputList   ->Write();
    fPhiCandidate   ->Write();
    fPhiEfficiency  ->Write();
    fKaonCandidate  ->Write();
    fKaonEfficiency ->Write();
    
    outFile     ->Close();
    return 0;
}
