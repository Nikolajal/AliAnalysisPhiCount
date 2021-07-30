
// Pythia
#include "Pythia8/Pythia.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TBenchmark.h"

using namespace std;
using namespace ROOT;
using namespace Pythia8;

TBenchmark *fBenchmark      =   new TBenchmark();
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
    printf(fMSG_PrintTimer.Data(),  fPrintEvt,  fSuffix.Data(), 100.*fFraction, fSpeedvsS,  fSuffix.Data(), fElapsedM,  fElapsedS -60.*fElapsedM,  fEta____M,  fEta____S -60.*fEta____M);
    fflush(stdout);
    
    // Resuming timer
    fBenchmark->Start(fTimerName.Data());
}
//
//_____________________________________________________________________________
//
Double_t    fGammaPhiValue  ( Double_t fYieldPhi, Double_t fYieldPhiPhi )  {
    return  2*fYieldPhiPhi/fYieldPhi -fYieldPhi;
}
//
//_____________________________________________________________________________
//
Double_t    fGammaPhiError  ( Double_t fYieldPhi, Double_t fYieldPhiPhi, Double_t fErrorPhi, Double_t fErrorPhiPhi)  {
    auto    fPar1   =   2*fErrorPhiPhi/fYieldPhi;
    auto    fPar2   =   (2*fYieldPhiPhi/(fYieldPhi*fYieldPhi)+1)*fErrorPhi;
    return  fPar1 + fPar2;
}

int main (int argc, char *argv[])
{
    // Check everything is good
    if (argc < 5)
    {
        cout << "ERROR: Insufficient parameters given!" << endl;
        cout << "Please use as: ./MCG_PhiQuickReference [filename] [nevents] [seed] [option]" << endl;
        return -1;
    }
    
    // Definition of number of events
    int   nEvents = atoi(argv[2]);
    
    /*
    int kMultBin    = 500;
    Double_t * fArrMult = new Double_t [kMultBin+1];
    for ( int i = 0; i <= kMultBin; i++ )
    {
        fArrMult[i] =   i;
    }
     */
    
    // Definition of option
    int   fOption = atoi(argv[4]);
     
    // Output File
    TFile * outFile     = new   TFile   (Form("%s.root",argv[1]),   "recreate", "", 101);
    /*
    // Event Count
    TH1D   *hEventCount     =   new TH1D    ("fQC_Event_Enumerate", "fQC_Event_Enumerate",      1,          -0.5,       0.5);
    TH1D   *hEventCountMult =   new TH1D    ("hEventCountMult",     "hEventCountMult",          kMultBin,   fArrMult);
    TH1D   *hPhiCount       =   new TH1D    ("hPhiCount",           "hPhiCount",                10,         -0.5,       9.5);
    TH2D   *hPhiCountMult   =   new TH2D    ("hPhiCountMult",       "hPhiCountMult",            10,         -0.5,       9.5,        kMultBin,    fArrMult);
    TH1D   *hPhiYield       =   new TH1D    ("hPhiYield",           "hPhiYield",                10,         -0.5,       9.5);
    TH2D   *hPhiYieldMult   =   new TH2D    ("hPhiYieldMult",       "hPhiYieldMult",            10,         -0.5,       9.5,        kMultBin,    fArrMult);
    TH1D   *hPhiGamma       =   new TH1D    ("hPhiGamma",           "hPhiGamma",                1,          -0.5,       0.5);
    TH1D   *hPhiGammaMult   =   new TH1D    ("hPhiGammaMult",       "hPhiGammaMult",            kMultBin,   fArrMult);
    */
    
    Float_t *fArrPT2D_Tst = new Float_t[12];
    
    fArrPT2D_Tst[0]     =   0.00; //0.4
    fArrPT2D_Tst[1]     =   0.40; //0.3
    fArrPT2D_Tst[2]     =   0.70; //0.2
    fArrPT2D_Tst[3]     =   0.90; //0.1
    fArrPT2D_Tst[4]     =   1.00; //0.2
    fArrPT2D_Tst[5]     =   1.20; //0.2
    fArrPT2D_Tst[6]     =   1.40; //0.2
    fArrPT2D_Tst[7]     =   1.60; //0.4
    fArrPT2D_Tst[8]     =   2.00; //0.8
    fArrPT2D_Tst[9]     =   2.80; //1.2
    fArrPT2D_Tst[10]    =   4.00; //6.0
    fArrPT2D_Tst[11]    =   10.0;
    
    
    
    TH2F       *hTST_2D;
    TH2D       *hTrue2D;
    //
    //  Defining Efficiency and check utilities
    //
    auto hName       =   Form("hTST_2D");
    auto hTitle      =   Form("hTST_2D");
    hTST_2D     =   new TH2F (hName,hTitle,11,fArrPT2D_Tst,11,fArrPT2D_Tst);
    //
    hName       =   Form("hTrue2D");
    hTitle      =   Form("hTrue2D");
    hTrue2D     =   new TH2D("hTrue2D","hTrue2D",   1200, 0., 12., 1200, 0., 12.);
    //
    // PYTHIA INITIALISATION
    Pythia8::Pythia pythia;
    //
    // Settings
    //
    // >-> Physics
    //
    pythia.readString("SoftQCD:all = on");
    pythia.readString("ParticleDecays:limitTau0 = on");
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 7000");
    pythia.readString(Form("333:mMin = %f",0.75));
    pythia.readString(Form("333:mMax = %f",1.25));
    switch (fOption)
    {
        case 1:
            pythia.readString("ColourReconnection:reconnect = on");
            pythia.readString("ColourReconnection:mode = 1");
            pythia.readString("BeamRemnants:remnantMode = 1");
            cout << "[INFO] Using Colour reconnection Mode 1" << endl;
            cout << "The new more QCD based scheme." << endl;
            break;
        case 2:
            pythia.readString("ColourReconnection:reconnect = on");
            pythia.readString("ColourReconnection:mode = 2");
            cout << "[INFO] Using Colour reconnection Mode 2" << endl;
            cout << "The new gluon-move model." << endl;
            break;
        case 3:
            pythia.readString("ColourReconnection:reconnect = on");
            pythia.readString("ColourReconnection:mode = 3");
            pythia.readString("ColourReconnection:forceResonance = on");
            pythia.readString("PartonLevel:earlyResDec = off");
            cout << "[INFO] Using Colour reconnection Mode 3" << endl;
            cout << "The SK I e^+ e^- CR model." << endl;
            break;
        case 4:
            pythia.readString("ColourReconnection:reconnect = on");
            pythia.readString("ColourReconnection:mode = 4");
            pythia.readString("ColourReconnection:forceResonance = on");
            pythia.readString("PartonLevel:earlyResDec = off");
            cout << "[INFO] Using Colour reconnection Mode 4" << endl;
            cout << "The SK II e^+ e^- CR model." << endl;
            break;
        case 5:
            pythia.readString("StringPT:sigma = 0.335");
            pythia.readString("StringZ:aLund = 0.36");
            pythia.readString("StringZ:aLund = 0.56");
            pythia.readString("StringFlav:probQQtoQ = 0.078");
            pythia.readString("StringFlav:probStoUD = 0.2");
            pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275, 0.0275, 0.0275, 0.0275");
            pythia.readString("MultiPartonInteractions:pT0Ref = 2.12");
            pythia.readString("BeamRemnants:remnantMode = 1");
            pythia.readString("BeamRemnants:saturation = 5");
            pythia.readString("ColourReconnection:mode = 1");
            pythia.readString("ColourReconnection:allowDoubleJunRem = off");
            pythia.readString("ColourReconnection:m0 = 2.9");
            pythia.readString("ColourReconnection:allowJunctions = on");
            pythia.readString("ColourReconnection:junctionCorrection = 1.43");
            pythia.readString("ColourReconnection:timeDilationMode = 0");
            cout << "[INFO] Using Colour reconnection Mode 5" << endl;
            cout << "Mode 0 from https://arxiv.org/pdf/1505.01681.pdf" << endl;
            break;
        case 6:
            pythia.readString("StringPT:sigma = 0.335");
            pythia.readString("StringZ:aLund = 0.36");
            pythia.readString("StringZ:aLund = 0.56");
            pythia.readString("StringFlav:probQQtoQ = 0.078");
            pythia.readString("StringFlav:probStoUD = 0.2");
            pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275, 0.0275, 0.0275, 0.0275");
            pythia.readString("MultiPartonInteractions:pT0Ref = 2.15");
            pythia.readString("BeamRemnants:remnantMode = 1");
            pythia.readString("BeamRemnants:saturation = 5");
            pythia.readString("ColourReconnection:mode = 1");
            pythia.readString("ColourReconnection:allowDoubleJunRem = off");
            pythia.readString("ColourReconnection:m0 = 0.3");
            pythia.readString("ColourReconnection:allowJunctions = on");
            pythia.readString("ColourReconnection:junctionCorrection = 1.20");
            pythia.readString("ColourReconnection:timeDilationMode = 2");
            pythia.readString("ColourReconnection:timeDilationPar = 0.18");
            cout << "[INFO] Using Colour reconnection Mode 6" << endl;
            cout << "Mode 2 from https://arxiv.org/pdf/1505.01681.pdf" << endl;
            break;
        case 7:
            pythia.readString("StringPT:sigma = 0.335");
            pythia.readString("StringZ:aLund = 0.36");
            pythia.readString("StringZ:aLund = 0.56");
            pythia.readString("StringFlav:probQQtoQ = 0.078");
            pythia.readString("StringFlav:probStoUD = 0.2");
            pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275, 0.0275, 0.0275, 0.0275");
            pythia.readString("MultiPartonInteractions:pT0Ref = 2.05");
            pythia.readString("BeamRemnants:remnantMode = 1");
            pythia.readString("BeamRemnants:saturation = 5");
            pythia.readString("ColourReconnection:mode = 1");
            pythia.readString("ColourReconnection:allowDoubleJunRem = off");
            pythia.readString("ColourReconnection:m0 = 0.3");
            pythia.readString("ColourReconnection:allowJunctions = on");
            pythia.readString("ColourReconnection:junctionCorrection = 1.15");
            pythia.readString("ColourReconnection:timeDilationMode = 3");
            pythia.readString("ColourReconnection:timeDilationPar = 0.073");
            cout << "[INFO] Using Colour reconnection Mode 7" << endl;
            cout << "Mode 3 from https://arxiv.org/pdf/1505.01681.pdf" << endl;
            break;
        case 8:
            pythia.readString("MultiPartonInteractions:pT0Ref = 2.15");
            pythia.readString("BeamRemnants:remnantMode = 1");
            pythia.readString("BeamRemnants:saturation = 5");
            pythia.readString("ColourReconnection:mode = 1");
            pythia.readString("ColourReconnection:allowDoubleJunRem = off");
            pythia.readString("ColourReconnection:m0 = 0.3");
            pythia.readString("ColourReconnection:allowJunctions = on");
            pythia.readString("ColourReconnection:junctionCorrection = 1.2");
            pythia.readString("ColourReconnection:timeDilationMode = 2");
            pythia.readString("ColourReconnection:timeDilationPar = 0.18");
            pythia.readString("Ropewalk:RopeHadronization = on");
            pythia.readString("Ropewalk:doShoving = on");
            pythia.readString("Ropewalk:tInit = 1.5");
            pythia.readString("Ropewalk:deltat = 0.05");
            pythia.readString("Ropewalk:tShove = 0.1");
            pythia.readString("Ropewalk:gAmplitude = 0");
            pythia.readString("Ropewalk:doFlavour = on");
            pythia.readString("Ropewalk:r0 = 0.5");
            pythia.readString("Ropewalk:m0 = 0.2");
            pythia.readString("Ropewalk:beta = 0.1");
            pythia.readString("PartonVertex:setVertex = on");
            pythia.readString("PartonVertex:protonRadius = 0.7");
            pythia.readString("PartonVertex:emissionWidth = 0.1");
            cout << "[INFO] Using Colour reconnection Mode 8" << endl;
            cout << "Custom mode 1: QCD based CR w/ ropes" << endl;
            break;
        case 9:
            pythia.readString("MultiPartonInteractions:pT0Ref = 2.15");
            pythia.readString("BeamRemnants:remnantMode = 1");
            pythia.readString("BeamRemnants:saturation = 5");
            pythia.readString("ColourReconnection:mode = 1");
            pythia.readString("ColourReconnection:allowDoubleJunRem = off");
            pythia.readString("ColourReconnection:allowJunctions = on");
            pythia.readString("ColourReconnection:m0 = 0.3");
            pythia.readString("ColourReconnection:junctionCorrection = 1.2");
            pythia.readString("ColourReconnection:timeDilationMode = 2");
            pythia.readString("ColourReconnection:timeDilationPar = 0.18");
            pythia.readString("Ropewalk:RopeHadronization = off");
            cout << "[INFO] Using Colour reconnection Mode 9" << endl;
            cout << "Custom mode 2: QCD based CR w/o ropes" << endl;
            break;
        default:
            pythia.readString("ColourReconnection:reconnect = on");
            pythia.readString("ColourReconnection:mode = 0");
            cout << "[INFO] Using default settings" << endl;
            cout << "To change settings choose from below:" << endl;
            cout << "0. The MPI-based original Pythia 8 scheme." << endl;
            cout << "1. The new more QCD based scheme." << endl;
            cout << "2. The new gluon-move model." << endl;
            cout << "3. The SK I e^+ e^- CR model." << endl;
            cout << "4. The SK II e^+ e^- CR model." << endl;
            cout << "5. Mode 0 from https://arxiv.org/pdf/1505.01681.pdf." << endl;
            cout << "6. Mode 2 from https://arxiv.org/pdf/1505.01681.pdf." << endl;
            cout << "7. Mode 3 from https://arxiv.org/pdf/1505.01681.pdf." << endl;
            cout << "8. Custom mode 1: QCD based CR w/ ropes" << endl;
            cout << "9. Custom mode 2: QCD based CR w/o ropes" << endl;
            break;
    }
    //
    // >-> General Settings
    //
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed = %i",atoi(argv[3])));
    pythia.readString("Print:quiet = on");
    //
    // Pythia event generation initialisation
    //
    pythia.init();
    //
    // Timer Start
    //
    fStartTimer("Production");
    
    for ( int iEvent = 0; iEvent < nEvents; iEvent++ )
    {
        // Next event
        pythia.next();
        fPrintLoopTimer("Production",iEvent,nEvents,10000);
        
        // Set Counters to 0
        Int_t   nPhiMesons      =   0;
        Float_t*fPTArray        =   new Float_t[100];
        
        // Starting cycling through event particles
        for ( int iParticle = 0; iParticle < pythia.event.size() ; iParticle++ )
        {
            // Saving particles
            const auto Current_Particle = pythia.event[iParticle];
            
            // Storing True Phis
            if ( Current_Particle.id() == 333 && (fabs(Current_Particle.p().rap()) <= 0.5) )
            {
                fPTArray[nPhiMesons]    =   Current_Particle.pT();
                nPhiMesons++;
            }
        }
        if ( nPhiMesons < 2 ) continue;
        for ( Int_t iTer = 0; iTer < nPhiMesons; iTer++ )   {
            for ( Int_t jTer = 0; jTer < nPhiMesons; jTer++ )   {
                if ( iTer == jTer ) continue;
                hTrue2D         ->  Fill( fPTArray[iTer], fPTArray[jTer], 0.5 );
                hTST_2D         ->  Fill( fPTArray[iTer], fPTArray[jTer], 0.5 );
            }
        }
    }
    
    /*
    // Cycling through events
    for ( int iEvent = 0; iEvent < nEvents; iEvent++ )
    {
        // Next event
        pythia.next();
        hEventCount->Fill(0);
        fPrintLoopTimer("Production",iEvent,nEvents,10000);
        
        // Set Counters to 0
        Int_t   nPhiMesons      =   0;
        Int_t   nMultiplicity   =   0;
        
        // Starting cycling through event particles
        for ( int iParticle = 0; iParticle < pythia.event.size() ; iParticle++ )
        {
            // Saving particles
            const auto Current_Particle = pythia.event[iParticle];
            
            // Storing True Phis
            if ( Current_Particle.id() == 333 && (fabs(Current_Particle.p().rap()) <= 0.5) )
            {
                nPhiMesons++;
            }
            
            //Skipping non-charged particles
            if ( !Current_Particle.isCharged() )    continue;
            
            //Skipping non-final particles
            if ( !Current_Particle.isFinal() )      continue;
            
            //Skipping non-charged particles
            if ( !(Current_Particle.eta() <= .5) )  continue;
            
            // Mutliplicity evaluation
            nMultiplicity++;
        }
        hPhiCount       ->  Fill(nPhiMesons);
        hPhiCountMult   ->  Fill(nPhiMesons,nMultiplicity);
        hEventCountMult ->  Fill(nMultiplicity);
    }
    
    hPhiCount       ->Write();
    hPhiCountMult   ->Write();
    
    hPhiCount       ->SetName("hPhiCount_norm");
    hPhiCountMult   ->SetName("hPhiCountMult_norm");
    
    hPhiCount       ->Scale(1./nEvents);
    hPhiCountMult   ->Scale(1./nEvents);
    
    for ( Int_t iBin = 1; iBin <= 10; iBin++ )
    {
        Double_t fYieldErr  =   0;
        for ( Int_t jBin = iBin; jBin < 10; jBin++ )
        {
            hPhiYield->Fill(iBin,TMath::Binomial(jBin,iBin)*hPhiCount->GetBinContent(jBin+1));
            fYieldErr += (TMath::Binomial(jBin,iBin)*(hPhiCount->GetBinError(jBin+1)))*(TMath::Binomial(jBin,iBin)*(hPhiCount->GetBinError(jBin+1)));
        }
        hPhiYield->SetBinError(iBin+1,sqrt(fYieldErr));
    }
    
    for ( Int_t kBin = 0; kBin < kMultBin; kBin++ )
    {
        for ( Int_t iBin = 1; iBin <= 10; iBin++ )
        {
            Double_t fYieldErrMult  =   0;
            for ( Int_t jBin = iBin; jBin < 10; jBin++ )
            {
                auto fMult = 0.5*(fArrMult[kBin+1]+fArrMult[kBin]);
                hPhiYieldMult->Fill(iBin,fMult,TMath::Binomial(jBin,iBin)*hPhiCountMult->GetBinContent(jBin+1,kBin+1));
                fYieldErrMult += (TMath::Binomial(jBin,iBin)*hPhiCountMult->GetBinError(jBin+1,kBin+1))*(TMath::Binomial(jBin,iBin)*hPhiCountMult->GetBinError(jBin+1,kBin+1));
            }
            hPhiYieldMult->SetBinError(iBin+1,kBin+1,sqrt(fYieldErrMult));
        }
    }

    auto fYphi      =   hPhiYield->GetBinContent(2);
    auto fYphiphi   =   hPhiYield->GetBinContent(3);
    auto fYEphi     =   hPhiYield->GetBinError(2);
    auto fYEphiphi  =   hPhiYield->GetBinError(3);
    hPhiGamma->SetBinContent(1,fGammaPhiValue(fYphi,fYphiphi));
    hPhiGamma->SetBinError(1,fGammaPhiError(fYphi,fYphiphi,fYEphi,fYEphiphi));
    
    for ( Int_t iBin = 1; iBin <= kMultBin; iBin++ )
    {
        auto fYphi      =   hPhiYieldMult->GetBinContent(2,iBin);
        auto fYphiphi   =   hPhiYieldMult->GetBinContent(3,iBin);
        auto fYEphi     =   hPhiYieldMult->GetBinError(2,iBin);
        auto fYEphiphi  =   hPhiYieldMult->GetBinError(3,iBin);
        
        if ( fYphi == 0 || fYphiphi == 0 ) continue;
        
        hPhiGammaMult->SetBinContent(iBin,fGammaPhiValue(fYphi,fYphiphi));
        hPhiGammaMult->SetBinError(iBin,fGammaPhiError(fYphi,fYphiphi,fYEphi,fYEphiphi));
    }
     */
    
    fStopTimer("Production");
    
    hTrue2D->Write();
    hTST_2D->Write();
    /*
    hEventCount     ->Write();
    hEventCountMult ->Write();
    hPhiCount       ->Write();
    hPhiCountMult   ->Write();
    hPhiYield       ->Write();
    hPhiYieldMult   ->Write();
    hPhiGamma       ->Write();
    hPhiGammaMult   ->Write();
    */
    
    outFile         ->Close();
    return 0;
}
