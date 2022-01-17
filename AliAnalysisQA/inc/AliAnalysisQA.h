# ifndef    ALIANALYSISQA_h
# define    ALIANALYSISQA_h
//
#include "../src/AliAnalysisQAStyle.cxx"
#include "../src/AliAnalysisQAEvent.cxx"
#include "../src/AliAnalysisQATrack.cxx"
#include "../src/AliAnalysisQAPID.cxx"
//
//>->->->->->->->->->->->->->->
//
//  Generation of Directories
//
//>->->->->->->->->->->->->->->
void
fGenerateOutputDirectories
()   {
    
    // Generate Output PID QA directory
    gROOT->ProcessLine(Form(".! mkdir -p %s/%s",fOutput_Directory.Data(),fOutput_Directory_PID.Data()));
    
    // Generate Output Track QA directory
    gROOT->ProcessLine(Form(".! mkdir -p %s/%s",fOutput_Directory.Data(),fOutput_Directory_Track.Data()));
    
    // Generate Output Multiplicity QA directory
    gROOT->ProcessLine(Form(".! mkdir -p %s/%s",fOutput_Directory.Data(),fOutput_Directory_Mult.Data()));
    
    // Generate Output Event QA directory
    gROOT->ProcessLine(Form(".! mkdir -p %s/%s",fOutput_Directory.Data(),fOutput_Directory_Mult.Data()));
    
    // Generate Output Event QA directory
    gROOT->ProcessLine(Form(".! mkdir -p %s/%s",fOutput_Directory.Data(),fOutput_Directory_Evnt.Data()));
    
}
//_____________________________________________________________________________
bool fSetInputFilesAndLists
()       {
    fInputFile_DT           =   TFile::Open(Form("%s/%s",fInputData_Directory.Data(),fInputData_FileName.Data()));
    if ( !fInputFile_DT )   { cout << "Data File not found! Aborting... " << endl;  return false;    }
    fInputList_DT           =   (TList*)fInputFile_DT->Get(fInputData_TList.Data());
    if ( !fInputFile_DT )   { cout << "Data Quality Control list not found! Aborting... " << endl;  return false;    }
    fInputFile_MC           =   TFile::Open(Form("%s/%s",fInputMntC_Directory.Data(),fInputMntC_FileName.Data()));
    if ( !fInputFile_DT )   { cout << "MC File not found! Aborting... " << endl;  return false;    }
    fInputList_MC           =   (TList*)fInputFile_MC->Get(fInputMntC_TList.Data());
    if ( !fInputFile_DT )   { cout << "MC Quality Control list not found! Aborting... " << endl;  return false;    }
    return true;
}
//_____________________________________________________________________________
bool
fSetInputHistograms
()  {
    //  Recovering Files and Lists
    if ( !fSetInputFilesAndLists() )    return false;
    //
    hQC_Event_Enumerate_DT =   (TH1D*)fInputList_DT   ->   FindObject(fQC_Event_Enumerate);
    if ( !hQC_Event_Enumerate_DT )   { cout << "Data Enumerate Histogram not found! Aborting... " << endl;  return false;    }
    hQC_Event_Enumerate_MC =   (TH1D*)fInputList_MC   ->   FindObject(fQC_Event_Enumerate);
    if ( !hQC_Event_Enumerate_MC )   { cout << "MC Enumerate Histogram not found! Aborting... " << endl;  return false;    }
    //
    return true;
}
//_____________________________________________________________________________
#endif
