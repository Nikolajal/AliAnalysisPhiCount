#include "../../inc/AliAnalysisPhiPair.h"
// !TODO: Split

void Check_sub_dataset ( TString fFileNames ... )
{
    //---------------------//
    //  Setting up input   //
    //---------------------//
    
    // >-> OPTIONS
    
    va_list args;
    const size_t fNumberOfFiles = sizeof...(fFileNames);
    va_start(args,);
    for ( Int_t iFile = 0; iFile <= fNumberOfFiles; iFile++ )   {
        if ( fFileNames.IsNull() )
        {
            cout << "[ERROR] Must Specify an input root file" << endl;
            cout << "[INFO] Usage PreProcessing(\"file1.root\",\"file2.root\",...)" << endl;
            return;
        }
    }
    
    std::vector<TString>
    
    return;
}

