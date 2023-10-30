#include "../../inc/AliAnalysisPhiPair.h"
#include "./SY_Analysis/SY_AN_SigExtraction.cxx"
#include "./SY_Analysis/SY_AN_PID.cxx"
//#include "./SY_AN_TRK.cxx"
#include "./SY_Analysis/SY_AN_XPU.cxx"
#include "./SY_Analysis/SY_AN_All.cxx"
#include "./SY_Analysis/SystematicsGeneralAnalysis.cxx"
// !TODO: Nothing

void
SY_Analysis
 ( TString fOption = "yield mult ", TString kFolder = "_p_p__5TeV", TString fType = " ALL" )    {
    //
    //  Signal Extraction Analysis
    if ( fType.Contains("SEX") ) SY_AN_SigExtraction( fOption, kFolder );
    //
    //  PID Analysis
    if ( fType.Contains("PID") ) SystematicsGeneralAnalysis( fOption, kFolder, "PID" );
    //
    //  TRK Analysis
    if ( fType.Contains("TRK") ) SystematicsGeneralAnalysis( fOption, kFolder, "TRK" );
    //
    //  Pile Up Analysis
    if ( fType.Contains("XPU") ) SY_AN_XPU( fOption, kFolder );
    //
    //  Full Analysis
    if ( fType.Contains("ALL") ) SY_AN_All( fOption, kFolder );
}
