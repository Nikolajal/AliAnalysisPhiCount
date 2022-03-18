#include "../../inc/AliAnalysisPhiPair.h"
#include "./SY_Analysis/SY_AN_SigExtraction.cxx"
#include "./SY_Analysis/SY_AN_PID.cxx"
#include "./SY_Analysis/SY_AN_TRK.cxx"
#include "./SY_Analysis/SY_AN_XPU.cxx"
#include "./SY_Analysis/SY_AN_All.cxx"
// !TODO: Nothing

void
SY_Analysis
 ( TString fOption = "mult", TString kFolder = "_p_p__5TeV", TString fType = "SEX" )    {
    //
    //  Signal Extraction Analysis
    if ( fType.Contains("SEX") ) SY_AN_SigExtraction( fOption, kFolder );
    //
    //  PID Analysis
    if ( fType.Contains("PID") ) SY_AN_PID( fOption, kFolder );
    //
    //  TRK Analysis
    if ( fType.Contains("TRK") ) SY_AN_TRK( fOption, kFolder );
    //
    //  Pile Up Analysis
    if ( fType.Contains("XPU") ) SY_AN_XPU( fOption, kFolder );
    //
    //  Full Analysis
    if ( fType.Contains("ALL") ) SY_AN_All( fOption, kFolder );
}
