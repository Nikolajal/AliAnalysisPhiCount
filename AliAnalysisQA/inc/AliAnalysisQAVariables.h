# ifndef    ALIANALYSISQAVARIABLES_h
# define    ALIANALYSISQAVARIABLES_h
//>->->->->->->->->->->->->->->
//
//  Definition of Directories
//
//>->->->->->->->->->->->->->->
//
//>->   Input Data Directory
//
TString fInputData_Directory    =   "/Volumes/NRUBINI_DATASTASH/Dataset/_Data/p_Pb_5TeV/";
TString fInputData_FileName     =   "LHC16qt_STD.root";
TString fInputData_TList        =   "fQCOutputList";
//
//>->->     Input Data Trees and Histograms
//
//>->->->       Event
//
TString fQC_Event_Enumerate     =   "fQC_Event_Enum_FLL";
TString fQC_Event_Vertex_Cut    =   "fQC_Event_Vertex_Cut";
TString fQC_Event_Enum_Mult     =   "fQC_Event_Enum_V0M";
//
//>->->->       Tracks
//
TString fQC_S_Momentum          =   "fQC_%s_Momentum";
TString fQC_S_TMomentum         =   "fQC_%s_TMomentum";
TString fQC_S_Eta               =   "fQC_%s_Eta";
TString fQC_S_Phi               =   "fQC_%s_Phi";
TString fQC_S_DCAXY_P           =   "fQC_%s_DCAXY_P";
TString fQC_S_DCAXY_PT          =   "fQC_%s_DCAXY_PT";
TString fQC_S_DCAZ_P            =   "fQC_%s_DCAZ_P";
TString fQC_S_DCAZ_PT           =   "fQC_%s_DCAZ_PT";
// !TODO: Below can be easily used for above
TString fQC_S_S_P               =   "fQC_%s_%s_P";
TString fQC_S_S_PT              =   "fQC_%s_%s_PT";
TString fQC_PID_S_S_P           =   "fQC_PID_%s_%s_P";
TString fQC_PID_S_S_PT          =   "fQC_PID_%s_%s_PT";
TString fQC_PID_S_S_SEL_S_P     =   "fQC_PID_%s_%s_SEL_%s_P";
TString fQC_PID_S_S_SEL_S_PT    =   "fQC_PID_%s_%s_SEL_%s_PT";
//
//>->->->       PID
//
TString fQC_PID_SignalDET_P     =   "fQC_PID_Signal%s_%s";
//
//>-> Input Monte Carlo Directory
//
TString fInputMntC_Directory    =   "/Volumes/NRUBINI_DATASTASH/Dataset/_Sim/p_Pb_5TeV/";
TString fInputMntC_FileName     =   "LHC17f2b_STD.root";
TString fInputMntC_TList        =   "fQCOutputList";
//
//>-> Output Data Directory
//
TString fOutput_Directory       =   "Result/QualityAssurance/p_Pb_5TeV/";
TString fOutput_Directory_PID   =   "PIDQA";
TString fOutput_Directory_Track =   "TRKQA";
TString fOutput_Directory_Mult  =   "MLTQA";
TString fOutput_Directory_Evnt  =   "EVTQA";
//
//>->-> Output Data Trees and Histograms
//
//>->->->->->->->->->->->->->->
//
//  Internal Variables
//
//>->->->->->->->->->->->->->->
//
//  Files and TLists
//
TFile  *fInputFile_DT;
TList  *fInputList_DT;
TFile  *fInputFile_MC;
TList  *fInputList_MC;
//
//  Histograms
//
TH1D  *hQC_Event_Enumerate_DT;
TH1D  *hQC_Event_Enumerate_MC;
//
#endif
