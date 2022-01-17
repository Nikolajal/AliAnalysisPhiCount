#include "../inc/AliAnalysisQAEvent.h"
//
//>->->->->->->->->->->->->->->->
//                              >
//  PID Checks on Selections    >
//                              >
//>->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
//
void    fCheckVertexPosition (  )  {
    // Speeding the process and avoid TCanvas popping up
    gROOT->SetBatch(kTRUE);
    //
    // Making the Vertex Directory
    gROOT->ProcessLine(Form(".! mkdir -p %s/%s/%s",fOutput_Directory.Data(),fOutput_Directory_Evnt.Data(),"VTX"));
    //
    // Recover the Vertex cut Histogram
    TH1F   *fQC_Event_Vertex_Cut_DT =   (TH1F*)fInputList_DT->FindObject(fQC_Event_Vertex_Cut);
    TH1F   *fQC_Event_Vertex_Cut_MC =   (TH1F*)fInputList_MC->FindObject(fQC_Event_Vertex_Cut);
    //
    //>->   Set Normalisation
    //
    auto    nEntries =  fQC_Event_Vertex_Cut_DT->GetEntries();
    fQC_Event_Vertex_Cut_DT->Scale(1./nEntries);
    fQC_Event_Vertex_Cut_DT->SetMinimum(1.e-4);
            nEntries =  fQC_Event_Vertex_Cut_MC->GetEntries();
    fQC_Event_Vertex_Cut_MC->Scale(1./nEntries);
    //
    //>->   Set Graphic Attributes
    //
    fQC_Event_Vertex_Cut_DT->SetTitle("");
    fSetVertexAxes(fQC_Event_Vertex_Cut_DT);
    fSetVertexAxes(fQC_Event_Vertex_Cut_MC);
    fSetColorScheme(fQC_Event_Vertex_Cut_DT,true);
    fSetColorScheme(fQC_Event_Vertex_Cut_MC,false);
    //
    //>->     Legend
    //
    TLegend    *fLegend             =   new TLegend(0.75,0.25,0.9,0.1);
    fLegend                         ->  AddEntry(fQC_Event_Vertex_Cut_DT,"Data","P");
    fLegend                         ->  AddEntry(fQC_Event_Vertex_Cut_MC,"MC","P");
    //
    //>->->     Canvas to save
    //
    TCanvas    *cCheck0             =   new TCanvas();
    gStyle                          ->  SetOptStat(0);
    gPad                            ->  SetLogy();
    fQC_Event_Vertex_Cut_DT         ->  Draw("P SAME");
    fQC_Event_Vertex_Cut_MC         ->  Draw("P SAME");
    fLegend                         ->  Draw("SAME");
    cCheck0                         ->  SaveAs(Form("./%s/%s/%s/%s.pdf",fOutput_Directory.Data(),fOutput_Directory_Evnt.Data(),"VTX","Vertex"));
    
    // Exiting Batch Mode
    gROOT->SetBatch(kFALSE);
}
//
void    fSetVertexAxes          ( TH1F* fTarget )   {
    fTarget->GetXaxis()->SetTitle("Vertex Z position (cm)");
    fTarget->GetYaxis()->SetTitle("Fraction of accepted events");
}
