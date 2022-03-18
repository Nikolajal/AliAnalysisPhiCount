// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../inc/AliAnalysisPhiPair.h"
#include "GeneralAnalysis.cxx"
#include "RooMsgService.h"

void
SystAnalysis
( TString fOption = "Yield", TString kFolder = "_p_p__7TeV", TString fType = "SEX", Int_t iTot = 0, Int_t iInit = 0 )    {
    //
    //  Generating the binning array--------------------------------------------------------------------------
    fSetAllBins();
    //
    if ( fType.Contains("SEX") ) {
        //  --- Recovering Standard Analysis
        TFile*  insFile_Data_YL     =   new TFile   (Form(kASigExtr_FitCheckRst,(TString("Yield")+kFolder).Data()));
        TFile*  insFile_Effc_YL     =   new TFile   (Form(kAnalysis_MCTruthHist,(TString("Yield")+kFolder).Data()));
        //
        auto        fHEventCount    =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "fQC_Event_Enum_FLL" );
        auto        h1D_Nraw        =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "h1D_Nraw_" );
        auto        h1D_Nrec        =   uLoadHistograms<0,TH1F> ( insFile_Effc_YL,  "h1D_Nrec" );
        auto        h1D_Ngen        =   uLoadHistograms<0,TH1F> ( insFile_Effc_YL,  "h1D_Ngen" );
        auto        h2D_Nraw        =   uLoadHistograms<0,TH2F> ( insFile_Data_YL,  "anSS2D_" );
        auto        h1D_Nrec_2Db    =   uLoadHistograms<0,TH1F> ( insFile_Effc_YL,  "h1D_Nrec_2Db" );
        auto        h1D_Ngen_2Db    =   uLoadHistograms<0,TH1F> ( insFile_Effc_YL,  "h1D_Ngen_2Db" );
        //      TODO: Make a separate function to calculate normalisation
        //  --- Minimum Bias Normalisation
        auto        kN_PU           =   (fHEventCount->GetBinContent(kEventCount::kPU_MB));
        auto        kN_Trg          =   -kN_PU +(fHEventCount->GetBinContent(kEventCount::kTrigger));
        auto        kN_Vtx          =   -kN_PU +(fHEventCount->GetBinContent(kEventCount::kVertex));
        auto        kN_MB           =   -kN_PU +(fHEventCount->GetBinContent(kEventCount::kVertex10));
        Double_t    f1DCorrection   =   (1./kBR)        *(1./kN_MB) *(kTriggerEff/1.)   *(kN_Vtx/kN_Trg);
        Double_t    f2DCorrection   =   (1./(kBR*kBR))  *(1./kN_MB) *(kTriggerEff/1.)   *(kN_Vtx/kN_Trg);
        if ( is_pp_anl && ( kEnergy == 5 ) )    {   f1DCorrection   *=  kTriggerEff15n17pq; f2DCorrection   *=  kTriggerEff15n17pq; }
        if ( is_pp_anl && ( kEnergy == 7 ) )    {   f1DCorrection   *=  kTriggerEff10bcdef; f2DCorrection   *=  kTriggerEff10bcdef; }
        //
                    h1D_Nraw        ->  Scale(1.,"width");
        //
        TH1F*       h1D_Nraw_stat   =   uEfficiencyCorrection1D ( h1D_Nraw, h1D_Nrec, h1D_Ngen, f1DCorrection );
        SetAxis(h1D_Nraw_stat,"PT1D");
        h1D_Nraw_stat   ->  SetName("h1D_Nraw_stat");
        //
        std::vector<TH1F*>  k1D_Variations;
        std::vector<TFile*> k1D_VarInFiles;
        for ( auto kCurrent_Syst : kSyst_SEX_1D_Options ) {
            hName           =   Form("h1D_Nraw_%s",kCurrent_Syst.Data());
            push_to_front( k1D_VarInFiles, new TFile ( Form(kASigExtr_FitChkRstSY,(TString("Yield")+kFolder+TString("/Systematics/")).Data(),(kCurrent_Syst).Data(),"1D") ) );
            auto    kCurrent_Variation  =   uEfficiencyCorrection1D ( (TH1F*)((k1D_VarInFiles.at(0))->Get(hName)), h1D_Nrec, h1D_Ngen, f1DCorrection );
            kCurrent_Variation  ->  Scale(1.,"width");
            kCurrent_Variation  ->  SetName(kCurrent_Syst.Data());
            k1D_Variations.push_back( kCurrent_Variation );
        }
        GeneralAnalysis( h1D_Nraw_stat, k1D_Variations, TString(Form(kAnalysis_Systemt_Dir,(TString("Yield")+kFolder).Data()))+TString("SignalExtraction") );
        //
                    h2D_Nraw        ->  Scale(1.,"width");
        //
        TH2F*       h2D_Nraw_stat   =   uEfficiencyCorrection2D ( h2D_Nraw, h1D_Nrec_2Db, h1D_Ngen_2Db, f2DCorrection );
        SetAxis(h2D_Nraw_stat,"PT2D");
        h2D_Nraw_stat   ->  SetName("h2D_Nraw_stat");
        //
        std::vector<TH2F*>  k2D_Variations;
        std::vector<TFile*> k2D_VarInFiles;
        for ( auto kCurrent_Syst : kSyst_SEX_2D_Options ) {
            hName           =   Form("anSS2D_%s",kCurrent_Syst.Data());
            push_to_front( k2D_VarInFiles, new TFile ( Form(kASigExtr_FitChkRstSY,(TString("Yield")+kFolder+TString("/Systematics/")).Data(),(kCurrent_Syst).Data(),"2D") ) );
            auto    kCurrent_Variation  =   uEfficiencyCorrection2D ( (TH2F*)((k2D_VarInFiles.at(0))->Get(hName)), h1D_Nrec_2Db, h1D_Ngen_2Db, f2DCorrection );
            kCurrent_Variation  ->  Scale(1.,"width");
            kCurrent_Variation  ->  SetName(kCurrent_Syst.Data());
            k2D_Variations.push_back( kCurrent_Variation );
        }
        GeneralAnalysis( h2D_Nraw_stat, k2D_Variations, TString(Form(kAnalysis_Systemt_Dir,(TString("Yield")+kFolder).Data()))+TString("SignalExtraction") );
        GeneralAnalysis( h2D_Nraw_stat, k2D_Variations, TString(Form(kAnalysis_Systemt_Dir,(TString("Yield")+kFolder).Data()))+TString("SignalExtraction") );
    }
}
        //
        //GeneralAnalysis(h1DStandard,f1DVars,TString(Form(kAnalysis_SgExSys_Dir,"yield/Systematics/Standard/")));
        
        /*
        
        //  YIELD Efficiency
        TFile      *insFile_EF_Yield    =   new TFile   ( Form(kAnalysis_MCTruthHist,"yield") ); // /Systematics/Standard/") );
        TFile      *insFile_DT_Yield    =   new TFile   ( Form(kASigExtr_FitCheckRst,"yield") ); // /Systematics/Standard/") );
        //
        TH1F       *hEFF_1D;
        hName                           =   Form("hEFF_1D");
        hEFF_1D                         =   (TH1F*)((insFile_EF_Yield)->Get(hName));
        //
        TH1F       *hEFF_2D_fr_1D;
        hName                           =   Form("hEFF_2D_fr_1D");
        hEFF_2D_fr_1D                   =   (TH1F*)((insFile_EF_Yield)->Get(hName));
        //
        TH1F       *h1DStandard;
        hName                           =   Form("hRAW_1D");
        h1DStandard                     =   (TH1F*)((insFile_DT_Yield)->Get(hName));
        //
        hName                           =   Form("fQC_Event_Enum_FLL");
        auto        kN_Trg          =   ((TH1F*)((insFile_DT_Yield)->Get(hName)))->GetBinContent(kEventCount::kTrigger);
        auto        kN_Vtx          =   ((TH1F*)((insFile_DT_Yield)->Get(hName)))->GetBinContent(kEventCount::kVertex);
        auto        kN_MB           =   ((TH1F*)((insFile_DT_Yield)->Get(hName)))->GetBinContent(kEventCount::kVertex10);
        Double_t    f1DCorrection   =   (1./kBR)        *(1./kN_MB) *(kTriggerEff/1.)   *(kN_Vtx/kN_Trg);
        Double_t    f2DCorrection   =   (1./(kBR*kBR))  *(1./kN_MB) *(kTriggerEff/1.)   *(kN_Vtx/kN_Trg);
        //
        h1DStandard->Scale(f1DCorrection);
        h1DStandard->Divide(h1DStandard,hEFF_1D,1.,1.,"b");
        //
        std::vector<TH1F*>  f1DVars;
        std::vector<TFile*> f1DInFiles;
        for ( auto kCurrent_Syst : kSyst_SEX_1D_Options ) {
            hName           =   Form("hRAW_1D");
            f1DInFiles.insert(f1DInFiles.begin() , new TFile   (Form("%s/ExtractionCheck/%s/1D/FitResults_%s.root",Form(kAnalysis_SgExSys_Dir,"yield"yield/Systematics/Standard/),kCurrent_Syst.Data(),kCurrent_Syst.Data())));
            auto    fTarget =   (TH1F*)((f1DInFiles.at(0))->Get(hName));
            fTarget->Scale(f1DCorrection);
            fTarget->Divide(fTarget,hEFF_1D,1.,1.,"b");
            fTarget->SetName(kCurrent_Syst.Data());
            f1DVars.push_back(fTarget);
        }
        //
        GeneralAnalysis(h1DStandard,f1DVars,TString(Form(kAnalysis_SgExSys_Dir,"yield"yield/Systematics/Standard/)));
        //
        TH2F       *h2DStandard;
        hName                           =   Form("hRAW_2D");
        h2DStandard                     =   (TH2F*)((insFile_DT_Yield)->Get(hName));
        //
        h2DStandard->Scale(f2DCorrection);
        h2DStandard->Divide(h2DStandard,hEFF_2D_fr_1D,1.,1.,"b");
        //
        std::vector<TH2F*>  f2DVars;
        std::vector<TFile*> f2DInFiles;
        auto kStop = 3;
        auto iFil = 0;
        for ( auto kCurrent_Syst : kSyst_SEX_2D_Options ) {
            iFil++;
            if ( iFil > kStop )  continue;
            hName           =   Form("hRAW_2D");
            f2DInFiles.insert(f2DInFiles.begin() , new TFile   (Form("%s/ExtractionCheck/%s/2D/FitResults_%s.root",Form(kAnalysis_SgExSys_Dir,"yield"yield/Systematics/Standard//),kCurrent_Syst.Data(),kCurrent_Syst.Data())));
            auto    fTarget =   (TH2F*)((f2DInFiles.at(0))->Get(hName));
            fTarget->Scale(f2DCorrection);
            fTarget->Divide(fTarget,hEFF_2D_fr_1D,1.,1.,"b");
            fTarget->SetName(kCurrent_Syst.Data());
            f2DVars.push_back(fTarget);
        }
        //
        GeneralAnalysis(h2DStandard,f2DVars,TString(Form(kAnalysis_SgExSys_Dir,"yield"/*yield/Systematics/Standard//)));
        //
        GeneralAnalysis(h1DStandard,f1DVars,h2DStandard,f2DVars,TString(Form(kAnalysis_SgExSys_Dir,"yield"/*yield/Systematics/Standard//)));
        //
        for ( auto CurrentFile : f1DInFiles ) CurrentFile->Close();
        for ( auto CurrentFile : f2DInFiles ) CurrentFile->Close();
        //
        insFile_DT_Yield->Close();
        insFile_EF_Yield->Close();
        //
        //GeneralAnalysis();
        //
        
        
        */
        
        /*
    } else if ( fType.Contains("GTK") ) {
        //
        fSetAllBins();
        //
        //  Standard Yield File
        TFile      *insFile_DT_Yield    =   new TFile   ( Form(kASigExtp_FitCheckRst,"yield") );
        //
        //  Standard Yield Spectrum
        //  -   1D
        TH1F       *h1DStandard;
        hName                           =   Form("hRES_1D_Stat");
        h1DStandard                     =   (TH1F*)((insFile_DT_Yield)->Get(hName));
        //  -   2D
        TH2F       *h2DStandard;
        hName                           =   Form("hRES_2D_Stat");
        h2DStandard                     =   new TH2F    ("h2DStandard","h2DStandard",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
        for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ ) {
            auto fInput = (TH1F*)(insFile_DT_Yield->Get(Form("hRES_2D_Cond1_Stat_%i",iPT2D)));
            for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ ) {
                h2DStandard->SetBinContent  (iPT2D+1,jPT2D+1,fInput->GetBinContent(jPT2D+1));
                h2DStandard->SetBinError    (iPT2D+1,jPT2D+1,fInput->GetBinError(jPT2D+1));
            }
        }
        //
        //  Building Variation Histograms
        //  -   1D
        std::vector<TH1F*>  f1DVariations;
        std::vector<TH2F*>  f2DVariations;
        for ( Int_t iBuild = 0; iBuild < 22; iBuild++ ) {
            auto    f1DTarget1  =   (TH1F*)(h1DStandard->Clone());
            auto    f1DTarget2  =   (TH1F*)(h1DStandard->Clone());
            auto    fScaleVal   =   uRandomGen->Gaus(1.,0.08);
                    f1DTarget1  ->  Scale(1+(1-fScaleVal));
                    f1DTarget1  ->  SetName(Form("%s%i",fType.Data(),2*iBuild));
            f1DVariations.push_back(f1DTarget1);
                    f1DTarget2  ->  Scale(1-(1-fScaleVal));
                    f1DTarget2  ->  SetName(Form("%s%i",fType.Data(),2*iBuild+1));
            f1DVariations.push_back(f1DTarget2);
            auto    f2DTarget1  =   (TH2F*)(h2DStandard->Clone());
            auto    f2DTarget2  =   (TH2F*)(h2DStandard->Clone());
                    f2DTarget1  ->  Scale(1+2*(1-fScaleVal));
                    f2DTarget1  ->  SetName(Form("%s%i",fType.Data(),2*iBuild));
            f2DVariations.push_back(f2DTarget1);
                    f2DTarget2  ->  Scale(1-2*(1-fScaleVal));
                    f2DTarget2  ->  SetName(Form("%s%i",fType.Data(),2*iBuild+1));
            f2DVariations.push_back(f2DTarget2);
        }
        GeneralAnalysis(h1DStandard,f1DVariations,TString(Form("%s/%s/",Form(kAnalysis_Systemt_Dir,"yield"),fType.Data())),true);
        GeneralAnalysis(h2DStandard,f2DVariations,TString(Form("%s/%s/",Form(kAnalysis_Systemt_Dir,"yield"),fType.Data())),true);
        GeneralAnalysis(h1DStandard,f1DVariations,h2DStandard,f2DVariations,TString(Form("%s/%s/",Form(kAnalysis_Systemt_Dir,"yield"),fType.Data())));
        GeneralAnalysis();
    } else if ( fType.Contains("TRK") || fType.Contains("PID") ) {
        //
        //  YIELD Efficiency
        TFile      *insFile_DT_Yield    =   new TFile   ( Form(kASigExtp_FitCheckRst,"yield") );
        //
        TH1F       *h1DStandard;
        hName                           =   Form("hRES_1D_Stat");
        h1DStandard                     =   (TH1F*)((insFile_DT_Yield)->Get(hName));
        //
        fSetAllBins();
        //
        auto iVec = 0;
        std::vector<TH1F*>  f1DVars;
        std::vector<TFile*> f1DInFiles;
        for ( Int_t iTer = iInit; iTer <  iTot; iTer++ )    {
            hName           =   Form("hRES_1D_Stat");
            f1DInFiles.push_back(new TFile  (Form("%s/%s/%s/SignalExtrapolation/FitResults.root",Form(kAnalysis_Systemt_Dir,"yield"),fType.Data(),Form("%s_%i",fType.Data(),iTer+1))));
            auto    fTarget =   (TH1F*)((f1DInFiles.at(iVec))->Get(hName));
            fTarget->SetName(Form("%s%i",fType.Data(),iTer+1));
            f1DVars.push_back(fTarget);
            iVec++;
        }
        //
        GeneralAnalysis(h1DStandard,f1DVars,TString(Form("%s/%s/",Form(kAnalysis_Systemt_Dir,"yield"),fType.Data())));
        //
        TH2F       *h2DStandard;
        hName                           =   Form("hRES_2D_Stat");
        h2DStandard                     =   new TH2F    ("h2DStandard","h2DStandard",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
        for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ ) {
            auto fInput = (TH1F*)(insFile_DT_Yield->Get(Form("hRES_2D_Cond1_Stat_%i",iPT2D)));
            for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ ) {
                h2DStandard->SetBinContent  (iPT2D+1,jPT2D+1,fInput->GetBinContent(jPT2D+1));
                h2DStandard->SetBinError    (iPT2D+1,jPT2D+1,fInput->GetBinError(jPT2D+1));
            }
        }
        //
        iVec = 0;
        std::vector<TH2F*>  f2DVars;
        std::vector<TFile*> f2DInFiles;
        for ( Int_t iTer = iInit; iTer < iTot; iTer++ )    {
            hName           =   Form("hRES_2D_Stat");
            f2DInFiles.push_back(new TFile   (Form("%s/%s/%s/SignalExtrapolation/FitResults.root",Form(kAnalysis_Systemt_Dir,"yield"),fType.Data(),Form("%s_%i",fType.Data(),iTer+1))));
            auto    fTarget =   new TH2F    (Form("%s%i",fType.Data(),iTer+1),Form("%s%i",fType.Data(),iTer),nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
            for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ ) {
                auto fInput = (TH1F*)(f2DInFiles.at(iVec)->Get(Form("hRES_2D_Cond1_Stat_%i",iPT2D)));
                for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ ) {
                    fTarget->SetBinContent  (iPT2D+1,jPT2D+1,fInput->GetBinContent(jPT2D+1));
                    fTarget->SetBinError    (iPT2D+1,jPT2D+1,fInput->GetBinError(jPT2D+1));
                }
            }
            f2DVars.push_back(fTarget);
            iVec++;
        }
        //
        GeneralAnalysis(h2DStandard,f2DVars,TString(Form("%s/%s/",Form(kAnalysis_Systemt_Dir,"yield"),fType.Data())));
        //
        GeneralAnalysis(h1DStandard,f1DVars,h2DStandard,f2DVars,TString(Form("%s/%s/",Form(kAnalysis_Systemt_Dir,"yield"),fType.Data())));
        //
        for ( auto CurrentFile : f1DInFiles ) CurrentFile->Close();
        for ( auto CurrentFile : f2DInFiles ) CurrentFile->Close();
        //
        insFile_DT_Yield->Close();
        //
        GeneralAnalysis();
        //
    } else if (  fType.Contains("FLL")  ) {
        SystAnalysis("GTK");
        SystAnalysis("SEX");
        SystAnalysis("PID",6,0);
        SystAnalysis("TRK",12,0);
    } else {
        cout << "[ERROR] The option " << fType.Data() << " is not recognized" << endl;
        cout << "Possible options:" << endl;
        cout << "(\"TRK\",12,0) Track Quality Cuts default " << endl;
        cout << "(\"PID\",6,0) Track Quality Cuts default " << endl;
        cout << "(\"SEX\") Signal Extraction default " << endl;
        cout << "(\"GTK\") Global Tracking default " << endl;
    }
         */
