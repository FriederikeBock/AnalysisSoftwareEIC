#include "../common/plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))
#include <fstream>

void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs);



void analyseTreeForTrackingResol(
                                    TString studyname           = "defaultLBL",
                                    TString singleFileName      = "/home/fbock/eic/Singularity/Fun4All_G4_LBLDetector_Org/10000/G4LBLDetector_g4tracking_eval.root",
                                    TString inputFileNameList   = ""
                                ){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString dateForOutput           = ReturnDateStringForOutput();

  TString outputDir               = "rootfiles/";
  gSystem->Exec("mkdir -p "+outputDir);

  //************************** Read data **************************************************
  Int_t nEne                      = 20;
  Double_t partE[21]              = { 0., 0.5, 1.0, 1.5, 2., 2.5, 3.0, 3.5, 4.0, 4.5, 
                                      5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5,
                                      10.};
  const Int_t nEta                = 7;                                        
  Double_t partEta[8]             = { 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0};

  TH2F* h_tracks_reso_pT[nEta+1]          = {NULL};
  TH2F* h_tracks_reso_p[nEta+1]           = {NULL};
  TH2F* h_tracks_resoEta_pT[nEta+1]       = {NULL};
  TH2F* h_tracks_resoPhi_pT[nEta+1]       = {NULL};
  TH2F* h_tracksTrue_Eta_pT               = new TH2F("h_tracksTrue_Eta_pT", "", 200, 0, 10, 320, 1, 4.);
  TH2F* h_tracksRec_Eta_pT                = new TH2F("h_tracksRec_Eta_pT", "", 200, 0, 10, 320, 1, 4.);
  TH2F* h_tracksTrue_Eta_p                = new TH2F("h_tracksTrue_Eta_p", "", 500, 0, 100, 320, 1, 4.);
  TH2F* h_tracksRec_Eta_p                 = new TH2F("h_tracksRec_Eta_p", "", 500, 0, 100, 320, 1, 4.);
  
  for (Int_t et = 0; et<nEta+1; et++){
    if (et < nEta){
      h_tracks_reso_pT[et]        = new TH2F(Form("h_tracks_reso_pT_%d",et), Form("%1.1f < #eta < %1.1f; #it{p}_{T,MC} (GeV/c); (#it{p}_{T,rec}-#it{p}_{T,MC})/#it{p}_{T,MC}",partEta[et], partEta[et+1]), 200, 0, 10, 1000, -0.5, 0.5);
      h_tracks_reso_p[et]         = new TH2F(Form("h_tracks_reso_p_%d",et), Form("%1.1f < #eta < %1.1f; #it{p}_{MC} (GeV/c); (#it{p}_{rec}-#it{p}_{MC})/#it{p}_{MC}",partEta[et], partEta[et+1]), 500, 0, 100, 1000, -0.5, 0.5);
      h_tracks_resoEta_pT[et]     = new TH2F(Form("h_tracks_reso_Eta_pT_%d",et), Form("%1.1f < #eta < %1.1f; #it{p}_{T,MC} (GeV/c); (#eta_{rec}-#eta_{MC})/#eta_{MC}",partEta[et], partEta[et+1]), 200, 0, 10, 1000, -0.1, 0.1);
      h_tracks_resoPhi_pT[et]     = new TH2F(Form("h_tracks_reso_Phi_pT_%d",et), Form("%1.1f < #eta < %1.1f; #it{p}_{T,MC} (GeV/c); (#phi_{rec}-#phi_{MC})/#phi_{MC}",partEta[et], partEta[et+1]), 200, 0, 10, 1000, -0.5, 0.5);
    } else {
      h_tracks_reso_pT[et]        = new TH2F(Form("h_tracks_reso_pT_%d",et), Form("%1.1f < #eta < %1.1f; #it{p}_{T,MC} (GeV/c); (#it{p}_{T,rec}-#it{p}_{T,MC})/#it{p}_{T,MC}",partEta[0], partEta[et]), 200, 0, 10, 1000, -0.5, 0.5);
      h_tracks_reso_p[et]         = new TH2F(Form("h_tracks_reso_p_%d",et), Form("%1.1f < #eta < %1.1f; #it{p}_{MC} (GeV/c); (#it{p}_{rec}-#it{p}_{MC})/#it{p}_{MC}",partEta[0], partEta[et]), 500, 0, 100, 1000, -0.5, 0.5);
      h_tracks_resoEta_pT[et]     = new TH2F(Form("h_tracks_reso_Eta_pT_%d",et), Form("%1.1f < #eta < %1.1f; #it{p}_{T,MC} (GeV/c); (#eta_{rec}-#eta_{MC})/#eta_{MC}",partEta[0], partEta[et]), 200, 0, 10, 1000, -0.1, 0.1);
      h_tracks_resoPhi_pT[et]     = new TH2F(Form("h_tracks_reso_Phi_pT_%d",et), Form("%1.1f < #eta < %1.1f; #it{p}_{T,MC} (GeV/c); (#phi_{rec}-#phi_{MC})/#phi_{MC}",partEta[0], partEta[et]), 200, 0, 10, 1000, -0.5, 0.5);
      
    }
  }
  
  

  TChain* t_tracks = new TChain("tracks", "tracks");
  if (singleFileName.CompareTo("") != 0){
    t_tracks->Add(singleFileName.Data());
  } else {
    ifstream in(inputFileNameList.Data());
    std::cout<<"File: "<<inputFileNameList.Data()<<std::endl;
    std::cout<<"reading input files "<<std::endl;
    // general number of triggers set
    while(!in.eof()){
      TString tempFilename = "";
      in >> tempFilename;
      if (tempFilename.CompareTo("") != 0 ){
        std::cout << "adding " << tempFilename.Data() << std::endl;
        t_tracks->Add(tempFilename.Data());
      }
    }      
  }

  if(t_tracks){
    Float_t trk_px,trk_py,trk_pz,trk_px_true,trk_py_true,trk_pz_true;
    Int_t trk_hits;
    t_tracks->SetBranchAddress("px", &trk_px);
    t_tracks->SetBranchAddress("py", &trk_py);
    t_tracks->SetBranchAddress("pz", &trk_pz);
    t_tracks->SetBranchAddress("gpx", &trk_px_true);
    t_tracks->SetBranchAddress("gpy", &trk_py_true);
    t_tracks->SetBranchAddress("gpz", &trk_pz_true);
    t_tracks->SetBranchAddress("nhits", &trk_hits);

    Int_t nentries = Int_t(t_tracks->GetEntries());
    Float_t lastEvt = 0;
    for (Int_t i = 0; i < nentries; ++i) {
      if (t_tracks->LoadTree(i) < 0)
        break;

      t_tracks->GetEntry(i);
      TVector3 vec_mom_true(trk_px, trk_py, trk_pz);
      TVector3 vec_mom(trk_px_true, trk_py_true, trk_pz_true);

      Double_t truePt     = TMath::Sqrt(TMath::Power(trk_px_true,2)+TMath::Power(trk_py_true,2));
      Double_t trueP      = TMath::Sqrt(TMath::Power(trk_px_true,2)+TMath::Power(trk_py_true,2)+TMath::Power(trk_pz_true,2));
      Double_t trueEta    = vec_mom_true.Eta();
      Double_t truePhi    = vec_mom_true.Phi();
      Double_t recPt      = TMath::Sqrt(TMath::Power(trk_px,2)+TMath::Power(trk_py,2));
      Double_t recP       = TMath::Sqrt(TMath::Power(trk_px,2)+TMath::Power(trk_py,2)+TMath::Power(trk_pz,2));
      Double_t recEta     = vec_mom.Eta();
      Double_t recPhi     = vec_mom.Phi();
      
      Int_t et = 0;
      while (partEta[et+1] < trueEta && et < nEta) et++;
      h_tracks_reso_pT[et]->Fill(truePt,(recPt-truePt)/truePt);
      h_tracks_reso_p[et]->Fill(trueP,(recP-trueP)/trueP);
      h_tracks_resoEta_pT[et]->Fill(truePt,(recEta-trueEta)/trueEta);
      h_tracks_resoPhi_pT[et]->Fill(truePt,(recPhi-truePhi)/truePhi);
      if (trueEta < partEta[nEta] && trueEta > partEta[0]){
        h_tracks_reso_pT[nEta]->Fill(truePt,(recPt-truePt)/truePt);
        h_tracks_reso_p[nEta]->Fill(trueP,(recP-trueP)/trueP);
        h_tracks_resoEta_pT[nEta]->Fill(truePt,(recEta-trueEta)/trueEta);
        h_tracks_resoPhi_pT[nEta]->Fill(truePt,(recPhi-truePhi)/truePhi);
      }
      h_tracksTrue_Eta_pT->Fill(truePt, trueEta);
      h_tracksRec_Eta_pT->Fill(recPt, recEta);
      h_tracksTrue_Eta_p->Fill(trueP, trueEta);
      h_tracksRec_Eta_p->Fill(recP, recEta);
    }
  }

  TString nameOutput      = Form("%s/histosTrackingResol_%s.root",outputDir.Data(), studyname.Data());
  TFile* fOutput          = new TFile(nameOutput,"RECREATE");

    for (Int_t et = 0; et < nEta+1; et++){
      h_tracks_reso_pT[et]->Write();
      h_tracks_reso_p[et]->Write();
      h_tracks_resoEta_pT[et]->Write();
      h_tracks_resoPhi_pT[et]->Write();
    }
    h_tracksTrue_Eta_pT->Write();
    h_tracksRec_Eta_pT->Write();
    h_tracksTrue_Eta_p->Write();
    h_tracksRec_Eta_p->Write();
  fOutput->Write();
  fOutput->Close();
  delete fOutput;
}

