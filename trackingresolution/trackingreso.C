#include "../common/plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))
void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs);

Bool_t checkPathFileAndAdd(TString pathtofile,TChain* &t_chaininp){
  if(!gSystem->AccessPathName(pathtofile.Data())){
    t_chaininp->Add(pathtofile.Data());
    return kTRUE;
  }else{
    cout << pathtofile.Data() << " not found!" << endl;
    return kFALSE;
  }
}

void trackingreso(
    TString suffix            = "pdf"
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString dateForOutput                       = ReturnDateStringForOutput();

  TString outputDir 						                  = Form("plots/%s",dateForOutput.Data());
  gSystem->Exec("mkdir -p "+outputDir);

  Double_t textSizeSinglePad                    = 0.05;
  Double_t textSizeLabelsPixel                    = 35;
  Double_t textSizeLabelsRel                    = 58./1300;


  //************************** Read data **************************************************
  const Int_t nEne = 23;
  const static Double_t partE[]   = {3.0, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 7.5, 8.0, 9.0,  10., 11., 13., 15., 20., 25., 30., 40., 50., 60., 70., 80., 150.};
  Int_t useE[nEne]                = { 0};
  Int_t activeE = {0};
  const Int_t nEta = 7;
  const static Double_t etaBins[]   = {1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0};

  TH1F* h_tracks_pT_spectrum[nEne] 	    = {NULL};
  TH1F* h_tracks_truepT_spectrum[nEne] 	    = {NULL};
  TH1F* h_tracks_mean_eta_reso[nEne] 	    = {NULL};
  TH1F* h_tracks_reso_EtaBin_Sigma[nEne] 	    = {NULL};
  TH1F* h_tracks_sigma_eta_reso[nEne] 	    = {NULL};
  TH2F* h_tracks_reso[nEne] 	    = {NULL};
  TH2F* h_tracks_reso_EtaBin[nEta] 	    = {NULL};
  TH2F* h_tracks_reso_Eta[nEne] 	    = {NULL};
  TH2F* h_tracks_nhits_Eta[nEne] 	    = {NULL};
  TH2F* h_tracks_nhits_pT[nEne] 	    = {NULL};
  Double_t etalowdef = 3.0;
  Double_t etahighdef = 4.0;

    for(Int_t etab=0; etab<NELEMS(etaBins); etab++){
        h_tracks_reso_EtaBin[etab] 	= new TH2F(Form("h_tracks_reso_EtaBin_%d",etab), "", 400, 0, 40, 100, -0.3, 0.3);
    }

    for(Int_t epart=0; epart<NELEMS(partE); epart++){
      cout << partE[epart] << endl;

      TChain* t_tracks = new TChain("tracks", "tracks");
      if(partE[epart]==10){
        for(Int_t ii=0; ii<27;ii++){
          checkPathFileAndAdd(Form("/media/nschmidt/external2/EICsimulationOutputs/%s/5000_%d/G4EICDetector_g4tracking_eval.root","Fun4All_G4_FullDetector_DefaultGranHCAL_PythiaJET_beampipe",ii),t_tracks);
          checkPathFileAndAdd(Form("/media/nschmidt/external2/EICsimulationOutputs/%s/5000_%d/G4EICDetector_g4tracking_eval.root","Fun4All_G4_FullDetector_HighGranHCAL_PythiaJET_beampipe",ii),t_tracks);
          checkPathFileAndAdd(Form("/media/nschmidt/external2/EICsimulationOutputs/%s/5000_%d/G4EICDetector_g4tracking_eval.root","Fun4All_G4_FullDetector_MaxHighGranHCAL_PythiaJET_beampipe",ii),t_tracks);
        }
      }
      checkPathFileAndAdd(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/800_%1.1f/G4EICDetector_g4tracking_eval.root","Fun4All_G4_FullDetector_Pythia_beampipe",partE[epart]),t_tracks);
      checkPathFileAndAdd(Form("/media/nschmidt/external2/EICsimulationOutputs/%s/10000_%1.1f/G4EICDetector_g4tracking_eval.root","Fun4All_G4_FullDetector_DefaultGranHCAL_beampipe",partE[epart]),t_tracks);
      for(Int_t ii=1; ii<4;ii++){
          checkPathFileAndAdd(Form("/media/nschmidt/external2/EICsimulationOutputs/%s/800_%d_%1.1f/G4EICDetector_g4tracking_eval.root","Fun4All_G4_FullDetector_Pythia_beampipe",ii,partE[epart]),t_tracks);
      }
      for(Int_t ii=2; ii<3;ii++){
          checkPathFileAndAdd(Form("/media/nschmidt/external2/EICsimulationOutputs/%s/10000_%d_%1.1f/G4EICDetector_g4tracking_eval.root","Fun4All_G4_FullDetector_Pythia_beampipe",ii,partE[epart]),t_tracks);
      }
    // if(partE[epart]!=13) continue;

      if(t_tracks){
        useE[epart] = 1;
        activeE++;
        h_tracks_pT_spectrum[epart] 	= new TH1F(Form("h_tracks_pT_spectrum%d",epart), "", 200, 0, 40);
        h_tracks_truepT_spectrum[epart] 	= new TH1F(Form("h_tracks_truepT_spectrum%d",epart), "", 200, 0, 40);
        h_tracks_reso[epart] 	= new TH2F(Form("h_tracks_reso_%d",epart), "", 200, 0, 40, 200, -1, 1);
        h_tracks_reso_Eta[epart] 	= new TH2F(Form("h_tracks_reso_%d",epart), "", 20, 0, 4, 200, -1, 1);
        h_tracks_nhits_Eta[epart] 	= new TH2F(Form("h_tracks_reso_%d",epart), "", 20, 0, 4, 10, -0.5, 9.5);
        h_tracks_nhits_pT[epart] 	= new TH2F(Form("h_tracks_reso_%d",epart), "", 200, 0, 40, 10, -0.5, 9.5);

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

          // if (((true_eta1 > 1.7 && true_eta1 < 3.0)
          //   || (true_eta2 > 1.7 && true_eta2 < 3.0))  //1.7 < 3
          // && (measured_energy1 > 0.1) // above MIP
          // && (measured_energy1 > 0.1) // above MIP
          // && nTowers1>=2 && nTowers2>=2 // actual cluster with >1 tower
          // ){
            h_tracks_reso[epart]->Fill(TMath::Sqrt(TMath::Power(trk_px_true,2)+TMath::Power(trk_py_true,2)),(TMath::Sqrt(TMath::Power(trk_px,2)+TMath::Power(trk_py,2))-TMath::Sqrt(TMath::Power(trk_px_true,2)+TMath::Power(trk_py_true,2)))/TMath::Sqrt(TMath::Power(trk_px_true,2)+TMath::Power(trk_py_true,2)));
            h_tracks_reso_Eta[epart]->Fill(vec_mom_true.Eta(),(TMath::Sqrt(TMath::Power(trk_px,2)+TMath::Power(trk_py,2))-TMath::Sqrt(TMath::Power(trk_px_true,2)+TMath::Power(trk_py_true,2)))/TMath::Sqrt(TMath::Power(trk_px_true,2)+TMath::Power(trk_py_true,2)));
            for(Int_t etabfill=0;etabfill<(nEta-1);etabfill++){
              if(vec_mom_true.Eta()>etaBins[etabfill] && vec_mom_true.Eta()<etaBins[etabfill+1]){
                // if(partE[epart]==10)
                h_tracks_reso_EtaBin[etabfill]->Fill(TMath::Sqrt(TMath::Power(trk_px_true,2)+TMath::Power(trk_py_true,2)),(TMath::Sqrt(TMath::Power(trk_px,2)+TMath::Power(trk_py,2))-TMath::Sqrt(TMath::Power(trk_px_true,2)+TMath::Power(trk_py_true,2)))/TMath::Sqrt(TMath::Power(trk_px_true,2)+TMath::Power(trk_py_true,2)));
              }
              else continue;
            }
            h_tracks_pT_spectrum[epart]->Fill(TMath::Sqrt(TMath::Power(trk_px,2)+TMath::Power(trk_py,2)));
            h_tracks_truepT_spectrum[epart]->Fill(TMath::Sqrt(TMath::Power(trk_px_true,2)+TMath::Power(trk_py_true,2)));
            h_tracks_nhits_Eta[epart]->Fill(TMath::Sqrt(TMath::Power(trk_px_true,2)+TMath::Power(trk_py_true,2)),trk_hits);
            h_tracks_nhits_pT[epart]->Fill(vec_mom_true.Eta(),trk_hits);
          // }
        }
        h_tracks_reso[epart]->Scale(1 / h_tracks_reso[epart]->GetEntries());
        h_tracks_reso_Eta[epart]->Scale(1 / h_tracks_reso_Eta[epart]->GetEntries());
        h_tracks_nhits_Eta[epart]->Scale(1 / h_tracks_nhits_Eta[epart]->GetEntries());
        h_tracks_nhits_pT[epart]->Scale(1 / h_tracks_nhits_pT[epart]->GetEntries());

        h_tracks_mean_eta_reso[epart] = new TH1F(Form("h_tracks_mean_eta_reso%d",epart),"",20, 0, 4);
        h_tracks_sigma_eta_reso[epart] = new TH1F(Form("h_tracks_sigma_eta_reso%d",epart),"",20, 0, 4);
        for (Int_t i=1; i < h_tracks_reso_Eta[epart]->GetNbinsX(); i++){
          TH1D* projectionYdummy = (TH1D*)h_tracks_reso_Eta[epart]->ProjectionY(Form("projectionYdummy%d%d",i,epart), i,i+1,"e");
          h_tracks_mean_eta_reso[epart]->SetBinContent(i,projectionYdummy->GetMean());
          h_tracks_mean_eta_reso[epart]->SetBinError(i,projectionYdummy->GetMeanError());
          h_tracks_sigma_eta_reso[epart]->SetBinContent(i,projectionYdummy->GetStdDev()); //GetMeanError()
          h_tracks_sigma_eta_reso[epart]->SetBinError(i,projectionYdummy->GetRMSError()); //GetMeanError()
        }
      }
    }
    for(Int_t etab=0; etab<NELEMS(etaBins); etab++){
      h_tracks_reso_EtaBin[etab]->Scale(1 / h_tracks_reso_EtaBin[etab]->GetEntries());

      h_tracks_reso_EtaBin_Sigma[etab] = new TH1F(Form("h_tracks_sigma_eta_reso%d",etab),"",400, 0, 40);
      for (Int_t i=1; i < h_tracks_reso_EtaBin[etab]->GetNbinsX(); i++){
        TH1D* projectionYdummy = (TH1D*)h_tracks_reso_EtaBin[etab]->ProjectionY(Form("projectionYdummy%d%d",i,etab), i,i+1,"e");
        h_tracks_reso_EtaBin_Sigma[etab]->SetBinContent(i,projectionYdummy->GetStdDev()); //GetMeanError()
        h_tracks_reso_EtaBin_Sigma[etab]->SetBinError(i,projectionYdummy->GetRMSError()); //GetMeanError()
      }
    }


  TCanvas* cPNG; // = new TCanvas("cPNG", "", gStyle->GetCanvasDefW() * 3,
	split_canvas(cPNG, "cPNG1", activeE);
  Int_t padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useE[epart]) continue;
    cPNG->cd(padnum+1);

    SetStyleHistoTH2ForGraphs(h_tracks_reso[epart], "#it{p}_{T}^{true}","(#it{p}_{T}^{rec}-#it{p}_{T}^{true})/#it{p}_{T}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);

    h_tracks_reso[epart]->Draw("colz");
    drawLatexAdd(Form("%1.1f GeV #pi^{-} in FST",partE[epart]),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    DrawGammaLines(0, 40, 0., 0., 1, kGray+2, 2);

    padnum++;
  }
  cPNG->Print(Form("%s/Tracking_Resolution.%s", outputDir.Data(), suffix.Data()));

	split_canvas(cPNG, "cPNGx", NELEMS(etaBins)-1);
  padnum =0;
  for(Int_t ietabl=0; ietabl<NELEMS(etaBins)-1; ietabl++){
    if(!useE[ietabl]) continue;
    cPNG->SetLeftMargin(0.1);
    cPNG->SetBottomMargin(0.1);
    cPNG->SetTopMargin(0.01);
    cPNG->SetRightMargin(0.01);
    cPNG->cd(padnum+1);
    (cPNG->cd(padnum+1))->SetLogx(1);
    (cPNG->cd(padnum+1))->SetLogz(1);
    // cPNG->SetLogx(1);
    h_tracks_reso_EtaBin[ietabl]->GetXaxis()->SetRangeUser(0.05,40);
    SetStyleHistoTH2ForGraphs(h_tracks_reso_EtaBin[ietabl], "#it{p}_{T}^{true} (GeV/#it{c})","(#it{p}_{T}^{rec}-#it{p}_{T}^{true})/#it{p}_{T}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);

    // h_tracks_reso_EtaBin[ietabl]->GetXaxis()->SetRangeUser(0.5,4);
    h_tracks_reso_EtaBin[ietabl]->Draw("col");

    DrawGammaSetMarker(h_tracks_reso_EtaBin_Sigma[ietabl], 20, 1, kBlack, kBlack);

    h_tracks_reso_EtaBin_Sigma[ietabl]->Draw("same");
    drawLatexAdd(Form("%1.1f<#eta<%1.1f in FST",etaBins[ietabl],etaBins[ietabl+1]),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("#sigma in black",etaBins[ietabl],etaBins[ietabl+1]),0.90,0.15,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    DrawGammaLines(0.1,40, 0., 0., 1, kGray+2, 2);

    padnum++;
  }
  cPNG->Print(Form("%s/Tracking_Resolution_EtaBins.%s", outputDir.Data(), suffix.Data()));

	split_canvas(cPNG, "cPNG1", activeE);
  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useE[epart]) continue;
    cPNG->SetLeftMargin(0.1);
    cPNG->SetBottomMargin(0.1);
    cPNG->SetTopMargin(0.01);
    cPNG->SetRightMargin(0.01);
    cPNG->cd(padnum+1);
    SetStyleHistoTH2ForGraphs(h_tracks_reso_Eta[epart], "#eta^{true}","(#it{p}_{T}^{rec}-#it{p}_{T}^{true})/#it{p}_{T}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);

    h_tracks_reso_Eta[epart]->GetXaxis()->SetRangeUser(0.5,4);
    h_tracks_reso_Eta[epart]->Draw("col,same");
    drawLatexAdd(Form("%1.1f GeV #pi^{-} in FST",partE[epart]),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    DrawGammaLines(0.5, 4, 0., 0., 1, kGray+2, 2);

    padnum++;
  }
  cPNG->Print(Form("%s/Tracking_Resolution_Eta.%s", outputDir.Data(), suffix.Data()));

	split_canvas(cPNG, "cPNG1", activeE);
  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useE[epart]) continue;
    cPNG->SetLeftMargin(0.1);
    cPNG->SetBottomMargin(0.1);
    cPNG->SetTopMargin(0.01);
    cPNG->SetRightMargin(0.01);
    cPNG->cd(padnum+1);
    SetStyleHistoTH2ForGraphs(h_tracks_nhits_Eta[epart], "#eta^{true}","num. hits", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);

    h_tracks_nhits_Eta[epart]->Draw("colz,same");
    drawLatexAdd(Form("%1.1f GeV #pi^{-} in FST",partE[epart]),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    DrawGammaLines(0, 40, 0., 0., 1, kGray+2, 2);

    padnum++;
  }
  cPNG->Print(Form("%s/Tracking_NHits_Eta.%s", outputDir.Data(), suffix.Data()));

  TH2F * histTrackResoMeanSigDummy                           = new TH2F("histTrackResoMeanSigDummy","histTrackResoMeanSigDummy",200, 0.5, 4,1000,-0.02, 0.39);
  SetStyleHistoTH2ForGraphs(histTrackResoMeanSigDummy, "#eta_{rec}","(#sigma)(#it{p}_{T}^{rec}-#it{p}_{T}^{true})/#it{p}_{T}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histTrackResoMeanSigDummy->GetXaxis()->SetNoExponent();
  histTrackResoMeanSigDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histTrackResoMeanSigDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
	split_canvas(cPNG, "cPNG1", activeE);
  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useE[epart]) continue;
    cPNG->SetLeftMargin(0.1);
    cPNG->SetBottomMargin(0.1);
    cPNG->SetTopMargin(0.01);
    cPNG->SetRightMargin(0.01);
    cPNG->cd(padnum+1);
    histTrackResoMeanSigDummy->DrawCopy();
    TLegend* legendMeanSig  = GetAndSetLegend2(0.17, 0.83-(3*textSizeLabelsRel), 0.6, 0.83,1.0*textSizeLabelsPixel, 1, "", 43, 0.15);

    DrawGammaSetMarker(h_tracks_mean_eta_reso[epart], 47, 2, kBlue+2, kBlue+2);
    DrawGammaSetMarker(h_tracks_sigma_eta_reso[epart], 20, 2, kOrange+2, kOrange+2);

    h_tracks_mean_eta_reso[epart]->Draw("colz,same");
    h_tracks_sigma_eta_reso[epart]->Draw("colz,same");
    drawLatexAdd(Form("%1.1f GeV #pi^{-} in FST",partE[epart]),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    DrawGammaLines(0.5, 4, 0., 0., 1, kGray+2, 2);

    legendMeanSig->AddEntry(h_tracks_mean_eta_reso[epart],Form("mean"),"p");
    legendMeanSig->AddEntry(h_tracks_sigma_eta_reso[epart],Form("sigma"),"p");
    legendMeanSig->Draw();

    padnum++;
  }
  cPNG->Print(Form("%s/Tracking_Resolution_Eta_MeanSig.%s", outputDir.Data(), suffix.Data()));

  TH1F* tracking_effi_pT[nEne] = {NULL};
  histTrackResoMeanSigDummy                           = new TH2F("histTrackResoMeanSigDummy","histTrackResoMeanSigDummy",200, 0, 40,1000,0.01, 1.09);
  SetStyleHistoTH2ForGraphs(histTrackResoMeanSigDummy, "#it{p}_{T} (GeV/#it{c})","#epsilon", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histTrackResoMeanSigDummy->GetXaxis()->SetNoExponent();
  histTrackResoMeanSigDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histTrackResoMeanSigDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
	split_canvas(cPNG, "cPNG1", activeE);
  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useE[epart]) continue;
    cPNG->SetLeftMargin(0.1);
    cPNG->SetBottomMargin(0.1);
    cPNG->SetTopMargin(0.01);
    cPNG->SetRightMargin(0.01);
    cPNG->cd(padnum+1);
    cPNG->SetLogy(1);
    histTrackResoMeanSigDummy->GetXaxis()->SetRangeUser(0,partE[epart]);
    histTrackResoMeanSigDummy->DrawCopy();
    // TLegend* legendMeanSig  = GetAndSetLegend2(0.17, 0.83-(3*textSizeLabelsRel), 0.6, 0.83,1.0*textSizeLabelsPixel, 1, "", 43, 0.15);
    tracking_effi_pT[epart] = (TH1F*)h_tracks_pT_spectrum[epart]->Clone(Form("efficiencyhisttrackspT%d",epart));
    tracking_effi_pT[epart]->Divide(h_tracks_truepT_spectrum[epart]);
    DrawGammaSetMarker(tracking_effi_pT[epart], 47, 2, kBlue+2, kBlue+2);
    // DrawGammaSetMarker(h_tracks_sigma_eta_reso[epart], 20, 2, kOrange+2, kOrange+2);

    tracking_effi_pT[epart]->Draw("same");
    // h_tracks_sigma_eta_reso[epart]->Draw("colz,same");
    drawLatexAdd(Form("%1.1f GeV #pi^{-} in FST",partE[epart]),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    DrawGammaLines(0.5, 4, 0., 0., 1, kGray+2, 2);

    // legendMeanSig->AddEntry(h_tracks_mean_eta_reso[epart],Form("mean"),"p");
    // legendMeanSig->AddEntry(h_tracks_sigma_eta_reso[epart],Form("sigma"),"p");
    // legendMeanSig->Draw();

    padnum++;
  }
  cPNG->Print(Form("%s/TrackingEfficiency.%s", outputDir.Data(), suffix.Data()));


}


void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs)
{
  if(numInputs<2){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 1, gStyle->GetCanvasDefH()*1);
	  // cPNG->Divide(2, 2);
  }else if(numInputs<5){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 2, gStyle->GetCanvasDefH()*2);
	  cPNG->Divide(2, 2);
  }else if(numInputs<7){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 3, gStyle->GetCanvasDefH()*2);
	  cPNG->Divide(3, 2);
  }else if(numInputs<10){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 3, gStyle->GetCanvasDefH()*3);
	  cPNG->Divide(3, 3);
  } else if(numInputs<13){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 4, gStyle->GetCanvasDefH()*3);
	  cPNG->Divide(4, 3);
  } else if(numInputs<17){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 4, gStyle->GetCanvasDefH()*4);
	  cPNG->Divide(4, 4);
  } else if(numInputs<21){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 5, gStyle->GetCanvasDefH()*4);
	  cPNG->Divide(5, 4);
  } else {
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 5, gStyle->GetCanvasDefH()*5);
	  cPNG->Divide(5, 5);
  }

}
