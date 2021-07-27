#include <TROOT.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TString.h>
#include <TFile.h>

#include <iostream>

// ANCHOR debug output verbosity
Int_t verbosityJRH = 1;

// ANCHOR create histograms globally
// Double_t xbinsDiff[16] = {0.,0.5, 1.2.,4.,6.,8.,10.,12.,14.,16.,18.,  20.,25.,30.,35.,40.,45.};
const int njettypes = 7;
TString jettype[njettypes] = {"track", "full","hcal", "calo","all", "nocluster", "emcal"};

const float min_eta[njettypes] = {-3.5, 0, 1.5, 1.5, 0, 1.5, 1.5};  // TODO Save this info as metadata...
const float max_eta[njettypes] = {3.5, 0, 3.5, 3.5, 0, 3.5, 3.5};

TH2F* h_jet_E_eta[njettypes] = {NULL};
TH2F* h_MCjet_E_eta[njettypes] = {NULL};
TH2F* h_jetscale_E[njettypes][nEta+1] = {NULL};
TH2F* h_jetscale_pT[njettypes][nEta+1] = {NULL};
TH2F* h_jetscale_eta[njettypes][nEta+1] = {NULL};
TH2F* h_jetscale_phi[njettypes][nEta+1] = {NULL};
TH2F* h_constituents_Eta[njettypes] = {NULL};
TH2F* h_constituents_true_Eta[njettypes] = {NULL};

TH1D *h_JES_pT[njettypes][nEta+1] = {NULL};
TH1D *h_JER_pT[njettypes][nEta+1] = {NULL};
TH1D *h_JES_E[njettypes][nEta+1] = {NULL};
TH1D *h_JER_E[njettypes][nEta+1] = {NULL};

TH2F* h2D_jet_EtaReso_E[njettypes][nEta+1] = {NULL};
TH1D *h_EtaReso_Mean_E[njettypes][nEta+1] = {NULL};
TH1D *h_EtaReso_Width_E[njettypes][nEta+1] = {NULL};
TH2F* h2D_jet_EtaReso_Eta[njettypes] = {NULL};
TH1D *h_EtaReso_Mean_Eta[njettypes] = {NULL};
TH1D *h_EtaReso_Width_Eta[njettypes] = {NULL};

// Study the scale and resolution in phi
TH2F* h2D_jet_PhiReso_E[njettypes][nEta+1] = {NULL};
TH1D *h_PhiReso_Mean_E[njettypes][nEta+1] = {NULL};
TH1D *h_PhiReso_Width_E[njettypes][nEta+1] = {NULL};
TH2F* h2D_jet_PhiReso_Eta[njettypes] = {NULL};
TH1D *h_PhiReso_Mean_Eta[njettypes] = {NULL};
TH1D *h_PhiReso_Width_Eta[njettypes] = {NULL};

// Plot the jet spectra in eta, phi, E, pT
TH2F *h2D_truth_reco_eta[njettypes] = {NULL};
TH1F *h_truth_eta[njettypes] = {NULL};
TH1F *h_reco_eta[njettypes] = {NULL};

TH2F *h2D_truth_reco_phi[njettypes] = {NULL};
TH1F *h_truth_phi[njettypes] = {NULL};
TH1F *h_reco_phi[njettypes] = {NULL};

TH2F *h2D_truth_reco_E[njettypes] = {NULL};
TH1F *h_truth_E[njettypes] = {NULL};
TH1F *h_reco_E[njettypes] = {NULL};

TH2F *h2D_truth_reco_pT[njettypes] = {NULL};
TH1F *h_truth_pT[njettypes] = {NULL};
TH1F *h_reco_pT[njettypes] = {NULL};

// For creating jet efficiency graphs
TH1F *h_truth_count[njettypes] = {NULL};
TH1F *h_matched_count[njettypes] = {NULL};

// Comparing eta-phi of truth jets and clusters
TH2F *h2D_truth_eta_phi[njettypes] = {NULL};
TH2F *h2D_cluster_eta_phi[njettypes] = {NULL};


// ANCHOR main function to be called in event loop
void jetresolutionhistos(std::tuple<std::shared_ptr<fastjet::ClusterSequenceArea>, std::vector<fastjet::PseudoJet>> recjets, std::tuple<std::shared_ptr<fastjet::ClusterSequenceArea>, std::vector<fastjet::PseudoJet>> truejets, int select){

  for( int isel=0;isel<njettypes;isel++){
    if(!h_MCjet_E_eta[isel])h_MCjet_E_eta[isel] = new TH2F(Form("h_MCjet_%s_E_eta",jettype[isel].Data()),"",40,0,200,80,-4,4);
    if(!h_jet_E_eta[isel])h_jet_E_eta[isel] = new TH2F(Form("h_jet_%s_E_eta",jettype[isel].Data()),"",40,0,200,80,-4,4);
    if(!h_constituents_Eta[isel])h_constituents_Eta[isel] = new TH2F(Form("h_constituents_Eta_%s",jettype[isel].Data()),"",80,-4,4,50,0,50);
    if(!h_constituents_true_Eta[isel])h_constituents_true_Eta[isel] = new TH2F(Form("h_constituents_true_Eta_%s",jettype[isel].Data()),"",80,-4,4,50,0,50);
    

    if(!h2D_jet_EtaReso_Eta[isel])h2D_jet_EtaReso_Eta[isel]  = new TH2F(Form("h2D_jet_EtaReso_%s_Eta",jettype[isel].Data()),"",40,0,4,200,-1,1);
    if(!h_EtaReso_Mean_Eta[isel])h_EtaReso_Mean_Eta[isel]    = new TH1D(Form("h_EtaReso_Mean_%s_Eta",jettype[isel].Data()),"",40,0,4);
    if(!h_EtaReso_Width_Eta[isel])h_EtaReso_Width_Eta[isel]  = new TH1D(Form("h_EtaReso_Width_%s_Eta",jettype[isel].Data()),"",40,0,4);

    if(!h2D_jet_PhiReso_Eta[isel])h2D_jet_PhiReso_Eta[isel]  = new TH2F(Form("h2D_jet_PhiReso_%s_Eta",jettype[isel].Data()),"",40, 0, 4,200,-1,1);
    if(!h_PhiReso_Mean_Eta[isel])h_PhiReso_Mean_Eta[isel]    = new TH1D(Form("h_PhiReso_Mean_%s_Eta",jettype[isel].Data()),"",40, 0, 4);
    if(!h_PhiReso_Width_Eta[isel])h_PhiReso_Width_Eta[isel]  = new TH1D(Form("h_PhiReso_Width_%s_Eta",jettype[isel].Data()),"",40, 0, 4);

    if (!h2D_truth_reco_eta[isel]) h2D_truth_reco_eta[isel] = new TH2F(Form("h2D_truth_reco_eta_%s", jettype[isel].Data()), "", 40, -7, 7, 40, -7, 7);
    if (!h_truth_eta[isel]) h_truth_eta[isel] = new TH1F(Form("h_truth_eta_%s", jettype[isel].Data()), "", 40, -7, 7);
    if (!h_reco_eta[isel]) h_reco_eta[isel] = new TH1F(Form("h_reco_eta_%s", jettype[isel].Data()), "", 40, -7, 7);
    
    if (!h2D_truth_reco_phi[isel]) h2D_truth_reco_phi[isel] = new TH2F(Form("h2D_truth_reco_phi_%s", jettype[isel].Data()), "", 40, 0, 7, 40, 0, 7);
    if (!h_truth_phi[isel]) h_truth_phi[isel] = new TH1F(Form("h_truth_phi_%s", jettype[isel].Data()), "", 40, 0, 7);
    if (!h_reco_phi[isel]) h_reco_phi[isel] = new TH1F(Form("h_reco_phi_%s", jettype[isel].Data()), "", 40, 0, 7);
    
    if (!h2D_truth_reco_E[isel]) h2D_truth_reco_E[isel] = new TH2F(Form("h2D_truth_reco_E_%s", jettype[isel].Data()), "", 40, 0, 150, 40, 0, 150);
    if (!h_truth_E[isel]) h_truth_E[isel] = new TH1F(Form("h_truth_E_%s", jettype[isel].Data()), "", 40, 0, 150);
    if (!h_reco_E[isel]) h_reco_E[isel] = new TH1F(Form("h_reco_E_%s", jettype[isel].Data()), "", 40, 0, 150);
    
    if (!h2D_truth_reco_pT[isel]) h2D_truth_reco_pT[isel] = new TH2F(Form("h2D_truth_reco_pT_%s", jettype[isel].Data()), "", 40, 0, 120, 40, 0, 120);
    if (!h_truth_pT[isel]) h_truth_pT[isel] = new TH1F(Form("h_truth_pT_%s", jettype[isel].Data()), "", 40, 0, 120);
    if (!h_reco_pT[isel]) h_reco_pT[isel] = new TH1F(Form("h_reco_pT_%s", jettype[isel].Data()), "", 40, 0, 120);


    if (!h_truth_count[isel]) h_truth_count[isel] = new TH1F(Form("h_truth_count_%s", jettype[isel].Data()), "", 40, 0, 200);
    if (!h_matched_count[isel]) h_matched_count[isel] = new TH1F(Form("h_matched_count_%s", jettype[isel].Data()), "", 40, 0, 200);

    if (!h2D_truth_eta_phi[isel]) h2D_truth_eta_phi[isel] = new TH2F(Form("h2D_truth_eta_phi_%s", jettype[isel].Data()), "", 40, -4, 4, 40, 0, 7);
    if (!h2D_cluster_eta_phi[isel]) h2D_cluster_eta_phi[isel] = new TH2F(Form("h2D_cluster_eta_phi_%s", jettype[isel].Data()), "", 40, -4, 4, 40, 0, 7);

    for (Int_t eT = 0; eT < nEta+1; eT++){
      if(!h_jetscale_pT[isel][eT])h_jetscale_pT[isel][eT] = new TH2F(Form("h_jetscale_%s_pT_%d",jettype[isel].Data(), eT),"",40,0,200,200,-1,1);
      if(!h_jetscale_E[isel][eT])h_jetscale_E[isel][eT]  = new TH2F(Form("h_jetscale_%s_E_%d",jettype[isel].Data(), eT),"",40,0,200,200,-1,1);
      if(!h_jetscale_eta[isel][eT])h_jetscale_eta[isel][eT]  = new TH2F(Form("h_jetscale_%s_eta_%d",jettype[isel].Data(), eT),"",40,0,200,200,-1,1);
      if(!h_jetscale_phi[isel][eT])h_jetscale_phi[isel][eT]  = new TH2F(Form("h_jetscale_%s_phi_%d",jettype[isel].Data(), eT),"",40,0,200,200,-1,1);

      if(!h_JES_pT[isel][eT])h_JES_pT[isel][eT]         = new TH1D(Form("h_JES_%s_pT_%d",jettype[isel].Data(), eT),"",40,0,200);
      if(!h_JER_pT[isel][eT])h_JER_pT[isel][eT]         = new TH1D(Form("h_JER_%s_pT_%d",jettype[isel].Data(), eT),"",40,0,200);
      if(!h_JES_E[isel][eT])h_JES_E[isel][eT]           = new TH1D(Form("h_JES_%s_E_%d",jettype[isel].Data(), eT),"",40,0,200);
      if(!h_JER_E[isel][eT])h_JER_E[isel][eT]           = new TH1D(Form("h_JER_%s_E_%d",jettype[isel].Data(), eT),"",40,0,200);

      if(!h2D_jet_EtaReso_E[isel][eT])h2D_jet_EtaReso_E[isel][eT]  = new TH2F(Form("h2D_jet_EtaReso_%s_E_%d",jettype[isel].Data(), eT),"",40,0,200,200,-1,1);
      if(!h_EtaReso_Mean_E[isel][eT])h_EtaReso_Mean_E[isel][eT]    = new TH1D(Form("h_EtaReso_Mean_%s_E_%d",jettype[isel].Data(), eT),"",40,0,200);
      if(!h_EtaReso_Width_E[isel][eT])h_EtaReso_Width_E[isel][eT]  = new TH1D(Form("h_EtaReso_Width_%s_E_%d",jettype[isel].Data(), eT),"",40,0,200);

      if(!h2D_jet_PhiReso_E[isel][eT])h2D_jet_PhiReso_E[isel][eT]  = new TH2F(Form("h2D_jet_PhiReso_%s_E_%d",jettype[isel].Data(), eT),"",40, 0, 200,200,-1,1);
      if(!h_PhiReso_Mean_E[isel][eT])h_PhiReso_Mean_E[isel][eT]    = new TH1D(Form("h_PhiReso_Mean_%s_E_%d",jettype[isel].Data(), eT),"",40, 0, 200);
      if(!h_PhiReso_Width_E[isel][eT])h_PhiReso_Width_E[isel][eT]  = new TH1D(Form("h_PhiReso_Width_%s_E_%d",jettype[isel].Data(), eT),"",40, 0, 200);

    }
  }
  int jets = 0;
  int matched = 0;
    for (std::size_t j = 0; j < std::get<1>(truejets).size(); j++) {
      h2D_truth_eta_phi[select]->Fill(std::get<1>(truejets)[j].eta(), std::get<1>(truejets)[j].phi());
      float eta = std::get<1>(truejets)[j].eta();
      if (eta < min_eta[select] || eta > max_eta[select]) {
        continue;
      }
      h_truth_phi[select]->Fill(std::get<1>(truejets)[j].phi());
      h_truth_eta[select]->Fill(eta);
      h_truth_E[select]->Fill(std::get<1>(truejets)[j].E());
      h_truth_pT[select]->Fill(std::get<1>(truejets)[j].pt());
      h_constituents_true_Eta[select]->Fill(std::get<1>(truejets)[j].eta(),(std::get<1>(truejets)[j].constituents()).size());
     
      jets++;
      h_truth_count[select]->Fill(std::get<1>(truejets)[j].E());    // Jet Efficiency
      for (std::size_t i = 0; i < std::get<1>(recjets).size(); i++) {
        eta = std::get<1>(recjets)[i].eta();
        if (eta < min_eta[select] || eta > max_eta[select]) {
          continue;
        }
        if (std::get<1>(truejets)[j].delta_R(std::get<1>(recjets)[i]) < 0.5) {
          h_matched_count[select]->Fill(std::get<1>(truejets)[j].E());
          matched++;
          break;
        }
      }
    }
    // std::cout << "jets: " << jets << " matched: " << matched << "\n";

  for (std::size_t i = 0; i < std::get<1>(recjets).size(); i++) {
    h2D_cluster_eta_phi[select]->Fill(std::get<1>(recjets)[i].eta(), std::get<1>(recjets)[i].phi());
    float eta = std::get<1>(recjets)[i].eta();
    if (eta < min_eta[select] || eta > max_eta[select]) {
      continue;
    }
    h_reco_phi[select]->Fill(std::get<1>(recjets)[i].phi());
    h_reco_eta[select]->Fill(std::get<1>(recjets)[i].eta());
    h_reco_E[select]->Fill(std::get<1>(recjets)[i].E());
    h_reco_pT[select]->Fill(std::get<1>(recjets)[i].pt());
    // if(verbosityJRH>1)std::cout << "jet " << i << ": "<< std::get<1>(recjets)[i].pt() << " " << std::get<1>(recjets)[i].rap() << " " << std::get<1>(recjets)[i].phi() << "  Ntrue: " << std::get<1>(truejets).size() << endl;
    h_constituents_Eta[select]->Fill(std::get<1>(recjets)[i].eta(),(std::get<1>(recjets)[i].constituents()).size());

    for (std::size_t j = 2; j < std::get<1>(truejets).size(); j++) {
      eta = std::get<1>(truejets)[j].eta();
      if (eta < min_eta[select] || eta > max_eta[select]) {
        continue;
      }
      Double_t deltaRTrueRec = std::get<1>(truejets)[j].delta_R(std::get<1>(recjets)[i]);
      // cout << deltaRTrueRec << endl;
      Int_t et = 0;
      while ( ( std::get<1>(truejets)[j].eta() > partEta[et+1] ) && ( et < nEta )) et++;
      
      if(deltaRTrueRec<0.25){
        // Spectra
        h2D_truth_reco_phi[select]->Fill(std::get<1>(truejets)[j].phi(), std::get<1>(recjets)[i].phi());
        h2D_truth_reco_eta[select]->Fill(std::get<1>(truejets)[j].eta(), std::get<1>(recjets)[i].eta());
        h2D_truth_reco_E[select]->Fill(std::get<1>(truejets)[j].E(), std::get<1>(recjets)[i].E());
        h2D_truth_reco_pT[select]->Fill(std::get<1>(truejets)[j].pt(), std::get<1>(recjets)[i].pt());

        h2D_jet_EtaReso_Eta[select]->Fill(std::get<1>(truejets)[j].eta(),(std::get<1>(recjets)[i].eta()-std::get<1>(truejets)[j].eta()));

        if(verbosityJRH>1) std::cout << "rec jet " << i << " (" << std::get<1>(recjets)[i].E() << "E)  matched with true jet " << j << " (" << std::get<1>(truejets)[j].E() << "E) with dR = " << deltaRTrueRec << std::endl;
        h_jetscale_E[select][et]->Fill(std::get<1>(truejets)[j].E(),(std::get<1>(recjets)[i].E()-std::get<1>(truejets)[j].E())/std::get<1>(truejets)[j].E());
        h_jetscale_pT[select][et]->Fill(std::get<1>(truejets)[j].pt(),(std::get<1>(recjets)[i].pt()-std::get<1>(truejets)[j].pt())/std::get<1>(truejets)[j].pt());
        h_jetscale_eta[select][et]->Fill(std::get<1>(truejets)[j].E(),(std::get<1>(recjets)[i].eta()-std::get<1>(truejets)[j].eta()));
        h_jetscale_phi[select][et]->Fill(std::get<1>(truejets)[j].E(),(std::get<1>(truejets)[j].delta_phi_to(std::get<1>(recjets)[i])));

        h2D_jet_EtaReso_E[select][et]->Fill(std::get<1>(truejets)[j].E(),(std::get<1>(recjets)[i].eta()-std::get<1>(truejets)[j].eta())/std::get<1>(truejets)[j].eta());
        h2D_jet_PhiReso_E[select][et]->Fill(std::get<1>(truejets)[j].E(),(std::get<1>(truejets)[j].delta_phi_to(std::get<1>(recjets)[i]))/std::get<1>(truejets)[j].phi());

        if(std::get<1>(truejets)[j].eta()>1.5){
          h_jetscale_E[select][nEta]->Fill(std::get<1>(truejets)[j].E(),(std::get<1>(recjets)[i].E()-std::get<1>(truejets)[j].E())/std::get<1>(truejets)[j].E());
          h_jetscale_pT[select][nEta]->Fill(std::get<1>(truejets)[j].pt(),(std::get<1>(recjets)[i].pt()-std::get<1>(truejets)[j].pt())/std::get<1>(truejets)[j].pt());

          h2D_jet_EtaReso_E[select][nEta]->Fill(std::get<1>(truejets)[j].E(),(std::get<1>(recjets)[i].eta()-std::get<1>(truejets)[j].eta()));
          h2D_jet_PhiReso_E[select][nEta]->Fill(std::get<1>(truejets)[j].E(),(std::get<1>(truejets)[j].delta_phi_to(std::get<1>(recjets)[i])));
        }
        h_MCjet_E_eta[select]->Fill(std::get<1>(truejets)[j].E(),std::get<1>(truejets)[j].eta());
        h_jet_E_eta[select]->Fill(std::get<1>(recjets)[i].E(),std::get<1>(recjets)[i].eta());
      }
    }
  } // jet loop end

}

// ANCHOR save function after event loop
void jetresolutionhistosSave(){
  for (Int_t eT = 0; eT < nEta+1; eT++){
    for( int isel=0;isel<njettypes;isel++){
      if(h_jetscale_pT[isel][eT]){
        // make projections for JES and JER
        for (Int_t i=1; i < h_jetscale_pT[isel][eT]->GetNbinsX()+1; i++){
          TH1D* projectionYdummy = (TH1D*)h_jetscale_pT[isel][eT]->ProjectionY(Form("projectionYdummy%d%d%d",isel,i,eT), i,i+1,"e");
          // projectionYdummy->GetXaxis()->SetRangeUser(-0.6,0);
          h_JES_pT[isel][eT]->SetBinContent(i,projectionYdummy->GetMean());
          h_JES_pT[isel][eT]->SetBinError(i,projectionYdummy->GetMeanError());
          h_JER_pT[isel][eT]->SetBinContent(i,projectionYdummy->GetStdDev()); //GetMeanError()
          h_JER_pT[isel][eT]->SetBinError(i,projectionYdummy->GetStdDevError()); //GetMeanError()
        }
      }
      if(h_jetscale_E[isel][eT]){
        for (Int_t i=1; i < h_jetscale_E[isel][eT]->GetNbinsX()+1; i++){
          // TH1D* projectionYdummy = (TH1D*)h_jetscale_E[isel][eT]->ProjectionY(Form("projectionYdummy%d%d",isel,i), h_jetscale_E[isel][eT]->GetXaxis()->FindBin(partP[i-1]+0.001), h_jetscale_E[isel][eT]->GetXaxis()->FindBin(partP[i]+0.001),"e");
          TH1D* projectionYdummy = (TH1D*)h_jetscale_E[isel][eT]->ProjectionY(Form("projectionYdummy2%d%d%d",isel,i,eT), i,i+1,"e");
          // projectionYdummy->GetXaxis()->SetRangeUser(-0.6,0);
          h_JES_E[isel][eT]->SetBinContent(i,projectionYdummy->GetMean());
          h_JES_E[isel][eT]->SetBinError(i,projectionYdummy->GetMeanError());
          h_JER_E[isel][eT]->SetBinContent(i,projectionYdummy->GetStdDev()); //GetMeanError()
          h_JER_E[isel][eT]->SetBinError(i,projectionYdummy->GetStdDevError()); //GetMeanError()
        }
      }
      if(h2D_jet_EtaReso_E[isel][eT]){
        for (Int_t i=1; i < h2D_jet_EtaReso_E[isel][eT]->GetNbinsX()+1; i++){
          TH1D* projectionYdummy = (TH1D*)h2D_jet_EtaReso_E[isel][eT]->ProjectionY(Form("projectionYdummy3%d%d%d",isel,i,eT), i,i+1,"e");
          h_EtaReso_Mean_E[isel][eT]->SetBinContent(i,projectionYdummy->GetMean());
          h_EtaReso_Mean_E[isel][eT]->SetBinError(i,projectionYdummy->GetMeanError());
          h_EtaReso_Width_E[isel][eT]->SetBinContent(i,projectionYdummy->GetStdDev()); //GetMeanError()
          h_EtaReso_Width_E[isel][eT]->SetBinError(i,projectionYdummy->GetStdDevError()); //GetMeanError()
        }
      }
      if (h2D_jet_PhiReso_E[isel][eT]) {
        for (Int_t i = 1; i < h2D_jet_PhiReso_E[isel][eT]->GetNbinsX() + 1; i++) {
          TH1D *projectionYdummy = (TH1D*)h2D_jet_PhiReso_E[isel][eT]->ProjectionY(Form("projectionYdummy3%d%d%d",isel,i,eT), i,i+1,"e");
          h_PhiReso_Mean_E[isel][eT]->SetBinContent(i,projectionYdummy->GetMean());
          h_PhiReso_Mean_E[isel][eT]->SetBinError(i,projectionYdummy->GetMeanError());
          h_PhiReso_Width_E[isel][eT]->SetBinContent(i,projectionYdummy->GetStdDev()); //GetMeanError()
          h_PhiReso_Width_E[isel][eT]->SetBinError(i,projectionYdummy->GetStdDevError()); //GetMeanError()
        }
      }
    }
  }
  for( int isel=0;isel<njettypes;isel++){
    if(h2D_jet_EtaReso_Eta[isel]){
      for (Int_t i=1; i < h2D_jet_EtaReso_Eta[isel]->GetNbinsX()+1; i++){
        TH1D* projectionYdummy = (TH1D*)h2D_jet_EtaReso_Eta[isel]->ProjectionY(Form("projectionYdummy4%d%d",isel,i), i,i+1,"e");
        h_EtaReso_Mean_Eta[isel]->SetBinContent(i,projectionYdummy->GetMean());
        h_EtaReso_Mean_Eta[isel]->SetBinError(i,projectionYdummy->GetMeanError());
        h_EtaReso_Width_Eta[isel]->SetBinContent(i,projectionYdummy->GetStdDev()); //GetMeanError()
        h_EtaReso_Width_Eta[isel]->SetBinError(i,projectionYdummy->GetStdDevError()); //GetMeanError()
      }
    }
    if (h2D_jet_PhiReso_Eta[isel]) {
      for (Int_t i = 1; i < h2D_jet_PhiReso_Eta[isel]->GetNbinsX() + 1; i++) {
        TH1D* projectionYdummy = (TH1D*)h2D_jet_PhiReso_Eta[isel]->ProjectionY(Form("projectionYdummy4%d%d",isel,i), i,i+1,"e");
        h_PhiReso_Mean_Eta[isel]->SetBinContent(i,projectionYdummy->GetMean());
        h_PhiReso_Mean_Eta[isel]->SetBinError(i,projectionYdummy->GetMeanError());
        h_PhiReso_Width_Eta[isel]->SetBinContent(i,projectionYdummy->GetStdDev()); //GetMeanError()
        h_PhiReso_Width_Eta[isel]->SetBinError(i,projectionYdummy->GetStdDevError()); //GetMeanError()
      }
    }
  }
  // define output file
  TFile* fileOutput = new TFile(Form("%s/output_JRH.root",outputDir.Data()),"RECREATE");

  // write histograms
  for( int isel=0;isel<njettypes;isel++){
    if(h_constituents_true_Eta[isel])h_constituents_true_Eta[isel]->Write();
    if(h_constituents_Eta[isel])h_constituents_Eta[isel]->Write();
    if(h_MCjet_E_eta[isel])h_MCjet_E_eta[isel]->Write();
    if(h_jet_E_eta[isel])h_jet_E_eta[isel]->Write();
    if(h2D_jet_EtaReso_Eta[isel])h2D_jet_EtaReso_Eta[isel]->Write();
    if(h_EtaReso_Mean_Eta[isel])h_EtaReso_Mean_Eta[isel]->Write();
    if(h_EtaReso_Width_Eta[isel])h_EtaReso_Width_Eta[isel]->Write();
    if(h2D_jet_PhiReso_Eta[isel])h2D_jet_PhiReso_Eta[isel]->Write();
    if(h_PhiReso_Mean_Eta[isel])h_PhiReso_Mean_Eta[isel]->Write();
    if(h_PhiReso_Width_Eta[isel])h_PhiReso_Width_Eta[isel]->Write();
    for (Int_t eT = 0; eT < nEta+1; eT++){
      if(h_jetscale_E[isel][eT])h_jetscale_E[isel][eT]->Write();
      if(h_jetscale_pT[isel][eT])h_jetscale_pT[isel][eT]->Write();
      if(h_jetscale_eta[isel][eT])h_jetscale_eta[isel][eT]->Write();
      if(h_jetscale_phi[isel][eT])h_jetscale_phi[isel][eT]->Write();
      if(h_JES_pT[isel][eT])h_JES_pT[isel][eT]->Write();
      if(h_JES_E[isel][eT])h_JES_E[isel][eT]->Write();
      if(h_JER_pT[isel][eT])h_JER_pT[isel][eT]->Write();
      if(h_JER_E[isel][eT])h_JER_E[isel][eT]->Write();
      if(h2D_jet_EtaReso_E[isel][eT])h2D_jet_EtaReso_E[isel][eT]->Write();
      if(h_EtaReso_Mean_E[isel][eT])h_EtaReso_Mean_E[isel][eT]->Write();
      if(h_EtaReso_Width_E[isel][eT])h_EtaReso_Width_E[isel][eT]->Write();
      if(h2D_jet_PhiReso_E[isel][eT])h2D_jet_PhiReso_E[isel][eT]->Write();
      if(h_PhiReso_Mean_E[isel][eT])h_PhiReso_Mean_E[isel][eT]->Write();
      if(h_PhiReso_Width_E[isel][eT])h_PhiReso_Width_E[isel][eT]->Write();
    }
    if (h2D_truth_reco_eta[isel]) h2D_truth_reco_eta[isel]->Write();
    if (h_truth_eta[isel]) h_truth_eta[isel]->Write();
    if (h_reco_eta[isel]) h_reco_eta[isel]->Write();

    if (h2D_truth_reco_phi[isel]) h2D_truth_reco_phi[isel]->Write();
    if (h_truth_phi[isel]) h_truth_phi[isel]->Write();
    if (h_reco_phi[isel]) h_reco_phi[isel]->Write();

    if (h2D_truth_reco_E[isel]) h2D_truth_reco_E[isel]->Write();
    if (h_truth_E[isel]) h_truth_E[isel]->Write();
    if (h_reco_E[isel]) h_reco_E[isel]->Write();

    if (h2D_truth_reco_pT[isel]) h2D_truth_reco_pT[isel]->Write();   
    if (h_truth_pT[isel]) h_truth_pT[isel]->Write();
    if (h_reco_pT[isel]) h_reco_pT[isel]->Write();

    if (h_truth_count[isel]) h_truth_count[isel]->Write(); 
    if (h_matched_count[isel]) h_matched_count[isel]->Write(); 

    if (h2D_truth_eta_phi[isel]) h2D_truth_eta_phi[isel]->Write();
    if (h2D_cluster_eta_phi[isel]) h2D_cluster_eta_phi[isel]->Write();

    std::cout << "writing..." << std::endl;
  }
  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
