
// ANCHOR debug output verbosity
Int_t verbosityJRH = 1;

// ANCHOR create histograms globally
// Double_t xbinsDiff[16] = {0.,0.5, 1.2.,4.,6.,8.,10.,12.,14.,16.,18.,  20.,25.,30.,35.,40.,45.};
const int njettypes = 5;
TString jettype[njettypes] = {"track", "full","hcal","calo","all"};

TH2F* h_jet_E_eta[njettypes] = {NULL};
TH2F* h_MCjet_E_eta[njettypes] = {NULL};
TH2F* h_jetscale_E[njettypes][nEta+1] = {NULL};
TH2F* h_jetscale_pT[njettypes][nEta+1] = {NULL};

TH1D *h_JES_pT[njettypes][nEta+1] = {NULL};
TH1D *h_JER_pT[njettypes][nEta+1] = {NULL};
TH1D *h_JES_E[njettypes][nEta+1] = {NULL};
TH1D *h_JER_E[njettypes][nEta+1] = {NULL};

// ANCHOR main function to be called in event loop
void jetresolutionhistos(std::tuple<std::shared_ptr<fastjet::ClusterSequenceArea>, std::vector<fastjet::PseudoJet>> recjets, std::tuple<std::shared_ptr<fastjet::ClusterSequenceArea>, std::vector<fastjet::PseudoJet>> truejets, int select){

  for( int isel=0;isel<njettypes;isel++){
    if(!h_MCjet_E_eta[isel])h_MCjet_E_eta[isel] = new TH2F(Form("h_MCjet_%s_E_eta",jettype[isel].Data()),"",nBinsP, binningP,80,-4,4);
    if(!h_jet_E_eta[isel])h_jet_E_eta[isel] = new TH2F(Form("h_jet_%s_E_eta",jettype[isel].Data()),"",nBinsP, binningP,80,-4,4);
    
    for (Int_t eT = 0; eT < nEta+1; eT++){
      if(!h_jetscale_pT[isel][eT])h_jetscale_pT[isel][eT] = new TH2F(Form("h_jetscale_%s_pT_%d",jettype[isel].Data(), eT),"",150,0,30,200,-1,1);
      if(!h_jetscale_E[isel][eT])h_jetscale_E[isel][eT]  = new TH2F(Form("h_jetscale_%s_E_%d",jettype[isel].Data(), eT),"",nBinsP, binningP,200,-1,1);

      if(!h_JES_pT[isel][eT])h_JES_pT[isel][eT]         = new TH1D(Form("h_JES_%s_pT_%d",jettype[isel].Data(), eT),"",150,0,30);
      if(!h_JER_pT[isel][eT])h_JER_pT[isel][eT]         = new TH1D(Form("h_JER_%s_pT_%d",jettype[isel].Data(), eT),"",150,0,30);
      if(!h_JES_E[isel][eT])h_JES_E[isel][eT]           = new TH1D(Form("h_JES_%s_E_%d",jettype[isel].Data(), eT),"",nP, partP);
      if(!h_JER_E[isel][eT])h_JER_E[isel][eT]           = new TH1D(Form("h_JER_%s_E_%d",jettype[isel].Data(), eT),"",nP, partP);
    }
  }

  for (std::size_t i = 0; i < std::get<1>(recjets).size(); i++) {
    // if(verbosityJRH>1)std::cout << "jet " << i << ": "<< std::get<1>(recjets)[i].pt() << " " << std::get<1>(recjets)[i].rap() << " " << std::get<1>(recjets)[i].phi() << "  Ntrue: " << std::get<1>(truejets).size() << endl;
    for (std::size_t j = 0; j < std::get<1>(truejets).size(); j++) {
      Double_t deltaRTrueRec = std::get<1>(truejets)[j].delta_R(std::get<1>(recjets)[i]);
      // cout << deltaRTrueRec << endl;
      Int_t et = 0;
      while (partEta[et+1] < std::get<1>(truejets)[j].eta() && et < nEta) et++;
      
      if(deltaRTrueRec<0.5){
        // if(verbosityJRH>1) cout << "rec jet " << i << " (" << std::get<1>(recjets)[i].pt() << "pT)  matched with true jet " << j << " (" << std::get<1>(truejets)[i].pt() << "pT) with dR = " << deltaRTrueRec << endl;
        h_jetscale_E[select][et]->Fill(std::get<1>(truejets)[j].E(),(std::get<1>(recjets)[i].E()-std::get<1>(truejets)[j].E())/std::get<1>(truejets)[j].E());
        h_jetscale_pT[select][et]->Fill(std::get<1>(truejets)[j].pt(),(std::get<1>(recjets)[i].pt()-std::get<1>(truejets)[j].pt())/std::get<1>(truejets)[j].pt());
        h_jetscale_E[select][nEta]->Fill(std::get<1>(truejets)[j].E(),(std::get<1>(recjets)[i].E()-std::get<1>(truejets)[j].E())/std::get<1>(truejets)[j].E());
        h_jetscale_pT[select][nEta]->Fill(std::get<1>(truejets)[j].pt(),(std::get<1>(recjets)[i].pt()-std::get<1>(truejets)[j].pt())/std::get<1>(truejets)[j].pt());
        
        h_MCjet_E_eta[select]->Fill(std::get<1>(truejets)[j].E(),std::get<1>(truejets)[j].eta());
        h_jet_E_eta[select]->Fill(std::get<1>(recjets)[j].E(),std::get<1>(recjets)[j].eta());
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
        for (Int_t i=1; i < nP+1; i++){
          TH1D* projectionYdummy = (TH1D*)h_jetscale_pT[isel][eT]->ProjectionY(Form("projectionYdummy%d%d",isel,i), i,i+1,"e");
          // projectionYdummy->GetXaxis()->SetRangeUser(-0.6,0);
          h_JES_pT[isel][eT]->SetBinContent(i,projectionYdummy->GetMean());
          h_JES_pT[isel][eT]->SetBinError(i,projectionYdummy->GetMeanError());
          h_JER_pT[isel][eT]->SetBinContent(i,projectionYdummy->GetStdDev()); //GetMeanError()
          h_JER_pT[isel][eT]->SetBinError(i,projectionYdummy->GetStdDevError()); //GetMeanError()
        }
      }
      if(h_jetscale_E[isel][eT]){
        for (Int_t i=1; i < nP+1; i++){
          TH1D* projectionYdummy = (TH1D*)h_jetscale_E[isel][eT]->ProjectionY(Form("projectionYdummy%d%d",isel,i), h_jetscale_E[isel][eT]->GetXaxis()->FindBin(partP[i-1]+0.001), h_jetscale_E[isel][eT]->GetXaxis()->FindBin(partP[i]+0.001),"e");
          // projectionYdummy->GetXaxis()->SetRangeUser(-0.6,0);
          h_JES_E[isel][eT]->SetBinContent(i,projectionYdummy->GetMean());
          h_JES_E[isel][eT]->SetBinError(i,projectionYdummy->GetMeanError());
          h_JER_E[isel][eT]->SetBinContent(i,projectionYdummy->GetStdDev()); //GetMeanError()
          h_JER_E[isel][eT]->SetBinError(i,projectionYdummy->GetStdDevError()); //GetMeanError()
        }
      }
    }
  }
  // define output file
  TFile* fileOutput = new TFile(Form("%s/output_JRH.root",outputDir.Data()),"RECREATE");

  // write histograms
  for( int isel=0;isel<njettypes;isel++){
    if(h_MCjet_E_eta[isel])h_MCjet_E_eta[isel]->Write();
    if(h_jet_E_eta[isel])h_jet_E_eta[isel]->Write();
    for (Int_t eT = 0; eT < nEta+1; eT++){
      if(h_jetscale_E[isel][eT])h_jetscale_E[isel][eT]->Write();
      if(h_jetscale_pT[isel][eT])h_jetscale_pT[isel][eT]->Write();
      if(h_JES_pT[isel][eT])h_JES_pT[isel][eT]->Write();
      if(h_JES_E[isel][eT])h_JES_E[isel][eT]->Write();
      if(h_JER_pT[isel][eT])h_JER_pT[isel][eT]->Write();
      if(h_JER_E[isel][eT])h_JER_E[isel][eT]->Write();
    }
  }
  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
