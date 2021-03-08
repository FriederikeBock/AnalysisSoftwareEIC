
// ANCHOR debug output verbosity
Int_t verbosityJRH = 1;

// ANCHOR create histograms globally
// Double_t xbinsDiff[16] = {0.,0.5, 1.2.,4.,6.,8.,10.,12.,14.,16.,18.,  20.,25.,30.,35.,40.,45.};
const int njettypes = 5;
TString jettype[njettypes] = {"track", "full","hcal","calo","all"};

TH2F* h_jetscale_pT[njettypes] = {NULL};
TH2F* h_jetscale_E[njettypes] = {NULL};

TH1D *h_JES_pT[njettypes] = {NULL};
TH1D *h_JER_pT[njettypes] = {NULL};
TH1D *h_JES_E[njettypes] = {NULL};
TH1D *h_JER_E[njettypes] = {NULL};

// ANCHOR main function to be called in event loop
void jetresolutionhistos(std::tuple<std::shared_ptr<fastjet::ClusterSequenceArea>, std::vector<fastjet::PseudoJet>> recjets, std::tuple<std::shared_ptr<fastjet::ClusterSequenceArea>, std::vector<fastjet::PseudoJet>> truejets, int select){

  for( int isel=0;isel<njettypes;isel++){
    if(!h_jetscale_pT[isel])h_jetscale_pT[isel] = new TH2F(Form("h_jetscale_%s_pT",jettype[isel].Data()),"",150,0,30,200,-1,1);
    if(!h_jetscale_E[isel])h_jetscale_E[isel]  = new TH2F(Form("h_jetscale_%s_E",jettype[isel].Data()),"",200,0,200,200,-1,1);

    if(!h_JES_pT[isel])h_JES_pT[isel]      = new TH1D(Form("h_JES_%s_pT",jettype[isel].Data()),"",150,0,30);
    if(!h_JER_pT[isel])h_JER_pT[isel]      = new TH1D(Form("h_JER_%s_pT",jettype[isel].Data()),"",150,0,30);
    if(!h_JES_E[isel])h_JES_E[isel]      = new TH1D(Form("h_JES_%s_E",jettype[isel].Data()),"",200,0,200);
    if(!h_JER_E[isel])h_JER_E[isel]      = new TH1D(Form("h_JER_%s_E",jettype[isel].Data()),"",200,0,200);
  }

  for (std::size_t i = 0; i < std::get<1>(recjets).size(); i++) {
    // if(verbosityJRH>1)std::cout << "jet " << i << ": "<< std::get<1>(recjets)[i].pt() << " " << std::get<1>(recjets)[i].rap() << " " << std::get<1>(recjets)[i].phi() << "  Ntrue: " << std::get<1>(truejets).size() << endl;
    for (std::size_t j = 0; j < std::get<1>(truejets).size(); j++) {
      Double_t deltaRTrueRec = std::get<1>(truejets)[j].delta_R(std::get<1>(recjets)[i]);
      // cout << deltaRTrueRec << endl;
      if(deltaRTrueRec<0.5){
        // if(verbosityJRH>1) cout << "rec jet " << i << " (" << std::get<1>(recjets)[i].pt() << "pT)  matched with true jet " << j << " (" << std::get<1>(truejets)[i].pt() << "pT) with dR = " << deltaRTrueRec << endl;
        h_jetscale_E[select]->Fill(std::get<1>(truejets)[j].E(),(std::get<1>(recjets)[i].E()-std::get<1>(truejets)[j].E())/std::get<1>(truejets)[j].E());
        h_jetscale_pT[select]->Fill(std::get<1>(truejets)[j].pt(),(std::get<1>(recjets)[i].pt()-std::get<1>(truejets)[j].pt())/std::get<1>(truejets)[j].pt());
      }
    }
  } // jet loop end

}

// ANCHOR save function after event loop
void jetresolutionhistosSave(){
  for( int isel=0;isel<njettypes;isel++){
    if(h_jetscale_pT[isel]){
      // make projections for JES and JER
      for (Int_t i=1; i < h_jetscale_pT[isel]->GetNbinsX(); i++){
        TH1D* projectionYdummy = (TH1D*)h_jetscale_pT[isel]->ProjectionY(Form("projectionYdummy%d%d",isel,i), i,i+1,"e");
        // projectionYdummy->GetXaxis()->SetRangeUser(-0.6,0);
        h_JES_pT[isel]->SetBinContent(i,projectionYdummy->GetMean());
        h_JES_pT[isel]->SetBinError(i,projectionYdummy->GetMeanError());
        h_JER_pT[isel]->SetBinContent(i,projectionYdummy->GetStdDev()); //GetMeanError()
        h_JER_pT[isel]->SetBinError(i,projectionYdummy->GetStdDevError()); //GetMeanError()
      }
    }
    if(h_jetscale_E[isel]){
      for (Int_t i=1; i < h_jetscale_E[isel]->GetNbinsX(); i++){
        TH1D* projectionYdummy = (TH1D*)h_jetscale_E[isel]->ProjectionY(Form("projectionYdummy%d%d",isel,i), i,i+1,"e");
        // projectionYdummy->GetXaxis()->SetRangeUser(-0.6,0);
        h_JES_E[isel]->SetBinContent(i,projectionYdummy->GetMean());
        h_JES_E[isel]->SetBinError(i,projectionYdummy->GetMeanError());
        h_JER_E[isel]->SetBinContent(i,projectionYdummy->GetStdDev()); //GetMeanError()
        h_JER_E[isel]->SetBinError(i,projectionYdummy->GetStdDevError()); //GetMeanError()
      }
    }
  }
  // make output directory
    gSystem->Exec("mkdir -p "+outputDir + "/jetresolutionhistosSave");
  // define output file
  TFile* fileOutput = new TFile(Form("%s/jetresolutionhistosSave/output_JRH.root",outputDir.Data()),"RECREATE");

  // write histograms
  for( int isel=0;isel<njettypes;isel++){
    if(h_jetscale_E[isel])h_jetscale_E[isel]->Write();
    if(h_jetscale_pT[isel])h_jetscale_pT[isel]->Write();
    if(h_JES_pT[isel])h_JES_pT[isel]->Write();
    if(h_JES_E[isel])h_JES_E[isel]->Write();
    if(h_JER_pT[isel])h_JER_pT[isel]->Write();
    if(h_JER_E[isel])h_JER_E[isel]->Write();
  }
  // write output file
  fileOutput->Write();
  fileOutput->Close();
}