
// ANCHOR debug output verbosity
Int_t verbosityJRH = 1;

// ANCHOR create histograms globally
// Double_t xbinsDiff[16] = {0.,0.5, 1.2.,4.,6.,8.,10.,12.,14.,16.,18.,  20.,25.,30.,35.,40.,45.};

TH2F* h_jetscale_pT = new TH2F(Form("h_jetscale_pT"),"",150,0,30,200,-1,1);
TH2F* h_jetscale_E  = new TH2F(Form("h_jetscale_E"),"",200,0,200,200,-1,1);

TH1D *h_JES_pT      = new TH1D(Form("h_JES_pT"),"",150,0,30);
TH1D *h_JER_pT      = new TH1D(Form("h_JER_pT"),"",150,0,30);
TH1D *h_JES_E      = new TH1D(Form("h_JES_E"),"",200,0,200);
TH1D *h_JER_E      = new TH1D(Form("h_JER_E"),"",200,0,200);

// ANCHOR main function to be called in event loop
void jetresolutionhistos(std::tuple<std::shared_ptr<fastjet::ClusterSequenceArea>, std::vector<fastjet::PseudoJet>> recjets, std::tuple<std::shared_ptr<fastjet::ClusterSequenceArea>, std::vector<fastjet::PseudoJet>> truejets){


  for (std::size_t i = 0; i < std::get<1>(recjets).size(); i++) {
    // if(verbosityJRH>1)std::cout << "jet " << i << ": "<< std::get<1>(recjets)[i].pt() << " " << std::get<1>(recjets)[i].rap() << " " << std::get<1>(recjets)[i].phi() << "  Ntrue: " << std::get<1>(truejets).size() << endl;
    for (std::size_t j = 0; j < std::get<1>(truejets).size(); j++) {
      Double_t deltaRTrueRec = std::get<1>(truejets)[i].delta_R(std::get<1>(recjets)[j]);
      // cout << deltaRTrueRec << endl;
      if(deltaRTrueRec<0.5){
        // if(verbosityJRH>1) cout << "rec jet " << i << " (" << std::get<1>(recjets)[i].pt() << "pT)  matched with true jet " << j << " (" << std::get<1>(truejets)[i].pt() << "pT) with dR = " << deltaRTrueRec << endl;
        h_jetscale_E->Fill(std::get<1>(truejets)[j].E(),(std::get<1>(recjets)[i].E()-std::get<1>(truejets)[j].E())/std::get<1>(truejets)[j].E());
        h_jetscale_pT->Fill(std::get<1>(truejets)[j].pt(),(std::get<1>(recjets)[i].pt()-std::get<1>(truejets)[j].pt())/std::get<1>(truejets)[j].pt());
      }
    }
  } // jet loop end

}

// ANCHOR save function after event loop
void jetresolutionhistosSave(){

  // make projections for JES and JER
  for (Int_t i=1; i < h_jetscale_pT->GetNbinsX(); i++){
    TH1D* projectionYdummy = (TH1D*)h_jetscale_pT->ProjectionY(Form("projectionYdummy%d",i), i,i+1,"e");
    // projectionYdummy->GetXaxis()->SetRangeUser(-0.6,0);
    h_JES_pT->SetBinContent(i,projectionYdummy->GetMean());
    h_JES_pT->SetBinError(i,projectionYdummy->GetMeanError());
    h_JER_pT->SetBinContent(i,projectionYdummy->GetStdDev()); //GetMeanError()
    h_JER_pT->SetBinError(i,projectionYdummy->GetStdDevError()); //GetMeanError()
  }
  for (Int_t i=1; i < h_jetscale_E->GetNbinsX(); i++){
    TH1D* projectionYdummy = (TH1D*)h_jetscale_E->ProjectionY(Form("projectionYdummy%d",i), i,i+1,"e");
    // projectionYdummy->GetXaxis()->SetRangeUser(-0.6,0);
    h_JES_E->SetBinContent(i,projectionYdummy->GetMean());
    h_JES_E->SetBinError(i,projectionYdummy->GetMeanError());
    h_JER_E->SetBinContent(i,projectionYdummy->GetStdDev()); //GetMeanError()
    h_JER_E->SetBinError(i,projectionYdummy->GetStdDevError()); //GetMeanError()
  }

  // make output directory
    gSystem->Exec("mkdir -p "+outputDir + "/jetresolutionhistosSave");
  // define output file
  TFile* fileOutput = new TFile(Form("%s/jetresolutionhistosSave/output_JRH.root",outputDir.Data()),"RECREATE");

  // write histograms
  h_jetscale_E->Write();
  h_jetscale_pT->Write();
  h_JES_pT->Write();
  h_JES_E->Write();
  h_JER_pT->Write();
  h_JER_E->Write();

  // write output file
  fileOutput->Write();
  fileOutput->Close();
}