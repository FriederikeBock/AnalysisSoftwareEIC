
// ANCHOR debug output verbosity
Int_t verbosityTRKEFF = 0;

// ANCHOR create histograms globally
TString str_TRKEFF_mcparticles[4] = {"electron", "cpion", "proton", "ckaon"};
int int_TRKEFF_mcparticles_PDG[4] = {11 /*e*/,211  /*pi*/,  2212/*p*/,  321/*K*/};

TH1F*  h_chargedpart_MC_pT 	= new TH1F("h_chargedpart_MC_pT", "", 60, 0,30);
TH1F*  h_particle_MC_pT[5];
TH1F*  h_chargedpart_pT 	= new TH1F("h_chargedpart_pT", "", 60, 0,30);
TH1F*  h_particle_pT[5];

// ANCHOR main function to be called in event loop
void trackingefficiency(){

  for(int imc=0; imc<_nMCPart; imc++){
    TVector3 mcpartvec(_mcpart_px[imc],_mcpart_py[imc],_mcpart_pz[imc]);

    if(mcpartvec.Eta()<3.0) continue;

    for(int ipart=0;ipart<4;ipart++){
      if(!h_particle_MC_pT[ipart]) h_particle_MC_pT[ipart] 	= new TH1F(Form("h_%s_MC_pT",str_TRKEFF_mcparticles[ipart].Data()), "", 60, 0,30);
      if(!h_particle_pT[ipart]) h_particle_pT[ipart] 	= new TH1F(Form("h_%s_pT",str_TRKEFF_mcparticles[ipart].Data()), "", 60, 0,30);

      if(_mcpart_PDG[imc] == int_TRKEFF_mcparticles_PDG[ipart]){
        h_particle_MC_pT[ipart]->Fill(TMath::Sqrt(TMath::Power(_mcpart_px[imc],2)+TMath::Power(_mcpart_py[imc],2)));
        h_chargedpart_MC_pT->Fill(TMath::Sqrt(TMath::Power(_mcpart_px[imc],2)+TMath::Power(_mcpart_py[imc],2)));

        // // see if true reco. track was found
        // for(Int_t itrk=0; itrk<_nTracks; itrk++){
        //   if(_track_trueID[itrk]==_mcpart_ID[imc]){
        //     TVector3 recpartvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
        //     if(recpartvec.Eta()<3.0) continue;
        //     h_particle_pT[ipart]->Fill(TMath::Sqrt(TMath::Power(_track_px[itrk],2)+TMath::Power(_track_py[itrk],2)));
        //     h_chargedpart_pT->Fill(TMath::Sqrt(TMath::Power(_track_px[itrk],2)+TMath::Power(_track_py[itrk],2)));
        //   }
        // }
      }
    }
  }

  // see if true reco. track was found
  for(Int_t itrk=0; itrk<_nTracks; itrk++){
    for(int ipart=0;ipart<4;ipart++){
      if(_mcpart_PDG[(int)_track_trueID[itrk]]==int_TRKEFF_mcparticles_PDG[ipart]){
        TVector3 recpartvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
        if(recpartvec.Eta()<3.0) continue;
        h_particle_pT[ipart]->Fill(TMath::Sqrt(TMath::Power(_track_px[itrk],2)+TMath::Power(_track_py[itrk],2)));
        h_chargedpart_pT->Fill(TMath::Sqrt(TMath::Power(_track_px[itrk],2)+TMath::Power(_track_py[itrk],2)));
      }
    }
  }

}

// ANCHOR save function after event loop
void trackingefficiencyhistosSave(){
  // make output directory
    gSystem->Exec("mkdir -p "+outputDir + "/trackingefficiency");

  // gSystem->Exec("mkdir -p treeProcessing/trackingefficiency");
  // define output file
  TFile* fileOutput = new TFile(Form("%s/trackingefficiency/output_TRKEFF.root",outputDir.Data()),"RECREATE");

  // write histograms
  if(h_chargedpart_MC_pT) h_chargedpart_MC_pT->Write();
  if(h_chargedpart_pT) h_chargedpart_pT->Write();
  if(h_chargedpart_MC_pT && h_chargedpart_pT) {
    TH1F*  h_chargedpart_eff_pT = (TH1F*)h_chargedpart_pT->Clone("h_chargedpart_eff_pT");
    h_chargedpart_eff_pT->Divide(h_chargedpart_MC_pT);
    h_chargedpart_eff_pT->Write();
  }
  TH1F*  h_particle_eff_pT[4];
  for(int ipart=0;ipart<4;ipart++){
    if(h_particle_MC_pT[ipart]) h_particle_MC_pT[ipart]->Write();
    if(h_particle_pT[ipart]) h_particle_pT[ipart]->Write();
    if(h_particle_MC_pT[ipart] && h_particle_pT[ipart]) {
      h_particle_eff_pT[ipart] = (TH1F*)h_particle_pT[ipart]->Clone(Form("h_%s_eff_pT",str_TRKEFF_mcparticles[ipart].Data()));
      h_particle_eff_pT[ipart]->Divide(h_particle_MC_pT[ipart]);
      h_particle_eff_pT[ipart]->Write();
    }
  }

  // write output file
  fileOutput->Write();
  fileOutput->Close();
}