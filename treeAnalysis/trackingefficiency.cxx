
// ANCHOR debug output verbosity
Int_t verbosityTRKEFF = 0;

// ANCHOR create histograms globally
TString str_TRKEFF_mcparticles[5] = {"electron", "cpion", "proton", "ckaon", "muon"};
int int_TRKEFF_mcparticles_PDG[5] = {11 /*e*/,211  /*pi*/,  2212/*p*/,  321/*K*/,  13/*mu*/};

TH2F*  h_chargedpart_MC_pT 	= new TH2F("h_chargedpart_MC_pT", "", 60, 0,30,80, -4.0,4.0);
TH2F*  h_particle_MC_pT[5];
TH2F*  h_chargedpart_rec_pT 	= new TH2F("h_chargedpart_rec_pT", "", 60, 0,30,80, -4.0,4.0);
TH2F*  h_chargedpart_rec_truepT 	= new TH2F("h_chargedpart_rec_truepT", "", 60, 0,30,80, -4.0,4.0);
TH2F*  h_particle_rec_pT[5];
TH2F*  h_particle_rec_truepT[5];

TH2F*  h_chargedpart_MC_p 	= new TH2F("h_chargedpart_MC_p", "", 200, 0,100,80, -4.0,4.0);
TH2F*  h_particle_MC_p[5];
TH2F*  h_chargedpart_rec_p 	= new TH2F("h_chargedpart_rec_p", "", 200, 0,100,80, -4.0,4.0);
TH2F*  h_chargedpart_rec_truep 	= new TH2F("h_chargedpart_rec_truep", "", 200, 0,100,80, -4.0,4.0);
TH2F*  h_particle_rec_p[5];
TH2F*  h_particle_rec_truep[5];

// ANCHOR main function to be called in event loop
void trackingefficiency(){
  // cout << "new evt" << endl;
  for(int imc=0; imc<_nMCPart; imc++){
    TVector3 mcpartvec(_mcpart_px[imc],_mcpart_py[imc],_mcpart_pz[imc]);
    float trueeta = mcpartvec.Eta();

    // cout << "\tMCID " << imc << "\tPDG: " << _mcpart_PDG[imc] << "\tpT: " << TMath::Sqrt(TMath::Power(_mcpart_px[imc],2)+TMath::Power(_mcpart_py[imc],2)) << "\tEta " << trueeta << endl;
    for(int ipart=0;ipart<5;ipart++){
      if(!h_particle_MC_pT[ipart]) h_particle_MC_pT[ipart] 	= new TH2F(Form("h_%s_MC_pT",str_TRKEFF_mcparticles[ipart].Data()), "", 60, 0,30,80, -4.0,4.0);
      if(!h_particle_rec_pT[ipart]) h_particle_rec_pT[ipart] 	= new TH2F(Form("h_%s_pT",str_TRKEFF_mcparticles[ipart].Data()), "", 60, 0,30,80, -4.0,4.0);
      if(!h_particle_rec_truepT[ipart]) h_particle_rec_truepT[ipart] 	= new TH2F(Form("h_%s_truepT",str_TRKEFF_mcparticles[ipart].Data()), "", 60, 0,30,80, -4.0,4.0);

      if(!h_particle_MC_p[ipart]) h_particle_MC_p[ipart] 	= new TH2F(Form("h_%s_MC_p",str_TRKEFF_mcparticles[ipart].Data()), "", 200, 0,100,80, -4.0,4.0);
      if(!h_particle_rec_p[ipart]) h_particle_rec_p[ipart] 	= new TH2F(Form("h_%s_p",str_TRKEFF_mcparticles[ipart].Data()), "", 200, 0,100,80, -4.0,4.0);
      if(!h_particle_rec_truep[ipart]) h_particle_rec_truep[ipart] 	= new TH2F(Form("h_%s_truep",str_TRKEFF_mcparticles[ipart].Data()), "", 200, 0,100,80, -4.0,4.0);

      if(abs(_mcpart_PDG[imc]) == int_TRKEFF_mcparticles_PDG[ipart]){
        h_particle_MC_pT[ipart]->Fill(TMath::Sqrt(TMath::Power(_mcpart_px[imc],2)+TMath::Power(_mcpart_py[imc],2)),trueeta);
        h_chargedpart_MC_pT->Fill(TMath::Sqrt(TMath::Power(_mcpart_px[imc],2)+TMath::Power(_mcpart_py[imc],2)),trueeta);

        h_particle_MC_p[ipart]->Fill(TMath::Sqrt(TMath::Power(_mcpart_px[imc],2)+TMath::Power(_mcpart_py[imc],2)+TMath::Power(_mcpart_pz[imc],2)),trueeta);
        h_chargedpart_MC_p->Fill(TMath::Sqrt(TMath::Power(_mcpart_px[imc],2)+TMath::Power(_mcpart_py[imc],2)+TMath::Power(_mcpart_pz[imc],2)),trueeta);
      }
    }
  }

  // see if true reco. track was found
  for(Int_t itrk=0; itrk<_nTracks; itrk++){
    TVector3 recpartvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
    // cout << itrk << endl;
    float receta = recpartvec.Eta();
    // cout << "\tTRKTRUEID " << _track_trueID[itrk]-1 << "\tPDG: " << _mcpart_PDG[(int)_track_trueID[itrk]-1] << "\tpT: " <<   TMath::Sqrt(TMath::Power(_track_px[itrk],2)+TMath::Power(_track_py[itrk],2)) << "\tEta " << receta << endl;
    float pt = TMath::Sqrt(TMath::Power(_track_px[itrk],2)+TMath::Power(_track_py[itrk],2));
    float truept = TMath::Sqrt(TMath::Power(_mcpart_px[(int)_track_trueID[itrk]-1],2)+TMath::Power(_mcpart_py[(int)_track_trueID[itrk]-1],2));
    float pmom = TMath::Sqrt(TMath::Power(_track_px[itrk],2)+TMath::Power(_track_py[itrk],2)+TMath::Power(_track_pz[itrk],2));
    float truepmom = TMath::Sqrt(TMath::Power(_mcpart_px[(int)_track_trueID[itrk]-1],2)+TMath::Power(_mcpart_py[(int)_track_trueID[itrk]],2)+TMath::Power(_mcpart_pz[(int)_track_trueID[itrk]-1],2));
    for(int ipart=0;ipart<5;ipart++){
      if(abs(_mcpart_PDG[(int)_track_trueID[itrk]-1])==int_TRKEFF_mcparticles_PDG[ipart]){
        h_particle_rec_pT[ipart]->Fill(pt,receta);
        h_particle_rec_truepT[ipart]->Fill(truept,receta);
        h_chargedpart_rec_pT->Fill(pt,receta);
        h_chargedpart_rec_truepT->Fill(truept,receta);

        h_particle_rec_p[ipart]->Fill(pmom,receta);
        h_particle_rec_truep[ipart]->Fill(truepmom,receta);
        h_chargedpart_rec_p->Fill(pmom,receta);
        h_chargedpart_rec_truep->Fill(truepmom,receta);
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
  if(h_chargedpart_rec_pT) h_chargedpart_rec_pT->Write();
  if(h_chargedpart_rec_truepT) h_chargedpart_rec_truepT->Write();
  if(h_chargedpart_MC_pT && h_chargedpart_rec_pT && h_chargedpart_rec_truepT) {
    TH2F*  h_chargedpart_eff_pT = (TH2F*)h_chargedpart_rec_pT->Clone("h_chargedpart_eff_pT");
    h_chargedpart_eff_pT->Divide(h_chargedpart_MC_pT);
    h_chargedpart_eff_pT->Write();
    TH2F*  h_chargedpart_eff_truepT = (TH2F*)h_chargedpart_rec_truepT->Clone("h_chargedpart_eff_truepT");
    h_chargedpart_eff_truepT->Divide(h_chargedpart_MC_pT);
    h_chargedpart_eff_truepT->Write();
  }
  if(h_chargedpart_MC_p) h_chargedpart_MC_p->Write();
  if(h_chargedpart_rec_p) h_chargedpart_rec_p->Write();
  if(h_chargedpart_rec_truep) h_chargedpart_rec_truep->Write();
  if(h_chargedpart_MC_p && h_chargedpart_rec_p && h_chargedpart_rec_truep) {
    TH2F*  h_chargedpart_eff_p = (TH2F*)h_chargedpart_rec_p->Clone("h_chargedpart_eff_p");
    h_chargedpart_eff_p->Divide(h_chargedpart_MC_p);
    h_chargedpart_eff_p->Write();
    TH2F*  h_chargedpart_eff_truep = (TH2F*)h_chargedpart_rec_truep->Clone("h_chargedpart_eff_truep");
    h_chargedpart_eff_truep->Divide(h_chargedpart_MC_p);
    h_chargedpart_eff_truep->Write();
  }
  TH2F*  h_particle_eff_pT[4];
  TH2F*  h_particle_eff_truepT[4];
  TH2F*  h_particle_eff_p[4];
  TH2F*  h_particle_eff_truep[4];
  for(int ipart=0;ipart<5;ipart++){
    if(h_particle_MC_pT[ipart]) h_particle_MC_pT[ipart]->Write();
    if(h_particle_rec_pT[ipart]) h_particle_rec_pT[ipart]->Write();
    if(h_particle_rec_truepT[ipart]) h_particle_rec_truepT[ipart]->Write();
    if(h_particle_MC_pT[ipart] && h_particle_rec_pT[ipart] && h_particle_rec_truepT[ipart]) {
      h_particle_eff_pT[ipart] = (TH2F*)h_particle_rec_pT[ipart]->Clone(Form("h_%s_eff_pT",str_TRKEFF_mcparticles[ipart].Data()));
      h_particle_eff_pT[ipart]->Divide(h_particle_MC_pT[ipart]);
      h_particle_eff_pT[ipart]->Write();
      h_particle_eff_truepT[ipart] = (TH2F*)h_particle_rec_truepT[ipart]->Clone(Form("h_%s_eff_truepT",str_TRKEFF_mcparticles[ipart].Data()));
      h_particle_eff_truepT[ipart]->Divide(h_particle_MC_pT[ipart]);
      h_particle_eff_truepT[ipart]->Write();
    }
    if(h_particle_MC_p[ipart]) h_particle_MC_p[ipart]->Write();
    if(h_particle_rec_p[ipart]) h_particle_rec_p[ipart]->Write();
    if(h_particle_rec_truep[ipart]) h_particle_rec_truep[ipart]->Write();
    if(h_particle_MC_p[ipart] && h_particle_rec_p[ipart] && h_particle_rec_truep[ipart]) {
      h_particle_eff_p[ipart] = (TH2F*)h_particle_rec_p[ipart]->Clone(Form("h_%s_eff_p",str_TRKEFF_mcparticles[ipart].Data()));
      h_particle_eff_p[ipart]->Divide(h_particle_MC_p[ipart]);
      h_particle_eff_p[ipart]->Write();
      h_particle_eff_truep[ipart] = (TH2F*)h_particle_rec_truep[ipart]->Clone(Form("h_%s_eff_truep",str_TRKEFF_mcparticles[ipart].Data()));
      h_particle_eff_truep[ipart]->Divide(h_particle_MC_p[ipart]);
      h_particle_eff_truep[ipart]->Write();
    }
  }
  // write output file
  fileOutput->Write();
  fileOutput->Close();
}