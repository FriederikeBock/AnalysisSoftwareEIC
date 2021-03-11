
// ANCHOR debug output verbosity
Int_t verbosityTRKEFF = 0;

// ANCHOR create histograms globally
const int nPart_TRKEFF = 5;
TString str_TRKEFF_mcparticles[nPart_TRKEFF] = {"electron", "cpion", "proton", "ckaon", "muon"};
int int_TRKEFF_mcparticles_PDG[nPart_TRKEFF] = {11 /*e*/,211  /*pi*/,  2212/*p*/,  321/*K*/,  13/*mu*/};

TH2F*  h_chargedpart_MC_pT 	= new TH2F("h_chargedpart_MC_pT", "", 60, 0,30,80, -4.0,4.0);
TH2F*  h_particle_MC_pT[nPart_TRKEFF];
TH2F*  h_chargedpart_rec_pT 	= new TH2F("h_chargedpart_rec_pT", "", 60, 0,30,80, -4.0,4.0);
TH2F*  h_chargedpart_rec_truepT 	= new TH2F("h_chargedpart_rec_truepT", "", 60, 0,30,80, -4.0,4.0);
TH2F*  h_particle_rec_pT[nPart_TRKEFF];
TH2F*  h_particle_rec_truepT[nPart_TRKEFF];

TH2F*  h_chargedpart_MC_p 	= new TH2F("h_chargedpart_MC_p", "", nBinsP, binningP, 80, -4.0,4.0);
TH2F*  h_particle_MC_p[nPart_TRKEFF];
TH2F*  h_chargedpart_rec_p 	= new TH2F("h_chargedpart_rec_p", "", nBinsP, binningP, 80, -4.0,4.0);
TH2F*  h_chargedpart_rec_truep 	= new TH2F("h_chargedpart_rec_truep", "", nBinsP, binningP, 80, -4.0,4.0);
TH2F*  h_particle_rec_p[nPart_TRKEFF];
TH2F*  h_particle_rec_truep[nPart_TRKEFF];

TH2F*  h_chargedpart_MC_E 	= new TH2F("h_chargedpart_MC_E", "", nBinsP, binningP, 80, -4.0,4.0);
TH2F*  h_particle_MC_E[nPart_TRKEFF];

// ANCHOR main function to be called in event loop
void trackingefficiency(){
  // cout << "new evt" << endl;
  for(int imc=0; imc<_nMCPart; imc++){
    TVector3 mcpartvec(_mcpart_px[imc],_mcpart_py[imc],_mcpart_pz[imc]);
    float trueeta = mcpartvec.Eta();
    float truept  = mcpartvec.Pt();
    float truep   = mcpartvec.Mag();

    // cout << "\tMCID " << imc << "\tPDG: " << _mcpart_PDG[imc] << "\tpT: " << TMath::Sqrt(TMath::Power(_mcpart_px[imc],2)+TMath::Power(_mcpart_py[imc],2)) << "\tEta " << trueeta << endl;
    for(int ipart=0;ipart<nPart_TRKEFF;ipart++){
      if(!h_particle_MC_pT[ipart]) h_particle_MC_pT[ipart] 	= new TH2F(Form("h_%s_MC_pT",str_TRKEFF_mcparticles[ipart].Data()), "", 60, 0,30,80, -4.0,4.0);
      if(!h_particle_rec_pT[ipart]) h_particle_rec_pT[ipart] 	= new TH2F(Form("h_%s_rec_pT",str_TRKEFF_mcparticles[ipart].Data()), "", 60, 0,30,80, -4.0,4.0);
      if(!h_particle_rec_truepT[ipart]) h_particle_rec_truepT[ipart] 	= new TH2F(Form("h_%s_rec_truepT",str_TRKEFF_mcparticles[ipart].Data()), "", 60, 0,30,80, -4.0,4.0);

      if(!h_particle_MC_p[ipart]) h_particle_MC_p[ipart] 	= new TH2F(Form("h_%s_MC_p",str_TRKEFF_mcparticles[ipart].Data()), "", nBinsP, binningP, 80, -4.0,4.0);
      if(!h_particle_rec_p[ipart]) h_particle_rec_p[ipart] 	= new TH2F(Form("h_%s_rec_p",str_TRKEFF_mcparticles[ipart].Data()), "", nBinsP, binningP, 80, -4.0,4.0);
      if(!h_particle_rec_truep[ipart]) h_particle_rec_truep[ipart] 	= new TH2F(Form("h_%s_rec_truep",str_TRKEFF_mcparticles[ipart].Data()), "", nBinsP, binningP, 80, -4.0,4.0);
      if(!h_particle_MC_E[ipart]) h_particle_MC_E[ipart] 	= new TH2F(Form("h_%s_MC_E",str_TRKEFF_mcparticles[ipart].Data()), "", nBinsP, binningP, 80, -4.0,4.0);
      
      if(abs(_mcpart_PDG[imc]) == int_TRKEFF_mcparticles_PDG[ipart]){
        h_particle_MC_pT[ipart]->Fill(truept,trueeta);
        h_chargedpart_MC_pT->Fill(truept,trueeta);

        h_particle_MC_p[ipart]->Fill(truep,trueeta);
        h_chargedpart_MC_p->Fill(truep,trueeta);
        h_particle_MC_E[ipart]->Fill(_mcpart_E[imc],trueeta);
        h_chargedpart_MC_E->Fill(_mcpart_E[imc],trueeta);
      }
    }
  }

  // see if true reco. track was found
  for(Int_t itrk=0; itrk<_nTracks; itrk++){
    TVector3 recpartvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
    TVector3 mcpartvec(_mcpart_px[(int)_track_trueID[itrk]-1],_mcpart_py[(int)_track_trueID[itrk]-1],_mcpart_pz[(int)_track_trueID[itrk]-1]);
    // cout << itrk << endl;
    float receta = recpartvec.Eta();
    float trueeta = mcpartvec.Eta();
    // cout << "\tTRKTRUEID " << _track_trueID[itrk]-1 << "\tPDG: " << _mcpart_PDG[(int)_track_trueID[itrk]-1] << "\tpT: " <<   TMath::Sqrt(TMath::Power(_track_px[itrk],2)+TMath::Power(_track_py[itrk],2)) << "\tEta " << receta << endl;
    float pt = recpartvec.Pt();
    float truept = mcpartvec.Pt();
    float pmom = recpartvec.Mag();
    float truepmom = mcpartvec.Mag();
    for(int ipart=0;ipart<nPart_TRKEFF;ipart++){
      if(abs(_mcpart_PDG[(int)_track_trueID[itrk]-1])==int_TRKEFF_mcparticles_PDG[ipart]){
        h_particle_rec_pT[ipart]->Fill(pt,trueeta);
        h_particle_rec_truepT[ipart]->Fill(truept,trueeta);
        h_chargedpart_rec_pT->Fill(pt,trueeta);
        h_chargedpart_rec_truepT->Fill(truept,trueeta);

        h_particle_rec_p[ipart]->Fill(pmom,trueeta);
        h_particle_rec_truep[ipart]->Fill(truepmom,trueeta);
        h_chargedpart_rec_p->Fill(pmom,trueeta);
        h_chargedpart_rec_truep->Fill(truepmom,trueeta);
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
    h_chargedpart_eff_pT->Write("h_chargedpart_eff_pT", TObject::kOverwrite);
    TH2F*  h_chargedpart_eff_truepT = (TH2F*)h_chargedpart_rec_truepT->Clone("h_chargedpart_eff_truepT");
    h_chargedpart_eff_truepT->Divide(h_chargedpart_MC_pT);
    h_chargedpart_eff_truepT->Write("h_chargedpart_eff_truepT", TObject::kOverwrite);
  }
  if(h_chargedpart_MC_p) h_chargedpart_MC_p->Write();
  if(h_chargedpart_rec_p) h_chargedpart_rec_p->Write();
  if(h_chargedpart_rec_truep) h_chargedpart_rec_truep->Write();
  if(h_chargedpart_MC_p && h_chargedpart_rec_p && h_chargedpart_rec_truep) {
    TH2F*  h_chargedpart_eff_p = (TH2F*)h_chargedpart_rec_p->Clone("h_chargedpart_eff_p");
    h_chargedpart_eff_p->Divide(h_chargedpart_MC_p);
    h_chargedpart_eff_p->Write("h_chargedpart_eff_p",TObject::kOverwrite);
    TH2F*  h_chargedpart_eff_truep = (TH2F*)h_chargedpart_rec_truep->Clone("h_chargedpart_eff_truep");
    h_chargedpart_eff_truep->Divide(h_chargedpart_MC_p);
    h_chargedpart_eff_truep->Write("h_chargedpart_eff_truep",TObject::kOverwrite);
  }
  if(h_chargedpart_MC_E) h_chargedpart_MC_E->Write();
  TH2F*  h_particle_eff_pT[nPart_TRKEFF];
  TH2F*  h_particle_eff_truepT[nPart_TRKEFF];
  TH2F*  h_particle_eff_p[nPart_TRKEFF];
  TH2F*  h_particle_eff_truep[nPart_TRKEFF];
  for(int ipart=0;ipart<nPart_TRKEFF;ipart++){
    if(h_particle_MC_pT[ipart]) h_particle_MC_pT[ipart]->Write();
    if(h_particle_rec_pT[ipart]) h_particle_rec_pT[ipart]->Write();
    if(h_particle_rec_truepT[ipart]) h_particle_rec_truepT[ipart]->Write();
    if(h_particle_MC_pT[ipart] && h_particle_rec_pT[ipart] && h_particle_rec_truepT[ipart]) {
      h_particle_eff_pT[ipart] = (TH2F*)h_particle_rec_pT[ipart]->Clone(Form("h_%s_eff_pT",str_TRKEFF_mcparticles[ipart].Data()));
      h_particle_eff_pT[ipart]->Divide(h_particle_MC_pT[ipart]);
      h_particle_eff_pT[ipart]->Write(Form("h_%s_eff_pT",str_TRKEFF_mcparticles[ipart].Data()),TObject::kOverwrite);
      h_particle_eff_truepT[ipart] = (TH2F*)h_particle_rec_truepT[ipart]->Clone(Form("h_%s_eff_truepT",str_TRKEFF_mcparticles[ipart].Data()));
      h_particle_eff_truepT[ipart]->Divide(h_particle_MC_pT[ipart]);
      h_particle_eff_truepT[ipart]->Write(Form("h_%s_eff_truepT",str_TRKEFF_mcparticles[ipart].Data()),TObject::kOverwrite);
    }
    if(h_particle_MC_p[ipart]) h_particle_MC_p[ipart]->Write();
    if(h_particle_rec_p[ipart]) h_particle_rec_p[ipart]->Write();
    if(h_particle_rec_truep[ipart]) h_particle_rec_truep[ipart]->Write();
    if(h_particle_MC_p[ipart] && h_particle_rec_p[ipart] && h_particle_rec_truep[ipart]) {
      h_particle_eff_p[ipart] = (TH2F*)h_particle_rec_p[ipart]->Clone(Form("h_%s_eff_p",str_TRKEFF_mcparticles[ipart].Data()));
      h_particle_eff_p[ipart]->Divide(h_particle_MC_p[ipart]);
      h_particle_eff_p[ipart]->Write(Form("h_%s_eff_p",str_TRKEFF_mcparticles[ipart].Data()),TObject::kOverwrite);
      h_particle_eff_truep[ipart] = (TH2F*)h_particle_rec_truep[ipart]->Clone(Form("h_%s_eff_truep",str_TRKEFF_mcparticles[ipart].Data()));
      h_particle_eff_truep[ipart]->Divide(h_particle_MC_p[ipart]);
      h_particle_eff_truep[ipart]->Write(Form("h_%s_eff_truep",str_TRKEFF_mcparticles[ipart].Data()),TObject::kOverwrite);
    }
    if(h_particle_MC_E[ipart]) h_particle_MC_E[ipart]->Write();
    
  }
  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
