// TOF pid analysis module

//----------------------------------------------------------------------------
// define histos

//----------------------------------------------------------------------------
bool bHistTOFCreated    = false;
TH2F* h_TPH_InvGenBeta_Eta_p[4][nPID][nEta+1]        = {{NULL}};
TH2F* h_TPH_InvBeta_Eta_p[4][nPID][nEta+1]           = {{NULL}};
TH2F* h_TPH_InvSmearBeta_Eta_p[4][nPID][nEta+1]      = {{NULL}};
TH2F* h_TPH_Res_InvBeta_Eta_p[4][nPID][nEta+1]       = {{NULL}};
TH2F* h_TPH_Res_InvSmearBeta_Eta_p[4][nPID][nEta+1]  = {{NULL}};

// random number generator
TRandom3 r3;
// define c
const double c                  = 29.9792458; // cm/ns
// set timing resolution
double sigmat             = 20e-3; // ns


void tofpidhistos(){
 
  if (!bHistTOFCreated){
    for (Int_t et = 0; et<nEta+1; et++){
      Double_t etaMin = partEta[0];
      Double_t etaMax = partEta[nEta];
      if (et < nEta){
        etaMin = partEta[et];
        etaMax = partEta[et+1];
      }      
      for (Int_t e = 0; e< 4; e++){
        for (Int_t id = 0; id < nPID; id++){
          h_InvGenBeta_Eta_p[e][id][et]      = new TH2F(Form("h_InvGenBeta%s_Eta_p_%s_%d", nameAddBeta[e].Data(), partName[id].Data(), et), 
                                                  Form("%1.1f < #eta < %1.1f; #it{p}_{MC} (GeV/c); 1/#beta_{MC}", etaMin, etaMax), 
                                                  nBinsPLow, binningP, 300, 0.9, 1.5);

          h_InvBeta_Eta_p[e][id][et]         = new TH2F(Form("h_InvBeta%s_Eta_p_%s_%d", nameAddBeta[e].Data(), partName[id].Data(), et), 
                                                  Form("%1.1f < #eta < %1.1f; #it{p}_{MC} (GeV/c); 1/#beta_{rec}", etaMin, etaMax), 
                                                  nBinsPLow, binningP, 300, 0.9, 1.5);
          h_InvSmearBeta_Eta_p[e][id][et]    = new TH2F(Form("h_InvSmearBeta%s_Eta_p_%s_%d", nameAddBeta[e].Data(),partName[id].Data(), et), 
                                                  Form("%1.1f < #eta < %1.1f; #it{p}_{MC} (GeV/c); 1/#beta_{smear}", etaMin, etaMax), 
                                                  nBinsPLow, binningP, 300, 0.9, 1.5);
          h_Res_InvBeta_Eta_p[e][id][et]     = new TH2F(Form("h_Res_InvBeta%s_Eta_p_%s_%d", nameAddBeta[e].Data(),partName[id].Data(), et), 
                                                  Form("%1.1f < #eta < %1.1f; #it{p}_{MC} (GeV/c); (1/#beta_{rec} - 1/#beta_{MC})/(1/#beta_{MC}) ", etaMin, etaMax), 
                                                  nBinsPLow, binningP, 500, -0.1, 0.1);
          h_Res_InvSmearBeta_Eta_p[e][id][et]= new TH2F(Form("h_Res_InvSmearBeta%s_Eta_p_%s_%d", nameAddBeta[e].Data(),partName[id].Data(), et), 
                                                  Form("%1.1f < #eta < %1.1f; #it{p}_{MC} (GeV/c); (1/#beta_{smear} - 1/#beta_{MC})/(1/#beta_{MC}) ", etaMin, etaMax), 
                                                  nBinsPLow, binningP, 500, -0.1, 0.1);
        }
      }
    }
    bHistTOFCreated = true;
  }
  
    
  for(Int_t itrk=0; itrk<(Int_t)_nTracks; itrk++){
    unsigned short trackSource = _track_source[itrk];
    TVector3 recpartvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
    TVector3 mcpartvec(_mcpart_px[(int)_track_trueID[itrk]],_mcpart_py[(int)_track_trueID[itrk]],_mcpart_pz[(int)_track_trueID[itrk]]);
    // cout << itrk << endl;
    float receta = recpartvec.Eta();
    float trueeta = mcpartvec.Eta();
    if (verbosityTRKEFF) cout << "\tTRKTRUEID " << (int)_track_trueID[itrk] << "\tPDG: " << _mcpart_PDG[(int)_track_trueID[itrk]] << "\tpT: " <<   TMath::Sqrt(TMath::Power(_track_px[itrk],2)+TMath::Power(_track_py[itrk],2)) << "\tEta " << receta << endl;
    float pt = recpartvec.Pt();
    float truept = mcpartvec.Pt();
    float pmom = recpartvec.Mag();
    float truepmom = mcpartvec.Mag();
    TParticlePDG *part = TDatabasePDG::Instance()->GetParticle(abs(_mcpart_PDG[(int)_track_trueID[itrk]]));
    Double_t mass = part->Mass();
    float trueE     = TMath::Sqrt(truepmom*truepmom+mass*mass);
    
  }
  
}
