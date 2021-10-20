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
  
   for(int imc=0; imc<_nMCPart; imc++){
    TVector3 mcpartvec(_mcpart_px[imc],_mcpart_py[imc],_mcpart_pz[imc]);
    float trueeta = _mcpart_Eta[imc];
    float truept  = mcpartvec.Pt();
    float truep   = mcpartvec.Mag();

    if (verbosityTRKEFF) cout << "\tMCID " << imc << "\tPDG: " << _mcpart_PDG[imc] << "\tpT: " << TMath::Sqrt(TMath::Power(_mcpart_px[imc],2)+TMath::Power(_mcpart_py[imc],2)) << "\tEta " << trueeta << endl;
    for(int ipart=0;ipart<nPart_TRKEFF;ipart++){
      // MC
      if(!h_particle_MC_pT[ipart]) h_particle_MC_pT[ipart] 	= new TH2F(Form("h_%s_MC_pT",str_TRKEFF_mcparticles[ipart].Data()), "", 60, 0,30,80, -4.0,4.0);
      if(!h_particle_MC_p[ipart]) h_particle_MC_p[ipart] 	= new TH2F(Form("h_%s_MC_p",str_TRKEFF_mcparticles[ipart].Data()), "", nBinsP, binningP, 80, -4.0,4.0);
      if(!h_particle_MC_E[ipart]) h_particle_MC_E[ipart] 	= new TH2F(Form("h_%s_MC_E",str_TRKEFF_mcparticles[ipart].Data()), "", nBinsP, binningP, 80, -4.0,4.0);
      // Rec
 
  
}
