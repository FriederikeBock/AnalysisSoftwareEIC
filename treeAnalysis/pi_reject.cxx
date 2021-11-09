


// pion rejection plots:
// Energy bin: 0, 0.5, 1, 2, 3, 5, 7, 10, 14, 20 => 9 bins { 0 < E <= 0.5, 0.5 < E <= 1. .........} 


const int energy_bin = 9;
const int N_EMCal = 3; // 0: EEMC, 1: BECAL, 2: FEMC

double e_range[10] = {0., 0.5, 1., 2., 3., 5., 7., 10., 14., 20.};
double e_bin[energy_bin] = {0.25, 0.75, 1.5, 2.5, 4., 6., 8.5, 12., 17.};
double e_bin_content_truthp[energy_bin][N_EMCal] = {}, e_bin_content_matchp[energy_bin][N_EMCal] = {};
double pion_e_like_truthp[energy_bin][N_EMCal] = {}, pion_e_like_matchp[energy_bin][N_EMCal] = {};
double rejection_factor_truthp[energy_bin] = {}, rejection_factor_matchp[energy_bin] = {};
double rej_erry_truthp[energy_bin] = {}, rej_erry_matchp[energy_bin] = {};
double rej_errx_truthp[energy_bin] = {}, rej_errx_matchp[energy_bin] = {};

double ep_electron[N_EMCal] = {};
double ep_pion[N_EMCal] = {};
double correction_EP[N_EMCal] = {0.897, 0.865, 1.};

int total_count[N_EMCal] = {};
int clus_2_count = 0, clus_3_count = 0;

TGraphErrors *pi_reject_ratio_truth_mom[N_EMCal], *pi_reject_ratio_match_mom[N_EMCal];

TH1F *h_E_over_truthP_per_bin[energy_bin][N_EMCal], *h_E_over_matchP_per_bin[energy_bin][N_EMCal];
TH1F *h_diff_mom_match_primary[N_EMCal][4];




void pi_reject()
{
  // CheckBarrelCaloVersion();

  // Loop all the track to find the primary particle in a single event
  //
  TLorentzVector PiGen;
  //  cout << endl << endl;
  //  cout << "Primary charge pion infomation: ";
  for( int nMC = 0 ; nMC < _nMCPart ; nMC++ )
    if( (_mcpart_PDG[nMC] == 211) && (_mcpart_ID[nMC] == 1 ) && ( _mcpart_ID_parent[nMC] == 0 ) )
      PiGen.SetPxPyPzE(_mcpart_px[nMC], _mcpart_py[nMC], _mcpart_pz[nMC], _mcpart_E[nMC]);
    //    cout << "Eta, Phi, P, E: [" << PiGen.Eta() << ", " << PiGen.Phi() << ", " << PiGen.P() << ", " << PiGen.E() << "]" << endl;
    //    cout << "Eta, Phi, E: [" << PiGen.Eta() << ", " << PiGen.Phi() << ", " << PiGen.E() << "]" << endl << endl;



  // Declare the histograms
  //
  for( int h_int = 0 ; h_int < energy_bin ; h_int++ )
  {
    for( int h_calo = 0 ; h_calo < N_EMCal ; h_calo++ )
    {
      if( !h_E_over_truthP_per_bin[h_int][h_calo] )
        h_E_over_truthP_per_bin[h_int][h_calo] = new TH1F(Form("h_E_over_truthP_per_bin_%d_%d", h_int, h_calo), "", 105, 0., 1.05);

      if( !h_E_over_matchP_per_bin[h_int][h_calo] )
        h_E_over_matchP_per_bin[h_int][h_calo] = new TH1F(Form("h_E_over_matchP_per_bin_%d_%d", h_int, h_calo), "", 105, 0., 1.05);
    }
  }

  for( int h_calo = 0 ; h_calo < N_EMCal ; h_calo++ )
      for( int h_p3 = 0 ; h_p3 < 4 ; h_p3++ )
        if( !h_diff_mom_match_primary[h_calo][h_p3] )
          h_diff_mom_match_primary[h_calo][h_p3] = new TH1F(Form("h_diff_mom_match_primary_%d_%d", h_calo, h_p3), "", 400, -1., 1.);

  // EP cut coefficient
  //
  const double EP_cut_strength = 1.6;

  // Looping all the calorimeter and the algorithm
  //
  for(int icalo = 0 ; icalo < _active_calo ; icalo++)
  {
    //===== Loop only over enabled calorimeters ===========//
    if (!caloEnabled[icalo]) continue;
    for(int ialgo = 0 ; ialgo < _active_algo ; ialgo++)
    {
      //if(verbosityERH) cout << "Calo-type: " << icalo << "\t Algorithm: " << ialgo << endl;

      if( !loadClusterizerInput(ialgo, icalo) ) continue;
      //    std::sort(_clusters_calo[ialgo][icalo].begin(), _clusters_calo[ialgo][icalo].end(), &acompareCl);

      // Pion rejection Analysis: EEMC, BECAL, FEMC
      // Only study with 1, 2, 3 cluster reconstructed cases
      //
      if( (str_calorimeter[icalo] == "EEMC") || (str_calorimeter[icalo] == "BECAL") || (str_calorimeter[icalo] == "FEMC") )
      {
        int clus_trk_ID = 0, match_trk_N = 0, N_cluster_EMC = _clusters_calo[ialgo][icalo].size();
        double E_over_truthP = 0., E_over_matchP = 0., EP_cut = 0., clus_E = 0.;
        TVector3 trk_vec3;

        // Single cluster with single track match case
        //
        if ( (_clusters_calo[ialgo][icalo].at(0)).cluster_isMatched == true && (N_cluster_EMC == 1) )
        {
          match_trk_N = ((_clusters_calo[ialgo][icalo].at(0)).cluster_matchedTrackIDs).size();

          if( match_trk_N == 1 )
          {
            clus_trk_ID = ((_clusters_calo[ialgo][icalo].at(0)).cluster_matchedTrackIDs).at(0);
            trk_vec3.SetXYZ(_track_px[clus_trk_ID], _track_py[clus_trk_ID], _track_pz[clus_trk_ID]);

            // Skip events with unreasonable energy reconstructed, I am not sure either the cut or the 1.05 is reasonable or not 
            //
            if( ((_clusters_calo[ialgo][icalo].at(0)).cluster_E) / trk_vec3.Mag() > 1.05 ) continue;

            // Find the corresponding calo, energy bin, and fill the data 
            //
            if( str_calorimeter[icalo] == "EEMC" ) 
            {
              //Calculate the E/P and the EP cut event by event
              clus_E = (_clusters_calo[ialgo][icalo].at(0)).cluster_E;
              E_over_truthP = clus_E / PiGen.P();
              E_over_matchP = clus_E / trk_vec3.Mag();

              EP_cut = 1. - EP_cut_strength * std::sqrt( 0.5 * 0.5 + 2.2 * 2.2 / PiGen.E() ) / 100.; //E_reso from Friedrike

              h_diff_mom_match_primary[0][0]->Fill( PiGen.P() - trk_vec3.Mag() );
              h_diff_mom_match_primary[0][1]->Fill( PiGen.Px() - trk_vec3.X() );
              h_diff_mom_match_primary[0][2]->Fill( PiGen.Py() - trk_vec3.Y() );
              h_diff_mom_match_primary[0][3]->Fill( PiGen.Pz() - trk_vec3.Z() );

              for( int i_e_bin = 0 ; i_e_bin < energy_bin ; i_e_bin++ )
              {
                if( e_range[i_e_bin] < PiGen.E() && e_range[i_e_bin+1] >= PiGen.E() )
                {
                  e_bin_content_truthp[i_e_bin][0]++;
                  e_bin_content_matchp[i_e_bin][0]++;
                  h_E_over_truthP_per_bin[i_e_bin][0]->Fill(E_over_truthP);
                  h_E_over_matchP_per_bin[i_e_bin][0]->Fill(E_over_matchP);

                  if( E_over_truthP > EP_cut )
                    pion_e_like_truthp[i_e_bin][0]++;
                  if( E_over_matchP > EP_cut )
                    pion_e_like_matchp[i_e_bin][0]++;
                }
              }
              total_count[0]++;
            }
            else if( str_calorimeter[icalo] == "BECAL" ) // BECAL
            {
              //Calculate the E/P and the EP cut event by event
              clus_E = (_clusters_calo[ialgo][icalo].at(0)).cluster_E;
              E_over_truthP = clus_E / PiGen.P();
              E_over_matchP = clus_E / trk_vec3.Mag();

              EP_cut = 1. - EP_cut_strength * std::sqrt( 0.8 * 0.8 + 1.7 * 1.7 / PiGen.E() ) / 100.; //E_reso from Friedrike

              h_diff_mom_match_primary[1][0]->Fill( PiGen.P() - trk_vec3.Mag() );
              h_diff_mom_match_primary[1][1]->Fill( PiGen.Px() - trk_vec3.X() );
              h_diff_mom_match_primary[1][2]->Fill( PiGen.Py() - trk_vec3.Y() );
              h_diff_mom_match_primary[1][3]->Fill( PiGen.Pz() - trk_vec3.Z() );

              for( int i_e_bin = 0 ; i_e_bin < energy_bin ; i_e_bin++ )
                {
                  if( e_range[i_e_bin] < PiGen.E() && e_range[i_e_bin+1] >= PiGen.E() )
              {
                e_bin_content_truthp[i_e_bin][1]++;
                e_bin_content_matchp[i_e_bin][1]++;
                h_E_over_truthP_per_bin[i_e_bin][1]->Fill(E_over_truthP);
                h_E_over_matchP_per_bin[i_e_bin][1]->Fill(E_over_matchP);

                if( E_over_truthP > EP_cut )
                  pion_e_like_truthp[i_e_bin][1]++;
                if( E_over_matchP > EP_cut )
                  pion_e_like_matchp[i_e_bin][1]++;
              }
                }
              total_count[1]++;
            }
              else if( str_calorimeter[icalo] == "FEMC" ) 
            {
              //Calculate the E/P and the EP cut event by event
              clus_E = (_clusters_calo[ialgo][icalo].at(0)).cluster_E;
              E_over_truthP = clus_E / PiGen.P();
              E_over_matchP = clus_E / trk_vec3.Mag();

              EP_cut = 1. - EP_cut_strength * std::sqrt( 0.1 * 0.1 + 7.3 * 7.3 / PiGen.E() ) / 100.; //E_reso from Friedrike

              h_diff_mom_match_primary[2][0]->Fill( PiGen.P() - trk_vec3.Mag() );
              h_diff_mom_match_primary[2][1]->Fill( PiGen.Px() - trk_vec3.X() );
              h_diff_mom_match_primary[2][2]->Fill( PiGen.Py() - trk_vec3.Y() );
              h_diff_mom_match_primary[2][3]->Fill( PiGen.Pz() - trk_vec3.Z() );

              for( int i_e_bin = 0 ; i_e_bin < energy_bin ; i_e_bin++ )
              {
                if( e_range[i_e_bin] < PiGen.E() && e_range[i_e_bin+1] >= PiGen.E() )
                {
                  e_bin_content_truthp[i_e_bin][2]++;
                  e_bin_content_matchp[i_e_bin][2]++;
                  h_E_over_truthP_per_bin[i_e_bin][2]->Fill(E_over_truthP);
                  h_E_over_matchP_per_bin[i_e_bin][2]->Fill(E_over_matchP);

                  if( E_over_truthP > EP_cut )
                    pion_e_like_truthp[i_e_bin][2]++;
                  if( E_over_matchP > EP_cut )
                    pion_e_like_matchp[i_e_bin][2]++;
                }
              }
              total_count[2]++;
            }
          }
        }// 1 cluster case done


          
        int N_clus_has_match = 0, N_trk_clus = 0;
        
        // Multiple clusters case
        //
        if( (N_cluster_EMC == 2) || (N_cluster_EMC == 3) )
        {

          for( int clus_idx = 0 ; clus_idx < N_cluster_EMC ; clus_idx++ )
          {
            if( (_clusters_calo[ialgo][icalo].at(clus_idx)).cluster_isMatched == true )
            {
              N_clus_has_match++;
              N_trk_clus = N_trk_clus + ((_clusters_calo[ialgo][icalo].at(clus_idx)).cluster_matchedTrackIDs).size();
            }
          }


          // Only focus on 1 trk matched to 1 cluster case
          if( (N_clus_has_match != 1) && (N_trk_clus != 1) )
            continue;

          for( int clus_idx = 0 ; clus_idx < N_cluster_EMC ; clus_idx++ )
          {
              // deal with the cluster has track match
              //
            if( (_clusters_calo[ialgo][icalo].at(clus_idx)).cluster_isMatched == true )
            {
              clus_trk_ID = ((_clusters_calo[ialgo][icalo].at(clus_idx)).cluster_matchedTrackIDs).at(0);
              trk_vec3.SetXYZ(_track_px[clus_trk_ID], _track_py[clus_trk_ID], _track_pz[clus_trk_ID]);

              // Again, Skip events with unreasonable energy reconstructed, I am not sure either the cut or the 1.05 is reasonable or not 
              //
              if( (_clusters_calo[ialgo][icalo].at(clus_idx)).cluster_E / trk_vec3.Mag() > 1.05 ) continue;


              // Find the corresponding energy bin, and fill the data 
              //
              if( str_calorimeter[icalo] == "EEMC" )
              {
                //Calculate the E/P and the EP cut event by event        
                clus_E = (_clusters_calo[ialgo][icalo].at(0)).cluster_E;
                E_over_truthP = clus_E / PiGen.P();
                E_over_matchP = clus_E / trk_vec3.Mag();

                EP_cut = 1. - EP_cut_strength * std::sqrt( 0.5 * 0.5 + 2.2 * 2.2 / PiGen.E() ) / 100.; //E_reso from Friedrike
                
                h_diff_mom_match_primary[0][0]->Fill( PiGen.P() - trk_vec3.Mag() );
                h_diff_mom_match_primary[0][1]->Fill( PiGen.Px() - trk_vec3.X() );
                h_diff_mom_match_primary[0][2]->Fill( PiGen.Py() - trk_vec3.Y() );
                h_diff_mom_match_primary[0][3]->Fill( PiGen.Pz() - trk_vec3.Z() );

                for( int i_e_bin = 0 ; i_e_bin < energy_bin ; i_e_bin++ )
                {
                  if( e_range[i_e_bin] < PiGen.E() && e_range[i_e_bin+1] >= PiGen.E() )
                  {
                    e_bin_content_truthp[i_e_bin][0]++;
                    e_bin_content_matchp[i_e_bin][0]++;
                    h_E_over_truthP_per_bin[i_e_bin][0]->Fill(E_over_truthP);
                    h_E_over_matchP_per_bin[i_e_bin][0]->Fill(E_over_matchP);

                    if( E_over_truthP > EP_cut )
                      pion_e_like_truthp[i_e_bin][0]++;
                    if( E_over_matchP > EP_cut )
                      pion_e_like_matchp[i_e_bin][0]++;
                  }
                }
                  total_count[0]++;
              }
              else if( str_calorimeter[icalo] == "BECAL" ) 
              {
                //Calculate the E/P and the EP cut event by event        
                clus_E = (_clusters_calo[ialgo][icalo].at(0)).cluster_E;
                E_over_truthP = clus_E / PiGen.P();
                E_over_matchP = clus_E / trk_vec3.Mag();

                EP_cut = 1. - EP_cut_strength * std::sqrt( 0.8 * 0.8 + 1.7 * 1.7 / PiGen.E() ) / 100.; //E_reso from Friedrike
                
                h_diff_mom_match_primary[1][0]->Fill( PiGen.P() - trk_vec3.Mag() );
                h_diff_mom_match_primary[1][1]->Fill( PiGen.Px() - trk_vec3.X() );
                h_diff_mom_match_primary[1][2]->Fill( PiGen.Py() - trk_vec3.Y() );
                h_diff_mom_match_primary[1][3]->Fill( PiGen.Pz() - trk_vec3.Z() );

                for( int i_e_bin = 0 ; i_e_bin < energy_bin ; i_e_bin++ )
                {
                  if( e_range[i_e_bin] < PiGen.E() && e_range[i_e_bin+1] >= PiGen.E() )
                  {
                    e_bin_content_truthp[i_e_bin][1]++;
                    e_bin_content_matchp[i_e_bin][1]++;
                    h_E_over_truthP_per_bin[i_e_bin][1]->Fill(E_over_truthP);
                    h_E_over_matchP_per_bin[i_e_bin][1]->Fill(E_over_matchP);
                
                    if( E_over_truthP > EP_cut )
                      pion_e_like_truthp[i_e_bin][1]++;
                    if( E_over_matchP > EP_cut )
                      pion_e_like_matchp[i_e_bin][1]++;
                  }
                }
                total_count[1]++;
              }
              else if( str_calorimeter[icalo] == "FEMC" ) 
              {
                //Calculate the E/P and the EP cut event by event        
                clus_E = (_clusters_calo[ialgo][icalo].at(0)).cluster_E;
                E_over_truthP = clus_E / PiGen.P();
                E_over_matchP = clus_E / trk_vec3.Mag();
                
                EP_cut = 1. - EP_cut_strength * std::sqrt( 0.1 * 0.1 + 7.3 * 7.3 / PiGen.E() ) / 100.; //E_reso from Friedrike

                h_diff_mom_match_primary[2][0]->Fill( PiGen.P() - trk_vec3.Mag() );
                h_diff_mom_match_primary[2][1]->Fill( PiGen.Px() - trk_vec3.X() );
                h_diff_mom_match_primary[2][2]->Fill( PiGen.Py() - trk_vec3.Y() );
                h_diff_mom_match_primary[2][3]->Fill( PiGen.Pz() - trk_vec3.Z() );

                for( int i_e_bin = 0 ; i_e_bin < energy_bin ; i_e_bin++ )
                {
                  if( e_range[i_e_bin] < PiGen.E() && e_range[i_e_bin+1] >= PiGen.E() )
                  {
                    e_bin_content_truthp[i_e_bin][2]++;
                    e_bin_content_matchp[i_e_bin][2]++;
                    h_E_over_truthP_per_bin[i_e_bin][2]->Fill(E_over_truthP);
                    h_E_over_matchP_per_bin[i_e_bin][2]->Fill(E_over_matchP);
                
                    if( E_over_truthP > EP_cut )
                      pion_e_like_truthp[i_e_bin][2]++;

                    if( E_over_matchP > EP_cut )
                      pion_e_like_matchp[i_e_bin][2]++;
                  }
                }
                total_count[2]++;
              }
            }// Find the cluster has the track match 
          }// loop multiple clusters
        }// multiple clusters case done !!!
      }// EEMC BECAL FEMC pi reject study done !!!
    }// loop algorithm      
  }// loop caloirmeter
} // End function




// ANCHOR save function after event loop
//
void pi_rejectSave()
{ 
  
  // define output file
  //
  TFile* fileOutput = new TFile(Form("%s/output_pi_reject.root", outputDir.Data()), "RECREATE");

  for( int hcal = 0 ; hcal < N_EMCal ; hcal++ )
  {
    if( hcal == 0 )
      cout << "EEMC pion rejection study: " << endl;
    else if( hcal == 1 )
      cout << "BECAL pion rejection study: " << endl;
    else if( hcal == 2 )
      cout << "FEMC pion rejection study: " << endl;

    cout << "We have total " << total_count[hcal] << " have the good track match." << endl;


    // All the information from the truth momentum
    //
    cout << "====================================== Truth ======================================" << endl;
    for( int i_e_bin = 0 ; i_e_bin < energy_bin ; i_e_bin++ )
    {
      cout << e_bin[i_e_bin] << "[GeV] has " << e_bin_content_truthp[i_e_bin][hcal] << " good counts";
      cout << " and " << pion_e_like_truthp[i_e_bin][hcal] << " are electron like !!!";

      if( e_bin_content_truthp[i_e_bin][hcal] > 0 && pion_e_like_truthp[i_e_bin][hcal] > 0 )
      {
        rejection_factor_truthp[i_e_bin] = e_bin_content_truthp[i_e_bin][hcal] / pion_e_like_truthp[i_e_bin][hcal];
        rej_erry_truthp[i_e_bin] = e_bin_content_truthp[i_e_bin][hcal] / pion_e_like_truthp[i_e_bin][hcal] * std::sqrt(pion_e_like_truthp[i_e_bin][hcal]) / pion_e_like_truthp[i_e_bin][hcal];

        cout << " Rejection factor: " << rejection_factor_truthp[i_e_bin] << " and the error: " << rej_erry_truthp[i_e_bin] << endl;
      }
      else
      {
        rejection_factor_truthp[i_e_bin] = 0.;
        rej_erry_truthp[i_e_bin] = 0;

        cout << " Rejection factor: " << rejection_factor_truthp[i_e_bin] << " and the error: 0" << endl;
      }

      if(h_E_over_truthP_per_bin[i_e_bin][hcal])
        h_E_over_truthP_per_bin[i_e_bin][hcal]->Write();
    }

    pi_reject_ratio_truth_mom[hcal] = new TGraphErrors(energy_bin, e_bin, rejection_factor_truthp, rej_errx_truthp, rej_erry_truthp);        
    pi_reject_ratio_truth_mom[hcal]->SetName(Form("TGE_rejection_factor_truthp_%d", hcal));
    if(pi_reject_ratio_truth_mom[hcal]) pi_reject_ratio_truth_mom[hcal]->Write();

    cout << "===================================================================================" << endl;



    cout << "====================================== Match ======================================" << endl;
    // All the information from the match momentum
    //
    for( int i_e_bin = 0 ; i_e_bin < energy_bin ; i_e_bin++ )
    {
      cout << e_bin[i_e_bin] << "[GeV] has " << e_bin_content_matchp[i_e_bin][hcal] << " good counts";
      cout << " and " << pion_e_like_matchp[i_e_bin][hcal] << " are electron like !!!";

      if( e_bin_content_matchp[i_e_bin][hcal] > 0 && pion_e_like_matchp[i_e_bin][hcal] > 0 )
      {
        rejection_factor_matchp[i_e_bin] = e_bin_content_matchp[i_e_bin][hcal] / pion_e_like_matchp[i_e_bin][hcal];
        rej_erry_matchp[i_e_bin] = e_bin_content_matchp[i_e_bin][hcal] / pion_e_like_matchp[i_e_bin][hcal] * std::sqrt(pion_e_like_matchp[i_e_bin][hcal]) / pion_e_like_matchp[i_e_bin][hcal];

        cout << " Rejection factor: " << rejection_factor_matchp[i_e_bin] << " and the error: " << rej_erry_matchp[i_e_bin] << endl;
      }
      else
      {
        rejection_factor_matchp[i_e_bin] = 0.;
        rej_erry_matchp[i_e_bin] = 0;

        cout << " Rejection factor: " << rejection_factor_matchp[i_e_bin] << " and the error: 0" << endl;
      }

      if(h_E_over_matchP_per_bin[i_e_bin][hcal])
        h_E_over_matchP_per_bin[i_e_bin][hcal]->Write();
    }

      pi_reject_ratio_match_mom[hcal] = new TGraphErrors(energy_bin, e_bin, rejection_factor_matchp, rej_errx_matchp, rej_erry_matchp);        
      pi_reject_ratio_match_mom[hcal]->SetName(Form("TGE_rejection_factor_matchp_%d", hcal));
      if(pi_reject_ratio_match_mom[hcal]) pi_reject_ratio_match_mom[hcal]->Write();

      cout << "===================================================================================" << endl;

      for(int p4 = 0 ; p4 < 4 ; p4++)
        if(h_diff_mom_match_primary[hcal][p4])
          h_diff_mom_match_primary[hcal][p4]->Write();

      cout << endl << endl << endl;
  }

  //  write output file
  fileOutput->Write();
  fileOutput->Close();

}
