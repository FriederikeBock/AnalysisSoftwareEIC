#include <algorithm>
// ANCHOR debug output verbosity
int verbosityTM = 1;

// ANCHOR main function to be called in event loop
void trackmatcher(){

}

// ANCHOR main function to be called in event loop
bool isClusterMatched(int clsID, string calo){
  int layerselect = -1;
  if(calo.find("FHCAL") != std::string::npos){
    layerselect = 5;
  } else if(calo.find("FEMC") != std::string::npos){
    layerselect = 6;
  } else {
    cout << "selected incorrect calorimeter: " << calo << endl;
    return false;
  }


  for(Int_t iproj=0; iproj<_nProjections; iproj++){
    if(_track_ProjLayer[iproj]!=layerselect) continue;
    if(verbosity>1) cout << "\tTrack: track " << iproj << "\twith true ID " << _track_Proj_x[iproj] << "\tand X = " << _track_Proj_y[iproj] << " cm" << endl;
    // check phi difference
    // if()
  }
  return false;
}


// int _nProjections;
// float* _track_ProjTrackID                 = new float[_maxNProjections];
// int* _track_ProjLayer                 = new int[_maxNProjections];
// float* _track_Proj_x             = new float[_maxNProjections];
// float* _track_Proj_y             = new float[_maxNProjections];
// float* _track_Proj_z             = new float[_maxNProjections];
// float* _track_Proj_t             = new float[_maxNProjections];
// float* _track_Proj_true_x             = new float[_maxNProjections];
// float* _track_Proj_true_y             = new float[_maxNProjections];
// float* _track_Proj_true_z             = new float[_maxNProjections];
// float* _track_Proj_true_t             = new float[_maxNProjections];


// ANCHOR save function after event loop
void saveHistosV3clusterizer(){
  // make output directory
  gSystem->Exec("mkdir -p treeProcessing/V3clusterizer");
  // define output file
  TFile* fileOutput = new TFile("treeProcessing/V3clusterizer/output_V3.root","RECREATE");

  // write histograms (e.g. cluster spectrum, M02 distribution, etc)
  // h_pion_fhcal_E->Write();

  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
