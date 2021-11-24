void CalculateDimensions(){
  Int_t nTowers         = 9040; //8370;
  
  Double_t dLoop  = 0.03;
  Double_t rArch  = 0.05;
  Double_t lTile  = 0.05;
  Double_t lTower = 1.40;
  Double_t tTileS = 0.004;
  Double_t tTileA = 0.016;
  Double_t tTile  = tTileS+tTileA;
  Double_t addWLS = 0.10;
  Double_t dFiber = 0.0005;
  Int_t layersTu  = 10;
  Int_t layers    = lTower/tTile;
  cout << lTower << "\t" << tTile << "\t"<< layers << endl;
  Double_t tlength = 0;
  for (Int_t l = 0; l < layers; l++){
    Double_t clength =  dLoop*TMath::Pi() + rArch/2 +0.25*2*rArch*TMath::Pi() + lTower-rArch-tTileA-tTile*l + addWLS;
    std::cout<< l << "\t"<< clength << std::endl;
    tlength+=clength;
  }
  std::cout << "total length per tower " << tlength << std::endl;
  
  Double_t densityA   = 8.05; //g/cm^3
  Double_t densityTu  = 19.3; //g/cm^3
  Double_t densityS   = 5; //g/cm^3
  Double_t densityWLS = 1.02; //g/cm^3

  Double_t mTocm    = 100;
  
  Double_t vTowerACm3   = lTile*lTile*tTileA*(layers-layersTu)*TMath::Power(mTocm,3); 
  Double_t vTowerAm3    = lTile*lTile*tTileA*(layers-layersTu); 
  Double_t vTowerACm3Tu = lTile*lTile*tTileA*layersTu*TMath::Power(mTocm,3); 
  Double_t vTowerAm3Tu  = lTile*lTile*tTileA*layersTu; 
  Double_t vTowerSCm3   = lTile*lTile*tTileS*layers*TMath::Power(mTocm,3); 
  Double_t vTowerSm3    = lTile*lTile*tTileS*layers; 
  Double_t vTowerWLSm3  = TMath::Pi()*(dFiber*dFiber)*tlength; //m3
  Double_t vTowerWLSCm3 = TMath::Pi()*(dFiber*dFiber*100*100)*tlength*100; //m3

  std::cout << "**************************************" << std::endl;
  std::cout << "LFHCAL" << std::endl;
  std::cout << "**************************************" << std::endl;

  std::cout << "V scint: " << vTowerAm3 << "m3  "<< vTowerSCm3   << "cm3 \t V abs steel: " << vTowerAm3 << "m3  "<< vTowerACm3   << "cm3 \t V abs tugsten: " << vTowerAm3Tu << "m3  "<< vTowerACm3Tu   << "cm3"<< std::endl;
  
  Double_t mTowerA      = vTowerACm3*densityA/1000;
  Double_t mTowerTu     = vTowerACm3Tu*densityTu/1000;
  Double_t mTowerS      = vTowerSCm3*densityS/1000;
  Double_t mTowerWLS    = vTowerWLSCm3*densityWLS/1000;

  Double_t weightLFHCALTower = (mTowerA+mTowerS+mTowerWLS+mTowerTu);
 
  std::cout << "m scint: " << mTowerS << " kg"<< "\t " << "m abs steel:" << mTowerA << " kg  \t m abs tungsten:" << mTowerTu << " kg \t m WLS: " << mTowerWLS << " kg total: "<< weightLFHCALTower  << " kg"<< endl; 
  std::cout << "total m scint: " << mTowerS*nTowers << " kg \t m abs tungsten:" << mTowerA*nTowers << " kg \t m abs tungsten:" << mTowerTu*nTowers << " kg\t m WLS:" << mTowerWLS*nTowers << " kg \t total:" << (mTowerA+mTowerS+mTowerWLS+mTowerTu)*nTowers << " kg"<< endl; 
  Double_t weightLFHCAL = (mTowerA+mTowerS+mTowerWLS+mTowerTu)*nTowers;
  
  Double_t pSteel       = 3.48 ;//($/kg)
  Double_t pSteelTower  = mTowerA*pSteel;
  
  cout<< "price steel/tower: " << pSteelTower << "\t total: " << pSteelTower*nTowers/1e3 << "K $" << endl;
  
  Double_t pWLS         = 2.57; //$/m
  Double_t pWLSTower    = pWLS*tlength; //$/m
  cout<< "price WLS/tower: " << pWLSTower << "\t total: " << pWLSTower*nTowers/1e3 << "K $"<< endl;
  
  
  
  cout << "total number of tiles/ tower: " <<  layers << endl;
  cout << "total number of tiles: " << layers*nTowers << endl;

  Int_t towersInner = 768;
  Int_t towersOuter = 3824;
  Int_t layersFEMC  = 66;
  Double_t zFEMC    = 0.375;
  Double_t tTilePb  = 0.0016;
  Double_t lAddFEMC = 0.05;
  Double_t lengthFEMCinner = 4 * 25 * (zFEMC+lAddFEMC) ;
  Double_t lengthFEMCouter = 4 * 9 * (zFEMC+lAddFEMC) ;
  Double_t dFiberFEMC = 0.0002;

  Double_t densityPb   = 11.34; //g/cm^3

  Double_t vTowerACm3FEMC = lTile*lTile*tTilePb*layersFEMC*TMath::Power(mTocm,3); 
  Double_t vTowerAm3FEMC  = lTile*lTile*tTilePb*layersFEMC; 
  Double_t vTowerSCm3FEMC = lTile*lTile*tTileS*layersFEMC*TMath::Power(mTocm,3); 
  Double_t vTowerSm3FEMC  = lTile*lTile*tTileS*layersFEMC; 
  Double_t vTowerWLSm3FEMCi  = TMath::Pi()*(dFiberFEMC*dFiberFEMC)*lengthFEMCinner; //m3
  Double_t vTowerWLSCm3FEMCi = TMath::Pi()*(dFiberFEMC*dFiberFEMC*100*100)*lengthFEMCinner*100; //m3
  Double_t vTowerWLSm3FEMCo  = TMath::Pi()*(dFiberFEMC*dFiberFEMC)*lengthFEMCouter; //m3
  Double_t vTowerWLSCm3FEMCo = TMath::Pi()*(dFiberFEMC*dFiberFEMC*100*100)*lengthFEMCouter*100; //m3
  
  std::cout << "**************************************" << std::endl;
  std::cout << "FEMC" << std::endl;
  std::cout << "**************************************" << std::endl;
  std::cout << "V scint: " << vTowerAm3FEMC << "m3  "<< vTowerSCm3FEMC   << "cm3 \t V abs: " << vTowerAm3FEMC << "m3  "<< vTowerACm3FEMC   << "cm3"<< std::endl;
  
  Double_t mTowerAFEMC      = vTowerACm3FEMC*densityPb/1000;
  Double_t mTowerSFEMC      = vTowerSCm3FEMC*densityS/1000;
  Double_t mTowerWLSFEMCi   = vTowerWLSCm3FEMCi*densityWLS/1000;
  Double_t mTowerWLSFEMCo   = vTowerWLSCm3FEMCo*densityWLS/1000;
  Double_t weightFEMC = (mTowerAFEMC+mTowerSFEMC+mTowerWLSFEMCi)*towersInner + (mTowerAFEMC+mTowerSFEMC+mTowerWLSFEMCo)*towersOuter;
  
  std::cout << "m scint: " << mTowerSFEMC << " kg"<< "\t " << "m abs:" << mTowerAFEMC << " kg \t m WLS: " << mTowerWLSFEMCi << "\t" << mTowerWLSFEMCo << endl; 
  std::cout << "total m scint: " << mTowerSFEMC*(towersInner+towersOuter) << " kg \t m abs:" << mTowerAFEMC*(towersInner+towersOuter) << " kg \t m WLS:" << mTowerWLSFEMCi*towersInner << " kg \t" << mTowerWLSFEMCo*towersOuter << " kg \ttotal:" << weightFEMC << " kg"<< endl; 
 
  std::cout << "**************************************" << std::endl;
  
  std::cout <<  "LFHCAL: " << weightLFHCAL/1000 << " t \t FEMC: " <<  weightFEMC/1000 << " t \t forward cal: " << (weightFEMC+weightLFHCAL)/1000 << " t" << std::endl;
  

  
}
