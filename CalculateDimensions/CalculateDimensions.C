void CalculateDimensions(){
  Int_t nTowers         = 8370;
  
  Double_t dLoop  = 0.03;
  Double_t rArch  = 0.05;
  Double_t lTile  = 0.05;
  Double_t lTower = 1.40;
  Double_t tTileS = 0.004;
  Double_t tTileA = 0.016;
  Double_t tTile  = tTileS+tTileA;
  Double_t addWLS = 0.10;
  Double_t dFiber = 0.0005;
  
  Int_t layers = lTower/tTile;
  cout << lTower << "\t" << tTile << "\t"<< layers << endl;
  Double_t tlength = 0;
  for (Int_t l = 0; l < layers; l++){
    Double_t clength =  dLoop*TMath::Pi() + rArch/2 +0.25*2*rArch*TMath::Pi() + lTower-rArch-tTileA-tTile*l + addWLS;
    std::cout<< l << "\t"<< clength << std::endl;
    tlength+=clength;
  }
  std::cout << "total length per tower " << tlength << std::endl;
  
  Double_t densityA   = 8.05; //g/cm^3
  Double_t densityS   = 5; //g/cm^3
  Double_t densityWLS = 1.02; //g/cm^3

  Double_t mTocm    = 100;
  
  Double_t vTowerACm3   = lTile*lTile*tTileA*layers*TMath::Power(mTocm,3); 
  Double_t vTowerAm3    = lTile*lTile*tTileA*layers; 
  Double_t vTowerSCm3   = lTile*lTile*tTileS*layers*TMath::Power(mTocm,3); 
  Double_t vTowerSm3    = lTile*lTile*tTileS*layers; 
  Double_t vTowerWLSm3  = TMath::Pi()*(dFiber*dFiber)*tlength; //m3
  Double_t vTowerWLSCm3 = TMath::Pi()*(dFiber*dFiber*100*100)*tlength*100; //m3
  
  std::cout << "V scint: " << vTowerAm3 << "m3  "<< vTowerSCm3   << "cm3 \t V abs: " << vTowerAm3 << "m3  "<< vTowerACm3   << "cm3"<< std::endl;
  
  Double_t mTowerA      = vTowerACm3*densityA/1000;
  Double_t mTowerS      = vTowerSCm3*densityS/1000;
  Double_t mTowerWLS    = vTowerWLSCm3*densityWLS/1000;
  
  std::cout << "m scint: " << mTowerS << " kg"<< "\t " << "m abs:" << mTowerA << " kg \t m WLS: " << mTowerWLS << endl; 
  std::cout << "total m scint: " << mTowerS*nTowers << " kg \t m abs:" << mTowerA*nTowers << " kg \t m WLS:" << mTowerWLS*nTowers << " kg \t total:" << (mTowerA+mTowerS+mTowerWLS)*nTowers << " kg"<< endl; 
 
  
  Double_t pSteel       = 3.48 ;//($/kg)
  Double_t pSteelTower  = mTowerA*pSteel;
  
  
  cout<< "price steel/tower: " << pSteelTower << "\t total: " << pSteelTower*nTowers/1e3 << "K $" << endl;
  
  Double_t pWLS         = 2.57; //$/m
  Double_t pWLSTower    = pWLS*tlength; //$/m
  cout<< "price WLS/tower: " << pWLSTower << "\t total: " << pWLSTower*nTowers/1e3 << "K $"<< endl;
  
  
  
  cout << "total number of tiles/ tower: " <<  layers << endl;
  cout << "total number of tiles: " << layers*nTowers << endl;
  
}
