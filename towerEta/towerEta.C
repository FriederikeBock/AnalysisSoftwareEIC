#include "../common/plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void towerEta(
    TString suffix            = "pdf"
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString dateForOutput                       = ReturnDateStringForOutput();

  TString outputDir 						                  = Form("plots/%s",dateForOutput.Data());
  gSystem->Exec("mkdir -p "+outputDir);

  Double_t textSizeSinglePad                    = 0.05;
  Double_t textSizeLabelsPixel                    = 35;
  Double_t textSizeLabelsRel                    = 58./1300;

// 2D PLOT
  TCanvas* cReso = new TCanvas("cReso","",0,0,1100,1000);
  DrawGammaCanvasSettings( cReso, 0.1, 0.1, 0.01, 0.105);
  // cReso->SetLogz();

  TH2F* histoTowerMapDummy   = new TH2F("histoTowerMapDummy","histoTowerMapDummy",276,0, 275.5,276,0, 275.5);
  SetStyleHistoTH2ForGraphs(histoTowerMapDummy, "#it{x} (cm)","#it{y} (cm)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.97);
  histoTowerMapDummy->GetXaxis()->SetNoExponent();
  histoTowerMapDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  histoTowerMapDummy->Draw();

  TEllipse *el4 = new TEllipse(0.0,0.0,250,250,0,90,0);
   el4->SetLineColor(4);
   el4->SetLineWidth(2);
  //  el4->Draw();
  Double_t zmid = 400;

  for(Int_t ix=0;ix<77;ix++){
    if((15+(10*ix))<265){
      DrawGammaLines(0, 5, 15+10*ix, 15+10*ix, 1, kGray+2, 1);
      DrawGammaLines(5, 15, 15+10*ix, 15+10*ix, 1, kGray+2, 1);
      drawLatexAdd(Form("%1.2f",-TMath::Log(TMath::Tan((15+10*ix)/zmid/2))),0.15,0.145+0.032*ix,0.6*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("#Delta#eta = %1.2f",TMath::Abs(-TMath::Log(TMath::Tan((15+10*(ix+1))/zmid/2))-(-TMath::Log(TMath::Tan((15+10*ix)/zmid/2))))),0.25,0.165+0.032*ix,0.6*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    }
  }
  DrawGammaLines(0, 100, 262, 262, 1, kGray+2, 1);
  DrawGammaLines(5, 5, 15, 262, 1, kGray+2, 1);
  DrawGammaLines(15, 15, 0, 262, 1, kGray+2, 1);

  Double_t yin  = 15;
  cout << "asin: " << TMath::ASin(yin/zmid) << endl;
  cout << "eta: " << -TMath::Log(TMath::Tan(yin/zmid/2)) << endl;

  Double_t zposFST[5]={35,53,77,101,125};
  Double_t radiusInFST[5]={4,4.5,5,6,6.5};
  Double_t radiusOutFST[5]={30,35,40,40,43};
  for(Int_t iFST=0;iFST<5;iFST++){
    cout << "FST layer " << iFST << "\tetaIn: " << -TMath::Log(TMath::Tan(radiusInFST[iFST]/zposFST[iFST]/2)) << "\tetaOut: " << -TMath::Log(TMath::Tan(radiusOutFST[iFST]/zposFST[iFST]/2))<< endl;

  }
  Double_t zposLGAD[6]={287,287,289,289,340,340};
  Double_t radiusInLGAD[6]={15,70.1,15,70.1,18,100.1};
  Double_t radiusOutLGAD[6]={70,170,70,170,100,240};
  for(Int_t iLGAD=0;iLGAD<6;iLGAD++){
    cout << "LGAD layer " << iLGAD << "\t" << zposLGAD[iLGAD] << "/" << radiusInLGAD[iLGAD] << "-" << radiusOutLGAD[iLGAD] << "\tetaIn: " << -TMath::Log(TMath::Tan(radiusInLGAD[iLGAD]/zposLGAD[iLGAD]/2)) << "\tetaOut: " << -TMath::Log(TMath::Tan(radiusOutLGAD[iLGAD]/zposLGAD[iLGAD]/2))<< endl;

  }
  Double_t zposLBL[5]={25,49,73,97,121};
  Double_t radiusInLBL[5]={3.18,3.18,3.5,4.7,5.91};
  Double_t radiusOutLBL[5]={18.5,36.26,43.23,43.23,43.23};
  for(Int_t iLBL=0;iLBL<5;iLBL++){
    cout << "LBL layer " << iLBL << "\t" << zposLBL[iLBL] << "/" << radiusInLBL[iLBL] << "-" << radiusOutLBL[iLBL] << "\tetaIn: " << -TMath::Log(TMath::Tan(radiusInLBL[iLBL]/zposLBL[iLBL]/2)) << "\tetaOut: " << -TMath::Log(TMath::Tan(radiusOutLBL[iLBL]/zposLBL[iLBL]/2))<< endl;

  }


  // drawLatexAdd(collisionSystem.Data(),0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  // drawLatexAdd(Form("Anti-#it{k}_{T}, HCAL+ECAL jets"),0.15,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  histoTowerMapDummy->Draw("axis,same");
  cReso->Print(Form("%s/towerGeometry.%s", outputDir.Data(), suffix.Data()));


  TH2F* histoTowerMapDummy2   = new TH2F("histoTowerMapDummy2","histoTowerMapDummy2",266,-132.5, 132.5,266,0, 266.5);
  SetStyleHistoTH2ForGraphs(histoTowerMapDummy2, "#it{x} (cm)","#it{y} (cm)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.97);
  histoTowerMapDummy2->GetXaxis()->SetNoExponent();
  histoTowerMapDummy2->GetXaxis()->SetMoreLogLabels(kTRUE);
  histoTowerMapDummy2->Draw();

  // TEllipse *el4 = new TEllipse(0.0,0.0,250,250,0,90,0);
  //  el4->SetLineColor(4);
  //  el4->SetLineWidth(2);
  // //  el4->Draw();
  zmid = 400;

  for(Int_t ix=0;ix<400;ix++){
    if(15+10*ix >150 && 15+10*ix <265){
      DrawGammaLines(0, -5, 15+10*ix, 15+10*ix, 1, kGray+2, 1);
      DrawGammaLines(-5, -15, 15+10*ix, 15+10*ix, 1, kGray+2, 1);
      drawLatexAdd(Form("%1.2f",-TMath::Log(TMath::Tan((15+10*ix)/zmid/2))),0.40,0.148+0.0345*ix,0.7*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("#Delta#eta = %1.2f",TMath::Abs(-TMath::Log(TMath::Tan((15+10*(ix+1))/zmid/2))-(-TMath::Log(TMath::Tan((15+10*ix)/zmid/2))))),0.25,0.168+0.0345*ix,0.7*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    }
    if(15+5*ix >70 && 15+5*ix <=150){
      DrawGammaLines(0, -5, 15+5*ix, 15+5*ix, 1, kGray+2, 1);
      DrawGammaLines(-5, -15, 15+5*ix, 15+5*ix, 1, kGray+2, 1);
      drawLatexAdd(Form("%1.2f",-TMath::Log(TMath::Tan((15+5*ix)/zmid/2))),0.42,0.148+0.01725*ix,0.35*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("#Delta#eta = %1.2f",TMath::Abs(-TMath::Log(TMath::Tan((15+5*(ix+1))/zmid/2))-(-TMath::Log(TMath::Tan((15+5*ix)/zmid/2))))),0.30,0.168+0.01725*ix,0.35*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    }
    if(15+2.5*ix >15 && 15+2.5*ix <=70){
      DrawGammaLines(0, -5, 15+2.5*ix, 15+2.5*ix, 1, kGray+2, 1);
      DrawGammaLines(-5, -15, 15+2.5*ix, 15+2.5*ix, 1, kGray+2, 1);
      drawLatexAdd(Form("%1.2f",-TMath::Log(TMath::Tan((15+2.5*ix)/zmid/2))),0.42,0.148+(0.01725/2*ix),0.175*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("#Delta#eta = %1.2f",TMath::Abs(-TMath::Log(TMath::Tan((15+2.5*(ix+1))/zmid/2))-(-TMath::Log(TMath::Tan((15+2.5*ix)/zmid/2))))),0.30,0.168+(0.01725/2*ix),0.175*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    }
  }
  DrawGammaLines(0, 0, 0, 150, 1, kGray+2, 1);
  DrawGammaLines(-10, -10, 70, 150, 1, kGray+2, 1);
  DrawGammaLines(-2.5, -2.5, 15, 70, 1, kGray+2, 1);
  DrawGammaLines(-7.5, -7.5, 15, 70, 1, kGray+2, 1);
  DrawGammaLines(-12.5, -12.5, 15, 70, 1, kGray+2, 1);
  DrawGammaLines(-10, -10, 15, 70, 1, kGray+2, 1);
  DrawGammaLines(-15, -0, 15, 15, 1, kGray+2, 1);

  for(Int_t ix=0;ix<24;ix++){
    DrawGammaLines(0, 5, 15+10*ix, 15+10*ix, 1, kGray+2, 1);
    DrawGammaLines(5, 15, 15+10*ix, 15+10*ix, 1, kGray+2, 1);
    drawLatexAdd(Form("%1.2f",-TMath::Log(TMath::Tan((15+10*ix)/zmid/2))),0.55,0.148+0.0345*ix,0.7*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    if(ix<23)
    drawLatexAdd(Form("#Delta#eta = %1.2f",TMath::Abs(-TMath::Log(TMath::Tan((15+10*(ix+1))/zmid/2))-(-TMath::Log(TMath::Tan((15+10*ix)/zmid/2))))),0.65,0.168+0.0345*ix,0.7*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  }
  DrawGammaLines(-5, -5, 15, 250, 1, kGray+2, 1);
  DrawGammaLines(5, 5, 15, 250, 1, kGray+2, 1);
  DrawGammaLines(15, 15, 0, 250, 1, kGray+2, 1);
  DrawGammaLines(-15, -15, 0, 250, 1, kGray+2, 1);

  // Double_t yin  = 15;
  // cout << "asin: " << TMath::ASin(yin/zmid) << endl;
  // cout << "eta: " << -TMath::Log(TMath::Tan(yin/zmid/2)) << endl;

  // // drawLatexAdd(collisionSystem.Data(),0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  // drawLatexAdd(Form("Anti-#it{k}_{T}, HCAL+ECAL jets"),0.15,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  histoTowerMapDummy2->Draw("axis,same");
  cReso->Print(Form("%s/towerGeometryChanged.%s", outputDir.Data(), suffix.Data()));


  histoTowerMapDummy   = new TH2F("histoTowerMapDummy","histoTowerMapDummy",191,-0.5, 190.5,191,-0.5, 190.5);
  SetStyleHistoTH2ForGraphs(histoTowerMapDummy, "#it{x} (cm)","#it{y} (cm)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.97);
  histoTowerMapDummy->GetXaxis()->SetNoExponent();
  histoTowerMapDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  histoTowerMapDummy->Draw();

  // TEllipse *el4 = new TEllipse(0.0,0.0,250,250,0,90,0);
  //  el4->SetLineColor(4);
  //  el4->SetLineWidth(2);
  //  el4->Draw();
  Double_t zmidEMC = 310;

  for(Int_t ix=0;ix<400;ix++){
    if(19.3725-5.535+5.535*ix <183){
      DrawGammaLines(0, 5.535/2, 19.3725-5.535+5.535*ix, 19.3725-5.535+5.535*ix, 1, kGray+2, 1);
      DrawGammaLines(5.535/2, 5.535/2+5.535, 19.3725-5.535+5.535*ix, 19.3725-5.535+5.535*ix, 1, kGray+2, 1);
      drawLatexAdd(Form("%1.2f",-TMath::Log(TMath::Tan((19.3725-5.535+5.535*ix)/zmidEMC/2))),0.15,0.160+0.0257*ix,0.6*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("#Delta#eta = %1.2f",TMath::Abs(-TMath::Log(TMath::Tan((19.3725-5.535+5.535*(ix+1))/zmidEMC/2))-(-TMath::Log(TMath::Tan((19.3725-5.535+5.535*ix)/zmidEMC/2))))),0.25,0.168+0.0257*ix,0.6*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    }
  }
  DrawGammaLines(5.535/2, 5.535/2, 19.3725-5.535,  190.5, 1, kGray+2, 1);
  DrawGammaLines(5.535/2+5.535, 5.535/2+5.535, 0,  190.5, 1, kGray+2, 1);

  Double_t yinEMC  = 19.3725-5.535;
  cout << "asin: " << TMath::ASin(yinEMC/zmidEMC) << endl;
  cout << "eta: " << -TMath::Log(TMath::Tan(yinEMC/zmidEMC/2)) << endl;

  // drawLatexAdd(collisionSystem.Data(),0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  // drawLatexAdd(Form("Anti-#it{k}_{T}, HCAL+ECAL jets"),0.15,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  histoTowerMapDummy->Draw("axis,same");
  cReso->Print(Form("%s/towerGeometryECAL.%s", outputDir.Data(), suffix.Data()));


}
