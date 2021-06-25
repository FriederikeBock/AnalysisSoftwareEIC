#include <TROOT.h>
#include <TString.h>
#include <TDatime.h>
#include <TF1.h>
#include <TF2.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <TLine.h>
#include <TLatex.h>
#include <TArrow.h>
#include <TGraphErrors.h>
#include <TGaxis.h>
#include <TLegend.h>

//__________________________________________________________________________________________________________
TString ReturnDateStringForOutput(){
  TDatime today;
  int iDate           = today.GetDate();
  int iYear           = iDate/10000;
  int iMonth          = (iDate%10000)/100;
  int iDay            = iDate%100;
  TString cMonth[12]  = {"Jan","Feb","Mar","Apr","May","Jun",
                      "Jul","Aug","Sep","Oct","Nov","Dec"};
  TString textDayth;
  if (iDay== 11){
    textDayth       = "th";
  } else if  (iDay== 12){
    textDayth       = "th";
  } else if  (iDay== 13){
    textDayth       = "th";
  } else if  (iDay%10 == 1){
    textDayth       = "st";
  } else if (iDay%10 == 2){
    textDayth       = "nd";
  } else if (iDay%10 == 3){
    textDayth       = "rd";
  } else {
    textDayth       = "th";
  }
  return Form("%i_%02d_%02d",iYear, iMonth, iDay);
}

//__________________________________________________________________________________________________________
void DrawGammaCanvasSettings( TCanvas* c1,
                            Double_t leftMargin,
                            Double_t rightMargin,
                            Double_t topMargin,
                            Double_t bottomMargin){
  c1->SetTickx();
  c1->SetTicky();
  c1->SetGridx(0);
  c1->SetGridy(0);
  c1->SetLogy(0);
  c1->SetLeftMargin(leftMargin);
  c1->SetRightMargin(rightMargin);
  c1->SetTopMargin(topMargin);
  c1->SetBottomMargin(bottomMargin);
  c1->SetFillColor(0);
}

//__________________________________________________________________________________________________________
void DrawVirtualPadSettings( TVirtualPad* c1,
                      Double_t leftMargin,
                      Double_t rightMargin,
                      Double_t topMargin,
                      Double_t bottomMargin){
  c1->SetTickx();
  c1->SetTicky();
  c1->SetGridx(0);
  c1->SetGridy(0);
  c1->SetLogy(0);
  c1->SetLeftMargin(leftMargin);
  c1->SetRightMargin(rightMargin);
  c1->SetTopMargin(topMargin);
  c1->SetBottomMargin(bottomMargin);
  c1->SetFillColor(0);
}
//__________________________________________________________________________________________________________
void DrawVirtualPadSettingsLogY( TVirtualPad* c1,
                      Double_t leftMargin,
                      Double_t rightMargin,
                      Double_t topMargin,
                      Double_t bottomMargin){
  c1->SetTickx();
  c1->SetTicky();
  c1->SetGridx(0);
  c1->SetGridy(0);
  c1->SetLogy(true);
  c1->SetLeftMargin(leftMargin);
  c1->SetRightMargin(rightMargin);
  c1->SetTopMargin(topMargin);
  c1->SetBottomMargin(bottomMargin);
  c1->SetFillColor(0);
}

//__________________________________________________________________________________________________________
TF1* DivideTF1(TF1* f1, TF1* f2, TString name) {

  if (!f1 || !f2) return NULL;

  Double_t xmin, xmax;
  f1->GetRange(xmin, xmax);
  Int_t nPar1                         = f1->GetNpar();
  Int_t nPar2                         = f2->GetNpar();
  TString formula1                    = f1->GetExpFormula();
  TString formula2                    = f2->GetExpFormula();

  for (Int_t i = 0; i< nPar2; i++){
    formula2.ReplaceAll(Form("[%d]",i), Form("[%d]",i+nPar1));
  }

  TF1* result = new TF1(name.Data(),Form("(%s)/(%s)",formula1.Data(), formula2.Data()), xmin, xmax);
  for (Int_t i = 0; i < nPar1; i++ ){
    result->SetParameter(i, f1->GetParameter(i));
  }
  for (Int_t j = 0; j < nPar2; j++ ){
    result->SetParameter(nPar1+j, f2->GetParameter(j));
  }

  return result;
}

//__________________________________________________________________________________________________________
void DrawGammaSetMarkerTGraph(  TGraph* graph,
                                Style_t markerStyle,
                                Size_t markerSize,
                                Color_t markerColor,
                                Color_t lineColor,
                                Width_t lineWidth       = 1,
                                Style_t lineStyle       = 1,
                                Bool_t boxes            = kFALSE,
                                Color_t fillColor       = 0,
                                Bool_t isHollow         = kFALSE
                              ) {
  graph->SetMarkerStyle(markerStyle);
  graph->SetMarkerSize(markerSize);
  graph->SetMarkerColor(markerColor);
  graph->SetLineColor(lineColor);
  graph->SetLineWidth(lineWidth);
  graph->SetLineWidth(lineStyle);
  if (boxes){
    graph->SetFillColor(fillColor);
    if (fillColor!=0){
      if (!isHollow){
        graph->SetFillStyle(1001);
      } else {
        graph->SetFillStyle(0);
      }
    } else {
      graph->SetFillStyle(0);
    }
  }
}

//__________________________________________________________________________________________________________
void DrawGammaSetMarkerTGraphAsym(  TGraphAsymmErrors* graph,
                                    Style_t markerStyle,
                                    Size_t markerSize,
                                    Color_t markerColor,
                                    Color_t lineColor,
                                    Width_t lineWidth   =1,
                                    Bool_t boxes        = kFALSE,
                                    Color_t fillColor   = 0,
                                    Bool_t isHollow     = kFALSE
                                  ) {
  graph->SetMarkerStyle(markerStyle);
  graph->SetMarkerSize(markerSize);
  graph->SetMarkerColor(markerColor);
  graph->SetLineColor(lineColor);
  graph->SetLineWidth(lineWidth);
  if (boxes){
    graph->SetFillColor(fillColor);
    if (fillColor!=0){
      if (!isHollow){
        graph->SetFillStyle(1001);
      } else {
        graph->SetFillStyle(0);
      }
    } else {
      graph->SetFillStyle(0);
    }
  }
}

//__________________________________________________________________________________________________________
TH1D* CalculateHistoRatioToFit (TH1D* histo, TF1* fit, Bool_t integrateFunction=kFALSE){
  TH1D* histo2                = (TH1D*)histo->Clone("Dummy");
  for( Int_t I = 1; I < histo->GetNbinsX() +1 ;I++){
    Double_t xValue         = histo2->GetBinCenter(I);
    Double_t yValue         = fit->Eval(xValue);
    if (integrateFunction){
      Double_t xMin       = histo2->GetXaxis()->GetBinLowEdge(I);
      Double_t xMax       = histo2->GetXaxis()->GetBinUpEdge(I);
      yValue              = fit->Integral(xMin,xMax)/(xMax-xMin);
    }
    Double_t formerYValue   = histo2->GetBinContent(I);
    if (yValue != 0){
      histo2->SetBinContent(I,formerYValue/yValue);
      histo2->SetBinError(I,histo2->GetBinError(I)/yValue);
    }
  }
  return histo2;
}
//__________________________________________________________________________________________________________
void SetStyleHistoTH2ForGraphs( TH2* histo,
                                TString XTitle,
                                TString YTitle,
                                Size_t xLableSize,
                                Size_t xTitleSize,
                                Size_t yLableSize,
                                Size_t yTitleSize,
                                Float_t xTitleOffset    = 1,
                                Float_t yTitleOffset    = 1,
                                Int_t xNDivisions       = 510,
                                Int_t yNDivisions       = 510,
                                Font_t textFontLabel    = 42,
                                Font_t textFontTitle    = 62
                              ){
  histo->SetXTitle(XTitle);
  histo->SetYTitle(YTitle);
  histo->SetTitle("");

  histo->GetXaxis()->SetLabelFont(textFontLabel);
  histo->GetYaxis()->SetLabelFont(textFontLabel);
  histo->GetZaxis()->SetLabelFont(textFontLabel);
  histo->GetXaxis()->SetTitleFont(textFontTitle);
  histo->GetYaxis()->SetTitleFont(textFontTitle);
  histo->GetZaxis()->SetTitleFont(textFontTitle);

  histo->GetXaxis()->SetLabelSize(xLableSize);
  histo->GetXaxis()->SetTitleSize(xTitleSize);
  histo->GetXaxis()->SetTitleOffset(xTitleOffset);
  histo->GetXaxis()->SetNdivisions(xNDivisions,kTRUE);

  histo->GetYaxis()->SetDecimals();
  histo->GetYaxis()->SetLabelSize(yLableSize);
  histo->GetYaxis()->SetTitleSize(yTitleSize);
  histo->GetYaxis()->SetTitleOffset(yTitleOffset);
  histo->GetYaxis()->SetNdivisions(yNDivisions,kTRUE);
  
  histo->GetZaxis()->SetDecimals();
  histo->GetZaxis()->SetLabelSize(yLableSize);
  histo->GetZaxis()->SetTitleSize(yTitleSize);
  
}

  //__________________________________________________________________________________________________________
void SetStyleHistoTH1ForGraphs( TH1* histo,
                                TString XTitle,
                                TString YTitle,
                                Size_t xLableSize,
                                Size_t xTitleSize,
                                Size_t yLableSize,
                                Size_t yTitleSize,
                                Float_t xTitleOffset    = 1,
                                Float_t yTitleOffset    = 1,
                                Int_t xNDivisions       = 510,
                                Int_t yNDivisions       = 510,
                                Font_t textFontLabel    = 42,
                                Font_t textFontTitle    = 62
                              ){
  histo->SetXTitle(XTitle);
  histo->SetYTitle(YTitle);
  histo->SetTitle("");

  histo->GetXaxis()->SetLabelFont(textFontLabel);
  histo->GetYaxis()->SetLabelFont(textFontLabel);
  histo->GetXaxis()->SetTitleFont(textFontTitle);
  histo->GetYaxis()->SetTitleFont(textFontTitle);

  histo->GetXaxis()->SetLabelSize(xLableSize);
  histo->GetXaxis()->SetTitleSize(xTitleSize);
  histo->GetXaxis()->SetTitleOffset(xTitleOffset);
  histo->GetXaxis()->SetNdivisions(xNDivisions,kTRUE);

  histo->GetYaxis()->SetDecimals();
  histo->GetYaxis()->SetLabelSize(yLableSize);
  histo->GetYaxis()->SetTitleSize(yTitleSize);
  histo->GetYaxis()->SetTitleOffset(yTitleOffset);
  histo->GetYaxis()->SetNdivisions(yNDivisions,kTRUE);
}


void DrawGammaLines(Float_t startX, Float_t endX,
                Float_t startY, Float_t endY,
                Float_t linew, Float_t lineColor = 4, Style_t lineStyle = 1){
    TLine * l1 = new TLine (startX,startY,endX,endY);
    l1->SetLineColor(lineColor);
    l1->SetLineWidth(linew);
    l1->SetLineStyle(lineStyle);
    l1->Draw("same");
}
//__________________________________________________________________________________________________________
void SetStyleTLatex( TLatex* text,
                    Size_t textSize,
                    Width_t lineWidth,
                    Color_t textColor = 1,
                    Font_t textFont = 42,
                    Bool_t kNDC = kTRUE,
                    Short_t align = 11
                ){
  if (kNDC) {text->SetNDC();}
  text->SetTextFont(textFont);
  text->SetTextColor(textColor);
  text->SetTextSize(textSize);
  text->SetLineWidth(lineWidth);
  text->SetTextAlign(align);
}

//__________________________________________________________________________________________________________
void drawLatexAdd(TString latextext, 
                  Double_t textcolumn, 
                  Double_t textrow, 
                  Double_t textSizePixel,
                  Bool_t setFont        = kFALSE, 
                  Bool_t setFont2       = kFALSE, 
                  Bool_t alignRight     = kFALSE, 
                  Color_t textcolor     = kBlack
                 ){
  TLatex *latexDummy                  = new TLatex(textcolumn ,textrow,latextext);
  SetStyleTLatex( latexDummy, textSizePixel,4);
  if(setFont)
    latexDummy->SetTextFont(62);
  if(setFont2)
    latexDummy->SetTextFont(43);
  if(alignRight)
    latexDummy->SetTextAlign(31);
  latexDummy->SetTextColor(textcolor);
  latexDummy->Draw();
}

//__________________________________________________________________________________________________________
void drawArrowAdd(Double_t x1, Double_t y1, Double_t x2, Double_t y2,Float_t size = 0.01, Color_t arrcolor = kBlack){
  TArrow *ar2 = new TArrow(x1,y1,x2,y2,size,"|>");
  ar2->SetAngle(60);
  ar2->SetLineColor(arrcolor);
  ar2->SetLineWidth(2);
  ar2->Draw();
}
//__________________________________________________________________________________________________________
void DrawGammaSetMarker(    TH1* histo1,
                            Style_t markerStyle,
                            Size_t markerSize,
                            Color_t markerColor,
                            Color_t lineColor ) {
  histo1->SetMarkerStyle(markerStyle);
  histo1->SetMarkerSize(markerSize);
  histo1->SetMarkerColor(markerColor);
  histo1->SetLineColor(lineColor);
  histo1->GetYaxis()->SetLabelFont(42);
  histo1->GetXaxis()->SetLabelFont(42);
  histo1->GetYaxis()->SetTitleFont(62);
  histo1->GetXaxis()->SetTitleFont(62);
}


    TF1* ScaleTF1(TF1* func, Double_t constant, TString name) {

        if (!func) return NULL;

        Double_t    xMin, xMax;
        TString     formula         = func->GetExpFormula();
        func->GetRange(xMin, xMax);
            #if !defined (__CINT__) || defined (__CLING__)
            for (Int_t i=0; i<func->GetNpar(); i++) {
                formula.ReplaceAll(Form("[p%d]", i), Form("[placeholder%d]",i+1));
            }
            for (Int_t i=1; i<func->GetNpar()+1; i++) {
                formula.ReplaceAll(Form("[placeholder%d]", i), Form("[p%d]",i));
            }
        #else
            for (Int_t i=0; i<func->GetNpar(); i++) {
                formula.ReplaceAll(Form("[%d]", i), Form("[placeholder%d]",i+1));
            }
            for (Int_t i=1; i<func->GetNpar()+1; i++) {
                formula.ReplaceAll(Form("[placeholder%d]", i), Form("[%d]",i));
            }
        #endif

        TF1* result                 = new TF1(name.Data(), Form("[0] * (%s)", formula.Data()), xMin, xMax);
        for (Int_t i=0; i<func->GetNpar()+1; i++) {
            if (i==0)   result->SetParameter(i, constant);
            else        result->SetParameter(i, func->GetParameter(i-1));
        }

        return result;
    }

//__________________________________________________________________________________________________________
void StyleSettingsThesis( TString format = ""){
  //gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptDate(0);   //show day and time
  gStyle->SetOptStat(0);  //show statistic
  gStyle->SetPalette(1,0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTextSize(0.5);
  gStyle->SetLabelSize(0.03,"xyz");
  gStyle->SetLabelOffset(0.002,"xyz");
  gStyle->SetTitleFontSize(0.04);
  gStyle->SetTitleOffset(1,"y");
  gStyle->SetTitleOffset(0.7,"x");
  gStyle->SetCanvasColor(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLineWidth(1);

  gStyle->SetPadTopMargin(0.03);
  gStyle->SetPadBottomMargin(0.09);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetPadLeftMargin(0.13);


  TGaxis::SetMaxDigits(3);
  gErrorIgnoreLevel=kError;

  if (format.CompareTo("eps") == 0 ||format.CompareTo("pdf") == 0  ) gStyle->SetLineScalePS(1);
}

//__________________________________________________________________________________________________________
void DrawGammaSetMarkerTGraphErr(   TGraphErrors* graph,
                                    Style_t markerStyle,
                                    Size_t markerSize,
                                    Color_t markerColor,
                                    Color_t lineColor,
                                    Width_t lineWidth       = 1,
                                    Bool_t boxes            = kFALSE,
                                    Color_t fillColor       = 0,
                                    Bool_t isHollow         = kFALSE) {
  graph->SetMarkerStyle(markerStyle);
  graph->SetMarkerSize(markerSize);
  graph->SetMarkerColor(markerColor);
  graph->SetLineColor(lineColor);
  graph->SetLineWidth(lineWidth);
  if (boxes){
    graph->SetFillColor(fillColor);
    if (fillColor!=0){
      if (!isHollow){
        graph->SetFillStyle(1001);
      } else {
        graph->SetFillStyle(0);
      }
    } else {
      graph->SetFillStyle(0);
    }
  }
}

//__________________________________________________________________________________________________________
TGraphErrors* CalculateGraphErrRatioToFit (TGraphErrors* graph_Org, TF1* fit){
  TGraphErrors* graph         = (TGraphErrors*)graph_Org->Clone(Form("%s_Dummy",graph_Org->GetName()));
  Double_t * xValue           = graph->GetX();
  Double_t * yValue           = graph->GetY();
  Double_t* xError            = graph->GetEX();
  Double_t* yError            = graph->GetEY();
  Int_t nPoints               = graph->GetN();
  for (Int_t i = 0; i < nPoints; i++){
    yValue[i]               = yValue[i]/fit->Eval(xValue[i]);
    yError[i]               = yError[i]/fit->Eval(xValue[i]);
  }
  TGraphErrors* returnGraph   = new TGraphErrors(nPoints,xValue,yValue,xError,yError);
  return returnGraph;
}

//__________________________________________________________________________________________________________
TGraph* CalculateGraphRatioToFit (TGraph* graph_Org, TF1* fit){
  TGraphErrors* graph         = (TGraphErrors*)graph_Org->Clone("Dummy2");
  Double_t * xValue       = graph->GetX();
  Double_t * yValue       = graph->GetY();
  Int_t nPoints           = graph->GetN();
  for (Int_t i = 0; i < nPoints; i++){
    yValue[i]           = yValue[i]/fit->Eval(xValue[i]);
  }
  TGraph* returnGraph     = new TGraph(nPoints,xValue,yValue);
  return returnGraph;
}

//__________________________________________________________________________________________________________
void DrawGammaSetMarkerTF1( TF1* fit1,
                            Style_t lineStyle,
                            Size_t lineWidth,
                            Color_t lineColor ) {
  fit1->SetLineColor(lineColor);
  fit1->SetLineStyle(lineStyle);
  fit1->SetLineWidth(lineWidth);
}

//__________________________________________________________________________________________________________
TLegend *GetAndSetLegend2(  Double_t positionX,
                            Double_t positionY,
                            Double_t positionXRight,
                            Double_t positionYUp,
                            Size_t textSize,
                            Int_t columns               = 1,
                            TString header              = "",
                            Font_t textFont             = 43,
                            Double_t margin             = 0
                          ){

  TLegend *legend = new TLegend(positionX,positionY,positionXRight,positionYUp);
  legend->SetNColumns(columns);
  legend->SetLineColor(0);
  legend->SetLineWidth(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetLineStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(textFont);
  legend->SetTextSize(textSize);
  if (margin != 0) legend->SetMargin(margin);
  if (header.CompareTo("")!= 0) legend->SetHeader(header);
  return legend;
}

//__________________________________________________________________________________________________________
void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs){
  
  if(numInputs<5){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 2, gStyle->GetCanvasDefH()*2);
    DrawGammaCanvasSettings( cPNG, 0, 0, 0, 0);
    cPNG->Divide(2, 2,0,0);
  }else if(numInputs<7){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 3, gStyle->GetCanvasDefH()*2);
    DrawGammaCanvasSettings( cPNG, 0, 0, 0, 0);
    cPNG->Divide(3, 2,0,0);
  }else if(numInputs<10){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 3, gStyle->GetCanvasDefH()*3);
    DrawGammaCanvasSettings( cPNG, 0, 0, 0, 0);
    cPNG->Divide(3, 3,0,0);
  } else if(numInputs<13){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 4, gStyle->GetCanvasDefH()*3);
    DrawGammaCanvasSettings( cPNG, 0, 0, 0, 0);
    cPNG->Divide(4, 3,0,0);
  } else if(numInputs<17){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 4, gStyle->GetCanvasDefH()*4);
    DrawGammaCanvasSettings( cPNG, 0, 0, 0, 0);
    cPNG->Divide(4, 4,0,0);
  } else if(numInputs<21){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 5, gStyle->GetCanvasDefH()*4);
    DrawGammaCanvasSettings( cPNG, 0, 0, 0, 0);
    cPNG->Divide(5, 4,0,0);
  } else {
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 5, gStyle->GetCanvasDefH()*5);
    DrawGammaCanvasSettings( cPNG, 0, 0, 0, 0);
    cPNG->Divide(5, 5,0,0);
  }

}

//__________________________________________________________________________________________________________
void SetPlotStyle() {
// 	const Int_t nRGBs = 7;
  const Int_t nRGBs = 5;
  const Int_t nCont = 255;

  Double_t stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[nRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[nRGBs] = { 0.31, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[nRGBs]  = { 0.51, 1., 0.12, 0.00, 0.00};

  TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);
  gStyle->SetNumberContours(nCont);
}


TString GetTrackerLabel(TString label){
  TString detLabel = "";
  if (label.Contains("defaultLBL") ){
    detLabel  = "ALL-SI";
  } else if (label.Contains("LBLwithLGAD") ){
    detLabel  = "ALL-SI+TTL";
  } else if (label.Contains("LBLwithACLGAD") ){
    detLabel  = "ALL-SI+TTL(AC-LGAD)";
  } else if (label.Contains("LBLwithFTTLS2LC-ETTL-CTTL") ){
    detLabel  = "ALL-SI+TTL(2 l's)";
  } else if (label.Contains("LBLwithFTTLS3LVC-ETTLLC-CTTLLC") ){
    detLabel  = "ALL-SI+TTL(1.3mm)";
  } else if (label.Contains("LBLwithFTTLSE1LC-ETTLSE1-CTTLSE1") ){
    detLabel  = "ALL-SI+TTL(1 l b. ECal)";
  } else if (label.Contains("LBLwithFTTLSE2LC-ETTL-CTTLSE1") ){
    detLabel  = "ALL-SI+TTL(1|2 l b. ECal)";
  } else if (label.Contains("defaultLANL") ){
    detLabel  = "LANL";
  } else if (label.Contains("LANLwithLGAD") ){
    detLabel  = "LANL+TTL";
  } else if (label.Contains("LANLwithACLGAD") ){
    detLabel  = "LANL+TTL(AC-LGAD)";
  } else if (label.Contains("LANLwithFTTLS3LVC-ETTLLC-CTTLLC") ){
    detLabel  = "LANL+TTL(1.3mm)";
  }
  return detLabel;
}

TString GetMagnetLabel(TString label){
  TString magnetLabel = "BABAR (B = 1.5T)";
  if (label.Contains("3T")){
    magnetLabel = "BEAST (B = 3T)";
  }
  return magnetLabel;
}

TString GetCollisionEnergy(TString label){
  TString collLabel = "Pythia 6, e+p, 10+250 GeV";
  if (label.Contains("e5p100") ){
    collLabel = "Pythia 6, e+p, 5+100 GeV";
  } else if (label.Contains("e18p275") ){
    collLabel = "Pythia 6, e+p, 18+275 GeV";
  } else if (label.Contains("EnergyDep") ){
    collLabel = "Pythia 6, e+p";
  }
  return collLabel;
}
