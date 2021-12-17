#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <stdlib.h>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TSpectrum.h"
#include "TLegend.h"


Double_t langaufun(Double_t *x, Double_t *par)
{
  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation),
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location

  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;

  // MP shift correction
  mpc = par[1] - mpshift * par[0];

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp-xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);

    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}

std::vector<int> ParseRunList(const std::string& str)
{
  std::vector<int> runList;

  std::stringstream ss(str);
  std::string token;
  while( std::getline(ss,token,',') )
  {
    std::stringstream ss2(token);
    std::string token2;
    int runMin = -1;
    int runMax = -1;
    while( std::getline(ss2,token2,'-') )
    {
      if( runMin != -1 && runMax == -1 ) runMax = atoi(token2.c_str());
      if( runMin == -1 ) runMin = atoi(token2.c_str());
    }
    if( runMax == -1 ) runMax = runMin;

    for(int run = runMin; run <= runMax; ++run)
    {
      runList.push_back(run);
    }
  }

  return runList;
}


void drawP1(TProfile* prof, std::string xtitle, std::string ytitle, std::string plotDir, std::string name="")
{
  if( prof == NULL ) return;

  gStyle -> SetOptStat(0);

  TCanvas* c = new TCanvas();
  prof -> SetMarkerStyle(20);
  prof -> SetMarkerColor(kBlack);
  prof -> GetXaxis() -> SetTitle(xtitle.c_str());
  prof -> GetYaxis() -> SetTitle(ytitle.c_str());
  prof -> Draw("P");
  if( name!="" ) c -> Print(Form("%s/%s.png",plotDir.c_str(),name.c_str()));
  else c -> Print(Form("%s/%s.png",plotDir.c_str(),prof->GetName()));

  delete c;
}

void drawP1_fitPol(TProfile* prof, std::string xtitle, std::string ytitle, std::string plotDir, TF1* f_pol, std::string name="")
{
  if( prof == NULL ) return;

  gStyle -> SetOptStat(0);

  TCanvas* c = new TCanvas();
  prof -> SetMarkerStyle(20);
  prof -> SetMarkerColor(kBlack);
  prof -> GetXaxis() -> SetTitle(xtitle.c_str());
  prof -> GetYaxis() -> SetTitle(ytitle.c_str());
  prof -> GetYaxis() -> SetRangeUser(-1000.,1000.);
  prof -> Draw("P");

  prof -> Fit(f_pol,"QNRS");
  f_pol -> Draw("same");

  if( name!="" ) c -> Print(Form("%s/%s.png",plotDir.c_str(),name.c_str()));
  else c -> Print(Form("%s/%s.png",plotDir.c_str(),prof->GetName()));

  delete c;
}


void drawH1(TH1F* h, std::string xtitle, std::string ytitle, std::string plotDir, bool logy=false, std::string name="", float minimum=-999., bool stat=true)
{
  if( h == NULL ) return;
  if( h->GetEntries() < 100 ) return;

  if( stat ) gStyle -> SetOptStat(1111);
  else gStyle -> SetOptStat(0);

  h -> SetMarkerStyle(20);
  h -> SetMarkerColor(kBlack);
  h -> SetLineColor(kBlack);
  h -> SetLineWidth(2);
  h -> GetXaxis() -> SetTitle(xtitle.c_str());
  h -> GetYaxis() -> SetTitle(ytitle.c_str());
  if( minimum != -999.) h -> SetMinimum( minimum );

  TCanvas* c = new TCanvas();
  if( logy ) c -> SetLogy();
  h -> Draw("Hist");
  if( name!="" ) c -> Print(Form("%s/%s.png",plotDir.c_str(),name.c_str()));
  else c -> Print(Form("%s/%s.png",plotDir.c_str(),h -> GetName()));

  delete c;
}


void drawH1_fitGaus(TH1F* h, float xMin, float xMax, std::string xtitle, std::string ytitle, std::string plotDir, TF1* f_gaus, bool logy=false, std::string name="")
{
  if( h == NULL ) return;
  if( h->GetEntries() < 100 ) return;

  gStyle -> SetOptStat(1111);

  h -> SetLineColor(kBlack);
  h -> SetLineWidth(2);
  h -> GetXaxis() -> SetTitle(xtitle.c_str());
  h -> GetYaxis() -> SetTitle(ytitle.c_str());
  h -> GetXaxis() -> SetRangeUser(xMin,xMax);

  f_gaus -> SetRange(xMin,xMax);
  f_gaus -> SetParameters(h->GetMaximum(),h->GetMean(),h->GetRMS());
  h -> Fit(f_gaus,"QNRLS");

  TCanvas* c = new TCanvas();
  if( logy ) c -> SetLogy();
  h -> Draw("Hist");
  f_gaus -> Draw("same");
  if( name!="" ) c -> Print(Form("%s/%s.png",plotDir.c_str(),name.c_str()));
  else c -> Print(Form("%s/%s.png",plotDir.c_str(),h -> GetName()));

  delete c;
}


void drawH1_arrays(TH1F* h1, TH1F* h2, std::string xtitle, std::string ytitle, std::string plotDir, bool logy=false, std::string name="", float minimum=-999., bool stat=true)
{
  if( h1 == NULL ) return;
  if( h2 == NULL ) return;

  if( stat ) gStyle -> SetOptStat(1111);
  else gStyle -> SetOptStat(0);

  h1 -> SetLineColor(kRed);
  h1 -> SetLineWidth(3);
  h1 -> GetXaxis() -> SetTitle(xtitle.c_str());
  h1 -> GetYaxis() -> SetTitle(ytitle.c_str());
  if( minimum != -999.) h1 -> SetMinimum( minimum );

  h2 -> SetLineColor(kBlue);
  h2 -> SetLineWidth(3);

  h1 -> SetMaximum(1.2*h1->GetMaximum());
  if( h2->GetMaximum() > h1->GetMaximum() )
    h1 -> SetMaximum(1.2*h2->GetMaximum());

  TCanvas* c = new TCanvas();
  if( logy ) c -> SetLogy();
  h1 -> Draw("Hist");
  h2 -> Draw("Hist,same");
  if( name!="" ) c -> Print(Form("%s/%s.png",plotDir.c_str(),name.c_str()));
  else c -> Print(Form("%s/%s.png",plotDir.c_str(),h1 -> GetName()));

  delete c;
}


void drawG_vector(std::vector<TGraphErrors*> g1s, std::string xtitle, std::string ytitle, float yMin, float yMax, std::string plotDir, bool logy=false, std::string name="", std::vector<std::string>* labels = NULL)
{

  TCanvas* c = new TCanvas();
  if( logy ) c -> SetLogy();
  gPad -> SetGridx();
  gPad -> SetGridy();

  TH1F* hPad = (TH1F*)( gPad->DrawFrame(-0.5,yMin,63.5,yMax) );
  hPad -> SetTitle(Form(";%s;%s",xtitle.c_str(),ytitle.c_str()));
  hPad -> Draw();

  TLegend* legend;
  if( labels )
  {
    legend = new TLegend(0.50,0.94-labels->size()*0.04,0.90,0.94);
    legend -> SetFillColor(0);
    legend -> SetFillStyle(1000);
    legend -> SetTextFont(42);
    legend -> SetTextSize(0.03);
  }

  int ii = 0;
  for(auto g1 : g1s)
  {
    g1 -> SetMarkerStyle(20);
    g1 -> SetMarkerColor(kBlack+ii);
    g1 -> SetLineColor(kBlack+ii);
    g1 -> GetXaxis() -> SetTitle(xtitle.c_str());
    g1 -> GetYaxis() -> SetTitle(ytitle.c_str());
    g1 -> Draw("PL,same");

    if( labels ) legend -> AddEntry(g1,labels->at(ii).c_str(),"PL");
    ++ii;
  }

  if( legend ) legend -> Draw("same");

  if( name!="" ) c -> Print(Form("%s/%s.png",plotDir.c_str(),name.c_str()));
  else c -> Print(Form("%s/%s.png",plotDir.c_str(),g1s.at(0)->GetName()));

  delete c;
}


void drawG_arrays(TGraphErrors* g1, TGraphErrors* g2, std::string xtitle, std::string ytitle, float yMin, float yMax, std::string plotDir, bool logy=false, std::string name="")
{
  if( g1 == NULL || g2 == NULL ) return;

  g1 -> SetMarkerStyle(20);
  g1 -> SetMarkerColor(kRed);
  g1 -> GetXaxis() -> SetTitle(xtitle.c_str());
  g1 -> GetYaxis() -> SetTitle(ytitle.c_str());

  g2 -> SetMarkerStyle(21);
  g2 -> SetMarkerColor(kBlue);

  TCanvas* c = new TCanvas();
  if( logy ) c -> SetLogy();
  gPad -> SetGridx();
  gPad -> SetGridy();

  TH1F* hPad = (TH1F*)( gPad->DrawFrame(-0.5,yMin,15.5,yMax) );
  hPad -> SetTitle(Form(";%s;%s",xtitle.c_str(),ytitle.c_str()));
  hPad -> Draw();

  g1 -> Draw("P,same");
  g2 -> Draw("P,same");

  if( name!="" ) c -> Print(Form("%s/%s.png",plotDir.c_str(),name.c_str()));
  else c -> Print(Form("%s/%s.png",plotDir.c_str(),g1->GetName()));

  delete c;
}


void drawH1_energy(TH1F* h_LR, TH1F* h_L, TH1F* h_R,
                   float xmin, float xmax, std::string xtitle, std::string ytitle,
                   std::string plotDir, TF1* f_langaus, bool logy=true, std::string name="")
{
  if( h_LR == NULL ) return;
  if( h_LR->GetEntries() < 100 ) return;

  gStyle -> SetOptStat(1111);

  h_LR -> SetLineColor(kBlack);
  h_LR -> SetLineWidth(2);
  h_LR -> GetXaxis() -> SetTitle(xtitle.c_str());
  h_LR -> GetYaxis() -> SetTitle(ytitle.c_str());
  h_LR -> GetXaxis() -> SetRangeUser(xmin,xmax);

  h_L -> SetLineColor(kRed);
  h_L -> SetLineWidth(1);
  h_L -> SetLineStyle(3);
  h_R -> SetLineColor(kBlue);
  h_R -> SetLineWidth(1);
  h_R -> SetLineStyle(3);

  TCanvas* c = new TCanvas();
  if( logy ) c -> SetLogy();
  h_LR -> Draw("Hist");
  h_L -> Draw("same");
  h_R -> Draw("same");

  float xMax = 0.;
  float max = 0;
  for(int bin = 1; bin <= h_LR->GetNbinsX(); ++bin)
  {
    float binCenter = h_LR -> GetBinCenter(bin);
    float binContent = h_LR -> GetBinContent(bin);
    if( binCenter < 200 ) continue;
    if( binContent > max )
    {
	  max = binContent;
	  xMax = binCenter;
	}
  }

  TF1* f_landau = new TF1("f_landau","[0]*TMath::Landau(x,[1],[2])",0.8*xMax,1000.);
  f_landau -> SetParameters(100.,xMax,10.);
  h_LR -> Fit(f_landau,"QNRS");

  f_langaus -> SetRange(0.8*f_landau->GetParameter(1),1000.);
  f_langaus -> SetParameters(f_landau->GetParameter(2),f_landau->GetParameter(1),h_LR->Integral(h_LR->FindBin(0.8*xMax),h_LR->FindBin(1000.))*h_LR->GetBinWidth(1),10.);
  h_LR -> Fit(f_langaus,"QRS+");
  f_langaus -> Draw("same");

  TLine* line = new TLine(0.80*f_langaus->GetParameter(1),h_LR->GetMinimum(),0.80*f_langaus->GetParameter(1),h_LR->GetMaximum());
  line -> SetLineStyle(7);
  line -> SetLineWidth(2);
  line -> SetLineColor(kBlack);
  line -> Draw("same");

  if( name!="" ) c -> Print(Form("%s/%s.png",plotDir.c_str(),name.c_str()));
  else c -> Print(Form("%s/%s.png",plotDir.c_str(),h_LR -> GetName()));

  //delete spectrum;
  delete f_landau;
  delete line;
  delete c;
}

void drawH1_deltaT(TH1F* h_raw, TH1F* h_corr,
                   std::string xtitle, std::string ytitle,
                   std::string plotDir, TF1* f_raw, TF1* f_corr, bool logy=false, std::string name="")
{
  if( h_raw == NULL ) return;
  if( h_raw -> GetEntries() < 100 ) return;
  if( h_corr == NULL ) return;
  if( h_corr -> GetEntries() < 100 ) return;

  gStyle -> SetOptStat(1111);

  h_corr -> SetLineColor(kBlue);
  h_corr -> SetLineWidth(2);
  h_corr -> GetXaxis() -> SetTitle(xtitle.c_str());
  h_corr -> GetYaxis() -> SetTitle(ytitle.c_str());

  h_raw -> SetLineColor(kRed);
  h_raw -> SetLineWidth(2);

  f_corr -> SetParameters(h_corr->GetMaximum(),h_corr->GetMean(),h_corr->GetRMS());
  h_corr -> Fit(f_corr,"QNLS","",h_corr->GetMean()-2.*h_corr->GetRMS(),h_corr->GetMean()+2.*h_corr->GetRMS());
  h_corr -> Fit(f_corr,"QLS+","",f_corr->GetParameter(1)-2.*f_corr->GetParameter(2),f_corr->GetParameter(1)+2.*f_corr->GetParameter(2));
  f_corr -> SetLineColor(kBlue);
  f_corr -> SetLineWidth(3);

  f_raw -> SetParameters(h_raw->GetMaximum(),h_raw->GetMean(),h_raw->GetRMS());
  h_raw -> Fit(f_raw,"QNLS","",h_raw->GetMean()-2.*h_raw->GetRMS(),h_raw->GetMean()+2.*h_raw->GetRMS());
  h_raw -> Fit(f_raw,"QLS+","",f_raw->GetParameter(1)-2.*f_raw->GetParameter(2),f_raw->GetParameter(1)+2.*f_raw->GetParameter(2));
  f_raw -> SetLineColor(kRed);
  f_raw -> SetLineWidth(3);

  h_corr -> GetXaxis() -> SetRangeUser(f_corr->GetParameter(1)-7.*f_corr->GetParameter(2),f_corr->GetParameter(1)+7.*f_corr->GetParameter(2));

  TCanvas* c = new TCanvas();
  // if( logy ) c -> SetLogy();
  h_corr -> Draw("Hist");
  h_raw -> Draw("same");
  f_corr -> Draw("same");
  f_raw -> Draw("same");

  if( name!="" ) c -> Print(Form("%s/%s.png",plotDir.c_str(),name.c_str()));
  else c -> Print(Form("%s/%s.png",plotDir.c_str(),h_corr -> GetName()));

  delete c;
}


void drawH2(TH2F* h, std::string xtitle, std::string ytitle, std::string ztitle, std::string plotDir, std::vector<std::string> labels={}, bool logz=false, std::string name="")
{
  gStyle -> SetOptStat(0);

  h -> SetMarkerStyle(20);
  h -> SetMarkerColor(kBlack);
  h -> SetLineColor(kBlack);
  h -> SetLineWidth(2);
  // h -> GetXaxis() -> SetNdivisions(2);
  h -> GetXaxis() -> SetTitle(xtitle.c_str());
  h -> GetYaxis() -> SetTitle(ytitle.c_str());
  h -> GetZaxis() -> SetTitle(ztitle.c_str());
  // h -> GetXaxis() -> SetBinLabel(1,"L");
  // h -> GetXaxis() -> SetBinLabel(2,"LR");
  // h -> GetXaxis() -> SetBinLabel(3,"R");
  for( unsigned int i=0; i<labels.size(); i++ )
        h -> GetYaxis() -> SetBinLabel(i+1,labels.at(i).c_str());

  TCanvas* c = new TCanvas("c","c",1000,700);
  if( logz ) c -> SetLogz();
  h -> Draw("COLZ");
  if( name!="" ) c -> Print(Form("%s/%s.png",plotDir.c_str(),name.c_str()));
  else c -> Print(Form("%s/%s.png",plotDir.c_str(),h -> GetName()));

  delete c;
}


int main(int argc, char** argv)
{
  gStyle -> SetOptStat(1111);

  int vth1Ref = 10;

  // parse arguments
  std::string runListStr(argv[1]);
  float timeLeng = 1.; if( argc > 2 ) timeLeng = atof(argv[2]);
  int prescale = 1;    if( argc > 3 ) prescale = atoi(argv[3]);
  int pedestals = 0;   if( argc > 4 ) pedestals = atoi(argv[4]);

  // open file
  TChain* data = new TChain("data","data");
  std::vector<int> runList = ParseRunList(runListStr);
  for(auto run : runList)
  {
    std::string inFileName = "/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_Jul2021/TOFHIR2/h8/reco/";
    if( pedestals ) inFileName.append(Form("%04d/*_ped_e.root",run));
    else            inFileName.append(Form("%04d/*_e.root",run));
    data -> Add(inFileName.c_str());
  }
  int nEntries = data -> GetEntries();

  std::cout << "\n\n============================================"  << std::endl;
  std::cout << "Runs:           " << runListStr << std::endl;
  std::cout << "Entries:        " << nEntries << std::endl;
  std::cout << "timeLeng:       " << timeLeng << " s" << std::endl;
  std::cout << "PreScale:       " << prescale << std::endl;
  std::cout << "Pedestals:      " << pedestals << std::endl;
  std::cout << "============================================\n\n"  << std::endl;

  // define branches
  float step1, step2;
  int channelIdx[128];
  std::vector<float>* tot = 0;
  std::vector<float>* energy = 0;
  std::vector<long long>* time = 0;
  std::vector<float>* qT1 = 0;
  std::vector<unsigned short>* t1fine = 0;
  data -> SetBranchStatus("*",0);
  data -> SetBranchStatus("step1",     1); data -> SetBranchAddress("step1",    &step1);
  data -> SetBranchStatus("step2",     1); data -> SetBranchAddress("step2",    &step2);
  data -> SetBranchStatus("channelIdx",1); data -> SetBranchAddress("channelIdx",channelIdx);
  data -> SetBranchStatus("tot",       1); data -> SetBranchAddress("tot",      &tot);
  data -> SetBranchStatus("energy",    1); data -> SetBranchAddress("energy",   &energy);
  data -> SetBranchStatus("time",      1); data -> SetBranchAddress("time",     &time);
  data -> SetBranchStatus("qT1",       1); data -> SetBranchAddress("qT1",      &qT1);
  data -> SetBranchStatus("t1fine",    1); data -> SetBranchAddress("t1fine",   &t1fine);

  // create plot folder
  std::string plotDir(Form("plots/",runListStr.c_str()));
  system(Form("mkdir -p %s",plotDir.c_str()));
  system(Form("mkdir -p %s/energy/",plotDir.c_str()));
  system(Form("mkdir -p %s/tot/",plotDir.c_str()));
  system(Form("mkdir -p %s/energyRatioCorr/",plotDir.c_str()));
  // system(Form("mkdir -p %s/totRatioCorr/",plotDir.c_str()));
  system(Form("mkdir -p %s/phaseCorr/",plotDir.c_str()));
  system(Form("mkdir -p %s/deltaT/",plotDir.c_str()));
  system(Form("mkdir -p %s/delta_tAvg/",plotDir.c_str()));
  system(Form("mkdir -p %s/timeResolution/",plotDir.c_str()));

  // define colors
  std::vector<int> colors;
  colors.push_back(kBlack);
  colors.push_back(kRed+1);
  colors.push_back(kBlue+1);
  colors.push_back(kOrange+1);
  colors.push_back(kGreen+1);
  colors.push_back(kYellow+1);
  colors.push_back(kCyan+1);
  colors.push_back(kViolet+1);
  colors.push_back(kGray+1);
  colors.push_back(kSpring+1);
  colors.push_back(kTeal+1);
  colors.push_back(kMagenta+1);
  colors.push_back(kAzure+1);

  // define channel mapping
  std::vector<int> ch_array_side1 = {14,12,10, 8, 6, 4, 2, 0, 1, 3, 5, 7, 9,11,13,15};
  std::vector<int> ch_array_side2 = {17,19,21,23,25,27,29,31,30,28,26,24,22,20,18,16};
  int num_bars = ch_array_side1.size();
  std::vector<std::string> ch_labels_array0;
  std::vector<std::string> ch_labels_array1;
  for(unsigned int i = 0; i < num_bars; ++i)
  {
    ch_labels_array0.push_back(std::string(Form("[%i,%i]",ch_array_side1.at(i),ch_array_side2.at(i))));
    ch_labels_array1.push_back(std::string(Form("[%i,%i]",ch_array_side1.at(i)+64,ch_array_side2.at(i)+64)));
  }

  // define histograms
  TFile* outFile = new TFile(Form("plots/simple_analysis_%s.root",runListStr.c_str()),"RECREATE");

  std::map<float, std::map<int, std::map<int, std::map<int,TH1F*> > > > h1_energy_L;
  std::map<float, std::map<int, std::map<int, std::map<int,TH1F*> > > > h1_energy_R;
  std::map<float, std::map<int, std::map<int, std::map<int,TH1F*> > > > h1_energy_LR;
  std::map<float, std::map<int, std::map<int, std::map<int,TF1*> > > > fit_energy_LR;

  std::map<float, std::map<int, std::map<int, std::map<int,TH1F*> > > > h1_energy_L_noShowerRejection;
  std::map<float, std::map<int, std::map<int, std::map<int,TH1F*> > > > h1_energy_R_noShowerRejection;
  std::map<float, std::map<int, std::map<int, std::map<int,TH1F*> > > > h1_energy_LR_noShowerRejection;
  std::map<float, std::map<int, std::map<int, std::map<int,TF1*> > > > fit_energy_LR_noShowerRejection;

  std::map<float, std::map<int, std::map<int, std::map<int,TH1F*> > > > h1_tot_L;
  std::map<float, std::map<int, std::map<int, std::map<int,TH1F*> > > > h1_tot_R;
  std::map<float, std::map<int, std::map<int, std::map<int,TH1F*> > > > h1_tot_LR;
  std::map<float, std::map<int, std::map<int, std::map<int,TF1*> > > > fit_tot_LR;

  std::map<float, std::map<int, std::map<int, std::map<int,TH1F*> > > > h1_tot_L_noShowerRejection;
  std::map<float, std::map<int, std::map<int, std::map<int,TH1F*> > > > h1_tot_R_noShowerRejection;
  std::map<float, std::map<int, std::map<int, std::map<int,TH1F*> > > > h1_tot_LR_noShowerRejection;
  std::map<float, std::map<int, std::map<int, std::map<int,TF1*> > > > fit_tot_LR_noShowerRejection;

  std::map<float, std::map<int, std::map<int, std::map<int,TH1F*> > > > h1_energyRatio;
  // std::map<float, std::map<int, std::map<int, std::map<int,TH1F*> > > > h1_totRatio;

  std::map<float, std::map<int, std::map<int, std::map<int,TF1*> > > > fit_energyRatio;
  // std::map<float, std::map<int, std::map<int, std::map<int,TF1*> > > > fit_totRatio;

  std::map<float, std::map<int, std::map<int, std::map<int,TProfile*> > > > p1_deltaT_vs_energyRatio;
  std::map<float, std::map<int, std::map<int, std::map<int,TF1*> > > > fit_energyRatioCorr;
  // std::map<float, std::map<int, std::map<int, std::map<int,TProfile*> > > > p1_deltaT_vs_totRatio;
  // std::map<float, std::map<int, std::map<int, std::map<int,TF1*> > > > fit_totRatioCorr;

  std::map<float, std::map<int, std::map<int, std::map<int,TH1F*> > > > h1_deltaT_raw;
  std::map<float, std::map<int, std::map<int, std::map<int,TH1F*> > > > h1_deltaT_energyRatioCorr;
  std::map<float, std::map<int, std::map<int, std::map<int,TH1F*> > > > h1_deltaT_energyRatioPhaseCorr;
  std::map<float, std::map<int, std::map<int, std::map<int,TF1*> > > > fit_deltaT_raw;
  std::map<float, std::map<int, std::map<int, std::map<int,TF1*> > > > fit_deltaT_energyRatioCorr;
  std::map<float, std::map<int, std::map<int, std::map<int,TF1*> > > > fit_deltaT_energyRatioPhaseCorr;

  std::map<float, std::map<int, std::map<int, std::map<int,TH2F*> > > > h2_deltaT_energyRatioCorr_vs_t1fine;
  std::map<float, std::map<int, std::map<int, std::map<int,TProfile*> > > > p1_deltaT_energyRatioCorr_vs_t1fine;

  std::map<float, std::map<int, std::map<int, TH1F*> > > h1_delta_tAvg_raw;
  std::map<float, std::map<int, std::map<int, TF1*> > > fit_delta_tAvg_raw;
  std::map<float, std::map<int, std::map<int, TH1F*> > > h1_delta_tAvg;
  std::map<float, std::map<int, std::map<int, TF1*> > > fit_delta_tAvg;
  std::map<float, std::map<int, std::map<int, TH1F*> > > h1_delta_tAvg_energyCorr;
  std::map<float, std::map<int, std::map<int, TF1*> > > fit_delta_tAvg_energyCorr;
  std::map<float, std::map<int, std::map<int, TH1F*> > > h1_delta_tAvg_energyAndPhaseCorr;
  std::map<float, std::map<int, std::map<int, TF1*> > >  fit_delta_tAvg_energyAndPhaseCorr;

  std::map<float, std::map<int, std::map<int, TH1F*> > > h1_delta_tFine;
  std::map<float, std::map<int, std::map<int, TF1*> > > fit_delta_tFine;
  std::map<float, std::map<int, std::map<int, TH2F*> > > h2_t1fine_array1_vs_array0;

  std::map<float, std::map<int, std::map<int, std::map<int,TProfile*> > > > p1_delta_tAvg_vs_energy;
  std::map<float, std::map<int, std::map<int, std::map<int,TF1*> > > > fit_tAvgEnergyCorr;
  std::map<float, std::map<int, std::map<int, std::map<int,TProfile*> > > > p1_delta_tAvgEnergyCorr_vs_t1fine;
  std::map<float, std::map<int, std::map<int, std::map<int,TH2F*> > > > h2_delta_tAvgEnergyCorr_vs_t1fine;

  std::map<int, bool> pass_shower_rejection;
  std::map<int, bool> pass_mip_selection;
  std::map<int, bool> pass_energyRatio_selection;


  // 1st loop over events: get energy and tot
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%100000 == 0 )
    {
      std::cout << ">>> 1st loop - reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
    }
    pass_shower_rejection[entry] = false;
    if( entry%prescale != 0 ) continue;
    data -> GetEntry(entry);

    float Vov = step1;
    int vth1 = int(step2/10000.)-1;
    int vth2 = int((step2-10000*(vth1+1))/100.)-1;

    std::vector<bool> hit_in_array = {false, false};
    std::vector<bool> hit_in_array_noShowerRejection = {false, false};

    // only look at bar 8 to restrict spatial distribution of hits
    int iBar = 8;

    for(int iArray = 0; iArray < 2; ++iArray)
    {
      if( iArray == 1 && !hit_in_array[0] ) continue;

      // get array energy for shower rejection
      float energy_array = 0.0;
      for(unsigned int iBar = 0; iBar < num_bars; ++iBar)
      {
        int chL = ch_array_side1[iBar];
        int chR = ch_array_side2[iBar];
        if( iArray == 1 ) chL += 64;
        if( iArray == 1 ) chR += 64;

        if( channelIdx[chL] < 0 ) continue;
        if( channelIdx[chR] < 0 ) continue;

        if( (*tot)[channelIdx[chL]]/1000. <  0. ) continue;
        if( (*tot)[channelIdx[chL]]/1000. > 20. ) continue;
        if( (*tot)[channelIdx[chR]]/1000. <  0. ) continue;
        if( (*tot)[channelIdx[chR]]/1000. > 20. ) continue;

        energy_array += 0.5*((*energy)[channelIdx[chL]]+(*energy)[channelIdx[chR]]);
      }

      int chL = ch_array_side1[iBar];
      int chR = ch_array_side2[iBar];
      if( iArray == 1 ) chL += 64;
      if( iArray == 1 ) chR += 64;

      if( channelIdx[chL] < 0 ) continue;
      if( channelIdx[chR] < 0 ) continue;

      if( (*tot)[channelIdx[chL]]/1000. <  0. ) continue;
      if( (*tot)[channelIdx[chL]]/1000. > 20. ) continue;
      if( (*tot)[channelIdx[chR]]/1000. <  0. ) continue;
      if( (*tot)[channelIdx[chR]]/1000. > 20. ) continue;

      hit_in_array_noShowerRejection[iArray] = true;

      // reject shower events
      if( 0.5*((*energy)[channelIdx[chL]]+(*energy)[channelIdx[chR]]) < 0.8*energy_array ) continue;

      hit_in_array[iArray] = true;
    }

    // fill noShowerRejection plots for each array only if both have hits
    if( !hit_in_array_noShowerRejection[0] || !hit_in_array_noShowerRejection[1] ) continue;
    for(int iArray = 0; iArray < 2; ++iArray)
    {
      int chL = ch_array_side1[iBar];
      int chR = ch_array_side2[iBar];
      if( iArray == 1 ) chL += 64;
      if( iArray == 1 ) chR += 64;

      if( !h1_energy_LR_noShowerRejection[Vov][vth1][vth2][iBar+num_bars*iArray] )
      {
        outFile -> cd();

        h1_energy_L_noShowerRejection[Vov][vth1][vth2][iBar+num_bars*iArray]  = new TH1F(Form("h1_energy_L_noShowerRejection_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d", iArray,iBar,Vov,vth1,vth2),"",200,-0.5,999.5);
        h1_energy_R_noShowerRejection[Vov][vth1][vth2][iBar+num_bars*iArray]  = new TH1F(Form("h1_energy_R_noShowerRejection_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d", iArray,iBar,Vov,vth1,vth2),"",200,-0.5,999.5);
        h1_energy_LR_noShowerRejection[Vov][vth1][vth2][iBar+num_bars*iArray] = new TH1F(Form("h1_energy_LR_noShowerRejection_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d",iArray,iBar,Vov,vth1,vth2),"",200,-0.5,999.5);
        fit_energy_LR_noShowerRejection[Vov][vth1][vth2][iBar+num_bars*iArray] = new TF1(Form("fit_energy_noShowerRejection_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d", iArray,iBar,Vov,vth1,vth2),langaufun,0.,1000.,4);

        h1_tot_L_noShowerRejection[Vov][vth1][vth2][iBar+num_bars*iArray]  = new TH1F(Form("h1_tot_L_noShowerRejection_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d", iArray,iBar,Vov,vth1,vth2),"",100,0.,20.);
        h1_tot_R_noShowerRejection[Vov][vth1][vth2][iBar+num_bars*iArray]  = new TH1F(Form("h1_tot_R_noShowerRejection_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d", iArray,iBar,Vov,vth1,vth2),"",100,0.,20.);
        h1_tot_LR_noShowerRejection[Vov][vth1][vth2][iBar+num_bars*iArray] = new TH1F(Form("h1_tot_LR_noShowerRejection_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d",iArray,iBar,Vov,vth1,vth2),"",100,0.,20.);
        fit_tot_LR_noShowerRejection[Vov][vth1][vth2][iBar+num_bars*iArray] = new TF1(Form("fit_tot_noShowerRejection_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d", iArray,iBar,Vov,vth1,vth2),"gaus(0)",0.,20.);
      }

      h1_energy_L_noShowerRejection[Vov][vth1][vth2][iBar+num_bars*iArray] -> Fill( (*energy)[channelIdx[chL]] );
      h1_energy_R_noShowerRejection[Vov][vth1][vth2][iBar+num_bars*iArray] -> Fill( (*energy)[channelIdx[chR]] );
      h1_energy_LR_noShowerRejection[Vov][vth1][vth2][iBar+num_bars*iArray] -> Fill( 0.5*((*energy)[channelIdx[chL]]+(*energy)[channelIdx[chR]]) );

      h1_tot_L_noShowerRejection[Vov][vth1][vth2][iBar+num_bars*iArray] -> Fill( (*tot)[channelIdx[chL]]/1000. );
      h1_tot_R_noShowerRejection[Vov][vth1][vth2][iBar+num_bars*iArray] -> Fill( (*tot)[channelIdx[chR]]/1000. );
      h1_tot_LR_noShowerRejection[Vov][vth1][vth2][iBar+num_bars*iArray] -> Fill( 0.5*((*tot)[channelIdx[chL]]+(*tot)[channelIdx[chR]])/1000. );
    }

    // fill plots for each array only if both have non-shower hits
    if( !hit_in_array[0] || !hit_in_array[1] ) continue;
    pass_shower_rejection[entry] = true;

    for(int iArray = 0; iArray < 2; ++iArray)
    {
      int chL = ch_array_side1[iBar];
      int chR = ch_array_side2[iBar];
      if( iArray == 1 ) chL += 64;
      if( iArray == 1 ) chR += 64;

      if( !h1_energy_LR[Vov][vth1][vth2][iBar+num_bars*iArray] )
      {
        outFile -> cd();

        h1_energy_L[Vov][vth1][vth2][iBar+num_bars*iArray]  = new TH1F(Form("h1_energy_L_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d", iArray,iBar,Vov,vth1,vth2),"",200,-0.5,999.5);
        h1_energy_R[Vov][vth1][vth2][iBar+num_bars*iArray]  = new TH1F(Form("h1_energy_R_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d", iArray,iBar,Vov,vth1,vth2),"",200,-0.5,999.5);
        h1_energy_LR[Vov][vth1][vth2][iBar+num_bars*iArray] = new TH1F(Form("h1_energy_LR_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d",iArray,iBar,Vov,vth1,vth2),"",200,-0.5,999.5);

        h1_tot_L[Vov][vth1][vth2][iBar+num_bars*iArray]  = new TH1F(Form("h1_tot_L_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d", iArray,iBar,Vov,vth1,vth2),"",100,0.,20.);
        h1_tot_R[Vov][vth1][vth2][iBar+num_bars*iArray]  = new TH1F(Form("h1_tot_R_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d", iArray,iBar,Vov,vth1,vth2),"",100,0.,20.);
        h1_tot_LR[Vov][vth1][vth2][iBar+num_bars*iArray] = new TH1F(Form("h1_tot_LR_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d",iArray,iBar,Vov,vth1,vth2),"",100,0.,20.);

        fit_energy_LR[Vov][vth1][vth2][iBar+num_bars*iArray] = new TF1(Form("fit_energy_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d", iArray,iBar,Vov,vth1,vth2),langaufun,0.,1000.,4);
        fit_tot_LR[Vov][vth1][vth2][iBar+num_bars*iArray] = new TF1(Form("fit_tot_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d", iArray,iBar,Vov,vth1,vth2),"gaus(0)",0.,20.);
      }

      h1_energy_L[Vov][vth1][vth2][iBar+num_bars*iArray] -> Fill( (*energy)[channelIdx[chL]] );
      h1_energy_R[Vov][vth1][vth2][iBar+num_bars*iArray] -> Fill( (*energy)[channelIdx[chR]] );
      h1_energy_LR[Vov][vth1][vth2][iBar+num_bars*iArray] -> Fill( 0.5*((*energy)[channelIdx[chL]]+(*energy)[channelIdx[chR]]) );

      h1_tot_L[Vov][vth1][vth2][iBar+num_bars*iArray] -> Fill( (*tot)[channelIdx[chL]]/1000. );
      h1_tot_R[Vov][vth1][vth2][iBar+num_bars*iArray] -> Fill( (*tot)[channelIdx[chR]]/1000. );
      h1_tot_LR[Vov][vth1][vth2][iBar+num_bars*iArray] -> Fill( 0.5*((*tot)[channelIdx[chL]]+(*tot)[channelIdx[chR]])/1000. );
    }
  }
  std::cout << std::endl;

  // draw 1st plots
  std::map<float, std::map<int, std::map<int, std::map<int, TGraphErrors*> > > >  g_energyPeak;

  for(auto mapIt : h1_energy_LR )
  {
    float Vov = mapIt.first;
    for(auto mapIt2 : mapIt.second)
    {
      int vth1 = mapIt2.first;

      for(auto mapIt3 : mapIt2.second)
      {
        int vth2 = mapIt3.first;

        for(int iArray = 0; iArray < 2; ++iArray)
        {
          for(unsigned int iBar = 8; iBar < 9; ++iBar)
          {
            drawH1_energy(h1_energy_LR[Vov][vth1][vth2][iBar+num_bars*iArray], h1_energy_L[Vov][vth1][vth2][iBar+num_bars*iArray], h1_energy_R[Vov][vth1][vth2][iBar+num_bars*iArray], 0., 1000, "energy [ADC]", "events", plotDir+"energy", fit_energy_LR[Vov][vth1][vth2][iBar+num_bars*iArray]);
            drawH1_energy(h1_tot_LR[Vov][vth1][vth2][iBar+num_bars*iArray], h1_tot_L[Vov][vth1][vth2][iBar+num_bars*iArray], h1_tot_R[Vov][vth1][vth2][iBar+num_bars*iArray], 0., 1000, "ToT [ns]", "events", plotDir+"tot", fit_tot_LR[Vov][vth1][vth2][iBar+num_bars*iArray]);

            drawH1_energy(h1_energy_LR_noShowerRejection[Vov][vth1][vth2][iBar+num_bars*iArray], h1_energy_L_noShowerRejection[Vov][vth1][vth2][iBar+num_bars*iArray], h1_energy_R_noShowerRejection[Vov][vth1][vth2][iBar+num_bars*iArray], 0., 1000, "energy [ADC]", "events", plotDir+"energy", fit_energy_LR_noShowerRejection[Vov][vth1][vth2][iBar+num_bars*iArray]);
            drawH1_energy(h1_tot_LR_noShowerRejection[Vov][vth1][vth2][iBar+num_bars*iArray], h1_tot_L_noShowerRejection[Vov][vth1][vth2][iBar+num_bars*iArray], h1_tot_R_noShowerRejection[Vov][vth1][vth2][iBar+num_bars*iArray], 0., 1000, "ToT [ns]", "events", plotDir+"tot", fit_tot_LR_noShowerRejection[Vov][vth1][vth2][iBar+num_bars*iArray]);

            if( fit_energy_LR[Vov][vth1][vth2][iBar+num_bars*iArray] == NULL ) continue;

            if( !g_energyPeak[Vov][vth1][vth2][iArray] )
            {
              g_energyPeak[Vov][vth1][vth2][iArray] = new TGraphErrors();
              g_energyPeak[Vov][vth1][vth2][iArray] -> SetName(Form("g_energyPeak_array%d_Vov%.1f_vth1_%02d_vth2_%02d",iArray,Vov,vth1,vth2));
            }
            g_energyPeak[Vov][vth1][vth2][iArray] -> SetPoint(g_energyPeak[Vov][vth1][vth2][iArray]->GetN(),iBar,fit_energy_LR[Vov][vth1][vth2][iBar+num_bars*iArray]->GetParameter(1));
          }
        }
      }
    }
  }

  for(auto mapIt : h1_energy_LR )
  {
    float Vov = mapIt.first;
    for(auto mapIt2 : mapIt.second)
    {
      int vth1 = mapIt2.first;

      for(auto mapIt3 : mapIt2.second)
      {
        int vth2 = mapIt3.first;

        drawG_arrays(g_energyPeak[Vov][vth1][vth2][0],g_energyPeak[Vov][vth1][vth2][1], "iBar", "MIP peak [ADC]", 0., 800., plotDir+"energy",false,Form("g_energyPeak_Vov%.1f_vth1_%02d_vth2_%02d",Vov,vth1,vth2));
      }
    }
  }

  // 2nd loop over events: get energy and tot ratios; calculate delta_tAvg_raw
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%100000 == 0 )
    {
      std::cout << ">>> 2nd loop - reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
    }
    if( entry%prescale != 0 ) continue;
    if ( !pass_shower_rejection[entry] ) continue;
    data -> GetEntry(entry);

    float Vov = step1;
    int vth1 = int(step2/10000.)-1;
    int vth2 = int((step2-10000*(vth1+1))/100.)-1;

    // only look at bar 8 to restrict spatial distribution of hits
    int iBar = 8;

    // apply MIP selection
    pass_mip_selection[entry] = false;
    std::vector<bool> pass_mip_selection_array = {false, false};
    for(int iArray = 0; iArray < 2; ++iArray)
    {
      int chL = ch_array_side1[iBar];
      int chR = ch_array_side2[iBar];
      if( iArray == 1 ) chL += 64;
      if( iArray == 1 ) chR += 64;

      float energyMean = 0.5*((*energy)[channelIdx[chL]]+(*energy)[channelIdx[chR]]);
      TF1* func = fit_energy_LR[Vov][vth1][vth2][iBar+num_bars*iArray];
      pass_mip_selection_array[iArray] = ( energyMean > 0.80*func->GetParameter(1) );
    }

    if( !pass_mip_selection_array[0] || !pass_mip_selection_array[1] ) continue;
    pass_mip_selection[entry] = true;
    
    // for delta_tAvg calculation
    // perhaps I should do something to ensure that these default values are never used
    std::vector<long double> tAvg = {0, 0};

    for(int iArray = 0; iArray < 2; ++iArray)
    {
      int chL = ch_array_side1[iBar];
      int chR = ch_array_side2[iBar];
      if( iArray == 1 ) chL += 64;
      if( iArray == 1 ) chR += 64;

      if( !h1_energyRatio[Vov][vth1][vth2][iBar+num_bars*iArray] )
      {
        outFile -> cd();

        h1_energyRatio[Vov][vth1][vth2][iBar+num_bars*iArray] = new TH1F(Form("h1_energyRatio_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d",iArray,iBar,Vov,vth1,vth2),"",500,0.5,1.5);
        // h1_totRatio[Vov][vth1][vth2][iBar+num_bars*iArray] = new TH1F(Form("h1_totRatio_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d",iArray,iBar,Vov,vth1,vth2),"",500,0.5,1.5);

        fit_energyRatio[Vov][vth1][vth2][iBar+num_bars*iArray] = new TF1(Form("fit_energyRatio_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d", iArray,iBar,Vov,vth1,vth2),"gaus(0)",0.,1000.);
        // fit_totRatio[Vov][vth1][vth2][iBar+num_bars*iArray] = new TF1(Form("fit_totRatio_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d", iArray,iBar,Vov,vth1,vth2),"gaus(0)",0.,20.);
      }

      h1_energyRatio[Vov][vth1][vth2][iBar+num_bars*iArray] -> Fill( (*energy)[channelIdx[chL]]/(*energy)[channelIdx[chR]] );
      // h1_totRatio[Vov][vth1][vth2][iBar+num_bars*iArray] -> Fill( (*tot)[channelIdx[chL]]/(*tot)[channelIdx[chR]] );

      tAvg[iArray] = 0.5*((*time)[channelIdx[chL]] + (*time)[channelIdx[chR]]);
    }

    if(!h1_delta_tAvg_raw[Vov][vth1][vth2])
    {
      h1_delta_tAvg_raw[Vov][vth1][vth2] = new TH1F(Form("h1_delta_tAvg_raw_Vov%.1f_vth1_%02d_vth2_%02d",Vov,vth1,vth2),"",100,-2000.,0.);
      fit_delta_tAvg_raw[Vov][vth1][vth2] = new TF1(Form("fit_delta_tAvg_raw_Vov%.1f_vth1_%02d_vth2_%02d", Vov,vth1,vth2),"gaus(0)",-2000.,0);
    }

    h1_delta_tAvg_raw[Vov][vth1][vth2] -> Fill( tAvg[1] - tAvg[0] );
  }
  std::cout << std::endl;

  // draw 2nd plots
  for(auto mapIt : h1_delta_tAvg_raw )
  {
    float Vov = mapIt.first;
    for(auto mapIt2 : mapIt.second)
    {
      int vth1 = mapIt2.first;

      for(auto mapIt3 : mapIt2.second)
      {
        int vth2 = mapIt3.first;

        drawH1_fitGaus(h1_delta_tAvg_raw[Vov][vth1][vth2], -2000., 0., "#Delta_{t_{avg}} [ps]", "events", plotDir+"delta_tAvg", fit_delta_tAvg_raw[Vov][vth1][vth2]);
      }
    }
  }

  for(auto mapIt : h1_energyRatio )
  {
    float Vov = mapIt.first;
    for(auto mapIt2 : mapIt.second)
    {
      int vth1 = mapIt2.first;

      for(auto mapIt3 : mapIt2.second)
      {
        int vth2 = mapIt3.first;

        for(int iArray = 0; iArray < 2; ++iArray)
        {
          for(unsigned int iBar = 8; iBar < 9; ++iBar)
          {
            drawH1_fitGaus(h1_energyRatio[Vov][vth1][vth2][iBar+num_bars*iArray], 0., 3., "energy ratio", "events", plotDir+"energyRatioCorr", fit_energyRatio[Vov][vth1][vth2][iBar+num_bars*iArray], false);

            // drawH1_fitGaus(h1_totRatio[Vov][vth1][vth2][iBar+num_bars*iArray], 0., 3., "ToT ratio", "events", plotDir+"totRatioCorr", fit_totRatio[Vov][vth1][vth2][iBar+num_bars*iArray], false);
          }
        }
      }
    }
  }

  // 3rd loop over events: derive energy ratio correction to deltaT
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%100000 == 0 )
    {
      std::cout << ">>> 3rd loop - reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
    }
    if( entry%prescale != 0 ) continue;
    if ( !pass_shower_rejection[entry] ) continue;
    if ( !pass_mip_selection[entry] ) continue;
    data -> GetEntry(entry);

    float Vov = step1;
    int vth1 = int(step2/10000.)-1;
    int vth2 = int((step2-10000*(vth1+1))/100.)-1;

    // only look at bar 8 to restrict spatial distribution of hits
    int iBar = 8;

    // for delta_tAvg calculation
    // perhaps I should do something to ensure that these default values are never used
    std::vector<long double> tAvg = {0, 0};
    std::vector<float> energyMean = {0., 0.};

    for(int iArray = 0; iArray < 2; ++iArray)
    {
      int chL = ch_array_side1[iBar];
      int chR = ch_array_side2[iBar];
      if( iArray == 1 ) chL += 64;
      if( iArray == 1 ) chR += 64;

      if( !p1_deltaT_vs_energyRatio[Vov][vth1][vth2][iBar+num_bars*iArray] )
      {
        outFile -> cd();

        float mean = fit_energyRatio[Vov][vth1][vth2][iBar+num_bars*iArray]->GetParameter(1);
        float sigma = fit_energyRatio[Vov][vth1][vth2][iBar+num_bars*iArray]->GetParameter(2);
        p1_deltaT_vs_energyRatio[Vov][vth1][vth2][iBar+num_bars*iArray] = new TProfile(Form("p1_deltaT_vs_energyRatio_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d",iArray,iBar,Vov,vth1,vth2),"",100,mean-3.*sigma,mean+3.*sigma);
        fit_energyRatioCorr[Vov][vth1][vth2][iBar+num_bars*iArray] = new TF1(Form("fit_energyRatioCorr_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d", iArray,iBar,Vov,vth1,vth2),"pol3",0.,5.);

        // mean = fit_totRatio[Vov][vth1][vth2][iBar+num_bars*iArray]->GetParameter(1);
        // sigma = fit_totRatio[Vov][vth1][vth2][iBar+num_bars*iArray]->GetParameter(2);
        // p1_deltaT_vs_totRatio[Vov][vth1][vth2][iBar+num_bars*iArray] = new TProfile(Form("p1_deltaT_vs_totRatio_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d",iArray,iBar,Vov,vth1,vth2),"",100,mean-3.*sigma,mean+3.*sigma);
        // fit_totRatioCorr[Vov][vth1][vth2][iBar+num_bars*iArray] = new TF1(Form("fit_totRatioCorr_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d", iArray,iBar,Vov,vth1,vth2),"pol3",0.,5.);
      }

      float energyRatio = (*energy)[channelIdx[chL]]/(*energy)[channelIdx[chR]];
      TF1* func = fit_energyRatio[Vov][vth1][vth2][iBar+num_bars*iArray];
      // if( fabs(energyRatio-func->GetParameter(1)) > 2.*func->GetParameter(2) ) continue;

      // float totRatio = (*tot)[channelIdx[chL]]/(*tot)[channelIdx[chR]];
      // func = fit_totRatio[Vov][vth1][vth2][iBar+num_bars*iArray];
      // if( fabs(totRatio-func->GetParameter(1)) > 2.*func->GetParameter(2) ) continue;

      float deltaT = (*time)[channelIdx[chL]] - (*time)[channelIdx[chR]];
      // if( fabs(deltaT) > 2000. ) continue;

      p1_deltaT_vs_energyRatio[Vov][vth1][vth2][iBar+num_bars*iArray] -> Fill( energyRatio,deltaT );
      // p1_deltaT_vs_totRatio[Vov][vth1][vth2][iBar+num_bars*iArray] -> Fill( totRatio,deltaT );

      energyMean[iArray] = 0.5*((*energy)[channelIdx[chL]]+(*energy)[channelIdx[chR]]);
      tAvg[iArray] = 0.5*((*time)[channelIdx[chL]] + (*time)[channelIdx[chR]]);
    }

    if(!p1_delta_tAvg_vs_energy[Vov][vth1][vth2][iBar+num_bars*0])
    {
      p1_delta_tAvg_vs_energy[Vov][vth1][vth2][iBar+num_bars*0] = new TProfile(Form("p1_delta_tAvg_vs_energy_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d",0,iBar,Vov,vth1,vth2),"",100,0.,1000.);
      p1_delta_tAvg_vs_energy[Vov][vth1][vth2][iBar+num_bars*1] = new TProfile(Form("p1_delta_tAvg_vs_energy_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d",1,iBar,Vov,vth1,vth2),"",100,0.,1000.);
      fit_tAvgEnergyCorr[Vov][vth1][vth2][iBar+num_bars*0] = new TF1(Form("fit_tAvgEnergyCorr_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d", 0,iBar,Vov,vth1,vth2),"pol3",0.,1000.);
      fit_tAvgEnergyCorr[Vov][vth1][vth2][iBar+num_bars*1] = new TF1(Form("fit_tAvgEnergyCorr_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d", 1,iBar,Vov,vth1,vth2),"pol3",0.,1000.);
    }

    TF1* func_delta_tAvg_raw = fit_delta_tAvg_raw[Vov][vth1][vth2];
    int delta_tAvg = tAvg[1] - tAvg[0] - func_delta_tAvg_raw->GetParameter(1);

    if( fabs(delta_tAvg) > 2000. ) continue;

    p1_delta_tAvg_vs_energy[Vov][vth1][vth2][iBar+num_bars*0] -> Fill( energyMean[0], delta_tAvg );
    p1_delta_tAvg_vs_energy[Vov][vth1][vth2][iBar+num_bars*1] -> Fill( energyMean[1], delta_tAvg );
  }
  std::cout << std::endl;

  // draw 3rd plots
  for(auto mapIt : p1_deltaT_vs_energyRatio )
  {
    float Vov = mapIt.first;
    for(auto mapIt2 : mapIt.second)
    {
      int vth1 = mapIt2.first;

      for(auto mapIt3 : mapIt2.second)
      {
        int vth2 = mapIt3.first;

        for(int iArray = 0; iArray < 2; ++iArray)
        {
          for(unsigned int iBar = 8; iBar < 9; ++iBar)
          {
            drawP1_fitPol(p1_deltaT_vs_energyRatio[Vov][vth1][vth2][iBar+num_bars*iArray], "energy ratio", "#DeltaT [ps]", plotDir+"energyRatioCorr", fit_energyRatioCorr[Vov][vth1][vth2][iBar+num_bars*iArray]);
            drawP1_fitPol(p1_delta_tAvg_vs_energy[Vov][vth1][vth2][iBar+num_bars*iArray], "Energy [ADC]", "#Delta_{t_{avg}} [ps]", plotDir+"energyRatioCorr", fit_tAvgEnergyCorr[Vov][vth1][vth2][iBar+num_bars*iArray]);

            // drawP1_fitPol(p1_deltaT_vs_totRatio[Vov][vth1][vth2][iBar+num_bars*iArray], "ToT ratio", "#DeltaT [ps]", plotDir+"totRatioCorr", fit_totRatioCorr[Vov][vth1][vth2][iBar+num_bars*iArray]);
          }
        }
      }
    }
  }

  // 4th loop over events: derive phase correction to deltaT
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%100000 == 0 )
    {
      std::cout << ">>> 4th loop - reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
    }
    if( entry%prescale != 0 ) continue;
    if ( !pass_shower_rejection[entry] ) continue;
    if ( !pass_mip_selection[entry] ) continue;
    data -> GetEntry(entry);

    float Vov = step1;
    int vth1 = int(step2/10000.)-1;
    int vth2 = int((step2-10000*(vth1+1))/100.)-1;

    // only look at bar 8 to restrict spatial distribution of hits
    int iBar = 8;

    // apply energyRatio selection
    pass_energyRatio_selection[entry] = false;
    std::vector<bool> pass_energyRatio_selection_array = {false, false};
    for(int iArray = 0; iArray < 2; ++iArray)
    {
      int chL = ch_array_side1[iBar];
      int chR = ch_array_side2[iBar];
      if( iArray == 1 ) chL += 64;
      if( iArray == 1 ) chR += 64;

      float energyRatio = (*energy)[channelIdx[chL]]/(*energy)[channelIdx[chR]];
      TF1* func_energyRatio = fit_energyRatio[Vov][vth1][vth2][iBar+num_bars*iArray];
      pass_energyRatio_selection_array[iArray] = ( fabs(energyRatio-func_energyRatio->GetParameter(1)) < 2.*func_energyRatio->GetParameter(2) );
    }

    if( !pass_energyRatio_selection_array[0] || !pass_energyRatio_selection_array[1] ) continue;
    pass_energyRatio_selection[entry] = true;

    std::vector<long double> tAvg = {0, 0};
    std::vector<long double> tAvgCorr = {0, 0};
    std::vector<float> tFine = {0., 0.};

    for(int iArray = 0; iArray < 2; ++iArray)
    {
      int chL = ch_array_side1[iBar];
      int chR = ch_array_side2[iBar];
      if( iArray == 1 ) chL += 64;
      if( iArray == 1 ) chR += 64;

      if( !h1_deltaT_raw[Vov][vth1][vth2][iBar+num_bars*iArray] )
      {
        outFile -> cd();

        h1_deltaT_raw[Vov][vth1][vth2][iBar+num_bars*iArray] = new TH1F(Form("h1_deltaT_raw_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d",iArray,iBar,Vov,vth1,vth2),"",500,-5000.,5000.);
        h1_deltaT_energyRatioCorr[Vov][vth1][vth2][iBar+num_bars*iArray] = new TH1F(Form("h1_deltaT_energyRatioCorr_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d",iArray,iBar,Vov,vth1,vth2),"",500,-5000.,5000.);

        fit_deltaT_raw[Vov][vth1][vth2][iBar+num_bars*iArray] = new TF1(Form("fit_deltaT_raw_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d", iArray,iBar,Vov,vth1,vth2),"gaus(0)",-5000.,5000.);
        fit_deltaT_energyRatioCorr[Vov][vth1][vth2][iBar+num_bars*iArray] = new TF1(Form("fit_deltaT_energyRatioCorr_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d", iArray,iBar,Vov,vth1,vth2),"gaus(0)",-5000.,5000.);

        p1_deltaT_energyRatioCorr_vs_t1fine[Vov][vth1][vth2][iBar+num_bars*iArray] = new TProfile(Form("p1_deltaT_energyRatioCorr_vs_t1fine_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d",iArray,iBar,Vov,vth1,vth2),"",128,-0.5,1023.5);
      }

      float deltaT = (*time)[channelIdx[chL]] - (*time)[channelIdx[chR]];
      float energyRatio = (*energy)[channelIdx[chL]]/(*energy)[channelIdx[chR]];
      TF1* func_energyRatio = fit_energyRatio[Vov][vth1][vth2][iBar+num_bars*iArray];
      TF1* func_energyRatioCorr = fit_energyRatioCorr[Vov][vth1][vth2][iBar+num_bars*iArray];
      float deltaT_energyRatioCorr = deltaT - func_energyRatioCorr->Eval(energyRatio) + func_energyRatioCorr->Eval(func_energyRatio->GetParameter(1));

      h1_deltaT_raw[Vov][vth1][vth2][iBar+num_bars*iArray] -> Fill( deltaT );
      h1_deltaT_energyRatioCorr[Vov][vth1][vth2][iBar+num_bars*iArray] -> Fill( deltaT_energyRatioCorr );

      tFine[iArray] =  0.5*((*t1fine)[channelIdx[chL]]+(*t1fine)[channelIdx[chR]]);
      if( fabs(deltaT_energyRatioCorr) < 1000 )
        p1_deltaT_energyRatioCorr_vs_t1fine[Vov][vth1][vth2][iBar+num_bars*iArray] -> Fill( tFine[iArray],deltaT_energyRatioCorr );

      float energyMean = 0.5*((*energy)[channelIdx[chL]] + (*energy)[channelIdx[chR]]);
      TF1* func_tAvgEnergyCorr = fit_tAvgEnergyCorr[Vov][vth1][vth2][iBar+num_bars*iArray];
      tAvg[iArray] = 0.5*((*time)[channelIdx[chL]] + (*time)[channelIdx[chR]]);
      tAvgCorr[iArray] = func_tAvgEnergyCorr->Eval(energyMean);
    }

    if(!h1_delta_tAvg[Vov][vth1][vth2])
    {
      h1_delta_tAvg[Vov][vth1][vth2] = new TH1F(Form("h1_delta_tAvg_Vov%.1f_vth1_%02d_vth2_%02d",Vov,vth1,vth2),"",100,-1000.,1000.);
      fit_delta_tAvg[Vov][vth1][vth2] = new TF1(Form("fit_delta_tAvg_Vov%.1f_vth1_%02d_vth2_%02d", Vov,vth1,vth2),"gaus(0)",-1000.,1000.);

      h1_delta_tAvg_energyCorr[Vov][vth1][vth2] = new TH1F(Form("h1_delta_tAvg_energyCorr_Vov%.1f_vth1_%02d_vth2_%02d",Vov,vth1,vth2),"",100,-1000.,1000.);
      fit_delta_tAvg_energyCorr[Vov][vth1][vth2] = new TF1(Form("fit_delta_tAvg_energyCorr_Vov%.1f_vth1_%02d_vth2_%02d", Vov,vth1,vth2),"gaus(0)",-1000.,1000.);

      h1_delta_tFine[Vov][vth1][vth2] = new TH1F(Form("h1_delta_tFine_Vov%.1f_vth1_%02d_vth2_%02d",Vov,vth1,vth2),"",100,-1000.,1000.);
      fit_delta_tFine[Vov][vth1][vth2] = new TF1(Form("fit_delta_tFine_Vov%.1f_vth1_%02d_vth2_%02d", Vov,vth1,vth2),"gaus(0)",-1000.,1000.);

      h2_delta_tAvgEnergyCorr_vs_t1fine[Vov][vth1][vth2][iBar+num_bars*0] = new TH2F(Form("h2_delta_tAvgEnergyCorr_vs_t1fine_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d",0,iBar,Vov,vth1,vth2),"",128,0,1024,100,-1000.,1000.);
      h2_delta_tAvgEnergyCorr_vs_t1fine[Vov][vth1][vth2][iBar+num_bars*1] = new TH2F(Form("h2_delta_tAvgEnergyCorr_vs_t1fine_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d",1,iBar,Vov,vth1,vth2),"",128,0,1024,100,-1000.,1000.);
      h2_t1fine_array1_vs_array0[Vov][vth1][vth2] = new TH2F(Form("h2_t1fine_array1_vs_array0_Vov%.1f_vth1_%02d_vth2_%02d",Vov,vth1,vth2),"",128,0,1024,128,0,1024);

      p1_delta_tAvgEnergyCorr_vs_t1fine[Vov][vth1][vth2][iBar+num_bars*0] = new TProfile(Form("p1_delta_tAvgEnergyCorr_vs_t1fine_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d",0,iBar,Vov,vth1,vth2),"",128,-0.5,1023.5);
      p1_delta_tAvgEnergyCorr_vs_t1fine[Vov][vth1][vth2][iBar+num_bars*1] = new TProfile(Form("p1_delta_tAvgEnergyCorr_vs_t1fine_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d",1,iBar,Vov,vth1,vth2),"",128,-0.5,1023.5);
    }

    TF1* func_delta_tAvg_raw = fit_delta_tAvg_raw[Vov][vth1][vth2];
    int delta_tAvg = tAvg[1] - tAvg[0] - func_delta_tAvg_raw->GetParameter(1);
    h1_delta_tAvg[Vov][vth1][vth2] -> Fill( delta_tAvg );
    float delta_tAvg_energyCorr =  delta_tAvg - tAvgCorr[0] - tAvgCorr[1];
    h1_delta_tAvg_energyCorr[Vov][vth1][vth2] -> Fill( delta_tAvg_energyCorr );

    h1_delta_tFine[Vov][vth1][vth2] -> Fill( tFine[1] - tFine[0] );

    h2_delta_tAvgEnergyCorr_vs_t1fine[Vov][vth1][vth2][iBar+num_bars*0] -> Fill( tFine[0], delta_tAvg_energyCorr );
    h2_delta_tAvgEnergyCorr_vs_t1fine[Vov][vth1][vth2][iBar+num_bars*1] -> Fill( tFine[1], delta_tAvg_energyCorr );
    h2_t1fine_array1_vs_array0[Vov][vth1][vth2] -> Fill( tFine[0], tFine[1] );

    if( fabs(delta_tAvg_energyCorr) < 1000 )
    {
        p1_delta_tAvgEnergyCorr_vs_t1fine[Vov][vth1][vth2][iBar+num_bars*0] -> Fill( tFine[0], delta_tAvg_energyCorr );
        p1_delta_tAvgEnergyCorr_vs_t1fine[Vov][vth1][vth2][iBar+num_bars*1] -> Fill( tFine[1], delta_tAvg_energyCorr );
    }
  }
  std::cout << std::endl;

  // draw 4th plots
  std::map<float, std::map<int, std::map<int, std::map<int, TGraphErrors*> > > >  g_tRes_energyRatioCorr;
  std::map<float, std::map<int, std::map<int, std::map<int, TGraphErrors*> > > > g_tRes_energyRatioCorr_vs_vth1;

  for(auto mapIt : h1_deltaT_raw )
  {
    float Vov = mapIt.first;
    for(auto mapIt2 : mapIt.second)
    {
      int vth1 = mapIt2.first;

      for(auto mapIt3 : mapIt2.second)
      {
        int vth2 = mapIt3.first;

        drawH1_fitGaus(h1_delta_tAvg[Vov][vth1][vth2], -1000., 1000., "#Delta_{t_{avg}} [ps]", "events", plotDir+"delta_tAvg", fit_delta_tAvg[Vov][vth1][vth2]);

        drawH1_fitGaus(h1_delta_tAvg_energyCorr[Vov][vth1][vth2], -1000., 1000., "#Delta_{t_{avg}} [ps]", "events", plotDir+"delta_tAvg", fit_delta_tAvg_energyCorr[Vov][vth1][vth2]);

        drawH1_deltaT(h1_delta_tAvg[Vov][vth1][vth2], h1_delta_tAvg_energyCorr[Vov][vth1][vth2], "#Delta_{t_{avg}} [ps]", "events", plotDir+"delta_tAvg", fit_delta_tAvg[Vov][vth1][vth2], fit_delta_tAvg_energyCorr[Vov][vth1][vth2]);

        drawH1_fitGaus(h1_delta_tFine[Vov][vth1][vth2], -1000., 1000., "#Delta_{t1fine}", "events", plotDir+"phaseCorr", fit_delta_tFine[Vov][vth1][vth2]);

        drawH2(h2_t1fine_array1_vs_array0[Vov][vth1][vth2], "t1fine, array 0", "t1fine, array 1", "events", plotDir+"phaseCorr");

        for(int iArray = 0; iArray < 2; ++iArray)
        {
          for(unsigned int iBar = 8; iBar < 9; ++iBar)
          {
            drawH1_deltaT(h1_deltaT_raw[Vov][vth1][vth2][iBar+num_bars*iArray],h1_deltaT_energyRatioCorr[Vov][vth1][vth2][iBar+num_bars*iArray], "#Deltat [ps]", "events", plotDir+"deltaT",fit_deltaT_raw[Vov][vth1][vth2][iBar+num_bars*iArray],fit_deltaT_energyRatioCorr[Vov][vth1][vth2][iBar+num_bars*iArray]);

            drawH2(h2_delta_tAvgEnergyCorr_vs_t1fine[Vov][vth1][vth2][iBar+num_bars*iArray], "t1fine", "#Delta_{t_{avg}} [ps]", "events", plotDir+"phaseCorr");

            drawP1(p1_deltaT_energyRatioCorr_vs_t1fine[Vov][vth1][vth2][iBar+num_bars*iArray], "t1fine", "#Deltat [ps]", plotDir+"phaseCorr");
            drawP1(p1_delta_tAvgEnergyCorr_vs_t1fine[Vov][vth1][vth2][iBar+num_bars*iArray], "t1fine", "#Delta_{t_{avg}} [ps]", plotDir+"phaseCorr");
          }
        }
      }
    }
  }

  // 5th loop over events: derive time resolution from width of fully corrected tDiff distribution
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%100000 == 0 )
    {
      std::cout << ">>> 5th loop - reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
    }
    if( entry%prescale != 0 ) continue;
    if ( !pass_shower_rejection[entry] ) continue;
    if ( !pass_mip_selection[entry] ) continue;
    if ( !pass_energyRatio_selection[entry] ) continue;
    data -> GetEntry(entry);

    float Vov = step1;
    int vth1 = int(step2/10000.)-1;
    int vth2 = int((step2-10000*(vth1+1))/100.)-1;

    // only look at bar 8 to restrict spatial distribution of hits
    int iBar = 8;

    std::vector<float> tFine = {0., 0.};
    std::vector<long double> tAvg = {0, 0};
    std::vector<long double> tAvgCorr = {0, 0};
    std::vector<long double> tAvgPhaseCorr = {0, 0};

    for(int iArray = 0; iArray < 2; ++iArray)
    {
      int chL = ch_array_side1[iBar];
      int chR = ch_array_side2[iBar];
      if( iArray == 1 ) chL += 64;
      if( iArray == 1 ) chR += 64;

      if( !h1_deltaT_energyRatioPhaseCorr[Vov][vth1][vth2][iBar+num_bars*iArray] )
      {
        outFile -> cd();

        h1_deltaT_energyRatioPhaseCorr[Vov][vth1][vth2][iBar+num_bars*iArray] = new TH1F(Form("h1_deltaT_energyRatioPhaseCorr_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d",iArray,iBar,Vov,vth1,vth2),"",500,-5000.,5000.);
        fit_deltaT_energyRatioPhaseCorr[Vov][vth1][vth2][iBar+num_bars*iArray] = new TF1(Form("fit_deltaT_energyRatioPhaseCorr_array%d_bar%02i_Vov%.1f_vth1_%02d_vth2_%02d", iArray,iBar,Vov,vth1,vth2),"gaus(0)",-5000.,5000.);
      }

      float energyRatio = (*energy)[channelIdx[chL]]/(*energy)[channelIdx[chR]];
      TF1* func_energyRatioCorr = fit_energyRatioCorr[Vov][vth1][vth2][iBar+num_bars*iArray];
      TF1* func_energyRatio = fit_energyRatio[Vov][vth1][vth2][iBar+num_bars*iArray];
      float deltaT = (*time)[channelIdx[chL]] - (*time)[channelIdx[chR]];
      float deltaT_energyRatioCorr = deltaT - func_energyRatioCorr->Eval(energyRatio) + func_energyRatioCorr->Eval(func_energyRatio->GetParameter(1));

      float energyMean = 0.5*((*energy)[channelIdx[chL]] + (*energy)[channelIdx[chR]]);
      TF1* func_tAvgEnergyCorr = fit_tAvgEnergyCorr[Vov][vth1][vth2][iBar+num_bars*iArray];
      tAvg[iArray] = 0.5*((*time)[channelIdx[chL]] + (*time)[channelIdx[chR]]);
      tAvgCorr[iArray] = func_tAvgEnergyCorr->Eval(energyMean);

      tFine[iArray] =  0.5*((*t1fine)[channelIdx[chL]]+(*t1fine)[channelIdx[chR]]);
      TProfile* prof = p1_deltaT_energyRatioCorr_vs_t1fine[Vov][vth1][vth2][iBar+num_bars*iArray];
      h1_deltaT_energyRatioPhaseCorr[Vov][vth1][vth2][iBar+num_bars*iArray] -> Fill( deltaT_energyRatioCorr - prof->GetBinContent(prof->FindBin(tFine[iArray])) );

      prof = p1_delta_tAvgEnergyCorr_vs_t1fine[Vov][vth1][vth2][iBar+num_bars*iArray];
      tAvgPhaseCorr[iArray] = prof->GetBinContent(prof->FindBin(tFine[iArray]));
    }

    if(!h1_delta_tAvg_energyAndPhaseCorr[Vov][vth1][vth2])
    {
      h1_delta_tAvg_energyAndPhaseCorr[Vov][vth1][vth2] = new TH1F(Form("h1_delta_tAvg_energyAndPhaseCorr_Vov%.1f_vth1_%02d_vth2_%02d",Vov,vth1,vth2),"",100,-1000.,1000.);
      fit_delta_tAvg_energyAndPhaseCorr[Vov][vth1][vth2] = new TF1(Form("fit_delta_tAvg_energyAndPhaseCorr_Vov%.1f_vth1_%02d_vth2_%02d", Vov,vth1,vth2),"gaus(0)",-1000.,1000.);
    }

    TF1* func_delta_tAvg_raw = fit_delta_tAvg_raw[Vov][vth1][vth2];
    float delta_tAvg = tAvg[1] - tAvg[0] - func_delta_tAvg_raw->GetParameter(1);
    float delta_tAvg_energyCorr =  delta_tAvg - tAvgCorr[0] - tAvgCorr[1];
    float tAvg_energyAndPhaseCorr =  delta_tAvg_energyCorr - tAvgPhaseCorr[0] - tAvgPhaseCorr[1];
    h1_delta_tAvg_energyAndPhaseCorr[Vov][vth1][vth2] -> Fill( tAvg_energyAndPhaseCorr );
  }
  std::cout << std::endl;

  // draw 5th plots
  std::map<float, std::map<int, std::map<int, std::map<int, TGraphErrors*> > > >  g_tRes_energyRatioPhaseCorr;
  std::map<float, std::map<int, std::map<int, std::map<int, TGraphErrors*> > > > g_tRes_energyRatioPhaseCorr_vs_vth1;
  std::map<float, std::map<int, std::map<int, std::map<int, TGraphErrors*> > > > g_deltaT_actualRes_vs_vth1;
  std::map<float, std::map<int, std::map<int, TGraphErrors*> > > g_tAvgRes_energyAndPhaseCorr_vs_vth1;
  std::map<float, std::map<int, std::map<int, TGraphErrors*> > > g_tAvg_actualRes_vs_vth1;

  for(auto mapIt : h1_deltaT_raw )
  {
    float Vov = mapIt.first;
    for(auto mapIt2 : mapIt.second)
    {
      int vth1 = mapIt2.first;

      for(auto mapIt3 : mapIt2.second)
      {
        int vth2 = mapIt3.first;

        for(unsigned int iBar = 8; iBar < 9; ++iBar)
        {

          drawH1_fitGaus(h1_delta_tAvg_energyAndPhaseCorr[Vov][vth1][vth2], -1000., 1000., "#Delta_{t_{avg}} [ps]", "events", plotDir+"delta_tAvg", fit_delta_tAvg_energyAndPhaseCorr[Vov][vth1][vth2]);

          drawH1_deltaT(h1_delta_tAvg_energyCorr[Vov][vth1][vth2], h1_delta_tAvg_energyAndPhaseCorr[Vov][vth1][vth2], "#Delta_{t_{avg}} [ps]", "events", plotDir+"delta_tAvg", fit_delta_tAvg_energyCorr[Vov][vth1][vth2], fit_delta_tAvg_energyAndPhaseCorr[Vov][vth1][vth2]);

          if( g_tAvgRes_energyAndPhaseCorr_vs_vth1[Vov][vth2][iBar] == NULL )
          {
            g_tAvgRes_energyAndPhaseCorr_vs_vth1[Vov][vth2][iBar] = new TGraphErrors();
            g_tAvgRes_energyAndPhaseCorr_vs_vth1[Vov][vth2][iBar] -> SetName(Form("g_tAvgRes_energyAndPhaseCorr_vs_vth1_bar%02d_Vov%.1f_vth2_%02d",iBar,Vov,vth2));
            g_tAvg_actualRes_vs_vth1[Vov][vth2][iBar] = new TGraphErrors();
            g_tAvg_actualRes_vs_vth1[Vov][vth2][iBar] -> SetName(Form("g_tAvg_actualRes_vs_vth1_bar%02d_Vov%.1f_vth2_%02d",iBar,Vov,vth2));
          }
          float delta_tAvg_energyAndPhaseCorr_sigma = fit_delta_tAvg_energyAndPhaseCorr[Vov][vth1][vth2]->GetParameter(2);
          g_tAvgRes_energyAndPhaseCorr_vs_vth1[Vov][vth2][iBar] -> SetPoint(g_tAvgRes_energyAndPhaseCorr_vs_vth1[Vov][vth2][iBar]->GetN(),vth1,delta_tAvg_energyAndPhaseCorr_sigma);
          g_tAvg_actualRes_vs_vth1[Vov][vth2][iBar] -> SetPoint(g_tAvg_actualRes_vs_vth1[Vov][vth2][iBar]->GetN(),vth1,delta_tAvg_energyAndPhaseCorr_sigma/1.414);

          for(int iArray = 0; iArray < 2; ++iArray)
          {
            drawH1_deltaT(h1_deltaT_energyRatioCorr[Vov][vth1][vth2][iBar+num_bars*iArray],h1_deltaT_energyRatioPhaseCorr[Vov][vth1][vth2][iBar+num_bars*iArray], "#Deltat [ps]", "events", plotDir+"deltaT",fit_deltaT_energyRatioCorr[Vov][vth1][vth2][iBar+num_bars*iArray],fit_deltaT_energyRatioPhaseCorr[Vov][vth1][vth2][iBar+num_bars*iArray]);

            if( fit_deltaT_energyRatioPhaseCorr[Vov][vth1][vth2][iBar+num_bars*iArray] == NULL ) continue;

            if( !g_tRes_energyRatioPhaseCorr[Vov][vth1][vth2][iArray] )
            {
              g_tRes_energyRatioPhaseCorr[Vov][vth1][vth2][iArray] = new TGraphErrors();
              g_tRes_energyRatioPhaseCorr[Vov][vth1][vth2][iArray] -> SetName(Form("g_tRes_energyRatioPhaseCorr_array%d_Vov%.1f_vth1_%02d_vth2_%02d",iArray,Vov,vth1,vth2));
            }
            g_tRes_energyRatioPhaseCorr[Vov][vth1][vth2][iArray] -> SetPoint(g_tRes_energyRatioPhaseCorr[Vov][vth1][vth2][iArray]->GetN(),iBar,fit_deltaT_energyRatioPhaseCorr[Vov][vth1][vth2][iBar+num_bars*iArray]->GetParameter(2));
            if( g_tRes_energyRatioPhaseCorr_vs_vth1[Vov][vth2][iBar][iArray] == NULL )
            {
              g_tRes_energyRatioPhaseCorr_vs_vth1[Vov][vth2][iBar][iArray] = new TGraphErrors();
              g_tRes_energyRatioPhaseCorr_vs_vth1[Vov][vth2][iBar][iArray] -> SetName(Form("g_tRes_energyRatioPhaseCorr_vs_vth1_bar%02d_array%d_Vov%.1f_vth2_%02d",iBar,iArray,Vov,vth2));
              g_deltaT_actualRes_vs_vth1[Vov][vth2][iBar][iArray] = new TGraphErrors();
              g_deltaT_actualRes_vs_vth1[Vov][vth2][iBar][iArray] -> SetName(Form("g_tRes_energyRatioPhaseCorr_vs_vth1_bar%02d_array%d_Vov%.1f_vth2_%02d",iBar,iArray,Vov,vth2));
            }
            float tDiff_sigma = fit_deltaT_energyRatioPhaseCorr[Vov][vth1][vth2][iBar+num_bars*iArray]->GetParameter(2);
            float delta_tAvg_energyCorr_sigma = fit_delta_tAvg_energyCorr[Vov][vth1][vth2]->GetParameter(2);
            g_tRes_energyRatioPhaseCorr_vs_vth1[Vov][vth2][iBar][iArray] -> SetPoint(g_tRes_energyRatioPhaseCorr_vs_vth1[Vov][vth2][iBar][iArray]->GetN(),vth1,tDiff_sigma);
            g_deltaT_actualRes_vs_vth1[Vov][vth2][iBar][iArray] -> SetPoint(g_deltaT_actualRes_vs_vth1[Vov][vth2][iBar][iArray]->GetN(),vth1,tDiff_sigma/2.);
          }
        }
      }
    }
  }

  // time resolution vs. bar
  for(auto mapIt : h1_deltaT_energyRatioPhaseCorr )
    {
    float Vov = mapIt.first;
    for(auto mapIt2 : mapIt.second)
    {
      int vth1 = mapIt2.first;
      if( vth1 != vth1Ref ) continue;

      for(auto mapIt3 : mapIt2.second)
      {
        int vth2 = mapIt3.first;

        drawG_arrays(g_tRes_energyRatioPhaseCorr[Vov][vth1][vth2][0],g_tRes_energyRatioPhaseCorr[Vov][vth1][vth2][1], "iBar", "#sigma_{t} [ps]", 0., 300., plotDir+"timeResolution",false,Form("g_tRes_energyRatioPhaseCorr_Vov%.1f_vth1_%02d_vth2_%02d",Vov,vth1,vth2));
      }
    }
  }

  // time resolution vs. threshold
  std::map<int, std::map<int, std::map<int, std::vector<TGraphErrors*> > > > g_tRes_energyRatioPhaseCorr_vs_vth1_vec;
  std::map<int, std::map<int, std::map<int, std::vector<std::string> > > > g_tRes_energyRatioPhaseCorr_vs_vth1_labels;
  std::map<int, std::map<int, std::vector<TGraphErrors*> > > g_tAvgRes_energyAndPhaseCorr_vs_vth1_vec;
  std::map<int, std::map<int, std::vector<std::string> > > g_tAvgRes_energyAndPhaseCorr_vs_vth1_labels;
  std::map<int, std::map<int, std::vector<TGraphErrors*> > > g_tResSummary_energyAndPhaseCorr_vs_vth1_vec;
  std::map<int, std::map<int, std::vector<std::string> > > g_tResSummary_energyAndPhaseCorr_vs_vth1_labels;
  for(auto mapIt : g_tRes_energyRatioPhaseCorr_vs_vth1 )
  {
    float Vov = mapIt.first;
    for(auto mapIt2 : mapIt.second)
    {
      int vth2 = mapIt2.first;

      for(auto mapIt3 : mapIt2.second)
      {
        int iBar = mapIt3.first;

        g_tAvgRes_energyAndPhaseCorr_vs_vth1_vec[vth2][iBar].push_back(g_tAvgRes_energyAndPhaseCorr_vs_vth1[Vov][vth2][iBar]);
        g_tAvgRes_energyAndPhaseCorr_vs_vth1_labels[vth2][iBar].push_back(Form("V_{OV} = %.1f V",Vov));

        for(auto mapIt4 : mapIt3.second)
        {
          int iArray = mapIt4.first;
          g_tRes_energyRatioPhaseCorr_vs_vth1_vec[vth2][iBar][iArray].push_back(g_tRes_energyRatioPhaseCorr_vs_vth1[Vov][vth2][iBar][iArray]);
          g_tRes_energyRatioPhaseCorr_vs_vth1_labels[vth2][iBar][iArray].push_back(Form("V_{OV} = %.1f V",Vov));

          g_tResSummary_energyAndPhaseCorr_vs_vth1_vec[vth2][iBar].push_back(g_deltaT_actualRes_vs_vth1[Vov][vth2][iBar][iArray]);
          g_tResSummary_energyAndPhaseCorr_vs_vth1_labels[vth2][iBar].push_back(Form("deltaT; array %i; V_{OV} = %.1f V",iArray,Vov));
        }

        g_tResSummary_energyAndPhaseCorr_vs_vth1_vec[vth2][iBar].push_back(g_tAvg_actualRes_vs_vth1[Vov][vth2][iBar]);
        g_tResSummary_energyAndPhaseCorr_vs_vth1_labels[vth2][iBar].push_back(Form("delta_tAvg; array 1; V_{OV} = %.1f V",Vov));
      }
    }
  }
  for(auto mapIt : g_tRes_energyRatioPhaseCorr_vs_vth1_vec)
  {
    int vth2 = mapIt.first;

    for(auto mapIt2 : mapIt.second)
    {
      int iBar = mapIt2.first;

      drawG_vector(g_tAvgRes_energyAndPhaseCorr_vs_vth1_vec[vth2][iBar], "vth_{1} [DAC]", "#sigma_{#Delta_{t_{avg}}} [ps]", 0., 300., plotDir+"timeResolution",false,Form("g_tAvgRes_energyAndPhaseCorr_vs_vth1_vth2_%02d_bar%02d",vth2,iBar),&g_tAvgRes_energyAndPhaseCorr_vs_vth1_labels[vth2][iBar]);

      drawG_vector(g_tResSummary_energyAndPhaseCorr_vs_vth1_vec[vth2][iBar], "vth_{1} [DAC]", "Time resolution [ps]", 0., 300., plotDir+"timeResolution",false,Form("g_tResSummary_energyAndPhaseCorr_vs_vth1_vec%02d_bar%02d",vth2,iBar),&g_tResSummary_energyAndPhaseCorr_vs_vth1_labels[vth2][iBar]);

      for(auto mapIt3 : mapIt2.second)
      {
        int iArray = mapIt3.first;
        drawG_vector(g_tRes_energyRatioPhaseCorr_vs_vth1_vec[vth2][iBar][iArray], "vth_{1} [DAC]", "#sigma_{t_{diff}} [ps]", 0., 300., plotDir+"timeResolution",false,Form("g_tRes_energyRatioPhaseCorr_vs_vth1_vth2_%02d_bar%02d_array%d",vth2,iBar,iArray),&g_tRes_energyRatioPhaseCorr_vs_vth1_labels[vth2][iBar][iArray]);
      }
    }
  }

  // write output file
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}
