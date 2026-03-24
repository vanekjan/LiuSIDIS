#ifndef _TENSOR_SIDIS_ND3_H_
#define _TENSOR_SIDIS_ND3_H_

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"
#include "TChain.h"
#include "TVector3.h"

#include "Lsidis.h"

using namespace std;

//set numbers of protons and neutrons of target nucleus

//from SoLID sidis 2020
//const double Np = 0.844 / 3.0 * 10.0 + 0.69 / 2.0 + 0.47 / 2.0;

//const double Nn = 0.844 / 3.0 * 7.0 + 0.69 / 2.0 + 0.47 / 2.0; //for NH3
//const double Nn = 0.844 / 3.0 * 10.0 + 0.69 / 2.0 + 0.47 / 2.0; //for ND3

//---------------------------

//from base SoLID folder
//const double Np = 0.334*10.0+0.593*2.0;

//const double Nn = 0.334*7.0+0.593*2.0; //for NH3
//const double Nn = 0.334*10.0+0.593*2.0;   //for ND3

//---------------------------
//pure deuteron, deal with dilution later

const double Np = 1;

const double Nn = 1;

//---------------------------

// Acceptance 
TFile * file_e = new TFile("Acceptance/acceptance_solid_SIDIS_NH3_electron_201710_1e7_output_final.root", "r");
TFile * file_pip = new TFile("Acceptance/acceptance_solid_SIDIS_NH3_pip_201710_1e7_output_final.root", "r");
TFile * file_pim = new TFile("Acceptance/acceptance_solid_SIDIS_NH3_pim_201710_1e7_output_final.root", "r");
TFile * file_kp = new TFile("Acceptance/acceptance_solid_SIDIS_NH3_kp_201710_1e7_output_final.root", "r");
TFile * file_km = new TFile("Acceptance/acceptance_solid_SIDIS_NH3_km_201710_1e7_output_final.root", "r");
TH3F * acc_FA_e = (TH3F *) file_e->Get("acceptance_ThetaPhiP_forwardangle");
TH3F * acc_LA_e = (TH3F *) file_e->Get("acceptance_ThetaPhiP_largeangle");
TH3F * acc_FA_pip = (TH3F *) file_pip->Get("acceptance_ThetaPhiP_forwardangle");
TH3F * acc_LA_pip = (TH3F *) file_pip->Get("acceptance_ThetaPhiP_largeangle");
TH3F * acc_FA_pim = (TH3F *) file_pim->Get("acceptance_ThetaPhiP_forwardangle");
TH3F * acc_LA_pim = (TH3F *) file_pim->Get("acceptance_ThetaPhiP_largeangle");
TH3F * acc_FA_kp = (TH3F *) file_kp->Get("acceptance_ThetaPhiP_forwardangle");
TH3F * acc_LA_kp = (TH3F *) file_kp->Get("acceptance_ThetaPhiP_largeangle");
TH3F * acc_FA_km = (TH3F *) file_km->Get("acceptance_ThetaPhiP_forwardangle");
TH3F * acc_LA_km = (TH3F *) file_km->Get("acceptance_ThetaPhiP_largeangle");

double Rfactor0 = 1.0e5;

double GetAcceptance_e(const TLorentzVector p, const char * detector = "all"){//Get electron acceptance
  double theta = p.Theta() / M_PI * 180.0;
  double phi = p.Phi() / M_PI * 180.0;
  if (theta > 50.0) return 0;
  double mom = p.P();
  double acc = 0;
  if (strcmp(detector, "FA") == 0 || strcmp(detector, "all") == 0)
    acc += acc_FA_e->GetBinContent(acc_FA_e->GetXaxis()->FindBin(theta), acc_FA_e->GetYaxis()->FindBin(phi), acc_FA_e->GetZaxis()->FindBin(mom));
  if (mom > 3.5 && (strcmp(detector, "LA") == 0 || strcmp(detector, "all") == 0))
    acc += acc_LA_e->GetBinContent(acc_LA_e->GetXaxis()->FindBin(theta), acc_LA_e->GetYaxis()->FindBin(phi), acc_LA_e->GetZaxis()->FindBin(mom));
  return acc;
}

double GetAcceptance_pip(const TLorentzVector p){//Get pi+ acceptance
  double theta = p.Theta() / M_PI * 180.0;
  double phi = p.Phi() / M_PI * 180.0;
  if (theta > 45.0) return 0;
  double mom = p.P();
  double acc = 0;
  acc += acc_FA_pip->GetBinContent(acc_FA_pip->GetXaxis()->FindBin(theta), acc_FA_pip->GetYaxis()->FindBin(phi), acc_FA_pip->GetZaxis()->FindBin(mom));
  return acc;
}

double GetAcceptance_pim(const TLorentzVector p){//Get pi+ acceptance
  double theta = p.Theta() / M_PI * 180.0;
  double phi = p.Phi() / M_PI * 180.0;
  if (theta > 45.0) return 0;
  double mom = p.P();
  double acc = 0;
  acc += acc_FA_pim->GetBinContent(acc_FA_pim->GetXaxis()->FindBin(theta), acc_FA_pim->GetYaxis()->FindBin(phi), acc_FA_pim->GetZaxis()->FindBin(mom));
  return acc;
}

double PKmax = 7.5;
double GetAcceptance_Kp(const TLorentzVector p){//Get K+ acceptance
  double theta = p.Theta() / M_PI * 180.0;
  double phi = p.Phi() / M_PI * 180.0;
  if (theta > 45.0) return 0;
  double mom = p.P();
  if (p.P() > PKmax) return 0;
  double acc = 0;
  acc += acc_FA_kp->GetBinContent(acc_FA_kp->GetXaxis()->FindBin(theta), acc_FA_kp->GetYaxis()->FindBin(phi), acc_FA_kp->GetZaxis()->FindBin(mom));
  return acc;
}

double GetAcceptance_Km(const TLorentzVector p){//Get K+ acceptance
  double theta = p.Theta() / M_PI * 180.0;
  double phi = p.Phi() / M_PI * 180.0;
  if (theta > 45.0) return 0;
  double mom = p.P();
  if (p.P() > PKmax) return 0;
  double acc = 0;
  acc += acc_FA_km->GetBinContent(acc_FA_km->GetXaxis()->FindBin(theta), acc_FA_km->GetYaxis()->FindBin(phi), acc_FA_km->GetZaxis()->FindBin(mom));
  return acc;
}

double GetAcceptance_hadron(const TLorentzVector p, const char * hadron = "pi+"){//
  if (strcmp(hadron, "pi+") == 0) return GetAcceptance_pip(p);
  else if (strcmp(hadron, "pi-") == 0) return GetAcceptance_pim(p);
  else if (strcmp(hadron, "K+") == 0) return GetAcceptance_Kp(p);
  else if (strcmp(hadron, "K-") == 0) return GetAcceptance_Km(p);
  else return 0;
}

int GetTotalRate(const double Ebeam, const char * hadron){//Estimate the total rate
  Lsidis sidis;
  TLorentzVector l(0, 0, Ebeam, Ebeam);
  TLorentzVector P(0, 0, 0, 0.938272);
  sidis.SetNucleus(Np,Nn);
  sidis.SetHadron(hadron);
  if (strcmp(hadron, "pi+") == 0 || strcmp(hadron, "pi-") == 0) sidis.ChangeTMDpars(0.604, 0.114);
  if (strcmp(hadron, "K+") == 0 || strcmp(hadron, "K-") == 0) sidis.ChangeTMDpars(0.604, 0.131);
  sidis.SetInitialState(l, P);
  sidis.SetPDFset("CJ15lo");
  sidis.SetFFset("DSSFFlo");
  double lumi = 1.0e+9 * pow(0.197327, 2);
  double Xmin[6] = {0.0, 1.0, 0.3, 0.0, -M_PI, -M_PI};
  double Xmax[6] = {0.7, 8.0, 0.7, 1.6, M_PI, M_PI};
  sidis.SetRange(Xmin, Xmax);
  double sum = 0.0;
  Long64_t Nsim = 100000000;
  double weight = 0.0;
  TLorentzVector lp(0, 0, 0, 0);
  TLorentzVector Ph(0, 0, 0, 0);
  for (Long64_t i = 0; i < Nsim; i++){
    //if (i%(Nsim/20) == 0) std::cout << i << std::endl;
    weight = sidis.GenerateEvent(0, 1);
    if (weight > 0){
      if (sidis.GetVariable("W") < 2.3) continue;
      if (sidis.GetVariable("Wp") < 1.6) continue;
      sidis.CalculateRfactor();
      if (sidis.GetVariable("Rfactor") > Rfactor0) continue;
      lp = sidis.GetLorentzVector("lp");
      Ph = sidis.GetLorentzVector("Ph");
      sum += weight * GetAcceptance_e(lp) * GetAcceptance_hadron(Ph, hadron);
    }
  }
  printf("\n");
  printf("Total rate: %.4E  (%.1f GeV %s)\n\n", sum * lumi / Nsim, Ebeam, hadron);
  std::cout << "Total rate: " << sum * lumi / Nsim << std::endl;
  return 0;
}

int MakeKinematicCoveragePlots(const int nEvents, const double Ebeam, const char * hadron, const char * savefile)
{
//cout<<"Test 1"<<endl;
  float lumi = 1.0e+9 * pow(0.197327, 2);
  float scale_lumi = lumi/nEvents;

  //--------------------------------------

  Lsidis sidis;

  TLorentzVector l(0, 0, Ebeam, Ebeam);
  //TLorentzVector P(0, 0, 0, 0.938272); //proton
  //TLorentzVector P(0, 0, 0, 0.939565); //neutron
  TLorentzVector P(0, 0, 0, 1.875612945); //deuteron

  sidis.SetNucleus(Np,Nn);
  sidis.SetHadron(hadron);
  if (strcmp(hadron, "pi+") == 0 || strcmp(hadron, "pi-") == 0) sidis.ChangeTMDpars(0.604, 0.114);
  if (strcmp(hadron, "K+") == 0 || strcmp(hadron, "K-") == 0) sidis.ChangeTMDpars(0.604, 0.131);
  sidis.SetInitialState(l, P); //Set incoming lepton and nucleon 4-momenta
  sidis.SetPDFset("CJ15lo");
  sidis.SetFFset("DSSFFlo");

  //set variable ranges for the simulation
  //two possibilities based on selection in GenerateEvent()
  //method == 0 - generate in x, y, z, Pt, phih, phiS
  //method == 1 - generate in x, Q2, z, Pt, phih, phiS - this one is used now
  double Xmin[6] = {0.05, 1.0, 0.3, 0.0, -M_PI, -M_PI}; //my values
  double Xmax[6] = {0.4, 3.0, 0.7, 2.0, M_PI, M_PI}; //my values
  //double Xmin[6] = {0.0, 1.0, 0.3, 0.0, -M_PI, -M_PI}; //dafault from SoLID
  //double Xmax[6] = {0.7, 10.0, 0.7, 2.0, M_PI, M_PI}; //default from SoLID
  sidis.SetRange(Xmin, Xmax);

//cout<<"Test 2"<<endl;
  TFile * fs = new TFile(savefile, "RECREATE");
  gStyle->SetOptStat(0);
  //(x, Q2)
  TH2D * xQ2_FA = new TH2D("xQ2_FA", "", 700, 0.0, 0.7, 900, 0.0, 9.0);
  xQ2_FA->GetXaxis()->SetTitle("x");
  xQ2_FA->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");

  TH2D * xQ2_new_bins = new TH2D("xQ2_new_bins", "", 50, 0.11, 0.26, 50, 1.4, 2.8);
  xQ2_new_bins->GetXaxis()->SetTitle("x");
  xQ2_new_bins->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");

  //TH2D * xQ2_LA = new TH2D("xQ2_LA", "", 700, 0.0, 0.7, 900, 0.0, 9.0);
  //xQ2_LA->GetXaxis()->SetTitle("x");
  //xQ2_LA->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");

  //(x, W)
  TH2D * xW_FA = new TH2D("xW_FA", "", 700, 0.0, 0.7, 500, 2.0, 4.5);
  xW_FA->GetXaxis()->SetTitle("x");
  xW_FA->GetYaxis()->SetTitle("W / GeV");

  //TH2D * xW_LA = new TH2D("xW_LA", "", 700, 0.0, 0.7, 500, 2.0, 4.5);
  //xW_LA->GetXaxis()->SetTitle("x");
  //xW_LA->GetYaxis()->SetTitle("W / GeV");

  //(x, Wp)
  TH2D * xWp_FA = new TH2D("xWp_FA", "", 700, 0.0, 0.7, 500, 1.5, 4.0);
  xWp_FA->GetXaxis()->SetTitle("x");
  xWp_FA->GetYaxis()->SetTitle("W' / GeV");

  //TH2D * xWp_LA = new TH2D("xWp_LA", "", 700, 0.0, 0.7, 500, 1.5, 4.0);
  //xWp_LA->GetXaxis()->SetTitle("x");
  //xWp_LA->GetYaxis()->SetTitle("W' / GeV");

  //(x, z)
  TH2D * xz_FA = new TH2D("xz_FA", "", 700, 0.0, 0.7, 600, 0.2, 0.8);
  xz_FA->GetXaxis()->SetTitle("x");
  xz_FA->GetYaxis()->SetTitle("z");

  //TH2D * xz_LA = new TH2D("xz_LA", "", 700, 0.0, 0.7, 600, 0.2, 0.8);
  //xz_LA->GetXaxis()->SetTitle("x");
  //xz_LA->GetYaxis()->SetTitle("z");

  //(x, Pt)
  TH2D * xPt_FA = new TH2D("xPt_FA", "", 700, 0.0, 0.7, 800, 0.0, 2.0);
  xPt_FA->GetXaxis()->SetTitle("x");
  xPt_FA->GetYaxis()->SetTitle("P_{T} / GeV");

  //TH2D * xPt_LA = new TH2D("xPt_LA", "", 700, 0.0, 0.7, 800, 0.0, 2.0);
  //xPt_LA->GetXaxis()->SetTitle("x");
  //xPt_LA->GetYaxis()->SetTitle("P_{T} / GeV");

  //(z, Pt)
  TH2D * zPt_FA = new TH2D("zPt_FA", "", 600, 0.2, 0.8, 800, 0.0, 2.0);
  zPt_FA->GetXaxis()->SetTitle("z");
  zPt_FA->GetYaxis()->SetTitle("P_{T} / GeV");

  //bins are inverted compared to default histogram above
  TH2D * zPt_new_bins = new TH2D("zPt_new_bins", "", 50, 0., 0.9, 50, 0.30, 0.70);
  zPt_new_bins->GetXaxis()->SetTitle("P_{T} / GeV");
  zPt_new_bins->GetYaxis()->SetTitle("z");


  //TH2D * zPt_LA = new TH2D("zPt_LA", "", 600, 0.2, 0.8, 800, 0.0, 2.0);
  //zPt_LA->GetXaxis()->SetTitle("z");
  //zPt_LA->GetYaxis()->SetTitle("P_{T} / GeV");

  //(z, Q2)
  TH2D * zQ2_FA = new TH2D("zQ2_FA", "", 600, 0.2, 0.8, 900, 0.0, 9.0);
  zQ2_FA->GetXaxis()->SetTitle("z");
  zQ2_FA->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");

  //TH2D * zQ2_LA = new TH2D("zQ2_LA", "", 600, 0.2, 0.8, 900, 0.0, 9.0);
  //zQ2_LA->GetXaxis()->SetTitle("z");
  //zQ2_LA->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");

  //(z, W)
  TH2D * zW_FA = new TH2D("zW_FA", "", 600, 0.2, 0.8, 500, 2.0, 4.5);
  zW_FA->GetXaxis()->SetTitle("z");
  zW_FA->GetYaxis()->SetTitle("W / GeV");

  //TH2D * zW_LA = new TH2D("zW_LA", "", 600, 0.2, 0.8, 500, 2.0, 4.5);
  //zW_LA->GetXaxis()->SetTitle("z");
  //zW_LA->GetYaxis()->SetTitle("W / GeV");

  //(z, Wp)
  TH2D * zWp_FA = new TH2D("zWp_FA", "", 600, 0.2, 0.8, 500, 1.5, 4.0);
  zWp_FA->GetXaxis()->SetTitle("z");
  zWp_FA->GetYaxis()->SetTitle("W' / GeV");

  //TH2D * zWp_LA = new TH2D("zWp_LA", "", 600, 0.2, 0.8, 500, 1.5, 4.0);
  //zWp_LA->GetXaxis()->SetTitle("z");
  //zWp_LA->GetYaxis()->SetTitle("W' / GeV");

  //(Pt, Q2)
  TH2D * PtQ2_FA = new TH2D("PtQ2_FA", "", 800, 0.0, 2.0, 900, 0.0, 9.0);
  PtQ2_FA->GetXaxis()->SetTitle("P_{T} / GeV");
  PtQ2_FA->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");

  //TH2D * PtQ2_LA = new TH2D("PtQ2_LA", "", 800, 0.0, 2.0, 900, 0.0, 9.0);
  //PtQ2_LA->GetXaxis()->SetTitle("P_{T} / GeV");
  //PtQ2_LA->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");

  //(Pt, W)
  TH2D * PtW_FA = new TH2D("PtW_FA", "", 800, 0.0, 2.0, 500, 2.0, 4.5);
  PtW_FA->GetXaxis()->SetTitle("P_{T} / GeV");
  PtW_FA->GetYaxis()->SetTitle("W / GeV");

  //TH2D * PtW_LA = new TH2D("PtW_LA", "", 800, 0.0, 2.0, 500, 2.0, 4.5);
  //PtW_LA->GetXaxis()->SetTitle("P_{T} / GeV");
  //PtW_LA->GetYaxis()->SetTitle("W / GeV");

  //(Pt, Wp)
  TH2D * PtWp_FA = new TH2D("PtWp_FA", "", 800, 0.0, 2.0, 500, 1.5, 4.0);
  PtWp_FA->GetXaxis()->SetTitle("P_{T} / GeV");
  PtWp_FA->GetYaxis()->SetTitle("W' / GeV");

  //TH2D * PtWp_LA = new TH2D("PtWp_LA", "", 800, 0.0, 2.0, 500, 1.5, 4.0);
  //PtWp_LA->GetXaxis()->SetTitle("P_{T} / GeV");
  //PtWp_LA->GetYaxis()->SetTitle("W' / GeV");

  //(W, Q2)
  TH2D * WQ2_FA = new TH2D("WQ2_FA", "", 500, 2.0, 4.5, 900, 0.0, 9.0);
  WQ2_FA->GetXaxis()->SetTitle("W / GeV");
  WQ2_FA->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");

  //TH2D * WQ2_LA = new TH2D("WQ2_LA", "", 500, 2.0, 4.5, 900, 0.0, 9.0);
  //WQ2_LA->GetXaxis()->SetTitle("W / GeV");
  //WQ2_LA->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");

  //(Wp, Q2)
  TH2D * WpQ2_FA = new TH2D("WpQ2_FA", "", 500, 1.5, 4.0, 900, 0.0, 9.0);
  WpQ2_FA->GetXaxis()->SetTitle("W' / GeV");
  WpQ2_FA->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");

  //TH2D * WpQ2_LA = new TH2D("WpQ2_LA", "", 500, 1.5, 4.0, 900, 0.0, 9.0);
  //WpQ2_LA->GetXaxis()->SetTitle("W' / GeV");
  //WpQ2_LA->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");
  
  //---------------------------------------------------------------------------
  //basic kinematic plots
  
  TH1D *lepton_P = new TH1D("lepton_P", "lepton_P", 100, 0, 10);
  TH1D *lepton_P_cuts = new TH1D("lepton_P_cuts", "lepton_P_cuts", 100, 0, 10);

  TH1D *lepton_phi = new TH1D("lepton_phi", "lepton_phi", 360, -TMath::Pi(), TMath::Pi());
  TH1D *lepton_phi_cuts = new TH1D("lepton_phi_cuts", "lepton_phi_cuts", 360, -TMath::Pi(), TMath::Pi());

  TH1D *lepton_theta_1 = new TH1D("lepton_theta_1", "lepton_theta_1", 180, 0, TMath::Pi());
  TH1D *lepton_theta_1_cuts = new TH1D("lepton_theta_1_cuts", "lepton_theta_1_cuts", 180, 0, TMath::Pi());

  TH1D *lepton_theta_2 = new TH1D("lepton_theta_2", "lepton_theta_2", 180, 0, TMath::Pi());
  TH1D *lepton_theta_2_cuts = new TH1D("lepton_theta_2_cuts", "lepton_theta_2_cuts", 180, 0, TMath::Pi());

  TH1D *lepton_eta = new TH1D("lepton_eta", "lepton_eta", 100, -5, 5);
  TH1D *lepton_eta_cuts = new TH1D("lepton_eta_cuts", "lepton_eta_cuts", 100, -5, 5);
  
  //---------------------------------------
  
  //histograms for compariosn of simc and SoLID SIDIS model
  
  TH2D *lepton_P_vs_hadron_z = new TH2D("lepton_P_vs_hadron_z", "lepton_P_vs_hadron_z", 100, 0, 10, 50, 0.30, 0.70);
  TH2D *lepton_P_vs_hadron_Pt = new TH2D("lepton_P_vs_hadron_z", "lepton_P_vs_hadron_z", 100, 0, 10, 50, 0., 0.9);
  
  TH2D *lepton_theta_vs_hadron_z = new TH2D("lepton_theta_vs_hadron_z", "lepton_theta_vs_hadron_z", 180, 0, TMath::Pi(), 0.30, 0.70);
  TH2D *lepton_theta_vs_hadron_Pt = new TH2D("lepton_theta_vs_hadron_Pt", "lepton_theta_vs_hadron_Pt", 180, 0, TMath::Pi(), 0., 0.9);
  
  TH2D *lepton_P_vs_hadron_z = new TH2D("lepton_P_vs_hadron_z", "lepton_P_vs_hadron_z", 360, -TMath::Pi(), TMath::Pi(), 0.30, 0.70);
  TH2D *lepton_P_vs_hadron_Pt = new TH2D("lepton_P_vs_hadron_z", "lepton_P_vs_hadron_z", 1360, -TMath::Pi(), TMath::Pi(), 0., 0.9);
  
  //-------
  
  TH2D *lepton_P_vs_hadron_z_cuts = new TH2D("lepton_P_vs_hadron_z_cuts", "lepton_P_vs_hadron_z_cuts", 100, 0, 10, 50, 0.30, 0.70);
  TH2D *lepton_P_vs_hadron_Pt_cuts = new TH2D("lepton_P_vs_hadron_Pt_cuts", "lepton_P_vs_hadron_Pt_cuts", 100, 0, 10, 50, 0., 0.9);
  
  TH2D *lepton_theta_vs_hadron_z_cuts = new TH2D("lepton_theta_vs_hadron_z_cuts", "lepton_theta_vs_hadron_z_cuts", 180, 0, TMath::Pi(), 0.30, 0.70);
  TH2D *lepton_theta_vs_hadron_Pt_cuts = new TH2D("lepton_theta_vs_hadron_Pt_cuts", "lepton_theta_vs_hadron_Pt_cuts", 180, 0, TMath::Pi(), 0., 0.9);
  
  TH2D *lepton_P_vs_hadron_z_cuts = new TH2D("lepton_P_vs_hadron_z_cuts", "lepton_P_vs_hadron_z_cuts", 360, -TMath::Pi(), TMath::Pi(), 0.30, 0.70);
  TH2D *lepton_P_vs_hadron_Pt_cuts = new TH2D("lepton_P_vs_hadron_Pt_cuts", "lepton_P_vs_hadron_Pt_cuts", 1360, -TMath::Pi(), TMath::Pi(), 0., 0.9);

  //---------------------------------------

  TH1D *hadron_P = new TH1D("hadron_P", "hadron_P", 100, 0, 10);
  TH1D *hadron_P_cuts = new TH1D("hadron_P_cuts", "hadron_P_cuts", 100, 0, 10);

  TH1D *hadron_phi = new TH1D("hadron_phi", "hadron_phi", 360, -TMath::Pi(), TMath::Pi());
  TH1D *hadron_phi_cuts = new TH1D("hadron_phi_cuts", "hadron_phi_cuts", 360, -TMath::Pi(), TMath::Pi());

  TH1D *hadron_theta_1 = new TH1D("hadron_theta_1", "hadron_theta_1", 360, -TMath::Pi(), TMath::Pi());
  TH1D *hadron_theta_1_cuts = new TH1D("hadron_theta_1_cuts", "hadron_theta_1_cuts", 360, -TMath::Pi(), TMath::Pi());

  TH1D *hadron_theta_2 = new TH1D("hadron_theta_2", "hadron_theta_2", 360, -TMath::Pi(), TMath::Pi());
  TH1D *hadron_theta_2_cuts = new TH1D("hadron_theta_2_cuts", "hadron_theta_2_cuts", 360, -TMath::Pi(), TMath::Pi());

  TH1D *hadron_eta = new TH1D("hadron_eta", "hadron_eta", 100, -5, 5);
  TH1D *hadron_eta_cuts = new TH1D("hadron_eta_cuts", "hadron_eta_cuts", 100, -5, 5);
  
  //---------------------------------------------------------------------------

  //F_U(LL),T histograms

  TH1D *F_UUT_hist = new TH1D("F_UUT_hist", "F_UUT_hist", 11, -0.1, 1);
  TH1D *F_UUT_weight_hist = new TH1D("F_UUT_weight_hist", "F_UUT_weight_hist", 11, -0.1, 1);

  TH1D *F_ULLT_hist = new TH1D("F_ULLT_hist", "F_ULLT_hist", 11, -0.1, 1);

  TH1D *F_UUT_hist_new_bins = new TH1D("F_UUT_hist_new_bins", "F_UUT_hist_new_bins", 10, 0, 0.7);
  TH1D *F_UUT_weight_hist_new_bins = new TH1D("F_UUT_weight_hist_new_bins", "F_UUT_weight_hist_new_bins", 10, 0., 0.7);
  
  TH1D *F_ULLT_hist_new_bins = new TH1D("F_ULLT_hist_new_bins", "F_ULLT_hist_new_bins", 10, 0, 0.7);

  //----

  TH1D *P_T_hist = new TH1D("P_T_hist", "P_T_hist", 11, -0.1, 1);
  TH1D *P_T_no_cut_hist = new TH1D("P_T_no_cut_hist", "P_T_no_cut_hist", 11, -0.1, 1);
  TH1D *P_T_no_kine_cut_hist = new TH1D("P_T_no_kine_cut_hist", "P_T_no_kine_cut_hist", 11, -0.1, 1);

  TH1D *P_T_weight_hist = new TH1D("P_T_weight_hist", "P_T_weight_hist", 11, -0.1, 1);
  TH1D *P_T_weight_no_cut_hist = new TH1D("P_T_weight_no_cut_hist", "P_T_weight_no_cut_hist", 11, -0.1, 1);
  TH1D *P_T_weight_no_kine_cut_hist = new TH1D("P_T_weight_no_kine_cut_hist", "P_T_weight_no_kine_cut_hist", 11, -0.1, 1);


  TH1D *P_T_hist_new_bins = new TH1D("P_T_hist_new_bins", "P_T_hist_new_bins", 10, 0.0, 0.7);
  TH1D *P_T_no_cut_hist_new_bins = new TH1D("P_T_no_cut_hist_new_bins", "P_T_no_cut_hist_new_bins", 10, 0.0, 0.7);
  TH1D *P_T_no_kine_cut_hist_new_bins = new TH1D("P_T_no_kine_cut_hist_new_bins", "P_T_no_kine_cut_hist_new_bins", 10, 0.0, 0.7);

  TH1D *P_T_weight_hist_new_bins = new TH1D("P_T_weight_hist_new_bins", "P_T_weight_hist_new_bins", 10, 0.0, 0.7);
  TH1D *P_T_weight_no_cut_hist_new_bins = new TH1D("P_T_weight_no_cut_hist_new_bins", "P_T_weight_no_cut_hist_new_bins", 10, 0.0, 0.7);
  TH1D *P_T_weight_no_kine_cut_hist_new_bins = new TH1D("P_T_weight_no_kine_cut_hist_new_bins", "P_T_weight_no_kine_cut_hist_new_bins", 10, 0.0, 0.7);

  //---------------------------------------------------------------------------
  
  TH1D *nEvents_hist = new TH1D("nEvents_hist", "nEvents_hist", 1, -0.5, 0.5);

  nEvents_hist->SetBinContent(1, nEvents);
  nEvents_hist->SetBinError(1, sqrt(nEvents));

  TH1D *nEvents_event_cuts_hist = new TH1D("nEvents_event_cuts_hist", "nEvents_event_cuts_hist", 1, -0.5, 0.5);

  TH1D *nEvents_all_cuts_hist = new TH1D("nEvents_all_cuts_hist", "nEvents_all_cuts_hist", 1, -0.5, 0.5);
  //--------------------------------------------------------------------------

  double x, Q2, z, z_my, Pt, Pt_my, W, Wp;
  double weight, acc_FA, acc_LA;
  TLorentzVector l_beam, lp, Ph;

      
  //-----------

  
  int nEvents_event_cuts = 0;
  int nEvents_all_cuts = 0;

  for (Long64_t i = 0; i < nEvents; i++)
  {
    //if (i % 1000000 == 0) std::cout << i*100/nEvents << " %" << std::endl; //original
    if (i % nEvents/10 == 0) std::cout << i*100/nEvents << " %" << std::endl;
    weight = sidis.GenerateEvent(0, 1);
    if (weight > 0)
    {
      //weight = 1.; //for testing
      
      P_T_no_cut_hist->Fill(Pt);
      P_T_no_cut_hist_new_bins->Fill(Pt);

      //these volumes neeed to be changed to the generator level ranges
      P_T_weight_no_cut_hist->Fill(Pt, weight);
      P_T_weight_no_cut_hist_new_bins->Fill(Pt, weight);

      z = sidis.GetVariable("z");
      if (z < 0.3 || z > 0.7) continue;

      Q2 = sidis.GetVariable("Q2");
      //if (Q2 < 1.0 || Q2 > 2.5) continue;

      W = sidis.GetVariable("W");
      if (W < 2.3) continue;

      Wp = sidis.GetVariable("Wp");
      if (Wp < 1.6) continue;

      x = sidis.GetVariable("x");
      //if (x < 0.1 || x > 0.35) continue;

      Pt = sidis.GetVariable("Pt");

      nEvents_event_cuts++;

      //-----------

      P_T_no_kine_cut_hist->Fill(Pt);
      P_T_no_kine_cut_hist_new_bins->Fill(Pt);

      P_T_weight_no_kine_cut_hist->Fill(Pt, weight);
      P_T_weight_no_kine_cut_hist_new_bins->Fill(Pt, weight);

      
      l_beam = sidis.GetLorentzVector("l");//beam lepton 4-momentum
      lp = sidis.GetLorentzVector("lp");//scattered lepton 4-momentum
      Ph = sidis.GetLorentzVector("Ph");//outgoing hadron 4-momentum

      TLorentzVector q_4mom = l_beam - lp;
           
      //cuts
      TVector3 z_axis(0,0,1);
      float theta_l = lp.Angle(z_axis);
      float theta_h = Ph.Angle(z_axis);

      //my calculation of pT (from LOI - from event cartoon)
      Pt_my = Ph.P() * sin(Ph.Angle(q_4mom.Vect()));

      //cuts QA histograms
      //electron
      lepton_P->Fill(lp.P(), weight);
      lepton_phi->Fill(lp.Phi(), weight);
      lepton_eta->Fill(lp.Eta(), weight);

      lepton_theta_1->Fill(theta_l, weight);
      lepton_theta_2->Fill(lp.Theta(), weight);
      //lepton_theta_2->Fill(2*atan(exp(-lp.Eta())), weight);

      //hadron
      hadron_P->Fill(Ph.P(), weight);
      hadron_phi->Fill(Ph.Phi(), weight);
      hadron_eta->Fill(Ph.Eta(), weight);
      
      hadron_theta_1->Fill(theta_h, weight);
      hadron_theta_2->Fill(Ph.Theta(), weight);
      //hadron_theta_2->Fill(2*atan(exp(-Ph.Eta())), weight);
      
      
      //2D distributions
      
      lepton_P_vs_hadron_z->Fill(lp.P(), z);
      lepton_P_vs_hadron_Pt->Fill(lp.P(), Pt);
      
      lepton_theta_vs_hadron_z->Fill(lp.Theta(), z);
      lepton_theta_vs_hadron_Pt->Fill(lp.Theta(), Pt);
      
      lepton_phi_vs_hadron_z->Fill(lp.Phi(), z);
      lepton_phi_vs_hadron_Pt->Fill(lp.Phi(), Pt);

      //-----------------------------------------
      //Q2 and x cuts 
      //if (Q2 < 1.0 || Q2 > 2.5) continue;
      //if (x < 0.1 || x > 0.35) continue;

      //electron
      //cuts from LOI: theta = (10.3 - 12.4) deg, phi = (-2.28 - 2.28) deg, P = (4.0 - 5.4) GeV
      if( lp.P() < 4.0 || lp.P() > 5.4 ) continue;
      //cuts from above converted to rad
      if( lp.Phi() < -0.03979 || lp.Phi() > 0.03979 ) continue;
      if( theta_l < 0.178 || theta_l > 0.2164208 ) continue; //theta calculated from pseudorapidity

      //hadron
      //cuts from LOI: theta = (5.0 - 15.0) deg, phi = (167 - 197) deg (this is +-13 deg around pi, i.e. directly back from electron), P = (2.0 - 4.0) GeV
      if( Ph.P() < 2.0 || Ph.P() > 4.0 ) continue;
      //cuts from above converted to rad
      if( Ph.Phi() < 2.9147 && Ph.Phi() > -2.9147 ) continue; //this is around pi (or -pi), i.e. in "backward" region 
      //if( Ph.Phi() < 2.9147 ) continue; //this is around pi (or -pi), i.e. in "backward" region 
      if( theta_h < 2.*0.0872665 || theta_h > 4.*0.0872665 ) continue; //theta = (10.0 - 20.0) deg

      nEvents_all_cuts++;

      //-----------------------------------------
      //
      //cuts QA histograms
      //electron
      lepton_P_cuts->Fill(lp.P(), weight);
      lepton_phi_cuts->Fill(lp.Phi(), weight);
      lepton_eta_cuts->Fill(lp.Eta(), weight);
 
      lepton_theta_1_cuts->Fill(theta_l, weight);
      lepton_theta_2_cuts->Fill(lp.Theta(), weight);
      //lepton_theta_2_cuts->Fill(2*atan(exp(-lp.Eta())), weight);
 
      //hadron
      hadron_P_cuts->Fill(Ph.P(), weight);
      hadron_phi_cuts->Fill(Ph.Phi(), weight);
      hadron_eta_cuts->Fill(Ph.Eta(), weight);
 
      hadron_theta_1_cuts->Fill(theta_h, weight);
      hadron_theta_2_cuts->Fill(Ph.Theta(), weight);
      //hadron_theta_2_cuts->Fill(2*atan(exp(-Ph.Eta())), weight);
      
      
      //2D distributions
      
      lepton_P_vs_hadron_z_cuts->Fill(lp.P(), z);
      lepton_P_vs_hadron_Pt_cuts->Fill(lp.P(), Pt);
      
      lepton_theta_vs_hadron_z_cuts->Fill(lp.Theta(), z);
      lepton_theta_vs_hadron_Pt_cuts->Fill(lp.Theta(), Pt);
      
      lepton_phi_vs_hadron_z_cuts->Fill(lp.Phi(), z);
      lepton_phi_vs_hadron_Pt_cuts->Fill(lp.Phi(), Pt);


      //acc_FA = GetAcceptance_e(lp, "FA") * GetAcceptance_hadron(Ph, hadron);
      acc_FA = 1.; //turn off SoLID acceptance x efficiency weights

      //acc_LA = GetAcceptance_e(lp, "LA") * GetAcceptance_hadron(Ph, hadron);
      acc_LA = 1.; //turn off SoLID acceptance x efficiency weights

      if (acc_FA > 0)
      {
        //kinematic histograms
	      xQ2_FA->Fill(x, Q2, acc_FA*weight);
	      xQ2_new_bins->Fill(x, Q2, acc_FA*weight);
	      xW_FA->Fill(x, W, acc_FA*weight);
      	xz_FA->Fill(x, z, acc_FA*weight);
      	xPt_FA->Fill(x, Pt, acc_FA*weight);
      	xWp_FA->Fill(x, Wp, acc_FA*weight);
      	zPt_FA->Fill(z, Pt, acc_FA*weight);
      	zPt_new_bins->Fill(Pt, z, acc_FA*weight);//this histogram has inverted axes compared to the default one
      	zQ2_FA->Fill(z, Q2, acc_FA*weight);
      	zW_FA->Fill(z, W, acc_FA*weight);
      	zWp_FA->Fill(z, Wp, acc_FA*weight);
      	PtQ2_FA->Fill(Pt, Q2, acc_FA*weight);
      	PtW_FA->Fill(Pt, W, acc_FA*weight);
      	PtWp_FA->Fill(Pt, Wp, acc_FA*weight);
      	WQ2_FA->Fill(W, Q2, acc_FA*weight);
      	WpQ2_FA->Fill(Wp, Q2, acc_FA*weight);


        //F_U(LL),T histograms
        F_UUT_hist->Fill(Pt, sidis.FUUT()*acc_FA);
        F_UUT_weight_hist->Fill(Pt, sidis.FUUT()*acc_FA*weight);
        //F_UUT_weight_hist->Fill(Pt, acc_FA*weight); //weight already contains FUUT()
        
        F_ULLT_hist->Fill(Pt, sidis.FUUT()*acc_FA*0.1); //fill F_UUT scaled by 1/10 to estimate F_ULLT

        F_UUT_hist_new_bins->Fill(Pt, sidis.FUUT()*acc_FA);
        F_UUT_weight_hist_new_bins->Fill(Pt, sidis.FUUT()*acc_FA);

        F_ULLT_hist_new_bins->Fill(Pt, sidis.FUUT()*acc_FA*0.1); //fill F_UUT scaled by 1/10 to estimate F_ULLT


        P_T_hist->Fill(Pt, acc_FA);
        P_T_weight_hist->Fill(Pt, acc_FA*weight);
        //P_T_weight_hist->Fill(Pt, acc_FA*weight/sidis.FUUT()); //get rid of FUUT in this distribution - for scaling of F_UUT_weight_hist
        
        P_T_hist_new_bins->Fill(Pt, acc_FA);
        P_T_weight_hist_new_bins->Fill(Pt, acc_FA*weight);

          
      }



    }
  }

  nEvents_event_cuts_hist->SetBinContent(1, nEvents_event_cuts);
  nEvents_event_cuts_hist->SetBinError(1, sqrt(nEvents_event_cuts));

  nEvents_all_cuts_hist->SetBinContent(1, nEvents_all_cuts);
  nEvents_all_cuts_hist->SetBinError(1, sqrt(nEvents_all_cuts));

  fs->Write();
  return 0;
}

int MakeRateDistributionPlots(const double Ebeam, const char * hadron, const char * savefile){
  double lumi = 1.0e+9 * pow(0.197327, 2);
  Lsidis sidis;
  TLorentzVector l(0, 0, Ebeam, Ebeam);
  TLorentzVector P(0, 0, 0, 0.938272);
  sidis.SetNucleus(Np,Nn);
  sidis.SetHadron(hadron);
  if (strcmp(hadron, "pi+") == 0 || strcmp(hadron, "pi-") == 0) sidis.ChangeTMDpars(0.604, 0.114);
  if (strcmp(hadron, "K+") == 0 || strcmp(hadron, "K-") == 0) sidis.ChangeTMDpars(0.604, 0.131);
  sidis.SetInitialState(l, P);
  sidis.SetPDFset("CJ15lo");
  sidis.SetFFset("DSSFFlo");
  double Xmin[6] = {0.0, 1.0, 0.3, 0.0, -M_PI, -M_PI};
  double Xmax[6] = {0.7, 10.0, 0.7, 2.0, M_PI, M_PI};
  sidis.SetRange(Xmin, Xmax);
  TFile * fs = new TFile(savefile, "RECREATE");
  gStyle->SetOptStat(0);
  //(x, Q2)
  TH2D * xQ2 = new TH2D("xQ2", "", 700, 0.0, 0.7, 900, 0.0, 9.0);
  xQ2->GetXaxis()->SetTitle("x");
  xQ2->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");
  //(x, W)
  TH2D * xW = new TH2D("xW", "", 700, 0.0, 0.7, 500, 2.0, 4.5);
  xW->GetXaxis()->SetTitle("x");
  xW->GetYaxis()->SetTitle("W / GeV");
  //(x, Wp)
  TH2D * xWp = new TH2D("xWp", "", 700, 0.0, 0.7, 500, 1.5, 4.0);
  xWp->GetXaxis()->SetTitle("x");
  xWp->GetYaxis()->SetTitle("W' / GeV");
  //(x, z)
  TH2D * xz = new TH2D("xz", "", 700, 0.0, 0.7, 600, 0.2, 0.8);
  xz->GetXaxis()->SetTitle("x");
  xz->GetYaxis()->SetTitle("z");
  //(x, Pt)
  TH2D * xPt = new TH2D("xPt", "", 700, 0.0, 0.7, 800, 0.0, 2.0);
  xPt->GetXaxis()->SetTitle("x");
  xPt->GetYaxis()->SetTitle("P_{T} / GeV");
  //(z, Pt)
  TH2D * zPt = new TH2D("zPt", "", 600, 0.2, 0.8, 800, 0.0, 2.0);
  zPt->GetXaxis()->SetTitle("z");
  zPt->GetYaxis()->SetTitle("P_{T} / GeV");
  //(z, Q2)
  TH2D * zQ2 = new TH2D("zQ2", "", 600, 0.2, 0.8, 900, 0.0, 9.0);
  zQ2->GetXaxis()->SetTitle("z");
  zQ2->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");
  //(z, W)
  TH2D * zW = new TH2D("zW", "", 600, 0.2, 0.8, 500, 2.0, 4.5);
  zW->GetXaxis()->SetTitle("z");
  zW->GetYaxis()->SetTitle("W / GeV");
  //(z, Wp)
  TH2D * zWp = new TH2D("zWp", "", 600, 0.2, 0.8, 500, 1.5, 4.0);
  zWp->GetXaxis()->SetTitle("z");
  zWp->GetYaxis()->SetTitle("W' / GeV");
  //(Pt, Q2)
  TH2D * PtQ2 = new TH2D("PtQ2", "", 800, 0.0, 2.0, 900, 0.0, 9.0);
  PtQ2->GetXaxis()->SetTitle("P_{T} / GeV");
  PtQ2->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");
  //(Pt, W)
  TH2D * PtW = new TH2D("PtW", "", 800, 0.0, 2.0, 500, 2.0, 4.5);
  PtW->GetXaxis()->SetTitle("P_{T} / GeV");
  PtW->GetYaxis()->SetTitle("W / GeV");
  //(Pt, Wp)
  TH2D * PtWp = new TH2D("PtWp", "", 800, 0.0, 2.0, 500, 1.5, 4.0);
  PtWp->GetXaxis()->SetTitle("P_{T} / GeV");
  PtWp->GetYaxis()->SetTitle("W' / GeV");
  //(W, Q2)
  TH2D * WQ2 = new TH2D("WQ2", "", 500, 2.0, 4.5, 900, 0.0, 9.0);
  WQ2->GetXaxis()->SetTitle("W / GeV");
  WQ2->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");
  //(Wp, Q2)
  TH2D * WpQ2 = new TH2D("WpQ2", "", 500, 1.5, 4.0, 900, 0.0, 9.0);
  WpQ2->GetXaxis()->SetTitle("W' / GeV");
  WpQ2->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");
  double x, Q2, z, Pt, W, Wp;
  double weight, acc;
  TLorentzVector lp, Ph;
  Long64_t Nsim = 100000000;
  printf("Generating rate plots: (%.1f GeV %s)\n", Ebeam, hadron);
  for (Long64_t i = 0; i < Nsim; i++){
    //if (i % 10000000 == 0) std::cout << i << " %" << std::endl;
    weight = sidis.GenerateEvent(0, 1);
    if (weight > 0){
      z = sidis.GetVariable("z");
      if (z < 0.3 || z > 0.7) continue;
      Q2 = sidis.GetVariable("Q2");
      if (Q2 < 1.0) continue;
      W = sidis.GetVariable("W");
      if (W < 2.3) continue;
      Wp = sidis.GetVariable("Wp");
      if (Wp < 1.6) continue;
      x = sidis.GetVariable("x");
      Pt = sidis.GetVariable("Pt");
      lp = sidis.GetLorentzVector("lp");
      Ph = sidis.GetLorentzVector("Ph");
      acc = GetAcceptance_e(lp, "all") * GetAcceptance_hadron(Ph, hadron);
      if (acc > 0){
	sidis.CalculateRfactor();
	if (sidis.GetVariable("Rfactor") > Rfactor0) continue;
	xQ2->Fill(x, Q2, weight * acc);
	xW->Fill(x, W, weight * acc);
	xz->Fill(x, z, weight * acc);
	xPt->Fill(x, Pt, weight * acc);
	xWp->Fill(x, Wp, weight * acc);
	zPt->Fill(z, Pt, weight * acc);
	zQ2->Fill(z, Q2, weight * acc);
	zW->Fill(z, W, weight * acc);
	zWp->Fill(z, Wp, weight * acc);
	PtQ2->Fill(Pt, Q2, weight * acc);
	PtW->Fill(Pt, W, weight * acc);
	PtWp->Fill(Pt, Wp, weight * acc);
	WQ2->Fill(W, Q2, weight * acc);
	WpQ2->Fill(Wp, Q2, weight * acc);
      }
    }
  }
  xQ2->Scale(lumi/Nsim);
  xW->Scale(lumi/Nsim);
  xz->Scale(lumi/Nsim);
  xPt->Scale(lumi/Nsim);
  xWp->Scale(lumi/Nsim);
  zW->Scale(lumi/Nsim);
  zQ2->Scale(lumi/Nsim);
  zWp->Scale(lumi/Nsim);
  zPt->Scale(lumi/Nsim);
  PtQ2->Scale(lumi/Nsim);
  WQ2->Scale(lumi/Nsim);
  PtW->Scale(lumi/Nsim);
  PtWp->Scale(lumi/Nsim);
  WpQ2->Scale(lumi/Nsim);
  fs->Write();
  return 0;
}

int GenerateBinInfoFile(const char * filename, const double Ebeam, const char * hadron){//Bin the data and create the bin info file
  FILE * fp = fopen(filename, "w");
  fprintf(fp, "Q2l\t Q2u\t zl\t zu\t Ptl\t Ptu\t xl\t xu\n");
  Lsidis sidis;
  TLorentzVector l(0, 0, Ebeam, Ebeam);
  TLorentzVector P(0, 0, 0, 0.938272);
  sidis.SetNucleus(Np,Nn);
  sidis.SetHadron(hadron);
  if (strcmp(hadron, "pi+") == 0 || strcmp(hadron, "pi-") == 0) sidis.ChangeTMDpars(0.604, 0.114);
  if (strcmp(hadron, "K+") == 0 || strcmp(hadron, "K-") == 0) sidis.ChangeTMDpars(0.604, 0.131);
  sidis.SetInitialState(l, P);
  sidis.SetPDFset("CJ15lo");
  sidis.SetFFset("DSSFFlo");
  double lumi = 1.0e+9 * pow(0.197327, 2);
  double eff = 0.85;
  double time = 55.0 * 24.0 * 3600.0;
  if (Ebeam < 10.0) time = 27.5 * 24.0 * 3600.0;
  double Nsim = 1.0e6;
  double Xmin[6] = {0.0, 0.0, 0.0, 0.0, -M_PI, -M_PI}; 
  double Xmax[6] = {0.7, 0.0, 0.0, 0.0, M_PI, M_PI};;//x, Q2, z, Pt, phih, phiS
  double weight = 0;
  double acc = 0;
  int Nx = 0;
  TLorentzVector lp(0, 0, 0, 0);
  TLorentzVector Ph(0, 0, 0, 0);
  double Q2list[7] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0};
  double statlist[6] = {1.0e7, 6.4e6, 3.2e6, 1.6e6, 1.2e6, 1.0e6};
  double zlist[9] = {0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7};
  double Ptlist[7] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.6};
  int xi = 1;
  for (int Qi = 0; Qi < 6; Qi++){//Q2 loop
    Xmin[1] = Q2list[Qi];
    Xmax[1] = Q2list[Qi+1];
    for (int zi = 0; zi < 8; zi++){//z loop
      Xmin[2] = zlist[zi];
      Xmax[2] = zlist[zi+1];
      Xmin[3] = Ptlist[0];
      for (int kj = 1; kj < 7;){//Pt loop
	Xmax[3] = Ptlist[kj];
	sidis.SetRange(Xmin, Xmax);
	TH1D * hx = new TH1D("hx", "hx", 7000, 0.0, 0.7);
	printf("Q2:%.1f-%.1f  z:%.2f-%.2f  Pt:%.1f-%.1f\n",
	       Xmin[1], Xmax[1], Xmin[2], Xmax[2], Xmin[3], Xmax[3]);
	for (Long64_t i = 0; i < Nsim; i++){//generate events
	  weight = sidis.GenerateEvent(0, 1);
	  if (weight > 0){
            if (sidis.GetVariable("W") < 2.3) continue;
            if (sidis.GetVariable("Wp") < 1.6) continue;
	    sidis.CalculateRfactor();
	    if (sidis.GetVariable("Rfactor") > Rfactor0) continue;
	    lp = sidis.GetLorentzVector("lp");
	    Ph = sidis.GetLorentzVector("Ph");
	    acc = GetAcceptance_e(lp) * GetAcceptance_hadron(Ph, hadron);
	    if (acc > 0)
	      hx->Fill(sidis.GetVariable("x"), weight * acc);
	  }
	}
	hx->Scale(lumi * time * eff / Nsim);
	if ((hx->Integral(1, -1) < statlist[Qi] && kj < 6) || (hx->Integral(1, -1) < 0.2 * statlist[Qi] && kj == 6)){
	  hx->Delete();
	  kj++;
	  continue;
	}
	Nx = 0;
	xi = 1;
	for (int xj = 1; xj <= 7000; xj++){
	  if (hx->Integral(xi, xj) > statlist[Qi] || xj == 7000){
	    fprintf(fp, "%.1f\t %.1f\t %.2f\t %.2f\t %.1f\t %.1f\t %.4f\t %.4f\n",
		    Xmin[1], Xmax[1], Xmin[2], Xmax[2], Xmin[3], Xmax[3],
		    hx->GetBinLowEdge(xi), hx->GetBinLowEdge(xj+1));
	    Nx++;
	    xi = xj + 1;
	  }
	}
	std::cout << Nx << std::endl;
	hx->Delete();
	Xmin[3] = Ptlist[kj];
	kj++;
      }
    }
  }
  fclose(fp);
  return 0;
}

int AnalyzeEstatUT3(const char * readfile, const char * savefile, const double Ebeam, const char * had){//bin analysis including stat. errors
  double Hadron = 0;
  if (strcmp(had, "pi+") == 0) Hadron = 0;
  else if (strcmp(had, "pi-") == 0) Hadron = 1;
  double Nucleon = 1;
  TFile * fs = new TFile(savefile, "RECREATE");
  TTree * Ts = new TTree("data", "data");
  Ts->SetDirectory(fs);
  double Eb = Ebeam;
  double x, y, z, Q2, Pt, phih, phiS;
  double dx, dy, dz, dQ2, dPt, dphih, dphiS, dv;
  double Nacc, fp;
  double Estatraw[3], Estat[3];
  Ts->Branch("Nucleon", &Nucleon, "Nucleon/D");
  Ts->Branch("Hadron", &Hadron, "Hadron/D");
  Ts->Branch("Ebeam", &Eb, "Ebeam/D");
  Ts->Branch("x", &x, "x/D");
  Ts->Branch("y", &y, "y/D");
  Ts->Branch("z", &z, "z/D");
  Ts->Branch("Q2", &Q2, "Q2/D");
  Ts->Branch("Pt", &Pt, "Pt/D");
  Ts->Branch("dx", &dx, "dx/D");
  Ts->Branch("dy", &dy, "dy/D");
  Ts->Branch("dz", &dz, "dz/D");
  Ts->Branch("dQ2", &dQ2, "dQ2/D");
  Ts->Branch("dPt", &dPt, "dPt/D");
  Ts->Branch("dphih", &dphih, "dphih/D");
  Ts->Branch("dphiS", &dphiS, "dphiS/D");
  Ts->Branch("dv", &dv, "dv/D");
  Ts->Branch("Nacc", &Nacc, "Nacc/D");
  Ts->Branch("fp", &fp, "fp/D");
  Ts->Branch("E0statraw", &Estatraw[0], "E0statraw/D");
  Ts->Branch("E1statraw", &Estatraw[1], "E1statraw/D");
  Ts->Branch("E2statraw", &Estatraw[2], "E2statraw/D");
  Ts->Branch("E0stat", &Estat[0], "E0stat/D");
  Ts->Branch("E1stat", &Estat[1], "E1stat/D");
  Ts->Branch("E2stat", &Estat[2], "E2stat/D");
  Lsidis sidis;
  TLorentzVector l(0, 0, Ebeam, Ebeam);
  TLorentzVector P(0, 0, 0, 0.938272);
  sidis.SetNucleus(Np,Nn);
  sidis.SetHadron(had);
  if (strcmp(had, "pi+") == 0 || strcmp(had, "pi-") == 0) sidis.ChangeTMDpars(0.604, 0.114);
  if (strcmp(had, "K+") == 0 || strcmp(had, "K-") == 0) sidis.ChangeTMDpars(0.604, 0.131);
  sidis.SetInitialState(l, P);
  sidis.SetPDFset("CJ15lo");
  sidis.SetFFset("DSSFFlo");
  double lumi = 1.0e+9 * pow(0.197327, 2);
  double eff = 0.85;
  double time = 55.0 * 24.0 * 3600.0;
  if (Ebeam < 10.0) time = 27.5 * 24.0 * 3600.0;
  Long64_t Nsim = 0;
  Long64_t Nrec = 0;
  double Xmin[6] = {0.0, 0.0, 0.0, 0.0, -M_PI, -M_PI}; 
  double Xmax[6] = {0.7, 0.0, 0.0, 0.0, M_PI, M_PI};;//x, Q2, z, Pt, phih, phiS
  double weight = 0;
  double weight_p = 0;
  double acc = 0;
  TLorentzVector lp(0, 0, 0, 0);
  TLorentzVector Ph(0, 0, 0, 0);
  Lsidis sidis_p;
  sidis_p.SetNucleus(1, 0);
  sidis_p.SetHadron(had);
  if (strcmp(had, "pi+") == 0 || strcmp(had, "pi-") == 0) sidis_p.ChangeTMDpars(0.604, 0.114);
  if (strcmp(had, "K+") == 0 || strcmp(had, "K-") == 0) sidis_p.ChangeTMDpars(0.604, 0.131);
  sidis_p.SetInitialState(l, P);
  sidis_p.SetPDFset("CJ15lo");
  sidis_p.SetFFset("DSSFFlo");
  ifstream infile(readfile);
  char tmp[300];
  infile.getline(tmp, 256);
  int Nt = 0;
  while (infile >> Xmin[1] >> Xmax[1] >> Xmin[2] >> Xmax[2] >> Xmin[3] >> Xmax[3] >> Xmin[0] >> Xmax[0]){
    printf("%.4d  Q2[%.1f,%.1f]  z[%.2f,%.2f]  Pt[%.1f,%.1f]  x[%.4f,%.4f]\n",
	   Nt++, Xmin[1], Xmax[1], Xmin[2], Xmax[2], Xmin[3], Xmax[3], Xmin[0], Xmax[0]);
    sidis.SetRange(Xmin, Xmax);
    sidis_p.SetRange(Xmin, Xmax);
    TH1D * hvar = new TH1D("hvar", "hvar", 7, -0.5, 6.5);
    TH2D * hs = new TH2D("hs", "hs", 36, -M_PI, M_PI, 18, 0, M_PI);
    Nsim = 0;
    Nrec = 0;
    for (Long64_t i = 0; i < 1.0e7; i++){
      Nsim++;
      weight = sidis.GenerateEvent(0, 1);
      if (weight > 0){
        if (sidis.GetVariable("W") < 2.3) continue;
        if (sidis.GetVariable("Wp") < 1.6) continue;
	sidis.CalculateRfactor();
	if (sidis.GetVariable("Rfactor") > Rfactor0) continue;
	lp = sidis.GetLorentzVector("lp");
	Ph = sidis.GetLorentzVector("Ph");
	acc = GetAcceptance_e(lp) * GetAcceptance_hadron(Ph, had);
	if (acc > 0){
	  sidis_p.SetFinalState(lp, Ph);
	  sidis_p.CalculateVariables();
	  weight_p = sidis_p.GetEventWeight(0, 1);
	  Nrec++;
	  hvar->Fill(0., weight_p * acc);
	  hvar->Fill(1., weight * acc);
	  hvar->Fill(2., weight * acc * sidis.GetVariable("x"));
	  hvar->Fill(3., weight * acc * sidis.GetVariable("y"));
	  hvar->Fill(4., weight * acc * sidis.GetVariable("z"));
	  hvar->Fill(5., weight * acc * sidis.GetVariable("Q2"));
	  hvar->Fill(6., weight * acc * sidis.GetVariable("Pt"));
	  hs->Fill(sidis.GetVariable("phih"), std::abs(sidis.GetVariable("phiS")), weight * acc);
	}
      }
      if (Nrec > 500000) break;
    }
    hvar->Scale(lumi * time * eff / Nsim);
    hs->Scale(lumi * time * eff / Nsim);
    Nacc = hvar->GetBinContent(2);
    fp = hvar->GetBinContent(1) / Nacc;
    x = hvar->GetBinContent(3) / Nacc;
    y = hvar->GetBinContent(4) / Nacc;
    z = hvar->GetBinContent(5) / Nacc;
    Q2 = hvar->GetBinContent(6) / Nacc;
    Pt = hvar->GetBinContent(7) / Nacc;
    TMatrixD MUT3(3,3);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
	MUT3(i,j) = 0.0;
    for (int i = 1; i <= 36; i++){
      for (int j = 1; j <= 18; j++){
	phih = hs->GetXaxis()->GetBinCenter(i);
	phiS = hs->GetYaxis()->GetBinCenter(j);
	MUT3(0,0) += sin(phih - phiS) * sin(phih - phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
	MUT3(0,1) += sin(phih - phiS) * sin(phih + phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
	MUT3(0,2) += sin(phih - phiS) * sin(3.0 * phih - phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
	MUT3(1,0) += sin(phih + phiS) * sin(phih - phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
	MUT3(1,1) += sin(phih + phiS) * sin(phih + phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
	MUT3(1,2) += sin(phih + phiS) * sin(3.0 * phih - phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
	MUT3(2,0) += sin(3.0 * phih - phiS) * sin(phih - phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
	MUT3(2,1) += sin(3.0 * phih - phiS) * sin(phih + phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
	MUT3(2,2) += sin(3.0 * phih - phiS) * sin(3.0 * phih - phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
      }
    }
    MUT3.Invert();
    for (int i = 0; i < 3; i++){
      Estatraw[i] = sqrt(2.0 * M_PI * M_PI / Nacc * (pow(MUT3(i,0),2) + pow(MUT3(i,1), 2) + pow(MUT3(i,2), 2)) * M_PI * M_PI);
      Estat[i] = Estatraw[i] / fp / 0.7;
      if (isnan(Estat[i]))
	std::cout << "NaN warning in Estat!" << std::endl;
    }
    Ts->Fill();
    hvar->Delete();
    hs->Delete();
  }
  fs->Write();
  infile.close();
  return 0;
}

  
double CheckCurrentCut(const double Ebeam, const char * hadron, const double kT2 = 0.5, const double MiT2 = 0.5, const double MfT2 = 0.5, const char * plotname = 0){

/*
  Lsidis sidis;
  TLorentzVector l(0, 0, Ebeam, Ebeam);
  TLorentzVector P(0, 0, 0, 0.938272);
  sidis.SetNucleus(Np,Nn);
  sidis.SetHadron(hadron);
  if (strcmp(hadron, "pi+") == 0 || strcmp(hadron, "pi-") == 0) sidis.ChangeTMDpars(0.604, 0.114);
  if (strcmp(hadron, "K+") == 0 || strcmp(hadron, "K-") == 0) sidis.ChangeTMDpars(0.604, 0.131);
  sidis.SetInitialState(l, P);
  sidis.SetPDFset("CJ15lo");
  sidis.SetFFset("DSSFFlo");
  double lumi = 1.0e+9 * pow(0.197327, 2);
  double Nsim = 1.0e7;
  TH2D * h0 = new TH2D("h0", "", 1, 0.2, 0.8, 1, 0.0, 1.6);
  h0->GetXaxis()->SetTitle("z");
  h0->GetXaxis()->CenterTitle(true);
  h0->GetXaxis()->SetTitleSize(0.05);
  h0->GetXaxis()->SetTitleOffset(1.15);
  h0->GetXaxis()->SetLabelSize(0.055);
  h0->GetYaxis()->SetTitle("P_{hT} / GeV");
  h0->GetYaxis()->CenterTitle(true);
  h0->GetYaxis()->SetTitleSize(0.05);
  h0->GetYaxis()->SetTitleOffset(1.15);
  h0->GetYaxis()->SetLabelSize(0.055);
  TH2D * hall = new TH2D("hall", "Before cut", 60, 0.2, 0.8, 160, 0.0, 1.6);
  TH2D * hcut = new TH2D("hcut", "After cut", 60, 0.2, 0.8, 160, 0.0, 1.6);
  double Xmin[6] = {0.0, 1.0, 0.3, 0.0, -M_PI, -M_PI};
  double Xmax[6] = {0.7, 8.0, 0.7, 1.6, M_PI, M_PI};
  sidis.SetRange(Xmin, Xmax);
  double weight = 0;
  double acc = 0;
  TLorentzVector lp, Ph;
  for (Long64_t i = 0; i < Nsim; i++){
    weight = sidis.GenerateEvent(0, 1);
    if (weight > 0){
      if (sidis.GetVariable("W") < 2.3) continue;
      if (sidis.GetVariable("Wp") < 1.6) continue;
      lp = sidis.GetLorentzVector("lp");
      Ph = sidis.GetLorentzVector("Ph");
      acc = GetAcceptance_e(lp) * GetAcceptance_hadron(Ph, hadron);
      if (acc > 0){
	hall->Fill(sidis.GetVariable("z"), sidis.GetVariable("Pt"), weight * acc);
	sidis.CalculateRfactor(kT2, MiT2, MfT2);
	if (sidis.GetVariable("Rfactor") < 0.4){
	  hcut->Fill(sidis.GetVariable("z"), sidis.GetVariable("Pt"), weight * acc);
	}
      }
    }
  }
  hall->Scale(lumi/Nsim);
  hcut->Scale(lumi/Nsim);
  double rate = hcut->Integral(1, -1);
  std::cout << "All: " << hall->Integral(1, -1) << "   Cut: " << hcut->Integral(1, -1) << std::endl;
  if (plotname != 0){
    gStyle->SetOptStat(0);
    //hall->GetZaxis()->SetRangeUser(0.01, hall->GetMaximum()/0.95);
    hcut->GetZaxis()->SetRangeUser(0.01, hall->GetMaximum());
    TCanvas * c0 = new TCanvas("c0", "", 1600, 600);
    c0->SetBorderMode(0);
    c0->SetBorderSize(2);
    c0->SetFrameBorderMode(0);
    c0->Divide(2, 1);
    c0->cd(1);
    c0->cd(1)->SetLeftMargin(0.15);
    c0->cd(1)->SetBottomMargin(0.15);
    h0->Draw();
    hall->Draw("samecolz");
    c0->cd(2);
    c0->cd(2)->SetLeftMargin(0.15);
    c0->cd(2)->SetBottomMargin(0.15);
    h0->Draw();
    hcut->Draw("samecolz");
    c0->Print(plotname);
    c0->Close();
  }
  h0->Delete();
  hall->Delete();
  hcut->Delete();
*/
  double rate = 0;
  return rate;

}

int CreateFileSivers(const char * rootfile1, const char * rootfile2, const char * csvfile){//Create file for Sivers analysis use
  TChain * Ts = new TChain("data", "data");
  Ts->Add(rootfile1);
  Ts->Add(rootfile2);
  double Nucleon, Hadron, Ebeam, x, y, z, Q2, Pt, stat, systrel, systabs, fp;
  Ts->SetBranchAddress("Nucleon", &Nucleon);
  Ts->SetBranchAddress("Hadron", &Hadron);
  Ts->SetBranchAddress("Ebeam", &Ebeam);
  Ts->SetBranchAddress("x", &x);
  Ts->SetBranchAddress("y", &y);
  Ts->SetBranchAddress("z", &z);
  Ts->SetBranchAddress("Q2", &Q2);
  Ts->SetBranchAddress("Pt", &Pt);
  Ts->SetBranchAddress("E1stat", &stat);
  Ts->SetBranchAddress("fp", &fp);
  FILE * file = fopen(csvfile, "w");
  fprintf(file, "i,Ebeam,x,y,z,Q2,pT,obs,value,stat,systrel,systabs,target,hadron,Experiment\n");
  for (int i = 0; i < Ts->GetEntries(); i++){
    std::cout << i << std::endl;
    Ts->GetEntry(i);
    systrel = 0.0;
    systabs = 0.0;
    systrel += pow(0.03, 2);//target polarization
    systrel += pow(0.05, 2);//dilution
    systrel += pow(0.025, 2);//radiative correction
    systrel += pow(0.03, 2);//diffractive meson
    systrel += pow(0.002, 2);//random coincidence
    if (Ebeam > 10.0)//raw asymmetry
      systabs += 7.78e-4 / 0.7 / fp;
    else
      systabs += 1.1e-3 / 0.7 / fp;
    systrel = sqrt(systrel);
    if (Hadron == 0)
      fprintf(file, "%d,%.1f,%.6f,%.6f,%.6f,%.6f,%.6f,%s,%.1f,%.6f,%.6f,%.6f,%s,%s,%s\n",
	      i, Ebeam, x, y, z, Q2, Pt, "AUT", 0.0, stat, systrel, systabs, "proton", "pi+", "solid");
    else if (Hadron == 1)
      fprintf(file, "%d,%.1f,%.6f,%.6f,%.6f,%.6f,%.6f,%s,%.1f,%.6f,%.6f,%.6f,%s,%s,%s\n",
	      i, Ebeam, x, y, z, Q2, Pt, "AUT", 0.0, stat, systrel, systabs, "proton", "pi-", "solid");
  }
  fclose(file);
  return 0;
}







#endif
