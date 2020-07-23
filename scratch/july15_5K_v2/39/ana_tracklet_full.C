#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"

using namespace std;
//#define PI 3.14159265359


void ana_tracklet_full(
    const char* var_name = "dimu_gphi",
    const char* title = ";dimu gphi;",
    const int nbin = 16,
    const double xmin = -3.14,
    const double xmax = 3.14
    ) {

//gStyle->SetOptStat(0);
gStyle->SetOptFit(1);
    
TChain *ttree;
ttree = new TChain("Truth");

char name_file[500];
bool KF_mode = false;

/*for (int ii=1; ii<=1000; ii++)
{
 sprintf(name_file,"/pnfs/e906/scratch/users/zakbar/SimChainDev/june9_1M/%i/out/trk_eval.root",ii);
ttree->Add(name_file);
cerr<<name_file<<endl;
}*/

//last files before magnet-bug fixed
/*
for (int ii=1; ii<=10; ii++)
{

sprintf(name_file,"/e906/app/users/zakbar/July2020/e1039-analysis/SimChainDev/scratch/july6_1K/%i/out/trk_eval.root",ii);
ttree->Add(name_file);
cerr<<name_file<<endl;
}

for (int ii=1; ii<=100; ii++)
{

sprintf(name_file,"/e906/app/users/zakbar/July2020/e1039-analysis/SimChainDev/scratch/july6_5K/%i/out/trk_eval.root",ii);
ttree->Add(name_file);
cerr<<name_file<<endl;
}
*/

for (int ii=1; ii<=100; ii++)
{

sprintf(name_file,"/e906/app/users/zakbar/July2020/e1039-analysis/SimChainDev/scratch/july7_5K/%i/out/trk_eval.root",ii);
ttree->Add(name_file);
cerr<<name_file<<endl;
}

for (int ii=1; ii<=1000; ii++)
{
 sprintf(name_file,"/pnfs/e906/persistent/users/zakbar/SimChainDev/july8_100K/%i/out/trk_eval.root",ii);
ttree->Add(name_file);
cerr<<name_file<<endl;
}

//ttree->Add("/e906/app/users/zakbar/July2020/e1039-analysis/SimChainDev/scratch/july6test/1/out/trk_eval.root");

cerr<<"additional_file"<<endl;

TTree *T = (TTree*)gROOT->FindObject("Truth");

unsigned short emu_trigger = 0;
  int n_tracks = 0;
  int gnhodo[100];
 // float dimu_gphi[100];
  float eval_var[100];

#define MAX_PARTICLES 100
bool draw=true;

float gpx[MAX_PARTICLES], gpy[MAX_PARTICLES], gpz[MAX_PARTICLES], gvz[MAX_PARTICLES], vz[MAX_PARTICLES];
int n_particles, krecstat;
float dimu_mass, dimu_gmass, dimu_gpz, dimu_pz, dimu_gpt, dimu_pt, dimu_gphi[100], dimu_phi[100];
float gpx_st1[100];
float gpy_st1[100];
float gpz_st1[100];
float gx_st1[100];
float gy_st1[100];
float gz_st1[100];
float px_st1[100];
float py_st1[100];
float pz_st1[100];
float x_st1[100];
float y_st1[100];
float z_st1[100];
int charge[100];
float chisquare[100];
float probability[100];
float covSt1_0[100];
float covSt1_1[100];
float covSt1_2[100];
float covSt1_3[100];
float covSt1_4[100];
float stateSt1_0[100];
float stateSt1_1[100];
float stateSt1_2[100];
float stateSt1_3[100];
float stateSt1_4[100];
float chisqSt1[100];
int nhitsSt1u[100];

float tracklet_tx[100];
float tracklet_ty[100];
float tracklet_x0[100];
float tracklet_y0[100];
float tracklet_invp[100];
float err_tracklet_tx[100];
float err_tracklet_ty[100];
float err_tracklet_x0[100];
float err_tracklet_y0[100];
float err_tracklet_invp[100];

float tracklet_mom_st1[100];
float tracklet_mom_st3[100];

float gpx_st3[100];
float gpy_st3[100];
float gpz_st3[100];

float gx_st3[100];
float gy_st3[100];
float gz_st3[100];
//
double charge_d;
double mom_truth;
     double mom_reco;
     double pull_mom_truth;
     double pull_mom_reco;

     double pull_tx_truth;
     double pull_tx_reco;

     double pull_ty_truth;
     double pull_ty_reco;

     double pull_x_truth;
     double pull_x_reco;

     double pull_y_truth;
     double pull_y_reco;

     double delta_mom;
     double delta_tx;
     double delta_ty;
     double delta_x0;
     double delta_y0;

     double mom_delta_st1;

//double d_tracklet_mom_

//

int nhits[100];
int gnhits[100];

double inv_sigma_mom = 1000.00;
double inv_sigma_tx = 3162.00;
double inv_sigma_ty = 3162.00;
double inv_sigma_x = 3162.00;
double inv_sigma_y = 3162.00;

  T->SetBranchAddress("emu_trigger", &emu_trigger);
  T->SetBranchAddress("n_tracks",    &n_tracks);
  T->SetBranchAddress("gnhodo",      &gnhodo);
  T->SetBranchAddress("dimu_gphi",   &dimu_gphi);
  T->SetBranchAddress("dimu_phi",   &dimu_phi);
  T->SetBranchAddress(var_name,      &eval_var);
  T->SetBranchAddress("gpx", &gpx);
  T->SetBranchAddress("gpy", &gpy);
  T->SetBranchAddress("gpz", &gpz);
  T->SetBranchAddress("gvz", &gvz);
  T->SetBranchAddress("vz", &vz);
  T->SetBranchAddress("krecstat", &krecstat);
  T->SetBranchAddress("dimu_mass", &dimu_mass);
  T->SetBranchAddress("dimu_gmass", &dimu_gmass);
  T->SetBranchAddress("dimu_gpz", &dimu_gpz);
  T->SetBranchAddress("dimu_pz", &dimu_pz);
  T->SetBranchAddress("dimu_gpt", &dimu_gpt);
  T->SetBranchAddress("dimu_pt", &dimu_pt);
//for tracking study
  T->SetBranchAddress("gpx_st1", &gpx_st1);
  T->SetBranchAddress("gpy_st1", &gpy_st1);
  T->SetBranchAddress("gpz_st1", &gpz_st1);
  T->SetBranchAddress("gx_st1", &gx_st1);
  T->SetBranchAddress("gy_st1", &gy_st1);
  T->SetBranchAddress("gz_st1", &gz_st1);
  T->SetBranchAddress("px_st1", &px_st1);
  T->SetBranchAddress("py_st1", &py_st1);
  T->SetBranchAddress("pz_st1", &pz_st1);
  T->SetBranchAddress("x_st1", &x_st1);
  T->SetBranchAddress("y_st1", &y_st1);
  T->SetBranchAddress("z_st1", &z_st1);
  T->SetBranchAddress("charge", &charge);
  T->SetBranchAddress("chisquare", &chisquare);
  T->SetBranchAddress("probability", &probability);
  T->SetBranchAddress("covSt1_0", &covSt1_0);
  T->SetBranchAddress("covSt1_1", &covSt1_1);
  T->SetBranchAddress("covSt1_2", &covSt1_2);
  T->SetBranchAddress("covSt1_3", &covSt1_3);
  T->SetBranchAddress("covSt1_4", &covSt1_4);
  T->SetBranchAddress("stateSt1_0", &stateSt1_0);
  T->SetBranchAddress("stateSt1_1", &stateSt1_1);
  T->SetBranchAddress("stateSt1_2", &stateSt1_2);
  T->SetBranchAddress("stateSt1_3", &stateSt1_3);
  T->SetBranchAddress("stateSt1_4", &stateSt1_4);
  T->SetBranchAddress("chisqSt1", &chisqSt1);
  T->SetBranchAddress("nhitsSt1u", &nhitsSt1u);
  T->SetBranchAddress("nhits", &nhits);
  T->SetBranchAddress("gnhits", &gnhits);
  T->SetBranchAddress("tracklet_tx", &tracklet_tx);
  T->SetBranchAddress("tracklet_ty", &tracklet_ty);
  T->SetBranchAddress("tracklet_x0", &tracklet_x0);
  T->SetBranchAddress("tracklet_y0", &tracklet_y0);
  T->SetBranchAddress("tracklet_invp", &tracklet_invp);

  T->SetBranchAddress("tracklet_mom_st1", &tracklet_mom_st1);
  T->SetBranchAddress("tracklet_mom_st3", &tracklet_mom_st3);

  T->SetBranchAddress("gpx_st3", &gpx_st3);
  T->SetBranchAddress("gpy_st3", &gpy_st3);
  T->SetBranchAddress("gpz_st3", &gpz_st3);
  T->SetBranchAddress("gx_st3", &gx_st3);
  T->SetBranchAddress("gy_st3", &gy_st3);
  T->SetBranchAddress("gz_st3", &gz_st3);

//for tracking study
  
//  Histogram definition here
//TH1D *pull_mom = new TH1D("pull_mom","pull_mom",50,-0.002,0.002);
TH1D *pull_mom = new TH1D("pull_mom","pull_mom",50,-5,5);
pull_mom->Sumw2();

TH1D *pull_x = new TH1D("pull_x","pull_x",50,-0.5,0.5);
pull_x->Sumw2();

TH1D *pull_y = new TH1D("pull_y","pull_y",50,-0.5,0.5);
pull_y->Sumw2();

TH1D *pull_tx = new TH1D("pull_tx","pull_tx",50,-0.02,0.02);
pull_tx->Sumw2();

TH1D *pull_ty = new TH1D("pull_ty","pull_ty",50,-0.02,0.02);
pull_ty->Sumw2();

TH1D *chisq = new TH1D("chisq","chisq",200,0,200);
chisq->Sumw2();

TH1D *prob = new TH1D("prob","prob",50,0,1);
prob->Sumw2();

TH1D *prob_st1 = new TH1D("prob_st1","prob_st1",50,0,1);
prob_st1->Sumw2();

TH1D *h_cov0 = new TH1D("h_cov0","h_cov0",200,-0.1,0.1);
h_cov0->Sumw2();

TH1D *h_cov1 = new TH1D("h_cov1","h_cov1",200,-0.1,0.1);
h_cov1->Sumw2();

TH1D *h_cov2 = new TH1D("h_cov2","h_cov2",200,-0.1,0.1);
h_cov2->Sumw2();

TH1D *dtx_st1 = new TH1D("dtx_st1","dtx_st1",40,-2,2);
dtx_st1->Sumw2();

////

TH1D *dmom_st1 = new TH1D("dmom_st1", "dmom_st1", 40,-5,5);
dmom_st1->Sumw2();

TH1D *dmom_st3 = new TH1D("dmom_st3", "dmom_st3", 40,-5,5);
dmom_st3->Sumw2();

TH1D *h_kick = new TH1D("h_kick","h_kick",40,0.35,0.45);
h_kick->Sumw2();

TH1D *dmom_st13 = new TH1D("dmom_st13", "dmom_st13", 40,-0.1,0.1);
dmom_st13->Sumw2();

//TH1D *



  for(int ientry=0; ientry<T->GetEntries(); ++ientry) 
{ //start of loop
 T->GetEntry(ientry);
 int ntrack_hodo = 0;
    for (int i=0; i<n_tracks; ++i) {
      if(gnhodo[i]>7) ++ntrack_hodo;
    }
 //cerr<<ientry<<" "<<T->GetEntries()<<" "<<" "<<krecstat<<endl;

 if(n_tracks > 0)
  {
   for (int i=0; i<n_tracks; i++)
    {
    
     charge_d = charge[i];
     mom_truth = TMath::Sqrt(gpx_st1[i]*gpx_st1[i]+gpy_st1[i]*gpy_st1[i]+gpz_st1[i]*gpz_st1[i]);
     mom_reco = TMath::Sqrt(px_st1[i]*px_st1[i]+py_st1[i]*py_st1[i]+pz_st1[i]*pz_st1[i]);
     pull_mom_truth = charge_d/mom_truth;
     pull_mom_reco = tracklet_invp[i];

     pull_tx_truth = gpx_st1[i]/gpz_st1[i];
     pull_tx_reco = tracklet_tx[i];

     pull_ty_truth = gpy_st1[i]/gpz_st1[i];
     pull_ty_reco = tracklet_ty[i];

     pull_x_truth = gx_st1[i] + pull_tx_truth*(0.0 - gz_st1[i]);
     pull_x_reco = tracklet_x0[i] + tracklet_tx[i]*617.274;

     pull_y_truth =  gy_st1[i] + pull_ty_truth*(0.0 - gz_st1[i]);
     pull_y_reco = tracklet_y0[i] + tracklet_ty[i]*617.274;

     delta_mom = pull_mom_truth - pull_mom_reco;
     delta_tx = pull_tx_truth - pull_tx_reco;
     delta_ty = pull_ty_truth - pull_ty_reco;
     delta_x0 = pull_x_truth - pull_x_reco;
     delta_y0 = pull_y_truth - pull_y_reco;

     mom_delta_st1 = mom_truth - tracklet_mom_st1[i];
     

     double cl_st1 = TMath::Prob(chisqSt1[i], nhitsSt1u[i]-5);

if (charge[i] < 2)
     {
      double err_c0 = TMath::Sqrt(fabs(covSt1_0[i]));
      double err_c1 = TMath::Sqrt(fabs(covSt1_1[i]));
      double err_c2 = TMath::Sqrt(fabs(covSt1_2[i]));
      double err_c3 = TMath::Sqrt(fabs(covSt1_3[i]));
      double err_c4 = TMath::Sqrt(fabs(covSt1_4[i]));


     /* inv_sigma_mom = 1.0/err_c0;
      inv_sigma_tx = 1.0/err_c1;
      inv_sigma_ty = 1.0/err_c2;
      inv_sigma_x = 1.0/err_c3;
      inv_sigma_y = 1.0/err_c4;*/

      inv_sigma_mom = 1.0;
      inv_sigma_tx = 1.0;
      inv_sigma_ty = 1.0;
      inv_sigma_x = 1.0;
      inv_sigma_y = 1.0;

      prob_st1->Fill(cl_st1);

      
     double truth_slope_st1 = gpx_st1[i]/gpz_st1[i];
     double truth_slope_st3 = gpx_st3[i]/gpz_st3[i];
     double truth_mom_st3 = TMath::Sqrt(gpx_st3[i]*gpx_st3[i]+gpy_st3[i]*gpy_st3[i]+gpz_st3[i]*gpz_st3[i]);
     double truth_mom_st1 = TMath::Sqrt(gpx_st1[i]*gpx_st1[i]+gpy_st1[i]*gpy_st1[i]+gpz_st1[i]*gpz_st1[i]);
     double kick_mom;
     double p_xz = TMath::Sqrt(gpx_st3[i]*gpx_st3[i]+gpz_st3[i]*gpz_st3[i]);
     if (charge_d > 0)
      {
        kick_mom = (truth_slope_st1 - truth_slope_st3)*p_xz;
      }
     else
      {
        kick_mom = (truth_slope_st3 - truth_slope_st1)*p_xz;
      }
   cerr<<charge_d<<" "<<krecstat<<" "<<ntrack_hodo<<"      "<<kick_mom<<endl;
cerr<<gpx_st3[i]<<" "<<gpy_st3[i]<<" "<<gpz_st3[i]<<"    "<<truth_slope_st3<<"    "<<truth_mom_st3<<"  "<<tracklet_mom_st3[i]<<endl;
cerr<<gpx_st1[i]<<" "<<gpy_st1[i]<<" "<<gpz_st1[i]<<"    "<<truth_slope_st1<<"  "<<truth_mom_st1<<"  "<<tracklet_mom_st1[i]<<endl;   

dmom_st1->Fill(mom_delta_st1);
dmom_st3->Fill(truth_mom_st3 - tracklet_mom_st3[i]);
h_kick->Fill(kick_mom);
dmom_st13->Fill(truth_mom_st1 - truth_mom_st3);


    }
    
  }
    
   } 
} //end of loop

TFile g ("recstat_track.root", "RECREATE", "Histograms from ntuples" );
if (draw)
{

/*TF1* gauss = new TF1("gauss","gaus",-3,3);
pull_mom->Fit("gauss","R");
pull_mom->GetFunction("gauss")->SetLineColor(2);
pull_tx->Fit("gauss","R");
pull_tx->GetFunction("gauss")->SetLineColor(2);
pull_ty->Fit("gauss","R");
pull_ty->GetFunction("gauss")->SetLineColor(2);
pull_x->Fit("gauss","R");
pull_x->GetFunction("gauss")->SetLineColor(2);
pull_y->Fit("gauss","R");
pull_y->GetFunction("gauss")->SetLineColor(2);*/


TCanvas *c1 = new TCanvas("c1","c1"); c1->SetGrid();
c1->cd();  
	dmom_st1->SetTitle("#Delta P at st1 ;#Delta P; N");
			dmom_st1->SetMarkerColor(kRed);
                        dmom_st1->SetMarkerStyle(20);
                        dmom_st1->SetMarkerSize(1);
			dmom_st1->SetLineColor(kRed);
			dmom_st1->Draw("e");
			c1->SaveAs("dmom_st1.png");

TCanvas *c2 = new TCanvas("c2","c2"); c2->SetGrid();
c2->cd();  
	dmom_st3->SetTitle("#Delta P at st3 ;#Delta P; N");
			dmom_st3->SetMarkerColor(kRed);
                        dmom_st3->SetMarkerStyle(20);
                        dmom_st3->SetMarkerSize(1);
			dmom_st3->SetLineColor(kRed);
			dmom_st3->Draw("e");
			c2->SaveAs("dmom_st3.png");

TCanvas *c3 = new TCanvas("c3","c3"); c3->SetGrid();
c3->cd();  
	h_kick->SetTitle("#Pt Kick KMag ;#Pt Kick; N");
			h_kick->SetMarkerColor(kRed);
                        h_kick->SetMarkerStyle(20);
                        h_kick->SetMarkerSize(1);
			h_kick->SetLineColor(kRed);
			h_kick->Draw("e");
			c3->SaveAs("KMag_kick.png");

TCanvas *c4 = new TCanvas("c4","c4"); c4->SetGrid();
c4->cd();  
	dmom_st13->SetTitle("#Delta P at st13 ;#Delta P; N");
			dmom_st13->SetMarkerColor(kRed);
                        dmom_st13->SetMarkerStyle(20);
                        dmom_st13->SetMarkerSize(1);
			dmom_st13->SetLineColor(kRed);
			dmom_st13->Draw("e");
			c4->SaveAs("dmom_st13b.png");


}

}
                      

