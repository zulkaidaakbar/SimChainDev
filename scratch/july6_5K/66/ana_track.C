#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include "TProfile.h"

using namespace std;
//#define PI 3.14159265359


void ana_track(
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

for (int ii=1; ii<=1000; ii++)
{
 sprintf(name_file,"/pnfs/e906/persistent/users/zakbar/SimChainDev/july2_100K/%i/out/trk_eval.root",ii);
ttree->Add(name_file);
cerr<<name_file<<endl;
}

for (int ii=1; ii<=80; ii++)
{
 sprintf(name_file,"/e906/app/users/zakbar/July2020/e1039-analysis/SimChainDev/scratch/july2_5K/%i/out/trk_eval.root",ii);
ttree->Add(name_file);
cerr<<name_file<<endl;
}

//ttree->Add("trk_aug29_8M.root");
//ttree->Add("old_trk.root");

//ttree->Add("june12_15K.root"); //last one before meeting june17
//ttree->Add("june18test.root");

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

//for tracking study
  
//  Histogram definition here
TH1D *pull_mom = new TH1D("pull_mom","pull_mom",100,-5,5);
pull_mom->Sumw2();

TH1D *pull_x = new TH1D("pull_x","pull_x",100,-5,5);
pull_x->Sumw2();

TH1D *pull_y = new TH1D("pull_y","pull_y",100,-5,5);
pull_y->Sumw2();

TH1D *pull_tx = new TH1D("pull_tx","pull_tx",100,-5,5);
pull_tx->Sumw2();

TH1D *pull_ty = new TH1D("pull_ty","pull_ty",100,-5,5);
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

TH2F *p_dist = new TH2F("p_dist","p_dist",16,10,90,200,-5,5);
p_dist->Sumw2();

TProfile* hprof;
hprof  = new TProfile("hprof","Profile of Delta P versus P",16,10,90,-5,5);

  for(int ientry=0; ientry<T->GetEntries(); ++ientry) 
{ //start of loop
 T->GetEntry(ientry);
 int ntrack_hodo = 0;
    for (int i=0; i<n_tracks; ++i) {
      if(gnhodo[i]>7) ++ntrack_hodo;
    }
 cerr<<ientry<<" "<<T->GetEntries()<<" "<<" "<<krecstat<<endl;

 if(n_tracks > 0)
  {
   for (int i=0; i<n_tracks; i++)
    {
     double charge_d;
     charge_d = charge[i];
     double mom_truth = TMath::Sqrt(gpx_st1[i]*gpx_st1[i]+gpy_st1[i]*gpy_st1[i]+gpz_st1[i]*gpz_st1[i]);
     double mom_reco = TMath::Sqrt(px_st1[i]*px_st1[i]+py_st1[i]*py_st1[i]+pz_st1[i]*pz_st1[i]);
     double pull_mom_truth = charge_d/mom_truth;
     double pull_mom_reco = stateSt1_0[i];

     double pull_tx_truth = gpx_st1[i]/gpz_st1[i];
     double pull_tx_reco = stateSt1_1[i];

     double pull_ty_truth = gpy_st1[i]/gpz_st1[i];
     double pull_ty_reco = stateSt1_2[i];

     double pull_x_truth = gx_st1[i] + pull_tx_truth*(0.0 - gz_st1[i]);
     double pull_x_reco = stateSt1_3[i];

     double pull_y_truth =  gy_st1[i] + pull_ty_truth*(0.0 - gz_st1[i]);
     double pull_y_reco = stateSt1_4[i];

     double delta_mom = pull_mom_truth - pull_mom_reco;
     double delta_tx = pull_tx_truth - pull_tx_reco;
     double delta_ty = pull_ty_truth - pull_ty_reco;
     double delta_x = pull_x_truth - pull_x_reco;
     double delta_y = pull_y_truth - pull_y_reco;

     double cl_st1 = TMath::Prob(chisqSt1[i], nhitsSt1u[i]-5);

if (charge[i] < 2)
     {
      double err_c0 = TMath::Sqrt(fabs(covSt1_0[i]));
      double err_c1 = TMath::Sqrt(fabs(covSt1_1[i]));
      double err_c2 = TMath::Sqrt(fabs(covSt1_2[i]));
      double err_c3 = TMath::Sqrt(fabs(covSt1_3[i]));
      double err_c4 = TMath::Sqrt(fabs(covSt1_4[i]));

      inv_sigma_mom = 1.0/err_c0;
      inv_sigma_tx = 1.0/err_c1;
      inv_sigma_ty = 1.0/err_c2;
      inv_sigma_x = 1.0/err_c3;
      inv_sigma_y = 1.0/err_c4;

      prob_st1->Fill(cl_st1);

      
     if (probability[i] > 0.0001)
      {
      prob->Fill(probability[i]);
      }
     if (delta_mom != 0)
      {
        //pull_mom->Fill(inv_sigma_mom*(pull_mom_truth - pull_mom_reco));
        pull_mom->Fill(mom_truth - fabs(1.0/pull_mom_reco));
        h_cov0->Fill(covSt1_0[i]);
        p_dist->Fill(mom_truth,(mom_truth - 1.0/fabs(pull_mom_reco)));
        hprof->Fill(mom_truth,(mom_truth - 1.0/fabs(pull_mom_reco)));
       }
     if (delta_tx != 0)
      {  
       pull_tx->Fill(inv_sigma_tx*(pull_tx_truth - pull_tx_reco));
       h_cov1->Fill(covSt1_1[i]);
      }
     if (delta_ty != 0)
      {
      pull_ty->Fill(inv_sigma_ty*(pull_ty_truth - pull_ty_reco));
      h_cov2->Fill(covSt1_2[i]);
      }
      if (delta_x != 0)
      {
      // pull_x->Fill(inv_sigma_x*(pull_x_truth - pull_x_reco));
        pull_x->Fill(inv_sigma_x*(gx_st1[i] - pull_x_reco));
       
      }
     if (delta_y != 0)
      {
     // pull_y->Fill(inv_sigma_y*(pull_y_truth - pull_y_reco));
        pull_y->Fill(inv_sigma_y*(gy_st1[i] - pull_y_reco));
      
      }
      // cerr<<gx_st1[i]<<" "<<x_st1[i]<<" "<<pull_x_reco<<endl;
      // cerr<<ientry<<" "<<chisquare[i]<<" "<<chisqSt1[i]<<" "<<err_c0<<" "<<err_c1<<" "<<err_c2<<" "<<err_c3<<" "<<err_c4<<endl; 
     // cerr<<ientry<<" "<<T->GetEntries()<<" "<<gnhits[i]<<" "<<nhits[i]<<" "<<nhitsSt1u[i]<<" "<<chisqSt1[i]<<" "<<TMath::Prob(chisqSt1[i], nhitsSt1u[i]-5)<<" "<<gx_st1[i]<<" "<<gy_st1[i]<<" "<<gz_st1[i]<<endl;    
    // cerr<<ientry<<" "<<T->GetEntries()<<" "<<ntrack_hodo<<" "<<krecstat<<" "<<delta_mom<<" "<<delta_tx<<" "<<delta_ty<<" "<<" "<<charge[i]<<" "<<covSt1_0[i]<<" "<<covSt1_1[i]<<" "<<covSt1_2[i]<<" "<<probability[i]<<" "<<chisquare[i]<<endl;
     }

    }
    cerr<<endl;
  }
    
    
} //end of loop

TFile g ("recstat_track.root", "RECREATE", "Histograms from ntuples" );
if (draw)
{

TF1* gauss = new TF1("gauss","gaus",-3,3);
pull_mom->Fit("gauss","R");
pull_mom->GetFunction("gauss")->SetLineColor(2);
pull_tx->Fit("gauss","R");
pull_tx->GetFunction("gauss")->SetLineColor(2);
pull_ty->Fit("gauss","R");
pull_ty->GetFunction("gauss")->SetLineColor(2);
pull_x->Fit("gauss","R");
pull_x->GetFunction("gauss")->SetLineColor(2);
pull_y->Fit("gauss","R");
pull_y->GetFunction("gauss")->SetLineColor(2);

TProfile *profile_p = p_dist->ProfileY("profile_p",0,15);


TCanvas *c1 = new TCanvas("c1","c1"); c1->SetGrid();
c1->cd();  
	pull_mom->SetTitle("q/p pull; q/p pull; N");
			/*pull_mom->SetMarkerColor(kRed);
                        pull_mom->SetMarkerStyle(20);
                        pull_mom->SetMarkerSize(1);
			pull_mom->SetLineColor(kRed);*/
			pull_mom->Draw("e");
			//recstat->SetStats(0);
			c1->SaveAs("pull_mom_covb.png");
TCanvas *c2 = new TCanvas("c2","c2"); c2->SetGrid();
c2->cd();  
	pull_tx->SetTitle("tx pull; tx pull; N");
			/*pull_tx->SetMarkerColor(kRed);
                        pull_tx->SetMarkerStyle(20);
                        pull_tx->SetMarkerSize(1);
			pull_tx->SetLineColor(kRed);*/
			/*pull_tx->Draw("e");
			//recstat->SetStats(0);
			c2->SaveAs("pull_tx_cov.png");*/
TCanvas *c3 = new TCanvas("c3","c3"); c3->SetGrid();
c3->cd();  
	pull_ty->SetTitle("ty pull; ty pull; N");
			/*pull_ty->SetMarkerColor(kRed);
                        pull_ty->SetMarkerStyle(20);
                        pull_ty->SetMarkerSize(1);
			pull_ty->SetLineColor(kRed);*/
			/*pull_ty->Draw("e");
			//recstat->SetStats(0);
			c3->SaveAs("pull_ty_cov.png");*/
/*
TCanvas *c4 = new TCanvas("c4","c4"); c4->SetGrid();
c4->cd();
        prob->SetTitle("Probability; P; N");
                        //pull_ty->SetMarkerColor(kRed);
                        //prob->SetMarkerStyle(20);
                        //prob->SetMarkerSize(1);
                        //prob->SetLineColor(kRed);
			prob->SetMaximum(200);
			prob->GetXaxis()->SetRangeUser(0.05,0.95);
			prob->SetFillColor(kBlue);
                        prob->Draw();
                        //recstat->SetStats(0);
                        c4->SaveAs("probability.png");

TCanvas *c5 = new TCanvas("c5","c5"); c5->SetGrid();
c5->cd();
        h_cov0->SetTitle("cov[0][0] distribution; cov[0][0]; N");
                        //pull_ty->SetMarkerColor(kRed);
                       // pull_ty->SetMarkerStyle(20);
                       // pull_ty->SetMarkerSize(1);
                       // pull_ty->SetLineColor(kRed);
                        h_cov0->Draw("e");
                        //recstat->SetStats(0);
                        c5->SaveAs("h_cov0.png");

TCanvas *c6 = new TCanvas("c6","c6"); c6->SetGrid();
c6->cd();
        h_cov1->SetTitle("cov[1][0] distribution; cov[1][0]; N");
                        //pull_ty->SetMarkerColor(kRed);
                       // pull_ty->SetMarkerStyle(20);
                       // pull_ty->SetMarkerSize(1);
                       // pull_ty->SetLineColor(kRed);
                        h_cov1->Draw("e");
                        //recstat->SetStats(0);
                        c6->SaveAs("h_cov1.png");

TCanvas *c7 = new TCanvas("c7","c7"); c7->SetGrid();
c7->cd();
        h_cov2->SetTitle("cov[2][0] distribution; cov[2][0]; N");
                        //pull_ty->SetMarkerColor(kRed);
                       // pull_ty->SetMarkerStyle(20);
                       // pull_ty->SetMarkerSize(1);
                       // pull_ty->SetLineColor(kRed);
                        h_cov2->Draw("e");
                        //recstat->SetStats(0);
                        c7->SaveAs("h_cov2.png");
*/
TCanvas *c8 = new TCanvas("c8","c8"); c8->SetGrid();
c8->cd();
        pull_x->SetTitle("x0 pull; x0 pull; N");
                        /*pull_tx->SetMarkerColor(kRed);
                        pull_tx->SetMarkerStyle(20);
                        pull_tx->SetMarkerSize(1);
                        pull_tx->SetLineColor(kRed);*/
                        /*pull_x->Draw("e");
                        //recstat->SetStats(0);
                        c8->SaveAs("pull_x_cov.png");*/
TCanvas *c9 = new TCanvas("c9","c9"); c9->SetGrid();
c9->cd();
        pull_y->SetTitle("y0 pull; y0 pull; N");
                        /*pull_ty->SetMarkerColor(kRed);
                        pull_ty->SetMarkerStyle(20);
                        pull_ty->SetMarkerSize(1);
                        pull_ty->SetLineColor(kRed);*/
                        /*pull_y->Draw("e");
                        //recstat->SetStats(0);
                        c9->SaveAs("pull_y_cov.png");*/

TCanvas *c10 = new TCanvas("c10","c10"); c10->SetGrid();
c10->cd();
        prob_st1->SetTitle("Probat St-1; P; N");
                        //pull_ty->SetMarkerColor(kRed);
                        //prob->SetMarkerStyle(20);
                        //prob->SetMarkerSize(1);
                        //prob->SetLineColor(kRed);
                       // prob->SetMaximum(200);
                      /* prob_st1->GetXaxis()->SetRangeUser(0.05,0.95);
                        prob_st1->SetFillColor(kBlue);
                        prob_st1->Draw();
                        //recstat->SetStats(0);
                        c10->SaveAs("probability_st1.png");*/

TCanvas *c11 = new TCanvas("c11","c11"); c11->SetGrid();
c11->cd();
        p_dist->SetTitle("; P; #Delta P");
                        /*pull_ty->SetMarkerColor(kRed);
                        pull_ty->SetMarkerStyle(20);
                        pull_ty->SetMarkerSize(1);
                        pull_ty->SetLineColor(kRed);*/
                        p_dist->Draw("COLZ");
                        //recstat->SetStats(0);
                        c11->SaveAs("p_dist.png");

TCanvas *c12 = new TCanvas("c12","c12"); c12->SetGrid();
c12->cd();
hprof->Draw();
c12->SaveAs("hprof.png");


}


}
                      

