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


void track_analysis() {

//gStyle->SetOptStat(0); //turn on/off the histogram statistic
gStyle->SetOptFit(1);
bool draw=true; //drawing option

bool KFREF=false;
bool DAFREF=true;
bool Test_file=true;
    
//add the file
TChain *ttree;
ttree = new TChain("Truth");
char name_file[500];

//KFREF:
if (KFREF)
{

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

}

//DAFREF:
if (DAFREF)
{

for (int ii=1; ii<=1000; ii++)
{
 sprintf(name_file,"/pnfs/e906/persistent/users/zakbar/SimChainDev/july17_100K/%i/out/trk_eval.root",ii);
ttree->Add(name_file);
cerr<<name_file<<endl;
}

}


//Test_file:
if (Test_file)
{

for (int ii=1; ii<=50; ii++)
{

sprintf(name_file,"/e906/app/users/zakbar/July2020/e1039-analysis/SimChainDev/scratch/july15_5K/%i/out/trk_eval.root",ii);
ttree->Add(name_file);
cerr<<name_file<<endl;
}

for (int ii=1; ii<=50; ii++)
{

sprintf(name_file,"/e906/app/users/zakbar/July2020/e1039-analysis/SimChainDev/scratch/july15_5K_v2/%i/out/trk_eval.root",ii);
ttree->Add(name_file);
cerr<<name_file<<endl;
}

}

TTree *T = (TTree*)gROOT->FindObject("Truth");

// variables
unsigned short emu_trigger = 0;
int n_tracks = 0;
int gnhodo[100];
//float eval_var[100];

#define MAX_PARTICLES 100
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

float uncSt1_0[100];
float uncSt1_1[100];
float uncSt1_2[100];
float uncSt1_3[100];
float uncSt1_4[100];
//


int nhits[100];
int gnhits[100];

TMatrixD covar_matrix[100];

//selected branch

  T->SetBranchAddress("emu_trigger", &emu_trigger);
  T->SetBranchAddress("n_tracks",    &n_tracks);
  T->SetBranchAddress("gnhodo",      &gnhodo);
  T->SetBranchAddress("dimu_gphi",   &dimu_gphi);
  T->SetBranchAddress("dimu_phi",   &dimu_phi);
  //T->SetBranchAddress(var_name,      &eval_var);
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
//  T->SetBranchAddress("z_st1", &z_st1);
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
  T->SetBranchAddress("covar_matrix", &covar_matrix);

  T->SetBranchAddress("uncSt1_0", &uncSt1_0);
  T->SetBranchAddress("uncSt1_1", &uncSt1_1);
  T->SetBranchAddress("uncSt1_2", &uncSt1_2);
  T->SetBranchAddress("uncSt1_3", &uncSt1_3);
  T->SetBranchAddress("uncSt1_4", &uncSt1_4);


  
//  Histogram definition here

// histogram keep for future study
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

TH1D *dmom_st1 = new TH1D("dmom_st1", "dmom_st1", 40,-5,5);
dmom_st1->Sumw2();

TH1D *dmom_st3 = new TH1D("dmom_st3", "dmom_st3", 40,-5,5);
dmom_st3->Sumw2();

TH1D *h_kick = new TH1D("h_kick","h_kick",40,0.35,0.45);
h_kick->Sumw2();

TH1D *dmom_st13 = new TH1D("dmom_st13", "dmom_st13", 40,-0.1,0.1);
dmom_st13->Sumw2();

//current histogram
TH1D *h_dreco_p_st1 = new TH1D("h_dreco_p_st1", "h_dreco_p_st1", 50,-1,1);
h_dreco_p_st1->Sumw2();

TH1D *h_dtracklet_p_st1 = new TH1D("h_dtracklet_p_st1", "h_dtracklet_p_st1", 50,-5,5);
h_dtracklet_p_st1->Sumw2();

TH1D *h_dreco_x_st1 = new TH1D("h_dreco_x_st1", "h_dreco_x_st1", 50,-0.002,0.002);
h_dreco_x_st1->Sumw2();

TH1D *h_dtracklet_x_st1 = new TH1D("h_dtracklet_x_st1", "h_dtracklet_x_st1", 50,-0.1,0.1);
h_dtracklet_x_st1->Sumw2();

TH1D *h_dreco_y_st1 = new TH1D("h_dreco_y_st1", "h_dreco_y_st1", 50,-0.1,0.1);
h_dreco_y_st1->Sumw2();

TH1D *h_dtracklet_y_st1 = new TH1D("h_dtracklet_y_st1", "h_dtracklet_y_st1", 50,-0.1,0.1);
h_dtracklet_y_st1->Sumw2();

//TH1D *h_dreco_x_st1 = new TH1D("h_dreco_x_st1", "h_dreco_x_st1", 50,-0.1,0.1);
//h_dreco_x_st1->Sumw2();

//
for(int ientry=0; ientry<T->GetEntries(); ++ientry) 
{ 
 T->GetEntry(ientry);
 cerr<<ientry<<" "<<krecstat<<endl;
 int ntrack_hodo = 0;
 for (int i=0; i<n_tracks; ++i) 
  {
    if(gnhodo[i]>7) ++ntrack_hodo; //ntrack in hodoscope, keep for future study, not used for now
  }

 
        

 if(n_tracks > 0)
  {
   for (int i=0; i<n_tracks; i++)
    {
     if (fabs(charge[i]) == 1) // if the track is not successfully reconstructed, charge show big number
     {
    
     double truth_p_st1 = TMath::Sqrt(gpx_st1[i]*gpx_st1[i]+gpy_st1[i]*gpy_st1[i]+gpz_st1[i]*gpz_st1[i]);
     double truth_x_st1 = gx_st1[i];
     double truth_y_st1 = gy_st1[i];
 
     double reco_p_st1 = 1.0/fabs(stateSt1_0[i]);
     double reco_x_st1 = stateSt1_3[i];
     double reco_y_st1 = stateSt1_4[i];

     //remember that the information (fit parameter) of tracklet are in St-3. but ty,y0 for st1 and st2 are the same
     double tracklet_p_st1 = tracklet_mom_st1[i];
     double PT_KICK_KMAG = 0.951*0.4016;
     double Z_KMAG_BEND = 1064.26;
     double tracklet_tx_st1 = tracklet_tx[i] + PT_KICK_KMAG*tracklet_invp[i]*charge[i];
     double tracklet_x0_st1 = tracklet_tx[i]*Z_KMAG_BEND + tracklet_x0[i] - tracklet_tx_st1*Z_KMAG_BEND;
     double tracklet_x_st1 = tracklet_x0_st1 + tracklet_tx_st1*617.274; //position of the D0X
     double tracklet_y_st1 = tracklet_y0[i] + tracklet_ty[i]*617.274; 

     h_dreco_p_st1->Fill((truth_p_st1 - reco_p_st1));
     h_dreco_x_st1->Fill((truth_x_st1 - reco_x_st1));
     h_dreco_y_st1->Fill((truth_y_st1 - reco_y_st1));

     h_dtracklet_p_st1->Fill(truth_p_st1 - tracklet_p_st1);
     h_dtracklet_x_st1->Fill(truth_x_st1 - tracklet_x_st1);
     h_dtracklet_y_st1->Fill(truth_y_st1 - tracklet_y_st1);
	cerr<<truth_p_st1<<" "<<reco_p_st1<<endl;
    

     


     }
    }
  } 
} 
//end of loop

TFile g ("histo_track.root", "RECREATE", "Histograms from ntuples" );
if (draw)
{

//reserve if we want to fit the histogram
/*TF1* gauss = new TF1("gauss","gaus",-3,3);
h_dreco_p_st1_mom->Fit("gauss","R");
h_dreco_p_st1->GetFunction("gauss")->SetLineColor(2);*/



/*TCanvas *c1 = new TCanvas("c1","c1"); c1->SetGrid();
c1->cd();  
	h_dreco_p_st1->SetTitle("#Delta P at st1 ;#Delta P; N");
			h_dreco_p_st1->SetMarkerColor(kRed);
                        h_dreco_p_st1->SetMarkerStyle(20);
                        h_dreco_p_st1->SetMarkerSize(1);
			h_dreco_p_st1->SetLineColor(kRed);
			h_dreco_p_st1->Draw("e");
			c1->SaveAs("hist_dreco_p_st1.png");*/

TCanvas *c2 = new TCanvas("c2","c2"); c2->SetGrid();
c2->cd();  
	h_dreco_x_st1->SetTitle("#Delta X at st1 ;#Delta X; N");
			h_dreco_x_st1->SetMarkerColor(kRed);
                        h_dreco_x_st1->SetMarkerStyle(20);
                        h_dreco_x_st1->SetMarkerSize(1);
			h_dreco_x_st1->SetLineColor(kRed);
			h_dreco_x_st1->Draw("e");
			c2->SaveAs("hist_dreco_x_st1.png");

TCanvas *c3 = new TCanvas("c3","c3"); c3->SetGrid();
c3->cd();  
	h_dreco_y_st1->SetTitle("#Delta Y at st1 ;#Delta Y; N");
			h_dreco_y_st1->SetMarkerColor(kRed);
                        h_dreco_y_st1->SetMarkerStyle(20);
                        h_dreco_y_st1->SetMarkerSize(1);
			h_dreco_y_st1->SetLineColor(kRed);
			h_dreco_y_st1->Draw("e");
			c3->SaveAs("hist_dreco_y_st1.png");

/*TCanvas *c4 = new TCanvas("c4","c4"); c4->SetGrid();
c4->cd();  
	h_dtracklet_p_st1->SetTitle("#Delta P at st1 ;#Delta P; N");
			h_dtracklet_p_st1->SetMarkerColor(kRed);
                        h_dtracklet_p_st1->SetMarkerStyle(20);
                        h_dtracklet_p_st1->SetMarkerSize(1);
			h_dtracklet_p_st1->SetLineColor(kRed);
			h_dtracklet_p_st1->Draw("e");
			c4->SaveAs("hist_dtracklet_p_st1.png");*/

TCanvas *c5 = new TCanvas("c5","c5"); c5->SetGrid();
c5->cd();  
	h_dtracklet_x_st1->SetTitle("#Delta X at st1 ;#Delta X; N");
			h_dtracklet_x_st1->SetMarkerColor(kRed);
                        h_dtracklet_x_st1->SetMarkerStyle(20);
                        h_dtracklet_x_st1->SetMarkerSize(1);
			h_dtracklet_x_st1->SetLineColor(kRed);
			h_dtracklet_x_st1->Draw("e");
			c5->SaveAs("hist_dtracklet_x_st1.png");

TCanvas *c6 = new TCanvas("c6","c6"); c6->SetGrid();
c6->cd();  
	h_dtracklet_y_st1->SetTitle("#Delta Y at st1 ;#Delta Y; N");
			h_dtracklet_y_st1->SetMarkerColor(kRed);
                        h_dtracklet_y_st1->SetMarkerStyle(20);
                        h_dtracklet_y_st1->SetMarkerSize(1);
			h_dtracklet_y_st1->SetLineColor(kRed);
			h_dtracklet_y_st1->Draw("e");
			c6->SaveAs("hist_dtracklet_y_st1.png");


}

}
                      

