#include "TEfficiency.h"
#include "TFile.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TTree.h"
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TProfile.h>
#include <TStyle.h>
void DrawSame(){
  TFile *fMC = TFile::Open("Mass_MCmatch_376.root");
  TFile *fKalman = TFile::Open("Mass_MCH_382.root");

  TH1D* hMC = (TH1D*)fMC->Get("LikeSign");
  TH1D* hKalman = (TH1D*)fKalman->Get("LikeSign");

  TH1D* Ratio = new TH1D("Ratio","Ratio",1000,0,10);

  hKalman->SetLineColor(2);

  hMC->Draw("Same");
  hKalman->Draw("Same");
  Ratio->Divide(hKalman, hMC);
  //Ratio->Draw();
}
