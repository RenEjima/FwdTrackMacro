#if !defined(__CLING__) || defined(__ROOTCLING__)

#ifdef __MAKECINT__
//#pragma link C++ class GlobalMuonTrack + ;
//#pragma link C++ class std::vector < GlobalMuonTrack> + ;
//#pragma link C++ class MatchingHelper + ;
#endif

#include "CommonConstants/MathConstants.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include "ReconstructionDataFormats/GlobalFwdTrack.h"
#include "DataFormatsMFT/TrackMFT.h"

#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TTree.h"
#include <TGeoGlobalMagField.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TProfile.h>
#include <TStyle.h>

#endif

using o2::MCTrackT;
using GlobalMuonTrack = o2::dataformats::GlobalFwdTrack;
using eventFoundTracks = std::vector<bool>;
using std::vector;

void LowMassVectorMesonPtEta(){
  const std::string o2sim_KineFile = "bkg_Kine.root";
  TFile *o2sim_KineFileIn = new TFile(o2sim_KineFile.c_str());
  TTree *o2SimKineTree = (TTree *)o2sim_KineFileIn->Get("o2sim");

  const std::string sgn_KineFile = "sgn_1_Kine.root";
  TFile *sgn_KineFileIn = new TFile(sgn_KineFile.c_str());
  TTree *sgnKineTree = (TTree *)sgn_KineFileIn->Get("o2sim");

  std::string outfilename = "LowMassPtEta.root";
  TFile outFile(outfilename.c_str(), "RECREATE");

  TH2F *pT_px_rho = new TH2F("pT_px_rho","#rho meson p_{T}-p_{x}",1000,0,10,1000,0,10);
  TH2F *Eta_py_rho = new TH2F("Eta_py_rho","#rho meson #eta-p_{y}",1000,0,10,1000,-4,4);
  TH2F *pT_px_omega = new TH2F("pT_px_omega","#omega meson p_{T}-p_{x}",1000,0,10,1000,0,10);
  TH2F *Eta_py_omega = new TH2F("Eta_py_omega","#omega meson #eta-p_{y}",1000,0,10,1000,-4,4);
  TH2F *pT_px_phi = new TH2F("pT_px_phi","#phi meson p_{T}-p_{x}",1000,0,10,1000,0,10);
  TH2F *Eta_py_phi = new TH2F("Eta_py_phi","#phi meson #eta-p_{y}",1000,0,10,1000,-4,4);
  TH2F *pT_px_a1 = new TH2F("pT_px_a1","a1 meson p_{T}-p_{x}",1000,0,10,1000,0,10);
  TH2F *Eta_py_a1 = new TH2F("Eta_py_a1","a1 meson #eta-p_{y}",1000,0,10,1000,-4,4);

  TH1F *pT_rho = new TH1F("pT_rho","#rho pT;p_{T}^{MC}[GeV/c];Entry",1000,0,10);
  TH1F *Eta_rho = new TH1F("Eta_rho","#rho #eta;#eta;Entry",1000,-10,10);
  TH1F *pT_omega = new TH1F("pT_omega","#omega pT;p_{T}^{MC}[GeV/c];Entry",1000,0,10);
  TH1F *Eta_omega = new TH1F("Eta_omega","#omega #eta;#eta;Entry",1000,-10,10);
  TH1F *pT_phi = new TH1F("pT_phi","#phi pT;p_{T}^{MC}[GeV/c];Entry",1000,0,10);
  TH1F *Eta_phi = new TH1F("Eta_phi","#phi #eta;#eta;Entry",1000,-10,10);
  TH1F *pT_Jpsi = new TH1F("pT_Jpsi","J/#psi pT;p_{T}^{MC}[GeV/c];Entry",1000,0,10);
  TH1F *Eta_Jpsi = new TH1F("Eta_Jpsi","J#psi #eta;#eta;Entry",1000,-10,10);
  TH1F *pT_a1 = new TH1F("pT_a1","a1 pT;p_{T}^{MC}[GeV/c];Entry",1000,0,10);
  TH1F *Eta_a1 = new TH1F("Eta_a1","a1 #eta;#eta;Entry",1000,-10,10);

  vector<MCTrackT<float>> *mcTr = nullptr;
  o2SimKineTree->SetBranchAddress("MCTrack", &mcTr);

  vector<MCTrackT<float>> *mcTrSgn = nullptr;
  sgnKineTree->SetBranchAddress("MCTrack", &mcTrSgn);

  Int_t numberOfEvents = o2SimKineTree->GetEntries();

  for (int iEvent = 0; iEvent < numberOfEvents; iEvent++) {
    std::cout<<"Start Event #"<<iEvent<<std::endl;
    o2SimKineTree->GetEntry(iEvent);
    int nMCTracks = (*mcTr).size();

    sgnKineTree->GetEntry(iEvent);
    int nMCTracksSgn = (*mcTrSgn).size();

    for (Int_t iMCTrack=0; iMCTrack<(*mcTr).size();++iMCTrack) {
        MCTrackT<float> *thisTrack = &(*mcTr).at(iMCTrack);
        if (thisTrack->GetPdgCode()==113) {
          pT_rho->Fill(thisTrack->GetPt());
          Eta_rho->Fill(thisTrack->GetEta());
          pT_px_rho->Fill(thisTrack->GetStartVertexMomentumX(),thisTrack->GetPt());
          Eta_py_rho->Fill(thisTrack->GetStartVertexMomentumY(),thisTrack->GetEta());
        }
        if (thisTrack->GetPdgCode()==223) {
          pT_omega->Fill(thisTrack->GetPt());
          Eta_omega->Fill(thisTrack->GetEta());
          pT_px_omega->Fill(thisTrack->GetStartVertexMomentumX(),thisTrack->GetPt());
          Eta_py_omega->Fill(thisTrack->GetStartVertexMomentumY(),thisTrack->GetEta());
        }
        if (thisTrack->GetPdgCode()==333) {
          pT_phi->Fill(thisTrack->GetPt());
          Eta_phi->Fill(thisTrack->GetEta());
          pT_px_phi->Fill(thisTrack->GetStartVertexMomentumX(),thisTrack->GetPt());
          Eta_py_phi->Fill(thisTrack->GetStartVertexMomentumY(),thisTrack->GetEta());
        }
        if (thisTrack->GetPdgCode()==443) {
          pT_Jpsi->Fill(thisTrack->GetPt());
          Eta_Jpsi->Fill(thisTrack->GetEta());
        }
        if (thisTrack->GetPdgCode()==20133) {
          pT_a1->Fill(thisTrack->GetPt());
          Eta_a1->Fill(thisTrack->GetEta());
          pT_px_a1->Fill(thisTrack->GetStartVertexMomentumX(),thisTrack->GetPt());
          Eta_py_a1->Fill(thisTrack->GetStartVertexMomentumY(),thisTrack->GetEta());
        }
    }
    for (Int_t iMCTrack=0; iMCTrack<(*mcTrSgn).size();++iMCTrack) {
      MCTrackT<float> *thisTrackSgn = &(*mcTrSgn).at(iMCTrack);
        if (thisTrackSgn->GetPdgCode()==113) {
          pT_rho->Fill(thisTrackSgn->GetPt());
          Eta_rho->Fill(thisTrackSgn->GetEta());
          pT_px_rho->Fill(thisTrackSgn->GetStartVertexMomentumX(),thisTrackSgn->GetPt());
          Eta_py_rho->Fill(thisTrackSgn->GetStartVertexMomentumY(),thisTrackSgn->GetEta());
        }
        if (thisTrackSgn->GetPdgCode()==223) {
          pT_omega->Fill(thisTrackSgn->GetPt());
          Eta_omega->Fill(thisTrackSgn->GetEta());
          pT_px_omega->Fill(thisTrackSgn->GetStartVertexMomentumX(),thisTrackSgn->GetPt());
          Eta_py_omega->Fill(thisTrackSgn->GetStartVertexMomentumY(),thisTrackSgn->GetEta());
        }
        if (thisTrackSgn->GetPdgCode()==333) {
          pT_phi->Fill(thisTrackSgn->GetPt());
          Eta_phi->Fill(thisTrackSgn->GetEta());
          pT_px_phi->Fill(thisTrackSgn->GetStartVertexMomentumX(),thisTrackSgn->GetPt());
          Eta_py_phi->Fill(thisTrackSgn->GetStartVertexMomentumY(),thisTrackSgn->GetEta());
        }
        if (thisTrackSgn->GetPdgCode()==443) {
          pT_Jpsi->Fill(thisTrackSgn->GetPt());
          Eta_Jpsi->Fill(thisTrackSgn->GetEta());
        }
        if (thisTrackSgn->GetPdgCode()==20133) {
          pT_a1->Fill(thisTrackSgn->GetPt());
          Eta_a1->Fill(thisTrackSgn->GetEta());
          pT_px_a1->Fill(thisTrackSgn->GetStartVertexMomentumX(),thisTrackSgn->GetPt());
          Eta_py_a1->Fill(thisTrackSgn->GetStartVertexMomentumY(),thisTrackSgn->GetEta());
        }
      }
    }
    pT_rho->Write();
    Eta_rho->Write();
    pT_omega->Write();
    Eta_omega->Write();
    pT_phi->Write();
    Eta_phi->Write();
    pT_Jpsi->Write();
    Eta_Jpsi->Write();
    pT_a1->Write();
    Eta_a1->Write();

    pT_px_rho->Write();
    Eta_py_rho->Write();
    pT_px_omega->Write();
    Eta_py_omega->Write();
    pT_px_phi->Write();
    Eta_py_phi->Write();
    pT_px_a1->Write();
    Eta_py_a1->Write();
    outFile.Close();
  }
