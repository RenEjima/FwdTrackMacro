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

void LowMassVectorMesonPtEta_DPG(){
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

  TH1F *pdg_list = new TH1F("pdg_list","PDG code distribution;PDG code;Entry",40000000,-10000000,10000000);
  TH1F *pdg_muon_mother = new TH1F("pdg_muon_mother","PDG code of muon's mother;PDG code;Entry",40000000,-10000000,10000000);
  TH1F *pdg_muon_mother_2 = new TH1F("pdg_muon_mother_2","PDG code of muon's second mother;PDG code;Entry",40000000,-10000000,10000000);
  TH1F *muFromA1Pt = new TH1F("muFromA1Pt","p_{T} of #mu from a_{1} meson;p_{T}[GeV/c];Entry",10000,0,10);
  TH2F *muFromA1PtPl = new TH2F("muFromA1PtPl","p_{L} - p_{T} correlation of #mu from a_{1} meson;p_{L}[GeV/c];p_{T}[GeV/c]",10000,0,10,10000,0,10);
  TH1F *a1_daughter_pdg = new TH1F("a1_daughter_pdg","PDG code of a1's daughter;PDG code;Entry",40000000,-10000000,10000000);

  TH2F *MuonThroughMFT = new TH2F("MuonThroughMFT","#mu in the acceptance of MFT",1000,0,10,1000,-4.5,-1.5);
  TH2F *RhoThroughMFT = new TH2F("RhoThroughMFT","#rho in the acceptance of MFT",1000,0,10,1000,-5.5,5.5);

  vector<MCTrackT<float>> *mcTr = nullptr;
  o2SimKineTree->SetBranchAddress("MCTrack", &mcTr);

  vector<MCTrackT<float>> *mcTrSgn = nullptr;
  sgnKineTree->SetBranchAddress("MCTrack", &mcTrSgn);

  Int_t numberOfEvents = o2SimKineTree->GetEntries();

  for (int iEvent = 0; iEvent < numberOfEvents; iEvent++) {

    o2SimKineTree->GetEntry(iEvent);
    int nMCTracks = (*mcTr).size();

    sgnKineTree->GetEntry(iEvent);
    int nMCTracksSgn = (*mcTrSgn).size();
    std::cout<<"Start Event #"<<iEvent<<" : nMCTracks_bkg = "<<nMCTracks<<" : nMCTracks_sgn = "<<nMCTracksSgn<<std::endl;

    for (Int_t iMCTrack=0; iMCTrack<(*mcTr).size();++iMCTrack) {
        MCTrackT<float> *thisTrack = &(*mcTr).at(iMCTrack);
        std::vector<Int_t>  TrackList;

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
        if (thisTrack->GetPdgCode()==20213) {
          pT_a1->Fill(thisTrack->GetPt());
          Eta_a1->Fill(thisTrack->GetEta());
          pT_px_a1->Fill(thisTrack->GetStartVertexMomentumX(),thisTrack->GetPt());
          Eta_py_a1->Fill(thisTrack->GetStartVertexMomentumY(),thisTrack->GetEta());
          Int_t a1_daughter_first_ID = thisTrack->getFirstDaughterTrackId();
          Int_t a1_daughter_last_ID = thisTrack->getLastDaughterTrackId();
          for (Int_t daughterID = a1_daughter_first_ID; daughterID< a1_daughter_last_ID;daughterID++){
            MCTrackT<float> *thisTrack_daughter = &(*mcTr).at(daughterID);
            a1_daughter_pdg->Fill(thisTrack_daughter->GetPdgCode());
            //std::cout<<"a1 #"<<iMCTrack<<" : daughter PDG = "<<thisTrack_daughter->GetPdgCode()<<std::endl;
          }
        }
        if (thisTrack->GetPdgCode()==-20213) {
          pT_a1->Fill(thisTrack->GetPt());
          Eta_a1->Fill(thisTrack->GetEta());
          pT_px_a1->Fill(thisTrack->GetStartVertexMomentumX(),thisTrack->GetPt());
          Eta_py_a1->Fill(thisTrack->GetStartVertexMomentumY(),thisTrack->GetEta());
          Int_t a1_daughter_first_ID = thisTrack->getFirstDaughterTrackId();
          Int_t a1_daughter_last_ID = thisTrack->getLastDaughterTrackId();
          for (Int_t daughterID = a1_daughter_first_ID; daughterID< a1_daughter_last_ID;daughterID++){
            MCTrackT<float> *thisTrack_daughter = &(*mcTr).at(daughterID);
            a1_daughter_pdg->Fill(thisTrack_daughter->GetPdgCode());
            //std::cout<<"a1 #"<<iMCTrack<<" : daughter PDG = "<<thisTrack_daughter->GetPdgCode()<<std::endl;
          }
        }

        if (thisTrack->GetPdgCode()==13||thisTrack->GetPdgCode()==-13) {
          Int_t motherID = thisTrack->getMotherTrackId();
          Int_t motherID_2 = thisTrack->getSecondMotherTrackId();
          if (thisTrack->GetPt()>0.5 && thisTrack->GetEta()>-3.6 && thisTrack->GetEta()<-2.5 && thisTrack->GetStartVertexCoordinatesZ()>-77.5){
            MuonThroughMFT->Fill(thisTrack->GetPt(),thisTrack->GetEta());
          }
          if (motherID != 0 && motherID !=-1) {
            //std::cout<<"iMCTrack = "<<iMCTrack<<" : motherID = "<<motherID<<" : motherID_2 = "<<motherID_2<<std::endl;
            MCTrackT<float> *thisTrack_mother = &(*mcTr).at(motherID);
            //std::cout<<"got mother track."<<std::endl;
            pdg_muon_mother->Fill(thisTrack_mother->GetPdgCode());
            if (thisTrack_mother->GetPdgCode() == 113 || thisTrack_mother->GetPdgCode() == -113) { //if thistrack's mother is rho,
              Int_t daughterOfRhoIDF = thisTrack_mother->getFirstDaughterTrackId();
              Int_t daughterOfRhoIDL = thisTrack_mother->getLastDaughterTrackId();
              auto search_result_F = std::find(TrackList.begin(),TrackList.end(),daughterOfRhoIDF);
              auto search_result_L = std::find(std::cbegin(TrackList),std::cend(TrackList),daughterOfRhoIDL);
              MCTrackT<float> *thisTrackF = &(*mcTr).at(daughterOfRhoIDF);
              MCTrackT<float> *thisTrackL = &(*mcTr).at(daughterOfRhoIDL);
              if ((thisTrackF->GetPt()>0.5)&&(thisTrackF->GetEta()>-3.6)&&(thisTrackF->GetEta()<-2.5)&&(thisTrackF->GetStartVertexCoordinatesZ()>-77.5)&&(thisTrackL->GetPt()>0.5)&&(thisTrackL->GetEta()>-3.6)&&(thisTrackL->GetEta()<-2.5)&&(thisTrackL->GetStartVertexCoordinatesZ()>-77.5)&&(search_result_F == std::cend(TrackList))&&(search_result_L == std::cend(TrackList))){
                RhoThroughMFT->Fill(thisTrack_mother->GetPt(),thisTrack_mother->GetEta());
                TrackList.push_back(daughterOfRhoIDF);
                TrackList.push_back(daughterOfRhoIDL);
                std::cout<<"IDF = "<<daughterOfRhoIDF<<" : IDL = "<<daughterOfRhoIDL<<std::endl;
                copy(TrackList.begin(), TrackList.end(), ostream_iterator<int>(cout, "; "));
                cout << endl;
              }
              Int_t motherOfRhoID = thisTrack_mother->getMotherTrackId();//get rho's mother ID
              if (motherOfRhoID<=(*mcTr).size()) {
                MCTrackT<float> *thisTrack_motherOfRho = &(*mcTr).at(motherOfRhoID); //get rho's mother track
                if (thisTrack_motherOfRho->GetPdgCode() == 20213 || thisTrack_motherOfRho->GetPdgCode() == -20213){ //if rho's mother is a1,
                  //std::cout<<"mu from a1!!"<<std::endl;
                  muFromA1Pt->Fill(thisTrack->GetPt());
                  muFromA1PtPl->Fill(thisTrack->GetStartVertexMomentumZ(),thisTrack->GetPt());
                }
              }
            }
          }
          if (motherID_2 != 0 && motherID_2 !=-1) {
            //std::cout<<"iMCTrack = "<<iMCTrack<<" : motherID = "<<motherID<<" : motherID_2 = "<<motherID_2<<std::endl;
            MCTrackT<float> *thisTrack_mother_2 = &(*mcTr).at(motherID_2);
            //std::cout<<"got mother track."<<std::endl;
            pdg_muon_mother->Fill(thisTrack_mother_2->GetPdgCode());
            if (thisTrack_mother_2->GetPdgCode() == 113 || thisTrack_mother_2->GetPdgCode() == -113) { //if thistrack's mother is rho,
              Int_t daughterOfRhoIDF = thisTrack_mother_2->getFirstDaughterTrackId();
              Int_t daughterOfRhoIDL = thisTrack_mother_2->getLastDaughterTrackId();
              auto search_result_F = std::find(TrackList.begin(),TrackList.end(),daughterOfRhoIDF);
              auto search_result_L = std::find(std::cbegin(TrackList),std::cend(TrackList),daughterOfRhoIDL);
              MCTrackT<float> *thisTrackF = &(*mcTr).at(daughterOfRhoIDF);
              MCTrackT<float> *thisTrackL = &(*mcTr).at(daughterOfRhoIDL);
              if ((thisTrackF->GetPt()>0.5)&&(thisTrackF->GetEta()>-3.6)&&(thisTrackF->GetEta()<-2.5)&&(thisTrackF->GetStartVertexCoordinatesZ()>-77.5)&&(thisTrackL->GetPt()>0.5)&&(thisTrackL->GetEta()>-3.6)&&(thisTrackL->GetEta()<-2.5)&&(thisTrackL->GetStartVertexCoordinatesZ()>-77.5)&&(search_result_F == std::cend(TrackList))&&(search_result_L == std::cend(TrackList))){
                RhoThroughMFT->Fill(thisTrack_mother_2->GetPt(),thisTrack_mother_2->GetEta());
                TrackList.push_back(daughterOfRhoIDF);
                TrackList.push_back(daughterOfRhoIDL);
              }
              Int_t motherOfRhoID_2 = thisTrack_mother_2->getMotherTrackId();//get rho's mother ID
              if (motherOfRhoID_2<=(*mcTr).size()){
                MCTrackT<float> *thisTrack_motherOfRho_2 = &(*mcTr).at(motherOfRhoID_2); //get rho's mother track
                if (thisTrack_motherOfRho_2->GetPdgCode() == 20213 || thisTrack_motherOfRho_2->GetPdgCode() == -20213){ //if rho's mother is a1,
                  //std::cout<<"mu from a1!!"<<std::endl;
                  muFromA1Pt->Fill(thisTrack->GetPt());
                  muFromA1PtPl->Fill(thisTrack->GetStartVertexMomentumZ(),thisTrack->GetPt());
                }
              }
            }
          }
        }

        pdg_list->Fill(thisTrack->GetPdgCode());
    }
    for (Int_t iMCTrack=0; iMCTrack<(*mcTrSgn).size();++iMCTrack) {
        MCTrackT<float> *thisTrack = &(*mcTrSgn).at(iMCTrack);
        std::vector<Int_t>  TrackList;

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
        if (thisTrack->GetPdgCode()==20213) {
          pT_a1->Fill(thisTrack->GetPt());
          Eta_a1->Fill(thisTrack->GetEta());
          pT_px_a1->Fill(thisTrack->GetStartVertexMomentumX(),thisTrack->GetPt());
          Eta_py_a1->Fill(thisTrack->GetStartVertexMomentumY(),thisTrack->GetEta());
          Int_t a1_daughter_first_ID = thisTrack->getFirstDaughterTrackId();
          Int_t a1_daughter_last_ID = thisTrack->getLastDaughterTrackId();
          for (Int_t daughterID = a1_daughter_first_ID; daughterID< a1_daughter_last_ID;daughterID++){
            MCTrackT<float> *thisTrack_daughter = &(*mcTrSgn).at(daughterID);
            a1_daughter_pdg->Fill(thisTrack_daughter->GetPdgCode());
            //std::cout<<"a1 #"<<iMCTrack<<" : daughter PDG = "<<thisTrack_daughter->GetPdgCode()<<std::endl;
          }
        }
        if (thisTrack->GetPdgCode()==-20213) {
          pT_a1->Fill(thisTrack->GetPt());
          Eta_a1->Fill(thisTrack->GetEta());
          pT_px_a1->Fill(thisTrack->GetStartVertexMomentumX(),thisTrack->GetPt());
          Eta_py_a1->Fill(thisTrack->GetStartVertexMomentumY(),thisTrack->GetEta());
          Int_t a1_daughter_first_ID = thisTrack->getFirstDaughterTrackId();
          Int_t a1_daughter_last_ID = thisTrack->getLastDaughterTrackId();
          for (Int_t daughterID = a1_daughter_first_ID; daughterID< a1_daughter_last_ID;daughterID++){
            MCTrackT<float> *thisTrack_daughter = &(*mcTrSgn).at(daughterID);
            a1_daughter_pdg->Fill(thisTrack_daughter->GetPdgCode());
            //std::cout<<"a1 #"<<iMCTrack<<" : daughter PDG = "<<thisTrack_daughter->GetPdgCode()<<std::endl;
          }
        }

        if (thisTrack->GetPdgCode()==13||thisTrack->GetPdgCode()==-13) {
          Int_t motherID = thisTrack->getMotherTrackId();
          Int_t motherID_2 = thisTrack->getSecondMotherTrackId();
          if (thisTrack->GetPt()>0.5 && thisTrack->GetEta()>-3.6 && thisTrack->GetEta()<-2.5 && thisTrack->GetStartVertexCoordinatesZ()>-77.5){
            MuonThroughMFT->Fill(thisTrack->GetPt(),thisTrack->GetEta());
          }
          if (motherID != 0 && motherID !=-1) {
            //std::cout<<"iMCTrack = "<<iMCTrack<<" : motherID = "<<motherID<<" : motherID_2 = "<<motherID_2<<std::endl;
            MCTrackT<float> *thisTrack_mother = &(*mcTrSgn).at(motherID);
            //std::cout<<"got mother track."<<std::endl;
            pdg_muon_mother->Fill(thisTrack_mother->GetPdgCode());
            if (thisTrack_mother->GetPdgCode() == 113 || thisTrack_mother->GetPdgCode() == -113) { //if thistrack's mother is rho,
              Int_t daughterOfRhoIDF = thisTrack_mother->getFirstDaughterTrackId();
              Int_t daughterOfRhoIDL = thisTrack_mother->getLastDaughterTrackId();
              auto search_result_F = std::find(TrackList.begin(),TrackList.end(),daughterOfRhoIDF);
              auto search_result_L = std::find(std::cbegin(TrackList),std::cend(TrackList),daughterOfRhoIDL);
              MCTrackT<float> *thisTrackF = &(*mcTrSgn).at(daughterOfRhoIDF);
              MCTrackT<float> *thisTrackL = &(*mcTrSgn).at(daughterOfRhoIDL);
              if ((thisTrackF->GetPt()>0.5)&&(thisTrackF->GetEta()>-3.6)&&(thisTrackF->GetEta()<-2.5)&&(thisTrackF->GetStartVertexCoordinatesZ()>-77.5)&&(thisTrackL->GetPt()>0.5)&&(thisTrackL->GetEta()>-3.6)&&(thisTrackL->GetEta()<-2.5)&&(thisTrackL->GetStartVertexCoordinatesZ()>-77.5)&&(search_result_F == std::cend(TrackList))&&(search_result_L == std::cend(TrackList))){
                RhoThroughMFT->Fill(thisTrack_mother->GetPt(),thisTrack_mother->GetEta());
                TrackList.push_back(daughterOfRhoIDF);
                TrackList.push_back(daughterOfRhoIDL);
                std::cout<<"IDF = "<<daughterOfRhoIDF<<" : IDL = "<<daughterOfRhoIDL<<std::endl;
                copy(TrackList.begin(), TrackList.end(), ostream_iterator<int>(cout, "; "));
                cout << endl;
              }
              Int_t motherOfRhoID = thisTrack_mother->getMotherTrackId();//get rho's mother ID
              if (motherOfRhoID<=(*mcTrSgn).size()) {
                MCTrackT<float> *thisTrack_motherOfRho = &(*mcTrSgn).at(motherOfRhoID); //get rho's mother track
                if (thisTrack_motherOfRho->GetPdgCode() == 20213 || thisTrack_motherOfRho->GetPdgCode() == -20213){ //if rho's mother is a1,
                  //std::cout<<"mu from a1!!"<<std::endl;
                  muFromA1Pt->Fill(thisTrack->GetPt());
                  muFromA1PtPl->Fill(thisTrack->GetStartVertexMomentumZ(),thisTrack->GetPt());
                }
              }
            }
          }
          if (motherID_2 != 0 && motherID_2 !=-1) {
            //std::cout<<"iMCTrack = "<<iMCTrack<<" : motherID = "<<motherID<<" : motherID_2 = "<<motherID_2<<std::endl;
            MCTrackT<float> *thisTrack_mother_2 = &(*mcTrSgn).at(motherID_2);
            //std::cout<<"got mother track."<<std::endl;
            pdg_muon_mother->Fill(thisTrack_mother_2->GetPdgCode());
            if (thisTrack_mother_2->GetPdgCode() == 113 || thisTrack_mother_2->GetPdgCode() == -113) { //if thistrack's mother is rho,
              Int_t daughterOfRhoIDF = thisTrack_mother_2->getFirstDaughterTrackId();
              Int_t daughterOfRhoIDL = thisTrack_mother_2->getLastDaughterTrackId();
              auto search_result_F = std::find(TrackList.begin(),TrackList.end(),daughterOfRhoIDF);
              auto search_result_L = std::find(std::cbegin(TrackList),std::cend(TrackList),daughterOfRhoIDL);
              MCTrackT<float> *thisTrackF = &(*mcTrSgn).at(daughterOfRhoIDF);
              MCTrackT<float> *thisTrackL = &(*mcTrSgn).at(daughterOfRhoIDL);
              if ((thisTrackF->GetPt()>0.5)&&(thisTrackF->GetEta()>-3.6)&&(thisTrackF->GetEta()<-2.5)&&(thisTrackF->GetStartVertexCoordinatesZ()>-77.5)&&(thisTrackL->GetPt()>0.5)&&(thisTrackL->GetEta()>-3.6)&&(thisTrackL->GetEta()<-2.5)&&(thisTrackL->GetStartVertexCoordinatesZ()>-77.5)&&(search_result_F == std::cend(TrackList))&&(search_result_L == std::cend(TrackList))){
                RhoThroughMFT->Fill(thisTrack_mother_2->GetPt(),thisTrack_mother_2->GetEta());
                TrackList.push_back(daughterOfRhoIDF);
                TrackList.push_back(daughterOfRhoIDL);
              }
              Int_t motherOfRhoID_2 = thisTrack_mother_2->getMotherTrackId();//get rho's mother ID
              if (motherOfRhoID_2<=(*mcTrSgn).size()) {
                MCTrackT<float> *thisTrack_motherOfRho_2 = &(*mcTrSgn).at(motherOfRhoID_2); //get rho's mother track
                if (thisTrack_motherOfRho_2->GetPdgCode() == 20213 || thisTrack_motherOfRho_2->GetPdgCode() == -20213){ //if rho's mother is a1,
                  //std::cout<<"mu from a1!!"<<std::endl;
                  muFromA1Pt->Fill(thisTrack->GetPt());
                  muFromA1PtPl->Fill(thisTrack->GetStartVertexMomentumZ(),thisTrack->GetPt());
                }
              }
            }
          }
        }

        pdg_list->Fill(thisTrack->GetPdgCode());
    }
    /*
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

        if (thisTrackSgn->GetPdgCode()==20233) {
          pT_a1->Fill(thisTrackSgn->GetPt());
          Eta_a1->Fill(thisTrackSgn->GetEta());
          pT_px_a1->Fill(thisTrackSgn->GetStartVertexMomentumX(),thisTrackSgn->GetPt());
          Eta_py_a1->Fill(thisTrackSgn->GetStartVertexMomentumY(),thisTrackSgn->GetEta());
        }

        pdg_list->Fill(thisTrackSgn->GetPdgCode());
      }
      */
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
    pdg_list->Write();
    pdg_muon_mother->Write();
    pdg_muon_mother_2->Write();
    muFromA1Pt->Write();
    muFromA1PtPl->Write();
    a1_daughter_pdg->Write();

    pT_px_rho->Write();
    Eta_py_rho->Write();
    pT_px_omega->Write();
    Eta_py_omega->Write();
    pT_px_phi->Write();
    Eta_py_phi->Write();
    pT_px_a1->Write();
    Eta_py_a1->Write();

    MuonThroughMFT->Write();
    RhoThroughMFT->Write();
    outFile.Close();
  }
