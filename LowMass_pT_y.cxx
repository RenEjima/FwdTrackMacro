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

bool DEBUG_VERBOSE = false;
bool EXPORT_HISTOS_IMAGES = false;

double getZField(double x, double y, double z)
{
  const auto grp = o2::parameters::GRPObject::loadFrom("o2sim_grp.root");
  std::unique_ptr<o2::parameters::GRPObject> mGRP = nullptr;
  mGRP.reset(grp);
  o2::base::Propagator::initFieldFromGRP(grp);
  auto field = static_cast<o2::field::MagneticField *>(TGeoGlobalMagField::Instance()->GetField());

  double position[3] = {x, y, z};
  auto Bz = field->getBz(position);
  printf("B field z = %f [kGauss]\n", Bz);
  return Bz;
}

//_________________________________________________________________________________________________
void InvariantMass(float chi2Cut = 5.f,
                      const std::string trkFile = "globalfwdtracks.root",
                      const std::string o2sim_KineFile = "bkg_Kine.root",
                      const std::string sig_KineFile = "sgn_1_Kine.root")
{
  // Files & Trees
  // MC
  TFile *o2sim_KineFileIn = new TFile(o2sim_KineFile.c_str());
  TTree *o2SimKineTree = (TTree *)o2sim_KineFileIn->Get("o2sim");
  TFile *sig_KineFileIn = new TFile(sig_KineFile.c_str());

  TTree *sigKineTree = (TTree *)sig_KineFileIn->Get("o2sim");

  vector<MCTrackT<float>> *mcTr = nullptr;
  o2SimKineTree->SetBranchAddress("MCTrack", &mcTr);
  o2::dataformats::MCEventHeader *eventHeader = nullptr;
  o2SimKineTree->SetBranchAddress("MCEventHeader.", &eventHeader);

  vector<MCTrackT<float>> *mcTrSig = nullptr;
  sigKineTree->SetBranchAddress("MCTrack", &mcTrSig);

  Int_t numberOfEvents = o2SimKineTree->GetEntries();

  // Global Muon Tracks
  TFile *trkFileIn = new TFile(trkFile.c_str());
  TTree *gmTrackTree = (TTree *)trkFileIn->Get("GlobalFwdTracks");
  std::vector<GlobalMuonTrack> trackGMVec, *trackGMVecP = &trackGMVec;
  gmTrackTree->SetBranchAddress("fwdtracks", &trackGMVecP);

  vector<o2::MCCompLabel> *mcLabels = nullptr;
  gmTrackTree->SetBranchAddress("MCTruth", &mcLabels);

  // MFT Tracks
  TFile *mfttrkFileIn = new TFile("mfttracks.root");
  TTree *mftTrackTree = (TTree *)mfttrkFileIn->Get("o2sim");
  std::vector<o2::mft::TrackMFT> trackMFTVec, *trackMFTVecP = &trackMFTVec;
  mftTrackTree->SetBranchAddress("MFTTrack", &trackMFTVecP);

  std::vector<o2::MCCompLabel> *mftMcLabels = nullptr;
  mftTrackTree->SetBranchAddress("MFTTrackMCTruth", &mftMcLabels);
  mftTrackTree->GetEntry(0);

  gmTrackTree->GetEntry(0);
  o2SimKineTree->GetEntry(0);
  sigKineTree->GetEntry(0);

  auto field_z = getZField(0, 0, -61.4); // Get field at Center of MFT

  std::string outfilename = "Mass.root";
  TFile outFile(outfilename.c_str(), "RECREATE");

  // Reconstructed Global Muon Tracks
  std::cout << "Loop over Global Muon Tracks!" << std::endl;
  auto iTrack = 0;
  Double_t m_mu = 0.106;
  Int_t Nbin = 1000;
  TH1F *SamePP = new TH1F("SamePP","Invariant Mass Spectrum (N_{++}^{Same});m^{#mu#mu}_{GM}[GeV/c^{2}]",Nbin,0,10);
  TH1F *SameMM = new TH1F("SameMM","Invariant Mass Spectrum (N_{--}^{Same});m^{#mu#mu}_{GM}[GeV/c^{2}]",Nbin,0,10);
  TH1F *SamePM = new TH1F("SamePM","Invariant Mass Spectrum (N_{+-}^{Same});m^{#mu#mu}_{GM}[GeV/c^{2}]",Nbin,0,10);
  TH1F *MixedPP = new TH1F("MixedPP","Invariant Mass Spectrum (N_{++}^{Mixed});m^{#mu#mu}_{GM}[GeV/c^{2}]",Nbin,0,10);
  TH1F *MixedMM = new TH1F("MixedMM","Invariant Mass Spectrum (N_{--}^{Mixed});m^{#mu#mu}_{GM}[GeV/c^{2}]",Nbin,0,10);
  TH1F *MixedPM = new TH1F("MixedPM","Invariant Mass Spectrum (N_{+-}^{Mixed});m^{#mu#mu}_{GM}[GeV/c^{2}]",Nbin,0,10);
  TH1F *LikeSign = new TH1F("LikeSign","Invariant Mass Spectrum (LikeSign-method);m^{#mu#mu}_{GM}[GeV/c^{2}]",Nbin,0,10);
  //TH1F *InvMassSpectrum_MC = new TH1F("InvMassSpectrum_MC","MC Invariant Mass Spectrum;m^{#mu#mu}_{GM}[GeV/c^{2}]",Nbin,0,10);
  //TH1F *PDGdist = new TH1F("PDGdist","PDGcode distribution;PDG code",400000,-200000,200000);
  for (auto &gmTrack : trackGMVec)
  {
    if (gmTrack.getMIDMatchingChi2() < 0 || gmTrack.getMFTMCHMatchingChi2() > chi2Cut) // Filter Muon Tracks and MFTMCH matching Chi2
    {
      iTrack++;
      continue;
    }
    auto bestMFTTrackMatchID = gmTrack.getMFTTrackID();
    auto &mftTrackMatch = trackMFTVec[bestMFTTrackMatchID];
    MCTrackT<float> *thisTrack;

    const auto &label = mcLabels->at(iTrack);

    auto thisTrkID = label.getTrackID();
    auto thisEvtID = label.getEventID();
    if (label.getSourceID() == 0)
    {
      o2SimKineTree->GetEntry(thisEvtID);
      thisTrack = &(mcTr->at(thisTrkID));
    }
    else if (label.getSourceID() == 1)
    {
      sigKineTree->GetEntry(thisEvtID);
      thisTrack = &(mcTrSig->at(thisTrkID));
    }
    else
    {
      iTrack++;
      continue;
    }

    auto iTrack2 = 0;
      for (auto &gmTrack2 : trackGMVec) {
        if (gmTrack2.getMIDMatchingChi2() < 0 || gmTrack2.getMFTMCHMatchingChi2() > chi2Cut) // Filter Muon Tracks and MFTMCH matching Chi2
        {
          iTrack2++;
          continue;
        }
        auto bestMFTTrackMatchID2 = gmTrack2.getMFTTrackID();
        auto &mftTrackMatch2 = trackMFTVec[bestMFTTrackMatchID2];
        MCTrackT<float> *thisTrack2;

        const auto &label2 = mcLabels->at(iTrack2);

        auto thisTrkID2 = label2.getTrackID();
        auto thisEvtID2 = label2.getEventID();
        if (label2.getSourceID() == 0)
        {
          o2SimKineTree->GetEntry(thisEvtID2);
          thisTrack2 = &(mcTr->at(thisTrkID2));
        }
        else if (label2.getSourceID() == 1)
        {
          sigKineTree->GetEntry(thisEvtID2);
          thisTrack2 = &(mcTrSig->at(thisTrkID2));
        }
        else
        {
          iTrack2++;
          continue;
        }
        if (label.getTrackID() >= label2.getTrackID()) { // to avoid duplication
          iTrack2++;
          continue;
        }

        // Invariant Mass
        gmTrack.propagateToZhelix(0., field_z);
        gmTrack2.propagateToZhelix(0., field_z);
        Double_t px_1 = gmTrack.getPx();
        Double_t py_1 = gmTrack.getPy();
        Double_t pz_1 = gmTrack.getPz();
        Double_t px_2 = gmTrack2.getPx();
        Double_t py_2 = gmTrack2.getPy();
        Double_t pz_2 = gmTrack2.getPz();
        Double_t p_1 = gmTrack.getP();
        Double_t p_2 = gmTrack2.getP();
        Double_t CosTheta = (px_1*px_2+py_1*py_2+pz_1*pz_2)/(p_1*p_2);
        Double_t InvMass = sqrt(2.*m_mu*m_mu+2.*(sqrt(m_mu*m_mu+p_1*p_1)*sqrt(m_mu*m_mu+p_2*p_2)-p_1*p_2*CosTheta));

        // Mixed event
        if (label.getEventID() != label2.getEventID()) {
          // charge
          if (gmTrack.getCharge() == gmTrack2.getCharge()) {
            if (gmTrack.getCharge() == 1.) { // charge = ++
              MixedPP->Fill(InvMass);
            }
            if (gmTrack.getCharge() == -1.) { // charge = --
              MixedMM->Fill(InvMass);
            }
            iTrack2++;
            continue;
          }
          MixedPM->Fill(InvMass); // charge = +-
          iTrack2++;
          continue;
        }

        // Same event
        // reconstruct mass from MC tracks assosiated with globalfwdtracks
        /*if (thisTrack->getMotherTrackId() == thisTrack2->getMotherTrackId()) { // from the same mother
          Double_t px_1_MC = thisTrack->GetStartVertexMomentumX();
          Double_t py_1_MC = thisTrack->GetStartVertexMomentumY();
          Double_t pz_1_MC = thisTrack->GetStartVertexMomentumZ();
          Double_t px_2_MC = thisTrack2->GetStartVertexMomentumX();
          Double_t py_2_MC = thisTrack2->GetStartVertexMomentumY();
          Double_t pz_2_MC = thisTrack2->GetStartVertexMomentumZ();
          Double_t p_1_MC = sqrt(px_1_MC*px_1_MC+py_1_MC*py_1_MC+pz_1_MC*pz_1_MC);
          Double_t p_2_MC = sqrt(px_2_MC*px_2_MC+py_2_MC*py_2_MC+pz_2_MC*pz_2_MC);
          Double_t CosTheta_MC = (px_1_MC*px_2_MC+py_1_MC*py_2_MC+pz_1_MC*pz_2_MC)/(p_1_MC*p_2_MC);
          Double_t InvMass_MC = sqrt(2.*m_mu*m_mu+2.*(sqrt(m_mu*m_mu+p_1_MC*p_1_MC)*sqrt(m_mu*m_mu+p_2_MC*p_2_MC)-p_1_MC*p_2_MC*CosTheta_MC));
          MCTrackT<float> *thisTrack_Mother;
          thisTrack_Mother = &(mcTr->at(thisTrack->getMotherTrackId()));
          Double_t PDGcode = thisTrack_Mother->GetPdgCode();
          InvMassSpectrum_MC->Fill(InvMass_MC);
          PDGdist->Fill(PDGcode);
        }*/
        // charge
        if (gmTrack.getCharge() == gmTrack2.getCharge()) {
          if (gmTrack.getCharge() == 1.) { // charge = ++
            SamePP->Fill(InvMass);
          }
          if (gmTrack.getCharge() == -1.) { // charge = --
            SameMM->Fill(InvMass);
          }
          iTrack2++;
          continue;
        }
        SamePM->Fill(InvMass); // charge = +-

        //std::cout<<"1 iTrack = "<<iTrack<<" : Event ID = "<<label.getEventID()<<" : MCHtrack ID = "<<label.getTrackID()<<" : MFTtrack ID = "<<gmTrack.getMFTTrackID()<<endl;
        //std::cout<<"2 iTrack = "<<iTrack2<<" : Event ID = "<<label2.getEventID()<<" : MCHtrack ID = "<<label2.getTrackID()<<" : MFTtrack ID = "<<gmTrack2.getMFTTrackID()<<endl;
        //std::cout<<"CosTheta = "<<CosTheta<<" : Invariant Mass = "<<InvMass<<"GeV/c^2"<<endl;
        //std::cout<<"-------------------------------------------------------------------------"<<endl;
        iTrack2++;
      } // Loop on GMtracks2

    iTrack++;

  } // Loop on GMTracks

  // LikeSign-method
  for (Int_t bin = 0; bin < Nbin; bin++) {
    Double_t NsamePP = SamePP->GetBinContent(bin);
    Double_t NsameMM = SameMM->GetBinContent(bin);
    Double_t NsamePM = SamePM->GetBinContent(bin);
    Double_t NmixedPP = MixedPP->GetBinContent(bin);
    Double_t NmixedMM = MixedMM->GetBinContent(bin);
    Double_t NmixedPM = MixedPM->GetBinContent(bin);
    if (NmixedPP == 0. || NmixedMM == 0.) {
      bin++;
      continue;
    }
    Double_t Rfactor = NmixedPM/(2.*sqrt(NmixedPP*NmixedMM));
    Double_t Nsignal = NsamePM-2.*Rfactor*sqrt(NsamePP*NsameMM);
    //std::cout<<"bin = "<<bin<<endl;
    //std::cout<<"Same PP = "<<NsamePP<<" : MM = "<<NsameMM<<" : PM = "<<NsamePM<<" : Mixed PP = "<<NmixedPP<<" : MM = "<<NmixedMM<<" : PM = "<<NmixedPM<<endl;
    //std::cout<<"Rfactor = "<<Rfactor<<" : Nsignal = "<<Nsignal<<endl;
    //std::cout<<"-------------------------------------------------------------------------"<<endl;
    LikeSign->SetBinContent(bin,Nsignal);
  }

  SamePM->Write();
  SamePP->Write();
  SameMM->Write();
  MixedPM->Write();
  MixedPP->Write();
  MixedMM->Write();
  LikeSign->Write();
  //InvMassSpectrum_MC->Write();
  //PDGdist->Write();
  outFile.Close();

  return;
}
