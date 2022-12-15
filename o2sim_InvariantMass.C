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
void o2sim_InvariantMass(float chi2Cut = 5.f,
                      const std::string trkFile = "globalfwdtracks.root",
                      const std::string o2sim_KineFile = "o2sim_Kine.root")
{
  //===================================================================================================================================
  //  Read Files and create TTree
  //===================================================================================================================================

  // MC
      TFile *o2sim_KineFileIn = new TFile(o2sim_KineFile.c_str());
      TTree *o2SimKineTree = (TTree *)o2sim_KineFileIn->Get("o2sim");

      vector<MCTrackT<float>> *mcTr = nullptr;
      o2SimKineTree->SetBranchAddress("MCTrack", &mcTr);

      o2::dataformats::MCEventHeader *eventHeader = nullptr;
      o2SimKineTree->SetBranchAddress("MCEventHeader.", &eventHeader);

      Int_t numberOfEvents = o2SimKineTree->GetEntries();

      std::cout<<"number of event = "<<numberOfEvents<<std::endl;

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

  // Set all of TTree's entry
      mftTrackTree->GetEntry(0);
      gmTrackTree->GetEntry(0);
      o2SimKineTree->GetEntry(0);

  // Get Magnetic Field at the cetner of MFT
      auto field_z = getZField(0, 0, -61.4);

  // Output file
      std::string outfilename = "Mass.root";
      TFile outFile(outfilename.c_str(), "RECREATE");

  // Histograms
      Int_t Nbin = 1000;

      TH1F *SamePP = new TH1F("SamePP","Invariant Mass Spectrum (N_{++}^{Same});m^{#mu#mu}_{GM}[GeV/c^{2}]",Nbin,0,10);
      TH1F *SameMM = new TH1F("SameMM","Invariant Mass Spectrum (N_{--}^{Same});m^{#mu#mu}_{GM}[GeV/c^{2}]",Nbin,0,10);
      TH1F *SamePM = new TH1F("SamePM","Invariant Mass Spectrum (N_{+-}^{Same});m^{#mu#mu}_{GM}[GeV/c^{2}]",Nbin,0,10);
      TH1F *MixedPP = new TH1F("MixedPP","Invariant Mass Spectrum (N_{++}^{Mixed});m^{#mu#mu}_{GM}[GeV/c^{2}]",Nbin,0,10);
      TH1F *MixedMM = new TH1F("MixedMM","Invariant Mass Spectrum (N_{--}^{Mixed});m^{#mu#mu}_{GM}[GeV/c^{2}]",Nbin,0,10);
      TH1F *MixedPM = new TH1F("MixedPM","Invariant Mass Spectrum (N_{+-}^{Mixed});m^{#mu#mu}_{GM}[GeV/c^{2}]",Nbin,0,10);
      TH1F *LikeSign = new TH1F("LikeSign","Invariant Mass Spectrum (LikeSign-method);m^{#mu#mu}_{GM}[GeV/c^{2}]",Nbin,0,10);
      TH1F *InvMassSpectrum_MC = new TH1F("InvMassSpectrum_MC","#rho MC Invariant Mass Spectrum;m^{#mu#mu}_{GM}[GeV/c^{2}]",Nbin,0,10);
      TH1F *InvMassDD = new TH1F("InvMassDD","Invariant Mass Spectrum from D^{0};m^{#mu#mu}_{GM}[GeV/c^{2}]",Nbin,0,10);
      TH1F *InvMassBB = new TH1F("InvMassBB","Invariant Mass Spectrum from B^{0};m^{#mu#mu}_{GM}[GeV/c^{2}]",Nbin,0,10);
      TH1F *InvMassChargedDD = new TH1F("InvMassChargedDD","Invariant Mass Spectrum from D^{+,-};m^{#mu#mu}_{GM}[GeV/c^{2}]",Nbin,0,10);
      TH1F *InvMassChargedBB = new TH1F("InvMassChargedBB","Invariant Mass Spectrum from B^{+,-};m^{#mu#mu}_{GM}[GeV/c^{2}]",Nbin,0,10);
      TH1F *InvMassRho = new TH1F("InvMassRho","Invariant Mass Spectrum from #rho^{0};m^{#mu#mu}_{GM}[GeV/c^{2}]",Nbin,0,10);
      TH1F *PtDist = new TH1F("PtDist","p_{T} distribution;p_{T}^{MC};Entry",10000,0,5);

  // Muon's mass (for calculation of invariant mass)
      Double_t m_mu = 0.106;


  //===================================================================================================================================
  //  Loop of iTrack1
  //===================================================================================================================================

    // First process
          auto iTrack = 0;
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
            auto thisSrcID = label.getSourceID();
            if (thisSrcID == 0)
            {
              o2SimKineTree->GetEntry(thisEvtID);
              thisTrack = &(mcTr->at(thisTrkID));
            }
            else
            {
              iTrack++;
              continue;
            }
            PtDist->Fill(gmTrack.getPt());

            Int_t Pdg1 = thisTrack->GetPdgCode();
            gmTrack.propagateToZhelix(0., field_z);
            Double_t px_1 = gmTrack.getPx();
            Double_t py_1 = gmTrack.getPy();
            Double_t pz_1 = gmTrack.getPz();
            Double_t p_1 = gmTrack.getP();
            Double_t charge_1 = gmTrack.getCharge();
            Int_t MotherTrackID_of_thisTrack = thisTrack->getMotherTrackId();
            MCTrackT<float> *thisTrack_Mother;
            thisTrack_Mother = &(mcTr->at(MotherTrackID_of_thisTrack));
            Int_t motherPdg = thisTrack_Mother->GetPdgCode();
            Int_t MotherTrackID_of_mother_of_thisTrack = thisTrack_Mother->getMotherTrackId();

        //===================================================================================================================================
        //  Loop of iTrack2
        //===================================================================================================================================

            // First process
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
                    auto thisSrcID2 = label2.getSourceID();

                    if (thisSrcID2 == 0)
                    {
                      o2SimKineTree->GetEntry(thisEvtID2);
                      thisTrack2 = &(mcTr->at(thisTrkID2));
                    }
                    else
                    {
                      iTrack2++;
                      continue;
                    }
                    if (thisTrkID >= thisTrkID2) { // to avoid duplication
                      iTrack2++;
                      continue;
                    }

                    Int_t Pdg2 = thisTrack2->GetPdgCode();
                    gmTrack2.propagateToZhelix(0., field_z);
                    Double_t px_2 = gmTrack2.getPx();
                    Double_t py_2 = gmTrack2.getPy();
                    Double_t pz_2 = gmTrack2.getPz();
                    Double_t p_2 = gmTrack2.getP();
                    Double_t charge_2 = gmTrack2.getCharge();
                    Int_t MotherTrackID_of_thisTrack2 = thisTrack2->getMotherTrackId();

                // Invariant Mass
                    Double_t CosTheta = (px_1*px_2+py_1*py_2+pz_1*pz_2)/(p_1*p_2);
                    Double_t InvMass = sqrt(2.*m_mu*m_mu+2.*(sqrt(m_mu*m_mu+p_1*p_1)*sqrt(m_mu*m_mu+p_2*p_2)-p_1*p_2*CosTheta));

                // Mixed event
                    if (thisEvtID != thisEvtID2) {
                          if (charge_1 == charge_2) {
                            if (charge_1 == 1.) { // charge = ++
                              MixedPP->Fill(InvMass);
                            }
                            if (charge_1 == -1.) { // charge = --
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
                    if (thisEvtID == thisEvtID2){

                      if (thisSrcID == 0 && thisSrcID2 == 0){
                        MCTrackT<float> *thisTrack_Mother2;
                        thisTrack_Mother2 = &(mcTr->at(MotherTrackID_of_thisTrack2));
                        Int_t motherPdg2 = thisTrack_Mother2->GetPdgCode();
                        Int_t MotherTrackID_of_mother_of_thisTrack2 = thisTrack_Mother2->getMotherTrackId();
                        std::cout<<"thisTrack[EventID = "<<thisEvtID<<" : TrackID = "<<thisTrkID<<" : PDG = "<<Pdg1<<" : MotherID = "<<MotherTrackID_of_thisTrack<<" : MotherPDG = "<<motherPdg<<"] thisTrack2[EventID = "<<thisEvtID2<<" : TrackID = "<<thisTrkID2<<" : PDG = "<<Pdg2<<" : MotherID = "<<MotherTrackID_of_thisTrack2<<" : MotherPDG = "<<motherPdg2<<"]"<<std::endl;
                        if (MotherTrackID_of_thisTrack != MotherTrackID_of_thisTrack2) {
                          if (((motherPdg == 421) && (motherPdg2 == -421)) || ((motherPdg == -421) && (motherPdg2 == 421))) {
                            if (MotherTrackID_of_mother_of_thisTrack == MotherTrackID_of_mother_of_thisTrack2) {
                              std::cout<<"Correlated DD decay."<<std::endl;
                              InvMassDD->Fill(InvMass);
                            }
                          }
                          if (((motherPdg == 511) && (motherPdg2 == -511)) || ((motherPdg == -511) && (motherPdg2 == 511))) {
                            if (MotherTrackID_of_mother_of_thisTrack == MotherTrackID_of_mother_of_thisTrack2) {
                              std::cout<<"Correlated BB decay."<<std::endl;
                              InvMassBB->Fill(InvMass);
                            }
                          }
                          if (((motherPdg == 411) && (motherPdg2 == -411)) || ((motherPdg == -411) && (motherPdg2 == 411))) {
                            if (MotherTrackID_of_mother_of_thisTrack == MotherTrackID_of_mother_of_thisTrack2) {
                              std::cout<<"Correlated charged(?) DD decay."<<std::endl;
                              InvMassChargedDD->Fill(InvMass);
                            }
                          }
                          if (((motherPdg == 521) && (motherPdg2 == -521)) || ((motherPdg == -521) && (motherPdg2 == 521))) {
                            if (MotherTrackID_of_mother_of_thisTrack == MotherTrackID_of_mother_of_thisTrack2) {
                              std::cout<<"Correlated charged(?) BB decay."<<std::endl;
                              InvMassChargedBB->Fill(InvMass);
                            }
                          }
                        }
                        if (MotherTrackID_of_thisTrack == MotherTrackID_of_thisTrack2) { // from the same mother
                          std::cout<<"dimuon mother's PDG code is "<<motherPdg<<std::endl;
                          if (motherPdg==113 || motherPdg==-113){
                            std::cout<<"Rho meson decay."<<std::endl;
                            InvMassRho->Fill(InvMass);
                          }
                        }
                      }

                      // charge
                          if (charge_1 == charge_2) {
                              if (charge_1 == 1.) { // charge = ++
                                SamePP->Fill(InvMass);
                              }
                              if (charge_1 == -1.) { // charge = --
                                SameMM->Fill(InvMass);
                              }
                              iTrack2++;
                              continue;
                            }
                            SamePM->Fill(InvMass); // charge = +-

                    }
                iTrack2++;
              } // Loop on GMtracks2
        iTrack++;
      } // Loop on GMTracks

  //===================================================================================================================================
  //  Like Sign Method
  //===================================================================================================================================
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
        LikeSign->SetBinContent(bin,Nsignal);
      }

  //===================================================================================================================================
  //  Write histograms on output file
  //===================================================================================================================================
      PtDist->Write();
      SamePM->Write();
      SamePP->Write();
      SameMM->Write();
      MixedPM->Write();
      MixedPP->Write();
      MixedMM->Write();
      LikeSign->Write();
      //InvMassSpectrum_MC->Write();
      InvMassDD->Write();
      InvMassBB->Write();
      InvMassChargedDD->Write();
      InvMassChargedBB->Write();
      InvMassRho->Write();
      //PDGdist->Write();
      outFile.Close();

  return;
}
