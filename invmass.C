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

#include <time.h>

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
void invmass(float chi2Cut = 5.f,
                      const std::string trkFile = "globalfwdtracks.root",
                      const std::string o2sim_KineFile = "bkg_Kine.root",
                      const std::string sig_KineFile = "sgn_1_Kine.root")
{
  //===================================================================================================================================
  //  Read Files and create TTree
  //===================================================================================================================================

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
            o2::dataformats::MCEventHeader *eventHeaderSig = nullptr;
            sigKineTree->SetBranchAddress("MCEventHeader.", &eventHeaderSig);

            Int_t numberOfEvents = o2SimKineTree->GetEntries();
            Int_t numberOfEventsSig = sigKineTree->GetEntries();

            std::cout<<"number of bkg event = "<<numberOfEvents<<" number of sgn event = "<<numberOfEventsSig<<std::endl;

        // Global Muon Tracks
            TFile *trkFileIn = new TFile(trkFile.c_str());
            TTree *gmTrackTree = (TTree *)trkFileIn->Get("GlobalFwdTracks");
            std::vector<GlobalMuonTrack> trackGMVec, *trackGMVecP = &trackGMVec;
            gmTrackTree->SetBranchAddress("fwdtracks", &trackGMVecP);

            vector<o2::MCCompLabel> *mcLabels = nullptr;
            gmTrackTree->SetBranchAddress("MCTruth", &mcLabels);

        // Set all of TTree's entry
            gmTrackTree->GetEntry(0);
            o2SimKineTree->GetEntry(0);
            sigKineTree->GetEntry(0);

        // Get Magnetic Field at the center of MFT
            auto field_z = getZField(0, 0, -61.4); // Get field at Center of MFT

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

        // Muon's Mass
            Double_t m_mu = 0.106;

  //===================================================================================================================================
  //  Loop of iTrack1
  //===================================================================================================================================

        std::cout << "Loop over Global Muon Tracks!" << std::endl;
        auto iTrack = 0;

        for (auto &gmTrack : trackGMVec)
        {

          // First process
              if (gmTrack.getMIDMatchingChi2() < 0 || gmTrack.getMFTMCHMatchingChi2() > chi2Cut) // Filter Muon Tracks and MFTMCH matching Chi2
              {
                iTrack++;
                continue;
              }

              MCTrackT<float> *thisTrack;
              const auto &label = mcLabels->at(iTrack);

              auto thisTrkID = label.getTrackID();
              auto thisEvtID = label.getEventID();
              //auto thisSrcID = label.getSourceID();

              /*if (thisSrcID == 0)
              {
                o2SimKineTree->GetEntry(thisEvtID);
                thisTrack = &(mcTr->at(thisTrkID));
              }
              else if (thisSrcID == 1)
              {
                sigKineTree->GetEntry(thisEvtID);
                thisTrack = &(mcTrSig->at(thisTrkID));
              }
              else
              {
                iTrack++;
                continue;
              }*/

          // Get iTrack1's parameters
              gmTrack.propagateToZhelix(0., field_z);
              Double_t px_1 = gmTrack.getPx();
              Double_t py_1 = gmTrack.getPy();
              Double_t pz_1 = gmTrack.getPz();
              Double_t p_1 = gmTrack.getP();
              Double_t charge_1 = gmTrack.getCharge();

            //===================================================================================================================================
            //  Loop of iTrack2
            //===================================================================================================================================

                      auto iTrack2 = 0;
                        for (auto &gmTrack2 : trackGMVec) {
                          std::chrono::system_clock::time_point  start, end; // 型は auto で可
                          start = std::chrono::system_clock::now(); // 計測開始時間
                                //First process
                                      if (gmTrack2.getMIDMatchingChi2() < 0 || gmTrack2.getMFTMCHMatchingChi2() > chi2Cut) // Filter Muon Tracks and MFTMCH matching Chi2
                                      {
                                        iTrack2++;
                                        continue;
                                      }

                                      MCTrackT<float> *thisTrack2;
                                      const auto &label2 = mcLabels->at(iTrack2);

                                      auto thisTrkID2 = label2.getTrackID();
                                      if (thisTrkID >= thisTrkID2) { // to avoid duplication
                                        iTrack2++;
                                        continue;
                                      }
                                      auto thisEvtID2 = label2.getEventID();
                                      /*auto thisSrcID2 = label2.getSourceID();
                                      if (thisSrcID2 == 0)
                                      {
                                        o2SimKineTree->GetEntry(thisEvtID2);
                                        thisTrack2 = &(mcTr->at(thisTrkID2));
                                      }
                                      else if (thisSrcID2 == 1)
                                      {
                                        sigKineTree->GetEntry(thisEvtID2);
                                        thisTrack2 = &(mcTrSig->at(thisTrkID2));
                                      }
                                      else
                                      {
                                        iTrack2++;
                                        continue;
                                      }*/

                                // Get iTrack2's parameters
                                      gmTrack2.propagateToZhelix(0., field_z);
                                      Double_t px_2 = gmTrack2.getPx();
                                      Double_t py_2 = gmTrack2.getPy();
                                      Double_t pz_2 = gmTrack2.getPz();
                                      Double_t p_2 = gmTrack2.getP();
                                      Double_t charge_2 = gmTrack2.getCharge();

                                // Invariant Mass
                                      Double_t CosTheta = (px_1*px_2+py_1*py_2+pz_1*pz_2)/(p_1*p_2);
                                      Double_t InvMass = sqrt(2.*m_mu*m_mu+2.*(sqrt(m_mu*m_mu+p_1*p_1)*sqrt(m_mu*m_mu+p_2*p_2)-p_1*p_2*CosTheta));

                                // Mixed event
                                      if (thisEvtID != thisEvtID2) {
                                            if (charge_1 == charge_2) {
                                                if (charge_1 == 1.) { // charge = ++
                                                    MixedPP->Fill(InvMass);
                                                    iTrack2++;
                                                    continue;
                                                }
                                                else { // charge = --
                                                    MixedMM->Fill(InvMass);
                                                    iTrack2++;
                                                    continue;
                                                }
                                            }
                                            else {  // charge = +-
                                                  MixedPM->Fill(InvMass);
                                                  iTrack2++;
                                                  continue;
                                            }

                                          }
                                // Same event
                                      else {
                                            if (charge_1 == charge_2) {
                                                  if (charge_1 == 1.) { // charge = ++
                                                      SamePP->Fill(InvMass);
                                                      iTrack2++;
                                                      continue;
                                                  }
                                                  else { // charge = --
                                                      SameMM->Fill(InvMass);
                                                      iTrack2++;
                                                      continue;
                                                  }
                                            }
                                            else {  // charge = +-
                                                SamePM->Fill(InvMass);
                                                iTrack2++;
                                                continue;
                                            }
                                      }
                                iTrack2++;
                        } // Loop on GMtracks2

          iTrack++;
        } // Loop on GMTracks
  std::cout<<"Next is Likesign Loop"<<std::endl;
  //===================================================================================================================================
  //  Like Sign Method
  //===================================================================================================================================
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
          LikeSign->SetBinContent(bin,Nsignal);
        }

  //===================================================================================================================================
  //  Write histograms on output file
  //===================================================================================================================================
        SamePM->Write();
        SamePP->Write();
        SameMM->Write();
        MixedPM->Write();
        MixedPP->Write();
        MixedMM->Write();
        LikeSign->Write();

        outFile.Close();

  return;
}
