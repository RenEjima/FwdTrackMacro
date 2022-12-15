#include "DataFormatsMCH/TrackMCH.h"

#include<iostream>
#include<memory>

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

void MCHTrackChecker(const std::string trkFile = "mchtracks.root",
                     const std::string muonFile = "muontracks.root",
                     const std::string SgntrkFile = "sgn_1_Kine.root",
                     const std::string BkgtrkFile = "bkg_Kine.root"){

  // Set MCH track information
  std::unique_ptr<TFile> trkFileIn(new TFile(trkFile.c_str()));
  std::unique_ptr<TTree> mchTrackTree((TTree*)trkFileIn->Get("o2sim"));
  std::vector<o2::mch::TrackMCH> trackMCHVec, *trackMCHVecP = &trackMCHVec;
  mchTrackTree->SetBranchAddress("tracks", &trackMCHVecP);
  std::vector<o2::MCCompLabel>* trkLabels = nullptr;
  if (mchTrackTree->GetBranch("tracklabels")) {
    mchTrackTree->SetBranchAddress("tracklabels", &trkLabels);
  } else {
    printf("No Monte-Carlo information in this file\n");
    return;
  }
  mchTrackTree->GetEntry(0);

  // Set MUON(MCH+MID) track information
  std::unique_ptr<TFile> muonFileIn(new TFile(muonFile.c_str()));
  std::unique_ptr<TTree> muonTrackTree((TTree*)muonFileIn->Get("o2sim"));
  std::vector<o2::dataformats::TrackMCHMID> trackMUONVec, *trackMUONVecP = &trackMUONVec;
  muonTrackTree->SetBranchAddress("tracks", &trackMUONVecP);
  std::vector<o2::MCCompLabel>* MUONtrkLabels = nullptr;
  if (muonTrackTree->GetBranch("tracklabels")) {
    muonTrackTree->SetBranchAddress("tracklabels", &MUONtrkLabels);
  } else {
    printf("No Monte-Carlo information in this file\n");
    return;
  }
  muonTrackTree->GetEntry(0);
  std::vector<int> MCHidMID;
  for(auto& muonTrack : trackMUONVec){
    Int_t MCHId = muonTrack.getMCHRef().getIndex();
    MCHidMID.push_back(MCHId);
  }

  /*std::cout<<"MUONtrack's MCHtrack ID list"<<std::endl;

  for (size_t i = 0; i < MCHidMID.size(); ++i) {
        cout << MCHidMID.at(i) << "; ";
    }
    cout << endl;
*/
  // Set MC track information
  std::unique_ptr<TFile> Sgn_KineFileIn(new TFile(SgntrkFile.c_str()));
  std::unique_ptr<TTree> SgnKineTree((TTree*)Sgn_KineFileIn->Get("o2sim"));
  vector<o2::MCTrackT<float>>* mcTrSgn = nullptr;
  SgnKineTree->SetBranchAddress("MCTrack", &mcTrSgn);
  //o2::dataformats::MCEventHeader* eventHeader = nullptr;
  //SgnKineTree->SetBranchAddress("MCEventHeader.", &eventHeader);

  std::unique_ptr<TFile> Bkg_KineFileIn(new TFile(BkgtrkFile.c_str()));
  std::unique_ptr<TTree> BkgKineTree((TTree*)Bkg_KineFileIn->Get("o2sim"));
  vector<o2::MCTrackT<float>>* mcTrBkg = nullptr;
  BkgKineTree->SetBranchAddress("MCTrack", &mcTrBkg);

  Int_t iTrack = 0;//N_track counter

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

  //Load MC track information
  o2::MCTrackT<float>* thisTrack;

  for (auto& mchTrack : trackMCHVec) {

    //Get label information
    const auto& label = trkLabels->at(iTrack);
    //Get MC track id
    auto trkID = label.getTrackID();
    //Get Event id
    auto evtID = label.getEventID();
    /*
    //Get Source id
    auto srcID = label.getSourceID();
    //select only comes from collision
    if(srcID==0){
      BkgKineTree->GetEntry(evtID);
      thisTrack = &(mcTrBkg->at(trkID));
    }
    if(srcID==1){
      SgnKineTree->GetEntry(evtID);
      thisTrack = &(mcTrSgn->at(trkID));
    }
    else{
      iTrack++;
      continue;
    }
    */
    //mchTrack.propagateToZhelix(0., field_z);

    bool found = std::find(MCHidMID.begin(), MCHidMID.end(), iTrack) != MCHidMID.end();
    //std::cout<<"iTrack = "<<iTrack<<" : thisEvent = "<<evtID<<" : thisMCHtrackID = "<<trkID<<" : found = "<<found<<std::endl;
    if (found){
      Double_t px_1 = mchTrack.getPx();
      Double_t py_1 = mchTrack.getPy();
      Double_t pz_1 = mchTrack.getPz();
      Double_t p_1 = mchTrack.getP();
      Double_t charge_1 = mchTrack.getSign();

      Int_t iTrack2 = 0;
      for(auto& mchTrack2 : trackMCHVec){
        //Get label information
        const auto& label2 = trkLabels->at(iTrack2);
        //Get MC track id
        auto trkID2 = label2.getTrackID();
        //Get Event id
        auto evtID2 = label2.getEventID();
        if(trkID>=trkID2){
          iTrack2++;
          continue;
        }
        /*
        //Get Source id
        auto srcID2 = label2.getSourceID();
        //select only comes from collision
        if(srcID2==0){
          BkgKineTree->GetEntry(evtID2);
          thisTrack2 = &(mcTrBkg->at(trkID2));
        }
        if(srcID2==1){
          SgnKineTree->GetEntry(evtID2);
          thisTrack2 = &(mcTrSgn->at(trkID2));
        }
        else{
          iTrack++;
          continue;
        }
*/
        //mchTrack2.propagateToZhelix(0., field_z);
        bool found2 = std::find(MCHidMID.begin(), MCHidMID.end(), iTrack2) != MCHidMID.end();
        if (found2){
          Double_t px_2 = mchTrack2.getPx();
          Double_t py_2 = mchTrack2.getPy();
          Double_t pz_2 = mchTrack2.getPz();
          Double_t p_2 = mchTrack2.getP();
          Double_t charge_2 = mchTrack2.getSign();

          Double_t CosTheta = (px_1*px_2+py_1*py_2+pz_1*pz_2)/(p_1*p_2);
          Double_t InvMass = sqrt(2.*m_mu*m_mu+2.*(sqrt(m_mu*m_mu+p_1*p_1)*sqrt(m_mu*m_mu+p_2*p_2)-p_1*p_2*CosTheta));

          // Mixed event
                if (evtID != evtID2) {
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
        }
        iTrack2++;
      }
    }
    iTrack++;
  }
  std::cout<<"Next is Likesign Loop"<<std::endl;

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
