#include "DataFormatsMCH/TrackMCH.h"

#include<iostream>
#include<memory>


void MCHTrackChecker(const std::string trkFile = "mchtracks.root"){

  // Set global muon track (MFT + MCH) information
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

  // Set MC track information
    std::unique_ptr<TFile> o2sim_KineFileIn(new TFile("o2sim_Kine.root"));
    std::unique_ptr<TTree> o2SimKineTree((TTree*)o2sim_KineFileIn->Get("o2sim"));
    vector<o2::MCTrackT<float>>* mcTr = nullptr;
    o2SimKineTree->SetBranchAddress("MCTrack", &mcTr);
    //o2::dataformats::MCEventHeader* eventHeader = nullptr;
    //o2SimKineTree->SetBranchAddress("MCEventHeader.", &eventHeader);

  auto iTrack = 0;
  int nPions = 0;
  int nMuons = 0;
  int nElectrons = 0;
  int nKaons = 0;
  int nProtons = 0;
  //Load MC track information
  o2::MCTrackT<float>* thisTrack;
  for (auto& mchTrack : trackMCHVec) {

    //Get label information
    const auto& label = trkLabels->at(iTrack);

    //Get MC track id
    auto trkID = label.getTrackID();
    //Get Event id
    auto evtID = label.getEventID();
    printf("EventID = %d : TrackID = %d \n ",evtID,trkID);

    if(label.getSourceID()==0){
    o2SimKineTree->GetEntry(evtID);
    thisTrack = &(mcTr->at(trkID));
    }else{
      iTrack++;
      continue;
    }
    //MC track information
    auto vx_MC = thisTrack->GetStartVertexCoordinatesX();//get particle production x-posiiton
    auto vy_MC = thisTrack->GetStartVertexCoordinatesY();//get particle production y-posiiton
    auto vz_MC = thisTrack->GetStartVertexCoordinatesZ();//get particle production z-posiiton
    auto Pt_MC = thisTrack->GetPt();//get particle pt
    auto p_MC = sqrt(pow((thisTrack->GetStartVertexMomentumX()),2.)+pow((thisTrack->GetStartVertexMomentumY()),2.)+pow( (thisTrack->GetStartVertexMomentumZ()),2.));//get particle total momentum
    auto pdgcode_MC = thisTrack->GetPdgCode();//get particle pdg code

    printf("EventID = %d : TrackID = %d : vx_MC = %4f : vy_MC = %4f : vz_MC = %4f : Pt_MC = %4f : vP_MC = %4f : PdgCode = %d\n ",evtID,trkID,vx_MC,vy_MC,vz_MC,Pt_MC,p_MC,pdgcode_MC);

    if(thisTrack->isPrimary()){
      if(fabs(pdgcode_MC)==211){
        nPions++;
      }
      if(fabs(pdgcode_MC)==13){
        nMuons++;
      }
      if(fabs(pdgcode_MC)==11){
        nElectrons++;
      }
      if(fabs(pdgcode_MC)==321){
        nKaons++;
      }
      if(fabs(pdgcode_MC)==2212){
        nProtons++;
      }
    }

    iTrack++;
  }
printf("iTrack = %7d \n",iTrack);
printf("N_Pion = %d : N_Muon = %d : N_Electron = %d : N_Kaons = %d : NProton = %d \n",nPions,nMuons,nElectrons,nKaons,nProtons);
}
