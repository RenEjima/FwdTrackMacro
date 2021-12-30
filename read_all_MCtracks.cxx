const std::string o2sim_KineFile = “o2sim_Kine.root”;
TFile *o2sim_KineFileIn = new TFile(o2sim_KineFile.c_str());
TTree *o2SimKineTree = (TTree *)o2sim_KineFileIn->Get("o2sim");

vector<MCTrackT<float>> *mcTr = nullptr;
o2SimKineTree->SetBranchAddress("MCTrack", &mcTr);

Int numberOfEvents = o2SimKineTree->GetEntries();

for (int iEvent = 0; iEvent < numberOfEvents; iEvent++) {
   o2SimKineTree->GetEntry(iEvent);
   int nMCTracks = (*mcTr).size();

   for (Int_t iMCTrack=0; iMCTrack<(*mcTr).size();++iMCTrack) {
       MCTrackT<float> *thisTrack = &(*mcTr).at(iMCTrack);
       float eta = thisTrack->GetEta()
       float pt = thisTrack->GetPt()
   }
}
