#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <TFile.h>
#include <TTree.h>
#include <TH2.h>

#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "MFTBase/GeometryTGeo.h"
#include "Math/GenVector/Cartesian3D.h"
#include "MathUtils/Utils.h"
#include "GPUCommonDef.h"
#include "SimulationDataFormat/MCTrack.h"
#include "TMatrixD.h"

#endif

void MFTTrackInOut() {

  // class CompClusterExt : public CompCluster
  // DataFormats/Detectors/ITSMFT/common/include/DataFormatsITSMFT/CompCluster.h
  using o2::itsmft::CompClusterExt;
  using o2::MCTrack;

  //Read MCtracks
  TFile fileK("o2sim_Kine.root");
  TTree* kineTree = (TTree*)fileK.Get("o2sim");
  std::vector<o2::MCTrack> mcTrkVec, *mcTrkVecP = &mcTrkVec;
  kineTree->SetBranchAddress("MCTrack",&mcTrkVecP);

  int nEventsKine = kineTree->GetEntries();
  //printf("Number of MCevents %d \n", nEventsKine);

  // Geometry and matrix transformations
  std::string inputGeom = "o2sim_geometry.root";
  o2::base::GeometryManager::loadGeometry(inputGeom);
  auto gman = o2::mft::GeometryTGeo::Instance();
  gman->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));

  // Cluster pattern dictionary
  std::string dictfile = "MFTdictionary.bin";
  o2::itsmft::TopologyDictionary dict;
  std::ifstream file(dictfile.c_str());
  if (file.good()) {
    //printf("Running with dictionary: %s \n", dictfile.c_str());
    dict.readBinaryFile(dictfile);
  } else {
    //printf("Can not run without dictionary !\n");
    return;
  }

  // Clusters

  TFile fileC("mftclusters.root");
  TTree *clsTree = (TTree*)fileC.Get("o2sim");
  std::vector<CompClusterExt> clsVec, *clsVecP = &clsVec;
  clsTree->SetBranchAddress("MFTClusterComp", &clsVecP);
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* clsLabels = nullptr;
  if (clsTree->GetBranch("MFTClusterMCTruth")) {
    clsTree->SetBranchAddress("MFTClusterMCTruth", &clsLabels);
  } else {
    //printf("No Monte-Carlo information in this file\n");
    return;
  }

  int nEntries = clsTree->GetEntries();
  //printf("Number of entries in clusters tree %d \n", nEntries);

  clsTree->GetEntry(0);

  int nClusters = clsVec.size();
  //printf("Number of clusters %d \n", nClusters);

  // Tracks
  // class TrackMFT : public o2::track::TrackParCovFwd
  // class TrackParCovFwd : public TrackParFwd
  // DataFormats/Detectors/ITSMFT/MFT/include/DataFormatsMFT/TrackMFT.h
  // DataFormats/Reconstruction/include/ReconstructionDataFormats/TrackFwd.h
  TFile fileT("mfttracks.root");
  TTree *trackTree = (TTree*)fileT.Get("o2sim");
  std::vector<o2::mft::TrackMFT> trackVec, *trackVecP = &trackVec;
  trackTree->SetBranchAddress("MFTTrack", &trackVecP);
  std::vector<o2::MCCompLabel>* trkLabels = nullptr;
  if (trackTree->GetBranch("MFTTrackMCTruth")) {
    trackTree->SetBranchAddress("MFTTrackMCTruth", &trkLabels);
  } else {
    //printf("No Monte-Carlo information in this file\n");
    return;
  }

  std::vector<int> trackExtClsVec, *trackExtClsVecP = &trackExtClsVec;
  trackTree->SetBranchAddress("MFTTrackClusIdx", &trackExtClsVecP);

  trackTree->GetEntry(0);

  int srcID, trkID, evnID;
  bool fake;
  int iTrack = 0;

  std::string outfilename = "DeltaZ.root";
  TFile outFile(outfilename.c_str(), "RECREATE");

  Double_t Delta_Z;
  Int_t trackID;
  Int_t eventID;
  TTree *DeltaZtree = new TTree("DeltaZtree", "DeltaZtree");
  DeltaZtree->Branch("Delta_Z", &Delta_Z, "Delta_Z/D");
  DeltaZtree->Branch("MFTTrackID", &trackID, "trackID/I");
  DeltaZtree->Branch("EventID", &eventID, "eventID/I");

  Double_t trkX;
  Double_t trkY;
  Double_t trkZ;
  //auto trkLabel = trkLabels->getLabels(iTrack);
  Double_t trkOutX;
  Double_t trkOutY;
  Double_t trkOutZ;

  Double_t trkXvtx;
  Double_t trkYvtx;
  Double_t trkZvtx;

  //Double_t mcVtxX, mcVtxY, mcVtxZ;

  TH1F *DeltaZ = new TH1F("DeltaZ","#DeltaZ;Entry",1000,20,31);
  TH1F *priDeltaZ = new TH1F("priDeltaZ","#DeltaZ of primary MFTtracks;Entry",1000,20,31);
  TH1F *secDeltaZ = new TH1F("secDeltaZ","#DeltaZ of secondary MFTtracks;Entry",1000,20,31);
  TH1F *DeltaZmu = new TH1F("DeltaZmu","#DeltaZ;Entry",1000,20,31);
  TH1F *DeltaZmupri = new TH1F("DeltaZmupri","#DeltaZ;Entry",1000,20,31);
  TH1F *DeltaZmusec = new TH1F("DeltaZmusec","#DeltaZ;Entry",1000,20,31);
  TH1F *DeltaZpi = new TH1F("DeltaZpi","#DeltaZ;Entry",1000,20,31);
  TH1F *DeltaZpipri = new TH1F("DeltaZpipri","#DeltaZ;Entry",1000,20,31);
  TH1F *DeltaZpisec = new TH1F("DeltaZpisec","#DeltaZ;Entry",1000,20,31);
  TH1F *DeltaZka = new TH1F("DeltaZka","#DeltaZ;Entry",1000,20,31);
  TH1F *DeltaZkapri = new TH1F("DeltaZkapri","#DeltaZ;Entry",1000,20,31);
  TH1F *DeltaZkasec = new TH1F("DeltaZkasec","#DeltaZ;Entry",1000,20,31);
  TH1F *DeltaZpro = new TH1F("DeltaZpro","#DeltaZ;Entry",1000,20,31);
  TH1F *DeltaZpropri = new TH1F("DeltaZpropri","#DeltaZ;Entry",1000,20,31);
  TH1F *DeltaZprosec = new TH1F("DeltaZprosec","#DeltaZ;Entry",1000,20,31);
  TH1F *DeltaZele = new TH1F("DeltaZele","#DeltaZ;Entry",1000,20,31);
  TH1F *DeltaZelepri = new TH1F("DeltaZelepri","#DeltaZ;Entry",1000,20,31);
  TH1F *DeltaZelesec = new TH1F("DeltaZelesec","#DeltaZ;Entry",1000,20,31);

  TH1F *StartVertexZpri = new TH1F("StartVertexZpri","primary track's StartVertexZ;Z[cm];Entry",1000,-90,10);
  TH1F *StartVertexZsec = new TH1F("StartVertexZsec","secondary track's StartVertexZ;Z[cm];Entry",1000,-90,10);

  TH2F *Vertex = new TH2F("Vertex","Vertexing of MFTtracks;X;Y",1000,-10,10,1000,-10,10);
  TH2F *priVertex = new TH2F("priVertex","Vertexing of primary MFTtracks;X;Y",1000,-10,10,1000,-10,10);
  TH2F *secVertex = new TH2F("secVertex","Vertexing of secondary MFTtracks;X;Y",1000,-10,10,1000,-10,10);
  //TH3F *VertexMC =new TH3F("Vertex3D","Vertex3D;x[cm];y[cm];z[cm]",1000,-50,50,1000,-50,50,1000,-90,10);
  //TH3F *priVertexMC =new TH3F("priVertex3D","primary Vertex3D;x[cm];y[cm];z[cm]",1000,-50,50,1000,-50,50,1000,-90,10);
  //TH3F *secVertexMC =new TH3F("secVertex3D","secondary Vertex3D;x[cm];y[cm];z[cm]",1000,-50,50,1000,-50,50,1000,-90,10);

  for (auto &track : trackVec) {
    trkX = track.getX();
    trkY = track.getY();
    trkZ = track.getZ();
    //auto trkLabel = trkLabels->getLabels(iTrack);
    auto pretrkLabels = *trkLabels;
    auto trkLabel = pretrkLabels.at(iTrack);
    //trkLabel[0].printf();
    eventID = trkLabel.getEventID();
    trackID = trkLabel.getTrackID();
    auto outParam = track.getOutParam();
    trkOutX = outParam.getX();
    trkOutY = outParam.getY();
    trkOutZ = outParam.getZ();
    Delta_Z = trkZ - trkOutZ;
    printf("iTrack %3d   isCA %1d   x,y,z-in  %7.3f  %7.3f  %7.3f  x,y,z-out  %7.3f  %7.3f  %7.3f   evID %2d  labels   trID %4d   Delta_Z  %7.3f \n" , iTrack, track.isCA(), trkX, trkY, trkZ, trkOutX, trkOutY, trkOutZ, eventID, trackID, Delta_Z);
    DeltaZtree->Fill();
    track.propagateToZlinear(0.);
    trkXvtx = track.getX();
    trkYvtx = track.getY();
    trkZvtx = track.getZ();
    //printf("x,y,z-priVtx  %7.3f  %7.3f  %7.3f\n",trkXvtx, trkYvtx, trkZvtx);

    /*
    printf("Now Starting the Loop of iLab\n");
    for (auto ilab = 0; ilab < trkLabel.maxTrackID(); ++ilab) {
      printf("Now In the Loop of iLab");
      auto trkID = trkLabel.getTrackID();

      //printf("getTrackID = %4d \n", trkID);
    }
    printf("\nNow End of the Loop of iLab\n");
    //printf("\n");
    */

    auto ncls = track.getNumberOfPoints();
    auto offset = track.getExternalClusterIndexOffset();

    for (int icls = 0; icls < ncls; ++icls) {
      //printf("Now In the Loop of icls");
      auto clsEntry = trackExtClsVec[offset + icls];
      //auto cluster = clsVec[clsEntry];
      auto& clsLabel = (clsLabels->getLabels(clsEntry))[0];
      if (!clsLabel.isNoise()) {
	clsLabel.get(trkID, evnID, srcID, fake);
	//auto chipID = cluster.getChipID();
	//auto pattID = cluster.getPatternID();
	//o2::math_utils::Point3D<float> locC;
	//int npix = 0;
  /*
	if (pattID == o2::itsmft::CompCluster::InvalidPatternID || dict.isGroup(pattID)) {
	  // temporary fix ...
	  //printf("temporary fix for group pattern ...\n");
	  locC = dict.getClusterCoordinates(cluster);

	  //o2::itsmft::ClusterPattern patt(pattIt);
	  //locC = dict.getClusterCoordinates(cluster, patt);
	} else {
	  locC = dict.getClusterCoordinates(cluster);
	  npix = dict.getNpixels(pattID);
	}
	// Transformation to the local --> global
	auto gloC = gman->getMatrixL2G(chipID) * locC;
	track.propagateToZlinear(gloC.Z());
	auto deltaX = gloC.X()-track.getX();
	auto deltaY = gloC.Y()-track.getY();
  */

      }//if (!clsLabel.isNoise())
    }//for (int icls = 0; icls < ncls; ++icls) {

    kineTree->GetEntry(evnID);
    int nMCTracks = mcTrkVec.size();
    //printf("Event %d has %d MC tracks\n", evnID, nMCTracks);
    int nt = 0;
    for (auto mcTrack : mcTrkVec) {
      //printf("Now In the Loop of mcTrack");
      if(nt == trkID){
        printf("\nEvent ID %4d   MCTrack ID %4d   PDG %4d   isSec %d   E %7.3f \n",
        evnID,
	nt,
	mcTrack.GetPdgCode(),
               //we can get particclename = TDatabasePDG::Instance()->GetParticle(mcTrack.GetPdgCode())->GetName(),
               mcTrack.isSecondary(),
               mcTrack.GetEnergy());
        auto mcVtxX = mcTrack.GetStartVertexCoordinatesX();
        auto mcVtxY = mcTrack.GetStartVertexCoordinatesY();
        auto mcVtxZ = mcTrack.GetStartVertexCoordinatesZ();
	//printf("Start Vertex Coordinate x,y,z   %7.3f   %7.3f   %7.3f\n", mcVtxX, mcVtxY,mcVtxZ);
    if(TMath::Abs(mcTrack.GetPdgCode())==13){
      DeltaZmu->Fill(Delta_Z);
      if(mcTrack.isSecondary()){
        DeltaZ->Fill(Delta_Z);
        Vertex->Fill(trkXvtx,trkYvtx);
        secDeltaZ->Fill(Delta_Z);
        secVertex->Fill(trkXvtx,trkYvtx);
        DeltaZmusec->Fill(Delta_Z);
	//VertexMC->Fill(mcVtxX,mcVtxY,mcVtxZ);
	//secVertexMC->Fill(mcVtxX,mcVtxY,mcVtxZ);
	StartVertexZsec->Fill(mcVtxZ);
      }
      if(!mcTrack.isSecondary()){
        DeltaZ->Fill(Delta_Z);
        Vertex->Fill(trkXvtx,trkYvtx);
        priDeltaZ->Fill(Delta_Z);
        priVertex->Fill(trkXvtx,trkYvtx);
        DeltaZmupri->Fill(Delta_Z);
	//VertexMC->Fill(mcVtxX,mcVtxY,mcVtxZ);
        //priVertexMC->Fill(mcVtxX,mcVtxY,mcVtxZ);
	StartVertexZpri->Fill(mcVtxZ);
      }
    }
      if(TMath::Abs(mcTrack.GetPdgCode())==211){
      DeltaZpi->Fill(Delta_Z);
      if(mcTrack.isSecondary()){
        DeltaZ->Fill(Delta_Z);
        Vertex->Fill(trkXvtx,trkYvtx);
        secDeltaZ->Fill(Delta_Z);
        secVertex->Fill(trkXvtx,trkYvtx);
        DeltaZpisec->Fill(Delta_Z);
	//VertexMC->Fill(mcVtxX,mcVtxY,mcVtxZ);
        //secVertexMC->Fill(mcVtxX,mcVtxY,mcVtxZ);
	StartVertexZsec->Fill(mcVtxZ);
      }
      if(!mcTrack.isSecondary()){
        DeltaZ->Fill(Delta_Z);
        Vertex->Fill(trkXvtx,trkYvtx);
        priDeltaZ->Fill(Delta_Z);
        priVertex->Fill(trkXvtx,trkYvtx);
        DeltaZpipri->Fill(Delta_Z);
	//VertexMC->Fill(mcVtxX,mcVtxY,mcVtxZ);
        //priVertexMC->Fill(mcVtxX,mcVtxY,mcVtxZ);
	StartVertexZpri->Fill(mcVtxZ);
      }
    }

    if(TMath::Abs(mcTrack.GetPdgCode())==321){
      DeltaZka->Fill(Delta_Z);
      if(mcTrack.isSecondary()){
        DeltaZ->Fill(Delta_Z);
        Vertex->Fill(trkXvtx,trkYvtx);
        secDeltaZ->Fill(Delta_Z);
        secVertex->Fill(trkXvtx,trkYvtx);
        DeltaZkasec->Fill(Delta_Z);
	//VertexMC->Fill(mcVtxX,mcVtxY,mcVtxZ);
        //secVertexMC->Fill(mcVtxX,mcVtxY,mcVtxZ);
	StartVertexZsec->Fill(mcVtxZ);
      }
      if(!mcTrack.isSecondary()){
        DeltaZ->Fill(Delta_Z);
        Vertex->Fill(trkXvtx,trkYvtx);
        priDeltaZ->Fill(Delta_Z);
        priVertex->Fill(trkXvtx,trkYvtx);
        DeltaZkapri->Fill(Delta_Z);
	//VertexMC->Fill(mcVtxX,mcVtxY,mcVtxZ);
        //priVertexMC->Fill(mcVtxX,mcVtxY,mcVtxZ);
	StartVertexZpri->Fill(mcVtxZ);
      }
    }

    if(TMath::Abs(mcTrack.GetPdgCode())==2212){
      DeltaZpro->Fill(Delta_Z);
      if(mcTrack.isSecondary()){
        DeltaZ->Fill(Delta_Z);
        Vertex->Fill(trkXvtx,trkYvtx);
        secDeltaZ->Fill(Delta_Z);
        secVertex->Fill(trkXvtx,trkYvtx);
        DeltaZprosec->Fill(Delta_Z);
	//VertexMC->Fill(mcVtxX,mcVtxY,mcVtxZ);
        //secVertexMC->Fill(mcVtxX,mcVtxY,mcVtxZ);
	StartVertexZsec->Fill(mcVtxZ);
      }
      if(!mcTrack.isSecondary()){
        DeltaZ->Fill(Delta_Z);
        Vertex->Fill(trkXvtx,trkYvtx);
        priDeltaZ->Fill(Delta_Z);
        priVertex->Fill(trkXvtx,trkYvtx);
        DeltaZpropri->Fill(Delta_Z);
	//VertexMC->Fill(mcVtxX,mcVtxY,mcVtxZ);
        //priVertexMC->Fill(mcVtxX,mcVtxY,mcVtxZ);
	StartVertexZpri->Fill(mcVtxZ);
      }
    }

    if(TMath::Abs(mcTrack.GetPdgCode())==11){
      DeltaZele->Fill(Delta_Z);
      if(mcTrack.isSecondary()){
        DeltaZ->Fill(Delta_Z);
        Vertex->Fill(trkXvtx,trkYvtx);
        secDeltaZ->Fill(Delta_Z);
        secVertex->Fill(trkXvtx,trkYvtx);
        DeltaZelesec->Fill(Delta_Z);
	//VertexMC->Fill(mcVtxX,mcVtxY,mcVtxZ);
        //secVertexMC->Fill(mcVtxX,mcVtxY,mcVtxZ);
	StartVertexZsec->Fill(mcVtxZ);
      }
      if(!mcTrack.isSecondary()){
        DeltaZ->Fill(Delta_Z);
        Vertex->Fill(trkXvtx,trkYvtx);
        priDeltaZ->Fill(Delta_Z);
        priVertex->Fill(trkXvtx,trkYvtx);
        DeltaZelepri->Fill(Delta_Z);
	//VertexMC->Fill(mcVtxX,mcVtxY,mcVtxZ);
        //priVertexMC->Fill(mcVtxX,mcVtxY,mcVtxZ);
	StartVertexZpri->Fill(mcVtxZ);
      }
    }
      }//if(nt == trkID)
      nt = nt+1;
    }//for{auto mcTrack : mcTrkVec}
    ++iTrack;
  }

  DeltaZtree->Write();
  DeltaZ->Write();
  priDeltaZ->Write();
  secDeltaZ->Write();
  Vertex->Write();
  priVertex->Write();
  secVertex->Write();
  //VertexMC->Write();
  StartVertexZpri->Write();
  StartVertexZsec->Write();
  //priVertexMC->Write();
  //secVertexMC->Write();
  DeltaZmu->Write();
  DeltaZpi->Write();
  DeltaZka->Write();
  DeltaZpro->Write();
  DeltaZele->Write();
  DeltaZmupri->Write();
  DeltaZpipri->Write();
  DeltaZkapri->Write();
  DeltaZpropri->Write();
  DeltaZelepri->Write();
  DeltaZmusec->Write();
  DeltaZpisec->Write();
  DeltaZkasec->Write();
  DeltaZprosec->Write();
  DeltaZelesec->Write();
  fileT.Close();
  fileC.Close();

}
