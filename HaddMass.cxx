void HaddMass() {
  TFile *fin = new TFile("InvMass.root", "read");

  TH1F *SamePP = (TH1F*)fin->Get("SamePP");
  TH1F *SameMM = (TH1F*)fin->Get("SameMM");
  TH1F *SamePM = (TH1F*)fin->Get("SamePM");
  TH1F *MixedPP = (TH1F*)fin->Get("MixedPP");
  TH1F *MixedMM = (TH1F*)fin->Get("MixedMM");
  TH1F *MixedPM = (TH1F*)fin->Get("MixedPM");
  TH1F *LikeSign = (TH1F*)fin->Get("LikeSign");
  Int_t Nbin = 1000;
  //TH1F *LikeSign = new TH1F("LikeSign","Invariant Mass Spectrum (LikeSign-method);m^{#mu#mu}_{GM}[GeV/c^{2}]",Nbin,0,10);
  TFile *fout = new TFile("InvariantMass.root", "recreate");

  Int_t bin;
  /*for (bin = 0; bin < Nbin; bin++) {
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
    Double_t Signal = NsamePM-2.*Rfactor*sqrt(NsamePP*NsameMM);
    //std::cout<<"bin = "<<bin<<endl;
    //std::cout<<"Same PP = "<<SamePP<<" : MM = "<<SameMM<<" : PM = "<<SamePM<<" : Mixed PP = "<<NmixedPP<<" : MM = "<<NmixedMM<<" : PM = "<<NmixedPM<<endl;
    //std::cout<<"Rfactor = "<<Rfactor<<" : Signal = "<<Signal<<endl;
    //std::cout<<"-------------------------------------------------------------------------"<<endl;
    LikeSign->SetBinContent(bin,Signal);
  }*/

  //LikeSign->Draw();

  Int_t bkg_range_min = 80;
  Int_t bkg_range_max = 999;
  Int_t bin_sig_range_min = 200;
  Int_t bin_sig_range_max = 450;

  TH1F *bkg_histo = new TH1F("bkg_histo","bkg_histo",1000,0,10);
  for (bin = bkg_range_min; bin<bkg_range_max; bin++) {
    bkg_histo->SetBinContent(bin+1,LikeSign->GetBinContent(bin+1));
  }
  //bkg_histo = (TH1F*)LikeSign->Clone("bkg_histo");
  bkg_histo->Draw();
  for(bin = bin_sig_range_min; bin < bin_sig_range_max; bin++) {
    bkg_histo->SetBinContent(bin+1,0);
  }
  //bkg_histo->Draw();
  TF1 *FitBkg = new TF1("FitBkg","pol5",bkg_range_min,bkg_range_max);
  //Double_t constant = 50.;
  //Double_t variable = 1.;
  //FitBkg->SetParameters(constant,variable);
  bkg_histo->Fit(FitBkg,"R","",bkg_range_min,bkg_range_max);
  //bkg_histo->Draw("same");

  TH1F *sig_histo = (TH1F*)LikeSign->Clone("sig_histo");
  for (bin = bkg_range_min; bin<bkg_range_max; bin++) {
    Double_t bin_center_val = LikeSign->GetXaxis()->GetBinCenter(bin+1);
    Double_t bkg = FitBkg->Eval(bin_center_val);
    Double_t e_bkg = 0.;
    Double_t data = sig_histo->GetBinContent(bin+1);
    Double_t e_data = sig_histo->GetBinError(bin+1);
    Double_t sig = data-bkg;
    Double_t e_sig = sqrt(pow(e_data,2)+pow(e_bkg,2));

    sig_histo->SetBinContent(bin+1,sig);
    sig_histo->SetBinError(bin+1,e_sig);
  }

  Double_t peak_mean_by_my_eye = 3.1;
  Double_t peak_width_by_my_eye = 0.2;

  TF1 *FitSig = new TF1("FitSig","gaus",bkg_range_min,bkg_range_max);
  FitSig->SetParameters(sig_histo->GetBinContent(sig_histo->GetMaximumBin()),peak_mean_by_my_eye,peak_width_by_my_eye);
  sig_histo->Fit(FitSig,"R","",bkg_range_min,bkg_range_max);
  //sig_histo->Draw("same");

  TH1F *sig_histo_2 = (TH1F*)LikeSign->Clone("sig_histo_2");
  for (bin = bkg_range_min; bin<bkg_range_max; bin++) {
    Double_t bin_center_val = LikeSign->GetXaxis()->GetBinCenter(bin+1);
    Double_t bkg = FitBkg->Eval(bin_center_val);
    Double_t e_bkg = 0.;
    Double_t sig_1 = FitSig->Eval(bin_center_val);
    Double_t e_sig_1 = 0.;
    Double_t data = sig_histo_2->GetBinContent(bin+1);
    Double_t e_data = sig_histo_2->GetBinError(bin+1);
    Double_t sig = data-bkg-sig_1;
    Double_t e_sig = sqrt(pow(e_data,2)+pow(e_bkg,2));

    sig_histo_2->SetBinContent(bin+1,sig);
    sig_histo_2->SetBinError(bin+1,e_sig);
  }

  Double_t peak_mean_by_my_eye_2 = 3.7;
  Double_t peak_width_by_my_eye_2 = 0.6;

  TF1 *FitSig_2 = new TF1("FitSig_2","gaus",bkg_range_min,bkg_range_max);
  FitSig_2->SetParameters(sig_histo_2->GetBinContent(sig_histo_2->GetMaximumBin()),peak_mean_by_my_eye,peak_width_by_my_eye);
  sig_histo_2->Fit(FitSig_2,"R","",bkg_range_min,bkg_range_max);
  //sig_histo_2->Draw("same");

  SamePM->Write();
  SamePP->Write();
  SameMM->Write();
  MixedPM->Write();
  MixedPP->Write();
  MixedMM->Write();
  LikeSign->Write();
  sig_histo->Write();
  sig_histo_2->Write();
  bkg_histo->Write();
  fout->Close();

  return;
}
