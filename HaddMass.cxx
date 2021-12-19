void HaddMass() {
  TFile *fin = new TFile("InvMass500x50.root", "read");

  TH1F *SamePP = (TH1F*)fin->Get("SamePP");
  TH1F *SameMM = (TH1F*)fin->Get("SameMM");
  TH1F *SamePM = (TH1F*)fin->Get("SamePM");
  TH1F *MixedPP = (TH1F*)fin->Get("MixedPP");
  TH1F *MixedMM = (TH1F*)fin->Get("MixedMM");
  TH1F *MixedPM = (TH1F*)fin->Get("MixedPM");
  TH1F *LikeSign = (TH1F*)fin->Get("LikeSign");

  TFile *fout = new TFile("InvariantMassFit.root", "recreate");

  Int_t bin;
  Int_t Nbin = 1000;
  Int_t bin_histo_min = 0; //GeV
  Int_t bin_histo_max = 10; //GeV
  Int_t bin_histo_range = bin_histo_max - bin_histo_min; //GeV
  
  Int_t bkg_range_min = 80;
  Int_t bkg_range_max = 999;
  
  Int_t sig_range_min_1 = 200;
  Int_t sig_range_max_1 = 330;
  
  Int_t sig_range_min_2 = 350;
  Int_t sig_range_max_2 = 450;
  

  TH1F *bkg_histo = new TH1F("bkg_histo","bkg_histo",Nbin,bin_histo_min,bin_histo_max);
  for (bin = bkg_range_min; bin<bkg_range_max; bin++) {
    bkg_histo->SetBinContent(bin,LikeSign->GetBinContent(bin));
  }
  //bkg_histo = (TH1F*)LikeSign->Clone("bkg_histo");
  bkg_histo->Draw();
  for(bin = sig_range_min_1; bin < sig_range_max_2; bin++) {
    bkg_histo->SetBinContent(bin,0);
  }
  TF1 *FitBkg = new TF1("FitBkg","pol5",bkg_range_min,bkg_range_max);
  bkg_histo->Fit(FitBkg,"R","",bkg_range_min,bkg_range_max);

  //TH1F *sig_histo_1 = (TH1F*)LikeSign->Clone("sig_histo_1");
  TH1F *sig_histo_1 = new TH1F("sig_histo_1","sig_histo_1",Nbin,bin_histo_min,bin_histo_max);
  /*for (bin = bkg_range_min; bin<bkg_range_max; bin++) {
    sig_histo_1->SetBinContent(bin,LikeSign->GetBinContent(bin));
  }*/
  for (bin = bkg_range_min; bin<bkg_range_max; bin++) {
    Double_t bin2GeV = bin_histo_range*bin/Nbin;
    Double_t bkg = FitBkg->Eval(bin2GeV);
    //Double_t e_bkg = 0.;
    Double_t data = LikeSign->GetBinContent(bin);
    //Double_t e_data = LikeSign->GetBinError(bin);
    Double_t sig = data-bkg;
    //Double_t e_sig = sqrt(pow(e_data,2)+pow(e_bkg,2));

    sig_histo_1->SetBinContent(bin,sig);
    //sig_histo_1->SetBinError(bin,e_sig);
  }

  Double_t peak_mean_by_my_eye_1 = 3.1;
  Double_t peak_width_by_my_eye_1 = 0.2;

  TF1 *FitSig_1 = new TF1("FitSig_1","gaus",bkg_range_min,bkg_range_max);
  FitSig_1->SetParameters(sig_histo_1->GetBinContent(sig_histo_1->GetMaximumBin()),peak_mean_by_my_eye_1,peak_width_by_my_eye_1);
  sig_histo_1->Fit(FitSig_1,"","",sig_range_min_1,sig_range_max_1);

  //TH1F *sig_histo_2 = (TH1F*)LikeSign->Clone("sig_histo_2");
  TH1F *sig_histo_2 = new TH1F("sig_histo_2","sig_histo_2",Nbin,bin_histo_min,bin_histo_max);
  /*for (bin = bkg_range_min; bin<bkg_range_max; bin++) {
    sig_histo_2->SetBinContent(bin,LikeSign->GetBinContent(bin));
  }*/
  for (bin = bkg_range_min; bin<bkg_range_max; bin++) {
    Double_t bin2GeV = bin_histo_range*bin/Nbin;
    Double_t bkg = FitBkg->Eval(bin2GeV);
    //Double_t e_bkg = 0.;
    Double_t sig_1 = FitSig_1->Eval(bin2GeV);
    //Double_t e_sig_1 = 0.;
    Double_t data = LikeSign->GetBinContent(bin);
    //Double_t e_data = LikeSign->GetBinError(bin);
    Double_t sig = data-bkg-sig_1;
    //Double_t e_sig = sqrt(pow(e_data,2)+pow(e_bkg,2));

    sig_histo_2->SetBinContent(bin,sig);
    //sig_histo_2->SetBinError(bin,e_sig);
  }

  Double_t peak_mean_by_my_eye_2 = 3.7;
  Double_t peak_width_by_my_eye_2 = 0.6;

  TF1 *FitSig_2 = new TF1("FitSig_2","gaus",bkg_range_min,bkg_range_max);
  FitSig_2->SetParameters(sig_histo_2->GetBinContent(sig_histo_2->GetMaximumBin()),peak_mean_by_my_eye_2,peak_width_by_my_eye_2);
  sig_histo_2->Fit(FitSig_2,"","",sig_range_min_2,sig_range_max_2);
  
  sig_histo_1->Draw();
  sig_histo_2->Draw("same");

  SamePM->Write();
  SamePP->Write();
  SameMM->Write();
  MixedPM->Write();
  MixedPP->Write();
  MixedMM->Write();
  LikeSign->Write();
  sig_histo_1->Write();
  sig_histo_2->Write();
  bkg_histo->Write();
  fout->Close();

  return;
}
