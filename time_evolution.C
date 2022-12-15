#include "TEfficiency.h"
#include "TFile.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TTree.h"
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TProfile.h>
#include <TStyle.h>
void time_evolution(){

  //衝突>>ハード散乱>>熱平衡化(QGP相に相転移,カイラル相転移)>>ハドロン化が始まる>>ケミカルフリーズアウト,ハドロン相に相転移>>カイラル相転移>>キネマチックフリーズアウト

  Double_t t_collision = 0.0;  //ビームが衝突した時刻
  Double_t t_netsu_heikou = t_collision + 1.0; //熱平衡化が始まる時刻
  Double_t t_chemi = 40.0; //ケミカルフリーズアウトする時刻
  Double_t t_kine = 60.0;  //キネマチックフリーズアウトする時刻

  Double_t t_hadronize = 7.0/200; //ハドロン化する時刻
  Double_t t_chiral = 16.0/200;  //カイラル対称性が自発的に破れ始める温度まで下がってきたときの時刻

  Double_t Atom = 207; //鉛
  Double_t v_QGP = 0.66;  //QGPのtransverse方向への空間発展の速度 光速c=1
  Double_t r_Pb = 1.3*TMath::Power(Atom,-1/3);  //鉛原子核の半径

  Double_t r_hadronize = r_Pb + v_QGP*(t_hadronize - t_collision); //ハドロン化が始まったときの系のtransverse方向の半径
  Double_t r_chiral = r_Pb + v_QGP*(t_chiral - t_collision); //カイラル対称性が自発的に破れ始めるときの系のtransverse方向の半径

  //寿命、生成点、生成時刻、pTを乱数振って決める。すると崩壊するときの時刻、座標が計算できて、その時刻での系のサイズと比較して中か外か判定する。中だったときのpTを粒子種ごとのヒストに詰めればいいってコト！？

  //ヒストグラム
  TH1F *a1_decay_inside_pT = new TH1F("a1_decay_inside_pT","pT distribution of a_{1} meson which decay inside;p_{T}[GeV/c];Entry",1000,0,10);
  TH1F *a1_generated_pT = new TH1F("a1_generated_pT","pT distribution of generated a1;p_{T}[GeV/c];Entry",1000,0,10);
  TH1F *a1_generated_time = new TH1F("a1_generated_time","when is a1 generated;t[fm/c];Entry",1000,0,100);
  TH1F *a1_generated_coordination = new TH1F("a1_generated_coordination","a1's generation point;r[fm];Entry",1000,-50,50);
  TH1F *a1_decay_coordination = new TH1F("a1_decay_coordination","a1's decay point;r[fm];Entry",1000,-100,100);
  TH1F *a1_generated_lifetime = new TH1F("a1_generated_lifetime","a1's lifetime distribution;t[fm/c];Entry",1000,0,1000);
  TH1F *a1_lifetime_distribution = new TH1F("a1_lifetime_distribution","rest a1's life time distribution;t[fm/c];Entry",1000,0,1000);

  TH1F *rho_decay_inside_pT = new TH1F("rho_decay_inside_pT","pT distribution of rho meson which decay inside;p_{T}[GeV/c];Entry",1000,0,10);
  TH1F *rho_generated_pT = new TH1F("rho_generated_pT","pT distribution of generated rho;p_{T}[GeV/c];Entry",1000,0,10);
  TH1F *rho_generated_time = new TH1F("rho_generated_time","when is rho generated;t[fm/c];Entry",1000,0,100);
  TH1F *rho_generated_coordination = new TH1F("rho_generated_coordination","rho's generation point;r[fm];Entry",1000,-50,50);
  TH1F *rho_decay_coordination = new TH1F("rho_decay_coordination","rho's decay point;r[fm];Entry",1000,-100,100);
  TH1F *rho_generated_lifetime = new TH1F("rho_generated_lifetime","rho's lifetime distribution;t[fm/c];Entry",1000,0,1000);
  TH1F *rho_lifetime_distribution = new TH1F("rho_lifetime_distribution","rest rho's life time distribution;t[fm/c];Entry",1000,0,1000);
  TH1F *rho_beta = new TH1F("rho_beta","rho's beta;#beta;Entry",250,0,1);
  TH1F *rho_inside_rate = new TH1F("rho_inside_rate","rho inside rate;p_{T}[GeV/c];decayed inside / all",1000,0,10);

  TH1F *omega_decay_inside_pT = new TH1F("omega_decay_inside_pT","pT distribution of omega meson which decay inside;p_{T}[GeV/c];Entry",1000,0,10);
  TH1F *omega_generated_pT = new TH1F("omega_generated_pT","pT distribution of generated omega;p_{T}[GeV/c];Entry",1000,0,10);
  TH1F *omega_generated_time = new TH1F("omega_generated_time","when is omega generated;t[fm/c];Entry",1000,0,100);
  TH1F *omega_generated_coordination = new TH1F("omega_generated_coordination","omega's generation point;r[fm];Entry",1000,-50,50);
  TH1F *omega_decay_coordination = new TH1F("omega_decay_coordination","omega's decay point;r[fm];Entry",1000,-100,100);
  TH1F *omega_generated_lifetime = new TH1F("omega_generated_lifetime","omega's lifetime distribution;t[fm/c];Entry",1000,0,1000);
  TH1F *omega_lifetime_distribution = new TH1F("omega_lifetime_distribution","rest omega's life time distribution;t[fm/c];Entry",1000,0,1000);
  TH1F *omega_beta = new TH1F("omega_beta","omega's beta;#beta;Entry",250,0,1);
  TH1F *omega_inside_rate = new TH1F("omega_inside_rate","omega inside rate;p_{T}[GeV/c];decayed inside / all",1000,0,10);

  TH1F *phi_decay_inside_pT = new TH1F("phi_decay_inside_pT","pT distribution of phi meson which decay inside;p_{T}[GeV/c];Entry",1000,0,10);
  TH1F *phi_generated_pT = new TH1F("phi_generated_pT","pT distribution of generated phi;p_{T}[GeV/c];Entry",1000,0,10);
  TH1F *phi_generated_time = new TH1F("phi_generated_time","when is phi generated;t[fm/c];Entry",1000,0,100);
  TH1F *phi_generated_coordination = new TH1F("phi_generated_coordination","phi's generation point;r[fm];Entry",1000,-50,50);
  TH1F *phi_decay_coordination = new TH1F("phi_decay_coordination","phi's decay point;r[fm];Entry",1000,-100,100);
  TH1F *phi_generated_lifetime = new TH1F("phi_generated_lifetime","phi's lifetime distribution;t[fm/c];Entry",1000,0,1000);
  TH1F *phi_lifetime_distribution = new TH1F("phi_lifetime_distribution","rest phi's life time distribution;t[fm/c];Entry",1000,0,1000);
  TH1F *phi_beta = new TH1F("phi_beta","phi's beta;#beta;Entry",250,0,1);
  TH1F *phi_inside_rate = new TH1F("phi_inside_rate","phi inside rate;p_{T}[GeV/c];decayed inside / all",1000,0,10);
  //出力ファイル
  std::string outfilename = "time_evolution.root";
  TFile outFile(outfilename.c_str(), "RECREATE");

  Int_t max = 1000000000;
/*
  //a1
  srand(time(NULL));  //乱数の種を初期化
  Int_t max = 1000000000;
  Double_t a1_pT, a1_lorentz_gamma, t_generated_a1, rRange_t_generated_a1, r_generated_a1, a1_lifetime, rRange_t_disappear_a1, r_decay_a1;
  for(Int_t ite=0;ite<max;ite++){
    Double_t x = (double)rand()/RAND_MAX*10;
    Double_t y = (double)rand()/RAND_MAX*2;
    Float_t p0, p1, p2, p3, p4, p5, p6;
    p0 = 1;
    p1 = 0.442572;
    p2 = 3.5986;
    p3 = 0.875139;
    p4 = 10.2394;
    p5 = 1.260;
    p6 = 0.770;
    Double_t a1_pT_dist = p0 * TMath::Sqrt(x*x + p5*p5 - p6*p6) * TMath::Power(x,2.) / TMath::Power(p1 + TMath::Power(x / p2, p3), p4);
    if(y<a1_pT_dist){
      a1_pT = x; //pT
      a1_lorentz_gamma = TMath::Sqrt(1.0 + a1_pT*a1_pT/(1.260*1.260));  //ローレンツガンマ
    }
    else{
      continue;
    }

    t_generated_a1 = (double)rand()/RAND_MAX*(t_chiral - t_hadronize) + t_hadronize; //生成時刻

    rRange_t_generated_a1 = r_Pb + v_QGP*t_generated_a1;  //生成時刻における系のサイズ

    if ((double)rand()/RAND_MAX<0.5){
      r_generated_a1 = (double)rand()/RAND_MAX*rRange_t_generated_a1; //生成点
    }
    else{
      r_generated_a1 = -1.0*(double)rand()/RAND_MAX*rRange_t_generated_a1;  //生成点
    }

    Double_t i = (double)rand()/RAND_MAX*1000;
    Double_t j = (double)rand()/RAND_MAX*2;
    Double_t a1_lifetime_dist = TMath::Exp(-i*0.25);  //a1の崩壊幅 = 250 to 600 MeV
    if(j<a1_lifetime_dist){
      a1_lifetime = i*a1_lorentz_gamma; //寿命
    }
    else{
      continue;
    }

    rRange_t_disappear_a1 = r_Pb + v_QGP*(t_generated_a1 + a1_lifetime);  //崩壊時の系のサイズ
    r_decay_a1 = a1_pT/(a1_lorentz_gamma*1.260)*(t_generated_a1 + a1_lifetime) + r_generated_a1;  //崩壊時の位置

    a1_generated_pT->Fill(a1_pT);
    a1_generated_time->Fill(t_generated_a1);
    a1_generated_coordination->Fill(r_generated_a1);
    a1_decay_coordination->Fill(r_decay_a1);
    a1_generated_lifetime->Fill(a1_lifetime);
    a1_lifetime_distribution->Fill(i);

    if (r_decay_a1 < rRange_t_disappear_a1){  //系の中で崩壊する場合
      a1_decay_inside_pT->Fill(a1_pT);
      std::cout<<"[a1][Inside] pT = "<<a1_pT<<" GeV/c : t_generated_a1 = "<<t_generated_a1<<" fm/c : r_generated_a1 = "<<r_generated_a1<<" fm: r_decay_a1 = "<<r_decay_a1<<" fm : lifetime = "<<a1_lifetime<<" fm/c"<<std::endl;
    }
    else{ //系の外で崩壊する場合
      std::cout<<"[a1][Outside] pT = "<<a1_pT<<" GeV/c : t_generated_a1 = "<<t_generated_a1<<" fm/c : r_generated_a1 = "<<r_generated_a1<<" fm : r_decay_a1 = "<<r_decay_a1<<" fm : lifetime = "<<a1_lifetime<<" fm/c"<<std::endl;
    }
  }

*/
  //Rho
  srand(time(NULL));  //乱数の種を初期化
  Double_t rho_pT, rho_lorentz_gamma, t_generated_rho, rRange_t_generated_rho, r_generated_rho, rho_lifetime, rRange_t_disappear_rho, r_decay_rho;
  for(Int_t ite=0;ite<max;ite++){
    Double_t x = (double)rand()/RAND_MAX*10;
    Double_t y = (double)rand()/RAND_MAX*2;
    Float_t p0,p1,p2,p3,p4;
  	p0 = 1.;
  	p1 = 0.512514;
  	p2 = 1.11739;
  	p3 = 1.88096;
    p4 = 2.43008;
    Double_t rho_pT_dist = p0 * TMath::Power(x,2.) / TMath::Power(p1 + TMath::Power(x / p2, p3), p4);
    if(y<rho_pT_dist){
      rho_pT = x; //pT
      rho_lorentz_gamma = TMath::Sqrt(1.0 + rho_pT*rho_pT/(0.770*0.770));  //ローレンツガンマ
    }
    else{
      continue;
    }

    t_generated_rho = (double)rand()/RAND_MAX*(t_chiral - t_hadronize) + t_hadronize; //生成時刻

    rRange_t_generated_rho = r_Pb + v_QGP*t_generated_rho;  //生成時刻における系のサイズ

    if ((double)rand()/RAND_MAX<0.5){
      r_generated_rho = (double)rand()/RAND_MAX*rRange_t_generated_rho; //生成点
    }
    else{
      r_generated_rho = -1.0*(double)rand()/RAND_MAX*rRange_t_generated_rho;  //生成点
    }

    Double_t i = (double)rand()/RAND_MAX*1000;
    Double_t j = (double)rand()/RAND_MAX*2;
    Double_t rho_lifetime_dist = TMath::Exp(-i*0.149);  //rhoの崩壊幅 = 149 MeV
    if(j<rho_lifetime_dist){
      rho_lifetime = i*rho_lorentz_gamma; //寿命
    }
    else{
      continue;
    }

    rRange_t_disappear_rho = r_Pb + v_QGP*(t_generated_rho + rho_lifetime);  //崩壊時の系のサイズ
    r_decay_rho = rho_pT/(rho_lorentz_gamma*0.770)*(t_generated_rho + rho_lifetime) + r_generated_rho;  //崩壊時の位置

    rho_generated_pT->Fill(rho_pT);
    rho_generated_time->Fill(t_generated_rho);
    rho_generated_coordination->Fill(r_generated_rho);
    rho_decay_coordination->Fill(r_decay_rho);
    rho_generated_lifetime->Fill(rho_lifetime);
    rho_lifetime_distribution->Fill(i);

    if (r_decay_rho < rRange_t_disappear_rho){  //系の中で崩壊する場合
      rho_decay_inside_pT->Fill(rho_pT);
      Double_t beta = rho_pT/(rho_lorentz_gamma*0.770);
      rho_beta->Fill(beta);
      std::cout<<"[rho][Inside] pT = "<<rho_pT<<" GeV/c : t_generated_rho = "<<t_generated_rho<<" fm/c : r_generated_rho = "<<r_generated_rho<<" fm: r_decay_rho = "<<r_decay_rho<<" fm : lifetime = "<<rho_lifetime<<" fm/c"<<std::endl;
    }
    else{ //系の外で崩壊する場合
      std::cout<<"[rho][Outside] pT = "<<rho_pT<<" GeV/c : t_generated_rho = "<<t_generated_rho<<" fm/c : r_generated_rho = "<<r_generated_rho<<" fm : r_decay_rho = "<<r_decay_rho<<" fm : lifetime = "<<rho_lifetime<<" fm/c"<<std::endl;
    }
  }

  //omega
  srand(time(NULL));  //乱数の種を初期化
  Double_t omega_pT, omega_lorentz_gamma, t_generated_omega, rRange_t_generated_omega, r_generated_omega, omega_lifetime, rRange_t_disappear_omega, r_decay_omega;
  for(Int_t ite=0;ite<max;ite++){
    Double_t x = (double)rand()/RAND_MAX*10;
    Double_t y = (double)rand()/RAND_MAX*2;
    Float_t p0,p1,p2,p3,p4;
  	p0 = 1.;
  	p1 = 0.826671;
  	p2 = 0.975843;
  	p3 = 1.77957;
    p4 = 2.85792;
    Double_t omega_pT_dist = p0 * TMath::Power(x,2.) / TMath::Power(p1 + TMath::Power(x / p2, p3), p4);
    if(y<omega_pT_dist){
      omega_pT = x; //pT
      omega_lorentz_gamma = TMath::Sqrt(1.0 + omega_pT*omega_pT/(0.782*0.782));  //ローレンツガンマ
    }
    else{
      continue;
    }

    t_generated_omega = (double)rand()/RAND_MAX*(t_chiral - t_hadronize) + t_hadronize; //生成時刻

    rRange_t_generated_omega = r_Pb + v_QGP*t_generated_omega;  //生成時刻における系のサイズ

    if ((double)rand()/RAND_MAX<0.5){
      r_generated_omega = (double)rand()/RAND_MAX*rRange_t_generated_omega; //生成点
    }
    else{
      r_generated_omega = -1.0*(double)rand()/RAND_MAX*rRange_t_generated_omega;  //生成点
    }

    Double_t i = (double)rand()/RAND_MAX*1000;
    Double_t j = (double)rand()/RAND_MAX*2;
    Double_t omega_lifetime_dist = TMath::Exp(-i*0.00868);  //omegaの崩壊幅 = 8.68 MeV
    if(j<omega_lifetime_dist){
      omega_lifetime = i*omega_lorentz_gamma; //寿命
    }
    else{
      continue;
    }

    rRange_t_disappear_omega = r_Pb + v_QGP*(t_generated_omega + omega_lifetime);  //崩壊時の系のサイズ
    r_decay_omega = omega_pT/(omega_lorentz_gamma*0.782)*(t_generated_omega + omega_lifetime) + r_generated_omega;  //崩壊時の位置

    omega_generated_pT->Fill(omega_pT);
    omega_generated_time->Fill(t_generated_omega);
    omega_generated_coordination->Fill(r_generated_omega);
    omega_decay_coordination->Fill(r_decay_omega);
    omega_generated_lifetime->Fill(omega_lifetime);
    omega_lifetime_distribution->Fill(i);

    if (r_decay_omega < rRange_t_disappear_omega){  //系の中で崩壊する場合
      omega_decay_inside_pT->Fill(omega_pT);
      Double_t beta = omega_pT/(omega_lorentz_gamma*0.782);
      omega_beta->Fill(beta);
      std::cout<<"[omega][Inside] pT = "<<omega_pT<<" GeV/c : t_generated_omega = "<<t_generated_omega<<" fm/c : r_generated_omega = "<<r_generated_omega<<" fm: r_decay_omega = "<<r_decay_omega<<" fm : lifetime = "<<omega_lifetime<<" fm/c"<<std::endl;
    }
    else{ //系の外で崩壊する場合
      std::cout<<"[omega][Outside] pT = "<<omega_pT<<" GeV/c : t_generated_omega = "<<t_generated_omega<<" fm/c : r_generated_omega = "<<r_generated_omega<<" fm : r_decay_omega = "<<r_decay_omega<<" fm : lifetime = "<<omega_lifetime<<" fm/c"<<std::endl;
    }
  }

    //phi
    srand(time(NULL));  //乱数の種を初期化
    Double_t phi_pT, phi_lorentz_gamma, t_generated_phi, rRange_t_generated_phi, r_generated_phi, phi_lifetime, rRange_t_disappear_phi, r_decay_phi;
    for(Int_t ite=0;ite<max;ite++){
      Double_t x = (double)rand()/RAND_MAX*10;
      Double_t y = (double)rand()/RAND_MAX*2;
      Float_t p0,p1,p2,p3,p4;
    	p0 = 1.;
    	p1 = 0.907455;
    	p2 = 0.945388;
    	p3 = 1.71594;
      p4 = 2.50996;
      Double_t phi_pT_dist = p0 * TMath::Power(x,2.) / TMath::Power(p1 + TMath::Power(x / p2, p3), p4);
      if(y<phi_pT_dist){
        phi_pT = x; //pT
        phi_lorentz_gamma = TMath::Sqrt(1.0 + phi_pT*phi_pT/(1.020*1.020));  //ローレンツガンマ
      }
      else{
        continue;
      }

      t_generated_phi = (double)rand()/RAND_MAX*(t_chiral - t_hadronize) + t_hadronize; //生成時刻

      rRange_t_generated_phi = r_Pb + v_QGP*t_generated_phi;  //生成時刻における系のサイズ

      if ((double)rand()/RAND_MAX<0.5){
        r_generated_phi = (double)rand()/RAND_MAX*rRange_t_generated_phi; //生成点
      }
      else{
        r_generated_phi = -1.0*(double)rand()/RAND_MAX*rRange_t_generated_phi;  //生成点
      }

      Double_t i = (double)rand()/RAND_MAX*1000;
      Double_t j = (double)rand()/RAND_MAX*2;
      Double_t phi_lifetime_dist = TMath::Exp(-i*0.004249);  //phiの崩壊幅 = 4.249 MeV
      if(j<phi_lifetime_dist){
        phi_lifetime = i*phi_lorentz_gamma; //寿命
      }
      else{
        continue;
      }

      rRange_t_disappear_phi = r_Pb + v_QGP*(t_generated_phi + phi_lifetime);  //崩壊時の系のサイズ
      r_decay_phi = phi_pT/(phi_lorentz_gamma*1.020)*(t_generated_phi + phi_lifetime) + r_generated_phi;  //崩壊時の位置

      phi_generated_pT->Fill(phi_pT);
      phi_generated_time->Fill(t_generated_phi);
      phi_generated_coordination->Fill(r_generated_phi);
      phi_decay_coordination->Fill(r_decay_phi);
      phi_generated_lifetime->Fill(phi_lifetime);
      phi_lifetime_distribution->Fill(i);

      if (r_decay_phi < rRange_t_disappear_phi){  //系の中で崩壊する場合
        phi_decay_inside_pT->Fill(phi_pT);
        Double_t beta = phi_pT/(phi_lorentz_gamma*1.020);
        phi_beta->Fill(beta);
        std::cout<<"[phi][Inside] pT = "<<phi_pT<<" GeV/c : t_generated_phi = "<<t_generated_phi<<" fm/c : r_generated_phi = "<<r_generated_phi<<" fm: r_decay_phi = "<<r_decay_phi<<" fm : lifetime = "<<phi_lifetime<<" fm/c"<<std::endl;
      }
      else{ //系の外で崩壊する場合
        std::cout<<"[phi][Outside] pT = "<<phi_pT<<" GeV/c : t_generated_phi = "<<t_generated_phi<<" fm/c : r_generated_phi = "<<r_generated_phi<<" fm : r_decay_phi = "<<r_decay_phi<<" fm : lifetime = "<<phi_lifetime<<" fm/c"<<std::endl;
      }
    }
/*
  a1_generated_pT->Write();
  a1_generated_time->Write();
  a1_generated_coordination->Write();
  a1_decay_coordination->Write();
  a1_generated_lifetime->Write();
  a1_lifetime_distribution->Write();
  a1_decay_inside_pT->Write();
*/
  rho_generated_pT->Write();
  rho_generated_time->Write();
  rho_generated_coordination->Write();
  rho_decay_coordination->Write();
  rho_generated_lifetime->Write();
  rho_lifetime_distribution->Write();
  rho_decay_inside_pT->Write();
  rho_beta->Write();
  omega_generated_pT->Write();
  omega_generated_time->Write();
  omega_generated_coordination->Write();
  omega_decay_coordination->Write();
  omega_generated_lifetime->Write();
  omega_lifetime_distribution->Write();
  omega_decay_inside_pT->Write();
  omega_beta->Write();
  phi_generated_pT->Write();
  phi_generated_time->Write();
  phi_generated_coordination->Write();
  phi_decay_coordination->Write();
  phi_generated_lifetime->Write();
  phi_lifetime_distribution->Write();
  phi_decay_inside_pT->Write();
  phi_beta->Write();
  rho_inside_rate->Divide(rho_decay_inside_pT,rho_generated_pT);
  omega_inside_rate->Divide(omega_decay_inside_pT,omega_generated_pT);
  phi_inside_rate->Divide(phi_decay_inside_pT,phi_generated_pT);
  rho_inside_rate->Write();
  omega_inside_rate->Write();
  phi_inside_rate->Write();
  outFile.Close();
}
