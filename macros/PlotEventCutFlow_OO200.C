// PlotEventCutFlow_OO200.C -- lightflavorspectra_OO200
// Plots the (hadd-merged) output of MakeEventCutFlow_OO200.C. Plain ROOT, no STAR libraries:
//   root -l -b -q 'macros/PlotEventCutFlow_OO200.C("cutflow_OO200_all.root")'
// Writes OO200_CutFlow.txt, OO200_Vz_before_after.png, OO200_TPC_vs_Vz.png and OO200_Luminosity.png.

#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <utility>
#include "TFile.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TBox.h"
#include "TStyle.h"
#include "TLatex.h"

void PlotEventCutFlow_OO200(const char* inName = "cutflow_OO200_all.root"){
  gStyle->SetOptStat(0); gStyle->SetPadLeftMargin(0.13); gStyle->SetPadBottomMargin(0.12);
  TFile* f = TFile::Open(inName);
  if(!f || f->IsZombie()){ std::cout << "cannot open " << inName << std::endl; return; }

  // ---- cut flow table
  TH1D* flow = (TH1D*)f->Get("hCutFlow");
  std::ofstream txt("OO200_CutFlow.txt");
  txt << "step\tevents\tfraction_of_previous\tfraction_of_trigger\n";
  for(int i=1;i<=flow->GetNbinsX();i++){
    double n = flow->GetBinContent(i), prev = (i>1) ? flow->GetBinContent(i-1) : n, trg = flow->GetBinContent(2);
    txt << flow->GetXaxis()->GetBinLabel(i) << "\t" << (long long)n << "\t" << (prev>0? n/prev:0) << "\t" << (i>=2 && trg>0 ? n/trg : 0) << "\n";
    std::cout << Form("%-32s %14.0f", flow->GetXaxis()->GetBinLabel(i), n) << std::endl;
  }
  // all trigger IDs, most frequent first
  TTree* tt = (TTree*)f->Get("triggers");
  if(tt){
    unsigned int id; double nev; tt->SetBranchAddress("trigId",&id); tt->SetBranchAddress("nEvents",&nev);
    std::vector<std::pair<double,unsigned int> > v;
    std::map<unsigned int,double> sum;
    for(Long64_t i=0;i<tt->GetEntries();i++){ tt->GetEntry(i); sum[id]+=nev; }   // hadd keeps one row per part
    for(std::map<unsigned int,double>::iterator it=sum.begin(); it!=sum.end(); ++it) v.push_back(std::make_pair(it->second,it->first));
    std::sort(v.rbegin(), v.rend());
    txt << "\ntrigger_id\tevents_with_this_id\n";
    for(size_t i=0;i<v.size();i++) txt << v[i].second << "\t" << (long long)v[i].first << "\n";
  }
  txt.close();

  // ---- Vz before / after
  TH1D* vzT = (TH1D*)f->Get("hVz_trigger"); TH1D* vzF = (TH1D*)f->Get("hVz_final");
  TCanvas* c1 = new TCanvas("c1","",1200,700); c1->SetLogy();
  vzT->SetLineColor(kAzure+3); vzT->SetLineWidth(2); vzT->SetTitle(""); vzT->GetXaxis()->SetRangeUser(-150,150);
  vzT->SetMinimum(0.5); vzT->Draw("hist");
  vzF->SetLineColor(kOrange+7); vzF->SetFillColor(kOrange+7); vzF->SetFillStyle(1001); vzF->Draw("hist same");
  double ymax = vzT->GetMaximum()*3;
  TLine* l1 = new TLine(-2,0.5,-2,ymax); TLine* l2 = new TLine(2,0.5,2,ymax);
  l1->SetLineStyle(2); l2->SetLineStyle(2); l1->Draw(); l2->Draw();
  TLegend* lg = new TLegend(0.58,0.74,0.89,0.89); lg->SetBorderSize(0);
  lg->AddEntry(vzT, Form("Trigger 860003: %.1f M", vzT->GetEntries()/1e6), "l");
  lg->AddEntry(vzF, Form("After all cuts: %.1f M", vzF->GetEntries()/1e6), "f");
  lg->Draw(); c1->SaveAs("OO200_Vz_before_after.png");

  // ---- TPC performance vs Vz
  TCanvas* c2 = new TCanvas("c2","",1400,600); c2->Divide(2,1);
  const char* pn[2] = {"pRefMult_vs_Vz","pNHitsFit_vs_Vz"};
  for(int k=0;k<2;k++){
    c2->cd(k+1); TProfile* p = (TProfile*)f->Get(pn[k]); if(!p) continue;
    p->SetMarkerStyle(20); p->SetMarkerColor(kAzure+3); p->SetLineColor(kAzure+3); p->GetXaxis()->SetRangeUser(-150,150);
    p->Draw();
    double lo = p->GetMinimum(), hi = p->GetMaximum();
    TBox* b = new TBox(-2, lo, 2, hi); b->SetFillColorAlpha(kOrange+7,0.35); b->Draw();
    p->Draw("same");
  }
  c2->SaveAs("OO200_TPC_vs_Vz.png");

  // ---- luminosity
  TCanvas* c3 = new TCanvas("c3","",1800,600); c3->Divide(3,1);
  c3->cd(1); TH1D* z = (TH1D*)f->Get("hZDCx_final"); z->SetLineColor(kAzure+3); z->SetLineWidth(2); z->Draw("hist");
  c3->cd(2);
  TTree* tr = (TTree*)f->Get("runs");
  if(tr){
    int run; double nEv, zs, nAll = 0, nTrig = 0;
    tr->SetBranchAddress("run",&run); tr->SetBranchAddress("nEvents",&nEv); tr->SetBranchAddress("zdcSum_kHz",&zs);
    const bool hasAll = tr->GetBranch("nAll") != 0;   // older outputs only have post-cut rows
    if(hasAll){ tr->SetBranchAddress("nAll",&nAll); tr->SetBranchAddress("nTrig",&nTrig); }
    std::map<int,std::pair<double,double> > m;        // hadd keeps one row per job: sum by run
    std::map<int,std::pair<double,double> > mAT;      // run -> (all, 860003)
    for(Long64_t i=0;i<tr->GetEntries();i++){
      tr->GetEntry(i); m[run].first+=nEv; m[run].second+=zs;
      if(hasAll){ mAT[run].first+=nAll; mAT[run].second+=nTrig; }
    }
    TGraph* g = new TGraph(); int k=0, nRunTrig=0; double sumAll=0, sumTrig=0, sumAllTrigRuns=0;
    std::ofstream rt("OO200_RunList.txt"); rt << "run\tall_events\ttrigger_860003\tevents_after_cuts\tmean_ZDCx_kHz\n";
    for(std::map<int,std::pair<double,double> >::iterator it=m.begin(); it!=m.end(); ++it){
      double mean = it->second.first>0 ? it->second.second/it->second.first : 0;
      double a = hasAll ? mAT[it->first].first : -1, t = hasAll ? mAT[it->first].second : -1;
      rt << it->first << "\t" << (long long)a << "\t" << (long long)t << "\t" << (long long)it->second.first << "\t" << mean << "\n";
      if(hasAll){ sumAll += a; sumTrig += t; if(t > 0){ nRunTrig++; sumAllTrigRuns += a; } }
      if(it->second.first > 0){ g->SetPoint(k, k, mean); k++; }
    }
    rt.close();
    if(hasAll){
      std::cout << Form("Runs read: %d; runs with trigger 860003: %d", (int)m.size(), nRunTrig) << std::endl;
      std::cout << Form("All events: %.0f; in runs with 860003: %.0f; with 860003: %.0f", sumAll, sumAllTrigRuns, sumTrig) << std::endl;
      txt.open("OO200_CutFlow.txt", std::ios::app);
      txt << "\nruns_read\t" << m.size() << "\nruns_with_860003\t" << nRunTrig
          << "\nall_events_in_runs_with_860003\t" << (long long)sumAllTrigRuns << "\n";
      txt.close();
    }
    g->SetTitle(Form("%d runs;run index;#LTZDC coincidence rate#GT (kHz)", k));
    g->SetMarkerStyle(20); g->SetMarkerSize(0.6); g->SetMarkerColor(kAzure+3); g->Draw("AP");
  }
  c3->cd(3); TProfile* pz = (TProfile*)f->Get("pRefMult_vs_ZDCx"); pz->SetMarkerStyle(20); pz->SetMarkerColor(kAzure+3); pz->Draw();
  c3->SaveAs("OO200_Luminosity.png");
  std::cout << "Wrote OO200_CutFlow.txt, OO200_RunList.txt and the three PNGs" << std::endl;
}
