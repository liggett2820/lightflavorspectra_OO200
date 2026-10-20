//____________________________________________________________________________________________________
// Step 2c for the EPD totnMIP centrality (O+O 200 GeV, trigger 860003): centrality cuts 0-100% from the
// best NBD x Glauber point, plus the 860003 efficiency (data/MC) vs totnMIP and per class.
//
// Replaces StNbdFitMaker::CalculateCentrality for totnMIP. That function only goes to 80% (we need 80-100%),
// switches data -> MC at a fixed 20% instead of at the fit range, and leaves unfilled entries of its
// centBin array uninitialized (garbage printed when a class is never reached).
//
// Input : RatioChi2Files/Ratio_npp..._k..._x..._eff....root of the best point (find_best_nbd_OO200.py prints it).
//         It holds "hRefMult" (data = corrected totnMIP, unit bins [n, n+1) from 0) and "hRefMultSim"
//         (MC, already normalized to data for totnMIP >= multCut by StNbdFitMaker::Fit).
// Method: N(>= c) = data above multCut (where 860003 is efficient) + MC between c and multCut;
//         N_tot = N(>= 0), the same hybrid (MC events with totnMIP = 0 included).
//         Each 5% boundary is the integer c whose N(>= c)/N_tot is closest to the target; the class is
//         totnMIP_corr >= c. Repeated with N_tot x 0.95 and x 1.05 (the usual +-5% total cross-section check).
// Check : data/MC vs totnMIP is the 860003 trigger+vertex+event-cut efficiency. It is compared with the
//         independent estimate  eff_MB(official StRefMultCorr totnmip weight) x R(860003/MB)/R_plateau,
//         R from macros/PlotTrigRatio_OO200.C. The official weight is for the official cuts and correction,
//         so agreement is a consistency check, not a proof.
// Usage : root -l -b -q 'CentCuts_OO200.C("RatioChi2Files/Ratio_npp24.000_k4.000_x0.050_eff0.000.root", 200)'
// Output: OO200_NBD_totnMIP_cuts.txt, OO200_NBD_totnMIP_fit.png, OO200_NBD_totnMIP_cuts.root
//____________________________________________________________________________________________________
#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"
#include "TSystem.h"
#include <cstdio>
#include <iostream>

// official StRefMultCorr("totnmip"), Param.h mParamStr_ref7 (2021:200:22130029,22144006:-30,30)
// lower edges of 75-80%, 70-75%, ..., 0-5% (peripheral -> central)
const int    kNOffCut = 16;
const double kOffCutPC[kNOffCut] = {37,46,55,66,77,91,105,122,141,162,186,213,242,276,315,367};
const double kOffP0 = 1.5443231276439915, kOffP1 = -0.021034900424299365, kOffP2 = -0.6284734564390636;
const double kOffNormStop = 355.;

double OffEffMB(double m){
  if(m >= kOffNormStop) return 1.0;
  return 1.0/(1.0 + kOffP0*TMath::Exp(kOffP1*m + kOffP2));
}

// official lower edge for the boundary at "pct" percent (5..80), -1 if none
double OffCut(int pct){
  if(pct < 5 || pct > 80) return -1;
  return kOffCutPC[(80 - pct)/5];
}

void CentCuts_OO200(const char* ratioFile = "RatioChi2Files/Ratio_best.root",
                    const double multCut = 200,
                    const char* trigRatioFile = "../../OO200_EPD_TrigRatio.root",
                    const char* tag = "OO200_NBD_totnMIP")
{
  gStyle->SetOptStat(0);
  TFile* in = TFile::Open(ratioFile);
  if(!in || in->IsZombie()){ std::cout << "cannot open " << ratioFile << std::endl; return; }
  TH1D* hD = (TH1D*) in->Get("hRefMult");
  TH1D* hM = (TH1D*) in->Get("hRefMultSim");
  if(!hD || !hM){ std::cout << "hRefMult / hRefMultSim missing in " << ratioFile << std::endl; return; }
  hD = (TH1D*) hD->Clone("hData"); hM = (TH1D*) hM->Clone("hMC");
  hD->SetDirectory(0); hM->SetDirectory(0);
  in->Close();

  const int nb = hD->GetNbinsX();
  if(hM->GetNbinsX() != nb || TMath::Abs(hD->GetBinWidth(1) - 1.0) > 1e-6 || TMath::Abs(hD->GetXaxis()->GetXmin()) > 1e-6){
    std::cout << "expected identical unit bins from 0 in data and MC" << std::endl; return;
  }
  const int mc = TMath::Nint(multCut);

  // cumulative N(>= c), c = 0..nb (c = nb: overflow only)
  double* nAbove = new double[nb+1];
  double* dAbove = new double[nb+1];
  double* mAbove = new double[nb+1];
  dAbove[nb] = hD->GetBinContent(nb+1); mAbove[nb] = hM->GetBinContent(nb+1);
  for(int c = nb-1; c >= 0; c--){
    dAbove[c] = dAbove[c+1] + hD->GetBinContent(c+1);
    mAbove[c] = mAbove[c+1] + hM->GetBinContent(c+1);
  }
  for(int c = 0; c <= nb; c++){
    nAbove[c] = (c >= mc) ? dAbove[c] : dAbove[mc] + (mAbove[c] - mAbove[mc]);
  }
  // N_tot is the same hybrid as N(>= c): data above multCut + MC below. Using the MC total instead would put
  // any MC/data difference above multCut (shape mismatch, MC beyond the data range) into the most peripheral class.
  const double nTot  = nAbove[0];
  const double nData = dAbove[0];
  const double fitMatch = (mAbove[mc] > 0) ? dAbove[mc]/mAbove[mc] : 0;

  // 5% boundaries 5..95 for the three total cross-section scales
  const int nB = 19;
  int    cut[3][nB];
  double frac[3][nB];
  const double scale[3] = {1.0, 0.95, 1.05};
  for(int s = 0; s < 3; s++){
    for(int j = 0; j < nB; j++){
      const double target = 0.05*(j+1);
      double best = 1e9; int bc = 0;
      for(int c = 0; c < nb; c++){
        const double f = nAbove[c]/(nTot*scale[s]);
        if(TMath::Abs(f - target) < best){ best = TMath::Abs(f - target); bc = c; }
      }
      cut[s][j] = bc; frac[s][j] = nAbove[bc]/(nTot*scale[s]);
    }
  }

  // 860003 / MB relative efficiency from PlotTrigRatio_OO200.C (optional)
  TGraphAsymmErrors* gR = 0;
  double rPlat = 0;
  TFile* fr = TFile::Open(trigRatioFile);
  if(fr && !fr->IsZombie()){
    gR = (TGraphAsymmErrors*) fr->Get("gRatio_860003_over_mb");
    if(gR){
      gR = (TGraphAsymmErrors*) gR->Clone("gR"); 
      int n = 0;
      for(int i = 0; i < gR->GetN(); i++){ if(gR->GetX()[i] >= 250 && gR->GetX()[i] < 400){ rPlat += gR->GetY()[i]; n++; } }
      if(n > 0) rPlat /= n; else gR = 0;
    }
    fr->Close();
  } else {
    std::cout << "no " << trigRatioFile << ": skipping the official-weight x ratio prediction" << std::endl;
  }

  // text output
  FILE* out = fopen(Form("%s_cuts.txt", tag), "w");
  FILE* fps[2] = {stdout, out};
  for(int k = 0; k < 2; k++){
    FILE* fp = fps[k];
    fprintf(fp, "EPD totnMIP centrality from %s\n", ratioFile);
    fprintf(fp, "data above multCut = %d, MC below. N_tot (hybrid) = %.4g, N_data(860003) = %.4g, overall 860003 efficiency = %.4f\n",
            mc, nTot, nData, nData/nTot);
    fprintf(fp, "MC events with totnMIP = 0: %.3f%% of N_tot\n", 100.*hM->GetBinContent(1)/nTot);
    fprintf(fp, "data/MC for totnMIP >= %d (should be ~1 for a good fit): %.4f\n", mc, fitMatch);
    if(gR) fprintf(fp, "860003/MB plateau (totnMIP 250-400) = %.4f\n", rPlat);
    fprintf(fp, "class = totnMIP_corr >= lower edge. Columns: class, edge, achieved %%, edge(-5%% xsec), edge(+5%% xsec), "
                "official edge, N_data, N_MC, data/MC%s\n", gR ? ", predicted (effMB_off x R/R_plat)" : "");
    for(int j = 0; j <= nB; j++){
      const int lo = (j < nB) ? cut[0][j] : 0;
      const int hi = (j == 0) ? nb : cut[0][j-1];
      double nd = 0, nm = 0, np = 0;
      for(int c = lo; c < hi; c++){
        nd += hD->GetBinContent(c+1); nm += hM->GetBinContent(c+1);
        if(gR){
          double r = gR->Eval(c + 0.5); if(c + 0.5 > 400) r = rPlat;
          np += hM->GetBinContent(c+1)*OffEffMB(c + 0.5)*r/rPlat;
        }
      }
      if(j == 0){ nd += hD->GetBinContent(nb+1); nm += hM->GetBinContent(nb+1); np += gR ? hM->GetBinContent(nb+1) : 0; }
      const int pct = 5*(j+1);
      const double off = OffCut(pct);
      if(j < nB)
        fprintf(fp, "%3d-%3d%%  %4d  %6.2f  %4d  %4d  %4s  %10.4g  %10.4g  %6.4f", 5*j, pct, lo, 100.*frac[0][j],
                cut[1][j], cut[2][j], off > 0 ? Form("%.0f", off) : "-", nd, nm, nm > 0 ? nd/nm : 0.);
      else
        fprintf(fp, "%3d-%3d%%  %4d  %6s  %4s  %4s  %4s  %10.4g  %10.4g  %6.4f", 5*j, pct, 0, "100", "-", "-", "-",
                nd, nm, nm > 0 ? nd/nm : 0.);
      if(gR) fprintf(fp, "  %6.4f", nm > 0 ? np/nm : 0.);
      fprintf(fp, "\n");
    }
    fprintf(fp, "edges, central -> peripheral (5..95%%): ");
    for(int j = 0; j < nB; j++) fprintf(fp, "%d%s", cut[0][j], j < nB-1 ? "," : "\n");
  }
  fclose(out);

  // plot: data and MC, data/MC with prediction and class edges
  TCanvas* cv = new TCanvas("cv", "", 900, 900);
  TPad* p1 = new TPad("p1", "", 0, 0.42, 1, 1); p1->SetBottomMargin(0.02); p1->SetLogy(); p1->Draw();
  TPad* p2 = new TPad("p2", "", 0, 0, 1, 0.42);  p2->SetTopMargin(0.02); p2->SetBottomMargin(0.22); p2->Draw();
  const double xmax = 650;
  p1->cd();
  hD->GetXaxis()->SetRangeUser(0, xmax); hD->SetLineColor(kBlack); hD->SetTitle("");
  hD->GetYaxis()->SetTitle("Events"); hD->GetXaxis()->SetLabelSize(0);
  hD->SetMinimum(TMath::Max(1.0, 1e-7*hD->GetMaximum()));
  hD->Draw("hist"); hM->SetLineColor(kRed+1); hM->Draw("hist same");
  TLine fl; fl.SetLineStyle(2); fl.SetLineColor(kBlue); fl.DrawLine(mc, hD->GetMinimum(), mc, hD->GetMaximum());
  TLegend* lg = new TLegend(0.55, 0.70, 0.88, 0.88); lg->SetBorderSize(0);
  lg->AddEntry(hD, "data, 860003 (corrected totnMIP)", "l"); lg->AddEntry(hM, "Glauber #times NBD", "l");
  lg->AddEntry(&fl, Form("fit range #geq %d", mc), "l"); lg->Draw();
  p2->cd();
  TH1D* hR = (TH1D*) hD->Clone("hDataOverMC"); hR->Divide(hM);
  hR->SetYTitle("data / MC"); hR->SetXTitle("totnMIP_{corr}");
  hR->GetXaxis()->SetLabelSize(0.07); hR->GetXaxis()->SetTitleSize(0.08);
  hR->GetYaxis()->SetLabelSize(0.07); hR->GetYaxis()->SetTitleSize(0.08); hR->GetYaxis()->SetTitleOffset(0.55);
  hR->SetMinimum(0.0); hR->SetMaximum(1.3); hR->SetMarkerStyle(20); hR->SetMarkerSize(0.4);
  hR->Draw("p");
  TGraph* gP = 0;
  if(gR){
    gP = new TGraph(); gP->SetName("gPredEff");
    for(int c = 0; c < (int)xmax; c++){
      double r = gR->Eval(c + 0.5); if(c + 0.5 > 400) r = rPlat;
      gP->SetPoint(gP->GetN(), c + 0.5, OffEffMB(c + 0.5)*r/rPlat);
    }
    gP->SetLineColor(kGreen+2); gP->SetLineWidth(2); gP->Draw("l same");
  }
  TLine one; one.SetLineStyle(2); one.DrawLine(0, 1, xmax, 1);
  TLine ed; ed.SetLineColor(kGray+1); ed.SetLineStyle(3);
  TLatex lab; lab.SetTextSize(0.05); lab.SetTextAlign(21);
  for(int j = 1; j < nB; j += 2){
    ed.DrawLine(cut[0][j], 0, cut[0][j], 1.3);
    lab.DrawLatex(cut[0][j], ((j/2) % 2 == 0) ? 0.05 : 0.13, Form("%d%%", 5*(j+1)));
  }
  if(gP){
    TLegend* l2 = new TLegend(0.45, 0.28, 0.88, 0.42); l2->SetBorderSize(0); l2->SetTextSize(0.05);
    l2->AddEntry(hR, "data / MC", "p");
    l2->AddEntry(gP, "official MB eff. #times R(860003/MB)/R_{plateau}", "l"); l2->Draw();
  }
  cv->SaveAs(Form("%s_fit.png", tag));

  TFile* fo = TFile::Open(Form("%s_cuts.root", tag), "RECREATE");
  hD->Write(); hM->Write(); hR->Write(); if(gP) gP->Write();
  fo->Close();
  delete[] nAbove; delete[] dAbove; delete[] mAbove;
}
