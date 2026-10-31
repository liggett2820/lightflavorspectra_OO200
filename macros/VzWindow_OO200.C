// VzWindow_OO200.C -- lightflavorspectra_OO200
//
// Is |Vz| <= 30 cm justified by the vertex distribution? Reads hVz_trigger (minbias-bbc_etof,
// before the Vz cut; 0.5 cm bins over +-200 cm) from the merged MakeEventCutFlow_OO200 output and
// prints, for windows |Vz| <= W:
//   N      events inside the window
//   frac   N / all triggered events
//   flat   events from the flat component (level from the sidebands 90 < |Vz| < 150 cm) inside W
//   purity 1 - flat/N
//   shell  for the 10 cm added on each side going from W-10 to W: events added and their flat share
// and the RMS of the peak (flat level subtracted, |Vz| < 60 cm).
// The flat level is an estimate: it assumes the plateau seen at large |Vz| continues under the peak.
//
// Usage: root -l -b -q 'macros/VzWindow_OO200.C("cutflow_OO200.root")'
// Written for ROOT 5 (C++98, no ternaries).

#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TString.h"

double sumRange(TH1D* h, double lo, double hi, int& nb){
  int b1 = h->FindBin(lo + 1e-6);
  int b2 = h->FindBin(hi - 1e-6);
  nb = b2 - b1 + 1;
  return h->Integral(b1, b2);
}

void VzWindow_OO200(const char* inName = "cutflow_OO200.root"){
  TFile* f = TFile::Open(inName);
  if(!f || f->IsZombie()){ std::cout << "cannot open " << inName << std::endl; return; }
  TH1D* h = (TH1D*) f->Get("hVz_trigger");
  if(!h){ std::cout << "no hVz_trigger in " << inName << std::endl; return; }

  int nb = 0;
  double all = h->Integral(0, h->GetNbinsX() + 1);
  int nL = 0, nR = 0;
  double sL = sumRange(h, -150., -90., nL);
  double sR = sumRange(h, 90., 150., nR);
  double flatPerBin = (sL + sR) / (nL + nR);
  std::cout << Form("All triggered events (incl. under/overflow): %.4g", all) << std::endl;
  std::cout << Form("Flat level from 90<|Vz|<150 cm: %.4g per %.2g cm bin (left %.4g, right %.4g)",
                    flatPerBin, h->GetBinWidth(1), sL / nL, sR / nR) << std::endl;

  // peak RMS with the flat level subtracted, |Vz| < 60 cm
  double s0 = 0, s1 = 0, s2 = 0;
  for(int b = h->FindBin(-60. + 1e-6); b <= h->FindBin(60. - 1e-6); b++){
    double w = h->GetBinContent(b) - flatPerBin;
    if(w < 0) w = 0;
    double x = h->GetBinCenter(b);
    s0 += w; s1 += w * x; s2 += w * x * x;
  }
  double mean = s1 / s0;
  std::cout << Form("Peak (flat subtracted, |Vz|<60): mean %.2f cm, RMS %.2f cm", mean, TMath::Sqrt(s2 / s0 - mean * mean)) << std::endl;

  std::cout << "\n    W     N         frac    flat      purity   | shell W-10..W: added   flat share" << std::endl;
  double prevN = 0, prevB = 0;
  double Ws[8] = {2., 10., 20., 30., 40., 50., 60., 80.};
  for(int i = 0; i < 8; i++){
    double W = Ws[i];
    double N = sumRange(h, -W, W, nb);
    double B = flatPerBin * nb;
    TString shell = "";
    if(W >= 20.){
      double dN = N - prevN, dB = B - prevB;
      double share = 0; if(dN > 0) share = dB / dN;
      shell = Form("%10.4g   %5.1f%%", dN, 100. * share);
    }
    std::cout << Form("  %4.0f  %10.4g  %5.1f%%  %10.4g  %5.1f%%   | ", W, N, 100. * N / all, B, 100. * (1. - B / N)) << shell << std::endl;
    prevN = N; prevB = B;
  }
  std::cout << "\nThe flat share in each shell is the fraction of the events added there that belong to the plateau"
            << " (beam-gas, pile-up or misreconstructed vertices, not yet identified)." << std::endl;
}
