//______________________________________________________________________________
// NpartPerClass_OO200.C
//
// <Npart>, <Ncoll> and <Nspec> = 32 - <Npart> in the six EPD totnMIP centrality classes
// (0-5, 5-10, 10-20, 20-40, 40-80, 80-100%), for the nominal Glauber and each Glauber
// variation, each with its own best NBD fit and its own CentCuts edges.
//
// The model is generated exactly as StNbdFitMaker::Fit does it for the O+O fits:
//   (Npart, Ncoll) from hNcoll_Npart with GetRandom2, floored to integers (SetIntegerGlauberSampling),
//   events with Npart < 2 or Ncoll < 1 skipped;
//   multiplicity as in StNegativeBinomial::GetMultiplicity: sum of Nint((1-x)Npart/2 + x Ncoll)
//   draws from the NBD(npp, k) histogram (100 unit bins), then 2*multIdeal trials, each kept when
//   a uniform number on (0, 2) is <= 0.98 (1 - d multIdeal/540) (StGlauberUtilities::GetUniform2;
//   multiplicity-dependent mode, trigger bias 1).
// A class is totnMIP >= lower edge (and below the next class's edge), as in CentCuts_OO200.C.
//
// Run from centrality/star_glauber (needs ncoll_npart.root and ncoll_npart_<type>.root), in the
// SL7 container after starver:
//   root -l -b -q 'NpartPerClass_OO200.C+'
// Output: printed tables and NpartPerClass_OO200.txt
//______________________________________________________________________________
// rootcint (ROOT 5 dictionary generation for ACLiC) cannot parse the variation table below,
// so it only sees the function declaration; the compiler sees everything.
#if !defined(__CINT__) || !defined(__MAKECINT__)
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
#include <cstdio>

namespace {
  const Int_t kNClass = 6;
  const char* kClassName[kNClass] = { "0-5%", "5-10%", "10-20%", "20-40%", "40-80%", "80-100%" };
  const Double_t kA2 = 32.0; // nucleons in O+O

  struct Variation {
    const char* name;
    const char* glauberFile;
    Double_t npp, k, x, d;
    Int_t edge[kNClass]; // lower totnMIP edge of each class, central -> peripheral; last is 0
  };

  // Best fits (find_best_nbd_OO200.py) and edges (CentCuts_OO200.C), Oct 23 2026
  const Int_t kNVar = 6;
  const Variation kVar[kNVar] = {
    { "nominal",   "ncoll_npart.root",           24.0, 6.0, 0.200, 0.00, { 364, 310, 234, 130, 36, 0 } },
    { "small",     "ncoll_npart_small.root",     26.0, 4.0, 0.175, 0.04, { 365, 311, 235, 133, 37, 0 } },
    { "large",     "ncoll_npart_large.root",     24.0, 4.0, 0.175, 0.00, { 365, 312, 236, 133, 37, 0 } },
    { "smallXsec", "ncoll_npart_smallXsec.root", 28.0, 4.0, 0.075, 0.00, { 368, 315, 240, 139, 41, 0 } },
    { "largeXsec", "ncoll_npart_largeXsec.root", 24.0, 6.0, 0.225, 0.04, { 366, 312, 236, 131, 36, 0 } },
    { "gauss",     "ncoll_npart_gauss.root",     26.0, 4.0, 0.125, 0.00, { 366, 312, 236, 134, 39, 0 } }
  };

  // NBD(npp, k; n) as StNegativeBinomial::GetNegativeBinomial
  Double_t Nbd(Int_t n, Double_t npp, Double_t k)
  {
    const Double_t c = TMath::Exp(TMath::LnGamma(n + k) - TMath::LnGamma(n + 1) - TMath::LnGamma(k));
    const Double_t t = n * TMath::Log(npp / k) - (n + k) * TMath::Log(npp / k + 1.0);
    return c * TMath::Exp(t);
  }

  Int_t ClassOf(Int_t mult, const Int_t* edge)
  {
    for (Int_t i = 0; i < kNClass; i++) if (mult >= edge[i]) return i;
    return kNClass - 1;
  }
}

void NpartPerClass_OO200(const Int_t nevents = 2000000, const UInt_t seed = 20261023)
{
  gRandom = new TRandom3(seed); // used by TH2::GetRandom2 and TH1::GetRandom
  TRandom3 rnd(seed + 1);

  Double_t meanNpart[kNVar][kNClass], rmsNpart[kNVar][kNClass];
  Double_t meanNcoll[kNVar][kNClass], rmsNcoll[kNVar][kNClass];
  Double_t frac[kNVar][kNClass];

  std::ofstream out("NpartPerClass_OO200.txt");

  for (Int_t iv = 0; iv < kNVar; iv++) {
    const Variation& v = kVar[iv];
    TFile* f = TFile::Open(v.glauberFile);
    if (!f || f->IsZombie()) { std::cout << "Cannot open " << v.glauberFile << std::endl; return; }
    TH2D* h = (TH2D*) f->Get("hNcoll_Npart");
    if (!h) { std::cout << "No hNcoll_Npart in " << v.glauberFile << std::endl; return; }

    TH1D hNbd(Form("hNbd_%d", iv), "", 100, 0, 100);
    for (Int_t i = 0; i < 100; i++) hNbd.SetBinContent(i + 1, Nbd(i, v.npp, v.k));

    Double_t n[kNClass], s1[kNClass], s2[kNClass], c1[kNClass], c2[kNClass];
    for (Int_t ic = 0; ic < kNClass; ic++) { n[ic] = s1[ic] = s2[ic] = c1[ic] = c2[ic] = 0; }

    Int_t ievent = 0;
    while (ievent < nevents) {
      Double_t npart, ncoll;
      h->GetRandom2(npart, ncoll);
      npart = TMath::Floor(npart);
      ncoll = TMath::Floor(ncoll);
      if (npart < 2 || ncoll < 1) continue;

      const Double_t nsrc = (v.x == 0.0) ? npart : 0.5 * (1.0 - v.x) * npart + v.x * ncoll;
      const Int_t nchInt = TMath::Nint(nsrc);
      Int_t multIdeal = 0;
      for (Int_t ich = 0; ich < nchInt; ich++) multIdeal += (Int_t) hNbd.GetRandom();
      const Double_t eff = 0.98 * (1.0 - multIdeal * v.d / 540.0);
      Int_t mult = 0;
      // StGlauberUtilities::GetUniform2() is uniform on (0, 2): 2*multIdeal trials each kept with
      // probability eff/2, i.e. on average eff * multIdeal
      for (Int_t im = 0; im < 2 * multIdeal; im++) if (rnd.Uniform(0, 2) <= eff) mult++;

      const Int_t ic = ClassOf(mult, v.edge);
      n[ic] += 1; s1[ic] += npart; s2[ic] += npart * npart; c1[ic] += ncoll; c2[ic] += ncoll * ncoll;
      ievent++;
    }

    for (Int_t ic = 0; ic < kNClass; ic++) {
      frac[iv][ic] = n[ic] / nevents;
      meanNpart[iv][ic] = (n[ic] > 0) ? s1[ic] / n[ic] : 0;
      rmsNpart[iv][ic]  = (n[ic] > 0) ? TMath::Sqrt(TMath::Max(0.0, s2[ic] / n[ic] - meanNpart[iv][ic] * meanNpart[iv][ic])) : 0;
      meanNcoll[iv][ic] = (n[ic] > 0) ? c1[ic] / n[ic] : 0;
      rmsNcoll[iv][ic]  = (n[ic] > 0) ? TMath::Sqrt(TMath::Max(0.0, c2[ic] / n[ic] - meanNcoll[iv][ic] * meanNcoll[iv][ic])) : 0;
    }

    TString head = Form("== %s (%s): npp = %.3f, k = %.3f, x = %.3f, d = %.3f, %d events\n",
                        v.name, v.glauberFile, v.npp, v.k, v.x, v.d, nevents);
    TString cols = "  class     edge   model %     <Npart>  rms    <Ncoll>  rms    <Nspec>\n";
    std::cout << head.Data() << cols.Data();
    out << head.Data() << cols.Data();
    for (Int_t ic = 0; ic < kNClass; ic++) {
      TString line = Form("  %-8s %5d   %6.2f   %8.2f %5.2f  %8.2f %5.2f  %8.2f\n",
                          kClassName[ic], v.edge[ic], 100.0 * frac[iv][ic],
                          meanNpart[iv][ic], rmsNpart[iv][ic], meanNcoll[iv][ic], rmsNcoll[iv][ic],
                          kA2 - meanNpart[iv][ic]);
      std::cout << line.Data();
      out << line.Data();
    }
    f->Close();
  }

  // Systematic: largest deviation of any variation from the nominal, per class
  TString head = "\n== Glauber systematic: nominal value and largest deviation of any variation (relative)\n"
                 "  class      <Npart>   max dev        <Ncoll>   max dev        <Nspec>   max dev\n";
  std::cout << head.Data();
  out << head.Data();
  for (Int_t ic = 0; ic < kNClass; ic++) {
    Double_t dNp = 0, dNc = 0;
    for (Int_t iv = 1; iv < kNVar; iv++) {
      if (TMath::Abs(meanNpart[iv][ic] - meanNpart[0][ic]) > TMath::Abs(dNp)) dNp = meanNpart[iv][ic] - meanNpart[0][ic];
      if (TMath::Abs(meanNcoll[iv][ic] - meanNcoll[0][ic]) > TMath::Abs(dNc)) dNc = meanNcoll[iv][ic] - meanNcoll[0][ic];
    }
    const Double_t ns = kA2 - meanNpart[0][ic];
    TString line = Form("  %-8s %8.2f  %+6.2f (%4.1f%%)  %8.2f  %+6.2f (%4.1f%%)  %8.2f  %+6.2f (%4.1f%%)\n",
                        kClassName[ic],
                        meanNpart[0][ic], dNp, 100.0 * TMath::Abs(dNp) / meanNpart[0][ic],
                        meanNcoll[0][ic], dNc, 100.0 * TMath::Abs(dNc) / meanNcoll[0][ic],
                        ns, -dNp, 100.0 * TMath::Abs(dNp) / ns);
    std::cout << line.Data();
    out << line.Data();
  }
  out.close();
  std::cout << "\nWritten NpartPerClass_OO200.txt" << std::endl;
}
#else
void NpartPerClass_OO200(const Int_t nevents = 2000000, const UInt_t seed = 20261023);
#endif
