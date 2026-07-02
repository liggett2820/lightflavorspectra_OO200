// DrawSpectraAndDnDy.C -- lightflavorspectra_OO200
//
// Reads SpectraFitter's output file (RunSpectraFitter.C's a_outputFile argument, e.g.
// "fit_output.root") and saves two families of diagnostic plots as .png files:
//
//   1. Invariant spectra: TotalSpectra/<ParticleCharge>/totalSpectra_<ParticleCharge>_
//      Cent_<CC>_yIndex<YY> (TGraphErrors, invariant yield vs m_T-m_0). One canvas per
//      (particle, charge, centrality bin), overlaying every rapidity bin found for that
//      centrality, offset by decades for readability (standard convention for this kind
//      of plot -- see a_offsetByDecade below).
//
//   2. dN/dy vs rapidity: ParticleYields/<ParticleCharge>/dNdy_<ParticleCharge>_Cent<CC>_
//      <Variant> (TH1D). One canvas per (particle, charge, fit variant), overlaying
//      every centrality bin found for that variant. Restricted to the "Nominal" and
//      "Systematic" variants -- SpectraFitter also writes several intermediate
//      extrapolation-piece histograms (LowMtM0Extrap/MidMtM0Extrap/HighMtM0Extrap/
//      Fiducial/Fiducial_Measured/Fiducial_And_Extrap, per source/SpectraFitter.cxx's
//      writeOutputs()) that are components of the dN/dy determination itself rather
//      than "the" dN/dy distribution, so those are deliberately not drawn here.
//
// Nothing here assumes a specific number of centrality/rapidity bins actually got
// filled -- RunSpectraFitter.C's a_numCent and CutClass's rapidity binning both affect
// that, and this macro doesn't load either. Instead it probes every (centrality,
// rapidity) slot up to the fixed array bounds declared in headers/SpectraFitter.h
// (9 particles / 16 centralities / g_SpecFit_MaxNumRapBins=35 rapidity bins -- only
// partIndex 0-2 are ever fit by RunSpectraFitter.C, so this only loops 0-2) via
// TDirectory::Get(), which returns NULL for anything that was never written
// (HistogramUtilities::ConditionalWrite skips NULL/empty objects at write time), and
// simply skips whatever isn't there.
//
// Object names are reconstructed with ParticleInfo::GetParticleName(), instantiated
// the same way SpectraFitter itself does (`new ParticleInfo()`, no args -- see
// source/SpectraFitter.cxx's constructor) so the naming matches exactly.
//
// PREREQUISITE: a_fitOutputFile is RunSpectraFitter.C's output (see RunSpectraFitter.C).
//
// USAGE (run from within macros/, loading LoadDrawingLibs.C first):
//   root -l -q -b LoadDrawingLibs.C 'DrawSpectraAndDnDy.C("fit_output.root")'
//   root -l -q -b LoadDrawingLibs.C 'DrawSpectraAndDnDy.C("fit_output.root","./myPlots",true,true,true)'

#include "../makefile_toggles.h"
// Needed so cling knows this type when the macro is interpreted -- gSystem->Load() in
// LoadDrawingLibs.C loads the compiled implementation, but doesn't by itself guarantee
// the class declaration is known to the interpreter (same reasoning as RunPicoBinner.C).
#include "../submodule/ParticleInfo/ParticleInfo/ParticleInfo.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TMath.h"
#include "TSystem.h"
#include <vector>
#include <string>

// Fixed array bounds copied from headers/SpectraFitter.h -- see header comment above.
static const int g_MaxCent    = 16;
static const int g_MaxRapBins = 35;

void DrawSpectraAndDnDy(string a_fitOutputFile,
                         string a_outputDir      = "./spectra_plots", // created if missing
                         bool   a_drawSpectra     = true,
                         bool   a_drawDnDy        = true,
                         bool   a_offsetByDecade  = true){  // multiply successive rapidity-bin spectra by 10x for readability (standard for this plot type; only affects the spectra plots, not dN/dy)

  gSystem->mkdir(a_outputDir.c_str(), true);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile* inFile = TFile::Open(a_fitOutputFile.c_str(), "READ");
  if(!inFile || inFile->IsZombie()){
    cout << "ERROR: could not open " << a_fitOutputFile << endl;
    return;
  }

  ParticleInfo* particleInfo = new ParticleInfo(); // default ctor -- matches SpectraFitter's own m_partInfo exactly

  const int nColors = 12;
  int colors[nColors] = {kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+1,
                          kCyan+2, kViolet+1, kSpring+4, kAzure+1, kPink+1, kGray+2};

  int nSpectraPlotsWritten = 0;
  int nDnDyPlotsWritten    = 0;

  //=========================================================================
  // 1. SPECTRA
  //=========================================================================
  if(a_drawSpectra && inFile->GetDirectory("TotalSpectra")){
    for(int partIndex = 0; partIndex < 3; partIndex++){ // pion/kaon/proton -- only species RunSpectraFitter.C ever fits
      for(int chargeSign = 0; chargeSign < 2; chargeSign++){
        int charge = (chargeSign == 0) ? 1 : -1;
        TString particleName = particleInfo->GetParticleName(partIndex, charge); // e.g. "PionPlus"
        TString particleSymbol = particleInfo->GetParticleSymbol(partIndex, charge);

        TDirectory* particleDir = inFile->GetDirectory(Form("TotalSpectra/%s", particleName.Data()));
        if(!particleDir) continue;

        for(int centIndex = 0; centIndex < g_MaxCent; centIndex++){
          std::vector<TGraphErrors*> graphs;
          std::vector<int> yIndices;

          for(int yIndex = 0; yIndex < g_MaxRapBins; yIndex++){
            TString objName = Form("totalSpectra_%s_Cent_%02d_yIndex%02d", particleName.Data(), centIndex, yIndex);
            TGraphErrors* g = (TGraphErrors*) particleDir->Get(objName.Data());
            if(g && g->GetN() > 0){
              graphs.push_back(g);
              yIndices.push_back(yIndex);
            }
          }
          if(graphs.empty()) continue;

          TCanvas* c = new TCanvas("c_spectra", "c_spectra", 50, 50, 800, 700);
          c->SetLogy();
          c->SetTicks();
          c->SetLeftMargin(0.13);
          c->SetRightMargin(0.04);
          c->SetTopMargin(0.06);
          c->SetBottomMargin(0.11);

          // Clone + optionally decade-offset so overlapping rapidity bins stay legible.
          std::vector<TGraphErrors*> drawGraphs;
          double maxY = -1e300, minY = 1e300, maxX = -1e300;
          for(unsigned int i = 0; i < graphs.size(); i++){
            TGraphErrors* g = (TGraphErrors*) graphs[i]->Clone(Form("%s_draw", graphs[i]->GetName()));
            double scale = a_offsetByDecade ? TMath::Power(10.0, (double) i) : 1.0;
            if(scale != 1.0){
              for(int p = 0; p < g->GetN(); p++){
                g->SetPoint(p, g->GetX()[p], g->GetY()[p]*scale);
                g->SetPointError(p, g->GetErrorX(p), g->GetErrorY(p)*scale);
              }
            }
            drawGraphs.push_back(g);
            if(g->GetN() > 0){
              double locMaxY = TMath::MaxElement(g->GetN(), g->GetY());
              double locMinY = TMath::MinElement(g->GetN(), g->GetY());
              double locMaxX = TMath::MaxElement(g->GetN(), g->GetX());
              if(locMaxY > maxY) maxY = locMaxY;
              if(locMinY < minY && locMinY > 0) minY = locMinY;
              if(locMaxX > maxX) maxX = locMaxX;
            }
          }
          if(minY <= 0 || minY >= 1e299) minY = 1e-8;
          if(maxX <= 0) maxX = 3.0;

          TH1F* frame = c->DrawFrame(0.0, minY*0.2, maxX*1.15, maxY*5.0);
          frame->GetXaxis()->SetTitle("m_{T} - m_{0} (GeV/c^{2})");
          frame->GetYaxis()->SetTitle("(1/N_{evt})(1/2#pi m_{T}) d^{2}N/dm_{T}dy (GeV/c^{2})^{-2}");
          frame->GetXaxis()->SetTitleSize(0.035);
          frame->GetYaxis()->SetTitleSize(0.035);
          frame->GetYaxis()->SetTitleOffset(1.5);

          TLegend* leg = new TLegend(0.62, 0.55, 0.94, 0.94);
          leg->SetBorderSize(0);
          leg->SetFillStyle(0);
          leg->SetTextSize(0.025);
          leg->SetHeader(a_offsetByDecade ? "yIndex (#times10^{n})" : "yIndex");

          for(unsigned int i = 0; i < drawGraphs.size(); i++){
            int color = colors[i % nColors];
            drawGraphs[i]->SetMarkerColor(color);
            drawGraphs[i]->SetLineColor(color);
            drawGraphs[i]->SetMarkerStyle(20);
            drawGraphs[i]->SetMarkerSize(0.7);
            drawGraphs[i]->Draw("PZ SAME");
            if(a_offsetByDecade) leg->AddEntry(drawGraphs[i], Form("%02d (#times10^{%d})", yIndices[i], (int) i), "P");
            else                 leg->AddEntry(drawGraphs[i], Form("%02d", yIndices[i]), "P");
          }
          leg->Draw("SAME");

          TPaveText* title = new TPaveText(0.15, 0.87, 0.62, 0.94, "brNDC");
          title->SetBorderSize(0);
          title->SetFillStyle(0);
          title->SetTextAlign(12);
          title->AddText(Form("%s spectra -- Cent bin %02d", particleSymbol.Data(), centIndex));
          title->Draw("SAME");

          TString outName = Form("%s/spectra_%s_Cent%02d.png", a_outputDir.c_str(), particleName.Data(), centIndex);
          c->SaveAs(outName.Data());
          nSpectraPlotsWritten++;
          delete c;
        } // centIndex
      } // chargeSign
    } // partIndex
  }

  //=========================================================================
  // 2. dN/dy vs RAPIDITY
  //=========================================================================
  if(a_drawDnDy && inFile->GetDirectory("ParticleYields")){
    const char* variants[2] = {"Nominal", "Systematic"};

    for(int partIndex = 0; partIndex < 3; partIndex++){
      for(int chargeSign = 0; chargeSign < 2; chargeSign++){
        int charge = (chargeSign == 0) ? 1 : -1;
        TString particleName   = particleInfo->GetParticleName(partIndex, charge);
        TString particleSymbol = particleInfo->GetParticleSymbol(partIndex, charge);

        TDirectory* particleDir = inFile->GetDirectory(Form("ParticleYields/%s", particleName.Data()));
        if(!particleDir) continue;

        for(int v = 0; v < 2; v++){
          std::vector<TH1D*> hists;
          std::vector<int> centIndices;

          for(int centIndex = 0; centIndex < g_MaxCent; centIndex++){
            TString objName = Form("dNdy_%s_Cent%02d_%s", particleName.Data(), centIndex, variants[v]);
            TH1D* h = (TH1D*) particleDir->Get(objName.Data());
            if(h){
              hists.push_back(h);
              centIndices.push_back(centIndex);
            }
          }
          if(hists.empty()) continue;

          TCanvas* c = new TCanvas("c_dndy", "c_dndy", 50, 50, 800, 700);
          c->SetTicks();
          c->SetLeftMargin(0.13);
          c->SetRightMargin(0.04);
          c->SetTopMargin(0.06);
          c->SetBottomMargin(0.11);

          double maxY = -1e300, minX = 1e300, maxX = -1e300;
          for(unsigned int i = 0; i < hists.size(); i++){
            if(hists[i]->GetMaximum() > maxY) maxY = hists[i]->GetMaximum();
            double lowEdge  = hists[i]->GetXaxis()->GetXmin();
            double highEdge = hists[i]->GetXaxis()->GetXmax();
            if(lowEdge  < minX) minX = lowEdge;
            if(highEdge > maxX) maxX = highEdge;
          }
          if(maxY <= 0) maxY = 1.0;

          TH1F* frame = c->DrawFrame(minX, 0.0, maxX, maxY*1.3);
          frame->GetXaxis()->SetTitle("rapidity y");
          frame->GetYaxis()->SetTitle("dN/dy");
          frame->GetXaxis()->SetTitleSize(0.035);
          frame->GetYaxis()->SetTitleSize(0.035);
          frame->GetYaxis()->SetTitleOffset(1.5);

          TLegend* leg = new TLegend(0.7, 0.6, 0.94, 0.94);
          leg->SetBorderSize(0);
          leg->SetFillStyle(0);
          leg->SetTextSize(0.028);
          leg->SetHeader("Centrality bin");

          for(unsigned int i = 0; i < hists.size(); i++){
            int color = colors[i % nColors];
            hists[i]->SetLineColor(color);
            hists[i]->SetMarkerColor(color);
            hists[i]->SetMarkerStyle(20);
            hists[i]->SetMarkerSize(0.8);
            hists[i]->SetLineWidth(2);
            hists[i]->Draw("PE SAME");
            leg->AddEntry(hists[i], Form("%02d", centIndices[i]), "PL");
          }
          leg->Draw("SAME");

          TPaveText* title = new TPaveText(0.15, 0.87, 0.62, 0.94, "brNDC");
          title->SetBorderSize(0);
          title->SetFillStyle(0);
          title->SetTextAlign(12);
          title->AddText(Form("%s dN/dy -- %s fit", particleSymbol.Data(), variants[v]));
          title->Draw("SAME");

          TString outName = Form("%s/dNdy_%s_%s.png", a_outputDir.c_str(), particleName.Data(), variants[v]);
          c->SaveAs(outName.Data());
          nDnDyPlotsWritten++;
          delete c;
        } // variant
      } // chargeSign
    } // partIndex
  }

  cout << "-------------  End of DrawSpectraAndDnDy Macro  -------------" << endl;
  cout << "Wrote " << nSpectraPlotsWritten << " spectra plot(s) and " << nDnDyPlotsWritten << " dN/dy plot(s) to " << a_outputDir << endl;
  if(nSpectraPlotsWritten == 0 && a_drawSpectra) cout << "WARNING: no spectra objects found under TotalSpectra/ -- check a_fitOutputFile is SpectraFitter's writeOutputs() output." << endl;
  if(nDnDyPlotsWritten == 0 && a_drawDnDy)       cout << "WARNING: no dN/dy objects found under ParticleYields/ -- check a_fitOutputFile is SpectraFitter's writeOutputs() output." << endl;

  inFile->Close();
}
