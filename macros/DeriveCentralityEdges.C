// DeriveCentralityEdges.C -- lightflavorspectra_OO200
//
// WHY THIS EXISTS: SetCutClass.C's setCentralities(5, centCutsArray, percents) call
// uses refMult thresholds {44,37,28,17,5} for the 5/10/20/40/80% bins -- these are the
// OFFICIAL cuts from a separate, already-validated analysis, not something to re-derive
// or second-guess. (A 2026-07-07 investigation into a reported population skew traced
// that symptom to a stale compiled bin/SetCutClass_C.so on SDCC still running an older
// cut scheme, not to these values being wrong -- see README's "Centrality binning"
// section, "UPDATE 2026-07-07".)
//
// This macro is an independent sanity-check tool: it derives the refMult value at each
// percentile boundary directly from a yield file's own "refMult" histogram (filled by
// PicoBinner, post event/track/vertex cuts, pre-centrality-cut -- see PicoBinner.cxx's
// a_makeBasicHistos block), by integrating the distribution from the HIGHEST
// multiplicity bin downward (centrality percentiles are defined most-central-first:
// 0-5% = highest-multiplicity 5% of events). Use it to cross-check that this analysis's
// own event sample is consistent with the official cuts -- not to replace them.
//
// Usage:
//   root -l -b -q 'DeriveCentralityEdges.C("/path/to/yieldHistos_OOGeV_proton_....root")'
//   root -l -b -q 'DeriveCentralityEdges.C("file.root","refMult",{5,10,20,40,80})'  // explicit percents
//
// Output: prints a ready-to-paste centCutsArray[]/percents[] pair for SetCutClass.C's
// setCentralities() call. Descending edges, ascending percents, matching the ordering
// CutClass::centralityIndex() requires (see that function's / README's "Centrality
// binning" section for why the order matters).
//
// NOTE: run this once per species yield file if their event samples differ (different
// triggers/pileup rejection could in principle shift the distribution) -- if pion/kaon/
// proton all share the exact same event-level cuts and trigger, one file's refMult
// histogram is representative of all three and only needs to be run once.

#include "TFile.h"
#include "TH1.h"
#include <iostream>
#include <vector>
using namespace std;

void DeriveCentralityEdges(string a_yieldFile,
                            string a_histName = "refMult",
                            vector<double> a_percents = {5,10,20,40,80}){

  TFile* inFile = new TFile(a_yieldFile.c_str(),"READ");
  if(!inFile || inFile->IsZombie()){
    cout << "!!! Could not open " << a_yieldFile << endl;
    return;
  }

  TH1* refMult = (TH1*) inFile->Get(a_histName.c_str());
  if(!refMult){
    cout << "!!! Histogram \"" << a_histName << "\" not found in " << a_yieldFile
         << " -- run FindHistoInFile.C to see what's actually in this file." << endl;
    return;
  }

  int nBins = refMult->GetNbinsX();
  double total = refMult->Integral(1, nBins);
  if(total <= 0){
    cout << "!!! " << a_histName << " has zero (or negative) integral -- nothing to derive from." << endl;
    return;
  }

  cout << "Total events in " << a_histName << ": " << total << endl;

  vector<double> edges(a_percents.size(), -1.0);
  double cumulative = 0.0;
  size_t percentIdx = 0;

  // integrate from the highest-multiplicity bin down to the lowest -- most central first
  for(int bin = nBins; bin >= 1; bin--){
    cumulative += refMult->GetBinContent(bin);
    double cumFracPercent = 100.0 * cumulative / total;
    while(percentIdx < a_percents.size() && cumFracPercent >= a_percents[percentIdx]){
      // this bin's LOW edge is the refMult value at/above which you get the most-central
      // cumFracPercent% of events -- i.e. the cut value for this percentile boundary
      edges[percentIdx] = refMult->GetBinLowEdge(bin);
      percentIdx++;
    }
    if(percentIdx >= a_percents.size()) break;
  }

  if(percentIdx < a_percents.size()){
    cout << "!!! WARNING: only resolved " << percentIdx << "/" << a_percents.size()
         << " percentile edges -- histogram may not extend low enough in multiplicity "
         << "to cover the largest requested percentile (" << a_percents.back() << "%)." << endl;
  }

  cout << "\nDerived from " << a_yieldFile << " (\"" << a_histName << "\"):\n" << endl;

  cout << "  Percentile boundary : refMult threshold" << endl;
  for(size_t iii = 0; iii < a_percents.size(); iii++){
    cout << "  " << a_percents[iii] << "% : " << edges[iii] << endl;
  }

  cout << "\n// paste into SetCutClass.C in place of the current centCutsArray/percents:" << endl;
  cout << "double centCutsArray[" << edges.size() << "] = {";
  for(size_t iii = 0; iii < edges.size(); iii++){
    cout << edges[iii];
    if(iii + 1 < edges.size()) cout << ",";
  }
  cout << "};" << endl;

  cout << "int percents[" << a_percents.size() << "] = {";
  for(size_t iii = 0; iii < a_percents.size(); iii++){
    cout << a_percents[iii];
    if(iii + 1 < a_percents.size()) cout << ",";
  }
  cout << "};" << endl;
  cout << "cuts->setCentralities(" << edges.size() << ", &centCutsArray[0], &percents[0]);" << endl;

  inFile->Close();
}
