// CheckBTOFCorrection.C
//
// Sanity check for spectra_OO200_Proton_Corrected.root: for every
// raw_InvBetaBTOFSpectra_<Species>_<Cent>_yIndex_NN graph found under
// Spectra_ProtonPlus / Spectra_ProtonMinus, finds its corrected_ counterpart and
// compares them by matching points on X (mT-m0) value -- NOT by array index, since the
// BTOF correction can (and does) drop points entirely (e.g. via a quality/error cut),
// so raw and corrected commonly have different lengths. Reports, per graph: how many
// raw points got dropped (no X match in corrected), and for the points that DID
// survive, the min/mean/max corrected/raw Y ratio, so you can tell whether the
// correction is scaling the surviving points up as expected (ratio == 1.0 everywhere
// would mean it didn't apply) and how much data is being lost to drops.
//
// Usage (no custom-class dictionaries are strictly required just to read TGraphErrors,
// but loading them first avoids "no dictionary" noise when ROOT scans the file's
// StreamerInfo -- run from macros/, same convention as everything else in this repo):
//   root -l -q LoadRawSpectraModifierLibs.C 'CheckBTOFCorrection.C("/Users/aliggett/data/OO/zfitter_output/spectra_OO200_Proton_Corrected.root")'
//
// Writes a plain-text summary to btof_correction_check.txt in the current directory
// (in addition to printing to stdout), so the results are easy to hand off/attach
// elsewhere.

#include "TFile.h"
#include "TDirectory.h"
#include "TGraphErrors.h"
#include "TKey.h"
#include "TList.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

// Compares one raw/corrected TGraphErrors pair by matching points on X (mT-m0) value
// rather than assuming identical indices/lengths -- BTOF correction can (and does) drop
// points entirely (e.g. via the maxFracError cut), so index-based comparison bails out
// the moment lengths differ and tells you nothing about the points that DID survive.
// Reports: how many raw points had no X match in corrected (dropped), how many raw
// points had no counterpart to compare against, and for the matched points, the
// corrected/raw Y ratio stats.
void compareOnePair(const string& a_label, TGraphErrors* a_raw, TGraphErrors* a_corrected, ofstream& a_out){
  int nRaw = a_raw->GetN();
  int nCorr = a_corrected->GetN();

  if(nRaw == 0 && nCorr == 0){
    string msg = "  " + a_label + ": both empty (0 points) -- nothing to compare";
    cout << msg << endl;
    a_out << msg << endl;
    return;
  }

  // xTol: tolerance for considering two X values "the same" point. mT-m0 bins in this
  // analysis are spaced far more than 1e-6 apart, so a small absolute tolerance is safe
  // and avoids needing to know the exact bin width here.
  const double xTol = 1e-6;

  double minRatio = 1e300, maxRatio = -1e300, sumRatio = 0.0;
  int nUnchanged = 0;      // ratio within 0.1% of 1.0
  int nComparable = 0;     // matched AND raw Y > 0
  int nRawZeroOrNeg = 0;   // matched but raw <= 0 -- ratio undefined, skipped from stats
  int nDropped = 0;        // raw point with no X match in corrected at all

  for(int i = 0; i < nRaw; i++){
    double xRaw, yRaw;
    a_raw->GetPoint(i, xRaw, yRaw);

    // Linear scan for the corrected point with matching X. nRaw/nCorr are small
    // (tens of points per graph), so O(n^2) here is negligible.
    int jMatch = -1;
    for(int j = 0; j < nCorr; j++){
      double xCorr, yCorr;
      a_corrected->GetPoint(j, xCorr, yCorr);
      if(fabs(xCorr - xRaw) < xTol){
        jMatch = j;
        break;
      }
    }

    if(jMatch < 0){
      nDropped++;
      continue;
    }

    double xCorr, yCorr;
    a_corrected->GetPoint(jMatch, xCorr, yCorr);

    if(yRaw <= 0){
      nRawZeroOrNeg++;
      continue;
    }
    double ratio = yCorr / yRaw;
    if(ratio < minRatio) minRatio = ratio;
    if(ratio > maxRatio) maxRatio = ratio;
    sumRatio += ratio;
    nComparable++;
    if(fabs(ratio - 1.0) < 0.001) nUnchanged++;
  }

  // Corrected points with no raw match at all would mean the correction invented new
  // X values, which shouldn't happen but is worth flagging if it does.
  int nExtraInCorrected = 0;
  for(int j = 0; j < nCorr; j++){
    double xCorr, yCorr;
    a_corrected->GetPoint(j, xCorr, yCorr);
    bool found = false;
    for(int i = 0; i < nRaw; i++){
      double xRaw, yRaw;
      a_raw->GetPoint(i, xRaw, yRaw);
      if(fabs(xRaw - xCorr) < xTol){ found = true; break; }
    }
    if(!found) nExtraInCorrected++;
  }

  string msg;
  if(nComparable == 0){
    msg = "  " + a_label + ": raw=" + to_string(nRaw) + " corrected=" + to_string(nCorr)
        + " -- dropped(no X match)=" + to_string(nDropped)
        + " | matched-but-raw<=0=" + to_string(nRawZeroOrNeg)
        + " | extra-in-corrected=" + to_string(nExtraInCorrected)
        + " -- no comparable (matched, raw>0) points";
  } else {
    double meanRatio = sumRatio / nComparable;
    msg = "  " + a_label + ": raw=" + to_string(nRaw) + " corrected=" + to_string(nCorr)
        + " -- " + to_string(nComparable) + " matched&comparable point(s) -- corrected/raw ratio: min=" + to_string(minRatio)
        + " mean=" + to_string(meanRatio) + " max=" + to_string(maxRatio)
        + " | unchanged(<0.1%% diff)=" + to_string(nUnchanged)
        + " | dropped(no X match)=" + to_string(nDropped)
        + " | matched-but-raw<=0=" + to_string(nRawZeroOrNeg)
        + " | extra-in-corrected=" + to_string(nExtraInCorrected);
  }
  cout << msg << endl;
  a_out << msg << endl;
}

void CheckBTOFCorrection(string a_correctedFile){
  ofstream out("btof_correction_check.txt");
  string header = "BTOF efficiency correction check -- file: " + a_correctedFile;
  cout << header << endl;
  out << header << endl;

  TFile* f = TFile::Open(a_correctedFile.c_str(), "READ");
  if(!f || f->IsZombie()){
    cout << "ERROR: could not open " << a_correctedFile << endl;
    out << "ERROR: could not open " << a_correctedFile << endl;
    return;
  }

  vector<string> speciesDirs = {"Spectra_ProtonPlus", "Spectra_ProtonMinus"};

  int totalPairsChecked = 0;
  int totalPairsMissingCorrected = 0;

  for(const string& dirName : speciesDirs){
    TDirectory* dir = (TDirectory*)f->Get(dirName.c_str());
    string dirMsg = "\n--- " + dirName + " ---";
    cout << dirMsg << endl;
    out << dirMsg << endl;
    if(!dir){
      string msg = "  (directory not found in file -- skipping)";
      cout << msg << endl;
      out << msg << endl;
      continue;
    }

    // Walk every key in this directory, find ones starting with "raw_", and look
    // for the matching "corrected_" key (same suffix after the prefix).
    TList* keys = dir->GetListOfKeys();
    if(!keys) continue;
    for(int i = 0; i < keys->GetEntries(); i++){
      TKey* key = (TKey*)keys->At(i);
      string keyName = key->GetName();
      string prefix = "raw_";
      if(keyName.rfind(prefix, 0) != 0) continue; // doesn't start with "raw_"

      string suffix = keyName.substr(prefix.size());
      string correctedName = "corrected_" + suffix;

      TGraphErrors* rawGraph = (TGraphErrors*)dir->Get(keyName.c_str());
      TGraphErrors* corrGraph = (TGraphErrors*)dir->Get(correctedName.c_str());

      if(!rawGraph){
        continue; // shouldn't happen since we got the name from the key list, but be safe
      }
      if(!corrGraph){
        string msg = "  " + suffix + ": no corrected_ counterpart found (looked for '" + correctedName + "')";
        cout << msg << endl;
        out << msg << endl;
        totalPairsMissingCorrected++;
        continue;
      }

      compareOnePair(suffix, rawGraph, corrGraph, out);
      totalPairsChecked++;
    }
  }

  string footer = "\nTOTAL: " + to_string(totalPairsChecked) + " raw/corrected pair(s) compared, "
                 + to_string(totalPairsMissingCorrected) + " raw graph(s) had no corrected_ counterpart.";
  cout << footer << endl;
  out << footer << endl;

  out.close();
  cout << "\nWrote btof_correction_check.txt in the current directory." << endl;

  f->Close();
}
