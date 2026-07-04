// DiffMTM0ZTPCDumps.C
//
// Library-free (no Histo2D/ZFitter dependency): reads the two CSVs written by
// DumpMTM0ZTPCProj2D.C (one from lightflavorspectra_etof, one from
// lightflavorspectra_OO200) and compares the RAW (pre-fit) per-(rapidity bin, mT-m0
// bin) dE/dx-amplitude counts bin-by-bin -- i.e. the actual input data each repo's TPC
// fits are drawn from, one level upstream of the dEdxSpectra_* fit-yield comparison
// already done by CompareSpectra_OO200_vs_etof.C.
//
// Usage:
//   root -l -q 'DiffMTM0ZTPCDumps.C("/Users/aliggett/data/OO/etof_compare/dump_etof_ProtonPlus_TPC.csv","/Users/aliggett/data/OO/etof_compare/dump_oo200_ProtonPlus_TPC.csv")'

#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <cmath>
using namespace std;

struct BinKey{
  int cent, rapIndex, mtm0Bin;
  bool operator<(const BinKey& o) const {
    if(cent != o.cent) return cent < o.cent;
    if(rapIndex != o.rapIndex) return rapIndex < o.rapIndex;
    return mtm0Bin < o.mtm0Bin;
  }
};

map<BinKey,double> readDump(string a_csvPath){
  map<BinKey,double> result;
  ifstream in(a_csvPath.c_str());
  if(!in.is_open()){ cout << "Could not open " << a_csvPath << endl; return result; }

  string line;
  getline(in, line); // header
  while(getline(in, line)){
    if(line.empty()) continue;
    stringstream ss(line);
    string field;
    vector<string> fields;
    while(getline(ss, field, ',')) fields.push_back(field);
    if(fields.size() < 6) continue;

    BinKey key;
    key.cent     = stoi(fields[0]);
    key.rapIndex = stoi(fields[1]);
    key.mtm0Bin  = stoi(fields[2]);
    double integral = stod(fields[5]);
    result[key] = integral;
  }
  return result;
}

void DiffMTM0ZTPCDumps(string a_etofCsv, string a_oo200Csv, string a_outDir = "./spectraComparison_OO200_vs_etof/", double a_relTolerance = 0.05){
  gStyle->SetOptStat(0);

  map<BinKey,double> etofData  = readDump(a_etofCsv);
  map<BinKey,double> oo200Data = readDump(a_oo200Csv);

  if(etofData.empty() || oo200Data.empty()){
    cout << "One or both dumps were empty -- nothing to compare." << endl;
    return;
  }

  // union of keys seen on either side
  map<BinKey,bool> allKeys;
  for(auto& kv : etofData)  allKeys[kv.first] = true;
  for(auto& kv : oo200Data) allKeys[kv.first] = true;

  int nCompared = 0, nBad = 0, nOneSided = 0;
  double worstDev = -1;
  BinKey worstKey{-1,-1,-1};

  // per (cent, rapIndex): accumulate ratios across mT-m0 bins
  map<pair<int,int>, vector<double>> ratiosByCentRap;

  for(auto& kv : allKeys){
    const BinKey& key = kv.first;
    double vEtof  = etofData.count(key)  ? etofData[key]  : 0;
    double vOO200 = oo200Data.count(key) ? oo200Data[key] : 0;
    if(vEtof == 0 && vOO200 == 0) continue;

    if(vEtof == 0 || vOO200 == 0){
      nOneSided++;
      nCompared++;
      continue;
    }

    double ratio = vEtof / vOO200;
    nCompared++;
    double dev = fabs(ratio - 1.0);
    if(dev > a_relTolerance) nBad++;
    if(dev > worstDev){ worstDev = dev; worstKey = key; }

    ratiosByCentRap[{key.cent, key.rapIndex}].push_back(ratio);
  }

  nBad += nOneSided;

  cout << "=== Raw mT-m0-vs-dE/dx bin comparison (pre-fit input data) ===" << endl;
  cout << nCompared << " populated bins compared, " << nOneSided << " one-sided-only, "
       << nBad << " outside " << (a_relTolerance*100) << "% tol, worst dev="
       << (worstDev < 0 ? 0 : worstDev) << " at Cent" << worstKey.cent
       << " rapIndex" << worstKey.rapIndex << " mtm0Bin" << worstKey.mtm0Bin
       << "  -> " << (nBad == 0 ? "PASS" : "FAIL") << endl;

  gSystem->mkdir(a_outDir.c_str(), true);

  // group by rapIndex -> one combined plot per rapidity bin (x=centrality, y=mean ratio)
  map<int, vector<int>> centsByRap;
  for(auto& kv : ratiosByCentRap) centsByRap[kv.first.second].push_back(kv.first.first);

  for(auto& rapEntry : centsByRap){
    int rapIndex = rapEntry.first;
    vector<int>& cents = rapEntry.second;

    vector<double> xVals, yVals, yErrs;
    for(int cent : cents){
      vector<double>& ratios = ratiosByCentRap[{cent, rapIndex}];
      double sum = 0, sum2 = 0;
      for(double r : ratios){ sum += r; sum2 += r*r; }
      double mean = sum / ratios.size();
      double rms  = sqrt(max(0.0, sum2/ratios.size() - mean*mean));
      double sem  = rms / sqrt((double) ratios.size());
      xVals.push_back(cent);
      yVals.push_back(mean);
      yErrs.push_back(sem);
    }

    TGraphErrors* g = new TGraphErrors((int) xVals.size(), xVals.data(), yVals.data(), nullptr, yErrs.data());
    g->SetName(Form("g_rawTPC_RapBin%02d", rapIndex));
    g->SetTitle(Form("Raw TPC input (pre-fit): etof / OO200, ProtonPlus, RapBin %02d;Centrality bin;etof / OO200", rapIndex));
    g->SetLineColor(kRed);
    g->SetMarkerColor(kRed);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.8);

    TCanvas* c = new TCanvas(Form("c_rawTPC_RapBin%02d", rapIndex), Form("Raw TPC RapBin %02d", rapIndex), 900, 700);
    c->SetGridy();

    g->Draw("AP");
    g->GetYaxis()->SetRangeUser(0.5, 1.5);
    if(xVals.size() > 0) g->GetXaxis()->SetLimits(xVals.front() - 0.5, xVals.back() + 0.5);
    gPad->Modified();
    gPad->Update();

    TLine* line = new TLine(g->GetXaxis()->GetXmin(), 1, g->GetXaxis()->GetXmax(), 1);
    line->SetLineStyle(2);
    line->Draw();

    string outName = Form("%s/rawTPC_ProtonPlus_RapBin%02d_ratio_summary.png", a_outDir.c_str(), rapIndex);
    c->SaveAs(outName.c_str());
  }

  cout << "Plots written to: " << a_outDir << " (rawTPC_ProtonPlus_RapBin##_ratio_summary.png)" << endl;
}
