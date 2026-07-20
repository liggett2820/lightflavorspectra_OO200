// CheckChiSqrValues.C -- lightflavorspectra_OO200
//
// Diagnostic: for a given ZFitter TPC/BTOF output file pair, print real (yield content,
// chi^2/ndf) pairs for a sample of bins that have positive yield, plus summary
// statistics (min/max/mean/median) of chi^2/ndf across ALL such bins. Written after
// Andrew reported that PresentZFitterSpectra.C's new chi^2/ndf cut (a_maxChiSqrNdf, see
// that macro's header comment) suppressed EVERY rapidity bin at its initial guessed
// default of 3.0. This answers two questions directly from the real fit output instead
// of guessing again: (1) is suppressBadFitQualityBins() actually reading the right
// (rapBin, mTm0Bin) cell of the chi^2 histogram -- sanity-checked here by printing
// chi^2/ndf right next to the yield it corresponds to, for real bins -- and (2) what
// does this analysis's real chi^2/ndf distribution actually look like, so a_maxChiSqrNdf
// can be picked empirically instead of assuming the textbook ~1-3 "good fit" range
// applies here (it may well not, for these background-subtracted dE/dx/1-beta template
// fits).
//
// Object paths/names match PresentZFitterSpectra.C's own reads exactly -- see that
// macro's PresentZFitterSpectra() (yield histos) and the chi^2 histogram block right
// after it (including the ETOF "ZbTOF_ChiSqr_..." naming quirk, not reproduced here
// since ETOF isn't populated for this analysis yet).

#include <vector>
#include <algorithm>

// 2026-07-20: shared by both the full-distribution report and the high-mTm0-tail-only
// report below (see reportChiSqrStats) -- factored out so the two don't drift apart.
void reportChiSqrDistribution(string a_label, vector<double> a_values){
  if(a_values.empty()){
    cout << a_label << ": no bins in this subset -- nothing to summarize." << endl;
    return;
  }
  sort(a_values.begin(),a_values.end());
  double sum = 0;
  for(size_t i = 0; i < a_values.size(); i++) sum += a_values[i];
  double mean = sum / a_values.size();
  double median = a_values[a_values.size()/2];
  cout << a_label << ": " << a_values.size() << " bins -- "
       << "chi2/ndf min=" << a_values.front() << " max=" << a_values.back()
       << " mean=" << mean << " median=" << median << endl;

  double percentiles[] = {0.50, 0.75, 0.90, 0.95, 0.99};
  cout << a_label << ": percentiles -- ";
  for(size_t pIdx = 0; pIdx < sizeof(percentiles)/sizeof(double); pIdx++){
    size_t idx = (size_t)(percentiles[pIdx] * (a_values.size()-1));
    cout << "p" << (int)(percentiles[pIdx]*100) << "=" << a_values[idx] << "  ";
  }
  cout << endl;

  double candidateCuts[] = {3.0, 5.0, 7.0, 10.0, 15.0, 20.0};
  cout << a_label << ": bins surviving candidate a_maxChiSqrNdf cuts -- ";
  for(size_t cIdx = 0; cIdx < sizeof(candidateCuts)/sizeof(double); cIdx++){
    int nSurviving = 0;
    for(size_t i = 0; i < a_values.size(); i++) if(a_values[i] <= candidateCuts[cIdx]) nSurviving++;
    cout << candidateCuts[cIdx] << "->" << nSurviving << "/" << a_values.size()
         << " (" << Form("%.0f",100.0*nSurviving/a_values.size()) << "%)  ";
  }
  cout << endl;
}

// 2026-07-20: significance (content/error) distribution, same shape as
// reportChiSqrDistribution above but for the OTHER cut PresentZFitterSpectra.C applies
// (a_minSignificance, suppressInsignificantBins() in that file -- default 5.0, already
// active). Added after the chi2/ndf tail check came back clean (tail bins already pass
// even a chi2/ndf<=3 cut at ~100%), to check whether the "lower right" messiness Andrew
// reported instead correlates with significance values sitting close to that 5-sigma
// floor -- i.e. bins that clear the cut but only just, and are still statistically noisy.
// Direction is reversed vs the chi2 helper: HIGHER significance is better, so "surviving"
// here means value >= candidate cut, not <=.
void reportSignificanceDistribution(string a_label, vector<double> a_values){
  if(a_values.empty()){
    cout << a_label << ": no bins in this subset -- nothing to summarize." << endl;
    return;
  }
  sort(a_values.begin(),a_values.end());
  double sum = 0;
  for(size_t i = 0; i < a_values.size(); i++) sum += a_values[i];
  double mean = sum / a_values.size();
  double median = a_values[a_values.size()/2];
  cout << a_label << ": " << a_values.size() << " bins -- "
       << "significance min=" << a_values.front() << " max=" << a_values.back()
       << " mean=" << mean << " median=" << median << endl;

  double percentiles[] = {0.05, 0.10, 0.25, 0.50, 0.75};
  cout << a_label << ": percentiles -- ";
  for(size_t pIdx = 0; pIdx < sizeof(percentiles)/sizeof(double); pIdx++){
    size_t idx = (size_t)(percentiles[pIdx] * (a_values.size()-1));
    cout << "p" << (int)(percentiles[pIdx]*100) << "=" << a_values[idx] << "  ";
  }
  cout << endl;

  double candidateCuts[] = {5.0, 7.0, 10.0, 15.0, 20.0};
  cout << a_label << ": bins surviving candidate a_minSignificance cuts -- ";
  for(size_t cIdx = 0; cIdx < sizeof(candidateCuts)/sizeof(double); cIdx++){
    int nSurviving = 0;
    for(size_t i = 0; i < a_values.size(); i++) if(a_values[i] >= candidateCuts[cIdx]) nSurviving++;
    cout << candidateCuts[cIdx] << "->" << nSurviving << "/" << a_values.size()
         << " (" << Form("%.0f",100.0*nSurviving/a_values.size()) << "%)  ";
  }
  cout << endl;
}

void reportChiSqrStats(TFile* a_file, string a_yieldPath, string a_chiSqrPath, string a_label, int a_maxPrint = 20){
  if(!a_file){ cout << a_label << ": file not open" << endl; return; }
  TH2* yieldHist = (TH2*) a_file->Get(a_yieldPath.c_str());
  if(!yieldHist){ cout << a_label << ": yield histogram not found at " << a_yieldPath << endl; return; }
  TH2* chiSqrHist = (TH2*) a_file->Get(a_chiSqrPath.c_str());
  if(!chiSqrHist){ cout << a_label << ": chi^2/ndf histogram not found at " << a_chiSqrPath << endl; return; }

  vector<double> chiSqrValues;
  // 2026-07-20: separately tracked for bins in the top 30% of the mT-m0 (Y) axis range --
  // added after Andrew reported the high-mTm0/low-yield corner of PresentZFitterSpectra's
  // plot ("lower right") still looked messy even after a first chi^2/ndf cut, to check
  // whether that corner's mess actually correlates with worse chi^2/ndf (in which case a
  // tighter global cut would help it) or not (in which case the messiness has some other
  // cause, e.g. low-statistics scatter that clears the significance cut anyway).
  vector<double> tailChiSqrValues;
  int tailMtBinStart = (int)(0.7 * yieldHist->GetNbinsY()) + 1;

  // 2026-07-20: significance (content/error) tracked alongside chi2/ndf, same
  // full-set/tail-only split -- see reportSignificanceDistribution's header comment.
  vector<double> significanceValues;
  vector<double> tailSignificanceValues;

  int nPrinted = 0;
  for(int rapBin = 1; rapBin <= yieldHist->GetNbinsX(); rapBin++){
    for(int mtBin = 1; mtBin <= yieldHist->GetNbinsY(); mtBin++){
      double yieldVal = yieldHist->GetBinContent(rapBin,mtBin);
      if(yieldVal <= 0) continue; // only bins that actually have a fit result worth checking
      double chiSqrNdf = chiSqrHist->GetBinContent(rapBin,mtBin);
      chiSqrValues.push_back(chiSqrNdf);
      if(mtBin >= tailMtBinStart) tailChiSqrValues.push_back(chiSqrNdf);

      // Matches suppressInsignificantBins()'s own convention in PresentZFitterSpectra.C:
      // significance = content/error, or 0 if error<=0 (rather than dividing by zero).
      double yieldErr = yieldHist->GetBinError(rapBin,mtBin);
      double significance = (yieldErr > 0) ? (yieldVal/yieldErr) : 0.0;
      significanceValues.push_back(significance);
      if(mtBin >= tailMtBinStart) tailSignificanceValues.push_back(significance);

      if(nPrinted < a_maxPrint){
        cout << a_label << ": rapBin=" << rapBin << " mtm0Bin=" << mtBin
             << "  yield=" << yieldVal << "  chi2/ndf=" << chiSqrNdf
             << "  significance=" << significance << endl;
        nPrinted++;
      }
    }
  }
  if(chiSqrValues.empty()){
    cout << a_label << ": no bins with positive yield found -- nothing to summarize." << endl;
    return;
  }
  reportChiSqrDistribution(a_label, chiSqrValues);
  cout << a_label << " -- TOP 30% mT-m0 BINS ONLY (mtm0Bin >= " << tailMtBinStart << ", the \"lower right\" tail):" << endl;
  reportChiSqrDistribution(a_label+" (tail)", tailChiSqrValues);
  reportSignificanceDistribution(a_label, significanceValues);
  reportSignificanceDistribution(a_label+" (tail)", tailSignificanceValues);
}

void CheckChiSqrValues(string a_tpcFile, string a_btofFile = "", string a_speciesName = "Pion", int a_centIndex = 0){
  TFile* tpcFile  = (a_tpcFile  != "") ? new TFile(a_tpcFile.c_str(), "READ")  : NULL;
  TFile* btofFile = (a_btofFile != "") ? new TFile(a_btofFile.c_str(),"READ")  : NULL;
  string cent = Form("Cent%02d",a_centIndex);

  if(tpcFile){
    reportChiSqrStats(tpcFile,
      Form("DeDx_FitData/%s/dEdxSpectra_%sPlus_%s", a_speciesName.c_str(),a_speciesName.c_str(),cent.c_str()),
      Form("DeDx_FitData/%s/ZTPC_ChiSqr_%sPlus_%s",  a_speciesName.c_str(),a_speciesName.c_str(),cent.c_str()),
      a_speciesName+" TPC Plus");
    reportChiSqrStats(tpcFile,
      Form("DeDx_FitData/%s/dEdxSpectra_%sMinus_%s", a_speciesName.c_str(),a_speciesName.c_str(),cent.c_str()),
      Form("DeDx_FitData/%s/ZTPC_ChiSqr_%sMinus_%s",  a_speciesName.c_str(),a_speciesName.c_str(),cent.c_str()),
      a_speciesName+" TPC Minus");
  }
  if(btofFile){
    reportChiSqrStats(btofFile,
      Form("BTOF_FitData/%s/InvBetaBTOFSpectra_%sPlus_%s", a_speciesName.c_str(),a_speciesName.c_str(),cent.c_str()),
      Form("BTOF_FitData/%s/ZbTOF_ChiSqr_%sPlus_%s",        a_speciesName.c_str(),a_speciesName.c_str(),cent.c_str()),
      a_speciesName+" BTOF Plus");
    reportChiSqrStats(btofFile,
      Form("BTOF_FitData/%s/InvBetaBTOFSpectra_%sMinus_%s", a_speciesName.c_str(),a_speciesName.c_str(),cent.c_str()),
      Form("BTOF_FitData/%s/ZbTOF_ChiSqr_%sMinus_%s",        a_speciesName.c_str(),a_speciesName.c_str(),cent.c_str()),
      a_speciesName+" BTOF Minus");
  }
}
