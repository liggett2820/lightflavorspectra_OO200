// CheckLastRapBin.C -- lightflavorspectra_OO200
//
// Quick diagnostic: for a given ZFitter TPC/BTOF output file pair, report the last
// (highest rapidity-bin index, i.e. highest |y|) bin that actually has any positive
// yield in it, for both charges. Answers "what's the last filled rapidity bin on the
// spectra?" directly from the real fit output, rather than assuming the analysis's
// nominal binning (41 bins, -2.05 to 2.05, so bin 41 = y=2.00 -- see SetCutClass.C's
// setVariableRapMtM0BinningInfo()) is actually populated with real fit results all the
// way out to that edge. A rapidity bin can easily be "booked" by the binning but never
// get a successful/nonzero fit if there's just not enough statistics that far out.

void reportLastFilledBin(TFile* a_file, string a_path, string a_label){
  if(!a_file){ cout << a_label << ": file not open" << endl; return; }
  TH2* h = (TH2*) a_file->Get(a_path.c_str());
  if(!h){ cout << a_label << ": histogram not found at " << a_path << endl; return; }
  // Sum bin content directly instead of calling ProjectionY() -- ProjectionY() attaches
  // its new histogram to whatever the CURRENT directory is (here, the still-open input
  // TFile), and an explicit delete right after fighting over ownership with ROOT's own
  // directory cleanup is exactly the kind of thing that segfaults (confirmed: this is
  // what crashed the first version of this macro). Summing GetBinContent() directly
  // never allocates anything, so there's nothing to manage or delete.
  int lastFilled = -1;
  for(int bx = 1; bx <= h->GetNbinsX(); bx++){
    double integral = 0;
    for(int by = 1; by <= h->GetNbinsY(); by++) integral += h->GetBinContent(bx,by);
    if(integral > 0) lastFilled = bx;
  }
  if(lastFilled < 0){
    cout << a_label << ": no filled rapidity bins at all (checked all " << h->GetNbinsX() << ")" << endl;
    return;
  }
  cout << a_label << ": last filled rapidity bin = " << lastFilled
       << "  (y=" << h->GetXaxis()->GetBinCenter(lastFilled)
       << ", edges [" << h->GetXaxis()->GetBinLowEdge(lastFilled) << ","
       << h->GetXaxis()->GetBinUpEdge(lastFilled) << "])"
       << "  out of " << h->GetNbinsX() << " total booked rapidity bins"
       << " (nominal last bin would be y=" << h->GetXaxis()->GetBinCenter(h->GetNbinsX()) << ")"
       << endl;
}

void CheckLastRapBin(string a_tpcFile, string a_btofFile = "", string a_speciesName = "Pion", int a_centIndex = 0){
  TFile* tpcFile  = (a_tpcFile  != "") ? new TFile(a_tpcFile.c_str(), "READ")  : NULL;
  TFile* btofFile = (a_btofFile != "") ? new TFile(a_btofFile.c_str(),"READ")  : NULL;
  string cent = Form("Cent%02d",a_centIndex);

  if(tpcFile){
    reportLastFilledBin(tpcFile, Form("DeDx_FitData/%s/dEdxSpectra_%sPlus_%s", a_speciesName.c_str(),a_speciesName.c_str(),cent.c_str()),  a_speciesName+" TPC Plus");
    reportLastFilledBin(tpcFile, Form("DeDx_FitData/%s/dEdxSpectra_%sMinus_%s",a_speciesName.c_str(),a_speciesName.c_str(),cent.c_str()),  a_speciesName+" TPC Minus");
  }
  if(btofFile){
    reportLastFilledBin(btofFile, Form("BTOF_FitData/%s/InvBetaBTOFSpectra_%sPlus_%s", a_speciesName.c_str(),a_speciesName.c_str(),cent.c_str()), a_speciesName+" BTOF Plus");
    reportLastFilledBin(btofFile, Form("BTOF_FitData/%s/InvBetaBTOFSpectra_%sMinus_%s",a_speciesName.c_str(),a_speciesName.c_str(),cent.c_str()), a_speciesName+" BTOF Minus");
  }
}
