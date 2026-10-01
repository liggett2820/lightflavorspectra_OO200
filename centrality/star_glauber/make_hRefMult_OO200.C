//____________________________________________________________________________________________________
// make_hRefMult_OO200.C -- data input for the STAR NBD x Glauber fit (StNbdFitMaker), O+O 200 GeV
//
// Reads the raw refMult histogram that PicoBinner books after the event cuts (trigger 860003,
// |Vz| <= 2 cm, Vr <= 1 cm, ...) -- TH1I "refMult", 1000 unit bins from 0 -- and writes it as
// TH1D "hRefMult", the name and type StNbdFitMaker::ReadData() expects.
//
// Binning requirement (StNbdFitMaker::CalculateCentrality): bin i+1 must hold refMult = i,
// i.e. unit-width bins starting at 0. Checked below.
//
// Plain ROOT, no STAR libraries:
//   root -l -b -q 'make_hRefMult_OO200.C("yieldHistos_OO200_pion.root")'
//____________________________________________________________________________________________________
void make_hRefMult_OO200(
    const char* inputFile  = "yieldHistos_OO200_pion.root",
    const char* inputName  = "refMult",
    const char* outputFile = "hRefMult_OO200.root"
)
{
  TFile* fin = TFile::Open(inputFile, "READ");
  if(!fin || fin->IsZombie()){ Error("make_hRefMult_OO200", "cannot open %s", inputFile); return; }

  TH1* hin = dynamic_cast<TH1*>(fin->Get(inputName));
  if(!hin){ Error("make_hRefMult_OO200", "no histogram '%s' in %s", inputName, inputFile); return; }

  const Int_t    nbins = hin->GetNbinsX();
  const Double_t xmin  = hin->GetXaxis()->GetXmin();
  const Double_t xmax  = hin->GetXaxis()->GetXmax();
  if( xmin != 0.0 || TMath::Abs((xmax - xmin)/nbins - 1.0) > 1e-9 ){
    Error("make_hRefMult_OO200", "'%s' must have unit bins starting at 0 (has %d bins in [%g, %g])",
          inputName, nbins, xmin, xmax);
    return;
  }

  TH1D* h = new TH1D("hRefMult", "O+O 200 GeV, raw refMult after event cuts;refMult;Events", nbins, xmin, xmax);
  for(Int_t i=1; i<=nbins; i++){
    h->SetBinContent(i, hin->GetBinContent(i));
    h->SetBinError(i, TMath::Sqrt(hin->GetBinContent(i))); // raw event counts
  }
  h->SetEntries(hin->GetEntries());

  Int_t lastFilled = 0;
  for(Int_t i=nbins; i>=1; i--) if(h->GetBinContent(i) > 0){ lastFilled = i; break; }
  cout << "hRefMult: " << (Long64_t)h->GetEntries() << " events, mean refMult = " << h->GetMean()
       << ", underflow = " << hin->GetBinContent(0) << ", overflow = " << hin->GetBinContent(nbins+1)
       << ", highest filled refMult = " << lastFilled-1 << endl;

  TFile* fout = TFile::Open(outputFile, "RECREATE");
  h->Write();
  fout->Close();
  fin->Close();
  cout << "Wrote " << outputFile << endl;
}
