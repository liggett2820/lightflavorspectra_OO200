// TestErrorBarDraw.C -- temporary diagnostic (2026-07-20), safe to delete after use.
//
// Isolates whether "PE" vs "P" draw options actually produce different PNG output in
// this ROOT environment, independent of PresentZFitterSpectra.C's much larger
// canvas/frame/pad/log-scale setup. Context: switching that macro's Draw("P SAME") to
// Draw("PE SAME") produced a PNG that's byte-for-byte identical to the "P SAME" version,
// even after confirming (via a temporary debug print right before the Draw call) that
// the histograms being drawn carry real, nonzero per-bin errors (several percent
// relative error) at exactly the point they're drawn. This rules out "the errors are
// getting lost somewhere upstream" and narrows it to something about how "PE" renders
// (or doesn't) in this specific setup. This macro reproduces the same Draw("P")/
// Draw("PE") contrast on a minimal single-bin histogram with an obviously large,
// unmissable error bar (content=10, error=2 -- 20%), so there's no ambiguity about
// whether an error bar SHOULD be visible.
void TestErrorBarDraw(){
  TH1D* h1 = new TH1D("h1","h1 -- Draw(\"P\")",5,0,5);
  h1->SetBinContent(3,10);
  h1->SetBinError(3,2);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1.2);
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  h1->Draw("P");
  c1->Print("test_draw_P.png");

  TH1D* h2 = new TH1D("h2","h2 -- Draw(\"PE\")",5,0,5);
  h2->SetBinContent(3,10);
  h2->SetBinError(3,2);
  h2->SetMarkerStyle(20);
  h2->SetMarkerSize(1.2);
  TCanvas* c2 = new TCanvas("c2","c2",600,600);
  h2->Draw("PE");
  c2->Print("test_draw_PE.png");

  cout << "h1 (P) bin 3: content=" << h1->GetBinContent(3) << " error=" << h1->GetBinError(3) << endl;
  cout << "h2 (PE) bin 3: content=" << h2->GetBinContent(3) << " error=" << h2->GetBinError(3) << endl;
  cout << "ROOT version: " << gROOT->GetVersion() << endl;
  cout << "Wrote test_draw_P.png and test_draw_PE.png in the current directory." << endl;
}
