// TestErrorBarDraw2.C -- temporary diagnostic (2026-07-20), safe to delete after use.
//
// Follow-up to TestErrorBarDraw.C: that test showed "P" and "PE" render IDENTICALLY
// (both draw the error bar) for a plain histogram on a plain linear-scale canvas with no
// pre-existing frame. But the real PresentZFitterSpectra.C draws onto a LOG-Y pad with a
// frame created via gPad->DrawFrame() first (see lines ~619, ~846), then Draw("...SAME")
// on top of that frame -- and in that real macro, NEITHER "P" nor "PE" show any error
// bar at all for the points actually checked. This test reproduces that exact
// combination (SetLogy + DrawFrame + SAME) to see if THAT specific combination is what
// suppresses the error bar, isolating it from the plain-canvas case that worked fine.
void TestErrorBarDraw2(){
  TCanvas* c3 = new TCanvas("c3","c3 -- SetLogy + DrawFrame + PE SAME",600,600);
  gPad->SetLogy();
  TH1F* frame = gPad->DrawFrame(0,0.1,5,100);
  frame->GetXaxis()->SetTitle("x");
  frame->GetYaxis()->SetTitle("y (log)");
  frame->Draw();

  TH1D* h3 = new TH1D("h3","h3",5,0,5);
  h3->SetBinContent(3,10);
  h3->SetBinError(3,2); // 20% error, same as before
  h3->SetMarkerStyle(20);
  h3->SetMarkerSize(1.2);
  h3->SetMarkerColor(kBlue);
  h3->SetLineColor(kBlue);
  h3->Draw("PE SAME");

  c3->Print("test_draw_PE_logy_frame_same.png");
  cout << "h3 bin 3: content=" << h3->GetBinContent(3) << " error=" << h3->GetBinError(3) << endl;
  cout << "Wrote test_draw_PE_logy_frame_same.png" << endl;
}
