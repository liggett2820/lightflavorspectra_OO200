// TestErrorBarDraw5.C -- temporary diagnostic (2026-07-20), safe to delete after use.
//
// TestErrorBarDraw4.C showed the deciding factor is value magnitude, not pipeline: small
// content (~0.004) never shows a bar regardless of pipeline, large content (~10) always
// does. But those two cases also used DIFFERENT axis spans (6 decades for the small
// case, matching the real plot's yFloor=1e-5; 3 decades for the large case) -- so a 9.4%
// relative error might just be genuinely too few pixels tall on a 6-decade axis to
// render (a resolution limit, not a bug), rather than small VALUES being broken outright.
// This test keeps the SAME 6-decade axis (1e-5 to 1) as the failing case, but uses a much
// larger, unmistakable 50% relative error at the same small absolute value (content=
// 0.004, error=0.002) -- if THIS shows a bar, small values render fine and the real
// plot's issue is just that most of its actual error bars are too thin to see at this
// axis's dynamic range, not a rendering bug.
void TestErrorBarDraw5(){
  TH1D* h = new TH1D("h","h",5,0,5);
  h->SetBinContent(3,0.004);
  h->SetBinError(3,0.002); // 50% relative error, same small absolute value as before
  TCanvas* c = new TCanvas("c","c -- small value, LARGE (50%) relative error, 6-decade axis",600,600);
  gPad->SetLogy();
  TH1F* fr = gPad->DrawFrame(0,1e-5,5,1);
  fr->Draw();
  h->SetMarkerStyle(24); h->SetMarkerSize(0.7); h->SetMarkerColor(kBlue); h->SetLineColor(kBlue);
  h->Draw("PE SAME");
  c->Print("test_case_E_small_bigrelerr.png");
  cout << "Wrote test_case_E_small_bigrelerr.png" << endl;
}
