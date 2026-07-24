// TestErrorBarDraw3.C -- temporary diagnostic (2026-07-20), safe to delete after use.
//
// Two simpler reproductions (TestErrorBarDraw.C: plain histogram, P vs PE;
// TestErrorBarDraw2.C: SetLogy + DrawFrame + SAME) both showed a real, clearly visible
// error bar. Yet the actual PresentZFitterSpectra_Pion_Cent00.png shows none on its
// sparse high-mT-m0 tail points, even for bins confirmed (via a debug print inside
// PresentZFitterSpectra.C itself) to carry ~9% relative error at the exact point of
// drawing. This test closes the gap: it loads PresentZFitterSpectra.C's REAL
// combineMirroredProjections() function (not a reimplementation) via ".L", runs a
// synthetic bin with the same content/error as a real BTOF tail bin through it (mimicking
// the actual pipeline: raw projection -> combineMirroredProjections -> Draw), and draws
// the result with the same SetLogy+DrawFrame+PE SAME pattern. If THIS shows a bar, the
// bug is somewhere else in the real macro (e.g. interaction between many overlaid
// SAME draws); if it doesn't, combineMirroredProjections() itself is the culprit.
void TestErrorBarDraw3(){
  gROOT->ProcessLine(".L PresentZFitterSpectra.C");

  // Mimic a real BTOF tail bin (grp18 BTOF from the 2026-07-20 debug run:
  // content=0.00408945, error=0.000382926, ~9.4% relative error).
  TH1D* projI = new TH1D("projI","projI",5,0,5);
  projI->SetBinContent(3,0.00408945);
  projI->SetBinError(3,0.000382926);
  TH1D* projJ = new TH1D("projJ","projJ",5,0,5); // empty mirror partner -- projI-only case

  TH1D* combined = (TH1D*) gROOT->ProcessLine(Form("combineMirroredProjections((TH1D*)%p,(TH1D*)%p,\"combinedTest\")",(void*)projI,(void*)projJ));

  cout << "combined bin 3: content=" << combined->GetBinContent(3) << " error=" << combined->GetBinError(3) << endl;

  TCanvas* c4 = new TCanvas("c4","c4 -- real combineMirroredProjections + PE SAME",600,600);
  gPad->SetLogy();
  TH1F* frame = gPad->DrawFrame(0,1e-5,5,1);
  frame->Draw();
  combined->SetMarkerStyle(24);
  combined->SetMarkerSize(0.7);
  combined->SetMarkerColor(kBlue);
  combined->SetLineColor(kBlue);
  combined->Draw("PE SAME");
  c4->Print("test_draw_real_pipeline.png");
  cout << "Wrote test_draw_real_pipeline.png" << endl;
}
