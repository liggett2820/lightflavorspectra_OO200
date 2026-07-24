// TestErrorBarDraw4.C -- temporary diagnostic (2026-07-20), safe to delete after use.
//
// TestErrorBarDraw2.C (fresh histogram, content=10/error=2, SetLogy+DrawFrame+PE SAME)
// showed a bar. TestErrorBarDraw3.C (same draw pattern, but content=0.004/error=0.0004
// AND passed through the real combineMirroredProjections() clone-then-overwrite
// pipeline) showed no bar, even though the printed content/error going into Draw() were
// confirmed correct. Two things changed at once between those two tests: the value's
// magnitude, and whether it went through combineMirroredProjections(). This isolates
// which one actually matters by testing all 4 combinations.
void TestErrorBarDraw4(){
  gROOT->ProcessLine(".L PresentZFitterSpectra.C");

  // Case A: small value (like the real BTOF tail bin), FRESH histogram, no clone/combine.
  TH1D* hA = new TH1D("hA","hA",5,0,5);
  hA->SetBinContent(3,0.00408945);
  hA->SetBinError(3,0.000382926);
  TCanvas* cA = new TCanvas("cA","cA -- small value, fresh histogram",600,600);
  gPad->SetLogy();
  TH1F* frA = gPad->DrawFrame(0,1e-5,5,1);
  frA->Draw();
  hA->SetMarkerStyle(24); hA->SetMarkerSize(0.7); hA->SetMarkerColor(kBlue); hA->SetLineColor(kBlue);
  hA->Draw("PE SAME");
  cA->Print("test_case_A_small_fresh.png");

  // Case B: large value (like the earlier working test), through combineMirroredProjections.
  TH1D* projI_B = new TH1D("projI_B","projI_B",5,0,5);
  projI_B->SetBinContent(3,10);
  projI_B->SetBinError(3,2);
  TH1D* projJ_B = new TH1D("projJ_B","projJ_B",5,0,5);
  TH1D* hB = (TH1D*) gROOT->ProcessLine(Form("combineMirroredProjections((TH1D*)%p,(TH1D*)%p,\"hB\")",(void*)projI_B,(void*)projJ_B));
  cout << "hB bin 3: content=" << hB->GetBinContent(3) << " error=" << hB->GetBinError(3) << endl;
  TCanvas* cB = new TCanvas("cB","cB -- large value, through combineMirroredProjections",600,600);
  gPad->SetLogy();
  TH1F* frB = gPad->DrawFrame(0,0.1,5,100);
  frB->Draw();
  hB->SetMarkerStyle(24); hB->SetMarkerSize(0.7); hB->SetMarkerColor(kBlue); hB->SetLineColor(kBlue);
  hB->Draw("PE SAME");
  cB->Print("test_case_B_large_combined.png");

  // Case C: small value, through combineMirroredProjections (repeat of TestErrorBarDraw3
  // for completeness, same-session so directly comparable).
  TH1D* projI_C = new TH1D("projI_C","projI_C",5,0,5);
  projI_C->SetBinContent(3,0.00408945);
  projI_C->SetBinError(3,0.000382926);
  TH1D* projJ_C = new TH1D("projJ_C","projJ_C",5,0,5);
  TH1D* hC = (TH1D*) gROOT->ProcessLine(Form("combineMirroredProjections((TH1D*)%p,(TH1D*)%p,\"hC\")",(void*)projI_C,(void*)projJ_C));
  cout << "hC bin 3: content=" << hC->GetBinContent(3) << " error=" << hC->GetBinError(3) << endl;
  TCanvas* cC = new TCanvas("cC","cC -- small value, through combineMirroredProjections",600,600);
  gPad->SetLogy();
  TH1F* frC = gPad->DrawFrame(0,1e-5,5,1);
  frC->Draw();
  hC->SetMarkerStyle(24); hC->SetMarkerSize(0.7); hC->SetMarkerColor(kBlue); hC->SetLineColor(kBlue);
  hC->Draw("PE SAME");
  cC->Print("test_case_C_small_combined.png");

  // Case D: large value, FRESH histogram, no clone/combine.
  TH1D* hD = new TH1D("hD","hD",5,0,5);
  hD->SetBinContent(3,10);
  hD->SetBinError(3,2);
  TCanvas* cD = new TCanvas("cD","cD -- large value, fresh histogram",600,600);
  gPad->SetLogy();
  TH1F* frD = gPad->DrawFrame(0,0.1,5,100);
  frD->Draw();
  hD->SetMarkerStyle(24); hD->SetMarkerSize(0.7); hD->SetMarkerColor(kBlue); hD->SetLineColor(kBlue);
  hD->Draw("PE SAME");
  cD->Print("test_case_D_large_fresh.png");

  cout << "Wrote test_case_A_small_fresh.png, test_case_B_large_combined.png, test_case_C_small_combined.png, test_case_D_large_fresh.png" << endl;
}
