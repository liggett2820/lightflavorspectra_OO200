//____________________________________________________________________________________________________
// NBD x Glauber fit for O+O 200 GeV (StNbdFitMaker), raw refMult.
// Same calls as doNbdFitMaker.C, plus SetIntegerGlauberSampling(kTRUE): (Npart, Ncoll) are drawn as
// the integer Glauber values instead of n + U(0,1) from the unit-width bins of hNcoll_Npart.
//   data histogram : "hRefMult" (made by make_hRefMult_OO200.C; unit bins [n, n+1) from 0)
//   isConstEfficiency = kFALSE : multiplicity-dependent efficiency eff = 0.98 (1 - d M/540), d = "efficiency"
//____________________________________________________________________________________________________
void scan(
    const Int_t    nevents         = 1000000,
    const Char_t*  realData        = "hRefMult_OO200.root",
    const Char_t*  glauber         = "ncoll_npart.root",
    const Double_t multiplicityCut = 15,
    const Int_t    nppbin = 1, const Double_t nppmin = 2.4,  const Double_t nppmax = 2.4,
    const Int_t    kbin   = 1, const Double_t kmin   = 2.0,  const Double_t kmax   = 2.0,
    const Int_t    xbin   = 1, const Double_t xmin   = 0.13, const Double_t xmax   = 0.13,
    const Double_t efficiency        = 0.00,
    const Double_t triggerbias       = 1.00,
    const Bool_t   isConstEfficiency = kFALSE
)
{
  gSystem->Load("St_base");
  gSystem->Load("StUtilities");
  gSystem->Load("StGlauberUtilities");
  gSystem->Load("StCentralityMaker");

  StNbdFitMaker* maker = new StNbdFitMaker();
  maker->SetMinimumMultiplicityCut(multiplicityCut);
  maker->SetIntegerGlauberSampling(kTRUE);
  maker->SetParameters(nppmin, kmin, xmin, efficiency, triggerbias, isConstEfficiency);
  maker->ReadData(realData, glauber, "hRefMult");

  maker->Scan(
      nevents,
      nppbin, nppmin, nppmax,
      kbin, kmin, kmax,
      xbin, xmin, xmax,
      efficiency, triggerbias, isConstEfficiency
      );
}

//____________________________________________________________________________________________________
// Single fit at fixed parameters (README "fix the parameters" step), writes MC/data ratio + histograms
void fit(
    const Char_t*  outputFileName  = "NbdFit_OO200_best.root",
    const Int_t    nevents         = 1000000,
    const Char_t*  realData        = "hRefMult_OO200.root",
    const Char_t*  glauber         = "ncoll_npart.root",
    const Double_t multiplicityCut = 15,
    const Double_t npp             = 2.4,
    const Double_t k               = 2.0,
    const Double_t x               = 0.13,
    const Double_t efficiency      = 0.00,
    const Double_t triggerbias     = 1.00,
    const Bool_t   isConstEfficiency = kFALSE
)
{
  gSystem->Load("St_base");
  gSystem->Load("StUtilities");
  gSystem->Load("StGlauberUtilities");
  gSystem->Load("StCentralityMaker");

  StNbdFitMaker* maker = new StNbdFitMaker();
  maker->SetMinimumMultiplicityCut(multiplicityCut);
  maker->SetIntegerGlauberSampling(kTRUE);
  maker->SetParameters(npp, k, x, efficiency, triggerbias, isConstEfficiency);
  maker->ReadData(realData, glauber, "hRefMult");
  maker->Fit(nevents, outputFileName);
}
