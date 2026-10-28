# STAR Glauber + NBD centrality for O+O 200 GeV (Run 21: refMult and EPD totnMIP)

This is STAR's own centrality package (`star-sw/StRoot/PWGTools/CentralityCalibration/Glauber/star_glauber`, last upstream change 9319c7f9be, Apr 2023) with O+O added. Everything else in the chain (StFastGlauberMcMaker, StNbdFitMaker, StGlauberAnalysisMaker, the condor scripts) is the upstream code and runs exactly as the upstream README describes. `star_glauber_OO200.patch` is the full diff against upstream.

## What was added or changed

**O+O system.** `StFastGlauberMcMaker` gets `InitOO()` (system string `OO`). The settings are:

- A = 16, spherical.
- σ_NN = 42 mb (upstream 200 GeV value; `largeXsec`/`smallXsec` give ±1 mb).
- Hard-core NN collisions, b sampled up to 20 fm, as upstream.

**¹⁶O density.** The density is a three-parameter Fermi (3pF): ρ(r) ∝ (1 + w r²/R²)/(1 + exp((r−R)/d)), with R = 2.608 fm, d = 0.513 fm and w = −0.051. These are the De Vries, De Jager & De Vries values (At. Data Nucl. Data Tables 36, 495 (1987)), the same ¹⁶O parameterization TGlauberMC uses (Loizides, Kamin & d'Enterria, PRC 97, 054910 (2018)). Its rms radius is 2.73 fm.

STAR's code only had 2pF Woods–Saxon, so `GlauberUtilities::WoodsSaxon3pF` was added. Setting w = 0 reproduces the old function, so no other system changes.

**Systematic types.** `large`/`small` follow the convention already documented in `Print("type")` and `submit_glauber.pl`: R ±2% with d ∓10%. For a nucleus this light, d dominates, so `large` is the more compact one (rms 2.64 fm) and `small` the more extended one (2.82 fm).

**StCentrality.** `Init_OO200GeV()` (system `OO_200GeV`) holds the raw-refMult NBD fit and edges from Step 3 (filled Oct 7 2026, commit 7ab3afa). The analysis now uses the EPD totnMIP classes of Step 5, which are applied in the main repo (`headers/EpdCentrality_OO200.h`, `macros/SetCutClass.C`), not through StCentrality.

**Centrality edge convention.** `StNbdFitMaker` prints, for each class, the lowest refMult *included* in it (class = refMult ≥ printed value). This is also what StRefMultCorr's `> cut` selects after its [0,1) smearing. `StCentrality::GetCentrality()` uses a strict `>` on the integer MC multiplicity, so `Init_OO200GeV()` subtracts 1 from the pasted values. This is the same "− 1" as the Run 4 Au+Au table in that file.

**What this means in your analysis.** Select class i as `refMult >= cut_i`.

**Npart sampling.** `StNbdFitMaker::SetIntegerGlauberSampling()` is new and off by default. `TH2::GetRandom2` on the unit-width `hNcoll_Npart` returns Npart + U(0,1). That is negligible for Au+Au but is up to +25% at Npart = 2. The O+O macro `doNbdFitMaker_OO200.C` turns it on.

**Script changes.**

- `doFastGlauberMcMaker.csh`, `doAnalysisMaker.csh` and `doPlotMaker.csh` use `$GLAUBER_STARVER` when set and default to SL24y (as PicoBinner). Upstream mixes SL19b, SL16d and stardev. `submit_condor.pl` runs each condor job inside the SDCC SL7 container.
- `submit_condor.pl` sends condor error mail to `$USER` instead of the original author.

**New O+O files:**

- `build_OO200.csh`
- `all_submit_OO200.csh`
- `createList_OO200.csh`
- `make_hRefMult_OO200.C`
- `doNbdFitMaker_OO200.C`
- `doScanX_OO200.csh`
- `submit_doScan_OO200.pl`
- `find_best_nbd_OO200.py`

`find_best_nbd_OO200.py` replaces `getBestChi2_fromCat/*.sh`, which read fixed log line numbers and character columns.

## Run it on RCF

There are two kinds of commands. Anything that runs `starver`, `cons` or `root4star` (the build, the interactive test, `addNcollVsNpart.C`) runs **inside** the SDCC SL7 container, the same one PicoBinner uses:

```bash
singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif /bin/tcsh
```

The submission scripts (`all_submit_OO200.csh`, `submit_doScan_OO200.pl`, `all_doAnalysisMaker.csh`) only call `condor_submit`, so run them **outside** the container, on the host. `submit_condor.pl` wraps every condor job in that same container image, so the jobs run inside it.

The STAR version defaults to SL24y, as for PicoBinner. To use another, `setenv GLAUBER_STARVER <version>` in both shells; the variable is passed into the jobs.

Inside the container:

```csh
cd <repo>/centrality/star_glauber
./prepare.sh                      # first set outDir to your PWG disk
./build_OO200.csh                 # syncs Makers/St* into StRoot/ and runs cons
./doFastGlauberMcMaker.csh output/test_OO.root 1000 OO 200 default kFALSE   # quick test
```

On the host, for the submissions:

```csh
cd <repo>/centrality/star_glauber
set path = ( . $path )            # submit_glauber.pl calls submit_condor.pl without ./
```

### Step 1: Glauber trees

```csh
./all_submit_OO200.csh 0 1 0      # dry run: prints the condor commands
./all_submit_OO200.csh 0 10 1     # 10 jobs x 200k per type (jobs 0..9)
```

Check one output first, `output/fastglaubermc_OO_200GeV_default_spherical_run0000.root`. `hWoodsSaxon_0` (the sampled density ρ(r)) should be flat in the core and drop to half its central value near r = 2.6 fm, and the tree header should show A = 16 and σ_NN = 42 mb.

```csh
./createList_OO200.csh
root4star -b -q -l addNcollVsNpart.C      # inside the container -> ncoll_npart.root
```

### Step 2: NBD fit to raw refMult

Make the data input from the PicoBinner output (plain ROOT; works on the laptop too):

```csh
root -l -b -q 'make_hRefMult_OO200.C("yieldHistos_OO200_pion.root")'   # -> hRefMult_OO200.root
```

The input is PicoBinner's `refMult` after the event cuts (trigger 860003, |Vz| ≤ 2 cm, Vr ≤ 1 cm). Use the histogram from your current production.

Scan, then rank:

```csh
./submit_doScan_OO200.pl          # prints the 180 jobs
./submit_doScan_OO200.pl -run     # submits them
python3 find_best_nbd_OO200.py LOG_Scan
```

The coarse grid is:

| Parameter | Range | Points |
|---|---|---|
| npp | 1.6–3.2 | 9 |
| k | 1–9 | 5 |
| x | 0.08–0.28 | 6 |
| d | 0–0.12 | 4 |

The fit uses refMult ≥ 15 (`multCut` in `doScanX_OO200.csh`). Repeat with 10 and 20 as a cross-check.

`find_best_nbd_OO200.py` flags a best point on the edge of the grid. Narrow the ranges around the minimum and rescan before trusting it.

d = 0 means a constant 98% efficiency. That is consistent with the flat TPC efficiency vs centrality seen in the embedding, but let the scan decide.

Then check the best point's `RatioChi2Files/Ratio_*.root`. `hRatio` (MC/data) should be flat at 1 above `multCut`, and the deviation below it is the trigger/vertex inefficiency. Follow the upstream README's reweight-and-refit check if you want the flat-ratio closure plot (`fit()` in `doNbdFitMaker_OO200.C` does a single fixed-parameter fit).

### Step 3: Centrality table, Npart/Ncoll and systematics

1. In `Makers/StCentralityMaker/StCentrality.cxx` → `Init_OO200GeV()`:
   - put in the best npp, k, x and d;
   - paste the 32 `mMultiplicityCut[...]` lines that `find_best_nbd_OO200.py` printed, between the PASTE markers, replacing the placeholder lines and the placeholder loop;
   - set `kIsPlaceholder = kFALSE`.
2. Run `./build_OO200.csh`, then regenerate the trees (Step 1), as the upstream README requires.
3. Run the analysis and tables:
   ```csh
   ./createList_OO200.csh
   ./all_doAnalysisMaker.csh OO_200GeV kFALSE kFALSE
   ./all_doPlotMaker.csh 200        # inside the container (runs root4star directly); tables in ./table, figures in ./figure
   ```

The tables cover 5% classes and the wider ones in `StGlauberConstUtilities`: 0–10, 10–20, 20–40, 40–60, 60–80, 40–80 and others. The systematic error per class is the maximum deviation over the types.

### Step 4: refMult corrections and trigger efficiency (Oct 2026)

The fits above use raw refMult with |Vz| <= 2 cm. The StRefMultCorr-style corrections for
trigger 860003 (Vz, luminosity, pile-up; see the main README, "StRefMultCorr-style
centrality") produce a corrected refMult histogram, `hRefMultCorr_OO200.root` (TH1D
`hRefMult`, unit bins from 0), for |Vz| <= 30 cm. To rescan on it, copy it to the name the
scan reads (`cp hRefMultCorr_OO200.root hRefMult_OO200.root`, after keeping the old file),
rerun `submit_doScan_OO200.pl` and `find_best_nbd_OO200.py`, and refill `Init_OO200GeV()`.

Trigger/vertex efficiency below the fit threshold, from the best point's Ratio file:
```csh
root -l -b -q 'FitTriggerEfficiency_OO200.C+("Ratio_best_mc15/Ratio_npp2.750_k2.500_x0.170_eff0.060.root", 15)'
```
It fits data/MC with 1 - exp(-p0 m^p1) (the form `StCentrality::GetReweighting` uses; the
lowrw/highrw +-2 sigma values are printed) and with 1 - p0 exp(-p1 m^p2), which allows a
nonzero efficiency at m = 0, and writes `OO200_TrigEff_params.txt` and `OO200_TrigEff.png`.
To use the weights in the Glauber tables, put the form-(1 - exp) parameters into
`Init_OO200GeV()` (`mParReweighting`) and run `all_doAnalysisMaker.csh OO_200GeV kTRUE kFALSE`
with `lowrw`/`highrw` added to its list.

### Step 5: EPD totnMIP centrality (Oct 2026, in use)

This is the centrality the spectra use. refMult cannot resolve O+O beyond about 75%, while the
EPD (2.1 < |η| < 5.1) can, and it is away from the mid-rapidity tracks the spectra are
measured with.

**Input.** The input is corrected totnMIP for trigger 860003 at |Vz| ≤ 30 cm, from the main
repo:
1. `macros/MakeEpdCentInputs_OO200.C` (SDCC job `xml/runEpdCentInputs_OO200_SDCC.xml`).
2. `macros/FitEpdCent_OO200.C`, which writes `hTotnMIPCorr_OO200.root`.

Copy that file here before scanning.

**1. Scan, on the host.** `-totnmip` reads `hTotnMIPCorr_OO200.root` and fits totnMIP ≥ 200;
logs go to `LOG_Scan_totnmip_mc200`. Run the coarse pass first, then the fine pass around its
minimum:
```csh
./submit_doScan_OO200.pl -totnmip -run
python3 find_best_nbd_OO200.py LOG_Scan_totnmip_mc200
./submit_doScan_OO200.pl -totnmip -fine -run
python3 find_best_nbd_OO200.py LOG_Scan_totnmip_mc200_fine
```

**2. Fit-range check.** Add `-multcut 150` or `-multcut 250` to the scan command.

**3. Glauber variations.** First make `ncoll_npart_<TYPE>.root` with `addNcollVsNpart.C`
for each TYPE (small, large, smallXsec, largeXsec, gauss). Then run
`./submit_doScan_OO200.pl -glauber TYPE -run` for each one. Each variation runs in its own
`glauber_TYPE/` directory.

**4. Cuts, 0–100%.** `CentCuts_OO200.C` replaces `StNbdFitMaker::CalculateCentrality`, which
stops at 80%:
```csh
root -l -b -q 'CentCuts_OO200.C("RatioChi2Files/Ratio_<best point>.root", 200)'
```
It takes the data above the fit range and the MC below it. It prints the class edges for
N_tot and for N_tot × 0.95 and × 1.05, and the 860003 efficiency (data/MC) vs totnMIP and per
class. Output: `OO200_NBD_totnMIP_cuts.txt`, `OO200_NBD_totnMIP_fit.png`.

**5. ⟨Npart⟩, ⟨Ncoll⟩, ⟨Nspec⟩ per class**, nominal and for each Glauber variation. This runs
inside the container:
```csh
root -l -b -q 'NpartPerClass_OO200.C+'
```

**Result.** The fine scan's best point is npp = 24, k = 6, x = 0.20, d = 0, with χ²/ndf = 0.94
for totnMIP ≥ 200. This gives edges 364/310/234/130/36/0 for 0–5, 5–10, 10–20, 20–40,
40–80 and 80–100%.

**Fit-range check result.** Refits to totnMIP ≥ 150 and ≥ 250 give 365/311/235/130/36 and
364/310/233/130/36. The per-class 860003 efficiencies agree within 0.02.

**Comparison with the official table.** The official StRefMultCorr totnMIP edges (triggers
860001/2/11/12) agree with this fit's −5% σ_tot edges within 3 units for 0–75%.

**Cross-checks** in the main repo: `macros/PlotTrigRatio_OO200.C` (860003 relative to 860001|860002
vs totnMIP) and `macros/CompareVzCent_OO200.C` (class fractions at |Vz| ≤ 30 vs ≤ 2 cm).

**Open:**
- Below 50%, the 860003 efficiency from the fit is 4–21% below the PlotTrigRatio-based
  estimate.
- The refit with the official normalisation is still to do.
- Approval by the centrality group.

**Note on the NBD model.** StNbdFitMaker generates integer NBD multiplicities, whereas totnMIP
is a continuous sum of truncated nMIP. The fit quality above suggests the NBD describes it
over the fit range.

## Things to keep in mind for O+O

**Multiplicity granularity.** refMult is a small integer here, so a 5% class cannot be exactly 5%. The printed edge is the integer whose cumulative fraction is closest. Quote the actual fraction each class selects, from the data histogram above the fit range and from the MC below it.

**refMult corrections.** STAR's procedure expects luminosity- and Vz-corrected refMult. The raw-refMult fit in Steps 2–3 uses |Vz| ≤ 2 cm, where the Vz dependence is small; Step 4 derives the corrections for |Vz| ≤ 30 cm. The analysis itself now uses the EPD totnMIP classes of Step 5.

**Nuclear structure.** ¹⁶O is generated spherical with uncorrelated nucleons. α-clustering and short-range correlations are not in STAR's model.
