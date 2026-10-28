# lightflavorspectra_OO200

π±, K± and p/p̄ spectra in O+O collisions at √s_NN = 200 GeV (STAR Run 21), ColliderCenter
geometry. This is a simplified fork of `lightflavorspectra_etof`, cut down to this one system.

The work is split between two places:

- **SDCC/RCF** runs everything that reads picoDsts as batch jobs: PicoBinner, and the
  event-level jobs for centrality, the cut flow and the ZDC.
- **A laptop** runs stages 2–4 (ZFitter, RawSpectraModifier, SpectraFitter) and the
  efficiency fits, on the merged PicoBinner output.

The internal PID-contamination modelling in the ported classes is left intact; it tracks all
9 particle species, not just π/K/p. Only the *driver* macros, the `Run*.C` entry points you
actually invoke, are restricted to π/K/p and to this one energy and geometry.

Physics and fit logic in every ported class (`PicoBinner`, `ZFitter`, `RawSpectraModifier`,
`SpectraFitter`, `EfficiencyFitter`, plus their utility dependencies) is byte-verbatim from
the original repo wherever the original code has no energy- or species-specific branching to
remove: it was copied with `cp` and checked with `md5sum`, not retyped. Where a run macro
*did* need to change, the change is documented in a header comment at the top of that macro,
including which values were kept verbatim, which became function arguments, and why.

## Current configuration (Oct 2026)

| | |
|---|---|
| Data | `production_OO_200GeV_2021`, **P24iy**, library **SL24y** (picoDsts) |
| Trigger | **860003** (minbias-bbc_etof), the only O+O trigger that reads out eTOF |
| Runs | 56: the runs containing 860003, minus the four on the official O+O bad-run list and 22134014 (no P24iy picoDsts); `runlists/OO200_860003_runs.txt` |
| Vertex | \|Vz\| ≤ 30 cm (2 cm until 2026-10-26), Vr ≤ 1 cm |
| Pile-up | refMult–nBTOFMatch pol4 band from the refMult-correction pass 1, plus nBTOFMatch ≤ refMult + 100 and refMult ≤ 260 |
| Centrality | corrected **EPD totnMIP**, six classes 0–5, 5–10, 10–20, 20–40, 40–80, 80–100% with edges 364/310/234/130/36/0 ("Centrality" below) |
| Track cuts | nHitsFit ≥ 10, nHitsDedx ≥ 10, nHitsFit/nHitsPoss ≥ 0.52, DCA ≤ 1 cm (`macros/SetCutClass.C`) |

All of this is in `macros/SetCutClass.C`, the one place the event and track cuts are set.

The **October 2026 PicoBinner production** with this configuration has 225,273,772 events
after all event cuts, the same in all three species files. Its merged yield files on RCF are
`/star/data03/pwg/liggett/PicoBinner_OO200/yieldHistos_OO200_{pion,kaon,proton}_Oct2026.root`.
In the |Vz| < 2 cm window it has 30.53M events, about 0.1% fewer than the July cut flow over
the same 56 runs (30.57M), which is what the pile-up band added in October removes. (The July
yield files hold 30.86M because they also include the four runs later dropped as bad.)

The July 2026 yield files, and anything made from them, used |Vz| ≤ 2 cm and raw-refMult
classes {44,37,28,17,5,0}. They are not compatible with the October files class by class.

## Centrality

Centrality is assigned once per event, inside PicoBinner, which calls
`CutClass::centralityIndex()` and bakes the class index into the raw yield histograms.
**Any change to the class scheme therefore needs PicoBinner re-run for all three species.**
`RunZFitter.C`'s `nCentToRun`, `RunRawSpectraModifier.C`'s `nCentBins`, `RunSpectraFitter.C`'s
`a_numCent` and `RunEfficiencyFitter.C`'s `a_numCentBins` all have to match the number of
classes (6). They are kept in sync by hand, not read from `CutClass`.

### EPD totnMIP classes (current)

`SetCutClass.C` calls `cuts->setUseEpdCentrality(true)` with edges `{364,310,234,130,36,0}` for
0–5, 5–10, 10–20, 20–40, 40–80 and 80–100%. An event is in class i if its corrected totnMIP
≥ edge_i (and below the edge of the class above).

Why the EPD and not refMult:

- **Resolution.** refMult cannot resolve centrality beyond about 75% in O+O. The 80% edge falls
  at refMult = 0, so 80–100% cannot be a refMult class.
- **Autocorrelation.** The EPD (2.1 < |η| < 5.1) is away from the mid-rapidity tracks the
  spectra are measured with.

**totnMIP definition.** totnMIP is the Run 21 O+O quantity of StRefMultCorr("totnmip"): the sum
over all EPD hits of nMIP, where tiles with nMIP < 0.3 are skipped and tiles with nMIP > 6
count as 6.

**Trigger.** The official O+O totnMIP calibration is for triggers 860001/2/11/12, so for 860003
the whole chain is redone here. The official triggers are used only for comparison.

| Step | What | Where |
|---|---|---|
| 1 | Inputs: totnMIP vs Vz, ZDCx, nBTOFMatch, refMult, run; 860001/2 comparison set | `macros/MakeEpdCentInputs_OO200.C`, `xml/runEpdCentInputs_OO200_SDCC.xml` → `hadd -f epd_OO200.root …/EpdCent_OO200/*_epd.root` |
| 2 | Vz and luminosity corrections, pile-up band, corrected totnMIP, comparison with the official table, run QA | `macros/FitEpdCent_OO200.C` → `OO200_EpdCent_params.txt`, `hTotnMIPCorr_OO200.root`, `OO200_EPD_*.png` |
| 3 | The corrections as used in PicoBinner | `headers/EpdCentrality_OO200.h` (c_vz per 1 cm, largest 0.35%; c_lumi = p0/(p0 + p1·ZDCx), −0.33% per kHz) |
| 4 | NBD × Glauber scan on corrected totnMIP | `centrality/star_glauber/submit_doScan_OO200.pl -totnmip` (fit totnMIP ≥ 200), then `-totnmip -fine`; see that directory's `README_OO200.md` |
| 5 | 0–100% cuts and the 860003 efficiency vs totnMIP from the best point | `centrality/star_glauber/CentCuts_OO200.C` → `OO200_NBD_totnMIP_cuts.txt` |
| 6 | ⟨Npart⟩, ⟨Ncoll⟩, ⟨Nspec⟩ per class, nominal and Glauber variations | `centrality/star_glauber/NpartPerClass_OO200.C` |
| 7 | Cross-checks | `macros/PlotTrigRatio_OO200.C` (860003 / (860001\|860002) vs totnMIP), `macros/CompareVzCent_OO200.C` (class fractions at \|Vz\| ≤ 30 vs ≤ 2 cm) |

**Best fit.** The fine scan gives npp = 24, k = 6, x = 0.20, d = 0, with χ²/ndf = 0.94 for
totnMIP ≥ 200.

**Fit-range check.** Refits to totnMIP ≥ 150 and ≥ 250 give edges 365/311/235/130/36 and
364/310/233/130/36. The per-class 860003 efficiencies agree within 0.02.

**Comparison with the official table.** The official StRefMultCorr totnMIP edges agree with
this fit's −5% σ_tot edges within 3 units for 0–75%. The edges in use are the nominal ones.

**Open items, before approval by the centrality group:**

- **Peripheral efficiency.** Below 50% the 860003 efficiency from the fit is 4–21% below the
  independent estimate from `PlotTrigRatio_OO200.C`.
- **Refit with the official normalisation.** Still to be done.

### History: raw refMult classes (until Oct 2026)

Until the October production the six classes were cut on raw refMult with edges
`{44,37,28,17,5,0}`, the cuts of an earlier, separately validated analysis. They are kept,
commented out, in `SetCutClass.C`. To go back, use `setUseEpdCentrality(false)` with those edges.

The refMult route was also taken through STAR's full StRefMultCorr procedure for 860003 at
\|Vz\| ≤ 30 cm, and those inputs are still in the repo:

- **Inputs.** `macros/MakeRefMultCorrInputs_OO200.C` with `xml/runRefMultCorrInputs_OO200_SDCC.xml`.
- **Vz, luminosity and pile-up corrections.** `macros/FitRefMultCorr_OO200.C`. Its pass-1
  pile-up band is the one now applied to every event.
- **NBD rescan on corrected refMult.** `submit_doScan_OO200.pl -wide`.
- **Trigger efficiency.** `centrality/star_glauber/FitTriggerEfficiency_OO200.C`.

`StCentrality::Init_OO200GeV()` holds the raw-refMult NBD fit and edges (commit 7ab3afa).

**Stale builds.** One failure from the refMult period still applies.
`makeLibs_SL24y.C` compiles `SetCutClass.C` through `gSystem->CompileMacro()`, which only
recompiles when the source is newer than `bin/SetCutClass_C.so`. If a new `SetCutClass.C`
reaches SDCC with an older mtime (e.g. `rsync`/`scp -p`), SDCC silently keeps the old cuts.
After any cut change, rebuild with `./scripts/rebuild_libs.sh`, not by running
`makeLibs_SL24y.C("clean")` and then the build by hand: "clean" also deletes
`libStPicoDst.so`, and the script restores it.

## PicoDst reader version

`submodule/PicoDstReader_SL24y/` vendors STAR's PicoDst class definitions (StPicoEvent,
StPicoTrack, StPicoArrays, etc. -- the classes that define the actual branch layout
written into `.picoDst.root` files) at STAR library version **SL24y** (as of
2026-07-02, replacing an earlier SL23a-based submodule, to match the real production
dataset's STAR library), re-vendored directly from STAR's own public source
(`https://github.com/star-bnl/star-sw`, tag `SL24y_3`, path `StRoot/StPicoEvent/`) into
this repo's existing standalone-package layout (Grigory Nigmatkulov's `StPicoEvent`
packaging conventions -- Makefile, `picodst_env.sh`, `macros/`, README all unchanged).

Diffing STAR's own SL23a vs SL24y source directly (not just trusting this repo's prior
SL23a-vendored copy) found that almost every class is byte-for-byte identical --
`StPicoTrack`, `StPicoBTofPidTraits`, `StPicoBEmcPidTraits`, `StPicoTrackCovMatrix`, and
every hit/trigger/PID class carried over unchanged. Three things actually differ:

- **`StPicoArrays`**: 4 new detector array types added (`StPicoFwdTrack`,
  `StPicoFwdVertex`, `StPicoFcsHit`, `StPicoFcsCluster` -- STAR's newer Forward/FCS
  detector upgrade), array count 20->24. These 4 new classes are vendored here too, for
  schema completeness, even though this analysis never reads Forward/FCS branches.
- **`StPicoMcTrack`**: 8 new MC-hit-count accessor methods added (unused by this
  analysis).
- **`StPicoEvent`**'s eTOF good-event-flag API: this one is NOT purely additive. The
  flag array grew from 108 entries (per-counter) to 1728 (per-"Get4" subcounter, 16x
  finer), and `eTofGoodEventFlag()` gained a required 4th argument (`iGet4`). This
  repo's only call site (`source/PicoBinner.cxx`, in the eTOF good-event-flag fraction
  loop) has been updated: a counter is now treated as "good" if ANY of its 16 Get4
  sub-flags are good, per explicit instruction -- the closest match to the old, coarser
  per-counter semantics. See that loop's inline comment and `PicoBinner.cxx`'s
  top-of-file header comment (item 9) for the full writeup.

All of `StPicoHelix`/`StPicoPhysicalHelix`/`PhysicalConstants.h`/`SystemOfUnits.h`/
`StPicoMessMgr.h` (generic StarClassLibrary math/utility headers, not PicoDst-specific)
were confirmed byte-identical between the two STAR tags and needed no changes.

**One gap had to be patched, not just ported.** STAR's own real SL24y_3 `StPicoEvent.h`
declares `eTofPulserStatus(iSector,iModule,iCounter)` and
`setETofHasPulsersFlag(vector<bool>)`, but STAR's own `StPicoEvent.cxx` never defines
either of them (confirmed by direct fetch -- a genuine upstream gap, not something lost
in this port). That's harmless for analysis code that never calls them, but this
repo's build (`rootcint` + `StPicoDstLinkDef.h`, which reflects the *entire* class for
ROOT's CINT dictionary) needs every declared method to have a real body just to link
`libStPicoDst.so` -- this only surfaced when building on an actual SDCC node (a Mac
build apparently never linked far enough to hit it). Minimal implementations were
added in `StPicoEvent.cxx`, clearly marked as not-from-upstream, following the same
index-formula convention as the neighboring `eTofGoodEventFlag()`/
`setETofGoodEventFlag()` pair. `PicoBinner.cxx` doesn't call either one, so their exact
behavior doesn't affect any physics result here -- they exist purely so the dictionary
links.

**SL24y is the only reader (since 2026-10-28).** The SL23c reader
(`submodule/PicoDstReader_SL23c`), `macros/makeLibs_SL23c.C` and the `_PICO_READER_SL23c_`
branches were removed in commit 89c7bd3. They had been added in September only to process
the P23ic embedding MuDsts under `starver SL23c`. `submodule/MuDstProcessEmbedding` now
builds against SL24y for a P24iy embedding request (its RequestID is still a placeholder,
`SET_P24IY_REQUEST_ID`). To reprocess the existing P23ic embedding, check out 89c7bd3^.

## One-time build setup

1. Build the PicoDst reader submodule and copy its library into `bin/`:
   ```
   cd submodule/PicoDstReader_SL24y
   make
   cp libStPicoDst.so ../../bin/
   cp StPicoDst_Dict_rdict.pcm ../../bin/ 2>/dev/null || true   # see note below
   cd ../..
   ```
   **On `StPicoDst_Dict_rdict.pcm`**: a local macOS build produces this dictionary
   sidecar alongside `libStPicoDst.so`, but SDCC's ROOT/rootcling does not generate one
   for this submodule at all (confirmed via an actual SDCC build directory listing --
   only `.o` files plus `StPicoDst_Dict.C/.h`, no `.pcm`), and it isn't needed there
   either (`xml/runPicoBinner_OO200_SDCC_template.xml` packages only `libStPicoDst.so`,
   matching how lightflavorspectra_etof's own working production job does it). Copy it
   opportunistically if it exists; don't treat its absence as an error, especially on
   SDCC. `scripts/rebuild_libs.sh` handles this automatically.
2. Compile everything else via ACLiC (run from the repo root):
   ```
   root -l -q -b macros/makeLibs_SL24y.C
   ```
   (Renamed 2026-08-24 from `macros/makeLibs.C` -- this is the SL24y-STAR_LEVEL build
   macro, used both here on a laptop and on SDCC for stage-1 PicoBinner batch jobs.)
   This compiles, in dependency order: `namespaces.cxx` (Helix/PhysMath/HistogramUtilities/
   PhysMathHistogramUtilities combined, matching the original's own build), `ParticleInfo.cxx`,
   `CutClass.cxx`, `CalibrationClass.cxx`, `PicoBinner.cxx`, `SetCutClass.C`, `GausMixModel.cxx`
   + `Histo2D.cxx` + `ZFitter.cxx`, `SpectraStructure.cxx` + `SplineFitter.cxx` +
   `RawSpectraModifier.cxx`, `FitUtilities.cxx` + `SpectraFitter.cxx`, and the efficiency
   tools: `GraphFitter.cxx` + `EfficiencyFitter.cxx` + `RunEfficiencyFitter.C`, then
   `MattMcEvent.cxx` + `MattMcTrack.cxx` + `EmbeddingBinner.cxx` + `RunEmbeddingBinner.C` and
   `HybridFeeddownMaker.cxx`. Expect this to take several minutes the first time (ZFitter.cxx
   and SpectraFitter.cxx are each tens of thousands of lines once their `#include` chains are
   counted).

   **ROOT 5 at RCF.** `starver SL24y` gives ROOT 5.34/38 and gcc 4.8.5, so anything compiled
   or run there must be C++98 and CINT-safe: no `auto`, `nullptr`, lambdas or range-for in
   anything rootcint sees. Several macros guard their bodies with
   `#if !defined(__CINT__) || !defined(__MAKECINT__)` for this reason. The laptop has ROOT 6,
   so code that builds locally can still fail at RCF.

   For a clean rebuild (e.g. after editing a header any class depends on), run
   `./scripts/rebuild_libs.sh` from the repo root (or `bash scripts/rebuild_libs.sh`). It rebuilds
   the PicoDst reader submodule, wipes `bin/` with `makeLibs_SL24y.C("clean")`, copies
   `libStPicoDst.so` back into `bin/`, and runs `makeLibs_SL24y.C`. Do not run the "clean" macro
   on its own followed by the build: "clean" also deletes `libStPicoDst.so`, and the build then
   fails to link PicoBinner (undefined references to `StPicoDst::picoArrays` etc.).

## Pipeline stages

Each stage is a `Load<Stage>Libs.C` + `Run<Stage>.C` pair. **Always load the libs macro
first, then the run macro, as two arguments to the same `root` invocation** -- this
two-step split avoids a cling JIT/link-resolution issue that shows up if everything is
loaded and executed inside one macro (see the comments in `RunPicoBinner.C` for the full
explanation). Run everything from inside `macros/`.

### 1. PicoBinner -- raw yield histograms from PicoDst files

Run once per particle species (three invocations), each producing its own output file:

```
root -l -q -b LoadPicoBinnerLibs.C 'RunPicoBinner.C("pionFileList.list","pionYield.root",9,"PION","true","false","false","false")'
root -l -q -b LoadPicoBinnerLibs.C 'RunPicoBinner.C("kaonFileList.list","kaonYield.root",9,"KAON","true","false","false","false")'
root -l -q -b LoadPicoBinnerLibs.C 'RunPicoBinner.C("protonFileList.list","protonYield.root",9,"PROTON","true","false","false","false")'
```

(`numParticles` should stay `9` regardless of which single species you're binning --
see the header comment in `source/PicoBinner.cxx`.) These three output files
(`pionYield.root`/`kaonYield.root`/`protonYield.root`) are the inputs every later stage
expects.

For a full production dataset, use `scripts/run_picobinner_resilient.sh` instead of
calling `RunPicoBinner.C` directly -- it chunks the input filelist and retries failed
chunks individually, working around a non-deterministic ROOT crash that shows up when
reading very large filelists in one process (see that script's own header comment for
the full story).

The production yield files come from SDCC instead: per-species PicoBinner jobs, a parallel
chunk merge and a final `hadd` ("SDCC batch submission" below). The October 2026 production
was 767 jobs per species. Copy the three merged files to the laptop with `sftp` (`get`) and
continue with stage 2.

### 2. ZFitter -- dE/dx (TPC) and 1/beta (BTOF) PID fits

Run once per (particle, detector) combination -- six invocations total:

```
root -l -q -b LoadZFitterLibs.C 'RunZFitter.C(0,0,"pionYield.root","kaonYield.root","protonYield.root","/path/to/outputDir/","etaBetaGammaDeDx.root")'   # pion,  TPC
root -l -q -b LoadZFitterLibs.C 'RunZFitter.C(0,1,"pionYield.root","kaonYield.root","protonYield.root","/path/to/outputDir/","etaBetaGammaDeDx.root")'   # pion,  BTOF
root -l -q -b LoadZFitterLibs.C 'RunZFitter.C(1,0,"pionYield.root","kaonYield.root","protonYield.root","/path/to/outputDir/","etaBetaGammaDeDx.root")'   # kaon,  TPC
root -l -q -b LoadZFitterLibs.C 'RunZFitter.C(1,1,"pionYield.root","kaonYield.root","protonYield.root","/path/to/outputDir/","etaBetaGammaDeDx.root")'   # kaon,  BTOF
root -l -q -b LoadZFitterLibs.C 'RunZFitter.C(2,0,"pionYield.root","kaonYield.root","protonYield.root","/path/to/outputDir/","etaBetaGammaDeDx.root")'   # proton,TPC
root -l -q -b LoadZFitterLibs.C 'RunZFitter.C(2,1,"pionYield.root","kaonYield.root","protonYield.root","/path/to/outputDir/","etaBetaGammaDeDx.root")'   # proton,BTOF
```

`a_partIndex`: 0=pion, 1=kaon, 2=proton. `a_detectorIndex`: 0=TPC, 1=BTOF (ETOF fit
drivers were not ported -- see `headers/ZFitter.h`). Each of the six invocations calls
`fitter->makeSpectra(...)` and writes its own spectra file into `a_outputDir`, named
`spectra_OO200_<Particle>_<Detector>_RapBin_-1_-1_-1_-1_-1_-1.root`. Keep the trailing `/`
on `a_outputDir`: the file names are appended to it directly.

**Arguments used for the production fits.** The fits pass four more arguments after the
η–βγ file: `"",false,false,false`, meaning no quick-start file, no quick load, no quick-start
creation, and refit the η–βγ dE/dx curve (written to the η–βγ file path given) instead of
loading it. Each run takes hours. Log every run to its own file.

**Merging TPC and BTOF.** `RawSpectraModifier` needs TPC and BTOF in one file per species.
Merge with `macros/MergeDetectorSpectraFiles.C`, **not** `hadd`:

```
for s in Pion Kaon Proton; do
  root -l -q -b "MergeDetectorSpectraFiles.C(\"out/spectra_OO200_${s}_TPC_RapBin_-1_-1_-1_-1_-1_-1.root\",\"out/spectra_OO200_${s}_BTOF_RapBin_-1_-1_-1_-1_-1_-1.root\",\"\",\"out/spectra_OO200_${s}_Merged.root\")"
done
```

The third argument is the eTOF spectra file. It stays `""` until eTOF is calibrated.

A bare `hadd` gets two things wrong:

- **Histo2D content is lost.** Without the Histo2D dictionary loaded, `hadd` reduces the
  `DeDx_FitData/.../mTm0_ZTPC_Proj2D` objects to empty `TNamed`.
- **The βγ histograms double.** Every run fills `betagamma_FitData/`, so `hadd` counts it twice.

`MergeDetectorSpectraFiles.C` avoids both: it copies only the directories `RawSpectraModifier`
reads, and loads the dictionary first. `macros/MergeZFitterSpectra.C` (a `hadd` with the
dictionary loaded) fixes the first problem only.

### 3. RawSpectraModifier -- efficiency/energy-loss/feed-down corrections

This stage reads the yield files, a merged spectra file from stage 2 and one efficiency
file. Any of the three yield-file arguments can be `""` to skip that species. The July
outputs were made one species at a time: `spectra_OO200_<S>_Merged.root` in,
`spectra_OO200_<S>_Modified.root` out.

```
root -l -q -b LoadRawSpectraModifierLibs.C 'RunRawSpectraModifier.C("pionYield.root","","","spectra_OO200_Pion_Merged.root","efficiency_combined.root","spectra_OO200_Pion_Modified.root")'
```

**The efficiency file** is a `hadd` of two independent pieces. Their object paths do not
overlap, so a plain `hadd` is safe here:

1. **The BTOF matching efficiency, from data.** `macros/RunEfficiencyFitter.C` reads the
   same PicoBinner yield files (`nTPCTracks*_ForBTOF` and `nBTOFTracks*`) and fits
   BTOF/TPC vs mT−m0. It writes `BTOFEfficiencyGraphs/` and `BTOFEfficiencyFits/`. It does not
   depend on ZFitter, so it can run before, after or alongside stage 2. Pass the proton yield
   file as the 4th argument; the default `""` skips protons, a leftover from a corrupt July
   proton file.
   ```
   root -l -q -b LoadEfficiencyFitterLibs.C 'RunEfficiencyFitter.C("pionYield.root","kaonYield.root","efficiency_BTOF.root","protonYield.root","./EfficiencyFitterImages/")'
   ```
   By default it makes one all-centrality fit per species, charge and rapidity bin, and copies
   it into each of the six classes (`a_fitAllCentBTOF`, `a_useAllCentBTOFFits`). Rapidity bins
   outside the BTOF acceptance have no BTOF tracks. They print
   `truncateHistogram … wouldn't have enough bins` followed by ROOT's
   `TAxis::Set: bins must be in increasing order`, which is expected. Check the log for
   `BAD FIT` instead.
2. **TPC efficiency and energy loss, from embedding.**
   - `submodule/MuDstProcessEmbedding` (RCF) turns the embedding MuDsts into trees.
   - `macros/RunEmbeddingBinner.C` (laptop) bins them.
   - `macros/RunEfficiencyFitter_Embedding.C` fits them and writes `efficiency_OO200_Embedding.root`.

   The only O+O 200 GeV embedding that exists is **P23ic** (SL23c). The current
   `efficiency_OO200_Embedding.root` (Sep 2026) was made from it. A P24iy request, generated
   flat in y to cover eTOF, is still to be submitted.

```
hadd efficiency_combined.root efficiency_BTOF.root efficiency_OO200_Embedding.root
```

`macros/CheckBTOFCorrection.C` compares raw and BTOF-corrected graphs in a `*_Modified.root`
file. It reports dropped points and the corrected/raw ratio for each graph, which confirms
that the correction was applied.

See the extensive header comment in `RunRawSpectraModifier.C` for which of its optional
arguments (momentum ranges, rapidity cuts, hybrid feed-down toggle, trim data) were fixed
verbatim from the original vs. exposed as arguments because they genuinely vary by
dataset and have no valid default for a brand-new O+O 200 GeV analysis.

### 4. SpectraFitter -- final invariant-yield spectra fits + dN/dy

One invocation, consuming stage 3's output:

```
root -l -q -b LoadSpectraFitterLibs.C 'RunSpectraFitter.C("spectra_modified.root","fit_output.root")'
```

Fits every (particle, centrality, charge, rapidity bin) twice -- once "Nominal"
(Tsallis for pion/kaon, BlastWave for proton) and once "Systematic" (BoseEinstein for
pion, MtM0Expo for kaon, Thermal for proton) -- matching the only fitting call sequence
confirmed live in the original Collider-mode macro (see `RunSpectraFitter.C`'s header
comment for how that was verified: a large fraction of the original macro turned out to
be dead, block-commented-out alternative fit sequences, confirmed both by a
character-level comment/string-state script and by direct reading).

### Convenience script for stages 2-4

**Out of date: do not use it as is.** The script still merges the six ZFitter files with a
bare `hadd` (see stage 2 for why that is wrong) and runs RawSpectraModifier once on all three
species. It also uses ZFitter's default arguments, not the production ones. Run stages 2–3 by
hand as described above until it is updated.

`scripts/run_stages_2_to_4.sh` wraps the six ZFitter invocations, the `hadd` merge step,
`RawSpectraModifier`, and `SpectraFitter` into one command, using default arguments for
everything not explicitly passed in:

```
./scripts/run_stages_2_to_4.sh pionYield.root kaonYield.root protonYield.root \
    etaBetaGammaDeDxFile.root efficiencyFile.root /path/to/workDir
```

Edit the script directly if you need to override any of `RawSpectraModifier`'s or
`SpectraFitter`'s non-default arguments (momentum ranges, rapidity cuts, fit seeds, etc.)
for your specific dataset -- see each stage's own `Run*.C` header comment for what those
defaults are and which ones are worth sanity-checking rather than trusting blindly.

### 5. Drawing spectra and dN/dy plots (optional)

`macros/DrawSpectraAndDnDy.C` reads stage 4's output (`fit_output.root`) and saves two
families of `.png` plots: invariant spectra (`TotalSpectra/...`, one .png per particle/
charge/centrality bin, overlaying every rapidity bin found, decade-offset for
readability) and dN/dy vs rapidity (`ParticleYields/.../dNdy_..._Nominal` and
`..._Systematic`, one .png per particle/charge/fit-variant, overlaying every centrality
bin found). It doesn't assume which bins actually got filled -- it probes every slot up
to the fixed array bounds in `headers/SpectraFitter.h` and just skips anything that
wasn't written. See the file's own header comment for details (including why the
Fiducial/extrapolation-piece dN/dy histograms SpectraFitter also writes are deliberately
not drawn here).

```
root -l -q -b LoadDrawingLibs.C 'DrawSpectraAndDnDy.C("fit_output.root")'
```

Optional arguments: output directory (default `./spectra_plots`), and booleans to draw
only spectra, only dN/dy, or turn off the decade-offsetting -- see the macro's function
signature.

## What's deliberately NOT in this repo

The 6-other-energy configuration data/paths, FixedTarget-mode code paths, and ZFitter's
`create2DProjectionFile`/deuteron-triton/electron-study branches (all dead or irrelevant
to this configuration) were dropped rather than carried forward as inert code. Each
`Run*.C` macro's header comment documents exactly what was dropped from its
corresponding original macro and why it was safe to drop.

Also not ported:

- **ZFitter's eTOF fit drivers.** `a_detectorIndex` is 0 (TPC) or 1 (BTOF) only.
- **EfficiencyFitter's feed-down and knockout fits.** Only the BTOF (data) and the
  embedding TPC/energy-loss fits were ported.
- **SDCC batch submission for stages 2–4.** Those stages run locally on PicoBinner's
  merged output.

The eTOF calibration itself (position, run-by-run Y, T0, and GET4 clock-jump states for
2021) is worked on in the separate `eTOFCalibrations` repo.

## SDCC batch jobs (real production data)

The build and pipeline stages 1–4 run the same way on a laptop and on an SDCC/RCF node:
same source, same macros, same physics. What runs on SDCC are the jobs that read the full
picoDst dataset, which is too large for one local `root` process. SUMS fans the file list
out across many independent worker jobs. Stages 2–4 then run on the laptop, on the merged
PicoBinner output.

### What's here

All jobs use the same catalog query (P24iy, SL24y, the 56-run whitelist) and the SL7
container `<shell>` line. Each XML's header comment lists its output directory, the `mkdir`
it needs and its merge command.

| Job (XML) | Runs | Output |
|---|---|---|
| `runPicoBinner_OO200_SDCC_template.xml` | PicoBinner, one species per job, via `scripts/run_picobinner_sdcc.bash` | `*_Processed.root` yield fragments |
| `haddPicoBinner_SDCC_template.xml` | stage 1 of the merge: `hadd` of 15 fragments per process, via `scripts/hadd_stage1_chunk_worker.sh` | `chunk_*.root` |
| (by hand) `scripts/hadd_final_merge.sh <chunkDir> <out.root> "chunk_*.root"` | stage 2 of the merge, inside the container | `yieldHistos_OO200_<sp>_<tag>.root` |
| `runCutFlow_OO200_SDCC.xml` | `MakeEventCutFlow_OO200.C` via `scripts/run_cutflow_sdcc.bash` (Event branch only) | per-run cut flow, Vz, luminosity |
| `runRefMultCorrInputs_OO200_SDCC.xml` | `MakeRefMultCorrInputs_OO200.C` | refMult-correction inputs |
| `runEpdCentInputs_OO200_SDCC.xml` | `MakeEpdCentInputs_OO200.C` | EPD totnMIP centrality inputs |
| `runZdcNeutrons_OO200_SDCC.xml` | `MakeZdcSpectatorNeutrons_OO200.C` | ZDC ADC by EPD class (`ZdcNeutrons_OO200_Oct2026/`) |
| `runPicoBinner_OO200_VzStudy.xml` | PicoBinner on 1500 files, for the Vz-window study | input to `macros/VzStudy_OO200.C` |

`scripts/cleanup_sdcc_scratch.sh` clears the scratch-area logs and scheduler files that pile
up from repeated submissions (see "Cleaning up scratch/logs" below). It does not touch
permanent output.

The XML template was restructured (2026-07) around the format of
`lightflavorspectra_etof`'s own real, working `pions_OO_200GeV` job -- including a
`<shell>` line that runs the job inside SDCC's required singularity container
(`/cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif`) and a `filesPerHour` throttle -- with
`library=`/`production=` kept at this repo's `SL24y`/`P24iy` port rather than reverting to
that job's `SL23c`/`P23ic`. Paths are filled in for `liggett2820`'s real SDCC account
(`/star/data03/scratch/liggett/...` for scratch logs, `/star/data03/pwg/liggett/...` for
permanent output) rather than left as `EDITME` placeholders. If someone else runs this
template, every `/liggett` path needs updating to their own account first.

### Before you can submit

1. **Push local commits, then `git pull` on SDCC, before rebuilding.** Any source edit
   made on a laptop (or anywhere else) only exists on SDCC once it's been pushed to
   GitHub *and* pulled into the SDCC checkout -- `bin/*.so` there is built from whatever
   `git pull` last brought in, not from whatever's newest on your laptop. It's easy to
   edit `PicoBinner.cxx`/`CutClass.cxx`/`SetCutClass.C` etc., rebuild and test locally,
   and then submit an SDCC job assuming it picked up the change, when SDCC is actually
   still running an older checkout -- the job won't error, it'll just silently produce
   output missing whatever the new code was supposed to add (this happened in practice:
   a new `eta_nHitsFit` histogram was added, committed, and even rebuilt locally, but the
   commit sat unpushed while a real production PicoBinner job ran on SDCC against the old
   code, so the resulting yield file simply didn't have the histogram -- nothing failed,
   it just wasn't there). Checklist, every time before submitting:
   ```
   # wherever you made the edit:
   git push origin master
   # on the SDCC/RCF interactive node:
   git pull
   ```
2. **Rebuild on SDCC, not just on your laptop.** `bin/*.so` are not committed to this
   repo (see `.gitignore`) because a laptop build produces macOS/arm64 binaries that
   will not load on BNL's Linux nodes. After `git pull` on an SDCC/RCF interactive node,
   run the "One-time build setup" steps above (build `submodule/PicoDstReader_SL24y`,
   then `root -l -q -b macros/makeLibs_SL24y.C`) there -- the `<SandBox>` packages whatever's
   actually in `bin/` at submission time, so it needs to already be genuine Linux
   binaries built from the code you just pulled. If in doubt whether a change actually
   recompiled, force it with `./scripts/rebuild_libs.sh` (a full clean rebuild, including
   `libStPicoDst.so`; see
   "Centrality" above for a case where a stale `.so` silently kept running old
   cut values even though the source was correct).
3. **If you're not liggett2820, update every `/liggett` path in the XML** -- the scratch
   log directory (`<stdout>`/`<stderr>`), the `<Generator><Location>`, and the permanent
   output directory (the `mv` line in `<command>` plus the `<output>` element). These are
   filled in for one specific SDCC account and disk allocation; they won't resolve for
   anyone else as-is.
4. **The catalog query.** `trgsetupname=production_OO_200GeV_2021,production=P24iy,library=SL24y`
   with the 56-run whitelist resolves to files: the October 2026 production ran on it. The
   whitelist is the same in every XML and in `runlists/OO200_860003_runs.txt`. It is the runs
   that contain trigger 860003 (22133032–22136010), minus:
   - 22134014, which has no P24iy picoDsts;
   - the four runs on the official O+O bad-run list in StRefMultCorr/BadRun.h: 22133032,
     22133037, 22133038 and 22135039.
5. **One copy of the template per species** (PION/KAON/PROTON). In each, edit the `PARTNAME`
   argument in the `run_picobinner_sdcc.bash` call, the job `name`, and the output paths, so
   the three jobs don't collide. The committed template writes to `SL24y_PION/`; the October
   2026 production used `PicoBinner_OO200/Oct2026_<SP>/` and `Oct2026_<SP>_chunks/`.
6. **Delete the sandbox zip after any code change** (`rm -f <PackageName>.zip` in the
   directory you submit from). SUMS reuses an existing zip without repacking it, so jobs would
   otherwise run the old `bin/` and `macros/`.

### Submitting and merging

Host vs container: `star-submit`, `condor_q` and `condor_rm` run **on the host**. `root`,
`hadd` and anything after `starver SL24y` run **inside the SL7 container**:

```
singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif /bin/tcsh
```

1. Submit each species from the repo root on the host:
   `star-submit xml/<your PION copy>.xml`, and the same for KAON and PROTON.
2. When `condor_q` is empty, check that every job wrote output and that none of the outputs
   are tiny (see "Practical notes at RCF" below).
3. Merge stage 1: list each species' `*_Processed.root` into a file list and submit
   `haddPicoBinner_SDCC_template.xml`, one copy per species.
4. Merge stage 2, inside the container:
   ```
   starver SL24y
   ./scripts/hadd_final_merge.sh <chunkDir>/ <pwg>/yieldHistos_OO200_<sp>_<tag>.root "chunk_*.root"
   ```
5. Copy the three merged files to the laptop with `sftp` and continue with stage 2 (ZFitter).
   Before fitting, check that the three files have the same number of events (`zVertex`
   entries; it is filled after all event cuts).

### Practical notes at RCF

**Jobs that hang after writing their output.** Macros that load the PicoDst libraries
through `LoadCutFlowLibs.C` can abort or hang during ROOT teardown, because the PicoDst
dictionaries are registered twice. SUMS then never copies the output back. The fix,
already in `MakeEpdCentInputs_OO200.C` and `MakeZdcSpectatorNeutrons_OO200.C`, is to close
the output file, flush stdout and stderr, and call `_exit(0)`. Do the same in any new macro
of this kind. To spot it, grep the job `.out` files for `root exited`: anything other than
`0` is this problem.

**Disk quota.** `/star/data03/pwg/liggett` has a quota. When it is full, jobs and `hadd`
write empty or truncated files without failing loudly. After every production or merge,
look for small files:
```
find <outDir> -name '*.root' -size -10k | wc -l
```
Before a large production, delete old per-job fragments and chunk files once their merged
file has been checked.

**csh.** The RCF login shell is tcsh, so:
- write `${var}:` rather than `$var:`, because `:` after a variable is a csh modifier;
- don't use `2>/dev/null`;
- avoid `!` in heredocs, because it triggers history expansion.

### Cleaning up scratch/logs

Repeated submissions pile up `.out`/`.err` logs under
`/star/data03/scratch/liggett/picoBinner_OO200/` and scheduler auxiliary files under
`/star/data03/scratch/liggett/schedFiles/picoBinner_OO200/`.
`scripts/cleanup_sdcc_scratch.sh` clears both (with a confirmation prompt), without
touching the permanent `pwg` output directory:

```
./scripts/cleanup_sdcc_scratch.sh --dry-run   # see what would be removed, delete nothing
./scripts/cleanup_sdcc_scratch.sh             # prompts before deleting
./scripts/cleanup_sdcc_scratch.sh -y          # skip the confirmation prompt
```

## Other studies in this repo

These macros run alongside the spectra pipeline. Each has a usage line in its header comment.

| Study | Macros | Notes |
|---|---|---|
| Event cut flow | `MakeEventCutFlow_OO200.C`, `PlotEventCutFlow_OO200.C` | Per-run counts (all / 860003 / after cuts), Vz, TPC vs Vz, luminosity. Still uses the July selection (\|Vz\| ≤ 2 cm, linear pile-up lines). |
| Vz window | `VzStudy_OO200.C` (reads PicoBinner's `VzStudy/` directory), `CompareVzCent_OO200.C` | Raw yields, BTOF matching and EPD class fractions per Vz slice |
| TPC–BTOF consistency | `PlotTpcBtofMatch_OO200.C` | BTOF matching efficiency, TPC vs BTOF/eff yields and h−/h+ at mid-rapidity, from a `*_Merged.root` and the yield file |
| ZDC spectator neutrons | `MakeZdcSpectatorNeutrons_OO200.C`, `PlotZdcSpectatorNeutrons_OO200.C`, `FitZdcSingleNeutronPeak_OO200.C`, `CheckZdcNeutronCapture_OO200.C` | Unattenuated ZDC ADC by EPD class, 1n-peak calibration. Since Oct 26 it uses the production selection and EPD classes (40–100% also in 10% sub-classes). The July plots used \|Vz\| ≤ 2 cm and refMult classes. |
| Detector read-out check | `CheckDetectorHits_OO200.C` | Which hit collections are filled for 860003 |
| Event display | `MakeEventDisplay_OO200.C` (`LoadEventDisplayLibs.C`) | One central event, track helices. Still uses refMult ≥ 44 for 0–5%. |
| Dose estimates | `dose/` | A separate side project (energy-dose fractions, Mars GCR and neutron Geant4 sims), not part of the spectra analysis |

Slides and the analysis note draft are in `docs/`.
