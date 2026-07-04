# lightflavorspectra_OO200

A simplified fork of `lightflavorspectra_etof`, restricted to a single physics
configuration: **Au+Au-style O+O collisions at 200 GeV, ColliderCenter geometry,
pion/kaon/proton spectra only**, running entirely on a local laptop (no SDCC/RCF batch
infrastructure).

The internal PID-contamination-modeling machinery in the ported classes (which tracks
all 9 particle species, not just pi/K/p, for cross-contamination corrections) is left
fully intact. Only the *driver* macros -- the `Run*.C` entry points you actually invoke
-- are restricted to pi/K/p and to this one energy/geometry.

Physics and fit logic in every ported class (`PicoBinner`, `ZFitter`, `RawSpectraModifier`,
`SpectraFitter`, plus their utility dependencies) is byte-verbatim from the original repo
(copied with `cp` and verified with `md5sum`, not retyped) wherever the original code has
no energy- or species-specific branching to remove. Where a run macro *did* need to
change (dropping the 7-energy scan, dropping SDCC paths, restricting to pi/K/p), every
change is documented in a header comment at the top of that macro, including exactly
which values were fixed verbatim vs. exposed as new function arguments, and why.

## Centrality binning

`macros/SetCutClass.C` configures a **5-bin** raw refMult centrality scheme via
`CutClass::setCentralities()`: 0-5%, 5-10%, 10-20%, 20-40%, 40-80% (as of 2026-07-02,
replacing an earlier 9-bin scheme -- see that call's header comment for the previous
9-bin cut values, kept commented out there for reference). The refMult cut values for
this 5-bin scheme are not newly derived -- they're the same cuts the 9-bin scheme
already used at the percentile boundaries that survive merging (5/10/20/40/80% were
already 9-bin edges; only the intermediate 30/50/60/70% edges are dropped).

Centrality classes are defined by integrating the Glauber-model multiplicity
distribution from most-central (highest multiplicity) to most-peripheral, so the
refMult cut must be *highest* for the 0-5% bin and decrease monotonically to the most
peripheral bin. The as-ported 9-bin cut array had these values ascending in parallel
with ascending percents (same magnitudes, wrong order) -- `CutClass::centralityIndex()`
only discriminates central from peripheral correctly with descending edges, so this was
corrected (edge order reversed, values not re-derived) as part of the 2026-07-02 change.

Centrality assignment happens once, per event, inside `PicoBinner` (it calls
`CutClass::centralityIndex()` and bakes the resulting centIndex into the raw yield
histograms) -- so **any change to this scheme requires re-running PicoBinner for all
three species**; existing yield files binned under a different scheme are not
compatible and must be regenerated. `RunZFitter.C`'s `nCentToRun`, `RunRawSpectraModifier.C`'s
`nCentBins`, and `RunSpectraFitter.C`'s `a_numCent` default all need to match whatever
scheme is actually configured here (currently 5) -- they're kept in sync manually, not
read from `CutClass` automatically, so double-check all three if you change this again.

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

## One-time build setup

1. Build the PicoDst reader submodule and copy its library and dictionary pcm into
   `bin/` (2026-07-04: this step used to only copy `libStPicoDst.so`, missing
   `StPicoDst_Dict_rdict.pcm` -- that's a separate file `make` also produces, and
   without it in `bin/`, RunPicoBinner.C's SDCC SandBox packaging silently ends up
   missing it, since the XML references `bin/StPicoDst_Dict_rdict.pcm` directly):
   ```
   cd submodule/PicoDstReader_SL24y
   make
   cp libStPicoDst.so ../../bin/
   cp StPicoDst_Dict_rdict.pcm ../../bin/
   cd ../..
   ```
2. Compile everything else via ACLiC (run from the repo root):
   ```
   root -l -q -b macros/makeLibs.C
   ```
   This compiles, in dependency order: `namespaces.cxx` (Helix/PhysMath/HistogramUtilities/
   PhysMathHistogramUtilities combined, matching the original's own build), `ParticleInfo.cxx`,
   `CutClass.cxx`, `CalibrationClass.cxx`, `PicoBinner.cxx`, `SetCutClass.C`, `GausMixModel.cxx`
   + `Histo2D.cxx` + `ZFitter.cxx`, `SpectraStructure.cxx` + `SplineFitter.cxx` +
   `RawSpectraModifier.cxx`, and `FitUtilities.cxx` + `SpectraFitter.cxx`. Expect this to take
   several minutes the first time (ZFitter.cxx and SpectraFitter.cxx are each tens of
   thousands of lines once their `#include` chains are counted).

   Re-run `root -l -q -b 'macros/makeLibs.C("clean")'` then re-run the build if you ever need
   a clean rebuild (e.g. after editing a header any class depends on).

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
`fitter->makeSpectra(...)` (see `RunZFitter.C` line ~501) and writes its OWN spectra file
into `a_outputDir`, named `spectra_OO200_<Particle>_<Detector>_RapBin_....root` -- there
is no automatic merge step across the six. Combine them into the single file
`RawSpectraModifier` expects with ROOT's `hadd` before moving to stage 3 (this matches
the original repo's own naming convention -- its per-energy spectra files are all named
with `_Combined_` in them, e.g. `spectra_27GeV_2023_08_12_Combined_modified.root` --
and its `SeparateAnalysis/combineAndTailorEfficiencies.bash` script uses the same
`hadd`-based approach for its efficiency files):
```
hadd spectra.root /path/to/outputDir/spectra_OO200_*.root
```

### 3. RawSpectraModifier -- efficiency/energy-loss/feed-down corrections

One invocation, consuming stage 1's three yield files and stage 2's spectra file, plus
an embedding-derived efficiency file you supply for this dataset:

```
root -l -q -b LoadRawSpectraModifierLibs.C 'RunRawSpectraModifier.C("pionYield.root","kaonYield.root","protonYield.root","spectra.root","efficiency.root","spectra_modified.root")'
```

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

Full SDCC/RCF batch submission for stages 2-4 (ZFitter/RawSpectraModifier/SpectraFitter)
is also not here -- only PicoBinner runs on SDCC (see below); the rest of the pipeline is
meant to run locally against PicoBinner's output, exactly as documented in the "Pipeline
stages" section above.

## SDCC batch submission (PicoBinner only, real production data)

Everything above (build, and pipeline stages 1-4) runs identically whether you're on a
laptop or an SDCC/RCF node -- same source, same macros, same physics. The only thing
added here is a way to run stage 1 (PicoBinner) as a proper SDCC batch job against a
full production dataset, since a real production filelist is too large to bin in one
long-lived local `root` process (see `scripts/run_picobinner_resilient.sh`'s own header
comment for the crash this works around locally; SDCC's scheduler sidesteps the same
problem differently, by fanning the filelist out across many independent worker jobs).

**Only PicoBinner is meant to run on SDCC.** The intended workflow is: run PicoBinner on
SDCC against the full production dataset, copy the resulting `*_Processed.root` yield
files back to a laptop, then run ZFitter/RawSpectraModifier/SpectraFitter locally exactly
as described above in "Pipeline stages" (2-4) -- nothing about those three stages changes
based on where PicoBinner ran.

### What's here

- `xml/runPicoBinner_OO200_SDCC_template.xml` -- a `star-submit`-style job description,
  one species per job (PION/KAON/PROTON -- make three copies, one per species, changing
  the `PARTNAME` line and job name/output paths in each).
- `scripts/run_picobinner_sdcc.bash` -- the wrapper the XML's `<command>` block invokes;
  it just resolves paths and calls this repo's own `LoadPicoBinnerLibs.C` +
  `RunPicoBinner.C` pair with the arguments the scheduler passes through, naming its
  output the way the XML's `<output fromScratch="*_Processed.root">` glob expects.
- `scripts/cleanup_sdcc_scratch.sh` -- housekeeping for the scratch-area log/scheduler
  clutter that piles up from repeated submissions (see "Cleaning up scratch/logs"
  below). Does not touch permanent output.

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

1. **Rebuild on SDCC, not just on your laptop.** `bin/*.so` are not committed to this
   repo (see `.gitignore`) because a laptop build produces macOS/arm64 binaries that
   will not load on BNL's Linux nodes. After cloning this repo onto an SDCC/RCF
   interactive node, run the "One-time build setup" steps above (build
   `submodule/PicoDstReader_SL24y`, then `root -l -q -b macros/makeLibs.C`) there first --
   the `<SandBox>` packages whatever's actually in `bin/` at submission time, so it needs
   to already be genuine Linux binaries.
2. **If you're not liggett2820, update every `/liggett` path in the XML** -- the scratch
   log directory (`<stdout>`/`<stderr>`), the `<Generator><Location>`, and the permanent
   output directory (the `mv` line in `<command>` plus the `<output>` element). These are
   filled in for one specific SDCC account and disk allocation; they won't resolve for
   anyone else as-is.
3. **Double-check the catalog query.** The `trgsetupname=production_OO_200GeV_2021` part
   of the `<input URL="catalog:...">` line, and its runnumber whitelist, are copied from
   `lightflavorspectra_etof`'s working O+O 200 GeV XML job; `production=`/`library=` have
   been changed to `P24iy`/`SL24y` per this repo's port (the reference job instead paired
   `production=P23ic` with `library=SL23c`). This `production=`/`library=` pairing -- with
   or without the runnumber whitelist carried over from the SL23c/P23ic job -- has not
   been independently re-verified against the current StarCatalog, so confirm it still
   resolves to files before relying on it (an unverified pairing may just return zero
   files).
4. Make one copy of the template per species (PION/KAON/PROTON), editing the `PARTNAME`
   argument in the `run_picobinner_sdcc.bash` call, the job `name`, and the output paths
   in each copy so the three jobs don't collide.

### Submitting

```
star-submit xml/runPicoBinner_OO200_SDCC_template.xml
```

(or whatever you renamed your per-species copies to). Once all three species' jobs
finish, copy the merged `*_Processed.root` output files back to your laptop and continue
with stage 2 (ZFitter) exactly as documented above.

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
