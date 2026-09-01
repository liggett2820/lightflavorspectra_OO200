#ifndef _MAKEFILE_TOGGLES_
#define _MAKEFILE_TOGGLES_

// =====================================================================================
// makefile_toggles_RCF.h -- lightflavorspectra_OO200/submodule/MuDstProcessEmbedding
//
// RCF-flavored companion to the top-level ../../makefile_toggles.h, for compiling ONLY
// this embedding submodule (CutClass, HistogramUtilities, PhysMath, ParticleInfo, Helix,
// MattMcEvent/Track, StMuAnalysisMaker) under `starver` at RCF/SDCC. compile.bash's own
// header comment explains why the top-level file can't be reused unmodified here: it
// hardcodes _MAC_OSX_ (= "not RCF") for the local/laptop build, which selects the WRONG
// branches in CutClass.h/CutClass.cxx for an RCF build (e.g. it would skip the
// StPicoEvent*/StPicoTrack*-pointer overloads this Maker needs). This file is that
// second, RCF-flavored version -- copied into place by compile.bash, NOT ../../makefile_toggles.h.
//
// Every toggle below was decided from real evidence, not guessed:
//
//   _STREFMULTCORR_ -- LEFT UNDEFINED, on purpose. Checked two independent sources:
//     (1) This repo's OWN local pipeline (source/PicoBinner.cxx, macros/SetCutClass.C)
//         explicitly never defines it, using raw refMult centrality bins instead
//         ("matching the original's _OO_200_COL_ config" -- see those files' own comments).
//     (2) The reference repo this was ported from (lightflavorspectra_etof)'s real,
//         working top-level makefile_toggles.h has an actual "#ifndef _MAC_OSX_ / #ifdef
//         _OO_200_COL_" (RCF) block -- and it ALSO never defines _STREFMULTCORR_ there,
//         unlike every other collision energy in that same file (7.7/9.2/11.5/14.6/
//         17.3/19.6/27 GeV all define _STREFMULTCORR_ + a matching SL/energy variant on
//         RCF; O+O 200 GeV's RCF block conspicuously does not).
//     So this is a real, deliberate, twice-confirmed choice for O+O 200 GeV specifically,
//     not an oversight -- and it must stay consistent between this embedding Maker and
//     the data-side PicoBinner, or embedding efficiency would be binned into centrality
//     classes differently than the spectra it's meant to correct.
//     (Separately: an actual O+O 200 GeV StRefMultCorr calibration DOES exist upstream --
//     STAR's public star-sw, Param.h mParamStr_ref6[0]/mParamStr_ref7[0], "Run 21 O+O 200
//     GeV", run range [22130029,22144006], which matches this repo's own picoBinner XML
//     run whitelist exactly. That's recorded in compile.bash's own header for whoever
//     needs it later -- it's just not what this build uses, per the two points above.)
//
//   _PICO_READER_SL23c_ -- matches compile.bash's actual `cp -r
//     ../../submodule/PicoDstReader_SL23c ./submodule/` (NOT PicoDstReader_SL24y, the
//     local/Mac build's reader -- see compile.bash's own header comment for why they
//     differ). NOTE: as of this writing, headers/CutClass.h has NO include branch for
//     _PICO_READER_SL23c_ at all -- its "#ifndef _MAC_OSX_" block only knows about
//     _PICO_READER_, _PICO_READER_SL22b_, _PICO_READER_SL22c_, _PICO_READER_SL23a_
//     (which, confusingly, actually points at PicoDstReader_SL24y/ -- a leftover from
//     when that macro was renamed but the #include paths under it were not), and
//     _PICO_READER_SL23d_/_PICO_READER_SL19b_. None of those point at
//     PicoDstReader_SL23c. Left defined here anyway (rather than silently substituting a
//     wrong existing toggle) so the mismatch is visible/greppable -- see compile.bash's
//     header for the flagged follow-up this needs before it will actually compile.
//
//   _ETOF_GOOD_EVENT_FLAG_, _HAS_ETOF_ -- kept on, matching both this repo's local build
//     and the reference repo's real _OO_200_COL_ RCF block (both need eTOF for O+O).
//
//   _PICO_READER_FINISH_BYPASS_, _CPP11_, _BEFORE_ROOT_6p24_ -- carried over unchanged
//     from ../../makefile_toggles.h; these aren't _MAC_OSX_-conditioned in either the
//     local file or the etof reference, so there's no RCF-specific reason to change them.
//
//   Everything _PICOBINNER_*-prefixed, the DEBUG block, and _LIGGETTDRAWEDITS_ are
//   OMITTED here (unlike the top-level file) -- they only gate source/PicoBinner.cxx and
//   analysis macros, none of which this submodule compiles (only CutClass/
//   HistogramUtilities/PhysMath/ParticleInfo/Helix/MattMcEvent/MattMcTrack/
//   StMuAnalysisMaker do, per compile.bash's file list). Add one back only if a build
//   error shows this submodule's sources actually need it.
// =====================================================================================

// _MAC_OSX_ intentionally NOT defined -- this is the RCF/SDCC branch.

#define _PICO_READER_SL23c_
#define _ETOF_GOOD_EVENT_FLAG_
#define _HAS_ETOF_
#define _PICO_READER_FINISH_BYPASS_

#define _CPP11_
#define _BEFORE_ROOT_6p24_

// _STREFMULTCORR_ -- deliberately not defined. See header comment above.

#endif
