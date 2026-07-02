#ifndef _MAKEFILE_TOGGLES_
#define _MAKEFILE_TOGGLES_

// =====================================================================================
// makefile_toggles.h -- lightflavorspectra_OO200
//
// This is a simplified, hardcoded version of the original repo's makefile_toggles.h,
// specialized for ONE configuration only: O+O at 200 GeV, ColliderCenter, running
// locally (not on SDCC/RCF). There is no energy switch and no _MAC_OSX_/RCF branch --
// this repo only ever builds this one way.
//
// Every toggle below is exactly what _OO_200_COL_ selected in the original repo's
// "Not on RCF" (_MAC_OSX_) branch, confirmed by reading the original makefile_toggles.h
// and source/PicoBinner.cxx line by line. Toggles that were never defined for this
// configuration (_STREFMULTCORR_, _JIA_CHEN_19p6_Modification_, _PICOBINNER_VZ_STUDY_,
// _PICOBINNER_2D_DEDX_BTOF_, _PICO_BINNER_DEBUG_, _USE_CALIB_CLASS_, _MATTDST_MODE_,
// FixedTarget-anything) have been removed rather than left commented out, since the
// corresponding dead code branches were also removed from the ported source files.
//
// _PICO_READER_SL23a_ was changed to _PICO_READER_SL24y_ on 2026-07-02, to match the
// production data's actual STAR library (SL24y). The PicoDst reader submodule itself
// was re-vendored from STAR's real SL24y source (submodule/PicoDstReader_SL24y/,
// replacing submodule/PicoDstReader_SL23a/) -- see that directory and
// source/PicoBinner.cxx's own header comment for what actually changed between the two
// versions (most classes are byte-identical; StPicoArrays/StPicoEvent/StPicoMcTrack
// are not, and StPicoEvent's eTOF good-event-flag API changed in a way that required a
// PicoBinner.cxx call-site update -- see PicoBinner.cxx). NOTE: unlike the toggles
// above, this one is only read in dead code paths in headers/CutClass.h (guarded by
// `#ifndef _MAC_OSX_`, which is always false here) -- every actual PicoDstReader
// include elsewhere in this repo is a hardcoded literal path, not conditioned on this
// macro. It's kept in sync anyway for clarity/consistency, not because anything
// actually branches on it.
// =====================================================================================

#define _PICO_READER_
#define _PICO_READER_SL24y_
#define _ETOF_GOOD_EVENT_FLAG_
#define _HAS_ETOF_
#define _PICO_READER_FINISH_BYPASS_

#define _CPP11_
// _MAC_OSX_ is NOT an "Apple hardware" flag -- in the original repo it means
// "building locally, not on RCF/SDCC" and selects the local code paths in
// CutClass.h/.cxx and HistogramUtilities.h/.cxx (e.g. it skips the StPicoEvent*/
// StPicoTrack*-pointer convenience overloads that PicoBinner never calls anyway,
// and enables calcGausData()). Since this repo only ever runs locally, this must
// stay defined -- an earlier draft of this file mistakenly dropped it as if it
// were RCF-only infrastructure; that was wrong and has been corrected.
#define _MAC_OSX_
#define _BEFORE_ROOT_6p24_

#define _PICOBINNER_BTOF_ISOLATED_HISTOS_
#define _PICOBINNER_SIMPLE_SPECTRA_
#define _PICOBINNER_VZ_CALIBRATION_
#define _PICOBINNER_ETA_PT_EFF_
#define _PICOBINNER_USE_BICHSEL_
#define _DRAWSIMPLE_FITNAME_MOD_

#define _SPECTRA_FITTER_NUMERIC_BLASTWAVE_INTEGRATION_

//########    DEBUG STATEMENTS  ##########
// Kept only the ones that were actually enabled upstream and are still relevant
// to the stages we've ported so far. Add more back only if a later stage needs one
// that's still commented out here.
#define _DNDY_ANALYSIS_DEBUG_
#define _EFFICIENCY_FITTER_DEBUG_
#define _RAWSPECMOD_DEBUG_
#define _SPECTRA_FITTER_DEBUG_
#define _SYSTEMATIC_ERROR_STUDIES_DEBUG_
#define _ZFITTER_DEBUG_

#define _LIGGETTDRAWEDITS_

#endif
