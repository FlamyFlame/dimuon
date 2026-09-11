#ifndef __PlotMCDataComprBaseClass_h__
#define __PlotMCDataComprBaseClass_h__

#include <TROOT.h>
#include "McDataComprColors.h"
#include "McDataComprConfig.h"
#include "McDataComprRatio.h"
#include <TFile.h>
#include <TSystem.h>
#include <string.h>
#include <stdlib.h>
#include "vector"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include <algorithm>
#include <cmath>
// <iostream> MUST come before helper_functions.c: that file uses std::cout but does not
// include <iostream> itself, so ACLiC compilation fails without this.
#include <iostream>
#include <string>
#include "../helper_functions.c"
#include "../DimuonPlottingBaseClass.cxx"
#include "../../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../../MuonObjectsParamsAndHelpers/PPBaseClass.h"

// =================================================================================================
// PlotMCDataComprBaseClass -- shared plumbing of the pp24 MC-vs-data comparison plot set.
//
// The set lives in TWO families, in two SUBDIRECTORIES of plots/mc_data_compr/ (see
// docs/tracking/mc_data_compr_signal_generic_split.md):
//
//   signal/   the data-like single-b SIGNAL REGION. Data and Pythia carry the SAME cuts
//             (m_uu in (1.08,2.9), pair pT > ParamsSet::signal_pair_pt_min (9 GeV since
//             2026-09-08), fiducial gap cut on both muons; the MC in
//             TRUTH quantities), so every efficiency correction on the data side is applied
//             inside the region in which it was measured. No POWHEG, no jacobian variants.
//   generic/  no signal-region cut on either side, gap cut only. Shapes over the full phase
//             space; the absolute level is NOT a cross-section (see generic/README.md).
//
// The Pythia partner is the pp24-condition FULLSIM FULL sample, NOT the private sample and NOT
// the Pythia TRUTH full sample. Physics reason: the truth full sample's weight carries the Pb
// isospin average 4:6:6:9 unconditionally, so its absolute sigma is a Pb-averaged NN cross
// section, not a pp one (16.6 % here) -- the wrong object to place beside pp24 data. The fullsim
// pp24 sample is pp beam only with isospin weight 1 and is AMI-weighted, so the ONLY conversion
// it needs is nb -> pb (the AMI `crossSection` is in nb, the data in pb). The hand-tuned
// `TRUTH_COMBINE_RENORM = 1e3` that the private sample required is therefore DELETED.
//
// Every 1D panel carries an MC/data RATIO PAD (McDataComprRatio.h). The shared registered binning
// (Design Decision D2) is what makes it possible at all: before it, data and MC sat on different
// axes and only their shapes could be compared.
// =================================================================================================

class PlotMCDataComprBaseClass : public DimuonPlottingBaseClass{
protected:
    // THREE input samples. There used to be four, with POWHEG entered TWICE (`powheg_bb` and
    // `powheg_cc` pointing at two 2025 legacy files). The current producer
    // (RDFBasedHistFillingPowhegTruth.cxx) writes ONE combined file, and the flavour summing that
    // the two entries existed for is what produced the 2.881x over-count -- so both collapse to a
    // single `powheg` entry reading the INCLUSIVE histogram.
    static const int s_nDtTypes = 3;

    enum DataType{
        powheg,
        pythia,
        pp_2024_2mu4
    };

    double pi = acos(-1.0);
    // The pp24 2mu4 luminosity is READ from PPBaseClass, the single source of truth the
    // cross-section itself uses. It was hard-coded here as 1/410.815 pb^-1 and was left behind
    // when the value was corrected to 400.412 pb^-1 on 2026-06-15 (commit 009f3b8) -- so the
    // data curve in every one of these plots was 2.6 % low relative to the cross-section it is
    // meant to be compared against. Retyping a normalization is exactly how that happens.
    // Filled in initialize(), because two entries are family-dependent:
    //   [pythia] = NB_TO_PB, the AMI weight being in nb;
    //   [powheg] = 1, and that is NOT an oversight -- see PowhegNormFactor();
    //   [data]   = 1/L_int, UNLESS the histogram is already d(sigma) (`data_already_dsigma`).
    std::array<float, s_nDtTypes> norm_factor = {1., 1., 1.};
    std::array<Color_t, s_nDtTypes> colors = {kRed, kBlue, kGreen+2};

    // Colour convention of the SIGNAL family lives at namespace scope in McDataComprColors.h,
    // because `plot_mc_data_pair_pt_in_eta.cxx` draws the same family without deriving from this
    // class. These aliases keep the existing call sites short.
    static constexpr Color_t kSignalDataColor  = McDataComprColors::kSignalData;
    static constexpr Color_t kSignalMcColor    = McDataComprColors::kSignalMc;
    static constexpr Color_t kSignalMcRefColor = McDataComprColors::kSignalMcRef;

    // The AMI `crossSection` carried by the Pythia fullsim weight is in nb; the data is
    // normalized to pb. This is the ONLY factor the fullsim Pythia needs -- it is physics, not a
    // tuned constant.
    static double NbToPb() { return 1.0e3; }

    // -----------------------------------------------------------------------------------------
    // POWHEG normalization is EXACTLY 1, and it must NOT be given Pythia's nb -> pb x1000.
    //
    // Established from the producer, not from a fit: `PowhegAlgCoreT.c:209-212` fills
    // `weight_norm = weight / N_gen` with `weight = EventWeights[0] * filter_effcy`, and the
    // POWHEG event weight is already a cross-section in **pb** -- not the nb of an AMI
    // `crossSection` field. There is therefore no unit conversion left to apply.
    //
    // Sanity check on the files as they stand: POWHEG OS total 14122 pb against pp24 data OS
    // 28055 pb. Multiplying by 1000 would put POWHEG three orders of magnitude above the data,
    // which is how anyone "fixing" this by symmetry with Pythia would find out -- but only after
    // the plot had been believed. Leave it at 1.
    // -----------------------------------------------------------------------------------------
    static double PowhegNormFactor() { return 1.0; }

    std::string pythia_with_data_resonance_cuts_suffix;

    std::array<TFile*, s_nDtTypes> f;
    std::array<bool, s_nDtTypes> is_data = {false, false, true};
    std::array<std::string, s_nDtTypes> dt_paths;
    std::array<std::string, s_nDtTypes> fnames;
    // "POWHEG" = bb + cc, as it must be. The two are POWHEG GENERATOR MODES -- each requires a
    // b-bbar (resp. c-cbar) pair in the NLO hard scattering -- so they are two exclusive
    // contributions to one process and their cross sections ADD. Every plotting code is REQUIRED
    // to sum them (.claude/conventions/atlas-plotting.md, MANDATORY; Analysis/docs/powheg.md).
    //
    // It read "POWHEG bb" from 2026-08-25 to 2026-09-11 because muon_pairs_powheg_cc_truth.root
    // did not exist (parts 2 and 6 had never been produced) -- an honest label on an incomplete
    // curve. All 6 cc parts exist since the 2026-09-10 rerun; cc is merged and
    // histograms_powheg_truth.root now carries bb + cc, each normalized to its OWN N_gen
    // (RDFBasedHistFillingPowheg::CreateBaseRDFsPowhegCommon, per-sample via DefinePerSample).
    // Measured effect on THIS curve: the generic family grows x1.065; the single_b family is
    // unchanged, because `from_same_b` is false for every cc pair -- charm events contain no
    // b-hadrons, so the charm sample contributes nothing to the single-b signal. That zero is a
    // physics result, and it is the reason the sample may not simply be dropped.
    std::array<std::string, s_nDtTypes> dtTitles = {"POWHEG", "Pythia", "pp data 2024"};

    // ---------------------------------------------------------------------------------------
    // THE ONE PLACE where the mc_data_compr histogram keys live.
    // A producer-side rename is a one-line edit in the table below (PlotMCDataComprBaseClass.c)
    // plus, at most, the key builders directly under it -- and nowhere else in the plot code.
    struct McDataObservable {
        std::string kin;         // plot id; also the PNG file-name stem
        std::string x_title;     // x axis title
        std::string y_symbol;    // the symbol under d in d(sigma)/d<symbol>
        std::string y_unit;      // "pb" for a dimensionless observable, "pb GeV^{-1}" otherwise
        std::string data_var;    // variable token inside the DATA histogram key
        std::string mc_var;      // variable token inside the Pythia fullsim histogram key
        std::string powheg_var;  // variable token inside the POWHEG histogram key ("" = no POWHEG)
        bool logx = false;
    };

    static const std::vector<McDataObservable>& Observables();
    static const McDataObservable& Observable(const std::string& kin);

    // The y axis title. `jacobian` is not cosmetic: a 1/dR-WEIGHTED spectrum is a different
    // quantity from the unweighted one, and the two panels sit side by side in generic/. They
    // used to carry byte-identical y titles, which made them indistinguishable from the image.
    static std::string YTitle(const McDataObservable& o, bool jacobian){
        if (jacobian)
            return "(1/" + o.y_symbol + ") d#sigma/d" + o.y_symbol + " [" + o.y_unit + "]";
        return "d#sigma/d" + o.y_symbol + " [" + o.y_unit + "]";
    }

    // --- key builders -----------------------------------------------------------------------
    // SIGNAL family. The data histograms are already d(sigma) in pb (their weight carries 1/L),
    // hence the `_dsigma` stem; do NOT re-apply the luminosity factor to them.
    static std::string SignalDataKey(const McDataObservable& o, int isign){
        return "h1d_crossx_" + o.data_var + "_w_signal_cuts_" + SignDataTag(isign) + "_dsigma";
    }
    // SIGNAL family, MC. The pair CATEGORY differs between the two pads and that is physics:
    //   OS pad -> `_single_b_pass_signal_truth`: truth OS pairs from the SAME b, the single-b
    //             signal the measurement is after. (`_single_b` is built in the producer as
    //             df_op.Filter("from_same_b"), so it is OS by construction.)
    //   SS pad -> `_ss_pass_signal_truth`: ALL truth same-sign pairs in the signal region, with
    //             no `from_same_b` -- `from_same_b` has no same-sign counterpart, and the data
    //             SS pad is the combinatorial-background estimate, whose MC partner is the
    //             inclusive SS yield.
    //   OS pad, second curve -> `_op_pass_signal_truth`: ALL truth OS pairs. This is the
    //             like-for-like partner of the data OS yield, which still contains
    //             gluon-splitting and combinatorial background.
    static const char* McSignalCatSingleB(){ return "_single_b_pass_signal_truth"; }
    static const char* McSignalCatAllOS()  { return "_op_pass_signal_truth"; }
    static const char* McSignalCatAllSS()  { return "_ss_pass_signal_truth"; }
    static std::string SignalMcKey(const McDataObservable& o, const std::string& cat){
        return "h_" + o.mc_var + cat;
    }

    // GENERIC family. The gap-cut generic filter is named identically in the Pythia fullsim and
    // in the POWHEG truth producer, so ONE pair of category tokens serves both.
    static const char* McGenericCatOS(){ return "_op_gapcut_truth"; }
    static const char* McGenericCatSS(){ return "_ss_gapcut_truth"; }
    static std::string GenericDataKey(const McDataObservable& o, int isign, bool jacobian){
        return "h_" + o.data_var + "_" + SignDataTag(isign)
             + (jacobian ? "_jacobian_corrected" : "");
    }
    static std::string GenericMcKey(const McDataObservable& o, int isign, bool jacobian){
        return "h_" + o.mc_var + (isign == 0 ? McGenericCatSS() : McGenericCatOS())
             + (jacobian ? "_jacobian_corrected" : "");
    }
    // ⚠ OPEN (2026-09-11, user): POWHEG here carries NEITHER `from_same_b` NOR the truth signal
    // cuts, while the Pythia signal family carries both. The requirement is that POWHEG and Pythia
    // both require `from_same_b` AND the data signal cuts evaluated on TRUTH quantities. As it
    // stands:
    //   Pythia  `_single_b_pass_signal_truth` = from_same_b && truth_minv in (1.08,2.9)
    //                                        && truth_pair_pt > 9 && gap cuts        <-- correct
    //   POWHEG  `_op_gapcut_truth`            = gap cuts ONLY                         <-- NOT
    // and POWHEG is drawn only in the GENERIC family: plot_mc_data_compr_signal.cxx pushes three
    // Pythia series and no POWHEG one at all. The POWHEG truth producer does build the full signal
    // selection (RDFBasedHistFillingPowhegTruth::FillHistogramsSignalAcceptance) but consumes it
    // internally for the acceptance RATIO and never writes it as a filter family, which is why
    // there is no key to read. Closing this needs a `_single_b_pass_signal_truth` family in the
    // POWHEG truth producer + a category argument here, mirroring SignalMcKey.
    //
    // POWHEG, 2026-08-25 rework. The INCLUSIVE histogram is read directly. The producer also
    // writes `..._flavor_binned_{single_b,both_from_b,both_from_c,others}`, which sums to this
    // same inclusive to <= 4e-12 -- but nothing here sums anything any more: the old code summed
    // three OVERLAPPING classification axes and over-counted by 2.881x. One key, no arithmetic.
    static std::string PowhegKey(const McDataObservable& o, int isign){
        if (o.powheg_var.empty()) return "";
        return "h_" + o.powheg_var + (isign == 0 ? McGenericCatSS() : McGenericCatOS());
    }

    static const char* SignDataTag(int isign){ return (isign == 0) ? "ss" : "op"; }

    // --- fetching ---------------------------------------------------------------------------
    TH1D* GetHist1D(int idt, const std::string& base_name) const;   // direct lookup, clones

    // --- axis guard -------------------------------------------------------------------------
    // Data, Pythia and POWHEG are now all booked from the SAME registered binning, so a
    // divergence is a BUG on the producer side, never something to silently plot. POWHEG was
    // exempt until 2026-08-25 (its 2025 file had dR on 100 bins against the data's 40); the
    // reworked producer puts it on the shared axes, verified edge by edge in all 12 comparisons,
    // so it now goes through this guard like every other curve. Compare with a RELATIVE
    // tolerance: one side is built as (nbins, min, max) and the other from a generated edge
    // vector, and those agree only to double rounding.
    static void AssertSameAxis1D(const TH1* a, const TH1* b, const std::string& ctx);

    virtual void initialize();
    virtual void set_legend_position();

    // --- output directory --------------------------------------------------------------------
    // Every save path used to be a string literal inside the SaveAs call, so one macro could
    // only ever write to one place. `output_subdir` ("" = the legacy top level) is prepended by
    // OutputPath(), which also creates the directory.
    static const char* PlotBaseDir(){
        return "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/mc_data_compr/";
    }
    std::string OutputPath(const std::string& filename) const;

public:
    bool pythia_with_data_resonance_cuts = false;

    // SIGNAL family: the data histograms already carry 1/L_int, so norm_factor[data] must be 1.
    bool data_already_dsigma = false;

    // MUON WORKING POINT. Nominal = Tight. See McDataComprConfig::MuonWP for exactly what it
    // switches (the DATA input file, selection and corrections together) and what it cannot
    // switch (the MC, whose histograms are truth quantities and have no reconstruction WP).
    McDataComprConfig::MuonWP muon_wp = McDataComprConfig::MuonWP::Tight;

    std::string output_subdir;   // "" | "signal" | "generic"

    // Log y for the 1D distributions. It lives HERE, not in the derived plotter, because
    // set_legend_position() below must know about it: the legend's default corner depends on it.
    bool logy = false;
    std::array<float,4> legend_position_same_sign = {s_NULL, s_NULL, s_NULL, s_NULL};
    std::array<float,4> legend_position_opp_sign = {s_NULL, s_NULL, s_NULL, s_NULL};
};

#endif
