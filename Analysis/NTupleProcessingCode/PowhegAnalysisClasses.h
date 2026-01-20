#include "PowhegAlgCoreT.c"
#include "PowhegTruthExtras.c"
#include "PowhegFullSimExtras.c"
#include "PowhegFullSimOverlayExtras.c"

class PowhegTruthAnalysis
  : public PowhegAlgCoreT<
        MuonPairPowhegTruth, MuonPowhegTruth,
        PowhegTruthAnalysis,
        PowhegTruthExtras<PowhegTruthAnalysis>
    >
  , public PowhegTruthExtras<PowhegTruthAnalysis>
{};

class PowhegFullSimAnalysisNoTruth
  : public PowhegAlgCoreT<
        MuonPairPowhegFullSimNoTruth, MuonPowhegFullSimNoTruth,
        PowhegFullSimAnalysisNoTruth,
        PowhegFullSimExtras<PowhegFullSimAnalysisNoTruth>
    >
  , public PowhegFullSimExtras<PowhegFullSimAnalysisNoTruth>
{};

class PowhegFullSimAnalysisWTruth
  : public PowhegAlgCoreT<
        MuonPairPowhegFullSimWTruth, MuonPowhegFullSimWTruth,
        PowhegFullSimAnalysisWTruth,
        PowhegFullSimExtras<PowhegFullSimAnalysisWTruth>,
        PowhegTruthExtras<PowhegFullSimAnalysisWTruth>
    >
  , public PowhegFullSimExtras<PowhegFullSimAnalysisWTruth>
  , public PowhegTruthExtras<PowhegFullSimAnalysisWTruth>
{};

class PowhegFullSimOverlayAnalysisNoTruth
  : public PowhegAlgCoreT<
        MuonPairPowhegFullSimOverlayNoTruth, MuonPowhegFullSimOverlayNoTruth,
        PowhegFullSimOverlayAnalysisNoTruth,
        PowhegFullSimExtras<PowhegFullSimOverlayAnalysisNoTruth>,
        PowhegFullSimOverlayExtras<PowhegFullSimOverlayAnalysisNoTruth>
    >
  , public PowhegFullSimExtras<PowhegFullSimOverlayAnalysisNoTruth>
  , public PowhegFullSimOverlayExtras<PowhegFullSimOverlayAnalysisNoTruth>
{};

class PowhegFullSimOverlayAnalysisWTruth
  : public PowhegAlgCoreT<
        MuonPairPowhegFullSimOverlayWTruth, MuonPowhegFullSimOverlayWTruth,
        PowhegFullSimOverlayAnalysisWTruth,
        PowhegFullSimExtras<PowhegFullSimOverlayAnalysisWTruth>,
        PowhegFullSimOverlayExtras<PowhegFullSimOverlayAnalysisWTruth>,
        PowhegTruthExtras<PowhegFullSimOverlayAnalysisWTruth>
    >
  , public PowhegFullSimExtras<PowhegFullSimOverlayAnalysisWTruth>
  , public PowhegFullSimOverlayExtras<PowhegFullSimOverlayAnalysisWTruth>
  , public PowhegTruthExtras<PowhegFullSimOverlayAnalysisWTruth>
{};


// class PowhegFullSimNTupleFirstPass : public virtual PowhegNTupleFirstPass{
// public: 
//     PowhegFullSimNTupleFirstPass(int file_batch_input, std::string mc_mode_input)
//         : PowhegNTupleFirstPass(file_batch_input, mc_mode_input, true){
//             output_single_muon_tree = false;
//         }

//     ~PowhegFullSimNTupleFirstPass(){}
// };


