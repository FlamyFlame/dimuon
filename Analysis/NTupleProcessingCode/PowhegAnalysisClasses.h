#include "PowhegAlgCoreT.c"
#include "PowhegTruthExtras.c"
#include "PowhegFullSimExtras.c"
#include "PowhegFullSimOverlayExtras.c"

class PowhegTruthAnalysis
  : public PowhegAlgCoreT<
        MuonPairPowhegTruth, MuonPowhegTruth,
        PowhegTruthAnalysis,
        PowhegTruthExtras<MuonPairPowhegTruth,PowhegTruthAnalysis>
    >
  , public PowhegTruthExtras<MuonPairPowhegTruth,PowhegTruthAnalysis>
{
public:
  PowhegTruthAnalysis(int file_batch_input, std::string mc_mode_input, bool useLocal_input = false)
    : PowhegAlgCoreT(file_batch_input, mc_mode_input, 0, useLocal_input){
        }
};

class PowhegFullSimAnalysisNoTruth
  : public PowhegAlgCoreT<
        MuonPairPowhegFullSimNoTruth, MuonPowhegFullSimNoTruth,
        PowhegFullSimAnalysisNoTruth,
        PowhegFullSimExtras<MuonPairPowhegFullSimNoTruth, MuonPowhegFullSimNoTruth, PowhegFullSimAnalysisNoTruth>
    >
  , public PowhegFullSimExtras<MuonPairPowhegFullSimNoTruth, MuonPowhegFullSimNoTruth, PowhegFullSimAnalysisNoTruth>
{
public:
  PowhegFullSimAnalysisNoTruth(int file_batch_input, std::string mc_mode_input, int run_year_input = 17, bool useLocal_input = false)
    : PowhegAlgCoreT(file_batch_input, mc_mode_input, run_year_input, useLocal_input){
        }
};

class PowhegFullSimAnalysisWTruth
  : public PowhegAlgCoreT<
        MuonPairPowhegFullSimWTruth, MuonPowhegFullSimWTruth,
        PowhegFullSimAnalysisWTruth,
        PowhegFullSimExtras<MuonPairPowhegFullSimWTruth, MuonPowhegFullSimWTruth, PowhegFullSimAnalysisWTruth>,
        PowhegTruthExtras<MuonPairPowhegFullSimWTruth, PowhegFullSimAnalysisWTruth>
    >
  , public PowhegFullSimExtras<MuonPairPowhegFullSimWTruth, MuonPowhegFullSimWTruth, PowhegFullSimAnalysisWTruth>
  , public PowhegTruthExtras<MuonPairPowhegFullSimWTruth, PowhegFullSimAnalysisWTruth>
{
public:
  PowhegFullSimAnalysisWTruth(int file_batch_input, std::string mc_mode_input, int run_year_input = 17, bool useLocal_input = false)
    : PowhegAlgCoreT(file_batch_input, mc_mode_input, run_year_input, useLocal_input){
        }
};

class PowhegFullSimOverlayAnalysisNoTruth
  : public PowhegAlgCoreT<
        MuonPairPowhegFullSimOverlayNoTruth, MuonPowhegFullSimOverlayNoTruth,
        PowhegFullSimOverlayAnalysisNoTruth,
        PowhegFullSimExtras<MuonPairPowhegFullSimOverlayNoTruth, MuonPowhegFullSimOverlayNoTruth, PowhegFullSimOverlayAnalysisNoTruth>,
        PowhegFullSimOverlayExtras<PowhegFullSimOverlayAnalysisNoTruth>
    >
  , public PowhegFullSimExtras<MuonPairPowhegFullSimOverlayNoTruth, MuonPowhegFullSimOverlayNoTruth, PowhegFullSimOverlayAnalysisNoTruth>
  , public PowhegFullSimOverlayExtras<PowhegFullSimOverlayAnalysisNoTruth>
{
public:
  PowhegFullSimOverlayAnalysisNoTruth(int file_batch_input, std::string mc_mode_input, int run_year_input = 24, bool useLocal_input = false)
    : PowhegAlgCoreT(file_batch_input, mc_mode_input, run_year_input, useLocal_input){
        }
};

class PowhegFullSimOverlayAnalysisWTruth
  : public PowhegAlgCoreT<
        MuonPairPowhegFullSimOverlayWTruth, MuonPowhegFullSimOverlayWTruth,
        PowhegFullSimOverlayAnalysisWTruth,
        PowhegFullSimExtras<MuonPairPowhegFullSimOverlayWTruth, MuonPowhegFullSimOverlayWTruth, PowhegFullSimOverlayAnalysisWTruth>,
        PowhegFullSimOverlayExtras<PowhegFullSimOverlayAnalysisWTruth>,
        PowhegTruthExtras<MuonPairPowhegFullSimOverlayWTruth,PowhegFullSimOverlayAnalysisWTruth>
    >
  , public PowhegFullSimExtras<MuonPairPowhegFullSimOverlayWTruth, MuonPowhegFullSimOverlayWTruth, PowhegFullSimOverlayAnalysisWTruth>
  , public PowhegFullSimOverlayExtras<PowhegFullSimOverlayAnalysisWTruth>
  , public PowhegTruthExtras<MuonPairPowhegFullSimOverlayWTruth,PowhegFullSimOverlayAnalysisWTruth>
{
public:
  PowhegFullSimOverlayAnalysisWTruth(int file_batch_input, std::string mc_mode_input, int run_year_input = 24, bool useLocal_input = false)
    : PowhegAlgCoreT(file_batch_input, mc_mode_input, run_year_input, useLocal_input){
        }
};
