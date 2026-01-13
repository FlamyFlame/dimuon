#include "PowhegAlgCoreT.c"
#include "PowhegTruthExtras.c"
#include "PowhegFullSimExtras.c"

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
