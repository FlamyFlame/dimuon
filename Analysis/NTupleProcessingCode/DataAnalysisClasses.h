#include "DimuonDataAlgCoreT.c"
#include "PPExtras.c"
#include "PbPbExtras.c"

class PPAnalysis
	: public DimuonDataAlgCoreT<
				MuonPairPP, MuonPP,
				PPAnalysis,
				PPExtras<PPAnalysis>
		>
	, public PPExtras<PPAnalysis>
{
public:
		PPAnalysis(int run_year_input, int file_batch_input)
				: DimuonDataAlgCoreT(run_year_input, file_batch_input){
				}
};

class PbPbAnalysis
	: public DimuonDataAlgCoreT<
				MuonPairPbPb, MuonPbPb,
				PbPbAnalysis,
				PbPbExtras<PbPbAnalysis>
		>
	, public PbPbExtras<PbPbAnalysis>
{
public:
		PbPbAnalysis(int run_year_input, int file_batch_input)
				: DimuonDataAlgCoreT(run_year_input, file_batch_input){
				}
};
