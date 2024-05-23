#include "Pythia8/Pythia.h"
#include "Pythia8/PartonDistributions.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TEnv.h"
#include <assert.h>

//#include "fastjet/PseudoJet.hh"
//#include "fastjet/ClusterSequence.hh"

// Nucleus is not thoroughly checked yet
// Define a nuclear EPA class to replace the muon EPA flux, implements the PDF interface
class Nucleus2gamma: public Pythia8::PDF 
{
  public:
    Nucleus2gamma(int idBeamIn): PDF(idBeamIn) {}

    void xfUpdate(int id, double x, double Q2)
    {
      double bmin = 2* 6.62; // Minimum impact parameter. Units: length

      double Z = 82.0; // Atomic number for lead
      double m_eff = 0.9315; // Effective nucleon mass. Units: mass (duh)
      double alphaEM = 0.007297353080; // Dimensionless electromagnetic coupling constant
      double hbarc = 0.197; // Dimensionless (in natural units) constant
      double xi = x * m_eff * bmin / hbarc; // Dimensionless effective momentum fraction
      double bK0 = TMath::BesselK0(xi);
      double bK1 = TMath::BesselK1(xi);
      // Equation 5 from Bertulani, Klein, Nystrad, "Physics of Ultra-Peripheral Nuclear Collisions" (2005)
      double intB = xi * bK1 * bK0 - 0.5 * xi * xi * (bK1*bK1-bK0*bK0);
      
      // Protected attribute inherited from PDF -- plugs into functions from PDF interface
      xgamma = 2. * Z * Z * alphaEM * intB / TMath::Pi(); // Finishes up eqn 5 from above taking beta=1.
    }
};

int main(int argc, char **argv) 
{

  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " <integer> <integer> <integer>" << std::endl;
      return 1;
  }

  int nEvent =100000; // The total number of events we will end up keeping
  int nTried =nEvent; // To be modified -- the total number of events generated without cuts
  int kinematicRange = atoi(argv[2]); // Configure following variables depending on the intended kinematic range (in terms of jet pT intervals). 0: low, 1: medium, 2: high, 3: very high  

  assert(kinematicRange >= 0 && kinematicRange <= 3);

  // define pTHat ranges to apply Phase Space Cuts on the hard process
  // pTHatMin, pTHatMax cuts are to be applied in the rest frame of the hard subprocess
  // with crossx adjusted accordingly
  // since event crossx tends to decrease with pT
  // if we want to study the kinematics at high pT without crazily large statistical error
  // we can generate events in different kinematic ranges
  // and be careful about overall weighting/normalization when combining the events

  double pTHatMin, pTHatMax; // Defining the parameters which need to change with pT cut interval
  if(kinematicRange == 0) // Low pT hard scattering -- filtering is harder here!
  {
  	pTHatMin=3.0;
  	pTHatMax=10.0;
  }
  else if(kinematicRange == 1) // Mid pT hard scattering
  {
  	pTHatMin=10.0;
  	pTHatMax=25.0;
  }
  else if(kinematicRange == 2) // High pT hard scattering
  {
  	pTHatMin=25.0;
  	pTHatMax=60.0;
  }
  else if(kinematicRange == 3) // Very high pT hard scattering
  {
  	pTHatMin=60.0;
  	pTHatMax=800.0;
  }
  // minimum jet pT
  //double jetCut = 15.0;
  // Center of mass energy range
  //double Wmax = 1142; 
  // Definitions to be changed later
  double efficiency=1;
  double totalSigma=0;

  // Initialize Pythia and its RNG
  Pythia8::Pythia pythia;
  pythia.readString("Random:setSeed=on");

  // Set the random seed to the number supplied as an argument
  char buffer[50];
  sprintf(buffer,"Random:seed=%s",argv[1]);
  //sprintf(buffer,"Random:seed=88003");
  pythia.readString(buffer);

  // Just some settings
  pythia.readString("Init:showChangedSettings = off");
  pythia.readString("Init:showChangedParticleData = off");
  pythia.readString("Next:numberCount = 0");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Setup objects to reference information about the simulation
  Pythia8::Settings& settings = pythia.settings;
  // Pythia8::Info& info         = pythia.info;
  Pythia8::Info info         = pythia.info;

  //A14 tune, base tune is Monash 2013
  pythia.readString("Tune:ee = 7");
  pythia.readString("Tune:pp = 14");
  pythia.readString("SpaceShower:rapidityOrder = on");
  pythia.readString("SigmaProcess:alphaSvalue = 0.140");
  pythia.readString("SpaceShower:pT0Ref = 1.56");
  pythia.readString("SpaceShower:pTmaxFudge = 0.91");
  pythia.readString("SpaceShower:pTdampFudge = 1.05");
  pythia.readString("SpaceShower:alphaSvalue = 0.127");
  pythia.readString("TimeShower:alphaSvalue = 0.127");
  pythia.readString("BeamRemnants:primordialKThard = 1.88");
  pythia.readString("MultipartonInteractions:pT0Ref = 2.09");
  pythia.readString("MultipartonInteractions:alphaSvalue = 0.126");
  //pythia.readString("BeamRemnants:reconnectRange  = 1.71");

  pythia.readString("Beams:eCM = 5020."); // CoM energy of the interaction (5.02 TeV)
  
  int beam_type = atoi(argv[3]);
  assert(beam_type >= 0 && beam_type <= 3);
  if(beam_type == 0){ // Runtime option: proton or neutron simulations
    pythia.readString("Beams:idA = 2212"); // Monte carlo code proton -- may change later
    pythia.readString("Beams:idB = 2212");
  }
  else if(beam_type == 1){ // pn
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2112");
  }
  else if(beam_type == 2){ // np
    pythia.readString("Beams:idA = 2112");
    pythia.readString("Beams:idB = 2212");
  }
  else{ // nn
    pythia.readString("Beams:idA = 2112");
    pythia.readString("Beams:idB = 2112");
  }
  pythia.readString("PDF:pSet = LHAPDF6:nCTEQ15npFullNuc_208_82/0001"); // Set PDFs
  // Enable photon flux sampling from a PDF. Currently set to muon EPA but will be re-pointed to custom nuclear
  //pythia.readString("PDF:lepton2gamma = on");

  // Use the nucleus distribution we defined
  // Following line removes photon virtuality test. Required for custom PDF. Small error since Q2 is small
  //pythia.readString("Photon:sampleQ2 = off"); 
  //Nucleus2gamma* photonFlux = new Nucleus2gamma(13);
  //pythia.setPhotonFluxPtr(photonFlux, 0);

  // Apply the above constraints
  pythia.readString("Photon:Q2max = 1.0"); // Upper limit on photon virtuality
  //pythia.readString("Photon:Wmin  = 10.0"); // Set lower limit for the invariant mass of 2-photon systems
  //pythia.readString("Photon:Wmax = 1142."); // Set upper limit for the invariant mass of 2-photon systems
  settings.parm("PhaseSpace:pTHatMin", pTHatMin); // Apply kinematics constraint on hard processes
  settings.parm("PhaseSpace:pTHatMax", pTHatMax); // Apply kinematics constraint on hard processes

  // Turn on all channels for processes we want
  //pythia.readString("HardQCD:hardccbar = on");
  //pythia.readString("HardQCD:hardbbbar = on");
  pythia.readString("HardQCD:all = on"); // Enable all hard QCD interactions (quark and gluon couplings)
  //if(kinematicRange == 0)
  //pythia.readString("SoftQCD:nonDiffractive = on"); // Enable soft QCD processes which are not diffractive
  //pythia.readString("PartonLevel:MPI = off");
  //pythia.readString("PhotonParton:all = on"); // Enable relevant interactions between photons and partons
  //settings.mode("Photon:ProcessType", 0); // Enables all 4 types of photon-photon interactions

  settings.listChanged(); // Print all settings modifications


  //settings.listAll();
  if(!pythia.init()) // Initialize pythia with all of these settings
  {
    std::cerr << "Failed to initialize" << std::endl;
    exit(0);
  }

  
  Pythia8::Event& event=pythia.event; // Make a reference to the event we are generating
  // Set up an anti-kT jet clustering algorithm with R=0.4 and energy recombination scheme  
  //fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, 0.4, fastjet::E_scheme, fastjet::Best);

  TFile *file = TFile::Open("pytreeEff.root","recreate"); // Place to write out
  TTree *T = new TTree("PyTree","ev1 Tree"); // Tree we will write data into
  std::vector<int>* b_barcode=new std::vector<int>(); // Individual identifiers for particles
  std::vector<int>* b_status=new std::vector<int>(); // What role it plays in calculation -- virtual, real, final state, input, etc.
  std::vector<int>* b_id=new std::vector<int>(); // Particle identity: p, n, etc
  std::vector<int>* b_mother1=new std::vector<int>(); // Particle first mother particle
  std::vector<int>* b_mother2=new std::vector<int>(); // Particle last mother particle
  std::vector<double>* b_pt=new std::vector<double>(); // Kinematics
  std::vector<double>* b_eta=new std::vector<double>(); // Kinematics
  std::vector<double>* b_phi=new std::vector<double>(); // Kinematics
  std::vector<double>* b_m=new std::vector<double>(); // Kinematics

  std::vector<int>* b_muon_inds=new std::vector<int>(); // muon indices - for easier pairings

  std::vector<double>* b_muon1_pt=new std::vector<double>();
  std::vector<double>* b_muon1_eta=new std::vector<double>();
  std::vector<double>* b_muon1_phi=new std::vector<double>();
  std::vector<int>* b_muon1_ch=new std::vector<int>();
  std::vector<int>* b_muon1_id=new std::vector<int>();
  std::vector<int>* b_muon1_bar=new std::vector<int>();
  std::vector<double>* b_muon2_pt=new std::vector<double>();
  std::vector<double>* b_muon2_eta=new std::vector<double>();
  std::vector<double>* b_muon2_phi=new std::vector<double>();
  std::vector<int>* b_muon2_ch=new std::vector<int>();
  std::vector<int>* b_muon2_id=new std::vector<int>();
  std::vector<int>* b_muon2_bar=new std::vector<int>();


  double b_weight=0.; // Event weighting, usually 0
  double b_QHard=0; // Momentum transfer
  int    b_nmuons=0;
  int    b_nmuonPairs=0;

  // Branch the tree for all simulation parameters we want to store
  T->Branch("event_weight",&b_weight);
  T->Branch("QHard",&b_QHard);
  T->Branch("nmuons",&b_nmuons);
  T->Branch("nmuonPairs",&b_nmuonPairs);

  T->Branch("truth_barcode",&b_barcode);
  T->Branch("truth_status",&b_status);
  T->Branch("truth_id",&b_id);
  T->Branch("truth_mother1",&b_mother1);
  T->Branch("truth_mother2",&b_mother2);
  T->Branch("truth_pt",&b_pt);
  T->Branch("truth_eta",&b_eta);
  T->Branch("truth_phi",&b_phi);
  T->Branch("truth_m",&b_m);

  T->Branch("truth_mupair_pt1",&b_muon1_pt);
  T->Branch("truth_mupair_eta1",&b_muon1_eta);
  T->Branch("truth_mupair_phi1",&b_muon1_phi);
  T->Branch("truth_mupair_ch1",&b_muon1_ch);
  T->Branch("truth_mupair_id1",&b_muon1_id);
  T->Branch("truth_mupair_bar1",&b_muon1_bar);
  T->Branch("truth_mupair_pt2",&b_muon2_pt);
  T->Branch("truth_mupair_eta2",&b_muon2_eta);
  T->Branch("truth_mupair_phi2",&b_muon2_phi);
  T->Branch("truth_mupair_ch2",&b_muon2_ch);
  T->Branch("truth_mupair_id2",&b_muon2_id);
  T->Branch("truth_mupair_bar2",&b_muon2_bar);

  
  // The event loop. A for loop without an incrementer because we will increment in the loop
  for (int iEvent = 0; iEvent < nEvent;){

    if (!pythia.next()) continue; // Generate next event

    // Restart muon & muon pair counting
    b_nmuons = 0;
    b_nmuonPairs = 0;

    // Empty out all of the vectors from the previous event
    b_barcode->clear();
    b_status->clear();
    b_id->clear();
    b_mother1->clear();
    b_mother2->clear();
    b_pt->clear();
    b_eta->clear();
    b_phi->clear();
    b_m->clear();

    b_muon_inds->clear();
    b_muon1_pt->clear();
    b_muon1_eta->clear();
    b_muon1_phi->clear();
    b_muon1_ch->clear();
    b_muon1_id->clear();
    b_muon1_bar->clear();
    b_muon2_pt->clear();
    b_muon2_eta->clear();
    b_muon2_phi->clear();
    b_muon2_ch->clear();
    b_muon2_id->clear();
    b_muon2_bar->clear();
    
    // Make sure all of our vectors have enough space to handle the outputs of the event.
    b_barcode->reserve(event.size());
    b_status->reserve(event.size());
    b_id->reserve(event.size());
    b_mother1->reserve(event.size());
    b_mother2->reserve(event.size());
    b_pt->reserve(event.size());
    b_eta->reserve(event.size());
    b_phi->reserve(event.size());
    b_m->reserve(event.size());
  
    // Set kinematic parameters from the event info
    b_weight=info.weight();
    b_QHard=info.QFac();

    for(int ip=1; ip<event.size(); ip++)
    {
      Pythia8::Particle& p=event[ip];
      if(p.status() > 0){
      	int aid=p.idAbs();
      	if (aid==13){
      	  if (std::abs(p.pT())>3.5){
      	    if (std::abs(p.eta())<2.5){
      	      b_nmuons++;
              b_muon_inds->push_back(ip);
      	    //std::cout << "muon found, pT:  "<< p.pT() << " eta: " << std::abs(p.eta()) << " b_nmuons: " << b_nmuons << std::endl;
      	    }
      	  }
      	}
      }
    }

    if (b_nmuons >= 2){
      iEvent++;
      //std::cout << "recording event" << std::endl;
      for(int ip=0; ip<event.size(); ip++){// Loop through all of the output particles in the event
    	  Pythia8::Particle& p=event[ip]; // Get the particle from the event
    	  // s > 0 => final state particle. 0 < s < 30 => Beam particles or intermediates for hardest subprocess
    	  //if(p.status() > 0 || std::abs(p.status()) < 30 )
    	  //{
    	  // Write all of the particle properties to the tree vectors
    	  b_barcode->push_back(ip);
    	  b_status->push_back(p.status());
    	  b_id->push_back(p.id());
    	  b_mother1->push_back(p.mother1());
    	  b_mother2->push_back(p.mother2());
    	  b_pt->push_back(p.pT());
    	  b_eta->push_back(p.eta());
    	  b_phi->push_back(p.phi());
    	  b_m->push_back(p.m());
    	}
      
      // Example cut! Re-start if there are no particles
      if(b_pt->size() == 0){
    	  iEvent--;
    	  continue;
    	}
      if(iEvent % 250 == 0){
	      printf("%i\n", iEvent);
      }

      // muon pairing
      for (int imu = 0; imu < b_nmuons - 1; imu++){
        for (int jmu = imu+1; jmu < b_nmuons; jmu++){
          b_muon1_pt->push_back(b_pt->at(b_muon_inds->at(imu)));
          b_muon1_eta->push_back(b_eta->at(b_muon_inds->at(imu)));
          b_muon1_phi->push_back(b_phi->at(b_muon_inds->at(imu)));
          b_muon1_ch->push_back((b_id->at(b_muon_inds->at(imu)) > 0)? -1 : 1); // +1 for mu+; -1 for mu-
          b_muon1_id->push_back(b_id->at(b_muon_inds->at(imu)));
          b_muon1_bar->push_back(b_muon_inds->at(imu));

          b_muon2_pt->push_back(b_pt->at(b_muon_inds->at(jmu)));
          b_muon2_eta->push_back(b_eta->at(b_muon_inds->at(jmu)));
          b_muon2_phi->push_back(b_phi->at(b_muon_inds->at(jmu)));
          b_muon2_ch->push_back((b_id->at(b_muon_inds->at(jmu)) > 0)? -1 : 1); // +1 for mu+; -1 for mu-
          b_muon2_id->push_back(b_id->at(b_muon_inds->at(jmu)));
          b_muon2_bar->push_back(b_muon_inds->at(jmu));
        }
      }

      T->Fill(); // Write this event to the root tree
    }else{
      continue;
    }
  }

  pythia.stat(); // Status results of the Pythia simulation
  settings.writeFile("settings.cfg");
  TEnv env("settings.cfg");

  // Now make another tree to include overall properties of this simulation
  TTree* meta_tree=new TTree("meta_tree","parameters");
  std::vector<int>* procCode=new std::vector<int>();
  std::vector<std::string>* procName=new std::vector<std::string>();
  std::vector<double>* sigma=new std::vector<double>();
  nTried=info.getCounter(4);  // Find the total number of events tried for jets
  efficiency=static_cast<double>(nEvent)/static_cast<double>(nTried); // Ratio of events w/ cuts to without  
  // Branch the tree, at first mostly with parameters of the simulation.
  meta_tree->Branch("nEvent",&nEvent);
  meta_tree->Branch("nTried",&nTried);
  // meta_tree->Branch("filterPtMin",&filterPtMin);
  // meta_tree->Branch("filterPtMax",&filterPtMax);
  meta_tree->Branch("pTHatMin",&pTHatMin);
  meta_tree->Branch("pTHatMax",&pTHatMax);
  meta_tree->Branch("totalSigma",&totalSigma); // Total cross-section, calculated below
  meta_tree->Branch("efficiency",&efficiency); // Efficiency, as defined above
  meta_tree->Branch("processCodes",&procCode); // All process codes
  meta_tree->Branch("processNames",&procName); // All process names
  meta_tree->Branch("processSigma",&sigma);    // All cross-sections, given by process

  // Make space for the cross-sections by process
  unsigned int nCodes=info.codesHard().size();
  procCode->reserve(nCodes);
  procName->reserve(nCodes);
  sigma->reserve(nCodes);

  // Save all of the process cross-sections
  for(unsigned int iproc=0; iproc<nCodes; iproc++){
    int jproc=info.codesHard()[iproc];
    procCode->push_back(jproc);
    procName->push_back(info.nameProc(jproc));
    double s1=info.sigmaGen(jproc);
    sigma->push_back(s1);
  }
  // Get total cross-section
  totalSigma=info.sigmaGen()*1e3 / info.weightSum();
  meta_tree->Fill(); // Write meta data

  T->Write(); // Write the event tree to the output file
  env.Write("pythia_params");
  meta_tree->Write(); // Write the meta data tree to the output file
  delete file; // Clean up the memory from the pointer
  return 0;
}
