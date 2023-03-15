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
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

// We generate events and record kinematics of general particles, muons and jets
// and event char's such as xGamma, yGamma, xA
// We require muons to satisfy the experimental constraints: abs(eta) < 2.7 and pT > 3 

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
  int nEvent = 1000; // The total number of events we will end up keeping
  int nTried =nEvent; // To be modified -- the total number of events generated without cuts
  int kinematicRange = 2; // Configure following variables depending on the intended kinematic range (in terms of jet pT intervals). 0: low, 1: medium, 2: high
  // Actual pT range of sim
  double filterPtMin, filterPtMax, pTHatMin, pTHatMax, Wmin; // Defining the parameters which need to change with pT cut interval
  if(kinematicRange == 0) // Low pT muons -- filtering is harder here!
    { // no pTHat filter due to conflict with soft processes (momentum transf not clearly defn) 
      pTHatMin=3.0;
      pTHatMax=8.0;
      Wmin = 5.0;
  }
  else if(kinematicRange == 1) // Mid-low pT muons
  {
  	pTHatMin=8.0;
  	pTHatMax=25.0;
  	Wmin = 10.0;
  }
  else if(kinematicRange == 2) // Mid pT muons
  {
  	pTHatMin=25.0;
  	pTHatMax=40.0;
  	Wmin = 10.0;
  }
  
  else if(kinematicRange == 3) // High pT jets
  {
  	pTHatMin=40.0;
  	pTHatMax=320.0;
  	Wmin = 10.0;
  }
  
  // minimum jet pT
  double jetCut = 15.0;
  // Center of mass energy range
  double Wmax = 1142; 
  // Definitions to be changed later
  double efficiency=1;
  double totalSigma=0;

  // Initialize Pythia and its RNG
  Pythia8::Pythia pythia;
  pythia.readString("Random:setSeed=on");

  // Set the random seed to the number supplied as an argument
  char buffer[50];
  sprintf(buffer,"Random:seed=%s",argv[1]);
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
  Pythia8::Info& info         = pythia.info;

  
  //A14 tune, base tune is Monash 2013
  pythia.readString("Tune:ee = 7");
  pythia.readString("Tune:pp = 14");
  pythia.readString("SpaceShower:rapidityOrder = on");
  pythia.readString("SigmaProcess:alphaSvalue = 0.140"); //coupling strength of strong intr
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
  pythia.readString("Beams:idA = 13"); // Monte carlo code muon
  if(atoi(argv[2]) == 0)  // Runtime option: proton or neutron simulations
    pythia.readString("Beams:idB = 2212"); // Monte carlo code proton
  else
    pythia.readString("Beams:idB = 2112"); // Monte carlo code neutron
  //pythia.readString("PDF:pSet = LHAPDF6:nCTEQ15npFullNuc_208_82/0001"); // Set PDFs
  // Enable photon flux sampling from a PDF. Currently set to muon EPA but will be re-pointed to custom nuclear
  pythia.readString("PDF:lepton2gamma = on");

  // Use the nucleus distribution we defined
  // Following line removes photon virtuality test. Required for custom PDF. Small error since Q2 is small
  pythia.readString("Photon:sampleQ2 = off"); 
  Nucleus2gamma* photonFlux = new Nucleus2gamma(13);
  pythia.setPhotonFluxPtr(photonFlux, 0);

  // Apply the above constraints
  pythia.readString("Photon:Q2max = 1.0"); // Upper limit on photon virtuality
  if(kinematicRange == 0)
  	pythia.readString("Photon:Wmin  = 5.0"); // Set lower limit for the invariant mass of 2-photon systems
  else
  	pythia.readString("Photon:Wmin  = 10.0"); // Set lower limit for the invariant mass of 2-photon systems
  pythia.readString("Photon:Wmax = 1142."); // Set upper limit for the invariant mass of 2-photon systems
  if (kinematicRange != 0){
    settings.parm("PhaseSpace:pTHatMin", pTHatMin); // Apply kinematics constraint on hard processes if not low pT
    settings.parm("PhaseSpace:pTHatMax", pTHatMax); // Apply kinematics constraint on hard processes if not low pT
  }
  
  // Turn on all channels for processes we want
  //pythia.readString("HardQCD:all = on"); // Enable all hard QCD interactions (quark and gluon couplings)
  pythia.readString("HardQCD:gg2ccbar = on");
  pythia.readString("HardQCD:qqbar2ccbar = on");
  pythia.readString("HardQCD:gg2bbbar = on");
  pythia.readString("HardQCD:qqbar2bbbar = on");

  //if(kinematicRange == 0)
  //pythia.readString("SoftQCD:nonDiffractive = on"); // Enable soft QCD processes which are not diffractive
  pythia.readString("PhotonParton:all = on"); // Enable relevant interactions between photons and partons
  settings.mode("Photon:ProcessType", 0); // Enables all 4 types of photon-photon interactions

  settings.listChanged(); // Print all settings modifications


  //settings.listAll();
  if(!pythia.init()) // Initialize pythia with all of these settings
  {
    std::cerr << "Failed to initialize" << std::endl;
    exit(0);
  }

  
  Pythia8::Event& event=pythia.event; // Make a reference to the event we are generating
  // event: a wrapper for a vector of the particles to be generated
  // Set up an anti-kT jet clustering algorithm with R=0.4 and energy recombination scheme  
  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, 0.4, fastjet::E_scheme, fastjet::Best);

  TFile *file = TFile::Open("pytreeEff.root","recreate"); // Place to write out
  TTree *T = new TTree("PyTree","ev1 Tree"); // Tree we will write data into
  std::vector<int>* b_barcode=new std::vector<int>();
  std::vector<int>* b_status=new std::vector<int>();
  std::vector<int>* b_id=new std::vector<int>(); // Particle identity: p, n, etc
  std::vector<double>* b_pt=new std::vector<double>(); // Kinematics
  std::vector<double>* b_eta=new std::vector<double>(); // Kinematics
  std::vector<double>* b_phi=new std::vector<double>(); // Kinematics
  std::vector<double>* b_m=new std::vector<double>(); // Kinematics
  std::vector<double>* b_jet_pt=new std::vector<double>(); // Jet Kinematics
  std::vector<double>* b_jet_eta=new std::vector<double>(); // Jet Kinematics
  std::vector<double>* b_jet_phi=new std::vector<double>(); // Jet Kinematics
  std::vector<double>* b_jet_m=new std::vector<double>(); // Jet Kinematics
  std::vector<int>* b_mother1=new std::vector<int>(); // Index of parent1
  std::vector<int>* b_mother2=new std::vector<int>(); // Index of parent2
  std::vector<double>* b_muon_pt=new std::vector<double>(); // Muon Kinematics
  std::vector<double>* b_muon_eta=new std::vector<double>(); // Muon Kinematics 
  std::vector<double>* b_muon_phi=new std::vector<double>(); // Muon Kinematics
  std::vector<double>* b_muon_m=new std::vector<double>(); // Muon Kinematics
  std::vector<int>* b_muon_mother1_id=new std::vector<int>(); // First mother
  //std::vector<int>* b_muon_mother2_id=new std::vector<int>(); // Last mother
  //std::vector<std::vector<int>>* b_muon_mother1_list=new std::vector<std::vector<int>>();
  std::vector<int>* b_muon_earliest_mother_id=new std::vector<int>();
  
  double b_weight=0.; // Event weighting, usually 0
  double b_yGamma=0.; // Fraction of photon energy that goes into the hard scattering
  double b_xGamma=0.; // Fraction of nucleus energy carried by the photon
  double b_xA=0.; // Fraction of nucleus energy carried by the parton
  double b_QGamma=0; // Momentum transfer
  double b_QHard=0; // Momentum transfer
  int b_muonCount=0; // Number of muons in the current event
  bool b_isNonDiffractive = false;
  
  // Branch the tree for all simulation parameters we want to store
  T->Branch("event_weight",&b_weight);
  T->Branch("yGamma",&b_yGamma);
  T->Branch("xGamma",&b_xGamma);
  T->Branch("xA",&b_xA);
  T->Branch("QGamma",&b_QGamma);
  T->Branch("QHard",&b_QHard);
  T->Branch("muonCount",&b_muonCount);
  T->Branch("isNonDiffractive",&b_isNonDiffractive);
  T->Branch("particle_barcode",&b_barcode);
  T->Branch("particle_status",&b_status);
  T->Branch("particle_id",&b_id);
  T->Branch("particle_pt",&b_pt);
  T->Branch("particle_eta",&b_eta);
  T->Branch("particle_phi",&b_phi);
  T->Branch("particle_m",&b_m);
  T->Branch("mother1_ind",&b_mother1);
  T->Branch("mother2_ind",&b_mother2);
  T->Branch("jet_pt",&b_jet_pt);
  T->Branch("jet_eta",&b_jet_eta);
  T->Branch("jet_phi",&b_jet_phi);
  T->Branch("jet_m",&b_jet_m);
  T->Branch("muon_pt",&b_muon_pt);
  T->Branch("muon_eta",&b_muon_eta);
  T->Branch("muon_phi",&b_muon_phi);
  T->Branch("muon_m",&b_muon_m);
  T->Branch("muon_mother1_id",&b_muon_mother1_id);
  //T->Branch("muon_mother2_id",&b_muon_mother2_id);
  //T->Branch("muon_mother1_list",&b_muon_mother1_list);
  T->Branch("muon_earliest_mother_id",&b_muon_earliest_mother_id);
  
  // The event loop. A for loop without an incrementer because we will increment in the loop
  for (int iEvent = 0; iEvent < nEvent;)
  {

    if (!pythia.next()) continue; // Generate next event

    // Empty out all of the vectors from the previous event
    b_barcode->clear();
    b_status->clear();
    b_id->clear();
    b_pt->clear();
    b_eta->clear();
    b_phi->clear();
    b_m->clear();
    b_mother1->clear();
    b_mother2->clear();
    b_muon_pt->clear();
    b_muon_eta->clear();
    b_muon_phi->clear();
    b_muon_m->clear();
    b_muon_mother1_id->clear();
    //b_muon_mother2_id->clear();
    //b_muon_mother1_list->clear();
    b_muon_earliest_mother_id->clear();
    b_muonCount	= 0;  // restart on muon counting 
    
    std::vector <fastjet::PseudoJet> fjInputs; // Create a pesudojet to run our clustering algorithm
    //std::vector<Pythia8::Particle>* muonArr=new std::vector<Pythia8::Particle>();
    //std::vector<double>* muon_pt=new std::vector<double>();
    
    // Make sure all of our vectors have enough space to handle the outputs of the event.
    fjInputs.reserve(event.size()); // event.size() ~ # of particles generated in the event
    b_barcode->reserve(event.size());
    b_status->reserve(event.size());
    b_id->reserve(event.size());
    b_pt->reserve(event.size());
    b_eta->reserve(event.size());
    b_phi->reserve(event.size());
    b_m->reserve(event.size());
    b_mother1->reserve(event.size());
    b_mother2->reserve(event.size());
    b_muon_pt->reserve(event.size());
    b_muon_eta->reserve(event.size());
    b_muon_phi->reserve(event.size());
    b_muon_m->reserve(event.size());
    b_muon_mother1_id->reserve(event.size());
    //b_muon_mother2_id->reserve(event.size());
    //b_muon_mother1_list->reserve(event.size());
    b_muon_earliest_mother_id->reserve(event.size());
    
    // Set kinematic parameters from the event info
    b_weight=info.weight();
    b_yGamma=info.xGammaA();
    b_xGamma=info.x1pdf()/info.xGammaA();
    b_xA=info.x2pdf();
    b_QGamma=std::sqrt(info.Q2GammaA());
    b_QHard=info.QFac();
    b_isNonDiffractive=info.isNonDiffractive();
    
    // No longer able to detect the event well for higher energy photons
    if (b_yGamma > 0.05) continue;

    int ind = 0; // index of the current parent in the c/c hadron parental chain
    int cid = 0;
    std::vector<int> curr_event_mother1_list; // c/c hadron parental chain of one muon
    
    for(int ip=1; ip<event.size(); ip++)// Loop through all of the output particles in the event
    {
      Pythia8::Particle& p=event[ip]; // Get the particle from the event
      // s > 0 => final state particle. 0 < s < 30 => Beam particles or intermediates for hardest subprocess
      // s < 0 -> has already decayed into other particles
      // s = 11-19: beam particles; s = 20-29: pts of the hardest subprocess

      //if(p.status() > 0 || std::abs(p.status()) < 30 )
      if(1)
      {
	// Write all of the particle properties to the tree vectors
        b_barcode->push_back(ip);
        b_status->push_back(p.status());
        b_id->push_back(p.id());
        b_pt->push_back(p.pT());
        b_eta->push_back(p.eta());
        b_phi->push_back(p.phi());
        b_m->push_back(p.m());
	      b_mother1->push_back(p.mother1());
	      b_mother2->push_back(p.mother2());
        if(p.status() > 0) // If it's a final state particle
        {
          int aid=p.idAbs();
	        // If the particle is not an electron/muon/tau neutrino (12/14/16) or a muon (13)
	        // Also cut extremely forward particles from the jets
          //if(aid!=12 && aid!=13 && aid!=14 && aid!=16 && (std::abs(p.eta()) < 4.9))  
      	  if(aid!=12 && aid!=13 && aid!=14 && aid!=16)
      	    fjInputs.push_back( fastjet::PseudoJet( p.px(), p.py(), p.pz(), p.e() ) ); // Add to jet clustering
      	  else if(aid==13 && abs(p.eta()) < 2.7 && p.pT() > 4){ // Muon that satisfies |eta| < 7.2 & pT > 4 GeV -> experimental constraints & get rid of fake beam remnant
      	    b_muon_pt->push_back(p.pT());
      	    b_muon_eta->push_back(p.eta());
      	    b_muon_phi->push_back(p.phi());
      	    b_muon_m->push_back(p.m());
            b_muon_mother1_id->push_back(event[p.mother1()].id() % 10000);
            //b_muon_mother2_id->push_back(event[p.mother2()].id() % 10000);
      	    ind = p.mother1(); // initiaze the index to be that of the muon's immediate parent
      	    cid = abs(event[ind].id()) % 10000;

      	    if (cid == 15){ // if immed parent is tau, trace back one step and record its immed parent instead
      	      ind = event[ind].mother1();
      	      cid = abs(event[ind].id()) % 10000;
      	    }

      	    while (cid == 13){ // if immed parent is muon (bremstralung), trace back until no longer is muon
      	      ind = event[ind].mother1();
      	      cid = abs(event[ind].id()) % 10000;
      	    }
      	    
      	    curr_event_mother1_list.clear(); // clear the vector that records the c/c hadron parental history of the muon

      	    while ((cid >= 400 && cid < 500) || (cid >= 4000 && cid < 5000)){ // if the current parent is a c/c hadron, then trace one step back
      	      curr_event_mother1_list.push_back(event[ind].id() % 10000); //note : youngest parent appears first in the c/c hadron mother list
      	      ind = event[ind].mother1(); // update ind to be index of the parent particle of the current particle
      	      cid = abs(event[ind].id()) % 10000; // |current id| % 10,000
      	    }
      	    // result: trace back until the parent particle is not a c hadron
      	    // so far: the list contains ONLY ancestors that are c hadrons

      	    curr_event_mother1_list.push_back(event[ind].id()); // record the youngest ancestor that is NOT a c hadron
      	    // Note: if muons's immediate parent is not a c hadron, only the immediate parent is recorded
      	    b_muon_earliest_mother_id->push_back(curr_event_mother1_list.at(curr_event_mother1_list.size()-1)); // the youngest ancestor of the muon that is not a c hadron
      	    //b_muon_mother1_list->push_back(curr_event_mother1_list);
      	    b_muonCount += 1;
      	  }
      	}
      }
    }

    fastjet::ClusterSequence clustSeq(fjInputs, jetDef); // Setup jet inputs and clustering algorithm
    std::vector<fastjet::PseudoJet> sortedJets=fastjet::sorted_by_pt(clustSeq.inclusive_jets(jetCut)); // Jets!
    double pTlead=0;
    //if(sortedJets.size() > 0) 
    //pTlead=sortedJets[0].perp(); // If there are any jets, get the leading pT
    if(b_muonCount >= 2)
      iEvent++;
    else
      continue;
    
    // Check to make sure that the groomer didn't remove all the jets
    /*    if(sortedJets.size() == 0 || muon_pt.size() != 0)
    {
      iEvent--;
      continue;
    }
    */

    if(iEvent % 50 == 0)
      printf("%i\n", iEvent);

    // Now, clear out the jet vectors
    b_jet_pt->clear();
    b_jet_eta->clear();
    b_jet_phi->clear();
    b_jet_m->clear();

    // Make space for all the jets we found
    b_jet_pt->reserve(sortedJets.size());
    b_jet_eta->reserve(sortedJets.size());
    b_jet_phi->reserve(sortedJets.size());
    b_jet_m->reserve(sortedJets.size());

    // Loop through and add all the jets to the tree by kinematics and partonic/hadronic kinematics
    double pzSum = 0.0, pxSum = 0.0, pySum = 0.0, ESum = 0.0;
    for(unsigned int ijet=0; ijet<sortedJets.size(); ijet++)
    {
      const fastjet::PseudoJet& p=sortedJets.at(ijet);
      b_jet_pt->push_back(p.perp());
      b_jet_eta->push_back(p.pseudorapidity());
      b_jet_phi->push_back(p.phi_std());
      b_jet_m->push_back(p.m());
      pzSum += p.pz();
      pxSum += p.px();
      pySum += p.py();
      ESum += p.E();
    }
    T->Fill(); // Write this event to the root tree
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
  meta_tree->Branch("filterPtMin",&filterPtMin);
  meta_tree->Branch("filterPtMax",&filterPtMax);
  if (kinematicRange != 0){
    meta_tree->Branch("pTHatMin",&pTHatMin);
    meta_tree->Branch("pTHatMax",&pTHatMax);
  }
  meta_tree->Branch("Wmin",&Wmin);
  meta_tree->Branch("Wmax",&Wmax);
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
  for(unsigned int iproc=0; iproc<nCodes; iproc++)
  {
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
