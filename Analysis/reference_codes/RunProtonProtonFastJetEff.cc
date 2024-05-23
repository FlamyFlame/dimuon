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
#include "TLorentzVector.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"


int main(int argc, char **argv)
{
  int nEvent = 100000; // The total number of events we will end up keeping
  int nTried = nEvent; // To be modified -- the total number of events generated without cuts
  int kinematicRange = 2; // Configure following variables depending on the intended kinematic range (in terms of jet pT intervals). 0: low, 1: medium, 2: high, 3: very high
  // Actual pT range of sim
  double filterPtMin, filterPtMax, pTHatMin, pTHatMax; // Defining the parameters which need to change with pT cut interval
  if (kinematicRange == 0) // Low pT jets -- filtering is harder here!
  {
    filterPtMin = 5;
    filterPtMax = 150;
    pTHatMin = 2.0;
    pTHatMax = 300.0;
  }
  else if (kinematicRange == 1) // Mid pT jets
  {
    filterPtMin = 10;
    filterPtMax = 150;
    pTHatMin = 7.0;
    pTHatMax = 300.0;
  }
  else if (kinematicRange == 2) // High pT jets
  {
    filterPtMin = 15;
    filterPtMax = 150;
    pTHatMin = 8.0;
    pTHatMax = 300.0;
  }
  else if (kinematicRange == 3) // High pT jets
  {
    filterPtMin = 40;
    filterPtMax = 150;
    pTHatMin = 10.0;
    pTHatMax = 320.0;
  }
  // else if (kinematicRange == 2) // High pT jets
  // {
  //   filterPtMin = 60;
  //   filterPtMax = 160;
  //   pTHatMin = 15.0;
  //   pTHatMax = 320.0;
  // }
  // else if (kinematicRange == 3) // Very high pT jets
  // {
  //   filterPtMin = 160;
  //   filterPtMax = 400;
  //   pTHatMin = 40.0;
  //   pTHatMax = 800.0;
  // }
  // minimum jet pT
  double jetCut = 5.0;
  // Definitions to be changed later
  double efficiency = 1;
  double totalSigma = 0;

  // Initialize Pythia and its RNG
  Pythia8::Pythia pythia;
  pythia.readString("Random:setSeed=on");

  // Set the random seed to the number supplied as an argument
  char buffer[50];
  // int largeSeed = 12345678 + atoi(argv[1]);
  int largeSeed = 45671234 + atoi(argv[1]);
  sprintf(buffer, "Random:seed=%d", largeSeed);
  // int largeSeed = 876543210 - atoi(argv[1]);
  // sprintf(buffer, "Random:seed=%d", largeSeed);
  // int largeSeed = atoi(argv[1]);
  // sprintf(buffer, "Random:seed=%s", argv[1]);
  pythia.readString(buffer);

  // Just some settings
  pythia.readString("Init:showChangedSettings     = off");
  pythia.readString("Init:showChangedParticleData = off");
  pythia.readString("Next:numberCount             = 0"  );
  pythia.readString("Next:numberShowInfo          = 0"  );
  pythia.readString("Next:numberShowProcess       = 0"  );
  pythia.readString("Next:numberShowEvent         = 0"  );

  // Setup objects to reference information about the simulation
  Pythia8::Settings& settings = pythia.settings;
  Pythia8::Info& info         = const_cast<Pythia8::Info&>(pythia.info);

  pythia.readString("Beams:eCM = 13000."); // CoM energy of the interaction (13 TeV)
  pythia.readString("Beams:idA = 2212"  ); // Monte carlo code proton
  pythia.readString("Beams:idB = 2212"  ); // Monte carlo code proton

  // Turn on all channels for processes we want
  pythia.readString("HardQCD:all               = on"); // Enable all hard QCD interactions (quark and gluon couplings)
  // pythia.readString("SoftQCD:inelastic         = on");
  // pythia.readString("SoftQCD:nonDiffractive    = on");
  // pythia.readString("SoftQCD:singleDiffractive = on");
  // pythia.readString("SoftQCD:doubleDiffractive = on");
  pythia.readString("SoftQCD:all               = off");
  pythia.readString("PartonLevel:MPI           = off");
  pythia.readString("PartonLevel:ISR           = on");
  // pythia.readString("PartonLevel:FSR           = on");

  settings.parm("PhaseSpace:pTHatMin", pTHatMin); // Apply kinematics constraint on hard processes
  settings.parm("PhaseSpace:pTHatMax", pTHatMax); // Apply kinematics constraint on hard processes

  // ------------------------------------------------------------------------
  // Using Monash 2013 now
  pythia.readString("Tune:ee = 7");
  pythia.readString("Tune:pp = 14");

  // Final-state radiation (FSR) parameters
  pythia.readString("TimeShower:alphaSvalue  = 0.1365");  // Effective alphaS(mZ) value
  pythia.readString("TimeShower:alphaSorder  = 1"     );  // Running order
  pythia.readString("TimeShower:alphaSuseCMW = false" );  // Translation from MS to CMW
  pythia.readString("TimeShower:pTmin        = 0.5"   );  // Cutoff for QCD radiation
  pythia.readString("TimeShower:pTminChgQ    = 0.5"   );  // Cutoff for QED radiation

  // String breaks: pT and z distributions
  pythia.readString("StringPT:sigma            = 0.335");  // Soft pT in string breaks (in GeV)
  pythia.readString("StringPT:enhancedFraction = 0.01" );  // Fraction of breakups with enhanced pT
  pythia.readString("StringPT:enhancedWidth    = 2.0"  );  // Enhancement factor
  pythia.readString("StringZ:aLund             = 0.68" );  // Lund FF a (hard fragmentation supp)
  pythia.readString("StringZ:bLund             = 0.98" );  // Lund FF b (soft fragmentation supp)
  pythia.readString("StringZ:aExtraSquark      = 0.0"  );  // Extra a when picking up an s quark
  pythia.readString("StringZ:aExtraDiquark     = 0.97" );  // Extra a when picking up a diquark
  pythia.readString("StringZ:rFactC            = 1.32" );  // Lund-Bowler c-quark parameter
  pythia.readString("StringZ:rFactB            = 0.855");  // Lund-Bowler b-quark parameter

  // Flavour composition: mesons
  pythia.readString("StringFlav:ProbStoUD     = 0.217");  // Strangeness-to-UD ratio
  pythia.readString("StringFlav:mesonUDvector = 0.5"  );  // Light-flavour vector suppression
  pythia.readString("StringFlav:mesonSvector  = 0.55" );  // Strange vector suppression
  pythia.readString("StringFlav:mesonCvector  = 0.88" );  // Charm vector suppression
  pythia.readString("StringFlav:mesonBvector  = 2.2"  );  // Bottom vector suppression
  pythia.readString("StringFlav:etaSup        = 0.60" );  // Suppression of eta mesons
  pythia.readString("StringFlav:etaPrimeSup   = 0.12" );  // Suppression of eta’ mesons

  // Flavour composition: baryons
  pythia.readString("StringFlav:probQQtoQ        = 0.081" );  // Diquark rate (for baryon production)
  pythia.readString("StringFlav:probSQtoQQ       = 0.915" );  // Strange-diquark suppression
  pythia.readString("StringFlav:probQQ1toQQ0     = 0.0275");  // Vector diquark suppression
  // pythia.readString("StringFlav:decupletSup      = 1.0"   );  // Spin-3/2 baryon suppression
  pythia.readString("StringFlav:suppressLeadingB = false" );  // Optional leading-baryon suppression
  pythia.readString("StringFlav:popcornSpair     = 0.9"   );  //
  pythia.readString("StringFlav:popcornSmeson    = 0.5"   );  //

  // PDF and matrix-element (ME) parameter
  pythia.readString("PDF:pSet                            = 13"   ); // PDF set for the proton
  pythia.readString("SigmaProcess:alphaSvalue            = 0.130"); // alphaS(MZ) for matrix elements
  pythia.readString("MultipartonInteractions:alphaSvalue = 0.130"); // alphaS(MZ) for MPI

  // Initial-state radiation (ISR) and primordial-kT parameter
  // All parameters for this section are the same as paper

  // Multi-Parton-Interaction (MPI), color-reconnection (CR), and diffractive parameters
  // All parameters for this section are the same as paper
  // pythia.readString("BeamRemnants:reconnectRange    = 1.8"   );  // Color reconnection
  // ------------------------------------------------------------------------


  // ------------------------------------------------------------------------
  // ussing A14
  // pythia.readString("Tune:ee = 7");
  // pythia.readString("Tune:pp = 19");
  // pythia.readString("PDF:pSet = 8"); // PDF set for the proton    CTEQ6L1, LO alpha_s(M_Z) = 0.1298
  // ------------------------------------------------------------------------



  settings.listChanged(); // Print all settings modifications


  //settings.listAll();
  if (!pythia.init()) // Initialize pythia with all of these settings
  {
    std::cerr << "Failed to initialize" << std::endl;
    exit(0);
  }


  Pythia8::Event& event = pythia.event; // Make a reference to the event we are generating
  // Set up an anti-kT jet clustering algorithm with R=0.4 and energy recombination scheme
  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, 0.4, fastjet::E_scheme, fastjet::Best);

  TFile *file = TFile::Open("pytreeEff.root", "recreate"); // Place to write out
  TTree *T = new TTree("PyTree", "ev1 Tree"); // Tree we will write data into
  // std::vector<int>*    b_barcode = new std::vector<int>();
  // std::vector<int>*    b_status  = new std::vector<int>();
  std::vector<int>*    b_id      = new std::vector<int>(); // Particle identity: p, n, etc
  std::vector<double>* b_pt      = new std::vector<double>(); // Kinematics
  std::vector<double>* b_eta     = new std::vector<double>(); // Kinematics
  std::vector<double>* b_phi     = new std::vector<double>(); // Kinematics
  std::vector<double>* b_m       = new std::vector<double>(); // Kinematics
  std::vector<double>* b_tau     = new std::vector<double>(); // Kinematics
  std::vector<double>* b_jet_pt  = new std::vector<double>(); // Jet Kinematics
  std::vector<double>* b_jet_eta = new std::vector<double>(); // Jet Kinematics
  std::vector<double>* b_jet_phi = new std::vector<double>(); // Jet Kinematics

  double b_weight = 0.; // Event weighting, usually 1
  double b_QHard  = 0.; // Momentum transfer
  // Branch the tree for all simulation parameters we want to store
  T->Branch("event_weight"    , &b_weight );
  T->Branch("QHard"           , &b_QHard  );
  // T->Branch("particle_barcode", &b_barcode);
  // T->Branch("particle_status" , &b_status );
  T->Branch("particle_id"     , &b_id     );
  T->Branch("particle_pt"     , &b_pt     );
  T->Branch("particle_eta"    , &b_eta    );
  T->Branch("particle_phi"    , &b_phi    );
  T->Branch("particle_m"      , &b_m      );
  T->Branch("particle_tau"    , &b_tau    );
  T->Branch("jet_pt"          , &b_jet_pt );
  T->Branch("jet_eta"         , &b_jet_eta);
  T->Branch("jet_phi"         , &b_jet_phi);

  unsigned int totalEvent = 0;
  // The event loop. A for loop without an incrementer because we will increment in the loop
  for (int iEvent = 0; iEvent < nEvent;)
  {
    if (!pythia.next()) continue; // Generate next event
    totalEvent++;
    if (totalEvent > 1) {
      float log10Entry = floor(log (totalEvent) / log(10));
      int   modEntry   = totalEvent % (int) pow(10, log10Entry);
      if (modEntry == 0) printf("                               PYTHIA8 event %i\n", totalEvent);
    }

    // Empty out all of the vectors from the previous event
    // b_barcode->clear();
    // b_status->clear();
    b_id    ->clear();
    b_pt    ->clear();
    b_eta   ->clear();
    b_phi   ->clear();
    b_m     ->clear();
    b_tau   ->clear();

    std::vector <fastjet::PseudoJet> fjInputs; // Create a pesudojet to run our clustering algorithm

    // Make sure all of our vectors have enough space to handle the outputs of the event.
    fjInputs.reserve(event.size());
    // b_barcode->reserve(event.size());
    // b_status->reserve(event.size());
    b_id    ->reserve(event.size());
    b_pt    ->reserve(event.size());
    b_eta   ->reserve(event.size());
    b_phi   ->reserve(event.size());
    b_m     ->reserve(event.size());
    b_tau   ->reserve(event.size());

    // Set kinematic parameters from the event info
    b_weight = info.weight();
    b_QHard = info.QFac();

    for (int ip = 0; ip < event.size(); ip++) {// Loop through all of the output particles in the event

      Pythia8::Particle& p = event[ip]; // Get the particle from the event
      // s > 0 => final state particle. 0 < s < 30 => Beam particles or intermediates for hardest subprocess
      // if (p.status() > 0 || std::abs(p.status()) < 30 )

      if (std::abs(p.eta()) > 4.9) {continue;}
      // if (p.tau() * 1E9 / 3E8 < 300) {continue;}

      int aid = p.idAbs();
      // if (p.isFinal() && (aid == 211 || aid == 321 || aid == 2212)) {
      if (p.isFinal() &&
          (aid == 211 || aid == 111 ||
           aid == 321 || aid == 311 || aid == 310 || aid == 130 ||
           aid == 2212 || aid == 2112 ||
           aid == 3122 || aid == 3222 || aid == 3212 || aid == 3112 ||
           aid == 11 || aid == 13 || aid == 22)) {
        // Write all of the particle properties to the tree vectors
        // b_barcode->push_back(ip        );
        // b_status ->push_back(p.status());
        b_id ->push_back(p.id() );
        b_pt ->push_back(p.pT() );
        b_eta->push_back(p.eta());
        b_phi->push_back(p.phi());
        b_m  ->push_back(p.m()  );
        b_tau->push_back(p.tau());

        TLorentzVector v;
        v.SetPtEtaPhiM(p.pT(), p.eta(), p.phi(), 0);

        // If the particle is not an electron/muon/tau neutrino (12/14/16) or a muon (13)
        // if (aid != 12 && aid != 13 && aid != 14 && aid != 16)
        // fjInputs.push_back( fastjet::PseudoJet( p.px(), p.py(), p.pz(), p.e() ) ); // Add to jet clustering
        fjInputs.push_back( fastjet::PseudoJet(v.Px(), v.Py(), v.Pz(), v.E()) ); // Add to jet clustering
      }
    }

    fastjet::ClusterSequence clustSeq(fjInputs, jetDef); // Setup jet inputs and clustering algorithm
    std::vector<fastjet::PseudoJet> sortedJets = fastjet::sorted_by_pt(clustSeq.inclusive_jets(jetCut)); // Jets!
    double pTlead = 0;


    // --------------------------------
    // only save jet events  -- default
    if (sortedJets.size() > 1) { // If there are any jets, get the leading pT
      pTlead = sortedJets[0].perp();
      if (pTlead < filterPtMin) {continue;}

      bool goodJet       = false;
      // bool haveJetInFcal = false;
      for (unsigned int ij = 0; ij < sortedJets.size(); ij++) {
        double pT  = sortedJets[ij].pt();
        double eta = sortedJets[ij].eta();
        // if (std::abs(eta) > 2.5) {haveJetInFcal = true;}
        if (pT > filterPtMin && std::abs(eta) < 2.1) {
          for (unsigned int jj = 0; jj < sortedJets.size(); jj++) {
            if (ij == jj) {continue;}
            if (sortedJets[jj].pt() < 10.) {continue;}
            if (std::abs(sortedJets[jj].eta()) > 4.5) {continue;}

            double dphi  = sortedJets[ij].phi_std() - sortedJets[jj].phi_std();
            dphi = std::atan2(std::sin(dphi), std::cos(dphi));
            if (dphi > 5 * TMath::Pi() / 6 || dphi < -5 * TMath::Pi() / 6) {goodJet = true;}
          }
        }
      }

      // if (goodJet && !haveJetInFcal) {iEvent++;}
      if (goodJet) {iEvent++;}
      else {continue;}
    }
    else {continue;}
    // --------------------------------

    // iEvent++;
    if (iEvent > 1) {
      // float log10Entry = floor(log (iEvent) / log(10));
      // int   modEntry   = iEvent % (int) pow(10, log10Entry);
      // if (modEntry == 0) printf("Processed (jet) event %i\n", iEvent);
      if (iEvent % 1000 == 0) printf("Processed (jet) event %i\n", iEvent);
    }

    if (sortedJets.size() != 0) {
      // Now, clear out the jet vectors
      b_jet_pt ->clear();
      b_jet_eta->clear();
      b_jet_phi->clear();

      // Make space for all the jets we found
      b_jet_pt ->reserve(sortedJets.size());
      b_jet_eta->reserve(sortedJets.size());
      b_jet_phi->reserve(sortedJets.size());

      // Loop through and add all the jets to the tree by kinematics
      for (unsigned int ijet = 0; ijet < sortedJets.size(); ijet++) {
        const fastjet::PseudoJet& p = sortedJets.at(ijet);
        if (p.pt() < 10.) continue;
        if (std::abs(p.eta()) > 4.5) {continue;}
        b_jet_pt ->push_back(p.pt     ());
        b_jet_eta->push_back(p.eta    ());
        b_jet_phi->push_back(p.phi_std());
      }
    }
    T->Fill(); // Write this event to the root tree
  }

  pythia.stat(); // Status results of the Pythia simulation
  settings.writeFile("settings.cfg");
  TEnv env("settings.cfg");

  // Now make another tree to include overall properties of this simulation
  TTree* meta_tree = new TTree("meta_tree", "parameters");
  std::vector<int>* procCode = new std::vector<int>();
  std::vector<std::string>* procName = new std::vector<std::string>();
  std::vector<double>* sigma = new std::vector<double>();
  nTried = info.getCounter(4); // Find the total number of events tried for jets
  efficiency = static_cast<double>(nEvent) / static_cast<double>(nTried); // Ratio of events w/ cuts to without
  // Branch the tree, at first mostly with parameters of the simulation.
  meta_tree->Branch("nEvent"     , &nEvent     );
  meta_tree->Branch("seed"       , &largeSeed  );
  meta_tree->Branch("nTried"     , &nTried     );
  meta_tree->Branch("filterPtMin", &filterPtMin);
  meta_tree->Branch("filterPtMax", &filterPtMax);
  meta_tree->Branch("pTHatMin"   , &pTHatMin   );
  meta_tree->Branch("pTHatMax"   , &pTHatMax   );
  meta_tree->Branch("totalSigma" , &totalSigma ); // Total cross-section, calculated below
  meta_tree->Branch("efficiency" , &efficiency ); // Efficiency, as defined above
  // meta_tree->Branch("processCodes", &procCode); // All process codes
  // meta_tree->Branch("processNames", &procName); // All process names
  meta_tree->Branch("processSigma", &sigma);   // All cross-sections, given by process

  // Make space for the cross-sections by process
  unsigned int nCodes = info.codesHard().size();
  procCode->reserve(nCodes);
  procName->reserve(nCodes);
  sigma->reserve(nCodes);

  // Save all of the process cross-sections
  for (unsigned int iproc = 0; iproc < nCodes; iproc++)
  {
    int jproc = info.codesHard()[iproc];
    procCode->push_back(jproc);
    procName->push_back(info.nameProc(jproc));
    double s1 = info.sigmaGen(jproc);
    sigma->push_back(s1);
  }
  // Get total cross-section
  totalSigma = info.sigmaGen() * 1e3 / info.weightSum();
  // totalSigma = info.sigmaGen() / info.weightSum();
  meta_tree->Fill(); // Write meta data

  T->Write(); // Write the event tree to the output file
  // env.Write("pythia_params");
  meta_tree->Write(); // Write the meta data tree to the output file
  delete file; // Clean up the memory from the pointer
  return 0;
}
