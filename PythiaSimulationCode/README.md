# Pythia Simulation Private Sample Submission Script

## Event generation principles
* We generate nucleon-nucleon ("pp") collisions at CoM energy of 5.02GeV, with four beam-particle combinations: pp, np, pn, nn at ratio 4 : 6 : 6 : 9
  * This mimics the neutron to proton ratio 1.5:1 in Pb nucleus (A = 207, Z = 82)
* We select on dimuon events with a slightly lower muon-pT requirement (say, 3.7GeV) than the experimental minimum-muon-pT requirement of 4GeV
* Because the pThat of event is power-law suppressed, we generate events in five kinematic ranges:
  * pTHat in [5,10], [10,25], [25,60], [60,120], [120,3200] GeV
  * The high-pTHat ranges have very low cross sections (generate more events to minimize statistical errors), but much higher efficiency (% of events with two muons passing minimum muon-pT requirement) than low-pTHat ranges, so when requiring N events passing dimuon selection per job, need a much lower #trials --> job finishes faster
  	* We hence run larger #jobs and require much lower #events (pass 3.7GeV-pT dimuon requirement) per job for lower-pTHat ranges
  	  * As an example, each dimuon event for the lowest pTHat range (5-10GeV) takes ~1 hour --> only requires 10 events (passing filter) each job
  	  * #jobs and #target-events/job need to be set so that the total #events give us sufficient statistics
  	  * This is separated into two places: the #jobs for each kinmetic range and beam-particle combination is in the job-submission scripts; #target-events/job for each kinematic range is in the .cc file



## Directory structure & basic workflow
* The script RunHardQCD_dimuon.cc specifies all event/kinmetic settings for pythia scripts, and #target-events for each job at each kinematic range. It then generates trial events, filters on two muons with minimum pT requirement, writes all truth particles (including mothers and daughers), and dimuon properties to an output TTree, repeating the process until #events passing the dimuon selection reaches #target-events. It finishes the job by writing meta-data to meta-data TTree, and prints out summary (metadata) information.

* The bash scripts in submission_scripts/ 
  * takes the job name as input
  * specifies code + output directories & job types (nucleon-nucleon combinations: pp, pn, np, nn)
  * g++ compiles the latest version of the .cc script to get executable
  * calculate #jobs for each (job type, kinematic range) combination 
  * makes (kinematic range, job type)-specific logfile + outputfile sub-directories in the output area
  * then for each job, copies over the executable, writes condor-submission script + bash script (sets up environment + runs the executable) , and condor_submit the condor script