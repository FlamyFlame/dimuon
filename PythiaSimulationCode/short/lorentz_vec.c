PyTree->Scan("truth_id : truth_mother1 : truth_daughter1 : truth_status : truth_px : truth_py : truth_pz : truth_e : truth_m","","",1,0)

PyTree->Scan("truth_id : truth_mother1 : truth_daughter1 : truth_status : truth_pt : truth_eta : truth_phi : truth_m","","",1,5)


TLorentzVector Mg, Mc, Mcbar;
Mg.SetPxPyPzE(7.993e-14, -1.82e-13, -194.5147, 194.51479)
Mc.SetPxPyPzE(-0.417258, 0.9552411, -39.11787, 39.160500)

SetPxPyPzE()

TLorentzVector Mg1, Mb1, Mbbar1;
Mbbar1.SetPtEtaPhiM(54.584226, -2.121826, 1.9364878, 4.8)
Mb1.SetPtEtaPhiM(45.731632, -0.755733, 0.9558686, 4.8)
Mg1 = Mb1 + Mbbar1

Mbbar1.SetPtEtaPhiM(66.913547, 1.8725390, -0.088124,       4.8)
Mb1.SetPtEtaPhiM(66.913547, 0.3325877, 3.0534684,       4.8)



// what follows is an attempt at Q^2 determination for the g --> b bbar splitting
// in pytree_trial1.root, event 0
// I tried both SetPtEtaPhiM and SetPxPyPzE to see if the result is similar
// this is an important sanity/consistency check since I'm not sure what the recording of the mass
// as being the on-shell mass (while the particles in an ISR/FSR are virtual) has on the Q^2 calculation

TLorentzVector Mg, Mb, Mbbar;
Mg.SetPtEtaPhiM(4.440e-16, 41.147049, -1.570796,         0)
Mb.SetPtEtaPhiM(16.514535, 0.9427780, 3.0208023,       4.8)
Mbbar = Mg - Mb
Mbbar.M() // -47.411811
Mbbar.Mt() // -44.442659

// below is a way to get kT^2 from Q^2, for either initial or final state radiation
TVector3 vg = Mg.Vect()
TVector3 vb = Mb.Vect()
TVector3 vbbar = Mbbar.Vect()
double z = vbbar.Dot(vg) / vg.Dot(vg)
double kT2 = -z * Mbbar.M2() // for final state radiation, kT^2 = z * (1-z) * Mbbar.M2()


TLorentzVector Mg2, Mb2, Mbbar2;
Mg2.SetPxPyPzE(0, -4.44e-16, 164.57990, 164.57990)
Mb2.SetPxPyPzE(-16.39420, 1.9899486, 17.980686, 24.881217)
Mbbar2 = Mg2 - Mb2
Mbbar2.M() // -47.416634
Mbbar2.Mt() // -44.447807



// the following is to test out what we get when we pretend we are ignorant of whether the gluon splitting
// happens as the primordial kT splitting right before the hard scatter
// or a hard FSR splitting right after
// so that in powheg this would be g X --> b bbar X
// and one thing that is sensible to do is the compare the minv of the b bbar pair
// with the energy scale of the hard scattering
// for instance by taking the ratio m_bbbar / m_gX (mHat in pythia)
// or m_bbbar / mT_gX (pTHat in pythia)

TLorentzVector Mg3, Mb3, Mbbar3;
Mbbar3.SetPxPyPzE(70.997387, 18.087201, 20.770439, 76.303519)
Mb3.SetPxPyPzE(-16.39420, 1.9899486, 17.980686, 24.881217)
Mg3 = Mb3 + Mbbar3
Mg3.M() // 73.158056
Mg3.Mt() // 93.470322



SetPtEtaPhiM()

M1.SetPtEtaPhiM(25.244, 0.431, 1.843);
M2.SetPtEtaPhiM(15.355, 2.113, -1.146);



TLorentzVector Mg1, Mg2, Ms, Msbar, Mb, Mbbar;
Mg1.SetPtEtaPhiM(7.324e-15, 37.77, -1.326, 0)
Mg2.SetPtEtaPhiM(2.232e-14, -36.73, 0.102, 0)
Mb.SetPtEtaPhiM(15.355, 2.11, -1.146, 4.18)
// Mb.SetPtEtaPhiM(13.591, 2.24, -1.18, 4.18)
Mbbar.SetPtEtaPhiM(36.55, 0.164, 1.89, 4.18)
Msbar.SetPtEtaPhiM(21.345, -0.35, -1.326, 0.095)
Ms.SetPtEtaPhiM(2.052, -4.170, -3.042, 0.095)
// Ms.SetPtEtaPhiM(2.491, -3.976, -2.463, 0.095)




