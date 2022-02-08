// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
//#include <Vector4.hh> //For sorting by pT: http://rivet.hepforge.org/code/2.2.0/a00121.html#details
#include <iostream>

#define MZ 91.188 

using namespace std;

namespace Rivet {

class WZ_ATLAS : public Analysis {

  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(WZ_ATLAS);

		float totalEvents=0;
		float vetoedEvents=0;
		float selectedEvents=0;	
		float jetsEvents=0;

    /// Book histograms and initialise projections before the run
    void init() {

			const FinalState fs(Cuts::abseta < 100.);
			declare(fs , "fs");
			
			const VisibleFinalState vfs(Cuts::abseta < 5.);
			declare(vfs , "vfs");

      FastJets jetfs(vfs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs, "jets");

			IdentifiedFinalState leadingLeptons(fs);
			leadingLeptons.acceptIdPair(PID::ELECTRON); // ** PID = Particle ID
			leadingLeptons.acceptIdPair(PID::MUON);
			leadingLeptons.acceptIdPair(12);
			leadingLeptons.acceptIdPair(14);
			declare(leadingLeptons, "leptons");

		 	book(h_ptz,8, 1,1);
		 	book(h_ptw,10, 1,1);
		 	book(h_mwz,12,1,1);
		 	book(h_phiwz,14,1,1);
		 	book(h_ptnu,16, 1,1);
		 	book(h_deltarap,18, 1,1);
		 	book(h_njet, 20,1,1);
		 	book(h_mjj, 22,1,1);
		 
		 	cout << "---------------- END INIT ------------------" << endl;
     
    }

   /// Perform the per-event analysis
	void analyze(const Event& event) {
	
		totalEvents++;

		Particles cand_electrons;
		Particles cand_muons;
		Particles cand_neutrinos;
		Jets cand_jets;
			
		double mw = 0;
		double mwz = 0; 
		unsigned int njet = 0;
		unsigned int n_charged_lep=0;

		//    you are getting a vector of the leptons for this event
		for( const Particle particle : apply<IdentifiedFinalState>(event,"leptons").particlesByPt() ){	 
		
				// Electron 1
				if ( cand_electrons.size() == 0 && abs(particle.eta())< 2.5 && abs(particle.pid()) == 11 && particle.pt() > 15*GeV) 	cand_electrons.push_back(particle);

				// Muon 
				else if ( abs(particle.eta())< 2.5 && abs(particle.pid()) == 13 && particle.pt() > 20*GeV ) 		cand_muons.push_back(particle);
				
				// Electron 2
				else if ( cand_electrons.size() == 1 && abs(particle.eta()) < 2.5 && abs(particle.pid())== 11 && particle.pt() > 15*GeV) 	cand_electrons.push_back(particle);
		
				//Neutrino
				else if ( abs(particle.pid())==14 ) cand_neutrinos.push_back(particle);
		}
		
		n_charged_lep = cand_muons.size() + cand_electrons.size();
		
		// ** here you are getting the jets with pT > 25 for this event
		//    then looping over them, and then requiring isolation. Nice :)
		for (const Jet& jet : apply<FastJets>(event, "jets").jetsByPt(20.0*GeV) ) {
			if ( abs( jet.eta() ) < 4.5 ) {  
				bool isolated = true;
				for (const Particle& neutrino : cand_neutrinos) {
					if (deltaR(neutrino, jet) < 0.1) {
						isolated = false;
						break;
				 	}
				} 
				if (isolated) {
					for (const Particle& lepton : cand_muons) {
						if (deltaR(lepton, jet) < 0.1) { 
							isolated = false;
							break;
						}
					}
					for (const Particle& lepton : cand_electrons) {
						if (deltaR(lepton, jet) < 0.1) { 
							isolated = false;
							break;
						}
					}
				}
				if (isolated && jet.momentum().pT()>25*GeV ){
					cand_jets.push_back(jet);
					njet++; 
				}
			}
		}	

		if ( n_charged_lep != 3 ) {
			//cout << "veto, event with selected leptons = " << n_charged_lep << endl;
			if ( cand_neutrinos.size() != 1 ) cout << "WARNING: Missing neutrino!! " << cand_neutrinos.size() << endl; 
			vetoedEvents++;			
			vetoEvent;
		}
		else if (cand_electrons.size()==2 && cand_muons.size() == 1 ){
		
			//Reconstruct missing energy				
			Particles vfs_particles = apply<VisibleFinalState>(event, "vfs").particles();

			// ** here you are calculating the missing transverse momentum (pTmiss)
			//    using charged particles. this could be done better ..
			FourMomentum temp;                                 
			for( const Particle & p : vfs_particles ) temp -= p.momentum();

			FourMomentum reco_nu, lepw, z, w;
			reco_nu.setEtaPhiME (0.,temp.phi(),0., temp.pT());
			lepw = cand_muons[0].momentum();
	  	z = cand_electrons[0].momentum() + cand_electrons[1].momentum();
	  	w = lepw + reco_nu;
	  				
			if (  2 * reco_nu.pT() * lepw.pT() * ( 1 - cos( deltaPhi(reco_nu.phi() , lepw.phi() ) ) ) < 0. ) {
				cout << "WARNING: something wrong in lepton selection !!!" << endl;
				vetoedEvents++;			
				vetoEvent;
			}
			mw = sqrt( 2 * reco_nu.pT() * lepw.pT() * ( 1 - cos( deltaPhi(reco_nu.phi() , lepw.phi() ) ) ) );
	  	
	   	//M WZ transverse
	   	double sum3pt = cand_muons[0].momentum().pt() + cand_electrons[0].momentum().pt() + cand_electrons[1].momentum().pt();
	   	double sum3px = cand_muons[0].momentum().px() + cand_electrons[0].momentum().px() + cand_electrons[1].momentum().px();
	   	double sum3py = cand_muons[0].momentum().py() + cand_electrons[0].momentum().py() + cand_electrons[1].momentum().py(); 
	  	mwz = sqrt( pow( sum3pt + reco_nu.Et(),2) - ( pow( sum3px + reco_nu.px(),2) + pow( sum3py + reco_nu.py(),2) ) );

	  	double pt_miss = cand_neutrinos[0].momentum().pT();
	  	//cout << "neutrino pt " << pt_miss << endl; 
	  	
	  	if ( abs( z.mass() - MZ ) > 10. ){
	  		//cout << "veto, phase space cuts not satisfied, mz " <<  mz << endl;
				vetoedEvents++;			
				vetoEvent;
	  	}
	  	else if( mw < 30.){
	  		//cout << "veto, phase space cuts not satisfied, mw " <<  mw << endl;
	  		vetoedEvents++;			
				vetoEvent;	
	  	}
  		else if( deltaR(cand_electrons[0].momentum(),  cand_electrons[1].momentum()) < 0.2 ||
	  					deltaR(cand_muons[0].momentum(),  cand_electrons[0].momentum()) < 0.3 ||
	  					deltaR(cand_muons[0].momentum(),  cand_electrons[1].momentum()) < 0.3  	){
	  		// cout << "veto, phase space cuts not satisfied,  deltaR " << endl;
	  		vetoedEvents++;			
				vetoEvent;	
	  	}
		  else{
	  		selectedEvents++; 
	  		h_ptz -> fill( z.pT() );
	  		h_mwz -> fill(mwz);
	  		h_phiwz -> fill( deltaPhi(z, w) );
	  		h_ptnu -> fill( reco_nu.pT() );
	  		h_deltarap -> fill( abs( z.rapidity() - lepw.rapidity()));
	  		h_ptw -> fill( w.pT() );
	  	}
		}else{
			cout << "WARNING: something wrong in lepton selection !!!" << endl;
			vetoedEvents++;			
			vetoEvent;
		}

		//selectet the leading jet 
		if ( njet > 2 ){
			jetsEvents++;
			h_mjj -> fill( (cand_jets[0].momentum() + cand_jets[1].momentum()).mass() );
		}
			
		h_njet->fill(njet);

	}

	void finalize() {

		std::cout << "---------------- FINALIZING ------------------" << endl;  
		double xsec = crossSection();
		cout << "xsec before cuts is: " << xsec << endl;
		cout << "selected Events = " << selectedEvents << ", vetoed = " << vetoedEvents << ", total= " << totalEvents << endl;
		
		float efficiency= selectedEvents/totalEvents;
		xsec = xsec*efficiency; 
		cout << "efficiency = " << efficiency << " , expected xsec after cuts = " << xsec << endl;	
		cout << "SumOfWeights() = " << sumOfWeights() << endl;
		
		const double sf(crossSection()/femtobarn/sumOfWeights());
  
    // scale to cross section		
    // this normalization is not correct to confront exp data with make-plots! It will give diff xs. 
    // The hepdata files contains Delta_xs 
    scale(h_njet,sf);
    scale(h_mwz, sf);
    scale(h_phiwz,sf);
    scale(h_ptnu,sf);
    scale(h_ptz, sf); 
    scale(h_ptw, sf);
    scale(h_deltarap, sf);
    scale(h_mjj, sf);    
    	

	}

private:

		Histo1DPtr h_njet;
		Histo1DPtr h_mwz;
		Histo1DPtr h_ptnu;
		Histo1DPtr h_phiwz;
		Histo1DPtr h_deltarap;
		Histo1DPtr h_ptz;
		Histo1DPtr h_ptw;
		Histo1DPtr h_mjj;

};

// The hook for the plugin system
DECLARE_RIVET_PLUGIN (WZ_ATLAS);

}
