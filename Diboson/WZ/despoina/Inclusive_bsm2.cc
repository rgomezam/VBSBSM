// -*- C++ -*-
//
#include "Rivet/Analysis.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/TauFinder.hh"
#include "Rivet/Projections/Smearing.hh"
#include "Rivet/Run.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {


  class TESTDET : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(TESTDET);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {



      Cut FS_Zlept = Cuts::abseta < 2.5 && Cuts::pT > 15*GeV;

      FinalState fs;
      Cut fs_z = Cuts::abseta < 2.5 && Cuts::pT > 15*GeV;
      Cut fs_j = Cuts::abseta < 4.5 && Cuts::pT > 25*GeV;

      PromptFinalState photons(Cuts::abspid == PID::PHOTON);

      // Electrons and muons in Fiducial PS
      PromptFinalState leptons(FinalState(fs_z && (Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON)));
      leptons.acceptTauDecays(false);
      DressedLeptons dressedleptons(photons, leptons, 0.1, FS_Zlept, true);
   //riv2   //addProjection(dressedleptons, "DressedLeptons");
        declare(dressedleptons, "DressedLeptons");
        
      // Electrons and muons in Total PS
      PromptFinalState leptons_total(Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON);
      leptons_total.acceptTauDecays(false);
      DressedLeptons dressedleptonsTotal(photons, leptons_total, 0.1, Cuts::open(), true);
      //addProjection(dressedleptonsTotal, "DressedLeptonsTotal");
      declare(dressedleptonsTotal, "DressedLeptonsTotal");


      // Promot neutrinos 
      IdentifiedFinalState nu_id;
      nu_id.acceptNeutrinos();
      PromptFinalState neutrinos(nu_id);
      neutrinos.acceptTauDecays(false);
      declare(neutrinos, "Neutrinos");
      MSG_WARNING("\033[91;1mLIMITED VALIDITY - check info file for details!\033[m");

      // Jets
      VetoedFinalState veto;
      veto.addVetoOnThisFinalState(dressedleptons);
      FastJets jets(veto, FastJets::ANTIKT, 0.4);
      declare(jets, "Jets");

/////////////////////////////////////////////////////////////////////

      // Initialise and register projections
       declare(FinalState(Cuts::abseta < 5 && Cuts::pT > 100*MeV), "FS");
//      declare(HeavyHadrons(Cuts::abseta < 2.5 && Cuts::pT > 5*GeV), "Bhadrons");

/*//riv2      _h_Sptlepton    = bookHisto1D("1_h_Sptlepton",pt_CUTS);
      _h_DeltaPhiWZ   = bookHisto1D("1_h_DeltaPhiWZ",DeltaPhiWZ_CUTS);
      _h_mtWZ  = bookHisto1D("1_h_mtWZ",mtWZ_CUTS);
      _h_Njets=bookHisto1D("1_h_Njets",Njets_CUTS);
      _h_Deltayjj=bookHisto1D("1_h_Deltayjj",Deltayjj_CUTS);
      _h_mjj=bookHisto1D("1_h_mjj",mjj_CUTS);
      _h_DeltaPhijj=bookHisto1D("1_h_DeltaPhijj",DeltaPhijj_CUTS);
      _h_Njgap=bookHisto1D("1_h_Njgap",Njgap_CUTS);
      _histo = bookHisto1D("SigmaFid",1,0.5,1.5);
*/

/*     book(_h_Sptlepton,pt_CUTS);
       book(_h_DeltaPhiWZ,DeltaPhiWZ_CUTS);
       book(_h_mtWZ,mtWZ_CUTS);
       book(_h_Njets,Njets_CUTS);
       book(_h_Deltayjj,Deltayjj_CUTS);
       book(_h_mjj,mjj_CUTS);
       book(_h_DeltaPhijj,DeltaPhijj_CUTS);
       book(_h_Njgap,Njgap_CUTS);
      book(_histo,1,0.5,1.5);
  */      
     
       book(_h_ptZ,8,1,1);  
       book(_h_ptW,10,1,1);  
       book(_h_mtWZ,12,1,1); 
       book(_h_DeltaPhiWZ,14,1,1);
       book(_h_ptv,16,1,1);
        book(_histo,"1_histo",sigma);

/*       book(_h_Sptlepton,8,1,1);
       book(_h_mtWZ,6,1,1);   
       book(_h_DeltaPhiWZ,10,1,1);
       book(_h_Njets,12,1,1);
       book(_h_Deltayjj,16,1,1);
       book(_h_mjj,14,1,1);        
       book(_h_DeltaPhijj,18,1,1);
       book(_h_Njgap,20,1,1);
       book(_histo,"1_histo",sigma);
      book(_h_mtWZ,"1_h_mtWZ",mtWZ_CUTS); 
      book(_h_mtWZ_2t_clip,"1_h_mtWZ_2t_clip",mtWZ_CUTS); 
      book(_h_mtWZ_1_5t_clip,"1_h_mtWZ_1_5t_clip",mtWZ_CUTS); 
      book(_h_mtWZ_1t_clip,"1_h_mtWZ_1t_clip",mtWZ_CUTS); 
  */
       cnt_fid=0;

    }

      Jet getLeadingJetOnOtherSide(const Jets& jets, bool& two_jets) {
          double eta1 = jets.at(0).eta();
        two_jets = false;
    Jet jet2;
    for (const Jet& jet : jets) {
      if(jet.eta()*eta1 < 0.) {
        jet2 = jet;
        two_jets = true;
        break;
      }
    }
    return jet2;
      }



    /// Perform the per-event analysis
    void analyze(const Event& event) {
//      const double weight = event.weight();
   
// ParticlePair beams(const Event& event);
// double sqrts = sqrtS(beams);
        const vector<DressedLepton>& dressedleptons = apply<DressedLeptons>(event, "DressedLeptons").dressedLeptons();
      const vector<DressedLepton>& dressedleptonsTotal = apply<DressedLeptons>(event, "DressedLeptonsTotal").dressedLeptons();
      const Particles& neutrinos = apply<PromptFinalState>(event, "Neutrinos").particlesByPt();
     Jets jets = apply<JetAlg>(event, "Jets").jetsByPt( (Cuts::abseta < 4.5) && (Cuts::pT > 25*GeV) );

      if (dressedleptonsTotal.size() < 3 || neutrinos.size() < 1) vetoEvent;
      //---Total PS: assign leptons to W and Z bosons using Resonant shape algorithm
      // NB: This resonant shape algorithm assumes the Standard Model and can therefore
      //     NOT be used for any kind of reinterpretation in terms of new-physics models..

      int i, j, k;
      double MassZ01 = 0., MassZ02 = 0., MassZ12 = 0.;
      double MassW0 = 0., MassW1 = 0., MassW2 = 0.;
      double WeightZ1, WeightZ2, WeightZ3;
      double WeightW1, WeightW2, WeightW3;
      double M1, M2, M3;
      double WeightTotal1, WeightTotal2, WeightTotal3;

      //try Z pair of leptons 01
      if ( (dressedleptonsTotal[0].pid() ==-(dressedleptonsTotal[1].pid())) && (dressedleptonsTotal[2].abspid()==neutrinos[0].abspid()-1)){
        MassZ01 = (dressedleptonsTotal[0].momentum()+dressedleptonsTotal[1].momentum()).mass();
        MassW2 = (dressedleptonsTotal[2].momentum()+neutrinos[0].momentum()).mass();
      }
      //try Z pair of leptons 02
      if ( (dressedleptonsTotal[0].pid()==-(dressedleptonsTotal[2].pid())) && (dressedleptonsTotal[1].abspid()==neutrinos[0].abspid()-1)){
        MassZ02 = (dressedleptonsTotal[0].momentum()+dressedleptonsTotal[2].momentum()).mass();
        MassW1 = (dressedleptonsTotal[1].momentum()+neutrinos[0].momentum()).mass();
      }
      //try Z pair of leptons 12
      if ( (dressedleptonsTotal[1].pid()==-(dressedleptonsTotal[2].pid())) && (dressedleptonsTotal[0].abspid()==neutrinos[0].abspid()-1)){
        MassZ12 = (dressedleptonsTotal[1].momentum()+dressedleptonsTotal[2].momentum()).mass();
        MassW0 = (dressedleptonsTotal[0].momentum()+neutrinos[0].momentum()).mass();
      }
      WeightZ1 = 1/(pow(MassZ01*MassZ01 - MZ_PDG*MZ_PDG,2) + pow(MZ_PDG*GammaZ_PDG,2));
      WeightW1 = 1/(pow(MassW2*MassW2 - MW_PDG*MW_PDG,2) + pow(MW_PDG*GammaW_PDG,2));
      WeightTotal1 = WeightZ1*WeightW1;
      M1 = -1*WeightTotal1;

      WeightZ2 = 1/(pow(MassZ02*MassZ02- MZ_PDG*MZ_PDG,2) + pow(MZ_PDG*GammaZ_PDG,2));
      WeightW2 = 1/(pow(MassW1*MassW1- MW_PDG*MW_PDG,2) + pow(MW_PDG*GammaW_PDG,2));
      WeightTotal2 = WeightZ2*WeightW2;
      M2 = -1*WeightTotal2;

      WeightZ3 = 1/(pow(MassZ12*MassZ12 - MZ_PDG*MZ_PDG,2) + pow(MZ_PDG*GammaZ_PDG,2));
      WeightW3 = 1/(pow(MassW0*MassW0 - MW_PDG*MW_PDG,2) + pow(MW_PDG*GammaW_PDG,2));
      WeightTotal3 = WeightZ3*WeightW3;
      M3 = -1*WeightTotal3;

      if( (M1 < M2 && M1 < M3) || (MassZ01 != 0 && MassW2 != 0 && MassZ02 == 0 && MassZ12 == 0) ){
        i = 0; j = 1; k = 2;
      }
      if((M2 < M1 && M2 < M3) || (MassZ02 != 0 && MassW1 != 0 && MassZ01 == 0 && MassZ12 == 0) ){
        i = 0; j = 2; k = 1;
      }
      if((M3 < M1 && M3 < M2) || (MassZ12 != 0 && MassW0 != 0 && MassZ01 == 0 && MassZ02 == 0) ){
        i = 1; j = 2; k = 0;
      }

      FourMomentum ZbosonTotal   = dressedleptonsTotal[i].momentum()+dressedleptonsTotal[j].momentum();
      if (!( ZbosonTotal.mass() >= 66*GeV && ZbosonTotal.mass() <= 116*GeV) ) vetoEvent;

      //---end Total PS


      //---Fiducial PS: assign leptons to W and Z bosons using Resonant shape algorithm
      if (dressedleptons.size() < 3)  vetoEvent;

      int EventType = -1;
      int Nel = 0, Nmu = 0;

      for (const DressedLepton& l : dressedleptons) {
        if (l.abspid() == 11)  ++Nel;
        if (l.abspid() == 13)  ++Nmu;
      }

      if ( (Nel == 3)  && (Nmu==0) ) { EventType = 3; //cout<<" eee"<<endl;
      }
      if ( (Nel == 2)  && (Nmu==1) ) { EventType = 2; //cout<<" eem"<<endl;
      }
      if ( (Nel == 1)  && (Nmu==2) ) { EventType = 1; //cout<<" mee"<<endl;
      }
      if ( (Nel == 0)  && (Nmu==3) ) { EventType = 0; //cout<<" mmm"<<endl;
      }


      int EventCharge = -dressedleptons[0].charge() * dressedleptons[1].charge() * dressedleptons[2].charge();

      MassZ01 = 0; MassZ02 = 0; MassZ12 = 0;
      MassW0 = 0;  MassW1 = 0;  MassW2 = 0;

      // try Z pair of leptons 01
      if (dressedleptons[0].pid() == -dressedleptons[1].pid()) {
        MassZ01 = (dressedleptons[0].momentum() + dressedleptons[1].momentum()).mass();
        MassW2 = (dressedleptons[2].momentum() + neutrinos[0].momentum()).mass();
      }
      // try Z pair of leptons 02
      if (dressedleptons[0].pid() == -dressedleptons[2].pid()) {
        MassZ02 = (dressedleptons[0].momentum() + dressedleptons[2].momentum()).mass();
        MassW1 = (dressedleptons[1].momentum() + neutrinos[0].momentum()).mass();
      }
      // try Z pair of leptons 12
      if (dressedleptons[1].pid() == -dressedleptons[2].pid()) {
        MassZ12 = (dressedleptons[1].momentum() + dressedleptons[2].momentum()).mass();
        MassW0 = (dressedleptons[0].momentum() + neutrinos[0].momentum()).mass();
      }
      WeightZ1 = 1/(pow(MassZ01*MassZ01 - MZ_PDG*MZ_PDG,2) + pow(MZ_PDG*GammaZ_PDG,2));
      WeightW1 = 1/(pow(MassW2*MassW2 - MW_PDG*MW_PDG,2) + pow(MW_PDG*GammaW_PDG,2));
      WeightTotal1 = WeightZ1*WeightW1;
      M1 = -1*WeightTotal1;

      WeightZ2 = 1/(pow(MassZ02*MassZ02- MZ_PDG*MZ_PDG,2) + pow(MZ_PDG*GammaZ_PDG,2));
      WeightW2 = 1/(pow(MassW1*MassW1- MW_PDG*MW_PDG,2) + pow(MW_PDG*GammaW_PDG,2));
      WeightTotal2 = WeightZ2*WeightW2;
      M2 = -1*WeightTotal2;

      WeightZ3 = 1/(pow(MassZ12*MassZ12 - MZ_PDG*MZ_PDG,2) + pow(MZ_PDG*GammaZ_PDG,2));
      WeightW3 = 1/(pow(MassW0*MassW0 - MW_PDG*MW_PDG,2) + pow(MW_PDG*GammaW_PDG,2));
      WeightTotal3 = WeightZ3*WeightW3;
      M3 = -1*WeightTotal3;

      if( (M1 < M2 && M1 < M3) || (MassZ01 != 0 && MassW2 != 0 && MassZ02 == 0 && MassZ12 == 0) ) {
        i = 0; j = 1; k = 2;
      }
      if((M2 < M1 && M2 < M3) || (MassZ02 != 0 && MassW1 != 0 && MassZ01 == 0 && MassZ12 == 0) ) {
        i = 0; j = 2; k = 1;
      }
      if((M3 < M1 && M3 < M2) || (MassZ12 != 0 && MassW0 != 0 && MassZ01 == 0 && MassZ02 == 0) ) {
        i = 1; j = 2; k = 0;
      }


      FourMomentum Zlepton1 = dressedleptons[i].momentum();
      FourMomentum Zlepton2 = dressedleptons[j].momentum();
      FourMomentum Wlepton  = dressedleptons[k].momentum();
      FourMomentum Zboson   = dressedleptons[i].momentum()+dressedleptons[j].momentum();
      FourMomentum Wboson   = dressedleptons[k].momentum()+neutrinos[0].momentum();
      FourMomentum WZ       = Zboson+Wboson;
      double mWZ            = WZ.mass()/GeV;

      double WZ_pt = Zlepton1.pt() + Zlepton2.pt() + Wlepton.pt() + neutrinos[0].pt();
      double WZ_px = Zlepton1.px() + Zlepton2.px() + Wlepton.px() + neutrinos[0].px();
      double WZ_py = Zlepton1.py() + Zlepton2.py() + Wlepton.py() + neutrinos[0].py();
      double mTWZ = sqrt( pow(WZ_pt, 2) - ( pow(WZ_px, 2) + pow(WZ_py,2) ) )/GeV;
      double Wboson_mT = sqrt( 2 * Wlepton.pT() * neutrinos[0].pt() * (1 - cos(deltaPhi(Wlepton, neutrinos[0]))) );

      if (fabs(Zboson.mass()/GeV - MZ_PDG) >= 10.) vetoEvent;
      if (Wboson_mT <= 30*GeV)                     vetoEvent;
      if (Wlepton.pT() <= 20*GeV)                  vetoEvent;

//      if (deltaR(Zlepton1, jet1) <= 0.3)        vetoEvent;
//      if (deltaR(Zlepton2, jet1) <= 0.3)        vetoEvent;
//      if (deltaR(Wlepton, jet1) <= 0.3)        vetoEvent;

//      if (deltaR(Zlepton1, jet2) <= 0.3)        vetoEvent;
//      if (deltaR(Zlepton2, jet2) <= 0.3)        vetoEvent;
//      if (deltaR(Wlepton, jet2) <= 0.3)        vetoEvent;

      if (deltaR(Zlepton1, Wlepton)  <= 0.3)        vetoEvent;
      if (deltaR(Zlepton2, Wlepton)  <= 0.3)        vetoEvent;
      if (deltaR(Zlepton1, Zlepton2) <= 0.2)        vetoEvent;
     // Jets alljets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 4.5);

/*        ifilter_discard(jets, [&](const Jet& j) {
        return (fabs(j.eta()) > 4.5 && j.pt()/GeV < 40);
          });
        for (const Particle& l : dressedleptons)
          ifilter_discard(jets, deltaRLess(l, 0.3));
      const Jets alljets = jets;

        // put in ZZ VBS analysis cuts here
        if (alljets.size() < 2)  vetoEvent;

      if (alljets.at(0).pt()< 40.*GeV) vetoEvent;
      if (alljets.at(1).pt()< 40.*GeV) vetoEvent;
      if (fabs(alljets.at(0).eta())>4.5) vetoEvent;
      if (fabs(alljets.at(1).eta())>4.5) vetoEvent;

        bool two_jets = false; 
        Jet jet2 = getLeadingJetOnOtherSide(jets,two_jets);   

        if (!two_jets) vetoEvent;

        double mjj = (jets.at(0).mom() + jet2.mom()).mass()/GeV;
        if (mjj < 500.) vetoEvent;
*/
        _histo->fill(1 );



      cnt_fid++;
////////////////////////////////////////////////////////////////////
      //cout<<"mWZ  "<<mWZ<<endl;
    //  cout<<"mTWZ  "<<mTWZ<<endl;
      
      
//      sptl=Zlepton1.pT()+Zlepton2.pT()+Wlepton.pT();


  //    DeltaYjj = fabs(jet1.rapidity()-jet2.rapidity());
 /*    DeltaYjj=1; 
      if(sptl/GeV<650){
          _h_Sptlepton->fill(sptl/GeV);
      }else{
          _h_Sptlepton->fill(600);
      }
   */

      if(mTWZ>=600.){  
          _h_ptZ->fill(Zboson.pt() );
          _h_ptW->fill(Wboson.pt() );

          _h_DeltaPhiWZ->fill(acos(cos(Zboson.phi()-Wboson.phi())));

          _h_ptv->fill(neutrinos[0].pt());
          if(mTWZ>=1000.){
              _h_mtWZ ->fill(999.);}
          else{
              _h_mtWZ ->fill(mTWZ);
          }
      }

  //    if (alljets.size() < 5)  _h_Njets->fill(alljets.size() );
//      else  _h_Njets->fill(5 );
//
    /*  if(DeltaYjj<7.0){
          _h_Deltayjj->fill(DeltaYjj );
      }else{
          _h_Deltayjj->fill(6 );
      }

      if(mjj/GeV <2000.) _h_mjj->fill(mjj/GeV );
      else  _h_mjj->fill(1700 );
*/
      /*      if(    ( (jet1+jet2).mass())/GeV<2000 ){
          _h_mjj->fill(((jet1+jet2).mass())/GeV );
      }  else{
          _h_mjj->fill(1700 );
      }
      _h_DeltaPhijj->fill(acos(cos(jet1.phi()-jet2.phi())) );

      NjGap = 0.;
      for(int i=1; i<jets.size(); i++){
          if(jets[i].eta()==jet1eta || (jets[i].eta()==jet2eta) )continue;         
          if(jets[i].eta()>0){ 
              if(jet1eta>0){
                  if(jets[i].eta()<jet1eta) NjGap++;
              }else if(jet2eta>0){            
                  if(jets[i].eta()<jet2eta) NjGap++;
              }
          }else{
              if(jet1eta<0){
                  if(jets[i].eta()>jet1eta) NjGap++;
              }else if(jet2eta<0){            
                  if(jets[i].eta()>jet2eta) NjGap++;
              }
          }
      }
    
      if (NjGap < 3)  _h_Njgap->fill(NjGap );
      else  _h_Njgap->fill(3 );
*/}
    /// Normalise histograms etc., after the run
    void finalize() {
     const double xs_pb(crossSection() / picobarn);
     const double xs_fb(crossSection() / femtobarn);
      const double sumw(sumOfWeights());
      const double sf_pb(xs_pb / sumw);
      const double sf_fb(xs_fb / sumw);
cout<<"XSECTION  in fb  "<<crossSection() / femtobarn<<endl;
cout<<"SUM of W  "<<sumw<<endl;
cout<<"sf  "<<sf_fb<<endl;
cout<<"EVENTS PASSED "<<cnt_fid<<endl;

      scale(_histo, sf_fb/4.);
      //scale(_h_Sptlepton, sf_fb/4.);
    //  scale(_h_DeltaPhiWZ, sf_fb/4.);
      scale(_h_mtWZ , sf_fb/4.);
      scale(_h_DeltaPhiWZ , sf_fb/4.);
      scale(_h_ptZ , sf_fb/4.);
      scale(_h_ptW , sf_fb/4.);
      scale(_h_ptv , sf_fb/4.);
  
     
      //   scale(_h_mtWZ_2t_clip , sf_fb/4.);
    //  scale(_h_mtWZ_1_5t_clip , sf_fb/4.);
     // scale(_h_mtWZ_1t_clip , sf_fb/4.);
     // scale(_h_Njets , sf_fb/4.);
   ///   scale(_h_Deltayjj , sf_fb/4.);
   //   scale(_h_mjj , sf_fb/4.);
   //   scale(_h_DeltaPhijj , sf_fb/4.);
   //   scale(_h_Njgap , sf_fb/4.);

    }

    //@}


  private:

    /// @name Histograms
    Histo1DPtr _histo,_h_mtWZ_clipped_1;
    Histo1DPtr _h_Sptlepton,_h_DeltaPhiWZ,_h_mtWZ,_h_Njets,_h_Deltayjj,_h_mjj,_h_DeltaPhijj,_h_Njgap;
   Histo1DPtr _h_ptv,_h_ptZ,_h_ptW; 
    double MZ_PDG = 91.1876;
    double MW_PDG = 80.385;
    double GammaZ_PDG = 2.4952;
    double GammaW_PDG = 2.085;
    const vector<double> pt_CUTS = {60,150,250,350,500,650};
    const vector<double> DeltaPhiWZ_CUTS = { 0.0, 0.6, 1.2,1.8,2.5,3.15 };
    const vector<double> mtWZ_CUTS = { 150,200,250,300,400,500};
    const vector<double> Njets_CUTS = { 1.5,2.5,3.5,4.5,5.5 };
    const vector<double> Deltayjj_CUTS = { 0.0, 2.5, 3.5,4.5,7.0 };
    const vector<double> mjj_CUTS = {500,700,1000,1500,2000};
    const vector<double> DeltaPhijj_CUTS = { 0.0, 1.0,1.8,2.3,2.8,3.15};
    const vector<double> Njgap_CUTS = {-0.5,0.5, 1.5,2.5,3.5};
    const vector<double> sigma ={-0.5,1.5};
    double sptl;
    double DeltaYjj;
    double NjGap;
    int cnt_fid;

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TESTDET);

}

