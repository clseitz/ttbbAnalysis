#include "PhysicsTools/KinFitter/interface/TFitConstraintM.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEtEtaPhi.h"
#include "PhysicsTools/KinFitter/interface/TKinFitter.h"

#include "TLorentzVector.h"

#include <iostream>
#include <algorithm>
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Modified version of the ExampleEtEtaPhi2CFit.C macro from CMS AN 2005/025.
// Jet error parametrization from CMS AN 2005/005.
//
// To run this macro in a root session do:
// root [0] gSystem->Load("libPhysicsToolsKinFitter.so");
// root [1] .x PhysicsTools/KinFitter/test/testKinFitFWLite.C+
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

TLorentzVector FittedValues(TFitParticleEtEtaPhi * part){
  TLorentzVector update_part;
  update_part.SetPtEtaPhiE(part->getCurr4Vec()->Pt(),part->getCurr4Vec()->Eta(),part->getCurr4Vec()->Phi(),part->getCurr4Vec()->E());
  return update_part;

}



Double_t ErrEt(Float_t Et, Float_t Eta) {
  Double_t InvPerr2, a, b, c;
  if(fabs(Eta) < 1.4){
    a = 5.6;
    b = 1.25;
    c = 0.033;
  }
  else{
    a = 4.8;
    b = 0.89;
    c = 0.043;
  }
  InvPerr2 = (a * a) + (b * b) * Et + (c * c) * Et * Et;
  return InvPerr2;
}

Double_t ErrEta(Float_t Et, Float_t Eta) {
  Double_t InvPerr2, a, b, c;
  if(fabs(Eta) < 1.4){
    a = 1.215;
    b = 0.037;
    c = 7.941 * 0.0001;
  }
  else{
    a = 1.773;
    b = 0.034;
    c = 3.56 * 0.0001;
  }
  InvPerr2 = a/(Et * Et) + b/Et + c;
  return InvPerr2;
}

Double_t ErrPhi(Float_t Et, Float_t Eta) {
  Double_t InvPerr2, a, b, c;
  if(fabs(Eta) < 1.4){
    a = 6.65;
    b = 0.04;
    c = 8.49 * 0.00001;
  }
  else{
    a = 2.908;
    b = 0.021;
    c = 2.59 * 0.0001;
  }
  InvPerr2 = a/(Et * Et) + b/Et + c;
  return InvPerr2;
}


Double_t*  testKinFitFWLite(TLorentzVector vA, TLorentzVector vB, TLorentzVector vC,TLorentzVector vD, TLorentzVector vE, TLorentzVector vF, int constrained = 0)
{

  std::vector<TLorentzVector >  vJets;
  //  std::vector<TLorentzVector >  fitJets;

  vJets.push_back(vA);
  vJets.push_back(vB);
  vJets.push_back(vC);
  vJets.push_back(vD);
  vJets.push_back(vE);
  vJets.push_back(vF);


  const unsigned short size = 6;
  int numbers[size] = {0, 1, 2, 3, 4, 5};

  unsigned int numJetForFit  = 6;
  int converged = 0;
  float minChi2 = 999;
  float corr_top1_m, corr_top1_pt, corr_top1_eta, corr_top1_phi;
  corr_top1_m  =corr_top1_pt = corr_top1_eta = corr_top1_phi = 0;
  float corr_top2_m, corr_top2_pt, corr_top2_eta, corr_top2_phi;
  corr_top2_m  =corr_top2_pt = corr_top2_eta = corr_top2_phi = 0;
  Double_t probChi2 = 0;
  int bestCombo[size] = {0,0,0,0,0,0};

  if (1==1){
    do { 	 
      if (numbers[0] == 0 && numbers[1] == 1 && numbers[2] < numbers[3] && numbers[4] < numbers[5]){

	/*
	if(tagged){
	  for(unsigned short i=0; i<size; ++i) {
	    std::cout << numbers[i] << (i+1!=size?" ":"\n");
	  }		
	}
	*/
	

	TLorentzVector v1;	TLorentzVector v2; 	TLorentzVector v3;
	TLorentzVector v4;	TLorentzVector v5; 	TLorentzVector v6;

	TLorentzVector fitted_b1;	TLorentzVector fitted_bbar; 	TLorentzVector fitted_wq;
	TLorentzVector fitted_wqbar;	TLorentzVector fitted_wp; 	TLorentzVector fitted_wpbar;

	v1=vJets[numbers[0]]; 
	v2=vJets[numbers[1]]; 
	v3=vJets[numbers[2]]; 

	v4=vJets[numbers[3]]; 
	v5=vJets[numbers[4]]; 
	v6=vJets[numbers[5]]; 


	TMatrixD m1(3,3);
	TMatrixD m2(3,3);
	TMatrixD m3(3,3);
	TMatrixD m4(3,3);
	TMatrixD m5(3,3);
	TMatrixD m6(3,3);
	m1.Zero();
	m2.Zero();
	m3.Zero();
	m4.Zero();
	m5.Zero();
	m6.Zero();

	//In this example the covariant matrix depends on the transverse energy and eta of the jets
	m1(0,0) = ErrEt (v1.Et(), v1.Eta()); // et
	m1(1,1) = ErrEta(v1.Et(), v1.Eta()); // eta
	m1(2,2) = ErrPhi(v1.Et(), v1.Eta()); // phi
	m2(0,0) = ErrEt (v2.Et(), v2.Eta()); // et
	m2(1,1) = ErrEta(v2.Et(), v2.Eta()); // eta
	m2(2,2) = ErrPhi(v2.Et(), v2.Eta()); // phi
	m3(0,0) = ErrEt (v3.Et(), v3.Eta()); // et
	m3(1,1) = ErrEta(v3.Et(), v3.Eta()); // eta
	m3(2,2) = ErrPhi(v3.Et(), v3.Eta()); // phi

	m4(0,0) = ErrEt (v4.Et(), v4.Eta()); // et
	m4(1,1) = ErrEta(v4.Et(), v4.Eta()); // eta
	m4(2,2) = ErrPhi(v4.Et(), v4.Eta()); // phi

	m5(0,0) = ErrEt (v5.Et(), v5.Eta()); // et
	m5(1,1) = ErrEta(v5.Et(), v5.Eta()); // eta
	m5(2,2) = ErrPhi(v5.Et(), v5.Eta()); // phi

	m6(0,0) = ErrEt (v6.Et(), v6.Eta()); // et
	m6(1,1) = ErrEta(v6.Et(), v6.Eta()); // eta
	m6(2,2) = ErrPhi(v6.Et(), v6.Eta()); // phi
	
	TFitParticleEtEtaPhi *jetB = new TFitParticleEtEtaPhi( "jetB", "jetB", &v1, &m1 );
	TFitParticleEtEtaPhi *jetBbar = new TFitParticleEtEtaPhi( "jetBbar", "jetBbar", &v2, &m2 );
	TFitParticleEtEtaPhi *jetWq = new TFitParticleEtEtaPhi( "jetWq", "jetWq", &v3, &m3 );
	TFitParticleEtEtaPhi *jetWqbar = new TFitParticleEtEtaPhi( "jetWqbar", "jetWqbar", &v4, &m4 );
	TFitParticleEtEtaPhi *jetWp = new TFitParticleEtEtaPhi( "jetWp", "jetWp", &v5, &m5 );
	TFitParticleEtEtaPhi *jetWpbar = new TFitParticleEtEtaPhi( "jetWpbar", "jetWpbar", &v6, &m6 );
	

	//vec3 and vec4 must make a W boson
	TFitConstraintM *mCons1 = new TFitConstraintM( "WMassConstraint", "WMass-Constraint", 0, 0 , 80.4);
	mCons1->addParticles1( jetWq, jetWqbar );
	//vec5 and vec6 must make a W boson
	TFitConstraintM *mCons2 = new TFitConstraintM( "WMassConstraint", "WMass-Constraint", 0, 0 , 80.4);
	mCons2->addParticles1( jetWp, jetWpbar );

	TFitConstraintM *mCons3 = new TFitConstraintM( "EqualMasses", "EqualMasses", 0, 0 , 0.0);
	mCons3->addParticles1( jetWp, jetWpbar, jetBbar );
	mCons3->addParticles2( jetWq, jetWqbar, jetB  );

	TFitConstraintM *mCons4 = new TFitConstraintM( "TopMassConstraint", "TopMass-Constraint", 0, 0 , 173.);
	mCons4->addParticles1( jetWq, jetWqbar, jetB );
	TFitConstraintM *mCons5 = new TFitConstraintM( "TopMassConstraint", "TopMass-Constraint", 0, 0 , 173.);
	mCons5->addParticles1( jetWp, jetWpbar, jetBbar );
	
	//Definition of the fitter
	//Add  measured particles(jets)
	//Add  constraints
	TKinFitter* fitter = new TKinFitter("fitter", "fitter");
	fitter->addMeasParticle( jetWq );
	fitter->addMeasParticle( jetWqbar );
	fitter->addMeasParticle( jetB );

	fitter->addMeasParticle( jetWp );
	fitter->addMeasParticle( jetWpbar );
	fitter->addMeasParticle( jetBbar );
	fitter->addConstraint( mCons1 );
	fitter->addConstraint( mCons2 );
	if(!constrained)fitter->addConstraint( mCons3 );
	
	  if(constrained){
	    fitter->addConstraint( mCons5 );
	    fitter->addConstraint( mCons4 );
	  }

	
	//Set convergence criteira
	fitter->setMaxNbIter(500);
	fitter->setMaxDeltaS( 5e-3 );
	fitter->setMaxF( 1e-4 );
	fitter->setVerbosity(1);

	converged = fitter->fit();
	
	float Chi2 = fitter->getS();

	if (Chi2 < minChi2){
	  minChi2 = Chi2;
	  probChi2 = TMath::Prob(fitter->getS(),fitter->getNDF());
	  
	  fitted_b1=FittedValues(jetB);
	  fitted_bbar=FittedValues(jetBbar);
	  fitted_wq=FittedValues(jetWq);
	  fitted_wqbar=FittedValues(jetWqbar);
	  fitted_wp=FittedValues(jetWp);
	  fitted_wpbar=FittedValues(jetWpbar);
	  corr_top1_m = (fitted_b1 + fitted_wq + fitted_wqbar).M();
	  corr_top1_pt = (fitted_b1 + fitted_wq + fitted_wqbar).Pt();
	  corr_top1_eta = (fitted_b1 + fitted_wq + fitted_wqbar).Eta();
	  corr_top1_phi = (fitted_b1 + fitted_wq + fitted_wqbar).Phi();
	  
	  corr_top2_m = (fitted_bbar + fitted_wp + fitted_wpbar).M();
	  corr_top2_pt = (fitted_bbar + fitted_wp + fitted_wpbar).Pt();
	  corr_top2_eta = (fitted_bbar + fitted_wp + fitted_wpbar).Eta();
	  corr_top2_phi = (fitted_bbar + fitted_wp + fitted_wpbar).Phi();
	  
	  for (unsigned int b=0; b < size; b++){ 
	    bestCombo[b] = numbers[b];
	  }
	}
	
  
	delete jetWq;
	delete jetWqbar;
	delete jetB;
	delete jetWp;
	delete jetWpbar;
	delete jetBbar;
	delete mCons1;
	delete mCons2;
	delete mCons3;
	delete mCons4;
	delete mCons5;


	delete fitter;
      }
    }while(std::next_permutation(numbers, numbers + size));
  }
  

  //  cout<<"minChi2 "<<minChi2<<endl;

  Int_t ncombo=0;

  for (unsigned int b=0; b < size; b++) ncombo += bestCombo[b]*pow(10,size-b -1);


  //cout << ncombo << " ncombo " << endl;
  static Double_t fitJets[12];
  fitJets[0] = converged;
  fitJets[1] = ncombo;
  fitJets[2] = minChi2;
  fitJets[3] = probChi2;
  fitJets[4] = corr_top1_pt;
  fitJets[5] = corr_top1_eta;
  fitJets[6] = corr_top1_phi;
  fitJets[7] = corr_top1_m;
  fitJets[8] = corr_top2_pt;
  fitJets[9] = corr_top2_eta;
  fitJets[10] = corr_top2_phi;
  fitJets[11] = corr_top2_m;

  return fitJets;
}
