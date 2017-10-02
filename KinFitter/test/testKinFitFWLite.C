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

void print(TKinFitter *fitter)
{
  std::cout << "=============================================" << std ::endl;
  std::cout << "-> Number of measured Particles  : " << fitter->nbMeasParticles() << std::endl;
  std::cout << "-> Number of unmeasured particles: " << fitter->nbUnmeasParticles() << std::endl;
  std::cout << "-> Number of constraints         : " << fitter->nbConstraints() << std::endl;
  std::cout << "-> Number of degrees of freedom  : " << fitter->getNDF() << std::endl;
  std::cout << "-> Number of parameters A        : " << fitter->getNParA() << std::endl;
  std::cout << "-> Number of parameters B        : " << fitter->getNParB() << std::endl;
  std::cout << "-> Maximum number of iterations  : " << fitter->getMaxNumberIter() << std::endl;
  std::cout << "-> Maximum deltaS                : " << fitter->getMaxDeltaS() << std::endl;
  std::cout << "-> Maximum F                     : " << fitter->getMaxF() << std::endl;
  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++" << std ::endl;
  std::cout << "-> Status                        : " << fitter->getStatus() << std::endl;
  std::cout << "-> Number of iterations          : " << fitter->getNbIter() << std::endl;
  std::cout << "-> S                             : " << fitter->getS() << std::endl;
  std::cout << "-> F                             : " << fitter->getF() << std::endl;
  std::cout << "=============================================" << std ::endl;
}

Double_t*  testKinFitFWLite(TLorentzVector vA, TLorentzVector vB, TLorentzVector vC,TLorentzVector vD, TLorentzVector vE, TLorentzVector vF, int tagged=0)
{

  std::vector<TLorentzVector >  vJets;
  //  std::vector<TLorentzVector >  fitJets;
  static Double_t fitJets[2];
  vJets.push_back(vA);
  vJets.push_back(vB);
  vJets.push_back(vC);
  vJets.push_back(vD);
  vJets.push_back(vE);
  vJets.push_back(vF);


  const unsigned short size = 6;
  int numbers[size] = {0, 1, 2, 3, 4, 5};

  unsigned int numJetForFit  = 6;
  int iteration = 0;
  float minChi2 = 999;
  float minChi2Iteration = 0;
  int bestCombo[size] = {0,0,0,0,0,0};

  if (1==1){
    do { 	 
      if (numbers[0] < numbers[1] && numbers[2] < numbers[3] && numbers[4] < numbers[5]){
	if(tagged && (numbers[0] != 0 || numbers[1] != 1)) continue;
	/*
	if(tagged){
	  for(unsigned short i=0; i<size; ++i) {
	    std::cout << numbers[i] << (i+1!=size?" ":"\n");
	  }
		
		
	}
	*/

	
	TLorentzVector v1;	TLorentzVector v2; 	TLorentzVector v3;
	TLorentzVector v4;	TLorentzVector v5; 	TLorentzVector v6;
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
	//vec3 and vec4 and vec1 must make a top quark
	TFitConstraintM *mCons2 = new TFitConstraintM( "TopMassConstraint", "TopMass-Constraint", 0, 0 , 175.);
	mCons2->addParticles1( jetWq, jetWqbar, jetB );

	//vec5 and vec6 must make a W boson
	TFitConstraintM *mCons3 = new TFitConstraintM( "WMassConstraint", "WMass-Constraint", 0, 0 , 80.4);
	mCons3->addParticles1( jetWp, jetWpbar );
	//vec5 and vec6 and vec2 must make a top quark
	TFitConstraintM *mCons4 = new TFitConstraintM( "TopMassConstraint", "TopMass-Constraint", 0, 0 , 175.);
	mCons4->addParticles1( jetWp, jetWpbar, jetBbar );

	//Definition of the fitter
	//Add three measured particles(jets)
	//Add two constraints
	TKinFitter* fitter = new TKinFitter("fitter", "fitter");
	fitter->addMeasParticle( jetWq );
	fitter->addMeasParticle( jetWqbar );
	fitter->addMeasParticle( jetB );

	fitter->addMeasParticle( jetWp );
	fitter->addMeasParticle( jetWpbar );
	fitter->addMeasParticle( jetBbar );
	fitter->addConstraint( mCons1 );
	fitter->addConstraint( mCons2 );
	fitter->addConstraint( mCons3 );
	fitter->addConstraint( mCons4 );

	//Set convergence criteira
	fitter->setMaxNbIter( 30 );
	fitter->setMaxDeltaS( 1e-2 );
	fitter->setMaxF( 1e-1 );
	fitter->setVerbosity(1);

	//Perform the fit
	//  std::cout << "Performing kinematic fit..." << std::endl;
	//print(fitter);

	fitter->fit();
	//std::cout << "Done." << std::endl;
	//  print(fitter);
	//std::cout << "Iteration: " << iteration<<std::endl;

	iteration++;
	float Chi2 = fitter->getS();
	if (Chi2 < minChi2){
	  minChi2 = Chi2;
	  minChi2Iteration = iteration;
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
	delete fitter;
      }
    }while(std::next_permutation(numbers, numbers + size));
  }
  

  //  cout<<"minChi2 "<<minChi2<<" minChi2Iteration: "<< minChi2Iteration<<endl;
  TLorentzVector chi2val; 
  Int_t ncombo=0;

  for (unsigned int b=0; b < size; b++) ncombo += bestCombo[b]*pow(10,size-b -1);


  //cout << ncombo << " ncombo " << endl;
  fitJets[0] = ncombo;
  fitJets[1] = minChi2;
  /*
  chi2val.SetPtEtaPhiE(ncombo,0,0,minChi2);
  fitJets.push_back(chi2val);
  */
  return fitJets;
}
