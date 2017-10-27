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



std::vector<Double_t > jet_resolutions(Double_t et,Double_t eta,TString partname){
  Double_t res_et,res_eta,res_phi;
  res_et=res_eta=res_phi=0;
  if(partname == "udsc"){

    if( 0.000<=abs(eta) && abs(eta)<0.087){

      res_et = et * (sqrt(pow(0.06,2) + pow((1.0023/sqrt(et)),2) + pow((0/et),2)));
      res_eta  = sqrt(pow(0.00901,2) + pow((1.5284/et),2));
      res_phi  = sqrt(pow(0.0104,2) + pow((1.6004/et),2));
    }
    if( 0.087<=abs(eta) && abs(eta)<0.174){
      res_et = et * (sqrt(pow(0.0633,2) + pow((0.964/sqrt(et)),2) + pow((1.49/et),2)));
      res_eta  = sqrt(pow(0.00927,2) + pow((1.486/et),2));
      res_phi  = sqrt(pow(0.01113,2) + pow((1.5354/et),2));
    }
  
    if( 0.174<=abs(eta) && abs(eta)<0.261){
      res_et = et * (sqrt(pow(0.0595,2) + pow((0.973/sqrt(et)),2) + pow((1.52/et),2)));
      res_eta  = sqrt(pow(0.00958,2) + pow((1.4794/et),2));
      res_phi  = sqrt(pow(0.01093,2) + pow((1.5387/et),2));
    }
    if( 0.261<=abs(eta) && abs(eta)<0.348){
      res_et = et * (sqrt(pow(0.058,2) + pow((0.991/sqrt(et)),2) + pow((1.3/et),2)));
      res_eta  = sqrt(pow(0.00884,2) + pow((1.5269/et),2));
      res_phi  = sqrt(pow(0.01107,2) + pow((1.5398/et),2));
    }
    if( 0.348<=abs(eta) && abs(eta)<0.435){
      res_et = et * (sqrt(pow(0.0597,2) + pow((0.96/sqrt(et)),2) + pow((1.82/et),2)));
      res_eta  = sqrt(pow(0.00915,2) + pow((1.5229/et),2));
      res_phi  = sqrt(pow(0.01079,2) + pow((1.557/et),2));
    }
    if( 0.435<=abs(eta) && abs(eta)<0.522){
      res_et = et * (sqrt(pow(0.0582,2) + pow((0.977/sqrt(et)),2) + pow((1.45/et),2)));
      res_eta  = sqrt(pow(0.00936,2) + pow((1.5322/et),2));
      res_phi  = sqrt(pow(0.01055,2) + pow((1.5636/et),2));
    }
    if( 0.522<=abs(eta) && abs(eta)<0.609){
      res_et = et * (sqrt(pow(0.0603,2) + pow((0.96/sqrt(et)),2) + pow((1.58/et),2)));
      res_eta  = sqrt(pow(0.00959,2) + pow((1.5176/et),2));
      res_phi  = sqrt(pow(0.01042,2) + pow((1.5547/et),2));
    }
    if( 0.609<=abs(eta) && abs(eta)<0.696){
      res_et = et * (sqrt(pow(0.0535,2) + pow((1/sqrt(et)),2) + pow((1.31/et),2)));
      res_eta  = sqrt(pow(0.00971,2) + pow((1.5233/et),2));
      res_phi  = sqrt(pow(0.01043,2) + pow((1.5674/et),2));
    }
    if( 0.696<=abs(eta) && abs(eta)<0.783){
      res_et = et * (sqrt(pow(0.0472,2) + pow((1.039/sqrt(et)),2) + pow((0.67/et),2)));
      res_eta  = sqrt(pow(0.00966,2) + pow((1.5239/et),2));
      res_phi  = sqrt(pow(0.01021,2) + pow((1.5725/et),2));
    }
    if( 0.783<=abs(eta) && abs(eta)<0.870){
      res_et = et * (sqrt(pow(0.0561,2) + pow((1.016/sqrt(et)),2) + pow((1.31/et),2)));
      res_eta  = sqrt(pow(0.00969,2) + pow((1.5407/et),2));
      res_phi  = sqrt(pow(0.00981,2) + pow((1.5962/et),2));
    }
    if( 0.870<=abs(eta) && abs(eta)<0.957){
      res_et = et * (sqrt(pow(0.0543,2) + pow((1.0701/sqrt(et)),2) + pow((0/et),2)));
      res_eta  = sqrt(pow(0.00976,2) + pow((1.5745/et),2));
      res_phi  = sqrt(pow(0.01039,2) + pow((1.6025/et),2));
    }
    if( 0.957<=abs(eta) && abs(eta)<1.044){
      res_et = et * (sqrt(pow(0.0544,2) + pow((1.071/sqrt(et)),2) + pow((1.2/et),2)));
      res_eta  = sqrt(pow(0.01025,2) + pow((1.5794/et),2));
      res_phi  = sqrt(pow(0.01002,2) + pow((1.6162/et),2));
    }
    if( 1.044<=abs(eta) && abs(eta)<1.131){
      res_et = et * (sqrt(pow(0.0537,2) + pow((1.1222/sqrt(et)),2) + pow((0/et),2)));
      res_eta  = sqrt(pow(0.01038,2) + pow((1.5695/et),2));
      res_phi  = sqrt(pow(0.01093,2) + pow((1.6176/et),2));
    }
    if( 1.131<=abs(eta) && abs(eta)<1.218){
      res_et = et* (sqrt(pow(0.0509,2) + pow((1.1539/sqrt(et)),2) + pow((0/et),2)));
      res_eta  = sqrt(pow(0.01099,2) + pow((1.6058/et),2));
      res_phi  = sqrt(pow(0.01018,2) + pow((1.6654/et),2));
    }
    if( 1.218<=abs(eta) && abs(eta)<1.305){
      res_et = et * (sqrt(pow(0.0444,2) + pow((1.2024/sqrt(et)),2) + pow((0/et),2)));
      res_eta  = sqrt(pow(0.01353,2) + pow((1.6026/et),2));
      res_phi  = sqrt(pow(0.0108,2) + pow((1.6873/et),2));
    }
    if( 1.305<=abs(eta) && abs(eta)<1.392){
      res_et = et * (sqrt(pow(0.0426,2) + pow((1.2401/sqrt(et)),2) + pow((0/et),2)));
      res_eta  = sqrt(pow(0.01863,2) + pow((1.5378/et),2));
      res_phi  = sqrt(pow(0.01338,2) + pow((1.7234/et),2));
    }
    if( 1.392<=abs(eta) && abs(eta)<1.479){
      res_et = et * (sqrt(pow(0.0435,2) + pow((1.2544/sqrt(et)),2) + pow((0/et),2)));
      res_eta  = sqrt(pow(0.01234,2) + pow((1.6169/et),2));
      res_phi  = sqrt(pow(0.01361,2) + pow((1.7677/et),2));
    }
    if( 1.479<=abs(eta) && abs(eta)<1.566){
      res_et = et * (sqrt( pow((1.2566/sqrt(et)),2) + pow((0/et),2)));
      res_eta  = sqrt(pow(0.01286,2) + pow((1.6066/et),2));
      res_phi  = sqrt(pow(0.01145,2) + pow((1.7966/et),2));
    }
    if( 1.566<=abs(eta) && abs(eta)<1.653){
      res_et = et * (sqrt( pow((1.1734/sqrt(et)),2) + pow((0/et),2)));
      res_eta  = sqrt(pow(0.01086,2) + pow((1.6108/et),2));
      res_phi  = sqrt(pow(0.01042,2) + pow((1.792/et),2));
    }
    if( 1.653<=abs(eta) && abs(eta)<1.740){
      res_et = et * (sqrt( pow((1.1259/sqrt(et)),2) + pow((0/et),2)));
      res_eta  = sqrt(pow(0.01056,2) + pow((1.626/et),2));
      res_phi  = sqrt(pow(0.00983,2) + pow((1.7587/et),2));
    }
    if( 1.740<=abs(eta) && abs(eta)<1.830){
      res_et = et * (sqrt( pow((1.0982/sqrt(et)),2) + pow((0/et),2)));
      res_eta  = sqrt(pow(0.00922,2) + pow((1.6977/et),2));
      res_phi  = sqrt(pow(0.0094,2) + pow((1.7323/et),2));
    }
    if( 1.830<=abs(eta) && abs(eta)<1.930){
      res_et = et * (sqrt( pow((0.988/sqrt(et)),2) + pow((2.678/et),2)));
      res_eta  = sqrt(pow(0.00982,2) + pow((1.6827/et),2));
      res_phi  = sqrt(pow(0.0085,2) + pow((1.7074/et),2));
    }
    if( 1.930<=abs(eta) && abs(eta)<2.043){
      res_et = et * (sqrt( pow((0.957/sqrt(et)),2) + pow((2.569/et),2)));
      res_eta  = sqrt(pow(0.01029,2) + pow((1.6801/et),2));
      res_phi  = sqrt(pow(0.00834,2) + pow((1.6954/et),2));
    }
    if( 2.043<=abs(eta) && abs(eta)<2.172){
      res_et = et * (sqrt( pow((0.9455/sqrt(et)),2) + pow((2.48/et),2)));
      res_eta  = sqrt(pow(0.01114,2) + pow((1.6469/et),2));
      res_phi  = sqrt(pow(0.0082,2) + pow((1.6705/et),2));
    }
    if( 2.172<=abs(eta) && abs(eta)<2.322){
      res_et = et * (sqrt( pow((0.9015/sqrt(et)),2) + pow((2.75/et),2)));
      res_eta  = sqrt(pow(0.0105,2) + pow((1.6086/et),2));
      res_phi  = sqrt(pow(0.00883,2) + pow((1.6729/et),2));
    }
    if( 2.322<=abs(eta) && abs(eta)<2.500){
      res_et = et * (sqrt( pow((0.9007/sqrt(et)),2) + pow((3.059/et),2)));
      res_eta  = sqrt(pow(0.01117,2) + pow((1.926/et),2));
      res_phi  = sqrt(pow(0.01045,2) + pow((1.7223/et),2));
    }
  }

  if(partname == "b"){

    if(0.000<=abs(eta) && abs(eta)<0.087){
      res_et= et * (sqrt(pow(0.0849,2) + pow((0.855/sqrt(et)),2) + pow((3.43/et),2)));
      res_eta= sqrt(pow(0.00672,2) + pow((1.5978/et),2));
      res_phi  = sqrt(pow(0.00842,2) + pow((1.6991/et),2));
    }
    
    if(0.087<=abs(eta) && abs(eta)<0.174){
      res_et= et * (sqrt(pow(0.08,2) + pow((0.959/sqrt(et)),2) + pow((2.15/et),2)));
      res_eta= sqrt(pow(0.00595,2) + pow((1.6273/et),2));
      res_phi  = sqrt(pow(0.00802,2) + pow((1.7211/et),2));
    }
    
    if(0.174<=abs(eta) && abs(eta)<0.261){
      res_et= et * (sqrt(pow(0.0734,2) + pow((0.998/sqrt(et)),2) + pow((1.67/et),2)));
      res_eta= sqrt(pow(0.00656,2) + pow((1.6177/et),2));
      res_phi  = sqrt(pow(0.00789,2) + pow((1.7235/et),2));
    }
    
    if(0.261<=abs(eta) && abs(eta)<0.348){
      res_et= et * (sqrt(pow(0.0786,2) + pow((0.935/sqrt(et)),2) + pow((2.23/et),2)));
      res_eta= sqrt(pow(0.00618,2) + pow((1.6293/et),2));
      res_phi  = sqrt(pow(0.00779,2) + pow((1.7145/et),2));
    }
    
    if(0.348<=abs(eta) && abs(eta)<0.435){
      res_et= et * (sqrt(pow(0.0721,2) + pow((1.001/sqrt(et)),2) + pow((1.42/et),2)));
      res_eta= sqrt(pow(0.00623,2) + pow((1.6427/et),2));
      res_phi  = sqrt(pow(0.00806,2) + pow((1.7107/et),2));
    }
    
    if(0.435<=abs(eta) && abs(eta)<0.522){
      res_et= et * (sqrt(pow(0.0682,2) + pow((1.011/sqrt(et)),2) + pow((1.37/et),2)));
      res_eta= sqrt(pow(0.00678,2) + pow((1.657/et),2));
      res_phi  = sqrt(pow(0.00806,2) + pow((1.7229/et),2));
    }
    
    if(0.522<=abs(eta) && abs(eta)<0.609){
      res_et= et * (sqrt(pow(0.0637,2) + pow((1.037/sqrt(et)),2) + pow((1.24/et),2)));
      res_eta= sqrt(pow(0.00633,2) + pow((1.6528/et),2));
      res_phi  = sqrt(pow(0.00785,2) + pow((1.722/et),2));
    }
    
    if(0.609<=abs(eta) && abs(eta)<0.696){
      res_et= et * (sqrt(pow(0.0658,2) + pow((1.032/sqrt(et)),2) + pow((0.83/et),2)));
      res_eta= sqrt(pow(0.00684,2) + pow((1.6606/et),2));
      res_phi  = sqrt(pow(0.00777,2) + pow((1.7348/et),2));
    }
    
    if(0.696<=abs(eta) && abs(eta)<0.783){
      res_et= et * (sqrt(pow(0.0661,2) + pow((1.0633/sqrt(et)),2) + pow((0/et),2)));
      res_eta= sqrt(pow(0.00722,2) + pow((1.6615/et),2));
      res_phi  = sqrt(pow(0.00786,2) + pow((1.7394/et),2));
    }
    
    if(0.783<=abs(eta) && abs(eta)<0.870){
      res_et= et * (sqrt(pow(0.0649,2) + pow((1.0755/sqrt(et)),2) + pow((0/et),2)));
      res_eta= sqrt(pow(0.00803,2) + pow((1.655/et),2));
      res_phi  = sqrt(pow(0.0077,2) + pow((1.7591/et),2));
    }
    
    if(0.870<=abs(eta) && abs(eta)<0.957){
      res_et= et * (sqrt(pow(0.0731,2) + pow((1.054/sqrt(et)),2) + pow((0.6/et),2)));
      res_eta= sqrt(pow(0.00817,2) + pow((1.678/et),2));
      res_phi  = sqrt(pow(0.00807,2) + pow((1.7585/et),2));
    }
    
    if(0.957<=abs(eta) && abs(eta)<1.044){
      res_et= et * (sqrt(pow(0.068,2) + pow((1.0925/sqrt(et)),2) + pow((0/et),2)));
      res_eta= sqrt(pow(0.0085,2) + pow((1.6774/et),2));
      res_phi  = sqrt(pow(0.00805,2) + pow((1.7778/et),2));
    }
    
    if(1.044<=abs(eta) && abs(eta)<1.131){
      res_et= et * (sqrt(pow(0.0662,2) + pow((1.1339/sqrt(et)),2) + pow((0/et),2)));
      res_eta= sqrt(pow(0.00858,2) + pow((1.7079/et),2));
      res_phi  = sqrt(pow(0.00852,2) + pow((1.7953/et),2));
    }
    
    if(1.131<=abs(eta) && abs(eta)<1.218){
      res_et= et * (sqrt(pow(0.064,2) + pow((1.1553/sqrt(et)),2) + pow((0/et),2)));
      res_eta= sqrt(pow(0.00927,2) + pow((1.7597/et),2));
      res_phi  = sqrt(pow(0.00856,2) + pow((1.8331/et),2));
    }
    
    if(1.218<=abs(eta) && abs(eta)<1.305){
      res_et= et * (sqrt(pow(0.0692,2) + pow((1.1655/sqrt(et)),2) + pow((0/et),2)));
      res_eta= sqrt(pow(0.01203,2) + pow((1.7792/et),2));
      res_phi  = sqrt(pow(0.00904,2) + pow((1.8867/et),2));
    }
    
    if(1.305<=abs(eta) && abs(eta)<1.392){
      res_et= et * (sqrt(pow(0.0756,2) + pow((1.1773/sqrt(et)),2) + pow((0/et),2)));
      res_eta= sqrt(pow(0.01479,2) + pow((1.774/et),2));
      res_phi  = sqrt(pow(0.0108,2) + pow((1.9241/et),2));
    }
    
    if(1.392<=abs(eta) && abs(eta)<1.479){
      res_et= et * (sqrt(pow(0.0761,2) + pow((1.1932/sqrt(et)),2) + pow((0/et),2)));
      res_eta= sqrt(pow(0.01191,2) + pow((1.791/et),2));
      res_phi  = sqrt(pow(0.011,2) + pow((2.006/et),2));
    }
    
    if(1.479<=abs(eta) && abs(eta)<1.566){
      res_et= et * (sqrt(pow(0.0631,2) + pow((1.2178/sqrt(et)),2) + pow((0/et),2)));
      res_eta= sqrt(pow(0.01287,2) + pow((1.764/et),2));
      res_phi  = sqrt(pow(0.01025,2) + pow((1.998/et),2));
    }
    
    if(1.566<=abs(eta) && abs(eta)<1.653){
      res_et= et * (sqrt(pow(0.0429,2) + pow((1.2103/sqrt(et)),2) + pow((0/et),2)));
      res_eta= sqrt(pow(0.01092,2) + pow((1.789/et),2));
      res_phi  = sqrt(pow(0.00904,2) + pow((1.99/et),2));
    }
    
    if(1.653<=abs(eta) && abs(eta)<1.740){
      res_et= et * (sqrt(pow(0,2) + pow((1.2206/sqrt(et)),2) + pow((0/et),2)));
      res_eta= sqrt(pow(0.00985,2) + pow((1.784/et),2));
      res_phi  = sqrt(pow(0.0083,2) + pow((1.954/et),2));
    }
    
    if(1.740<=abs(eta) && abs(eta)<1.830){
      res_et= et * (sqrt(pow(0,2) + pow((1.1902/sqrt(et)),2) + pow((0/et),2)));
      res_eta= sqrt(pow(0.00916,2) + pow((1.881/et),2));
      res_phi  = sqrt(pow(0.00728,2) + pow((1.937/et),2));
    }
    
    if(1.830<=abs(eta) && abs(eta)<1.930){
      res_et= et * (sqrt(pow(0,2) + pow((1.1441/sqrt(et)),2) + pow((0/et),2)));
      res_eta= sqrt(pow(0.01144,2) + pow((1.811/et),2));
      res_phi  = sqrt(pow(0.00664,2) + pow((1.9/et),2));
    }
    
    if(1.930<=abs(eta) && abs(eta)<2.043){
      res_et= et * (sqrt(pow(0,2) + pow((1.1221/sqrt(et)),2) + pow((0/et),2)));
      res_eta= sqrt(pow(0.01113,2) + pow((1.864/et),2));
      res_phi  = sqrt(pow(0.00618,2) + pow((1.857/et),2));
    }
    
    if(2.043<=abs(eta) && abs(eta)<2.172){
      res_et= et * (sqrt(pow(0,2) + pow((1.0843/sqrt(et)),2) + pow((1.73/et),2)));
      res_eta= sqrt(pow(0.01176,2) + pow((1.844/et),2));
      res_phi  = sqrt(pow(0.00624,2) + pow((1.855/et),2));
    }
    
    if(2.172<=abs(eta) && abs(eta)<2.322){
      res_et= et * (sqrt(pow(0,2) + pow((1.0579/sqrt(et)),2) + pow((1.78/et),2)));
      res_eta= sqrt(pow(0.01076,2) + pow((1.821/et),2));
      res_phi  = sqrt(pow(0.00543,2) + pow((1.884/et),2));
    }
   
    if(2.322<=abs(eta) && abs(eta)<2.500){
      res_et= et * (sqrt(pow(0,2) + pow((1.1037/sqrt(et)),2) + pow((1.62/et),2)));
      res_eta= sqrt(pow(0.00883,2) + pow((2.189/et),2));
      res_phi  = sqrt(pow(0.00836,2) + pow((1.959/et),2));
    }



  }
  std::vector<Double_t> resolutions;
  resolutions.push_back(res_et);
  resolutions.push_back(res_eta);
  resolutions.push_back(res_phi);
  return resolutions;
  
}



TLorentzVector FittedValues(TFitParticleEtEtaPhi * part){
  TLorentzVector update_part;
  update_part.SetPtEtaPhiE(part->getCurr4Vec()->Pt(),part->getCurr4Vec()->Eta(),part->getCurr4Vec()->Phi(),part->getCurr4Vec()->E());
  return update_part;

}



Double_t*  KinFitFWLite( std::vector<TLorentzVector >  &vJets, const  int constrained)

{

  std::vector<TLorentzVector >  vFit;
  const unsigned short size = 6;
  unsigned int numJetForFit  = 6;
  int converged = 0;
  float minChi2 = 999;
  Double_t probChi2 = 0;

 


  TMatrixD m1(3,3); m1.Zero();
  TMatrixD m2(3,3); m2.Zero();
  TMatrixD m3(3,3); m3.Zero();
  TMatrixD m4(3,3); m4.Zero();
  TMatrixD m5(3,3); m5.Zero();
  TMatrixD m6(3,3); m6.Zero();
  
 	
  std::vector<std::vector<Double_t> > res;

  for(int i = 0; i < size; i++){
    TString partname = "b";
    if(i > 2) partname = "udsc";
    res.push_back(jet_resolutions(vJets[i].Et(), vJets[i].Eta(),partname));
  }


  for(int i =0; i<3;i++){
    m1(i,i) = res[0][i];
    m2(i,i) = res[1][i];
    m3(i,i) = res[2][i];
    m4(i,i) = res[3][i];
    m5(i,i) = res[4][i];
    m6(i,i) = res[5][i];

  }


  TFitParticleEtEtaPhi *jetB = new TFitParticleEtEtaPhi( "jetB", "jetB", &vJets[0], &m1 );
  TFitParticleEtEtaPhi *jetBbar = new TFitParticleEtEtaPhi( "jetBbar", "jetBbar", &vJets[1], &m2 );
  TFitParticleEtEtaPhi *jetWq = new TFitParticleEtEtaPhi( "jetWq", "jetWq", &vJets[2], &m3 );
  TFitParticleEtEtaPhi *jetWqbar = new TFitParticleEtEtaPhi( "jetWqbar", "jetWqbar", &vJets[3], &m4 );
  TFitParticleEtEtaPhi *jetWp = new TFitParticleEtEtaPhi( "jetWp", "jetWp", &vJets[4], &m5 );
  TFitParticleEtEtaPhi *jetWpbar = new TFitParticleEtEtaPhi( "jetWpbar", "jetWpbar", &vJets[5], &m6 );


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
	
  minChi2 = fitter->getS();

	

  probChi2 = TMath::Prob(fitter->getS(),fitter->getNDF());
  vFit.push_back(FittedValues(jetB));
  vFit.push_back(FittedValues(jetBbar));
  vFit.push_back(FittedValues(jetWq));
  vFit.push_back(FittedValues(jetWqbar));
  vFit.push_back(FittedValues(jetWp));
  vFit.push_back(FittedValues(jetWpbar));
  
	  
	
	
  
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
    
  
  static Double_t fitJets[11];

  fitJets[0] = converged;
  fitJets[1] = minChi2;
  fitJets[2] = probChi2;
  fitJets[3] = (vFit[0] + vFit[2] + vFit[3]).Pt();
  fitJets[4] = (vFit[0] + vFit[2] + vFit[3]).Eta();
  fitJets[5] = (vFit[0] + vFit[2] + vFit[3]).Phi();
  fitJets[6] = (vFit[0] + vFit[2] + vFit[3]).M();
  fitJets[7] = (vFit[1] + vFit[4] + vFit[5]).Pt();
  fitJets[8] = (vFit[1] + vFit[4] + vFit[5]).Eta();
  fitJets[9] = (vFit[1] + vFit[4] + vFit[5]).Phi();
  fitJets[10] = (vFit[1] + vFit[4] + vFit[5]).M();

  return fitJets;
}
