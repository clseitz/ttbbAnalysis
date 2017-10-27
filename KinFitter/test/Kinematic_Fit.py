# import ROOT in batch mode
import sys
import ROOT
from array import array
from itertools import permutations, combinations
import numpy as np
import re
ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load("libPhysicsToolsKinFitter.so")
ROOT.gROOT.ProcessLine(".L Kinfit.C+")

from ROOT import KinFitFWLite,TFile,TTree,TString

class Particle:
      'Defines a particle with its properties and create a TTree branch associated'
      def __init__(self,name,tree,mc = 0):
            self.pt = array( 'f', [ 0. ] )
            self.mass = array( 'f', [ 0. ] )
            self.eta = array( 'f', [ 0. ] )
            self.phi = array( 'f', [ 0. ] )
            if mc:
                  self.MCid = array( 'i', [ 0 ] )
                  tree.Branch(name + '_MCid',self.MCid,name + '_MCid/I')
                  #self.BMCmatch = array( 'i', [ 0 ] )
                  #tree.Branch(name + '_BMCmatch',self.BMCmatch,name + '_BMCmatch/I')
                  self.TopMCmatch = array( 'i', [ 0 ] )
                  tree.Branch(name + '_TopMCmatch',self.TopMCmatch,name + '_TopMCmatch/I')
                  
            tree.Branch(name + '_m',self.mass,name + '_m/F')
            tree.Branch(name + '_pt',self.pt,name + '_pt/F')
            tree.Branch(name + '_phi',self.phi,name + '_phi/F')
            tree.Branch(name + '_eta',self.eta,name + '_eta/F')

      def UpdateParam(self,vec,MCid=0,MCmatch=0):
            self.pt[0] = vec.Pt()
            self.eta[0] = vec.Eta()
            self.phi[0] = vec.Phi()
            self.mass[0] = vec.M()
            if hasattr(self, 'MCid'):
                  self.MCid[0] = MCid
                  self.TopMCmatch[0] = MCmatch

def JetComb(njets, nbjets =0,bpos_=[]):
      '''return all the combinations for jet and b-jet  positioning'''
      if nbjets == 0: return list(combinations(range(njets),6))
      else:
            bpermut_ = list(combinations(bpos_,2))

#            print bpermut_, ' bpermut '
            jetpermut_ = list(permutations((np.delete(np.array(range(njets)),bpos_)),4))#permutations of all jets - bjets
            completelist_ = []
            for i in range(len(bpermut_)):
                  for j in range(len(jetpermut_)):
                        if jetpermut_[j][0] > jetpermut_[j][1] or jetpermut_[j][2] > jetpermut_[j][3]: continue
                        completelist_.append(bpermut_[i] + jetpermut_[j])
#            print completelist_, ' complete list ', nbjets, ' nbjets ', njets, ' njets'
            return completelist_





#def KinAnalysis(tagged,filename,nmaxentries = -1):
if __name__ == "__main__":
 
      if len(sys.argv) > 2:
            filename = sys.argv[1]
            constrained = int(sys.argv[2])
            nmaxentries = -1
      if len(sys.argv) > 3:
            nmaxentries = int(sys.argv[3])
      

      t = ROOT.TChain('vhbb/tree')
#      t = ROOT.TChain('tree')



      hCount =ROOT.TH1F()
      hCount.SetName('hCount')
      with open(filename) as f:
            lines = f.read().splitlines()
            foutnumber = map(int, re.findall(r'\d+', lines[0]))[-1]           
            for line in lines:
                  t.Add(line)
                  f2 = ROOT.TFile.Open(line,'READ')
                  if f2:
                        hist = f2.Get('vhbb/Count')
                        if hist:
                              hCount.Add(hist)
                              
                        f2.Close()

                  
#      maxcomb = 210
      topjets = 6
#      btagbench = 0.8 #benchmark for btag cut
      foutname = 'Kinematic_Results_2011'
      if constrained:foutname += '_Constrained'
      foutname += str(foutnumber)+'.root'
      fout = TFile(foutname,'recreate')
      
      tkin =  TTree("tkin","tkin")
      if constrained == 1: print "Using top mass constraint"
      w1 = Particle('w1',tkin)
      w2 = Particle('w2',tkin)
      top1 = Particle('top1',tkin)
      top2 = Particle('top2',tkin)
      b1 = Particle('b1',tkin,1)
      b2 = Particle('b2',tkin,1)            
      tt = Particle('tt',tkin)
      
#      lightq1 = Particle('lightq1',tkin,1)
#      lightqb1 = Particle('lightqb1',tkin,1)
#      lightq2 = Particle('lightq2',tkin,1)
#      lightqb2 = Particle('lightqb2',tkin,1)
      

      corr_top1 = Particle('correct_top1',tkin)
      corr_top2 = Particle('correct_top2',tkin)
      corr_fitted_top1 = Particle('correct_fitted_top1',tkin)
      corr_fitted_top2 = Particle('correct_fitted_top2',tkin)
      fitted_top1 = Particle('fitted_top1',tkin)
      fitted_top2 = Particle('fitted_top2',tkin)
      fitted_tt = Particle('fitted_tt',tkin)
      corr_tt = Particle('correct_tt',tkin)
      corr_fitted_tt = Particle('correct_fitted_tt',tkin)
      ttmc = Particle('ttMC',tkin)
     
      chi2 = array( 'f', [ 0. ] )
      tkin.Branch('chi2',chi2,'chi2/F')
#      all_chi2 = array( 'f', maxcomb*[ 0. ] )
#      tkin.Branch('all_chi2',all_chi2,'all_chi2[n_combinations]/F')
#      combmasses = array( 'f', maxcomb*[ 0. ] )
#      tkin.Branch('comb_masses',combmasses,'comb_masses[n_combinations]/F')
      correct_chi2 = array( 'f', [ 0. ] )
      tkin.Branch('correct_chi2',correct_chi2,'correct_chi2/F')
      prob_chi2 = array( 'f', [ 0. ] )
      tkin.Branch('prob_chi2',prob_chi2,'prob_chi2/F')
#      all_prob_chi2 = array( 'f', maxcomb*[ 0. ] )
#      tkin.Branch('all_prob_chi2',all_prob_chi2,'all_prob_chi2[n_combinations]/F')
      correct_prob_chi2 = array( 'f', [ 0. ] )
      tkin.Branch('correct_prob_chi2',correct_prob_chi2,'correct_prob_chi2/F')
      n_bjets = array( 'I', [ 0 ] )
      tkin.Branch('n_bjets',n_bjets,'n_bjets/I')
      n_jets = array( 'i', [ 0 ] )
      tkin.Branch('n_jets',n_jets,'n_jets/I')
      n_topjets = array( 'i', [ 0 ] )
      tkin.Branch('n_topjets',n_topjets,'n_topjets/I')
      noWrong = array( 'i', [ 0 ] )
      tkin.Branch('no_Wrong',noWrong,'no_Wrong/I')
      noCorrect = array( 'i', [ 0 ] )
      tkin.Branch('no_Correct',noCorrect,'no_Correct/I')
      notConverged = array( 'i', [ 0 ] )
      tkin.Branch('not_Converged',notConverged,'not_Converged/I')
      n_sumIDtop = array( 'i', topjets*[ 0 ] )
      tkin.Branch('n_sumIDtop',n_sumIDtop,'n_sumIDtop[6]/I')
#      n_sumIDtop2 = array( 'i', topjets*[ 0 ] )
#      tkin.Branch('n_sumIDtop2',n_sumIDtop2,'n_sumIDtop2[n_jets]/I')
            
      




      
      for e,event in enumerate(t) :
#            print "processing from event number: ", nevtstart

            bpos_ = []
            jets_ = ROOT.vector('TLorentzVector')()
            ordered_jets_ = ROOT.vector('TLorentzVector')()
            bestcomb_ = []
            corrbestcomb_ = []

#            if t.njets >= 6 and t.n_bjets >=2 and t.njets - t.n_bjets >=4:
            if t.nJet >= 6 and t.nGenLepFromTop == 0 and t.nGenTop == 2:
                  ntopjets = 0
#                  if t.Jet_pt[5]<30 or not (t.HLT_BIT_HLT_PFHT450_SixJet40_BTagCSV_p056_v or t.HLT_BIT_HLT_PFHT400_SixJet30_DoubleBTagCSV_p056_v): continue
                  if t.Jet_pt[t.nJet-1]<40  or not (t.HLT_BIT_HLT_PFHT450_SixJet40_BTagCSV_p056_v or t.HLT_BIT_HLT_PFHT400_SixJet30_DoubleBTagCSV_p056_v): continue
                  nbjets = 0
                  etamin =0
                  for i in range(0,t.nJet):#find the b-tagged jets
                        jet = ROOT.TLorentzVector()
                        jet.SetPtEtaPhiM(t.Jet_pt[i],t.Jet_eta[i],t.Jet_phi[i],t.Jet_mass[i])
                        #jet.SetPtEtaPhiM(t.jets_pt[i],t.jets_eta[i],t.jets_phi[i],t.jets_mass[i])
                        jets_.push_back(jet)
                        if abs(t.Jet_eta[i]) > 2.4: etamin =1
                        if abs(t.Jet_mcMatchId[i]) == 6: ntopjets += 1
                        
                        if t.Jet_btagCSV[t.Jet_btagIdx[i]] > 0.679 :
                              bpos_.append(t.Jet_btagIdx[i])
                              nbjets +=1
                  
                  if nbjets < 2 or t.nJet - nbjets <4 or etamin > 0:continue
                  bestchi2 = corrbestchi2 =  999999999
                  probchi2 = corrprobchi2 = -10
                  n_bjets[0] = nbjets;
                  n_jets[0] = t.nJet;
                  n_topjets[0] = ntopjets
                  noWrong[0] = noCorrect[0] = 0

                  fit_top2 = ROOT.TLorentzVector()
                  fit_top1 = ROOT.TLorentzVector()
                  corr_fit_top2 = ROOT.TLorentzVector()
                  corr_fit_top1 = ROOT.TLorentzVector()


#                  print t.njets, ' njets ',nbjets, ' nbjets ', ' bpos ', bpos_ 
                  comb_ = JetComb(t.nJet,nbjets,bpos_)
                  
                  for i in range(len(comb_)):
                        ncorrcomb =0
                        for j in comb_[i]:ordered_jets_.push_back(jets_[j])
                        fitJets = KinFitFWLite(ordered_jets_,constrained)
                        ordered_jets_.clear()
                        
                        notConverged[0] = int(fitJets[0])
                        
                        ncorr =0
                        sumtop1 = sumtop2=0
                              
                        sumtop1 = t.Jet_mcMatchId[comb_[i][0]] + t.Jet_mcMatchId[comb_[i][2]] + t.Jet_mcMatchId[comb_[i][3]]
                        sumtop2 = t.Jet_mcMatchId[comb_[i][1]] + t.Jet_mcMatchId[comb_[i][4]] + t.Jet_mcMatchId[comb_[i][5]]
                        if 6 == ntopjets and abs(sumtop1) == 18 and abs(sumtop2) == 18:ncorrcomb+=1

                       
                        if (bestchi2 > fitJets[1]) and ncorrcomb == 0 : 
                              bestchi2 = fitJets[1]
                              probchi2 = fitJets[2]
#                              print bestpermut_, ' best permut '
                              fit_top1.SetPtEtaPhiM(fitJets[3],fitJets[4],fitJets[5],fitJets[6])
                              fit_top2.SetPtEtaPhiM(fitJets[7],fitJets[8],fitJets[9],fitJets[10])
                              
                              bestcomb_=comb_[i]
                              n_sumIDtop[0] = t.Jet_mcMatchId[comb_[i][0]]
                              n_sumIDtop[1] = t.Jet_mcMatchId[comb_[i][2]]
                              n_sumIDtop[2] = t.Jet_mcMatchId[comb_[i][3]]
                              n_sumIDtop[3] = t.Jet_mcMatchId[comb_[i][1]]
                              n_sumIDtop[4] = t.Jet_mcMatchId[comb_[i][4]]
                              n_sumIDtop[5] = t.Jet_mcMatchId[comb_[i][5]]
                              
                        if corrbestchi2 > fitJets[1] and ncorrcomb == 1:

                              corrbestchi2 = fitJets[1]
                              corrprobchi2 = fitJets[2]
                              corr_fit_top1.SetPtEtaPhiM(fitJets[3],fitJets[4],fitJets[5],fitJets[6])
                              corr_fit_top2.SetPtEtaPhiM(fitJets[7],fitJets[8],fitJets[9],fitJets[10])
#                              print bestpermut_, ' best permut '

                              corrbestcomb_=comb_[i]
                                    

                  
                  if len(corrbestcomb_) == 0: #There were no correct combinations...
                        noCorrect[0] = 1
                        corrbestchi2 = -10
                        corrprobchi2 = -10
                        corrbestcomb_ = bestcomb_
                        corr_fit_top1 = fit_top1
                        corr_fit_top2 = fit_top2
                        
#                  if (len(corrbestcomb_) == 0) or (len(bestcomb_) == 0): continue
#                  print bestcomb_, ' best comb ', corrbestcomb_, ' best corrcomb \n'

                  corr_fitted_top1.UpdateParam(corr_fit_top1)
                  corr_fitted_top2.UpdateParam(corr_fit_top2)
                  corr_fitted_tt.UpdateParam(corr_fit_top1+corr_fit_top2)
                  fitted_top1.UpdateParam(fit_top1)
                  fitted_top2.UpdateParam(fit_top2)            
                  fitted_tt.UpdateParam(fit_top1+fit_top2)
                  
                  chi2[0] = bestchi2
                  prob_chi2[0] = probchi2
                  correct_chi2[0] = corrbestchi2
                  correct_prob_chi2[0] = corrprobchi2

#                  lightq1.UpdateParam(jets_[bestcomb_[2]])
#                  lightqb1.UpdateParam(jets_[bestcomb_[3]])
#                  lightq2.UpdateParam(jets_[bestcomb_[4]])
#                  lightqb2.UpdateParam(jets_[bestcomb_[5]])
#                  lightq1.TopMCmatch[0] =t.jets_mcMatchId[bestcomb_[2]] 
#                  lightqb1.TopMCmatch[0] =t.jets_mcMatchId[bestcomb_[3]]
#                  lightq2.TopMCmatch[0] =t.jets_mcMatchId[bestcomb_[4]] 
#                  lightqb2.TopMCmatch[0] =t.jets_mcMatchId[bestcomb_[5]] 
#                  lightq1.MCid[0] =t.jets_mcFlavour[bestcomb_[2]] 
#                  lightqb1.MCid[0] =t.jets_mcFlavour[bestcomb_[3]]
#                  lightq2.MCid[0] =t.jets_mcFlavour[bestcomb_[4]] 
#                  lightqb2.MCid[0] =t.jets_mcFlavour[bestcomb_[5]] 
                  
                  b1.UpdateParam(jets_[bestcomb_[0]],t.Jet_mcFlavour[bestcomb_[0]],t.Jet_mcMatchId[bestcomb_[0]])
                  b2.UpdateParam(jets_[bestcomb_[1]],t.Jet_mcFlavour[bestcomb_[1]],t.Jet_mcMatchId[bestcomb_[1]])
                  w1.UpdateParam(jets_[bestcomb_[2]]+jets_[bestcomb_[3]])
                  w2.UpdateParam(jets_[bestcomb_[5]]+jets_[bestcomb_[4]])
                  top1.UpdateParam(jets_[bestcomb_[0]]+jets_[bestcomb_[3]]+jets_[bestcomb_[2]])
                  top2.UpdateParam(jets_[bestcomb_[1]]+jets_[bestcomb_[4]]+jets_[bestcomb_[5]])
                  tt.UpdateParam(jets_[bestcomb_[0]]+jets_[bestcomb_[1]]+jets_[bestcomb_[2]]+jets_[bestcomb_[3]]+jets_[bestcomb_[4]]+jets_[bestcomb_[5]])

                  
                  corr_top1.UpdateParam(jets_[corrbestcomb_[0]]+jets_[corrbestcomb_[3]]+jets_[corrbestcomb_[2]])
                  corr_top2.UpdateParam(jets_[corrbestcomb_[1]]+jets_[corrbestcomb_[4]]+jets_[corrbestcomb_[5]])
                  corr_tt.UpdateParam(jets_[corrbestcomb_[0]]+jets_[corrbestcomb_[1]]+jets_[corrbestcomb_[2]]+jets_[corrbestcomb_[3]]+jets_[corrbestcomb_[4]]+jets_[corrbestcomb_[5]])

                  t1 = ROOT.TLorentzVector()
                  t2 = ROOT.TLorentzVector()
#                        t1.SetPtEtaPhiM(t.genTopHad_pt[0],t.genTopHad_eta[0],t.genTopHad_phi[0],t.genTopHad_mass[0])
#                        t2.SetPtEtaPhiM(t.genTopHad_pt[1],t.genTopHad_eta[1],t.genTopHad_phi[1],t.genTopHad_mass[1])

                  t1.SetPtEtaPhiM(t.GenTop_pt[0],t.GenTop_eta[0],t.GenTop_phi[0],t.GenTop_mass[0])
                  t2.SetPtEtaPhiM(t.GenTop_pt[1],t.GenTop_eta[1],t.GenTop_phi[1],t.GenTop_mass[1])
                  
                  ttmc.UpdateParam(t1+t2)
                  
                        

                  tkin.Fill()
                  if nmaxentries > 0:
                        if e > nmaxentries: break


      hCount.Write()
      fout.Write()
      fout.Close()
#      f.Close()
