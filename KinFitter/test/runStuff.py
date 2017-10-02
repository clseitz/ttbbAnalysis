# import ROOT in batch mode
import sys
import ROOT
from array import array
from itertools import combinations
import numpy as np
ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load("libPhysicsToolsKinFitter.so")
ROOT.gROOT.ProcessLine(".L testKinFitFWLite.C+")


from ROOT import testKinFitFWLite,TFile,TTree,TString

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
                  self.MCmatch = array( 'i', [ 0 ] )
                  tree.Branch(name + '_MCmatch',self.MCmatch,name + '_MCmatch/I')

            tree.Branch(name + '_m',self.mass,name + '_m/F')
            tree.Branch(name + '_pt',self.pt,name + '_pt/F')
            tree.Branch(name + '_phi',self.phi,name + '_phi/F')
            tree.Branch(name + '_eta',self.eta,name + '_eta/F')

      def UpdateParam(self,vec):
            self.pt[0] = vec.Pt()
            self.eta[0] = vec.Eta()
            self.phi[0] = vec.Phi()
            self.mass[0] = vec.M()



def JetComb(njets, nbjets =0,bpos_=[]):
      '''return all the combinations for jet and b-jet  positioning'''
      if nbjets == 0: return list(combinations(range(njets),6))
      else:
            bpermut_ = list(combinations(bpos_,2))
#            print bpermut, ' bpermut '
            jetpermut_ = list(combinations(np.delete(np.array(range(njets)),bpos_),4))#permutations of all jets - bjets
            completelist_ = []
            for i in range(len(bpermut_)):
                  for j in range(len(jetpermut_)):
                        completelist_.append(bpermut_[i] + jetpermut_[j])
#            print completelist_, ' complete list '
            return completelist_





def KinAnalysis(tagged,nmaxentries):
      '''Runs the entire code for jet combinationfor the option of only btagging (tagged == True) or full'''
      if tagged: print 'Running tagged option'
      else: print 'Running full option'

      btagbench = 0.8 #benchmark for btag cut
      f = ROOT.TFile.Open("/shome/clseitz/Files/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root", "read")
      t = f.Get("tree")
      if tagged: fout = TFile('Kinematic_Results_Tagged.root','recreate')
      else: fout = TFile('Kinematic_Results_Full.root','recreate')
      
      tkin =  TTree("tkin","tkin")
      if not tagged:
            w1 = Particle('w1',tkin)
            w2 = Particle('w2',tkin)
            top1 = Particle('top1',tkin)
            top2 = Particle('top2',tkin)
            b1 = Particle('b1',tkin,1)
            b2 = Particle('b2',tkin,1)            
            tt = Particle('tt',tkin)
            
            chi2 = array( 'f', [ 0. ] )
            tkin.Branch('chi2',chi2,'chi2/F')



      w1tag = Particle('w1_tag',tkin)
      w2tag = Particle('w2_tag',tkin)
      top1tag = Particle('top1_tag',tkin)
      top2tag = Particle('top2_tag',tkin)
      b1tag = Particle('b1_tag',tkin,1)
      b2tag = Particle('b2_tag',tkin,1)            
      tttag = Particle('tt_tag',tkin)
      ttmc = Particle('ttMC',tkin)
      
      chi2tag = array( 'f', [ 0. ] )
      tkin.Branch('chi2_tag',chi2tag,'chi2_tag/F')
      n_bjets = array( 'I', [ 0 ] )
      tkin.Branch('n_bjets',n_bjets,'n_bjets/I')
      n_jets = array( 'i', [ 0 ] )
      tkin.Branch('n_jets',n_jets,'n_jets/I')
            
      




      
      for e,event in enumerate(t) : 
            bpos_ = []
            jets_ = []
            bestcomb_ = []
            
            if t.njets >= 6:
                  nbjets = 0
                  bestchi2 = 99999
                  
                  for i in range(0,t.njets):#find the b-tagged jets
                        jet = ROOT.TLorentzVector()
                        jet.SetPtEtaPhiM(t.jets_pt[i],t.jets_eta[i],t.jets_phi[i],t.jets_mass[i])
                        jets_.append(jet)
                        if (t.jets_btagCSV[i] > btagbench):
                              bpos_.append(i)
                              nbjets +=1
                  if nbjets < 2: continue
                  
                  if tagged and t.njets - nbjets <4: continue
                  n_bjets[0] = nbjets;
                  n_jets[0] = t.njets;
                  if not tagged:
                        notagcomb_ = JetComb(t.njets)                                    
                        for i in range(len(notagcomb_)):
                              fitJets = testKinFitFWLite(jets_[notagcomb_[i][0]],jets_[notagcomb_[i][1]],jets_[notagcomb_[i][2]],jets_[notagcomb_[i][3]],jets_[notagcomb_[i][4]], jets_[notagcomb_[i][5]])

                              if bestchi2 > fitJets[1]: # fitJets[6].E() is the chi2 of the combination
                                    del bestcomb_[:]
                                    bestchi2 = fitJets[1]
                                    bestpermut_ = [int(d) for d in str(int(fitJets[0]))] # fitJets[6].Pt() is the permutation that results in the chi2 value
                                    if len(bestpermut_) == 5: bestpermut_ = [0] + bestpermut_
                                    for  j in  bestpermut_:
                                          bestcomb_.append(notagcomb_[i][j]) 

#                        print bestcomb_, ' bestcomb from text '
                        chi2[0] = bestchi2
                        b1.UpdateParam(jets_[bestcomb_[0]])
                        b2.UpdateParam(jets_[bestcomb_[1]])
                        w1.UpdateParam(jets_[bestcomb_[2]]+jets_[bestcomb_[3]])
                        w2.UpdateParam(jets_[bestcomb_[5]]+jets_[bestcomb_[4]])
                        top1.UpdateParam(jets_[bestcomb_[0]]+jets_[bestcomb_[3]]+jets_[bestcomb_[2]])
                        top2.UpdateParam(jets_[bestcomb_[1]]+jets_[bestcomb_[4]]+jets_[bestcomb_[5]])
                        tt.UpdateParam(jets_[bestcomb_[0]]+jets_[bestcomb_[1]]+jets_[bestcomb_[2]]+jets_[bestcomb_[3]]+jets_[bestcomb_[4]]+jets_[bestcomb_[5]])
                        
                        b1.MCid[0] =t.jets_mcFlavour[bestcomb_[0]] 
                        b2.MCid[0] = t.jets_mcFlavour[bestcomb_[1]]
                        b1.MCmatch[0] =t.jets_matchBfromHadT[bestcomb_[0]] 
                        b2.MCmatch[0] = t.jets_matchBfromHadT[bestcomb_[1]]

                                    
                  bestchi2 = 999999
#                  print t.njets, ' njets ',nbjets, ' nbjets ', ' bpos ', bpos_ 
                  tagcomb_ = JetComb(t.njets,nbjets,bpos_)
#                  print tagcomb_, ' tagcomb '
                  for i in range(len(tagcomb_)):
                        fitJets = testKinFitFWLite(jets_[tagcomb_[i][0]],jets_[tagcomb_[i][1]],jets_[tagcomb_[i][2]],jets_[tagcomb_[i][3]],jets_[tagcomb_[i][4]], jets_[tagcomb_[i][5]],1)
                        if bestchi2 > fitJets[1]:
                              del bestcomb_[:]
                              bestchi2 = fitJets[1]
                              bestpermut_ = [0]+[int(d) for d in str(int(fitJets[0]))]
#                              print bestpermut_, ' best permut '
                              for  j in  bestpermut_:
                                    bestcomb_.append(tagcomb_[i][j]) 

#                  print bestcomb_, ' bestcomb from tag '
            
                  chi2tag[0] = bestchi2
                  b1tag.UpdateParam(jets_[bestcomb_[0]])
                  b2tag.UpdateParam(jets_[bestcomb_[1]])
                  w1tag.UpdateParam(jets_[bestcomb_[2]]+jets_[bestcomb_[3]])
                  w2tag.UpdateParam(jets_[bestcomb_[5]]+jets_[bestcomb_[4]])
                  top1tag.UpdateParam(jets_[bestcomb_[0]]+jets_[bestcomb_[3]]+jets_[bestcomb_[2]])
                  top2tag.UpdateParam(jets_[bestcomb_[1]]+jets_[bestcomb_[4]]+jets_[bestcomb_[5]])
                  tttag.UpdateParam(jets_[bestcomb_[0]]+jets_[bestcomb_[1]]+jets_[bestcomb_[2]]+jets_[bestcomb_[3]]+jets_[bestcomb_[4]]+jets_[bestcomb_[5]])
                  b1tag.MCid[0] =t.jets_mcFlavour[bestcomb_[0]] 
                  b2tag.MCid[0] = t.jets_mcFlavour[bestcomb_[1]]
                  b1tag.MCmatch[0] =t.jets_matchBfromHadT[bestcomb_[0]] 
                  b2tag.MCmatch[0] = t.jets_matchBfromHadT[bestcomb_[1]]
                  
            
                  if t.ngenTopHad == 2:
                        t1 = ROOT.TLorentzVector()
                        t2 = ROOT.TLorentzVector()
                        t1.SetPtEtaPhiM(t.genTopHad_pt[0],t.genTopHad_eta[0],t.genTopHad_phi[0],t.genTopHad_mass[0])
                        t2.SetPtEtaPhiM(t.genTopHad_pt[1],t.genTopHad_eta[1],t.genTopHad_phi[1],t.genTopHad_mass[1])
                        
                        ttmc.UpdateParam(t1+t2)
                  else:
                        ttmc.mass[0] = -10.0
                        ttmc.pt[0] =-10.0
                        ttmc.phi[0] =-10.0
                        ttmc.eta[0] =-10.0
                        

                  tkin.Fill()

            if e > nmaxentries: break

      fout.Write()
      f.Close()
