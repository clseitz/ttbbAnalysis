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
            w1m = array( 'f', [ 0. ] )
            w1pt = array( 'f', [ 0. ] )
            w1phi = array( 'f', [ 0. ] )
            w1eta = array( 'f', [ 0. ] )
            tkin.Branch('w1m',w1m,'w1m/F')
            tkin.Branch('w1pt',w1pt,'w1pt/F')
            tkin.Branch('w1phi',w1phi,'w1phi/F')
            tkin.Branch('w1eta',w1eta,'w1eta/F')
      
      
            w2m = array( 'f', [ 0. ] )
            w2pt = array( 'f', [ 0. ] )
            w2phi = array( 'f', [ 0. ] )
            w2eta = array( 'f', [ 0. ] )
            tkin.Branch('w2m',w2m,'w2m/F')
            tkin.Branch('w2pt',w2pt,'w2pt/F')
            tkin.Branch('w2phi',w2phi,'w2phi/F')
            tkin.Branch('w2eta',w2eta,'w2eta/F')


            top1m = array( 'f', [ 0. ] )
            top1pt = array( 'f', [ 0. ] )
            top1phi = array( 'f', [ 0. ] )
            top1eta = array( 'f', [ 0. ] )
            tkin.Branch('top1m',top1m,'top1m/F')
            tkin.Branch('top1pt',top1pt,'top1pt/F')
            tkin.Branch('top1phi',top1phi,'top1phi/F')
            tkin.Branch('top1eta',top1eta,'top1eta/F')


            top2m = array( 'f', [ 0. ] )
            top2pt = array( 'f', [ 0. ] )
            top2phi = array( 'f', [ 0. ] )
            top2eta = array( 'f', [ 0. ] )
            tkin.Branch('top2m',top2m,'top2m/F')
            tkin.Branch('top2pt',top2pt,'top2pt/F')
            tkin.Branch('top2phi',top2phi,'top2phi/F')
            tkin.Branch('top2eta',top2eta,'top2eta/F')


            b1m = array( 'f', [ 0. ] )
            b1pt = array( 'f', [ 0. ] )
            b1phi = array( 'f', [ 0. ] )
            b1eta = array( 'f', [ 0. ] )
            tkin.Branch('b1m',b1m,'b1m/F')
            tkin.Branch('b1pt',b1pt,'b1pt/F')
            tkin.Branch('b1phi',b1phi,'b1phi/F')
            tkin.Branch('b1eta',b1eta,'b1eta/F')
            
            
            b2m = array( 'f', [ 0. ] )
            b2pt = array( 'f', [ 0. ] )
            b2phi = array( 'f', [ 0. ] )
            b2eta = array( 'f', [ 0. ] )
            tkin.Branch('b2m',b2m,'b2m/F')
            tkin.Branch('b2pt',b2pt,'b2pt/F')
            tkin.Branch('b2phi',b2phi,'b2phi/F')
            tkin.Branch('b2eta',b2eta,'b2eta/F')
            


            ttm = array( 'f', [ 0. ] )
            ttpt = array( 'f', [ 0. ] )
            ttphi = array( 'f', [ 0. ] )
            tteta = array( 'f', [ 0. ] )
            chi2 = array( 'f', [ 0. ] )
            tkin.Branch('ttm',ttm,'ttm/F')
            tkin.Branch('ttpt',ttpt,'ttpt/F')
            tkin.Branch('ttphi',ttphi,'ttphi/F')
            tkin.Branch('tteta',tteta,'tteta/F')
            tkin.Branch('chi2',chi2,'chi2/F')

            b1mcid = array( 'i', [ 0 ] )
            b2mcid = array( 'i', [ 0 ] )
      
            b1mcid[0] = 0
            b2mcid[0] = 0
      
            tkin.Branch('b1mcid',b1mcid,'b1mcid/I')
            tkin.Branch('b2mcid',b2mcid,'b2mcid/I')



      w1mtag = array( 'f', [ 0. ] )
      w1pttag = array( 'f', [ 0. ] )
      w1phitag = array( 'f', [ 0. ] )
      w1etatag = array( 'f', [ 0. ] )
      tkin.Branch('w1mtag',w1mtag,'w1mtag/F')
      tkin.Branch('w1pttag',w1pttag,'w1pttag/F')
      tkin.Branch('w1phitag',w1phitag,'w1phitag/F')
      tkin.Branch('w1etatag',w1etatag,'w1etatag/F')
      
      
      w2mtag = array( 'f', [ 0. ] )
      w2pttag = array( 'f', [ 0. ] )
      w2phitag = array( 'f', [ 0. ] )
      w2etatag = array( 'f', [ 0. ] )
      tkin.Branch('w2mtag',w2mtag,'w2mtag/F')
      tkin.Branch('w2pttag',w2pttag,'w2pttag/F')
      tkin.Branch('w2phitag',w2phitag,'w2phitag/F')
      tkin.Branch('w2etatag',w2etatag,'w2etatag/F')
      
      
      top1mtag = array( 'f', [ 0. ] )
      top1pttag = array( 'f', [ 0. ] )
      top1phitag = array( 'f', [ 0. ] )
      top1etatag = array( 'f', [ 0. ] )
      tkin.Branch('top1mtag',top1mtag,'top1mtag/F')
      tkin.Branch('top1pttag',top1pttag,'top1pttag/F')
      tkin.Branch('top1phitag',top1phitag,'top1phitag/F')
      tkin.Branch('top1etatag',top1etatag,'top1etatag/F')
      
      
      top2mtag = array( 'f', [ 0. ] )
      top2pttag = array( 'f', [ 0. ] )
      top2phitag = array( 'f', [ 0. ] )
      top2etatag = array( 'f', [ 0. ] )
      tkin.Branch('top2mtag',top2mtag,'top2mtag/F')
      tkin.Branch('top2pttag',top2pttag,'top2pttag/F')
      tkin.Branch('top2phitag',top2phitag,'top2phitag/F')
      tkin.Branch('top2etatag',top2etatag,'top2etatag/F')
      
      
      b1mtag = array( 'f', [ 0. ] )
      b1pttag = array( 'f', [ 0. ] )
      b1phitag = array( 'f', [ 0. ] )
      b1etatag = array( 'f', [ 0. ] )
      tkin.Branch('b1mtag',b1mtag,'b1mtag/F')
      tkin.Branch('b1pttag',b1pttag,'b1pttag/F')
      tkin.Branch('b1phitag',b1phitag,'b1phitag/F')
      tkin.Branch('b1etatag',b1etatag,'b1etatag/F')
      
      
      b2mtag = array( 'f', [ 0. ] )
      b2pttag = array( 'f', [ 0. ] )
      b2phitag = array( 'f', [ 0. ] )
      b2etatag = array( 'f', [ 0. ] )
      b2mtag = array( 'f', [ 0. ] )
      tkin.Branch('b2mtag',b2mtag,'b2mtag/F')
      tkin.Branch('b2pttag',b2pttag,'b2pttag/F')
      tkin.Branch('b2phitag',b2phitag,'b2phitag/F')
      tkin.Branch('b2etatag',b2etatag,'b2etatag/F')
      
      
      
      ttmtag = array( 'f', [ 0. ] )
      ttpttag = array( 'f', [ 0. ] )
      ttphitag = array( 'f', [ 0. ] )
      ttetatag = array( 'f', [ 0. ] )
      chi2tag = array( 'f', [ 0. ] )
      tkin.Branch('ttmtag',ttmtag,'ttmtag/F')
      tkin.Branch('ttpttag',ttpttag,'ttpttag/F')
      tkin.Branch('ttphitag',ttphitag,'ttphitag/F')
      tkin.Branch('ttetatag',ttetatag,'ttetatag/F')
      tkin.Branch('chi2tag',chi2tag,'chi2tag/F')
      
      
      ttmcm = array( 'f', [ 0. ] )
      ttmcpt = array( 'f', [ 0. ] )
      ttmcphi = array( 'f', [ 0. ] )
      ttmceta = array( 'f', [ 0. ] )
      b1tagmcid = array( 'i', [ 0 ] )
      b2tagmcid = array( 'i', [ 0 ] )
      b1tagmcmatchtop = array( 'i', [ 0 ])
      b2tagmcmatchtop = array( 'i', [ 0 ])
      
      ttmcm[0] = 0
      ttmcpt[0] =0
      ttmcphi[0] =0
      ttmceta[0] =0
      b1tagmcid[0] = 0
      b2tagmcid[0] = 0
      b1tagmcmatchtop[0] = 0
      b2tagmcmatchtop[0] = 0
      
      tkin.Branch('ttmcm',ttmcm,'ttmcm/F')
      tkin.Branch('ttmcpt',ttmcpt,'ttmcpt/F')
      tkin.Branch('ttmcphi',ttmcphi,'ttmcphi/F')
      tkin.Branch('ttmceta',ttmceta,'ttmceta/F')
      tkin.Branch('ttmceta',ttmceta,'ttmceta/F')
      tkin.Branch('b1tagmcid',b1tagmcid,'b1tagmcid/I')
      tkin.Branch('b2tagmcid',b2tagmcid,'b2tagmcid/I')
      tkin.Branch('b1tagmcmatchtop',b1tagmcmatchtop,'b1tagmcmatchtop/I')
      tkin.Branch('b2tagmcmatchtop',b2tagmcmatchtop,'b2tagmcmatchtop/I')





      
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
                  if not tagged:
                        notagcomb_ = JetComb(t.njets)
                                    
                        for i in range(len(notagcomb_)):
                              fitJets = testKinFitFWLite(jets_[notagcomb_[i][0]],jets_[notagcomb_[i][1]],jets_[notagcomb_[i][2]],jets_[notagcomb_[i][3]],jets_[notagcomb_[i][4]], jets_[notagcomb_[i][5]]) 
                              if bestchi2 > fitJets[6].E(): # fitJets[6].E() is the chi2 of the combination
                                    del bestcomb_[:]
                                    bestchi2 = fitJets[6].E()
                                    bestpermut_ = [int(d) for d in str(int(fitJets[6].Pt()))] # fitJets[6].Pt() is the permutation that results in the chi2 value
                                    if len(bestpermut_) == 5: bestpermut_ = [0] + bestpermut_
                                    for  j in  bestpermut_:
                                          bestcomb_.append(notagcomb_[i][j]) 

#                        print bestcomb_, ' bestcomb from text '
                        chi2[0] = bestchi2            
                        w1m[0] = (jets_[bestcomb_[2]]+jets_[bestcomb_[3]]).M()
                        w1pt[0] = (jets_[bestcomb_[2]]+jets_[bestcomb_[3]]).Pt()
                        w1phi[0] = (jets_[bestcomb_[2]]+jets_[bestcomb_[3]]).Phi()
                        w1eta[0] = (jets_[bestcomb_[2]]+jets_[bestcomb_[3]]).Eta()
                        
                        ttm[0] = (jets_[bestcomb_[0]]+jets_[bestcomb_[1]]+jets_[bestcomb_[2]]+jets_[bestcomb_[3]]+jets_[bestcomb_[4]]+jets_[bestcomb_[5]]).M()
                        ttpt[0] = (jets_[bestcomb_[0]]+jets_[bestcomb_[1]]+jets_[bestcomb_[2]]+jets_[bestcomb_[3]]+jets_[bestcomb_[4]]+jets_[bestcomb_[5]]).Pt()
                        ttphi[0] = (jets_[bestcomb_[0]]+jets_[bestcomb_[1]]+jets_[bestcomb_[2]]+jets_[bestcomb_[3]]+jets_[bestcomb_[4]]+jets_[bestcomb_[5]]).Phi()
                        tteta[0] = (jets_[bestcomb_[0]]+jets_[bestcomb_[1]]+jets_[bestcomb_[2]]+jets_[bestcomb_[3]]+jets_[bestcomb_[4]]+jets_[bestcomb_[5]]).Eta()
                        
                        w2m[0] = (jets_[bestcomb_[5]]+jets_[bestcomb_[4]]).M()
                        w2pt[0] = (jets_[bestcomb_[5]]+jets_[bestcomb_[4]]).Pt()
                        w2phi[0] = (jets_[bestcomb_[5]]+jets_[bestcomb_[4]]).Phi()
                        w2eta[0] =(jets_[bestcomb_[5]]+jets_[bestcomb_[4]]).Eta()
                        
                        top1m[0] = (jets_[bestcomb_[0]]+jets_[bestcomb_[3]]+jets_[bestcomb_[2]]).M()
                        top1pt[0] = (jets_[bestcomb_[0]]+jets_[bestcomb_[3]]+jets_[bestcomb_[2]]).Pt()
                        top1phi[0] = (jets_[bestcomb_[0]]+jets_[bestcomb_[3]]+jets_[bestcomb_[2]]).Phi()
                        top1eta[0] = (jets_[bestcomb_[0]]+jets_[bestcomb_[3]]+jets_[bestcomb_[2]]).Eta()
                        
                        top2m[0] = (jets_[bestcomb_[5]]+jets_[bestcomb_[4]]+jets_[bestcomb_[1]]).M()
                        top2pt[0] = (jets_[bestcomb_[5]]+jets_[bestcomb_[4]]+jets_[bestcomb_[1]]).Pt()
                        top2phi[0] = (jets_[bestcomb_[5]]+jets_[bestcomb_[4]]+jets_[bestcomb_[1]]).Phi()
                        top2eta[0] =(jets_[bestcomb_[5]]+jets_[bestcomb_[4]]+jets_[bestcomb_[1]]).Eta()
                        
                        
                        b1m[0] = (jets_[bestcomb_[0]]).M()
                        b1pt[0] = (jets_[bestcomb_[0]]).Pt()
                        b1phi[0] = (jets_[bestcomb_[0]]).Phi()
                        b1eta[0] = (jets_[bestcomb_[0]]).Eta()
                        
                        b2m[0] = (jets_[bestcomb_[1]]).M()
                        b2pt[0] = (jets_[bestcomb_[1]]).Pt()
                        b2phi[0] = (jets_[bestcomb_[1]]).Phi()
                        b2eta[0] = (jets_[bestcomb_[1]]).Eta()


                        b1mcid[0] =t.jets_mcFlavour[bestcomb_[0]] 
                        b2mcid[0] = t.jets_mcFlavour[bestcomb_[1]]
                                    
                  bestchi2 = 9999
#                  print t.njets, ' njets ',nbjets, ' nbjets ', ' bpos ', bpos_ 
                  tagcomb_ = JetComb(t.njets,nbjets,bpos_)
#                  print tagcomb_, ' tagcomb '
                  for i in range(len(tagcomb_)):
                        fitJets = testKinFitFWLite(jets_[tagcomb_[i][0]],jets_[tagcomb_[i][1]],jets_[tagcomb_[i][2]],jets_[tagcomb_[i][3]],jets_[tagcomb_[i][4]], jets_[tagcomb_[i][5]],1)
                        if bestchi2 > fitJets[6].E():
                              del bestcomb_[:]
                              bestchi2 = fitJets[6].E()
                              bestpermut_ = [0]+[int(d) for d in str(int(fitJets[6].Pt()))]
#                              print bestpermut_, ' best permut '
                              for  j in  bestpermut_:
                                    bestcomb_.append(tagcomb_[i][j]) 

#                  print bestcomb_, ' bestcomb from tag '
            
                  chi2tag[0] = bestchi2

                  w1mtag[0] = (jets_[bestcomb_[2]]+jets_[bestcomb_[3]]).M()
                  w1pttag[0] = (jets_[bestcomb_[2]]+jets_[bestcomb_[3]]).Pt()
                  w1phitag[0] = (jets_[bestcomb_[2]]+jets_[bestcomb_[3]]).Phi()
                  w1etatag[0] = (jets_[bestcomb_[2]]+jets_[bestcomb_[3]]).Eta()
                  
                  ttmtag[0] = (jets_[bestcomb_[0]]+jets_[bestcomb_[1]]+jets_[bestcomb_[2]]+jets_[bestcomb_[3]]+jets_[bestcomb_[4]]+jets_[bestcomb_[5]]).M()
                  ttpttag[0] = (jets_[bestcomb_[0]]+jets_[bestcomb_[1]]+jets_[bestcomb_[2]]+jets_[bestcomb_[3]]+jets_[bestcomb_[4]]+jets_[bestcomb_[5]]).Pt()
                  ttphitag[0] = (jets_[bestcomb_[0]]+jets_[bestcomb_[1]]+jets_[bestcomb_[2]]+jets_[bestcomb_[3]]+jets_[bestcomb_[4]]+jets_[bestcomb_[5]]).Phi()
                  ttetatag[0] = (jets_[bestcomb_[0]]+jets_[bestcomb_[1]]+jets_[bestcomb_[2]]+jets_[bestcomb_[3]]+jets_[bestcomb_[4]]+jets_[bestcomb_[5]]).Eta()
                  
                  w2mtag[0] = (jets_[bestcomb_[5]]+jets_[bestcomb_[4]]).M()
                  w2pttag[0] = (jets_[bestcomb_[5]]+jets_[bestcomb_[4]]).Pt()
                  w2phitag[0] = (jets_[bestcomb_[5]]+jets_[bestcomb_[4]]).Phi()
                  w2etatag[0] =(jets_[bestcomb_[5]]+jets_[bestcomb_[4]]).Eta()
                  
                  top1mtag[0] = (jets_[bestcomb_[0]]+jets_[bestcomb_[3]]+jets_[bestcomb_[2]]).M()
                  top1pttag[0] = (jets_[bestcomb_[0]]+jets_[bestcomb_[3]]+jets_[bestcomb_[2]]).Pt()
                  top1phitag[0] = (jets_[bestcomb_[0]]+jets_[bestcomb_[3]]+jets_[bestcomb_[2]]).Phi()
                  top1etatag[0] = (jets_[bestcomb_[0]]+jets_[bestcomb_[3]]+jets_[bestcomb_[2]]).Eta()
                  
                  top2mtag[0] = (jets_[bestcomb_[5]]+jets_[bestcomb_[4]]+jets_[bestcomb_[1]]).M()
                  top2pttag[0] = (jets_[bestcomb_[5]]+jets_[bestcomb_[4]]+jets_[bestcomb_[1]]).Pt()
                  top2phitag[0] = (jets_[bestcomb_[5]]+jets_[bestcomb_[4]]+jets_[bestcomb_[1]]).Phi()
                  top2etatag[0] =(jets_[bestcomb_[5]]+jets_[bestcomb_[4]]+jets_[bestcomb_[1]]).Eta()
                  
                  
                  b1mtag[0] = (jets_[bestcomb_[0]]).M()
                  b1pttag[0] = (jets_[bestcomb_[0]]).Pt()
                  b1phitag[0] = (jets_[bestcomb_[0]]).Phi()
                  b1etatag[0] = (jets_[bestcomb_[0]]).Eta()
                  
                  b2mtag[0] = (jets_[bestcomb_[1]]).M()
                  b2pttag[0] = (jets_[bestcomb_[1]]).Pt()
                  b2phitag[0] = (jets_[bestcomb_[1]]).Phi()
                  b2etatag[0] = (jets_[bestcomb_[1]]).Eta()
                  
            
                  if t.ngenTopHad == 2:
                        top1 = ROOT.TLorentzVector()
                        top2 = ROOT.TLorentzVector()
                        top1.SetPtEtaPhiM(t.genTopHad_pt[0],t.genTopHad_phi[0],t.genTopHad_eta[0],t.genTopHad_mass[0])
                        top2.SetPtEtaPhiM(t.genTopHad_pt[1],t.genTopHad_phi[1],t.genTopHad_eta[1],t.genTopHad_mass[1])
                        ttmcm[0] = (top1+top2).M()
                        ttmcpt[0] = (top1+top2).Pt()
                        ttmcphi[0] = (top1+top2).Phi()
                        ttmceta[0] = (top1+top2).Eta()
                  else:
                        ttmcm[0] = -10.0
                        ttmcpt[0] =-10.0
                        ttmcphi[0] =-10.0
                        ttmceta[0] =-10.0
                        
                        b1tagmcid[0] =t.jets_mcFlavour[bestcomb_[0]] 
                        b2tagmcid[0] = t.jets_mcFlavour[bestcomb_[1]]
                        b1tagmcmatchtop[0] =t.jets_matchBfromHadT[bestcomb_[0]] 
                        b2tagmcmatchtop[0] = t.jets_matchBfromHadT[bestcomb_[1]]

                  tkin.Fill()

            

#            print "W1", (fitJets[0]+fitJets[1]).M()
#            print "Top1", (fitJets[0]+fitJets[1]+fitJets[2]).M()

#            print "W2", (fitJets[3]+fitJets[4]).M()
#            print "Top2", (fitJets[3]+fitJets[4]+fitJets[5]).M()

            if e > nmaxentries: break

      fout.Write()
      f.Close()
