# import ROOT in batch mode
import sys
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load("libPhysicsToolsKinFitter.so")
ROOT.gROOT.ProcessLine(".L testKinFitFWLite.C+")


from ROOT import testKinFitFWLite
f = ROOT.TFile.Open("/shome/clseitz/Files/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root", "read")
t = f.Get("tree")

for e,event in enumerate(t) : 
      jets = []
      for i in range(0,t.njets):
            
            jet = ROOT.TLorentzVector()
            jet.SetPtEtaPhiM(t.jets_pt[i],t.jets_eta[i],t.jets_phi[i],t.jets_mass[i])
            jets.append(jet)


      if t.njets >= 6:
            fitJets = testKinFitFWLite(jets[0],jets[1],jets[2],jets[3],jets[4], jets[5]) 

            print "W1", (fitJets[0]+fitJets[1]).M()
            print "Top1", (fitJets[0]+fitJets[1]+fitJets[2]).M()

            print "W2", (fitJets[3]+fitJets[4]).M()
            print "Top2", (fitJets[3]+fitJets[4]+fitJets[5]).M()

      if e > 10: break

