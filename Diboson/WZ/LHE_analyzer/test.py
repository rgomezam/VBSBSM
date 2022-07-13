from __future__ import division

import ROOT
import lheanalyzer
import math

import sys
#######
#
# usage: python test.py input.lhe output.root
#
#######

#Create lhe analysis from lhe file
analysis = lheanalyzer.LHEAnalysis(sys.argv[1])
ROOT.TH1.SetDefaultSumw2()

outfile = ROOT.TFile(sys.argv[2] ,"RECREATE")

#h_lead_pt = ROOT.TH1F("l0_pt","p_{T}^{l0}",50,0,500)
#h_ee_pt = ROOT.TH1F("ee_pt","p_{T}^{ee}",50,0,500)
#h_Z_pt = ROOT.TH1F("Z_pt","p_{T}^{Z}",100,0,1000)
#h_Wlep_pt = ROOT.TH1F("Wlep_pt","p_{T}^{Wlep}",100,0,1000)
#h_mu_pt = ROOT.TH1F("mu_pt","p_{T}^{Wlep}",100,0,1000)


h_Mjj = ROOT.TH1F("Mjj","m_{jj}",100,0,1000)
h_MTWZ = ROOT.TH1F("MTWZ","m_{T}^{WZ}",70,0.0,3500.0)
h_mWZ = ROOT.TH1F("mWZ","m^{WZ}",100,0,10000)

h_aZylW = ROOT.TH1F("aZylW","aZylW",15,0,5)
h_centrality = ROOT.TH1F("centrality","centrality",20,-4,4)


h_dpiWZ = ROOT.TH1F("dphiWZ","dPhi^{WZ}",5,0,3.142)
h_dpiWZ_2G = ROOT.TH1F("dphiWZ_2G","dPhi^{WZ}",5,0,3.142)
h_dpiWZ_1_5G = ROOT.TH1F("dphiWZ_1_5G","dPhi^{WZ}",5,0,3.142)
h_dpiWZ_1G = ROOT.TH1F("dphiWZ_1G","dPhi^{WZ}",5,0,3.142)
h_pt3l = ROOT.TH1F("pt3l","pt^{3l}",5000,0,10000)
h_pt3l_2G = ROOT.TH1F("pt3l_2G","pt^{3l}",5000,0,10000)
h_pt3l_1_5G = ROOT.TH1F("pt3l_1_5G","pt^{3l}",5000,0,10000)
h_pt3l_1G = ROOT.TH1F("pt3l_1G","pt^{3l}",5000,0,10000)
h_MTWZ_2G = ROOT.TH1F("MTWZ_2G","m_{T}^{WZ}}",125,0,10000)
h_MTWZ_1_5G = ROOT.TH1F("MTWZ_1_5G","m_{T}^{WZ}",125,0,10000)
h_MTWZ_1G = ROOT.TH1F("MTWZ_1G","m_{T}^{WZ}",125,0,10000)
h_mWZ_beforecuts = ROOT.TH1F("mWZ_beforecuts","m^{WZ}",100,0,10000)
nevents =0
for event in analysis:
    nevents = nevents+1
    #if nevents%10000==0:
#    print "EVENT "
#    print nevents
    jets=[21,2,4,1,3,-2,-4,-1,-3,5, -5] 
 
    ls_WZ = filter(lambda particle: abs(particle.pdgId) in [23,24], event.particles)
    ls_leps = filter(lambda particle: abs(particle.pdgId) in [11,13], event.particles)
##################VBS#################
    ls_jets = filter(lambda particle: abs(particle.pdgId) in jets, event.particles)
    ls_finaljets = filter(lambda particle: particle.status in [1], ls_jets)
    ls_finalBjets = filter(lambda particle: abs(particle.pdgId) in[5,-5],ls_finaljets)
    ls_Z=filter(lambda particle: abs(particle.pdgId) in [23],event.particles) 
    ls_Ze=filter(lambda particle: abs(particle.pdgId) in [11],event.particles) 
    ls_Wl=filter(lambda particle: abs(particle.pdgId) in [13],event.particles) 
    if len(ls_WZ)==2 and len(ls_leps)==3: 
        h_mWZ_beforecuts.Fill((ls_WZ[0].tLorentzVector+ls_WZ[1].tLorentzVector).M(),event.weight)
        if len(ls_finaljets)==2 and len(ls_finalBjets)==0:
            if (ls_finaljets[0].tLorentzVector+ls_finaljets[1].tLorentzVector).M()>=500. and  ls_finaljets[0].eta*ls_finaljets[1].eta<0:
                deltaR=math.sqrt((ls_finaljets[0].eta-ls_finaljets[1].eta)*(ls_finaljets[0].eta-ls_finaljets[1].eta)+ (ls_finaljets[0].phi-ls_finaljets[1].phi)*(ls_finaljets[0].phi-ls_finaljets[1].phi))   
                if deltaR>0.3:
                    m_WZ=(ls_WZ[0].tLorentzVector+ls_WZ[1].tLorentzVector).M() 
                    h_mWZ.Fill((ls_WZ[0].tLorentzVector+ls_WZ[1].tLorentzVector).M(),event.weight)
#                    h_MTWZ.Fill((ls_WZ[0].tLorentzVector+ls_WZ[1].tLorentzVector).Mt(),event.weight)
                    mtwz_timi=(ls_WZ[0].tLorentzVector+ls_WZ[1].tLorentzVector).Mt()
                    pt3l_timi=ls_leps[0].pt+ls_leps[1].pt+ls_leps[2].pt
                    dpiwz_timi=abs(ls_WZ[0].phi-ls_WZ[1].phi)
                    h_MTWZ.Fill(mtwz_timi,event.weight)
                    h_dpiWZ.Fill(dpiwz_timi,event.weight)
                    h_pt3l.Fill(pt3l_timi,event.weight)
#                    azylw = abs(ls_Z[0].y-ls_Wl[0].y)
                    #print azylw;
#                    h_aZylW.Fill(azylw,event.weight)
#                    eta_W=ls_Wl[0].eta
#                    eta_Z1=ls_Ze[0].eta
#                    eta_Z2=ls_Ze[1].eta
#                    eta_j1=ls_finaljets[0].eta
#                    eta_j2=ls_finaljets[1].eta

#                    EtaMinLept = min(eta_W,eta_Z1,eta_Z2)
#                    EtaMaxLept = max(eta_W,eta_Z1,eta_Z2)
#                    EtaMinJets = min(eta_j1,eta_j2)
#                    EtaMaxJets = max(eta_j1,eta_j2)
#                    DeltaMinus = EtaMinLept - EtaMinJets
#                    DeltaPlus = EtaMaxJets - EtaMaxLept
#                    ZetaLept = min(DeltaMinus,DeltaPlus )
                    #print ZetaLept
#                    h_centrality.Fill(ZetaLept,event.weight)
                    if m_WZ<2000.0:
                        h_dpiWZ_2G.Fill(dpiwz_timi,event.weight)
                        h_pt3l_2G.Fill(pt3l_timi,event.weight)
                        h_MTWZ_2G.Fill(mtwz_timi,event.weight)
                    if m_WZ<1500.0:
                        h_pt3l_1_5G.Fill(pt3l_timi,event.weight)
                        h_dpiWZ_1_5G.Fill(dpiwz_timi,event.weight)
                        h_MTWZ_1_5G.Fill(mtwz_timi,event.weight)
                    if m_WZ<1000.0:
                        h_MTWZ_1G.Fill(mtwz_timi,event.weight)
                        h_pt3l_1G.Fill(pt3l_timi,event.weight)
                        h_dpiWZ_1G.Fill(dpiwz_timi,event.weight)
                        


    #h_lead_pt.Fill(max(ls[0].pt,ls[1].pt),event.weight)
    #h_ee_pt.Fill((ls[0].tLorentzVector+ls[1].tLorentzVector).Pt(),event.weight)
    #h_mu_pt.Fill((ls3[0].pt),event.weight)
#    if ls2:
#        h_Z_pt.Fill((ls2[0].pt),event.weight)
#    if ls4:
#        h_Wlep_pt.Fill((ls3[0].pt),event.weight)
#    if ls6:
#    h_Mjj.Fill((ls6[0].tLorentzVector+ls6[1].tLorentzVector).M(),event.weight)
#    if len(ls7)==2:
#        h_MTWZ.Fill((ls7[0].tLorentzVector+ls7[1].tLorentzVector).Mt(),event.weight)
#        h_mWZ.Fill((ls7[0].tLorentzVector+ls7[1].tLorentzVector).M(),event.weight)
#    if len(ls7)>2:

outfile.Write()
outfile.Close()
