import os,sys
from array import array
import ROOT as rt

rt.gStyle.SetOptStat(0)

inputfile = rt.TFile("out_recoE_v_time_ana_dm2_Ue4sq.root")

"""
EXPECTED CONTENTS
TFile**		out_recoE_v_time_ana_dm2_Ue4sq.root	
 TFile*		out_recoE_v_time_ana_dm2_Ue4sq.root	
  KEY: TH2F	hnull_cc_smeared;1	
  KEY: TH2F	hnull_nc_smeared;1	
"""


# NULL
hcc2d = inputfile.Get("hnull_cc_smeared")
hnc2d = inputfile.Get("hnull_nc_smeared")
hcc = hcc2d.ProjectionX()
hnc = hnc2d.ProjectionX()

# OSC
hcc_osc2d = inputfile.Get("hosc_cc_smeared")
hnc_osc2d = inputfile.Get("hosc_nc_smeared")
hcc_osc = hcc_osc2d.ProjectionX()
hnc_osc = hnc_osc2d.ProjectionX()

c = rt.TCanvas("c","c",800,500)
c.Draw()

# NULL STACK
hstack = rt.THStack()
hstack.Add( hnc )
hstack.Add( hcc )
hstack.Draw()
hcc.SetLineColor( rt.kBlue )
hnc.SetLineColor( rt.kRed )
hcc.SetLineWidth( 1 )
hnc.SetLineWidth( 1 )

# OSC STACK
hstack_osc = rt.THStack()
hcc_osc.SetLineColor( rt.kBlue )
hnc_osc.SetLineColor( rt.kRed )
hcc_osc.SetLineWidth( 1 )
hnc_osc.SetLineWidth( 1 )
hcc_osc.SetLineStyle( 2 )
hnc_osc.SetLineStyle( 2 )
hcc_osc.SetFillStyle( 0 )
hnc_osc.SetFillStyle( 0 )
hstack_osc.Add( hnc_osc )
hstack_osc.Add( hcc_osc )
hstack_osc.Draw()


hstack.Draw("hist")
hstack_osc.Draw("histsame")
hstack.SetTitle(";reconstructed visible energy (MeVee);")

"""
fX1NDC                        0.674185            X1 point in NDC coordinates
fY1NDC                        0.593684            Y1 point in NDC coordinates
fX2NDC                        0.818296            X2 point in NDC coordinates
fY2NDC                        0.848421            Y2 point in NDC coordinates
"""

tlen = rt.TLegend( 0.67, 0.59, 0.818, 0.848 )
tlen.AddEntry( hcc, "CC #nu_{e}", "L" )
tlen.AddEntry( hnc, "NC #nu",   "L" )
tlen.AddEntry( hcc_osc, "osc. CC #nu_{e}", "L" )
tlen.AddEntry( hnc_osc, "osc. NC #nu",   "L" )
tlen.SetBorderSize(0)
tlen.Draw()

c.Update()

c.SaveAs("reco_spectrum_wosc.pdf")

raw_input()
