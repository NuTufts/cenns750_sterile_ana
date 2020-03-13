import os,sys
from array import array
import ROOT as rt

rt.gStyle.SetOptStat(0)
rt.gStyle.SetPalette(54)
rt.gStyle.SetNumberContours(2)

rt.gROOT.ProcessLine("void plot_magenta() {\
  const UInt_t Number = 4;\
  Double_t Red[Number]    = {0.30, 1.00, 1.00, 1.00};\
  Double_t Green[Number]  = {0.00, 0.00, 0.40, 0.90};\
  Double_t Blue[Number]   = {0.10, 1.00, 1.00, 1.00};\
  Double_t Length[Number] = {0.00, 0.30, 0.60, 1.00};\
  const Int_t nb=50;\
  TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);\
  gStyle->SetNumberContours(nb);\
  return;\
}")

ex1 = rt.TExec("ex1","plot_magenta();");

inputfile19 = rt.TFile("out_ana_dm2_Ue4sq_L19.root")
inputfile29 = rt.TFile("out_ana_dm2_Ue4sq_L29.root")

"""
EXPECTED CONTENTS
TFile**		out_ana_dm2_Ue4sq_L29.root	
 TFile*		out_ana_dm2_Ue4sq_L29.root	
  KEY: TH2F	h2d_nosys;1	
  KEY: TH2F	h2d_sys;1	
  KEY: TH2F	hdms2_nosys;1	
  KEY: TH2F	hdms2_sys;1	
  KEY: TH2F	hdms2_llr_nosys;1	
  KEY: TH2F	hdms2_llr_sys;1	
"""


hsys19   = inputfile19.Get("hdms2_sys")
hnosys19 = inputfile19.Get("hdms2_nosys")

hsys29   = inputfile29.Get("hdms2_sys")
hnosys29 = inputfile29.Get("hdms2_nosys")


contlevels = array('d',[5.0,200.0])
for h in [hsys19,hnosys19,hsys29,hnosys29]:
    h.SetContour(2,contlevels)

sin2_em_bf = 4.0*0.163*0.163*0.117*0.117
dm2_bf = 1.75
bf_sin2 = rt.TMarker( sin2_em_bf , dm2_bf, 20 )
bf_sin2.SetMarkerColor( rt.kRed )

c = rt.TCanvas("c","c",800,600)

# DRAW L19 first
c.Draw()
c.SetLogx(1)
c.SetLogy(1)


hnosys19.SetTitle(";sin^{2}(2#theta_{e#mu});#Delta m^{2}_{41} (eV^{2})")
hnosys19.Draw("Cont1")
hnosys19.SetLineColor(rt.kBlack)
hsys19.SetLineColor(rt.kBlack)
hnosys19.GetXaxis().SetRangeUser(3.0e-4,6.0e-3)
hsys19.SetLineStyle(2)
hsys19.Draw("Cont1SAME")

ex1.Draw()
hnosys29.SetLineColor(rt.kMagenta)
hsys29.SetLineColor(rt.kMagenta)
hnosys29.Draw("Cont1SAME")
hsys29.SetLineStyle(2)
hsys29.Draw("Cont1SAME")

l = rt.TLegend(0.15,0.15,0.6,0.3)
l.AddEntry(hnosys19,"5 sig. excl, no sys. @ 19.5 m","L")
l.AddEntry(hsys19,"5 sig. excl, w. 10% cc/nc sys. @ 19.5 m","L")
l.AddEntry(hnosys29,"5 sig. excl, no sys. @ 29 m","L")
l.AddEntry(hsys29,"5 sig. excl, w. 10% cc/nc sys. @ 29 m","L")
l.AddEntry(bf_sin2,"Global best fit @ (0.00145,1.75)","P")
l.SetBorderSize(0)
l.Draw()

bf_sin2.Draw()
c.SaveAs("inelastic_sterile_sens_chi2.pdf")

raw_input()
