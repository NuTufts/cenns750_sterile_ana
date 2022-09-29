import os,sys
import ROOT as rt

from apply_smearing_matrix import get_time_projections

tnumu, tnue, tnumubar = get_time_projections()

print tnumu, tnue, tnumubar
