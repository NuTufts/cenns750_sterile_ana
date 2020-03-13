import os,sys
from array import array
import numpy as np
import ROOT as rt
from math import sqrt,log

rt.gStyle.SetOptStat(0)

snowglobes_dir = "../snowglobes"
sterile_dir    = snowglobes_dir + "/coherent/sterile"
sys.path.append(snowglobes_dir+"/coherent/sterile")

from oscillate_flux_sterile import read_flux
from get_channel_data import get_channel_data


def read_output_data( outdir, fluxname, chanfile, exptconfig, data_dir ):

    channeldata = get_channel_data( chanfile, data_dir=data_dir )

    # Now get the info from the files.  Assumption is that all have the same binning.
    nchannels = channeldata["numchans"]
    maxpoints = 1000

    # setup data arrays
    import numpy as np
    data = { "total_events":np.zeros( maxpoints ),
             "total_events_smeared":np.zeros( maxpoints ),
             "en":np.zeros( maxpoints ),
             "eng":np.zeros( maxpoints ),
             "events":np.zeros( (nchannels,maxpoints) ),
             "events_smeared":np.zeros( (nchannels,maxpoints) ),
             "nue_es_smeared":np.zeros( maxpoints ),
             "nuebar_es_smeared":np.zeros( maxpoints ),
             "nux_es_smeared":np.zeros( maxpoints ),
             "channeldata":channeldata }

    # read in data
    for ifile in xrange(nchannels):

        # read unsmeared data        
        filename = outdir+"/"+fluxname+"_"+channeldata["channame"][ifile]+"_"+exptconfig+"_events.dat"
        f = open(filename,'r')
        ll = f.readlines()
        npoints = 0
        #print("parsing ",filename)
        for i,l in enumerate(ll):
            l = l.strip()

            if "---" in l or "Total" in l:
                continue
            
            info = l.split()            
            data["en"][i]  = float(info[0])*1000.0
            data["eng"][i] = float(info[0])*1000.0
            nevents = float(info[1])
            if abs(nevents)>1e9 or info[1].strip()=="NaN":
                nevents = 0.0
            data["events"][ ifile, i ] = nevents

            # Account for the number of targets relative to reference target-- for unweighted file only
            #data["events"][ ifile, i ] *= channeldata["num_target_factor"][ifile]

            data["total_events"][i] += data["events"][ifile,i]
            npoints += 1
        #print(" found %d points for channel %d. integral=%f"%(npoints,ifile,np.sum( data["events"][ifile,:] ) ))
        f.close()

        # read smeared data
        filename_smeared = outdir+"/"+fluxname+"_"+channeldata["channame"][ifile]+"_"+exptconfig+"_events_smeared.dat";
        f = open(filename_smeared,'r')
        ll = f.readlines()
        npoints_smeared = 0        
        for i,l in enumerate(ll):
            l = l.strip()

            if "---" in l or "Total" in l:
                continue
            
            info = l.split()
            nevents = float(info[1])
            if abs(nevents)>1e9 or info[1].strip()=="NaN":
                data["events_smeared"][ ifile, i ] = 0.0
            else:
                data["events_smeared"][ ifile, i ] = nevents

            # Account for the number of targets relative to reference target-- for unweighted file only
            #data["events_smeared"][ ifile, i ] *= channeldata["num_target_factor"][ifile]

            data["total_events_smeared"][i] += data["events_smeared"][ifile,i]
            npoints_smeared += 1

            for es_chans in ["nue_e","nuebar_e","numu_e","numubar_e","nutau_e","nutaubar_e"]:
                if channeldata["channame"][ifile] == es_chans:
                    if es_chans in ["nue_e","nuebar_e"]:
                        data["%ss_smeared"%(es_chans)][i] += data["events_smeared"][ifile,i]
                    else:
                        data["nux_es_smeared"][i] += data["events_smeared"][ifile,i]
            
        #print(" found %d points for smeared channel %d; integral=%f"%(npoints_smeared,ifile, np.sum( data["events_smeared"][ifile,:] ) ))

    return data

def hist2npy_2d( hist ):

    hnpy = np.zeros( (hist.GetXaxis().GetNbins(),hist.GetYaxis().GetNbins()) )
    for xbin in xrange( 1,hist.GetXaxis().GetNbins()+1):
        for ybin in xrange( 1, hist.GetYaxis().GetNbins()+1 ):
            hnpy[xbin-1,ybin-1] = hist.GetBinContent( xbin, ybin )
        print "xbin[",xbin-1,"]: ",hnpy[xbin-1,:].sum()

    return hnpy
    

def apply_smearing_matrix( data, ccsmear, ncsmear, hcc, hnc ):
    """
    apply smearing matrix to calculation
    """
    # get cc events
    cc = data["events"][6,:] + data["events"][7,:] # nue + nuebar
    # bin edges
    en = data["en"]
    #print "cc true E bins: ",len(cc)

    outcc = np.zeros( hcc.GetYaxis().GetNbins() )
    outen = np.zeros( hcc.GetYaxis().GetNbins()+1 )
    for ibin in xrange( hcc.GetYaxis().GetNbins()+1 ):
        outen[ibin] = hcc.GetYaxis().GetBinLowEdge(ibin+1)

    # loop over energy bins
    nlost = 0.0
    for ibin in xrange( len(en)-1 ):
        if en[ibin]<=0.0 or cc[ibin]<=0.0:
            #if cc[ibin]>0:
            #print "events lost: ",cc[ibin]            
            nlost += cc[ibin]            
            continue
        if en[ibin]>=hcc.GetXaxis().GetXmax():
            #if cc[ibin]>0:            
            #    #print "events lost: ",cc[ibin]
            nlost += cc[ibin]            
            continue

        ebin1 = en[ibin]
        ebin2 = en[ibin+1]
        xbin  = hcc.GetXaxis().FindBin( en[ibin] )
        lobin = hcc.GetXaxis().GetBinLowEdge( xbin )
        mibin = hcc.GetXaxis().GetBinLowEdge( xbin+1 )
        hibin = hcc.GetXaxis().GetBinLowEdge( xbin+2 )

        bincc = np.zeros( hcc.GetYaxis().GetNbins() )        
        if ebin1<mibin and ebin2<=mibin:
            # event true e bin, inside smearing matrix true e bin
            bincc[:] += cc[ibin]*ccsmear[xbin,:]
        elif ebin1<mibin and ebin2>mibin:
            # event true e bin straddles smearing matrix true e bin boundary
            bincc[:] += (mibin-ebin1)/(ebin2-ebin1)*cc[ibin]*ccsmear[xbin,:]
            bincc[:] += (ebin2-mibin)/(ebin2-ebin1)*cc[ibin]*ccsmear[xbin+1,:]
        #print "ibin[",ibin,"] ebin=",[ebin1,ebin2],": hist true bin=",xbin," edges: ",[lobin,mibin,hibin]," inputcc=",cc[ibin]," outbincc=",bincc.sum(),"diff=",cc[ibin]-bincc.sum()
        outcc += bincc

    #print "input cc: ",cc.sum()
    #print "output cc: ",outcc.sum()
    #print "reco cc: ",outcc
    #print "reco lost: ",nlost
    #print "bins: ",outen
    data["ccreco"]    = outcc
    data["en_ccreco"] = outen
    return outcc.sum()
    

def ana_data( oscdata, nulldata ):
    """
    [ 6   nue_Ar40_marley1 ]  327.01578253498656
    [ 7   nuebar_Ar40 ]  0.0
    [ 8   nc_nue_Ar40 ]  26.45700936979586
    [ 9   nc_numu_Ar40 ]  10.90995486894403
    [ 10   nc_nutau_Ar40 ]  0.006368001369250145
    [ 11   nc_nuebar_Ar40 ]  0.009812826296037637
    [ 12   nc_numubar_Ar40 ]  44.057154025162575
    [ 13   nc_nutaubar_Ar40 ]  0.009812826296037637
    """
    # single bin analysis
    totnc_null = 0.0
    totnc_osc  = 0.0
    for fid in xrange(8,14):
        totnc_null += np.sum(nulldata["events"][fid,:])
        totnc_osc  += np.sum(oscdata["events"][fid,:])
    diff_nc = totnc_osc-totnc_null
    onebin_nc_chi2 = diff_nc*diff_nc/totnc_null

    cc_null = nulldata["events"][6,:] + nulldata["events"][7,:]
    cc_osc  = oscdata["events"][6,:]  + oscdata["events"][7,:]
    diff_cc = cc_osc - cc_null

    nonzero = np.argwhere( cc_null>0 )
    # rebin [5,55] MeV, in 10 bins
    cc_null_rebin = np.zeros( (10) )
    cc_osc_rebin  = np.zeros( (10) )    
    for idx in nonzero:
        en = nulldata["en"][idx[0]]
        ebin = int(en)/5-1        
        #print "diff[cc][",idx[0],", bin=",ebin,"]: ",nulldata["en"][idx[0]]," ",diff_cc[idx[0]]
        cc_null_rebin[ebin] += cc_null[idx[0]]
        cc_osc_rebin[ebin]  += cc_osc[idx[0]]        

    for ebin in xrange(len(cc_null_rebin)):
        print "[ebin ",ebin," ",5+ebin*5," MeV ] null=",cc_null_rebin[ebin]," ",cc_osc_rebin[ebin]
    
    tot_cc_osc  = np.sum(cc_osc)
    tot_cc_null = np.sum(cc_null)
    tot_cc_diff = tot_cc_osc - tot_cc_null
    tot_nc_diff = totnc_osc - totnc_null
    onebin_cc_chi2 = tot_cc_diff*tot_cc_diff/tot_cc_null
    onebin_nc_chi2 = tot_nc_diff*tot_nc_diff/totnc_null

    R_null = tot_cc_null/totnc_null
    R_osc  = tot_cc_osc/totnc_osc

    R_null_sig = R_null*sqrt( 1.0/tot_cc_null + 1.0/totnc_null )
    R_osc_sig  = R_osc*sqrt( 1.0/tot_cc_osc + 1.0/totnc_osc )
    
    R_chi2 = (R_null-R_osc)*(R_null-R_osc)/(R_null_sig*R_null_sig)

    print "tot_cc(null)=",tot_cc_null
    print "tot_cc(osc)=",tot_cc_osc," tot_cc_diff(osc-null)=",tot_cc_diff
    print "tot_nc(null)=",totnc_null
    print "tot_nc(osc)=",totnc_osc," tot_nc_diff(osc-null)=",tot_nc_diff
    print "R(null)=",R_null," +/- ",R_null_sig
    print "R(osc)=",R_osc," +/- ",R_osc_sig
    print "chi2(one-bin) = ",onebin_cc_chi2," + ",onebin_nc_chi2
    print "chi2(R) = ",R_chi2

    # constrain uncertainty with NC
    nc_sig = sqrt(totnc_osc)
    ccnc_ratio_sig = 0.05 # percent uncertainty in ratio

    # calculate the log-likelihood ratio, binned-chi2, binned-chi2 w/ sys
    llr =  0.
    chi2 = 0.0
    chi2_nosys = 0.0
    for ebin in xrange(len(cc_null_rebin)):
        # neg log-likelihood
        e = cc_null_rebin[ebin]
        o = cc_osc_rebin[ebin]
        llr += 2.0*( e-o )
        if o>0:
            llr += 2.0*(log(o) - log(e))
        bin_sig = sqrt( e*(1.0+ccnc_ratio_sig) + totnc_osc )
        # bin chi2
        binchi2 = (e-o)*(e-0)/(bin_sig*bin_sig)
        chi2 += binchi2
        chi2_nosys += (e-o)*(e-0)/e
    print "-2llr: ",llr
    print "nbins: ",len(cc_null_rebin)
    print "-2llr/NDF: ",llr/float(len(cc_null_rebin))
    print "chi2(stat): ",chi2_nosys," chi2/ndf=",chi2_nosys/float(len(cc_null_rebin))
    print "chi2: ",chi2," chi2/NDF=",chi2/float(len(cc_null_rebin))
    ndf = float(len(cc_null_rebin))
    
    # Binned chi2
    #diff    = cc_osc-cc_null
    #diff2   = np.power(diff,2)
    #chi2    = diff2[cc_null>0]/cc_null[cc_null>0]
    #ndf     = len(cc_null[cc_null>0])

    #print "cc_osc: ",cc_osc[cc_null>0]
    #print "cc_null: ",cc_null[cc_null>0]
    #print "diff: ",diff[cc_null>0]
    #print "chi2: ",chi2,ndf
    #schi2 = (np.sum(chi2) + (totnc_osc-totnc_null)*(totnc_osc-totnc_null)/totnc_null)/float(ndf+1-1)
    #print "sum(diff): ",np.abs(diff[cc_null>0]).sum()
    #print "sum(diff2):",np.sum(diff2[cc_null>0])
    #print "sum(chi2):",np.sum(chi2)
    #print "chi2/ndf: ",schi2

        
    return chi2/ndf,chi2_nosys/ndf,llr/ndf

def ana_data_reco( oscdata, nulldata ):
    """
    [ 6   nue_Ar40_marley1 ]  327.01578253498656
    [ 7   nuebar_Ar40 ]  0.0
    [ 8   nc_nue_Ar40 ]  26.45700936979586
    [ 9   nc_numu_Ar40 ]  10.90995486894403
    [ 10   nc_nutau_Ar40 ]  0.006368001369250145
    [ 11   nc_nuebar_Ar40 ]  0.009812826296037637
    [ 12   nc_numubar_Ar40 ]  44.057154025162575
    [ 13   nc_nutaubar_Ar40 ]  0.009812826296037637
    """
    # single bin analysis
    totnc_null = 0.0
    totnc_osc  = 0.0
    for fid in xrange(8,14):
        totnc_null += np.sum(nulldata["events"][fid,:])
        totnc_osc  += np.sum(oscdata["events"][fid,:])
    diff_nc = totnc_osc-totnc_null
    onebin_nc_chi2 = diff_nc*diff_nc/totnc_null

    cc_null = nulldata["ccreco"]
    cc_osc  = oscdata["ccreco"]
    enreco  = nulldata["en_ccreco"]
    diff_cc = cc_osc - cc_null


    for ebin in xrange(len(cc_null)):
        print "[ebin ",ebin," ",enreco[ebin]," MeV ] null=",cc_null[ebin]," ",cc_osc[ebin]
    
    tot_cc_osc  = np.sum(cc_osc)
    tot_cc_null = np.sum(cc_null)
    tot_cc_diff = tot_cc_osc - tot_cc_null
    tot_nc_diff = totnc_osc - totnc_null
    onebin_cc_chi2 = tot_cc_diff*tot_cc_diff/tot_cc_null
    onebin_nc_chi2 = tot_nc_diff*tot_nc_diff/totnc_null

    R_null = tot_cc_null/totnc_null
    R_osc  = tot_cc_osc/totnc_osc

    R_null_sig = R_null*sqrt( 1.0/tot_cc_null + 1.0/totnc_null )
    R_osc_sig  = R_osc*sqrt( 1.0/tot_cc_osc + 1.0/totnc_osc )
    
    R_chi2 = (R_null-R_osc)*(R_null-R_osc)/(R_null_sig*R_null_sig)

    print "tot_cc(null)=",tot_cc_null
    print "tot_cc(osc)=",tot_cc_osc," tot_cc_diff(osc-null)=",tot_cc_diff
    print "tot_nc(null)=",totnc_null
    print "tot_nc(osc)=",totnc_osc," tot_nc_diff(osc-null)=",tot_nc_diff
    print "R(null)=",R_null," +/- ",R_null_sig
    print "R(osc)=",R_osc," +/- ",R_osc_sig
    print "chi2(one-bin) = ",onebin_cc_chi2," + ",onebin_nc_chi2
    print "chi2(R) = ",R_chi2

    # constrain uncertainty with NC
    nc_sig = sqrt(totnc_osc)
    ccnc_ratio_sig = 0.05 # percent uncertainty in ratio

    # calculate the log-likelihood ratio, binned-chi2, binned-chi2 w/ sys
    llr =  0.
    chi2 = 0.0
    chi2_nosys = 0.0
    ndf = 0.0    
    for ebin in xrange(len(cc_null)):
        # neg log-likelihood
        e = cc_null[ebin]
        o = cc_osc[ebin]

        if e<=0:
            continue
        
        llr += 2.0*( e-o )
        if o>0:
            llr += 2.0*(log(o) - log(e))
        bin_sig = sqrt( e*(1.0+ccnc_ratio_sig) + totnc_osc )
        # bin chi2
        binchi2 = (e-o)*(e-0)/(bin_sig*bin_sig)
        chi2 += binchi2
        chi2_nosys += (e-o)*(e-0)/e
        ndf += 1.0
    print "-2llr: ",llr
    print "nbins: ",ndf
    print "-2llr/NDF: ",llr/ndf
    print "chi2(stat): ",chi2_nosys," chi2/ndf=",chi2_nosys/ndf
    print "chi2: ",chi2," chi2/NDF=",chi2/ndf
            
    return chi2/ndf,chi2_nosys/ndf,llr/ndf

if __name__ == "__main__":

    # GET SMEARING MATRICES
    sminfile = rt.TFile("out_smearing_matrix.root","read")
    hcc = sminfile.Get("hcc_pmt")
    hnc = sminfile.Get("hnc_pmt")

    ccnpy = hist2npy_2d( hcc )
    ncnpy = hist2npy_2d( hnc )
    
    tfout = rt.TFile("out_null_reco_spectra.root","recreate")


    # get null hypothesis
    flux_unosc = read_flux( snowglobes_dir+"/coherent/sterile/fluxes/stpi.dat" )

    channame    = "argon_marley1"
    expt_config = "ar40kt"
    L = 29.0
    if False:
        pnull = os.popen( sterile_dir+"/./supernova.pl stpi %s %s 0 %s"%(channame,expt_config,snowglobes_dir) )
        for p in pnull:
            print p.strip()

    nulldata = read_output_data("out","stpi","argon_marley1","ar40kt",snowglobes_dir)
    # scale null data to 3 years
    for ic,chan in enumerate(nulldata["channeldata"]["channame"]):
        nulldata["events"][ic,:] *= 3.0*3.14e7*0.000612/40.0*(20.0*20.0)/(L*L)*(5000.0/8766.0)
        print "[ ",ic," ",chan," ] ",nulldata["events"][ic,:].sum()

    # apply smearing matrix
    apply_smearing_matrix( nulldata, ccnpy, ncnpy, hcc, hnc )


    if True:
        sys.exit(-1)

    # chi2 for dm2 vs. Ue4
    c = rt.TCanvas("osc","c",1400,600)
    c.Divide(2,1)

    c.cd(1).SetLogy(1)
    c.cd(1).SetLogx(1)
    h2d_nosys.Draw("colz")
    h2d_nosys.GetZaxis().SetRangeUser(0,5)
    
    c.cd(2).SetLogy(1)
    c.cd(2).SetLogx(1)
    h2d_sys.Draw("colz")
    h2d_sys.GetZaxis().SetRangeUser(0,5)
    c.Update()


    # vs. sin^2
    csin2 = rt.TCanvas("osc_vs_sin2","dm2 vs sin2",1400,600)
    csin2.Divide(2,1)
    
    csin2.cd(1).SetLogy(1)
    csin2.cd(1).SetLogx(1)
    h2dms2_nosys.Draw("colz")
    h2dms2_nosys.GetZaxis().SetRangeUser(0,5)
    bf_sin2.Draw()
    
    csin2.cd(2).SetLogy(1)
    csin2.cd(2).SetLogx(1)
    h2dms2_sys.Draw("colz")
    h2dms2_sys.GetZaxis().SetRangeUser(0,5)
    bf_sin2.Draw()
    csin2.Update()
    
    # LLR for dm2 vs. sin^2
    cllr_sin2 = rt.TCanvas("llr_vs_sin2","nLLR for dm2 vs sin2",800,600)
    
    cllr_sin2.cd(1).SetLogy(1)
    cllr_sin2.cd(1).SetLogx(1)
    h2dms2_llr_nosys.Draw("colz")
    h2dms2_llr_nosys.GetZaxis().SetRangeUser(0,5)
    bf_sin2.Draw()

    raw_input()
    
    tfout.Write()
    tfout.Close()
