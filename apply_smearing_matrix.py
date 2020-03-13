import os,sys
import numpy as np
import ROOT as rt


def hist2npy_2d( hist ):
    """ convert smearing matrix histogram into numpy array """
    
    hnpy = np.zeros( (hist.GetXaxis().GetNbins(),hist.GetYaxis().GetNbins()) )
    for xbin in xrange( 1,hist.GetXaxis().GetNbins()+1):
        for ybin in xrange( 1, hist.GetYaxis().GetNbins()+1 ):
            hnpy[xbin-1,ybin-1] = hist.GetBinContent( xbin, ybin )
        #print "xbin[",xbin-1,"]: ",hnpy[xbin-1,:].sum()

    return hnpy
    

def get_time_projections( rootfile="sns_out_BERT_convolved.root", tbins=[1000,1826,2201,20001] ):
    """ get 1D PDF distributions for time for each flavor """
    rfile = rt.TFile( rootfile )
    nue2d     = rfile.Get("convolved_energy_time_of_nu_e")
    numu2d    = rfile.Get("convolved_energy_time_of_nu_mu")
    numubar2d = rfile.Get("convolved_energy_time_of_anti_nu_mu")

    px_nue = nue2d.ProjectionX()
    px_numu = numu2d.ProjectionX()
    px_numubar = numubar2d.ProjectionX()
    
    tnue  = np.zeros( 3 )
    tnumu = np.zeros( 3 )
    tnumubar = np.zeros( 3 )
    for ibin in xrange(3):
        tnue[ibin]     = px_nue.Integral(     tbins[ibin]+1, tbins[ibin+1] )
        tnumu[ibin]    = px_numu.Integral(    tbins[ibin]+1, tbins[ibin+1] )
        tnumubar[ibin] = px_numubar.Integral( tbins[ibin]+1, tbins[ibin+1] )        
    tnue[:] /= tnue.sum()
    tnumu[:] /= tnumu.sum()
    tnumubar[:] /= tnumubar.sum()
    return tnumu, tnue, tnumubar

def apply_smearing_matrix( data, ccsmear, ncsmear, hcc, hnc ):
    """
    apply smearing matrix to calculation

    index of data["events"] numpy array containing events
    [ 6   nue_Ar40_marley1 ]  327.01578253498656
    [ 7   nuebar_Ar40 ]  0.0
    [ 8   nc_nue_Ar40 ]  26.45700936979586
    [ 9   nc_numu_Ar40 ]  10.90995486894403
    [ 10   nc_nutau_Ar40 ]  0.006368001369250145
    [ 11   nc_nuebar_Ar40 ]  0.009812826296037637
    [ 12   nc_numubar_Ar40 ]  44.057154025162575
    [ 13   nc_nutaubar_Ar40 ]  0.009812826296037637
    """

    # energy bins
    en = data["en"]    
    # get cc events
    cc = data["events"][6,:] + data["events"][7,:]
    # get nc events
    nc = np.zeros( en.shape[0] )
    for i in xrange(8,14):
        nc += data["events"][i,:]

    if hcc.GetYaxis().GetNbins()!=hnc.GetYaxis().GetNbins():
        raise ValueError("NC and CC smearing matrix has different number of reco bins. not supported")
        
    outen = np.zeros( hcc.GetYaxis().GetNbins()+1 )
    for ibin in xrange( hcc.GetYaxis().GetNbins()+1 ):
        outen[ibin] = hcc.GetYaxis().GetBinLowEdge(ibin+1)

    outcc = np.zeros( hcc.GetYaxis().GetNbins() )
    outnc = np.zeros( hnc.GetYaxis().GetNbins() )    
        
    # loop over energy bins
    nlost_cc = 0.0
    nlost_nc = 0.0    
    for ibin in xrange( len(en)-1 ):
        if en[ibin]<=0.0 or cc[ibin]<=0.0:
            #if cc[ibin]>0:
            #print "events lost: ",cc[ibin]            
            nlost_cc += cc[ibin]            
            continue
        if en[ibin]>=hcc.GetXaxis().GetXmax():
            #if cc[ibin]>0:            
            #    #print "events lost: ",cc[ibin]
            nlost_cc += cc[ibin]            
            continue

        ebin1 = en[ibin]
        ebin2 = en[ibin+1]
        xbin  = hcc.GetXaxis().FindBin( en[ibin] )
        lobin = hcc.GetXaxis().GetBinLowEdge( xbin )
        mibin = hcc.GetXaxis().GetBinLowEdge( xbin+1 )
        hibin = hcc.GetXaxis().GetBinLowEdge( xbin+2 )

        bincc = np.zeros( hcc.GetYaxis().GetNbins() )
        binnc = np.zeros( hnc.GetYaxis().GetNbins() )
        
        if ebin1<mibin and ebin2<=mibin:
            # event true e bin, inside smearing matrix true e bin

            # cc events
            bincc[:] += cc[ibin]*ccsmear[xbin,:]
            # nc events
            binnc[:] += nc[ibin]*ncsmear[xbin,:]
            
        elif ebin1<mibin and ebin2>mibin:
            # event true e bin straddles smearing matrix true e bin boundary

            # cc events
            bincc[:] += (mibin-ebin1)/(ebin2-ebin1)*cc[ibin]*ccsmear[xbin,:]
            bincc[:] += (ebin2-mibin)/(ebin2-ebin1)*cc[ibin]*ccsmear[xbin+1,:]
            # nc events
            binnc[:] += (mibin-ebin1)/(ebin2-ebin1)*nc[ibin]*ncsmear[xbin,:]
            binnc[:] += (ebin2-mibin)/(ebin2-ebin1)*nc[ibin]*ncsmear[xbin+1,:]            
            #print "ibin[",ibin,"] ebin=",[ebin1,ebin2],": hist true bin=",xbin," edges: ",[lobin,mibin,hibin]," inputcc=",cc[ibin]," outbincc=",bincc.sum(),"diff=",cc[ibin]-bincc.sum()
        outcc += bincc
        outnc += binnc

    print "cc: input=",cc.sum()," out=",outcc.sum()
    print "nc: input=",nc.sum()," out=",outnc.sum()
    #print "reco cc: ",outcc
    #print "reco lost: ",nlost
    #print "bins: ",outen
    data["ccreco"]    = outcc
    data["ncreco"]    = outnc
    data["en_ccreco"] = outen
    return outcc.sum(),outnc.sum()

def apply_smearing_and_time_matrix( data_prompt, data_delayed, ccsmear, ncsmear, hcc, hnc, tnumu, tnue, tnumubar ):
    """
    apply smearing matrix to calculation

    index of data["events"] numpy array containing events
    [ 6   nue_Ar40_marley1 ]  327.01578253498656
    [ 7   nuebar_Ar40 ]  0.0
    [ 8   nc_nue_Ar40 ]  26.45700936979586
    [ 9   nc_numu_Ar40 ]  10.90995486894403
    [ 10   nc_nutau_Ar40 ]  0.006368001369250145
    [ 11   nc_nuebar_Ar40 ]  0.009812826296037637
    [ 12   nc_numubar_Ar40 ]  44.057154025162575
    [ 13   nc_nutaubar_Ar40 ]  0.009812826296037637
    """

    # true energy bins
    en = data_prompt["en"]    
    # get cc events, in true energy bins
    cc_prompt  = data_prompt["events"][6,:]  + data_prompt["events"][7,:]
    cc_delayed = data_delayed["events"][6,:] + data_delayed["events"][7,:]    
    # get nc events, in true energy bins
    nc_prompt  = np.zeros( en.shape[0] )
    nc_delayed = np.zeros( en.shape[0] )    
    for i in xrange(8,14):
        nc_prompt  += data_prompt["events"][i,:]
        nc_delayed += data_delayed["events"][i,:]

    if hcc.GetYaxis().GetNbins()!=hnc.GetYaxis().GetNbins():
        raise ValueError("NC and CC smearing matrix has different number of reco bins. not supported")
        
    outen = np.zeros( hcc.GetYaxis().GetNbins()+1 )
    for ibin in xrange( hcc.GetYaxis().GetNbins()+1 ):
        outen[ibin] = hcc.GetYaxis().GetBinLowEdge(ibin+1)

    # output in (recoE,time)
    outcc = np.zeros( (hcc.GetYaxis().GetNbins(),3) )
    outnc = np.zeros( (hnc.GetYaxis().GetNbins(),3) )    
        
    # loop over true energy bins
    nlost_cc = 0.0
    nlost_nc = 0.0    
    for ibin in xrange( len(en)-1 ):
        if en[ibin]<=0.0 or ibin>175:
            if (cc_prompt[ibin]+cc_delayed[ibin])>0:
                print "events lost from true bin[",ibin,"] en=",en[ibin],": ",cc_prompt[ibin]+cc_delayed[ibin]
            nlost_cc += cc_prompt[ibin]+cc_delayed[ibin]
            continue
        if en[ibin]>=hcc.GetXaxis().GetXmax():
            if cc_prompt[ibin]+cc_delayed[bin]>0:            
                print "events lost: ",cc_prompt[ibin]+cc_delayed[ibin]
            nlost_cc += cc_prompt[ibin]+cc_delayed[ibin]
            continue

        # get data matrix true bin edges
        ebin1 = en[ibin]
        ebin2 = en[ibin+1]
        # get smear matrix true bin index
        xbin  = hcc.GetXaxis().FindBin( en[ibin] )
        # get smear matrix true bin edges (and next bin edge)
        lobin = hcc.GetXaxis().GetBinLowEdge( xbin )
        mibin = hcc.GetXaxis().GetBinLowEdge( xbin+1 )
        hibin = hcc.GetXaxis().GetBinLowEdge( xbin+2 )

        # get bin contributions for this true energy bin
        bincc_prompt  = np.zeros( (hcc.GetYaxis().GetNbins(),3) )
        bincc_delayed = np.zeros( (hcc.GetYaxis().GetNbins(),3) )        
        binnc_prompt  = np.zeros( (hnc.GetYaxis().GetNbins(),3) )
        binnc_delayed = np.zeros( (hnc.GetYaxis().GetNbins(),3) )        
        
        if ebin1<mibin and ebin2<=mibin:
            # event true e bin, inside smearing matrix true e bin and over time bins
            for tbin in xrange(3):
                # cc events
                bincc_prompt[:,tbin]  += cc_prompt[ibin]*ccsmear[xbin,:]*tnumu[tbin]
                bincc_delayed[:,tbin] += cc_delayed[ibin]*ccsmear[xbin,:]*tnue[tbin]
                # nc events
                binnc_prompt[:,tbin]  += nc_prompt[ibin]*ncsmear[xbin,:]*tnumu[tbin]
                binnc_delayed[:,tbin] += nc_delayed[ibin]*ncsmear[xbin,:]*tnue[tbin]
            
        elif ebin1<mibin and ebin2>mibin:
            # event true e bin straddles smearing matrix true e bin boundary

            # cc events
            for tbin in xrange(3):
                bincc_prompt[:,tbin] += (mibin-ebin1)/(ebin2-ebin1)*cc_prompt[ibin]*ccsmear[xbin,:]*tnumu[tbin]
                bincc_prompt[:,tbin] += (ebin2-mibin)/(ebin2-ebin1)*cc_prompt[ibin]*ccsmear[xbin+1,:]*tnumu[tbin]

                bincc_delayed[:,tbin] += (mibin-ebin1)/(ebin2-ebin1)*cc_delayed[ibin]*ccsmear[xbin,:]*tnue[tbin]
                bincc_delayed[:,tbin] += (ebin2-mibin)/(ebin2-ebin1)*cc_delayed[ibin]*ccsmear[xbin+1,:]*tnue[tbin]

                # nc events
                binnc_prompt[:,tbin] += (mibin-ebin1)/(ebin2-ebin1)*nc_prompt[ibin]*ncsmear[xbin,:]*tnumu[tbin]
                binnc_prompt[:,tbin] += (ebin2-mibin)/(ebin2-ebin1)*nc_prompt[ibin]*ncsmear[xbin+1,:]*tnumu[tbin]  

                binnc_delayed[:,tbin] += (mibin-ebin1)/(ebin2-ebin1)*nc_delayed[ibin]*ncsmear[xbin,:]*tnue[tbin]
                binnc_delayed[:,tbin] += (ebin2-mibin)/(ebin2-ebin1)*nc_delayed[ibin]*ncsmear[xbin+1,:]*tnue[tbin]

        #print "ibin[",ibin,"] ebin=",[ebin1,ebin2],": hist true bin=",xbin," edges: ",[lobin,mibin,hibin]
        #print "  inputcc: prompt=",cc_prompt[ibin]," delayed=",cc_delayed[ibin]
        #print "  outbincc prompt=",bincc_prompt.sum()," delayed=",bincc_delayed.sum()

        outcc += bincc_prompt
        outcc += bincc_delayed
        
        outnc += binnc_prompt
        outnc += binnc_delayed

    print "cc: input=",cc_prompt.sum()+cc_delayed.sum()," out=",outcc.sum()
    print "nc: input=",nc_prompt.sum()+nc_delayed.sum()," out=",outnc.sum()
    #print "reco cc: ",outcc
    #print "reco lost: ",nlost
    #print "bins: ",outen
    data = {"ccreco":outcc,
            "ncreco":outnc,
            "en_ccreco":outen}
    return outcc.sum(),outnc.sum(),data
