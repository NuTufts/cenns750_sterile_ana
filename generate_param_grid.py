import os,sys
from math import sqrt,log,exp

"""
Generate list of points in (dm^2,Ux4)-space.
"""


def generate_param_dm2_Ue4sq_file( param_file, dm2_nbins, dm2_min, dm2_max, dU_nbins, Ue4sq_min, Ue4sq_max ):
    """ 
    generate parameter values to calculate event rate 
    generates Y, the relic abundance, and M, the mass of the DM particle.
    In the file is also epsilon, the coupling parameter put into the run config file.
    """
    
    Um4sq = 0.014
    Ut4sq = 0.2

    lndm2_min = log(dm2_min)
    lndm2_max = log(dm2_max)

    ln_U_min  = log(Ue4sq_min)
    ln_U_max  = log(Ue4sq_max)

    ddm2 = (lndm2_max-lndm2_min)/float(dm2_nbins)
    dU   = (ln_U_max-ln_U_min)/float(dU_nbins)

    parlist = []
    for idm in xrange(dm2_nbins+1):
        for iu in xrange(dU_nbins+1):

            lndm2 = lndm2_min + ddm2*idm
            lnU   = ln_U_min + dU*iu
            
            parlist.append( (exp(lnU), exp(lndm2)) )

    
    
    parfile = open(param_file,'w')
    print>>parfile,"dm2",'\t',"Ue4^2",'\t',"Um4^2",'\t',"Ut4^2"
    for (U,dm) in parlist:
        print>>parfile,dm,'\t',U,'\t',Um4sq,'\t',Ut4sq
    print "generated ",len(parlist)," parameters into ",param_file
    parfile.close()


if __name__ == "__main__":

    param_file = "param_file_dm2_Ue4sq.dat"
    generate_param_dm2_Ue4sq_file( param_file, 10, 1.0e-2, 1.0e2, 10, 1.0e-3, 1.0e-1 )


