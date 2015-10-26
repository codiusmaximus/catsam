import numpy as np

def makepopii(snapshot,mvir,vmax,rvir,mvir_prog,coldgas,hotgas,blowout,mstar):

#make the timesteps
    print " snapshot: %3.2e" % (snapshot)
    print "     mvir: %3.2e" % (mvir)
    print "mvir_prog: %3.2e" % (mvir_prog)
    print "     vmax: %3.2f" % (vmax)
    print "     rvir: %3.2f" % (rvir)
    print "  coldgas: %3.2e" % (coldgas)
    print "   hotgas: %3.2e" % (hotgas)
    print "  blowout: %3.2e" % (blowout)
    print "    mstar: %3.2e" % (mstar)

    nsteps = 50.
    timestep_mini = np.array(np.range(nsteps))*(t_current - t_prev)/nsteps + t_prev

#initialise some variables to return here
    coldgas_now = 0.
    hotgas_now = mvir*bac.fb
    mstar_now = 0.
    blowout_now = 0.
    macc_now = 0.

#check the DM mass difference for Macc
    delta_mdm = mvir - mvir_prog
    if delta_mdm > 0:

        for ti in range(0,nsteps):
            mstar_now += (hotgas_now/tdyn*timestep_mini[ti])

    return "return your output fields here"
    #coldgas, hotgas, blowout, mstar
