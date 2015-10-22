def popii_makestars(snapshot,mvir,vmax,rvir,mvir_prog,coldgas,hotgas,blowout,mstar)

#make the timesteps

    nsteps = 50.
    timestep_mini = np.array(range(nsteps))*(t_current - t_prev)/nsteps + t_prev

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
