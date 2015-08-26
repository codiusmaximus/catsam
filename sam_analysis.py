#THIS GETS CALLED FROM THE sam_master.py
#this is the analysis part of the re-newed SAM of Agarwal et al. 12 written in python
#handles the methods that do something with the halo depending on its history

#useful definitions of variables
#------------------------------
#determines what we want to do with halo
#star_key = 3 (Pop III) , 2 (Pop II) , 1 (DCBH) , 0 (empty)
#LW flux would be computed only if halo is empty, i.e. star_key = 0

#depending on the key, the actual sub-method gets called
#in case for Pop II stars: star_key = 2, the info about the progenitors will be accessed
#the progenitor info is stored in the dictionary 'prog_info'

def main_worker(snapshot, tree_id, id, pid, origid, desc_id, scale, phantom, mvir, rvir, rs, vrms, mmp,
                scale_of_last_MM, vmax, posx, posy, posz, spin, mvir_prog, key, bh_switch, coldgas, hotgas, blowout,
                mstar)

#---# case empty
    if key == 0 :
        
        # need to make a list of active stars here based on cell_info, t_form and t_death (or t_contrib)
  
        #step1: compute the local LW radiation
        local_LW = computelocal_LW(snapshot, posx, posy, posz)        
        total_LW = local_LW + global_LW

        #step2: Pop III or no Pop III
        make_PopIII = ifpopiii(total_LW, mvir)

        #step3: pass final output
        if ifpopiii == 'yes':
            smass_iii = random_popiii_mass()
            key_update == 3
            

#---# case DCBH
    elif key == 1 :

        #do nothing here
        bh_switch = 1

#---# case Pop II stars already exist
    # here is where the prog_info comes in handy
    elif star_key == 2 :


#---# case Pop III formed here at some point , make Pop II now
    elif star_key == 3 : 

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

        macc_now = bac.fb*delta_mdm - blowout - mstar
        mhot_now = macc_now - hotgas_now/tdyn*timestep_mini[0]
        coldgas_now = mhot_now/tdyn*timestep_mini[0]

        for ti in range(1,nsteps):

            macc_step = bac.fb*delta_mdm - blowout - mstar

            mcold_now = (mhot_now/tdyn - mstar_now/timestep_mini[ti] - blowout_now/timestep_mini[ti])\
                        *timestep_mini[ti]
    else:

    return "return your output fields here"



        
