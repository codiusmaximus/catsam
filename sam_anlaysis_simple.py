"""
THIS GETS CALLED FROM THE sam_master.py
this is the analysis part of the re-newed SAM of Agarwal et al. 12 written in python
handles the methods that do something with the halo depending on its history

useful definitions of variables
------------------------------
determines what we want to do with halo
--star_key = 3 (Pop III) , 2 (Pop II) , 1 (DCBH) , 0 (empty)

--LW flux would be computed only if halo is empty, i.e. star_key = 0

depending on the key, the actual sub-method gets called
-- in case for Pop II stars: star_key = 2, the info about the progenitors will be accessed
   the progenitor info is stored in the dictionary 'prog_info'
"""

# import
import makelwgalaxy as makelw
import lwlib
import zlookup
import popiiimodule as popiii

#def main_worker(snapshot, tree_id, id, pid, origid, desc_id, scale, phantom, mvir, rvir, rs, vrms, mmp,
#                scale_of_last_MM, vmax, posx, posy, posz, spin, mvir_prog, key, bh_switch, coldgas, hotgas, blowout,
#                mstar, min_star_cut_id, stardata):

def main_worker(haloi,ffo, stardata,snapshot):
    z_current = htils.get_z_snap(hpath,snapshot)[0]
#---# case empty
    if ffo['key'] == 0 :

        # LW to be computed here
        # Way it is done now is to lookback one snapshot, then find galaxies with stars in them, track them to rhe start
        # and make arrays of star mass and ages for that 'one' position at the last snapshot. Thus you create an SED
        # accordingly and then compute ht ekde/kdi/Jlw

        gal_mass, gal_age, pos_gal  = makelw.lwgal(snapshot, stardata)
        kde , kdi , jlw = lwlib.sedcompute(gal_mass,gal_age,pos_gal)

        jlw_global = lwlib.lwglobal(z_current)
        # Global J is computed as a fit

        total_lw = jlw + jlw_global

        #step2: Pop III or no Pop III
        make_PopIII = popiii.checkpopiii_lw(mvir, total_lw)

        #step3: pass final output
        if make_PopIII == 'yes':
            smass_iii = popiii.makepopiii()
            # NEED TO PRINT THIS TO A FILE: GLOBALSTARCAT
            key_update = 3


#---# case DCBH
    elif ffo['key'] == 1 :

        #do nothing here
        bh_switch = 1

#---# case Pop II stars already exist
    # here is where the prog_info comes in handy

    elif ffo['key'] == 2 :
        popii_makestars(snapshot,mvir,vmax,rvir,mvir_prog,coldgas,hotgas,blowout,mstar)
        # NEED TO PRINT THIS TO A FILE: GLOBALSTARCAT

#---# case Pop III formed here at some point , make Pop II now
    elif ffo['key'] == 3 :
        popiii_makestars(snapshot,mvir,vmax,rvir,mvir_prog,coldgas,hotgas,blowout,mstar)
        # NEED TO PRINT THIS TO A FILE: GLOBALSTARCAT



