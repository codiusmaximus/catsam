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
import numpy as np
import makelwgalaxy as makelw
import lwlib
import zlookup
import popiiimodule as popiii
import popiimodule as popii
import haloutils as htils

#number of popiii stars populated in minihaloes
niii = 1

#initialise the bh_switch
bh_switch = 0

#this gets called only once at the start for initialising the star catalogue
def main_worker_ini(haloi,ffo,snapshot,min_snap):

    total_lw = 0.
    make_PopIII = popiii.checkpopiii_lw(mvir, total_lw)

    #step3: pass final output
    if make_PopIII == 'yes':
        smass_iii = popiii.makepopiii(niii)
        mstar = smass_iii

    # NEED TO PRINT THIS TO A FILE: GLOBALSTARCAT
        print "MADE IT"
        key_update = 3

def main_worker(haloi,ffo, stardata,snapshot,min_snap):
    hpath = "/bigbang/data/AnnaGroup/caterpillar/halos/middle_mass_halos/H1387186/H1387186_EB_Z127_P7_LN7_LX14_O4_NV4/"
    z_current = htils.get_z_snap(hpath,snapshot)[0]
    #print "Current redshift",z_current
#---# case empty
    if ffo['key'] == 0 :

        if snapshot != min_snap:
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
                smass_iii = popiii.makepopiii(niii)
                mstar = smass_iii
                # NEED TO PRINT THIS TO A FILE: GLOBALSTARCAT

                print "MADE IT"
                key_update = 3

#---# case DCBH
    elif ffo['key'] == 1 :

        #do nothing here
        ffo['bh_switch'] = 1

#---# case Pop II stars already exist
    # here is where the prog_info comes in handy

    elif ffo['key'] == 2 :
        popii.makepopii(haloi['snapshot'],haloi['mvir'],haloi['vmax'],haloi['rvir'],ffo['mvir_prog'],ffo['coldgas'],ffo['hotgas'],ffo['blowout'],ffo['mstar'])
        # NEED TO PRINT THIS TO A FILE: GLOBALSTARCAT
        ffo['key'] = 2
#---# case Pop III formed here at some point , make Pop II now
    elif ffo['key'] == 3 :
        popiii.makepopiii(haloi['snapshot'],haloi['mvir'],haloi['vmax'],haloi['rvir'],ffo['mvir_prog'],ffo['coldgas'],ffo['hotgas'],ffo['blowout'],ffo['mstar'])
        # NEED TO PRINT THIS TO A FILE: GLOBALSTARCAT
        ffo['key'] = 2

    if ffo['key'] != 1:
        ffo['bh_switch'] = 0

    #return the output fields
    return ffo['key'], ffo['bh_switch'], ffo['coldgas'], ffo['hotgas'], ffo['blowout'], ffo['mstar']