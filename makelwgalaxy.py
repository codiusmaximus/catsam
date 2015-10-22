"""
Code to loop back over the progs in the star catalogue, so we can use just 1 galaxy at a time from the previous snapshot

INPUT
param snapnow: snapshot that you are at now
param stardata: stardata is the input form the file that contains all the star fields

OUTPUT

"""

import numpy as np
import ba_constants as bac

def lwgal(snapnow, stardata):

    minsnap =bac.minsnap

    # galaxies at the prev snapshot that will be passed to the lw module, with the pos at the prev snapshot
    galaxies = np.array(np.where((stardata['snapshot'] == snapnow - 1)))
    stardata_galaxies = stardata[galaxies]

    # all stars that formed before snapnow-1 to be connected to 'galaxies'
    prev_gal = np.array(np.where((stardata['snapshot'] < snapnow - 1)))
    stardata_prev_gal = stardata[prev_gal]

    M_gal = [np.array([i]) for i in np.array([stardata_galaxies['mstar']])]
    age_gal = [np.array([j]) for j in np.array([stardata_galaxies['agestar']])] # agestar needs to be defined wrt t_form
    pos_gal = [np.array([k]) for k in np.array([stardata_galaxies['postion']])]
    # correct this later to xpos,ypos,zpos

    loopback_snaps = list(reversed(range(minsnap,snapnow-2)))

    for gal_index, gali in stardata_galaxies():

        # just look for the progenitors one step back from your current galaxy giving out LW, which is at snapnow-1
        step1_ids = np.array(np.where((stardata_prev_gal['snap'] == snapnow - 2) &
                                      (stardata_prev_gal['treeid'] == stardata_galaxies['treeid'][gali]) &
                                      ((stardata_prev_gal['descid'] == stardata_galaxies['haloid'][gali]))))

        if step1_ids.size != 0:

            M_gal[gal_index] = np.append(M_gal[gal_index],stardata_prev_gal['Mstar'][step1_ids])
            age_gal[gal_index] = np.append(age_gal[gal_index],stardata_prev_gal['age'][step1_ids])

            # store the info of these progenitors at snapnow-2 into an array, and we will make the connection for each
            # of these till we find no descendants

            buffer_array = stardata_prev_gal[step1_ids]

            # looping over each of the entries in the prog list
            for prevgal_index, prev_gali in buffer_array():

                #looping back in time till no progs are found
                # append M_gal, age_gal as you go back and find entries

                for subi, subsnap in enumerate(loopback_snaps):
                    step2_ids = np.array(np.where((buffer_array['snap'] == subsnap) &
                                                  (buffer_array['treeid'] == buffer_array['treeid'][prev_gali]) &
                                                  (buffer_array['descid'] == buffer_array['haloid'][prev_gali])))

                    M_gal[gal_index] = np.append(M_gal[gal_index],buffer_array['Mstar'][step2_ids])
                    age_gal[gal_index] = np.append(age_gal[gal_index],buffer_array['age'][step2_ids])

    return M_gal,age_gal,pos_gal