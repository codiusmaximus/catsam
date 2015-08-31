# THIS IS THE MASTER CODE THAT BRINGS/CALLS EVERYTHING TOGETHER

# ---initiating block---
import sam_analysis_simple as sa
import numpy as np
import pandas as pd
import zlookup as zl
import lwmodule as lwm

# ---define constants etc.---

# ---reading the catalogues created by you---
input_cat_path = "/bigbang/data/bagarwal/halocatalogues/input/"
input_cat_file_prefix = "snaspshot_"

# Since consistent tree finder inserts haloes between snaps i-1 and i+1, if the halo was skipped at i,
# the link file is not needed. All descendants of a halo at snap i will be found at i+1

output_cat_path = "/bigbang/data/bagarwal/halocatalogues/output/"
output_cat_file_prefix = "snaspshot_out_"

min_snap = 3
max_snap = 97

input_cols = ['snapshot', 'tree_id', 'id', 'pid', 'origid', 'desc_id', 'scale', 'phantom', 'mvir', 'rvir', 'rs', 'vrms',
              'mmp', 'scale_of_last_MM', 'vmax', 'posx', 'posy', 'posz', 'spin']

output_cols = ['snapshot','tree_id', 'id', 'pid', 'origid', 'desc_id', 'scale', 'phantom', 'mvir', 'rvir', 'rs', 'vrms',
               'mmp','scale_of_last_MM','vmax','posx', 'posy', 'posz', 'spin', 'key', 'bh_switch', 'Jiii', 'Jii', 'Jbg',
               'coldgas', 'hotgas', 'blowout', 'mstar']

for snapi in np.arange(min_snap, max_snap):

    input_file_name = input_cat_path + input_cat_file_prefix + snapi
    input_data = pd.read_csv(input_file_name, delim_whitespace=True, names=input_cols)

    snapi_prev = snapi - 1

    output_cat_file_name = output_cat_path + output_cat_file_prefix + snapi_prev
    output_cat_data = pd.read_csv(output_cat_file_name, delim_whitespace=True, names=output_cols)

    star_data = pd.reac_csv()

    z_current = zl.zreturn(snapi)

    # matching has to happen here for each halo in the current input file
    for index, haloi in input_data.iterrows():
        match_id = np.array(np.where((output_cat_data['tree_id'] == haloi['tree_id']) & (
            output_cat_data['desc_id'] == haloi['id'])))

        # you need to add these halo fields and pass them to the main_worker, as all these matched haloes
        # from previous redshifts descend into this current guy
        if match_id.size != 0:
            fields_from_output = [sum(output_cat_data['mvir'][match_id]),max(output_cat_data['key'][match_id]),
                                  max(output_cat_data['bh_switch'][match_id]),sum(output_cat_data['coldgas'][match_id]),
                                  sum(output_cat_data['hotgas'][match_id]),sum(output_cat_data['blowout'][match_id]),
                                  sum(output_cat_data['mstar'][match_id])]
        else:
            fields_from_output = [0., 0., 0., 0., 0., 0.]

        # make the redshift cut in the star catalogue
        # needs to be done once per snapshot
        z_cut = lwm.lwlookback(z_current)
        star_ids_cut = np.where(star_data['z_form'] <= z_cut)
        # NOTE: IF YOU ARE USING THE LOOPBACK MODULE, DO NOT CUT IT
        # READ IN THE STARS IN THE LWMODULE THEN

        # pass the min of these ids to the workers, so incase LW needs to be computed, you can lookup in the star
        # catalogue there and only consider the star_data from the [min:*]
        haloi_analysed = sa.main_worker(haloi, fields_from_output, star_ids_cut[0])

        #>>> PRINT this analysed halo into the file