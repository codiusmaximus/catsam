# THIS IS THE MASTER CODE THAT BRINGS/CALLS EVERYTHING TOGETHER

# ---initiating block---
import sam_anlaysis_simple as sa
import numpy as np
import pandas as pd
import lwmodule as lwm
import haloutils as htils
import numpy as np
import sys,os

hpath = "/bigbang/data/AnnaGroup/caterpillar/halos/middle_mass_halos/H1387186/H1387186_EB_Z127_P7_LN7_LX14_O4_NV4/"

# ---define constants etc.---
min_snap = 3
max_snap = 97

# --- reading the catalogues created by you---
input_cat_path = "/bigbang/data/bagarwal/halocatalogues/input/"
input_cat_file_prefix = "snapshot_"

# Since consistent tree finder inserts haloes between snaps i-1 and i+1, if the halo was skipped at i,
# the link file is not needed. All descendants of a halo at snap i will be found at i+1

output_cat_path = "/bigbang/data/bagarwal/halocatalogues/output/"
output_cat_file_prefix = "snapshot_out_"

# ---- reading in the stardata-----
print "> Initializing star data..."
star_file_name = "/bigbang/data/bagarwal/halocatalogues/globalstarcat"
star_cols = ['snapshot', 'mstar' , 'posx','posy','posz', 'type', 'tform', 'tdeath']

stardata_df = pd.DataFrame({'snapshot':min_snap,
                            'mstar':[],
                            'posx':[],
                            'posy':[],
                            'posz':[],
                            'type':[],
                            'tform':[],
                            'tdeath':[]},index=[0])

stardata_df.to_csv(star_file_name,sep=' ')
stardata = pd.read_csv(star_file_name, delim_whitespace=True)

# ---- inital output for min_snap -----
print "> Initializing output data..."
input_cols = ['snapshot', 'tree_id', 'id', 'pid', 'origid', 'desc_id', 'scale', 'phantom', 'mvir', 'rvir', 'rs', 'vrms',
              'mmp', 'scale_of_last_MM', 'vmax', 'posx', 'posy', 'posz', 'spin']

fof_cols = ['key', 'bh_switch', 'Jiii', 'Jii', 'Jbg','coldgas', 'hotgas', 'blowout', 'mstar']

output_cols = input_cols + fof_cols

first_in = input_cat_path + input_cat_file_prefix + str(min_snap)
first_data = pd.read_csv(first_in, delim_whitespace=True, names=input_cols)

ffo = {key:0. for key in fof_cols}

for ini_index, halo_first in first_data.iterrows():

    o_key , o_bh_swtich, o_coldgas, o_hotgas, o_blowout, o_mstar = sa.main_worker_ini(haloi,ffo, stardata,min_snap)

    first_out_filename = output_cat_path + output_cat_file_prefix + str(min_snap)

    ffo_data = pd.DataFrame({'key':o_key, 
                             'bh_switch':o_bh_swtich,
                             'coldgas':o_coldgas, 'hotgas':o_hotgas, 
                             'blowout':0, 
                             'mstar':o_mstar},index=first_data.index)
    
    first_data_with_ffo = first_data.join(ffo_data)
    first_data_with_ffo.to_csv(first_out_filename,sep=' ')

print "> Now looping through snapshots (%i > %i)" % (min_snap+1,max_snap)

for snapi in np.arange(min_snap + 1, max_snap):
    z_current = htils.get_z_snap(hpath,snapi)[0]

    print " >> Snapshot: %i | %3.2f" % (snapi,z_current)
    print " >> Reading input_data..."
    input_file_name = input_cat_path + input_cat_file_prefix + str(snapi)
    input_data = pd.read_csv(input_file_name, delim_whitespace=True, names=input_cols)

    snapi_prev = snapi - 1
    print " >> Reading output_catalogue..."
    output_cat_file_name = output_cat_path + output_cat_file_prefix + str(snapi_prev)
    output_cat_data = pd.read_csv(output_cat_file_name, delim_whitespace=True)

    #sys.exit()
    #z_current = zl.zreturn(snapi)
    
    # matching has to happen here for each halo in the current input file
    print " >> Looping through halos..."
    for index, haloi in input_data.iterrows():
        match_id = np.array(np.where((output_cat_data['tree_id'] == haloi['tree_id']) & (
            output_cat_data['desc_id'] == haloi['id'])))

        # you need to add these halo fields and pass them to the main_worker, as all these matched haloes
        # from previous redshifts descend into this current guy
        if match_id.size != 0:
            ffo = {'key':np.max(output_cat_data['key'][match_id]), 
                   'bh_switch':np.max(max(output_cat_data['bh_switch'][match_id])),\
                   'mvir_prog':np.sum(output_cat_data['mvir'][match_id]),\
                   'Jiii':0,
                   'Jii':0, \
                   'Jbg':0,\
                   'coldgas':np.sum(output_cat_data['coldgas'][match_id]), \
                   'hotgas':np.sum(output_cat_data['hotgas'][match_id]), \
                   'blowout':np.sum(output_cat_data['blowout'][match_id]), \
                   'mstar':np.sum(output_cat_data['mstar'][match_id])}

        else:
            ffo = {'key':0, 'bh_switch':0, 'Jiii':0,'mvir_prog':0,'Jii':0, 'Jbg':0,'coldgas':0, 'hotgas':0, 'blowout':0, 'mstar':0}

        # make the redshift cut in the star catalogue
        # needs to be done once per snapshot
        # z_cut = lwm.lwlookback(z_current)
        # star_ids_cut = np.where(star_data['z_form'] <= z_cut)
        # NOTE: IF YOU ARE USING THE LOOPBACK MODULE, DO NOT CUT IT
        # READ IN THE STARS IN THE LWMODULE THEN

        # pass the min of these ids to the workers, so incase LW needs to be computed, you can lookup in the star
        # catalogue there and only consider the star_data from the [min:*]
        # haloi_analysed = sa.main_worker(haloi, ffo, star_ids_cut[0])

        haloi_analysed = sa.main_worker(haloi, ffo, stardata, snapi,min_snap)

        #>>> PRINT this analysed halo into the file