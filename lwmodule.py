"""
# module to compute
# A
# 1. periodic coords
# 2. light cone analysis for radiation
# 3. make selection
# B
# 1. distance
# 2. alpha/beta/rates/Jlw
"""
import numpy as np
import ba_constants as bac
import pandas as pd

def lwlookback(z_current):
    z_lookback = 12.1/11.2*(1 + z_current) - 1
    return z_lookback

def lwglobal(z):
    lw_backgorund = 1.6*np.exp(-((z - 10.9)/5.)**2 / 5)
    return lw_backgorund

def lwmodule_master(z_current, t_current, r_arr, box_l, t_contrib, t_form, z_form):
    """
    input:

    :param z_current: current redshift
    :param t_current: age of the universe now

    :param r_arr: in Mpc - halo pos - star array pos
    the x,y,z array of co-ords that represent the
    del_r = x_o - x_i , y_o - y_i , z_o - z_i
    for all i stars/objects (i) measured from point 'r_o'

    :param box_l: in co-moving Mpc

    :t_contrib: for stars

    :t_form: for stars

    :param alpha, beta: rate params for each stellar population

    :param L_lw: the 13.6 output in erg/s/Hz to be converted to specific intensity J21

    output:
    k_de, k_di, J21
    """

    j21 = lwglobal(z_current)
    k_de = 0.
    k_di = 0.

    # computing the periodic bit
    n_rows = len(r_arr)
    coord_arr_periodic = np.zeros((n_rows, 3), float)
    dist_arr_comov = np.zeros((n_rows), float)
    cols = 3

    for ri in range(0, n_rows):
        for icord in range(0, cols):
            if r_arr[ri][icord] > box_l * 0.5:
                coord_arr_periodic[ri][icord] = r_arr[ri][icord] - box_l * 0.5
            elif r_arr[ri][icord] <= (-1.) * box_l * 0.5:
                coord_arr_periodic[ri][icord] = r_arr[ri][icord] + box_l * 0.5
            else:
                coord_arr_periodic[ri][icord] = r_arr[ri][icord]

        # computing distances here
        dist_arr_comov[ri] = np.sum(coord_arr_periodic[ri]**2)  # comoving distance

        # light travel time beteeen now amd the formation time of the star
        d_lt_if = bac.year * bac.c * (t_current - t_form[ri]) / bac.parsec

        # light travel time between now and the death time of the star
        d_lt_ic = bac.parsec * bac.c * (t_current - t_contrib[ri]) / bac.parsec

        # selecting the ones that satisfy the light cone analysis
        if dist_arr_comov[ri] / (1 + z_form) <= d_lt_if[ri] & dist_arr_comov[ri] / (1 + z_form) > d_lt_ic:

            j_buffer = l_lw / (dist_arr_comov[ri] * bac.parsec * bac.mtocm * (1 + z_form) / (1 + z_current) ** 2) ** 2

        """ADD A LINE HERE FOR THE LOOKUP TABLE"""

            k_de += bac.kappa_de * alpha[ri] * j_buffer
            k_di += bac.kappa_de * beta[ri] * j_buffer
            j21 += j_buffer


    return k_de, k_di, j21


# noinspection PyUnreachableCode
def sedcompute(mass, age, pos):

    file_path = '/Users/BAMAC/GitStuff/BaAstro/Erik/'
    file_name = file_path + 'SED_2'
    size = 54 # tied to file name, check if you change filename

    age_index = [str(x) for x in range(0,size)]
    sed_data = pd.read_csv(file_name,delim_whitespace=True,names=age_index,skiprows=1)

    age_data = []
    age_data = [np.append(age_data,sed_data[str(x)][0]) for x in range(0,size)]
    age_data_log = np.log10(age_data)

    wvlt_file = file_path+'wvlt_A_Erik'
    wvlt_data = pd.read_csv(file_path,names='A',skiprows=1) # wvlt in angstroms
    wvlt_A = np.array(wvlt_data['A'])

    n_gal = len(pos)

    # define the LW parameters that you will use
    kde = 0.
    kdi =0.
    jlw =0.

    # each galaxy here
    for gi in range(0,n_gal):

        # you could modify it to store an SED for *each* of the galaxies, but i dont see a point
        SED_gal = 0.

        for gi_i in range(0,len(mass[gi])):
            agetomatch = np.log10(age[gi][gi_i])
            age_match_id = min(enumerate(age_data_log), key=lambda age_log: abs(age_data_log[1] - agetomatch))

            SED_gal += sed_data[str(int(age_match_id[0]))][1:]*mass[gi][gi_i]/1e6

        # pass this SED to the alpha/beta/jlw module
        lw_quantities = np.array(lwlib.computerates(SED_gal,wvlt_A), pos[gi]) # check the pos bit
        kde += lw_quantities[0]
        kdi += lw_quantities[1]
        jlw += lw_quantities[2]

    return kde,kdi,jlw