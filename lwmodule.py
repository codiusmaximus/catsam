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
import lwlib

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
