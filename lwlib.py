import numpy as np
from scipy import integrate
import ba_constants as bac

def computerates(SED, wvlt, pos, z_current):
    """
    code to compute the rates and jlw output from an input SED
    note that jlwhere is defined at 13.6 eV and in units of 1e-21 erg/s/Hz/sr/cm^2

    :param SED: in erg/s/A
    :param wvlt: in Angstroms
    :param pos: the r_arr of the [halo - star] to be used for distance
    :return: kde, kdi, lwoutput
    """

    wvlt_13 = 911.64 # 13.6 ev
    wvlt_11 = 1107 # 11.2 ev
    wvlt_1 = 16313.717 # 0.76 eV
    kappa_de= 1.1e-10
    kappa_di= 1.38e-12

    dist = np.sqrt(sum(pos**2))*1e6*bac.pctom*bac.mtocm/(1 + z_current) # physical distance in cm

    # computing the jlw

    match_13 = min(enumerate(wvlt), key=lambda x: abs(x[1]-wvlt_13))
    SED_sel_13 = SED[match_13-1 : match_13+1]
    wvlt_sel_13 = wvlt[match_13-1 : match_13+1]
    jlw = abs(integrate.simps(SED_sel_13,wvlt_sel_13)) # er

    SED_norm = SED/SED[match_13] # SED normalised to it's own value at 13.6 eV in erg/s now
    freq = np.array([bac.c/(wvlt*1e-10)])

    # computing beta, kdi

    match_11 = min(enumerate(wvlt), key=lambda x: abs(x[1]-wvlt_11))
    match_12 = min(enumerate(wvlt), key=lambda x: abs(x[1]-wvlt_13))

    beta = abs(integrate.simps(SED_norm[match_11:match_12],freq[match_11:match_12])/
               (freq[match_11]-freq[match_12]))/kappa_de

    kdi = kappa_di * beta * jlw/dist**2 * bac.c21

    # computing alpha, kde

    # cross section for 1eV photons
    wvlt_cross = 1.6419*1e-6*1e10 # in angstroms
    # freq_cross = c/(wvlt_cross*1e-10)
    l_o= 1.6419 # in microns

    wvlt_cross_ids = np.array(np.where((wvlt >= wvlt_13) & (wvlt < wvlt_cross)))
    count_sel = len(wvlt_cross_ids)

    f_l_john = np.zeros(count_sel)
    l = wvlt[wvlt_cross_ids]*1e-10 # wvlt in m

    n = np.array([1.,2.,3.,4.,5.,6.])
    Cn = np.array([152.519, 49.534, -118.858, 92.536 , -34.195, 4.982])

    for li in range(0,count_sel):
         for ni in range(0,6):
            f_l_john[li] = f_l_john[li] + (Cn[ni]*(1/(l[li]*1e6) - 1./(l_o))**((n[ni]-1)/2.))

    sigma_pd =  1e-18*(l*1e6)^3 * (1/(l*1e6) - 1/l_o)**1.5 * f_l_john

    yfunc =4*3.14*SED[wvlt_cross_ids]*(wvlt[wvlt_cross_ids]**2 * 1e-10**2/bac.c * 1e10)/(SED[match_13]*wvlt[match_13]**2
            * 1e-10**2/bac.c *1e10)*1e-21 * sigma_pd/(bac.h*bac.joulestoerg*freq[wvlt_cross_ids])

    alpha = integrate.simps(SED[wvlt_cross_ids],freq[wvlt_cross_ids])/kappa_di

    kde = kappa_de * alpha * jlw/dist**2 * bac.c21

    return kde, kdi, jlw
