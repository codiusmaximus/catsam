import numpy as np


def checkpopiii_lw(mass,LW):
    """
    Code tells you if halo can host Pop III or not depending on the LW specific intensity and halo mass

    input:
    mass = virial mass of the halo in Msun
    LW = net LW spec intensity from all stars(units of 1e-21 erg/s/cm^2/Hz/sr)
    this might get tricky later on, leave as is for now*

    output:
    halo mass required for Pop III formation given the irraditaing spec inten.
    """""

    #Machachek et al 2001 * 4 (as reported by Oshea and Norman 08)
    mass_req = 1.25e5 + 8.7e5*(4*3.14*LW)**0.47
    if mass_req <= mass
        make_popiii = "yes"
    else:
        make_popiii = "no"
    return make_popiii

def makepopiii(num):
    """
    code to generate Pop III stars randomly between 100-1000 solar masses

    :param num: number of stars to be made
    :return: array of size num, with the masses of stars generated
    """

    min_iii = 100
    max_iii = 1000
    star_masses = np.array(np.random.random(num)*(max_iii-min_iii) + min_iii)

    return star_masses