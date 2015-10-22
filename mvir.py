"""
code that returns Mvir at a given z for a given Tvir
"""

def mvir(T,z):
    h =0.6711
    omega_m = 0.3125
    M = 1e8/h * (T/(1+z))**1.5 * (0.6*10/1.22/1.98e4)**1.5 * (18*3.14*3.14/178/omega_m)**0.5 #in solar masses
    return M
