def pollute_rad(z, mvir, starmass):
    """

    :param z: redshift that you're at
    :param mvir: virial mass of the halo
    :param starmass: starmass array, optional
    :return: radius of pollution by metals in phy. kpc
    """

    vwind = 100 # km/s
    omegam=0.265
    omegal=0.735
    omegak=0.0
    h=0.71
    delc=178

    if starmass == 0 :

        omegamz= omegam*(1+z)**3 / ( omegam*(1+z)**3 + omegal + omegak*(1+z)**2 )

        Vc= 23.4 * (mvir*h/1e8)**(1./3.) * ((omegam/omegamz)*(delc/(18*3.14*3.14)))**(1./6.) * ((1.+ z)/10.)**0.5  #km/s

        rvir = 0.784* ( (mvir/(1e8 * 1./h))**(1./3) ) *( ( (omegam/omegamz)*(delc/(18*3.14*3.14)) )**(-1./3.) ) \
           *( ( (1 + z)/10.)**(-1) )* (1./h) # phy kpc


        pollrad = rvir * (vwind / Vc)# in phy kpc

    return pollrad