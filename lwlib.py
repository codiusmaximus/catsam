#LW modules go in here


#this one computes the LW flux in J_21 units, alpha and beta
def computeLW(input_dist, mass_stellar, age, sk):

    #load the SED tables
    popII_sed_block  = 
    #all ages loaded in here, extrapolation will happen here

    popIII_sed =
    
    #initialise the arrays
    n_points = len(input_dist)
    alpha = np.ndarray((1,n_points),float)
    beta = np.ndarray((1,n_points),float)
    LW_21 = np.ndarray((1,n_points),float)

    for i in range(0, n_points): 
        alpha[i],beta[i],LW_21[i] = computeparam(SED, wvlt, sk[i])

    return alpha, beta, LW_21

def computeparam(SED, wvlt, param):
    #CONSTANTS for stuff
    h_cgs =4.13e-15 # ev s
    jtoerg = 1e7
    evtoerg = 1/6.24e11
    c= 3e8 # m/s
    h = 6.626e-34 # J s
    kb = 1.3806e-23 # m^2 kg /s^2/K
    kb_cgs = 8.61e-5 # cm^2 kg/s^2 /K
    
    #define constants etc. here
    
    #compute LW here

    #compute alpha here
    
    #reaction constant
    k_const = 1.1d-10

    #base wvlt in microns corresponds to 0.76 eV
    l_o = 1.6419
    nu_o = c/(1.6419 * 1d-6) # in Hz


    #checking the xvalues
    sel_ids = where(xvalues*1d10 >= 911.6 and xvalues*1d6 <= l_o)
    if sel_ids :
      print 'check your wavelength array!'
      
   xvalues_sel = xvalues[sel_ids]
   yvalues_sel = yvalues[sel_ids]

   l = xvalues_sel
   bufferspec = yvalues_sel
   freq = c/l #Hz

   #the lyman limit, needed for integration limit and cross section--
   lylimit = closest(911.6,l*1d10)


    #compute LW here
