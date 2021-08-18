#!/usr/bin/python3
# Filename: bn.py
# Aim: to calculate the bn

import argparse, os, time
import numpy as np
import scipy.constants as spc
from scipy import linalg

class Constants_IS():
    def __init__(self):
        import scipy.constants as spc
        self.e         = spc.e
        self.c         = spc.c 
        self.h         = spc.h 
        self.k         = spc.k 
        self.pi        = spc.pi
        self.m_e       = spc.m_e
        self.m_p       = spc.m_p
        self.hbar      = spc.hbar
        self.alpha     = spc.alpha
        self.epsilon_0 = spc.epsilon_0

        self.atomic_mass = spc.physical_constants['atomic mass constant'][0]
        self.amu_1hydrogen = 1.007825032231  
        self.amu_4helium = 4.002603254130
        self.amu_12carbon = 12

        self.R_inf  = self.m_e*self.e**4/(8*self.h**3*self.c*self.epsilon_0**2)
        self.R_H = self.R_inf * self.m_p/(self.m_p+self.m_e)
        self.E1_h = 1 * self.R_H * self.h * self.c  # E = Rhc, 13.6 eV in J
        self.a0 = 4*self.pi*self.epsilon_0*self.hbar**2/(self.m_e*self.e**2)

    def En_h(self,n):
        return self.E1_h/(n**2)

    def hv_kt(self,n,T):
        return self.En_h(n)/(self.k*T)

    def an(self,n):
        return self.a0*n**2

class Constants_CGS():
    def __init__(self):
        import scipy.constants as spc
        self.e    = 4.80320425e-10 # wikipedia
        self.c    = spc.c    * 1e2
        self.h    = spc.h    * 1e7
        self.k    = spc.k    * 1e7
        self.pi   = spc.pi
        self.m_e  = spc.m_e  * 1e3
        self.m_p  = spc.m_p  * 1e3
        self.hbar = spc.hbar * 1e7
        self.alpha     = spc.alpha

        self.atomic_mass = spc.physical_constants['atomic mass constant'][0]*1e3
        self.amu_1hydrogen = 1.007825032231  
        self.amu_4helium = 4.002603254130
        self.amu_12carbon = 12

        self.R_inf = 2*self.pi**2*self.m_e*self.e**4/(self.h**3*self.c)
        self.R_H = self.R_inf * self.m_p/(self.m_p+self.m_e)
        self.E1_h = 1 * self.R_H * self.h * self.c  # E = Rhc, 13.6 eV in erg
        self.a0 = self.hbar**2/(self.m_e*self.e**2)

    def En_h(self,n):
        return self.E1_h/(n**2)

    def hv_kt(self,n,T):
        return self.En_h(n)/(self.k*T)

    def an(self,n):
        return self.a0*n**2

def Te(t_c,t_l,fwhm,freq=1.3,p_he=0.1):
    """To calculate the RMS electron temperature.

    Parameters:
    -----------
    t_c: float or 1D numpy array.
         The continuum intensity.
    t_l: float or 1D numpy array.
         The line peak intensity.
    fwhm: float or 1D numpy array.
          The line width in km/s.
    freq: float.
          The line rest frequency in GHz.

    Returns:
    --------
    t_e: float or 1D numpy array.
    """

    c2l = t_c/(t_l*fwhm)

    t_e = np.power(7103.3* np.power(freq,1.1)*c2l/(1+p_he),0.87)

    return t_e

def luminosity(Sp, fwhm, v_res, D,unit=None):
    """To calculate the Luminosity of radio source

    Parameters:
    -----------
    Sp: float or 1D numpy array.
         The line peak flux intensity in Jy.
    D: float or 1D numpy array.
         The distance in kpc.
    fwhm: float or 1D numpy array.
          The line width in km/s.
    v_res: float.
          The velocity resolution in km/s.

    Returns:
    --------
    L: float or 1D numpy array.
        the luminosity in Watt/Hz.
        L = Sv * (4*np.pi) * D**2
    """
    L_solar = 3.86e26        # W -> L_solar
    D_meter = D*3.086e+19    # kpc -> m

    Sv = np.multiply(Sp, fwhm)/ v_res * 1e-26 # W m^-2 Hz^-1.

    L = np.multiply(Sv, np.power(D_meter,2)) * (4*np.pi)  

    if unit == 'solar':
        return L/L_solar # unit: L_solar
    else:
        return L

def eV2J(x):
    return  x * spc.e

def recombination_line_frequency(n_l, n_u=None,delta=1,unit=1,element='H'):

    if n_u is None:
        n_u = n_l+delta

    const = Constants_IS()

    if element == 'H':
        u = const.amu_1hydrogen
    elif element == 'He':
        u = const.amu_4helium
    elif element == 'C':
        u = const.amu_12carbon
    else:
        print('Element: {} is not known in the program'.format(element))
        return None
    m = const.atomic_mass * u
    R_m = const.R_inf * (m-const.m_e)/m
    freq = R_m * const.c * (1/n_l**2 - 1/n_u**2)

    return freq / unit # Hz or MHz

def recombination_line_frequency_H(n_l, n_u=None, delta=1,unit=1):

    if n_u is None:
        n_u = n_l+delta
    R_inf = spc.m_e*spc.e**4/(8*spc.h**3*spc.c*spc.epsilon_0**2)
    R_H = R_inf * spc.m_p/(spc.m_p+spc.m_e)
    freq = R_H * spc.c * ( 1/n_l**2 - 1/n_u**2) 

    return freq / unit # Hz or MHz

def recombination_line_width_thermal(Tk,freq=None,element='H'):
    # NRAO ERA Equ 7.35 & 7.34
    cgs = Constants_CGS()

    if element == 'H':
        u = cgs.amu_1hydrogen
    elif element == 'He':
        u = cgs.amu_4helium
    elif element == 'C':
        u = cgs.amu_12carbon
    else:
        print('Element: {} is not known in the program'.format(element))
        return None
    m = cgs.atomic_mass * u

    delta_v = (8*np.log(2)*cgs.k*Tk/m)**0.5
    if freq is not None:
        delta_f = delta_v/cgs.c*freq
        return delta_f
    else:
        return delta_v
def recombination_line_name(n,delta=1,element='H'):
    return '{element}{level}{type}'.format(element=element,level=n,type=chr(944+delta))

def recombine_rate_H(n,T=1e4,method='YJH'):
    from scipy.special import expn
    """To calculate the Hydrogen recombination rate

    TODO: YJH and SW methods give different results. need to check.

    Parameters:
    ----------
    n: primary quantum number
    T: LTE temperature

    Returns:
    -------
    rate
    """
    cgs = Constants_CGS()
    hv_kt = cgs.hv_kt(n,T)
    
    if method == 'YJH':
        # -----------------------------------------------------------
        #Equation 7.12: Chinese Radiative process in Astrophysics. You JH.

        c1 = 2**9 * cgs.pi**5 / (6 * cgs.pi)**(3/2) * cgs.e**10 \
                / (cgs.m_e**2 * cgs.c**3 * cgs.h**3)
        c2 = np.power((cgs.m_e/(cgs.k*T)),3/2) /(n**3)
        c3 = np.exp(hv_kt)

        if len(n) == 1:
            c4 = expn(1,hv_kt)
        else:
            c4 = np.empty(len(n))
            for i,t in enumerate(hv_kt):
                c4[i] = expn(1,t)

        rate = c1*c2*c3*c4 

    elif method == 'SK':
        ## -----------------------------------------------------------
        ## Equation 4.42, Physics and Chemistry of the ISM by Sun Kwok

        g_bf = 0.9
        c1 = 8/spc.pi**0.5 * 16/(3*3**0.5) * cgs.e**2 * cgs.h \
                / (cgs.m_e**2*cgs.c**3)
        c2 = (cgs.m_e/(2*cgs.k*T))**(3/2) * (cgs.En_h(n)/cgs.m_e)**2 /n**3 * g_bf
        c3 = np.exp(hv_kt)

        if len(n) == 1:
            c4 = expn(1,hv_kt)
        else:
            c4 = np.empty(len(n))
            for i,t in enumerate(hv_kt):
                c4[i] = expn(1,t)
                
        rate = c1*c2*c3*c4 
        ## -----------------------------------------------------------
    else:
        print('Method is not correct: YJH or SW.')
        rate = None

    return rate

def spontaneous_emission_rate_H(n,delta=1,unit='CGS'):
    n_u = n
    n_l = n_u - delta

    freq = recombination_line_frequency_H(n_l,delta=delta)

    if unit == 'CGS':
    # -------------------------------
        cgs = Constants_CGS()
        power = (2*cgs.pi)**4/3 * cgs.e**2/cgs.c**3 * freq**4 * cgs.an(n)**2
        A_ul = power/(cgs.h*freq)
    # -------------------------------
    elif unit =='IS':
        cis = Constants_IS()
        power = 1./(6*cis.pi*cis.epsilon_0) * cis.e**2/cis.c**3 \
                * (2*cis.pi*freq)**4 * cis.an(n)**2/2
        A_ul = power/(cis.h*freq)
    # -------------------------------
    else:
        print('UNIT not correct: CGS or IS')
        A_ul = None

    return A_ul

def energy_level_lifetime_H(n,delta=1):
    return 1./spontaneous_emission_rate_H(n,delta=1)

def Einstein_Coefficients_H(n,delta=1,unit='CGS'):

    n_u = n
    n_l = n_u - delta
    g_l = 2*n_l**2
    g_u = 2*n_u**2

    freq = recombination_line_frequency_H(n_l,delta=delta)
    A_ul = spontaneous_emission_rate_H(n,delta=delta,unit=unit)


    if unit == 'CGS':
        # NRAO ERA
        cgs = Constants_CGS()
        B_ul = A_ul * cgs.c**3/(8*cgs.pi*cgs.h*freq**3) 
    elif unit == 'IS':
        # wikipedia
        cis = Constants_IS()
        B_ul = A_ul * cis.c**3/(2*cis.h*freq**3)
    else:
        print('UNIT not correct: CGS or IS')
        B_ul = None

    B_lu = g_u/g_l * B_ul

    return A_ul, B_ul, B_lu

def Boltzmann_equation(n,delta=1,T=1e4):
    '''
    https://www.cv.nrao.edu/~sransom/web/Ch7.html
    Equation 7.43
    '''
    n_u = n
    n_l = n_u - delta

    g_l = 2*n_l**2
    g_u = 2*n_u**2

    R_inf = spc.m_e*spc.e**4/(8*spc.h**3*spc.c*spc.epsilon_0**2)
    R_H = R_inf * spc.m_p/(spc.m_p+spc.m_e)
    delta_E = R_H * spc.c * spc.h * ( 1/n_l**2 - 1/n_u**2) # hv

    r_Nu2Nl = g_u/g_l * np.exp(-1*delta_E/(spc.k*T))

    return r_Nu2Nl

def Saha_equation(n,Ne=1e3,Te=1e4,Np=None,bn=1):
    '''
    https://www.cv.nrao.edu/~sransom/web/Ch7.html
    Equation 7.92
    '''
    if Np is None:
        Np = Ne

    cgs = Constants_CGS()
    hv_kt = cgs.hv_kt(n,Te)
    gn = 2*n**2
    N_n = bn * Np * Ne * gn/2 *(cgs.h**2/(2*cgs.pi*cgs.m_e*cgs.k*Te))**(3/2)*np.exp(hv_kt)

    return N_n
def absorption_coefficient_cont(freq,Te=1e4,Ne=1e3,Ni=None,bn=1): #\kappa
    # YJH: EQ 7.42
    if Ni is None:
        Ni = Ne
    kc = 0.01 * Te**(-1.5)*freq**(-2)*Ne*Ni*np.log10(5e7*Te**(1.5)/freq)
    return kc

def absorption_coefficient_line(n,delta=1,Te=1e4,Ne=1e3,Np=None,bn=1): #\kappa
    # NRAO ESA: EQ 7.87
    n_l = n
    n_u = n_l + delta

    g_l = 2*n_l**2
    g_u = 2*n_u**2
    
    cgs = Constants_CGS()

    #freq = recombination_line_frequency(n_l,delta=delta,element='H')
    freq = recombination_line_frequency_H(n_l,delta=delta)
    delta_f = recombination_line_width_thermal(Te,freq=freq,element='H')
    phi_v0 = (np.log(2)/cgs.pi)**(0.5) * 2 / delta_f
    Nn = Saha_equation(n,Ne=Ne,Te=Te,Np=Np,bn=bn)
    A_ul = spontaneous_emission_rate_H(n,delta=delta,unit='CGS')

    kappa = cgs.c**2/(8*cgs.pi*freq**2) * g_u/g_l * Nn * A_ul * phi_v0 \
            * (1 - np.exp(-1*cgs.h*freq/(cgs.k*Te)))
    return kappa

def optical_depth(k,L):
    return k*L*3.08567758128e+18

def departure_coefficients(n):
    return 0

#==============================================================================
def main_bn():

    H = colrat(0,0,0,0)

    Te = 1000
    Ne = 1000
    n_min = 2
    icyc = 5

    N_low = 40
    N_high = 300
    Label = 'a'


    if n_min == 1:
        icase = 'Case A'
    elif n_min == 2:
        icase = 'Case B'


    icyc = max(1, icyc)
    npage = 1
    nline = 0
    nd = N_hi - N_low + 1
    itm = 1
    if (Te > 1000):
        itm = 3

    Te_12 = np.sqrt(Te)
    Te_32 = np.power(T3,3/2)

    for i in range(20,707):
        cx = 0
        cte = 15.778/Te
        arg = cte/(i**2)
        if arg < 165:
            cx = np.exp(-arg)
        cxp[i] = cx
    
# --------------------------------
# Call COLION(N, IONZ, T, QI)
# computes the collisional ionization rate from level N for ionis of effective charge ionz at the electron temperature T.
# When called with N=0, COLION computes and stores quantities which depend only upon temperature and effective charge. It is assumed that T and IONZ remain constant until the next call with N=0.

def colion(n, ionz, Te):

    cte = 15.778/Te

    cons = (ionz*ionz) * cte
    ind = np.arange(10,507)
    expx = np.exp(-1*cons/ind**2)

    x = cons/(n*n)
    dxp = expx[int(n)]

    if n > 507:
        dxp = np.exp(-1*x)

    if x <= 1:
        e1 = -0.57721566+x*(0.9999193+x*(-0.24991055+x*(0.5519968E-1+x*(-0.9760041E-2+x*0.107857E-2))))-np.log(x)
    else:
        e1 = ((0.250621+x*(2.334733+x))/(1.681534+x*(3.330657+x)))*dxp/x

    ei = dxp*(1.666667-.6666667*x)/x+e1*(.6666667*x-1.)-0.5*e1*e1/dxp
    qi = 5.444089*ei/np.power(Te,3/2)

    return qi


# --------------------------------
# Call RADCOL(T, MVAL, IC, NMIN)
# computes and stores in COMMON block RCRATS the total rates of radiative and collisional transitions from each energy level.
# NMIN is 1 for case A, 2 for case B.

def radcol(Te, Mval, ic, N_min):
    # len(Mval): 75
    # N_min: 1 for case A, 2 for case B.

    # ---------------------------------
    # Radiative cascade coefficients
    # Sum from N to all levels down to N_min (= 1 ro 2)

    for i in range(len(Mval)):
        N = Mval[i]
        tot = 0
        if N > N_min:
            K = N - 1
            for j in range(N_min, K):
                A = rad(j,N) # Einstein A, and Gaunt factor G.
                tot += A
        rad_tot[i] = tot
    # ---------------------------------
    # Collisional Rate totals on to level N
    # Summed from N_min to Infinity
        tot = 0
        for j in range(N_min, K):
            tot +=colrat(j,N,Te)
        L = N + 1
        N_max = N + 40
        for j in range(L, N_max):
            c = colrat(N,j,Te)
            cx = cxp(j)
            cxn = cxp(N)
            if cx < 1.0e-30 or cxn < 1.0e-30:
                Q = np.exp(15.778*(1/j**2-1/N**2)/Te)
            else:
                Q = cxn/cx
            c = c*(j/N)**2*Q
            tot += c

        col_tot[i] = tot

    return rad_tot, col_tot

def oscillator_strength(n_u, n_l):
    # Equation (1.15) from Menzel and Pekeris, 1935MNRAS..96...77M
    from scipy.special import hyp2f1

    D = 1 # for bond-bond transitions
    w_l = 2*n_l**2

    Delta = (hyp2f1(-1*n_u+1, -1*n_l, 1, -4*n_u*n_l/(n_u-n_l)**2))**2 - \
            (hyp2f1(-1*n_l+1, -1*n_u, 1, -4*n_u*n_l/(n_u-n_l)**2))**2
    f = 2**6/3* D/w_l * \
        np.abs(((n_u-n_l)/(n_u+n_l))**(2*(n_u+n_l)) / \
                (n_u**2*n_l**2*(1/n_l**2 - 1/n_u**2)**3) * \
               Delta/(n_u-n_l))

    return f

def gaunt_factor(n_u,n_l):
    # Equation (1.31 & 1.33) from Menzel and Pekeris, 1935MNRAS..96...77M
    w_l = 2*n_l**2
    f = oscillator_strength(n_u,n_l)
    f_ = 2**6/(3*np.sqrt(3)*np.pi*w_l) / \
            (1/n_l**2 - 1/n_u**2)**3 * \
            (1/n_u**3) * (1/n_l**3)
    g = f / f_
    return g

def spontaneous_coefficient(n_u,n_l):
    # Computes Einstein radiative rates A and Gaunt factors G for transitions from N_ to N
    # Equation (1.3 & 1.15) from Menzel and Pekeris, 1935MNRAS..96...77M
    cgs = Constants_CGS()
    f = oscillator_strength(n_u,n_l)
    v = recombination_line_frequency(n_l, n_u=n_u)
    A = 2*n_l**2/(2*n_u**2)*(8*np.pi**2*cgs.e**2*v**2)*f/(cgs.m_e*cgs.c**3)
    return A

def rad(N, N_):
    # Computes Einstein radiative rates A and Gaunt factors G for transitions from N_ to N
    # Equation (3.1 & 3.2) from Brocklehurst & Seaton, 1970MNRAS.148..417B
    n_u = N_
    n_l = N
    const = Constants_IS()
    g = gaunt_factor(n_u, n_l)
    # gamma = 15.7457e9 
    gamma = 16*const.alpha**4*const.c/(3*const.pi*np.sqrt(3)*const.a0)
    A = g * gamma /(n_l*n_u**3*(n_u**2-n_l**2))
    return A

def colrat(N,Np,Te):
    # Calculates rate of collisions from level N to higher level Np at electron temperature Te.
    # Set rate = 0 for Np-N > 40.
    # GPLR: Gee, Percival, Lodge and Richards, MNRAS,1976, 175, 209-215
    # GPLR need: 10^6/N^2 < Te << 3e9; Good for n > 40 and Te > 1000.

    if Te < 1e6/N**2:
        return colgl(N,Np,Te)

    beta = 1.58e5/Te
    beta1 = 1.4*np.sqrt(N*Np)

    if 0.2*(Np-N)/(N*Np) > 0.02:
        F1 = (1. - 0.3*(Np-N)/(N*Np))**(1+2*(Np-N))
    else:
        F1 = 1. -(1+2*(Np-N)) * 0.3*(Np-N)/(N*Np)


    A = (2.666667/(Np-N))*(Np/((Np-N)*N))**3*s23trm[int(Np-N)] * F1
    L = 0.85/beta
    L = np.log((1.+0.53*L**2*N*Np)/(1.+0.4*L))

    J1 = 1.333333 * A * L * (beta1/beta)/(beta1+beta)
    
    if 0.3*(Np-N)/(N*Np) > 0.02:
        F1 = (1. - 0.3*(Np-N)/(N*Np))**(1+2*(Np-N))
    else:
        F1 = 1. -(1+2*(Np-N)) * 0.3*(Np-N)/(N*Np)

    J2 = 1.777778 * F1 * (Np*(np.sqrt(2.-N**2/Np**2) +1.)/((N+Np)*(Np-N)))* \
            1*3*np.exp(-1*beta/beta1)/(beta/(1.-al18s4[int(Np-N)]))

    xi = 2./(N**2*(np.sqrt(2.-N**2/Np**2) -1))
    z = 0.75 * xi * (beta1+beta)

    J4 = 2./(z*(2.+z*(1.+np.exp(-1*z))))
    J3 = 0.25 * (N**2*xi/Np)**3 * J4 * np.log(1.+0.5*beta*xi)/(beta1+beta)
    rate = N**4*(J1 + J2 + J3)/np.power(Te,1.5)

    return rate

def colgl(N,N_,Te):
    # calculates collision rates from level N to higher level N_ at Te.
    # Uses Gauss-Laguerre integration of cross-sections (function cross) over maxwell distribution.
    # This function is used for values outside the region of validity of function colrat.

    ngl = 3
    xgl = np.array([.4157746,2.294280,6.289945])
    wgl = np.array([.7110930,.2785177,1.038926E-2])
    beta = 1.58e5/Te
    de = 1./N**2-1./N_**2

    r = 0
    for i in range(ngl):
        e = xgl[i]/beta+de
        r = r + wgl[i] * cross(N,N_,e)*e*beta

    r = r*6.21241e5*np.sqrt(Te)*N**2/N_**2

    return r

def cross(N,Np,E):

    # computes cross section for transition from level N to higher level Np due to collision with election of energy E.
    # The formula is valid for energies in the range 4/N**2 < E << 137**2
    # This program does not check that E is within this range
    # GPLR Theory: Gee, Percival, Lodge and Richards, MNRAS, 1976, 175, 209-215

    z = np.ones(506)+1

    al18s4 = np.log(18*z/(4.*z)
    s23trm = (0.184 - 0.04*z**(-0.6666667))

    def c2(x,y):
        return x*x*np.log(1.+0.6666667*x)/(y+y+1.5*x)
    
    if 0.2*(Np-N)/(N*Np) > 0.02:
        D = (1. -0.2*(Np-N)/(N*Np) )**(1+2*(Np-N))
    else:
        D = 1. - (1+2*(Np-N))*0.2*(Np-N)/(N*Np)

    A = (2.666667/(Np-N))*(Np/((Np-N)*N))**3*s23trm[int(Np-N)]*D

    if 1./(E*E*N*Np) < 150:
        D = np.exp(-1./(E*E*N*Np))
    else:
        D = 0
    L = np.log((1.+0.53*E*E*N*Np)/(1.+0.4*E))
    F = (1. - 0.3*(Np-N)*D/(N*Np))**(1+2*(Np-N))
    G = 0.5*(E*N*N/Np)**3 

    y = 1./(1. - D*al18s4[int(Np-N)])
    xp = 2./(E*N*N*(np.sqrt(2.-N**2/Np**2) +1.))
    xm = 2./(E*N*N*(np.sqrt(2.-N**2/Np**2) -1.))
    H = c2(xm,y) - c2(xp,y)

    cs = 8.797016e-17 * (N**4/E) * (A*D*L + F*G*H)

    return cs

# --------------------------------
# Call REDUCE(MVAL, IC, IR, SK)
# computes the values of the elements of the condensed matrix SK(of dimension IC*IC) by using lagrangian interpolation of order 2*(IR+1).

def reduce(M, ic, ir, sk):

    # Given a set of integers: M[it], it = 1 -> ic, such that:
    # 1) M[it+1] = M[it] + 1 for it < ia, where ia >=1, and
    # 2) (M[it+1] - M[it])) > 1 for it >= ia,
    # And given a function subprogram bk
    # which calculates the elements of a large M[ic]*M[ic] matrix,
    # This subroutine uses largrange interpolation of order 2(ir+1) to calculate a
    # smaller ic*ic matrix sk.
    # requires a function subprogram: phi
    # ir must be <= (ia - 1)

    lg = 2*(ir + 1)
    ib = ic - ir

    sk = np.empty([ic,ic])
    for i in range(ic):
        for j in range(ic):
            sk[i,j] = bk(M[i],M[j],i)
    for i in range(ic):
        ia = i
        if M[i+1] - M[i] > 1:
            if ia = ic:
                return sk
            if ia = ib:
                if ir == 0:
                    return sk
    for it in range(ia,ib-1):
        n1 = M[it] +1
        n2 = M[it+1] -1 
        for itau in range(1, lg):
            ind = it - ir -1 + itau
            iq[itau] = M[ind]
        for itau in range(1, lg):
            store1[itau] = phi(iq, lg, itau, iq[itau])
        for n in range(n1, n2):
            for itau in range(1, lg):
                store2[itau] = phi(iq, lg, itau, n)
            for si in range(1, ic):
                duckit = bk(M[si], n, si)
                for itau in range(1, lg):
                    fl = store2[itau]/store1[itau]
                    ind = it-ir-1+itau
                    sk[si, ind] = sk[si, ind] + duckit*fl

    if ir == 0:
        return sk

    for itau in range(1, lg):
        ind = ic-lg+itau
        iq[itau] = M[ind]

    for itau in range(1, lg):
        phitau=1./phi(iq,lg,itau,iq[itau])
        for it in range(ib, ic-1):
            n1 = M[it] +1
            n2 = M[it+1] -1
            for n in range(n1,n2):
                fl = phi(iq,lg,itau,n)*phiitau
                for si in range(ic):
                    ind = ic-lg+itau
                    sk[si,ind] = sk[si,ind]+bk(M[si],n,si)*fl
    return sk

def phi(iq,lg,itau,n):
    # lagrangian interpolation
    phi0 = 1.0
    for l in range(lg):
        if l == itau:
            continue
        phi0 = phi0 * (n-iq[l]) 

    return phi0

def bk(N,N_,si):

    # calls appropriate routines for calculation of atomic data for array sk
    # N = initial level, N_ = final level, 
    # si = subscript which identifies value of N in condensed matrix

    if N == N_:
        rt = colion(n,1,Te)
        bk = -1*radtot[si] - (coltot[si]+rt)*Ne
        if N > 20:
            bk = bk+cor(N,3)

    if N > N_:
        c = colrat(N, N_, Te)
        cx = cxp(N_)
        cxn = cxp(N)
        if cx < 1e-30 or cxn < 1e-30 or N_ > 707:
            Temp = np.exp(-cte*(1./N**2-1./N_**2))
        else:
            Temp = cxn/cx*(N_/N)**2

        if N > 20 and N_ == N+1:
            bk = bk + cor(N,1)

    if N < N_:
        c = colrat(N_,N,Te)
        bk = c*Ne
        if N_ == N-1 and N>20:
            bk = bk+cor(N,2)

    return bk

# --------------------------------
# Call RHS(CO, MVAL, IC)
# calculates the right-hand side, CO, of eq. (I.2.7) in the condensed system.

def rhs(co,Mval,ic):
    # equation 2.7 of Brocklehurst, MNRAS, 1970, 148, 417

    co = np.empty([ic])
    for i in range(1,ic):
        j = Mval[i]
        rt = colion(j,i,Te)
        alfa = recomb(1,Te,j)
        co[i] = -alfa*cxp[j]*np.power(Te,3/2)*0.24146879e16/(j*j)-rt*Ne
    return co

# --------------------------------
# Call JMD
# calculate corrections to the condensed matrix (SK in the program) to allow for levels above the highest level explicitly included in the equations.
def jmd(sk,co,Mval,ic,nfit,ival,itm):
    limit = 200
    ng = 24
    def sos(i,a):
        return np.power(np.sqrt(-1.*A),(2*i+itm)/np.log(-1.*A))
    for j in range(1,nfit):
        k = ival[j]
        aj = -1./float(mval[k])**2
        for i in range(1,nfit):
            afit[j,i] = sos(i,aj)

    
    #TODO: matinv(afit,nfit,az,0,d,irror,4,ipiv,ind)
    az = linalg.inv(afit)
    b = 1./float(mval[ic]+limit)**2
    a = -0.5*b
    nh = ng/2
    for k in range(1,nh):
        val[2k-1] = a+value[k]*b
        val[2k] = a - value[k]*b
    for k in range(1, limit):
        a = mval[ic]+k
        ac = -1./a**2
        for j in range(1,nfit):
            store1[k,j] = sos(j,ac)
    for k in range(1,ng):
        for j in range(1,nfit):
            inp = k+limit
            store[inp,j] = sos(j,val[k])
    kk = limit + ng
    for k in range(1, kk):
        for j in range(1,nfit):
            store3[k,j] = dmj1(k,j)
        store3[k,nfit+1]= dmj2(k)

    for j in range(1,ic):
        i = mval[j]
        for k in range(1,limit):
            kk = mval[ic]+k
            akk = kk
            store2[k] = bk(i,kk,0)
        for k in range(1,ng):
            ak = np.sqrt(-1./val[k])
            rtval[k] = ak**3*0.5.
            kk = ak
            inp = k + limit
            store2[inp] = bk(i,kk,0)
        aid = helpme(nfit+1,b,limit,rtval,store2,store3)
        co[j] = co[j] - aid
        for km in range(1, nfit):
            aid = helpme(km,b,limit,rtval,store2,store3)
            l = ival[km]
            sk[j,l] = sk[j,l] + aid

    return sk, az

def pol(k,km,limit,rtval,store2,store3):
    ind = k + limit
    return rtval[k]*store3[ind,km]*store2[ind]

def helpme(km,b,limit,rtval,store2,store3):
    s = 0.
    for k in range(1,limit):
        c = store2[k]*store3[k,km]
        s += c
    s = s - 0.5*c
    y =     0.6170614899993600 - 2.0*(pol( 1,km,limit,rtval,store2,store3)+pol( 2,km,limit,rtval,store2,store3))
    y = y + 0.1426569431446683 - 1.0*(pol( 3,km,limit,rtval,store2,store3)+pol( 4,km,limit,rtval,store2,store3))
    y = y + 0.2213871940870990 - 1.0*(pol( 5,km,limit,rtval,store2,store3)+pol( 6,km,limit,rtval,store2,store3))
    y = y + 0.2964929245771839 - 1.0*(pol( 7,km,limit,rtval,store2,store3)+pol( 8,km,limit,rtval,store2,store3))
    y = y + 0.3667324070554015 - 1.0*(pol( 9,km,limit,rtval,store2,store3)+pol(10,km,limit,rtval,store2,store3))
    y = y + 0.4309508076597664 - 1.0*(pol(11,km,limit,rtval,store2,store3)+pol(12,km,limit,rtval,store2,store3))
    y = y + 0.4880932605205694 - 1.0*(pol(13,km,limit,rtval,store2,store3)+pol(14,km,limit,rtval,store2,store3))
    y = y + 0.5372213505798282 - 1.0*(pol(15,km,limit,rtval,store2,store3)+pol(16,km,limit,rtval,store2,store3))
    y = y + 0.5775283402686280 - 1.0*(pol(17,km,limit,rtval,store2,store3)+pol(18,km,limit,rtval,store2,store3))
    y = y + 0.6083523646390170 - 1.0*(pol(19,km,limit,rtval,store2,store3)+pol(20,km,limit,rtval,store2,store3))
    y = y + 0.6291872817341415 - 1.0*(pol(21,km,limit,rtval,store2,store3)+pol(22,km,limit,rtval,store2,store3))
    y = y + 0.6396909767337608 - 1.0*(pol(23,km,limit,rtval,store2,store3)+pol(24,km,limit,rtval,store2,store3))

    return s+b*y
# --------------------------------
# Call MATINV
# is a standard routine which solves simultaneous equations or inverts a matrix
# LB: I use scipy.linalg.inv and scipy.linalg.solve to replace this function.
# --------------------------------
# Call INTERP(MVAL, CO, VAL, DVAL, IC, IR)
#   computes the bn values (VAL) and their derivatives (DVAL) for all n values
#   from the solutions CO at the IC condensed points defined by vector MVAL.
#   LAGRANGIAN INTERPOLATION OF ORDER 2*(IR+1) is used.
# LB: I use scipy.interpolate.lagrange to replace this function.
# --------------------------------

if __name__ == '__main__':


    parser = argparse.ArgumentParser()
    parser.add_argument('--Te',   type=float, default=10000, help='Electron Temperature in K, default 10000 K')
    parser.add_argument('--Ne',   type=float, default=10000, help='Electron Density in cm-3, default 10000 cm-3')
    parser.add_argument('--N',    type=int,   default=40,    help='Primary quantum level, default 40')
    parser.add_argument('--N_min', type=int,   default=10,    help='lower limit of Primary quantum level, default 10')
    parser.add_argument('--N_max', type=int,   default=507,   help='upper limit Primary quantum level, default 507')
    parser.add_argument('--N_l',  type=int,   default=10,   help='lower level, default 10')
    parser.add_argument('--N_u',  type=int,   default=11,   help='upper level, default 11')
    parser.add_argument('--case', type=int,   default=2,     help='1 for Case A, 2 for Case B, default 2')

    args = parser.parse_args()

    start_time = time.time()
    # --------------------------
    # processing
    #print(colion(40,1,args.Te))
    print(spontaneous_coefficient(args.N_u, args.N_l))
    print(rad(args.N_l, args.N_u))
    # --------------------------
    print("--- {:.3f} seconds ---".format(time.time() - start_time))
