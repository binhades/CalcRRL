#! /home/bliu/anaconda3/bin/python
#
# hiirrl.py
# start-up script 

import argparse, os, sys
import numpy as np
import matplotlib.pyplot as plt

light_speed = 3e5 #km/s

def kms2khz(d_v,f0):
    # input: v should be in km/s

#    f = f0 * (1-v/light_speed)

    d_f = f0*d_v/light_speed

    return d_f

def khz2kms(d_f,f0):

    d_v = d_f/f0*light_speed

    return d_v

def calc_EM(n_e, D):
    # EM in unit pc cm-6
    # ne in unit cm-3
    # D  in unit pc

    EM = n_e*n_e*D

    return EM

def calc_taul(T_e, EM,delta_f):
    tau_L = 1.92e3*np.power(T_e,-2.5)*EM*(1./delta_f)

    return tau_L

def calc_Tl(T_e, tau_L):

    T_L = T_e * tau_L

    return T_L

def make_hii(n_e, D, d):
    return 0

def make_telescope():
    return 0

def ghz2khz(f):
    return f*1e6

def calc_departure():

    return 0

def main(args):

    f0 = ghz2khz(args.frequency)
    T_e = args.temperature

    M = 1836*9.11e-28

    lw_v = np.power(8.*np.log(2)*1.38e-16,0.5)*np.power(T_e/M,0.5) * 1e-5 # cm/s to km/s

    print('T_e = {:d}, Zero turb Line width: {:4.1f} km/s'.format(int(T_e),lw_v))
    print('T_e = {:d}, Equa turb Line width: {:4.1f} km/s'.format(int(T_e),lw_v*2))

    lw_v = []
    te = []

    for T_e in range(2000,15000,200):

        te.append(T_e)
        lw_v.append(np.power(8.*np.log(2)*1.38e-16,0.5)*np.power(T_e/M,0.5) * 1e-5)

    plt.plot(te, lw_v)
    plt.xlabel('Electron Temperature (K)')
    plt.ylabel('Line width without Turb (km/s)')
    plt.show()


    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--temperature', '-t', type=float,default=10000, help='The electron temperature in K')
    parser.add_argument('--frequency', '-f', type=float,default=5, help='Frequency in GHz')

    args = parser.parse_args()
    main(args)
