#!/usr/bin/python3

# line_intensity.py
# start-up script 

import argparse, os, sys
import numpy as np
import astropy 

light_speed = 3e5 #km/s

def kms2khz(d_v,f0):
    # input: v should be in km/s

#    f = f0 * (1-v/light_speed)

    d_f = f0*d_v/light_speed

    return d_f

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

    delta_v = 20 #km/s
    T_e = 10000 # K
    n_e = 20 # cm-3
    D = 5 # pc

    delta_freq = kms2khz(delta_v,f0)

    EM = calc_EM(n_e, D)
    tau_L = calc_taul(T_e, EM, delta_freq)
    T_L = calc_Tl(T_e,tau_L)

    print('Line Temperature:{:6.3f} K'.format(T_L))

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--diameter', '-d', type=float,default=10, help='the diameter of the HII region')
    parser.add_argument('--density', '-ne', type=float,default=10, help='the diameter of the HII region')
    parser.add_argument('--frequency', '-f', type=float,default=5, help='the diameter of the HII region')

    args = parser.parse_args()
    main(args)
