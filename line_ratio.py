#! /home/bliu/anaconda3/bin/python
# line_ratio.py

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

def main(args):

    f0 = ghz2khz(args.frequency0)
    f1 = ghz2khz(args.frequency1)

    delta_v0 = args.linewidth0
    delta_v1 = args.linewidth1

    delta_freq0 = kms2khz(delta_v0,f0)
    delta_freq1 = kms2khz(delta_v1,f1)

    ratio_0_1 = delta_freq1/delta_freq0

    print('Line ratio: {:3.1f} '.format(ratio_0_1))

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--frequency0', '-f0', type=float,default=1.3, help='the line rest frequency 0 (GHz)')
    parser.add_argument('--frequency1', '-f1', type=float,default=6.0, help='the line rest frequency 1 (GHz)')
    parser.add_argument('--linewidth0', '-lw0', type=float,default=25, help='the 0 linewidth in km/s')
    parser.add_argument('--linewidth1', '-lw1', type=float,default=25, help='the 1 linewidth in km/s')

    args = parser.parse_args()
    main(args)
