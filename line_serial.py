#! /home/bliu/anaconda3/bin/python
#

import argparse, os, sys
import numpy as np
import matplotlib.pyplot as plt

def plot_lines(lines):
    for line in lines:
        plt.axvline(x=line)
    plt.show()
def calc_line_freq(n,type='alpha'):
    # n2 -> n1
    n2 = n+1
    n1 = n
    if type == 'alpha':
        freq = 6.58e15/(n*n*n)/1e6 # MHz
        #============================
        # need to be checked
        #RH = 109677.57
        #c = 2.99792458e10
        #Z = 1.0
        #freq = RH*c*Z*Z*(1/(n1*n1*n1) - 1/(n2*n2*n2)) /1e6 #MHz
        #============================
    return freq
def main(args):

    n0 = args.energy_level_start
    n1 = args.energy_level_stop

    lines = []

    for n in range(n0,n1):

        line_freq = calc_line_freq(n)
        #print('Line freq: {:8.3f} MHz at n = {:3d}'.format(line_freq,n))
        print('{:03d} {:08.3f}'.format(n,line_freq))
        lines.append(line_freq)

    plot_lines(lines)

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--energy_level_start', '-n0', type=int,default=100, help='start energy level')
    parser.add_argument('--energy_level_stop', '-n1', type=int,default=200, help='stop energy level')

    args = parser.parse_args()
    main(args)
