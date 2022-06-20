#! /usr/bin/python3

import argparse, os, sys
import numpy as np
import pyrrl as rrl
import matplotlib.pyplot as plt

def plot_lines(lines):
    for line in lines:
        plt.axvline(x=line)
    plt.show()

def main(args):

    n_arr = np.arange(args.n_start, args.n_stop+1)
    freq = rrl.calc.recombination_line_frequency(n_arr,delta=args.delta_n,element=args.element,unit=1e6)

    for i,n in enumerate(n_arr):
        name = rrl.calc.recombination_line_name(n,delta=args.delta_n,element=args.element)
        print('{} n = {:3d} -> Line freq: {:8.3f} MHz'.format(name, n, freq[i]))

    if args.plot:
        plot_lines(freq)

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--n_start', type=int,default=100, help='energy level to start')
    parser.add_argument('--n_stop', type=int,default=200, help='energy level to stop')
    parser.add_argument('--delta_n', type=int,default=1, help='number of energy levels to jump')
    parser.add_argument('--element', type=str,default='H', help='H, He, or C')
    parser.add_argument('--plot', action='store_true', help='set to plot')


    args = parser.parse_args()
    main(args)
