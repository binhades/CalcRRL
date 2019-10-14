#
# hiirrl.py
# start-up script 

import argparse, os, sys
import numpy as np
import matplotlib.pyplot as plt

def plot_lines(lines):
    for line in lines:
        plt.axvline(x=line)
    plt.show()

def main(args):

    n0 = args.energy_level_start
    n1 = args.energy_level_stop

    lines = []

    for n in range(n0,n1):

        line_freq = calc_line_freq(n)
        print('Line freq: {:8.3f} MHz at n = {:3d}'.format(line_freq,n))
        lines.append(line_freq)

    plot_lines(lines)

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--energy_level_start', '-n0', type=float,default=100, help='start energy level')
    parser.add_argument('--energy_level_stop', '-n1', type=float,default=200, help='stop energy level')

    args = parser.parse_args()
    main(args)
