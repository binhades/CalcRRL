#!/usr/bin/python3
# Filename: radiative recombination

import argparse, os, time
import numpy as np
import matplotlib.pyplot as plt
import pyrrl as rrl

def plot(n_arr,y_arr,labels,format='png',issave=False,isshow=True):
    fig = plt.figure(figsize=(10,6))
    ax0 = plt.subplot(111)
    mark = ['-','-.','--']

    for i,y in enumerate(y_arr):
        ax0.plot(n_arr,y,mark[i],label=labels[i])
    ax0.set_title('Rate Coefficient of Radiative Recombination')
    ax0.set_xlabel('Primary level $n$')
    ax0.set_ylabel(r'$\alpha_{n,rr}$ ($cm^{3}s^{-1}$)')
    ax0.set_yscale('log')
    ax0.grid(True)
    ax0.legend()

    plt.tight_layout()
    if issave:
        plt.savefig('radia_recomb.png',format=format,dpi=300,bbox_inches='tight')
    if isshow:
        plt.show()
    plt.close(fig)
    return 0

def main(args):

    n_arr = np.arange(args.N_max) + 1 # ind: 0 -> N_max-1; n: 1 -> N_max
    if args.compare == 'M':
        a_1   = np.zeros(args.N_max)
        a_2   = np.zeros(args.N_max)
        methods = ['S59','YJH']
        labels = methods
        for i, n in enumerate(n_arr):
            a_1[i] = rrl.atom.recombination_rate(n,args.Te,method=methods[0])
            a_2[i] = rrl.atom.recombination_rate(n,args.Te,method=methods[1])
        y_arr = [a_1,a_2]
    elif args.compare == 'T':
        a_1   = np.zeros(args.N_max)
        a_2   = np.zeros(args.N_max)
        a_3   = np.zeros(args.N_max)
        Te = [5000, 10000, 20000]
        labels = ['Te = 5e3 K', 'Te = 1e4 K', 'Te = 2e4 K']
        for i, n in enumerate(n_arr):
            a_1[i] = rrl.atom.recombination_rate(n,Te[0],method=args.method)
            a_2[i] = rrl.atom.recombination_rate(n,Te[1],method=args.method)
            a_3[i] = rrl.atom.recombination_rate(n,Te[2],method=args.method)
        y_arr = [a_1,a_2,a_3]

    plot(n_arr,y_arr,labels)

    return 0

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--N_max', type=int,   default=300, \
            help='upper limit Primary quantum level, default 300')
    parser.add_argument('--Te',   type=float, default=10000, \
            help='Electron Temperature in K, default 10000 K')
    parser.add_argument('--method', type=str,   default='S59', \
    parser.add_argument('--compare', type=str,   default='T', \
        help='Choose to compare the result for different: method(M) or Te(T)')
    args = parser.parse_args()
    start_time = time.time()
    main(args)
    print("--- {:.3f} seconds ---".format(time.time() - start_time))
