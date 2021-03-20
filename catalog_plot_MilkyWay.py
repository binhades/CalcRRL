#!/usr/bin/python3

# catalog_plot_MilkyWay.py
# to plot the source in catalog onto the pane of MW. 

import argparse, os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

def load_catalog_csv(file_csv):
    import csv
    catalog = []
    with open(file_csv,mode='rt') as fcsv:
        reader = csv.DictReader(fcsv)
        for row in reader:
            catalog.append(row)
    return catalog 

def plot_MilkyWay(ax,R_sun=8.34,xlim=[-5,13],ylim=[-5,13],csize=4):
    c1 = patches.Circle((0,0), csize, facecolor='darkgray')
    c2 = patches.Circle((0,0), 2*csize, facecolor='silver')
    c3 = patches.Circle((0,0), 3*csize, facecolor='lightgray')

    ax.add_patch(c3)
    ax.add_patch(c2)
    ax.add_patch(c1)
    ax.set_aspect('equal')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.tick_params(labelsize='xx-large')
    ax.set_xlabel('X (kpc)',fontsize='xx-large')
    ax.set_ylabel('Y (kpc)',fontsize='xx-large')
    ax.plot(0,0,'*',markersize=30,color='k')
    ax.plot(0,R_sun,'o',markersize=15, markerfacecolor='none', markeredgewidth=2,color='yellow')
    ax.plot(0,R_sun,'.',color='yellow')
    ax.text(-0.3,-0.8,'GC',fontsize='xx-large')
    ax.text(0.4,R_sun,'Sun',fontsize='xx-large')

    return ax

def main(args):

    if not os.path.isfile(args.file_csv):
        print("CSV file not found:",args.file_csv)
        return 0

    catalog = load_catalog_csv(args.file_csv)
    R_sun = 8.34 # kpc

    green = [14,15]
    blue  = [2,4,8,18,19,20,21]
    red   = [1,3,5,6,7,9,10,11,12,13,16,17]
    colors= ['red','green','blue']

    fig,ax = plt.subplots(figsize=(6,8))
    plot_MilkyWay(ax,xlim=[-3,7],ylim=[-2,10],csize=3)

    for i, item in enumerate(catalog):
        d_n = float(item['D_near'])
        d_f = float(item['D_far'])
        l = np.deg2rad(float(item['GLon']))

        if i+1 in green:
            c = 'green'
        elif i+1 in blue:
            c = 'blue'
        elif i+1 in red:
            c = 'red'

        x_n = d_n*np.sin(l)
        x_f = d_f*np.sin(l)
        y_n = R_sun - d_n*np.cos(l)
        y_f = R_sun - d_f*np.cos(l)
        #ax.plot(x,y,'o',markersize=10,alpha=0.3,markerfacecolor=c,markeredgecolor='k')
        ax.plot(x_n,y_n,'+',markersize=16,color=c)
        ax.plot(x_f,y_f,'x',markersize=8,color=c)

    fig.tight_layout()
    if args.file_png is not None:
        plt.savefig(args.file_png,dpi=300,bbox_inches='tight')
    plt.show()
    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file_csv', type=str, help='the source catalog')
    parser.add_argument('--file_png', type=str, help='output fig file')
    args = parser.parse_args()
    main(args)

