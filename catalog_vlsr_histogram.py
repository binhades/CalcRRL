#!/usr/bin/python3

# catalog_vlsr_histogram.py

import argparse, os
import matplotlib.pyplot as plt

def load_catalog_csv(file_csv):
    import csv
    with open(file_csv,mode='rt') as fcsv:
        reader = csv.DictReader(fcsv)
        return list(reader)

def main(args):

    if not os.path.isfile(args.file_csv):
        print("CSV file not found:",args.file_csv)
        return 0

    catalog = load_catalog_csv(args.file_csv)
    
    index = []
    vlsr = []
    green = [14,15]
    blue  = [2,4,8,18,19,20,21]
    red   = [1,3,5,6,7,9,10,11,12,13,16,17]

    colors = ['red','green','blue']


    for item in catalog:

        index.append(int(item['Index'] ))
        vlsr.append(float(item['VLSR']))

    n_bins = 10
    bins = [*range(35,95,5)]
    fig = plt.figure(figsize=(8,6))
    plt.grid(True,axis='y',linewidth=1,linestyle='dotted')
    plt.hist([[vlsr[i-1] for i in green],[vlsr[i-1] for i in red],[vlsr[i-1] for i in blue]],\
            bins,histtype='bar',color=['green','red','blue'],alpha=0.7,stacked=True,\
            edgecolor='black', linewidth=1.2)
    plt.xlabel('LSR Velocity (km/s)',fontsize='xx-large')
    plt.ylabel('Counts',fontsize='xx-large')
    plt.xticks(fontsize='xx-large')
    plt.yticks(fontsize='xx-large')

    fig.tight_layout()
    if args.file_out is not None:
        plt.savefig(args.file_out,dpi=300,bbox_inches='tight')
    plt.show()
    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file_csv', type=str, help='the source catalog')
    parser.add_argument('--file_out', type=str, help='the output figure name')
    args = parser.parse_args()
    main(args)

