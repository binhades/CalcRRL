#!/usr/bin/python3

# obs_calc_luminosity.py
# to Calculate luminosity of RRL emitting cloud. 

import argparse, os
from pyrrl.calc import luminosity

def load_catalog_csv(file_csv):
    import csv
    catalog = []
    with open(file_csv,mode='rt') as fcsv:
        reader = csv.DictReader(fcsv)
        for row in reader:
            catalog.append(row)
    return catalog 

def main(args):

    if not os.path.isfile(args.file_csv):
        print("CSV file not found:",args.file_csv)
        return 0

    catalog = load_catalog_csv(args.file_csv)

    for item in catalog:
        #print(item)
        d =    float(item['D_near'])
        sp =   float(item['Peak']  ) /1000. # mJy -> Jy
        vlsr = float(item['VLSR']  )
        fwhm = float(item['FWHM']  )
        v_res = 0.5

        L = luminosity(sp,fwhm,0.5,d,unit='solar')
        print('{} {:4.1f}'.format(item['Index'],L*1e12))

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file_csv', type=str, help='the source catalog')
    args = parser.parse_args()
    main(args)

