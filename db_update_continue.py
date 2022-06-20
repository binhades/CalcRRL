#!/usr/bin/python3

# Filename : db_update_continue.py
# Aim: Load fits image, print the flux of given ds9 regions.

import argparse, sqlite3
from astropy import wcs
from astropy.io import fits
from astropy.coordinates import Angle, SkyCoord
import numpy as np
import regions
from fastspec import db

def regions_from_sources(source_list):
    region_list = []
    for source in source_list:
        l = source['Glon']
        b = source['Glat']
        r = Angle(source['radius'],'arcsec')
        center = SkyCoord(l,b,unit='deg',frame='galactic')
        region_list.append(regions.CircleSkyRegion(center, r))
    return region_list
        
def main(args):
    hdul = fits.open(args.file_fits,'readonly')
    hdr = hdul[0].header
    img = np.squeeze(hdul[0].data)
    hdul.close()

    db_conn = db.create_connection(args.file_db)
    #db.create_table(db_conn,'PhyPara')
    #db.add_new_column(db_conn,'PhyPara','Tb_vgps','FLOAT')
    #db.add_new_column(db_conn,'PhyPara','Tb_ao','FLOAT')

    source_list = db.read_from_table(db_conn,'Source')
    img_wcs = wcs.WCS(hdr)

    reg_list = regions_from_sources(source_list)

    with db_conn:

        for i, region in enumerate(reg_list):
            print(i)
            reg_pix = region.to_pixel(img_wcs.celestial)
            mask = reg_pix.to_mask()
            cutout = mask.cutout(img, fill_value=0.0,copy=True)
            if cutout is None:
                print('cutout is None')
                continue
            gname = 'G{:07.3f}{:+07.3f}'.format(region.center.l.value,region.center.b.value)

            cutout = cutout.astype(np.float64)

#            data = {'gname':gname,'Tb_vgps':np.mean(cutout)}
#            db.add_to_phypara(db_conn,data)

#            data = {'gname':gname,'Tb_ao':np.mean(cutout)}
#            db.update_phypara(db_conn,data)
            print(data)

    return 0

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('file_fits', type=str, help='the input file name')
    parser.add_argument('--file_db', type=str, default=':memort:', help='data base name')
    args = parser.parse_args()
    main(args)
