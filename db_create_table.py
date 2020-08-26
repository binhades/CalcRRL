#! /home/bliu/anaconda3/bin/python

# Filename: db_create_table.py
# Aim: to load the data cube, extract spectra.

import argparse, sqlite3
import astropy.units as u
from spectral_cube import SpectralCube
from fastspec import db
import regions

def db_add_source(conn,file_cata,ftype='reg'):

    if ftype == 'reg':
        reg_list = regions.read_ds9(file_cata)
        source_list = []
        for region in reg_list:
            source_list.append(db.region_to_source(region))
        db.add_to_source(conn,source_list)
    if ftype == 'db':
        #TODO
        1

    return 0

def db_add_spectrum(conn,file_cata,file_cube,ftype='reg'):

    cube = SpectralCube.read(file_cube)
    cube = cube.with_spectral_unit(u.km/u.s)
    cube = cube.spectral_slab(-200 * u.km/u.s, +200*u.km/u.s)
    velo = cube.spectral_axis.value

    if ftype == 'reg':
        reg_list = regions.read_ds9(file_cata)
        spec_list = []
        for region in reg_list:
            gname = db.get_gname(region=region)
            try:
                sub_cube = cube.subcube_from_regions([region])
            except ValueError as err:
                print(region)
                print(err)
                continue
            spec = sub_cube.mean(axis=(1,2))
            spec_list.append({'gname':gname,'spec':spec,'velo':velo})
        db.add_to_spectrum(conn,spec_list)

    if ftype == 'db':
        #TODO
        1
    return 0
def main(args):

    if not os.path.isfile(args.file_cube):
        print('File not found: ', args.file_cube)
        return 0

    db_conn = db.create_connection(args.file_db)
    with db_conn:
        db.create_table(db_conn,'Source')
        db.create_table(db_conn,'Spectrum')
        db.create_table(db_conn,'SpFitHa')

        if args.file_reg is not None:
            if not os.path.isfile(args.file_cata):
                print('File not found: ', args.file_cube)
                return 0
            db_add_source(db_conn,args.file_reg,ftype='reg')
            db_add_spectrum(db_conn,args.file_reg,args.file_cube,ftype='reg')

        if args.file_dbi is not None:
            if not os.path.isfile(args.file_reg):
                print('File not found: ', args.file_reg)
                return 0
            db_add_source(db_conn,args.file_reg,ftype='db')
            db_add_spectrum(db_conn,args.file_reg,args.file_cube,ftype='db')
    return 0

#----------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--file_cube',type=str, default='',required=True,help='The datacube file')
# TODO: Group
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--file_dbi', type=str, required=True,help='The input catalog db file')
    group.add_argument('--file_reg', type=str,  required=True,help='The ds9 region file')

    parser.add_argument('--file_db',  type=str, default=':memory:',help='The output sqlite database name')
    args = parser.parse_args()
    main(args)
