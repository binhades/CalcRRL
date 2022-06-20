#!/usr/bin/python3

# Filename: db_table_fromCSV.py
# Aim: to load the WISE HII catalog from CSV into SQLite3 database.

import argparse, sqlite3, csv

def table_format():
    table="wisehii"
    tab_hdr = '''
                WISE_Name,
                Catalog,
                GLong,
                GLat,
                Radius,
                HII_Region,
                Membership,
                VLSR,
                Author,
                VLSR_Mol,
                Molecule,
                KDAR,
                Dist,
                Err_Dist,
                Dist_Method,
                RGal,
                z,
                GLIMPSE_8um,
                WISE_12um,
                WISE_22um,
                MIPSGAL_24um,
                Hi_Gal_70um,
                Hi_Gal_160um,
                HRDS_3cm,
                MAGPIS_20cm,
                VGPS_21cm
    '''
    tab_command = "INSERT INTO wisehii VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"

    return table, tab_hdr, tab_command
def main(args):

    db_conn = sqlite3.connect(args.file_db)
    db_cur  = db_conn.cursor()
    table, tab_hdr, tab_command = table_format()
    db_cur.execute("CREATE TABLE if not exists {} ({})".format(table,tab_hdr))

    with open(args.file_csv, 'r') as fcsv:
        rows = csv.reader(fcsv,delimiter=',')
        next(rows,None)
        db_cur.executemany(tab_command, rows)

    db_cur.execute("SELECT * FROM wisehii")
    print(db_cur.fetchall())

    db_conn.commit()
    db_conn.close()

    return 0

#----------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--file_db',  type=str, default=':memory:',help='The sqlite database name')
    parser.add_argument('--file_csv', type=str, default='',required=True,help='The ds9 region file')
    args = parser.parse_args()
    main(args)
