#!/bin/bash

dir=/home/bliu/data/fast/G034.39+00.22/
fdb=${dir}/db/hii-wise_20240710.db
fcsv=${dir}/catalog/wise_hii_V2.2.csv

db_table_fromCSV.py --file_csv ${fcsv}  --file_db ${fdb}
