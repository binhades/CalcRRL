#!/bin/bash

cd /home/bliu/data/fast/G034.39+00.22/
fdb=./db/hii-wise.db
fcsv=./catalog/wise_hii_V2.2.csv

db_table_fromCSV.py --file_csv ${fcsv}  --file_db ${fdb}
