#!/bin/csh -f
mkdir data
rm -f data/top.x?f
max -p3 /u1b/remko/bury/data/*.adr data/top
rm -f xorms.grd
inimini data/top
