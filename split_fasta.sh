#!/bin/bash
mkdir new
cd new
csplit -s -z /home/yp11/Desktop/genomes/hg19/chromFa_all/hg19.fa '/>/' '{*}'
for i in xx* ; do \
  n=$(sed 's/>// ; s/ .*// ; 1q' "$i") ; \
  mv "$i" "$n.fa" ; \
 done
