#!/bin/bash
#
# $1  path to Finngen file in gcloud
# $2  path to store GCTA-COJO formatted file, excludes empty RSIDs and MAF < 0.01 
#     MUST end in .ma

ma_file="${2%.*}.ma"

gcloud storage cat $1 \
    | gzip -dc \
    | awk '(NR == 1) {print "SNP\tA1\tA2\tfreq\tb\t se\tp\tn"} \
        (NR > 1 && $5 != "" && $11 >= 0.01) {sub(/,.*/, "", $5); print $5"\t"$3"\t"$4"\t" 1 - $11 "\t"$9"\t"$10"\t"$7"\tNA"}' \
    > $ma_file
