#!/bin/bash

# Estimate LD scores from the hap550 data set, separately for each
# chromosome. The LD scores are stored in files ldsc_*.l2.ldscore.gz.

# These settings should work on the PPS cluster (for an updated node
# such as spudling21). Please adjust these settings as needed to suit
# your computing environment.
PLINK=/mnt/gluster/data/software/plink-1.90b5.2/plink
PYTHON=/opt/anaconda/anaconda2/bin/python
LDSC=/mnt/gluster/data/software/ldsc-1.0.0/ldsc.py

# Repeat for chromosomes 1 through 22.
for (( i=1 ; i <= 22 ; i++ )); do
  echo "chromosome ${i}"

  # Extract the genotypes on the chromosome.
  $PLINK --bfile hap550 --chr $i --make-bed > plink.out

  # Estimate univariate LD scores for all SNPs on the chromosome.
  $PYTHON $LDSC --bfile plink --l2 --ld-wind-kb 100 \
    --out ldsc_${i} > ldsc_${i}.out
done
