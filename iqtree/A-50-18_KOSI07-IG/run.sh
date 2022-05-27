#!/bin/bash -l

#SBATCH -A snic2017-7-283
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 24:00:00
#SBATCH -J iqtree
#SBATCH -o Andricus-50more-all18-codon_at1.log

echo $(date)

../iqtree-1.6.12-Linux/bin/iqtree -s ../data/50-18.phy -nt 19 -m KOSI07+G+I -pre A-50-18codGIb --runs 2 -st CODON -bb 1000

wait
echo $(date)

