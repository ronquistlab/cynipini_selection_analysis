#!/bin/bash -l

#SBATCH -A snic2017-7-283
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 12:00:00
#SBATCH -J iqtree
#SBATCH -o Andricus-C60-275-18-50_at4.log

echo $(date)

../iqtree-1.6.12-Linux/bin/iqtree -s ../../data/A275r50-18-1.phy -nt 19 -m C60 -pre C60-A275-50-18-1 --runs 2 -bb 1000
../iqtree-1.6.12-Linux/bin/iqtree -s ../../data/A275r50-18-2.phy -nt 19 -m C60 -pre C60-A275-50-18-2 --runs 2 -bb 1000
../iqtree-1.6.12-Linux/bin/iqtree -s ../../data/A275r50-18-3.phy -nt 19 -m C60 -pre C60-A275-50-18-3 --runs 2 -bb 1000
../iqtree-1.6.12-Linux/bin/iqtree -s ../../data/A275r50-18-4.phy -nt 19 -m C60 -pre C60-A275-50-18-4 --runs 2 -bb 1000

wait
echo $(date)
