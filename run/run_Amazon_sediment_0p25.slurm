#!/bin/bash
#SBATCH -J Amazon_ctl_0p5
#SBATCH -p second
#SBATCH -N 1
#SBATCH -n 160
#SBATCH -w tms12
##SBATCH --ntasks-per-node=40
#SBATCH --mem=100G
#SBATCH -o logs/Amazon_ctl_0p5.o
#SBATCH -e logs/Amazon_ctl_0p5.e
#SBATCH -t 2-00:00:00
#SBATCH --exclusive

set -e

mkdir -p logs
{
  echo "job start: $(date)"
  echo "host: $(hostname)"
  echo "mpirun: $MPI"
  ldd ./colm.x | egrep 'libnetcdf|libnetcdff|libhdf5|libmpi|libmpifort|libiomp5|libifcore|not found' || true
} > logs/Amazon_sed_0p25_env.txt

#"$MPI" -np 240 ./mksrfdata.x "${NML}" > logs/Amazon_sed_0p25_srf.txt 2>&1
#mpirun -np 128 ./mksrfdata.x Amazon_sediment_0p25.nml > logs/Amazon_sed_0p25_srf.txt 2>&1
#mpirun -np 128 ./mkinidata.x Amazon_sediment_0p25.nml > logs/Amazon_sed_0p25_ini.txt 2>&1
mpirun -np 160 ./colm.x  Amazon_sediment_0p25.nml > logs/Amazon_sed_0p25_colm_test.txt 2>&1


exit $?
#echo job end time is `date`
