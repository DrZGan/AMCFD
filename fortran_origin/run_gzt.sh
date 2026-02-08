#! /bin/bash
#PBS -l nodes=1:ppn=24
#PBS -l walltime=24:00:00
#PBS -N abqDB
#PBS -o cluster-log.log
#PBS -e cluster-error.log

echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR 
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo The following processors are allocated to this job:
echo `cat $PBS_NODEFILE`
NP=`wc -l < $PBS_NODEFILE`

./cluster_main

