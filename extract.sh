#!/bin/bash
#PBS -Pw85
#PBS -qnormal
#PBS -N extractts
#PBS -m ae
#PBS -M craig.arthur@ga.gov.au
#PBS -lwalltime=06:00:00
#PBS -lmem=32GB,ncpus=16,jobfs=4000MB
#PBS -W umask=0022
#PBS -joe
#PBS -o /home/547/cxa547/extremes/logs/extractts.out.log
#PBS -e /home/547/cxa547/extremes/logs/extractts.err.log
#PBS -lstorage=gdata/w85+gdata/rt52+scratch/w85

# Run this with the following command line:
#
# qsub -v NJOBS=41,YEAR=1979 extract.sh
#
# This will run the process for 41 years, starting 1979

module purge
module load pbs
module load dot

module load python3/3.7.4
module load netcdf/4.6.3
export PYTHONPATH=/g/data/w85/.local/lib/python3.7/site-packages:$PYTHONPATH
export PYTHONPATH=$PYTHONPATH:$HOME/extremes
# Suppresses an error related to HDF5 libraries:
export HDF5_DISABLE_VERSION_CHECK=2
umask 0022

ECHO=/bin/echo

cd $HOME/extremes

DATAPATH=/g/data/rt52/era5/single-levels/reanalysis/i10fg
OUTPUTPATH=/scratch/w85/cxa547/reanalysis/timeseries
FINALPATH=/scratch/w85/cxa547/reanalysis/timeseries/station

for YEAR in {1979..2020..1}; do
    $ECHO "Extracting data for $YEAR"
    python3 extract.py -p $DATAPATH -y $YEAR -o $OUTPUTPATH > $HOME/extremes/logs/extract.stdout.$YEAR 2>&1
    if [ -f STOP_SEQUENCE ] ; then
        $ECHO  "Terminating sequence at $YEAR"
        exit 0
    fi
done


