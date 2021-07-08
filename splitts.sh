#!/bin/bash
#PBS -Pw85
#PBS -qnormal
#PBS -N splitts
#PBS -m ae
#PBS -M craig.arthur@ga.gov.au
#PBS -lwalltime=06:00:00
#PBS -lmem=32GB,ncpus=16,jobfs=4000MB
#PBS -W umask=0022
#PBS -joe
#PBS -o /home/547/cxa547/extremes/logs/splitts.out.log
#PBS -e /home/547/cxa547/extremes/logs/splitts.err.log
#PBS -lstorage=scratch/w85

# Run this with the following command line:
#
# qsub splitts.sh
#
# This will run the process for 41 years, starting 1979
module purge
module load pbs
module load dot

module load netcdf/4.6.3
module load nco/4.9.2
module load cdo/1.9.8
# Suppresses an error related to HDF5 libraries:
export HDF5_DISABLE_VERSION_CHECK=2
umask 0022

ECHO=/bin/echo

cd $HOME/extremes

DATAPATH=/g/data/rt52/era5/single-levels/reanalysis/i10fg
OUTPUTPATH=/scratch/w85/cxa547/reanalysis/timeseries
FINALPATH=/scratch/w85/cxa547/reanalysis/timeseries/station

STNLIST=`ncks --trd -s '%d ' -C -v stnid $OUTPUTPATH/ts.i10fg.2020.nc | grep "^[^stnid]"`

for YEAR in {1979..2020..1}; do
    STNNUM=0
    $ECHO "Splitting station data for $YEAR"
    for STN in $STNLIST; do
        ncks -d stnid,$STNNUM $OUTPUTPATH/ts.i10fg.$YEAR.nc -O $OUTPUTPATH/ts.i10fg.$YEAR.$STN.nc
        STNNUM=$(($STNNUM+1))
    done
done

for STN in $STNLIST; do
    $ECHO "Merging all data for station $STN"
    cdo mergetime $OUTPUTPATH/ts.i10fg.*.$STN.nc $FINALPATH/ts.$STN.nc
done

