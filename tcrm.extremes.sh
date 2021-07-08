#!/bin/bash
#PBS -Pw85
#PBS -qexpress
#PBS -N fitgpd
#PBS -m ae
#PBS -M craig.arthur@ga.gov.au
#PBS -lwalltime=2:00:00
#PBS -lmem=32GB,ncpus=32,jobfs=1000MB
#PBS -joe
#PBS -lstorage=gdata/w85+scratch/w85

# BEFORE RUNNING:
# Check there is a 'fit' folder under 'plots' in the designated output folder
# e.g. $OUTPUTPATH/plots/fit

module purge
module load pbs
module load dot

module load python3/3.7.4
module load netcdf/4.6.3
module load hdf5/1.10.5
module load geos/3.8.0
module load proj/6.2.1
module load gdal/3.0.2
module load openmpi/4.0.3

# Need to ensure we get the correct paths to access the local version of gdal bindings. 
# The module versions are compiled against Python3.6
export PYTHONPATH=/g/data/w85/.local/lib/python3.7/site-packages:$PYTHONPATH

# Add the local Python-based scripts to the path:
export PATH=/g/data/w85/.local/bin:$PATH

# Needs to be resolved, but this suppresses an error related to HDF5 libs
export HDF5_DISABLE_VERSION_CHECK=2


module list
DATE=`date +%Y%m%d%H%M`


# Add path to where TCRM is installed. 
SOFTWARE=/g/data/w85/software

# Add to the Python path we need to ensure we set the paths in the correct order
# to access the locally installed version of the GDAL bindings
export PYTHONPATH=$PYTHONPATH:$SOFTWARE/tcrm/master:$SOFTWARE/tcrm/master/Utilities
export PYTHONPATH=$PYTHONPATH:$HOME/extremes

# Suppresses an error related to HDF5 libraries:
export HDF5_DISABLE_VERSION_CHECK=2


# This can essentially be a normal TCRM configuration file
# (for a hazard simulation):

CONFIGFILE=$HOME/extremes/tcrmextremes.ini

cd $HOME/extremes


mpirun -mca mpi_warn_on_fork 0 -n $PBS_NCPUS python3 tcrmextremes.py -c $CONFIGFILE > tcrmextremes.stdout.$DATE 2>&1

