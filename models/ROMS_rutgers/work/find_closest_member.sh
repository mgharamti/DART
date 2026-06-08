#!/bin/bash

#PBS -N CMT
#PBS -e cmt.err
#PBS -o cmt.out
#PBS -l select=5:ncpus=128
#PBS -l walltime=00:10:00
#PBS -A NHAP0012 
#PBS -q main
#PBS -l job_priority=premium

work_dir=/glade/derecho/scratch/${USER}/inacawo/DART_training/models/ROMS_rutgers/work

# Cleanup
rm -f $work_dir/dart_log.out       \
      $owrk_dir/closest_restart.nc \
      $work_dir/core 

# Run closest_member_tool
mpirun -n 500 $work_dir/closest_member_tool
