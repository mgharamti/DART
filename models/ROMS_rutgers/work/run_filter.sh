#!/bin/bash

#PBS -N roms_dart
#PBS -e roms_dart.err
#PBS -o roms_dart.out
#PBS -l select=40:ncpus=128
#PBS -l walltime=00:20:00
#PBS -A NHAP0012 
#PBS -q main
#PBS -l job_priority=premium

work_dir=/glade/derecho/scratch/${USER}/inacawo/DART_training/models/ROMS_rutgers/work

source $work_dir/../../../scripts/load_modules.sh

output_dir=$work_dir/filter_output

mkdir -p $output_dir

# Cleanup
rm -f $work_dir/dart_log.out           \
      $work_dir/filter_out.log         \
      $work_dir/core                   \
      $work_dir/roms_dart.*            \
      $work_dir/preassim_*.nc          \
      $work_dir/LargeInnov.txt         \
      $work_dir/filter_done            \
      $work_dir/filter_failed  

# Run filter
mpirun -n 5000 $work_dir/filter &> filter_out.log
filter_status=$?

# Check status and do diagnostics
if [ $filter_status -eq 0 ]; then 
   ncdiff $work_dir/output_mean.nc $work_dir/preassim_mean.nc \
          $work_dir/increment.nc
   touch  $work_dir/filter_done

   # Prepare observation space diagnostics
   cd $work_dir
   ./obs_diag
   ./obs_seq_to_netcdf

   # Clean up the directory
   mv output_mean.nc     $output_dir
   mv output_sd.nc       $output_dir
   mv preassim_*nc       $output_dir
   mv output_priorinf*nc $output_dir
   mv obs_seq.final      $output_dir
   mv increment.nc       $output_dir
   mv obs_diag_output.nc $output_dir
   mv obs_epoch_*nc      $output_dir
else 
   touch $work_dir/filter_failed
   echo  "Exit code: ${filter_status}" > $work_dir/filter_failed
fi


