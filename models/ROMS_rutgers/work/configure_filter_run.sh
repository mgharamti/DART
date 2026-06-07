#! /bin/bash

work_dir=/glade/derecho/scratch/${USER}/inacawo/DART_training/models/ROMS_rutgers/work

ens_size=$(grep -m 1 ens_size $work_dir/input.nml | sed 's/.*=//')

# cleanup
rm -f $work_dir/dart_log.* 

# Prior ensemble file list
ls -d /glade/derecho/scratch/gharamti/inacawo/roms_rst/roms_ens80/roms_mem_*nc > $work_dir/filter_input_list.txt

# Posterior ensemble 
post_dir=post_ens
mkdir -p $post_dir

for iens in $(seq -w 01 ${ens_size}); do
    echo $work_dir/$post_dir/roms_mem_out_${iens}.nc
done   > $work_dir/filter_output_list.txt
