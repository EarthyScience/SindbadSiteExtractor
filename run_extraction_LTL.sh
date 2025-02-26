#!/bin/bash
output_path='/Net/Groups/BGI/scratch/skoirala/v202312_wroasted/fluxNet_0.04_CLIFF/'
pythonpath='/Net/Groups/BGI/scratch/skoirala/con3/envs/sinproc004/bin/python'
configuration_file='/Net/Groups/BGI/scratch/skoirala/toolDev/sindbad-preproc/configuration/sindbad/sindbad_extraction_main.json'
# fn_version='FLUXNET2015'
fn_version='LaThuile'
time_scale='daily'
echo ${pythonpath} ${PWD}/run_extraction.py ${configuration_file} ${fn_version} ${time_scale} ${output_path}
${pythonpath} ${PWD}/run_extraction.py ${configuration_file} ${fn_version} ${time_scale} ${output_path}
exit 0