## Install conda environment
conda env create -f /Net/Groups/BGI/work_3/sindbad/sindbad-preproc/fluxnet_workflow/sindbad_preprocess.yml

## Install the FLUXCOM repo:
cd /Net/Groups/BGI/work_3/biomass/BIOMASCAT/scripts/data/fluxcom
pip install -e .

## Run matlab wrapper
matlab /Net/Groups/BGI/work_3/sindbad/sindbad-preproc/matlab_wrapper/wrapper.m

## Note
The python script being called in the matlab wrapper is:
/Net/Groups/BGI/work_3/sindbad/sindbad-preproc/matlab_wrapper/cliff4matlab.py

