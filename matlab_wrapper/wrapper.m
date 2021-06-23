pyversion('/User/homes/sbesnard/miniconda3/condabin/conda')
cd /Net/Groups/BGI/work_3/sindbad/sindbad-preproc/matlab_wrapper
flag = int32 (bitor (2, 8));
py.sys.setdlopenflags (flag);
py.importlib.import_module('cliff4matlab');
cubepath = '/Net/Groups/BGI/scratch/FLUXCOM/v0.01/site_cube';
version = 'FLUXNET2015';
site =  'DE-Hai';
varname = {['TA'],['VPD'], ['P'],['SW_IN'], ['NETRAD'],['LE'],['H']};
resolution = 'hourly';
output = struct(py.cliff4matlab.FLUXCOMdata4cliff(cubepath, version, site, varname, resolution));

