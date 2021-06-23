pyversion('/User/homes/sbesnard/miniconda3/condabin/conda')
cd /Net/Groups/BGI/work_3/sindbad/sindbad-preproc/matlab_wrapper
flag = int32 (bitor (2, 8));
py.sys.setdlopenflags (flag);
py.importlib.import_module('getSiteId');
cubepath = '/Net/Groups/BGI/scratch/FLUXCOM/v0.01/site_cube';
version = 'LaThuile';
varname = {['TA'], ['TA_QC'], ['P'], ['P_QC']};
resolution = 'hourly';
output = struct(py.getSiteId.getSiteID(cubepath, version, varname));

