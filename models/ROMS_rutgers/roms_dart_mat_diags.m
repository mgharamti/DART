close all
clc

% Observation Space Diagnostics in MATLAB:
% ****************************************
%
% For a long time, DART used MATLAB for its observation space diagnostics. 
% So far in this training, we've used "pydartdiags" to explore DART's output and closely look at the obs_seq.final file. 
%
% Diagnostic Tools:
% =================
%
% In MATLAB, we provide several diagnostic functions and tools. 
% Some of these are similar to the ones we used in pydartdiags. 
% - plot_evolution           : plots the temporal evolution of any of the quantities above for each variable for specified levels. 
% - plot_profile             : plots the spatial and temporal average of any specified quantity as a function of height. 
% - plot_rmse_xxx_evolution  : same as plot_evolution but will overlay rmse on the same axis
% - plot_rmse_xxx_profile    : same as plot_profile with an overlay of rmse
% - plot_bias_xxx_profile    : same as plot_profile with an overlay of bias
% - two_experiments_evolution: same as plot_evolution.m but will overlay multiple (more than two, actually) experiments
% - two_experiments_profile  : same as plot_profile.m but will overlay multiple (more than two, actually) experiments 
% - plot_rank_histogram      : will create rank histograms for any variable that has that information present in obs_diag_output.nc
% - read_obs_netcdf          : reads a particular variable and copy from a netCDF-format observation sequence file and returns a single structure
% - plot_obs_netcdf          : creates a 3D scatterplot of the observation locations, color-coded to the observation values
% - plot_obs_netcdf_diffs    : creates a 3D scatterplot of the difference between two ‘copies’ of an observation
% - plot_wind_vectors        : reates a 2D ‘quiver’ plot of a wind field
% - link_obs                 : creates multiple figures that have linked attributes

% Load the obs_diag file
obdiag_loc500 = 'filter_output/obs_diag_output.nc';
obdiag_loc250 = 'filter_output_300_250/obs_diag_output.nc';
obdiag_loc100 = 'filter_output_300_100/obs_diag_output.nc';
obdiag_locno  = 'filter_output_300_no/obs_diag_output.nc';

obseq_loc500 = 'filter_output/obs_epoch_001.nc';
obseq_loc250 = 'filter_output_300_250/obs_epoch_001.nc';
obseq_loc100 = 'filter_output_300_100/obs_epoch_001.nc';
obseq_locno  = 'filter_output_300_no/obs_epoch_001.nc';

%% Compare Profiles
files    = {obdiag_locno, obdiag_loc500, obdiag_loc250, obdiag_loc100};
titles   = {'No V-Loc', 'V-Loc: 500m', 'V-Loc: 250m','V-Loc: 100m'};
varnames = {'FLOAT_TEMPERATURE', 'FLOAT_SALINITY'};
qtty     = 'rmse';
prpo     = 'posterior';

two_experiments_profile(files, titles, varnames, qtty, prpo, 'pause', true)



%% LINK OBS
global obsmat

% Read the observations 
fname         = obseq_loc250;
ObsTypeString = 'SATELLITE_BLENDED_SST';
ObsCopyString = 'observation';
CopyString    = 'prior ensemble mean';
QCString      = 'DART quality control';
region        = [85 150 -20 20 -Inf Inf];


link_obs(fname, ObsTypeString, ObsCopyString, CopyString, QCString, region)
