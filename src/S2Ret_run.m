% This is a small helper script to run the Leaf Area Index (LAI) model
% developed by Amin et al (2021, https://doi.org/10.1016/j.rse.2020.112168)
% on a directory with one or more Sentinel-2 scenes in .SAFE archive
% structure (L2A processing level).
%
% The LAI model produces LAI maps in 10m spatial resolution. It uses
% information from the scene classification layer (SCL) to mask out clouds,
% cloud shadows, 
function [] = S2Ret_run(model_path, in_dir_safe, out_dir)

% instanciate a new object
s2ret = S2Ret(model_path);

% run model
s2ret.retrieval(in_dir_safe, out_dir);

% terminate Matlab
exit;

end
