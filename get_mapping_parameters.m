function P = get_mapping_parameters

% GET MAPPING PARAMETERS FROM FILE
[file, path, filterindex] = ...
    uigetfile('*.mat', 'Select a Parameters File', 'MultiSelect', 'off','C:\data\sample.mat');
filename = [path,file];
% load(sprintf('%s',filename), 'P')
P = load(filename)

end