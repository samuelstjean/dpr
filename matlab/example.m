% Example on how to call the main realignment function
% We assume the files have been previously resampled in EDTI so that each subject is represented by one streamline
% per bundle and that the distance between each point (roughly the stepsize, after metric extraction) is the same.

% All datasets are in the same folder for simplicity, but this can likely
% be adapted by changing the paths

folder = 'datasets';
metric = 'TractFA';

% Everything down here should not need to be modified
files = dir(fullfile(folder,'*.mat'));
nfiles = length(files);
lengths = zeros(nfiles, 1);

for k = 1:nfiles
    filename = fullfile(folder, files(k).name);
    lengths(k) = length(load(filename, metric));
end

max_length = max(lengths);
bundles = zeros(nfiles, max_length);

% Read all the files for the metric we want to realign.
% metrics have all different length for every subject
% so we read them one by one and put them in a huge array we measured before
for k = 1:nfiles
    filename = fullfile(folder, files(k).name);
    values = load(filename, metric);
    values = cell2mat(getfield(values, metric));
    bundles(k, 1:length(values)) = values;
end

% We can now call the main realignment function with our pre zero-padded data
[realigned, final_shifts] = diffusion_profile_realignment(bundles, 'default');
