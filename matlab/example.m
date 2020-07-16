% Example on how to call the main realignment function
% We assume the files have been previously resampled in EDTI so that each subject is represented by one streamline
% per bundle and that the distance between each point (roughly the stepsize, after metric extraction) is the same.

% All datasets are in the same folder for simplicity, but this can
% be adapted by changing the paths for your own data

folder = 'datasets';
metric = 'TractFA';
addpath('diffusion_profile_realignment')

disp('Now loading datasets')
tic;

% Everything down here should not need to be modified
files = dir(fullfile(folder,'*_ur_sp_ri.mat'));
nfiles = length(files);
lengths = zeros(nfiles, 1);

for k = 1:nfiles
    filename = fullfile(folder, files(k).name);
    file = load(filename, metric);
    lengths(k) = length(file.(metric){:});
end

max_length = max(lengths);
bundles = zeros(nfiles, max_length);

% Read all the files for the metric we want to realign.
% metrics have all different length for every subject
% so we read them one by one and put them in a huge array we preallocated before
for k = 1:nfiles
    filename = fullfile(folder, files(k).name);
    values = load(filename, metric);
    values = cell2mat(values.(metric));
    bundles(k, 1:length(values)) = values;
end

toc;

% We can now call the main realignment function with our pre zero-padded data
disp('Now performing diffusion profile realignment')
tic;
[realigned, final_shifts] = diffusion_profile_realignment(bundles, 'default');
toc;

disp('Now drawing the results')
tic;

% Note that the realigned values have an insane padding, so we truncate it for display purpose
% Here I chose to keep the length where at least 95% of the subjects are present
% Other options are 'shortest' (remove all padding, even if it truncates some data)
% or 'longest' (keep everything where at least one subject is present)
truncated = truncate_bundles(realigned, 95);

% We can also resample to N points (e.g. one per voxel) if we want.
% Here I pick 50 points because I want to
resampled = resample_bundles_to_same(truncated, 50);

% We can also compute statistics and then overlay them to draw fancy bundles
% This is really just an example and the test values themselves make no sense in this example context
% We just need some pvals to draw for the example.
% Do note that matlab ttest2 correctly handles nans, but for octave this is in the nan package
[h, pvals] = ttest2(truncated(1:50, :), truncated(51:100, :));


% Load the coordinates for the background gray bundle
coords = load('datasets/100307_DWI_b3000_af_left_AFD.mat', 'Tracts').Tracts;

% Load the coordinates for the truncated between ROIs blue bundle
coords_truncated = load('datasets/100307_DWI_b3000_af_left_AFD_ur.mat', 'Tracts').Tracts;

% Load the coordinates for the representative subject (roughly the middle) in green
coords_representative = load('datasets/100307_DWI_b3000_af_left_AFD_ur_sp.mat', 'Tracts').Tracts;

% Finally draw everything with that information.
% See the function for the default parameters and how to change them to achieve different effects such as labeling of the axes.
ax = [1, 3]; % Draw using the X and Z coordinates, so we effectively remove the Ys in the 2D projection
draw_fancy_graph(pvals, coords,  coords_truncated,  coords_representative, ax, 'default')

toc;
