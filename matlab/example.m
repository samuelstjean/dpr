% Example on how to call the main realignment function
% We assume the files have been previously resampled in EDTI so that each subject is represented by one streamline
% per bundle and that the distance between each point (roughly the stepsize, after metric extraction) is the same.

% All datasets are in the same folder for simplicity, but this can likely
% be adapted by changing the paths

folder = 'datasets';
metric = 'TractFA';
addpath('diffusion_profile_realignment')

disp('Now loading datasets')
tic;

% Everything down here should not need to be modified
files = dir(fullfile(folder,'*.mat'));
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
% so we read them one by one and put them in a huge array we measured before
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
writematrix(final_shifts,'final_shifts.csv');
writematrix(bundles,'bundles.csv');
writematrix(realigned,'realigned.csv')
toc;

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
% Do note that matlab correctly handles nans, but for octave this is in the nan package
[h,pvals] = ttest2(truncated(1:50, :), truncated(51:100, :));

% Since we also know the shift for each subject, we can also apply it elsewhere
% For example on the xyz coordinates so we can draw stuff now
% We use nans to initialize since we will also truncate them later

x = nan(nfiles, max_length);
y = nan(nfiles, max_length);
z = nan(nfiles, max_length);

% Read all the files for the xyz coordinates
for k = 1:nfiles
    filename = fullfile(folder, files(k).name);
    values = load(filename, 'Tracts');
    xyz = values.Tracts{:};
    len = size(xyz, 1);

    x(k, 1:len) = xyz(:, 1);
    y(k, 1:len) = xyz(:, 2);
    z(k, 1:len) = xyz(:, 3);
end

% Apply truncation on the coordinates
x_truncated = truncate_bundles(x, 95);
y_truncated = truncate_bundles(y, 95);
z_truncated = truncate_bundles(z, 95);

% Get the average pathway for drawing
x_average = mean(x_truncated, 1);
y_average = mean(y_truncated, 1);
z_average = mean(z_truncated, 1);

% Draw everything. This doesn't look super good in octave since it lacks transparency rendering though.
draw_fancy_graph(pvals, x, z, x_truncated, z_truncated, x_average, z_average, 'default')
