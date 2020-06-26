% Resample every line of a 2D matrix of metrics to have the same number of points
% accounting for padding during the 1D resampling

function [resampled] = resample_bundles_to_same(bundles, num_points)

    ndim = size(bundles);

    if length(ndim) ~= 2
        error('bundles needs to be a 2D array, but is %iD', length(ndim));
    end

    resampled = zeros(ndim(1), num_points);
    strip = zeros(ndim(1));

    for i = 1:ndim(1)
        strip(:) = bundles(i, :);
        len = size(strip);
        old_coord = linspace(1, len, len) / len;
        new_coord = linspace(1, len, num_points) / len;
        interpolated = interp1(new_coord, old_coord, strip);
        resampled(i, :) = interpolated;
    end
end
