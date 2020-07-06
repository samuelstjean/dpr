% Resample every line of a 2D matrix of metrics to have the same number of points
% accounting for padding during the 1D resampling

function [resampled] = resample_bundles_to_same(bundles, num_points)

    ndim = size(bundles);

    if length(ndim) ~= 2
        error('bundles needs to be a 2D array, but is %iD', length(ndim));
    end

    resampled = zeros(ndim(1), num_points);
    strip = zeros(ndim(2), 1);

    for k = 1:ndim(1)
        strip(:) = bundles(k, :);
        len = length(strip);
        old_coord = linspace(1, len, len) / len;
        new_coord = linspace(1, len, num_points) / len;
        
        % We only interpolate on non nans, the padding is magically kept, don't ask why
        % It mimics np.interp though so that's good enough 
        valid = find(~isnan(strip));
        interpolated = interp1(old_coord(valid), strip(valid), new_coord);
        resampled(k, :) = interpolated;
    end
end
