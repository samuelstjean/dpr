% Function to truncate a 2D array of extracted metrics to a common number of points, removing nans
% If truncmode is a number between 1 and 100, refers ot the percentage of overlapping to keep
% If truncmode = 'shortest', only keep segments common to all bundles
% If truncmode = 'longest', only remove the nan padding and keep all segments

function [truncated] = truncate_bundles(bundles, truncmode)

    if ndims(bundles) > 2
        error('Number of dimension must be equal or lower than 2, but was %i', ndims(bundles));
    end

    if isnumeric(truncmode)
        if truncmode > 100 || truncmode < 1
            error('truncmode must be between 1 and 100, but is %f', truncmode);
        end
        threshold = floor(truncmode * length(bundles, 1) / 100);
    elseif strcmp(truncmode, 'shortest')
        threshold = length(bundle, 1);
    elseif strcmp(truncmode, 'longest')
        threshold = 1;
    else
        error('Unrecognized truncation truncmode %s',truncmode);
    end

    indexes = sum(isfinite(bundles), 1) >= threshold;
    truncated = bundles(:, indexes);
end
