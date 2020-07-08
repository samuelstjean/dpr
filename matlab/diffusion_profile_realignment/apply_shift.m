% Shifts each line of the 2D array bundles by the amount indicated by the vector shifts for each position
% This is already applied as the last line of diffusion_profile_realignment, but is also useful
% to draw new coordinates systems or p-values that would also need to be moved

function [shifted_bundles] = apply_shift(bundles, shifts)

    padding = nan;
    shape = size(bundles);
    shifted_bundles = zeros(shape(1), 3 * shape(2));

    for idx = 1:shape(1)
        bundle = bundles(idx, :);
        current_shift = shifts(1, idx);

        % if a shift is nan, it's an outlier anyway
        if isnan(current_shift)
            current_shift = 0;
        end

        integer_shift = fix(current_shift);
        invoxel_shift = current_shift - integer_shift;
        pad = zeros(1, shape(2)) * padding;

        new = [pad bundle pad];
        new = circshift(new, int32(integer_shift));
        xs = 1:length(new);
        new = interp1(xs, new, invoxel_shift + xs, 'linear');

        shifted_bundles(idx, :) = new(:);
    end
end
