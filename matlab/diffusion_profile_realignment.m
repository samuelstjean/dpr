function [ys, final_shifts] = diffusion_profile_realignment(bundles, params)

    if exist(params, 'var') == 1
        percent = params.percent;
        padding = params.padding;
        epsilon = params.epsilon;
        remove_baseline = params.remove_baseline;
        whiten = params.whiten;
        rematch_outliers = params.rematch_outliers;
    %% Defaults
    else
        params.percent = 15;
        params.padding = 0;
        params.epsilon = 1e-5;
        params.remove_baseline = true;
        params.whiten = true;
        params.rematch_outliers = true;
    end

    % if we zero padded, put them to nan for now
    bundles(bundles == padding) = nan;

    npairs = size(bundles, 1);
    size_overlaps = sum(isfinite(bundles), 2);
    maxoverlaps = ceil(percent * size_overlaps / 100);
    shifts = zeros([npairs, npairs]);

    % pairs = filter_pairs(npairs);
    ffts = get_ffts(bundles, whiten, remove_baseline);
    ffta = zeros(size(ffts, 1));
    fftb = zeros(size(ffts, 1));

    for a = 1:npairs
        for b = 1:npairs
            ffta(:) = ffts(a);
            fftb(:) = ffts(b);
            shifts(a, b) = get_shift_from_fft(ffta, fftb);
        end
    end

    shifts(abs(shifts) < epsilon) = 0;  % force symmetry at ~zero shift

    % this is how many bundles are inside the overlapping threshold for bundle i
    sumstuff = sum(abs(shifts) < maxoverlaps, 2);

    % if we have more than one template, we pick the one with the min maximum shift
    % this does not seem to happen on real data, only on synthetic
    original_templater = find(sumstuff == max(sumstuff));

    if (size(original_templater) > 1)
        % abs(shifts) will be a symmetric matrix
        largest_shifts = max(abs(shifts), 2);
        largest_shifts = largest_shifts(original_templater);
        original_templater = original_templater(find(min(largest_shifts)));
    end

    if rematch_outliers
        condition = true;
        current_outliers = [];

        % while loop ends when outliers will be empty
        while condition

            % this is the indices of the outliers for the best alignment template bundle we just found
            outliers = abs(shifts(original_templater)) > maxoverlaps;
            outliers = find(outliers);

            candidates = abs(shifts(original_templater)) <= maxoverlaps;
            candidates = find(candidates);

            for idx = 1:size(outliers)
                outlier = outliers(idx);
                % find the match that minimize distance to original template without having very large displacement
                new_shifts = shifts(original_templater, :) + shifts(:, outlier);

                % remove candidate templates which are also outliers themselves
                [~, indexes] = sort(abs(new_shifts));
                indexes = setdiff(indexes, candidates);

                if ~isempty(indexes)
                    new_shifter = indexes(1);
                    % this is the shift required to realign bundle -> new_template -> original template
                    shifts(outlier, original_templater) = shifts(outlier, new_shifter) + shifts(new_shifter, original_templater);
                    shifts(original_templater, outlier) = shifts(original_templater, new_shifter) + shifts(new_shifter, outlier);
                end
            end
            % We have some outliers which do not overlap between the threshold
            % with the others, so we must break out of an infinite loop.
            % We use a set since order is not important and we replace them outside the loop.
            if isempty(setxor(current_outliers, outliers))
                rematch_outliers = false;
                break
            end

            current_outliers = outliers;
            condition = ~isempty(outliers);
        end
    else
        % no rematch outliers? put them to zero/nan then
        outliers = abs(shifts(original_templater)) > maxoverlaps;
        all_bundles = (1:npairs);
        outliers = all_bundles(outliers);

        if remove_outliers
            shifts(outliers) = nan;
            shifts(:, outliers) = nan;
        else
            shifts(outliers) = 0;
            shifts(:, outliers) = 0;
        end
        disp(['Possible outliers found : ', num2str(outliers)]);
    end

    % use the custom shift patterns
    ys = apply_shift(bundles, shifts(original_templater));
    final_shifts = shifts(original_templater);
end


% function [truncated] = truncate(bundles, truncmode='shortest')

%     if ndims(bundles) > 2:
%         error('Number of dimension must be equal or lower than 2, but was %i', ndims(bundles));
%     end

%     if isnumeric(truncmode)
%         if truncmode > 100 || truncmode < 1
%             error('truncmode must be between 1 and 100, but is %f', truncmode);
%         end
%         threshold = floor(truncmode * bundles.shape[0] / 100);
%     elseif truncmode == 'shortest'
%         threshold = bundles.shape[0];
%     elseif truncmode == 'longest';
%         threshold = 1;
%     else
%         error('Unrecognized truncation truncmode %s',truncmode);
%     end

%     indexes = sum(isfinite(bundles), 1) >= threshold;
%     truncated = bundles(:, indexes);
% end


% function [output] = filter_pairs(bundles)
%     how_much = size(bundles, 1);
%     output = cartesian(1:how_much, 1:how_much);
% end


% %% Stolen from https://stackoverflow.com/questions/9834254/cartesian-product-in-matlab
% function C = cartesian(varargin)
%     args = varargin;
%     n = nargin;

%     [F{1:n}] = ndgrid(args{:});

%     for i=n:-1:1
%         G(:,i) = F{i}(:);
%     end

%     C = unique(G , 'rows');
% end


function [shift] = get_shift_from_fft(x, y)

    spectras = crosscorrelation(x, y);

    peak = argmax(spectras);
    max_value = spectras(peak);
    % this works because we pad up to 2**N, so it's always even.
    half_length = length(spectras) / 2;

    m0 = spectras(peak - 1);
    m1 = spectras(peak + 1);
    p0 = peak - 1;
    p1 = peak + 1;

    extrapol = extrapolate([p0, peak, p1], [m0, max_value, m1]);
    shift = extrapol - half_length;
end


function [ffts] = get_ffts(bundles, whiten, remove_baseline)

    shape = size(bundles);
    finite = isfinite(bundles);

    if remove_baseline
        for i = 1:shape(1)
            bundle = bundles(i, finite(i));
            lin_poly = linspace(0., sqrt(max(bundle(:))), length(bundle));
            lin_fit = polyfit(lin_poly, bundle, 1);
            drift = polyval(lin_fit, lin_poly);
            bundles(i, finite(i)) = bundles(i, finite(i)) - drift;
        end
    end

    if whiten
        for i = 1:shape(1)
            mean_val = mean(bundles(i, finite(i)));
            std_val = std(bundles(i, finite(i)));
            bundles(i, finite(i)) = (bundles(i, finite(i)) - mean_val) / std_val;
        end
    end

    N = 2 * shape(2) - 1;
    pad = 2.^ceil(log2(N));
    ffts = zeros(shape(1), pad/2 + 1, 'complex128');

    for i = 1:shape(1)
        bundle = bundles(i, finite(i));
        ffts(i, :) = rfft(bundle, pad);
    end
end


function [result] = crosscorrelation(ffta, fftb)
    A = ffta;
    B = conj(fftb);

    result = irfft(A .* B);
    result = fftshift(result);
end


%% I'm too lazy to convert to fft and account for padding and everything, so here are rfft and irfft functions instead
%% all stolen from https://stackoverflow.com/questions/45778504/what-is-numpy-fft-rfft-and-numpy-fft-irfft-and-its-equivalent-code-in-matlab
function rfft = rfft(a, pad)
     ffta = fft(a, pad);
     rfft = ffta(1:(floor(length(ffta)/2)+1));
end


function irfft = irfft(x, even=true)
     n = 0; % the output length
     s = 0; % the variable that will hold the index of the highest
            % frequency below N/2, s = floor((n+1)/2)
     % Even case
     if (rem(length(x), 2)==0)
        n = 2 * (length(x) - 1);
        s = length(x) - 1;
    % Odd case
    else
        n = 2 * (length(x) - 1) + 1;
        s = length(x);
     end
     xn = zeros(1,n);
     xn(1:length(x)) = x;
     xn(length(x)+1:n) = conj(x(s:-1:2));
     irfft  = ifft(xn);
end


function [peak] = extrapolate(x, y)

    parabola = polyfit(x, y, 2);
    a, b, c = parabola;

    if a == 0 && b == 0 && c == 0
        peak = nan;
        % value = nan;
    else
        peak = -b / (2 * a);
        % value = a * peak.^2 + b * peak + c;
    end
end


function [shifted_bundles] = apply_shift(bundles, shifts)

    padding = nan;
    shape = size(bundles);
    shifted_bundles = zeros(shape(1), 3 * shape(2));

    for idx = 1:shape(1)
        bundle = bundles(idx,:);
        shift = shifts(idx,:);

        % if a shift is nan, it's an outlier anyway
        if isnan(shift)
            shift = 0;
        end

        integer_shift = fix(shift);
        invoxel_shift = shift - integer_shift;
        pad = zeros(shape(1)) * padding;

        new = [pad bundle pad];
        new = shift(new, integer_shift);
        new = interp1(new, invoxel_shift + (1:length(new)), 'linear', nan);

        shifted_bundles(idx, :) = new;
    end

end


% function [resampled] = resample_bundles_to_same(bundles, num_points)

%     ndim = size(bundles);

%     if length(ndim) ~= 2:
%         error('bundles needs to be a 2D array, but is %iD', length(ndim));

%     resampled = zeros(ndim(1), num_points);
%     strip = zeros(ndim(1));

%     for i = 1:ndim(1)
%         strip(:) = bundles(i, :);
%         len = size(strip);
%         old_coord = linspace(1, len, len) / len;
%         new_coord = linspace(1, len, num_points) / len;
%         interpolated = interp1(new_coord, old_coord, strip);
%         resampled(i, :) = interpolated;
%     end
% end
