% Main function for diffusion profile realignment and its helper subfunctions

function [realigned, final_shifts] = diffusion_profile_realignment(bundles, params)

    if strcmp(params, 'default')
        percent = 15;
        padding = 0;
        epsilon = 1e-5;
        remove_baseline = true;
        whiten = true;
        rematch_outliers = true;
    else
        percent = params.percent;
        padding = params.padding;
        epsilon = params.epsilon;
        remove_baseline = params.remove_baseline;
        whiten = params.whiten;
        rematch_outliers = params.rematch_outliers;
    end

    % if we zero padded, put them to nan for now
    bundles(bundles == padding) = nan;

    npairs = size(bundles, 1);
    size_overlaps = sum(isfinite(bundles), 2);
    maxoverlaps = ceil(percent * size_overlaps / 100);
    shifts = zeros([npairs, npairs]);

    % pairs = filter_pairs(npairs);
    ffts = get_ffts(bundles, whiten, remove_baseline);
    ffta = zeros(size(ffts, 2), 1);
    fftb = zeros(size(ffts, 2), 1);

    for a = 1:npairs
        ffta(:) = ffts(a, :);
        for b = 1:npairs
            fftb(:) = ffts(b, :);
            shifts(a, b) = get_shift_from_fft(ffta, fftb);
        end
    end

    shifts(abs(shifts) < epsilon) = 0;  % force symmetry at ~zero shift

    % this is how many bundles are inside the overlapping threshold for bundle i
    sumstuff = sum(abs(shifts) < maxoverlaps, 2);

    % if we have more than one template, we pick the one with the min maximum shift
    % this does not seem to happen on real data, only on synthetic
    % original_templater = find(sumstuff == max(sumstuff));
    [~, original_templater] = max(sumstuff(:));

    if (length(original_templater) > 1)
        % abs(shifts) will be a symmetric matrix
        largest_shifts = max(abs(shifts), [], 2);
        largest_shifts = largest_shifts(original_templater);
        % argmin = find(min(largest_shifts(:)) == largest_shifts);
        [~, argmin] = min(largest_shifts(:));
        original_templater = original_templater(argmin);
    end

    if rematch_outliers
        condition = true;
        current_outliers = [];

        % while loop ends when outliers will be empty
        while condition

            % this is the indices of the outliers for the best alignment template bundle we just found
            abs_shifts = abs(shifts(original_templater, :));
            outliers = find(abs_shifts(:) > maxoverlaps);

            for idx = 1:length(outliers)
                outlier = outliers(idx);
                % find the match that minimize distance to original template without having very large displacement
                new_shifts = shifts(original_templater, :) + shifts(:, outlier)';

                % remove candidate templates which are also outliers themselves
                % We can use a setdiff in matlab since it preserves ordering
                [~, indexes] = sort(abs(new_shifts));
                indexes = setdiff(indexes, outliers, 'stable');

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
            if isequal(intersect(current_outliers, outliers), union(current_outliers, outliers))
                rematch_outliers = false;
                break
            end

            current_outliers = outliers;
            condition = ~isempty(outliers);
        end
    end

    % This check (instead of if/else) allows us to enter this condition if the previous matching
    % was looping around if an outlier could not be realigned at all
    if (~rematch_outliers)
        % no rematch outliers? put them to zero/nan then
        outliers = abs(shifts(original_templater)) > maxoverlaps;
        all_bundles = (1:npairs);
        outliers = all_bundles(outliers);

        if remove_outliers
            shifts(outliers, :) = nan;
            shifts(:, outliers) = nan;
        else
            shifts(outliers, :) = 0;
            shifts(:, outliers) = 0;
        end
        disp(['Possible outliers found : ', num2str(outliers)]);
    end

    % use the custom shift patterns
    realigned = apply_shift(bundles, shifts(original_templater, :));
    final_shifts = shifts(original_templater, :);
end


function [shift] = get_shift_from_fft(x, y)

    spectras = crosscorrelation(x, y);
    [max_value, peak] = max(spectras);

    % this works because we pad up to 2**N + 1 for rfft, so it's always even.
    half_length = floor(length(spectras) / 2);

    m0 = spectras(peak - 1);
    m1 = spectras(peak + 1);
    p0 = peak - 1;
    p1 = peak + 1;

    % In the matlab version, we preshift the values for scaling the parabola fit
    shift = extrapolate([p0, peak, p1] - (half_length + 1), [m0, max_value, m1]);
end


function [ffts] = get_ffts(bundles, whiten, remove_baseline)

    shape = size(bundles);
    finite = isfinite(bundles);

    if remove_baseline
        for k = 1:shape(1)
            idx = finite(k, :);
            bundle = bundles(k, idx);
            lin_poly = linspace(0., sqrt(max(bundle(:))), length(bundle));
            lin_fit = polyfit(lin_poly, bundle, 1);
            drift = polyval(lin_fit, lin_poly);
            bundles(k, idx) = bundles(k, idx) - drift;
        end
    end

    if whiten
        for k = 1:shape(1)
            idx = finite(k, :);
            mean_val = mean(bundles(k, idx));
            std_val = std(bundles(k, idx));
            bundles(k, idx) = (bundles(k, idx) - mean_val) / std_val;
        end
    end

    N = 2 * shape(2) - 1;
    pad = 2.^ceil(log2(N));

    % Initialize a complex valued array
    ffts = complex(zeros(shape(1), pad/2 + 1), 0);

    for k = 1:shape(1)
        idx = finite(k, :);
        bundle = bundles(k, idx);
        ffts(k, :) = rfft(bundle, pad);
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


function irfft = irfft(x)
    % n: the output length
    % s: the variable that will hold the index of the highest
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
    irfft = real(ifft(xn));
end


function [peak] = extrapolate(x, y)
    p = polyfit(x, y, 2);

    a = p(1);
    b = p(2);
    c = p(3);

    if (a == 0) && (b == 0) && (c == 0)
        peak = nan;
    else
        peak = -b / (2 * a);
    end
end

