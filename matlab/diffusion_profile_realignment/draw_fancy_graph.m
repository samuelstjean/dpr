% Draws every streamline of a given bundle and overlay a map (such as pvals)
% This makes the fancy figures used in the manuscript where all subjects are plotted
% with a statistical test for visualization

function draw_fancy_graph(pval, coords, truncated_coords, average, ax, params)

    if isOctave()
       disp('Octave does not support all the drawing options, so for now we have to exit.')
       disp('This could change in a future version though, so try again later.')
       fprintf('At least it does not work (yet) in 5.2.0, and you have %s.\n', OCTAVE_VERSION);
       return
    end

    if strcmp(params, 'default')
        coord1_label = 'X';
        coord2_label = 'Y';
        pval_threshold = 1.;
        pval_cmap = 'hot';
        mean_fiber_cmap = [0.22953434 0.57685998 0.42976558]; % green
        bundle_cmap = [0.6093924 0.73212757 0.82249106]; %blue
        shadow_cmap = [0.75 0.75 0.75]; % gray;
        figtitle = 'p-values after realignment';
        draw_colorbar = true;
        filename = 'dpr_fancy_graph.png';
    else
        coord1_label = params.coord1_label;
        coord2_label = params.coord2_label;
        pval_threshold = params.pval_threshold;
        pval_cmap = params.pval_cmap;
        mean_fiber_cmap = params.mean_fiber_cmap;
        bundle_cmap = params.bundle_cmap;
        shadow_cmap = params.shadow_cmap;
        figtitle = params.title;
        draw_colorbar = params.draw_colorbar;
        filename = params.filename;
    end

    hold on;
    ax1 = ax(1);
    ax2 = ax(2);
    set(gca, 'Visible', 'off');

    % Draw the full length shadow bundle, each coordinates is a cell with a matrix of different size
    for idx = 1:length(coords)
        x = coords{idx}(:, ax1);
        z = coords{idx}(:, ax2);

        p1 = plot(x, z);%, 'LineWidth', 5);
        set(p1, 'Color', [shadow_cmap 0.1]) % transparency is the 4th channel
    end

    % Draw the original coord, but truncated between rois
    for idx = 1:length(truncated_coords)
        x = truncated_coords{idx}(:, ax1);
        z = truncated_coords{idx}(:, ax2);

        p2 = plot(x, z);%, 'LineWidth', 5);
        set(p2, 'Color', [bundle_cmap 0.3]) % transparency is the 4th channel
    end

    % Draw the mean coord
    x = average{1}(:, ax1);
    z = average{1}(:, ax2);
    plot(x, z, 'color', mean_fiber_cmap, 'LineWidth', 5) %, 'zorder', 5)

    % Again, we only interpolate on non nans, the padding is magically kept, don't ask why
    valid = find(~isnan(pval));
    new_coord = (1:length(x)) / length(x);
    old_coord = (1:length(pval)) / length(pval);

    % We resample the pvals because the coords may be at the tracking stepsize and the final results at the voxel resolution
    pval_resampled = interp1(old_coord(valid), pval(valid), new_coord);

    % This makes everything above the threshold invisible on the final graph
    pval_resampled(pval_resampled > pval_threshold) = nan;
    scatter(x, z, 'CData', pval_resampled, 'marker', 'o')
    colormap(pval_cmap)

    if draw_colorbar
        colorbar
        caxis([0, pval_threshold])
    end

    axis('equal')
    xlabel([coord1_label " coordinates (mm)"], 'fontsize', 12)
    ylabel([coord2_label " coordinates (mm)"], 'fontsize', 12)

    title(figtitle, 'fontsize', 20)
    hold off
    set(gca, 'Visible', 'on');
    saveas(gca, filename)
end


function [output] = isOctave ()
    output = exist('OCTAVE_VERSION', 'builtin') ~= 0;
end
