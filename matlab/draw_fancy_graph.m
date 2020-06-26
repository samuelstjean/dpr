% Draws every streamline of a given bundle and overlay a map (such as pvals)
% This makes the fancy figures used in the manuscript where all subjects are plotted
% with a statistical test for visualization

function [fig] = draw_fancy_graph(pval, coords1, coords2, truncated_coords1, truncated_coords2, average1, average2, params)

    if strcmp(params, 'default')
        coord1_label = 'X';
        coord2_label = 'Y';
        pval_threshold = 1.;
        pval_cmap = 'hot';
        mean_fiber_cmap = 'green';
        bundle_cmap = 'blue';
        shadow_cmap = [0.75 0.75 0.75] % gray;
        figtitle = 'p-values after realignment';
        draw_colorbar = true;
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
    end

    % Draw the full length shadow bundle
    for idx = 1:size(coords1, 1)
        x = coords1(idx);
        z = coords2(idx);
        x = x(isfinite(x));
        z = z(isfinite(z));
        plot(x, z, 'color', shadow_cmap, 'alpha', 0.1, 'zorder', 1);
    end

    % Draw the original coord, but truncated between rois
    for idx = 1:size(truncated_coords1, 1)
        x = truncated_coords1(idx);
        z = truncated_coords2(idx);
        x = x(isfinite(x));
        z = z(isfinite(z));
        plot(x, z, 'color', bundle_cmap, 'alpha', 0.3, 'zorder', 2)
    end

    % Draw the mean coord
    plot(average1, average2, 'color', mean_fiber_cmap, 'zorder', 5)

    % We resample the pvals because the coords may be at the tracking stepsize and the final results at the voxel resolution
    lenx = length(x);
    lenp = length(pval);
    pval_resampled = interp1(1:lenx, 1:lenp, pval);

    % This makes everything above the threshold invisible on the final graph
    pval_resampled(pval_resampled > pval_threshold) = nan;
    scatter(x, z, 10, pval_resampled, 'cmap', pval_cmap, 'zorder', 10, 'marker', '.', 'vmin', 0, 'vmax', pval_threshold)

    if draw_colorbar
        colorbar()
    end

    axis('equal')
    xlabel([coord1_label " coordinates (mm)"], 'fontsize', 12)
    ylabel([coord2_label " coordinates (mm)"], 'fontsize', 12)

    title(figtitle, 'fontsize', 20)
