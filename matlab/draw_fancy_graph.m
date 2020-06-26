function [fig] = draw_fancy_graph(pval, coords1, coords2, truncated_coords1, truncated_coords2, average1, average2, params)

    if exist(params, 'var') == 1
        coord1_label = params.coord1_label
        coord2_label = params.coord2_label
        pval_threshold = params.pval_threshold
        pval_cmap = params.pval_cmap
        mean_fiber_cmap = params.mean_fiber_cmap
        bundle_cmap = params.bundle_cmap
        shadow_cmap = params.shadow_cmap
        title = params.title
        draw_colorbar = params.draw_colorbar
    %% Defaults
    else
        coord1_label = 'X'
        coord2_label = 'Y'
        pval_threshold = 1.
        pval_cmap = 'plt.cm.hot'
        mean_fiber_cmap = 'green'
        bundle_cmap = 'blue'
        shadow_cmap = 'gray'
        title = 'p-values after realignment'
        draw_colorbar = true
    end

    fig, ax = subplots(1, 1, sharex='col', sharey='row', figsize=(8,8))

    % Draw the full length shadow bundle
    for x, z in zip(coords1, coords2):
        x = x(isfinite(x))
        z = z(isfinite(z))
        ax.plot(x, z, color=shadow_cmap, alpha=0.1, zorder=1)

    % Draw the original coord, but truncated between rois
    for x, z in zip(truncated_coords1, truncated_coords2):
        x = x(isfinite(x))
        z = z(isfinite(z))
        ax.plot(x, z, color=bundle_cmap, alpha=0.3, zorder=2)

    % Draw the mean coord, yes we now use x and y because reasons
    x = average1
    y = average2
    ax.plot(x, y, color=mean_fiber_cmap, zorder=5)

    % We resample the pvals because the coords may be at the tracking stepsize and the final results at the voxel resolution
    size = len(x)
    pval_resampled = np.interp(np.arange(size), np.arange(len(pval)), pval)

    % This makes everything above the threshold invisible on the final graph
    pval_resampled[pval_resampled > pval_threshold] = np.nan
    cmap_axis = ax.scatter(x, y, c=pval_resampled, cmap=pval_cmap, zorder=10, marker='.', vmin=0, vmax=pval_threshold)

    if draw_colorbar
        colorbar(cmap_axis)
    end

    axis('equal')
    set_xlabel("{} coordinates (mm)".format(coord1_label), fontsize=12)
    set_ylabel("{} coordinates (mm)".format(coord2_label), fontsize=12)

    set_title(title, fontsize=20)
