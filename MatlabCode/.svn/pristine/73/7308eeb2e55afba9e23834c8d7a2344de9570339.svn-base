function plot_fixation_positions(filepaths, titles, limits, save_plots)
PLOTS_PER_FIGURE = 25;
r = ceil(sqrt(PLOTS_PER_FIGURE));
c = ceil(PLOTS_PER_FIGURE/r);

if save_plots
    if ~exist('export_fig.m', 'file')
        error('save_plots == true requires export_fig in the path');
    end
    save_folder = ['fixpos-' datestr(now, 30)];
    mkdir(save_folder);
end

plot_counter = 0;
fig_counter = 0;
for ii = 1:length(filepaths)
    curr_files = filepaths{ii};
    stros = cell(length(curr_files), 1);
    for jj = 1:length(curr_files)
        try
            stro = nex2stro(curr_files{jj});
            assert(stro.sum.paradigmID ~= 102, 'Ignoring SMurray file');
            assert(~isempty(stro.ras(:,strcmp(stro.sum.rasterCells(1,:), 'AD11'))), ...
                'No eye position data');
        catch e
            fprintf('Skipping %s (%s)\n', curr_files{jj}, strtrim(e.message));
            continue
        end
        stros{jj} = stro;
    end
    stros = [stros{:}];
    if isempty(stros), continue; end
    
    [medH, medV] = get_median_eye_pos(stros);
    
    if ~mod(plot_counter, PLOTS_PER_FIGURE)
        if save_plots
            if fig_counter > 0
                save_figure(fig_counter, save_folder);
            end
            figure(1905);
        else
            figure();
        end
        fig_counter = fig_counter + 1;
    end
    
    subplot(r, c, mod(plot_counter, PLOTS_PER_FIGURE)+1); hold on;
    scattercloud(medH,medV,40,.5,limits,[],pmkmp([],'linlhot'));
    plot([limits(1) mean(limits(3:4)); limits(2) mean(limits(3:4))], ...
        [mean(limits(1:2)) limits(3); mean(limits(1:2)) limits(4)], '-w');
    axis equal; axis(limits);
    set(gca, 'visible', 'off');
    set(gca, 'color', 'none');
    
    if ~isempty(titles), text(mean(limits(1:2)), limits(4), titles{ii}, 'horiz', 'center', 'vert', 'bottom'); end
    
    drawnow();
    plot_counter = plot_counter + 1;
end
if save_plots, save_figure(fig_counter, save_folder); end

function save_figure(fig_num, save_folder)
export_fig(1905, [save_folder filesep sprintf('eyepos%03d', fig_num)], '-m1.5', '-pdf', '-transparent');
close(1905);

% use the median eye position as a proxy for fixation position
function [medH, medV] = get_median_eye_pos(stros)
eyeH = cell(length(stros), 1);
eyeV = cell(length(stros), 1);
for ii = 1:length(stros)
    eyeH{ii} = stros(ii).ras(:,strcmp(stros(ii).sum.rasterCells(1,:), 'AD11'));
    eyeV{ii} = stros(ii).ras(:,strcmp(stros(ii).sum.rasterCells(1,:), 'AD12'));
end
medH = 4096 / 400 * cellfun(@median, cat(1, eyeH{:}));
medV = 4096 / 400 * cellfun(@median, cat(1, eyeV{:}));
