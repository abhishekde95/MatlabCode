% Take thresholds from a set of models, weight those thresholds by the
% distance from the model to a point of interest (i.e., from each model's behavioral
% eccentricity to a neuron's receptive field), and return a set of points that
% live near the surface of a model that describes those weighted thresholds.
%
% This function doesn't do any model fitting. We hijack the face-vertex
% representation of an isosurface passing through the weighted thresholds. The
% mean of each vertex provides us a point in (L,M,TF) that "nicely" tiles the
% underlying model. This model is unknown since we haven't collected behavioral
% data at that location. Instead, we use this function and others in a
% semi-principled attempt to estimate that model.

function [subsampled,hp] = LMTF_subsample(target_nstims, fpars, domain, hplot)
if nargin < 5 || isempty(hplot)
    hplot = false;
    hp = [];
end

npoints = 50;
[L,M,TF] = meshgrid(linspace(domain(1,1), domain(1,2), npoints), ...
    linspace(domain(2,1), domain(2,2), npoints), ...
    logspace(log10(domain(3,1)), log10(domain(3,2)), npoints));

try
r = LMTF_thresh_from_model(L, M, TF, fpars);

LM_r = sqrt(L.^2+M.^2);
V = LM_r - r;
fv = isosurface(L, M, TF, V, 0);

if ishghandle(hplot) || hplot
    if islogical(hplot)
        figure; hold on; hplot = gca;
    end
    hp = patch('faces',fv.faces,'vertices',fv.vertices,'facealpha',.15,'edgecolor','none','parent',hplot);
    set(hplot,'Xlim',[min(L(:)) max(L(:))],'Ylim',[min(M(:)) max(M(:))],'View',[-10 30]);
end

% reducepatch destroyed the preferential sampling of small temporal frequencies
% (via the log-spaced lattice above). We've replaced it with a simpler heuristic
% that _seems_ to work well.
% [faces,verts] = reducepatch(fv.faces, fv.vertices, target_nstims/size(fv.faces, 1));

keep = round(linspace(1, size(fv.faces,1), target_nstims));
fv.faces = fv.faces(keep,:);

facerows = num2cell(fv.faces(:,1:3), 2);
subsampled = cell2mat(cellfun(@(x) {mean(fv.vertices(x,:), 1)}, facerows)); % use vertices instead of mean of vertices?
catch
    keyboard % GDLH
end
