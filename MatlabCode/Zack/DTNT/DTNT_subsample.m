function subsampled = DTNT_subsample(target_nstims, fpar, domain, plots)
if nargin < 4 || isempty(plots)
    plots = false;
end

npoints = 50;
[L,M,S] = meshgrid(linspace(domain(1,1), domain(1,2), npoints), ...
    linspace(domain(2,1), domain(2,2), npoints), ...
    linspace(domain(3,1), domain(3,2), npoints));
v = sum(abs([L(:) M(:) S(:)] * reshape(fpar(2:end), [3 3])) .^ fpar(1), 2);

fv = isosurface(L, M, S, reshape(v, size(L)), 1);

if plots
    figure; hold on;
    patch('faces',fv.faces,'vertices',fv.vertices,'facecolor',[.6 .6 .6],'facealpha',.15,'edgecolor','none');
    set(gca,'Xlim',[min(L(:)) max(L(:))],'Ylim',[min(M(:)) max(M(:))],'View',[-10 30]);
end

[faces,verts] = reducepatch(fv.faces, fv.vertices, target_nstims/size(fv.faces, 1));

facerows = num2cell(faces(:,1:3), 2);
subsampled = cell2mat(cellfun(@(x) {mean(verts(x,:), 1)}, facerows));

if plots
    plot3(subsampled(:,1), subsampled(:,2), subsampled(:,3), 'k.');
end
