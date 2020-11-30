% This function was a working sandbox for a clever but ultimately non-functional
% way of picking new stimuli for the LMTF paradigm. This script is left here for
% historical and pedagogical reasons of what _not_ to do.

function LMTF_simulation_old()
global DEBUG done_with_initial_ts
done_with_initial_ts = false;
DEBUG = true;
PLOT_TRIANGLES = true;
PLOT_NORMALS = false;

% Trying to sample a "twisty funnel" by a NeuroThresh type design.
% First making a Twisty funnel.
NPOINTS = 30;
minTF = 1;
maxTF = 25;
lum = linspace(-.2,.2,NPOINTS);
rg = linspace(-.2,.2,NPOINTS);
tf = linspace(minTF,maxTF,NPOINTS);
[x,y,z] = meshgrid(lum,rg,tf);

% Crude approximation of Charlie's TF sensitivity functions
pplum = pchip([0 .02 1 5 25],1./[10 10 10 10 10]);
pprg = pchip([0 .02 1 5 25],1./[10 10 10 10 10]);
lumthreshs = fnval(pplum,tf);
rgthreshs = fnval(pprg,tf);
% figure; axes; hold on;
% plot(log10(tf),1./rgthreshs,'r-');
% plot(log10(tf),1./lumthreshs,'k-');

lumthreshmat = permute(repmat(repmat(lumthreshs,NPOINTS,1),[1 1 NPOINTS]),[3 1 2]);
rgthreshmat = permute(repmat(repmat(rgthreshs,NPOINTS,1),[1 1 NPOINTS]),[3 1 2]);
fval = (x./lumthreshmat).^2+(y./rgthreshmat).^2;

% Now trying to rotate it so that the axes are L and M, not LUM and RG
rotmat = [1/sqrt(2) -1/sqrt(2); 1/sqrt(2) 1/sqrt(2)];
newxy = [x(:) y(:)]*rotmat';

figure; axes; hold on;
fv = isosurface(reshape(newxy(:,1),size(x)),reshape(newxy(:,2),size(y)),z,fval,1);
patch(fv,'facecolor',[.6 .6 .6],'facealpha',.15,'edgecolor','none','hittest','off');
set(gca,'Xlim',[min(lum) max(lum)],'Ylim',[min(rg) max(rg)],'View',[-10 30]);

% OK, now we're going to try an simulate a NeuroThresh-style experiment
% on our new twisty funnel

N_INTERLEAVED_TS = 4; % test at most this number of TS per round

% starting with four stimuli (in L, M, and TF coordinates)
starting_stimuli = [1 1 minTF; 1 -1 minTF; 1 1 maxTF; 1 -1 maxTF];
N_INITIAL_TS = size(starting_stimuli, 1);
trialspecs = setup_starting_stimuli(starting_stimuli);

while true
    probed_idxs = next_directions_to_test(trialspecs, N_INTERLEAVED_TS, minTF, maxTF);
    if isempty(probed_idxs), break; end
    
    % simulate a measurement along one or more directions
    [measured_thresholds,oog] = simulate_noisy_observer(cat(1,trialspecs(probed_idxs).colordir), [trialspecs(probed_idxs).tf], pplum, pprg);
    [trialspecs(oog).oog] = deal(true);
    [trialspecs(probed_idxs).measuredthreshold] = measured_thresholds{:};
    
    % plot the result of the measurement (magenta within tolerance; blue is scaled to represent the
    % error between prediction and measurement)
    if DEBUG, trialspecs = plot_progress(trialspecs, probed_idxs); end
    
    % subdivide triangles to make new directions to probe
    trialspecs = setup_next_rounds(trialspecs, probed_idxs);
    
    if length(trialspecs) == N_INITIAL_TS && sum([trialspecs.done]) == N_INITIAL_TS
        trialspecs = [trialspecs setup_second_round(trialspecs)]; %#ok<AGROW>
    end
    
    % black points show queued TS
    if DEBUG, trialspecs = plot_predictions_queue(trialspecs, PLOT_TRIANGLES, PLOT_NORMALS); end
    drawnow();
end

% return the N trialspecs indices of the next directions to test
function idxs = next_directions_to_test(trialspecs, N, minTF, maxTF)
global done_with_initial_ts

idxs = zeros(1, N);
not_done = ~[trialspecs.done];
if ~done_with_initial_ts
    initial_stimuli = all(isnan(cat(1, trialspecs.parentvertices)), 2)';
    initial_round_stim = find(initial_stimuli & not_done, N);
    idxs(1:length(initial_round_stim)) = initial_round_stim;
    if isempty(initial_round_stim), done_with_initial_ts = true; end
end

if all(idxs > 0), return; end

eligible_ts = not_done;
eligible_ts(idxs(idxs > 0)) = false;
% set your desired ordering of the trialspec indices in `your_order`
% your_order = 1:length(trialspecs);
[~,your_order] = sort([trialspecs.parent_error], 'descend');
% your_order = sort_by_curvature(trialspecs, not_done, minTF, maxTF);
eligible_idxs = your_order(eligible_ts(your_order));
fprintf('#%02d %g\n', [eligible_idxs; [trialspecs(eligible_idxs).parent_error]]);

if isempty(eligible_idxs), return; end

start_here = find(idxs == 0, 1);
idxs(start_here) = eligible_idxs(1);
% `idxs` must be trialspec indices that all lie on mutually-exclusive triangles. We flag `false` in
% `eligible_ts` for trialspec indices that lie on or within a triangle also shared by the TS indexed
% by `seed_ts_idx`. We repeat this process until we have up to `N` trialspec indices that lie on
% mutually-exclusive triangles.
for k = start_here+1:length(idxs)
    seed_ts_idx = idxs(k-1);
    eligible_ts(seed_ts_idx) = false;
    ts_parents = trialspecs(seed_ts_idx).parentvertices;
    % Construct a set of parents whose children are ineligible for testing
    if any(isinf(ts_parents))
        ts_parents(isinf(ts_parents)) = [];
        [other_parents, descendents] = find_relatives(trialspecs, ts_parents, seed_ts_idx);
        [other_parents, descendents] = process_relatives(trialspecs, other_parents, descendents);
        eligible_ts(descendents) = false;
        vert_design = fullfact([length(other_parents) 2]);
        parents_of_ineligible_ts = sort([other_parents(vert_design(:,1)) ts_parents(vert_design(:,2))' inf(size(vert_design, 1), 1)], 2);
    else
        parents_of_ineligible_ts = sort([ts_parents([1 2; 1 3; 2 3]) inf(3, 1)], 2);
    end
    ineligible_ts = find_ts_locations(cat(1, trialspecs.parentvertices), parents_of_ineligible_ts);
    eligible_ts(ineligible_ts) = false;
    next_idxs = your_order(eligible_ts(your_order));
    if isempty(next_idxs), break; end
    idxs(k) = next_idxs(1);
end
idxs = idxs(idxs > 0);

function trialspecs = setup_next_rounds(trialspecs, TS_idxs)
global DEBUG
for measured_vert = TS_idxs
    trialspecs(measured_vert).done = true;
    this_vert_error = (log(trialspecs(measured_vert).measuredthreshold) ...
        - log(trialspecs(measured_vert).predictedthreshold))^2;
    if ~isnan(trialspecs(measured_vert).predictedthreshold)
        parent_verts = trialspecs(measured_vert).parentvertices;
        if any(isinf(parent_verts)) % `measured_vert` lies on an edge
            parent_verts(isinf(parent_verts)) = [];
            % Find vertices connected to `parent_verts` and the descendents (children) of those
            % triangles; i.e., the triangle derived from TS with indices `[other_parents(i)
            % parent_verts]` and the midpoint trialspec index `descendents(i)`.
            [other_parents, descendents] = find_relatives(trialspecs, parent_verts, measured_vert);
            [other_parents, descendents] = process_relatives(trialspecs, other_parents, descendents);
            n = length(other_parents);
            vert_design = fullfact([n 2]);
            
            vert_list = [ ...
                % bisect between our immediate parents
                parent_verts(1) measured_vert inf
                parent_verts(2) measured_vert inf
                % subdivide adjacent triangles
                other_parents(vert_design(:,1)) parent_verts(vert_design(:,2))' measured_vert+zeros(size(vert_design, 1), 1)
                ];
            
            descendents_error = [trialspecs(descendents).parent_error]';
            parental_error = [this_vert_error([1; 1]); descendents_error(vert_design(:,1))];
            
            % We replace TS indexed by `descendents` with new TS defined as midpoints along the
            % edges between TS vertices indexed by `measured_vert` and `other_parents`.
            new_dirs = make_new_directions(trialspecs, ...
                [measured_vert+zeros(n, 1) other_parents inf(n, 1)], ...
                descendents_error); % use existing parent error for new directions
            if DEBUG
                delete([[trialspecs(descendents).htri] ...
                    [trialspecs(descendents).hpoint] ...
                    [trialspecs(descendents).hquiv]]);
            end
            trialspecs(descendents) = new_dirs;
        else % face
            vert_list = [ ...
                parent_verts(1) measured_vert inf
                parent_verts(2) measured_vert inf
                parent_verts(3) measured_vert inf
                parent_verts([1 2]) measured_vert
                parent_verts([1 3]) measured_vert
                parent_verts([2 3]) measured_vert
                ];
            parental_error = this_vert_error/6; % propagate this vertex's error to all children
        end
        trialspecs = [trialspecs make_new_directions(trialspecs, vert_list, parental_error)]; %#ok<AGROW>
    end
end

function new_ts = make_new_directions(trialspecs, parent_verts, parental_error)
global DEBUG
if length(parental_error) < size(parent_verts,1)
    parental_error = parental_error(ones(size(parent_verts, 1), 1));
end
parent_verts = sort(parent_verts, 2);
% don't make duplicate TS
Ldupe = find_ts_locations(parent_verts, cat(1, trialspecs.parentvertices));
parent_verts(Ldupe,:) = [];
parental_error(Ldupe) = [];
% ignore requested directions if all parents are oog
Lall_oog = all_parents_oog(trialspecs, parent_verts);
parent_verts(Lall_oog,:) = [];
parental_error(Lall_oog) = [];
ndirs = size(parent_verts, 1);
new_ts(1:ndirs) = trialspec_template();
for ii = 1:ndirs
    valid_idxs = parent_verts(ii,:);
    valid_idxs(isinf(valid_idxs)) = [];
    signs = sign(valid_idxs)';
    valid_idxs = abs(valid_idxs);
    
    unit_cdirs = cat(1, trialspecs(valid_idxs).colordir); % each row is a stimulus
    tfs = [trialspecs(valid_idxs).tf]';
    threshs = bsxfun(@times, unit_cdirs, signs.*cat(1, trialspecs(valid_idxs).measuredthreshold));
    pred_threshs = mean(threshs);
    
    new_ts(ii).parent_error = parental_error(ii);
    new_ts(ii).colordir = pred_threshs./norm(pred_threshs); % L, M
    new_ts(ii).tf = mean(tfs); % Hz
    new_ts(ii).predictedthreshold = norm(pred_threshs);
    new_ts(ii).parentvertices = parent_verts(ii,:);
    
    % computing normals
    vdiff = diff([threshs tfs]);
    if size(unit_cdirs, 1) < 3 % i.e., less than 3 parents (edge point)
        new_ts(ii).normal = [-vdiff(2) vdiff(1) vdiff(3)];
    else
        new_ts(ii).normal = cross(vdiff(1,:),vdiff(2,:));
    end
    new_ts(ii).normal = sign([new_ts(ii).colordir new_ts(ii).tf]*new_ts(ii).normal')*new_ts(ii).normal;
    new_ts(ii).normal = new_ts(ii).normal/norm(new_ts(ii).normal);
    
    if DEBUG, fprintf('made new direction from vertices %d %d %d\n', parent_verts(ii,:)); end
end

function tf = all_parents_oog(trialspecs, parent_verts)
tf = false(size(parent_verts,1),1);
for ii = 1:length(tf)
    tf(ii) = all([trialspecs(abs(parent_verts(ii,isfinite(parent_verts(ii,:))))).oog]);
end

function ts_loc = find_ts_locations(A, B) % identical behavior to C = ismember(A,B,'rows');
ts_loc = any(all(bsxfun(@eq, reshape(B', 1, 3, []), A), 2), 3);

% Find "descendents" with .parentvertices containing shared_verts or -shared_verts
function [third_parents,descendents] = find_relatives(trialspecs, shared_verts, skip)
% ismembc's second input must be sorted!!
if isempty(shared_verts), return; end
% ismembc puts 1s where elements of shared_verts are in all_parents. We have to also consider
% "matches" with negative vertex indices (a convention applied in the second round). This means,
% for example, that shared_verts = [1 3] will match to the following set of parents: [1 3 4],
% [-3 -1 7], and [1 3 Inf].
flipped_verts = -shared_verts(end:-1:1); % this could make -Infs (forced to +Infs on next line)
flipped_verts = [flipped_verts(isfinite(flipped_verts)) inf(1, sum(isinf(flipped_verts)))];
all_parents = cat(1, trialspecs.parentvertices);
unshared_verts_pos = ~ismembc(all_parents, shared_verts);
unshared_verts_pos(skip,:) = false;
unshared_verts_neg = ~ismembc(all_parents, flipped_verts);
ts_matched_parents = find(sum(unshared_verts_pos, 2) == 1);
ts_flipped_parents = find(sum(unshared_verts_neg, 2) == 1);
descendents = [ts_matched_parents; ts_flipped_parents];
% We want the logical indexing below to go row-wise so transpose the parents and the masks (this is
% crucial to match the vertex order in descendents and third_parents).
parents_matched = all_parents(ts_matched_parents,:)';
parents_flipped = -all_parents(ts_flipped_parents,:)';
third_parents = [parents_matched(unshared_verts_pos(ts_matched_parents,:)'); ...
    parents_flipped(unshared_verts_neg(ts_flipped_parents,:)')];

function [other_parents,descendents] = process_relatives(trialspecs, other_parents, descendents)
% We need to make new TS within the 1 or 2 adjacent triangles that share the edge `measured_vert`
% lies on. Only in these adjacent regions is it safe to make TS without crossing edges from previous
% rounds. We use the TS indices we found from `find_relatives(...)` to find the applicable adjacent
% regions. If `descendents(i) == other_parents(j)`, then this vertex makes a triangle with
% `parent_verts` _inside_ a triangle with vertices `other_parents(i)` and `parent_verts`. Thus, we
% don't want to subdivide triangles containing vertices `measured_vert` and `other_parents(i)`
% because we'll overlap the enclosed triangle. Therefore, we remove these i-th entries, called here
% `incestuous_verts`. Unfortunately, this heuristic isn't enough to ensure non-overlapping regions;
% we must also test for doneness.

% Get the elements of `descendents` that appear in `other_parents`. We have to abs `other_parents`
% because of a sign convention we use in the initial round to disambiguate vertex indices.
incestuous_verts = ismembc(descendents, sort(abs(other_parents)));
% We don't want to replace measured TS so we also check for doneness.
remove_verts = incestuous_verts | [trialspecs(descendents).done]';
descendents(remove_verts) = [];
other_parents(remove_verts) = [];

function ordered_idxs = sort_by_curvature(ts, not_done, minTF, maxTF)
% Some normals could be on the opposite sides of the surface; this happens because of the
% negative vertex index convention established on round 2. If the siblings span the "seam"
% between positive and negative indices, then we could erroneously calculate a large angle.
% We detect siblings spanning the seam by testing for any sign flips in the siblings' parent
% indices compared to the parents of the current vertex.
global DEBUG
labels = cell(size(not_done));
for curr_idx = find(not_done)
    parents = ts(curr_idx).parentvertices;
    is_face = all(isfinite(parents));
    if ~is_face
        parents = parents(isfinite(parents));
        parents_tfs = [ts(abs(parents)).tf];
        % exterior edge:
        if all(softEq(parents_tfs, minTF)) || all(softEq(parents_tfs, maxTF))
            [~,sibs1] = find_relatives(ts, [parents(1) inf], curr_idx);
            [~,sibs2] = find_relatives(ts, [parents(2) inf], curr_idx);
            siblings = [sibs1; sibs2];
            L = softEq([ts(siblings).tf], ts(curr_idx).tf) & ~[ts(siblings).done];
            siblings = siblings(L);
            labels{curr_idx} = 'ext';
        else % interior edge:
            [other_parents,siblings] = find_relatives(ts, parents, curr_idx);
            if length(siblings) > 2
                [~,siblings] = process_relatives(ts, other_parents, siblings);
            end
            labels{curr_idx} = 'int';
        end
        assert(length(siblings)==2, 'Expecting only two siblings for edge points');
        normals = cat(1, ts(siblings).normal);
    else % face:
        other_parents = inf(50, 1);
        siblings = inf(size(other_parents));
        end_ptr = 0;
        for perm = [1 2; 1 3; 2 3]'
            [otherstmp, sibstmp] = find_relatives(ts, parents(perm), curr_idx);
            [otherstmp, sibstmp] = process_relatives(ts, otherstmp, sibstmp);
            copy_range = end_ptr+1:end_ptr+length(otherstmp);
            other_parents(copy_range) = otherstmp;
            siblings(copy_range) = sibstmp;
            end_ptr = copy_range(end);
        end
        siblings = siblings(isfinite(other_parents));
        assert(2<=length(siblings)&&length(siblings)<=3, 'Expecting 2-3 siblings for face points');
        normals = cat(1, ts(curr_idx).normal, ts(siblings).normal); % 1st row is this point's normal
        labels{curr_idx} = 'fce';
    end
    % Check for seam crossing (i.e., when a subset of parents are -parents in siblings)
    flip_normal = sum(ismembc(-cat(1, ts(siblings).parentvertices), parents), 2) > 0;
    % The first row of `normals` is the normal of the current vertex ONLY if it's a face point. In
    % this case we need to prepend a false to `flip_normal` so we don't flip the sign of the current
    % vertex's normal.
    normals([false(is_face); flip_normal],:) = -normals([false(is_face); flip_normal],:);
    thetas = acos(normals(1,:)*normals(2:end,:)');
    if siblings(1) ~= siblings(2) % this needs to be more robust
        ts(curr_idx).curvature = max(thetas);
    else % degenerate case of parallel normal vectors (the exterior edges after round 1)
        ts(curr_idx).curvature = pi/2;
    end
end
[ordered_vals,ordered_idxs] = sort([ts.curvature], 'descend');
ordered_idxs = ordered_idxs(~isnan(ordered_vals));
labels = labels(ordered_idxs);
if DEBUG
    print_me = [labels; num2cell(ordered_idxs(:)'); num2cell([ts(ordered_idxs).curvature])];
    fprintf('%s #%02d -> %g\n', print_me{:});
end

function trialspecs = setup_starting_stimuli(stimuli)
trialspecs(1:size(stimuli, 1)) = trialspec_template();
for curr_dir = 1:length(trialspecs)
    trialspecs(curr_dir).colordir = stimuli(curr_dir,[1 2])./norm(stimuli(curr_dir,[1 2])); % L, M
    trialspecs(curr_dir).tf = stimuli(curr_dir,3); % Hz
end

function new_specs = setup_second_round(trialspecs)
% 3        4       -3       -4
% +---//---+---//---+---//---+ max temporal freq
% |        |        |        |
% |        |        |        |
% +--------+--------+--------+ min temporal freq
% 1        2       -1       -2
% Cutting each face on the same diagonal seems to be important so that you don't end up with points
% that are very close to each other (and areas that are unsampled).
perms = [1 2 4; 1 3 4; -3 -1 2; -3 2 4];
perms = [perms; 1 2 inf; -1 2 inf; 3 4 inf; -3 4 inf];
perms = [perms; 1 3 inf; 2 4 inf; 1 4 inf; -3 2 inf];
new_specs = make_new_directions(trialspecs, perms, nan);

function [measurements,truncated] = simulate_noisy_observer(LM_dirs, tfs, pplum, pprg)
noise_lim = 0;
lumthresh = fnval(pplum,tfs);
rgthresh = fnval(pprg,tfs);
% rotating colordir from (L, M) to (LUM, RG).
LUMRG_dirs = LM_dirs/sqrt(2)*[1 -1; 1 1];
surface_points = dot(LUMRG_dirs.^2, [1./lumthresh; 1./rgthresh]'.^2, 2)'.^-.5;
% The line above was a pain to derive. The point is that we want to
% find a scale factor, 'c', such that
% (c*cos(theta)).^2/a.^2 + (c*sin(theta)).^2/b.^2 = 1;
% where cdir = [cos(theta) sin(theta)];
mult_noise = exp(unifrnd(-noise_lim, noise_lim, size(surface_points)));
noisy_surface_points = surface_points.*mult_noise;
% truncate at the edge of the gamut
[measurements,truncated] = truncate_at_gamut_edge(LM_dirs, noisy_surface_points);
measurements = num2cell(measurements);

function [points,oog] = truncate_at_gamut_edge(dirs, points)
% M and bkgndrgb are from Dell4's calibration data
M = [...
    0.050668871815624   0.100663099362598   0.017291476480446
    0.018212331805184   0.105784629645293   0.026412529028706
    0.001638607887914   0.007307556884068   0.078692093480351];
bkgndrgb = [...
    0.563717835642795   0.468688047715644   0.496065369042598];

ccs = [dirs.*points(ones(size(dirs,2),1),:)' zeros(size(dirs,1),1)];
[in_gamut,gamut_scalars] = gamutCheck(ccs', bkgndrgb, M, 'both');
oog = ~in_gamut;
points(oog) = points(oog).*gamut_scalars(oog);

function yn = within_tolerance(thresh_ratio)
tol = 0.2;
yn = abs(log(thresh_ratio)) < abs(log(1+tol));

function trialspec = trialspec_template()
global DEBUG
trialspec = struct('colordir', [], 'tf', [], 'measuredthreshold', [], 'parent_error', nan, ...,
    'curvature', nan, 'predictedthreshold', nan, 'parentvertices', nan(1, 3), 'oog', false, ...
    'normal', nan(1,3), 'done', false);
if DEBUG, trialspec.hpoint = []; trialspec.htri = []; trialspec.hquiv = []; end

function trialspecs = plot_progress(trialspecs, TS_idxs)
for curr_dir = TS_idxs
    thresh_ratio = abs(log(trialspecs(curr_dir).measuredthreshold)-log(trialspecs(curr_dir).predictedthreshold));
    if isnan(thresh_ratio) || within_tolerance(thresh_ratio)
        plot_opts = {'marker' 'o' 'color' 'm' 'markerfacecolor' 'm'};
    else
        plot_opts = {'marker' 'o' 'color' 'b' 'markerfacecolor' 'b' 'markersize' 20*thresh_ratio+eps};
    end
    xy = trialspecs(curr_dir).measuredthreshold*trialspecs(curr_dir).colordir;
    if isempty(trialspecs(curr_dir).hpoint)
        trialspecs(curr_dir).hpoint(1) = plot3(xy(1),xy(2),trialspecs(curr_dir).tf,plot_opts{:});
        trialspecs(curr_dir).hpoint(2) = plot3(-xy(1),-xy(2),trialspecs(curr_dir).tf,plot_opts{:});
        set([trialspecs(curr_dir).hpoint], 'ButtonDownFcn', @(s,e) fprintf('\t-> (#%d) %g %g %g\n', ...
            curr_dir, xy, trialspecs(curr_dir).tf));
    else
        set(trialspecs(curr_dir).hpoint(1),'xdata',xy(1),'ydata',xy(2),'zdata',trialspecs(curr_dir).tf,plot_opts{:});
        set(trialspecs(curr_dir).hpoint(2),'xdata',-xy(1),'ydata',-xy(2),'zdata',trialspecs(curr_dir).tf,plot_opts{:});
    end
end

function trialspecs = plot_predictions_queue(trialspecs, plot_tris, plot_quivs)
queue_idxs = find(~[trialspecs.done] & cellfun('isempty', {trialspecs.measuredthreshold}));
if plot_tris
    set([trialspecs.htri], 'visible', 'off');
    set([trialspecs.hquiv], 'visible', 'off');
end
for curr_dir = queue_idxs
    if isnan(trialspecs(curr_dir).predictedthreshold), continue; end
    xy = trialspecs(curr_dir).predictedthreshold*trialspecs(curr_dir).colordir;
%     if isempty(trialspecs(curr_dir).hpoint)
%         trialspecs(curr_dir).hpoint(1) = plot3(xy(1),xy(2),trialspecs(curr_dir).tf,'k.');
%         trialspecs(curr_dir).hpoint(2) = plot3(-xy(1),-xy(2),trialspecs(curr_dir).tf,'k.');
%         set([trialspecs(curr_dir).hpoint], 'ButtonDownFcn', @(s,e) fprintf('\t-> (#%d) %g %g %g\n', ...
%             curr_dir, xy, trialspecs(curr_dir).tf));
%     end
    if plot_quivs
        if isempty(trialspecs(curr_dir).hquiv)
            normal = trialspecs(curr_dir).normal;
            trialspecs(curr_dir).hquiv = quiver3([-xy(1) xy(1)], [-xy(2) xy(2)], ...
                trialspecs(curr_dir).tf([1 1]), [-normal(1) normal(1)], [-normal(2) normal(2)], ...
                [-normal(3) normal(3)], 'hittest', 'off');
        else
            set([trialspecs(curr_dir).hquiv], 'visible', 'on');
        end
    end
    if plot_tris
        parent_idxs = [trialspecs(curr_dir).parentvertices];
        if isempty(parent_idxs) || any(isinf(parent_idxs)), continue; end
        if ~isempty(trialspecs(curr_dir).htri), set([trialspecs(curr_dir).htri], 'visible', 'on'); continue; end
        parent_signs = sign(parent_idxs)'; parent_idxs = abs(parent_idxs);
        parent_cdirs = cat(1, trialspecs(parent_idxs).colordir);
        parent_tfs = [trialspecs(parent_idxs).tf]';
        parent_xys = bsxfun(@times, parent_cdirs, parent_signs.*[trialspecs(parent_idxs).measuredthreshold]');
        trialspecs(curr_dir).htri(1) = patch(parent_xys(:,1),parent_xys(:,2),parent_tfs,'g','facealpha',.25,'edgealpha',.25,'hittest','off');
        trialspecs(curr_dir).htri(2) = patch(-parent_xys(:,1),-parent_xys(:,2),parent_tfs,'g','facealpha',.25,'edgealpha',.25,'hittest','off');
    end
end
