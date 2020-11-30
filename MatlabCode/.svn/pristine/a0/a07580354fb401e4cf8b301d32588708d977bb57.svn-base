function inRF = getRFfromSTA(STA, mode, thresh, gauss_var)
% inRF = getRFfromSTA(STA, mode, thresh, [sigma2])
%
% Support function for pulling out an spatial RF from a spatiotemporal
% color STA.
%
% STA shoud be n x maxT (the standard return format from
% getWhtnsStats) where n is 3*nstixperside^2
% mode 0 = highest energy stixel only
% mode 1 = consider only stixels with energy in the top 'thresh' percent
% relative to the noise.
% mode 2 = consider only stixels with energy above 'thresh' relative to
% a chi-squared CDF based on the noise.
%
% 'thresh' should be between 0 and 1
%
% The estimate of noise comes from the first frame.

DEBUGGING = 0;

if nargin == 1
    mode = 1;
    thresh = 1;
end
if ~ismatrix(STA)
    error('getRFfromSTA expects n x m input');
end

if mode ~= 0 && (thresh < 0 || thresh > 1)
    error('threshold should be between 0 and 1');
end
maxT = size(STA,2);
nstixperside = sqrt(size(STA,1)/3);
if floor(nstixperside) ~= nstixperside
    error('number of row must be 3*nstixperside^2');
end
STA = reshape(STA,[nstixperside.^2  3 maxT]);
if nargin < 4 % estimating variance from STA itself
    gauss_var = mean(var(STA(:,:,1)));
end
noise = sum((STA(:,:,1)./sqrt(gauss_var)).^2,2);

% Debugging
if DEBUGGING
    figure('Position',[6 385 560 420]); subplot(2,2,1); hold on;
    [n,x] = hist(noise);
    bar(x,n./sum(n));
    y = chi2pdf(x,3);
    plot(x,y./sum(y));
end

energy = sum(sum((STA./sqrt(gauss_var)).^2,2),3); % summing across guns and time
whichpix_center = energy == max(energy);

if mode == 0
    return
elseif mode == 1
    whichpix = whichpix_center | energy > maxT*prctile(noise,thresh*100);
elseif mode == 2
    whichpix = whichpix_center | energy > maxT*(chi2inv(thresh,3));
else
    error('unknown mode');
end

% Getting rid of points that are outside the convex hull of the main cluster
if sum(whichpix) > 1
    [tmpi,tmpj] = ind2sub([nstixperside nstixperside],find(whichpix));
    ij = [tmpi,tmpj];
    [center_i, center_j] = ind2sub([nstixperside nstixperside],find(whichpix_center));
    center_ij = [center_i, center_j];
    T = clusterdata(ij,'linkage','single','distance','euclidean','criterion','distance','cutoff',sqrt(2));
    
    %clusternmembers = [];
    %for k =1:max(T)
    %    clusternmembers(k) = sum(T == k);
    %end
    %dominantcluster = find(clusternmembers == max(clusternmembers));
    dominantcluster = T(all(ij == repmat(center_ij,size(ij,1),1),2)); % dominant cluster must contain the hottest stixel
    
    clustermat = zeros(nstixperside, nstixperside);
    clustermat(sub2ind(size(clustermat),ij(T==dominantcluster,1),ij(T==dominantcluster,2))) = 1;
    
    if sum(clustermat(:)) > 2 & length(unique(ij(T==dominantcluster,1))) > 1 & length(unique(ij(T==dominantcluster,2))) > 1
        % Then get convex hull
        dominantclusteridxs = ij(T==dominantcluster,:);
        K=convhull(dominantclusteridxs(:,1), dominantclusteridxs(:,2));
        tmp = dominantclusteridxs(K,:);
        [x,y] = meshgrid(1:nstixperside,1:nstixperside);
        inRF = reshape(inpolygon(x(:),y(:),dominantclusteridxs(K,2),dominantclusteridxs(K,1)),[nstixperside, nstixperside]);
    else
        inRF = clustermat;
    end
else
    inRF = reshape(whichpix,nstixperside,nstixperside);
end
if DEBUGGING
    subplot(2,2,3);
    imagesc(inRF);
    axis square;
    subplot(2,2,4);
    imagesc(reshape(energy,nstixperside,nstixperside));
    axis square;
    colormap(gray)
end
end