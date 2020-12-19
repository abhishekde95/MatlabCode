function [] = WN_plotSTAest(imagevectors, nstixperside, figtitle, save_figs, p2s)
%%
% This function plots the STA over time given the input of STAs for
% different offsets. 
%
% Inputs:
%   'imagevectors' [ncovariates, noffset] - 2D matrix where each column is 
%       vector of weights for a given offset.
%   'nstixperside' [int] - Number of stix per side of stimulus
%   'figtitle' [string] - String to title figure and name saved file.
%   'save_figs' [bool] - Save figure or don't
%   'p2s' [string] - path to save data
%
% Outputs:
%   None returned. Displays and can save figure. 
%%
rowidxs = reshape([1:3*nstixperside^2],[nstixperside^2 3]);
maxes = []; mins = []; means = [];
for i = 1:3
    maxes = [maxes; max(max(imagevectors(rowidxs(:,i),:)))];
    mins = [mins; min(min(imagevectors(rowidxs(:,i),:)))];
    means = [means; mean(mean(imagevectors(rowidxs(:,i),:)))];
end
potentialnormfactors = [([.5; .5; .5]-eps)./(maxes-means); ...
    (-[.5; .5; .5]+eps)./(mins-means)];    
% 'eps' in above line is a kludge that is required for avoiding
% out of bounds errors.

potentialnormfactors(potentialnormfactors < 0) = []; % if min > mu or max < mu
normfactor = min(potentialnormfactors);
muvect = reshape(repmat([.5 .5 .5],nstixperside^2,1),nstixperside^2*3,1);

% Plot each vector as an image
figure('Renderer', 'painters', 'Position', [100 100 275*size(imagevectors, 2) 200])
for jj = 1:size(imagevectors,2)
    subplot(1,size(imagevectors,2),jj)
    STAlars = normfactor*(imagevectors(:,jj)-mean(imagevectors(:,jj)))+muvect;
    STAlars = reshape(STAlars,[nstixperside nstixperside 3]);
    image(STAlars);
    title([figtitle ' - Offset: ' num2str(jj-1)],'FontSize', 12);
    axis square
end
if save_figs; saveas(gcf, [p2s figtitle '.png']); end