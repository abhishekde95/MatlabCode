
% How much variance in STA RGBs do we capture with a single line through
% RGB space versus two vectors (not necesarily 180° apart) that meet at 
% the origin? 

close all; clearvars;
filenames = {'K040708003.nex','M121417002.nex','K110408003.nex','M082217003.nex','K052308001.nex','S071510004.nex'};

Separability_index = []; % Adding a separability index

for fileidx = 1:length(filenames)
    if filenames{fileidx}(1) == 'K'
        minipath = ['Greg',filesep,'Kali',filesep,'2008'];
    elseif filenames{fileidx}(1) == 'S'
        minipath = ['Greg',filesep,'Sedna',filesep,'2010'];
    elseif filenames{fileidx}(1) == 'M'
        minipath = ['Abhishek',filesep,'Maui'];
    else
        minipath = [];
    end
    WN=nex2stro([nexfilepath,filesep,minipath,filesep,filenames{fileidx}]);
    
    maxT = 9;
    spikename = 'sig001a'; % A hack. Might not be correct for all cells.
    framerate = WN.sum.exptParams.framerate;
    nstixperside = WN.sum.exptParams.nstixperside;
    ntrials = length(WN.sum.absTrialNum);
    stimonidx = find(strcmp(WN.sum.trialFields(1,:),'stim_on'));
    stimoffidx = find(strcmp(WN.sum.trialFields(1,:),'all_off'));
    nframesidx = find(strcmp(WN.sum.trialFields(1,:),'num_frames'));
    noisetypeidx = find(strcmp(WN.sum.trialFields(1,:),'noise_type'));
    sigmaidxs = strmatch('sigma',WN.sum.trialFields(1,:));
    
    gammaTable = WN.sum.exptParams.gamma_table;
    gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
    
    % Reconstructing the M matrix
    fundamentals = WN.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = WN.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    
    % Getting the background rgb/lms
    ridx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_r'));
    gidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_g'));
    bidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_b'));
    bkgndRGB = [mode(WN.trial(:,ridx)), mode(WN.trial(:,gidx)), mode(WN.trial(:,bidx))];
    bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
    bkgndlms = M*bkgndrgb;
    L = WN.trial(:,noisetypeidx) == 1; % gun noise only
    WN.ras(~L,:) = [];
    WN.trial(~L,:) = [];
    out = getWhtnsStats(WN,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);
    STAs = out{1}./out{3};
    sqrtenergy = sqrt(sum(STAs.^2));
    peakframe = find(sqrtenergy == max(sqrtenergy),1,'first');
    whichframes = peakframe-1:peakframe+1;
    weightingcoefs = sqrtenergy(whichframes);
    weightingcoefs = weightingcoefs./norm(weightingcoefs);
    %weightingcoefs = [0 1 0]; % Just pulling out peak frame
    
    weightedSTA = STAs(:,whichframes)*weightingcoefs';
    reshapedweightedSTA = reshape(weightedSTA,[nstixperside.^2 3]);
    
    figure; axes; hold on;
    title(filenames{fileidx});
    plot3(reshapedweightedSTA(:,1),reshapedweightedSTA(:,2),reshapedweightedSTA(:,3),'ko');
    plot3(0,0,0,'m*');
    coef = pca(reshapedweightedSTA); 
    plotscalefactor = max(sqrt(sum(reshapedweightedSTA.^2,2))); % give PC1 a reasonable length for plotting
    plot3([-1 1]*coef(1,1)*plotscalefactor,[-1 1]*coef(2,1)*plotscalefactor,[-1 1]*coef(3,1)*plotscalefactor,'k-');
    
    v1 = coef(:,1);
    v2 = -coef(:,1);
    
    myvariance(reshapedweightedSTA,[v1 v2]); % The mean squared projection onto the PC1
    f = @(x)-1*myvariance(reshapedweightedSTA, x); % Setting up anonymous function for fmincon
    [vs,fval] = fmincon(f,[v1 v2],[],[],[],[],[],[],@unitnormcon);
     
    plot3([0 1]*vs(1,1)*plotscalefactor,[0 1]*vs(2,1)*plotscalefactor,[-0 1]*vs(3,1)*plotscalefactor,'r-');
    plot3([0 1]*vs(1,2)*plotscalefactor,[0 1]*vs(2,2)*plotscalefactor,[-0 1]*vs(3,2)*plotscalefactor,'r-');
    axis equal;
    
    vec1_pca = coef(:,1)/norm(coef(:,1));
    vec1_proj = vs(:,1)/norm(vs(:,1));
    vec2_proj = vs(:,2)/norm(vs(:,2));
    
    tmp_SI = abs(dot(vec1_proj,vec2_proj))*max([abs(dot(vec1_pca,vec1_proj)) abs(dot(vec1_pca,vec2_proj))]);
    Separability_index = [Separability_index; tmp_SI];
    
    drawnow;
%     keyboard
    disp(['Original variance: ',num2str(myvariance(reshapedweightedSTA,[v1 v2]))])
    disp(['New variance: ',num2str(-fval)])
    disp(['Ratio of variances: ',num2str(-fval/myvariance(reshapedweightedSTA,[v1 v2]))]);
   
    angles = min(acos(vs'*[v1 v2])*180/pi);
    disp(['angles between PC1 and new vectors (°): ',num2str(angles)]);
    disp(' '); % Skipping a line
end

% Note: this function does not subtract the mean from the data, which makes
% sense in this context, so the "variance" it returns is slightly different 
% from the variance of projections onto the PC1, even when "basis" = [PC1
% -PC1].
function v = myvariance(data, basis)
projs = data*basis;
projs(projs<0) = 0;
projs = max(projs,[],2);
v = mean(projs.^2);
end

function [c,ceq] = unitnormcon(x)
c = [];
ceq = sum(x.^2)-1;
end