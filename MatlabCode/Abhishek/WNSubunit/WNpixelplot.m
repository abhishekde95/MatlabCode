function [plot_counter,STAs] = WNpixelplot(stro,plot_counter,mask_changes,scriptname)
% Abhishek's copy of WNSubunitPLot.m

global maskidx spikeidx nstixperside seedidx nframesidx stimonidx muidxs sigmaidxs
global msperframe maxT xx yy M 

%disp(mask_changes);
% STCOVmex - This function calculates the covariance matrices, efficiently
% and quickly. Function takes 2 arguments. The string 'init' as the first
% argument and the number of stimulus dimensions as the second argument in
% the format [npics, ncols, nframes].
% In subsequent cells, send as the first argument a vector of RGBs (the
% total number of elements should be xpix*ncols*nframes, and the order
% shoud be first all the elements for frame 1, then all the elements for
% frame 2, and so on). The second argument should be an nframesx1 matrix of
% spike counts.
% Send the 'return' string as a solo argument to get the output - a
% structure containing {STS, STCross, and n}
if isfield(stro.sum.exptParams,'nrepframes')
    if ~isnan(stro.sum.exptParams.nrepframes)
        nvblsperstimupdate = stro.sum.exptParams.nrepframes;
    else
        nvblsperstimupdate = 1;
    end
else
    nvblsperstimupdate = 1;
end
for mask_span = mask_changes
    STCOVmex('init', {nstixperside^2 3 maxT});
    for k = mask_span(1):mask_span(2)
        nframes = stro.trial(k,nframesidx);
        if nframes == 0, continue; end
        
        seed = stro.trial(k,seedidx);
        mu = stro.trial(k,muidxs)/1000;
        sigma = stro.trial(k,sigmaidxs)/1000;
        
        org_mask = stro.ras{k,maskidx}; % useful for subunits computation
        nrandnums_perchannel = nstixperside^2; % In this case, it is the 100 pixels which are flickering when no subunits are selected
        % assuming Gaussian gun noise only, random number generator
        % routine as a mexfile (getEJrandnums.mexw64)
        invnormcdf = bsxfun(@plus, bsxfun(@times, yy, sigma), mu);
        %         figure(3),plot(invnormcdf);
        randnums = getEJrandnums(3*nrandnums_perchannel*nframes, seed); % random numbers for 300 x 9 pixels
        randnums = reshape(randnums, 3*nrandnums_perchannel, nframes);
        for gun = 1:3
            idxs = (1:nrandnums_perchannel)+nrandnums_perchannel*(gun-1);
            randnums(idxs,:) = reshape(invnormcdf(randnums(idxs,:)+1,gun),[length(idxs) nframes]);
        end
        
        rgbs = randnums;
        t_stimon = stro.trial(k, stimonidx);
        spiketimes = (stro.ras{k,spikeidx}-t_stimon)*1000; % observing spiketimes in milliseconds
%         frametimes = linspace(0, nframes*msperframe, nframes)+(msperframe/2)';
        frametimes = linspace(0, nframes*msperframe*nvblsperstimupdate, nframes)+(msperframe*nvblsperstimupdate/2)';
        spiketimes(spiketimes < maxT*msperframe) = [];
        spiketimes(spiketimes > frametimes(end)) = [];
        n = hist(spiketimes, frametimes);
        STCOVmex(rgbs(:),n);
    end
    out = STCOVmex('return');
    STS = out{1};
    STCross = out{2};
    nspikes = out{3};
    clear STCOVmex;
    clear out;
    STAs = STS/nspikes;
    STCs = zeros(size(STCross));
    for i = 1:maxT
        tmp = STS(:,i)*STS(:,i)';
        STCs(:,i) = (nspikes.*STCross(:,i)-tmp(:))/(nspikes*(nspikes-1));
    end
    an_mask = [];
    an_mask = zeros(nstixperside, nstixperside);   % Pixel mask.  Someday make a tool to make this non-zero.
    Lmask = logical(repmat(~an_mask(:),[3 1]));
    PCs = [];
    for i = 1:size(STCs,2)
        STC = reshape(STCs(:,i), 3*nstixperside^2, 3*nstixperside^2);
        subSTC = STC(Lmask, Lmask);
        subSTA = STAs(Lmask,i);
        % subtracting the samples from the STA to ensure the PCs are
        % orthogonal to the the STA.
        P = eye(size(subSTC)) - subSTA*inv(subSTA'*subSTA)*subSTA';
        subSTC = P*subSTC*P';
        [tmp,d] = eig(subSTC);
        v = repmat(double(Lmask),[1 size(tmp,2)]);
        v(Lmask,:) = tmp;
        [~, idxs] = sort(diag(d));
        v = v(:,idxs);
        v = v(:,[end end-1 end-2]);  % Collecting the first 3 principle components
        PCs = cat(3, PCs, v);
    end
    muvect = reshape(repmat([.5 .5 .5],nstixperside^2,1),nstixperside^2*3,1); % creating a 300 x 1 array with each entry as 0.5
    
    % The plotting begins here
    FigHandle = figure(plot_counter);
    set(FigHandle, 'Position', [50, 70, 1449, 705]);
 
    for i = 1:size(STAs,2) % evaluates for each frame
        normfactor = 0.5/((max(abs(STAs(:))))*1.05); % This step is necessary to constrain the values within [-0.5, 0.5]
        STA = normfactor*(STAs(:,i))+muvect; % This makes the values fall back within a range of 0 and 1.
        STA = reshape(STA,[nstixperside nstixperside 3]); % Decomposing the STA into R, G and B plane
        subplot(5,size(STAs,2),i);
        image(STA); % for viewing the image
        set(gca,'XTick',[],'YTick',[]); axis square;
        if (i == 1)
            ylabel('STA');
        end
        for j = 1:size(v,2) % size(v,2) denotes the number of principle components we are interested in looking at
            PC_int = PCs(:,j,i);
            PC_int = 0.5*(PC_int)/(max(abs(PC_int)));
            normfactor = 1;
            PC = normfactor*PC_int+muvect;
            PC = reshape(PC,[nstixperside nstixperside 3]);
            subplot(5,size(STAs,2), j*size(STAs,2)+i);
            image(PC);
            if (i == 1)
                ylabel(strcat('PC',num2str(j)));
            end
            set(gca,'XTick',[],'YTick',[]); axis square;
        end
        STC = reshape(STCs(:,i),[sqrt(length(STCs(:,i))),sqrt(length(STCs(:,i)))]);
        STV = reshape(diag(STC),[nstixperside nstixperside 3]);
        STV = (STV-mean(STV(:)))./(2*range(STV(:)))+.5;
        subplot(6,size(STAs,2),5*size(STAs,2)+i);
        image(STV);
        set(gca,'XTick',[],'YTick',[]); axis square;
        if (i == 1)
            ylabel('STV');
        end
    end
    plot_counter = plot_counter + 1;
    
    %     The code below was provided by GH to enhance the contrast of the
    %     mask. Used it for my rotation talk, 2015.
    %     timecourse = sum(STAs.^2);
    %     weightedSTA = STAs*timecourse';
    %     weightedSTA = (1./(1+exp(-weightedSTA/500000)))-.5;
    %     normfactor = 0.5/((max(abs(weightedSTA)))*1.05);
    %     figure(plot_counter); axes;
    %     image(reshape(normfactor*weightedSTA+muvect,10,10,3))
    %     axis square; plot_counter = plot_counter + 1;
end

%%
% Significance testing on the STAs, STVs, broken down by color
noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
L = stro.trial(:,noisetypeidx)==1;
mu1idx = find(strcmp(stro.sum.trialFields(1,:),'mu1'));
mu2idx = find(strcmp(stro.sum.trialFields(1,:),'mu2'));
mu3idx = find(strcmp(stro.sum.trialFields(1,:),'mu3'));
sigma1idx = find(strcmp(stro.sum.trialFields(1,:),'sigma1'));
sigma2idx = find(strcmp(stro.sum.trialFields(1,:),'sigma2'));
sigma3idx = find(strcmp(stro.sum.trialFields(1,:),'sigma3'));
muvect = unique(stro.trial(L,[mu1idx mu2idx mu3idx]),'rows')/1000;
sigmavect = unique(stro.trial(L,[sigma1idx sigma2idx sigma3idx]),'rows')/1000;
sigmavect(all(any(sigmavect == 0),2),:) = [];
gausslims = [stro.sum.exptParams.gauss_locut stro.sum.exptParams.gauss_hicut]/1000;
mumat = repmat(reshape(repmat(muvect,nstixperside^2,1),[nstixperside^2*3, 1]),[1,size(STS,2)]);
sigmamat = repmat(reshape(repmat(sigmavect,nstixperside^2,1),[nstixperside^2* 3, 1]),[1,size(STS,2)]);
zscoremeans = (STAs-mumat)./(sigmamat/sqrt(nspikes));
% Doing the calculations for the variances
% Only considering one correction factor per dimension (assuming variances on green and blue guns are same as red gun)

NPOINTS = 65536;
x = linspace(gausslims(1),gausslims(2),NPOINTS);
Fx = norminv(x)*sigmavect(1);
sigmacorrectionfactor = std(Fx)./sigmavect(1);
for i = 1:size(STCs,2)
    STC = reshape(STCs(:,i),[sqrt(length(STCs(:,i))),sqrt(length(STCs(:,i)))]);
    STVs(:,i) = diag(STC);
end
muvar = (sigmavect(1)*sigmacorrectionfactor)^2;
sigmavar = muvar*sqrt(2/nspikes); % Check out this page for more info https://groups.google.com/forum/#!topic/sci.stat.math/dsgmWBLJoHc
zscorevars = (STVs-muvar)./sigmavar;
maxzscoremeans = max(abs(zscoremeans(:)));
maxzscorevars = max(abs(zscorevars(:)));
alpha = 0.025;
crit = norminv(1-alpha,0,1);

% Plotting
FigHandle = figure(plot_counter);
set(FigHandle, 'Position', [50, 70, 1449, 705]);
plot_counter = plot_counter + 1;
for i = 1:size(STAs,2)
    % First STAs
    zmat = reshape(zscoremeans(:,i),[nstixperside nstixperside 3]);
    for j = 1:3
        pmat = logical(abs(zmat(:,:,j))>crit);
        im = zmat(:,:,j)./(2*maxzscoremeans)+.5; % ensures that all values lie between 0 and 1
        im = repmat(im,[1 1 3]);
        sigidxs = find(pmat);
        im(sigidxs) = .5;  % red to .5 where sig.  Looks red on dark and cyan on bright.
        subplot(6,size(STS,2),(j-1)*size(STS,2)+i);
        image(im);
        axis image;
        if (i == 1)
            if (j == 1)
                ylabel('STA-R');
            elseif (j == 2)
                ylabel('STA-G');
            elseif (j == 3)
                ylabel('STA-B');
            end
        end
        set(gca,'XTick',[],'YTick',[]);
    end
    % Then STVs
    zmat = reshape(zscorevars(:,i),[nstixperside nstixperside 3]);
    for j = 4:6
        pmat = logical(abs(zmat(:,:,j-3))>crit);
        im = zmat(:,:,j-3)./(2*maxzscorevars)+.5; % ensures that all values lie between 0 and 1
        im = repmat(im,[1 1 3]);
        sigidxs = find(pmat);
        im(sigidxs) = .5;
        subplot(6,size(STCross,2),(j-1)*size(STCross,2)+i);
        image(im);
        axis image;
        if (i == 1)
            if (j == 4)
                ylabel('STV-R');
            elseif (j == 5)
                ylabel('STV-G');
            elseif (j == 6)
                ylabel('STV-B');
            end
        end
        set(gca,'XTick',[],'YTick',[]);
    end
end

%% Derived from STCGUI.m in Greg's Whitenoise folder

if (strcmp(scriptname,'WN_nosubunit')==1 || strcmp(scriptname,'WNSubunit_ST')==1)
    normfactor = 0.5/((max(abs(STAs(:))))*1.05);
    muvect = reshape(repmat([.5 .5 .5],nstixperside^2,1),nstixperside^2*3,1);
    FigHandle = figure(plot_counter); plot_counter = plot_counter + 1;
    set(FigHandle, 'Position', [50, 250, 1449, 350]);
    for j = 1:maxT % evaluates for each frame
        STA1 = normfactor*(STAs(:,j))+muvect; % This makes the values fall back within a range of 0 and 1.
        STA1 = reshape(STA1,[nstixperside nstixperside 3]);
        subplot(1,maxT,j);
        image(STA1); % for viewing the image
        set(gca,'XTick',[],'YTick',[]); axis square;
        if (j == 1)
            ylabel('STA');
        end
    end
    [plot_counter] = Selectpixel(plot_counter,STAs); % function call
    plot_counter = plot_counter + 1;
    
end
    
end



