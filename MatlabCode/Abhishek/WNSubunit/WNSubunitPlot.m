function [STAs, plot_counter] = WNSubunitPlot(stro,plot_counter,mask_changes)
% Abhishek's copy of WNSubunitPLot.m

global spikename maskidx spikeidx nstixperside ngammasteps seedidx nframesidx stimonidx muidxs sigmaidxs
global msperframe ntrials maxT xx yy

%disp(mask_changes);
% STCOVmex - This function calculates the covariance matrices, efficiently
% and quickly. Function takes 2 arguments. The string 'init' as the first
% argument and the number of stimulus dimensions as the second argument in
% the format [npics, ncols, nframes].
% In subsequent cells, send as the first argument a vector of RGBs (the
% total number of elements should be xpix*ncols*nframes, and the order
% shoud be first all the elemets for frame 1, then all the elements for
% frame 2, and so on). The second argument should be an nframesx1 matrix of
% spike counts.
% Send the 'return' string as a solo argument to get the output - a
% structure containing {STS, STCross, and n}

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
        frametimes = linspace(0, nframes*msperframe, nframes)+(msperframe/2)';
        spiketimes(spiketimes < maxT*msperframe) = [];
        spiketimes(spiketimes > frametimes(end)) = [];
        n = hist(spiketimes, frametimes);
        STCOVmex(rgbs(:),n);
    end
    out = STCOVmex('return');
    STAs = out{1};
    STCs = out{2};
    nspikes = out{3};
    clear STCOVmex;
    clear out;
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
        %         figure(5); plot(subSTA); pause(0.2);
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
    figure(plot_counter);
    for i = 1:size(STAs,2) % evaluates for each frame
        normfactor = 0.5/((max(abs(STAs(:))))*1.01); % This step is necessary to constrain the values within [-0.5, 0.5]
        STA = normfactor*(STAs(:,i)-muvect)+muvect; % This makes the values fall back within a range of 0 and 1.
        STA = reshape(STA,[nstixperside nstixperside 3]); % Decomposing the STA into R, G and B plane
        subplot(5,size(STAs,2),i);
        image(STA); % for viewing the image
        disp(i);
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
end