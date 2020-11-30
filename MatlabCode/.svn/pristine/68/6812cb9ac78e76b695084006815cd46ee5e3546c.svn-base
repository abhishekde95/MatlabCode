function GridLMSubunits_ParadigmTest()
global gl

%This code is for modeling sparse to fine sampling of RF structure.

% 2/12/13   Created.    JPW
% 3/19/13   Adapted to test actual online paradigm



%% Set Up Variables

% User-defined variables
gl.CurrentQueue = [];
gl.epoch = 1; % 1 = getting timing/threshold/square responses, 2 = stixel mapping, 3 = gridLMplane
gl.rnd = 1; % for "phases" within a given epoch (epoch 1 will have only 1 round, epoch 2 will have ~3 rounds, and epoch 3 has an indefinite number)
gl.initNPts = 5;
gl.nRnds = 3;
gl.nStixelGrid = (gl.initNPts - 1) * 2 .^ (gl.nRnds - 1) + 1;
gl.DVAPerStixel = .1;
gl.StimDuration = .1; %In second
gl.nPres = 1;
gl.square = [];
gl.stixels = [];
gl.gridX = [];
gl.gridY = [];
%abortflag = 0;
gl.pixelmask = zeros(gl.nStixelGrid, gl.nStixelGrid);


% Grid for probing
%masterGrid1D_C = -(nPts1D_MG_NC-1)/2:(nPts1D_MG_NC-1)/2;
%[masterGrid_C_c masterGrid_C_r] = meshgrid(masterGrid1D_C,fliplr(masterGrid1D_C));
%responses = nan(size(masterGrid_C_c));
%ptsInRnd_test_MG_C_c = [];
%ptsInRnd_test_MG_C_r = [];

%Set Xlim and Ylim to values that can be called for all figures
%xLim = [min(masterGrid1D_C)*1.1 max(masterGrid1D_C)*1.1];
%yLim = [min(masterGrid1D_C)*1.1 max(masterGrid1D_C)*1.1];

%% Actual receptive field (delete this in real experiment)

% Super fine grid just for reference
X = linspace(-(gl.nStixelGrid-1)/2,(gl.nStixelGrid-1)/2,100);
[xN yN] = meshgrid(X,X);
ZZ = RFModel(yN,xN,gl.nStixelGrid,gl.DVAPerStixel,0);

masterGrid1D_C = -(gl.nStixelGrid-1)/2:(gl.nStixelGrid-1)/2;
xLim = [min(masterGrid1D_C)*1.1 max(masterGrid1D_C)*1.1];
yLim = [min(masterGrid1D_C)*1.1 max(masterGrid1D_C)*1.1];


% Plot receptive field for reference
figure(1); clf; hold on; grid on;
xlim(xLim)
ylim(yLim)
surf(xN,yN,ZZ)
set(gca,'XTick',linspace(min(X),max(X),5))
set(gca,'YTick',linspace(min(X),max(X),5))
set(gca,'XTickLabel',linspace(-.5,.5,numel(get(gca,'XTick'))))
set(gca,'YTickLabel',linspace(-.5,.5,numel(get(gca,'YTick'))))
xlabel('DVA from Mapped RF center')
ylabel('DVA from Mapped RF center')
title('"Real" Underlying Receptive Field')

% if stops == 1
%     keyboard
% end

%% Pre round screening (large squares)


    disp('In Epoch 1: Testing Squares')
    LM = [.05 .05; .05 -.05; -.05 -.05; -.05 .05];


    gl.masterGrid1D_C = -(gl.nStixelGrid-1)/2:(gl.nStixelGrid-1)/2;
    [gl.masterGrid_C_c gl.masterGrid_C_r] = meshgrid(gl.masterGrid1D_C,fliplr(gl.masterGrid1D_C));
    gl.responses = nan(size(gl.masterGrid_C_c));    
    if isempty(gl.square)
        gl.square = SetUpSquares(gl.square);
    end
    % Should the stuff above be in "SetUpSquares"?

    % Large Square Pre-Screening
    if isempty(gl.CurrentQueue) && gl.epoch == 1
        
        
        % Set up queue ([L M],[x y])
        LMxyIdx = fullfact([size(LM,1) numel(gl.square)]);
        
        % Randomize order for presentation
        LMxyIdx_r = LMxyIdx(randperm(size(LMxyIdx,1)),:);        
        
        % Convert parameters that are not changing from single values to vectors
        S = zeros(size(LMxyIdx,1),1);
        param_nstixgrid = repmat(gl.nStixelGrid,size(LMxyIdx,1),1);
        param_dvaperstix = repmat(gl.DVAPerStixel,size(LMxyIdx,1),1);
        param_stimdur = repmat(gl.StimDuration,size(LMxyIdx,1),1);
        
        % Assign values to CurrentQueue
        % ([L,M,S,nstixgrid,dvaperstix,stimdur])
        gl.CurrentQueue = [LM(LMxyIdx_r(:,1),:) S param_nstixgrid param_dvaperstix param_stimdur];
        for n = 1:size(LMxyIdx_r,1)
            gl.gridX(n,:) = gl.square{LMxyIdx_r(n,2)}.c(:);
            gl.gridY(n,:) = gl.square{LMxyIdx_r(n,2)}.r(:);
        end
        
        % Lastly, set gl.epoch to 2 (assuming the user wants to advance to round 1)
        % to begin looping through the rounds. At the moment, only checking
        % when the CurrentQueue is empty. 
        a = get(gcf,'UserData');
        if (get(a.uicontrols.startmapping,'Value') == 1)
            gl.square = SummarizeSquareResponses(gl.square);
            gl.CurrentQueue = [];
            gl.epoch = 2;
            gl.rnd = 1;
        end
        
    end



squareEdges = linspace(min(masterGrid1D_C),max(masterGrid1D_C),gl.initNPts)...
    +max(masterGrid1D_C)+1;
sqnum = 1;
for sqrow = 1:numel(squareEdges)-1
    for sqcol = 1:numel(squareEdges)-1
        square(sqnum).r = masterGrid_C_r(squareEdges(sqrow):squareEdges(sqrow+1),...
            squareEdges(sqrow):squareEdges(sqrow+1));
        square(sqnum).c = masterGrid_C_c(squareEdges(sqcol):squareEdges(sqcol+1),...
            squareEdges(sqcol):squareEdges(sqcol+1));
        sqnum = sqnum+1;
    end
end

% Set up figure
figure; hold on; grid on;

% Go through squares and test for responses.
tempCounter = 0;
tempResp = nan(nPres,1);
for sq = 1:numel(square)
    
    % Responses
    % **NOTE: I am currently taking the mean of the response at each point.
    % This may not reflect the properties of real nuerons, so consider
    % alternate strategies.  Using the sum of large areas will yield a
    % response at all locations, even those without an underlying subunit.
    for rep = 1:nPres
        tempResp(rep) = mean(mean(RFModel(square(sq).r,square(sq).c,nPts1D_MG_NC,DVAperStix,addnoise)));
    end
    sqResp = mean(tempResp);
    %sqResp = RFModel(square(sq).c,square(sq).r,0);
    
    if sum(sum(sqResp)) > threshold
        
        square(sq).test = 1;
        tempCounter = tempCounter + 1;
        
        % Plot square boundaries
        %         upperLeft_rc(tempCounter,:) = [square(sq).r(1,1) square(sq).c(1,1)] + [1 -1];
        %         lowerLeft_rc(tempCounter,:) = [square(sq).r(end,1) square(sq).c(end,1)] + [-1 -1];
        %         upperRight_rc(tempCounter,:) = [square(sq).r(1,end) square(sq).c(1,end)] + [1 1];
        %         lowerRight_rc(tempCounter,:) = [square(sq).r(end,end) square(sq).c(end,end)] + [-1 1];
        upperLeft_rc = [square(sq).r(1,1) square(sq).c(1,1)] + [1 -1];
        lowerLeft_rc = [square(sq).r(end,1) square(sq).c(end,1)] + [-1 -1];
        upperRight_rc = [square(sq).r(1,end) square(sq).c(1,end)] + [1 1];
        lowerRight_rc = [square(sq).r(end,end) square(sq).c(end,end)] + [-1 1];
        
        square(sq).topline = [upperLeft_rc; upperRight_rc]';
        square(sq).rightline = [upperRight_rc; lowerRight_rc]';
        square(sq).bottomline = [lowerRight_rc; lowerLeft_rc]';
        square(sq).leftline = [lowerLeft_rc; upperLeft_rc]';
        
        square(sq).perimeter = cat(2,square(sq).topline,square(sq).rightline,...
            square(sq).bottomline,square(sq).leftline);
        
        plot(square(sq).perimeter(2,:),square(sq).perimeter(1,:),'k--')
        plot(square(sq).c,square(sq).r,'co')
        legend('Area spanned by squares that evoked a response')
        
    else
        square(sq).test = 0;
    end
    
end

% Plot "real" underlying RF contours with boxes over it
contour(xN,yN,ZZ)
xlim(xLim)
ylim(yLim)
set(gca,'XTick',linspace(min(X),max(X),5))
set(gca,'YTick',linspace(min(X),max(X),5))
set(gca,'XTickLabel',linspace(-.5,.5,numel(get(gca,'XTick'))))
set(gca,'YTickLabel',linspace(-.5,.5,numel(get(gca,'YTick'))))
xlabel('DVA from Mapped RF center')
ylabel('DVA from Mapped RF center')
title('RF Pre-Screening Using Large Partially Overlapping Squares')

if stops == 1
    keyboard
end

%% Simulate the Procession through Rounds


disp('Selecting Next Stimulus...')
LM = [.05 .05; .05 -.05; -.05 -.05; -.05 .05];

% Set up squares
% Taking this outside the
% "if isempty(gl.CurrentQueue) && gl.epoch == 1"
% clause so that we can set up the square structure artificially for
% debugging. This change should not affect the normal function of the
% script. GDLH 3/17/13
% Set up grid for probing and tracking responses
gl.masterGrid1D_C = -(gl.nStixelGrid-1)/2:(gl.nStixelGrid-1)/2;
[gl.masterGrid_C_c gl.masterGrid_C_r] = meshgrid(gl.masterGrid1D_C,fliplr(gl.masterGrid1D_C));
gl.responses = nan(size(gl.masterGrid_C_c));    
if isempty(gl.square)
    gl.square = SetUpSquares(gl.square);
end
% Should the stuff above be in "SetUpSquares"?

% Large Square Pre-Screening
if isempty(gl.CurrentQueue) && gl.epoch == 1
    
    
    % Set up queue ([L M],[x y])
    LMxyIdx = fullfact([size(LM,1) numel(gl.square)]);
    
    % Randomize order for presentation
    LMxyIdx_r = LMxyIdx(randperm(size(LMxyIdx,1)),:);
    
    % Convert parameters that are not changing from single values to vectors
    S = zeros(size(LMxyIdx,1),1);
    param_nstixgrid = repmat(gl.nStixelGrid,size(LMxyIdx,1),1);
    param_dvaperstix = repmat(gl.DVAPerStixel,size(LMxyIdx,1),1);
    param_stimdur = repmat(gl.StimDuration,size(LMxyIdx,1),1);
    
    % Assign values to CurrentQueue
    % ([L,M,S,nstixgrid,dvaperstix,stimdur])
    gl.CurrentQueue = [LM(LMxyIdx_r(:,1),:) S param_nstixgrid param_dvaperstix param_stimdur];
    for n = 1:size(LMxyIdx_r,1)
        gl.gridX(n,:) = gl.square{LMxyIdx_r(n,2)}.c(:);
        gl.gridY(n,:) = gl.square{LMxyIdx_r(n,2)}.r(:);
    end
    
    % Lastly, set gl.epoch to 2 (assuming the user wants to advance to round 1)
    % to begin looping through the rounds. At the moment, only checking
    % when the CurrentQueue is empty.
    a = get(gcf,'UserData');
    if (get(a.uicontrols.startmapping,'Value') == 1)
        gl.square = SummarizeSquareResponses(gl.square);
        gl.CurrentQueue = [];
        gl.epoch = 2;
        gl.rnd = 1;
    end
    
end
end

%% Epoch 2
if isempty(gl.CurrentQueue) && gl.epoch == 2
    
    % Set up structure for stixels
    if isempty(gl.stixels)
        
        gl.stixels = SetUpStixels(gl.stixels);
        LMrcQueue = [];
        
    end
    
    fieldnames = {'pLpM','pLmM','mLpM','mLmM'};
    if gl.rnd == 1
        
        for f = 1:numel(fieldnames)
            
            for sq = 1:numel(gl.square)
                
                if gl.square{sq}.(fieldnames{f}).subsample == 1
                    
                    % Use only corners of the squares as the initial grid.
                    upperLeft = [gl.square{sq}.r(1,1) gl.square{sq}.c(1,1)];
                    lowerLeft = [gl.square{sq}.r(end,1) gl.square{sq}.c(end,1)];
                    upperRight = [gl.square{sq}.r(1,end) gl.square{sq}.c(1,end)];
                    lowerRight = [gl.square{sq}.r(end,end) gl.square{sq}.c(end,end)];
                    
                    gl.stixels.(fieldnames{f}).testing.r = cat(1,...
                        gl.stixels.(fieldnames{f}).testing.r,...
                        [upperLeft(1); lowerLeft(1); upperRight(1); lowerRight(1)]);
                    gl.stixels.(fieldnames{f}).testing.c = cat(1,...
                        gl.stixels.(fieldnames{f}).testing.c,...
                        [upperLeft(2); lowerLeft(2); upperRight(2); lowerRight(2)]);
                    
                    % Keep track of all eligible points (points comprised
                    % by responsive squares)
                    gl.stixels.(fieldnames{f}).eligible.r = gl.square{sq}.r;
                    gl.stixels.(fieldnames{f}).eligible.c = gl.square{sq}.c;
                    
                end
                
            end
            
            % Delete repeated stimuli
            try
                rc = unique([gl.stixels.(fieldnames{f}).testing.r gl.stixels.(fieldnames{f}).testing.c],'rows');
                gl.stixels.(fieldnames{f}).testing.r = rc(:,1);
                gl.stixels.(fieldnames{f}).testing.c = rc(:,2);
            catch
                disp('Something wrong with indexing here...')
                keyboard
            end
            
            % Set up queue ([L M],[x y])
            if strcmp(fieldnames{f},'pLpM')
                LM = [.05 .05];
            elseif strcmp(fieldnames{f},'mLmM')
                LM = [-.05 -.05];
            elseif strcmp(fieldnames{f},'pLmM')
                LM = [.05 -.05];
            elseif strcmp(fieldnames{f},'mLpM')
                LM = [-.05 .05];
            end
            
            LMrcQueue = cat(1,LMrcQueue,...
                [repmat(LM,size(gl.stixels.(fieldnames{f}).testing.r,1),1)...
                gl.stixels.(fieldnames{f}).testing.r...
                gl.stixels.(fieldnames{f}).testing.c]);
            
        end
        
        % Repeat stimuli nPres times
        LMrcQueue = repmat(LMrcQueue,gl.nPres,1);
        LMrcQueue = LMrcQueue(randperm(size(LMrcQueue,1)),:);
        
        % Set up Current Queue, gridX, and gridY
        gl.gridX = LMrcQueue(:,4);
        gl.gridY = LMrcQueue(:,3);
        param_S = zeros(size(LMrcQueue,1),1);
        param_nstixelgrid = repmat(gl.nStixelGrid,size(LMrcQueue,1),1);
        param_dvaperstixel = repmat(gl.DVAPerStixel,size(LMrcQueue,1),1);
        param_stimdur = repmat(gl.StimDuration,size(LMrcQueue,1),1);
        gl.CurrentQueue = [LMrcQueue(:,1:2) param_S param_nstixelgrid...
            param_dvaperstixel param_stimdur];
        
        % Supplement the round ticker
        gl.rnd = gl.rnd+1;
        
    elseif gl.rnd > 1
        
        ptsInRnd = (gl.initNPts - 1) * 2 .^ (gl.rnd - 1) + 1;
        ptsInterval = (gl.nStixelGrid-1)/(ptsInRnd-1);
        crossring = [ptsInterval 0; 0 ptsInterval; -ptsInterval 0; 0 -ptsInterval];
        diagring = [ptsInterval ptsInterval; -ptsInterval ptsInterval;...
            ptsInterval -ptsInterval; -ptsInterval -ptsInterval];
        LMrcQueue = [];
        outofrangeI = [];
        
        % Transfer response information from figure to stixel structure
        a = get(gcf,'UserData');
        gl.threshold = a.uicontrols.threshslider;
        
        % Pull out stixels that have been tested
        b = get(a.axes.lm(1),'UserData');
        
        %pLpM
        stixIdx = find(b.nsp(:,:,1)>=0);
        tempfield = 'pLpM';
        [tempr tempc] = ind2sub(size(b.nsp(:,:,1)),stixIdx);
        gl.stixels.(tempfield).tested.r = tempr + min(gl.masterGrid1D_C) - 1;
        gl.stixels.(tempfield).tested.c = tempc + min(gl.masterGrid1D_C) - 1;
        tempIdx = sub2ind(size(b.nsp),tempr,tempc,ones(size(tempr)));
        gl.stixels.(tempfield).tested.fr = b.nsp(tempIdx)./b.n(tempIdx);
        
        %mLmM
        stixIdx = find(b.nsp(:,:,2)>=0);
        tempfield = 'mLmM';
        [tempr tempc] = ind2sub(size(b.nsp(:,:,2)),stixIdx);
        gl.stixels.(tempfield).tested.r = tempr + min(gl.masterGrid1D_C) - 1;
        gl.stixels.(tempfield).tested.c = tempc + min(gl.masterGrid1D_C) - 1;
        tempIdx = sub2ind(size(b.nsp),tempr,tempc,repmat(2,size(tempr)));
        gl.stixels.(tempfield).tested.fr = b.nsp(tempIdx)./b.n(tempIdx);
        
        %pLmM
        b = get(a.axes.lVm(1),'UserData');
        stixIdx = find(b.nsp(:,:,1)>=0);
        tempfield = 'pLmM';
        [tempr tempc] = ind2sub(size(b.nsp(:,:,1)),stixIdx);
        gl.stixels.(tempfield).tested.r = tempr + min(gl.masterGrid1D_C) - 1;
        gl.stixels.(tempfield).tested.c = tempc + min(gl.masterGrid1D_C) - 1;
        tempIdx = sub2ind(size(b.nsp),tempr,tempc,ones(size(tempr)));
        gl.stixels.(tempfield).tested.fr = b.nsp(tempIdx)./b.n(tempIdx);
        
        %mLpM
        stixIdx = find(b.nsp(:,:,2)>=0);
        tempfield = 'mLpM';
        [tempr tempc] = ind2sub(size(b.nsp(:,:,2)),stixIdx);
        gl.stixels.(tempfield).tested.r = tempr + min(gl.masterGrid1D_C) - 1;
        gl.stixels.(tempfield).tested.c = tempc + min(gl.masterGrid1D_C) - 1;
        tempIdx = sub2ind(size(b.nsp),tempr,tempc,repmat(2,size(tempr)));
        gl.stixels.(tempfield).tested.fr = b.nsp(tempIdx)./b.n(tempIdx);
        
        for f = 1:numel(fieldnames)
            
            ptsInRnd_noTest_rc = [];
            
            for st = 1:size(gl.stixels.(fieldnames{f}).tested.r,1)
                
                % For any stixel that prevously responded, entertain testing surrounding stixels
                if gl.stixels.(fieldnames{f}).tested.fr(st) > gl.threshold
                    
                    % Quick reference for current stixel (bc the
                    % strucure is large)
                    currentstixel = [gl.stixels.(fieldnames{f}).tested.r(st) gl.stixels.(fieldnames{f}).tested.c(st)];
                    
                    % Add "ring" to responsive points - these will be
                    % considered for the upcoming round.
                    temptest = repmat(currentstixel,4,1) + crossring;
                    
                    % If any point appears twice on this list, it
                    % appears between two responsive points, and we
                    % now assume that it will also be responsive.  This
                    % stixel is taken off the list of stixels to be
                    % tested.
                    tempnotest = intersect(temptest,[gl.stixels.(fieldnames{f}).testing.r gl.stixels.(fieldnames{f}).testing.c],'rows');
                    
                    if ~isempty(tempnotest)
                        
                        
                        ptsInRnd_noTest_rc = cat(1,ptsInRnd_noTest_rc,tempnotest);
                        
                        % Add only points that don't occur twice to the
                        % upcoming queue.
                        temptest = setxor(temptest,[gl.stixels.(fieldnames{f}).testing.r(st) gl.stixels.(fieldnames{f}).testing.c(st)],'rows');
                        
                    end
                    
                    % Eliminate points outside of the stixel range
                    outofrangeI = find(temptest(:,1) > max(gl.masterGrid1D_C));
                    outofrangeI = cat(1,outofrangeI,find(temptest(:,1) < min(gl.masterGrid1D_C)));
                    outofrangeI = cat(1,outofrangeI,find(temptest(:,2) > max(gl.masterGrid1D_C)));
                    outofrangeI = cat(1,outofrangeI,find(temptest(:,2) < min(gl.masterGrid1D_C)));
                    if any(outofrangeI)
                        disp('WARNING: responsive stixels at the endge of the RF')
                        try
                            temptest(outofrangeI,:) = [];
                        catch
                            keyboard
                        end
                    end
                    
                    % Add points the don't occur twice into the "testing" field.
                    gl.stixels.(fieldnames{f}).testing.r = cat(1,gl.stixels.(fieldnames{f}).testing.r,temptest(:,1));
                    gl.stixels.(fieldnames{f}).testing.c = cat(1,gl.stixels.(fieldnames{f}).testing.c,temptest(:,2));
                    
                end
                
            end
            
            % Clear out some fields
            temptest = [];
            outofrangeI = [];
            
            % To all points that responded previously, add diagonal points
            % as potential candidates for next round
            for st = 1:size(gl.stixels.(fieldnames{f}).tested.r,1)
                
                % Quick reference for current stixel (bc the strucure is large)
                currentstixel = [gl.stixels.(fieldnames{f}).tested.r(st) gl.stixels.(fieldnames{f}).tested.c(st)];
                
                % Add "ring" to responsive points - these will be considered for the upcoming round.
                temptest = cat(1,temptest,repmat(currentstixel,4,1) + diagring);
                
            end
            
            % If a diagonal point appears 4 times in the queue, it is
            % completely surrounded by 4 responsive points, and is
            % assumed to also be responsive.
            if ~isempty(temptest)
                [ptsInRnd_test_rc ia ib] = unique(temptest,'rows');
                tempIB = unique(ib);
                getridof = [];
                for u = 1:numel(tempIB)
                    nrep = sum(tempIB(u)==ib);
                    if nrep >= 4
                        ptsInRnd_noTest_rc = cat(1,ptsInRnd_noTest_rc,...
                            ptsInRnd_test_rc(tempIB(u),:));
                        getridof = cat(1,getridof,tempIB(u));
                    end
                end
                ptsInRnd_test_rc(getridof,:) = [];
            end
            
            
            % Eliminate points outside of the stixel range
            outofrangeI = cat(1,outofrangeI,find(ptsInRnd_test_rc(:,1) > max(gl.masterGrid1D_C)));
            outofrangeI = cat(1,outofrangeI,find(ptsInRnd_test_rc(:,1) < min(gl.masterGrid1D_C)));
            outofrangeI = cat(1,outofrangeI,find(ptsInRnd_test_rc(:,2) > max(gl.masterGrid1D_C)));
            outofrangeI = cat(1,outofrangeI,find(ptsInRnd_test_rc(:,2) < min(gl.masterGrid1D_C)));
            if any(outofrangeI)
                disp('WARNING: responsive stixels at the endge of the RF')
            end
            ptsInRnd_test_rc(outofrangeI,:) = [];
            
            gl.stixels.(fieldnames{f}).testing.r = ptsInRnd_test_rc(:,1);
            gl.stixels.(fieldnames{f}).testing.c = ptsInRnd_test_rc(:,2);
            
            % Put "no test" stixels into the stixel structure, and
            % assign an indicative fr and nspikes
            ptsInRnd_noTest_rc = unique(ptsInRnd_noTest_rc,'rows');
            tempnspikes = ones(size(ptsInRnd_noTest_rc,1),1) * 1000;
            if any(ptsInRnd_noTest_rc)
                tempfr = ones(size(tempnspikes)) * 1000;
                gl.stixels.(fieldnames{f}).tested.r = cat(1,gl.stixels.(fieldnames{f}).tested.r,ptsInRnd_noTest_rc(:,1));
                gl.stixels.(fieldnames{f}).tested.c = cat(1,gl.stixels.(fieldnames{f}).tested.c,ptsInRnd_noTest_rc(:,2));
                gl.stixels.(fieldnames{f}).tested.fr = cat(1,gl.stixels.(fieldnames{f}).tested.fr,tempfr);
                
                %Must put "no test" points into figure 'user data'
                assumeResp_r = ptsInRnd_noTest_rc(:,1) + max(gl.masterGrid1D_C) + 1;
                assumeResp_c = ptsInRnd_noTest_rc(:,2) + max(gl.masterGrid1D_C) + 1;
                assumeRespIdx = sub2ind([size(b.nsp,1),size(b.nsp,2),1],assumeResp_r,assumeResp_c)
                
                if f == 1
                    handle = 'lm';
                    fig = 1;
                elseif f == 2
                    handle = 'lm';
                    fig = 2;
                elseif f == 3
                    handle = 'lVm'
                    fig = 1;
                else
                    handle = 'lVm'
                    fig = 2;
                end
                
                a = get(gcf,'UserData');
                b = get(a.axes.(handle)(1),'UserData');
                tempresps = b.nsp(:,:,fig);
                tempresps(assumeRespIdx) = 1000;
                b.nsp(:,:,fig) = tempresps;
                set(a.axes.(handle)(1),'UserData',b);
                
            end
            
            
            % Put stixels into a temp queue to be randomized later ([L M],[x y])
            if strcmp(fieldnames{f},'pLpM')
                LM = [.05 .05];
            elseif strcmp(fieldnames{f},'mLmM')
                LM = [-.05 -.05];
            elseif strcmp(fieldnames{f},'pLmM')
                LM = [.05 -.05];
            elseif strcmp(fieldnames{f},'mLpM')
                LM = [-.05 .05];
            end
            
            LMrcQueue = cat(1,LMrcQueue,...
                [repmat(LM,size(gl.stixels.(fieldnames{f}).testing.r,1),1)...
                gl.stixels.(fieldnames{f}).testing.r...
                gl.stixels.(fieldnames{f}).testing.c]);
            
        end
        
        % Repeat stimuli nPres times
        LMrcQueue = repmat(LMrcQueue,gl.nPres,1);
        LMrcQueue = LMrcQueue(randperm(size(LMrcQueue,1)),:);
        
        % Set up Current Queue, gridX, and gridY
        gl.gridX = LMrcQueue(:,4);
        gl.gridY = LMrcQueue(:,3);
        param_S = zeros(size(LMrcQueue,1),1);
        param_nstixelgrid = repmat(gl.nStixelGrid,size(LMrcQueue,1),1);
        param_dvaperstixel = repmat(gl.DVAPerStixel,size(LMrcQueue,1),1);
        param_stimdur = repmat(gl.StimDuration,size(LMrcQueue,1),1);
        gl.CurrentQueue = [LMrcQueue(:,1:2) param_S param_nstixelgrid...
            param_dvaperstixel param_stimdur];
        
        
        if gl.rnd == gl.nRnds
            gl.epoch = 3;
            gl.rnd = 1;
        else
            gl.rnd = gl.rnd + 1;
        end
        
        keyboard
        
    end
    
end


% Send to REX to be tested.
tempResp = nan(numel(ptsInRnd_test_MG_C_c),nPres);
for rep = 1:nPres
    tempResp(:,rep) = RFModel(ptsInRnd_test_MG_C_r,ptsInRnd_test_MG_C_c,nPts1D_MG_NC,DVAperStix,addnoise);
end
responses(ptsInRnd_test_MG_I) = mean(tempResp,2);
%responses(ptsInRnd_test_MG_I) = RFModel(ptsInRnd_test_MG_C_c,ptsInRnd_test_MG_C_r);





% Plot the "real" underlying RF overlayed by tested and untested
% points.
figure; clf; hold on; grid on;
leg = [];
plot(ptsInRnd_test_MG_C_c,ptsInRnd_test_MG_C_r,'*r')
xlim(xLim);
ylim(yLim);
leg{1} = 'Current Stimuli';
forgonePts = find(responses==1000);

% Plot points that are presumed to be responsive
if ~isempty(forgonePts)
    plot(masterGrid_C_c(forgonePts),masterGrid_C_r(forgonePts),'*')
    leg{end+1} = 'Presumed Responsive';
end

% Plot points that previously gave a response
if r>1 && nansum(nansum(ptsInRnd_gaveResp_L))
    ptsInRnd_gaveResp_MG_I = ptsInRnd_MG_I(ptsInRnd_gaveResp_L);
    ptsInRnd_gaveResp_MG_C_c = masterGrid_C_c(ptsInRnd_gaveResp_MG_I);
    ptsInRnd_gaveResp_MG_C_r = masterGrid_C_r(ptsInRnd_gaveResp_MG_I);
    plot(ptsInRnd_gaveResp_MG_C_c,ptsInRnd_gaveResp_MG_C_r,'g*')
    leg{end+1} = 'Previously Responded';
end


% Plot squares from pre-screening
for sq = 1:numel(square)
    if square(sq).test == 1
        plot(square(sq).perimeter(2,:),square(sq).perimeter(1,:),'k--')
    end
end
leg{end+1} = 'Square Perimeters';

% Plot current grid
plot(ptsInRnd_MG_C_c,ptsInRnd_MG_C_r,'ok')
leg{end+1} = 'Current grid';

set(gca,'XTick',linspace(min(X),max(X),5))
set(gca,'YTick',linspace(min(X),max(X),5))
set(gca,'XTickLabel',linspace(-.5,.5,numel(get(gca,'XTick'))))
set(gca,'YTickLabel',linspace(-.5,.5,numel(get(gca,'YTick'))))
xlabel('DVA from Mapped RF center')
ylabel('DVA from Mapped RF center')
title(['Round #' num2str(r)])

%legend(leg)
contour(xN,yN,ZZ)

if stops == 1
    keyboard
end

end




function [response] = RFModel(testy, testx, nSizeGrid, DVAperStix, addnoise) % Model Neuron

if nargin == 4
    addnoise = 1;
end

cellBaseline = 1;


% Add noise to eye position (Gaussian distributed)
if addnoise == 1
    eyePosNoise = randn(1,2)./12;
    eyePosNoise(eyePosNoise>.2) = .2;
    eyePosNoise(eyePosNoise<-.2) = -.2;
    eyePosNoise = eyePosNoise./DVAperStix;
    testx = testx + eyePosNoise(1);
    testy = testy + eyePosNoise(2);
end

% I'm thinking of this area as roughly 1DVA.^2
OnOffUnit = (cos(testy./(nSizeGrid).*pi).*sin(testx./(nSizeGrid/2).*pi));
height = 50;
mean = 3;
std = (nSizeGrid-1)/6;
Gaussian = height.*exp(-((testx-mean).^2 ./ (2*std.^2)) - ((testy-mean).^2 ./ (2*std.^2)));

if addnoise == 1
    response = poissrnd(OnOffUnit .* Gaussian + cellBaseline);
else
    response = OnOffUnit .* Gaussian;
end


response(isnan(response)) = 0;

end


function square = SetUpSquares(square)
global gl        

    if (isempty(square))
        squareEdges = linspace(min(gl.masterGrid1D_C),max(gl.masterGrid1D_C),gl.initNPts)...
            +max(gl.masterGrid1D_C)+1;
        sqnum = 1;
        square = cell((numel(squareEdges)-1).^2,1);
        for sqcol = 1:numel(squareEdges)-1
            for sqrow = 1:numel(squareEdges)-1
                
                % Set up [r,c] values for squares
                square{sqnum}.r = gl.masterGrid_C_r(squareEdges(sqrow):squareEdges(sqrow+1),...
                    squareEdges(sqrow):squareEdges(sqrow+1));
                square{sqnum}.c = gl.masterGrid_C_c(squareEdges(sqcol):squareEdges(sqcol+1),...
                    squareEdges(sqcol):squareEdges(sqcol+1));
                
                % Set up fields for responses
                square{sqnum}.pLpM.spiketimes = cell([0,0]);
                square{sqnum}.mLmM.spiketimes = cell([0,0]);
                square{sqnum}.pLmM.spiketimes = cell([0,0]);
                square{sqnum}.mLpM.spiketimes = cell([0,0]);
                
                square{sqnum}.pLpM.fr = [];
                square{sqnum}.mLmM.fr = [];
                square{sqnum}.pLmM.fr = [];
                square{sqnum}.mLpM.fr = [];
                
                
                sqnum = sqnum+1;
            end
        end
    end
end







