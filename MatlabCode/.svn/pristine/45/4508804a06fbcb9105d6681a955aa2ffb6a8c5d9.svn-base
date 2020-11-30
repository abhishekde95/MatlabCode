function GridLMSubunits_ParadigmTest()

%This code is for modeling sparse to fine sampling of RF structure.

% 2/12/13   Created.    JPW
clear all


%% Set Up Variables

% User-defined variables
initNPts = 5;
threshold = 5;
nPres = 3;
DVAperStix = .1;
stops = 0;
nPts1D_MG_NC = 17;
addNoise = 0;
S = 0;


% Calculating total number of rounds
pts = initNPts;
nRnds = 1;
while pts < nPts1D_MG_NC
    pts = (pts)*2-1;
    nRnds = nRnds+1;
end
if pts ~= nPts1D_MG_NC
    disp('Warning: initial number of points and final number of points do not agree!')
    keyboard
end


% Grid for probing
masterGrid1D_C = -(nPts1D_MG_NC-1)/2:(nPts1D_MG_NC-1)/2;
[masterGrid_C_c masterGrid_C_r] = meshgrid(masterGrid1D_C,fliplr(masterGrid1D_C));
responses.poslum = nan(size(masterGrid_C_c));
responses.neglum = nan(size(masterGrid_C_c));
ptsInRnd_test_MG_C_c = [];
ptsInRnd_test_MG_C_r = [];

%Set Xlim and Ylim to values that can be called for all figures
xLim = [min(masterGrid1D_C)*1.1 max(masterGrid1D_C)*1.1];
yLim = [min(masterGrid1D_C)*1.1 max(masterGrid1D_C)*1.1];



%% "Real" receptive field (delete this in real experiment)

% Super fine grid just for reference
X = linspace(-(nPts1D_MG_NC-1)/2,(nPts1D_MG_NC-1)/2,100);
[xN yN] = meshgrid(X,X);
ZZ = RFModel(xN,yN,nPts1D_MG_NC,DVAperStix,1,1,S,0);

% Plot receptive field for reference
figure(10); clf; hold on; grid on;
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

if stops == 1
    keyboard
end


%% Pre-round screening (large squares)
fig = 11;

squareEdges = linspace(min(masterGrid1D_C),max(masterGrid1D_C),initNPts)...
    +max(masterGrid1D_C)+1;
sqnum = 1;
for sqrow = 1:numel(squareEdges)-1
    for sqcol = 1:numel(squareEdges)-1
        
        % Set up squares that test +luminance
        square(sqnum).poslum.r = masterGrid_C_r(squareEdges(sqrow):squareEdges(sqrow+1),...
            squareEdges(sqrow):squareEdges(sqrow+1));
        square(sqnum).poslum.c = masterGrid_C_c(squareEdges(sqcol):squareEdges(sqcol+1),...
            squareEdges(sqcol):squareEdges(sqcol+1));
        square(sqnum).poslum.L = 1;
        square(sqnum).poslum.M = 1;
        
        % Set up squares that test -luminance
        square(sqnum).neglum.r = masterGrid_C_r(squareEdges(sqrow):squareEdges(sqrow+1),...
            squareEdges(sqrow):squareEdges(sqrow+1));
        square(sqnum).neglum.c = masterGrid_C_c(squareEdges(sqcol):squareEdges(sqcol+1),...
            squareEdges(sqcol):squareEdges(sqcol+1));
        square(sqnum).neglum.L = -1;
        square(sqnum).neglum.M = -1;
        
        sqnum = sqnum+1;
    end
end


% Go through squares and test for responses.
tempResp = nan(nPres,1);
fn = fieldnames(square);

% Loop through each stimulus type
for st = 1:numel(fn)
    
    % Specify stimulus
    stim = cell2mat(fn(st));
    
    % Set up figure
    figure(fig); clf; hold on; grid on;
    fig = fig+1;
    
    % Loop through each square
    for sq = 1:numel(square)
        
        % Responses
        % **NOTE: I am currently taking the mean of the response at each point.
        % This may not reflect the properties of real nuerons, so consider
        % alternate strategies.  Using the sum of large areas will yield a
        % response at all locations, even those without an underlying subunit.
        
        % Send stimulus to REX
        for rep = 1:nPres
            tempResp(rep) = mean(mean(RFModel(square(sq).(stim).c,square(sq).(stim).r,nPts1D_MG_NC,DVAperStix,...
                square(sq).(stim).L,square(sq).(stim).M,S,addNoise)));
        end
        sqResp = mean(tempResp);
        
        if sum(sum(sqResp)) > threshold
            
            square(sq).(stim).test = 1;
            
            % Plot square boundaries
            upperLeft_rc = [square(sq).(stim).r(1,1) square(sq).(stim).c(1,1)] + [.25 -.25];
            lowerLeft_rc = [square(sq).(stim).r(end,1) square(sq).(stim).c(end,1)] + [-.25 -.25];
            upperRight_rc = [square(sq).(stim).r(1,end) square(sq).(stim).c(1,end)] + [.25 .25];
            lowerRight_rc = [square(sq).(stim).r(end,end) square(sq).(stim).c(end,end)] + [-.25 .25];
            
            square(sq).(stim).topline = [upperLeft_rc; upperRight_rc]';
            square(sq).(stim).rightline = [upperRight_rc; lowerRight_rc]';
            square(sq).(stim).bottomline = [lowerRight_rc; lowerLeft_rc]';
            square(sq).(stim).leftline = [lowerLeft_rc; upperLeft_rc]';
            
            square(sq).(stim).perimeter = cat(2,square(sq).(stim).topline,square(sq).(stim).rightline,...
                square(sq).(stim).bottomline,square(sq).(stim).leftline);
            
            plot(square(sq).(stim).perimeter(2,:),square(sq).(stim).perimeter(1,:),'k--')
            plot(square(sq).(stim).c,square(sq).(stim).r,'co')
            legend('Area spanned by squares that evoked a response')
            
        else
            
            square(sq).(stim).test = 0;
            
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
    
end

if stops == 1
    keyboard
end


%% Simulate the Procession through Rounds

for r = 1:nRnds
    
    % Round 1 is special because it is using the squares from the
    % pre-screening to guide point placement.  Subsequent rounds look for
    % boarders between points that evoke responses and points that do not
    % evoke responses.
    if r == 1
        
        rndNPts1D = initNPts;
        
        for sttype = 1:numel(fn)
            
            % Specify stimulus
            stim = cell2mat(fn(sttype));
            
            % Look only at points in squares that evoked a response.
            sqIdx = [];
            for t = 1:numel(square)
                sqIdx(t) = square(t).(stim).test;
            end
            sqIdx = find(sqIdx);
            
            for sq = 1:numel(sqIdx)
                
                % Use only corners of the squares as the initial grid.
                upperLeft = [square(sqIdx(sq)).(stim).r(1,1) square(sqIdx(sq)).(stim).c(1,1)];
                lowerLeft = [square(sqIdx(sq)).(stim).r(end,1) square(sqIdx(sq)).(stim).c(end,1)];
                upperRight = [square(sqIdx(sq)).(stim).r(1,end) square(sqIdx(sq)).(stim).c(1,end)];
                lowerRight = [square(sqIdx(sq)).(stim).r(end,end) square(sqIdx(sq)).(stim).c(end,end)];
                
                ptsInRnd_test_MG_C_r = cat(1,ptsInRnd_test_MG_C_r,...
                    [upperLeft(1); lowerLeft(1); upperRight(1); lowerRight(1)]);
                ptsInRnd_test_MG_C_c = cat(1,ptsInRnd_test_MG_C_c,...
                    [upperLeft(2); lowerLeft(2); upperRight(2); lowerRight(2)]);
                
            end
            
            % Housekeeping for subsequent rounds...
            
            % Set the response of all points outside of the squares that
            % evoked a response to 0.
            allSquares_MG_C_r = [];
            allSquares_MG_C_c = [];
            for t = 1:numel(sqIdx)
                allSquares_MG_C_r = cat(1,allSquares_MG_C_r,square(sqIdx(t)).(stim).r);
                allSquares_MG_C_c = cat(1,allSquares_MG_C_c,square(sqIdx(t)).(stim).c);
            end
            allSquares_MG_r = allSquares_MG_C_r + max(masterGrid1D_C)+1;
            allSquares_MG_c = allSquares_MG_C_c + max(masterGrid1D_C)+1;
            allSquares_I = sub2ind(size(responses.poslum),allSquares_MG_r(:),allSquares_MG_c(:));
            tempRespL = ones(size(responses.poslum));
            tempRespL(allSquares_I) = 0;
            responses.(stim)(logical(tempRespL)) = 0;
            
            ptsInRnd_MG_C_r = ptsInRnd_test_MG_C_r;
            ptsInRnd_MG_C_c = ptsInRnd_test_MG_C_c;
            
            % Indexing of points on larger grid (to keep track of responses)
            ptsInRnd_MG_I = sub2ind(size(masterGrid_C_c),...
                ptsInRnd_test_MG_C_r+max(masterGrid1D_C)+1,...
                ptsInRnd_test_MG_C_c+max(masterGrid1D_C)+1);
            
            % Send stimuli to rex
            for st = 1:numel(ptsInRnd_test_MG_C_c)
                tempResp = nan(nPres,1);
                stim_c = ptsInRnd_MG_C_c(st);
                stim_r = ptsInRnd_MG_C_r(st);
                for rep = 1:nPres
                    tempResp(rep) = RFModel(stim_c,stim_r,nPts1D_MG_NC,DVAperStix,...
                        square(sq).(stim).L,square(sq).(stim).M,S,addNoise);
                end
                
                responses.(stim)(ptsInRnd_MG_I(st)) = mean(tempResp);
            end
            
        end
        
        
    else %additional rounds
        
        % Establish a new grid
        rndNPts1D = (rndNPts1D-1)*2+1;
        ptsInRnd1D_MG_C = linspace(min(masterGrid1D_C),max(masterGrid1D_C),rndNPts1D);
        [ptsInRnd_MG_C_c ptsInRnd_MG_C_r] = meshgrid(ptsInRnd1D_MG_C,ptsInRnd1D_MG_C);
        ptsInRnd_MG_I = sub2ind(size(masterGrid_C_c),...
            ptsInRnd_MG_C_r+max(masterGrid1D_C)+1,...
            ptsInRnd_MG_C_c+max(masterGrid1D_C)+1);
        
        % Identify the points that have envoked a (suprathreshold) response
        ptsInRnd_gaveResp_L = sign(responses(ptsInRnd_MG_I)-threshold)==1;
        [ptsInRnd_gaveResp_r ptsInRnd_gaveResp_c] = ind2sub(size(ptsInRnd_gaveResp_L),find(ptsInRnd_gaveResp_L));
        
        % If no points (on the first round) evoked responses, subsample
        % everywhere.
        if nansum(nansum(ptsInRnd_gaveResp_L)) == 0
            
            % Test all points except previously tested points
            ptsInRnd_test_MG_C_c = ptsInRnd_MG_C_c(isnan(responses(ptsInRnd_MG_I)));
            ptsInRnd_test_MG_C_r = ptsInRnd_MG_C_r(isnan(responses(ptsInRnd_MG_I)));
            ptsInRnd_test_MG_I = ptsInRnd_MG_I(isnan(responses(ptsInRnd_MG_I)));
            
            
        else % Otherwise, use previous responses to establish new candidate
            % points to be tested.
            
            % For all points that responded, consider immediately neighboring
            % points (up, down, left, right) as candidates.
            crossring = [1 0; 0 1; -1 0; 0 -1];
            diagring = [1 1; -1 1; 1 -1; -1 -1];
            ptsInRnd_noTest_rc = [];
            ptsInRnd_test_rc = [];
            
            for t = 1:numel(ptsInRnd_gaveResp_c)
                
                %Points that fall between two like-signed points don't get
                %tested.  These points are assigned an arbitrarily high
                %response so that they will be counted as 'tested' in
                %subsequent rounds.
                tempTestPts = repmat([ptsInRnd_gaveResp_r(t) ptsInRnd_gaveResp_c(t)],4,1) + crossring;
                if t == 1
                    ptsInRnd_test_rc = tempTestPts;
                else
                    tempnotest = intersect(tempTestPts,ptsInRnd_test_rc,'rows');
                    if ~isempty(tempnotest)
                        if isempty(ptsInRnd_noTest_rc)
                            ptsInRnd_noTest_rc = tempnotest;
                        else
                            ptsInRnd_noTest_rc = cat(1,ptsInRnd_noTest_rc,tempnotest);
                        end
                    end
                    ptsInRnd_test_rc = setxor(tempTestPts,ptsInRnd_test_rc,'rows');
                end
            end
            
            
            % To all points that responded previously, add diagonal points
            % as potential candidates for next round
            for t = 1:numel(ptsInRnd_gaveResp_c)
                
                tempTestPts = repmat([ptsInRnd_gaveResp_r(t) ptsInRnd_gaveResp_c(t)],4,1) + diagring;
                ptsInRnd_test_rc = cat(1,tempTestPts,ptsInRnd_test_rc);
                
            end
            
            % If a diagonal point appears 4 times in the queue, it is
            % completely surrounded by tested points, and does not, itself,
            % need to be tested.
            if ~isempty(ptsInRnd_test_rc)
                [ptsInRnd_test_rc ia ib] = unique(ptsInRnd_test_rc,'rows');
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
            
            
            % Set up points to test (indexing)
            try
                ptsInRnd_test_I = sub2ind(size(ptsInRnd_MG_I),...
                    ptsInRnd_test_rc(:,1),ptsInRnd_test_rc(:,2));
            catch err
                disp('Stimulus larger than initial grid!')
                keyboard
            end
            ptsInRnd_test_MG_I = ptsInRnd_MG_I(ptsInRnd_test_I);
            ptsInRnd_test_MG_C_c = masterGrid_C_c(ptsInRnd_test_MG_I);
            ptsInRnd_test_MG_C_r = masterGrid_C_r(ptsInRnd_test_MG_I);
            
            % Assign arbitrily large responses to noTest points (to designate
            % them at points which will almost surely give a response, but we
            % don't know what that repsonse is.
            if any(ptsInRnd_noTest_rc)
                ptsInRnd_noTest_I = sub2ind(size(ptsInRnd_MG_I),...
                    ptsInRnd_noTest_rc(:,1),ptsInRnd_noTest_rc(:,2));
                ptsInRnd_noTest_MG_I = ptsInRnd_MG_I(ptsInRnd_noTest_I);
                responses(ptsInRnd_noTest_MG_I) = 1000;
            end
            
        end
        
        % Send to REX to be tested.
        for st = 1:numel(ptsInRnd_test_MG_C_c)
            tempResp = nan(nPres,1);
            stim_c = ptsInRnd_test_MG_C_c(st);
            stim_r = ptsInRnd_test_MG_C_r(st);
            for rep = 1:nPres
                tempResp(rep) = RFModel(stim_c,stim_r,nPts1D_MG_NC,DVAperStix,addNoise);
            end
            responses(ptsInRnd_test_MG_I(st)) = mean(tempResp);
        end
        
    end
    
    % Plot the "real" underlying RF overlayed by tested and untested
    % points.
    figure(r); clf; hold on; grid on;
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
            plot(square(sq).(stim).perimeter(2,:),square(sq).(stim).perimeter(1,:),'k--')
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
    
    legend(leg)
    contour(xN,yN,ZZ)
    
    if stops == 1 && r < nRnds
        keyboard
    end
    
end




function [response] = RFModel(testx, testy, nSizeGrid, DVAperStix, L, M, S, addnoise) % Model Neuron

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
mean = 0;
std = (nSizeGrid-1)/6;
Gaussian = height.*exp(-((testx-mean).^2 ./ (2*std.^2)) - ((testy-mean).^2 ./ (2*std.^2)));

% DIfferentiate responses for +lum and -lum stimuli
if L == 1 && M == 1
    if addnoise == 1
        response = poissrnd(OnOffUnit .* Gaussian + cellBaseline);
    else
        response = OnOffUnit .* Gaussian;
    end
elseif L == -1 && M == -1
    if addnoise == 1
        response = poissrnd(-OnOffUnit .* Gaussian + cellBaseline);
    else
        response = -OnOffUnit .* Gaussian;
    end
end

response(isnan(response)) = 0;







