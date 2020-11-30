function [] = GLMS_StimSelection()
global gl

% This code is an attempt to divide GLMS into subfunctions for more
% efficient calculation.

% 3/4/13    Created.    JPW


%% First, set up the space.

% User-Defined Variables (some (all?) of these will be passed in)
gl.initNPts = 5;
threshold = 5;
nPres = 3;
DVAperStix = .1;
stops = 1;
nPts1D_MG_NC = 17;
addNoise = 1;
S = 0;


% Calculating total number of rounds
pts = gl.initNPts;
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
gl.masterGrid1D_C = -(nPts1D_MG_NC-1)/2:(nPts1D_MG_NC-1)/2;
[gl.masterGrid_C_c gl.masterGrid_C_r] = meshgrid(gl.masterGrid1D_C,fliplr(gl.masterGrid1D_C));
responses = nan(size(gl.masterGrid_C_c));
ptsInRnd_test_MG_C_c = [];
ptsInRnd_test_MG_C_r = [];

%Set Xlim and Ylim to values that can be called for all figures
xLim = [min(gl.masterGrid1D_C)*1.1 max(gl.masterGrid1D_C)*1.1];
yLim = [min(gl.masterGrid1D_C)*1.1 max(gl.masterGrid1D_C)*1.1];


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

%% Pre-Screening
if isempty(currentQueue) && isempty(rnd)
    
    [square] = SetUpSquares()
    
    

end




%% Proceed through rounds

    
    
if rnd == 1
    
    
end
    
    


%% Test stimuli

currentQueue = []



end


function [square] = SetUpSquares() % Round 1 (Large Squares Pre-Screening)
    global gl

    
    squareEdges = linspace(min(gl.masterGrid1D_C),max(gl.masterGrid1D_C),gl.initNPts)...
        +max(gl.masterGrid1D_C)+1;
    sqnum = 1;
    for sqrow = 1:numel(squareEdges)-1
        for sqcol = 1:numel(squareEdges)-1
        
            % Set up squares that test +luminance
            square(sqnum).r = gl.masterGrid_C_r(squareEdges(sqrow):squareEdges(sqrow+1),...
                squareEdges(sqrow):squareEdges(sqrow+1));
            square(sqnum).c = gl.masterGrid_C_c(squareEdges(sqcol):squareEdges(sqcol+1),...
                squareEdges(sqcol):squareEdges(sqcol+1));
        
            sqnum = sqnum+1;
        end
    end
end


function [response] = RFModel(testx, testy, DVAperStix, stimDur, L, M, S, nSizeGrid, addnoise) % Model Neuron

% stimDur is ignored here, but is not ignored in the real paradigm.

% nSizeGrid is not passed into the real slave program, but is necessary
% here only to create the RF.

% addnoise is for the model only.


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
    

end


