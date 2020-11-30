function GridLMSubunitOnline_workingJPW()

% 2/10/13  New paradigm created.    JPW

global udpCom gl

disp('In GridLMSubunitOnline')

gl.CurrentQueue = [];
gl.rnd = [];
gl.nPres = 3;

% These things will be sent in by rex?
gl.nStixelGrid = 17;
gl.DVAperStix = .1;
gl.StimDur = .2;
gl.initNPts = 5;
gl.threshold = 5;



udpCom.sock = [];
udpCom.port = 6665;
udpCom.rexip = '192.168.1.120';
udpCom.plexip = '192.168.1.121';
udpCom.mac = '192.168.1.122';

C = GridLMSubunitCodes();% defined in a separate m file

p = InitPStruct(0, C.HDRCOMPLETECD);
s = InitPlex();
initGUI();
initTrialStruct();

% Some initializations
TS = TrialSpecCodes;
currentspiketally = [];
stimon_t = [];
fix_t = [];


% The main loop
socketOpen = 0;
while ~socketOpen
    [udpCom.sock, socketOpen] = pnetStart(udpCom.port);
end
plxServer = InitPlex();

bounceOut = false;
while ~bounceOut
    if CheckForESCKey() || dealWithMsgs(udpCom.sock)
        bounceOut = true;
    end
    
    [n, eventList] = PL_GetTS(s);
    if n 
        p = ProcessEventList(p, eventList);
    end

    % do your magic hither
    if (any(p.events == C.FPACQCD))
        fix_t = p.times(find(p.events == C.FPACQCD,1));
        abortflag = 0;
    end
    if (any(p.events == C.STIMONCD))
        stimon_t = p.times(find(p.events == C.STIMONCD,1));
        nspikes = sum(currentspiketally > 0); % Spike times are relative to fix_t
        p.lastprocessed_t = p.times(find(p.events == C.STIMONCD,1));
    end
    if (any(p.spikes{1}) && ~isempty(fix_t) && ~abortflag)
        currentspiketally = [currentspiketally; p.spikes{1}-fix_t]; %Spike times are relative to fix_t
        %PlotRaster(currentspiketally);
        p.spikes{1} = []; 
    end
    if (any(p.events == C.ABORTCD))
        abortflag = 1;
        currentspiketally = [];
        %ClearRaster;
        p.lastprocessed_t = p.times(find(p.events == C.ABORTCD,1,'last'));
    end
    if (any(p.events == C.FPOFFCD))
        stimoff_t = p.times(find(p.events == C.FPOFFCD,1));

%         trial.fpacq(tcounter) = fix_t;
%         trial.stimon(tcounter) = stimon_t;
%         trial.stimoff(tcounter) = p.times(find(p.events == C.FPOFFCD,1));    
%         trial.duration(tcounter) = stimoff_t-stimon_t;
%         trial.tspikes{tcounter} = currentspiketally + stimon_t;
%         trial.normtspikes{tcounter} = currentspiketally;
        
        % Spike Historgram
%         PlotSpikeHist
        
        % Clear out fields
        currentspiketally = [];
        stimon_t = [];
        fix_t = [];        
        p.lastprocessed_t = stimoff_t;
        
        % Organize and Plot data
%         organizeParStruct();
%         plat = organizePlatStruct();
%         plotGUI(plat);
        
        disp('Sending plexdone message...')
        pnet(udpCom.sock, 'write', 'plexdone>> >>');
        pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
        disp('plexdone message sent')
        
    end
    p = CleanUpEvents(p);
end

PL_Close(plxServer);

end


function PlotRaster(spiketally)
global trial

    if any(~isnan(trial.duration))
        xmax = max(trial.duration);
    else
        xmax = 3000;
    end
    
    figure(1);
    subplot(2,5,[1 2]); cla; hold on;
    axis([-500 xmax -.5 .5])
    plot([spiketally, spiketally]',[-.5; .5]*ones(1,length(spiketally)),'k-');
    title('Incoming Spikes')
    
end

function PlotSpikeHist
global trial tcounter

%     % Set up variables
%     bins = linspace(-500,max(trial.duration),75);
%     PSTH = zeros(1,length(bins));
%     circDelay = .05;
% 
%     % Generate spiking profile across all trials
%     for i = 1:tcounter
%         if ~isempty(trial.normtspikes{i})
%             PSTH = PSTH + histc(trial.normtspikes{i}',bins);
%         end
%     end
% 
%     % Find baseline spikerate
%     nsp =  sum(trial.tspikes{tcounter} > trial.fpacq(tcounter) + circDelay & trial.tspikes{tcounter} < trial.stimon(tcounter));
%     dt = trial.stimon(tcounter) - trial.fpacq(tcounter);
%     trial.baselinefr(tcounter) = nsp./(dt - circDelay);
% 
%     % Fit Gaussian
%     gauss = @(x,params)(params(2)+params(1)*exp((-(x-params(3)).^2)/params(4).^2));
%     bguess = mean(PSTH(bins<.25));
%     aguess = max(PSTH)-bguess;
%     muguess = bins(find(PSTH == max(PSTH),1,'first'));
%     sigmaguess = .1; % terrible. Need something better.
%     fittedparams = fminsearch(@(params)sum((gauss(bins,params)-PSTH).^2),[aguess,bguess,muguess,sigmaguess]);
%     offset = [fittedparams(3)-fittedparams(4) fittedparams(3)+fittedparams(4)];
% 
%     % Plot spiking profile
%     figure(1);
%     subplot(2,5,[6 7]); cla; hold on;
%     xlim([-500 max(trial.duration)])
%     bar(bins,PSTH);
%     plot(bins,gauss(bins,fittedparams),'b-.','LineWidth',3);
%     plot(offset,[fittedparams(1) fittedparams(1)],'k-','linewidth',3);
%     xlabel('Normalized time (ms)')
%     ylabel('Spike Count')
%     title('Spike Histogram')
% 
%     % Now counting up spikes in a window
%     for i = 1:tcounter
%         trial.nspikes(i) = sum(trial.normtspikes{i} > offset(1) & trial.normtspikes{i} < offset(2));
%         dt = offset(2) - offset(1);
%         trial.fr(i) = trial.nspikes(i)./dt;
%     end

end


function allDone = dealWithMsgs(socket)
    global udpCom
    allDone = false;
    msgSize = pnet(socket, 'readpacket', 250, 'noblock');
    if ~msgSize, return; end

    message = pnet(socket, 'read', msgSize, 'char');
    if ~isempty(message)
        [message,allDone] = parseMsg(message,socket);
        evalMsg(message);
    else
        allDone = true;
    end
end


% here we implement the logic/function calls associated with the requested
% variables within the message
function [message,doneFlag] = parseMsg(message, socket)
    global udpCom
    %disp('In parseMsg')
    doneFlag = false;


    matches = strtrim(regexp(message, ',', 'split'));
    if strncmp(matches{1}, 'sendToRex', 9)
        varToSend = matches{2}; % the second match contains the requested variable(s)
        if strfind(varToSend, 'stimColor')
            keyboard;
        elseif strfind(varToSend, 'some.other.variable')
            % do something else
    %     else
    %         eval(message)
        end
    elseif strncmp(message, 'return', 6) % not sure what catastrophe this prevents
        stk = dbstack();
        if ~strcmp(stk(end).name, mfilename)
            doneFlag = true;
        end
    end
end


function evalMsg(message)
    global udpCom gl %#ok<NUSED> eval needs access to gl
    %disp('In evalMsg')
    try
        eval(message);
    catch exception
        fprintf('Trouble with message: "%s"\n', message);
        disp(getReport(exception));
    end
end
     

function [stim] = getStimParams()
global gl

    disp('Selecting Next Stimulus...')
    
    % Large Square Pre-Screening
    if isempty(gl.CurrentQueue) && isempty(gl.rnd)
        
        % Set up grid for probing and tracking responses
        nPts1D_MG_NC = gl.nStixelGrid;
        gl.masterGrid1D_C = -(nPts1D_MG_NC-1)/2:(nPts1D_MG_NC-1)/2;
        [gl.masterGrid_C_c gl.masterGrid_C_r] = meshgrid(gl.masterGrid1D_C,fliplr(gl.masterGrid1D_C));
        gl.responses = nan(size(gl.masterGrid_C_c));
    
        % Set up squares
        [square] = SetUpSquares();
        
        % Set up queue ([L M],[x y])
        LM = [1 1; 1 -1; -1 -1; -1 1];
        LMxyIdx = fullfact([size(LM,1) numel(square)]);
        
        % Randomize order for presentation
        LMxyIdx_r = LMxyIdx(randperm(size(LMxyIdx,1)),:);        
        
        % Convert parameters that are not changing from single values to vectors
        S = zeros(size(LMxyIdx,1),1);
        param_nstixgrid = repmat(gl.nStixelGrid,size(LMxyIdx,1),1);
        param_dvaperstix = repmat(gl.DVAperStix,size(LMxyIdx,1),1);
        param_stimdur = repmat(gl.StimDur,size(LMxyIdx,1),1);
        
        % Assign values to CurrentQueue ([L,M,S,nstixgrid,dvaperstix,stimdur])

        %%%%% PLAN: Seperate [x y] coordinates from the other parameters.
        %%%%% This will allow us to send overy as many [x y] pairs as we
        %%%%% like.  gl.CurrentQueue is now divoid of xy pairs...
        
        gl.CurrentQueue = [LM(LMxyIdx_r(:,1),:) S param_nstixgrid param_dvaperstix param_stimdur];
        for n = 1:size(LMxyIdx_r,1)
            gl.x(n,:) = square(LMxyIdx_r(n,2)).c(:);
            gl.y(n,:) = square(LMxyIdx_r(n,2)).r(:);
        end
        
        
        % Lastly, set gl.rnd to 1 to begin looping through the rounds
        %gl.rnd = 1;
        
    
    elseif isempty(gl.CurrentQueue) && gl.epoch == 2
        
        
        % Look only at points in squares that evoked a response.
        sqIdx = find(cat(1,gl.square.test)==1);
        for sq = 1:numel(sqIdx)
            
            % Use only corners of the squares as the initial grid.
            upperLeft = [gl.square(sqIdx(sq)).r(1,1) gl.square(sqIdx(sq)).c(1,1)];
            lowerLeft = [gl.square(sqIdx(sq)).r(end,1) gl.square(sqIdx(sq)).c(end,1)];
            upperRight = [gl.square(sqIdx(sq)).r(1,end) gl.square(sqIdx(sq)).c(1,end)];
            lowerRight = [gl.square(sqIdx(sq)).r(end,end) gl.square(sqIdx(sq)).c(end,end)];
            
            ptsInRnd_test_MG_C_r = cat(1,ptsInRnd_test_MG_C_r,...
                [upperLeft(1); lowerLeft(1); upperRight(1); lowerRight(1)]);
            ptsInRnd_test_MG_C_c = cat(1,ptsInRnd_test_MG_C_c,...
                [upperLeft(2); lowerLeft(2); upperRight(2); lowerRight(2)]);
            
        end
        
        % Housekeeping for subsequent rounds...
        
        % Set the response of all points outside of the squares that
        % evoked a response to 0.
        allSquares_MG_C_r = cat(1,square(sqIdx).r);
        allSquares_MG_C_c = cat(1,square(sqIdx).c);
        allSquares_MG_r = abs(allSquares_MG_C_r - max(masterGrid1D_C)-1);
        allSquares_MG_c = allSquares_MG_C_c + max(masterGrid1D_C)+1;
        allSquares_I = sub2ind(size(responses),allSquares_MG_r(:),allSquares_MG_c(:));
        tempRespL = ones(size(responses));
        tempRespL(allSquares_I) = 0;
        responses(logical(tempRespL)) = 0;
               
        %ptsInRnd_MG_C_r = ptsInRnd_test_MG_C_r;
        %ptsInRnd_MG_C_c = ptsInRnd_test_MG_C_c;
        
        % Indexing of points on larger grid (to keep track of responses)
        ptsInRnd_MG_I = sub2ind(size(masterGrid_C_c),...
            abs(ptsInRnd_test_MG_C_r-max(masterGrid1D_C)-1),...
            ptsInRnd_test_MG_C_c+max(masterGrid1D_C)+1);
        
        % Set up CurrentQueue
        % ([L,M,S,nstixgrid,dvaperstix,stimdur])
        for rep = 1:gl.nPres
            
            %Randomize order
            randIdx = randperm(numel(ptsInRnd_MG_C_c));
            gl.CurrentQueue = cat(1,gl.CurrentQueue,[L M S pre_nstixgrid pres_dvaperstix pres_stimdur]);
            gl.x = cat(1,gl.x,ptsInRnd_test_MG_C_c(randIdx));
            gl.y = cat(1,gl.y,ptsInRnd_test_MG_C_r(randIdx));
            
        end
        
        
    end

    
    % Send top stimulus in the Queue to REX, then delete it from the queue.
    %stim = [L M S nStixelGrid dvaPerStixel stimDur];
    stim = gl.CurrentQueue(1,:);
    gl.CurrentQueue(1,:) = [];
    
    %%%% TO DO: Send a plexdone message, as well...
    disp('Sending stimulus parameters to REX...')
    
    
end

function x_coordinates = getXCoordinates()
global gl

    % Return the first row vector, then delete it (like the CurrentQueue)
    x_coordinates = gl.x(1,:);
    gl.x(1,:) = [];

end

function y_coordinates = getYCoordinates()
global gl

    % Return the first row vector, then delete it (like the CurrentQueue)
    y_coordinates = gl.y(1,:);
    gl.y(1,:) = [];

end


function out = TrialSpecCodes


end


function p = CleanUpEvents(p)

    if (p.processnowcode ~= 0)
        pncodeidx = find(p.events == p.processnowcode,1);
        p.lastprocessed_t = p.times(pncodeidx);
    end

    L = p.times <= p.lastprocessed_t;
    
    p.times(L) = [];
    p.events(L) = [];
    for i = 1:length(p.spikes)
        spiketimevect = p.spikes{i};
        p.spikes{i} = [spiketimevect(spiketimevect > p.lastprocessed_t)];
    end
    p.processnowflag = 0;
end


function initGUI()

    close all

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



