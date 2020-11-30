%% This code is a uniform sampling sceme in the LM plane.
% Created   10/31/11    JPW
% Revised   12/7/11     JPW    
%       -Added Polar Grid
% Revised   12/19/11    JPW    
%       -Added S-Cone Plateau
% Revised   2/19/12     JPW
%       -Changed to logarithmic spacing
% Revised   7/18/12     JPW
%       -Adding plateau for each variable

clear all
close all

%% User-Defined Variables

% Number of samples at each stim
nTestFr = 2;

% Number of rounds
nRounds = 1;

% Sent into getStimParams.m
grate.theta = pi;
grate.sigx = .5;
grate.sigy = .1;
grate.DR = 3;
grate.nstd = 3;
grate.nexp = 6;

% Plateau Variables
plat.Scc = [-1 0 1]'; %Scc amounts
plat.theta = [0 pi]'; %Rotations?
plat.sigX = 1; %Multiples?
plat.sigY = 1; %Multiples?
plat.driftRate = 1; %Multiples?
plat.nstd = 1; %Multiples? 
plat.nexp = 1; %Multiples?                                                                                                                                                                                                                                                                                                                                        

% General variables
testedptMat = [];
rnd = [];
currentQueue = [];
testedptMat = [];
gridspace = [];
sisterAngs = [0 pi/2 pi 3*pi/2]';

% Init par struct
par.polstim = [];
par.cartstim = [];
par.frs = [];
par.meanfr = [];

% Init trial struct
trial.polstim = [];
trial.fr = [];
trial.parIdx = [];
tcounter = 0;


%% Present Stim

for i = 1:nRounds

    if isempty(rnd)
        rnd = 1;
    else
        rnd = rnd + 1;
    end

    % Construct Polar Grid
    if rnd == 1
        thetaspace = pi/4;
        rhospace = .5;
    elseif rnd>1 && mod(rnd,2)~=0
        thetaspace = thetaspace * .5;
    elseif rnd>1 && mod(rnd,2)==0
        rhospace = rhospace * .5;
    end
        
    theta = 0:thetaspace:pi/2;
    rho = 0:rhospace:1;
    polTheta = repmat(theta,numel(rho),1);
    polRho = repmat(rho,numel(theta),1)';
    
    % Polar Grid Organization
    allptMat = [polTheta(:) polRho(:)];
    quad1 = allptMat((allptMat(:,1) > 0) & (allptMat(:,2) > 0),:);
    unqueuedptMat = quad1(~ismember(quad1,testedptMat,'rows'),:);
    
    for n = 1:size(unqueuedptMat,1)
        
        choice = randi(size(unqueuedptMat,1));
        choiceTheta = unqueuedptMat(choice,1);
        choiceRho = unqueuedptMat(choice,2);

        sisterStim = [repmat(choiceTheta,4,1) + sisterAngs repmat(choiceRho,4,1)];        

        % For plateaus
        PlatIdx = fullfact([size(sisterStim,1) numel(plat.Scc) numel(plat.theta) numel(plat.sigX)...
            numel(plat.sigY) numel(plat.driftRate) numel(plat.nstd) numel(plat.nexp)]);
        nplat = size(PlatIdx,1);
        
        out(:,1) = plat.Scc(PlatIdx(:,2));
        out(:,2) = repmat(grate.theta,nplat,1) + plat.theta(PlatIdx(:,3));
        out(:,3) = repmat(grate.sigx,nplat,1) .* plat.sigX(PlatIdx(:,4));
        out(:,4) = repmat(grate.sigy,nplat,1) .* plat.sigY(PlatIdx(:,5));
        out(:,5) = repmat(grate.DR,nplat,1) .* plat.driftRate(PlatIdx(:,6));
        out(:,6) = repmat(grate.nstd,nplat,1) .* plat.nexp(PlatIdx(:,7));
        out(:,7) = repmat(grate.nexp,nplat,1) .* plat.nstd(PlatIdx(:,8));

        stimFam = [repmat(sisterStim,nplat/size(sisterStim,1),1) out];
        
        % WORKING: a reference into stim type (does not yet work for plateaus!!)
        par.polLM = cat(1,par.polstim,stimFam(:,1:2));
        [L,M] = pol2cart(stimFam(:,1),stimFam(:,2));
        L = rndofferr(L,3);
        M = rndofferr(M,3);
        par.cartLM = cat(1,par.cartstim,[L M]);
        stimFam(:,1:2) = [L M];
        
        tempParIdx = [];
        
        % Randomize Order of Rows from StimFam, and repeat nTestFr times
        for q = 1:nTestFr
            order = randperm(size(stimFam,1))';
            tempParIdx = cat(1,tempParIdx,order);
            currentQueue = cat(1,currentQueue,stimFam(order,:));
        end
        
        trial.parIdx = cat(1,trial.parIdx,max([0; trial.parIdx]) + tempParIdx);
        unqueuedptMat(choice,:) = [];
    end
    testedptMat = unique(quad1,'rows');
    
%     % Plot Demonstration
%     reallyAllPts = [quad1; [quad1(:,1)+pi/2 quad1(:,2)]; [quad1(:,1)+pi quad1(:,2)]; [quad1(:,1)+3*pi/2 quad1(:,2)]];
%     reallyAllPts2 = cat(2,repmat(reallyAllPts,numel(splat),1),sort(repmat(splat(:),numel(reallyAllPts)/2,1)));
%     [x,y,z] = pol2cart(reallyAllPts2(:,1),reallyAllPts2(:,2),reallyAllPts2(:,3));
%     [a,b,c] = pol2cart(currentQueue(:,1),currentQueue(:,2),currentQueue(:,3));
%     for n = 1:size(currentQueue,1)
%         figure(1); clf; hold on, grid on; axis([-1.1 1.1 -1.1 1.1 -1.1 1.1]);
%         set(gca,'view',[-49 66]);
%         xlabel('L + M'); ylabel('L - M'); zlabel('S'); title('Polar Grid Stimulus Selection');
%         plot3(x,y,z,'kx')
%         plot3(a,b,c,'*b')
%         h = plot3(a(n),b(n),c(n),'ro');
%         set(h,{'LineWidth'},{1.5}')
%         pause(.25)
%     end
    
end
    
    
    