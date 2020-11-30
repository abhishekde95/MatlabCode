% PhaseNTsim

% 3D simulation revisited with better isoresponse measurements (equating
% responses for the 0° and 180° phase-shifted stimuli.

% 1) add fminsearch to find the contrast producing the criterion response.
% Assuming linear contrast-response function for now...

coneweights = [10 20 0; 0 0 0];
stimcc = [1 1 1;.1 -.1 0];
%stimcc = unifrnd(-1,1,2,3);

% whichmodel = 'NONLINEARENERGY';
whichmodel = 'ENERGY';
GAMUTEDGE = 100;
CRITFR = 20;

TIMEBINS = linspace(0,6*pi,100);
RFtemplate = sin(TIMEBINS).*normpdf(TIMEBINS,mean(TIMEBINS),2); 
RFtemplate(2,:) = cos(TIMEBINS).*normpdf(TIMEBINS,mean(TIMEBINS),2);

RFs = zeros(size(coneweights,2),length(TIMEBINS), 2, size(coneweights,1)); % cones x space x sin/cos x mechanism
for i = 1:size(coneweights,1) % mechanisms
    for j = 1:size(coneweights,2) % cones
        for k = 1:size(RFtemplate,1) % sin/cos
            RFs(j,:,k,i) = RFtemplate(k,:).*coneweights(i,j);
        end
    end
end

% Simulated experiment
% First, finding the contrast of the two components that produce the
% criterion firing rate.
componentcontrasts = [nan nan];

for whichcomponent = 1:2
    componentcontrasts(whichcomponent) = findContrast([stimcc(whichcomponent,:); 0 0 0], [0 0], TIMEBINS, RFs, whichmodel, CRITFR, GAMUTEDGE);
     
    %contrast = 1;
    %stim = MkStim(contrast, [stimcc(whichcomponent,:); 0 0 0],[0 0], TIMEBINS);
    %r = SimResp(stim, RFs, whichmodel);
   % componentcontrasts(whichcomponent) = 1/r;
   % stim = MkStim(1/r*contrasts, stimcc,[0 0], TIMEBINS);
end
if any(componentcontrasts == GAMUTEDGE)
    disp('One of the components does not drive a criterion response');
    keyboard
end
stimcc = diag(componentcontrasts)*stimcc; % Destructively modifying stimcc so that both components produce same response

% Now finding contrast of the compound grating in the 0 and pi phases that
% produce the criterion response.
phases = [0 pi];
startingcontrasts = [nan nan];
for phase = phases
    startingcontrasts(phases == phase) = findContrast(stimcc, [0 phase], TIMEBINS, RFs, whichmodel, CRITFR, GAMUTEDGE);

    %stim = MkStim(1, stimcc, [0 phase], TIMEBINS);
    %r = SimResp(stim,RFs, whichmodel);
    %startingcontrasts(phases == phase) = min(1/r,GAMUTEDGE); % Yuck
end

% figuring out which phases to test based on the relative sensitivity to
% the 0° and 180° phase-shifted stimuli.
phases = linspace(0,pi,21);
[x,y] = pol2cart(phases,1); % get a bunch of points on the unit circle
x = x.*startingcontrasts(1); y = y.*startingcontrasts(2); % stretch % This is wrong
phases = atan2(y,x); phases = phases*2; phases(phases<0) = 2*pi+phases(phases<0);

% Now the actual experiment
isorespcontrasts = zeros(size(phases));
for phase = phases
    isorespcontrasts(phases == phase) = findContrast(stimcc, [0 phase], TIMEBINS, RFs, whichmodel, CRITFR, GAMUTEDGE);

    
    %stim = MkStim(componentcontrasts, stimcc,[0 phase], TIMEBINS);
    %r = SimResp(stim,RFs, whichmodel);
    %isorespcontrasts(phases == phase) = min(1/r, GAMUTEDGE);
    
end
figure; subplot(2,1,1); 
LOOG = isorespcontrasts >= GAMUTEDGE;
polar(phases(~LOOG)/2,isorespcontrasts(~LOOG),'k-o'); hold on;
polar(phases(LOOG)/2,isorespcontrasts(LOOG),'ro');

subplot(2,1,2); hold on;
[x,y] = pol2cart(phases/2, isorespcontrasts);
plot(x(~LOOG)./max(x),y(~LOOG)./max(y),'k-o');
plot(x(LOOG)./max(x),y(LOOG)./max(y),'ro');
axis equal


% Nested functions
function c = findContrast(stimcc, phases, TIMEBINS, RFs, whichmodel, criterion, GAMUTEDGE)
   
    try
         x0 = 0;
         [c,fval,exitflag] = fminsearch(@nestedfun, x0); 
         if (c == x0 & exitflag == 1) | c > GAMUTEDGE
             c = GAMUTEDGE;
         end
     catch
         keyboard
     end
     
     function err = nestedfun(x0)
         err = (SimResp(MkStim(x0, stimcc, phases, TIMEBINS), RFs, whichmodel)-criterion).^2;
     end
end

function stim = MkStim(contrast, stimcc, phases, x)
     stim = zeros(size(stimcc,2),length(x),size(stimcc,1));
     stim(:,:,1) = contrast*stimcc(1,:)'*sin(x+phases(1));
     stim(:,:,2) = contrast*stimcc(2,:)'*sin(x+phases(2));
     stim = sum(stim,3);
end

function response = SimResp(stim, RFs, whichmodel)
    NAKARUSHTON = true;
    tmp =[];
    for k = 1:size(RFs,3) % sin/cos
        for l = 1:size(RFs,4) % mechanisms
            tmp = [tmp, sum(sum(RFs(:,:,k,l).*stim))];
        end
    end
    if strcmp(whichmodel, 'ENERGY')
        response = sqrt(sum(tmp.^2));
    elseif strcmp(whichmodel, 'NONLINEARENERGY')
        response = sqrt((abs(tmp(1))+abs(tmp(2)))^2+(abs(tmp(3))+abs(tmp(4)))^2);
    else
        error('Unknown model type');
    end
    if NAKARUSHTON
        response = 200*response/(response+50);
    end
end
