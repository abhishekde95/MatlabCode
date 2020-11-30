% LMTF population analyses
%
% Section 1) Summary statistics for fixational saccades
%
% Section 1.1) relationship between microsaccades and correct and
% incorrect (as a function of  TF and LM). Pretty rudamentary - this was
% hacked together for the grant. Requires that section 1 be run first.
%
% Section 1.2) Percent correct as a function of saccade direction. 
% Requires that section 1 be run first.
%
% Section 2) Looking at the residuals as a function of retinal position.
% Assumes that the variable "A" is in the workspace (structure of model
% parameters)
%
% Section 3) Comparing detection thresholds between upper and lower visual
% fields.
%
% Section 4) Cross-validated model comparisons.
%
% **** Section 4.1) Creating the LMTF model structure for Isosamp.d.
% The .mat file gets stored in MatlabCode/Zack/IsoSamp/private/data/LMTF.mat
% Trying also to get rid of the "domain" parameter.
%
% Section 4.2) Analysis of residuals from the model fit in the previous
% section (section 4.1).
%
% Section 4.3) Comparing two model fits to try to find stimuli for which
% the detection threshold is maximally different.
%
% Section 5) F-tests to compare the unconstrained model to the 11+2n model.
%
% Section 6) Just visually comparing TCSFs from humans and monkeys
%
% Section 7) Comparing sensitivity for different stimlus durations (looking
% at psychophysical temporal integration time).
%
% Section 8) Making a plot of screen locations tested broken down by
% observer (pie charts?)
%
% Section 9) Comparing human button box data to eye tracker data.
% 
% Section 10) Slices through a LMTF surface fit at a single spatial location + data.

%%
% Section 1) Summary statistics of fixational saccades. Note: getSacData.m
% will probably have to be changed to accept an optional veolcity threshold
% of saccades as measured with the eye tracker since it's (presumably)
% noisier than the eye coil.

SUBJECTINITIAL = 'U';
if SUBJECTINITIAL == 'G' || SUBJECTINITIAL == 'E'
    [filenames,spikecds] = fnamesFromTxt('LMTF','subjID',{SUBJECTINITIAL}, 'notes',{'eye tracker', 'eyetracker'});
    THRESHOLD = 15; % optical tracker
else
    [filenames,spikecds] = fnamesFromTxt('LMTF','subjID',{SUBJECTINITIAL}, 'monID',{'ProPixx'},'Nonstandard',{'0'});
    THRESHOLD = 8; % eye coil in rig 2
end

amplitudes = [];
directions = [];
peakv = [];
pathlengths = [];
durations = [];
trialparams = [];
sacduringstim = [];
BINWIDTH = .01;
offset = -0.5; % offset in sec from stim on (added, so make it negative)
timebins = [offset:BINWIDTH:.767];
PSTH = zeros(2,length(timebins));  % histogram of saccade initiations, not spikes

totalnumtrials = [0 0];
for a = 1:size(filenames,1)
    stro = nex2stro(findfile(char(filenames{a})));
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        continue;
    end
    ntrials = size(stro.trial,1);
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
    fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));

    tf = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
    lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'lcc'));
    mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'mcc'));
    correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    sacstats = getSacData(stro,Inf,THRESHOLD);
    close;
    
    for i = 1:ntrials
        st = sacstats.starttimes{i};
        Lsac = (sacstats.starttimes{i} > stimon_t(i)+offset) & (sacstats.endtimes{i} < fpoff_t(i));
        if any(Lsac)
            if any(sacstats.amplitudes{i}(Lsac) > 2)
                keyboard
            end
            amplitudes = [amplitudes; sacstats.amplitudes{i}(Lsac)];
            directions = [directions; sacstats.directions{i}(Lsac)];
            peakv = [peakv; sacstats.peakv{i}(Lsac)];
            pathlengths = [pathlengths; sacstats.pathlengths{i}(Lsac)];
            durations = [durations; sacstats.durations{i}(Lsac)];
            trialparams = [trialparams; repmat(correct(i),sum(Lsac),1) lcc(sacstats.trialnums{i}(Lsac)) mcc(sacstats.trialnums{i}(Lsac)) tf(sacstats.trialnums{i}(Lsac))];
            sactimes = [];
            for j = find(Lsac')
                tmp = sacstats.starttimes{i}(j)-stimon_t(i);
                sactimes = [sactimes; tmp((tmp > timebins(1)-BINWIDTH/2) & (tmp < timebins(end)+BINWIDTH/2))];
            end
            sacduringstim = [sacduringstim; sacstats.starttimes{i}(Lsac) > stimon_t(i)];
            [n,x] = hist(sactimes, timebins);
            PSTH(correct(i)+1,:) = PSTH(correct(i)+1,:) + n;
        end
        if length(sacduringstim) ~= length(peakv)
            keyboard
        end
        totalnumtrials(correct(i)+1) = totalnumtrials(correct(i)+1) + 1;
    end
end

figure;
subplot(3,2,1);
hist(amplitudes,50);
xlabel('amplitude'); ylabel('count');

subplot(3,2,2);
[rho, theta]= hist(directions,20);
polar([theta theta(1)],[rho rho(1)],'k-');

subplot(3,2,3);
[n,x] = hist2([amplitudes,peakv],[50 40]);
imagesc(flipud(n)); colormap(gray); axis xy;
m = (x{1}(end)-x{1}(1))./diff(get(gca,'XLim')+[.5 -.5]);
b = x{1}(1)-m;  % assuming first bin is at '1'
set(gca,'XTick',([0 .25 .5 .75 1]-b)./m,'XTickLabel',[0 .25 .5 .75 1]);
m = (x{2}(end)-x{2}(1))./diff(get(gca,'YLim')+[.5 -.5]);
b = x{2}(1)-m;  % assuming first bin is at '1'
bins = 25:50:150;
set(gca,'YTick',(bins-b)./m,'YTickLabel',bins);
xlabel('amplitude (deg)'); ylabel('peak vel. (deg/sec)');

subplot(3,2,4);
[n,x] = hist2([amplitudes,pathlengths],50);
imagesc(flipud(n)); colormap(gray); axis xy;
m = (x{1}(end)-x{1}(1))./diff(get(gca,'XLim')+[.5 -.5]);
b = x{1}(1)-m;  % assuming first bin is at '1'
set(gca,'XTick',([0 .25 .5 .75 1]-b)./m,'XTickLabel',[0 .25 .5 .75 1]);

m = (x{2}(end)-x{2}(1))./diff(get(gca,'YLim')+[.5 -.5]);
b = x{2}(1)-m;  % assuming first bin is at '1'
set(gca,'YTick',([0:.5:1.5]-b)./m,'YTickLabel',[0:.5:1.5]);
xlabel('amplitude (deg)'); ylabel('traj. length (deg)');

subplot(3,2,5); hold on;
plot(timebins, sum(PSTH)./(BINWIDTH*sum(totalnumtrials)),'k-');
plot([0 0],[0 3],'b:');
set(gca,'Xlim',[min(timebins) max(timebins)]);
xlabel('time wrt stimulus onset (s)'); ylabel('saccades/sec');

subplot(3,2,6); hold on;
plot(timebins, PSTH(1,:)./(BINWIDTH*totalnumtrials(1)),'k-');
plot(timebins, PSTH(2,:)./(BINWIDTH*totalnumtrials(2)),'b-');
legend('inc','cor');
ylims = get(gca,'Ylim');
plot([0 0],[0 ylims(2)],'b:');
set(gca,'Xlim',[min(timebins) max(timebins)]);
xlabel('time wrt stimulus onset (s)'); ylabel('saccades/sec');

% number of saccades, number of trials
L = timebins >= 0;
sum(PSTH(:,L),2)'
totalnumtrials
sum(PSTH(:,L),2)'./totalnumtrials
% For Apollo, 1 to 2% of trials have a microsaccade

%%
% Section 1.1
% Looking at early microsaccades and their relationhips with correct and
% incorrect (and TF and LM)

varnames = {'Up','Correct','L','M','TF'};
L = sacduringstim == 1; % Only looking at trials with microsaccades during the stimulus presentation (could be up or down)
Lup = directions > 0 & directions < pi;
%prettycorr([Lup(L) trialparams(L,:)],varnames);
mean(amplitudes(Lup)./durations(Lup)) % Speed of upward saccades
mean(peakv(Lup))

% Contingency table:
%     Up Down
%     -------
% Cor |
% Inc |
%[h,p] = fisherExact(sum(Lup(L)&trialparams(L,1)), sum(~Lup(L)&trialparams(L,1)),...
%                    sum(Lup(L)&~trialparams(L,1)), sum(~Lup(L)&~trialparams(L,1)))
% Here's the contingency table
[sum(Lup(L)&trialparams(L,1)), sum(~Lup(L)&trialparams(L,1));...
 sum(Lup(L)&~trialparams(L,1)), sum(~Lup(L)&~trialparams(L,1))]
% Looks like there are lots of "saccade down & incorrect" trials
[orat,~,~,p] = oddsratio(sum(Lup(L)&trialparams(L,1)), sum(~Lup(L)&trialparams(L,1)),...
                    sum(Lup(L)&~trialparams(L,1)), sum(~Lup(L)&~trialparams(L,1)),0.05)
% A little dodgy because each count is a microsaccade, not a trial.
                
% --------------------------------------------------
% L/M ratio
% --------------------------------------------------
% Only looking at LUM trials with microsaccades during the stimulus presentation (could be up or down)
LTF = trialparams(:,4) >= 5;  % High speed
Llum = trialparams(:,2) > 0 & trialparams(:,3) > 0;
L = sacduringstim == 1 & Llum & LTF;
    
%[h,p] = fisherExact(sum(Lup(L)&trialparams(L,1)), sum(~Lup(L)&trialparams(L,1)),...
%                    sum(Lup(L)&~trialparams(L,1)), sum(~Lup(L)&~trialparams(L,1)))

% Here's the contingency table
[sum(Lup(L)&trialparams(L,1)), sum(~Lup(L)&trialparams(L,1));...
 sum(Lup(L)&~trialparams(L,1)), sum(~Lup(L)&~trialparams(L,1))]
[orat,~,~,p] = oddsratio(sum(Lup(L)&trialparams(L,1)), sum(~Lup(L)&trialparams(L,1)),...
                    sum(Lup(L)&~trialparams(L,1)), sum(~Lup(L)&~trialparams(L,1)),0.05)
% An odds ratio >1 means up/correct occurs more often than down/correct
% --------------------------------------------------
% Only looking at RG trials with microsaccades during the stimulus presentation (could be up or down)
LTF = trialparams(:,4) >= 5;  % High speed
LRG = trialparams(:,2) < 0 & trialparams(:,3) > 0;
L = sacduringstim == 1 & LRG & LTF;
    
%[h,p] = fisherExact(sum(Lup(L)&trialparams(L,1)), sum(~Lup(L)&trialparams(L,1)),...
%                    sum(Lup(L)&~trialparams(L,1)), sum(~Lup(L)&~trialparams(L,1)))
% Here's the continegncy table
[sum(Lup(L)&trialparams(L,1)), sum(~Lup(L)&trialparams(L,1));...
 sum(Lup(L)&~trialparams(L,1)), sum(~Lup(L)&~trialparams(L,1))]

[orat,~,~,p] = oddsratio(sum(Lup(L)&trialparams(L,1)), sum(~Lup(L)&trialparams(L,1)),...
    sum(Lup(L)&~trialparams(L,1)), sum(~Lup(L)&~trialparams(L,1)),0.05)

% Sedna:
% At high TFs, down saccades precede incorrect choices. Speeding the
% stimulus on the retina.

% Apollo:
% Apollo's down saccades tend to precede incorrect choices. He doesn't make
% enough saccades to break this down by color/TF. He makes very few
% microsaccades.
% Freya shows nothing consistent.

%%
% Section 1.2
% Polar plot of correct probability as a function of saccade direction
% Probably not great that I'm not taking contrast into account
%
nbins=8;
colordir = atan2(trialparams(:,3),trialparams(:,2));
Lcol = true(length(colordir),1); % All color directions
%Lcol = Lcol & amplitudes > .3; % Adding an amplitude constraint
%Lcol = colordir>pi/6 & colordir<2*pi/6; % Lum
%Lcol = colordir> 4*pi/6 & colordir< 5*pi/6; % RG

tfbound = [5 25];
LTF = trialparams(:,4)>=tfbound(1) & trialparams(:,4)<=tfbound(2);
bins = linspace(0,2*pi,nbins+1);
bins(end) = [];
[t_all,r_all] = rose(directions(LTF&Lcol),bins);
counts_all = r_all([3:4:end]);
[t_cor,r_cor] = rose(directions(trialparams(:,1)==1&LTF&Lcol),bins);
counts_cor = r_cor([3:4:end]);
[t_sac,r_sac] = rose(directions(logical(sacduringstim)&LTF&Lcol),bins);
counts_sac = r_sac([3:4:end]);
[t_saccor,r_saccor] = rose(directions(logical(sacduringstim)&trialparams(:,1)==1&LTF&Lcol),bins);
counts_saccor = r_saccor([3:4:end]);

figure;
% All saccades
subplot(2,3,1); polar([bins bins(1)],[counts_all counts_all(1)]);
subplot(2,3,2); polar([bins bins(1)],[counts_cor counts_cor(1)]);
subplot(2,3,3); polar([bins bins(1)],[counts_cor counts_cor(1)]./[counts_all counts_all(1)],'k.-');

% Only saccades occuring during stimulus period
subplot(2,3,4); polar([bins bins(1)],[counts_sac counts_sac(1)]);
subplot(2,3,5); polar([bins bins(1)],[counts_saccor counts_saccor(1)]);
subplot(2,3,6);
polar([bins bins(1)],[counts_saccor counts_saccor(1)]./[counts_sac counts_sac(1)],'r.-');
hold on;
polar([bins bins(1)],[counts_cor counts_cor(1)]./[counts_all counts_all(1)],'k.--');

%%
% Section 2
% Looking at the residuals as a function of retinal position for the model
% in which the 10 filter-shape parameters are fixed and the other three
% parameters are allowed to vary unconstrained across positions.
% This script requires the LMTF "module" that is used by IsoSamp and
% contains all the model parameters.
% Black symbols are used for the FIXTHETA color direction and red symbols
% are used for the orthogonal direction.

FIXTHETA = 0; % <0 use theta from model fit, otherwise set to angle of interest (e.g. pi/4)
WEDGEWIDTH = 20; % degrees in the LM plane
PLOTSCALEFACTOR = 7;
FIGWIDTHINPIX = 900;
AXWIDTH = 70;
AXHEIGHT = 40;
isosamppath = which('IsoSampOnline');
isosamppath(find(isosamppath==filesep,1,'last'):end) = [];
load ([isosamppath,filesep,'private',filesep','data',filesep,'LMTF.mat']);

eccs = A.eccs;
models = A.legacy.secondroundmodels; % Shape parameters fixed, xi_lum, xi_rg, theta vary over retinal location
figure('Position',[500 200 FIGWIDTHINPIX FIGWIDTHINPIX]);
for i = 1:size(eccs,1)
   model = models(:,i);
   data = A.raw{i};
   if (data(1,5) ~= eccs(i,1) || data(1,6) ~= eccs(i,2))
       error ('Mismatch in retinal location between raw and ecc fields.');
   end
   xaxpos = PLOTSCALEFACTOR*eccs(i,2)+FIGWIDTHINPIX/2;
   axes('units','pixels','Position',[PLOTSCALEFACTOR*eccs(i,1) xaxpos AXWIDTH AXHEIGHT]); hold on;
   plot([1 25],[0 0],'k:');
   for j = 1:2 % LUM/RG
       if (j == 1)
           centertheta = model(end);
           if (FIXTHETA >= 0)
               centertheta = FIXTHETA;
           end
       else % The orthogonal direction
           centertheta = model(end)+pi/2;
           if (FIXTHETA >= 0)
               centertheta = FIXTHETA-pi/2;
           end
       end
       [th,r] = cart2pol(data(:,1),data(:,2));
     %  err = mod(th-centertheta,pi);
       err = mod(atan2(sin(th-centertheta),cos(th-centertheta)),pi);
       Langle = abs(err) < WEDGEWIDTH/2*(pi/180)
       if any(Langle)
           % Using the full 13-parameter model, not just a 1-D slice.
           % But taking a narrow wedge of data
           pred = LMTF_thresh_from_model(data(Langle,1:4),model);
           r = sqrt(data(Langle,1).^2+data(Langle,2).^2);
           resid = log10(r)-log10(pred);
           h = plot(data(Langle,3),resid,'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',2);
           if (j == 2)
               set(h,'MarkerFaceColor','red','MarkerEdgeColor','red')
           end
       end
   end
   set(gca,'YTick',[-.2 0 .2],'XLim',[1 25],'XScale','log');
end
disp('negative residuals mean that the model underestimates sensitivity');
disp('positive residuals mean that the model overestimates sensitivity');

%%
% Section 3
% Comparing detection thresholds between the upper and lower visual fields.
% Is it true the luminance sensitivity is slightly higher in the lower
% visual fields but that color sensitivity is the same in upper and lower
% VFs?
PLOTFUNNEL = false;
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
% filenames = fetch(conn, 'SELECT fileID FROM LMTF WHERE subjID=''P'' AND notes IN(''eye tracker'', ''eyetracker'')');
filenames = fetch(conn, 'CALL propixxStandardFilenames(''A'')');

[data, check] = getLMTFrawdata([], filenames, nexfilepath, []); % In Emily\Data Plotting\getLMTFrawdata.m
% Columns of "data"
% 1) L cone contrast
% 2) M cone contrast
% 3) Temporal frequency
% 4) Loog (1 = threshold was not measured; threshold point is placed on the edge of the display gamut)
% 5) X position of stimulus (1/10 deg)
% 6) Y position of stimulus (1/10 deg)
% 7) File identification number

uniqueXYs = unique(data(:,[5,6]),'rows');
% Finding and labeling uniqueXYs that do have matched upper and lower visual field
% pairs
L = false(size(uniqueXYs,1),1);
for i = 1:size(uniqueXYs,1)
    if (uniqueXYs(i,2) < 0)
       if any(uniqueXYs(i,1) == uniqueXYs(:,1) & abs(uniqueXYs(i,2)) == uniqueXYs(:,2))
           L(i) = 1;
       end
    end 
end
uniqueXYs = uniqueXYs(L,:); % Removing stimulus locations that aren't members of upper/lower VF pairs.
% By convention, all of the Y positions are negative at this point.

binedges = [1,2,4,8,16]; % [,)
binedges = [1,3,10];
nbins = length(binedges)-1;
Loog = logical(data(:,4));
% Settiong up for ellipse fits
pts = linspace(0, 2*pi, 100);
opt = optimoptions(@fmincon, 'Display', 'off', 'Algorithm','active-set');
for i = 1:size(uniqueXYs,1)
    L_lower = data(:,5) == uniqueXYs(i,1) & data(:,6) == uniqueXYs(i,2);
    L_upper = data(:,5) == uniqueXYs(i,1) & data(:,6) == -1*uniqueXYs(i,2);
    if (PLOTFUNNEL)
        figure; axes; hold on; set(gcf,'Name',num2str(uniqueXYs(i,:)));
        % First plotting lower visual field thresholds
        h = plot3(data(L_lower&~Loog,1),data(L_lower&~Loog,2),data(L_lower&~Loog,3),'kv','MarkerFaceColor','black');
        plot3(-data(L_lower&~Loog,1),-data(L_lower&~Loog,2),data(L_lower&~Loog,3),'kv','MarkerFaceColor','black');
        plot3(data(L_lower&Loog,1),data(L_lower&Loog,2),data(L_lower&Loog,3),'rv','MarkerFaceColor','black');
        plot3(-data(L_lower&Loog,1),-data(L_lower&Loog,2),data(L_lower&Loog,3),'rv','MarkerFaceColor','black');
        % Second plotting upper visual field thresholds
        plot3(data(L_upper&~Loog,1),data(L_upper&~Loog,2),data(L_upper&~Loog,3),'k^','MarkerFaceColor','blue');
        plot3(-data(L_upper&~Loog,1),-data(L_upper&~Loog,2),data(L_upper&~Loog,3),'k^','MarkerFaceColor','blue');
        plot3(data(L_upper&Loog,1),data(L_upper&Loog,2),data(L_upper&Loog,3),'r^','MarkerFaceColor','blue');
        plot3(-data(L_upper&Loog,1),-data(L_upper&Loog,2),data(L_upper&Loog,3),'r^','MarkerFaceColor','blue');
        set(gca,'Zscale','log');
    end
    figure; set(gcf,'Name',num2str(uniqueXYs(i,:)));
    for j = 1:nbins
        subplot(ceil(sqrt(nbins)),ceil(sqrt(nbins)),j); hold on;
        L_TF = data(:,3) >= binedges(j) & data(:,3) < binedges(j+1);
        maxcc = 0; % For equating axis limits
        for k = 1:2  % Not plotting OOG points for slice plots right now
            if (k == 1) % lower
                L = [data(L_TF&L_lower&~Loog,1); -data(L_TF&L_lower&~Loog,1)];
                M = [data(L_TF&L_lower&~Loog,2); -data(L_TF&L_lower&~Loog,2)];
                plot(L,M,'kv','MarkerFaceColor','black');
                color = 'k-';
            else
                L = [data(L_TF&L_upper&~Loog,1); -data(L_TF&L_upper&~Loog,1)];
                M = [data(L_TF&L_upper&~Loog,2); -data(L_TF&L_upper&~Loog,2)];
                plot(L,M,'b^','MarkerFaceColor','blue');
                color = 'b-';
            end
            maxcc = max([maxcc; abs([L;M])]);
            % Ellipse fitting
            if (length(L) >= 4) % constraint on minimal # of data points for ellipse fitting
                [v,d] = eig(cov([L M]));
                fitErrFn = @(params)ellipsefiterr(params,[L M],zeros(size(L))); % Assuming no OOGs
                alpha_b = fmincon(fitErrFn,[sqrt(d(2,2)), sqrt(d(1,1)),atan2(v(2,2), v(1,2))],[],[],[],[],[.01 .01 0],[10 10 2*pi], [], opt);
                theta = pts+alpha_b(3); rho = (alpha_b(1)*alpha_b(2))./sqrt((alpha_b(2)*cos(pts)).^2+((alpha_b(1)*sin(pts)).^2));
                [x,y] = pol2cart(theta,rho);
                plot(x,y,color);
            end
        end
        axis square; set(gca,'Xlim',[-1 1]*1.1*maxcc);  set(gca,'Ylim',[-1 1]*1.1*maxcc); 
        title([num2str(binedges(j)),' Hz to ',num2str(binedges(j+1)),' Hz'])
    end
end
% Black, downward pointing triangles are data from the lower hemifield
% Blue, upward pointing triangles are data from the upper hemifield

%%
% Section 4)
% Calling LMTF_generate_module_data.m with a data matrix as a step
% towards making a wrapper that will do cross validated model comparisons
% and create the model structures that are used in IsoSamp experiments.

% Currently assumes that there are 6 model to be tested:
% 1) firstroundmodels: Every retinal location is fitted independently
% 2) mode0: 11 parameters are fixed across retinal positions. LUM and RG gains only are allow to vary.
% 3) mode1: 10 parameters are fixed across retinal positions. LUM and RG gains and theta are allow to vary.
% 4) mode2: The rampy trough (LUM gain and RG gain are linear functions of r and phi^2)
% 5) mode3: Tilted rampy trough. (LUM gain is a linear function of r, phi, and phi^2)
% 6) mode4: Double tilted rampy trough (LUM and RG gains are linear
% functions of r, phi, and phi^2)
%
% Known error: For some datasets, there are too few thresholds at
% individual screen locations to fit all the models. Sometimes those points
% get selected for the leave-one-out. This causes a crash because some
% models don't make predictions at that location.
SAVEFILES = true;
niter = 1000;

conn = database('Nex_Paradigm_Sort','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
%filenames = fetch(conn, 'SELECT fileID FROM LMTF WHERE subjID=''E'' AND notes IN(''eye tracker'', ''eyetracker'', ''propixx'', ''human button box'')');
whichsubject = 'U';
%filenames = fetch(conn, ['CALL postPropixxFilenames(''',whichsubject,''')']);
filenames = fetch(conn, 'SELECT fileID FROM LMTF WHERE subjID=''U'' AND monID =''ProPixx'' and quality=1 AND nonstandard = 0');

close(conn);
[data, ~] = getLMTFrawdata(filenames);
CVerrors = [];
% Uncomment out to eliminate retinal locations (10,10) and (10,-10)
%L = abs(data(:,5)) == 100 & abs(data(:,6)) == 100;
%data(L,:) = [];

nuniqueXYs = size(unique(data(:,[5 6]),'rows','stable'),1);
NOOGidxs = Shuffle(find(data(:,4) == 0)); % "Not out of gamut indexes". Don't leave out an OOG point during CV; these rarely change the error
input_struct = {};
WARNINGS = {};
niter = min(niter,length(NOOGidxs));

fvs = [];
DEBUGGINGSTUFF = [];
for i = 1:niter
    i
    L = false(size(data,1),1); L(NOOGidxs(i)) = true;
    [output_struct, warning] = LMTF_generate_module_data(data(~L,:), false, input_struct); % Here's the workhorse
    if (sum(output_struct.legacy.mode4fvs) > sum(output_struct.legacy.mode5fvs))
        disp('Double tilted RT has more error than yoked');
        keyboard
    end
    WARNINGS{length(WARNINGS)+1} = warning;
    if (isempty(input_struct))
        input_struct = output_struct;
    end

    LXY = all(output_struct.eccs == repmat(data(L,[5 6]),size(output_struct.eccs,1),1),2);
    if sum(LXY) == 0 % If there is not enough data to fit the model at this location, skip this data point
        if isempty(CVerrors)
            CVerrors(i,:) = nan*ones(1,8);
            fvs(i,:) = nan*ones(1,8);
        else
            CVerrors(i,:) = nan*ones(1,size(CVerrors,2));
            fvs(i,:) = nan*ones(1,size(fvs,2));
        end
    else
        fvs(i,:) = [sum(output_struct.legacy.firstroundfvs),...
            sum(output_struct.legacy.mode0fvs),...
            sum(output_struct.legacy.mode1fvs),...
            sum(output_struct.legacy.mode1p1fvs),...
            sum(output_struct.legacy.mode1p2fvs),...
            sum(output_struct.legacy.mode2fvs),...
            sum(output_struct.legacy.mode3fvs),...
            sum(output_struct.legacy.mode4fvs),...
            sum(output_struct.legacy.mode5fvs)];
        CVerrors(i,:) = [tf_fiterr2(output_struct.legacy.firstroundmodels(:,LXY),data(L,:)),...
            tf_fiterr2(output_struct.legacy.mode0models(:,LXY),data(L,:)),...
            tf_fiterr2(output_struct.legacy.mode1models(:,LXY),data(L,:)),...
            tf_fiterr2(output_struct.legacy.mode1p1models(:,LXY),data(L,:)),...
            tf_fiterr2(output_struct.legacy.mode1p2models(:,LXY),data(L,:)),...
            tf_fiterr2(output_struct.legacy.mode2models(:,LXY),data(L,:)),...
            tf_fiterr2(output_struct.legacy.mode3models(:,LXY),data(L,:)),...
            tf_fiterr2(output_struct.legacy.mode4models(:,LXY),data(L,:)),...
            tf_fiterr2(output_struct.legacy.mode5models(:,LXY),data(L,:))];
    end
end
if SAVEFILES
    outfilename = [whichsubject,'_CVout'];
    eval(['save ',outfilename,' fvs CVerrors']);
end
[diffs, p] = prettypaircomp(CVerrors,{'unconstrained','11+2*n','10+3*n(theta)','10+3*n(1p1)', '10+3*n(1p2)', 'rampy trough','tilted rampy trough','double TRT', 'yoked dTRT'},'t-test');
title(['Subject: ',whichsubject,', Negative diff means column is superior to row']);

% Sanity checks on fitting error (fv)
L = ~all(isnan(fvs),2);
checks = [];
checks(1) = all(all(repmat(fvs(L,1),1,size(fvs,2)-1)<=fvs(L,[2:end]))); % firstround models should always have less error than any other fit
checks(2) = all(fvs(L,2)>=fvs(L,3)); % mode0_err > mode1_err. (Freeing theta should decrease error)
checks(3) = all(fvs(L,2)>=fvs(L,4)); % mode0_err  > mode1p1_err. (Freeing n_lum should decrease error)
checks(4) = all(fvs(L,2)>=fvs(L,5)); % mode0_err  > mode1p2_err. (Freeing n_rg should decrease error)
checks(5) = all(fvs(L,6)>=fvs(L,7)); % tilted rampy trough > rampy trough. (More df = less error)
checks(6) = all(fvs(L,7)>=fvs(L,8)); % double tilted rampy trough > tilted rampy trough. (More df = less error)
checks(7) = all(fvs(L,9)>=fvs(L,8)); % yoked double tilted rampy trough > free double tilted rampy trough. (More df = less error)

if (~all(checks))
    error('A constrained fit has lower *fitting error* than a flexible fit!');
end

%%
% Section 4.1
% Creating a model that IsoSamp.d can use.

SIDs = {  
    'P'};
SIDs = {'U', 'A'};

for i = 1:length(SIDs)
    conn = database('Nex_Paradigm_Sort','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
    command = sp
    rintf('CALL propixxStandardFilenames(''%s'')', SIDs{i});

    filenames = fetch(conn, command);
    if isempty(filenames)
        error = sprintf('User ID %s has no post-propixx files in the database', SIDs{i});
        disp(error);
        break;
    end
    close(conn);
    [data, check] = getLMTFrawdata(filenames);
    [output_struct, warning] = LMTF_generate_module_data(data, false);
    master_struct.(SIDs{i}) = output_struct;
end
if length(SIDs) == 1
    save(fullfile(fileparts(which('IsoSampOnline.m')), 'private', 'data', ['LMTF_',SIDs{1},'.mat']), ...
        '-struct', 'master_struct');
else
    save(fullfile(fileparts(which('IsoSampOnline.m')), 'private', 'data', 'LMTF.mat'), ...
        '-struct', 'master_struct');
end
%%
% Section 4.2
% Analysis of residuals. Assumed that section 4.1, above, has already been
% run so the models have already been fit.

% Labels for the columns of "data"
Lidx = 1;
Midx = 2;
TFidx = 3;
Xidx = 4;
Yidx = 5;
threshidx = 6;
predidx = 7;

load(fullfile(fileparts(which('IsoSampOnline.m')), 'private', 'data', 'LMTF.mat'))
SIDs = {'E','G','U','A'};
for i = 1:length(SIDs)
    s = eval(SIDs{i});
    neccs = length(s.raw);
    data = [];
    for j = 1:neccs
        Loog = logical(s.raw{j}(:,4));
        [~,r,f] = tf_fiterr2(s.legacy.mode5models(:,j), s.raw{j});
        data = [data; s.raw{j}(~Loog,[1:3,5:6]) r(~Loog) f(~Loog)];
    end  
    % Statistics we're going to need later
    resid = log10(data(:,threshidx))-log10(data(:,predidx)); % + means real threshold is too high (model overestimates senstivity)
    colordir = cart2pol(data(:,Lidx),data(:,Midx));
    
    % Standard plot of residuals against predicted valuesx
    figure; set(gcf,'Name',SIDs{i},'position',[70 590 1300 530]);
    subplot(2,3,1); hold on; 
    plot(log10(data(:,predidx)), resid,'.')
    xlabel('log10(pred threshold)');
    ylabel('log10(residual)');
    title(['Subject ',SIDs{i}]);

%     % Scatter plot of residual as a function of color direction and TF
%     subplot(2,3,2); hold on; 
%     plot3(colordir, log10(data(:,TFidx)),resid,'.');
%     set(gca,'Ylim',[0 max(log10(data(:,TFidx)))]);
%     set(gca,'Xtick',[0 pi/4 pi/2 3*pi/4],'Xticklabel',{'L', 'L+M', 'M' ,'L-M'});
%     %set(gca,'Zlim',[-.5 .5]); % Don't forget we're cuting off data here
%     xlabel('Color direction');
%     ylabel('log TF');
%     zlabel('residual');
%     axis vis3d
    subplot(2,3,2); hold on; 
    % Autocorrelation of the residuals as a function of color direction and TF
    thetabinedges = linspace(min(colordir),max(colordir),8);
    thetabinedges = 0:pi/8:pi;
    logTFbinedges = linspace(log10(1),log10(30),10);
    theta_TF_matrix = zeros(length(thetabinedges)-1,length(logTFbinedges)-1);
    % First calculating mean residuals in each bin
    for j = 1:length(thetabinedges)-1
        for k = 1:length(logTFbinedges)-1
            L = colordir>=thetabinedges(j) & colordir<thetabinedges(j+1) & log10(data(:,TFidx))>=logTFbinedges(k) & log10(data(:,TFidx))<logTFbinedges(k+1);
            theta_TF_matrix(j,k) = nanmean(resid(L));
        end
    end % nans where there are no residuals
    
    % Now calculating cross-correlation
    tmpmatrix = theta_TF_matrix;
    rs = zeros(length(thetabinedges)-1,length(logTFbinedges)-1);
    for j = 1:length(thetabinedges)-1
        for k = 1:length(logTFbinedges)-1
            tmpmatrix = theta_TF_matrix;
            if j > 1 % circular scrolling
                tmpmatrix = tmpmatrix([j:end,1:(j-1)],:);
            end
            if k > 1
                tmpmatrix = [nan(length(thetabinedges)-1,k-1),tmpmatrix(:,1:end-(k-1))];
            end
            r = nansum(nansum(theta_TF_matrix.*tmpmatrix))./sum(~isnan(theta_TF_matrix(:).*tmpmatrix(:)));
            rs(j,k) = r;
            if j == 1 & k == 1
                rs(j,k) = nan; % autocorrelation at zero lag is trivially large
            end
        end
    end
    imagesc(rs');
    axis tight
    set(gca,'Xtick',[0:2:length(thetabinedges)],'Xticklabel',round(10*thetabinedges(1:2:end))/10);
    set(gca,'Ytick',[0:2:length(logTFbinedges)],'Yticklabel',round(10.^logTFbinedges(1:2:end)));
    ylabel('TF'); xlabel('theta');

    % Bubble plot of mean (abs) residual as a function of retinal location
    subplot(2,3,3); hold on; set(gcf,'Name',SIDs{i});
    for j = 1:neccs
        L = all(data(:,[Xidx Yidx]) == repmat(s.eccs(j,:),size(data,1),1),2);
        resids = log10(data(L,threshidx))-log10(data(L,predidx));
        netresid = sign(median(resids)).*median(abs(resids));
        h = plot(s.eccs(j,1)/10,s.eccs(j,2)/10,'o');
        if (netresid > 0)
            set(h,'MarkerSize',50*netresid,'Color','black','MarkerFaceColor','black');
        else
            set(h,'MarkerSize',-50*netresid,'Color','red','MarkerFaceColor','red');
        end
    end
	title('Mean residuals. Red = underestimating sensitivity')
    xlabel('X');
    ylabel('Y');
    
    % regressing the residual against the predictors (reusing variable 'r')
    [phi,r] = cart2pol(data(:,Xidx)/10,data(:,Yidx)/10);
    [b,bint,~,~,stats] = regress(resid,[ones(length(colordir),1) cos(colordir) sin(colordir) log10(data(:,TFidx)) r r.*cos(2.*phi) r.*sin(2.*phi)]);
    labels = {'intercept','cos(colordir)','sin(colordir)','TF','R','R*cos(2phi)','R*sin(2*phi)'};
    
    %[b,bint,~,~,stats] = regress(resid,[ones(length(colordir),1) cos(colordir) sin(colordir) log10(data(:,TFidx))]);
    %labels = {'intercept','cos(colordir)','sin(colordir)','TF'};
    
    
    subplot(2,3,4);
    for j = 1:length(labels)
        h = text(0,j-1,[labels{j},': ',num2str([bint(j,1) bint(j,2)])]); % colordir, TF, X, Y
        if sign(bint(j,1)) == sign(bint(j,2))
            set(h,'FontWeight','bold')
        end
    end
    text(0,j,['F-test p-value ',num2str(stats(3))],'FontWeight','bold')

    text(0,j+1,'Confidence intervals for regressions on residuals','FontWeight','bold')
    set(gca,'Ylim',[0 j+1],'visible','off');
    
    % Comparing TCS peaks (not gain parameters) between retinal locations.
    lumpeakTFs = [];
    rgpeakTFs = [];
    options = optimset('Display','off');
    for i = 1:size(s.eccs,1)
        model = s.legacy.mode0models(:,i);
        lumf1 = @(omega)-1*(model(1)*abs(((1i*2*pi*10^model(5).*omega+1).^-model(3))-model(2)*((1i*2*pi*10^(model(5)+model(6)).*omega+1).^-(model(3)+model(4)))));
        rgf1 = @(omega)-1*(model(1+6)*abs(((1i*2*pi*10^model(5+6).*omega+1).^-model(3+6))-model(2+6)*((1i*2*pi*10^(model(5+6)+model(6+6)).*omega+1).^-(model(3+6)+model(4+6)))));
        model = s.legacy.mode3models(:,i);
        lumf3 = @(omega)-1*(model(1)*abs(((1i*2*pi*10^model(5).*omega+1).^-model(3))-model(2)*((1i*2*pi*10^(model(5)+model(6)).*omega+1).^-(model(3)+model(4)))));
        rgf3 = @(omega)-1*(model(1+6)*abs(((1i*2*pi*10^model(5+6).*omega+1).^-model(3+6))-model(2+6)*((1i*2*pi*10^(model(5+6)+model(6+6)).*omega+1).^-(model(3+6)+model(4+6)))));
        [~,fv1] = fmincon(lumf1,2,[],[],[],[],1,60,[],options); [~,fv3] = fmincon(lumf3,5,[],[],[],[],1,60,[],options);
        lumpeakTFs(i,:) = [-1*fv1, -1*fv3];
        [~,fv1] = fmincon(rgf1,2,[],[],[],[],1,60,[],options); [~,fv3] = fmincon(rgf3,5,[],[],[],[],1,60,[],options);
        rgpeakTFs(i,:) = [-1*fv1, -1*fv3];
    end
    
    % Comparing LUM peak sensitivity between mode1 and mode3 model (in which is must follow the tilted rampy trough)
    subplot(2,3,5); hold on;
    L = sign(lumpeakTFs(:,1)-lumpeakTFs(:,2)) == 1; % comparing gains, not thresholds. 1 = we're underestimating sensitivity with mode3 mdel
    plot3(s.eccs(:,1),s.eccs(:,2),lumpeakTFs(:,1),'rs','MarkerFaceColor','red','MarkerSize',4);
    plot3(s.eccs(:,1),s.eccs(:,2),lumpeakTFs(:,2),'ks','MarkerFaceColor','black','MarkerSize',4);
    plot3([s.eccs(L,1) s.eccs(L,1)]',[s.eccs(L,2) s.eccs(L,2)]',[lumpeakTFs(L,1) lumpeakTFs(L,2)]','r-','Linewidth',2)
    plot3([s.eccs(~L,1) s.eccs(~L,1)]',[s.eccs(~L,2) s.eccs(~L,2)]',[lumpeakTFs(~L,1) lumpeakTFs(~L,2)]','b-','Linewidth',2)
    zlabel('LUM peak sensitivity'); axis vis3d; set(gca,'Xlim',[0 200],'Ylim',[-100 100])
    title('Red: LUM gain is too low (underestimates sens.)');
    set(gca,'Zlim',[0 ceil(max(lumpeakTFs(:)/10))*10],'Zscale','log');
    
    % Comparing RG gain between mode1 model (in which it's free to vary)
    % and mode3 model (in which is must follow a rampy trough)
    subplot(2,3,6); hold on;
    L = sign(rgpeakTFs(:,1)-rgpeakTFs(:,2)) == 1; % comparing gains, not thresholds. 1 = we're underestimating sensitivity with mode3 mdel
    plot3(s.eccs(:,1),s.eccs(:,2),rgpeakTFs(:,1),'rs','MarkerFaceColor','red','MarkerSize',4);
    plot3(s.eccs(:,1),s.eccs(:,2),rgpeakTFs(:,2),'ks','MarkerFaceColor','black','MarkerSize',4);
    plot3([s.eccs(L,1) s.eccs(L,1)]',[s.eccs(L,2) s.eccs(L,2)]',[rgpeakTFs(L,1) rgpeakTFs(L,2)]','r-','Linewidth',2)
    plot3([s.eccs(~L,1) s.eccs(~L,1)]',[s.eccs(~L,2) s.eccs(~L,2)]',[rgpeakTFs(~L,1) rgpeakTFs(~L,2)]','b-','Linewidth',2)
    zlabel('RG peak sensitivity'); axis vis3d; set(gca,'Xlim',[0 200],'Ylim',[-100 100]);
    set(gca,'Zlim',[0 ceil(max(rgpeakTFs(:)/10))*10],'Zscale','log');
    title('Red: RG gain is too low');
end


%%
% Section 4.3
% Finding stimuli that one observer can see relatively easily and another
% cannot (under their respective models). Testing these stimuli empircally
% might be useful.

MODE = 5;
nseeds = 10; % Number of random initial guesses

load(fullfile(fileparts(which('IsoSampOnline.m')), 'private', 'data', 'LMTF.mat'))
SIDs = {'E','G','A','U'}; % Order must be human, human, monkey, monkey
% we're minimizing (-thresh(obs 1))-(-thresh(obs 2))
% or equivalently thresh(obs 2)-thresh(obs 1). Trying to a condition where
% the 1st observer is more sensitive than the 2nd (without taking abs())
for i = 1:length(SIDs)
    models(:,i) = eval([SIDs{i},'.legacy.mode',num2str(MODE),'params']); % global models
end
%compareLMTFmodelpreds([coldir TF X Y], models)
compareLMTFmodelpreds = @(data, model1, model2)(tf_fiterr3(model1,data,MODE)-tf_fiterr3(model2,data,MODE)); % Errors are already in log space
compareLMTFmodelpreds2 = @(data, model1_h, model2_h, model1_m, model2_m)...
    ((tf_fiterr3(model1_h,data,MODE)+tf_fiterr3(model2_h,data,MODE))-...
     (tf_fiterr3(model1_m,data,MODE)+tf_fiterr3(model2_m,data,MODE))); % humans vs monkeys

% Logic: tf_fiterr3 actually returns the *error*, but the "true"
% threshold in this script is always "1" because L^2+M^2 is constrained to be 1.
% So the error returned by tf_fiterr3 is log10(data point)-log10(pred) = log10(1)-log10(pred)
% -log10(pred), where pred is distance to surface
% So compareLMTFmodelpreds returns the difference in predictions between
% two models for a single stimulus direction. We want to maximize the
% absolute value of this difference (minimize -1*diff).

%compareLMTFmodelpreds([20 20 3 0 30 30],G.model,E.model) % Testing

nlincon = @(x) deal([],(x(1)^2+x(2)^2)-1); % output: c, ceq
max_diff = 0;
bestvect = [];
LB = [-1 -1 1  0 20 -50];% L, M, TF, LOOG, X, Y
UB = [ 1  1 10 0 50  50];
for i = 1:nseeds
    theta = unifrnd(-pi,pi); % random direction in LM plane
    initialguess = [cos(theta), sin(theta), unifrnd(LB(:,[3:end]), UB(:,[3:end]))]; % L, M, TF, LOOG, X, Y
    options = optimset('Algorithm', 'interior-point', 'MaxFunEvals', 5e4, ...
        'MaxIter', 5e4, 'TolFun', 1e-8, 'Display', 'off','Display','iter','MaxFunEvals',10^6);
    %[maxdiffvect,fv] = fmincon(@(params) -1*abs(compareLMTFmodelpreds(params, models(:,1),models(:,2))), initialguess,...
    %    [],[],[],[],LB,UB,nlincon,options);
    [maxdiffvect,fv] = fmincon(@(params) -1*abs(compareLMTFmodelpreds2(params, models(:,1),models(:,2),models(:,3),models(:,4))),...
        initialguess,[],[],[],[],LB,UB,nlincon,options);
    if fv < max_diff
        max_diff = fv;
        bestvect = maxdiffvect
    end
end

10^max_diff % Fold difference in contrast sensitivity
bestvect

%%
% Section 5 
% F-tests comparing unconstrained model to 11+2n model
% Currently just using a standard F-test, but should be 
% changed to use a bootstrap F-test.

SID = 'A';

% Loading the LMTF model structure
isosamppath = which('IsoSampOnline');
isosamppath(find(isosamppath==filesep,1,'last'):end) = [];
load ([isosamppath,filesep,'private',filesep','data',filesep,'LMTF.mat']);

conn = database('Nex_Paradigm_Sort','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filenames = fetch(conn, sprintf('CALL postPropixxFilenames(''%s'')',SID));
close(conn);
[data, check] = getLMTFrawdata(filenames);
uniqueXYs = eval([SID,'.eccs']);
nuniqueXYs = size(uniqueXYs,1);

unconstrained_params = eval([SID,'.legacy.firstroundmodels']);
constrained_params = eval([SID,'.legacy.mode2models']);

% We assume that the simple model is true and ask whether the more complex
% model accounts for more variance than expected

% First, getting residuals for each data point under each model
resids = zeros(size(data,1),2);
for i = 1:nuniqueXYs
    L = data(:,5) == uniqueXYs(i,1) & data(:,6) == uniqueXYs(i,2);
    pred1 = LMTF_thresh_from_model(data(L,[1 2 3]),unconstrained_params(:,i));
    pred2 = LMTF_thresh_from_model(data(L,[1 2 3]),constrained_params(:,i));
    resids(L,1) = log10(sqrt(sum(data(L,[1 2]).^2,2)))-log10(pred1);
    resids(L,2) = log10(sqrt(sum(data(L,[1 2]).^2,2)))-log10(pred2);
end
Loog = data(:,4) == 1; % OOG
resids(Loog,:) = nan; % not counting OOG points in this analysis

% Sanity checking
figure; axes; hold on;
plot(abs(resids(:,1)),abs(resids(:,2)),'.');
plot([0 1],[0 1],'k-');
xlabel('abs(unconstrained residuals)');
ylabel('abs(constrained residuals)');
nansum(abs(resids))
[sum(eval([SID,'.legacy.firstroundfvs'])) sum(eval([SID,'.legacy.mode2fvs']))] % includes OOGs

% Classical F-test
n = sum(~Loog);
RSS = nansum(resids.^2);
p1 = numel(unconstrained_params) % # parameters of unconstrained model
p2 = 11+2*size(constrained_params,2);
num = (RSS(2)-RSS(1))/(p1-p2);
den = RSS(1)/(n-p1);
F = num/den;
p = 1-fcdf(F,p1-p2,n-p1)


%%
% Section 5.1
% As above but now a bootstrap F-test
% Not resampling OOG points, but still using them in the fits.

SID = 'E';

% Loading the LMTF model structure
isosamppath = which('IsoSampOnline');
isosamppath(find(isosamppath==filesep,1,'last'):end) = [];
load ([isosamppath,filesep,'private',filesep','data',filesep,'LMTF.mat']);

conn = database('Nex_Paradigm_Sort','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filenames = fetch(conn, sprintf('CALL postPropixxFilenames(''%s'')',SID));
close(conn);
[data, check] = getLMTFrawdata(filenames);
uniqueXYs = eval([SID,'.eccs']);
nuniqueXYs = size(uniqueXYs,1);
% sorting the data
tmp = [];
for i = 1:nuniqueXYs
    L = data(:,5) == uniqueXYs(i,1) & data(:,6) == uniqueXYs(i,2);
    tmp = [tmp; data(L,:)];
end
data = tmp;
clear tmp;

unconstrained_params = eval([SID,'.legacy.firstroundmodels']); % unconstrained
constrained_params = eval([SID,'.legacy.mode0models']); % 11+2n

% We assume that the simple model is true and ask whether the more complex
% model accounts for more variance than expected

% First, getting residuals for each data point under the simpler (mode 2)
% model
LM_TF_Loog_pred_resid = {};
residuals = []; % One gigantic vector of residuals (not OOGs) that we will use in resampling
for i = 1:nuniqueXYs
    L = data(:,5) == uniqueXYs(i,1) & data(:,6) == uniqueXYs(i,2);
    pred = LMTF_thresh_from_model(data(L,[1 2 3]),constrained_params(:,i));
    resid = log10(sqrt(sum(data(L,[1 2]).^2,2)))-log10(pred);
    LM_TF_Loog_pred_resid{i} = [data(L,[1:4]) log10(pred), resid];
    residuals = [residuals; resid(logical(~data(L,4)))]; % ~OOG residuals for resampling
end
% getting the F-statistic from the real data
realF = [sum(eval([SID,'.legacy.mode0fvs'])) sum(eval([SID,'.legacy.firstroundfvs']))];

niter = 200;
% First make a fake data set. Then fit mode 2 model using the "truth"
% as an initial guess. Then fit the unconstrained model from the best

% mode 2 fit (and the original fit?)
% Setting up for fitting mode0 model to resampled data
mode0models = eval([SID,'.legacy.mode0models']);
initialguess = [mode0models([2:6, 8:13],2); reshape(mode0models([1 7],:),2*nuniqueXYs,1)];
options = optimset('Display','final','MaxFunEvals',10^5,'TolFun',10^-4,'MaxIter',10^5);
Fs = zeros(niter+1,2);
Fs = [];
parfor i = 0:niter
    i
    fakedata = []; % [L, M, TF, Loog, RFX, RFY]
    if i == 0
        resids = residuals;
    else
        resids = residuals(unidrnd(length(residuals),size(residuals)));
    end

    for j = 1:nuniqueXYs
        tmp = LM_TF_Loog_pred_resid{j}(:,[1:4]);
        Loog = tmp(:,4) == 1;
        preds = LM_TF_Loog_pred_resid{j}(~Loog,5); % log scaled
        modifiedpreds = preds+resids(1:length(preds));
        resids(1:length(preds)) = []; % popping residuals off the stack
        unitvectors = mkbasis(tmp(~Loog,[1 2])')';
        tmp(~Loog,[1 2]) = unitvectors.*repmat(10.^modifiedpreds,1,2);
        tmp(:,[5 6]) = repmat(uniqueXYs(j,:),size(tmp,1),1);
        fakedata = [fakedata; tmp];
    end
    if length(resids) ~= 0
        error('too many resids');
    end
    [fpar,fv_constrained,exitflag] = fminsearch(@(params) tf_fiterr3(params,fakedata,0),initialguess,options);
    if (exitflag ~= 1)
        error('Fitting did not converge');
    end
    % Now the unconstrained models
    ximat = reshape(fpar(12:end),2,nuniqueXYs);
    inmodels = [];
    fvs = [];
    LB = [0 0 1 0 -3 0 0 0 1 0 -3 0 0];
    UB = [500 1 40 10 -1 0.3010 500 1 40 10 -1 0.3010 1.5708];
    
    for j = 1:nuniqueXYs
        L = fakedata(:,5) == uniqueXYs(j,1) & fakedata(:,6) == uniqueXYs(j,2);
        ig = [ximat(1,j); fpar(1:5); ximat(2,j); fpar(6:11)];
        [modelparams,fv] = fmincon(@(params) tf_fiterr2(params,fakedata(L,[1:4])),ig, [],[],[],[],LB,UB,[],options);
        inmodels = [inmodels, modelparams];
        fvs = [fvs,fv];
        %original_unconstrained_fits = eval([SID,'.legacy.firstroundmodels']);
        %ig2 = original_unconstrained_fits(:,j);
        %[~,fv2] = fminsearch(@(params) tf_fiterr2(params,fakedata(L,[1:4])),ig2,options);
         % The original model still has an unfair advantage because we
         % worked harder to fit it. Need to incorporate
         % "CleanupFirstRoundFits" from generate module data
    end
    [outmodels, fvs] = CleanupFirstRoundFits(uniqueXYs, fakedata, inmodels, fvs)

    Fs(i+1,:) = [fv_constrained sum(fvs)]
end

figure; subplot(2,1,1); hold on;
plot(Fs([2:end],1),Fs([2:end],2),'bo')
plot(Fs(1,1),Fs(1,2),'m*') % First row is real data (residuals are not randomized)
xlabel('Constrained error');
ylabel('Unconstrained error');
axis square
title(SID);
subplot(2,1,2); hold on;
hist(Fs([2:end],1)./Fs([2:end],2),50);
plot(Fs(1,1)/Fs(1,2),0,'m*')
p = sum(Fs([2:end],1)./Fs([2:end],2) > Fs(1,1)/Fs(1,2))/(size(Fs,1)-1);
title(['p = ',num2str(p)]);
%%
% Section 5.2
% Comparing distributions of residuals across retinal locations
% Emily did this. No significant differences by Kruskal-Wallis for any
% observer (as long as OOG points are not considered which only makes
% sense).
%%
% Section 6
% Comparing 1-D temporal contrast senstivity functions across observers

mode = 5;
RFs = [5 0; 10 -10];
load /Users/greghorwitz/Desktop/MatlabCode/Zack/IsoSamp/private/data/LMTF.mat
SIDs = {'A','E','G','U'};
colors = [111 96 170; 224 126 39; 252 192 127; 175 180 219]/255;
tfs = logspace(log10(1),log10(40),100);
figure; 
h = [];
for i = 1:size(RFs,1)
    subplot(1, size(RFs,1),i); hold on;
    for j = 1:length(SIDs)
        s = eval(SIDs{j});
        globalparams = getfield(s.legacy,['mode',num2str(mode),'params']);
        model = LMTF_global_to_local_model(globalparams, RFs(i,1), RFs(i,2), mode);
        lumf1 = @(omega)(model(1)*abs(((1i*2*pi*10^model(5).*omega+1).^-model(3))-model(2)*((1i*2*pi*10^(model(5)+model(6)).*omega+1).^-(model(3)+model(4)))));
        rgf1 = @(omega)(model(1+6)*abs(((1i*2*pi*10^model(5+6).*omega+1).^-model(3+6))-model(2+6)*((1i*2*pi*10^(model(5+6)+model(6+6)).*omega+1).^-(model(3+6)+model(4+6)))));
        h(j) = plot(tfs,lumf1(tfs),'-','LineWidth',2,'Color',colors(j,:));
        plot(tfs,rgf1(tfs),':','LineWidth',2,'Color',colors(j,:));
    end
    title(['RF @ (',num2str(RFs(i,:)),')']);
    set(gca,'xscale','log','Yscale','log');
    set(gca,'Xlim',[tfs(1) tfs(end)]);
end
equatesubplotaxeslims;
legend(h,SIDs)
% sensitivity = 1/sqrt(Lcc^2+Mcc^2)
% so at threshold, assuming Lcc = Mcc, Lcc = 1/(sqrt(2)*sensitivity)
s = load('ProPixx.mat');
cal = s.cals{end};
s = load('T_cones_smj10');
fundamentals = s.T_cones_smj10;
M = fundamentals*SplineSpd([380:4:780]', cal.P_device, [380:5:780]');
% Ignoring gamma functions since this is ProPixx
unitvector = [1/sqrt(2) -1/sqrt(2) 0];
[in_gamut,gamut_scalars] = gamutCheck(unitvector, cal.bgColor, M,'both');
highest_contrast_chromatic_stim = unitvector*gamut_scalars;
plot([tfs(1) tfs(end)],[1 1]./norm(highest_contrast_chromatic_stim),'r:'); % Lower testable bound for chromatic sensitivity

%%
% Section 7
% Temporal integration times
STIMTYPE = 'RG'; % 'LUM' or 'RG'
MONKEY = 'A';
conn = database('Nex_Paradigm_Sort','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
if strcmp(STIMTYPE, 'LUM')
    longcomment = 'std inputs 666 lum';
    shortcomment = 'std inputs 333 lum';
else
    longcomment = 'std inputs 666 chr';
    shortcomment = 'std inputs 333 chr';
end
query = sprintf('SELECT fileID FROM LMTF WHERE subjID = ''%s'' AND rfX = 50 AND rfY = 0 AND quality = 1 AND notes LIKE (''%s'')',MONKEY,longcomment);
flist_long = fetch(conn, query);
query = sprintf('SELECT fileID FROM LMTF WHERE subjID = ''%s'' AND rfX = 50 AND rfY = 0 AND quality = 1 AND notes LIKE (''%s'')',MONKEY,shortcomment);
flist_short = fetch(conn, query);
close(conn)
if ~isempty(flist_long)
    data_long = getLMTFrawdata(flist_long);
else
    data_long= [];
end
if ~isempty(flist_short)
    data_short = getLMTFrawdata(flist_short);
else
    data_short= [];
end

data = [sqrt(data_long(:,1).^2+data_long(:,2).^2) data_long(:,3) repmat(666,size(data_long,1),1)];
if ~isempty(data_short)
    data = [data; sqrt(data_short(:,1).^2+data_short(:,2).^2) data_short(:,3) repmat(333,size(data_short,1),1)];
end
% Loading the model
isosamppath = which('IsoSampOnline');
isosamppath(find(isosamppath==filesep,1,'last'):end) = [];
load ([isosamppath,filesep,'private',filesep','data',filesep,'LMTF.mat']);
rfx = unique(data_long(:,5))/10;
rfy = unique(data_long(:,6))/10;

if MONKEY == 'A'
    model = LMTF_global_to_local_model(A.legacy.mode5params, rfx, rfy, 5);
else
    model = LMTF_global_to_local_model(U.legacy.mode5params, rfx, rfy, 5);
end
if strcmp(STIMTYPE,'RG')
    pred = @(omega)1./(model(1+6)*abs(((1i*2*pi*10^model(5+6).*omega+1).^-model(3+6))-model(2+6)*((1i*2*pi*10^(model(5+6)+model(6+6)).*omega+1).^-(model(3+6)+model(4+6)))));
else
    pred = @(omega)1./(model(1)*abs(((1i*2*pi*10^model(5).*omega+1).^-model(3))-model(2)*((1i*2*pi*10^(model(5)+model(6)).*omega+1).^-(model(3)+model(4)))));
end

figure; axes; hold on;
tmp_tfs = logspace(log10(1),log10(20),100);
plot(tmp_tfs,pred(tmp_tfs),'k-');
plot(tmp_tfs,pred(tmp_tfs)*sqrt(2),'r:');
TFs = unique(data(:,2),'sorted');
durations = unique(data(:,3),'sorted');
colors = [1 0 0; 0 0 0];
set(gca,'Yscale','log','Xscale','log');
for i = 1:length(TFs)
    tmp = [];
    for j = 1:length(durations)
        L = data(:,2) == TFs(i) & data(:,3) == durations(j);
        mn = 10.^mean(log10(data(L,1)));
        sem = sqrt(var(log10(data(L,1)))/sum(L));
        %plot(TFs(i),data(L,1),'x','MarkerFaceColor',colors(j,:),'MarkerEdgeColor',colors(j,:),'MarkerSize',2); % Individual data points
        plot(TFs(i)+.02*(j-1),median(data(L,1)),'s','MarkerFaceColor',colors(j,:),'MarkerEdgeColor',colors(j,:),'MarkerSize',10);
        plot(TFs(i)*[1 1]+.02*(j-1),[prctile(data(L,1),25) prctile(data(L,1),75)],'-','Color',colors(j,:),'Linewidth',1);
        plot(TFs(i)*[.95 1.05]+.02*(j-1),[prctile(data(L,1),25) prctile(data(L,1),25)],'-','Color',colors(j,:),'Linewidth',1);
        plot(TFs(i)*[.95 1.05]+.02*(j-1),[prctile(data(L,1),75) prctile(data(L,1),75)],'-','Color',colors(j,:),'Linewidth',1);

        tmp = [tmp; j*ones(sum(L),1) data(L,1)];
        %plot(TFs(i),mn,'o','MarkerFaceColor',colors(j,:),'MarkerEdgeColor',colors(j,:));
        %plot([TFs(i) TFs(i)],10.^[log10(mn)+sem log10(mn)-sem],'-','Color',colors(j,:),'LineWidth',2);
    end
    L1 = tmp(:,1) == 1;
    L2 = tmp(:,1) == 2;
 
    p = ranksum(tmp(L1,2),tmp(L2,2),'method','exact','tail','right'); % One-tailed assuming lower duration is first (true because of sorting)
    p
    if (p < 0.05)
       plot(TFs(i),median(tmp(L1,2)),'ys');
       plot(TFs(i),median(tmp(L2,2)),'ys');
    end
end
title(MONKEY);
set(gca,'Xlim',[TFs(1)*.9 TFs(end)*1.1]);
if strcmp(STIMTYPE, 'LUM')
    set(gca,'Ylim',[.04 .24]);
else
    set(gca,'Ylim',[.02 .2]);
end


% tmp = [];
% CIs = [];
% for i = 1:length(TFs)
%     L_short = data(:,2) == TFs(i) & data(:,3) == durations(1);
%     L_long = data(:,2) == TFs(i) & data(:,3) == durations(2);
%     tmp(i) = mean(log10(data(L_short,1)))-mean(log10(data(L_long,1))) % 333 threshold-666 threshold
%     [~,~,CI] = ttest2(log10(data(L_short,1)),log10(data(L_long,1)))
%     CIs(:,i) = CI;
% end
% figure; axes; hold on;
% plot(TFs,10.^tmp);
% CImat = 10.^CIs-repmat(10.^tmp,2,1);
% errorbar(TFs,10.^tmp,CImat(1,:),CImat(2,:))
% set(gca,'Xscale','log')
% plot([TFs(1) TFs(end)],sqrt(2)*[1 1],'-');

%%
% Section 8
% Pie charts of number of threshold points per subject across visual field
% locations
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
location_query = 'SELECT rfX, rfY, subjID FROM LMTF WHERE recDate > ''2016-05-27'' AND recDate < ''2017-02-03'' AND quality = 1 ORDER BY `LMTF`.`subjID` ASC';
alldata = fetch(conn, location_query);
close(conn);
allsubjects = cell2mat(alldata(:,3));
rfxy = cell2mat(alldata(:,[1 2]));
uniquexys = unique(rfxy,'rows');

% Setting up dummy variables for subjects
subjects = {'A','U','E','G'};
colors = [.9 .4 .2;0 0 1; 1 0 1; 0 0 0];
n_subjects = length(subjects);
n_files = zeros(size(uniquexys,1),length(subjects))

for i = 1:size(uniquexys,1)
    L = rfxy(:,1) == uniquexys(i,1) & rfxy(:,2) == uniquexys(i,2)
    for j = 1:length(subjects)
        n_files(i,j) = sum(L & allsubjects == subjects{j});
    end
end

% Now plotting based on number of files/number of subjects
% Ugly plot
% figure; axes; hold on;
% tmp = [cos(linspace(0,2*pi,100))' sin(linspace(0,2*pi,100))'];
% for i = 1:size(uniquexys,1)
%     r = sum(n_files(i,:),2);
%     plot(uniquexys(i,1)+tmp(:,1)*r, uniquexys(i,2)+tmp(:,2)*r,'-')
% end
% axis equal

% Pie chart all same size
figure; axes; hold on;
r = 5;
tmp = [cos(linspace(0,2*pi,100))' sin(linspace(0,2*pi,100))'];
for i = 1:size(uniquexys,1)
    Lsubject = n_files(i,:)> 0;
    nsubjects = sum(Lsubject);
    %plot(uniquexys(i,1)+tmp(:,1)*r, uniquexys(i,2)+tmp(:,2)*r,'k-');
    for j = 1:nsubjects
        subjectidxs = find(Lsubject);
        tmp = linspace((j-1)/nsubjects,j/nsubjects,50);
        tmpx = [0 r*cos(2*pi*tmp)];
        tmpy = [0 r*sin(2*pi*tmp)];
        h = patch(uniquexys(i,1)+tmpx,uniquexys(i,2)+tmpy,colors(subjectidxs(j),:));
        set(h,'EdgeColor',colors(subjectidxs(j),:))
    end
end
axis equal

sum(n_files > 0) % Number of locations tested per subject

%%
% Section 9
% Comparing eye tracker vs. human button box

conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
E_hbb_flist = fetch(conn, 'SELECT fileID FROM LMTF WHERE subjID = ''E'' AND monID = ''ProPixx'' AND notes like ''%human button box%'' AND quality = 1'); 
E_eyetracker_flist = fetch(conn, 'SELECT fileID FROM LMTF WHERE subjID = ''E'' AND monID = ''ProPixx'' AND notes like ''%eye tracker%'' AND quality = 1'); 
G_hbb_flist = fetch(conn, 'SELECT fileID FROM LMTF WHERE subjID = ''G'' AND monID = ''ProPixx'' AND notes like ''%human button box%'' AND quality = 1'); 
G_eyetracker_flist = fetch(conn, 'SELECT fileID FROM LMTF WHERE subjID = ''G'' AND monID = ''ProPixx'' AND notes like ''%eye tracker%'' AND quality = 1'); 

% counting up the file to make sure we have the right number
length(E_hbb_flist)+length(E_eyetracker_flist)
length(G_hbb_flist)+length(G_eyetracker_flist)

% seeing how many files we have at individual screen locations with HBB and
% eye tracker for each subject.
subjIDs = ['E','G'];
data = [];
for i = 1:2
    HBB_rf_locations = cell2mat(fetch(conn, ['SELECT rfX, rfY FROM LMTF WHERE subjID = ''',subjIDs(i),''' AND monID = ''ProPixx'' AND notes like ''%human button box%'' AND quality = 1']));
    EyeTracker_rf_locations = cell2mat(fetch(conn, ['SELECT rfX, rfY FROM LMTF WHERE subjID = ''',subjIDs(i),''' AND monID = ''ProPixx'' AND notes like ''%eye tracker%'' AND quality = 1']));
    uniquelocations = unique([HBB_rf_locations; EyeTracker_rf_locations],'rows');
    for j = 1:size(uniquelocations)
        hbb_n = sum(all(HBB_rf_locations == repmat(uniquelocations(j,:),size(HBB_rf_locations,1),1),2));
        et_n = sum(all(EyeTracker_rf_locations == repmat(uniquelocations(j,:),size(EyeTracker_rf_locations,1),1),2));
        data(j,:,i) = [uniquelocations(j,:), hbb_n, et_n];
    end
end

% For subject E, the only location where we have both eye tracker and button
% box data is 50,0. For subject G, we have both eye tracker and button
% box data from (50, 0), (50, -60), (50, 60) and (80, 0).

data = []; % reusing this variable
conditions = [{'E', 50, 0}; {'G', 50, 0}; {'G', 50, -60}; {'G', 50, 60}; {'G', 80, 0}];
for i = 1:size(conditions,1)
    str = ['SELECT fileID FROM LMTF WHERE subjID = ''',conditions{i,1},''' AND monID = ''ProPixx'' AND notes like ''%eye tracker%'' AND rfX = ',num2str(conditions{i,2}),' AND rfY = ',num2str(conditions{i,3}),' AND quality = 1'];
    eyetracker_flist = fetch(conn, str);
    str = ['SELECT fileID FROM LMTF WHERE subjID = ''',conditions{i,1},''' AND monID = ''ProPixx'' AND notes like ''%human button box%'' AND rfX = ',num2str(conditions{i,2}),' AND rfY = ',num2str(conditions{i,3}),' AND quality = 1'];
    HBB_flist = fetch(conn, str);
    et_data = getLMTFrawdata(eyetracker_flist);
    hbb_data = getLMTFrawdata(HBB_flist);
    
    % First just plotting as a sanity check
%     figure; axes; hold on;
%     plot3(et_data(:,1),et_data(:,2),et_data(:,3),'ko');
%     plot3(-et_data(:,1),-et_data(:,2),et_data(:,3),'ko');
%     plot3(hbb_data(:,1),hbb_data(:,2),hbb_data(:,3),'mo');
%     plot3(-hbb_data(:,1),-hbb_data(:,2),hbb_data(:,3),'mo');
%     set(gca,'Zscale','log')
    
    % How many exact repeats do we have?
    [et_theta, et_r] = cart2pol(et_data(:,1),et_data(:,2));
    [hbb_theta, hbb_r] = cart2pol(hbb_data(:,1),hbb_data(:,2));
    thresholds = [];
    for i = 1:length(et_theta)
        L = et_theta(i) == hbb_theta & et_data(i,3) == hbb_data(:,3) & ~hbb_data(:,4) & ~et_data(i,4);
        if any(L)
            thresholds = [thresholds; et_r(i) hbb_r(L)];
        end
    end
    figure;
    plot(log10(thresholds(:,1)),log10(thresholds(:,2)),'o');
    [~,p] = ttest(log10(thresholds(:,1))-log10(thresholds(:,2)));
    r = corrcoef([log10(thresholds(:,1)) log10(thresholds(:,2))]);
    data = [data; p r(1,2) size(thresholds,1)];
end
data
close(conn);
% Looks like only the 50,0 condition had enough data.

% Figuring out which data were collected when
dates_hbb = []; dates_et = [];
hbb_flist = E_hbb_flist;
eyetracker_flist = E_eyetracker_flist;
for i = 1:size(hbb_flist,1)
    tmp = hbb_flist{i}(2:7);
    dates_hbb = [dates_hbb; datetime(2000+str2num(tmp(5:6)),str2num(tmp(1:2)),str2num(tmp(3:4)))];
end
for i = 1:size(eyetracker_flist,1)
    tmp = eyetracker_flist{i}(2:7);
    dates_et = [dates_et; datetime(2000+str2num(tmp(5:6)),str2num(tmp(1:2)),str2num(tmp(3:4)))];
end
    
sum(dates_hbb<min(dates_et))
sum(dates_hbb>max(dates_et))
figure; axes; hold on;
plot(dates_hbb,'.-');
plot(dates_et,'.');


%%
% Section 10
% Slices through fitted surface at one location

load /Users/greghorwitz/Desktop/MatlabCode/Zack/IsoSamp/private/data/LMTF.mat;
SID = 'U';
RF = [50 0];

s = eval(SID);
Lecc = s.eccs(:,1) == RF(1) & s.eccs(:,2) == RF(2);
raw = s.raw{Lecc};
modelparams = s.legacy.firstroundmodels(:,Lecc);
%lumf1 = @(omega)(model(1)*abs(((1i*2*pi*10^model(5).*omega+1).^-model(3))-model(2)*((1i*2*pi*10^(model(5)+model(6)).*omega+1).^-(model(3)+model(4)))));
%rgf1 = @(omega)(model(1+6)*abs(((1i*2*pi*10^model(5+6).*omega+1).^-model(3+6))-model(2+6)*((1i*2*pi*10^(model(5+6)+model(6+6)).*omega+1).^-(model(3+6)+model(4+6)))));
%theta = modelparams(end);

% First detection ellipses in LM plane
TFs = [1 5 10 20];
margin = 1.5;
lm = [cos(linspace(0,2*pi,100))', sin(linspace(0,2*pi,100))'];
colors = linspace(.8,0,length(TFs))'*[1 1 1];
colors = [1 0 0; 0 1 0; 0 .5 1; 0 0 0];

symbols = ['o','+','^','x'];
figure; axes; hold on;
for i = 1:length(TFs)
    pred = LMTF_thresh_from_model([lm, TFs(i)*ones(length(lm),1)],modelparams);
    plot(pred.*lm(:,1), pred.*lm(:,2),'-','Color',colors(i,:),'LineWidth',2);
    Loog = raw(:,4);
    L = (raw(:,3) > TFs(i)/margin) & (raw(:,3) < TFs(i)*margin) & ~Loog;
    if any(L)
       % plot(raw(L,1),raw(L,2),'ko','MarkerSize',8,'MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:));
       % plot(-raw(L,1),-raw(L,2),'ko','MarkerSize',8,'MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:));
        plot(raw(L,1),raw(L,2),symbols(i),'MarkerSize',8,'MarkerEdgeColor',colors(i,:));
        plot(-raw(L,1),-raw(L,2),symbols(i),'MarkerSize',8,'MarkerEdgeColor',colors(i,:));
    end
end
axis square
set(gca,'Xlim',[-.3 .3],'Ylim',[-.3 .3],'TickDir','out');

% Second contrast sensitivity functions
thetas = [0 pi/4 pi/2 3*pi/4];
TFs = logspace(log10(1),log10(40),100)';
margin = 10*pi/180;
colors = [1 0 0; 0 0 0; 0 .75 .75; .5 .5 .5];
figure; axes; hold on;
for i = 1:length(thetas)
    pred = LMTF_thresh_from_model([cos(thetas(i))*ones(length(TFs),1),sin(thetas(i))*ones(length(TFs),1),TFs],modelparams);
    plot(TFs, 1./pred,'-','Color',colors(i,:),'LineWidth',2);
    t = cart2pol(raw(:,1),raw(:,2));
    L = abs(angdiff(t,repmat(thetas(i),size(t)))) < margin & ~raw(:,4);
    if any(L)
        plot(raw(L,3),1./sqrt(raw(L,1).^2+raw(L,2).^2),'ko','MarkerSize',5,'MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:));
    end
end
set(gca,'Yscale','log','Xscale','log','Xlim',[min(TFs) max(TFs)],'Ylim',[1 100],'TickDir','out');
axis square;
set(gcf,'renderer','painters');

% Double checking Propixx gamut to make sure all of these measured thresholds are inside. Checks out.
load('ProPixx')
cal = cals{end};
bkgndrgb = cal.bgColor;
fundamentals = load('T_cones_smj10');
funds = fundamentals.T_cones_smj10;
spds = SplineSpd([380:4:780]',cal.P_device,[380:5:780]')
M = funds*spds;
bkgndlms = M*bkgndrgb;
cc = [-0.0545 0.2859 0]';
lms = cc.*bkgndlms+bkgndlms;
rgb = inv(M)*lms
