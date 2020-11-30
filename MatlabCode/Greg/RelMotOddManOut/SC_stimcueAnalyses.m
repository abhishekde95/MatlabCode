%analysis of data from SC_stimcue paradigm
%adapted from "FixStimAnalyses.m"

%%
% Plotting eye trajectories from fpoff to fpoff+300 ms
% This only plots the short latency saccades

stro = nex2stro;
targ_x = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_x'));
targ_y = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_y'));
uniquetargxy = unique([targ_x targ_y],'rows');
ntrials = size(stro.trial,1);
Lcatchtrials = targ_x == 0 & targ_y == 0;
laseron_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
laseroff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
fix_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpacq_t'));
targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targon_t'));
targoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targoff_t'));
sacinit_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'saccinit_t'));
sacmade_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'saccmade_t'));

Llaser = ~isnan(laseron_t);
samplerate = stro.sum.analog.storeRates{1};
figure; axes; hold on;
for i = 1:size(stro.trial,1)
    x = stro.ras{i,2}*4096/400;
    y = stro.ras{i,3}*4096/400;
    t = stro.ras{i,5}+[0:1:length(x)-1]/samplerate;
    %Lt = t>fpoff_t(i) & t < sacmade_t(i)+.1;
    Lt = t>fpoff_t(i) & t < fpoff_t(i) +.3;

    h = plot(x(Lt),y(Lt),'k');
    if Lcatchtrials(i)
        set(h,'LineStyle','--','color',[0 .5 0]);
    end
    if Llaser(i)
        set(h,'color',[1 0 0]);
    end
end
axis square
set(gca,'Xlim',[-12 12],'Ylim',[-12 12]);


%%
% % Looking at consistency of events across trials
% % First, catch trials
% Lcatchtrial = targ_x == 0 & targ_y == 0;
% figure; subplot(2,1,1); hold on;
% plotcounter = 0;
% for i = find(Lcatchtrial')
%    plot(laseron_t(i)-fpoff_t(i),plotcounter,'m*')
%    plotcounter = plotcounter + 1; 
% end
% 
% % Now saccade trials
% % green line is target on to target off
% % red star is saccade initiation
% % magenta star is laser on, magenta triangle is laser off
% % 0 is fpoff
% 
% subplot(2,1,2); hold on;
% plotcounter = 0;
% for i = find(~Lcatchtrial')
%    plot(laseron_t(i)-fpoff_t(i),plotcounter,'m*');
%    plot(laseroff_t(i)-fpoff_t(i),plotcounter,'m^');
%    plot([targon_t(i) targoff_t(i)]-fpoff_t(i),[1 1]*plotcounter,'g-');
%    plot(sacinit_t(i)-fpoff_t(i),plotcounter,'r*');
%    plotcounter = plotcounter + 1; 
% end
%%
% Plotting X,Y eye position plots on a per condition basis
H = stro.ras(:,strcmp(stro.sum.rasterCells,'AD11'));
V = stro.ras(:,strcmp(stro.sum.rasterCells,'AD12'));
ep_t = stro.ras(:,strcmp(stro.sum.rasterCells,'anlgStartTime'));
ADfreq = stro.sum.analog.storeRates{1};

for i = 1:size(uniquetargxy,1)
    Ltarg = targ_x == uniquetargxy(i,1) & targ_y == uniquetargxy(i,2);
    trialsid = find(Ltarg);
    figure; axes; hold on;
    
    for counter = 1:sum(Ltarg)
        trialid = trialsid(counter);
        h = H{trialid}*4096/400;
        v = V{trialid}*4096/400;
        t = [0:length(h)-1]./ADfreq+ep_t{trialid};
        t = t-targon_t(trialid);
        Lt = t>0.0 & t<.33;
        if Llaser(trialid)
            plot(t(Lt),h(Lt),'g-','Color',[0 1 0]);
            plot(t(Lt),v(Lt),'r-','Color',[1 0 0]);
        else
            plot(t(Lt),h(Lt),'g--','Color',[.5 1 .5]);
            plot(t(Lt),v(Lt),'r--','Color',[1 .5 .5]);
        end
    end
    title(num2str(uniquetargxy(i,:)/10));
end
% Gren is horizontal, red is vertical
% solid is laser, dashed is no laser.
%% 
% Looking at potential effects of laser on latency

CULLBASEDONLATENCY = 0; % a total hack for a grant figure. Do NOT fool yourself.
data = [];
uniquetargxy = unique([targ_x targ_y],'rows');
for i = 1:size(uniquetargxy,1)
   Ltarg = targ_x == uniquetargxy(i,1) & targ_y == uniquetargxy(i,2);
   lat_laser = sacinit_t(Ltarg&Llaser)-fpoff_t(Ltarg&Llaser);
   lat_nolaser = sacinit_t(Ltarg&~Llaser)-fpoff_t(Ltarg&~Llaser);
   if CULLBASEDONLATENCY
       lat_laser(lat_laser > .3) = nan;
       lat_nolaser(lat_nolaser > .3) = nan;
   end
   
   data(i,:) = [nanmean(lat_laser) sum(~isnan(lat_laser)) nanmean(lat_nolaser) sum(~isnan(lat_nolaser)) ];
end

figure; subplot(2,2,1); hold on;
for i = 1:size(uniquetargxy,1)
    tmp = data(i,1)-data(i,3);
    if data(i,2) > 1 & data(i,4) > 1
        if tmp == 0 | isnan(tmp)
            h = plot(uniquetargxy(i,1),uniquetargxy(i,2),'k.','MarkerSize',1);
        else
            h = plot(uniquetargxy(i,1),uniquetargxy(i,2),'ko','MarkerSize',500*abs(tmp),'Linewidth',2);
            if tmp < 0
                set(h,'Color','red');
            end
        end
    else
        h = plot(uniquetargxy(i,1),uniquetargxy(i,2),'b*','MarkerSize',1);
    end
end
axis square;
axis equal;
title('Latencies. Red = no laser latency is longer');

% Number of successful saccades
subplot(2,2,2); hold on;
for i = 1:size(uniquetargxy,1)
    if data(i,2) > 0
        h = plot(uniquetargxy(i,1),uniquetargxy(i,2),'bo','MarkerSize',5*data(i,2));
    else
        plot(uniquetargxy(i,1),uniquetargxy(i,2),'b.','MarkerSize',2);
    end
    if data(i,4) > 0
        h = plot(uniquetargxy(i,1),uniquetargxy(i,2),'k*','MarkerSize',5*data(i,4));
    else
        plot(uniquetargxy(i,1),uniquetargxy(i,2),'k.','MarkerSize',2);
    end
end
axis square;
axis equal;
