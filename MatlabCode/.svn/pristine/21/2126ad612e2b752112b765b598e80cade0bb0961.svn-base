                  function GLMSGUI_Overview(~,~)
global DN GLMP
%DN = dn;
%GLMP = glmp;

%% Plot DN
figure(279); clf;
set(gcf,'pos',[100 200 1200 725],'numbertitle','off','name',['Overview: ' GLMP.datafile])
overviewfig.dnpanel = uipanel(gcf,'units','normalized','pos',[.275 .5 .715 .49]);
overviewfig.glmppanel = uipanel(gcf,'units','normalized','pos',[.275 .01 .715 .48]);
overviewfig.conpanel = uipanel(gcf,'units','normalized','pos',[.01 .01 .255 .24]);
overviewfig.PSTHpanel = uipanel(gcf,'units','normalized','pos',[.01 .26 .255 .73]);

conpanel.uicontrols.subonoff = uicontrol('parent',overviewfig.conpanel,...
    'Style','Checkbox','units','normalized','pos',[.1 .7 .8 .1],'Value',1,...
    'String','Subunit Visible','Callback',@subunitonoff);
conpanel.uicontrols.gaussianonoff = uicontrol('parent',overviewfig.conpanel,...
    'Style','Checkbox','units','normalized','pos',[.1 .5 .8 .1],'Value',1,...
    'String','Counting Window Visibile','Callback',@GaussianOnOff);
conpanel.uicontrols.loadconpan = uicontrol('parent',overviewfig.conpanel,...
    'Style','Pushbutton','units','normalized','pos',[.1 .1 .8 .2],...
    'String','Load Control Panel','Callback',@LoadConPanel);

% Set up some variables
nframesback = 10;
maxSTAval = nan(1,4);
minSTAval = nan(1,4);
fieldnames = {'pLpM' 'mLmM' 'pLmM' 'mLpM'};
labels = {'+L+M' '-L-M' '+L-M' '-L+M'};

% Pull out STAs and highest/lowest values to scale
for fn = 1:numel(fieldnames)
    STAs = DN.stats.(fieldnames{fn}).STA;
    maxSTAval(fn) = max(max(abs(STAs)));
    minSTAval(fn) = min(min(abs(STAs(1:DN.NStixGrid(1)^2,:))));
end
        
% Set up figure
maxval = max(maxSTAval);
minval = min(minSTAval);
panelXpos = linspace(35,750,nframesback);
panelYpos = linspace(250,nframesback,numel(labels));
MSperStix = 1000./DN.framerate;
dispTimes = linspace(0,-(nframesback-1)*MSperStix,nframesback);
panelsize = 75;
DN.times = dispTimes;


% Plot DN data
for ypos = 1:numel(panelYpos)
    field = fieldnames{ypos};
    im = reshape(DN.stats.(field).STA,[DN.NStixGrid(1) DN.NStixGrid(1) 3 nframesback]);
    for xpos = 1:numel(panelXpos)
        dnpanel.axes(ypos,xpos) = axes('Parent',overviewfig.dnpanel,'Units','pixels',...
            'Xtick',[],'Ytick',[],'box','on',...
            'Position',[panelXpos(xpos) panelYpos(ypos) panelsize panelsize]);
        tempim = im(:,:,:,xpos);
        tempim(:,:,3) = tempim(:,:,1);
        dnpanel.origSTAvals{ypos,xpos} = tempim(:,:,1);
        tempim = ((abs(tempim)- minval)./(maxval-minval));% * .9;
        dnpanel.scaledSTAvals{ypos,xpos} = tempim(:,:,1);
        h = image(tempim,'parent',dnpanel.axes(ypos,xpos)); hold on;
        set(h,'ButtonDownFcn',@StixelSelection)
        dnpanel.dnimage(ypos,xpos) = h;
        set(dnpanel.axes(ypos,xpos),'Xtick',[],'Ytick',[]);
        if ypos == 1
            dnpanel.labels.time(xpos) = text(.5,1.1,[num2str(dispTimes(xpos),'%.1f') 'ms'],'Units','normalized',...
                'HorizontalAlignment','center');
        end
        if xpos == 1
            dnpanel.labels.(field) = text(-.2,.5,(labels{ypos}),'Rotation',90,'Units','normalized',...
                'HorizontalAlignment','center');
        end
        if isfield(GLMP,'subunit')
            for sub = 1:2
                if isempty(GLMP.subunit{sub})
                    continue
                end
                cols = GLMP.subunit{sub}.gridX{1} + ceil(DN.NStixGrid(1)/2);
                rows = -GLMP.subunit{sub}.gridY{1} + ceil(DN.NStixGrid(1)/2);
                if GLMP.subunit{sub}.gridX{1} ~= 100
                    Idx = sub2ind(size(tempim),rows,cols,ones(size(rows)).*sub);
                    tempsub{sub}.cols = cols;
                    tempsub{sub}.rows = rows;
                    
                    if sub == 1
                        plotcolor = 'r';
                    elseif sub == 2
                        plotcolor = 'g';
                    end
                    k = plot(tempsub{sub}.cols,tempsub{sub}.rows,[plotcolor 'o']);
                    dnpanel.subunits(ypos,xpos,sub) = k;
                    set(k,'ButtonDownFcn',@StixelSelection)
                end
            end
        end
    end
end

%% Plot GLMP
maxcc = max(GLMP.rho);
scalefac = 30;

% Set up axes
glmppanel.axes.sub(1) = axes('Parent',overviewfig.glmppanel,'Units','normalized',...
    'OuterPosition',[.05 .05 .3 .9],...
    'XGrid','on','YGrid','on','ZGrid','on','box','on');hold on;
xlabel('Lcc'); ylabel('Mcc'); zlabel('Firing Rate (sp/s)')
title('Subunit 1','FontSize',12,'FontWeight','bold')

glmppanel.axes.sub(2) = axes('Parent',overviewfig.glmppanel,'Units','normalized',...
    'OuterPosition',[.35 .05 .3 .9],...
    'XGrid','on','YGrid','on','ZGrid','on','box','on'); hold on;
xlabel('Lcc'); ylabel('Mcc'); zlabel('Firing Rate (sp/s)') 
title('Subunit 2','FontSize',12,'FontWeight','bold')

glmppanel.axes.sub(3) = axes('Parent',overviewfig.glmppanel,'Units','normalized',...
    'OuterPosition',[.65 .05 .3 .9],...
    'XLim',[-maxcc maxcc],'YLim',[-maxcc maxcc],'ZLim',[0 25],...
    'XGrid','on','YGrid','on','ZGrid','on','box','on'); hold on;
xlabel('Lcc'); ylabel('Mcc'); zlabel('Firing Rate (sp/s)')

for n = 1:3
    axes(glmppanel.axes.sub(n))
    if ~isempty(GLMP.subunit{n})
        Lcc = GLMP.subunit{n}.uniqueLcc;
        Mcc = GLMP.subunit{n}.uniqueMcc;
        %nsp = GLMP.subunit{n}.meannspikes;
        nsp = GLMP.subunit{n}.meanfr;
        [maxnsp,maxidx] = max(nsp);
        [minnsp,minidx] = min(nsp);
        for i = 1:size(Lcc,1)
            h(i) = plot3(Lcc(i),Mcc(i),nsp(i),'ko');
            set(h(i),'MarkerFaceColor','black','MarkerSize',(((nsp(i)-minnsp)/maxnsp)+.2).*scalefac,'MarkerEdgeColor','white')
        end
        axis square tight;
        set(gca,'xtick',[-.4 -.2 0 .2 .4],'ytick',[-.4 -.2 0 .2 .4])
        legend([h(maxidx),h(minidx)],num2str(maxnsp),num2str(minnsp),'location','northwest')
        if GLMP.subunit{n}.gridX{1} == 100
            title('Gabors','FontSize',12,'FontWeight','bold')
        end
    else
        set(glmppanel.axes.sub(n),'color',[.3 .3 .3])
    end

end



%% PSTH

% % Set up a few variables
% binsize = .01; % in seconds
% bins = -.1:binsize:.5;
% PSTH = histc(cat(1,GLMP.normspiketimes{GLMP.flashL}),bins);
% PSTH = (PSTH./sum(GLMP.flashL)) ./ binsize;
% psthpanel.axes = axes('parent',overviewfig.PSTHpanel,'units','normalized',...
%     'pos',[.1 .1 .8 .8]); cla; hold on; box on;
% % 
% % plot hist
% psthpanel.psth = bar(bins,PSTH,'facecolor',[.5 0 .9],'edgecolor','k');
% set(gca,'xlim',[min(bins) max(bins)])
% title('PSTH')
% xlabel('Time from Stim Onset (ms)')
% ylabel('Mean Firing Rate (sp/s)')
% 
% % Set up gaussian filter
% whichbins = bins >= GLMP.countingwin(1) & bins <= GLMP.countingwin(2);
% gaussfilter = gausswin(sum(whichbins),1.5); % temportal form is a gaussian
% gf = (gaussfilter-min(gaussfilter)) ./ (max(gaussfilter)-min(gaussfilter)) .* max(PSTH);
% %psthpanel.tempwin = plot(bins(whichbins),gf,'r--');
% plot([GLMP.countingwin(1) GLMP.countingwin(1)],[0 max(gf)],'r--')
% plot([GLMP.countingwin(2) GLMP.countingwin(2)],[0 max(gf)],'r--')

%% Alternate PSTH

% Bin and smooth baseline 
binsize = .001; % in seconds
bins = -.2:binsize:.5;
gaussSize = .01;
gaussSize = ceil(gaussSize/binsize);
gaussfilter = gausswin(gaussSize,3);
gaussfilter = gaussfilter./sum(gaussfilter); %normalizing
PSTH = histc(cat(1,GLMP.normspiketimes{GLMP.flashL}),bins);
PSTH = (PSTH./sum(GLMP.flashL)) ./ binsize;
smoothPSTH = conv(PSTH,gaussfilter,'same');
nrows = sum(GLMP.flashL);
rowcoords = linspace(max(smoothPSTH),0,nrows*2+1);
rowcents = rowcoords(2:2:end);
rowspace = rowcoords(2)-rowcoords(1);

% plot hist
psthpanel.axes = axes('parent',overviewfig.PSTHpanel,'units','normalized',...
    'pos',[.1 .1 .8 .8]); cla; hold on; box on;
plot(bins,smoothPSTH,'linestyle','-','color',[.5 .5 .5])
plot([bins(1) bins(end)],[GLMP.blfrthresh GLMP.blfrthresh],'r--')
plot([GLMP.countingwin(1) GLMP.countingwin(1)],[0 max(smoothPSTH)],'r--')
plot([GLMP.countingwin(2) GLMP.countingwin(2)],[0 max(smoothPSTH)],'r--')
for n = 1:nrows
    tpts = GLMP.normspiketimes{n};
    if ~isempty(tpts)
        plot(repmat([tpts tpts],numel(tpts),1),repmat([rowcents(n)-rowspace rowcents(n)+rowspace],numel(tpts,1)),'k-')
    end
end
xlim([bins(1) bins(end)])
ylim([0 max(smoothPSTH)])
set(gca,'xlim',[min(bins) max(bins)])
title('PSTH')
xlabel('Time from Stim Onset (ms)')
ylabel('Mean Firing Rate (sp/s)')

% Save data
set(overviewfig.conpanel,'userdata',conpanel)
set(overviewfig.dnpanel,'userdata',dnpanel)
set(overviewfig.glmppanel,'userdata',glmppanel)
set(overviewfig.PSTHpanel,'userdata',psthpanel)
set(gcf,'userdata',overviewfig)

% For printing figures
%name = strcat('/Users/jpatrickweller/Documents/Matlab Working Directory/Patrick/GridLMSubunit/Figures/',GLMP.datafile,'overview1');
%name = strcat('/Users/jpatrickweller/Documents/Matlab Working Directory/Patrick/GridLMSubunit/Figures/diffnormLLscat');
%figure(1); clf;
%s = copyobj(gca,1)
%export_fig((name),'-eps','-depsc','-transparent')
%export_fig((name),'-eps','-depsc','-transparent')
%export_fig((name),'-tif','-transparent')


end

function GaussianOnOff(~,~)
%Simple button for turning the subunit selction visibility on andoff.

% Load figure variables
overviewfig = get(gcf,'UserData');
dnpanel = get(overviewfig.dnpanel,'UserData');
conpanel = get(overviewfig.conpanel,'UserData');
psthpanel = get(overviewfig.PSTHpanel,'userdata');

if get(conpanel.uicontrols.gaussianonoff,'Value') == 1
    set(psthpanel.tempwin,'Visible','on')
else
    set(psthpanel.tempwin,'Visible','off')
end

% Save figure variables
set(overviewfig.dnpanel,'UserData',dnpanel);
set(overviewfig.conpanel,'UserData',conpanel);
set(overviewfig.PSTHpanel,'userdata',psthpanel);

set(gcf,'UserData',overviewfig);


end

function subunitonoff(~,~)
%Simple button for turning gaussian counting window visibility on and off.

% Load figure variables
overviewfig = get(gcf,'UserData');
dnpanel = get(overviewfig.dnpanel,'UserData');
conpanel = get(overviewfig.conpanel,'UserData');

if get(conpanel.uicontrols.subonoff,'Value') == 1
    set(dnpanel.subunits,'Visible','on')
else
    set(dnpanel.subunits,'Visible','off')
end

% Save figure variables
set(overviewfig.dnpanel,'UserData',dnpanel);
set(overviewfig.conpanel,'UserData',conpanel);
set(gcf,'UserData',overviewfig);


end


function LoadConPanel(~,~)
global GLMP DN

% Analysis GUI (interactive data display and analysis)
GLMS_AnalysisGUI(GLMP,DN);

end

