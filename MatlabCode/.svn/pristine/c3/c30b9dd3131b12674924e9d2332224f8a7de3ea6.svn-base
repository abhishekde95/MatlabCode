function GLMSGUI_Overview(dn,glmp)
global DN GLMP
DN = dn;
GLMP = glmp;

%% Plot DN

figure(279); clf;
set(gcf,'pos',[100 100 1200 725],'numbertitle','off','name',['Overview: ' GLMP.datafile])
overviewfig.dnpanel = uipanel(gcf,'units','normalized','pos',[.275 .5 .715 .49]);
overviewfig.glmppanel = uipanel(gcf,'units','normalized','pos',[.275 .01 .715 .48]);
overviewfig.conpanel = uipanel(gcf,'units','normalized','pos',[.01 .01 .255 .24]);
overviewfig.PSTHpanel = uipanel(gcf,'units','normalized','pos',[.01 .26 .255 .73]);

conpanel.uicontrols.subonoff = uicontrol('parent',overviewfig.conpanel,...
    'Style','Checkbox','units','normalized','pos',[.1 .6 .8 .2],'Value',1,...
    'String','Turn On Subunit Visibility','Callback',@subunitonoff);
conpanel.uicontrols.loadconpan = uicontrol('parent',overviewfig.conpanel,...
    'Style','Pushbutton','units','normalized','pos',[.1 .3 .8 .2],...
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

glmppanel.axes.sub1 = axes('Parent',overviewfig.glmppanel,'Units','normalized',...
    'OuterPosition',[.05 .05 .3 .9],...
    'XGrid','on','YGrid','on','ZGrid','on','box','on');
    %'XLim',[-maxcc maxcc],'YLim',[-maxcc maxcc],'ZLim',[0 25],...
xlabel('Lcc'); ylabel('Mcc'); hold on;
title('Subunit 1','FontSize',12,'FontWeight','bold')
if ~isempty(GLMP.subunit{1})
    Lcc = GLMP.subunit{1}.uniqueLcc;
    Mcc = GLMP.subunit{1}.uniqueMcc;
    nsp = GLMP.subunit{1}.meannspikes;
    [maxnsp,maxidx] = max(nsp);
    [minnsp,minidx] = min(nsp);
    for i = 1:size(Lcc,1)
        h(i) = plot(Lcc(i),Mcc(i),'ko');
        set(h(i),'MarkerFaceColor','black','MarkerSize',(((nsp(i)-minnsp)/maxnsp)+.2).*30,'MarkerEdgeColor','white')
    end
    axis square equal tight;
    set(gca,'xtick',[-.4 -.2 0 .2 .4],'ytick',[-.4 -.2 0 .2 .4])
    legend([h(maxidx),h(minidx)],num2str(maxnsp),num2str(minnsp),'location','northwest')
end


glmppanel.axes.sub2 = axes('Parent',overviewfig.glmppanel,'Units','normalized',...
    'OuterPosition',[.35 .05 .3 .9],...
    'XGrid','on','YGrid','on','ZGrid','on','box','on');
xlabel('Lcc'); ylabel('Mcc'); hold on;
title('Subunit 2','FontSize',12,'FontWeight','bold')
if ~isempty(GLMP.subunit{2})
    Lcc = GLMP.subunit{2}.uniqueLcc;
    Mcc = GLMP.subunit{2}.uniqueMcc;
    nsp = GLMP.subunit{2}.meannspikes;
    [maxnsp,maxidx] = max(nsp);
    [minnsp,minidx] = min(nsp);
    for i = 1:size(Lcc,1)
        h(i) = plot(Lcc(i),Mcc(i),'ko');
        set(h(i),'MarkerFaceColor','black','MarkerSize',(((nsp(i)-minnsp)/maxnsp)+.2).*30,'MarkerEdgeColor','white')
    end
    axis square equal tight;
    set(gca,'xtick',[-.4 -.2 0 .2 .4],'ytick',[-.4 -.2 0 .2 .4])
    legend([h(maxidx),h(minidx)],num2str(maxnsp),num2str(minnsp),'location','northwest')
else
    set(glmppanel.axes.sub2,'color',[.3 .3 .3])
end

glmppanel.axes.sub3 = axes('Parent',overviewfig.glmppanel,'Units','normalized',...
    'OuterPosition',[.65 .05 .3 .9],...
    'XLim',[-maxcc maxcc],'YLim',[-maxcc maxcc],'ZLim',[0 25],...
    'XGrid','on','YGrid','on','ZGrid','on','box','on'); hold on;
xlabel('Lcc'); ylabel('Mcc')
if ~isempty(GLMP.subunit{3})
    if GLMP.subunit{3}.gridX{1} == 100
        title('Gabors','FontSize',12,'FontWeight','bold')
    else
        title('Subunit 2','Fontsize',12,'FontWeight','bold')
    end
    Lcc = GLMP.subunit{3}.uniqueLcc;
    Mcc = GLMP.subunit{3}.uniqueMcc;
    nsp = GLMP.subunit{3}.meannspikes;
    [maxnsp,maxidx] = max(nsp);
    [minnsp,minidx] = min(nsp);
    for i = 1:size(Lcc,1)
        h(i) = plot(Lcc(i),Mcc(i),'ko');
        set(h(i),'MarkerFaceColor','black','MarkerSize',(((nsp(i)-minnsp)/maxnsp)+.2).*30,'MarkerEdgeColor','white')
    end
    axis square equal tight;
    set(gca,'xtick',[-.4 -.2 0 .2 .4],'ytick',[-.4 -.2 0 .2 .4])
    legend([h(maxidx),h(minidx)],num2str(maxnsp),num2str(minnsp),'location','northwest')
else
    set(glmppanel.axes.sub3,'color',[.3 .3 .3])
end

%% PSTH

binsize = .02; % in seconds
bins = 0:binsize:.5;
PSTH = histc(cat(1,GLMP.normspiketimes{:}),bins);
PSTH = (PSTH./numel(GLMP.normspiketimes)) ./ binsize;

psthpanel.axes = axes('parent',overviewfig.PSTHpanel,'units','normalized',...
    'pos',[.1 .1 .8 .8]);
psthpanel.psth = bar(bins,PSTH,'facecolor',[.5 0 .9],'edgecolor','k');
%psthpanel.psth = plot(bins,PSTH,'k-');
set(gca,'xlim',[min(bins) max(bins)])
title('PSTH')
xlabel('Time from Stim Onset (ms)')
ylabel('Mean Firing Rate (sp/s)')


%% Save data

set(overviewfig.conpanel,'userdata',conpanel)
set(overviewfig.dnpanel,'userdata',dnpanel)
set(overviewfig.glmppanel,'userdata',glmppanel)
set(overviewfig.PSTHpanel,'userdata',psthpanel)
set(gcf,'userdata',overviewfig)

% Printing for sfn poster
name = strcat('/Users/jpatrickweller/Documents/Matlab Working Directory/Patrick/GridLMSubunit/Figures/',GLMP.datafile,'overview');
export_fig((name),'-eps','-depsc','-transparent')


end

function subunitonoff(~,~)
%Simple button for turning the subunit selction visibility on andoff.

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

