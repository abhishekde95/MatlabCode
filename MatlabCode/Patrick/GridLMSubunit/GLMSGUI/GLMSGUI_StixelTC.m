function [] = GLMSGUI_StixelTC(~,~)
global DN GLMP

% Highlights statistically significant stixels and shows time courses.

% Create time course figure
figure(50); clf;
set(gcf,'units','pixels','pos',[300 300 1100 550],'NumberTitle','off',...
    'Name',['Stixel Time Course (' DN.datafile ')']);
stixtc = get(50,'UserData');
stixtc.dnpanel = uipanel('Pos',[.025 .325 .775 .65],'Parent',gcf);
stixtc.tcpanel = uipanel('pos',[.025 .025 .775 .275],'parent',gcf);
stixtc.conpanel = uipanel('pos',[.81 .025 .16 .95],'parent',gcf);
dnpanel = get(stixtc.dnpanel,'UserData');
tcpanel = get(stixtc.tcpanel,'UserData');

% Set up some variables
nframesback = 10;
maxSTAval = nan(1,4);
minSTAval = nan(1,4);
fieldnames = {'pLpM' 'mLmM' 'pLmM' 'mLpM'};
labels = {'+L+M' '-L-M' '+L-M' '-L+M'};
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
    im = reshape(DN.stats.(field).STA,[DN.NStixGrid(1) DN.NStixGrid(1) 3 10]);
    for xpos = 1:numel(panelXpos)
        dnpanel.axes(ypos,xpos) = axes('Parent',stixtc.dnpanel,'Units','pixels',...
            'Xtick',[],'Ytick',[],'box','on',...
            'Position',[panelXpos(xpos) panelYpos(ypos) panelsize panelsize]);
        tempim = im(:,:,:,xpos);
        tempim(:,:,3) = tempim(:,:,1);
        dnpanel.origSTAvals{ypos,xpos} = tempim(:,:,1);
        tempim = ((abs(tempim)- minval)./(maxval-minval));
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
                if ~isempty(GLMP.subunit{sub})
                    cols = GLMP.subunit{sub}.gridX{1} + ceil(DN.NStixGrid(1)/2);
                    rows = -GLMP.subunit{sub}.gridY{1} + ceil(DN.NStixGrid(1)/2);
                    tempsub{sub}.cols = cols;
                    tempsub{sub}.rows = rows;
                    
                    if sub == 1
                        plotcolor = 'r';
                    elseif sub == 2
                        plotcolor = 'g';
                    else
                        plotcolor = 'b';
                    end
                    k = plot(tempsub{sub}.cols,tempsub{sub}.rows,[plotcolor 'o']);
                    hold on;
                    dnpanel.subunits(ypos,xpos,sub) = k;
                    set(k,'ButtonDownFcn',@StixelSelection)
                end
            end
        end
    end
end

% Printing figures for comittee meeting
% set(gcf,'paperpositionmode','auto')
% print(gcf,'-djpeg','lumcell1b.jpeg')
% keyboard

% Display stim #1 by default
set(dnpanel.axes,'visible','off')
set(dnpanel.dnimage,'visible','off')
set(dnpanel.axes(:,:,1),'visible','on')
set(dnpanel.dnimage(:,:,1),'visible','on')


% Set up buttons
conpanel.uicontrols.nspikes = uicontrol('style','edit','units','normalized',...
    'pos',[.55 .75 .35 .05],'parent',stixtc.conpanel,'fontsize',12,...
    'string',num2str(DN.stats.pLpM.nspikes));
conpanel.uicontros.nspikeslabel = uicontrol('style','text','units','normalized',...
    'pos',[.1 .75 .4 .04],'parent',stixtc.conpanel,'fontsize',12,...
    'string','# of spikes = ');
conpanel.uicontrols.sigstixpanel = uipanel('units','normalized',...
    'pos',[.05 .45 .9 .2],'parent',stixtc.conpanel);
conpanel.uicontrols.sigstixbutton = uicontrol('style','pushbutton','units','pixels',...
    'parent',conpanel.uicontrols.sigstixpanel,'pos',[5 20 140 40],'FontSize',12,...
    'string','Show Significant Stixels','callback',@ShowSigStix);
conpanel.labels.pval = uicontrol('style','text','units','normalized',...
    'parent',conpanel.uicontrols.sigstixpanel,'pos',[.2 .65 .2 .2],...
    'string','p = ','fontsize',12);
conpanel.uicontrols.sigval = uicontrol('style','edit','units','normalized',...
    'pos',[.4 .7 .4 .2],'parent',conpanel.uicontrols.sigstixpanel,'string','.05',...
    'callback',@SigVal);
conpanel.uicontrols.subonoff = uicontrol('Style','Checkbox','String','Turn On Subunit Visibility',...
    'pos',[20 350 150 30],'parent',stixtc.conpanel,'Value',1,'Callback',@subunitonoff);
conpanel.uicontrols.TCgrp = uibuttongroup('units','normalized',...
    'parent',stixtc.conpanel,'pos',[.05 .025 .9 .4]);
conpanel.uicontrols.TC1 = uicontrol('style','radiobutton',...
    'parent',conpanel.uicontrols.TCgrp,'pos',[15 165 100 25],...
    'string','Stixel TC #1','fontsize',12);
conpanel.uicontrols.TC2 = uicontrol('style','radiobutton',...
    'parent',conpanel.uicontrols.TCgrp,'pos',[15 135 100 25],...
    'string','Stixel TC #2','fontsize',12);
conpanel.uicontrols.TC3 = uicontrol('style','radiobutton',...
    'parent',conpanel.uicontrols.TCgrp,'pos',[15 105 100 25],...
    'string','Stixel TC #3','fontsize',12);
conpanel.uicontrols.TC4 = uicontrol('style','radiobutton',...
    'parent',conpanel.uicontrols.TCgrp,'pos',[15 75 100 25],...
    'string','Stixel TC #4','fontsize',12);
conpanel.uicontrols.TC5 = uicontrol('style','radiobutton',...
    'parent',conpanel.uicontrols.TCgrp,'pos',[15 45 100 25],...
    'string','Stixel TC #5','fontsize',12);
conpanel.uicontrols.TC6 = uicontrol('style','radiobutton',...
    'parent',conpanel.uicontrols.TCgrp,'pos',[15 15 100 25],...
    'string','Stixel TC #6','fontsize',12);


% Set up time course axes
pLpMmax = max(max(abs(DN.stats.pLpM.STS)));
mLmMmax = max(max(abs(DN.stats.mLmM.STS)));
pLmMmax = max(max(abs(DN.stats.pLmM.STS)));
mLpMmax = max(max(abs(DN.stats.mLpM.STS)));
pLpMmin = min(min(abs(DN.stats.pLpM.STS(1:DN.NStixGrid(1)^2,:))));
mLmMmin = min(min(abs(DN.stats.mLmM.STS(1:DN.NStixGrid(1)^2,:))));
pLmMmin = min(min(abs(DN.stats.pLmM.STS(1:DN.NStixGrid(1)^2,:))));
mLpMmin = min(min(abs(DN.stats.mLpM.STS(1:DN.NStixGrid(1)^2,:))));
maxval = max([pLpMmax mLmMmax pLmMmax mLpMmax]);
minval = min([pLpMmin mLmMmin pLmMmin mLpMmin]);
tcpanel.axes = axes('parent',stixtc.tcpanel,'units','normalized',...
    'pos',[.075 .1 .85 .8],'YLim',[minval maxval],'XLim',[1 nframesback],...
    'XTickLabel',DN.times,'box','on'); hold on
ylabel('STS')
theoMean = DN.stats.pLpM.nspikes/4;
theo2STD = 2*sqrt((.25*.75) * DN.stats.pLpM.nspikes);
tcpanel.meanln = plot([1 nframesback],[theoMean theoMean],'k');
tcpanel.STD(1) = plot([1 nframesback],[theoMean-theo2STD theoMean-theo2STD],'k--');
tcpanel.STD(2) = plot([1 nframesback],[theoMean+theo2STD theoMean+theo2STD],'k--');

% Save variables
set(stixtc.dnpanel,'UserData',dnpanel)
set(stixtc.tcpanel,'UserData',tcpanel)
set(stixtc.conpanel,'UserData',conpanel)
set(gcf,'UserData',stixtc)

end

function subunitonoff(~,~)
global GLMP
%Simple button for turning the subunit selction visibility on andoff.

% Load figure variables
stixtc = get(50,'UserData');
dnpanel = get(stixtc.dnpanel,'UserData');
conpanel = get(stixtc.conpanel,'UserData');
tcpanel = get(stixtc.tcpanel,'UserData');

if isfield(GLMP,'subunit')
    if get(conpanel.uicontrols.subonoff,'Value') == 1
        set(dnpanel.subunits,'Visible','on')
    else
        set(dnpanel.subunits,'Visible','off')
    end
else
    disp('No subunits to be displayed!')
end

% Save variables
set(stixtc.dnpanel,'UserData',dnpanel);
set(stixtc.conpanel,'UserData',conpanel);
set(stixtc.tcpanel,'UserData',tcpanel);
set(gcf,'UserData',stixtc);


end


function StixelSelection(~,~)
global DN
% Selection of an individual stixel for time course

% Identify the tagged stixel
h = gca;
whichpt = get(h,'CurrentPoint');
whichpt = round(whichpt(1,[1 2]));
whichpt = min([whichpt; [11 11]]);

% Load figure variables
stixtc = get(50,'UserData');
dnpanel = get(stixtc.dnpanel,'UserData');
conpanel = get(stixtc.conpanel,'UserData');
tcpanel = get(stixtc.tcpanel,'UserData');
whichstim = dnpanel.axes == h;
whichstim = mod(find(whichstim),4);
if whichstim == 0
    whichstim = 4;
end
%whichcc = get(conpanel.uicontrols.stimsel,'value');

labels = {'pLpM','mLmM','pLmM','mLpM'};
skipplot = 0;

% Set up some variables particular to this function
if ~isfield(dnpanel,'selstixim')
    dnpanel.selstixim = nan(6,10); %Hard coded, make more flexible
    tcpanel.TC = nan(6,1);
end

% Set up variables particular to each TC
whichTC = get(conpanel.uicontrols.TCgrp,'selectedObject');
if whichTC == conpanel.uicontrols.TC1
    col = [1 0 0];
    whichTC = 1;
    oldPt = get(conpanel.uicontrols.TC1,'UserData');
    set(conpanel.uicontrols.TC1,'UserData',[whichstim whichpt]);
elseif whichTC == conpanel.uicontrols.TC2
    col = [0 1 0];
    whichTC = 2;
    oldPt = get(conpanel.uicontrols.TC2,'UserData');
    set(conpanel.uicontrols.TC2,'UserData',[whichstim whichpt]);
elseif whichTC == conpanel.uicontrols.TC3
    col = [0 0 1];
    whichTC = 3;
    oldPt = get(conpanel.uicontrols.TC3,'UserData');
    set(conpanel.uicontrols.TC3,'UserData',[whichstim whichpt]);
elseif whichTC == conpanel.uicontrols.TC4
    col = [1 0 1];
    whichTC = 4;
    oldPt = get(conpanel.uicontrols.TC4,'UserData');
    set(conpanel.uicontrols.TC4,'UserData',[whichstim whichpt]);
elseif whichTC == conpanel.uicontrols.TC5
    col = [1 1 0];
    whichTC = 5;
    oldPt = get(conpanel.uicontrols.TC5,'UserData');
    set(conpanel.uicontrols.TC5,'UserData',[whichstim whichpt]);
elseif whichTC == conpanel.uicontrols.TC6
    col = [0 1 1];
    whichTC = 6;
    oldPt = get(conpanel.uicontrols.TC6,'UserData');
    set(conpanel.uicontrols.TC6,'UserData',[whichstim whichpt]);
else
    disp('Problem with TC identifier')
    keyboard
end

% Draw * on selected point, delete old *, and turn off * if same point
% selected twice
if isempty(oldPt)
    for n = 1:size(dnpanel.axes,2)
        axes(dnpanel.axes(whichstim,n))
        dnpanel.selstixim(whichTC,n) = plot(whichpt(1),whichpt(2),'*','color',col,'ButtonDownFcn',@StixelSelection);
    end
elseif oldPt == [whichstim whichpt]
    if strcmp(get(dnpanel.selstixim(whichTC,:),'Visible'),'on')
        set(dnpanel.selstixim(whichTC,:),'Visible','off')
        set(tcpanel.TC(whichTC),'Visible','off')
    else
        set(dnpanel.selstixim(whichTC,:),'Visible','on')
        set(tcpanel.TC(whichTC),'Visible','on')
    end
    skipplot = 1;
else
    delete(dnpanel.selstixim(whichTC,:))
    delete(tcpanel.TC(whichTC))
    for n = 1:size(dnpanel.axes,2)
        axes(dnpanel.axes(whichstim,n))
        dnpanel.selstixim(whichTC,n) = plot(whichpt(1),whichpt(2),'*','color',col,'ButtonDownFcn',@StixelSelection);
    end
end

if skipplot == 0
    % Draw time course
    stixIdx = sub2ind([DN.NStixGrid(1) DN.NStixGrid(1)],whichpt(2),whichpt(1));
    stixelTC = DN.stats.(labels{whichstim}).STS(stixIdx,:);
    axes(tcpanel.axes);
    tcpanel.TC(whichTC) = plot(abs(stixelTC),'-','color',col);
end

% Save variables
set(stixtc.dnpanel,'UserData',dnpanel);
set(stixtc.conpanel,'UserData',conpanel);
set(stixtc.tcpanel,'UserData',tcpanel);
set(gcf,'UserData',stixtc);

end

function ShowSigStix(~,~)
global DN

% Load figure variables
stixtc = get(50,'UserData');
dnpanel = get(stixtc.dnpanel,'UserData');
conpanel = get(stixtc.conpanel,'UserData');
tcpanel = get(stixtc.tcpanel,'UserData');

if isfield(dnpanel,'sigmask') && any(any(dnpanel.sigmask))
    delete(dnpanel.sigmask)
    dnpanel.sigmask = nan(size(dnpanel.dnimage));
end

% Get stats
theoMean = .25;
pval = str2num(get(conpanel.uicontrols.sigval,'string'));
nspikes = DN.stats.pLpM.nspikes;
sigma = sqrt((theomean*(1-theomean))/nspikes);
lb = norminv(pval, theoMean, sigma);
ub = norminv(1-pval, theoMean, sigma);

% fill the 'mask' field.  1 = mask it, 0 = don't.
mask = ones(size(dnpanel.origSTAvals{1,1}));
dnpanel.sigmaskvals = cell(size(dnpanel.dnimage));
[dnpanel.sigmaskvals{:,:}] = deal(mask);
dnpanel.sigmask = nan(size(dnpanel.dnimage));
set(dnpanel.dnimage,'visible','off')
for n = 1:numel(dnpanel.dnimage)
    dnpanel.sigmaskvals{n}(abs(dnpanel.origSTAvals{n}) >= ub) = 0;
    dnpanel.sigmaskvals{n}(abs(dnpanel.origSTAvals{n}) <= lb) = 0;
    scaledSTAvals = dnpanel.scaledSTAvals{n};
    maskoverlay = dnpanel.sigmaskvals{n};
    maskoverlay(logical(dnpanel.sigmaskvals{n})) = scaledSTAvals(logical(dnpanel.sigmaskvals{n}));
    maskoverlay = reshape(maskoverlay,size(scaledSTAvals));
    maskoverlay = repmat(maskoverlay,[1 1 3]);
    scaledSTAvals = repmat(scaledSTAvals,[1 1 3]);
    combImg = scaledSTAvals - maskoverlay;
    combImg(combImg==0) = .5;
    axes(dnpanel.axes(n));
    h = image(combImg);
    dnpanel.sigmask(n) = h;
    set(h,'ButtonDownFcn',@StixelSelection)
    
    %Replot subunit
    if get(conpanel.uicontrols.subonoff,'value') == 1
        tempchildren = get(dnpanel.axes(n),'children');
        whichisdnim = tempchildren(tempchildren == dnpanel.dnimage(n));
        restofchildren = tempchildren(tempchildren ~= dnpanel.dnimage(n) & tempchildren ~= h);
        set(dnpanel.axes(n),'children',[restofchildren;h;whichisdnim])
    end
    
end

set(conpanel.uicontrols.sigstixbutton,'string','Show All Stixels',...
    'callback',@HideSigStix);


% Save variables
set(stixtc.dnpanel,'UserData',dnpanel);
set(stixtc.conpanel,'UserData',conpanel);
set(stixtc.tcpanel,'UserData',tcpanel);
set(gcf,'UserData',stixtc);

end


function HideSigStix(~,~)

% Load figure variables
stixtc = get(50,'UserData');
dnpanel = get(stixtc.dnpanel,'UserData');
conpanel = get(stixtc.conpanel,'UserData');
tcpanel = get(stixtc.tcpanel,'UserData');

delete(dnpanel.sigmask)
dnpanel.sigmask = nan(size(dnpanel.dnimage));
set(dnpanel.dnimage,'visible','on')
set(conpanel.uicontrols.sigstixbutton,'string','Show Significant Stixels',...
    'callback',@ShowSigStix);

% Save variables
set(stixtc.dnpanel,'UserData',dnpanel);
set(stixtc.conpanel,'UserData',conpanel);
set(stixtc.tcpanel,'UserData',tcpanel);
set(gcf,'UserData',stixtc);

end


function SigVal(~,~)

% Load figure variables
stixtc = get(50,'UserData');
dnpanel = get(stixtc.dnpanel,'UserData');
conpanel = get(stixtc.conpanel,'UserData');
tcpanel = get(stixtc.tcpanel,'UserData');

if strcmp(get(conpanel.uicontrols.sigstixbutton,'string'),'Show All Stixels')
    
    ShowSigStix
    
end

end


function stimsel(hObj,~)
global DN
% For selecting between different lum/col Dense Noise conditions

% Load figure variables
stixtc = get(50,'UserData');
dnpanel = get(stixtc.dnpanel,'UserData');
conpanel = get(stixtc.conpanel,'UserData');
tcpanel = get(stixtc.tcpanel,'UserData');
val = get(hObj,'Value');
nframesback = 10;

set(dnpanel.axes,'Visible','off');
set(dnpanel.dnimage,'Visible','off')
set(dnpanel.axes(:,:,val),'Visible','on')
set(dnpanel.dnimage(:,:,val),'Visible','on')
set(conpanel.uicontrols.nspikes,'string',num2str(DN.stim{val}.stats.pLpM.nspikes));

try
    if isfield(tcpanel,'TC') && any(tcpanel.TC)
        delete(dnpanel.selstixim(~isnan(dnpanel.selstixim)))
        delete(tcpanel.TC(~isnan(tcpanel.TC)))
        tcpanel.TC = nan(size(tcpanel.TC));
        dnpanel.selstixim = nan(size(dnpanel.selstixim));
        set(conpanel.uicontrols.TC1,'UserData',[])
        set(conpanel.uicontrols.TC2,'UserData',[])
        set(conpanel.uicontrols.TC3,'UserData',[])
        set(conpanel.uicontrols.TC4,'UserData',[])
        set(conpanel.uicontrols.TC5,'UserData',[])
        set(conpanel.uicontrols.TC6,'UserData',[])
    end
    delete(tcpanel.axes)
catch
    keyboard
end

% Reset time course axes
%whichcc = get(conpanel.uicontrols.stimsel,'value');
pLpMmax = max(max(abs(DN.stats.pLpM.STS)));
mLmMmax = max(max(abs(DN.stats.mLmM.STS)));
pLmMmax = max(max(abs(DN.stats.pLmM.STS)));
mLpMmax = max(max(abs(DN.stats.mLpM.STS)));
pLpMmin = min(min(abs(DN.stats.pLpM.STS(1:DN.NStixGrid(1)^2,:))));
mLmMmin = min(min(abs(DN.stats.mLmM.STS(1:DN.NStixGrid(1)^2,:))));
pLmMmin = min(min(abs(DN.stats.pLmM.STS(1:DN.NStixGrid(1)^2,:))));
mLpMmin = min(min(abs(DN.stats.mLpM.STS(1:DN.NStixGrid(1)^2,:))));
maxval = max([pLpMmax mLmMmax pLmMmax mLpMmax]);
minval = min([pLpMmin mLmMmin pLmMmin mLpMmin]);
tcpanel.axes = axes('parent',stixtc.tcpanel,'units','normalized',...
    'pos',[.075 .1 .85 .8],'YLim',[minval maxval],'XLim',[1 nframesback],...
    'XTickLabel',DN.times,'box','on'); hold on
ylabel('STS')
theoMean = DN.stats.pLpM.nspikes/4;
theo2STD = 2*sqrt(.1875 * DN.stats.pLpM.nspikes);
tcpanel.meanln = plot([1 nframesback],[theoMean theoMean],'k');
tcpanel.STD(1) = plot([1 nframesback],[theoMean-theo2STD theoMean-theo2STD],'k--');
tcpanel.STD(2) = plot([1 nframesback],[theoMean+theo2STD theoMean+theo2STD],'k--');


% Save variables
set(stixtc.dnpanel,'UserData',dnpanel);
set(stixtc.conpanel,'UserData',conpanel);
set(stixtc.tcpanel,'UserData',tcpanel);
set(gcf,'UserData',stixtc);

end

