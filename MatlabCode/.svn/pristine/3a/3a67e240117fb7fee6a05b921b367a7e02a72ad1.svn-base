function [] = GLMSGUI_ClusterAnalysis(~,~)
global DN GLMP

% Set up figure
figure(150); clf;
set(gcf,'units','pixels','pos',[50 300 900 700],'NumberTitle','off',...
    'Name',['Cluster Analysis (' DN.datafile ')']);
clusterfig.dnpanel = uipanel('Pos',[.025 .025 .95 .5],'Parent',gcf);
clusterfig.clusterpanel = uipanel('pos',[.025 .55 .6 .425],'parent',gcf);
clusterfig.conpanel = uipanel('pos',[.65 .55 .326 .425],'parent',gcf);

% Set up control panel
conpanel.uicontrols.nspikes = uicontrol('style','edit','units','pixels',...
    'pos',[90 175 70 40],'parent',clusterfig.conpanel,'fontsize',12,...
    'string',num2str(DN.stats.pLpM.nspikes));
conpanel.uicontros.nspikeslabel = uicontrol('style','text','units','pixels',...
    'pos',[10 165 80 40],'parent',clusterfig.conpanel,'fontsize',12,...
    'string','# of spikes = ');
conpanel.uicontrols.subonoff = uicontrol('Style','Checkbox',...
    'String','Turn On Subunit Visibility','pos',[10 145 150 30],...
    'parent',clusterfig.conpanel,'Value',1,'Callback',@subunitonoff);
conpanel.nframesback = 10;

% DN panel
dnpanel.seltptim = nan(1,conpanel.nframesback);


% Calculate statistics and plot DN stimuli
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
panelXpos = linspace(35,750,conpanel.nframesback);
panelYpos = linspace(250,conpanel.nframesback,numel(labels));
MSperStix = 1000./DN.framerate;
dispTimes = linspace(0,-(conpanel.nframesback-1)*MSperStix,conpanel.nframesback);
panelsize = 75;
DN.times = dispTimes;

% Plot DN data
for ypos = 1:numel(panelYpos)
    field = fieldnames{ypos};
    im = reshape(DN.stats.(field).STA,[DN.NStixGrid(1) DN.NStixGrid(1) 3 10]);
    for xpos = 1:numel(panelXpos)
        dnpanel.axes(ypos,xpos) = axes('Parent',clusterfig.dnpanel,'Units','pixels',...
            'Xtick',[],'Ytick',[],'box','on',...
            'Position',[panelXpos(xpos) panelYpos(ypos) panelsize panelsize]);
        tempim = im(:,:,:,xpos);
        tempim(:,:,3) = tempim(:,:,1);
        dnpanel.origSTAvals{ypos,xpos} = tempim(:,:,1);
        tempim = ((abs(tempim)- minval)./(maxval-minval));% * .9;
        dnpanel.scaledSTAvals{ypos,xpos} = tempim(:,:,1);
        h = image(tempim,'parent',dnpanel.axes(ypos,xpos)); hold on;
        set(h,'ButtonDownFcn',@TimePtSelection)
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
                    Idx = sub2ind(size(tempim),rows,cols,ones(size(rows)).*sub);
                    DN.subunit{sub}.idx = Idx;
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
                    dnpanel.subunits(ypos,xpos,sub) = k;
                    set(k,'ButtonDownFcn',@TimePtSelection)
                end
            end
        end
    end
end

% Set up cluster panel
h = axes('parent',clusterfig.clusterpanel,'pos',[.15 .15 .75 .75]);
clusterpanel.axes = h;

% Save variables
set(clusterfig.dnpanel,'UserData',dnpanel)
set(clusterfig.clusterpanel,'UserData',clusterpanel)
set(clusterfig.conpanel,'UserData',conpanel)
set(gcf,'UserData',clusterfig)

end

function subunitonoff(~,~)
global GLMP
%Simple button for turning the subunit selction visibility on andoff.

% Load figure variables
clusterfig = get(gcf,'UserData');
dnpanel = get(clusterfig.dnpanel,'UserData');
conpanel = get(clusterfig.conpanel,'UserData');
clusterpanel = get(clusterfig.clusterpanel,'UserData');

if isfield(GLMP,'subunit')
    if get(conpanel.uicontrols.subonoff,'Value') == 1
        set(dnpanel.subunits,'Visible','on')
        set(clusterpanel.substix,'visible','on')
    else
        set(dnpanel.subunits,'Visible','off')
        set(clusterpanel.substix,'visible','off')
    end
else
    disp('No subunits to be displayed!')
end

% Save variables
set(clusterfig.dnpanel,'UserData',dnpanel);
set(clusterfig.conpanel,'UserData',conpanel);
set(clusterfig.clusterpanel,'UserData',clusterpanel);
set(gcf,'UserData',clusterfig);


end



function TimePtSelection(~,~)
global DN GLMP
% Selection of an time point for 3d analysis

% Identify the tagged stixel
h = gca;

% Load figure variables
clusterfig = get(gcf,'UserData');
dnpanel = get(clusterfig.dnpanel,'UserData');
conpanel = get(clusterfig.conpanel,'UserData');
clusterpanel = get(clusterfig.clusterpanel,'UserData');
whichtpt = dnpanel.axes == h;
if isempty(whichtpt)
    whichtpt = dnpanel.subunits == h;
end
[~,whichtpt] = find(whichtpt);


% Set up some variables particular to this function
if ~isfield(dnpanel,'seltptim')
end

set(dnpanel.labels.time,'backgroundcolor',[.9 .9 .9])
set(dnpanel.labels.time(whichtpt),'backgroundcolor',[1 .7 .7])

% Plot cluster
pLpMsta = abs(DN.stats.pLpM.STA(1:DN.NStixGrid(1)^2*2,whichtpt));
mLmMsta = abs(DN.stats.mLmM.STA(1:DN.NStixGrid(1)^2*2,whichtpt));
pLmMsta = abs(DN.stats.pLmM.STA(1:DN.NStixGrid(1)^2*2,whichtpt));
axes(clusterpanel.axes); cla;
clusterpanel.allstix = plot3(pLpMsta,mLmMsta,pLmMsta,'ko');
hold on; grid on;
set(clusterpanel.axes,'box','on')

for sub = 1:2
    if ~isempty(GLMP.subunit{sub})
        idx = DN.subunit{sub}.idx;
        if sub == 1
            col = '*r';
        elseif sub == 2
            col = 'g*';
        end
        h = plot3(pLpMsta(idx),mLmMsta(idx),pLmMsta(idx),col);
        clusterpanel.substix(sub) = h;
    end
end
xlabel('+L+M')
ylabel('-L-M')
zlabel('+L-M')
plot3(.25,.25,.25,'b*')
axis equal

% Save variables
set(clusterfig.dnpanel,'UserData',dnpanel);
set(clusterfig.conpanel,'UserData',conpanel);
set(clusterfig.clusterpanel,'UserData',clusterpanel);
set(gcf,'UserData',clusterfig);

end
