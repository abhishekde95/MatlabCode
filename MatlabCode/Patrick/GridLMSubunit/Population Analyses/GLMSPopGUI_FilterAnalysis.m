function GLMSPopGUI_FilterAnalysis(varargin)
global GLMSPopData
% This analysis is supposed to make a linear prediction

% Grab saved population data
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
load([library 'GLMSPopData.mat'])
%load([library 'GLMSPopData_0to350.mat'])
%load([library 'GLMSPopData_50to250.mat'])
%load([library 'GLMSPopData_80to280.mat'])

SetUpFig()
if ~isempty(varargin)
    for n = 1:numel(varargin)
        eval(varargin{n})
    end
end


end

function subunitonoff(~,~)
global GLMP
%Simple button for turning the subunit selction visibility on andoff.

% Load figure variables
filterfig = get(gcf,'UserData');
dnpanel = get(filterfig.dnpanel,'UserData');
conpanel = get(filterfig.conpanel,'UserData');
filterpanel = get(filterfig.filterpanel,'UserData');

if isfield(GLMP,'subunit')
    if get(conpanel.general.uicontrols.subonoff,'Value') == 1
        
        % Turn on Dense Noise subunits
        set(dnpanel.subunits,'visible','on')
        set(filterpanel.linfilter.subunits,'visible','on')
        set(filterpanel.diff.subunits,'visible','on')
       
    else
        
        % Turn all subunits off
        set(dnpanel.subunits,'Visible','off')
        set(filterpanel.linfilter.subunits,'visible','off')
        set(filterpanel.diff.subunits,'visible','off')
        
    end
end

% Save variables
set(filterfig.dnpanel,'UserData',dnpanel);
set(filterfig.conpanel,'UserData',conpanel);
set(filterfig.filterpanel,'UserData',filterpanel);
set(gcf,'UserData',filterfig);


end
 
function plotdiff(~,~)

% Load figure variables
filterfig = get(gcf,'UserData');
dnpanel = get(filterfig.dnpanel,'UserData');
conpanel = get(filterfig.conpanel,'UserData');
filterpanel = get(filterfig.filterpanel,'UserData');

val1 = get(conpanel.diff.uicontrols.stimsel1,'value');
val2 = get(conpanel.diff.uicontrols.stimsel2,'value');
STA1 = dnpanel.scaledSTAvals(val1,:);
STA2 = dnpanel.scaledSTAvals(val2,:);

temp = cat(1,STA1{:}) - cat(1,STA2{:});
minval = min(temp(:));
maxval = max(temp(:));

imgs = findall(filterpanel.diff.axes,'type','image');
delete(imgs);
for n = 1:numel(filterpanel.diff.axes)
    filterpanel.diff.vals{n} = (STA1{n}-STA2{n}-minval)./(maxval-minval);
    tempim = repmat(filterpanel.diff.vals{n},[1 1 3]);
    filterpanel.diff.diffimage(n) = image(tempim,'parent',filterpanel.diff.axes(n)); hold on;
    set(filterpanel.diff.axes(n),'Xlim',[1 size(tempim,1)],'Ylim',[1 size(tempim,2)],...
       'Xtick',[],'Ytick',[]);
   if strcmp(get(filterpanel.diff.axes(n).Children(1),'type'),'image')
       set(filterpanel.diff.axes(n),'Children',flipud(filterpanel.diff.axes(n).Children))
   end
end
           
% Save figure variables
set(filterfig.dnpanel,'UserData',dnpanel);
set(filterfig.conpanel,'UserData',conpanel);
set(filterfig.filterpanel,'UserData',filterpanel);
set(gcf,'UserData',filterfig);

end

function SetUpFig()
global GLMSPopData

% Set up figure
figure(160); clf;
set(gcf,'units','pixels','pos',[50 300 1350 700],'NumberTitle','off',...
    'Name','Population White Noise Analysis');
filterfig.dnpanel = uipanel('Pos',[.3 .475 .69 .5],'Parent',gcf,...
    'title','White Noise');
filterfig.filterpanel = uipanel('pos',[.3 .025 .69 .435],'parent',gcf,...
    'title','Filters','titleposition','centertop','fontsize',12);
filterfig.conpanel = uipanel('pos',[.01 .01 .28 .98],'parent',gcf,...
    'title','Control Panel');

% Set up control panel
conpanel.general.uicontrols.nspikes = uicontrol('style','edit','units','normalized',...
    'pos',[.6 .4 .25 .04],'parent',filterfig.conpanel,'fontsize',12,...
    'string',[]);
conpanel.general.uicontrols.nspikeslabel = uicontrol('style','text','units','normalized',...
    'pos',[.6 .45 .25 .03],'parent',filterfig.conpanel,'fontsize',12,...
    'string','# of spikes = ');
conpanel.general.uicontrols.subonoff = uicontrol('parent',filterfig.conpanel,'Style','Checkbox',...
    'units','normalized','pos',[.1 .4 .4 .05],'value',0,...
    'String','Turn On Subunit Visibility','Callback',@subunitonoff);
conpanel.general.uicontrols.reanalAll = uicontrol('parent',filterfig.conpanel,...
    'style','pushbutton','units','normalized','pos',[.1 .05 .4 .05],...
    'backgroundcolor',[1 0 0],'string','Reanalyze All','callback',@reanalAll);
conpanel.nframesback = 10;

% For difference plots
conpanel.fieldnames = {'pLpM' 'mLmM' 'pLmM' 'mLpM'};
conpanel.labels = {'+L+M' '-L-M' '+L-M' '-L+M'};
conpanel.diff.diffpanel = uipanel('parent',filterfig.conpanel,'units','normalized',...
    'pos',[.05 .15 .9 .2],'title','Difference Filter Controls');
conpanel.diff.uicontrols.stimsel1 = uicontrol('style','popupmenu','parent',conpanel.diff.diffpanel,...
    'units','normalized','pos',[.1 .3 .35 .05],'string',conpanel.labels,'fontsize',12,...
    'callback',@plotdiff);
conpanel.diff.uicontrols.stimsel2 = uicontrol('style','popupmenu','parent',conpanel.diff.diffpanel,...
    'units','normalized','pos',[.6 .3 .35 .05],'string',conpanel.labels,'fontsize',12,...
    'callback',@plotdiff);
conpanel.diff.conpanel.labels.stimsel1 = uicontrol('style','text','string','Stim Type #1',...
    'parent',conpanel.diff.diffpanel,'units','normalized','pos',[.1 .5 .35 .2],...
    'fontsize',12);
conpanel.diff.conpanel.labels.stimsel2 = uicontrol('style','text','string','Stim Type #2',...
    'parent',conpanel.diff.diffpanel,'units','normalized','pos',[.6 .5 .35 .2],...
    'fontsize',12);
conpanel.selectedidx = [];

%%% Display Table %%%
datatypes = GLMSPopData(1,:);
data = cell(size(GLMSPopData,1)-1,3);
data(:,1) = GLMSPopData(2:end,strcmp(datatypes,'Datafile'));
data(:,2) = GLMSPopData(2:end,strcmp(datatypes,'Subunit'));
data(:,3) = GLMSPopData(2:end,strcmp(datatypes,'normLLdiff'));
colname = {'Datafile' 'Sub' 'Norm LL Diff'};
colformat = cell(1,size(data,2));
colformat(:) = {'char'};
conpanel.table = uitable('parent',filterfig.conpanel,...
    'units','normalized','pos',[.01 .6  .98 .39],...
    'data',data,'columnname',colname,...
    'columnformat',colformat,...
    'CellSelectionCallback',@cellSelect);

%%% DN panel %%%
panelXpos = linspace(.05,.95,conpanel.nframesback+1);
panelYpos = linspace(.95,.05,numel(conpanel.labels)+1);
panelXpos = panelXpos(1:end-1);
panelYpos = panelYpos(2:end);
panelxsize = mean(diff(panelXpos))*.9;
panelysize = abs(mean(diff(panelYpos))*.9);

% Draw DN axes
for ypos = 1:(numel(panelYpos))
    field = conpanel.fieldnames{ypos};
    for xpos = 1:(numel(panelXpos))
        dnpanel.axes(ypos,xpos) = axes('parent',filterfig.dnpanel,'Units','normalized',...
            'pos',[panelXpos(xpos) panelYpos(ypos) panelxsize panelysize],...
            'Xtick',[],'Ytick',[],'box','on'); hold on; axis tight
        if xpos == 1
            dnpanel.labels.(field) = text(-.2,.5,(conpanel.labels{ypos}),'Rotation',90,...
                'Units','normalized','HorizontalAlignment','center');
        end
    end
end

%%% Linear Filter %%%
panelXpos = linspace(.05,.95,conpanel.nframesback+1);
panelYpos = linspace(.95,.05,4); % hard coded, only 3 analyses
panelXpos = panelXpos(1:end-1);
panelYpos = panelYpos(2:end);
panelxsize = mean(diff(panelXpos))*.9;
panelysize = abs(mean(diff(panelYpos))*.9);
    
% Draw Axes
ypos = 1;
for xpos = 1:numel(panelXpos)
    filterpanel.linfilter.axes(xpos) = axes('Parent',filterfig.filterpanel,'Units','normalized',...
        'box','on','Position',[panelXpos(xpos) panelYpos(ypos) panelxsize panelysize],...
        'xtick',[],'ytick',[]); hold on; axis tight
end

%%% Difference Panel %%%
ypos = 2;
for xpos = 1:numel(panelXpos)
    filterpanel.diff.axes(xpos) = axes('Parent',filterfig.filterpanel,'Units','normalized',...
        'box','on','Position',[panelXpos(xpos) panelYpos(ypos) panelxsize panelysize],...
        'xtick',[],'ytick',[]); hold on; axis tight
end

%%% Variance Panel %%%
ypos = 3;
filterpanel.var.axes = axes('Parent',filterfig.filterpanel,'Units','normalized',...
    'box','on','Position',[.05 panelYpos(ypos) .9 panelysize],...
    'XTick',[]); hold on; ylabel('Dev from Mean')


% Save variables
set(filterfig.dnpanel,'UserData',dnpanel);
set(filterfig.conpanel,'UserData',conpanel);
set(filterfig.filterpanel,'UserData',filterpanel);
set(gcf,'UserData',filterfig);

end

function cellSelect(a,b)
global GLMSPopData DN GLMP

% Load figure variables
filterfig = get(gcf,'UserData');
dnpanel = get(filterfig.dnpanel,'UserData');
conpanel = get(filterfig.conpanel,'UserData');
filterpanel = get(filterfig.filterpanel,'UserData');

% Set aside the index (+1 for referencing GLMSPopData)
if isnumeric(a) % a call from outside this function
    conpanel.selectedidx = a;
elseif strcmp(a.Type,'uitable') || strcmp(a.Type,'reanal')
    if isempty(b.Indices)
        return
    elseif conpanel.selectedidx == b.Indices(1)
        return
    end
    conpanel.selectedidx = b.Indices(1);
end
popidx = conpanel.selectedidx + 1;

% Unpack variables
datatypes = GLMSPopData(1,:);
DN = GLMSPopData{popidx,strcmp(datatypes,'DN')};
GLMP = GLMSPopData{popidx,strcmp(datatypes,'GLMP')};

% Fill in spike counter
set(conpanel.general.uicontrols.nspikes,'string',num2str(DN.stats.pLpM.nspikes));

% Clear out all old data
conpanel.times = 0:(1/DN.framerate):(conpanel.nframesback-1)*(1/DN.framerate);
imgs = findall(dnpanel.axes,'type','Image');
lines = findall(dnpanel.axes,'type','line');
delete(imgs)
delete(lines)
imgs = findall(filterpanel.linfilter.axes,'type','image');
lines = findall(filterpanel.linfilter.axes,'type','line');
delete(imgs)
delete(lines)
imgs = findall(filterpanel.diff.axes,'type','image');
lines = findall(filterpanel.diff.axes,'type','line');
delete(imgs)
delete(lines)
imgs = findall(filterpanel.var.axes,'type','image');
lines = findall(filterpanel.var.axes,'type','line');
delete(imgs)
delete(lines)

% Index the location of the subunits
for sub = 1:3
    if ~isempty(GLMP.subunit{sub})
        if GLMP.subunit{sub}.gridX{1} ~= 100;
            cols = GLMP.subunit{sub}.gridX{1} + ceil(DN.NStixGrid(1)/2);
            rows = -GLMP.subunit{sub}.gridY{1} + ceil(DN.NStixGrid(1)/2);
            conpanel.subxy{sub} = [cols rows];
            conpanel.subidx{sub} = sub2ind([DN.NStixGrid(1) DN.NStixGrid(1)],rows,cols);
        end
    end
end


%%% DN panel %%%
maxSTAval = nan(1,4);
minSTAval = nan(1,4);
for fn = 1:numel(conpanel.fieldnames)
    STAs = DN.stats.(conpanel.fieldnames{fn}).STA;
    maxSTAval(fn) = max(max(abs(STAs)));
    minSTAval(fn) = min(min(abs(STAs(1:DN.NStixGrid(1)^2,:))));
end
maxval = max(maxSTAval);
minval = min(minSTAval);
nstixgrid = DN.NStixGrid(1);


% Plot DN data
dnpanel.subunits = nan(4,conpanel.nframesback);
for ypos = 1:(size(dnpanel.axes,1))
    field = conpanel.fieldnames{ypos};
    im = reshape(DN.stats.(field).STA,[DN.NStixGrid(1) DN.NStixGrid(1) 3 10]);
    dnpanel.STA.(field) = squeeze(im(:,:,1,:));
    for xpos = 1:(size(dnpanel.axes,2))

        tempim = im(:,:,:,xpos);
        tempim(:,:,3) = tempim(:,:,1);
        tempim = ((abs(tempim)- minval)./(maxval-minval));
        dnpanel.scaledSTAvals{ypos,xpos} = abs(tempim(:,:,1));
        h = image(tempim,'parent',dnpanel.axes(ypos,xpos)); hold on;
        set(h,'ButtonDownFcn',@StixelSelection)
        set(dnpanel.axes(ypos,xpos),'Xtick',[],'Ytick',[]);
        if ypos == 1
            dnpanel.labels.time(xpos) = text(.5,1.1,[num2str(conpanel.times(xpos),'%.1f') 'ms'],...
                'Units','normalized','HorizontalAlignment','center');
        end
        if xpos == 1
            dnpanel.labels.(field) = text(-.2,.5,(conpanel.labels{ypos}),'Rotation',90,...
                'Units','normalized','HorizontalAlignment','center');
        end
        for sub = 1:2
            if ~isempty(GLMP.subunit{sub})
                if sub == 1
                    plotcolor = 'r';
                elseif sub == 2
                    plotcolor = 'g';
                end
                dnpanel.subunits(ypos,xpos,sub) = plot(conpanel.subxy{sub}(:,1),...
                    conpanel.subxy{sub}(:,2),[plotcolor 'o'],'parent',dnpanel.axes(ypos,xpos));
            end
        end
    end
end

%%% Linear filter %%%
linearSTA = DN.stats.pLpM.STA + DN.stats.mLmM.STA + DN.stats.pLmM.STA + DN.stats.mLpM.STA;
filterpanel.normlinearSTA = (linearSTA - min(linearSTA(:))) ./ (max(linearSTA(:))-min(linearSTA(:)));
ypos = 1;
filterpanel.linfilter.subunits = nan(conpanel.nframesback,1);
for xpos = 1:numel(filterpanel.linfilter.axes)
    
    % Plot linear filter
    tempim = reshape(filterpanel.normlinearSTA(:,xpos),[nstixgrid nstixgrid 3]);
    filterpanel.linfilter.image(ypos,xpos) = image(tempim,...
        'parent',filterpanel.linfilter.axes(ypos,xpos)); hold on;
    set(filterpanel.linfilter.axes(ypos,xpos),'Xtick',[],'Ytick',[]);
    
    % Draw subunits
    for sub = 1:2
        if ~isempty(GLMP.subunit{sub})
            if sub == 1
                plotcolor = 'r';
            elseif sub == 2
                plotcolor = 'g';
            end
            filterpanel.linfilter.subunits(xpos,sub) = plot(conpanel.subxy{sub}(:,1),...
                conpanel.subxy{sub}(:,2),[plotcolor 'o'],'parent',filterpanel.linfilter.axes(xpos));
        end
    end
end


%%% Difference Panel %%%
filterpanel.diff.subunits = nan(conpanel.nframesback,1);
for xpos = 1:numel(filterpanel.diff.axes)
    
    % Draw subunits
    if isfield(GLMP,'subunit')
        for sub = 1:2
            if ~isempty(GLMP.subunit{sub})
                if sub == 1
                    plotcolor = 'r';
                elseif sub == 2
                    plotcolor = 'g';
                end
                filterpanel.diff.subunits(xpos,sub) = plot(conpanel.subxy{sub}(:,1),...
                    conpanel.subxy{sub}(:,2),[plotcolor 'o'],'parent',filterpanel.diff.axes(xpos));
            end
        end
    end
end

%%% Deviation plot %%%
filterpanel.var.subunits = nan(conpanel.nframesback,1);
tempmat = nan(4,conpanel.nframesback,3);
for sub = 1:3
    if ~isempty(GLMP.subunit{sub})
        if GLMP.subunit{sub}.gridX{1} ~= 100
            for n = 1:conpanel.nframesback
                for m = 1:numel(conpanel.fieldnames)
                    temp = dnpanel.STA.(conpanel.fieldnames{m})(:,:,n);
                    filterpanel.var.vals.sub(sub).(conpanel.fieldnames{m})(n) = abs(mean(temp(conpanel.subidx{sub})));
                    tempmat(m,n,sub) = abs(mean(temp(conpanel.subidx{sub})));
                end
            end
        end
    end
end

% Claculate the distances
meanmat = repmat(.25,4,conpanel.nframesback,3);
diffmat = tempmat - meanmat;
for n = 1:3
    filterpanel.var.sub(n).deviation = sqrt(sum(diffmat(:,:,n).^2,1));
end

% Plot
axes(filterpanel.var.axes); cla; hold on; grid on;
plot(filterpanel.var.sub(1).deviation,'ro--')
plot(filterpanel.var.sub(2).deviation,'go--')
xlim([.5 conpanel.nframesback+.5])
ylim([0 max(cat(2,filterpanel.var.sub(1).deviation,filterpanel.var.sub(2).deviation))*1.1])


% Save variables
set(filterfig.dnpanel,'UserData',dnpanel);
set(filterfig.conpanel,'UserData',conpanel);
set(filterfig.filterpanel,'UserData',filterpanel);
set(gcf,'UserData',filterfig);

plotdiff()
subunitonoff()

end

function reanalAll(~,~)
global GLMSPopData

% Load figure variables
filterfig = get(gcf,'UserData');
dnpanel = get(filterfig.dnpanel,'UserData');
conpanel = get(filterfig.conpanel,'UserData');
filterpanel = get(filterfig.filterpanel,'UserData');

if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
datatypes = GLMSPopData(1,:);
for n = 1:(size(GLMSPopData,1)-1)

    b.Indices = n;
    a.Type = 'reanal';
    cellSelect(a,b);
    
    % Load figure variables
    filterfig = get(gcf,'UserData');
    dnpanel = get(filterfig.dnpanel,'UserData');
    conpanel = get(filterfig.conpanel,'UserData');
    filterpanel = get(filterfig.filterpanel,'UserData');

    % Pull out various filters
    filters.norm4filterSTAs = dnpanel.scaledSTAvals;
    nstix = sqrt(size(filterpanel.normlinearSTA,1)/3);
    nframes = size(filterpanel.normlinearSTA,2);
    filters.norm1filterSTA = reshape(filterpanel.normlinearSTA,[nstix nstix 3 nframes]);
    filters.subfilter = filterpanel.var.sub;
    
    % Insert filter data into pop structure
    GLMSPopData{n+1,strcmp(datatypes,'DN Filters')} = filters;

end

% Save pop data
save([library 'GLMSPopData.mat'],'GLMSPopData');


% Save variables
set(filterfig.dnpanel,'UserData',dnpanel);
set(filterfig.conpanel,'UserData',conpanel);
set(filterfig.filterpanel,'UserData',filterpanel);
set(gcf,'UserData',filterfig);

end


