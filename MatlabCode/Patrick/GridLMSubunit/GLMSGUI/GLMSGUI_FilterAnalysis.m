function GLMSGUI_FilterAnalysis(~,~)
global DN GLMP
% This analysis is supposed to make a linear prediction

% Set up figure
figure(160); clf;
set(gcf,'units','pixels','pos',[50 300 1350 700],'NumberTitle','off',...
    'Name',['Filter Analysis (' DN.datafile ')']);
filterfig.dnpanel = uipanel('Pos',[.3 .475 .69 .5],'Parent',gcf,...
    'title','White Noise');
filterfig.filterpanel = uipanel('pos',[.3 .025 .69 .435],'parent',gcf,...
    'title','Filters','titleposition','centertop','fontsize',12);
filterfig.conpanel = uipanel('pos',[.01 .01 .28 .98],'parent',gcf,...
    'title','Control Panel');

% Set up control panel
conpanel.general.uicontrols.nspikes = uicontrol('style','edit','units','normalized',...
    'pos',[.55 .8 .35 .05],'parent',filterfig.conpanel,'fontsize',12,...
    'string',num2str(DN.stats.pLpM.nspikes));
conpanel.general.uicontrols.nspikeslabel = uicontrol('style','text','units','normalized',...
    'pos',[.1 .8 .4 .04],'parent',filterfig.conpanel,'fontsize',12,...
    'string','# of spikes = ');
conpanel.general.uicontrols.subonoff = uicontrol('parent',filterfig.conpanel,'Style','Checkbox',...
    'units','normalized','pos',[.1 .7 .4 .05],'value',1,...
    'String','Turn On Subunit Visibility','Callback',@subunitonoff);
conpanel.nframesback = 10;
conpanel.times = 0:(1/DN.framerate):(conpanel.nframesback-1)*(1/DN.framerate);

% For difference plots
conpanel.fieldnames = {'pLpM' 'mLmM' 'pLmM' 'mLpM'};
conpanel.labels = {'+L+M' '-L-M' '+L-M' '-L+M'};
conpanel.diff.diffpanel = uipanel('parent',filterfig.conpanel,'units','normalized',...
    'pos',[.05 .05 .9 .2],'title','Difference Filter Controls');
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
for sub = 1:2
    if ~isempty(GLMP.subunit{sub})
        cols = GLMP.subunit{sub}.gridX{1} + ceil(DN.NStixGrid(1)/2);
        rows = -GLMP.subunit{sub}.gridY{1} + ceil(DN.NStixGrid(1)/2);
        conpanel.subxy{sub} = [cols rows];
        conpanel.subidx{sub} = sub2ind([DN.NStixGrid(1) DN.NStixGrid(1)],rows,cols);
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
nstixgrid = DN.NStixGrid(1);

% Set up figure
maxval = max(maxSTAval);
minval = min(minSTAval);
panelXpos = linspace(.05,.95,conpanel.nframesback+1);
panelYpos = linspace(.95,.05,numel(conpanel.labels)+1);
panelXpos = panelXpos(1:end-1);
panelYpos = panelYpos(2:end);
panelxsize = mean(diff(panelXpos))*.9;
panelysize = abs(mean(diff(panelYpos))*.9);

% Plot DN data
for ypos = 1:(numel(panelYpos))
    field = conpanel.fieldnames{ypos};
    im = reshape(DN.stats.(field).STA,[DN.NStixGrid(1) DN.NStixGrid(1) 3 10]);
    dnpanel.STA.(field) = squeeze(im(:,:,1,:));
    for xpos = 1:(numel(panelXpos))
        dnpanel.axes(ypos,xpos) = axes('parent',filterfig.dnpanel,'Units','normalized',...
            'pos',[panelXpos(xpos) panelYpos(ypos) panelxsize panelysize],...
            'Xtick',[],'Ytick',[],'box','on');
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
                k = plot(conpanel.subxy{sub}(:,1),conpanel.subxy{sub}(:,2),[plotcolor 'o']);
                dnpanel.subunits(ypos,xpos,sub) = k;
                set(k,'ButtonDownFcn',@StixelSelection)
            end
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

linearSTA = DN.stats.pLpM.STA + DN.stats.mLmM.STA + DN.stats.pLmM.STA + DN.stats.mLpM.STA;
filterpanel.normlinearSTA = (linearSTA - min(linearSTA(:))) ./ (max(linearSTA(:))-min(linearSTA(:)));
ypos = 1;
for xpos = 1:numel(panelXpos)
    
    % Draw Axes
    filterpanel.linfilter.axes(xpos) = axes('Parent',filterfig.filterpanel,'Units','normalized',...
        'box','on','Position',[panelXpos(xpos) panelYpos(ypos) panelxsize panelysize]);
    
    % Plot All Stim Filter
    tempim = reshape(filterpanel.normlinearSTA(:,xpos),[nstixgrid nstixgrid 3]);
    h = image(tempim,'parent',filterpanel.linfilter.axes(ypos,xpos)); hold on;
    filterpanel.linfilter.image(ypos,xpos) = h;
    set(filterpanel.linfilter.axes(ypos,xpos),'Xtick',[],'Ytick',[]);
    
    % Draw subunits
    if isfield(GLMP,'subunit')
        for sub = 1:2
            if ~isempty(GLMP.subunit{sub})
                if sub == 1
                    plotcolor = 'r';
                elseif sub == 2
                    plotcolor = 'g';
                end
                k = plot(conpanel.subxy{sub}(:,1),conpanel.subxy{sub}(:,2),[plotcolor 'o']);
                filterpanel.linfilter.subunits(xpos,sub) = k;
            end
        end
    end
end

%%% Difference Panel %%%
ypos = 2;
for xpos = 1:numel(panelXpos)
    
    % Draw Axes
    filterpanel.diff.axes(xpos) = axes('Parent',filterfig.filterpanel,'Units','normalized',...
        'box','on','Position',[panelXpos(xpos) panelYpos(ypos) panelxsize panelysize],...
        'ydir','reverse'); hold on;
    
    % Draw subunits
    if isfield(GLMP,'subunit')
        for sub = 1:2
            if ~isempty(GLMP.subunit{sub})
                if sub == 1
                    plotcolor = 'r';
                elseif sub == 2
                    plotcolor = 'g';
                end
                k = plot(conpanel.subxy{sub}(:,1),conpanel.subxy{sub}(:,2),[plotcolor 'o']);
                filterpanel.diff.subunits(xpos,sub) = k;
            end
        end
    end
end


%%% Variance Panel %%%
ypos = 3;

% Draw Axes
filterpanel.var.axes = axes('Parent',filterfig.filterpanel,'Units','normalized',...
    'box','on','Position',[.05 panelYpos(ypos) .9 panelysize],...
    'XTick',[]); hold on; ylabel('Dev from Mean')
    

% Save variables
set(filterfig.dnpanel,'UserData',dnpanel);
set(filterfig.conpanel,'UserData',conpanel);
set(filterfig.filterpanel,'UserData',filterpanel);
set(gcf,'UserData',filterfig);

plotdiff()
plotVar()

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
else
    disp('No subunits to be displayed!')
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

function plotVar(~,~)

% Load figure variables
filterfig = get(gcf,'UserData');
dnpanel = get(filterfig.dnpanel,'UserData');
conpanel = get(filterfig.conpanel,'UserData');
filterpanel = get(filterfig.filterpanel,'UserData');

tempmat = nan(4,conpanel.nframesback,2);
for sub = 1:2
    for n = 1:conpanel.nframesback
        for m = 1:numel(conpanel.fieldnames)
            temp = dnpanel.STA.(conpanel.fieldnames{m})(:,:,n);
            filterpanel.var.vals.sub(sub).(conpanel.fieldnames{m})(n) = abs(mean(temp(conpanel.subidx{sub})));
            tempmat(m,n,sub) = abs(mean(temp(conpanel.subidx{sub})));
        end
    end
end

meanmat = repmat(.25,4,conpanel.nframesback,2);
diffmat = tempmat - meanmat;
filterpanel.var.sub(1).deviation = sqrt(sum(diffmat(:,:,1).^2,1));
filterpanel.var.sub(2).deviation = sqrt(sum(diffmat(:,:,2).^2,1));

axes(filterpanel.var.axes); cla; hold on; grid on;
plot(filterpanel.var.sub(1).deviation,'ro--')
plot(filterpanel.var.sub(2).deviation,'go--')
xlim([.5 conpanel.nframesback+.5])
ylim([0 max(cat(2,filterpanel.var.sub(1).deviation,filterpanel.var.sub(2).deviation))*1.1])

% Save figure variables
set(filterfig.dnpanel,'UserData',dnpanel);
set(filterfig.conpanel,'UserData',conpanel);
set(filterfig.filterpanel,'UserData',filterpanel);
set(gcf,'UserData',filterfig);

end
