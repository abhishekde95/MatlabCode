function GLMSGUI_DiffMaps(~,~)
global DN GLMP

% Create time course figure
figure(80); clf;
set(gcf,'units','pixels','pos',[300 300 1100 550],'NumberTitle','off',...
    'Name',['Stixel Difference Maps (' DN.datafile ')']);
difftc = get(gcf,'UserData');
difftc.dnpanel = uipanel('Pos',[.025 .325 .775 .65],'Parent',gcf);
difftc.diffpanel = uipanel('pos',[.025 .025 .775 .275],'parent',gcf);
difftc.conpanel = uipanel('pos',[.81 .025 .16 .95],'parent',gcf);
conpanel = get(difftc.conpanel,'UserData');
dnpanel = get(difftc.dnpanel,'UserData');
diffpanel = get(difftc.diffpanel,'UserData');

% Set up control panel
labels = {'pLpM'; 'mLmM';'pLmM';'mLpM'};
conpanel.uicontrols.stimsel1 = uicontrol('style','popupmenu','parent',difftc.conpanel,...
    'units','normalized','pos',[.25 .8 .5 .05],'string',labels,'fontsize',12,...
    'callback',@plotdiff);
conpanel.labels.stimsel1 = uicontrol('style','text','string','Stim Type #1',...
    'parent',difftc.conpanel,'units','normalized','pos',[.25 .85 .5 .05],...
    'fontsize',12);
conpanel.uicontrols.stimsel2 = uicontrol('style','popupmenu','parent',difftc.conpanel,...
    'units','normalized','pos',[.25 .65 .5 .05],'string',labels,'fontsize',12,...
    'callback',@plotdiff);
conpanel.labels.stimsel2 = uicontrol('style','text','string','Stim Type #2',...
    'parent',difftc.conpanel,'units','normalized','pos',[.25 .7 .5 .05],...
    'fontsize',12);
conpanel.uicontrols.subonoff = uicontrol('Style','Checkbox','String','Turn On Subunit Visibility',...
    'pos',[20 250 150 30],'parent',difftc.conpanel,'Value',1,'Callback',@subunitonoff);

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
        dnpanel.axes(ypos,xpos) = axes('Parent',difftc.dnpanel,'Units','pixels',...
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
                    else
                        plotcolor = 'b';
                    end
                    k = plot(tempsub{sub}.cols,tempsub{sub}.rows,[plotcolor 'o']);
                    dnpanel.subunits(ypos,xpos,sub) = k;
                    set(k,'ButtonDownFcn',@StixelSelection)
                end
            end
        end
    end
end

% Set up diffpanel display
nstixgrid = DN.NStixGrid(1);
tempim = zeros(nstixgrid,nstixgrid,3);
for xpos = 1:numel(panelXpos)
    diffpanel.axes(xpos) = axes('Parent',difftc.diffpanel,'Units','pixels',...
        'box','on','Position',[panelXpos(xpos) 35 panelsize panelsize]);
    h = image(tempim,'parent',diffpanel.axes(xpos)); hold on;
    diffpanel.diffimage(xpos) = h;
    set(diffpanel.axes(xpos),'Xtick',[],'Ytick',[])
    if isfield(GLMP,'subunit')
        %nsubs = numel(GLMP.subunit);
        %if nsubs == 3
        %    nsubs = 2;
        %end
        for sub = 1:2
            if isempty(GLMP.subunit{sub})
                continue
            end
            if GLMP.subunit{sub}.gridX{1} ~= 100

                cols = GLMP.subunit{sub}.gridX{1} + ceil(DN.NStixGrid(1)/2);
            
                rows = -GLMP.subunit{sub}.gridY{1} + ceil(DN.NStixGrid(1)/2);
                %Idx = sub2ind(size(tempim),rows,cols,ones(size(rows)).*sub);
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
                diffpanel.subunits(xpos,sub) = k;
            end
        end
    end
end

% Display stim #1 by default
set(dnpanel.axes,'visible','off')
set(dnpanel.dnimage,'visible','off')
set(dnpanel.axes(:,:,1),'visible','on')
set(dnpanel.dnimage(:,:,1),'visible','on')

% Save figure variables
set(difftc.dnpanel,'UserData',dnpanel);
set(difftc.conpanel,'UserData',conpanel);
set(difftc.diffpanel,'UserData',diffpanel);
set(gcf,'UserData',difftc);

end

function subunitonoff(~,~)
global GLMP
%Simple button for turning the subunit selction visibility on andoff.

% Load figure variables
difftc = get(gcf,'UserData');
dnpanel = get(difftc.dnpanel,'UserData');
conpanel = get(difftc.conpanel,'UserData');
diffpanel = get(difftc.diffpanel,'UserData');

if isfield(GLMP,'subunit')
    if get(conpanel.uicontrols.subonoff,'Value') == 1
        set(dnpanel.subunits,'Visible','on')
        set(diffpanel.subunits,'visible','on')
        for n = 1:numel(diffpanel.axes)
            try
            uistack(shiftdim(dnpanel.subunits(1,n,:)),'top')
            uistack(shiftdim(dnpanel.subunits(2,n,:)),'top')
            uistack(shiftdim(dnpanel.subunits(3,n,:)),'top')
            uistack(shiftdim(dnpanel.subunits(4,n,:)),'top')
            uistack(diffpanel.subunits(n,:),'top')
            catch
                keyboard
            end
        end
    else
        set(dnpanel.subunits,'Visible','off')
        set(diffpanel.subunits,'visible','off')
    end
else
    disp('No subunits to be displayed!')
end

% Save figure variables
set(difftc.dnpanel,'UserData',dnpanel);
set(difftc.conpanel,'UserData',conpanel);
set(difftc.diffpanel,'UserData',diffpanel);
set(gcf,'UserData',difftc);


end

function plotdiff(~,~)
global DN

% Load figure variables
difftc = get(gcf,'UserData');
dnpanel = get(difftc.dnpanel,'UserData');
conpanel = get(difftc.conpanel,'UserData');
diffpanel = get(difftc.diffpanel,'UserData');

labels = {'pLpM'; 'mLmM';'pLmM';'mLpM'};
val1 = get(conpanel.uicontrols.stimsel1,'value');
val2 = get(conpanel.uicontrols.stimsel2,'value');
nstixgrid = DN.NStixGrid(1);
%STA1 = abs(DN.stim{1}.stats.(labels{val1}).STA(1:nstixgrid^2,:));
%STA2 = abs(DN.stim{1}.stats.(labels{val2}).STA(1:nstixgrid^2,:));
STA1 = abs(DN.stats.(labels{val1}).STA(1:nstixgrid^2,:));
STA2 = abs(DN.stats.(labels{val2}).STA(1:nstixgrid^2,:));

diffSTA = STA1 - STA2;
diffSTA = diffSTA - min(min(diffSTA));
diffSTA = diffSTA ./ max(max(diffSTA));

for n = 1:numel(diffpanel.axes)
    delete(diffpanel.diffimage(n))
    tempim = repmat(reshape(diffSTA(:,n),[nstixgrid nstixgrid]),[1 1 3]);
    h = image(tempim,'parent',diffpanel.axes(n)); hold on;
    diffpanel.diffimage(n) = h;
    set(diffpanel.axes(n),'Xtick',[],'Ytick',[]);
end
           
% Save figure variables
set(difftc.dnpanel,'UserData',dnpanel);
set(difftc.conpanel,'UserData',conpanel);
set(difftc.diffpanel,'UserData',diffpanel);
set(gcf,'UserData',difftc);

end
