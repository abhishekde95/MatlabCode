function GLMSGUI_LinPred(~,~)
global DN GLMP

% This analysis is supposed to make a linear prediction


% Set up figure
figure(100); clf;
set(gcf,'units','pixels','pos',[50 300 1350 700],'NumberTitle','off',...
    'Name',['Linear Analysis (' DN.datafile ')']);
filterfig = get(gcf,'UserData');
filterfig.dnpanel = uipanel('Pos',[.175 .45 .625 .525],'Parent',gcf);
filterfig.filterpanel = uipanel('pos',[.175 .025 .625 .4],'parent',gcf);
filterfig.conpanel = uipanel('pos',[.81 .025 .16 .95],'parent',gcf);
filterfig.specpanel = uipanel('pos',[.015 .025 .15 .95],'parent',gcf);
dnpanel = get(filterfig.dnpanel,'UserData');
filterpanel = get(filterfig.filterpanel,'UserData');
specpanel = get(filterfig.specpanel,'UserData');

% Set up some variables
nframesback = 10;
allStimCond = cat(2,DN.lumCC,DN.colCC);
uniqueStim = unique(allStimCond,'rows');

% Calculate statistics and plot DN stimuli
for u = 1:size(uniqueStim,1)    

    maxSTAval = nan(1,4);
    minSTAval = nan(1,4);
    fieldnames = {'pLpM' 'mLmM' 'pLmM' 'mLpM'};
    labels = {'+L+M' '-L-M' '+L-M' '-L+M'};
    for fn = 1:numel(fieldnames)
        STS = DN.stim{u}.stats.(fieldnames{fn}).STS;
        n = DN.stim{u}.stats.(fieldnames{fn}).nspikes;
        STAs = DN.stim{u}.stats.(fieldnames{fn}).STA;
        maxSTAval(fn) = max(max(abs(STAs)));
        minSTAval(fn) = min(min(abs(STAs(1:DN.NStixGrid(1)^2,:))));
    end
        
    % Set up figure
    if u == 1
        maxval = max(maxSTAval);
        minval = min(minSTAval);
        panelXpos = linspace(35,750,nframesback);
        panelYpos = linspace(250,nframesback,numel(labels));
        MSperStix = 1000./DN.framerate;
        dispTimes = linspace(0,-(nframesback-1)*MSperStix,nframesback);
        panelsize = 75;
        DN.times = dispTimes;
    else
        set(dnpanel.axes(:,:,u-1),'Visible','off')
    end

    % Plot DN data
    for ypos = 1:numel(panelYpos)
        field = fieldnames{ypos};
        im = reshape(DN.stim{u}.stats.(field).STA,[DN.NStixGrid(1) DN.NStixGrid(1) 3 10]);
        for xpos = 1:numel(panelXpos)
            
            dnpanel.axes(ypos,xpos,u) = axes('Parent',filterfig.dnpanel,'Units','pixels',...
                'Xtick',[],'Ytick',[],'box','on',...
                'Position',[panelXpos(xpos) panelYpos(ypos) panelsize panelsize]);
            tempim = im(:,:,:,xpos);
            tempim(:,:,3) = tempim(:,:,1);
            dnpanel.origSTAvals{ypos,xpos,u} = tempim(:,:,1);
            tempim = ((abs(tempim)- minval)./(maxval-minval));% * .9;
            dnpanel.scaledSTAvals{ypos,xpos} = tempim(:,:,1);
            h = image(tempim,'parent',dnpanel.axes(ypos,xpos,u)); hold on;
            set(h,'ButtonDownFcn',@StixelSelection)
            dnpanel.dnimage(ypos,xpos,u) = h;
            set(dnpanel.axes(ypos,xpos,u),'Xtick',[],'Ytick',[]);
            if ypos == 1
                dnpanel.labels.time(xpos) = text(.5,1.1,[num2str(dispTimes(xpos),'%.1f') 'ms'],'Units','normalized',...
                    'HorizontalAlignment','center');
            end
            if xpos == 1
                dnpanel.labels.(field) = text(-.2,.5,(labels{ypos}),'Rotation',90,'Units','normalized',...
                    'HorizontalAlignment','center');
            end
            if isfield(GLMP,'subunit')
                nsubs = numel(GLMP.subunit);
                if nsubs == 3
                    nsubs = 2;
                end
                for sub = 1:nsubs
                    cols = GLMP.subunit{sub}.gridX(1,:) + ceil(DN.NStixGrid(1)/2);
                    rows = -GLMP.subunit{sub}.gridY(1,:) + ceil(DN.NStixGrid(1)/2);
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
    
    % Set up linear filters
    pLpMSTA = DN.stim{u}.stats.pLpM.STA;
    mLmMSTA = DN.stim{u}.stats.mLmM.STA;
    pLmMSTA = DN.stim{u}.stats.pLmM.STA;
    mLpMSTA = DN.stim{u}.stats.mLpM.STA;
    
    % Over all stixels
    linearSTA = pLpMSTA + mLmMSTA + pLmMSTA + mLpMSTA;
    linearSTA = linearSTA - min(min(linearSTA));
    linearSTA = linearSTA ./ max(max(linearSTA));
    DN.stim{u}.stats.filters.all = linearSTA;
    
    % For subunits
    nstixgrid = DN.NStixGrid(1);
    for n = 1:numel(tempsub)
        rows = repmat(tempsub{n}.rows,[1 3]);
        cols = repmat(tempsub{n}.cols,[1 3]);
        temp = ones(size(tempsub{n}.rows));
        rgb = [temp*1 temp*2 temp*3];
        idx = sub2ind([nstixgrid nstixgrid 3],rows,cols,rgb);
        DN.stim{u}.stats.filters.subunits{n} = DN.stim{u}.stats.filters.all(idx,:);
    end

    
    
end

% Display stim #1 by default
set(dnpanel.axes,'visible','off')
set(dnpanel.dnimage,'visible','off')
set(dnpanel.axes(:,:,1),'visible','on')
set(dnpanel.dnimage(:,:,1),'visible','on')


% Set up diffpanel display
panelYpos = [25 125 225];
linearSTA = DN.stim{1}.stats.filters.all;
for xpos = 1:numel(panelXpos)
    filterpanel.single.axes(xpos) = axes('Parent',filterfig.filterpanel,'Units','pixels',...
        'box','on','Position',[panelXpos(xpos) 85 panelsize panelsize]);
    tempim = reshape(linearSTA(:,xpos),[nstixgrid nstixgrid 3]);
    h = image(tempim,'parent',filterpanel.single.axes(xpos)); hold on;
    filterpanel.linimage(xpos) = h;
    set(filterpanel.single.axes(xpos),'Xtick',[],'Ytick',[])
    if isfield(GLMP,'subunit')
        nsubs = numel(GLMP.subunit);
        if nsubs == 3
            nsubs = 2;
        end
        for sub = 1:nsubs
            cols = GLMP.subunit{sub}.gridX(1,:) + ceil(DN.NStixGrid(1)/2);
            rows = -GLMP.subunit{sub}.gridY(1,:) + ceil(DN.NStixGrid(1)/2);
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
            filterpanel.subunits(ypos,xpos,sub) = k;
        end
    end
end


% Set up spectral panel
Lfilter = sum(DN.stim{1}.stats.filters.all(1:nstixgrid^2,:));
Mfilter = sum(DN.stim{1}.stats.filters.all(nstixgrid^2+1:nstixgrid^2*2,:));
specpanel.axes.filter1 = axes('parent',filterfig.specpanel,'units','pixels','pos',[10 400 175 175]);
k = polar(Lfilter,Mfilter,'ok-');
specpanel.filter1 = k;

% Set up controls
conpanel.uicontrols.stimsel = uicontrol('Style','popupmenu','Parent',filterfig.conpanel,...
    'units','normalized','pos',[.1 .85 .8 .05],...
    'String',num2str(uniqueStim),'Value',1,'Callback',@stimsel);
conpanel.labels.stimsel = uicontrol('Parent',filterfig.conpanel,'Units','normalized',...
    'style','edit','string','Lum/Col Cone Contrast','FontSize',10,...
    'HorizontalAlignment','center','Position',get(conpanel.uicontrols.stimsel,'Position')+[0 .07 0 0]);
conpanel.uicontrols.nspikes = uicontrol('style','edit','units','normalized',...
    'pos',[.55 .8 .35 .05],'parent',filterfig.conpanel,'fontsize',12,...
    'string',num2str(DN.stim{1}.stats.pLpM.nspikes));
conpanel.uicontros.nspikeslabel = uicontrol('style','text','units','normalized',...
    'pos',[.1 .8 .4 .04],'parent',filterfig.conpanel,'fontsize',12,...
    'string','# of spikes = ');
conpanel.uicontrols.subonoff = uicontrol('Style','Checkbox','String','Turn On Subunit Visibility',...
    'pos',[20 485 150 30],'parent',filterfig.conpanel,'Value',1,'Callback',@subunitonoff);
conpanel.uicontrols.nfilters = uibuttongroup('units','normalized',...
    'parent',filterfig.conpanel,'pos',[.05 .5 .9 .2],'title','# of Filters');
conpanel.uicontrols.nfilter1 = uicontrol('style','radiobutton',...
    'parent',conpanel.uicontrols.nfilters,'pos',[15 85 100 25],...
    'string','1 Filter','fontsize',12);
conpanel.uicontrols.nfilter2 = uicontrol('style','radiobutton',...
    'parent',conpanel.uicontrols.nfilters,'pos',[15 60 100 25],...
    'string','2 Filters','fontsize',12);
conpanel.uicontrols.nfilter3 = uicontrol('style','radiobutton',...
    'parent',conpanel.uicontrols.nfilters,'pos',[15 35 100 25],...
    'string','3 Filters','fontsize',12);
conpanel.uicontrols.nfilter4 = uicontrol('style','radiobutton',...
    'parent',conpanel.uicontrols.nfilters,'pos',[15 10 100 25],...
    'string','4 Filters','fontsize',12);
conpanel.uicontrols.filtspecprops = uipanel('parent',filterfig.conpanel,...
    'units','normalized','pos',[.05 .275 .9 .2],'title','Filter Spectral Properties');
conpanel.uicontrols.colinfilter = uicontrol('parent',conpanel.uicontrols.filtspecprops,...
    'style','checkbox','pos',[10 80 200 30],'string','Colinear Filters',...
    'fontsize',12,'enable','off');
conpanel.uicontrols.orthogfilter = uicontrol('parent',conpanel.uicontrols.filtspecprops,...
    'style','checkbox','pos',[10 50 200 30],'string','Orthogonal Filters',...
    'fontsize',12,'enable','off');
conpanel.uicontrols.filtspatprops = uibuttongroup('parent',filterfig.conpanel,...
    'units','normalized','pos',[.05 .05 .9 .2],'title','Filter Spatial Properties');
conpanel.uicontrols.allsitxfilter = uicontrol('parent',conpanel.uicontrols.filtspatprops,...
    'style','radiobutton','pos',[10 80 200 30],'string','All Stixels',...
    'fontsize',12);
conpanel.uicontrols.subunitfilter = uicontrol('parent',conpanel.uicontrols.filtspatprops,...
    'style','radiobutton','pos',[10 50 200 30],'string','Subunit Only',...
    'fontsize',12);
conpanel.uicontrols.otherstixfilter = uicontrol('parent',conpanel.uicontrols.filtspatprops,...
    'style','radiobutton','pos',[10 20 200 30],'string','Select Stixels',...
    'fontsize',12);

% Save variables
set(filterfig.dnpanel,'UserData',dnpanel);
set(filterfig.conpanel,'UserData',conpanel);
set(filterfig.filterpanel,'UserData',filterpanel);
set(filterfig.specpanel,'UserData',specpanel);
set(gcf,'UserData',filterfig);


end



function stimsel(hObj,~)
global DN
% For selecting between different lum/col Dense Noise conditions

% Load figure variables
filterfig = get(gcf,'UserData');
dnpanel = get(filterfig.dnpanel,'UserData');
conpanel = get(filterfig.conpanel,'UserData');
filterpanel = get(filterfig.filterpanel,'UserData');
val = get(hObj,'Value');
%nframesback = 10;

set(dnpanel.axes,'Visible','off');
set(dnpanel.dnimage,'Visible','off');
set(dnpanel.axes(:,:,val),'Visible','on');
set(dnpanel.dnimage(:,:,val),'Visible','on');
set(conpanel.uicontrols.nspikes,'string',num2str(DN.stim{val}.stats.pLpM.nspikes));
set(filterpanel.axes,'Visible','off');
set(filterpanel.linimage,'visible','off');
set(filterpanel.axes(:,:,val),'visible','on');
set(filterpanel.linimage(:,:,val),'visible','on');


% Save variables
set(filterfig.dnpanel,'UserData',dnpanel);
set(filterfig.conpanel,'UserData',conpanel);
set(filterfig.filterpanel,'UserData',filterpanel);
set(gcf,'UserData',filterfig);

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
    if get(conpanel.uicontrols.subonoff,'Value') == 1
        set(dnpanel.subunits,'Visible','on')
        set(filterpanel.subunits,'visible','on')
    else
        set(dnpanel.subunits,'Visible','off')
        set(filterpanel.subunits,'visible','off')
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
 
