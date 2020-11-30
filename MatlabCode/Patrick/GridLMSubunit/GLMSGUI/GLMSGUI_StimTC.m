function GLMSGUI_StimTC(~,~)
global GLMP


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  GLMP Time Course Analysis functions  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(41); clf;
set(gcf,'numbertitle','off','name',['GLMP Time Course of Responses (' GLMP.datafile ')'],'pos',[50 100 800 500])
StimTC = get(gcf,'UserData');
StimTC.controls = uipanel('parent',gcf,'units','normalized',...
    'pos',[.025 .025 .35 .45]);
StimTC.stimmap = uipanel('parent',gcf,'units','normalized',...
    'pos',[.025 .5 .35 .475]);
StimTC.movie = uipanel('parent',gcf,'units','normalized',...
    'pos',[.4 .025 .575 .95]);

% Load variables
cpanel = get(StimTC.controls,'UserData');
mappanel = get(StimTC.movie,'UserData');
movpanel = get(StimTC.stimmap,'UserData');

% Subunit Selection
cpanel.uicontrols.subbuttons = uibuttongroup('Parent',StimTC.controls,...
    'units','normalized','pos',[.025 .525 .5 .45],...
    'title','Subunit Selection','SelectionChangeFcn',@StimTCsubsel);
cpanel.uicontrols.sub1 = uicontrol('style','radiobutton',...
    'parent',cpanel.uicontrols.subbuttons,'pos',[15 58 100 25],...
    'string','Subunit #1','fontsize',12);
cpanel.uicontrols.sub2 = uicontrol('style','radiobutton',...
    'parent',cpanel.uicontrols.subbuttons,'pos',[15 33 100 25],...
    'string','Subunit #2','fontsize',12);
cpanel.uicontrols.sub3 = uicontrol('style','radiobutton',...
    'parent',cpanel.uicontrols.subbuttons,'pos',[15 8 100 25],...
    'string','Subunit #3','fontsize',12);
cpanel.subselect = 1;
if isempty(GLMP.subunit{2})
    set(cpanel.uicontrols.sub2,'enable','off')
end
if isempty(GLMP.subunit{3})
    set(cpanel.uicontrols.sub3,'enable','off')
end
cpanel.uicontrols.playmov = uicontrol('style','pushbutton',...
    'parent',StimTC.controls,'units','normalized','pos',[.6 .05 .35 .2],...
    'string','Play Movie','fontsize',12,'Callback',@GenerateTC);
cpanel.uicontrols.gaussfilter = uipanel('parent',StimTC.controls,...
    'units','normalized','pos',[.025 .025 .5 .45],'title','Gaussian Filter');
cpanel.uicontrols.sigma = uicontrol('parent',cpanel.uicontrols.gaussfilter,...
    'style','edit','units','normalized','pos',[.5 .25 .45 .45],'string',.02);
cpanel.labels.sigma = uicontrol('parent',cpanel.uicontrols.gaussfilter,...
    'style','text','units','normalized','pos',[.05 .35 .45 .2],'string','Sigma = ');

% Draw stimulus map for stim #1 (as an initial condition)
mappanel.axes.map = axes('parent',StimTC.stimmap,...
    'units','normalized','pos',[.1 .1 .8 .8]);
grid on; hold on;
markSize = (((GLMP.subunit{1}.meanfr) /(max(GLMP.subunit{1}.meanfr)))+.1)*50;
for t = 1:numel(GLMP.subunit{1}.uniqueLcc)
    plot(GLMP.subunit{1}.uniqueLcc(t),GLMP.subunit{1}.uniqueMcc(t),'ro','Markersize',markSize(t));
end
axis equal tight;
xlabel('Lcc')
ylabel('Mcc')

% Set up movie panel
maxcc = max(max(GLMP.subunit{1}.Lcc),max(GLMP.subunit{1}.Mcc));
movpanel.axes.movie = axes('parent',StimTC.movie,'units','normalized',...
    'pos',[.1 .15 .8 .8],...
    'XLim',[-maxcc maxcc],'YLim',[-maxcc maxcc],...
    'ZLim',[0 max(GLMP.subunit{1}.nspikes)*1.1],...
    'XGrid','on','YGrid','on','ZGrid','on','box','on');
movpanel.uicontrols.tcslider = uicontrol('parent',StimTC.movie,...
    'style','slider','min',0,'max',.5,'units','normalized',...
    'pos',[.1 .05 .8 .05],'callback',@ChooseTP);


% Save variables
set(StimTC.controls,'UserData',cpanel);
set(StimTC.stimmap,'UserData',mappanel);
set(StimTC.movie,'UserData',movpanel);
set(gcf,'UserData',StimTC);


end


function StimTCsubsel(~,eventdata)
global GLMP
% Set up subunit selection

% Load variables
StimTC = get(gcf,'UserData');
cpanel = get(StimTC.controls,'UserData');
mappanel = get(StimTC.stimmap,'UserData');

subval = get(eventdata.NewValue,'string');
axes(mappanel.axes.map); cla;
if strcmp(subval,'Subunit #1')
    cpanel.subselect = 1;
    mfc = 'ro';
elseif strcmp(subval,'Subunit #2')
    cpanel.subselect = 2;
    mfc = 'go';
elseif strcmp(subval,'Subunit #3')
    cpanel.subselect = 3;
    mfc = 'bo';
end
%set(cpanel.axes.stim,'ButtonDownFcn',@StimTCStimulusSelection)

markSize = (((GLMP.subunit{cpanel.subselect}.meanfr) /(max(GLMP.subunit{cpanel.subselect}.meanfr)))+.1)*50;
grid on; hold on;
for t = 1:numel(GLMP.subunit{cpanel.subselect}.uniqueLcc)
    h = plot(GLMP.subunit{cpanel.subselect}.uniqueLcc(t),...
        GLMP.subunit{cpanel.subselect}.uniqueMcc(t),...
        mfc,'Markersize',markSize(t));
end
axis equal tight;
xlabel('Lcc')
ylabel('Mcc')

% Save variables
set(StimTC.controls,'UserData',cpanel);
set(StimTC.stimmap,'UserData',mappanel);
set(gcf,'UserData',StimTC);

end

function GenerateTC(~,~)
global GLMP
% Does the actual time course calculation for each time point

% Load variables
StimTC = get(gcf,'UserData');
cpanel = get(StimTC.controls,'UserData');
mappanel = get(StimTC.stimmap,'UserData');
movpanel = get(StimTC.movie,'UserData');

s = cpanel.subselect;
maxcc = max(max(GLMP.subunit{s}.Lcc),max(GLMP.subunit{s}.Mcc));

% Convolve with Gaussian profile
mintpt = -.2;
%maxtpt = mean(GLMP.subunit{s}.stimDur)+.2;
maxtpt = .6;
tptstep = .01;
edges = mintpt:tptstep:maxtpt;
movpanel.deets.edges = edges;
x = 0:.01:.1;
sigma = get(cpanel.uicontrols.sigma,'string');
sigma = str2double(sigma);
gauss = normpdf(x,mean(x),sigma);
gauss = gauss./max(gauss);
for n = 1:numel(GLMP.subunit{s}.meanfr)
    rast = histc(GLMP.subunit{s}.spiketimes_cat{n},edges);
    gauss = gauss/max(gauss);
    spksmooth(n,:) = conv(rast,gauss,'same');
end

% Configure slider bar
sliderstring = rnddec(linspace(mintpt,maxtpt,9),1);
set(movpanel.uicontrols.tcslider,'min',mintpt,'max',maxtpt,...
    'sliderstep',[tptstep tptstep*10],'value',mintpt,...
    'string',sliderstring,'Callback',@sliderSurfSel)
movpanel.labels.sliderTP = axes('parent',StimTC.movie,...
    'units','normalized','pos',[.125 .05 .75 .001],...
    'XTickLabel',sliderstring,'XLim',[mintpt maxtpt],'YTickLabel',[]);

% Create a linear interpolation surface for each time point
%M = cell(size(GLMP.subunit{s}.spksmooth,2),1);
for n = 1:size(spksmooth,2)
    x = GLMP.subunit{s}.uniqueLcc;
    y = GLMP.subunit{s}.uniqueMcc;
    z =  spksmooth(:,n);
    F = TriScatteredInterp(x,y,z);
    [qx,qy] = meshgrid(min(x):.01:max(x),min(y):.01:max(y));
    
    % Save x,y,z
    GLMP.subunit{s}.spikesurf{n,1} = qx;
    GLMP.subunit{s}.spikesurf{n,2} = qy;
    GLMP.subunit{s}.spikesurf{n,3} = F(qx,qy);
    
    % Concatinate frames into a movie
    axes(movpanel.axes.movie);
    view = get(gca,'view');
    surfc(qx,qy,F(qx,qy))
    set(gca,'View',view,'xlim',[-maxcc maxcc],'ylim',[-maxcc maxcc],...
        'zlim',[0 max(GLMP.subunit{s}.nspikes)*1.1])
    %M{n} = getframe(gca);
    drawnow;
    
    % Adjust slider
    set(movpanel.uicontrols.tcslider,'Value',edges(n))
    
end

%GLMP.subunit{s}.movie = M;

% Save variables
set(StimTC.controls,'UserData',cpanel);
set(StimTC.stimmap,'UserData',mappanel);
set(StimTC.movie,'UserData',movpanel);
set(gcf,'UserData',StimTC);


end


function sliderSurfSel(~,~)
global GLMP

% Load variables
StimTC = get(gcf,'UserData');
cpanel = get(StimTC.controls,'UserData');
mappanel = get(StimTC.stimmap,'UserData');
movpanel = get(StimTC.movie,'UserData');

s = cpanel.subselect;
maxcc = max(max(GLMP.subunit{s}.Lcc),max(GLMP.subunit{s}.Mcc));

% Get frame
frame = get(movpanel.uicontrols.tcslider,'Value');
[junk,idx] = min(abs(frame - movpanel.deets.edges));

% Draw surface
view = get(gca,'view');
axes(movpanel.axes.movie);
qx = GLMP.subunit{s}.spikesurf{idx,1};
qy = GLMP.subunit{s}.spikesurf{idx,2};
qz = GLMP.subunit{s}.spikesurf{idx,3};
surfc(qx,qy,qz)
set(gca,'View',view,'xlim',[-maxcc maxcc],'ylim',[-maxcc maxcc],...
    'zlim',[0 max(GLMP.subunit{s}.nspikes)*1.1]);
    
% Save variables
set(StimTC.controls,'UserData',cpanel);
set(StimTC.stimmap,'UserData',mappanel);
set(StimTC.movie,'UserData',movpanel);
set(gcf,'UserData',StimTC);
    
end

