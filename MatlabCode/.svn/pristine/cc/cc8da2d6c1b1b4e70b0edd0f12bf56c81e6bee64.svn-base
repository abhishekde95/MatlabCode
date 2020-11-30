function GLMSGUI_Mean2Var(a,b)
global GLMP

figure(65); clf;
set(gcf,'units','pixels','pos',[300 100 1000 400],'NumberTitle','off',...
    'Name',['Mean to Variance Analysis (' GLMP.datafile ')']);

% Set up panels
m2vfig = get(gcf,'UserData');
m2vfig.conpanel = uipanel('pos',[.01 .01 .33 .98],'parent',gcf,...
    'title','Control Panel');
m2vfig.fitspanel = uipanel('Pos',[.35 .01 .69 .98],'Parent',gcf,...
    'title','Fits to Data');

% Display stimuli
conpanel.stimmap = axes('parent',m2vfig.conpanel,'units','normalized',...
    'pos',[.1 .55 .8 .4]);
conpanel.stim = polar(GLMP.subunit{1}.uniquetheta,GLMP.subunit{1}.uniquerho,'ro',...
    'parent',conpanel.stimmap);

% Display info
fitspanel.axes = axes('parent',m2vfig.fitspanel,'pos',[.1 .1 .8 .8],'box','on');

% Set up Control Panel
strs = {'All Directions' 'Select Directions'};
conpanel.whichdirectionsmenu = uicontrol(m2vfig.conpanel,'style','popupmenu',...
    'units','normalized','pos',[.1 .05 .8 .1],'string',strs,...
    'fontsize',12,'callback',@TuneDirSelection);

% Set up Subunit Selection
conpanel.subbuttons = uibuttongroup('Parent',m2vfig.conpanel,...
    'units','normalized','pos',[.6 .2 .35 .25],...
    'SelectionChangeFcn',@subsel);
conpanel.sub1 = uicontrol('style','radiobutton','units','normalized',...
    'parent',conpanel.subbuttons,'pos',[.1 .7 .8 .25],...
    'string','Subunit #1','fontsize',12);
conpanel.sub2 = uicontrol('style','radiobutton','units','normalized',...
    'parent',conpanel.subbuttons,'pos',[.1 .4 .8 .25],...
    'string','Subunit #2','fontsize',12);
conpanel.sub3 = uicontrol('style','radiobutton','units','normalized',...
    'parent',conpanel.subbuttons,'pos',[.1 .1 .8 .25],...
    'string','Subunit #3','fontsize',12);
conpanel.subselect = 1;
if isempty(GLMP.subunit{2})
    set(conpanel.sub2,'enable','off')
end
if isempty(GLMP.subunit{3})
    set(conpanel.sub3,'enable','off')
end

% Fit which function
conpanel.fitfuncspanel = uipanel(m2vfig.conpanel,'units','normalized',...
    'pos',[.05 .2 .5 .25]);
conpanel.poisson = uicontrol(conpanel.fitfuncspanel,'style','checkbox','fontsize',12,...
    'units','normalized','pos',[.05 .1 .8 .15],'string','Poisson','callback',@plotfig);
conpanel.negbinomial = uicontrol(conpanel.fitfuncspanel,'style','checkbox','fontsize',12,....
    'units','normalized','pos',[.05 .4 .8 .15],'string','Negative Binomail','callback',@plotfig);
conpanel.quasipoisson = uicontrol(conpanel.fitfuncspanel,'style','checkbox','fontsize',12,...
    'units','normalized','pos',[.05 .7 .8 .15],'string','Quasi-Poisson','callback',@plotfig);

conpanel.means = GLMP.subunit{1}.meannspikes;
conpanel.vars = GLMP.subunit{1}.varnspikes;

% Save figure variables
set(m2vfig.conpanel,'userdata',conpanel)
set(m2vfig.fitspanel,'userdata',fitspanel)
set(gcf,'userdata',m2vfig)

plotfig()

end

function subsel(~,b)
global GLMP

% Load variables
m2vfig = get(gcf,'UserData');
conpanel = get(m2vfig.conpanel,'UserData');
fitspanel = get(m2vfig.fitspanel,'UserData');

subval = get(b.NewValue,'string');
axes(conpanel.stimmap); cla;
if strcmp(subval,'Subunit #1')
    conpanel.stim = polar(GLMP.subunit{1}.uniquetheta,GLMP.subunit{1}.uniquerho,'ro',...
        'parent',conpanel.stimmap);
    conpanel.subselect = 1;
elseif strcmp(subval,'Subunit #2')
    conpanel.stim = polar(GLMP.subunit{2}.uniquetheta,GLMP.subunit{2}.uniquerho,'go');
    conpanel.subselect = 2;
elseif strcmp(subval,'Subunit #3')
    conpanel.stim = polar(GLMP.subunit{3}.uniquetheta,GLMP.subunit{3}.uniquerho,'bo',...
        'parent',conpanel.stimmap);
    conpanel.subselect = 3;
end
conpanel.whichdirectionsmenu.Value = 1;

% Save figure variables
set(m2vfig.conpanel,'userdata',conpanel)
set(m2vfig.fitspanel,'userdata',fitspanel)
set(gcf,'userdata',m2vfig)

TuneDirSelection()
plotfig()

end


function TuneDirSelection(~,~)
global GLMP

% Load figure variables
m2vfig = get(gcf,'userdata');
conpanel = get(m2vfig.conpanel,'userdata');
fitspanel = get(m2vfig.fitspanel,'userdata');

% Turn on direciton selection and clear plot
if conpanel.whichdirectionsmenu.Value == 1
    whichsub = conpanel.subselect;
    conpanel.means = GLMP.subunit{whichsub}.meannspikes;
    conpanel.vars = GLMP.subunit{whichsub}.varnspikes;
elseif conpanel.whichdirectionsmenu.Value == 2
    set(conpanel.stim,'ButtonDownFcn',@SelectDirection)
    axes(fitspanel.axes); cla;
end

% Save figure variables
set(m2vfig.conpanel,'userdata',conpanel)
set(gcf,'userdata',m2vfig);

end

function SelectDirection(a,b)
global GLMP

% Grab current point before anything else
h = gca;
whichpt = get(h,'CurrentPoint');
whichpt = whichpt(1,[1 2]);

% Load variables
m2vfig = get(gcf,'UserData');
conpanel = get(m2vfig.conpanel,'UserData');
fitspanel = get(m2vfig.fitspanel,'UserData');

% Grab other variables
whichsub = conpanel.subselect;
[theta,~] = cart2pol(whichpt(1),whichpt(2));
[~,idx] = min(abs(GLMP.subunit{whichsub}.theta - theta));
idx = find(GLMP.subunit{whichsub}.theta==GLMP.subunit{whichsub}.theta(idx));
theta = GLMP.subunit{whichsub}.theta(idx(1));
L = GLMP.subunit{whichsub}.uniquetheta == theta;
maxrho = max(GLMP.subunit{whichsub}.uniquerho(L));

% Save stimuli inside conpanel
conpanel.means = GLMP.subunit{whichsub}.meannspikes(L);
conpanel.vars = GLMP.subunit{whichsub}.varnspikes(L);

% Plot stim 
cols = {'r' 'g' 'b'};
axes(conpanel.stimmap); cla; hold on;
conpanel.stim = polar(GLMP.subunit{whichsub}.uniquetheta,GLMP.subunit{whichsub}.uniquerho,[cols{whichsub} 'o'],...
        'parent',conpanel.stimmap);
polar([theta theta],[0 maxrho],'k-')
set(conpanel.stim,'ButtonDownFcn',@SelectDirection)

% Plot mean points
axes(fitspanel.axes); cla; hold on;
plot(conpanel.means,conpanel.vars,'k*')

% Save figure variables
set(m2vfig.conpanel,'userdata',conpanel)
set(m2vfig.fitspanel,'userdata',fitspanel)
set(gcf,'userdata',m2vfig)

plotfig()

end


function plotfig(~,~)

% Load figure variables
m2vfig = get(gcf,'userdata');
conpanel = get(m2vfig.conpanel,'userdata');
fitspanel = get(m2vfig.fitspanel,'userdata');

% Set up some variables
means = conpanel.means;
vars = conpanel.vars;
L = vars == 0;
vars(L) = [];
means(L) = [];
if isempty(means)
    disp('Data has no variance...')
    return
end
x = linspace(0,max(means),100);
leg = {};

% Set up axes (clear and replot every call)
axes(fitspanel.axes); cla; hold on;
xlabel('Measured Mean (#spikes)')
ylabel('Variance (#spikes)')
if conpanel.poisson.Value == 1
    %poisson_neg_LL = sum(GLMP.meannspikes)-(log(GLMP.meannspikes'+.001)*GLMP.varnspikes);
    plot(x,x,'r--')
    leg{numel(leg)+1} = 'Var = Mean';
end
if conpanel.quasipoisson.Value == 1
    data = means;
    b = regress(vars,data);
    y = b(end).*x;
    plot(x,y,'g--')
    leg{numel(leg)+1} = ['Var = ' num2str(b(end)) ' * Mean'];
end
if conpanel.negbinomial.Value == 1
    data = means.^2;
    b = regress(vars-means,data);
    y = x + b(end)*x.^2;
    plot(x,y,'b--')
    leg{numel(leg)+1} = ['Var = Mean + ' num2str(b(end)) ' * Mean^2'];
end
if ~isempty(leg)
    legend(leg,'fontsize',12,'location','northwest')
end
plot(means,vars,'k*')

% Save figure variables
set(m2vfig.conpanel,'userdata',conpanel)
set(m2vfig.fitspanel,'userdata',fitspanel)
set(gcf,'userdata',m2vfig)

end