function GLMSGUI_CompareDNGLMP(~,~)
global GLMP

% Set up figure
figure(91); clf;
set(gcf,'pos',[100 200 1200 500],'NumberTitle','off','Name',['Correlation Analysis (' GLMP.datafile ')'])
CompareFig.TCpanel = uipanel('parent',gcf,'units','normalized','pos',[.01 .01 .79 .38]);
CompareFig.GLMPpanel = uipanel('parent',gcf,'units','normalized','pos',[.01 .4 .24 .59]);
CompareFig.corrpanelT = uipanel('parent',gcf,'units','normalized','pos',[.26 .4 .2 .59]);
CompareFig.corrpanelP = uipanel('parent',gcf,'units','normalized','pos',[.47 .4 .33 .59]);
CompareFig.STApanel = uipanel('parent',gcf,'units','normalized','pos',[.81 .01 .18 .98]);

% Set default variables
CompareFig.colors = [.8 .8 1; .3 .3 .5; .8 .5 .5; .5 .8 .5];
CompareFig.stim = {'pLpM' 'mLmM' 'pLmM' 'mLpM'};
CompareFig.labels = {'+L+M' '-L-M' '+L-M' '-L+M'};
CompareFig.angs = [pi/4 -3*pi/4 -pi/4 3*pi/4];

% Set up corrT (table) panel
corrpanelT.table = uitable('parent',CompareFig.corrpanelT,'units','normalized','pos',[.1 .4 .8 .5]);
corrpanelT.corrval = uicontrol('parent',CompareFig.corrpanelT,'units','normalized',...
    'pos',[.2 .25 .6 .1],'style','edit');
subidx = find(~cellfun(@isempty,GLMP.subunit));
subsel = cell(numel(subidx),1);
for n = 1:numel(subidx) %This is so sloppy! Can cellfun fix this?
    subsel{n} = ['Subunit #' num2str(subidx(n))];
end
corrpanelT.subunit = uicontrol('style','popupmenu','parent',CompareFig.corrpanelT,'units','normalized',...
    'pos',[.1 .05 .8 .1],'value',1,'string',subsel,'callback',@switchsub);

% Set up CorrP (plot) panel
corrpanelP.axes = axes('parent',CompareFig.corrpanelP,'units','normalized','pos',[.2 .2 .7 .7],'box','on');
xlabel('WN STA'); ylabel('GLMP STA');

% Set up GLMP panel
glmppanel.axes = axes('parent',CompareFig.GLMPpanel,'units','normalized','pos',[.15 .15 .75 .75],'box','on');
xlabel('Lcc'); ylabel('Mcc');
glmppanel.sub = 1;

% Set up Time Course panel
tcpanel.axes = axes('parent',CompareFig.TCpanel,'units','normalized','pos',[.05 .15 .9 .75],'box','on',...
    'ButtonDownFcn',@switchframe);
tcpanel.frame = 6;

% Set up STA panel
yvals = linspace(.96,.04,5);
for n = 1:numel(CompareFig.stim)
    stapanel.axes.(CompareFig.stim{n}) = axes('parent',CompareFig.STApanel,'pos',[.3 yvals(n+1) .6 .2]);
    set(stapanel.axes.(CompareFig.stim{n}),'box','on','XTick',[],'YTick',[],'ydir','reverse')
    stapanel.labels.(CompareFig.stim{n}) = text(-.2,.5,(CompareFig.labels{n}),'Rotation',90,'Units','normalized',...
        'HorizontalAlignment','center','fontsize',15);
end

% Save figure variables
set(CompareFig.TCpanel,'userdata',tcpanel)
set(CompareFig.STApanel,'userdata',stapanel)
set(CompareFig.GLMPpanel,'userdata',glmppanel)
set(CompareFig.corrpanelP,'userdata',corrpanelP)
set(CompareFig.corrpanelT,'userdata',corrpanelT)
set(gcf,'userdata',CompareFig)

% Fill in figure panel by panel
ExecuteAnal()

end

function ExecuteAnal()

PlotTC()
PlotSTA()
PlotGLMP()
CorrelationTable()
CorrelationFigure()

end
function switchsub()

% Load figure variables
CompareFig = get(gcf,'UserData');
corrpanelT = get(CompareFig.corrpanelT,'userdata');

% Lookup which subunit and save idx
str = get(corrpanelT.subunit,'String');
sel = get(corrpanelT.subunit,'Value');
glmppanel.sub = str2double(str{sel}(end));

% Save figure variables
set(CompareFig.corrpanelT,'userdata',corrpanelT)
set(gcf,'userData',CompareFig)

end

function switchframe(~,b)

% Load figure variables
CompareFig = get(gcf,'UserData');
tcpanel = get(CompareFig.TCpanel,'userdata');

tcpanel.frame = round(b.IntersectionPoint(1));

% Save figure variables
set(CompareFig.TCpanel,'userdata',tcpanel)
set(gcf,'userData',CompareFig)

% Fill in figure panel by panel
ExecuteAnal()

end

function PlotTC()
global DN GLMP

% Load figure variables
CompareFig = get(gcf,'UserData');
tcpanel = get(CompareFig.TCpanel,'UserData');
glmppanel = get(CompareFig.GLMPpanel,'UserData');
stapanel = get(CompareFig.STApanel,'UserData');

% Load some variables from figure
nstix = DN.NStixGrid(1);
sub = glmppanel.sub;
cols = GLMP.subunit{sub}.gridX{1} + ceil(nstix/2);
rows = -GLMP.subunit{sub}.gridY{1} + ceil(nstix/2);
tcpanel.idx = sub2ind([nstix nstix],rows,cols);
tcpanel.STAs = nan(nstix,nstix,numel(CompareFig.stim),10);

% Plot White Noise STA over frames
axes(tcpanel.axes); cla; hold on; ylabel('STA'); xlabel('Frames')
for n = 1:numel(CompareFig.stim)
    tcpanel.STAs(:,:,n,:) = abs(reshape(DN.stats.(CompareFig.stim{n}).STA(1:nstix^2,:),nstix,nstix,1,[]));
    tcpanel.(CompareFig.stim{n}).stixSTA = abs(DN.stats.(CompareFig.stim{n}).STA(tcpanel.idx,:));
    tcpanel.(CompareFig.stim{n}).meanSTA = abs(mean(tcpanel.(CompareFig.stim{n}).stixSTA));
    plot(tcpanel.(CompareFig.stim{n}).meanSTA,'col',CompareFig.colors(n,:))
end

% Plot frame marker
yval = get(gca,'ylim');
theo2STD = 2*sqrt((.25*.75)/DN.stats.pLpM.nspikes);
plot([tcpanel.frame tcpanel.frame],yval,'r-');
plot([1 numel(tcpanel.pLpM.meanSTA)],[.25+theo2STD .25+theo2STD],'r--')
plot([1 numel(tcpanel.pLpM.meanSTA)],[.25-theo2STD .25-theo2STD],'r--')
set(gca,'ylim',yval)

% Save figure variables
set(CompareFig.TCpanel,'userdata',tcpanel)
set(gcf,'userData',CompareFig)

end

function PlotSTA()
global DN

% Load figure variables
CompareFig = get(gcf,'UserData');
tcpanel = get(CompareFig.TCpanel,'UserData');
stapanel = get(CompareFig.STApanel,'UserData');

% Initialize some variables
[row,col] = ind2sub([DN.NStixGrid(1) DN.NStixGrid(1)],tcpanel.idx);

% Load whole STAs for the selected frame%
STAs = (tcpanel.STAs-min(tcpanel.STAs(:)))./(max(tcpanel.STAs(:))-min(tcpanel.STAs(:))); % normalize over range (so min=0, max=1)
for n = 1:numel(CompareFig.stim) % and plot
    axes(stapanel.axes.(CompareFig.stim{n})); hold on;
    tempSTA = repmat(STAs(:,:,n,tcpanel.frame),1,1,3);
    image(tempSTA); axis tight;
    plot(col,row,'ro')
end

% Save figure variables
set(CompareFig.STApanel,'userdata',stapanel)
set(gcf,'userData',CompareFig)

end



function PlotGLMP()
global GLMP

% Load figure variables
CompareFig = get(gcf,'UserData');
glmppanel = get(CompareFig.GLMPpanel,'UserData');

sub = glmppanel.sub;

% Do indexing of WN stimuli/stixels
normtot = 0;
for n = 1:numel(CompareFig.stim)
    idx1 = find(GLMP.subunit{sub}.uniquetheta == CompareFig.angs(n));
    [~,idx2] = max(GLMP.subunit{sub}.uniquerho(idx1));
    glmppanel.(CompareFig.stim{n}).idx = idx1(idx2);
    glmppanel.(CompareFig.stim{n}).Lcc = GLMP.subunit{sub}.uniqueLcc(idx1(idx2));
    glmppanel.(CompareFig.stim{n}).Mcc = GLMP.subunit{sub}.uniqueMcc(idx1(idx2));
    glmppanel.(CompareFig.stim{n}).meannsp = GLMP.subunit{sub}.meannspikes(idx1(idx2));
    normtot = normtot + glmppanel.(CompareFig.stim{n}).meannsp;
end

% Plot stimuli
axes(glmppanel.axes); cla; axis equal; hold on; box on;
for n = 1:numel(CompareFig.stim)
    glmppanel.(CompareFig.stim{n}).normresp = glmppanel.(CompareFig.stim{n}).meannsp/normtot;
    p = plot(glmppanel.(CompareFig.stim{n}).Lcc,glmppanel.(CompareFig.stim{n}).Mcc,...
        'o','color',CompareFig.colors(n,:),'MarkerSize',(glmppanel.(CompareFig.stim{n}).normresp+1)*10);
    set(p,'MarkerFaceColor',CompareFig.colors(n,:));
end
xlim([-max(GLMP.subunit{sub}.Lcc)*1.1 max(GLMP.subunit{sub}.Lcc)*1.1]);
ylim([-max(GLMP.subunit{sub}.Mcc)*1.1 max(GLMP.subunit{sub}.Mcc)*1.1]);

% Save figure variables
set(CompareFig.GLMPpanel,'userdata',glmppanel)
set(gcf,'userdata',CompareFig)

end

function CorrelationTable()
global DN GLMP

% Load figure variables
CompareFig = get(gcf,'userdata');
tcpanel = get(CompareFig.TCpanel,'userdata');
glmppanel = get(CompareFig.GLMPpanel,'userdata');
corrpanelT = get(CompareFig.corrpanelT,'userdata');

% Initialize/load/set variables
data = nan(4,2);
frame = tcpanel.frame;
sub = glmppanel.sub;

% Construct 2x4 table
for n = 1:numel(CompareFig.stim);
    data(n,1) = tcpanel.(CompareFig.stim{n}).meanSTA(frame);
    data(n,2) = glmppanel.(CompareFig.stim{n}).normresp;
end

% Plot table
set(corrpanelT.table,'data',data,'RowName',[],'fontsize',15,'columnname',{'WN' 'GLMP'},...
    'unit','normalized','BackgroundColor',CompareFig.colors,'pos',[.1 .4 .8 .5]);
set(corrpanelT.table,'units','pixels')
width = corrpanelT.table.Position(3)/2;
set(corrpanelT.table,'ColumnWidth',{width width})
corrpanelT.table.Position(3:4) = corrpanelT.table.Extent(3:4);

% Set up multinomial covariance matrix
W = ones(3,3);
p = .25;
W = W * p.^2;
for n = 1:size(W,1)
    W(n,n) = p * (1-p);
end
W = W ./ numel(GLMP.subunit{sub}.gridX{1});

% find p-value for the given STA
xbar = data(1:size(W,1),1);
pmat = repmat(p,size(W,1),1);
n = numel(CompareFig.stim);
nsp = DN.stats.pLpM.nspikes;
%tsq = nsp*(xbar - pmat)' * inv(W) * (xbar-pmat);
tsq = nsp*(xbar - pmat)' * (W\(xbar-pmat));
f = ((nsp-n)/(n*(nsp-1))) * tsq;
pval = 1 - fcdf(f,n,nsp-n);
set(corrpanelT.corrval,'string',['p-val = ' num2str(pval)],'fontsize',15)

% Calculate and Display correlation value
%corrval = rnddec(corr(data(:,1),data(:,2)),4);
%set(corrpanelT.corrval,'string',['rho = ' num2str(corrval)],'fontsize',15)

% Save figure variables
set(CompareFig.TCpanel,'userdata',tcpanel)
set(CompareFig.GLMPpanel,'userdata',glmppanel)
set(CompareFig.corrpanelT,'userdata',corrpanelT)
set(gcf,'userdata',CompareFig)

end

function CorrelationFigure()

% Load figure variables
CompareFig = get(gcf,'userdata');
corrpanelT = get(CompareFig.corrpanelT,'userdata');
corrpanelP = get(CompareFig.corrpanelP,'userdata');

% Pull out data from table
data = get(corrpanelT.table,'data');
%onemat = ones(size(data(:,1)));
%b = regress(data(:,1),cat(2,onemat,data(:,2)));

% Plot points
axes(corrpanelP.axes); cla; hold on;
for n = 1:numel(CompareFig.stim)
    plot(data(n,1),data(n,2),'*','color',CompareFig.colors(n,:))
end
%xvals = linspace(min(data(:,1)),max(data(:,1)),50);
%yvals = (xvals * b(1)) + b(2);
%plot(xvals,yvals,'k--')

% Save figure variables
set(CompareFig.corrpanelT,'userdata',corrpanelT)
set(CompareFig.corrpanelP,'userdata',corrpanelP)
set(gcf,'userdata',CompareFig)

end
