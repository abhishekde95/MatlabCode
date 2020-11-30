function GLMSPopGUI_PredDiff
global GLMSPopData

% Grab saved population data
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
load([library 'GLMSPopData.mat'])


% Set up figure and panels
figure(6458); clf; 
set(gcf,'NumberTitle','off','Name','Population Cell Selection','units','normalized',...
    'pos',[.1 .1 .8 .8])
prediffig.table = uitable('units','normalized','pos',[.01 .51 .48 .48],...
    'BackgroundColor',[1 1 1]);%,'cellselectioncallback',@cellselect);
prediffig.poppanel = uipanel('units','normalized','pos',[.01 .01 .48 .48]);
prediffig.indivpanel = uipanel('units','normalized','pos',[.5 .01 .49 .98]);


% Individual panel
indpan.onedsurf = axes('parent',prediffig.indivpanel,'units','normalized',...
    'pos',[.1 .6 .35 .35]); cla; hold on; grid on; box on;
title('1D Surface')
indpan.twodsurf = axes('parent',prediffig.indivpanel,'units','normalized',...
    'pos',[.6 .6 .35 .35]); cla; hold on; grid on; box on;
title('2D Surface')
indpan.bubbleplot = axes('parent',prediffig.indivpanel,'units','normalized',...
    'pos',[.1 .1 .35 .35]); cla; hold on; grid on; box on;
indpan.diffplot = axes('parent',prediffig.indivpanel,'units','normalized',...
    'pos',[.6 .1 .35 .35]); cla; hold on; grid on; box on;

% Population panel
popan.axes = axes('parent',prediffig.poppanel,'units','normalized',...
    'pos',[.1 .1 .8 .8]); hold on; grid on; box on;

%Save variables
set(prediffig.indivpanel,'UserData',indpan);
set(prediffig.poppanel,'userdata',popan);
set(gcf,'UserData',prediffig);

UnpackPopData()
CalcMaxDiff()

PopulateTable()

PlotPopData()

end

function PopulateTable()
global GLMSPopData

% Load user data
prediffig = get(gcf,'userdata');
popan = get(prediffig.poppanel,'UserData');

% Populate table
datatypes = GLMSPopData(1,:);
data(:,1) = GLMSPopData(2:end,strcmp(datatypes,'Datafile'));
data(:,2) = GLMSPopData(2:end,strcmp(datatypes,'Subunit'));
data(:,3) = GLMSPopData(2:end,strcmp(datatypes,'nD'));
data(:,4) = num2cell(popan.vars.maxdiffnsp);
data(:,5) = num2cell(popan.vars.diffnormLL);
data(:,6) = num2cell(([GLMSPopData{2:end,strcmp(datatypes,'Tuning')}]'));
data(:,7) = num2cell(popan.vars.PDstd);
data(:,8) = num2cell(popan.vars.mod);
colname = {'Datafile' 'Sub' 'nD' 'MaxDiffnsp' 'NormLLDiff' 'Tuning' 'Tun Std' 'Mod'};
colformat = cell(1,size(data,2));
colformat(:) = {'char'};
set(prediffig.table,'data',data,'columnname',colname,...
    'columnformat',colformat,'CellSelectionCallback',@cellselect);

end



function PlotPopData(~,~)

% Load user data
prediffig = get(gcf,'userdata');
popan = get(prediffig.poppanel,'UserData');

% Plot population max differences vs norm-LL difference (LL improvement from 2D fit)
axes(popan.axes); cla; hold on; grid on; box on;
h = plot(popan.vars.diffnormLL,popan.vars.maxdiffnsp,'ko');
xlabel('2D Norm LL - 1D Norm LL')
ylabel('Maximum Difference Between Fits (spike count)')
set(popan.axes,'ButtonDownFcn',@cellselect)
set(h,'ButtonDownFcn',@cellselect)

% % Set table colors
% PopulateTable()
% cols = ones(numel(goodL),3);
% cols(~goodL,:) = repmat([1 0 0],sum(~goodL),1);
% set(celselfig.table,'BackgroundColor',cols,'CellSelectionCallback',@cellselect);

end

function cellselect(a,b)
global GLMSPopData

% Load user data
prediffig = get(gcf,'userdata');
popan = get(prediffig.poppanel,'UserData');
indpan = get(prediffig.indivpanel,'userdata');

if  strcmp(get(a,'type'),'line') || strcmp(get(a,'type'),'axes')
    normLLdist = popan.vars.diffnormLL - b.IntersectionPoint(1);
    diffnspdist = popan.vars.maxdiffnsp - b.IntersectionPoint(2);
    [~,idx] = min(sqrt(normLLdist.^2 + diffnspdist.^2));
elseif strcmp(get(a,'type'),'uitable')
    if isempty(b.Indices)
        return
    end
    idx = b.Indices(1);
end

% Plot population data again and highlight selected point
PlotPopData()
plot(popan.vars.diffnormLL(idx),popan.vars.maxdiffnsp(idx),'b*')
cols = ones(size(GLMSPopData,1)-1,3);
cols(idx,:) = [.5 .5 1];
set(prediffig.table,'backgroundcolor',cols)

% Unpack parameters, etc.
datatypes = GLMSPopData(1,:);
GLMP = GLMSPopData{idx+1,strcmp(datatypes,'GLMP')};
sub = GLMSPopData{idx+1,strcmp(datatypes,'Subunit')};

% Plot 1D Surface
Lcc = GLMP.subunit{sub}.Lcc;
Mcc = GLMP.subunit{sub}.Mcc;
nsp = GLMP.subunit{sub}.nspikes;
x = linspace(-max(GLMP.subunit{sub}.rho),max(GLMP.subunit{sub}.rho),50);
[xx,yy] = meshgrid(x,x);
surface1d = ComputeNakaRushtonJPW(popan.vars.fitparams.oneD(idx,:),[xx(:) yy(:)],'conicsection');
surface1d = reshape(surface1d,size(xx));
axes(indpan.onedsurf); cla; hold on; grid on; axis square;
ticks = linspace(min(Lcc),max(Lcc),5);
set(gca,'XTick',ticks,'YTick',ticks);
p = surf(xx,yy,surface1d);
set(p,'edgecolor','none')
alpha(.3);
contour(xx,yy,surface1d);
h = plot3(Lcc,Mcc,nsp,'ko');
set(h,'markersize',2,'markerfacecolor','k')
xlabel('Lcc');
ylabel('Mcc');
zlabel('# of spikes')
title('1D Fit')

% Plot 2D Surface
x = linspace(-max(GLMP.subunit{sub}.rho),max(GLMP.subunit{sub}.rho),50);
[xx,yy] = meshgrid(x,x);
surface2d = ComputeNakaRushtonJPW(popan.vars.fitparams.twoD(idx,:),[xx(:) yy(:)],'conicsection');
surface2d = reshape(surface2d,size(xx));
axes(indpan.twodsurf); cla; hold on; grid on; axis square;
ticks = linspace(min(Lcc),max(Lcc),5);
set(gca,'XTick',ticks,'YTick',ticks);
p = surf(xx,yy,surface2d);
set(p,'edgecolor','none')
alpha(.3);
contour(xx,yy,surface2d);
h = plot3(Lcc,Mcc,nsp,'ko');
set(h,'markersize',2,'markerfacecolor','k')
xlabel('Lcc');
ylabel('Mcc');
zlabel('# of spikes')
title('2D Fit')

% Bubble plot
axes(indpan.bubbleplot); cla; hold on; grid on; axis square;
maxnsp = max(GLMP.subunit{sub}.meannspikes);
Lcc = GLMP.subunit{sub}.uniqueLcc;
Mcc = GLMP.subunit{sub}.uniqueMcc;
nsp = GLMP.subunit{sub}.meannspikes;
for i = 1:numel(Lcc)
    mn = nsp(i)/maxnsp*10;
    h = plot3(Lcc(i),Mcc(i),mn,'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
xlabel('Lcc');
ylabel('Mcc');
zlabel('# of spikes')
title('Bubble Plot')


% Plot bubble plot of differences between stimulus predictions
Lcc = GLMP.subunit{sub}.uniqueLcc;
Mcc = GLMP.subunit{sub}.uniqueMcc;
pred1d = ComputeNakaRushtonJPW(popan.vars.fitparams.oneD(idx,:),...
    [Lcc Mcc],'conicsection');
pred2d = ComputeNakaRushtonJPW(popan.vars.fitparams.twoD(idx,:),...
    [Lcc Mcc],'conicsection');
diff1d2d = pred2d - pred1d;
axes(indpan.diffplot); cla; hold on; grid on; axis square;
ticks = linspace(min(Lcc),max(Lcc),5);
set(gca,'XTick',ticks,'YTick',ticks);
p = surf(xx,yy,surface2d-surface1d);
set(p,'edgecolor','none')
alpha(.3);
contour(xx,yy,surface2d-surface1d);
%maxnsp = max(diff1d2d);
for i = 1:numel(Lcc)
    mn = abs(diff1d2d(i))/maxnsp*10;
    h = plot3(Lcc(i),Mcc(i),diff1d2d(i),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
xlabel('Lcc');
ylabel('Mcc');
zlabel('Difference between model predictions')
title('Difference Plot')

end



function UnpackPopData()
global GLMSPopData

% Load user data
prediffig = get(gcf,'userdata');
popan = get(prediffig.poppanel,'UserData');

% Load data
datatypes = GLMSPopData(1,:);

% Create nd index
popan.vars.nds = GLMSPopData(2:end,strcmp(datatypes,'nD'));
popan.vars.onedL = strcmp(popan.vars.nds,'1D');

% Fit Parameters
poparams = [GLMSPopData{2:end,strcmp(datatypes,'Surface Parameters')}];
onedpars = [poparams.oneD];
twodpars = [poparams.twoD];
popan.vars.fitparams.oneD = cat(1,onedpars.parvals);
popan.vars.fitparams.twoD = cat(1,twodpars.parvals);

% Modulation
popan.vars.mod = nan(size(popan.vars.onedL));
popan.vars.mod(popan.vars.onedL) = popan.vars.fitparams.oneD(popan.vars.onedL,1)...
    - popan.vars.fitparams.twoD(popan.vars.onedL,end-2);
popan.vars.mod(~popan.vars.onedL) = popan.vars.fitparams.oneD(~popan.vars.onedL,1)...
    - popan.vars.fitparams.twoD(~popan.vars.onedL,end-2);

% Tuning Std
popan.vars.PDstd = [GLMSPopData{2:end,strcmp(datatypes,'Tun Std')}]';

% Diff NormLL
popan.vars.diffnormLL = cat(1,poparams.diffnormLL);

% Maximum difference in surface prediction
popan.vars.maxdiffnsp = cat(1,poparams.maxdiffnsp);

%Save variables
set(prediffig.poppanel,'UserData',popan);
set(gcf,'UserData',prediffig);

end

function CalcMaxDiff()
global GLMSPopData

% Load user data
prediffig = get(gcf,'userdata');
popan = get(prediffig.poppanel,'UserData');

% Save maxdiffnsp for each unit
datatypes = GLMSPopData(1,:);
popan.vars.maxdiffnsp = nan(size(GLMSPopData,1)-1,1);
for n = 2:size(GLMSPopData,1)
   
    GLMP = GLMSPopData{n,strcmp(datatypes,'GLMP')};
    sub = GLMSPopData{n,strcmp(datatypes,'Subunit')};
    Lcc = GLMP.subunit{sub}.uniqueLcc;
    Mcc = GLMP.subunit{sub}.uniqueMcc;
    pred1d = ComputeNakaRushtonJPW(popan.vars.fitparams.oneD(n-1,:),...
        [Lcc Mcc],'conicsection');
    pred2d = ComputeNakaRushtonJPW(popan.vars.fitparams.twoD(n-1,:),...
        [Lcc Mcc],'conicsection');
    diff1d2d = diff([pred2d pred1d],[],2);
    %diff1d2d = diff1d2d ./ sum([pred1d pred2d],2);
    popan.vars.maxdiffnsp(n-1) = max(abs(diff1d2d));
    surfparams = GLMSPopData{n,strcmp(datatypes,'Surface Parameters')};
    surfparams.maxdiffnsp = max(abs(diff1d2d));
    GLMSPopData{n,strcmp(datatypes,'Surface Parameters')} = surfparams;
    
end

% Grab saved population data
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
save([library 'GLMSPopData.mat'],'GLMSPopData')

%Save variables
set(prediffig.poppanel,'UserData',popan);
set(gcf,'UserData',prediffig);

end

