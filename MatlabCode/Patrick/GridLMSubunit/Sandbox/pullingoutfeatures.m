function pullingoutfeatures()
global GLMSPopData

% Grab saved population data
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
load([library 'GLMSPopData.mat'])




% Which features are you looking for?
%DNgridsizes()
%TestedContrasts()
%WNDVAs()
%RFlocations()
%eachmonk()
nstimsamples()

end

function nstimsamples()
global GLMSPopData

datatypes = GLMSPopData(1,:);
GLMPs = GLMSPopData(2:end,strcmp(datatypes,'GLMP'));
subs = [GLMSPopData{2:end,strcmp(datatypes,'Subunit')}]';
nobs = nan(size(subs));

for n = 1:size(GLMPs,1)
    nobs(n) = size(GLMPs{n}.subunit{subs(n)}.spiketimes_col,2);
end

end

function eachmonk()
global GLMSPopData

datatypes = GLMSPopData(1,:);
dfs = cat(1,GLMSPopData{2:end,strcmp(datatypes,'Datafile')});

[a,b,c] = unique(dfs,'rows');
idx = false(size(c));
idx(b) = true;
twosubs = dfs(~idx,:);

% number of 2 sub neurons for each monkey
ndubs = sum(twosubs(:,1) == 'N')
mdubs = sum(twosubs(:,1) == 'M')

% number of 1 sub neurons for each monkey
sum(a(:,1) == 'N') - ndubs
sum(a(:,1) == 'M') - mdubs

end

function RFlocations()
global GLMSPopData

datatypes = GLMSPopData(1,:);

rfs = [GLMSPopData{2:end,strcmp(datatypes,'RF Params')}]';
rhos = [rfs.rfrho]';

figure(5); clf; hold on;
set(gcf,'name','DVAs')
plot(rhos,'k*')
xlabel('Datafile')
ylabel('DVA per stixel')

end


function WNDVAs()
global GLMSPopData

datatypes = GLMSPopData(1,:);
datafiles = GLMSPopData(2:end,strcmp(datatypes,'Datafile'));
[uniquedfs,uniqueidx] = unique(datafiles);
dvas = nan(numel(uniquedfs),1);
for n = 1:numel(uniqueidx)
    DN = GLMSPopData{uniqueidx(n)+1,strcmp(datatypes,'DN')};    
    dvas(n) = DN.DVAPerStix(1);
end

figure(4); clf; hold on;
set(gcf,'name','DVAs')
plot(uniqueidx,dvas,'k*')
xlabel('Datafile')
ylabel('DVA per stixel')

uniconts = unique(dvas);
for n = 1:size(uniconts,1)
    total = sum(dvas == uniconts(n));
    disp(['Neurons tested with DN DVA/stixel = ' num2str(uniconts(n)) ' = ' num2str(total)])
end

end



function TestedContrasts()
global GLMSPopData

datatypes = GLMSPopData(1,:);
datafiles = GLMSPopData(2:end,strcmp(datatypes,'Datafile'));
[uniquedfs,uniqueidx] = unique(datafiles);
maxcolDN = nan(numel(uniquedfs),1);
maxlumDN = nan(numel(uniquedfs),1);
maxcolGLMP = nan(numel(uniquedfs),1);
maxlumGLMP = nan(numel(uniquedfs),1);
for n = 1:numel(uniqueidx)
    DN = GLMSPopData{uniqueidx(n)+1,strcmp(datatypes,'DN')};    
    maxcolDN(n) = max(DN.colCC);
    maxlumDN(n) = max(DN.lumCC);
    GLMP = GLMSPopData{uniqueidx(n)+1,strcmp(datatypes,'GLMP')};
    maxcolGLMP(n) = max(GLMP.rho(GLMP.theta == -pi/4));
    maxlumGLMP(n) = max(GLMP.rho(GLMP.theta == pi/4));
end


figure(2); clf; hold on; box on;
set(gcf,'name','Max contrasts used in DN/GLMP')
plot(maxcolDN,maxcolGLMP,'r*')
plot(maxlumDN,maxlumGLMP,'ko')
plot([0 1],[0 1],'k--')
legend('Chromatic','Luminance')
xlabel('DN Contrast')
ylabel('GLMP Contrast')

uniconts = unique([maxcolDN maxlumDN],'rows');
for n = 1:size(uniconts,1)
    total = sum(maxcolDN == uniconts(n,1) & maxlumDN == uniconts(n,2));
    disp(['Neurons tested with DN lum = ' num2str(uniconts(n,1)) ' and col = ' num2str(uniconts(n,2)) ' = ' num2str(total)])
end
uniconts = unique([maxcolGLMP maxlumGLMP],'rows');
for n = 1:size(uniconts,1)
    total = sum(maxcolGLMP == uniconts(n,1) & maxlumGLMP == uniconts(n,2));
    disp(['Neurons tested with GLMP lum = ' num2str(uniconts(n,1)) ' and col = ' num2str(uniconts(n,2)) ' = ' num2str(total)])
end


end

function DNgridsizes()
global GLMSPopData

datatypes = GLMSPopData(1,:);
datafiles = GLMSPopData(2:end,strcmp(datatypes,'Datafile'));
[uniquedfs,uniqueidx] = unique(datafiles);
dnsize = nan(numel(uniquedfs),1);
for n = 1:numel(uniqueidx)
    %GLMP = GLMSPopData{uniqueidx(n)+1,strcmp(datatypes,'GLMP')};
    DN = GLMSPopData{uniqueidx(n)+1,strcmp(datatypes,'DN')};
    dnsize(n) = DN.NStixGrid(1);
end

figure(1); clf; hold on;
set(gcf,'name','Number of stixels in WN grid')
plot(uniqueidx,dnsize,'ko')
xlabel('Datafile number')
ylabel('NStixelGrid')

uniquesizes = unique(dnsize);
for n = 1:numel(uniquesizes)
    disp(['Neurons tested with grid size ' num2str(uniquesizes(n)) 'x' num2str(uniquesizes(n)) ' = ' num2str(sum(dnsize == uniquesizes(n)))])
end

end
