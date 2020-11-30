% Building figures for a presentation


%% Gabor movie
x = linspace(0,pi*6,15);
[xx,yy] = meshgrid(x,x);
ntpts = 100;
z = cos(xx);
sigma = [2*pi 0; 0 2*pi];
gauss = mvnpdf([xx(:) yy(:)],mean(x),sigma);
gauss = gauss./max(gauss);
gauss = reshape(gauss,size(xx));
RGBs = nan(size(xx,1),size(xx,2),3,ntpts);

for n = 1:ntpts
    gabor = z.*gauss;
    temp = repmat(gabor,1,1,3);
    temp = (temp+1)/2;
    RGBs(:,:,:,n) = temp;
    z = circshift(z,1,2);    
end
    
mov = immovie(RGBs);
implay(mov)

gabmov = VideoWriter('gabormovie');
open(gabmov)
writeVideo(gabmov,mov)
close(gabmov)

%% White nosie movie
global DN

if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/nex files/';
elseif ispc
    library = 'C:\Documents and Settings\JPatrickWeller\My Documents\Dropbox\Patrick\GLMS Data\nex files\';
end
rawdata = nex2stro([char(library) DN.datafile '.nex']);


lumCC = sqrt((max(DN.lumCC)^2)/2);
colCC = sqrt((max(DN.colCC)^2)/2);
ccs = [lumCC lumCC 0; -lumCC -lumCC 0; colCC -colCC 0; -colCC colCC 0];
nstix = max(DN.NStixGrid);
pixperstix = 15;
nframes = 1;
randorder = randi(4,nstix^2,nframes);
bkgndrgb = rawdata.sum.exptParams.bkgndrgb';
bkgndlms =  bkgndrgb * DN.M;
stimlms = ccs .* repmat(bkgndlms,size(ccs,1),1) + repmat(bkgndlms,size(ccs,1),1);
stimrgb = (inv(DN.M)*stimlms')'% * ;(inv(M)*LMSimage')'

% from slave
% LccLum = sqrt((max(DN.lumCC)^2)/2);
% MccLum = LccLum;
% LccCol = sqrt((colcc^2)/2);
% MccCol = LccCol;
% cc = [LccLum MccLum Scc; -LccLum -MccLum Scc; LccCol -MccCol Scc; -LccCol MccCol Scc]
% gammaTable = reshape(DN.gammatable,[numel(DN.gammatable)/3 3])
% 
% gl.stim.mu = [gammaTable(bkgndrgb(1)+1,1),...
%     gl.cal.gammaTable(gl.bkgndRGB(2)+1,2),...
%     gl.cal.gammaTable(gl.bkgndRGB(3)+1,3)];
% mulms = gl.cal.M*gl.stim.mu';
% stimlms = cc.*repmat(mulms',size(cc,1),1) + repmat(mulms',size(cc,1),1);
% rgbmat = gl.cal.invM*stimlms'



WNmov = nan(nstix*pixperstix,nstix*pixperstix,3,nframes);
for n = 1:nframes
    temp = stimrgb(randorder(:,n),:);
    temp = reshape(temp,[nstix nstix 3]);
    for m = 1:nstix
        for p = 1:nstix
            Rframe = repmat(temp(m,p,1),[pixperstix pixperstix]);
            Gframe = repmat(temp(m,p,2),[pixperstix pixperstix]);
            Bframe = repmat(temp(m,p,3),[pixperstix pixperstix]);
            WNmov((pixperstix*(m-1)+1):(pixperstix*m),(pixperstix*(p-1)+1):(pixperstix*p),:,n) =...
                cat(3,Rframe,Gframe,Bframe);
        end
    end
    %WNmov(:,:,:,n) = reshape(temp,[nstix,nstix,3]);
end

figure(1); clf; image(WNmov(:,:,:,n));

% mov = immovie(WNmov);
% implay(mov);

% whitenoisemov = VideoWriter('whitenoisemovie');
% open(whitenoisemov);
% writeVideo(whitenoisemov,mov);
% close(whitenoisemov);


%% Pulling out WN images
global DN


filterfig = get(160,'userdata');
fp = get(filterfig.filterpanel,'userdata');
dnp = get(filterfig.dnpanel,'userdata')';
    
% Formatting
f = figure(1); clf;
set(f,'units','normalized','pos',[.1 .5 .6 .4])
%p = copyobj(dnp.axes(:,:),f);
p = copyobj(fp.axes(2,:),f);
set(gcf,'PaperPositionMode','auto')
name = [DN.datafile ' PSTH'];
export_fig(['/Users/jpatrickweller/Documents/Committee Meetings/Lab Presentation/' name],'-depsc')

%% Pull out frames from overview
global DN

cf = get(279,'userdata');
dn = get(cf.dnpanel,'userdata');
ax = 8;
figure(1); clf;
p = copyobj(dn.axes(1,ax),1);
set(p,'xlabel',[],'ylabel',[],'title',[],'units','normalized','pos',[.05 .55 .4 .4])
axis(p,'equal','tight')
p = copyobj(dn.axes(2,ax),1);
set(p,'xlabel',[],'ylabel',[],'title',[],'units','normalized','pos',[.55 .05 .4 .4])
axis(p,'equal','tight')
p = copyobj(dn.axes(3,ax),1);
set(p,'xlabel',[],'ylabel',[],'title',[],'units','normalized','pos',[.55 .55 .4 .4])
axis(p,'equal','tight')
p = copyobj(dn.axes(4,ax),1);
set(p,'xlabel',[],'ylabel',[],'title',[],'units','normalized','pos',[.05 .05 .4 .4])
axis(p,'equal','tight')
set(gcf,'color',[1 1 1])
name = [DN.datafile ' WN sta'];
export_fig(['/Users/jpatrickweller/Documents/Committee Meetings/Lab Presentation/' name],'-depsc')


%% Pull out Bubbleplot

cf = get(279,'userdata');
bp = get(cf.glmppanel,'userdata');
sub = 1;
f = figure(1); clf;
p = copyobj(bp.axes.sub(sub),1);
set(p,'xlabel',[],'ylabel',[],'title',[],'units','normalized','pos',[.1 .1 .8 .8])
set(gcf,'color',[1 1 1])
name = [DN.datafile ' s' num2str(sub) ' bubbleplot'];
export_fig(['/Users/jpatrickweller/Documents/Committee Meetings/Lab Presentation/' name],'-depsc')


%% Pull out pop histograms

tunfig = get(15,'userdata');
popan = get(tunfig.popanel,'userdata');
figure(1); clf;
p = copyobj(popan.axes.oneD,1);
set(p,'xlabel',[],'ylabel',[],'title',[],'units','normalized','pos',[.1 .1 .8 .8])
name = 'poptun_1d_unique'
print(name,'-depsc')




