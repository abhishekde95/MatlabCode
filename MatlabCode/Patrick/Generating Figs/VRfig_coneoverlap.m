% Making an additional figure for the VR paper

freq = 0:1:1000;
lim = [-.5 .5];
gun1freq = round(max(freq)/2.4);
gun2freq = max(freq) - round(max(freq)/2.4);
sig = max(freq)/10;
gun1idx = find(freq == gun1freq);
gun2idx = find(freq == gun2freq);
gun1act = normrnd(0,.2,[1000 1]);
gun2act = normrnd(0,.2,[1000 1]);
while any(gun1act < lim(1)) || any(gun1act > lim(2))
    L = gun1act < lim(1) | gun1act > lim(2);
    gun1act(L) = normrnd(0,.2,[sum(L) 1]);
end
while any(gun2act < lim(1)) || any(gun2act > lim(2))
    L = gun2act < lim(1) | gun2act > lim(2);
    gun2act(L) = normrnd(.5,.2,[sum(L) 1]);
end
boxedges = [lim(1) lim(1); lim(1) lim(2); lim(2) lim(2); lim(2) lim(1); lim(1) lim(1)];



%% Nonoverlapping cone absorption spectra
mu1 = max(freq)/4;
mu2 = max(freq) - max(freq)/4;
cone1 = normpdf(freq,mu1,sig);
cone1 = cone1 ./ max(cone1);
cone2 = normpdf(freq,mu2,sig);
cone2 = cone2 ./ max(cone2);

% Plot cone absorption spectra
figure(10); clf; hold on; box on;
xlabel('Wavelength')
ylabel('Probability of Absorption')
plot(freq,cone1,'g')
plot(freq,cone2,'r')
plot([gun1idx gun1idx],[0 1],'g--')
plot([gun2idx gun2idx],[0 1],'r--')

%Format and save figure
set(gcf,'PaperPositionMode','auto')
if ismac
    print('/Users/jpatrickweller/Dropbox/VisionResearchPaper/6a','-depsc')
else ispc
    print('C:\Users\jpweller\Dropbox\VisionResearchPaper\6a','-depsc');
end
disp('Fig 6a done.')


% Build gamut box
cone1box = cone1(gun1idx) .* boxedges(:,1) + cone1(gun2idx) .* boxedges(:,2);
cone2box = cone2(gun1idx) .* boxedges(:,1) + cone2(gun2idx) .* boxedges(:,2);

% Gaussian white noise distribution
conemod1 = (cone1(gun1idx) .* gun1act) + (cone1(gun2idx) .* gun2act);
conemod2 = (cone2(gun1idx) .* gun1act) + (cone2(gun2idx) .* gun2act);

% White Noise cone mod distribution
figure(11); clf; hold on; box on; axis square;
xlabel('Cone 2 Activity')
ylabel('Cone 1 Activity')
h = plot(conemod1,conemod2,'o');
%set(h,'MarkerFaceColor','k','MarkerEdgeColor','w')
plot(cone1box,cone2box,'r--')

%Format and save figure
set(gcf,'PaperPositionMode','auto')
if ismac
    print('/Users/jpatrickweller/Dropbox/VisionResearchPaper/6b','-depsc')
else ispc
    print('C:\Users\jpweller\Dropbox\VisionResearchPaper\6b','-depsc');
end
disp('Fig 6b done.')


%% Partially overlapping cone absorption spectra
mu1 = max(freq)/2.4;
mu2 = max(freq) - max(freq)/2.4;
cone1 = normpdf(freq,mu1,sig);
cone1 = cone1 ./ max(cone1);
cone2 = normpdf(freq,mu2,sig);
cone2 = cone2 ./ max(cone2);

% Plot cone absorption spectra
figure(20); clf; hold on; box on;
xlabel('Wavelength')
ylabel('Probability of Absorption')
plot(freq,cone1,'g')
plot(freq,cone2,'r')
plot([gun1idx gun1idx],[0 1],'g--')
plot([gun2idx gun2idx],[0 1],'r--')

%Format and save figure
set(gcf,'PaperPositionMode','auto')
if ismac
    print('/Users/jpatrickweller/Dropbox/VisionResearchPaper/6c','-depsc')
else ispc
    print('C:\Users\jpweller\Dropbox\VisionResearchPaper\6c','-depsc');
end
disp('Fig 6c done.')

% Build gamut box
cone1box = cone1(gun1idx) .* boxedges(:,1) + cone1(gun2idx) .* boxedges(:,2);
cone2box = cone2(gun1idx) .* boxedges(:,1) + cone2(gun2idx) .* boxedges(:,2);

% Gaussian white noise distribution
conemod1 = (cone1(gun1idx) .* gun1act) + (cone1(gun2idx) .* gun2act);
conemod2 = (cone2(gun1idx) .* gun1act) + (cone2(gun2idx) .* gun2act);

% White Noise cone mod distribution
figure(21); clf; hold on; box on; axis square;
xlabel('Cone 2 Activity')
ylabel('Cone 1 Activity')
h = plot(conemod1,conemod2,'o');
%set(h,'MarkerFaceColor','k','MarkerEdgeColor','w')
plot(cone1box,cone2box,'r--')

%Format and save figure
set(gcf,'PaperPositionMode','auto')
if ismac
    print('/Users/jpatrickweller/Dropbox/VisionResearchPaper/6d','-depsc')
else ispc
    print('C:\Users\jpweller\Dropbox\VisionResearchPaper\6d','-depsc');
end
disp('Fig 6d done.')


%% Entirely overlapping cone absorption spectra
mu1 = max(freq)/2.03;
mu2 = max(freq) - max(freq)/2.03;
cone1 = normpdf(freq,mu1,sig);
cone1 = cone1 ./ max(cone1);
cone2 = normpdf(freq,mu2,sig);
cone2 = cone2 ./ max(cone2);

% Plot cone absorption spectra
figure(30); clf; hold on; box on;
xlabel('Wavelength')
ylabel('Probability of Absorption')
plot(freq,cone1,'g')
plot(freq,cone2,'r')
plot([gun1idx gun1idx],[0 1],'g--')
plot([gun2idx gun2idx],[0 1],'r--')

%Format and save figure
set(gcf,'PaperPositionMode','auto')
if ismac
    print('/Users/jpatrickweller/Dropbox/VisionResearchPaper/6e','-depsc')
else ispc
    print('C:\Users\jpweller\Dropbox\VisionResearchPaper\6e','-depsc');
end
disp('Fig 6e done.')

% Gaussian white noise distribution
conemod1 = (cone1(gun1idx) .* gun1act) + (cone1(gun2idx) .* gun2act);
conemod2 = (cone2(gun1idx) .* gun1act) + (cone2(gun2idx) .* gun2act);

% Build gamut box
cone1box = cone1(gun1idx) .* boxedges(:,1) + cone1(gun2idx) .* boxedges(:,2);
cone2box = cone2(gun1idx) .* boxedges(:,1) + cone2(gun2idx) .* boxedges(:,2);

% White Noise cone mod distribution
figure(31); clf; hold on; axis square; box on;
xlabel('Cone 2 Activity')
ylabel('Cone 1 Activity')
h = plot(conemod1,conemod2,'o');
%set(h,'MarkerFaceColor','k','MarkerEdgeColor','w')
plot(cone1box,cone2box,'r--')

%Format and save figure
set(gcf,'PaperPositionMode','auto')
if ismac
    print('/Users/jpatrickweller/Dropbox/VisionResearchPaper/6f','-depsc')
else ispc
    print('C:\Users\jpweller\Dropbox\VisionResearchPaper\6f','-depsc');
end
disp('Fig 6f done.')


