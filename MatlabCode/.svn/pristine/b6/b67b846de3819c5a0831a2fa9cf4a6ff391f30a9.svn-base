%% Sigma

x = 0:.01:1;
y = zeros(size(x));
sig1 = .3;
sig2 = .5;
sig3 = .7;
exp = 5;
A = 40;
bl = 0.1;

params1 = [A sig1 sig1 exp bl 0];
params2 = [A sig2 sig2 exp bl 0];
params3 = [A sig3 sig3 exp bl 0];


z1 = ComputeNakaRushtonJPW(params1,[x' y'],'surface7');
z2 = ComputeNakaRushtonJPW(params2,[x' y'],'surface7');
z3 = ComputeNakaRushtonJPW(params3,[x' y'],'surface7');

figure(50); clf; hold on;
plot(x,z1,'r-')
plot(x,z2,'g-')
plot(x,z3,'b-')

set(gcf,'PaperPositionMode','auto')
print('-depsc','NRSigmas')

%% Exponent

x = 0:.01:1;
y = zeros(size(x));
sig = .5;
exp1 = 3;
exp2 = 5;
exp3 = 8;
A = 40;
bl = 0.1;

params1 = [A sig sig exp1 bl 0];
params2 = [A sig sig exp2 bl 0];
params3 = [A sig sig exp3 bl 0];


z1 = ComputeNakaRushtonJPW(params1,[x' y'],'surface7');
z2 = ComputeNakaRushtonJPW(params2,[x' y'],'surface7');
z3 = ComputeNakaRushtonJPW(params3,[x' y'],'surface7');

figure(50); clf; hold on;
plot(x,z1,'r-')
plot(x,z2,'g-')
plot(x,z3,'b-')

set(gcf,'PaperPositionMode','auto')
print('-depsc','NRExp')


%% Upper Asymptote

x = 0:.01:1;
y = zeros(size(x));
sig = .5;
exp = 5;
A1 = 20;
A2 = 40;
A3 = 60;
bl = 0.1;

params1 = [A1 sig sig exp bl 0];
params2 = [A2 sig sig exp bl 0];
params3 = [A3 sig sig exp bl 0];


z1 = ComputeNakaRushtonJPW(params1,[x' y'],'surface7');
z2 = ComputeNakaRushtonJPW(params2,[x' y'],'surface7');
z3 = ComputeNakaRushtonJPW(params3,[x' y'],'surface7');

figure(50); clf; hold on;
plot(x,z1,'r-')
plot(x,z2,'g-')
plot(x,z3,'b-')

set(gcf,'PaperPositionMode','auto')
print('-depsc','NRUpperA')

%% Baseline

x = 0:.01:1;
y = zeros(size(x));
sig = .5;
exp = 5;
A = 40;
bl1 = 0.1;
bl2 = 5;
bl3 = 10;
params1 = [A sig sig exp bl1 0];
params2 = [A sig sig exp bl2 0];
params3 = [A sig sig exp bl3 0];


z1 = ComputeNakaRushtonJPW(params1,[x' y'],'surface7');
z2 = ComputeNakaRushtonJPW(params2,[x' y'],'surface7');
z3 = ComputeNakaRushtonJPW(params3,[x' y'],'surface7');

figure(50); clf; hold on;
plot(x,z1,'r-')
plot(x,z2,'g-')
plot(x,z3,'b-')

set(gcf,'PaperPositionMode','auto')
print('-depsc','NRBaseline')


%% Plotting DN difference maps

figure(50); clf;
a = get(60,'children');
b = get(a(3),'children');
copyobj(b(1),50);
set(gcf,'PaperPositionMode','auto','pos',[300 300 800 400])
name = get(60,'Name');
axis equal
width = 5.5;         % Initialize a variable for width.
height = 3;
papersize = get(gcf,'papersize')
left = (papersize(1)-width)/2
bottom = (papersize(2)-height)/2
myfiguresize = [left,bottom,width,height];
set(gcf,'PaperPosition',myfiguresize);
print('-depsc',name)

