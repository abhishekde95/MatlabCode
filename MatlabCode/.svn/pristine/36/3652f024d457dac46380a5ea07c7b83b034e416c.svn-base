function Compare1Dand2DSurf()
global GLMP

% Load Figure Variables
surffig = get(60,'UserData');
conpanel = get(surffig.conpanel,'UserData');
fitspanel = get(surffig.fitspanel,'UserData');

% Load data parameters
p = conpanel.subselect;
Lcc = GLMP.subunit{p}.Lcc;
Mcc = GLMP.subunit{p}.Mcc;
nsp = GLMP.subunit{p}.nspikes;
params1D = fitspanel.params.oneD;
params2D = fitspanel.params.twoD;

% Calculate predicted responses
predresp1D = ComputeNakaRushtonJPW(params1D,[Lcc Mcc],'surface7');
predresp2D = ComputeNakaRushtonJPW(params2D,[Lcc Mcc],'surface8');

% Find error with each fit
error1D = nsp - predresp1D;
error2D = nsp - predresp2D;

% plot
figure(61); clf;
set(gcf,'position',[100 100 1200 600])
ax1 = axes('parent',gcf,'units','normalized','pos',[.05 .05 .26 .9]);
stem3(Lcc,Mcc,predresp2D-predresp1D)
xlabel('Lcc'); ylabel('Mcc'); zlabel('2D-1D')
ax2 = axes('parent',gcf,'units','normalized','pos',[.36 .05 .26 .9]);
stem3(Lcc,Mcc,error1D)
xlabel('Lcc'); ylabel('Mcc'); zlabel('Error in 1D Fit')
ax3 = axes('parent',gcf,'units','normalized','pos',[.67 .05 .26 .9]);
stem3(Lcc,Mcc,error2D)
xlabel('Lcc'); ylabel('Mcc'); zlabel('Error in 2D Fit')
keyboard
end