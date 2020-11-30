% script for generating projection figure

global GLMP
datafile = GLMP.datafile;


sub = 1;
Lcc = GLMP.subunit{sub}.uniqueLcc;
Mcc = GLMP.subunit{sub}.uniqueMcc;
meannsp = GLMP.subunit{sub}.meannspikes;
stderr = GLMP.subunit{sub}.stderr;

surffig = get(60,'userdata');
fitspanel = get(surffig.fitspanel,'userdata');
params1 = fitspanel.params.twoD;

rot = params1(end-1);
rotMat = [cos(rot) sin(rot); -sin(rot) cos(rot)];
proj0 = [Lcc Mcc] * rotMat;


% Axis 1
proj = proj0(:,1);
params2 = [params1(1) params1(1) params1(3) params1(5) params1(6) params1(6) params1(7)];
fitpts = linspace(min(proj),max(proj),100);
fit = ComputeNakaRushtonJPW(params2,fitpts,'asymmetric');

figure(2563); clf; hold on; grid on
plot(-proj,meannsp,'k*')
plot(fitpts,fit,'r-')
errorbar(-proj,meannsp,stderr,'*k')
ylim([0 max(meannsp)*1.01])
box on;
name = strcat('/Users/jpatrickweller/Documents/Matlab Working Directory/Patrick/GridLMSubunit/Figures/',GLMP.datafile,'proj1');
export_fig((name),'-eps','-depsc','-transparent')

% Axis 2
proj = proj0(:,2);
params3 = [params1(1) params1(1) params1(4) params1(4) params1(6) params1(6) params1(7)];
fitpts = linspace(min(proj),max(proj),100);
fit = ComputeNakaRushtonJPW(params3,fitpts,'asymmetric');

figure(2564); clf; hold on; grid on
plot(proj,meannsp,'k*')
plot(-fitpts,fit,'r-')
errorbar(proj,meannsp,stderr,'*k')
ylim([0 max(meannsp)*1.01])
box on;
name = strcat('/Users/jpatrickweller/Documents/Matlab Working Directory/Patrick/GridLMSubunit/Figures/',GLMP.datafile,'proj2');
export_fig((name),'-eps','-depsc','-transparent')