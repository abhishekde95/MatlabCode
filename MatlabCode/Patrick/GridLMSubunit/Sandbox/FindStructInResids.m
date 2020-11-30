function FindStructInResids()
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
params1 = fitspanel.params.oneD;

% Calculate predicted responses
predresp = ComputeNakaRushtonJPW(params1,[Lcc Mcc],'surface7');
resids = nsp-predresp;

% Try using only the points within 20 deg of the 2 directions
prefdir1 = params1(end);
if prefdir1 > 0
    prefdir2 = prefdir1 - pi;
else
    prefdir2 = prefdir1 + pi;
end
orthodir1 = prefdir1 - pi/2;
orthodir2 = prefdir1 + pi/2;
tol = 20/180*pi;
colinptsL1 = GLMP.subunit{p}.theta > (prefdir1-tol) & GLMP.subunit{p}.theta < (prefdir1+tol);
colinptsL2 = GLMP.subunit{p}.theta > (prefdir2-tol) & GLMP.subunit{p}.theta < (prefdir2+tol);
colinptsL = colinptsL1 | colinptsL2;
orthoptsL1 = GLMP.subunit{p}.theta > (orthodir1-tol) & GLMP.subunit{p}.theta < (orthodir1+tol);
orthoptsL2 = GLMP.subunit{p}.theta > (orthodir2-tol) & GLMP.subunit{p}.theta < (orthodir2+tol);
orthoptsL = orthoptsL1 | orthoptsL2;
colinptsL = true(size(colinptsL));
orthoptsL = colinptsL;

% Plot for sanity check
figure(937690); clf; hold on;
plot(GLMP.subunit{p}.Lcc,GLMP.subunit{p}.Mcc,'ok')
plot(GLMP.subunit{p}.Lcc(colinptsL),GLMP.subunit{p}.Mcc(colinptsL),'r*')
plot(GLMP.subunit{p}.Lcc(orthoptsL),GLMP.subunit{p}.Mcc(orthoptsL),'g*')

% Project data onto unit vector orthogonal to preferred direction
prefdir = params1(end);
[x,y] = pol2cart(prefdir,1);
unitvect = [x y];
%colinpts = [Lcc Mcc] * unitvect';
colinpts = [Lcc(colinptsL) Mcc(colinptsL)] * unitvect';
orthunitvect = unitvect * [cos(pi/2) sin(pi/2); -sin(pi/2) cos(pi/2)];
%orthopts = [Lcc Mcc] * orthunitvect';
orthopts = [Lcc(orthoptsL) Mcc(orthoptsL)] * orthunitvect';

% Create baseline surf
x = linspace(min(Lcc),max(Lcc),100);
[xx,yy] = meshgrid(x,x);
zerosurf = zeros(size(xx));

% Plot data
figure(5439406); clf; hold on; grid on;
plot3(Lcc,Mcc,nsp,'k*')
plot3(Lcc,Mcc,predresp,'ro')
xlabel('Lcc'); ylabel('Mcc'); zlabel('Response (#spikes)')

figure(83758); clf; hold on; grid on;
stem3(Lcc,Mcc,resids)
surf(xx,yy,zerosurf)
xlabel('Lcc'); ylabel('Mcc'); zlabel('Residuals')

figure(59493); clf;
set(gcf,'position',[200 100 1000 700])
colinaxes = axes('parent',gcf,'units','normalized','pos',[.1 .55 .8 .35],'box','on');
hold on; grid on;
stem(colinpts,resids(colinptsL))
[curve, goodness, output] = fit(colinpts,resids(colinptsL),'smoothingspline');
[curve, goodness, output] = fit(colinpts,resids(colinptsL),'poly2');
%[curve, goodness, output] = fit(colinpts,resids(colinptsL),'exp1');
plot(curve,colinpts,resids(colinptsL));
xlabel('L/M Stim prjected onto Preferred Axis')
ylabel('Residuals')


orthoaxes = axes('parent',gcf,'units','normalized','pos',[.1 .1 .8 .35],'box','on');
hold on; grid on;
stem(orthopts,resids(orthoptsL))
[curve, goodness, output] = fit(orthopts,resids(orthoptsL),'smoothingspline');
[curve, goodness, output] = fit(orthopts,resids(orthoptsL),'poly2');
%[curve, goodness, output] = fit(orthopts,resids(orthoptsL),'exp1');
plot(curve,orthopts,resids(orthoptsL));
xlabel('L/M Stim prjected onto Orthogonal Axis')
ylabel('Residuals')



figure(83586); clf; hold on; grid on;
F = TriScatteredInterp(Lcc,Mcc,resids);
[qx,qy] = meshgrid(min(Lcc):.01:max(Lcc),min(Mcc):.01:max(Mcc));
plot3(Lcc,Mcc,resids,'k*')
surfc(qx,qy,F(qx,qy))
xlabel('Lcc'); ylabel('Mcc'); zlabel('Residuals')
% 
% % Fit residuals with two functions: a line and an exponential (all Poisson error)
% figure(5436543); clf; hold on; grid on;
% markSize = (((resids - min(resids)) ./ max(resids-min(resids)))) + .1;
% markSize = markSize * 25;
% for t = 1:numel(Lcc)
%     plot(Lcc(t),Mcc(t),'ko','MarkerSize',markSize(t))
% end
keyboard






