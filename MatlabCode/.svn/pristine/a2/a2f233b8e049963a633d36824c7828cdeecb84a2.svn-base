% This script is for generating figure for a Vision Research paper

% Pre-Defined Variables
thetaspace = pi/4;
rhospace = .5;
gl.nRnds = 3;
gl.nPres = 5;
gl.lummaxcc = .9;
gl.colmaxcc = .09;

% Construct Polar Grid
if mod(gl.nRnds,2) == 0
    thetaspace = thetaspace / 2^((gl.nRnds-2)/2);
    rhospace = rhospace / 2^(gl.nRnds/2);
elseif mod(gl.nRnds,2) == 1
    thetaspace = thetaspace / 2^((gl.nRnds-1)/2);
    rhospace = rhospace / 2^((gl.nRnds-1)/2);
end
thetas = shiftdim(0:thetaspace:2*pi-thetaspace);
rhos = shiftdim(rhospace:rhospace:1);

% Ennumerate all conditions
PolRhoIdx = fullfact([numel(thetas) numel(rhos)]);
rhothetalist = repmat([thetas(PolRhoIdx(:,1)) rhos(PolRhoIdx(:,2))],gl.nPres,1);
thetas = rhothetalist(:,1);
rhos = rhothetalist(:,2);

% Transform from polar to cartesian coordinates
[tempLcc,tempMcc] = pol2cart(thetas,rhos);

% Scale Cone Contrast Units for Monitor
scale = gl.lummaxcc*gl.colmaxcc./sqrt((gl.colmaxcc.*cos(thetas-pi/4)).^2 ...
    +(gl.lummaxcc.*sin(thetas-pi/4)).^2);
Lcc = tempLcc .* scale;
Mcc = tempMcc .* scale;

% Retrieve stimuli and parameters
ang = gl.allAngs(gl.currentAng);
%L = surfpanel.angle(gl.currentAng).samp(gl.currentSamp).Lcc;
%M = surfpanel.angle(gl.currentAng).samp(gl.currentSamp).Mcc;
L = surfpanel.Lcc;
M = surfpanel.Mcc;

%rotMat = [cos(ang) sin(ang); -sin(ang) cos(ang)];

% Set up some variables
params = [50 100 .3 3 .1 ang];

% Responses to stimuli
nsp = ComputeNakaRushtonJPW(params,[L M],'surface7');
nsp = poissrnd(nsp);

%surfpanel.angle(gl.currentAng).samp(gl.currentSamp).Lcc = L;
%surfpanel.angle(gl.currentAng).samp(gl.currentSamp).Mcc = M;
%surfpanel.angle(gl.currentAng).samp(gl.currentSamp).responses = nsp;
surfpanel.nsp = nsp;

% Plot figure
axes(surfpanel.axes); cla; hold on;
plot3(L,M,nsp,'ko')
xlabel('Lcc'); ylabel('Mcc'); title('Model Surface')
axis equal square
