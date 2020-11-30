% This script emulates the stimulus selection proceedure from GLMS.

% User defined variables
colmaxcc = .09;
lummaxcc = .7;
nRnds = 3;
nRpt = 5;

% Non-user-defined variables
thetaspace = pi/4;
rhospace = .5;

% Construct Polar Grid
n = floor((nRnds)/2);
rhospace = rhospace/(2^n);
if n>0 
    thetaspace = thetaspace/(2^(n-1));
    if mod(nRnds,2)
        thetaspace = thetaspace/2;
    end
end
thetavals = shiftdim(thetaspace:thetaspace:2*pi);
rhovals = shiftdim(rhospace:rhospace:1);

% Ennumerate all conditions
polrhoIdx = fullfact([numel(thetavals) numel(rhovals)]);
thetarholist = [thetavals(polrhoIdx(:,1)) rhovals(polrhoIdx(:,2))];
thetas = thetarholist(:,1);
rhos = thetarholist(:,2);

% Transform from polar to cartesian coordinates
[tempLcc,tempMcc] = pol2cart(thetas,rhos);

% Scale Cone Contrast Units for Monitor
scale = lummaxcc*colmaxcc./sqrt((colmaxcc.*cos(thetas-pi/4)).^2 ...
    +(lummaxcc.*sin(thetas-pi/4)).^2);
Lcc = tempLcc .* scale;
Mcc = tempMcc .* scale;
stim = [Lcc Mcc];
stim = repmat(stim,nRpt,1);

% Plot figure
figure(1); clf;
plot(Lcc,Mcc,'ko')
axis square

% Save data
%save('StimListJPW','stim')
