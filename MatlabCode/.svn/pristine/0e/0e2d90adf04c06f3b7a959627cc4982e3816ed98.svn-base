% Another attempt to model a 2d surface


% NakaRushton variables
topfr1 = 50;
topfr2 = 50;
topfr3 = 50;
topfr4 = 50;
baseline = .0001;
exp1 = 2;
exp2 = 2;
exp3 = 2;
exp4 = 2;
sigma1 = .5;
sigma2 = .50;
sigma3 = .5;
sigma4 = .50;
params0 = [topfr1, topfr3, sigma1, sigma3, exp1, exp3, baseline];
params1 = [topfr2, topfr4, sigma2, sigma4, exp2, exp4, baseline];
x = -1:.05:1;

fitPts0 = ComputeNakaRushtonJPW(params0,x,'asymmetric');
fitPts1 = ComputeNakaRushtonJPW(params1,x,'asymmetric');


figure(1); clf; hold on; grid on;
plot(x,fitPts0,'r-')
plot(x,fitPts1,'b-')

[xx,yy] = meshgrid(x,x);

[theta,rho] = cart2pol(xx,yy);
theta = theta(:);
rho = rho(:);

% Vary parameters
sigmasQ1 = (sigma1*sigma2) ./ ((sigma1*sin(theta)).^2 + (sigma2*cos(theta)).^2).^.5;
sigmasQ2 = (sigma2*sigma3) ./ ((sigma2*cos(theta)).^2 + (sigma3*sin(theta)).^2).^.5;
sigmasQ3 = (sigma3*sigma4) ./ ((sigma3*sin(theta)).^2 + (sigma4*cos(theta)).^2).^.5;
sigmasQ4 = (sigma4*sigma1) ./ ((sigma4*cos(theta)).^2 + (sigma1*sin(theta)).^2).^.5;

topfrsQ1 = (topfr1*topfr2) ./ ((topfr1*sin(theta)).^2 + (topfr2*cos(theta)).^2).^.5;
topfrsQ2 = (topfr2*topfr3) ./ ((topfr2*cos(theta)).^2 + (topfr3*sin(theta)).^2).^.5;
topfrsQ3 = (topfr3*topfr4) ./ ((topfr3*sin(theta)).^2 + (topfr4*cos(theta)).^2).^.5;
topfrsQ4 = (topfr4*topfr1) ./ ((topfr4*cos(theta)).^2 + (topfr1*sin(theta)).^2).^.5;

expsQ1 = (exp1*exp2) ./ ((exp1*sin(theta)).^2 + (exp2*cos(theta)).^2).^.5;
expsQ2 = (exp2*exp3) ./ ((exp2*cos(theta)).^2 + (exp3*sin(theta)).^2).^.5;
expsQ3 = (exp3*exp4) ./ ((exp3*sin(theta)).^2 + (exp4*cos(theta)).^2).^.5;
expsQ4 = (exp4*exp1) ./ ((exp4*cos(theta)).^2 + (exp1*sin(theta)).^2).^.5;

fitPts2 = nan(numel(xx),1);

for n = 1:numel(xx)
    
    if theta(n) > pi/2 
        params2 = [topfrsQ1(n), sigmasQ1(n), expsQ1(n), baseline];
    elseif theta(n) > 0 && theta(n) <= pi/2
        params2 = [topfrsQ2(n), sigmasQ2(n), expsQ2(n), baseline];
    elseif theta(n) <= 0 && theta(n) > -pi/2
        params2 = [topfrsQ3(n), sigmasQ3(n), expsQ3(n), baseline];
    elseif theta(n) <= -pi/2
        params2 = [topfrsQ4(n), sigmasQ4(n), expsQ4(n), baseline];
    end
    
    fitPts2(n) = ComputeNakaRushtonJPW(params2,rho(n),'symmetric');
end

fitPts2 = reshape(fitPts2,size(xx));

figure(2); clf; hold on; grid on;
surface(xx,yy,fitPts2)

figure(3); clf; hold on; grid on;
contour(xx,yy,fitPts2)
