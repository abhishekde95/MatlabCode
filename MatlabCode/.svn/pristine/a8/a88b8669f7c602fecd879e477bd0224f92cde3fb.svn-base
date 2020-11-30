% Creating a new script for fitting a spline in theta-log(rho) plane and
% converting the fit into the cartesian plane
% Data obtained from M030117001.nex THETA, RHO
% Data obtained from M122816003.nex THETA1, RHO1
% Data obtained from M030717003.nex THETA2, RHO2
close all; clearvars;
load THETA.mat
load RHO.mat 
load oog_idx.mat

% Fitting 
fitmodel1 = fit(THETA,log(RHO),'cubicinterp'); % Not the best strategy as this method just interpolates the points within
fitmodel2 = fit(THETA,log(RHO),'poly2');
fitmodel3 = fit(THETA,log(RHO),'smoothingspline');
fitmodel4 = fit(THETA,log(RHO),'fourier1');
fitmodel5 = fit(THETA,log(RHO),'fourier2');
outofgamut = zeros(size(THETA));
outofgamut(oog_idx) = 1;
outofgamut = logical(outofgamut);

% % Plotting 
% figure(1),subplot(121);plot(THETA,log(RHO),'bo'); hold on;
% plot(fitmodel1,'r'); plot(fitmodel2,'g'); plot(fitmodel3,'k');
% xlabel('THETA'), ylabel('log(RHO)'); legend('data','cubic','poly2','smoothingspline');
% hold off;

% Trying the cubic and smoothing spline interpolation with the coefficients
% cubic
breaks_cubic = fitmodel1.p.breaks;
coefs_cubic = fitmodel1.p.coefs;
num_cubic = coefs_cubic(:,end);
tmp_x = breaks_cubic(end)-breaks_cubic(end-1);
num_cubic = [num_cubic; coefs_cubic(end,1)*tmp_x.^3 + coefs_cubic(end,2)*tmp_x.^2 + coefs_cubic(end,3)*tmp_x + coefs_cubic(end,4)];

% smooth - spline
breaks_sspline = fitmodel3.p.breaks;
coefs_sspline = fitmodel3.p.coefs;
num_sspline = coefs_sspline(:,end);
tmp_x = breaks_sspline(end)-breaks_sspline(end-1);
num_sspline = [num_sspline; coefs_sspline(end,1)*tmp_x.^3 + coefs_sspline(end,2)*tmp_x.^2 + coefs_sspline(end,3)*tmp_x + coefs_sspline(end,4)];

% 2nd order polynomial
num_poly2 = fitmodel2.p1*(THETA.*THETA) + fitmodel2.p2*THETA + fitmodel2.p3;

% 1st order fourier 
num_fourier1 = fitmodel4.a0 + fitmodel4.a1*cos(THETA.*fitmodel4.w) + fitmodel4.b1*sin(THETA.*fitmodel4.w);  

% 4th order fourier
num_fourier2 = fitmodel5.a0 + fitmodel5.a1*cos(THETA.*fitmodel5.w) + fitmodel5.b1*sin(THETA.*fitmodel5.w) + ...
     fitmodel5.a2*cos(2*(THETA.*fitmodel5.w)) + fitmodel5.b2*sin(2*(THETA.*fitmodel5.w));
     
[~,I] = sort(THETA);
figure(1),subplot(231),plot(THETA(I),log(RHO(I)),'o','MarkerFaceColor',[0 0 1]); hold on; plot(THETA(I),num_cubic,'r','Linewidth',2);
plot(THETA(I),num_sspline,'k','Linewidth',2); plot(THETA(I),num_poly2(I),'g','Linewidth',2); plot(THETA(I),num_fourier1(I),'m','Linewidth',2);
plot(THETA(I),num_fourier2(I),'b','Linewidth',2); legend('Data','cubic','smooth','poly2','fourier1','fourier2'); xlabel('THETA'), ylabel('log(RHO)'); title('R-theta');hold off;

% Reconstructing the points in cartesian coordinates
% Trying out for the second order polynomial model
R_poly2 = exp(num_poly2); x_poly2 = R_poly2(I).*cos(THETA(I)*pi/180); y_poly2 = R_poly2(I).*sin(THETA(I)*pi/180);
R_sspline = exp(num_sspline); x_sspline = R_sspline.*cos(THETA(I)*pi/180); y_sspline = R_sspline.*sin(THETA(I)*pi/180);
R_cubic = exp(num_cubic); x_cubic = R_cubic.*cos(THETA(I)*pi/180); y_cubic = R_cubic.*sin(THETA(I)*pi/180);
R_fourier1 = exp(num_fourier1); x_fourier1 = R_fourier1(I).*cos(THETA(I)*pi/180); y_fourier1 = R_fourier1(I).*sin(THETA(I)*pi/180);
R_fourier2 = exp(num_fourier2); x_fourier2 = R_fourier2(I).*cos(THETA(I)*pi/180); y_fourier2 = R_fourier2(I).*sin(THETA(I)*pi/180);
x_orig = RHO.*cos(THETA*pi/180); y_orig = RHO.*sin(THETA*pi/180);
figure(1),subplot(232);plot(x_orig,y_orig,'o','MarkerFaceColor',[0 0 1]); hold on; plot(x_poly2,y_poly2,'g','Linewidth',2); 
plot(x_cubic,y_cubic,'r','Linewidth',2); plot(x_sspline,y_sspline,'k','Linewidth',2); plot(x_fourier1,y_fourier1,'m','Linewidth',2); 
plot(x_fourier2,y_fourier2,'b','Linewidth',2); xlabel('X'); ylabel('Y'); title('Cartesian');hold off;

% For just visualising the fourier2 fit
figure(1), subplot(233),plot(x_orig(~outofgamut),y_orig(~outofgamut),'o','MarkerFaceColor',[0 0 1]); hold on;
plot(x_orig(outofgamut),y_orig(outofgamut),'o','MarkerFaceColor',[0 1 0]);
plot([x_fourier2;x_fourier2(1)],[y_fourier2;y_fourier2(1)],'b','Linewidth',2); xlabel('X'), ylabel('Y'); title('Fourier2 fit'); hold off;

% Fitting a line a pair of parallel lines to the STAvsPC Neurothresh Data
x = [x_orig(~outofgamut) ones(size(x_orig(~outofgamut)))];
params = x\y_orig(~outofgamut); % obtained by linear regression in log(R)-theta plane
initial_guess = [-params(1);-1;-params(2)];
[final_model,fval] = linefit_AD(RHO, THETA, outofgamut,initial_guess); % Model paramas and error values from the line fit using fminsearch
x_line = linspace(min(x_orig),max(x_orig),101);
Slope = -1*final_model(1)/final_model(2);
Intercept = -1*final_model(3)/final_model(2);
y_line1 = Slope*x_line + Intercept;
y_line2 = Slope*x_line - Intercept;

% Fitting an ellipse to the data, currently not working that well
[final_modele,fvale] = ellipsefit_AD(RHO, THETA, outofgamut); % Model params and error values from the ellipse fit using fminsearch
A = cos(THETA*pi/180);
B = sin(THETA*pi/180);
r_ellipse2 = 1./(((A+final_modele(3)*B)/final_modele(1)).^2 + ((final_modele(3)*A-B)/final_modele(2)).^2);
r_ellipse = sqrt(r_ellipse2);
fact = 1/(1+final_modele(3).^2);
x_ellipse = fact*(final_modele(1)*A(I)+final_modele(3)*final_modele(2)*B(I));
y_ellipse = fact*(final_modele(3)*final_modele(1)*A(I)-final_modele(2)*B(I));

figure(1),subplot(234), plot(x_orig(~outofgamut),y_orig(~outofgamut),'o','MarkerFaceColor',[0 0 1]); hold on;
plot(x_orig(outofgamut),y_orig(outofgamut),'o','MarkerFaceColor',[0 1 0]);
plot(x_line,y_line1,'k','Linewidth',2); plot(x_line,y_line2,'k','Linewidth',2);
plot([x_ellipse; x_ellipse(1)],[y_ellipse; y_ellipse(1)],'g','Linewidth',2);
xlabel('X'), ylabel('Y'); title('Line & Ellipse Fit'); hold off;

% Fitting in the cartesian plane
figure(1), subplot(235),plot(x_orig(~outofgamut),y_orig(~outofgamut),'o','MarkerFaceColor',[0 0 1]); hold on;
plot(x_orig(outofgamut),y_orig(outofgamut),'o','MarkerFaceColor',[0 1 0]);
plot(x_line,params(1)*x_line+params(2),'r','Linewidth',2);
plot(x_line,params(1)*x_line-params(2),'r','Linewidth',2);
xlabel('X'), ylabel('Y'),title('Regression');hold off;

% Plotting the the line fit in log(R)-Theta domain

pred_RHO = -1*final_model(3)./(final_model(2)*sin(THETA*pi/180)+final_model(1)*cos(THETA*pi/180));
pred_regression = initial_guess(2)./(sin(THETA*pi/180)-initial_guess(1)*cos(THETA*pi/180));
figure(1), subplot(236), plot(THETA(I),log(abs(pred_RHO(I))),'k','Linewidth',2), hold on;
% plot(THETA(I),log(abs(pred_regression(I))),'r','Linewidth',2),
plot(THETA(I),log(r_ellipse(I)),'g','Linewidth',2), % Ellipse fit 
plot(THETA(~outofgamut),log(RHO(~outofgamut)),'o','MarkerFaceColor',[0 0 1]);
plot(THETA(outofgamut),log(RHO(outofgamut)),'o','MarkerFaceColor',[0 1 0]);xlabel('THETA'); ylabel('log(RHO');
title('Fit & Regression');hold off;

%% Lets try out grid search in the space spanned by Intercept and Slope
Intercept  = linspace(-2,2,51);
Slope  = linspace(-50,50,101);
error = zeros(numel(Intercept),numel(Slope));
for ii = 1:numel(Intercept)
    for jj = 1:numel(Slope)
        pred = (sin(THETA*pi/180)-(Slope(jj)*cos(THETA*pi/180)))/Intercept(ii);
        resid = log(abs(pred))-log(1./RHO);
        resid(outofgamut & (abs(pred)>1./RHO))=0;
        error(ii,jj) = sum(resid.^2);
    end
end
figure(2),surfl(error); xlabel('Intercept'); ylabel('Slope'); zlabel('Error');
title('Error w/ grid search');


