% Am doing some pruning to the fit_rho_theta.m
% Removing some of the unnecessary stuff
% Author - Abhishek De, 3/17
% Data obtained from M030117001.nex THETA, RHO
% Data obtained from M122816003.nex THETA1, RHO1
% Data obtained from M030717003.nex THETA2, RHO2
% Data obtained from M010317003.nex THETA3, RHO3
close all; clearvars;
load THETA.mat
load RHO.mat 
load oog_idx.mat

outofgamut = zeros(size(THETA));
outofgamut(oog_idx) = 1;
outofgamut = logical(outofgamut);
[~,I] = sort(THETA);

x_orig = RHO.*cos(THETA*pi/180); y_orig = RHO.*sin(THETA*pi/180);
% Fitting a pair of parallel lines to the STAvsPC Neurothresh Data
x = [x_orig(~outofgamut) ones(size(x_orig(~outofgamut)))];
params = x\y_orig(~outofgamut); % obtained by linear regression in log(R)-theta plane
initial_guess = [-RHO(find(THETA==90));-1;-params(2)];
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

figure(1),subplot(141), plot(x_orig(~outofgamut),y_orig(~outofgamut),'o','MarkerFaceColor',[0 0 1]); hold on;
plot(x_orig(outofgamut),y_orig(outofgamut),'o','MarkerFaceColor',[0 1 0]);
plot(x_line,y_line1,'k','Linewidth',2); plot(x_line,y_line2,'k','Linewidth',2);
plot([x_ellipse; x_ellipse(1)],[y_ellipse; y_ellipse(1)],'g','Linewidth',2);
xlabel('X'), ylabel('Y'); title('Line & Ellipse Fit'); grid on,
set(gca,'XLim',[-2 2],'Ylim',[-2 2]);hold off;

% Fitting in the cartesian plane
figure(1), subplot(142),plot(x_orig(~outofgamut),y_orig(~outofgamut),'o','MarkerFaceColor',[0 0 1]); hold on;
plot(x_orig(outofgamut),y_orig(outofgamut),'o','MarkerFaceColor',[0 1 0]);
plot(x_line,params(1)*x_line+params(2),'k','Linewidth',2);
plot(x_line,params(1)*x_line-params(2),'k','Linewidth',2);
xlabel('X'), ylabel('Y'),title('Regression'); grid on; hold off;

% Plotting the the line fit in log(R)-Theta domain
pred_RHO = -1*final_model(3)./(final_model(2)*sin(THETA*pi/180)+final_model(1)*cos(THETA*pi/180));
pred_regression = initial_guess(2)./(sin(THETA*pi/180)-initial_guess(1)*cos(THETA*pi/180));
figure(1), subplot(143), plot(THETA(I),log(abs(pred_RHO(I))),'k','Linewidth',2), hold on;
plot(THETA(I),log(r_ellipse(I)),'g','Linewidth',2), % Ellipse fit 
plot(THETA(~outofgamut),log(RHO(~outofgamut)),'o','MarkerFaceColor',[0 0 1]);
plot(THETA(outofgamut),log(RHO(outofgamut)),'o','MarkerFaceColor',[0 1 0]);xlabel('THETA'); ylabel('log(RHO');
title('Fit & Regression');hold off;

strcat('Line fit= ',num2str(fval))
strcat('Ellipse fit= ',num2str(fvale))

% Next experiment will be to just fit one line instead of a pair of parallel lines 
a  = linspace(-20,20,100);
b  = linspace(-20,20,100);
c  = linspace(-20,20,100);
error = zeros(numel(a),numel(b),numel(c));
valid_indices = ones(numel(a),numel(b),numel(c));
ref_error = 1000000000000;
for ii = 1:numel(a)
    for jj = 1:numel(b)
        for kk = 1:numel(c)
            pred = -1*(b(jj)*sin(THETA*pi/180)+a(ii)*cos(THETA*pi/180))/c(kk);
            resid = log(abs(pred))-log(1./RHO);
            resid(outofgamut & (abs(pred)<1./RHO))=0;
            error(ii,jj) = sum(resid.^2);
            if any(pred<0)
                valid_indices(ii,jj,kk) = 0;
            else
                if error(ii,jj,kk) < ref_error
                    ref_error = error(ii,jj,kk);
                    idx1 = ii;
                    idx2 = jj;
                    idx3 = kk;
                end
            end
        end
    end
end

if exist('idx1') % Only if a solution for a single line fitting exists
    Slope = -1*a(idx1)/b(idx2);
    Intercept = -1*c(idx3)/b(idx2);
    y_line = Slope*x_line + Intercept;
    figure(1),subplot(144),plot(x_line,y_line,'r','Linewidth',2); hold on;
    plot(x_orig(~outofgamut),y_orig(~outofgamut),'o','MarkerFaceColor',[0 0 1]);
    plot(x_orig(outofgamut),y_orig(outofgamut),'o','MarkerFaceColor',[0 1 0]);
    xlabel('X'); ylabel('Y'); title('Single Line fit'); grid on;  hold off;
    
    r_singleline = (-1*c(idx3))./(a(idx1)*cos(THETA*pi/180) + b(idx2)*sin(THETA*pi/180));
    figure(1),subplot(143),hold on; plot(THETA(I),log(r_singleline(I)),'r','Linewidth',2); % Single line fit
    set(gca,'YLim',[-4 3]); grid on; hold off;
end



