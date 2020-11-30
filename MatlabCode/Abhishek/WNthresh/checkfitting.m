% Author - Abhishek De, 9/18
% Control experiment to test if the fitting routine quadfit_AD2 can fit a
% narrower than linear curve or not
clearvars; close all;
plot_counter = 1;
x = logspace(-1,1,51)';
y = 1./x; % creating a rectangular hyperbola
y = x.^2; % creating a broader than linear arc
not_oog_idx = 1:1:numel(x);
outofgamut = logical(zeros(size(x)));
allthetas = linspace(0,pi/2,100);
% Converting x and y into polar coordinates 
[THETA,RHO] = cart2pol(x,y);
THETA = THETA*180/pi;

% Fitting a line
initguess1 = [x(not_oog_idx) y(not_oog_idx)]\ones(numel(x(not_oog_idx)),1);
[final_model1,fval_linR,fval_linOLS] = linefit_AD2(RHO, THETA,not_oog_idx,outofgamut,initguess1');
rho1 = 1./(final_model1*[cos(allthetas); sin(allthetas)]);
LOOGtmp1= rho1<0;
[x_lin,y_lin] = pol2cart(allthetas(~LOOGtmp1),rho1(~LOOGtmp1));

% Fitting a conic section
initguess3 = [0 0 0 final_model1];
%     initguess3 = [x(not_oog_idx).^2 y(not_oog_idx).^2 x(not_oog_idx).*y(not_oog_idx) x(not_oog_idx) y(not_oog_idx)]\ones(numel(x(not_oog_idx)),1); initguess3 = initguess3';
[final_model3,fval_quadR,fval_quadOLS] = quadfit_AD2(RHO, THETA, not_oog_idx,outofgamut,initguess3);
A = final_model3(1); B = final_model3(2); C = final_model3(3); D = final_model3(4); E = final_model3(5);
rho3 = [];
p = [A*cos(allthetas').^2+B*sin(allthetas').^2+C*(cos(allthetas').*sin(allthetas')) D*cos(allthetas')+E*sin(allthetas') -1*ones(numel(allthetas'),1)];
for kk = 1:size(p,1)
    rts = roots(p(kk,:));
    disp(rts);
%     if all(rts>0)
%         r = min(rts); % if both the roots are positive
%     else
        r = max(rts);
%     end
    rho3 = [rho3; r];
end
L = rho3>0 & rho3==real(rho3);
[x_quad,y_quad] = pol2cart(allthetas(L),rho3(L)');

figure(plot_counter); subplot(121); plot(x,y,'o','MarkerFacecolor',[0 0 1]); hold on; plot(x_lin,y_lin,'k','Linewidth',2); plot(x_quad,y_quad,'g','Linewidth',2);xlabel('x'); ylabel('y'); hold off; 
subplot(122); plot(THETA,log10(RHO),'o','MarkerFacecolor',[0 0 1]); hold on; plot(allthetas*180/pi,log10(rho1),'k','Linewidth',2); plot(allthetas*180/pi,log10(rho3),'g','Linewidth',2); set(gca,'Xlim',[-45 135]); xlabel('theta'); ylabel('logR'); title('Linear fit');hold off;

% As of now the fitting routine is doing a good job of fitting a broader 
% than linear contour but not narrower than linear contour

