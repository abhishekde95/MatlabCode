% Author - Abhishek De, 9/18
% Control experiment to test if the fitting routine quadfit_AD2 can fit a
% narrower than linear curve or not
clearvars; close all;
plot_counter = 1;
x = logspace(-1,1,21)';
% y = 12-x; % creating a line
% y = x.^2; % parabola
% y = 1./x; % creating a rectangular hyperbola
y = sqrt(120-x.^2); % creating a broader than linear arc
not_oog_idx = 1:1:numel(x);
outofgamut = logical(zeros(size(x)));
allthetas = linspace(0,pi/2,100);
% Converting x and y into polar coordinates
[THETA,RHO] = cart2pol(x,y);

% Fitting a line
initguess1 = [x(not_oog_idx) y(not_oog_idx)]\ones(numel(x(not_oog_idx)),1);
[final_model1] = tmp_linefit(RHO, THETA,not_oog_idx,outofgamut,initguess1');
rho1 = 1./(final_model1*[cos(allthetas); sin(allthetas)]);
LOOGtmp1= rho1<0;
[x_lin,y_lin] = pol2cart(allthetas(~LOOGtmp1),rho1(~LOOGtmp1));

% Fitting a conic section using quadfit_AD2
initguess3 = [0 0 0 final_model1];
[final_model3] = tmp_quadfit(RHO, THETA, not_oog_idx,outofgamut,initguess3);
[x_quad,y_quad,rho3] = tmp_calc_xyvalues(allthetas, final_model3);

% Fitting a conic section using conicsectionfit
[model_conic,fval] = tmp_conicsectionfit(x,y,final_model3);
[x_conic,y_conic,rho_conic] = tmp_calc_xyvalues(THETA', model_conic);
[model_conicquad] = tmp_quadfit(RHO, THETA, not_oog_idx,outofgamut,model_conic);
[x_cquad,y_cquad,rho_cquad] = tmp_calc_xyvalues(THETA', model_conicquad);

figure(plot_counter); subplot(421); plot(x,y,'o','MarkerFacecolor',[0 0 1]); hold on; plot(x_lin,y_lin,'k','Linewidth',2); xlabel('x'); ylabel('y'); hold off;
subplot(422); plot(THETA*180/pi,log10(RHO),'o','MarkerFacecolor',[0 0 1]); hold on; plot(allthetas*180/pi,log10(rho1),'k','Linewidth',2); set(gca,'Xlim',[-45 135]); xlabel('theta'); ylabel('logR'); title('Linear fit');hold off;
subplot(423); plot(x,y,'o','MarkerFacecolor',[0 0 1]); hold on; plot(x_quad,y_quad,'g','Linewidth',2); xlabel('x'); ylabel('y'); hold off;
subplot(424); plot(THETA*180/pi,log10(RHO),'o','MarkerFacecolor',[0 0 1]); hold on; plot(allthetas*180/pi,log10(rho3),'k','Linewidth',2); set(gca,'Xlim',[-45 135]); xlabel('theta'); ylabel('logR'); title('Quad fit');hold off;
subplot(425); plot(x,y,'o','MarkerFacecolor',[0 0 1]); hold on; plot(x_conic,y_conic,'r','Linewidth',2); xlabel('x'); ylabel('y'); hold off;
subplot(426); plot(THETA*180/pi,log10(RHO),'o','MarkerFacecolor',[0 0 1]); hold on; plot(THETA*180/pi,log10(rho_conic),'r','Linewidth',2); set(gca,'Xlim',[-45 135]); xlabel('theta'); ylabel('logR'); title('Conic section fit');hold off;
subplot(427); plot(x,y,'o','MarkerFacecolor',[0 0 1]); hold on; plot(x_cquad,y_cquad,'g','Linewidth',2); xlabel('x'); ylabel('y'); hold off;
subplot(428); plot(THETA*180/pi,log10(RHO),'o','MarkerFacecolor',[0 0 1]); hold on; plot(THETA*180/pi,log10(rho_cquad),'g','Linewidth',2); set(gca,'Xlim',[-45 135]); xlabel('theta'); ylabel('logR'); title('Conic-quad section fit');hold off;
plot_counter = plot_counter + 1;
%% Having tried and tested some simple conic sections, I want to test how well this fits the actual raw data
close all; clearvars;
plot_counter = 1;
load RHO_all.mat
load THETA_all.mat
load not_oog_idx_all.mat
load oog_idx_all.mat
load S1LMS.mat
load S2LMS.mat
[OCidx, LMidx, LUMidx, SOidx, hardtoclassifyidx] = classifycells(S1LMS,S2LMS);
idx = [OCidx LMidx LUMidx SOidx hardtoclassifyidx];
tmp_linmodelparams = []; fval_linmodel = [];
tmp_quadmodelparams = []; fval_quadmodel = [];
tmp_conicquadmodelparams = []; fval_conicquadmodel = [];
runsp = [];
for ii = 1:41
    RHO = RHO_all{1,idx(ii)};
    THETA = THETA_all{1,idx(ii)}*pi/180;
    [~,a] = sort(THETA); 
    not_oog_idx = not_oog_idx_all{1,idx(ii)};
    oog_idx = oog_idx_all{1,idx(ii)};
    outofgamut = zeros(size(THETA));
    outofgamut(oog_idx) = 1;
    outofgamut = logical(outofgamut); 
    allthetas = linspace(-pi/4,3*pi/4,100);
    [x,y] = pol2cart(THETA,RHO);
    % Fitting a line
    initguess1 = [x(not_oog_idx) y(not_oog_idx)]\ones(numel(x(not_oog_idx)),1);
    [final_model1,fval1] = tmp_linefit(RHO,THETA,not_oog_idx,outofgamut,initguess1');
    rho1 = 1./(final_model1*[cos(allthetas); sin(allthetas)]);
    LOOGtmp1= rho1<0;
    [x_lin,y_lin] = pol2cart(allthetas(~LOOGtmp1),rho1(~LOOGtmp1));
    rho1idx = imag(rho1)==0 & rho1>0;
    [~,p1] = tmp_calclinSSE(final_model1,RHO,THETA,not_oog_idx);
    runsp = [runsp; p1];
    
    % Fitting a conic section using quadfit_AD2
    initguess3 = [0 0 0 final_model1];
    [final_model3,fval3] = tmp_quadfit(RHO,THETA,not_oog_idx,outofgamut,initguess3);
    [x_quad,y_quad,rho3] = tmp_calc_xyvalues(allthetas, final_model3);
    equation1 = strcat(num2str(final_model3(1)),'*x.^2+',num2str(final_model3(2)),'*y.^2+',num2str(final_model3(3)),'*x.*y+',num2str(final_model3(4)),'*x+',num2str(final_model3(5)),'*y-1');
    rho3idx = imag(rho3)==0 & rho3>0;
%     fval3 = tmp_calcquadSSE(final_model3,RHO,THETA,not_oog_idx);
    
    % Fitting a conic section using conicsectionfit
    [model_conic,fval_conic] = tmp_conicsectionfit(x(not_oog_idx),y(not_oog_idx),final_model3);
    [x_conic,y_conic,rho_conic] = tmp_calc_xyvalues(allthetas, model_conic);
    rho_conicidx = imag(rho_conic)==0 & rho_conic>0;
    equation2 = strcat(num2str(model_conic(1)),'*x.^2+',num2str(model_conic(2)),'*y.^2+',num2str(model_conic(3)),'*x.*y+',num2str(model_conic(4)),'*x+',num2str(model_conic(5)),'*y-1');
    % Re-fitting the conic section using the estimates from conicsectionfit
    [model_conicquad,fval_conicquad] = tmp_quadfit(RHO, THETA, not_oog_idx,outofgamut,model_conic);
    [x_cquad,y_cquad,rho_cquad] = tmp_calc_xyvalues(allthetas, model_conicquad);
    rho_cquadidx = imag(rho_cquad)==0 & rho_cquad>0;
    equation3 = strcat(num2str(model_conicquad(1)),'*x.^2+',num2str(model_conicquad(2)),'*y.^2+',num2str(model_conicquad(3)),'*x.*y+',num2str(model_conicquad(4)),'*x+',num2str(model_conicquad(5)),'*y-1');
%     [~,p2] = tmp_calcquadSSE(model_conicquad,RHO,THETA,not_oog_idx);
    
    % Just plot the non-linear cells 
    figure(plot_counter); subplot(321); plot(x(not_oog_idx),y(not_oog_idx),'o','MarkerFacecolor',[0 0 1]); hold on; plot(x(outofgamut),y(outofgamut),'go','MarkerFacecolor',[0 1 0]); plot(x_lin,y_lin,'k','Linewidth',2); xlabel('x'); ylabel('y'); hold off;
    subplot(322); plot(THETA(not_oog_idx)*180/pi,log10(RHO(not_oog_idx)),'o','MarkerFacecolor',[0 0 1]); hold on; plot(THETA(outofgamut)*180/pi,log10(RHO(outofgamut)),'o','MarkerFacecolor',[0 1 0]); plot(allthetas(rho1idx)*180/pi,log10(rho1(rho1idx)),'ko','Linewidth',2); set(gca,'Xlim',[-45 135]); xlabel('theta'); ylabel('logR'); title('Linear fit');hold off;
    subplot(323); plot(x(not_oog_idx),y(not_oog_idx),'o','MarkerFacecolor',[0 0 1]); hold on; plot(x(outofgamut),y(outofgamut),'go','MarkerFacecolor',[0 1 0]); plot(x_quad,y_quad,'Linewidth',2);  xlabel('x'); ylabel('y'); hold off;
    subplot(324); plot(THETA(not_oog_idx)*180/pi,log10(RHO(not_oog_idx)),'o','MarkerFacecolor',[0 0 1]); hold on; plot(THETA(outofgamut)*180/pi,log10(RHO(outofgamut)),'o','MarkerFacecolor',[0 1 0]); plot(allthetas(rho3idx)*180/pi,log10(rho3(rho3idx)),'ko','Linewidth',2); set(gca,'Xlim',[-45 135]); xlabel('theta'); ylabel('logR'); title('Quad fit');hold off;
    subplot(325); plot(x(not_oog_idx),y(not_oog_idx),'o','MarkerFacecolor',[0 0 1]); hold on; plot(x(outofgamut),y(outofgamut),'go','MarkerFacecolor',[0 1 0]); ezplot(equation3); xlabel('x'); ylabel('y'); hold off;
    subplot(326); plot(THETA(not_oog_idx)*180/pi,log10(RHO(not_oog_idx)),'o','MarkerFacecolor',[0 0 1]); hold on; plot(THETA(outofgamut)*180/pi,log10(RHO(outofgamut)),'o','MarkerFacecolor',[0 1 0]); plot(allthetas(rho_cquadidx)*180/pi,log10(rho_cquad(rho_cquadidx)),'go','Linewidth',2); set(gca,'Xlim',[-45 135]); xlabel('theta'); ylabel('logR'); title('Conic-quad section fit');hold off;
    plot_counter = plot_counter + 1;

    tmp_linmodelparams = [tmp_linmodelparams; final_model1]; fval_linmodel = [fval_linmodel; fval1];
    tmp_quadmodelparams = [tmp_quadmodelparams; final_model3]; fval_quadmodel = [fval_quadmodel; fval3];
    tmp_conicquadmodelparams = [tmp_conicquadmodelparams; model_conicquad]; fval_conicquadmodel = [fval_conicquadmodel; fval_conicquad];
end

%%
% I am trying to come up with random conics and observing their solution in
% rho -theta space, not a very helpful exercise
newtheta = linspace(-pi/4,3*pi/4,51);
figure(plot_counter);
for jj = 1:100
    [X,Y,R] = tmp_calc_xyvalues(newtheta,randn(1,5));
    plot(newtheta*180/pi,log10(R),'Linewidth',2); hold on;
end
set(gca,'Xlim',[-45 135]); xlabel('theta'); ylabel('logR'); title('Linear fit');hold off;
plot_counter = plot_counter + 1;





